"""rabbit ParamModel for a fully differentiable SCETlib prediction.

The prediction is assembled in three steps::

    1  SCETlib cached rule replay  ->  sigma(p; g)       boson level, gen grid
    2  fold through the response R ->  sigma_reco(p; b)  gen -> reco
    3  ratio to the reference      ->  rnorm(b, proc)    handed to rabbit

``p`` is SCETlib's own differentiation vector, so every theory parameter the
calculation exposes is a continuous fit parameter with exact derivatives:

* ``alphaS``                     -- the strong coupling;
* the nonperturbative lambdas    -- the Collins-Soper and TMD form factors;
* the theory nuisance parameters -- ``gamma_cusp``, ``gamma_mu_q``,
  ``gamma_nu``, ``s``, ``h_qqV`` and the five beam-function TNPs, present
  whenever the runcard declares a ``[TNPs]`` block (which is also what makes the
  prediction N^{3+0}LL rather than N3LL);
* PDF eigenvector coefficients   -- carried as additional differentiable columns
  by a cache built with PDF variations, exact at ``c_e = 0, +-1``. Supplying the
  PDF set's alpha_s member pair alongside them folds the
  ``dsigma/dPDF . dPDF/dalpha_s`` piece into the ``alphaS`` slot, so ``alphaS``
  becomes the PDF-consistent coupling rather than alpha_s at fixed PDF.

Which of these are present is a property of the cache, not of this file: the
model reads ``gradient_param_names()`` and registers what it finds. Only the
profile-scale parameters (``kappaFO``, ``kappaf``, ``muf``, the transition
points) are outside SCETlib's autodiff and still need template nuisances.

Why the derivatives are injected rather than traced
---------------------------------------------------
SCETlib returns the value, the full Jacobian and the exact per-bin Hessian at the
current point, all from C++. Rather than differentiating THROUGH the py_function
boundary, we evaluate once and hand autodiff an exact local quadratic::

    d          = p - stop_gradient(p)                 # value 0, unit gradient
    sigma_gen  = stop_gradient(val) + J.d + 0.5 d^T K d

At the evaluation point the value is exact, the first derivative is exactly J and
the second exactly K, while everything downstream (fold, ratio, floor) stays pure
TF and is differentiated analytically. Two practical consequences:

* the ``PyFunc`` op sits behind ``stop_gradient``, so ``GradientTape.jacobian``
  (whose pfor vectorisation has no converter for a ``PyFunc``) never re-enters
  C++ once per fit parameter;
* J and K come from the C++ call for free, with no ``ForwardAccumulator`` loop.

``differentiate=through`` selects the alternative -- let TF drive the C++
callbacks through the bridge's nested ``custom_gradient`` wrappers -- so this is a
checkable claim rather than an assertion. Both paths run a full fit to the same
answer, and agree to 1.3e-16 on the gradient, 4e-15 on HVPs and 2.4e-17 on the
Hessian (``scripts/rabbit/scetlib_ad/check_differentiate_modes.py``).

``through`` is the DEFAULT: it is the ordinary TF idiom, needs no reasoning about
what autodiff can see, and is what the SCETlib example does. It is also FASTER at
small parameter counts -- measured on a 6-parameter gen-level card, 15.8 s of fit
time against 40.9 s for straightthrough.

``straightthrough`` is kept because the two scale oppositely in the number of FIT
parameters, counting every datacard nuisance and not just ours. rabbit builds the
postfit Hessian as ``t2.jacobian(grad, self.x)``, and pfor has no converter for a
``PyFunc``, so ``through`` costs one C++ HVP sweep per fit parameter (the bridge
caches values and Jacobians, but not HVPs). ``straightthrough`` instead pays one
value+Jacobian and one full Hessian per distinct parameter point, whatever the
parameter count. Since an HVP is ~2.5x a gradient and a materialised Hessian
~40x, the crossover is a few tens of parameters: ``through`` wins on a
handful, ``straightthrough`` on a card carrying hundreds of nuisances.

One requirement that ``through`` imposes and is easy to break: map rabbit's fit
vector into SCETlib's layout with a CONSTANT 0/1 MATRIX MULTIPLY, never
``tensor_scatter_nd_update``. A scatter's backward pass contains a gather, whose
gradient TF represents as ``tf.IndexedSlices``, and the bridge's second-order
py_function payloads call ``.numpy()`` on the incoming cotangent and fail on it.
The matmul is bit-identical (entries are exactly 0 and 1) and free at these sizes.

K is always included, so the composite Hessian is exact and the fit and the
postfit covariance are ONE rabbit job.

XLA cannot compile a ``PyFunc``, so the fit MUST run with ``--jitCompile off``.
The model checks this at construction.
"""

import copy

import numpy as np
import tensorflow as tf
from rabbit.param_models.param_model import ParamModel

from wremnants.postprocessing.scetlib_ad import params as adp
from wremnants.postprocessing.scetlib_ad.response import (
    DEFAULT_RESPONSE_GROUP,
    RATIO_FLOOR_MIN,
    RATIO_FLOOR_SCALE,
    R_info_from_auxiliary,
    crop_R_to_fit,
    marginalize_R_reco,
    np_anchor_from_meta,
)
from wremnants.postprocessing.scetlib_ad.xsec_backend import ScetlibADXsec

DTYPE = tf.float64

# Substring -> the registered parameter that makes a card syst a double count.
# Checked case-insensitively against indata.systs; only the entries whose model
# parameter is actually fitted are enforced, so a lambda-only run still tolerates
# a card carrying pdfAlphaS.
_CONFLICTS = (
    (
        "scetlibnp",
        "any NP lambda",
        lambda names: any(n.startswith("lambda") for n in names),
    ),
    ("pdfalphas", "alphaS", lambda names: "alphaS" in names),
    (
        "resumtnp",
        "a resummation TNP",
        lambda names: any(n.startswith(adp.TNP_PREFIX_OUT) for n in names),
    ),
)


def _as_name_tuple(value):
    """Spec tokens arrive as strings; accept ``a,b`` as well as a real tuple."""
    if value is None:
        return ()
    if isinstance(value, str):
        return tuple(s.strip() for s in value.split(",") if s.strip())
    return tuple(value)


class SCETlibADParamModel(ParamModel):
    """Fit SCETlib's own differentiable parameters directly.

    Usage::

        --paramModel wremnants.postprocessing.scetlib_ad.SCETlibADParamModel \\
            cache=<cache>.npz conf=<runcard>.conf gen_level=1 [key=value ...]
    """

    @classmethod
    def parse_args(cls, indata, *args, **kwargs):
        """``key=value`` spec tokens, typed off the ``__init__`` default.

        A later duplicate key wins, which is what lets a driver append overrides
        to a spec it inherited from a previous step's recorded arguments.
        """
        import inspect

        sig = inspect.signature(cls.__init__)
        valid = {
            n: p
            for n, p in sig.parameters.items()
            if n not in ("self", "indata")
            and p.kind is not inspect.Parameter.VAR_KEYWORD
        }
        positional = []
        for tok in args:
            key = tok.split("=", 1)[0] if isinstance(tok, str) and "=" in tok else None
            if key is not None and key not in valid:
                # Do NOT fall through to positional: every argument of this model
                # is keyword-with-default, so an unknown key=value is a typo (or a
                # token that has been removed), and silently treating it as
                # positional produces "got multiple values for argument 'cache'".
                raise TypeError(
                    f"{cls.__name__}: unknown spec token {key!r}. Valid tokens: "
                    + ", ".join(sorted(valid))
                )
            if key is not None:
                val = tok.split("=", 1)[1]
                default = valid[key].default
                if isinstance(default, bool):
                    val = str(val).strip().lower() in ("1", "true", "yes", "on")
                elif isinstance(default, float):
                    val = float(val)
                elif isinstance(default, int):
                    val = int(val)
                kwargs[key] = val
            else:
                positional.append(tok)
        return cls(indata, *positional, **kwargs)

    def __init__(
        self,
        indata,
        cache=None,
        conf=None,
        gen_level=False,
        signal_proc="Zmumu",
        Q_lo=60.0,
        Q_hi=120.0,
        fit_params=None,
        poi_params="alphaS",
        differentiate="through",
        threads=0,
        priors=False,
        prior_sigmas=None,
        xparam_default=None,
        check_anchor=True,
        anchor_tol=1e-6,
        response_group=DEFAULT_RESPONSE_GROUP,
        **kwargs,
    ):
        """
        Parameters
        ----------
        cache, conf
            The ``.npz`` written by ``ScetlibCachedXsecTF.save`` and the SCETlib
            runcard it was built from. Both are required: the cache holds the
            compressed rules and the frozen fixed-order grid, and the runcard is
            what rebuilds the identical calculation they attach to.
        gen_level
            Gen-level sigmaUL mode. The fit channel IS the gen (qT, |Y|) binning,
            so there is no response matrix and no fold: ``compute()`` returns the
            per-gen-bin ratio ``sigma_gen(p) / sigma_gen(p_anchor)``. Otherwise
            the reco fold reads R from the datacard's response auxiliary
            (see ``response_group``).
        signal_proc
            Process whose column carries the ratio; the others stay at 1.
        Q_lo, Q_hi
            The mass window the cache's single Q bin must span.
        fit_params
            Comma-separated rabbit-facing names to expose to the fit. Default:
            every parameter the cache carries except ``params.DEFAULT_FROZEN``
            (the tanh saturation scales and the b* convention, which are shape
            constants). Parameters not listed are held at their cache anchor and
            never reach rabbit, so they cannot contribute a zero-derivative
            (singular) Hessian row.
        poi_params
            Subset of ``fit_params`` reported as POIs (they must come first in
            the fitter's layout). Default ``alphaS``.
        differentiate
            How TensorFlow gets derivatives of sigma_gen.

            ``through`` (default) calls ``ScetlibCachedXsecTF`` inside the graph
            and lets TF differentiate it, via the nested ``tf.custom_gradient``
            wrappers that call back into C++. This is the ordinary TF idiom and
            the same one ``examples/matched_ad/tf_gradients.py`` uses.

            ``straightthrough`` instead evaluates SCETlib once per parameter point
            and hands autodiff an exact local quadratic (see the module docstring).

            Both give the same value, gradient and Hessian -- measured agreement
            1.3e-16, 4e-15 and 2.4e-17, and a full fit either way returns the same
            alphaS and the same uncertainties. They differ only in WHO drives the
            C++ calls, and therefore in how the cost scales: see the module
            docstring. Switch to ``straightthrough`` if the postfit Hessian
            becomes the bottleneck on a card with many nuisances.
        threads
            SCETlib worker threads for the batch replay (0 = one per hardware
            thread).
        priors
            Declare Gaussian priors. rabbit applies priors whenever a model
            declares ``prior_sigmas``, so this token IS the decision. Off by
            default -- everything floats free. Note the TNP defaults are sigma=1
            (they are genuine nuisances), unlike the lambdas.
        prior_sigmas
            Per-name override, ``name=value,...`` or a Mapping. ``nan`` frees a
            parameter.
        xparam_default
            ``name=value,...`` shifting the fit START (and the prior mean) off
            the cache anchor, for injection / closure tests. The ratio
            DENOMINATOR is not moved -- it always stays the anchor.
        check_anchor
            Cross-check the cache anchor against the nonperturbative values the
            card records for its own prediction. An anchor that disagrees
            is the silent-wrong-answer trap documented in
            ``knowledge/20_frameworks/gen_level_sigmaul_fit.md``: the ratio is
            still 1 at the start, so nothing looks broken, but the response is
            evaluated at the wrong point.
        """
        self.indata = indata
        if cache is None or conf is None:
            raise ValueError(
                "SCETlibADParamModel needs both cache=<cache>.npz and "
                "conf=<runcard>.conf spec tokens."
            )
        self._require_no_xla(kwargs)

        self._response_group = str(response_group)
        self.gen_level = bool(gen_level)
        if differentiate not in ("straightthrough", "through"):
            raise ValueError(
                f"differentiate must be 'straightthrough' or 'through', got "
                f"{differentiate!r}"
            )
        self._through = differentiate == "through"
        # The exact second-derivative term is always included: it is what makes
        # the composite Hessian exact, and it lets the fit and the postfit
        # covariance be ONE rabbit job. Dropping it (Gauss-Newton) was measured to
        # give identical numbers on Asimov and ~5x the speed, but that is a
        # speed/exactness trade we do not need at the binning we fit on (~1 s/bin
        # of serial Hessian work, scaling as 1 + P(P+1)/2), so it is not offered.
        # Reintroduce a switch here if a full-correction-grid fit ever needs it.
        self._curvature = not self._through

        # ---- Backend: rebuild the calculation and load the cache.
        self.core = ScetlibADXsec(conf, cache, threads=threads)
        self.scetlib_names = list(self.core.param_names)
        self.rabbit_names = [adp.rabbit_name(n) for n in self.scetlib_names]
        self._anchor = np.asarray(self.core.anchor, dtype=np.float64)

        # ---- Gen binning, and (reco path) the response matrix.
        self._setup_binning(indata, Q_lo, Q_hi)

        # ---- Parameter registration. Everything the fit does NOT expose stays
        # pinned at the anchor, so the SCETlib vector is always complete.
        self._register_params(fit_params, poi_params, xparam_default)

        # ---- Central: the ratio denominator, evaluated by the model itself at
        # the anchor so the ratio is exactly 1 at the start whatever the card's
        # own template looks like.
        self._eval_cache = None
        sigma_gen_anchor = self._sigma_gen_np(self._p_base_anchor)
        self.sigma_gen_central_flat = tf.constant(sigma_gen_anchor, dtype=DTYPE)
        if self.gen_level:
            self.sigma_reco_central = None
        else:
            self.sigma_reco_central = tf.linalg.matvec(
                self.R, self.sigma_gen_central_flat
            )
            n_bad = int(tf.reduce_sum(tf.cast(self.sigma_reco_central <= 0, tf.int32)))
            if n_bad:
                raise ValueError(
                    f"SCETlibADParamModel: {n_bad} reco bins have non-positive "
                    f"sigma_reco at the anchor. Likely a binning mismatch between "
                    f"R and the fit-tensor reco axes."
                )

        if check_anchor:
            self._check_anchor_against_card(anchor_tol)
        self._check_double_counting()
        self._check_no_inert_params()

        # ---- Process column.
        procs = [p.decode() if isinstance(p, bytes) else str(p) for p in indata.procs]
        if signal_proc not in procs:
            raise ValueError(
                f"SCETlibADParamModel: signal_proc={signal_proc!r} not in "
                f"indata.procs={procs[:10]}..."
            )
        self.signal_proc_idx = procs.index(signal_proc)
        self.nproc = indata.nproc
        self._signal_col_mask = tf.reshape(
            tf.one_hot(self.signal_proc_idx, self.nproc, dtype=indata.dtype),
            [1, self.nproc],
        )

        self._setup_priors(priors, prior_sigmas)

        print(
            f"[SCETlibADParamModel] {self.core} | {self._fold.describe()} | "
            f"{'gen-level' if self.gen_level else 'reco'} | "
            f"fitting {self.nparams} of {self.core.n_params} "
            f"({self.npoi} POI: {[n for n in self._param_order[: self.npoi]]}) | "
            f"differentiate={'through' if self._through else 'straightthrough'}",
            flush=True,
        )

    # =========================================================================
    # construction helpers
    # =========================================================================

    def _require_no_xla(self, kwargs):
        """Refuse to build under XLA -- a PyFunc has no XLA lowering.

        rabbit resolves ``--jitCompile auto`` to True in dense mode
        (Fitter.__init__), and the failure is an opaque compile error deep in the
        first loss evaluation, so trip here with the fix in the message.
        """
        opt = str(kwargs.get("jitCompile", "auto")).lower()
        sparse = bool(getattr(self.indata, "sparse", False))
        # --eager turns every tf.function into eager execution, so jit_compile
        # never applies and a PyFunc is fine.
        if kwargs.get("eager") or opt == "off" or (opt == "auto" and sparse):
            return
        raise ValueError(
            "SCETlibADParamModel calls into SCETlib through tf.py_function, "
            "which XLA cannot compile. Re-run with --jitCompile off "
            f"(got --jitCompile {opt}" + (", dense input)." if not sparse else ").")
        )

    def _setup_binning(self, indata, Q_lo, Q_hi):
        """Resolve the gen grid (and R, in the reco path) and map it onto the cache."""
        if self.gen_level:
            gen_axes = self._fit_axes(indata)
            if len(gen_axes) != 2:
                raise NotImplementedError(
                    "gen_level SCETlibADParamModel expects a single fit channel "
                    "with 2 gen axes (qT, |Y|); got "
                    f"{[n for n, _ in gen_axes]}"
                )
            self.R = None
            self.reco_shape = None
        else:
            R_info = R_info_from_auxiliary(indata, self._response_group)
            fit_reco_axes = self._fit_axes(indata)
            R_full, R_reco_axes = marginalize_R_reco(
                R_info["R"], R_info["reco_axes"], [n for n, _ in fit_reco_axes]
            )
            R_arr = crop_R_to_fit(R_full, R_reco_axes, fit_reco_axes)
            self.reco_shape = R_arr.shape[: len(fit_reco_axes)]
            gen_axes = R_info["gen_axes"]
            if R_info.get("N_gen") is None:
                raise ValueError(
                    f"SCETlibADParamModel: the {self._response_group!r} "
                    "auxiliary has no N_gen "
                    "(gen-total). Rebuild the datacard from a histmaker output "
                    "that carries the 'prefsr' xnorm hist."
                )
            n_reco = int(np.prod(self.reco_shape))
            n_gen = int(np.prod([len(e) - 1 for _, e in gen_axes]))
            R_raw = tf.constant(R_arr.reshape(n_reco, n_gen), dtype=DTYPE)
            n_gen_flat = tf.constant(
                np.asarray(R_info["N_gen"]).reshape(-1), dtype=DTYPE
            )
            # R must encode only the gen->reco mapping, not the MC's absolute gen
            # spectrum: normalise each gen column by the gen-total (empty gen bins
            # keep a zero column).
            safe = tf.where(n_gen_flat > 0, n_gen_flat, tf.ones_like(n_gen_flat))
            self.R = R_raw / safe[tf.newaxis, :]

        self.gen_axes = [
            (name, np.asarray(edges, dtype=np.float64)) for name, edges in gen_axes
        ]
        self.gen_shape = tuple(len(e) - 1 for _, e in self.gen_axes)
        self.Q_lo, self.Q_hi = float(Q_lo), float(Q_hi)

        # compute() returns one row per fit bin, so the shape it builds must be
        # the card's. A mismatch here would surface as an opaque broadcast error
        # inside the first loss evaluation.
        n_rows = int(np.prod(self.gen_shape if self.gen_level else self.reco_shape))
        if n_rows != int(indata.nbins):
            raise ValueError(
                f"SCETlibADParamModel: the model produces {n_rows} bins "
                f"({'gen' if self.gen_level else 'reco'} shape "
                f"{self.gen_shape if self.gen_level else self.reco_shape}) but the "
                f"card has {int(indata.nbins)}."
            )
        # Exact sum of cache bins onto the gen grid: handles a different nesting
        # order, a signed-Y cache folded onto |Y|, and a cache finer than the fit's
        # gen binning. Coverage is verified, so a cache that does not tile this
        # card's gen bins raises here rather than integrating over less phase space.
        self._fold = self.core.fold_for(self.gen_axes, self.Q_lo, self.Q_hi)

    def _fit_axes(self, indata):
        """(name, edges) of each axis of the single non-masked channel."""
        non_masked = [
            (name, info)
            for name, info in indata.channel_info.items()
            if not info.get("masked", False)
        ]
        if len(non_masked) != 1:
            raise NotImplementedError(
                f"SCETlibADParamModel supports a single non-masked channel; got "
                f"{len(non_masked)}: {[n for n, _ in non_masked]}"
            )
        _, info = non_masked[0]
        return [
            (ax.name, np.asarray(ax.edges, dtype=np.float64)) for ax in info["axes"]
        ]

    def _register_params(self, fit_params, poi_params, xparam_default):
        """Decide which SCETlib parameters rabbit sees, and their start values."""
        available = list(self.rabbit_names)
        # Default: alpha_s and the NP lambdas, minus the shape constants. TNPs are
        # NOT included by default even when the cache carries them -- an
        # analysis-faithful runcard registers all ten (theta=0 'level0' IS the
        # N^{3+0}LL prescription), and silently floating ten unconstrained theory
        # nuisances is not what "fit the NP model" should mean. Ask for them
        # explicitly via fit_params, with priors.
        requested = _as_name_tuple(fit_params) or tuple(
            n
            for n in available
            if n not in adp.DEFAULT_FROZEN and not n.startswith(adp.TNP_PREFIX_OUT)
        )
        unknown = [n for n in requested if n not in available]
        if unknown:
            raise ValueError(
                f"SCETlibADParamModel: fit_params {unknown} are not in this "
                f"cache's parameter set {available}. A parameter can only be "
                f"fitted if the runcard declared it before the rules were built."
            )
        pois = _as_name_tuple(poi_params)
        bad_pois = [n for n in pois if n not in requested]
        if bad_pois:
            raise ValueError(
                f"SCETlibADParamModel: poi_params {bad_pois} are not in "
                f"fit_params {list(requested)}."
            )
        # rabbit's layout contract: all POIs first, then the POUs.
        nou = tuple(n for n in requested if n not in pois)
        self._param_order = tuple(pois) + nou
        self.npoi = len(pois)
        self.npou = len(nou)
        self.params = np.array([p.encode() for p in self._param_order])

        # Position of each fitted parameter inside SCETlib's own vector. rabbit's
        # vector is NOT SCETlib's: it holds only what we fit, POIs first, while
        # SCETlib's has every registered parameter in its registry order.
        self._fit_idx = np.array(
            [available.index(n) for n in self._param_order], dtype=np.int64
        )
        # The map from rabbit's vector to SCETlib's, as a constant 0/1 matrix.
        # Deliberately NOT tensor_scatter_nd_update: a scatter's backward pass
        # contains a gather, whose gradient TF represents as tf.IndexedSlices, and
        # the SCETlib bridge's second-order py_function payloads call .numpy() on
        # the incoming cotangent and die on it. A matmul against a constant gives
        # a dense gradient, so the differentiate=through path survives at second
        # order. Bit-identical to the scatter (the entries are exactly 0 and 1)
        # and negligible in cost: (n_scetlib, n_fit) is at most ~25 x 25.
        self._select = np.zeros((len(available), len(self._param_order)))
        self._select[self._fit_idx, np.arange(len(self._param_order))] = 1.0

        # Start values: the cache anchor, optionally shifted for injection tests.
        # _p_base_anchor is the UNSHIFTED full vector and stays the ratio
        # denominator; _p_base carries the shift for the non-fitted slots only
        # (fitted slots are overwritten from the fit vector on every call).
        self._p_base_anchor = self._anchor.copy()
        self._p_base = self._anchor.copy()
        defaults = self._anchor[self._fit_idx].copy()
        for name, val in _parse_kv(xparam_default).items():
            if name not in available:
                raise KeyError(f"xparam_default: unknown parameter {name!r}")
            if name in self._param_order:
                defaults[self._param_order.index(name)] = val
            else:
                # not fitted: pin the held value at the shifted point
                self._p_base[available.index(name)] = val
                print(
                    f"[SCETlibADParamModel] xparam_default {name}={val:g} applies "
                    f"to a NON-fitted parameter; it is pinned there, not floated.",
                    flush=True,
                )
        if xparam_default:
            print(
                "[SCETlibADParamModel] start shifted: "
                f"{dict(zip(self._param_order, defaults))}",
                flush=True,
            )

        # lambdas can be legitimately zero or negative (delta_lambda2), so store
        # POIs directly rather than as sqrt(value).
        self.allowNegativeParam = True
        self.is_linear = False
        self.xparamdefault = tf.constant(defaults, dtype=self.indata.dtype)

        active = set(self._param_order)
        groups = {
            label: tuple(p for p in members if p in active)
            for label, members in adp.IMPACT_GROUP_MEMBERS.items()
        }
        tnps = adp.tnp_group(self._param_order)
        if tnps:
            groups["resumTNP"] = tnps
        self.param_impact_groups = {k: v for k, v in groups.items() if v}

    def _setup_priors(self, priors, prior_sigmas):
        """Declare ``prior_sigmas`` only when asked; rabbit's Fitter keys off it."""
        tnps = adp.tnp_group(self._param_order)
        if tnps and not priors:
            # theta is normalised upstream so |theta| = 1 IS the recommended
            # variation; floating a TNP free discards that, and the fit will
            # happily absorb a physical effect into an unconstrained nuisance.
            raise ValueError(
                f"SCETlibADParamModel: {len(tnps)} theory nuisance parameter(s) "
                f"are being fitted ({', '.join(tnps[:4])}"
                f"{', ...' if len(tnps) > 4 else ''}) but priors are off. TNPs "
                f"carry an N(0,1) constraint by construction. Pass priors=1 (the "
                f"registry gives every TNP sigma=1 and leaves the lambdas free), "
                f"or drop them from fit_params."
            )
        if not priors:
            return
        overrides = (
            _parse_kv(prior_sigmas)
            if not isinstance(prior_sigmas, dict)
            else dict(prior_sigmas)
        )
        unknown = [k for k in overrides if k not in self._param_order]
        if unknown:
            raise KeyError(
                f"prior_sigmas: {unknown} are not fitted parameters "
                f"({list(self._param_order)})"
            )
        sigmas = np.empty(self.nparams, dtype=np.float64)
        for i, name in enumerate(self._param_order):
            s = overrides.get(name, adp.prior_sigma(name))
            sigmas[i] = np.nan if s is None else float(s)
        self.prior_sigmas = sigmas
        constrained = {
            n: s for n, s in zip(self._param_order, sigmas) if np.isfinite(s) and s > 0
        }
        print(
            f"[SCETlibADParamModel] Gaussian priors on {len(constrained)} "
            f"parameter(s): {constrained}",
            flush=True,
        )

    def _check_anchor_against_card(self, tol):
        """Compare the cache anchor with the card's propagated lambda_central."""
        meta = getattr(self.indata, "metadata", None)
        if not meta:
            print(
                "[SCETlibADParamModel] WARNING: the card carries no metadata, so "
                "the cache anchor could not be cross-checked against the "
                "nonperturbative values the card was produced with. Pass "
                "check_anchor=0 to silence, but verify by hand.",
                flush=True,
            )
            return
        card = np_anchor_from_meta(meta)
        if not card:
            print(
                "[SCETlibADParamModel] WARNING: the card records no "
                "nonperturbative anchor, so the cache anchor was NOT "
                "cross-checked.",
                flush=True,
            )
            return
        mismatched = []
        for card_key, rabbit in adp.LAMBDA_CENTRAL_KEYS.items():
            if card_key not in card or rabbit not in self.rabbit_names:
                continue
            want = float(card[card_key])
            have = float(self._anchor[self.rabbit_names.index(rabbit)])
            if abs(have - want) > tol * max(1.0, abs(want)):
                mismatched.append((rabbit, want, have))
        if mismatched:
            raise ValueError(
                "SCETlibADParamModel: the cache anchor does not match the card's "
                "NP central. The ratio would still be 1 at the start, so this "
                "fails silently -- rebuild the cache at the card's runcard "
                "values, or the card at the cache's.\n"
                + "\n".join(
                    f"    {n}: card {w:.6g} vs cache {h:.6g}" for n, w, h in mismatched
                )
            )

    def _check_no_inert_params(self):
        """Refuse a fitted parameter the prediction does not depend on.

        Not hypothetical: with the analysis runcard ``tnp_b_qqDS`` has an
        identically zero gradient for the Z (the channel it scales does not
        contribute), and ``lambda_inf`` / ``lambda_inf_nu`` are nearly inert at
        the anchor. A zero Jacobian column is a zero row and column of the NLL
        Hessian, i.e. a singular covariance and a meaningless "uncertainty".
        """
        p_start = np.asarray(self.xparamdefault.numpy(), dtype=np.float64)
        _, jac = self.core.values_and_jacobian(self._full_vector(p_start))
        J = self._fold(np.asarray(jac, dtype=np.float64))[:, self._fit_idx]
        scale = np.max(np.abs(J)) or 1.0
        dead = [
            n
            for i, n in enumerate(self._param_order)
            if np.max(np.abs(J[:, i])) <= 1e-12 * scale
        ]
        if dead:
            raise ValueError(
                f"SCETlibADParamModel: {len(dead)} fitted parameter(s) have an "
                f"identically zero derivative at the start point and would make "
                f"the covariance singular: {', '.join(dead)}. Drop them from "
                f"fit_params."
            )

    def _check_double_counting(self):
        """Refuse a card that still carries the templates our parameters replace."""
        systs = getattr(self.indata, "systs", None)
        if systs is None or len(systs) == 0:
            return
        names = [s.decode() if isinstance(s, bytes) else str(s) for s in systs]
        lowered = [(s, s.lower()) for s in names]
        for substring, what, applies in _CONFLICTS:
            if not applies(self._param_order):
                continue
            clash = [s for s, low in lowered if substring in low]
            if clash:
                raise ValueError(
                    f"[SCETlibADParamModel] {len(clash)} card syst(s) matching "
                    f"{substring!r} describe the same physics as the fitted "
                    f"{what}; running both double-counts. Remake the datacard "
                    f"with setupRabbit --excludeNuisances '.*{substring}.*' "
                    f"(case as in the card). Conflicting systs:\n"
                    + "\n".join(f"    {s}" for s in clash[:20])
                )

    # =========================================================================
    # evaluation
    # =========================================================================

    def _full_vector(self, fit_values):
        """Fitted values -> the complete SCETlib parameter vector."""
        p = self._p_base.copy()
        p[self._fit_idx] = np.asarray(fit_values, dtype=np.float64)
        return p

    def _sigma_gen_np(self, p_full):
        """sigma_gen on the gen grid (flattened) at a full SCETlib vector."""
        vals, _ = self.core.values_and_jacobian(p_full)
        return self._fold(np.asarray(vals, dtype=np.float64))

    def _eval_np(self, fit_values):
        """(value, J[, K]) on the gen grid, restricted to the fitted columns.

        Single-entry cached on the parameter bytes: rabbit evaluates the loss and
        the gradient through two separate ``tf.function``s at the same point, and
        the Hessian path adds a third, so without this the C++ replay would run
        several times per minimizer step.
        """
        p_fit = np.ascontiguousarray(fit_values, dtype=np.float64)
        key = p_fit.tobytes()
        if self._eval_cache is not None and self._eval_cache[0] == key:
            return self._eval_cache[1]

        p_full = self._full_vector(p_fit)
        vals, jac = self.core.values_and_jacobian(p_full)
        val = self._fold(np.asarray(vals, dtype=np.float64))
        J = self._fold(np.asarray(jac, dtype=np.float64))[:, self._fit_idx]
        out = (val, J)
        if self._curvature:
            H = self._fold(self.core.hessian(p_full))
            K = H[:, self._fit_idx][:, :, self._fit_idx]
            out = (val, J, np.ascontiguousarray(K))
        self._eval_cache = (key, out)
        return out

    def _sigma_gen_straightthrough(self, param):
        """Exact value with an exact local quadratic exposed to autodiff."""
        n_gen = int(np.prod(self.gen_shape))
        n_fit = self.nparams
        p = tf.cast(param, DTYPE)

        if self._curvature:
            val, J, K = tf.py_function(
                lambda t: self._eval_np(t.numpy()),
                [p],
                Tout=[DTYPE, DTYPE, DTYPE],
            )
            K.set_shape((n_gen, n_fit, n_fit))
        else:
            val, J = tf.py_function(
                lambda t: self._eval_np(t.numpy()), [p], Tout=[DTYPE, DTYPE]
            )
            K = None
        val.set_shape((n_gen,))
        J.set_shape((n_gen, n_fit))

        # d == 0 numerically, but carries a unit gradient, so the expression below
        # reproduces (value, J, K) exactly at the evaluation point while keeping
        # the py_function itself off the differentiated graph.
        d = p - tf.stop_gradient(p)
        out = tf.stop_gradient(val) + tf.linalg.matvec(tf.stop_gradient(J), d)
        if K is not None:
            out = out + 0.5 * tf.einsum("nij,i,j->n", tf.stop_gradient(K), d, d)
        return out

    def _sigma_gen_through(self, param):
        """sigma_gen with the SCETlib bridge INSIDE the differentiated graph.

        The reference implementation: ``ScetlibCachedXsecTF.__call__`` is an
        ordinary TF-differentiable function (its backward pass is itself a
        ``custom_gradient`` whose gradient contracts Hessian-vector products), so
        nested tapes work and TF drives every C++ call itself. Kept so the
        straight-through default is a checkable claim rather than an assertion.
        """
        p = tf.cast(param, DTYPE)
        # held = the non-fitted slots at their anchor, zero where we fit, so
        # held + S.p reconstructs the full vector. See _select on why this is a
        # matmul and not a scatter.
        held = self._p_base.copy()
        held[self._fit_idx] = 0.0
        p_full = tf.constant(held, dtype=DTYPE) + tf.linalg.matvec(
            tf.constant(self._select, dtype=DTYPE), p
        )
        return self._fold.fold_tf(self.core.tf_fn(p_full))

    def _ratio_from_param(self, param):
        """Per-fit-bin ratio to the anchor prediction, softly floored positive."""
        sigma_gen = (
            self._sigma_gen_through(param)
            if self._through
            else self._sigma_gen_straightthrough(param)
        )
        if self.gen_level:
            ratio = sigma_gen / self.sigma_gen_central_flat
        else:
            ratio = tf.linalg.matvec(self.R, sigma_gen) / self.sigma_reco_central
        # The rules put no wall in front of pathological parameters: a bad point
        # can drive sigma negative, which is a NaN Poisson NLL and a dead
        # gradient. Soft-floor so such a point is a large-but-finite penalty with
        # a usable gradient. The scale is far below any physical response, so the
        # anchor closure is untouched.
        scale = tf.constant(RATIO_FLOOR_SCALE, dtype=ratio.dtype)
        return tf.maximum(
            scale * tf.math.softplus(ratio / scale),
            tf.constant(RATIO_FLOOR_MIN, dtype=ratio.dtype),
        )

    def compute(self, param, full=False):
        """(N_bins, N_proc) multiplicative scaling; only the signal column moves."""
        ratio = self._ratio_from_param(param)
        ratio_col = tf.cast(tf.reshape(ratio, [-1, 1]), self.indata.dtype)
        rnorm = 1.0 + (ratio_col - 1.0) * self._signal_col_mask

        n_masked = int(getattr(self.indata, "nbinsmasked", 0) or 0)
        if full and n_masked:
            # Masked channels sit after the fit bins in the full tensor and carry
            # no model prediction, so they scale by 1.
            rnorm = tf.concat(
                [rnorm, tf.ones([n_masked, self.nproc], dtype=rnorm.dtype)], axis=0
            )
        return rnorm

    # =========================================================================
    # introspection (validation scripts)
    # =========================================================================

    def sigma_gen_at(self, **overrides):
        """sigma_gen on the gen grid at the anchor with named overrides applied.

        ``model.sigma_gen_at(lambda2_nu=0.12)`` -- used by the validation
        scripts to compare against a native SCETlib run.
        """
        p = self._p_base_anchor.copy()
        for name, val in overrides.items():
            if name not in self.rabbit_names:
                raise KeyError(f"sigma_gen_at: unknown parameter {name!r}")
            p[self.rabbit_names.index(name)] = float(val)
        return self._sigma_gen_np(p).reshape(self.gen_shape)

    def __deepcopy__(self, memo):
        """Copy the model but share the (immutable, expensive) SCETlib backend."""
        cls = self.__class__
        new = cls.__new__(cls)
        memo[id(self)] = new
        for k, v in self.__dict__.items():
            if k == "core":
                new.core = v
            elif k == "_eval_cache":
                new._eval_cache = None
            else:
                new.__dict__[k] = copy.deepcopy(v, memo)
        return new


def _parse_kv(spec):
    """``"a=1,b=2"`` (or a Mapping, or None) -> ``{"a": 1.0, "b": 2.0}``."""
    if not spec:
        return {}
    if not isinstance(spec, str):
        return {str(k): float(v) for k, v in dict(spec).items()}
    out = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"expected name=value, got {item!r}")
        k, v = item.split("=", 1)
        out[k.strip()] = float(v)
    return out
