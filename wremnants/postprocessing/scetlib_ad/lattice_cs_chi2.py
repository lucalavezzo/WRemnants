"""Exact ASWZ lattice CS-kernel chi2 as a rabbit ``Regularizer``.

Adds 1/2 * (chi2_lat - offset) to the NLL, where

    chi2_lat = min_k1  r^T C^-1 r ,
    r_i = pert_i(alpha_s) + 1/2 gamma_nu^NP(b_i; lambda_inf_nu, lambda2_nu, lambda4_nu[, lambda6_nu])
          + k1 a_i / b_i - y_i

over the 21 ASWZ per-ensemble lattice points (arXiv:2402.06725): a SIMULTANEOUS fit of the lattice data with the Z
data. Conventions are those of WRemnantsHelpers studies/lattice-cs-kernel (260923-conventions-map,
260923-scetlib-kernel-fit, 261006-lattice-chi2-in-fit):

* The lattice gamma_q (MSbar, mu = 2 GeV) is SCETlib's gamma_zeta = gamma_nu / 2.
* pert_i: our SCETlib perturbative CS kernel at the lattice points (n_f = 5, mu = 2 GeV, N3LL, alpha_s(mZ) = 0.118,
  mu0 = ((b0/b_T)^4 + 1 GeV^4)^(1/4), sextic b* with b0/bmax_nu = 1 GeV). It is a frozen table in the inputs file;
  the evaluator it comes from is validated against SCETlib's Gamma_nu to 5e-11.
* gamma_nu^NP = -linf tanh(P(u)/linf), P(u) = l2 u + l4 u^2 (+ l6 u^3 for tanh_6), u = b_T^2 at the BARE b_T
  (Gamma_nu.cpp).
* k1 (lattice-spacing artefact) is PROFILED ANALYTICALLY. The term is quadratic in k1, so
  chi2 = r0^T M r0 with M = W - (W v)(W v)^T / (v^T W v), W = C^-1, v = a/b, r0 = r(k1 = 0), and
  k1_hat = -(v^T W r0) / (v^T W v). This is exactly equivalent to a free k1 parameter: same minimum, and the
  curvature in the other parameters is the k1-marginal (Schur complement). A rabbit regularizer cannot own fit
  parameters.

Covariance, ``syst=`` (components joined by '+'; each adds sum delta delta^T to the block-diagonal lattice stat
covariance; a linear point shift in the covariance is identical to a profiled unit-Gaussian nuisance):

    Jnf, Jbt   (default: ``Jnf+Jbt``) the n_f-scheme (n_f = 5 matched at mu = 1 GeV) and b_T-window (drop
               b_T < 0.2 fm) groups of the 2D summary term: each chosen parameter shift d_g = (dl2, dl4) mapped onto
               the points as J_NP d_g at the lattice best fit. By Woodbury the profiled (l2, l4) covariance is then
               stat + sum d_g d_g^T, and the minimum and chi2_min do not move.
    Jk         the k-form group (k2-only instead of k1). NOT USED by default (Luca, 2026-10-06): the k1 a/b_T form is
               the lattice authors' prescription, not a choice of ours to vary. Opt-in only.
    J          all three (Jnf+Jk+Jbt), i.e. exactly the 2D summary card term's covariance (pre-2026-10-06 default).
    direct_nf  the n_f group as the direct point shift pert(n_f = 5 matched at mu = 1 GeV) - pert(n_f = 5).
    pert       perturbative-kernel scale variation: the b-space boundary scale mu0 -> kappa mu0, kappa = 2 or 1/2 (the
               larger of the two under the stat metric), as a direct point shift.
    none       stat only.

``alphas=``:

    frozen  (default) pert at alpha_s(mZ) = 0.118.
    live    pert + (alpha_s - 0.118) d pert / d alpha_s(mZ), linear (exact to 2.2e-3 sigma_lat per point over
            +-0.002). The physical alpha_s is computed with the PARAM MODEL'S OWN reparametrisation coefficients for
            ``alphaS`` (``SCETlibADParamModel._rp_*``, anchor included), applied to the regularizer's ``params``
            argument. That argument is rabbit's ``Fitter.get_x()`` = [get_poi(), get_model_nui(), get_theta()], and
            ``get_poi()`` is exactly what the fitter hands the param model (``_compute_yields_noBBB``). It is the
            offset-applied, PHYSICAL frame, never the blinded internal ``Fitter.x``. Nothing alpha_s-dependent is
            printed or stored by this class.

``ydata=``:

    lattice  (default) the ASWZ points.
    asimov   the points replaced by the model at the truth: pert(0.118) + NP(lambda at the correction anchors) +
             k1_asimov a/b (``k1asimov=``, default 0.2). Used for Asimov closure (``-t -1``); use with ``offset=0``.

``offset=``: ``min`` (default) is the lattice-only chi2_min at lambda_inf_nu = 2, alpha_s frozen, same covariance, so
the term is 1/2 Delta chi2 (zero at the lattice best fit, like a Gaussian summary). Or a float.

THE exp(2 tau) COMPENSATION. rabbit multiplies EVERY regularizer penalty by exp(2 tau) (one common ``fitter.tau`` =
``--regularizationStrength``, meant for NPDampingWall). A likelihood term must not be scaled, so this term divides it
back out with the LIVE ``fitter.tau`` variable, found on the call stack at ``set_expectations`` (rabbit calls it from
``Fitter.arm_regularizers``). Offline (no fitter), ``tau=<float>`` is used, else 1. Both present and different raises.

theta -> physical for the CS lambdas: the same map as NPDampingWall (``resolve_wall_inputs``), resolved by name;
held lambdas (lambda_inf_nu) at the correction anchor.

Invoke (composes with the wall; each on its own -r):

    rabbit_fit.py ... --regularizationStrength 8 \\
      -r wremnants.postprocessing.scetlib_ad.np_damping_wall.NPDampingWall \\
         wremnants.postprocessing.scetlib_ad.np_damping_wall.NPDampingMapping margin=0 \\
      -r wremnants.postprocessing.scetlib_ad.lattice_cs_chi2.LatticeCSChi2 \\
         wremnants.postprocessing.scetlib_ad.lattice_cs_chi2.LatticeCSMapping \\
         [syst=Jnf+Jbt|Jnf+Jbt+pert|...] [alphas=frozen|live] [offset=min|<float>] [ydata=lattice|asimov] [inputs=<npz>] [tau=<f>]

Inputs file: ``data/lattice_aswz_inputs.npz`` next to this module (built by
WRemnantsHelpers studies/lattice-cs-kernel/261006-lattice-chi2-in-fit/scripts/build_inputs.py from the ASWZ files and
the validated evaluator; its .json sidecar records the provenance).
"""

import inspect
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_CANDIDATES = (
    os.path.join(_HERE, "data", "lattice_aswz_inputs.npz"),
    os.path.join(os.path.dirname(_HERE), "lattice_aswz_inputs.npz"),
)
DEFAULT_INPUTS = next((p for p in _CANDIDATES if os.path.exists(p)), _CANDIDATES[0])
CS_NAMES = ("lambda_inf_nu", "lambda2_nu", "lambda4_nu", "lambda6_nu")
DEFAULT_SYST = (
    "Jnf+Jbt"  # Luca 2026-10-06: no k-form systematic (the theorists' prescription)
)
SYST_COMPONENTS = ("J", "Jnf", "Jk", "Jbt", "direct_nf", "pert", "none")
_J_ROWS = {"Jnf": 0, "Jk": 1, "Jbt": 2}
ALPHAS_TABLE = (
    0.118  # alpha_s(mZ) of the pert table (= the AD cache / correction alphas_mu0)
)
LINF_REF = 2.0  # lambda_inf_nu of the lattice-only reference minimum (offset=min)


def _syst_rows(d, syst):
    comps = [c for c in str(syst).split("+") if c]
    bad = [c for c in comps if c not in SYST_COMPONENTS]
    if bad or not comps:
        raise ValueError(f"syst={syst!r}: components must be from {SYST_COMPONENTS}")
    rows = []
    for c in comps:
        if c == "J":
            rows += list(d["syst_J"])
        elif c in _J_ROWS:
            rows.append(d["syst_J"][_J_ROWS[c]])
        elif c == "direct_nf":
            rows.append(d["syst_direct_nf"][0])
        elif c == "pert":
            rows += list(d["syst_pert"])
    return rows


class LatticeCSCore:
    """Backend-agnostic exact lattice chi2 (numpy here; the TF path uses the same constants)."""

    def __init__(self, inputs=DEFAULT_INPUTS, syst=DEFAULT_SYST, y=None):
        d = np.load(inputs)
        self.inputs, self.syst = inputs, syst
        self.b_fm, self.a_fm = d["b_fm"], d["a_fm"]
        self.y = np.array(d["y"] if y is None else y, float)
        self.u = (self.b_fm * float(d["fm_to_gevinv"])) ** 2  # GeV^-2
        self.v = self.a_fm / self.b_fm
        self.pert = d["pert"]
        self.dpert_dalphas = d["dpert_dalphas"]
        cov = np.array(d["cov_stat"], float)
        self.csyst = sum(
            (np.outer(s, s) for s in _syst_rows(d, syst)), np.zeros_like(cov)
        )
        self.cov = cov + self.csyst
        W = np.linalg.inv(self.cov)
        self.W = 0.5 * (W + W.T)
        self.Wv = self.W @ self.v
        self.vWv = float(self.v @ self.Wv)
        M = self.W - np.outer(self.Wv, self.Wv) / self.vWv
        self.M = 0.5 * (M + M.T)
        self.c = self.pert - self.y  # r0 = c + np(lambda) + dalphas * dpert

    @staticmethod
    def np_zeta(u, linf, l2, l4, l6=0.0, tanh=np.tanh):
        return -0.5 * linf * tanh((l2 * u + l4 * u * u + l6 * u * u * u) / linf)

    def model(self, lam, k1=0.0, dalphas=0.0):
        """gamma_zeta at the points (k1 term included)."""
        g = self.np_zeta(
            self.u,
            lam["lambda_inf_nu"],
            lam["lambda2_nu"],
            lam.get("lambda4_nu", 0.0),
            lam.get("lambda6_nu", 0.0),
        )
        return self.pert + dalphas * self.dpert_dalphas + g + k1 * self.v

    def r0(self, lam, dalphas=0.0):
        return self.model(lam, 0.0, dalphas) - self.y

    def chi2(self, lam, dalphas=0.0):
        r = self.r0(lam, dalphas)
        return float(r @ self.M @ r)

    def dchi2_dalphas(self, lam, dalphas=0.0):
        """Analytic d chi2 / d alpha_s(mZ) (linear pert response)."""
        return float(2.0 * self.r0(lam, dalphas) @ self.M @ self.dpert_dalphas)

    def k1hat(self, lam, dalphas=0.0):
        return float(-(self.Wv @ self.r0(lam, dalphas)) / self.vWv)

    def fit(self, free=("lambda2_nu", "lambda4_nu"), fixed=None, x0=None):
        """Lattice-only minimum (scipy, alpha_s frozen), for the offset and the closure."""
        from scipy.optimize import minimize

        fixed = dict(
            dict(lambda_inf_nu=LINF_REF, lambda4_nu=0.0, lambda6_nu=0.0),
            **(fixed or {}),
        )
        x0 = np.array(x0 if x0 is not None else [0.15, 0.0][: len(free)], float)

        def f(x):
            return self.chi2(dict(fixed, **dict(zip(free, x))))

        best = None
        for s in [x0, x0 * 0.5 + 0.05, x0 + 0.01]:
            r = minimize(
                f,
                s,
                method="Nelder-Mead",
                options=dict(xatol=1e-12, fatol=1e-13, maxiter=40000, maxfev=80000),
            )
            r = minimize(f, r.x, method="BFGS", options=dict(gtol=1e-11))
            if best is None or r.fun < best.fun:
                best = r
        return best, dict(fixed, **dict(zip(free, best.x)))


_KEYS = ("inputs", "syst", "offset", "tau", "alphas", "ydata", "k1asimov")


def _make_mapping_class():
    from rabbit.mappings.mapping import BaseMapping

    class LatticeCSMapping(BaseMapping):
        """Carries the -r options (key=value) to LatticeCSChi2; see the module docstring."""

        def __init__(self, indata, key, **kw):
            super().__init__(indata, key)
            self.indata = indata
            self.inputs = kw.get("inputs", DEFAULT_INPUTS)
            self.syst = kw.get("syst", DEFAULT_SYST)
            self.offset = kw.get("offset", "min")
            self.tau = kw.get("tau")
            self.alphas = kw.get("alphas", "frozen")
            self.ydata = kw.get("ydata", "lattice")
            self.k1asimov = float(kw.get("k1asimov", 0.2))

        @classmethod
        def parse_args(cls, indata, *args):
            kw = {}
            for a in args:
                if "=" not in a:
                    raise ValueError(f"LatticeCSMapping: args are key=value, got {a!r}")
                k, v = a.split("=", 1)
                if k not in _KEYS:
                    raise ValueError(
                        f"LatticeCSMapping: unknown key {k!r}; use {_KEYS}"
                    )
                kw[k] = float(v) if k == "tau" else v
            return cls(indata, " ".join([cls.__name__, *args]), **kw)

    return LatticeCSMapping


def _find_fitter():
    """The rabbit Fitter calling set_expectations (Fitter.arm_regularizers), or None."""
    fr = inspect.currentframe()
    try:
        f = fr.f_back
        for _ in range(8):
            if f is None:
                return None
            obj = f.f_locals.get("self")
            if obj is not None and hasattr(obj, "tau") and hasattr(obj, "regularizers"):
                return obj
            f = f.f_back
    finally:
        del fr
    return None


def model_alphas_map(param_model):
    """(c0, c1) of the param model's OWN linear theta -> physical map for alphaS: physical = c0 + c1 * theta.

    Read from SCETlibADParamModel's reparametrisation arrays (``_scetlib_order``, ``_rp_quad``, ``_rp_c``), i.e. the
    exact coefficients ``_physical_tf`` applies before SCETlib sees the value."""
    order = list(getattr(param_model, "_scetlib_order", ()))
    if "alphaS" not in order:
        raise ValueError(
            "LatticeCSChi2(alphas=live): the fitter's param model has no SCETlib alphaS "
            f"({type(param_model).__name__}); refusing to guess its map"
        )
    i = order.index("alphaS")
    if not bool(param_model._rp_quad[i]) or float(param_model._rp_c[2, i]) != 0.0:
        raise ValueError(
            "LatticeCSChi2(alphas=live): alphaS map is not linear (unit) in the param model"
        )
    return float(param_model._rp_c[0, i]), float(param_model._rp_c[1, i])


def _make_regularizer_class():
    import tensorflow as tf

    from rabbit.regularization.regularizer import Regularizer
    from wremnants.postprocessing.scetlib_ad import np_damping_wall as wall

    class LatticeCSChi2(Regularizer):
        """1/2 (chi2_lat - offset) / exp(2 tau): exact ASWZ lattice chi2, k1 profiled; see the module docstring."""

        needs_observables = False

        def __init__(self, mapping, dtype):
            super().__init__(mapping, dtype)
            self.dtype = dtype
            self.mapping = mapping
            inp = wall.resolve_wall_inputs(mapping.indata)
            if inp["np_model_nu"] not in ("tanh_2", "tanh_6"):
                raise ValueError(
                    f"LatticeCSChi2: np_model_nu={inp['np_model_nu']!r} (need tanh_2 / tanh_6)"
                )
            self.form = inp["np_model_nu"]
            self.specs = {n: inp["specs"][n] for n in CS_NAMES if n in inp["specs"]}
            self.anchors = {
                n: inp["anchors"][n] for n in CS_NAMES if n in inp["anchors"]
            }
            self.inputs = getattr(mapping, "inputs", DEFAULT_INPUTS)
            self.syst = getattr(mapping, "syst", DEFAULT_SYST)
            self.ydata = getattr(mapping, "ydata", "lattice")
            if self.ydata not in ("lattice", "asimov"):
                raise ValueError(
                    f"LatticeCSChi2: ydata={self.ydata!r} (lattice|asimov)"
                )
            self.core = LatticeCSCore(self.inputs, self.syst)
            if self.ydata == "asimov":
                truth = {
                    n: float(self.anchors[n])
                    for n in CS_NAMES
                    if self.anchors.get(n) is not None
                }
                truth.setdefault("lambda4_nu", 0.0)
                y = self.core.model(truth, k1=float(getattr(mapping, "k1asimov", 0.2)))
                self.core = LatticeCSCore(self.inputs, self.syst, y=y)
            off = getattr(mapping, "offset", "min")
            if off == "min":
                res, lam = self.core.fit()
                self.offset, self.offset_point = float(res.fun), lam
            else:
                self.offset, self.offset_point = float(off), None
            self.tau_arg = getattr(mapping, "tau", None)
            self.alphas_mode = getattr(mapping, "alphas", "frozen")
            if self.alphas_mode not in ("frozen", "live"):
                raise ValueError(
                    f"LatticeCSChi2: alphas={self.alphas_mode!r} (frozen|live)"
                )
            # alphas=live without a fitter (offline harness): the correction-runcard map, which is what the param
            # model uses with its default anchor_source=correction; set_expectations replaces it by the model's own.
            self.as_map = None
            if self.alphas_mode == "live":
                from wremnants.postprocessing.scetlib_ad.response import (
                    corr_config_from_meta,
                )

                spec = wall.physical_spec(
                    corr_config_from_meta(mapping.indata.metadata)["config"], "alphaS"
                )
                if spec[0] != "quad" or spec[1][2] != 0.0:
                    raise ValueError(f"LatticeCSChi2: unexpected alphaS map {spec}")
                self.as_map_corr = (spec[1][0], spec[1][1])
                self._check_anchor(self.as_map_corr)
                self.as_map = self.as_map_corr
            self._ias, self._tau = None, None
            self._idx, self._held = {}, {}
            c = lambda a: tf.constant(np.asarray(a, float), dtype=dtype)  # noqa: E731
            self.tf_u, self.tf_c, self.tf_M = (
                c(self.core.u),
                c(self.core.c),
                c(self.core.M),
            )
            self.tf_das = c(self.core.dpert_dalphas)
            # Only alpha_s-INDEPENDENT quantities are printed (the offset is the frozen-table lattice-only minimum).
            print(
                f"[LatticeCSChi2] {len(self.core.y)} ASWZ points ({self.ydata}), inputs {self.inputs}, syst={self.syst}, "
                f"form {self.form}, k1 profiled analytically, alpha_s {self.alphas_mode} (table at {ALPHAS_TABLE}). "
                f"offset = {self.offset:.6f}"
                + (
                    f" (lattice-only min at lambda_inf_nu={LINF_REF}, {self.offset_point})"
                    if self.offset_point
                    else ""
                )
                + ". theta->physical: "
                + ", ".join(f"{n}: {self.specs[n]}" for n in self.specs),
                flush=True,
            )

        @staticmethod
        def _check_anchor(m):
            if abs(m[0] - ALPHAS_TABLE) > 1e-12:
                raise ValueError(
                    f"LatticeCSChi2(alphas=live): the alphaS anchor {m[0]} is not the pert table's {ALPHAS_TABLE}"
                )

        def set_expectations(self, initial_params, initial_observables, parms=None):
            if parms is None:
                raise ValueError(
                    "LatticeCSChi2 needs parms= (names) to resolve the CS lambdas"
                )
            names = list(np.asarray(parms).astype(str))
            self._idx, self._held = {}, {}
            for n in CS_NAMES:
                if n == "lambda6_nu" and self.form != "tanh_6":
                    continue
                if n in names:
                    self._idx[n] = names.index(n)
                else:
                    a = self.anchors.get(n)
                    if a is None:
                        raise ValueError(
                            f"LatticeCSChi2: {n} neither fitted nor carried by the correction runcard"
                        )
                    self._held[n] = float(a)
            if "lambda2_nu" not in self._idx:
                raise ValueError(
                    "LatticeCSChi2: lambda2_nu is not a fit parameter -- the term would be a constant"
                )
            fitter = _find_fitter()
            msg_as = ""
            self._ias = None
            if self.alphas_mode == "live":
                if "alphaS" not in names:
                    raise ValueError(
                        "LatticeCSChi2: alphas=live but alphaS is not a fit parameter"
                    )
                self._ias = names.index("alphaS")
                if fitter is not None:
                    pm = fitter.param_model
                    m = model_alphas_map(pm)
                    self._check_anchor(m)
                    # the regularizer's params = get_x(); its alphaS slot is get_poi()'s, the vector the model gets
                    if (
                        self._ias >= int(pm.npoi) + int(pm.npou)
                        or list(np.asarray(pm.params).astype(str)).index("alphaS")
                        != self._ias
                    ):
                        raise ValueError(
                            "LatticeCSChi2: alphaS position in get_x() differs from the param model's"
                        )
                    self.as_map = m
                    same = (
                        max(
                            abs(m[0] - self.as_map_corr[0]),
                            abs(m[1] - self.as_map_corr[1]),
                        )
                        < 1e-15
                    )
                    msg_as = (
                        f"; alpha_s = the param model's own map ({m[0]:g} + {m[1]:g}*theta, "
                        f"{'==' if same else '!='} the correction-runcard map) on get_x()[{self._ias}] = get_poi()"
                    )
                else:
                    msg_as = f"; alpha_s from the correction-runcard map (no fitter on the stack) on params[{self._ias}]"
            tau = fitter.tau if fitter is not None else None
            if tau is not None:
                if (
                    self.tau_arg is not None
                    and abs(float(tau.numpy()) - self.tau_arg) > 1e-12
                ):
                    raise ValueError(
                        f"LatticeCSChi2: tau={self.tau_arg} on the -r line but fitter.tau = {float(tau.numpy())}"
                    )
                self._tau = tau
                src = f"fitter.tau (live variable, now {float(tau.numpy()):g})"
            else:
                self._tau = None
                src = (
                    f"tau={self.tau_arg} from the -r line"
                    if self.tau_arg is not None
                    else "NONE (scale 1)"
                )
            print(
                f"[LatticeCSChi2] armed: fitted {sorted(self._idx)}, held {self._held}; exp(2 tau) compensation from "
                f"{src}{msg_as}",
                flush=True,
            )

        def physical(self, params):
            vals = {n: tf.constant(v, dtype=self.dtype) for n, v in self._held.items()}
            for n, i in self._idx.items():
                vals[n] = wall.physical_from_theta(self.specs[n], params[i], exp=tf.exp)
            return vals

        def dalphas_tf(self, params):
            """alpha_s(mZ) - 0.118 from params (= get_x(), physical frame); never printed."""
            c0, c1 = self.as_map
            return (
                tf.constant(c0, dtype=self.dtype)
                + tf.constant(c1, dtype=self.dtype) * params[self._ias]
            ) - (tf.constant(ALPHAS_TABLE, dtype=self.dtype))

        def chi2_tf(self, params):
            lam = self.physical(params)
            z = tf.constant(0.0, dtype=self.dtype)
            g = LatticeCSCore.np_zeta(
                self.tf_u,
                lam["lambda_inf_nu"],
                lam["lambda2_nu"],
                lam.get("lambda4_nu", z),
                lam.get("lambda6_nu", z),
                tanh=tf.tanh,
            )
            r = self.tf_c + g
            if self._ias is not None:
                r = r + self.dalphas_tf(params) * self.tf_das
            return tf.tensordot(r, tf.linalg.matvec(self.tf_M, r), 1)

        def compute_nll_penalty(self, params, observables):
            half = 0.5 * (
                self.chi2_tf(params) - tf.constant(self.offset, dtype=self.dtype)
            )
            if self._tau is not None:
                return half * tf.exp(-2.0 * tf.cast(self._tau, self.dtype))
            if self.tau_arg is not None:
                return half * tf.constant(np.exp(-2.0 * self.tau_arg), dtype=self.dtype)
            return half

    return LatticeCSChi2


def __getattr__(name):
    if name == "LatticeCSMapping":
        cls = _make_mapping_class()
        globals()[name] = cls
        return cls
    if name == "LatticeCSChi2":
        cls = _make_regularizer_class()
        globals()[name] = cls
        return cls
    raise AttributeError(name)
