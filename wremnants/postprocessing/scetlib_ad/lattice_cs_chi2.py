"""Exact ASWZ lattice CS-kernel chi2 as a rabbit ``Regularizer``.

Adds 1/2 * (chi2_lat - offset) to the NLL, where

    chi2_lat = min_k1  r^T C^-1 r ,
    r_i = pert_i(alpha_s, TNPs) + 1/2 gamma_nu^NP(b_i; lambda_inf_nu, lambda2_nu, lambda4_nu[, lambda6_nu])
          + k1 a_i / b_i - y_i

over the 21 ASWZ per-ensemble lattice points (arXiv:2402.06725): a SIMULTANEOUS fit of the lattice data with the Z
data. Conventions are those of WRemnantsHelpers studies/lattice-cs-kernel (260923-conventions-map,
260923-scetlib-kernel-fit, 261006-lattice-chi2-in-fit):

* The lattice gamma_q (MSbar, mu = 2 GeV) is SCETlib's gamma_zeta = gamma_nu / 2.
* pert_i: our SCETlib perturbative CS kernel at the lattice points (n_f = 5, mu = 2 GeV, N3LL, alpha_s(mZ) = 0.118,
  TNPs = 0, mu0 = ((b0/b_T)^4 + 1 GeV^4)^(1/4), sextic b* with b0/bmax_nu = 1 GeV). It is a table in the inputs file,
  with its derivatives; the evaluator it comes from is validated against SCETlib's Gamma_nu to 5e-11 (TNPs included).
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
    mu0scale   (alias ``pert``) b-space boundary-scale variation mu0 -> kappa mu0 (kappa = 1/2), a direct point shift.
               NOT USED by default (Luca, 2026-10-07): missing higher orders of the SCETlib kernel are the TNPs'
               job, and with ``pert=live`` the fitted TNPs move the kernel themselves. Opt-in only.
    none       stat only.

``pert=`` (how pert_i follows the fit; ``alphas=live`` is the older spelling of ``pert=alphas``):

    frozen  (default) the table: alpha_s(mZ) = 0.118, TNPs = 0.
    alphas  + (alpha_s - 0.118) d pert/d alpha_s (linear).
    live    every fitted parameter that enters the perturbative CS kernel at fixed (b_T, mu = 2 GeV):
              alpha_s (to second order), resumTNP_gamma_nu (3-loop gamma_nu boundary constant) and
              resumTNP_gamma_cusp (4-loop cusp), plus the alpha_s x TNP cross terms:
              pert + da dP_a + da^2/2 dP_aa + sum_t [ t dP_t + da t dP_at ],  da = alpha_s - 0.118, t = TNP value.
            The TNPs enter the kernel exactly linearly; the expansion is good to 4e-3 sigma_lat per point for
            |da| <= 0.002 and |t| <= 2 (2e-2 at |da| = 0.004). These are the only TNPs in SCETlib's
            ``qT::Gamma_nu::TNPs`` (gamma_cusp, gamma_nu); SCETlib has no beta-function TNP. Nothing else that is
            fitted enters gamma_zeta(b_T, mu): resumTransition2 is a qT-space profile transition point
            (Calculation_settings.transition_points), mu0_min and b0/bmax_nu are fixed settings, and the other TNPs
            (gamma_mu_q, s, b_*, h_*) sit in the hard/beam/soft boundary terms, not in the rapidity anomalous dimension.
    In ``alphas`` and ``live`` every physical value is computed with the PARAM MODEL'S OWN theta -> physical
    coefficients (``SCETlibADParamModel._rp_*``), applied to the regularizer's ``params`` argument = rabbit's
    ``Fitter.get_x()`` = [get_poi(), get_model_nui(), get_theta()]; ``get_poi()`` / ``get_model_nui()`` are exactly
    what the fitter hands the param model (``_compute_yields_noBBB``): the offset-applied PHYSICAL frame, never the
    blinded internal ``Fitter.x``. A TNP not in the fit is held at the correction runcard's value. Nothing alpha_s-
    or TNP-dependent is printed or stored by this class.

``ydata=``:

    lattice  (default) the ASWZ points.
    asimov   the points replaced by the model at the truth: pert(table point) + NP(lambda at the correction anchors) +
             k1_asimov a/b (``k1asimov=``, default 0.2). For Asimov closure (``-t 1 --toysDataMode expected
             --toysDataRandomize none --toysSystRandomize none``; rabbit never minimises ``-t -1``); use ``offset=0``.

``offset=``: ``min`` (default) is the lattice-only chi2_min at lambda_inf_nu = 2 (pert frozen, same covariance), so
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
         [syst=Jnf+Jbt|...] [pert=frozen|alphas|live] [offset=min|<float>] [ydata=lattice|asimov] [inputs=<npz>]

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
# Luca 2026-10-06: no k-form systematic; 2026-10-07: no mu0 scale variation (TNPs carry missing higher orders).
DEFAULT_SYST = "Jnf+Jbt"
SYST_COMPONENTS = ("J", "Jnf", "Jk", "Jbt", "direct_nf", "mu0scale", "pert", "none")
_J_ROWS = {"Jnf": 0, "Jk": 1, "Jbt": 2}
PERT_MODES = ("frozen", "alphas", "live")
# rabbit name -> (d pert / d t, d^2 pert / d alpha_s d t) keys in the inputs file
PERT_TNPS = {
    "resumTNP_gamma_nu": ("dpert_dtnp_nu1", "dpert_das_tnp_nu"),
    "resumTNP_gamma_cusp": ("dpert_dtnp_cusp1", "dpert_das_tnp_cusp"),
}
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
        elif c in ("mu0scale", "pert"):
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
        self.dpert_dalphas2 = d["dpert_dalphas2"] if "dpert_dalphas2" in d else None
        self.dpert_tnp = {
            n: (d[k1], d[k2])
            for n, (k1, k2) in PERT_TNPS.items()
            if k1 in d and k2 in d
        }
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
        self.c = self.pert - self.y  # r0 = c + np(lambda) + pert shift

    @staticmethod
    def np_zeta(u, linf, l2, l4, l6=0.0, tanh=np.tanh):
        return -0.5 * linf * tanh((l2 * u + l4 * u * u + l6 * u * u * u) / linf)

    def pert_shift(self, dalphas=0.0, tnps=None, mode="alphas"):
        """pert(alpha_s, TNPs) - table. mode 'alphas': linear in alpha_s only; 'live': the full expansion."""
        out = dalphas * self.dpert_dalphas
        if mode == "live":
            out = out + 0.5 * dalphas * dalphas * self.dpert_dalphas2
            for n, t in (tnps or {}).items():
                d1, dx = self.dpert_tnp[n]
                out = out + t * d1 + dalphas * t * dx
        return out

    def model(self, lam, k1=0.0, dalphas=0.0, tnps=None, mode="alphas"):
        """gamma_zeta at the points (k1 term included)."""
        g = self.np_zeta(
            self.u,
            lam["lambda_inf_nu"],
            lam["lambda2_nu"],
            lam.get("lambda4_nu", 0.0),
            lam.get("lambda6_nu", 0.0),
        )
        return self.pert + self.pert_shift(dalphas, tnps, mode) + g + k1 * self.v

    def r0(self, lam, dalphas=0.0, tnps=None, mode="alphas"):
        return self.model(lam, 0.0, dalphas, tnps, mode) - self.y

    def chi2(self, lam, dalphas=0.0, tnps=None, mode="alphas"):
        r = self.r0(lam, dalphas, tnps, mode)
        return float(r @ self.M @ r)

    def dchi2_dalphas(self, lam, dalphas=0.0):
        """Analytic d chi2 / d alpha_s(mZ) (linear pert response)."""
        return float(2.0 * self.r0(lam, dalphas) @ self.M @ self.dpert_dalphas)

    def k1hat(self, lam, dalphas=0.0, tnps=None, mode="alphas"):
        return float(-(self.Wv @ self.r0(lam, dalphas, tnps, mode)) / self.vWv)

    def fit(self, free=("lambda2_nu", "lambda4_nu"), fixed=None, x0=None):
        """Lattice-only minimum (scipy, pert frozen), for the offset and the closure."""
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


_KEYS = ("inputs", "syst", "offset", "tau", "alphas", "pert", "ydata", "k1asimov")


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
            alphas = kw.get("alphas", "frozen")
            if alphas not in ("frozen", "live"):
                raise ValueError(f"LatticeCSMapping: alphas={alphas!r} (frozen|live)")
            pert = kw.get("pert")
            if pert is not None and alphas == "live" and pert != "alphas":
                raise ValueError("LatticeCSMapping: give either pert= or alphas=live")
            self.pert = pert or ("alphas" if alphas == "live" else "frozen")
            if self.pert not in PERT_MODES:
                raise ValueError(f"LatticeCSMapping: pert={self.pert!r} {PERT_MODES}")
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


def model_param_map(param_model, name):
    """(c0, c1, c2) of the param model's OWN theta -> physical map for ``name``: physical = c0 + c1 t + c2 t^2.

    Read from SCETlibADParamModel's reparametrisation arrays (``_scetlib_order``, ``_rp_quad``, ``_rp_c``,
    ``_rp_id``, ``_rp_scale``), i.e. the exact coefficients ``_physical_tf`` applies before SCETlib sees the value.
    The identity branch is (0, scale, 0); a log map is refused."""
    order = list(getattr(param_model, "_scetlib_order", ()))
    if name not in order:
        raise ValueError(
            f"LatticeCSChi2: the fitter's param model has no SCETlib {name} "
            f"({type(param_model).__name__}); refusing to guess its map"
        )
    i = order.index(name)
    if bool(param_model._rp_quad[i]):
        return tuple(float(c) for c in param_model._rp_c[:, i])
    if bool(getattr(param_model, "_rp_id", np.ones(len(order), bool))[i]):
        scale = float(getattr(param_model, "_rp_scale", np.ones(len(order)))[i])
        return (0.0, scale, 0.0)
    raise ValueError(
        f"LatticeCSChi2: {name} has a non-polynomial map in the param model"
    )


def model_alphas_map(param_model):
    """(c0, c1) of the param model's linear alphaS map (kept for the alpha_s-only callers)."""
    c0, c1, c2 = model_param_map(param_model, "alphaS")
    if c2 != 0.0:
        raise ValueError(
            "LatticeCSChi2: alphaS map is not linear (unit) in the param model"
        )
    return c0, c1


def _make_regularizer_class():
    import tensorflow as tf

    from rabbit.regularization.regularizer import Regularizer
    from wremnants.postprocessing.scetlib_ad import np_damping_wall as wall
    from wremnants.postprocessing.scetlib_ad import params as adp
    from wremnants.postprocessing.scetlib_ad.response import corr_config_from_meta

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
            self.cfg = corr_config_from_meta(mapping.indata.metadata)["config"]
            self.inputs = getattr(mapping, "inputs", DEFAULT_INPUTS)
            self.syst = getattr(mapping, "syst", DEFAULT_SYST)
            self.ydata = getattr(mapping, "ydata", "lattice")
            if self.ydata not in ("lattice", "asimov"):
                raise ValueError(
                    f"LatticeCSChi2: ydata={self.ydata!r} (lattice|asimov)"
                )
            self.pert_mode = getattr(mapping, "pert", "frozen")
            if self.pert_mode not in PERT_MODES:
                raise ValueError(f"LatticeCSChi2: pert={self.pert_mode!r}")
            self.alphas_mode = "frozen" if self.pert_mode == "frozen" else "live"
            self.core = LatticeCSCore(self.inputs, self.syst)
            if self.pert_mode == "live" and (
                self.core.dpert_dalphas2 is None
                or set(self.core.dpert_tnp) != set(PERT_TNPS)
            ):
                raise ValueError(
                    f"LatticeCSChi2(pert=live): {self.inputs} lacks the second-order / TNP derivative tables"
                )
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
            # Offline (no fitter) maps: the correction-runcard ones, which is what the param model uses with its
            # default anchor_source=correction; set_expectations replaces them by the model's own.
            self.as_map = None
            if self.pert_mode != "frozen":
                spec = wall.physical_spec(self.cfg, "alphaS")
                if spec[0] != "quad" or spec[1][2] != 0.0:
                    raise ValueError(f"LatticeCSChi2: unexpected alphaS map {spec}")
                self.as_map_corr = (spec[1][0], spec[1][1])
                self._check_anchor(self.as_map_corr)
                self.as_map = self.as_map_corr
            self._ias, self._tau = None, None
            self._itnp, self._tnp_map, self._tnp_held = {}, {}, {}
            self._idx, self._held = {}, {}
            c = lambda a: tf.constant(np.asarray(a, float), dtype=dtype)  # noqa: E731
            self.tf_u, self.tf_c, self.tf_M = (
                c(self.core.u),
                c(self.core.c),
                c(self.core.M),
            )
            self.tf_das = c(self.core.dpert_dalphas)
            if self.pert_mode == "live":
                self.tf_daa = c(self.core.dpert_dalphas2)
                self.tf_dtnp = {
                    n: (c(d1), c(dx)) for n, (d1, dx) in self.core.dpert_tnp.items()
                }
            # Only alpha_s-INDEPENDENT quantities are printed (the offset is the frozen-table lattice-only minimum).
            print(
                f"[LatticeCSChi2] {len(self.core.y)} ASWZ points ({self.ydata}), inputs {self.inputs}, syst={self.syst}, "
                f"form {self.form}, k1 profiled analytically, pert {self.pert_mode} (table at alpha_s {ALPHAS_TABLE}"
                f", TNPs 0). offset = {self.offset:.6f}"
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
                    f"LatticeCSChi2(pert live): the alphaS anchor {m[0]} is not the pert table's {ALPHAS_TABLE}"
                )

        def _check_slot(self, pm, name, idx):
            """The regularizer's params = get_x(); its slot for a model parameter is the vector the model gets."""
            if (
                idx >= int(pm.npoi) + int(pm.npou)
                or list(np.asarray(pm.params).astype(str)).index(name) != idx
            ):
                raise ValueError(
                    f"LatticeCSChi2: {name} position in get_x() differs from the param model's"
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
            pm = getattr(fitter, "param_model", None)
            msg = ""
            self._ias = None
            self._itnp, self._tnp_map, self._tnp_held = {}, {}, {}
            if self.pert_mode != "frozen":
                if "alphaS" not in names:
                    raise ValueError(
                        f"LatticeCSChi2: pert={self.pert_mode} but alphaS is not a fit parameter"
                    )
                self._ias = names.index("alphaS")
                if pm is not None:
                    m = model_alphas_map(pm)
                    self._check_anchor(m)
                    self._check_slot(pm, "alphaS", self._ias)
                    self.as_map = m
                    same = (
                        max(
                            abs(m[0] - self.as_map_corr[0]),
                            abs(m[1] - self.as_map_corr[1]),
                        )
                        < 1e-15
                    )
                    msg = (
                        f"; alpha_s = the param model's own map ({m[0]:g} + {m[1]:g}*theta, "
                        f"{'==' if same else '!='} the correction-runcard map) on get_x()[{self._ias}] = get_poi()"
                    )
                else:
                    msg = f"; alpha_s from the correction-runcard map (no fitter on the stack) on params[{self._ias}]"
            if self.pert_mode == "live":
                parts = []
                for n in PERT_TNPS:
                    if n in names:
                        self._itnp[n] = names.index(n)
                        if pm is not None:
                            self._check_slot(pm, n, self._itnp[n])
                            self._tnp_map[n] = model_param_map(pm, n)
                            src = "param model"
                        else:
                            spec = adp.reparam(n)
                            if spec is not None:
                                raise ValueError(
                                    f"LatticeCSChi2: {n} has a REPARAM entry {spec}; offline map unknown"
                                )
                            self._tnp_map[n] = (0.0, 1.0, 0.0)
                            src = "identity (offline)"
                        parts.append(
                            f"{n} fitted, map {self._tnp_map[n]} ({src}), get_x()[{self._itnp[n]}]"
                        )
                    else:
                        v = adp.corr_anchor_value(self.cfg, n)
                        self._tnp_held[n] = float(v or 0.0)
                        parts.append(
                            f"{n} held at the correction's {self._tnp_held[n]:g}"
                        )
                msg += "; pert live in: alphaS, " + "; ".join(parts)
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
                f"{src}{msg}",
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

        def tnp_tf(self, params, name):
            """Physical TNP value (the model's map on get_x(), or the held correction value); never printed."""
            if name in self._itnp:
                c0, c1, c2 = (
                    tf.constant(v, dtype=self.dtype) for v in self._tnp_map[name]
                )
                t = params[self._itnp[name]]
                return c0 + c1 * t + c2 * t * t
            return tf.constant(self._tnp_held[name], dtype=self.dtype)

        def pert_shift_tf(self, params):
            da = self.dalphas_tf(params)
            out = da * self.tf_das
            if self.pert_mode == "live":
                out = out + 0.5 * da * da * self.tf_daa
                for n, (d1, dx) in self.tf_dtnp.items():
                    t = self.tnp_tf(params, n)
                    out = out + t * d1 + da * t * dx
            return out

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
                r = r + self.pert_shift_tf(params)
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
