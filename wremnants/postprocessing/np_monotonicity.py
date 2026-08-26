"""
Tanh-monotonicity wall for the SCETlib NP model in rabbit fits.

Adds a custom rabbit Regularizer + Mapping that, at fit time, maps the
six floating NP nuisance pulls back to physical lambda_i / Lambda_i values
and adds a hinge-loss penalty enforcing the tanh-monotonicity walls derived
in agents/knowledge/30_physics_global/np_parametrization_constraints.md.

The walls (criterion (b), monotonicity) are:
  CS:  lambda_4 >= -sqrt(3 * lambda_2 * lambda_6)
  TMD: 3*Lambda_4 + L_2(y)^3 >= -sqrt(20 * Lambda_6 * L_2(y))
       evaluated at y=0 and y=y_max and summed (safe under floating
       Delta_Lambda_2 of either sign); L_2(y) = Lambda_2 + Delta_Lambda_2 * y^2.

Criterion (a) (sign-preservation) replaces sqrt(3)->sqrt(4), sqrt(20)->sqrt(36).

The physical-parameter values that anchor the linearization theta -> param
are hard-coded below in PARAM_MAP. An optional verification helper
verify_param_map_against_hist(corr_hist) cross-checks them against the
histogram syst-axis labels written by the histmaker; setupRabbit calls
this so any drift in the templates surfaces as a loud error rather than
silent regularizer mismatch.

Supported NP card:  CS  = LatticeNoConstraints
                    TMD = Delta_Lambda

References:
  AN-25-085 theory.tex Eqs. eq:npgamma, eq:npf.
"""

import numpy as np

# rabbit/TF imports are deferred to inside the lazy class factories so the
# module is importable on the histmaker side without rabbit available.


# ---------------------------------------------------------------------------
# Hard-coded physics constants and per-nuisance template map.
#
# Sources for the template Up/Down labels and conventions:
#   CS  (LatticeNoConstraints): rabbit_theory_helper.py:686-711
#   TMD (Delta_Lambda):         rabbit_theory_helper.py:827-882
#
# Note the CS convention: var_names iterates ["Up","Down"] *outside* per-name,
# so the FIRST entry of lattice_vals becomes the "Up" template -- which is the
# *smaller* physical value. The piecewise param(theta) formula handles this via
# negative delta_up/delta_down without special-casing.
# ---------------------------------------------------------------------------

# Fixed (non-floating) higher-order coefficients added to the AN parametrization
# to enforce the correct large-bT asymptote.
LAMBDA_6 = 0.0007       # coefficient of b_T^6 / lambda_inf  (CS)
BIG_LAMBDA_6 = 0.016    # coefficient of b_T^5 / Lambda_inf  (TMD)
Y_MAX = 2.5             # operational |y_ll| range for the binding-y evaluation

# Default criterion for the wall: "b" = monotonicity (stricter), "a" = sign-preservation.
DEFAULT_CRITERION = "b"

# rabbit nuisance names produced by the supported NP card.
#
# CS side: add_gamma_np_uncertainties() with np_model='LatticeNoConstraints' adds 3
# correlated nuisances with the "scetlibNPgamma..." prefix and Lambda<N> uppercase.
#
# TMD side: add_uncorrelated_np_uncertainties() with np_model='LatticeNoConstraints'
# adds 3 nuisances decorrelated across chargeVgenNP and integrated over absYVgenNP,
# producing names "chargeVgenNP<bin>scetlibNP<sample_label><nuisance>" where
# <sample_label>='Z' for dilepton Z-only fits and <nuisance> is lowercase
# (lambda2 / delta_lambda2 / lambda4). For Z-only there is one charge bin, so the
# nuisance names contain "chargeVgenNP0" exactly.
NUIS_LAMBDA_2_GAMMA    = "scetlibNPgammaLambda2"
NUIS_LAMBDA_4_GAMMA    = "scetlibNPgammaLambda4"
NUIS_LAMBDA_INF_GAMMA  = "scetlibNPgammaLambdaInf"
NUIS_LAMBDA_2_TMD      = "chargeVgenNP0scetlibNPZlambda2"
NUIS_DELTA_LAMBDA_2    = "chargeVgenNP0scetlibNPZdelta_lambda2"
NUIS_LAMBDA_4_TMD      = "chargeVgenNP0scetlibNPZlambda4"


# Per-nuisance entry layout:
#   "param":     human-readable parameter name
#   "nominal":   physical value at theta=0 (AN-25-085 central)
#   "up_value":  physical value at theta=+1 (the "Up" template)
#   "down_value":physical value at theta=-1 (the "Down" template)
#   "hist_up_label", "hist_down_label": labels in the histmaker syst axis, used
#                                       only by verify_param_map_against_hist.
#
# CS-SIDE Up/Down convention (post-2026-05-11 fix):
# rabbit_theory_helper.py:706-715 now iterates `for direction in ["Down","Up"]`
# in the LatticeNoConstraints branch, so the lattice_vals (which are listed
# [smaller, larger]) are paired:
#     scetlibNPgammaLambda2Down <-> lambda2_nu0.0538   (smaller, ≈ −1σ_lat)
#     scetlibNPgammaLambda2Up   <-> lambda2_nu0.1202   (larger,  ≈ +1σ_lat)
# Positive rabbit pull → larger physical λ_2 (standard convention).
# IMPORTANT: HDF5s produced BEFORE this fix used the inverted convention.
# If verifying old HDF5s, swap the hist_up_label/hist_down_label entries.
#
# Linear interpolation:
#     param(theta) = nominal + max(theta,0)*(up_value - nominal)
#                            - max(-theta,0)*(nominal - down_value)
# Works for both symmetric and asymmetric Up/Down ranges, and naturally handles
# the CS inverted Up/Down convention (where up_value < nominal < down_value).
#
# IMPORTANT: these up/down physical values assume the histmaker's scale factor is 1.
# add_gamma_np_uncertainties() currently applies scale=10 internally for
# LatticeNoConstraints (rabbit_theory_helper.py:738-739) which inflates the CS-side
# template magnitude by 10x. To use this PARAM_MAP unchanged, re-run setupRabbit
# with `--scaleParams 'scetlibNPgamma=0.1'` to undo the inflation (multiplies the
# kfactor by 0.1, which composes with the helper's 10x to give net 1x).
PARAM_MAP = {
    # CS side: standard convention (Up = larger physical value) after the
    # 2026-05-11 fix in rabbit_theory_helper.py (iteration order ["Down","Up"]).
    NUIS_LAMBDA_2_GAMMA: {
        "param": "lambda_2",
        "nominal": 0.0870, "up_value": 0.1202, "down_value": 0.0538,
        "hist_up_label": "lambda2_nu0.1202",
        "hist_down_label": "lambda2_nu0.0538",
    },
    NUIS_LAMBDA_4_GAMMA: {
        "param": "lambda_4",
        "nominal": 0.0074, "up_value": 0.014, "down_value": 0.0008,
        "hist_up_label": "lambda4_nu0.014",
        "hist_down_label": "lambda4_nu0.0008",
    },
    NUIS_LAMBDA_INF_GAMMA: {
        "param": "lambda_inf",
        "nominal": 1.6853, "up_value": 2.1922, "down_value": 1.1784,
        "hist_up_label": "lambda_inf_nu2.1922",
        "hist_down_label": "lambda_inf_nu1.1784",
    },
    # TMD side (LatticeNoConstraints uncorrelated path): Up = larger physical value.
    # Histmaker label convention is mixed:
    #   lambda2: tail is the absolute value (so down=0.0, up=0.5)
    #   delta_lambda2: tail is the variation from the AN nominal (so down=0.105, up=0.145)
    #   lambda4: tail is the absolute value, asymmetric (down=0.01, up=0.12)
    NUIS_LAMBDA_2_TMD: {
        "param": "Lambda_2",
        "nominal": 0.25, "up_value": 0.50, "down_value": 0.00,
        "hist_up_label": "lambda20.5",
        "hist_down_label": "lambda20.0",
    },
    NUIS_DELTA_LAMBDA_2: {
        "param": "Delta_Lambda_2",
        # NB: In the LatticeNoConstraints branch the histmaker uses nominal=0
        # for delta_lambda2 (no y-dependence in central prediction); templates
        # are at ±0.02 absolute. This is NOT the AN-quoted 0.125 baseline,
        # which applies to the uppercase "Delta_Lambda" branch only.
        "nominal": 0.0, "up_value": 0.02, "down_value": -0.02,
        "hist_up_label": "delta_lambda20.02",
        "hist_down_label": "delta_lambda2-0.02",
    },
    NUIS_LAMBDA_4_TMD: {
        "param": "Lambda_4",
        "nominal": 0.06, "up_value": 0.12, "down_value": 0.01,   # asymmetric, range [0.01,0.12]
        "hist_up_label": "lambda40.12",
        "hist_down_label": "lambda40.01",
    },
}


# ---------------------------------------------------------------------------
# Optional: cross-check the hard-coded values against the histogram labels.
# Call this from setupRabbit / rabbit_theory_helper.py at corr_hist time to
# fail loudly if histmaker templates ever drift away from what PARAM_MAP assumes.
# ---------------------------------------------------------------------------

def _parse_hist_label(label, prefix):
    s = str(label)
    if not s.startswith(prefix):
        raise ValueError(f"Label '{label}' does not start with prefix '{prefix}'")
    return float(s[len(prefix):])


def verify_param_map_against_hist(corr_hist, syst_ax="vars", tolerance=1e-3):
    """Cross-check PARAM_MAP physical values against the syst-axis labels of corr_hist.

    For symmetric cases (CS lambdas, TMD Lambda_2, Delta_Lambda_2), checks that the
    label-encoded value matches the hard-coded up/down value. For absolute-value
    labels (CS '<name>_nu<value>', TMD 'Lambda4<value>') the comparison is direct;
    for delta-style TMD labels ('Lambda2<delta>', 'Delta_Lambda2<delta>') the
    comparison is against (value - nominal).

    Raises ValueError on any mismatch; returns silently otherwise.
    """
    label_set = {str(x) for x in corr_hist.axes[syst_ax]}

    # Prefix used to strip leading alpha; the value remainder is float()-parsable.
    # Pairs are (prefix, kind) where kind is 'abs' (tail is absolute value) or
    # 'delta' (tail is variation from nominal).
    label_kinds = {
        NUIS_LAMBDA_2_GAMMA:    ("lambda2_nu",     "abs"),
        NUIS_LAMBDA_4_GAMMA:    ("lambda4_nu",     "abs"),
        NUIS_LAMBDA_INF_GAMMA:  ("lambda_inf_nu",  "abs"),
        # TMD lowercase paths (LatticeNoConstraints uncorrelated branch):
        NUIS_LAMBDA_2_TMD:      ("lambda2",        "abs"),    # tail is absolute value
        NUIS_DELTA_LAMBDA_2:    ("delta_lambda2",  "abs"),    # tail is absolute (nominal=0 in this branch)
        NUIS_LAMBDA_4_TMD:      ("lambda4",        "abs"),    # tail is absolute value
    }

    mismatches = []
    for nuis_name, info in PARAM_MAP.items():
        prefix, kind = label_kinds[nuis_name]
        up_label = info["hist_up_label"]
        down_label = info["hist_down_label"]
        if up_label not in label_set:
            mismatches.append(f"{nuis_name}: Up label '{up_label}' missing from syst axis")
            continue
        if down_label not in label_set:
            mismatches.append(f"{nuis_name}: Down label '{down_label}' missing from syst axis")
            continue

        up_tail   = _parse_hist_label(up_label,   prefix)
        down_tail = _parse_hist_label(down_label, prefix)
        if kind == "abs":
            up_expected = info["up_value"]
            down_expected = info["down_value"]
        else:  # delta: label encodes (value - nominal)
            up_expected   = info["up_value"]   - info["nominal"]
            down_expected = info["down_value"] - info["nominal"]
        if abs(up_tail - up_expected) > tolerance:
            mismatches.append(
                f"{nuis_name}: hist Up value {up_tail} != PARAM_MAP {up_expected} (label '{up_label}')"
            )
        if abs(down_tail - down_expected) > tolerance:
            mismatches.append(
                f"{nuis_name}: hist Down value {down_tail} != PARAM_MAP {down_expected} (label '{down_label}')"
            )

    if mismatches:
        raise ValueError(
            "PARAM_MAP cross-check failed against corr_hist syst axis:\n  "
            + "\n  ".join(mismatches)
        )


# ---------------------------------------------------------------------------
# Rabbit-side Mapping and Regularizer (lazily constructed).
# ---------------------------------------------------------------------------

def _make_mapping_class():
    from rabbit.mappings.mapping import BaseMapping

    # Allowed physical-parameter names for kfactor overrides; must match the
    # "param" entries in PARAM_MAP above.
    _ALLOWED_KFACTOR_PARAMS = {info["param"] for info in PARAM_MAP.values()}

    class NPMonotonicityMapping(BaseMapping):
        """Vestigial BaseMapping that hangs onto indata, an optional
        criterion override, and optional per-nuisance kfactor overrides.

        Per-nuisance kfactors model the effect of rabbit `--scaleParams`:
        a kfactor k means a unit theta-pull moves the physical NP value by
        k * (Up_template - nominal) instead of the bare PARAM_MAP step.
        Without these the regularizer's view of the physical lambda becomes
        decoupled from the fit's (the wall stops firing).

        Invoked via:
            -r ...NPMonotonicityWall ...NPMonotonicityMapping
            -r ...NPMonotonicityWall ...NPMonotonicityMapping a
            -r ...NPMonotonicityWall ...NPMonotonicityMapping b lambda_2=5.24 lambda_4=2.24 Lambda_2=2.0 Delta_Lambda_2=2.0 Lambda_4=2.0
        Positional args after the criterion are <param_name>=<kfactor> tokens
        where <param_name> is one of: lambda_2, lambda_4, Lambda_2,
        Delta_Lambda_2, Lambda_4. Unspecified kfactors default to 1.
        """

        def __init__(self, indata, key, criterion=DEFAULT_CRITERION, kfactors=None):
            super().__init__(indata, key)
            self.indata = indata
            self.criterion = criterion
            self.kfactors = dict(kfactors) if kfactors else {}

        @classmethod
        def parse_args(cls, indata, *args):
            criterion = DEFAULT_CRITERION
            kfactors = {}
            for a in args:
                if "=" in a:
                    pname, kval = a.split("=", 1)
                    if pname not in _ALLOWED_KFACTOR_PARAMS:
                        raise ValueError(
                            f"NPMonotonicityMapping: unknown kfactor param '{pname}'; "
                            f"allowed: {sorted(_ALLOWED_KFACTOR_PARAMS)}"
                        )
                    try:
                        kfactors[pname] = float(kval)
                    except ValueError:
                        raise ValueError(
                            f"NPMonotonicityMapping: kfactor for '{pname}' not a float: '{kval}'"
                        )
                else:
                    if a not in ("a", "b"):
                        raise ValueError(
                            f"NPMonotonicityMapping: positional arg must be 'a', 'b', or "
                            f"'<param>=<kfactor>', got '{a}'"
                        )
                    criterion = a
            key_parts = [cls.__name__, criterion]
            for pname in sorted(kfactors):
                key_parts.append(f"{pname}={kfactors[pname]:g}")
            key = " ".join(key_parts)
            return cls(indata, key, criterion=criterion, kfactors=kfactors)

    return NPMonotonicityMapping


def _make_regularizer_class():
    from rabbit.regularization.regularizer import Regularizer
    import tensorflow as tf

    class NPMonotonicityWall(Regularizer):
        """Hinge-loss penalty enforcing tanh-monotonicity (or sign-preservation)
        of the SCETlib NP CS kernel and TMD boundary condition.

        Construct via:
            -r wremnants.postprocessing.np_monotonicity.NPMonotonicityWall \\
               wremnants.postprocessing.np_monotonicity.NPMonotonicityMapping [criterion]

        Wall hardness is controlled at fit time by rabbit's
        --regularizationStrength (the penalty is multiplied by exp(2*tau)
        inside fitter.py:2491).
        """

        # Pull criterion choice from the param map; default to (b) monotonicity.
        # Criterion (a): sign-preservation, factors sqrt(4) and sqrt(36).
        # Criterion (b): monotonicity,      factors sqrt(3) and sqrt(20).
        _K_CS = {"a": 4.0, "b": 3.0}
        _K_TMD = {"a": 36.0, "b": 20.0}

        def __init__(self, mapping, dtype):
            super().__init__(mapping, dtype)
            self.dtype = dtype
            self.mapping = mapping
            self.indata = mapping.indata

            self.criterion = getattr(mapping, "criterion", DEFAULT_CRITERION)
            if self.criterion not in self._K_CS:
                raise ValueError(
                    f"NPMonotonicityWall: criterion must be 'a' or 'b', got '{self.criterion}'"
                )
            self.k_cs = self._K_CS[self.criterion]
            self.k_tmd = self._K_TMD[self.criterion]

            # Locate each NP nuisance in indata.systs and stash linearization
            # coefficients for param(theta). The mapping's kfactors dict maps
            # physical-parameter names (e.g. "lambda_2") to a multiplicative
            # scale on the template half-range, matching the rabbit
            # --scaleParams used to build the workspace.
            kfactors = getattr(mapping, "kfactors", {}) or {}
            systs = [s.decode() if isinstance(s, (bytes, bytearray)) else s
                     for s in self.indata.systs]
            self._nuis_lookup = {}
            for nuis_name, info in PARAM_MAP.items():
                if nuis_name not in systs:
                    raise ValueError(
                        f"NPMonotonicityWall: expected nuisance '{nuis_name}' "
                        f"not found in indata.systs"
                    )
                pname = info["param"]
                k = float(kfactors.get(pname, 1.0))
                self._nuis_lookup[pname] = {
                    "syst_index": systs.index(nuis_name),
                    "nominal":  float(info["nominal"]),
                    "up_value": float(info["up_value"]),
                    "down_value": float(info["down_value"]),
                    "kfactor":  k,
                }

            # x = concat([poi, model_nui, theta]); theta_start resolved at first use.
            self._theta_start = None
            self._cast = lambda v: tf.constant(v, dtype=self.dtype)

        def set_expectations(self, initial_params, initial_observables, parms=None):
            # parms accepted for interface compatibility (rabbit passes the
            # current parameter names so regularizers can resolve by name).
            # DEPRECATED module -- if it is ever revived, resolve positions
            # from parms here rather than caching them at construction.
            nsyst = len(self.indata.systs)
            self._theta_start = int(initial_params.shape[0]) - nsyst

        def _physical_value(self, params, param_name):
            info = self._nuis_lookup[param_name]
            idx = info["syst_index"]
            theta = params[self._theta_start + idx]
            nominal = self._cast(info["nominal"])
            k = self._cast(info["kfactor"])
            delta_up = self._cast(info["up_value"] - info["nominal"]) * k
            delta_down = self._cast(info["nominal"] - info["down_value"]) * k
            theta_pos = tf.maximum(theta, self._cast(0.0))
            theta_neg = tf.maximum(-theta, self._cast(0.0))
            return nominal + theta_pos * delta_up - theta_neg * delta_down

        def compute_nll_penalty(self, params, observables):
            if self._theta_start is None:
                nsyst = len(self.indata.systs)
                self._theta_start = int(params.shape[0]) - nsyst

            zero = self._cast(0.0)

            # CS side
            l2 = self._physical_value(params, "lambda_2")
            l4 = self._physical_value(params, "lambda_4")
            l2_safe = tf.maximum(l2, zero)
            cs_floor = -tf.sqrt(self._cast(self.k_cs) * l2_safe * self._cast(LAMBDA_6))
            cs_pen = tf.square(tf.maximum(zero, cs_floor - l4))
            cs_l2_pos_pen = tf.square(tf.maximum(zero, -l2))  # small-bT positivity

            # TMD side: evaluate at y=0 and y=y_max, sum
            L2 = self._physical_value(params, "Lambda_2")
            DL2 = self._physical_value(params, "Delta_Lambda_2")
            L4 = self._physical_value(params, "Lambda_4")

            pens = [cs_pen, cs_l2_pos_pen]
            for y_sq in (0.0, Y_MAX * Y_MAX):
                L2y = L2 + DL2 * self._cast(y_sq)
                pens.append(tf.square(tf.maximum(zero, -L2y)))  # L_2(y) >= 0
                L2y_safe = tf.maximum(L2y, zero)
                c1 = self._cast(3.0) * L4 + L2y_safe ** 3
                floor = -tf.sqrt(self._cast(self.k_tmd) * self._cast(BIG_LAMBDA_6) * L2y_safe)
                pens.append(tf.square(tf.maximum(zero, floor - c1)))

            return tf.add_n(pens)

    return NPMonotonicityWall


# PEP-562 lazy class resolution: rabbit's loader does
#     module = importlib.import_module(...); cls = getattr(module, class_name)
# so we synthesize the classes on first attribute access. Keeps the module
# importable on the histmaker side without rabbit/TF.
def __getattr__(name):
    if name == "NPMonotonicityMapping":
        cls = _make_mapping_class()
        globals()["NPMonotonicityMapping"] = cls
        return cls
    if name == "NPMonotonicityWall":
        cls = _make_regularizer_class()
        globals()["NPMonotonicityWall"] = cls
        return cls
    raise AttributeError(name)
