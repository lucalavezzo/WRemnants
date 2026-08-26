"""Tanh-monotonicity wall for the FranksVals NP model.

FranksVals tanh_2 functional form (no b_T^6 / b_T^5 tails, lambda_4_nu
and lambda_inf_nu fixed on the CS side, Lambda_inf fixed on the TMD
side):

  CS:  gamma_NP(b_T)   = -(lambda_inf_nu/2) * tanh(A/lambda_inf_nu)
       A(b_T)          = lambda_2_nu * b_T^2     (lambda_4_nu = 0 fixed)
  TMD: f_NP(b_T,y)     = exp[-2 * Lambda_inf * b_T * tanh(B/Lambda_inf)]
       B(b_T,y)        = L_2(y) * b_T + c_1(y)/3 * b_T^3
       L_2(y)          = Lambda_2 + Delta_Lambda_2 * y^2
       c_1(y)          = 3 * Lambda_4 + L_2(y)^3

With lambda_6_CS = Lambda_6_TMD = 0 (no tail terms), the
criterion (a) sign-preservation and criterion (b) monotonicity walls
collapse to the same simple set:

  CS:   A(b_T) >= 0  for all b_T  <=>  lambda_2_nu >= 0
                                       (lambda_4_nu = 0 fixed)
  TMD:  B(b_T,y) >= 0 monotone in b_T <=>
            L_2(y) >= 0  AND  c_1(y) >= 0   at every y

Evaluate the TMD walls at y in {0, Y_MAX} and sum (safe under floating
Delta_Lambda_2 of either sign).

Floating nuisances (4):
  scetlibNPgammaLambda2                 -> lambda_2_nu
  chargeVgenNP0scetlibNPZlambda2        -> Lambda_2
  chargeVgenNP0scetlibNPZdelta_lambda2  -> Delta_Lambda_2
  chargeVgenNP0scetlibNPZlambda4        -> Lambda_4

Invoke via:
  rabbit_fit.py <indata>.hdf5 ...
    -r wremnants.postprocessing.np_monotonicity_franks.NPMonotonicityFranksWall \\
       wremnants.postprocessing.np_monotonicity_franks.NPMonotonicityFranksMapping \\
    --regularizationStrength <tau> \\
    --noConstrainParams 'scetlibNPgamma.*|scetlibNPZ.*'

Per-nuisance kfactor overrides (use when the workspace was built with
--scaleParams) follow the same protocol as np_monotonicity.py:
positional `<param_name>=<k>` tokens after the mapping class name. Valid
param names: lambda_2_nu, Lambda_2, Delta_Lambda_2, Lambda_4.
"""

import numpy as np


# FranksVals fixed (non-floating) values; not part of PARAM_MAP because
# they have no nuisance.  Kept here as documentation.
LAMBDA_INF_NU_FIXED = 2.0   # CS lambda_inf_nu
LAMBDA_4_NU_FIXED   = 0.0   # CS lambda_4_nu
LAMBDA_INF_FIXED    = 1.0   # TMD Lambda_inf

# Operational |y_ll| range for the TMD binding-y evaluation.
Y_MAX = 2.5

# Per-nuisance template map.  Centrals/Up/Down match np_param_map_franks.json
# (single source of truth for the plot_np_kernel_franks.py + this module).
PARAM_MAP = {
    # CS side (only lambda_2_nu varies in FranksVals).
    "scetlibNPgammaLambda2": {
        "param": "lambda_2_nu",
        "nominal": 0.15, "up_value": 0.25, "down_value": 0.05,
        "hist_up_label": "lambda2_nu0.25",
        "hist_down_label": "lambda2_nu0.05",
    },
    # TMD side (3 nuisances).
    "chargeVgenNP0scetlibNPZlambda2": {
        "param": "Lambda_2",
        "nominal": 0.40, "up_value": 1.00, "down_value": 0.00,
        "hist_up_label": "lambda21.0",
        "hist_down_label": "lambda20.0",
    },
    "chargeVgenNP0scetlibNPZdelta_lambda2": {
        "param": "Delta_Lambda_2",
        "nominal": 0.00, "up_value": 0.02, "down_value": -0.02,
        "hist_up_label": "delta_lambda20.02",
        "hist_down_label": "delta_lambda2-0.02",
    },
    "chargeVgenNP0scetlibNPZlambda4": {
        "param": "Lambda_4",
        "nominal": 0.40, "up_value": 1.00, "down_value": 0.00,
        "hist_up_label": "lambda41.0",
        "hist_down_label": "lambda40.0",
    },
}


def _make_mapping_class():
    from rabbit.mappings.mapping import BaseMapping

    _ALLOWED = {info["param"] for info in PARAM_MAP.values()}

    class NPMonotonicityFranksMapping(BaseMapping):
        """Holds optional per-nuisance kfactor overrides for the FranksVals
        wall.  No 'criterion' switch because the (a) and (b) walls collapse
        to the same set of constraints when the b_T^6/b_T^5 tails are zero.
        """

        def __init__(self, indata, key, kfactors=None):
            super().__init__(indata, key)
            self.indata = indata
            self.kfactors = dict(kfactors) if kfactors else {}

        @classmethod
        def parse_args(cls, indata, *args):
            kfactors = {}
            for a in args:
                if "=" not in a:
                    raise ValueError(
                        f"NPMonotonicityFranksMapping: positional arg must be "
                        f"'<param>=<kfactor>', got '{a}'"
                    )
                pname, kval = a.split("=", 1)
                if pname not in _ALLOWED:
                    raise ValueError(
                        f"NPMonotonicityFranksMapping: unknown kfactor param "
                        f"'{pname}'; allowed: {sorted(_ALLOWED)}"
                    )
                try:
                    kfactors[pname] = float(kval)
                except ValueError:
                    raise ValueError(
                        f"NPMonotonicityFranksMapping: kfactor for '{pname}' "
                        f"not a float: '{kval}'"
                    )
            key_parts = [cls.__name__]
            for pname in sorted(kfactors):
                key_parts.append(f"{pname}={kfactors[pname]:g}")
            key = " ".join(key_parts)
            return cls(indata, key, kfactors=kfactors)

    return NPMonotonicityFranksMapping


def _make_regularizer_class():
    from rabbit.regularization.regularizer import Regularizer
    import tensorflow as tf

    class NPMonotonicityFranksWall(Regularizer):
        """Hinge-loss penalty enforcing the FranksVals tanh_2 walls.

        Penalty = sum of:
          (1) [max(0, -lambda_2_nu)]^2
          (2) [max(0, -L_2(y))]^2     at y = 0 and y = Y_MAX
          (3) [max(0, -c_1(y))]^2     at y = 0 and y = Y_MAX

        Strength controlled at fit time by --regularizationStrength tau
        (rabbit multiplies by exp(2 tau) inside fitter.py).
        """

        def __init__(self, mapping, dtype):
            super().__init__(mapping, dtype)
            self.dtype = dtype
            self.mapping = mapping
            self.indata = mapping.indata

            kfactors = getattr(mapping, "kfactors", {}) or {}
            systs = [s.decode() if isinstance(s, (bytes, bytearray)) else s
                     for s in self.indata.systs]
            self._nuis_lookup = {}
            for nuis_name, info in PARAM_MAP.items():
                if nuis_name not in systs:
                    raise ValueError(
                        f"NPMonotonicityFranksWall: expected nuisance "
                        f"'{nuis_name}' not found in indata.systs"
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

            # CS wall: lambda_2_nu >= 0.
            l2n = self._physical_value(params, "lambda_2_nu")
            pens = [tf.square(tf.maximum(zero, -l2n))]

            # TMD walls at y = 0 and y = Y_MAX.
            L2  = self._physical_value(params, "Lambda_2")
            DL2 = self._physical_value(params, "Delta_Lambda_2")
            L4  = self._physical_value(params, "Lambda_4")
            for y_sq in (0.0, Y_MAX * Y_MAX):
                L2y = L2 + DL2 * self._cast(y_sq)
                pens.append(tf.square(tf.maximum(zero, -L2y)))  # L_2(y) >= 0
                # c_1(y) = 3*Lambda_4 + L_2(y)^3 >= 0
                L2y_safe = tf.maximum(L2y, zero)
                c1 = self._cast(3.0) * L4 + L2y_safe ** 3
                pens.append(tf.square(tf.maximum(zero, -c1)))

            return tf.add_n(pens)

    return NPMonotonicityFranksWall


# PEP-562 lazy class resolution so the module is importable on the
# histmaker side without rabbit/TF available.
def __getattr__(name):
    if name == "NPMonotonicityFranksMapping":
        cls = _make_mapping_class()
        globals()["NPMonotonicityFranksMapping"] = cls
        return cls
    if name == "NPMonotonicityFranksWall":
        cls = _make_regularizer_class()
        globals()["NPMonotonicityFranksWall"] = cls
        return cls
    raise AttributeError(name)
