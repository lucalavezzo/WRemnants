"""Parameter registry for the SCETlib autodiff param model.

SCETlib owns the parameter vector: :meth:`DrellYan.gradient_param_names` returns
the differentiable parameters in a FIXED order that is baked into the cached bin
rules (the rule fingerprint hashes the names in order, so any addition or
reordering invalidates a cache). This module is the translation layer between
that vector and the names rabbit sees.

Rabbit-facing names use the spelling the analysis tooling already reads
(``lambda2``, ``lambda2_nu``, …) so the postfit readers, the cross-run
fit-summary tools and the impact-group labels work unchanged. ``alphas`` becomes
``alphaS`` in PHYSICAL units (0.118-ish), NOT the ``pdfAlphaS`` template's
Delta(alpha_s) = 0.002-per-theta convention.
"""

# --- SCETlib gradient name -> rabbit-facing name -----------------------------
#
# Every entry here is an exact-match rename. Names not listed fall through
# :func:`rabbit_name`, which handles the two open-ended families (TNPs and, in
# a later phase, PDF eigenvector coefficients).
EXPLICIT_NAMES = {
    "alphas": "alphaS",
    "np_eff_lambda_inf": "lambda_inf",
    "np_eff_lambda2": "lambda2",
    "np_eff_lambda4": "lambda4",
    "np_eff_lambda6": "lambda6",
    "np_eff_delta_lambda2": "delta_lambda2",
    "np_gnu_lambda_inf": "lambda_inf_nu",
    "np_gnu_lambda2": "lambda2_nu",
    "np_gnu_lambda4": "lambda4_nu",
    "np_gnu_lambda6": "lambda6_nu",
    "np_gnu_b0_bmax": "b0_over_bmax_nu",
}

# TNPs: ``tnp_gamma_cusp`` -> ``resumTNP_gamma_cusp``. The prefix matches the
# group setupRabbit gives the discrete TNP templates (``resumTNP``), so a
# grouped-impact bar stays comparable between the template and model treatments.
TNP_PREFIX_IN = "tnp_"
TNP_PREFIX_OUT = "resumTNP_"

# PDF eigenvector coefficients, appended after the physics parameters when the
# cache carries PDF variations (phase 4). ``c_e`` are standard N(0,1) Hessian
# coefficients.
PDF_COEFF_FMT = "pdfEig{:d}"


def rabbit_name(scetlib_name):
    """SCETlib gradient-parameter name -> the name rabbit reports."""
    if scetlib_name in EXPLICIT_NAMES:
        return EXPLICIT_NAMES[scetlib_name]
    if scetlib_name.startswith(TNP_PREFIX_IN):
        return TNP_PREFIX_OUT + scetlib_name[len(TNP_PREFIX_IN) :]
    raise KeyError(
        f"scetlib_ad.params: no rabbit name for SCETlib parameter "
        f"{scetlib_name!r}. Add it to EXPLICIT_NAMES (and give it a prior in "
        f"PRIOR_SIGMAS / a group in IMPACT_GROUP_MEMBERS)."
    )


def scetlib_name(rabbit):
    """Inverse of :func:`rabbit_name` (raises on an unknown name)."""
    for k, v in EXPLICIT_NAMES.items():
        if v == rabbit:
            return k
    if rabbit.startswith(TNP_PREFIX_OUT):
        return TNP_PREFIX_IN + rabbit[len(TNP_PREFIX_OUT) :]
    raise KeyError(f"scetlib_ad.params: unknown rabbit parameter {rabbit!r}")


# --- Priors ------------------------------------------------------------------
#
# Only consulted when the model is constructed with ``priors=1``; otherwise every
# parameter floats free. ``None`` means "free even when priors are on".
#
# The lambda sigmas are the ones the analysis uses for the corresponding template
# nuisances, so the nonperturbative sector is constrained the same way whether it
# is fitted continuously or morphed. TNPs are genuine N(0,1) nuisances --
# theta is normalised upstream so |theta|=1 IS the recommended variation
# (prod/scetlib_run/examples/theory_nuisance_parameters/*.conf) -- and get
# sigma = 1 by default, unlike the free lambdas.
PRIOR_SIGMAS = {
    "alphaS": None,
    "lambda2": 0.50,
    "lambda4": 0.50,
    "lambda6": 0.10,
    "delta_lambda2": 0.50,
    "lambda_inf": None,
    "lambda2_nu": 0.10,
    "lambda4_nu": 0.50,
    "lambda6_nu": 0.10,
    "lambda_inf_nu": None,
    "b0_over_bmax_nu": None,
}
TNP_PRIOR_SIGMA = 1.0


def prior_sigma(rabbit):
    """Default Gaussian prior sigma for a rabbit-facing name (None = free)."""
    if rabbit.startswith(TNP_PREFIX_OUT):
        return TNP_PRIOR_SIGMA
    if rabbit.startswith("pdfEig"):
        return 1.0
    return PRIOR_SIGMAS.get(rabbit, None)


# --- Defaults ----------------------------------------------------------------
#
# Parameters frozen unless the user asks otherwise. These are shape constants of
# the SCETlib nonperturbative forms, not physics we fit: lambda_inf sets the
# saturation of the tanh forms and b0_over_bmax_nu the b* convention.
DEFAULT_FROZEN = ("lambda_inf", "lambda_inf_nu", "b0_over_bmax_nu")

# Grouped impacts over the model's own parameters (rabbit resolves these labels
# to floating x-indices; see Fitter._resolved_param_impact_groups). Membership is
# intersected with the parameters actually registered.
IMPACT_GROUP_MEMBERS = {
    "resumNonpert": (
        "lambda2",
        "lambda4",
        "lambda6",
        "delta_lambda2",
        "lambda_inf",
        "lambda2_nu",
        "lambda4_nu",
        "lambda6_nu",
        "lambda_inf_nu",
        "b0_over_bmax_nu",
    ),
    "scetlibNPFeff": (
        "lambda2",
        "lambda4",
        "lambda6",
        "delta_lambda2",
        "lambda_inf",
    ),
    "scetlibNPgammaNu": (
        "lambda2_nu",
        "lambda4_nu",
        "lambda6_nu",
        "lambda_inf_nu",
        "b0_over_bmax_nu",
    ),
}


def tnp_group(names):
    """The ``resumTNP`` impact group for whichever TNPs are registered."""
    return tuple(n for n in names if n.startswith(TNP_PREFIX_OUT))


# --- lambda_central cross-check ----------------------------------------------
#
# A histmaker output records the nonperturbative values its theory correction was
# generated at, under two sub-dicts using the histmaker's own spelling, and that
# is propagated into the datacard. Map those names onto rabbit-facing ones so the
# card's anchor can be compared against the cache's.
LAMBDA_CENTRAL_KEYS = {
    # card metadata key -> rabbit-facing name (identical here, but
    # spelled out so a future divergence is a one-line fix rather than a silent
    # mismatch)
    "lambda2": "lambda2",
    "lambda4": "lambda4",
    "lambda6": "lambda6",
    "delta_lambda2": "delta_lambda2",
    "lambda_inf": "lambda_inf",
    "lambda2_nu": "lambda2_nu",
    "lambda4_nu": "lambda4_nu",
    "lambda6_nu": "lambda6_nu",
    "lambda_inf_nu": "lambda_inf_nu",
    "b0_over_bmax_nu": "b0_over_bmax_nu",
}
