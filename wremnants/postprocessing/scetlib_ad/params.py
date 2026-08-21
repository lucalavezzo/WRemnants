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

import math

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
    # Profile scales and matching transition points, registered by
    # set_diff_scales(1). scale_kappa_F is inert in the kernel -- the slot exists
    # only for build_pdf_variations to tie the muF member pair to -- so it does
    # nothing unless the cache was built with has_muf.
    "scale_kappa_R": "resumScaleMuR",
    "scale_kappa_F": "resumScaleMuF",
    "scale_x1": "resumTransition1",
    "scale_x2": "resumTransition2",
    "scale_x3": "resumTransition3",
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
PDF_PREFIX_IN = "pdf_eig"
PDF_PREFIX_OUT = "pdfEig"


def rabbit_name(scetlib_name):
    """SCETlib gradient-parameter name -> the name rabbit reports."""
    if scetlib_name in EXPLICIT_NAMES:
        return EXPLICIT_NAMES[scetlib_name]
    if scetlib_name.startswith(TNP_PREFIX_IN):
        return TNP_PREFIX_OUT + scetlib_name[len(TNP_PREFIX_IN) :]
    if scetlib_name.startswith(PDF_PREFIX_IN):
        return PDF_PREFIX_OUT + scetlib_name[len(PDF_PREFIX_IN) :]
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
    if rabbit.startswith(PDF_PREFIX_OUT):
        return PDF_PREFIX_IN + rabbit[len(PDF_PREFIX_OUT) :]
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
    # --- profile scales and transition points -------------------------------
    # These REPLACE card nuisances, so their priors have to reproduce the
    # variation the card encodes, and two of the three only do so approximately.
    #
    # resumScaleMuF: exact. The muF pair is built at kappa_F = 0.5 / 2.0 and
    #   interpolated in t = ln(kappa_F)/ln(muf_hi), so the direction is already
    #   unit-normalised and |theta| = 1 IS the card's variation.
    # resumScaleMuR: APPROXIMATE. The parameter is kappa_R itself, central 1 and
    #   linear, while the card varies kappaFO by x2 and /2. sigma = 0.5 gives
    #   +-1 sigma = [0.5, 1.5], so the up-variation is understated ([0.5, 2.0]).
    #   A log-parametrised kappa_R would fix this properly; until then this is a
    #   deliberate approximation, not an equivalence.
    # resumTransition2: APPROXIMATE. Central 0.6, the card's variations are
    #   0.35 and 0.75 (variations_resummed.conf [35]/[36]), i.e. -0.25/+0.15;
    #   0.2 is the symmetric stand-in.
    # resumTransition1/3 are FROZEN by default -- the analysis varies only the
    #   CENTRAL transition point ("new recommendation from Frank for variation of
    #   central transition parameter only"), so floating the outer two would ADD
    #   uncertainty the card does not carry.
    # resumScaleMuR / MuF / Transition2 are reparametrised unit nuisances --
    # see REPARAM below; prior_sigma() returns 1.0 for them.
    "resumTransition1": None,
    "resumTransition3": None,
}
TNP_PRIOR_SIGMA = 1.0


def prior_sigma(rabbit):
    """Default Gaussian prior sigma for a rabbit-facing name (None = free)."""
    if rabbit in REPARAM:
        # Unit nuisance by construction: the map carries the physical range, so
        # theta = +-1 IS the variation the replaced template encoded.
        return 1.0
    if rabbit.startswith(TNP_PREFIX_OUT):
        return TNP_PRIOR_SIGMA
    if rabbit.startswith(PDF_PREFIX_OUT):
        # Hessian eigenvector coefficients: c = +-1 IS the member 1 sigma, by
        # construction of build_pdf_variations (exact at c = 0, +-1).
        return 1.0
    return PRIOR_SIGMAS.get(rabbit, None)


# --- Reparametrisation: unit nuisances for the profile scales -----------------
#
# SCETlib registers the profile scales as the PHYSICAL quantities:
# ``scale_kappa_R`` and ``scale_kappa_F`` are kappa itself with central 1, and
# ``scale_x1..x3`` are the transition points themselves. That is the right
# interface for a calculation, but it is the wrong one for a nuisance, because
# the template variations these REPLACE are not symmetric in the physical
# variable:
#
#   kappaFO   x2 and /2          -> symmetric in ln(kappa), not in kappa
#   x2        0.6 -> 0.35, 0.75  -> genuinely asymmetric, -0.25 / +0.15
#
# and rabbit's ParamModel priors are a single symmetric Gaussian per parameter
# (fitter.py: cw = 1/sigma^2, one scalar, no up/down hook). Tuning sigma cannot
# reproduce either variation: sigma = 0.5 on a linear kappa_R gives [0.5, 1.5],
# understating the up side.
#
# So the fitted parameter is a UNIT nuisance theta and the model maps it to the
# physical value, exactly as SCETlib itself does for the PDF eigenvectors and
# the muF pair ("exact at 0, +-1, quadratic in between"). Every replaced-template
# direction is then sigma = 1, the same convention as the TNPs and pdfEig*.
#
#   "log"  : value = exp(theta * L)          theta = +-1 -> exp(+-L)
#   "quad" : value = c0 + c1*theta + c2*theta^2
#
# The log form has a second benefit: exp() is positive by construction, so it
# cannot trip SCETlib's silent `p[_muf_index] > 0. ? ... : 1.` fallback, which
# would drop the muF variation with no error at all.
LN2 = math.log(2.0)
REPARAM = {
    # kappa_R: theta = +-1 -> kappa_R = 2 / 0.5, matching kappaFO x2 and /2.
    "resumScaleMuR": ("log", (LN2,)),
    # kappa_F: same. SCETlib converts internally to t = ln(kappa_F)/ln(2), so
    # theta = +-1 lands exactly on the two members that were built (0.5, 2.0).
    "resumScaleMuF": ("log", (LN2,)),
    # x2: the quadratic through the three points the analysis actually uses --
    # theta = -1 -> 0.35, theta = 0 -> 0.6, theta = +1 -> 0.75. Monotone for
    # |theta| < 2 (the derivative 0.20 - 0.10*theta vanishes at theta = 2).
    "resumTransition2": ("quad", (0.6, 0.20, -0.05)),
    # resumTransition1/3 are deliberately NOT reparametrised: they are frozen by
    # default and no reference variation exists for them, so a study that floats
    # them should do so in the physical variable and choose its own range.
}


def reparam(rabbit):
    """``(kind, coeffs)`` for a reparametrised name, else ``None``."""
    return REPARAM.get(rabbit)


# --- Defaults ----------------------------------------------------------------
#
# Parameters frozen unless the user asks otherwise. These are shape constants of
# the SCETlib nonperturbative forms, not physics we fit: lambda_inf sets the
# saturation of the tanh forms and b0_over_bmax_nu the b* convention.
DEFAULT_FROZEN = (
    "lambda_inf",
    "lambda_inf_nu",
    "b0_over_bmax_nu",
    # Only the CENTRAL matching transition point is varied in the analysis; see
    # the PRIOR_SIGMAS comment. Float these two only as a deliberate study.
    "resumTransition1",
    "resumTransition3",
)

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
    # Named to line up with the card groups they replace, so a grouped-impact
    # bar stays comparable between the template and model treatments.
    "resumScale": ("resumScaleMuR", "resumScaleMuF"),
    "resumTransition": (
        "resumTransition1",
        "resumTransition2",
        "resumTransition3",
    ),
}


def tnp_group(names):
    """The ``resumTNP`` impact group for whichever TNPs are registered."""
    return tuple(n for n in names if n.startswith(TNP_PREFIX_OUT))


def pdf_group(names):
    """The ``pdf`` impact group for whichever eigenvector coefficients exist."""
    return tuple(n for n in names if n.startswith(PDF_PREFIX_OUT))


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
