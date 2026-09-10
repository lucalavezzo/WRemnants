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
        f"{scetlib_name!r}. Add it to EXPLICIT_NAMES, decide whether it is FREE "
        f"or constrained (FREE_PARAMS), and give it a group in "
        f"IMPACT_GROUP_MEMBERS."
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
# NOTE the units. These sigmas are in RABBIT's parameter units, which for a
# REPARAM'd name is theta, not the physical variable. So a reparametrised
# parameter whose map already carries the physical scale wants sigma = 1.0 --
# |theta| = 1 IS 1 sigma, the same convention as the TNPs and pdfEig* -- and
# leaving the old PHYSICAL width here instead would silently tighten the prior
# by a factor of that width (lambda2 would go from 0.4 +- 0.5 to 0.4 +- 0.25).
# The physical widths the analysis chose now live in REPARAM as the map widths.
# Which parameters are FREE. Everything else is constrained at sigma = 1.
#
# That is the whole declaration, and it is deliberately a SET rather than a
# {name: sigma} table. A table invites the units bug: sigma is in RABBIT's
# units, which for a reparametrised name is theta, so a physical width left
# there gets multiplied by the map's width and the prior silently shrinks
# (lambda2 went from 0.4 +- 0.5 to 0.4 +- 0.25 exactly this way). With every
# constrained parameter at sigma = 1 there is no number here to get wrong: the
# physical 1 sigma is the REPARAM width, in one place, and
#
#     physical 1 sigma = width x sigma = width.
#
# A parameter with no REPARAM entry is already unit-normalised upstream -- the
# TNPs carry an N(0,1) constraint by construction, and pdfEig* have their CL
# convention in the coefficient map -- so sigma = 1 is right for them too, and
# they need no entry anywhere.
#
# Notes on the free ones:
#   alphaS            the POI.
#   lambda_inf, lambda_inf_nu, b0_over_bmax_nu
#                     shape constants of the NP form, frozen by default (see
#                     DEFAULT_FROZEN); free rather than constrained so a study
#                     that deliberately floats one is not fighting a prior it
#                     did not choose.
#   resumTransition1/3
#                     the analysis varies only the CENTRAL transition point
#                     ("new recommendation from Frank for variation of central
#                     transition parameter only"), so floating the outer two
#                     would ADD uncertainty the card does not carry. Frozen by
#                     default; no reference variation exists to normalise them
#                     against, which is also why they have no REPARAM map.
FREE_PARAMS = frozenset(
    {
        "alphaS",
        "lambda_inf",
        "lambda_inf_nu",
        "b0_over_bmax_nu",
        "resumTransition1",
        "resumTransition3",
    }
)


def prior_sigma(rabbit):
    """Gaussian prior sigma in RABBIT's units: None if free, else 1.0.

    There is no per-parameter width here by design -- see FREE_PARAMS. The
    physical width lives in REPARAM, so this answers only "is it constrained".
    """
    return None if rabbit in FREE_PARAMS else 1.0


# --- PDF confidence-level convention -----------------------------------------
#
# A Hessian PDF set's members are not all 1 sigma. CT18Z's are 90% CL
# displacements, and some sets additionally carry an analysis-chosen inflation.
# The analysis keeps both numbers in one place, ``theory_utils.pdfMap``, and
# composes them the same way everywhere:
#
#     scale = pdf_inflation_factor(pdfMap[set], noi) * pdfMap[set]["scale"]
#
# (``rabbit_theory_helper.add_pdf_uncertainty`` when it builds the templates;
# ``postfit_pdf_helper`` when it reads them back). For CT18Z with noi=alphaS
# that is 1.0 * 1/1.645 = 0.60790, and MEASURED off the 2D card's own logk it is
# 0.6025 +- 0.0055 over the 29 eigenvectors -- so the templates really do carry
# it and the model must too, or the swap inflates the PDF uncertainty by 1.645.
#
# WHERE it goes is the one place the model can do better than a template. The
# template route has no choice but to scale the RESPONSE, because a template is
# a fixed shape that can only morph linearly. The model evaluates SCETlib at
# whatever coefficient it is handed, and
#
#     I(c) = I_0 + c (I_+ - I_-)/2 + c^2 (I_+ + I_- - 2 I_0)/2
#
# is exactly quadratic in c (DrellYan.hpp, exact at c = 0, +-1). So the model
# scales the COEFFICIENT instead -- theta = +-1 evaluates the calculation at the
# 68% CL point in eigenvector space, which is what a 1 sigma PDF displacement
# physically is. The two agree to O(c) and differ by (scale^2 - scale) times the
# quadratic part, i.e. wherever the up and down members are not mirror images.
PDF_COEFF_SCALE_NOI_DEFAULT = ("alphaS",)


def pdf_set_key(lha_name, pdf_map=None):
    """``theory_utils.pdfMap`` key whose ``lha_name`` is *lha_name*.

    The cache runcard names the set the way LHAPDF does (``CT18ZNNLO``); the map
    is keyed by the analysis' short name (``ct18z``). Raises if the set is not
    in the map -- returning 1.0 for an unknown set would silently drop a
    convention that is a factor 1.645 for the one set we use.
    """
    if pdf_map is None:
        from wremnants.utilities import theory_utils

        pdf_map = theory_utils.pdfMap
    want = str(lha_name).strip().lower()
    for key, info in pdf_map.items():
        if str(info.get("lha_name", "")).lower() == want:
            return key
    raise KeyError(
        f"params.pdf_set_key: LHAPDF set {lha_name!r} is not in theory_utils."
        f"pdfMap, so its confidence-level convention is unknown. Add it there, "
        f"or pass pdf_coeff_scale=<float> explicitly."
    )


def pdf_coeff_scale(lha_name, noi=None, pdf_map=None):
    """Coefficient scale that makes ``pdfEig{i} = +-1`` the analysis' 1 sigma.

    Mirrors ``postfit_pdf_helper.PostfitPdfHelper`` (which does the same product
    the other way round) so the model and the templates cannot drift apart:
    the per-set 90%->68% ``scale`` times the nuisance-of-interest inflation.
    """
    from wremnants.utilities import theory_utils

    pdf_map = theory_utils.pdfMap if pdf_map is None else pdf_map
    info = pdf_map[pdf_set_key(lha_name, pdf_map)]
    noi = list(PDF_COEFF_SCALE_NOI_DEFAULT if noi is None else noi)
    return float(theory_utils.pdf_inflation_factor(info, noi)) * float(
        info.get("scale", 1)
    )


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
    #
    # --- "unit": value = <cache anchor> + width * theta -----------------------
    #
    # The remaining physical parameters. Before this, alphaS and the five NP
    # lambdas were the only ones rabbit saw in physical units, which is what
    # made them the odd ones out: 41 of 47 fitted parameters were already unit
    # nuisances. Normalising them keeps rabbit generic (it only ever sees
    # theta ~ O(1)) and collapses the curvature spread that made the
    # preconditioner necessary -- one block had max|diag| = 3.3e+09 against
    # singletons at exactly 1, and that nine-order range IS this units mismatch.
    #
    # c0 is NOT written here: the model fills it from the CACHE ANCHOR, so
    # theta = 0 reproduces the anchor by construction and cannot drift when a
    # cache is built at a different tune. Only the width is a choice, and it is
    # the parameter's own natural step:
    #   alphaS        the PDF set's alphasRange (0.002 for CT18Z) -- exactly the
    #                 Delta(alpha_s)-per-theta convention the pdfAlphaS template
    #                 used, so theta here means what it meant before.
    #   the lambdas   their PRIOR_SIGMAS widths, so |theta| = 1 is 1 sigma of the
    #                 prior the analysis already chose.
    "alphaS": ("unit", (0.002,)),
    "lambda2": ("unit", (0.50,)),
    "lambda4": ("unit", (0.50,)),
    "delta_lambda2": ("unit", (0.50,)),
    "lambda2_nu": ("unit", (0.10,)),
    # lambda6 / lambda6_nu: the tanh_6 form's third coefficient. Only active
    # when np_model / np_model_nu selects tanh_6, but reparametrised anyway so
    # that EVERY constrained parameter obeys the same rule -- sigma = 1 in
    # theta, physical width in the map. Leaving one physical is how a prior gets
    # silently rescaled the next time someone adds a map to it.
    "lambda6": ("unit", (0.10,)),
    "lambda6_nu": ("unit", (0.10,)),
    "lambda4_nu": ("unit", (0.50,)),
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
    # Inert for the Z, not a choice: b_qqDS scales a channel that does not
    # contribute, so its whole response is O(1e-16) and its Jacobian column is
    # identically zero. A zero column is a zero row+column of the NLL Hessian,
    # i.e. a singular covariance -- `_check_no_inert_params` refuses it, so this
    # name has to be frozen rather than "fitted with a prior". Measured
    # 2026-08-27 (studies/.../260827-authoritative-validation, 06_model_defaults).
    "resumTNP_b_qqDS",
    # Only the CENTRAL matching transition point is varied in the analysis; see
    # the PRIOR_SIGMAS comment. Float these two only as a deliberate study.
    "resumTransition1",
    "resumTransition3",
)

# resumTransition2 was frozen here from 2026-08-21 to 2026-08-25 because the
# derivative SCETlib handed us for it was sign-inverted -- moving the transition
# points moves muF while the per-node beam convolutions stayed at the config's
# muF. UNFROZEN: scetlib-cms bfc6be6 feeds the induced per-node muF shift into
# the muF member interpolation, and validate_variations now gives 1.1e-03 ..
# 3.4e-03 against the production templates where it was 1.1e-01 .. 1.99e-01.
#
# Two caveats that survive the fix:
#   * the residual is the interpolation's own limit, not a bug -- the induced
#     shift is carried by a quadratic through three knots (kappa_F = 0.5/1/2),
#     exact only AT the knots, and a transition variation lands between them.
#     Upstream quotes +7.8e-04 against an independent runcard route.
#   * the fix works THROUGH the muF member pair, so it does nothing without one.
#     Upstream guards this rather than leaving it silent: moving scale_x1..x3 off
#     the anchor on a cache with no muF pair RAISES ("the result would be WRONG,
#     sign included"), and says which call builds the pair. So a --no-pdf cache
#     cannot float the transition points at all -- which is also why
#     backend_check fails on one.

# Directions whose response has been MEASURED to disagree with the template it
# replaces, keyed by rabbit-facing name -> why. These are not frozen for physics
# reasons, so anyone who floats one anyway deserves to be told once, loudly,
# rather than to find out from a pull. Keep the strings short; the detail lives in
# studies/scetlib-ad-param-model/.
# Empty since 2026-08-25: every direction the model registers now agrees with its
# production template. The worst of 39 is 1.4e-02 (mufup, qT [0,1], which is
# template precision rather than a model error -- model/direct = 1.0000 there),
# and 37 of them are <= 7.5e-03. Kept as the hook, because "this direction's
# response is known-wrong, do not profile it" is worth being able to say loudly
# rather than in a comment.
KNOWN_BAD_RESPONSE = {}


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
