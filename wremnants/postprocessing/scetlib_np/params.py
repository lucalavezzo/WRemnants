"""Shared TF-free vocabulary and helpers for the SCETlib-NP param model.

Single source of truth for the λ parameter names, the np_model selector keys,
the reco/gen axis names, and the few numpy helpers used by both the TF core
(:mod:`sigma_gen`, :mod:`param_model`) and the TF-free tools
(:mod:`sigma_gen_at_lambda`, :mod:`np_function_plots`, :mod:`fitresult_lambdas`,
:mod:`lambda_central`). Import-light (numpy only) so the lightweight tools reach
these without pulling in the TF core.

Also home to :class:`NPTune`, the bundle every CLI passes around: the two form
names TOGETHER with the λ those forms use. Resolving a tune from its SOURCES (a
fitresults' postfit λ, a card's λ_central) needs hdf5/pickle readers and so lives
in :mod:`np_tune`, which keeps this module import-light.
"""

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

# The λ vocabulary (GNU_PARAMS / EFF_PARAMS / ALL_PARAMS) is DERIVED from the
# model registry further down — see :func:`_derive_params`. Adding a model that
# uses a new λ needs no edit here.

# np_model selector keys carried alongside the numeric λ in a tune dict.
EFF_MODEL_KEY = "np_model"
GNU_MODEL_KEY = "np_model_nu"

# Fit observable axis names.
RECO_AXES = ("ptll", "yll", "cosThetaStarll_quantile", "phiStarll_quantile")
GEN_AXES = ("ptVGen", "absYVGen")

# np_model selector aliases -> canonical model name. The form-factor branches in
# btgrid_tf and the registry below key on the canonical name.
_EFF_MODEL_ALIASES = {"hyp_tangent": "tanh_2", "square_root": "frac_2"}
_GNU_MODEL_ALIASES = {"hyp_tangent": "tanh_2", "linear": "frac_1"}

# ---- Model → λ registry (single source of truth) -----------------------------
# Per model: the λ it uses and their fit defaults (value = neutral start fallback,
# sigma = default prior width, None = free). Must stay in sync with the btgrid_tf
# form branches, which read exactly these λ by name.
EFF_MODEL_PARAMS = {
    "identity": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
    },
    "signed_lambda": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
    },
    "tanh_2": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
    # Structurally-damping TMD form: tanh_2 with BOTH coefficients of its tanh
    # argument pushed onto (0, ∞), so F_eff ≤ 1 is an IDENTITY for every λ and the
    # wall's whole TMD block becomes unnecessary. What gets mapped is NOT the raw
    # λ — tanh_2's damping condition is λ2_Y ≥ 0 and B = λ4 + λ2_Y³/(3λ∞²) ≥ 0,
    # so mapping λ4 positive would be STRICTER than physical and mapping δλ2
    # positive would be plain wrong (δλ2 < 0 is what the external tunes have):
    #   * λ2_Y = λ2 + δλ2·Y²  via ``btgrid_tf.pos_anchor_lambda2_Y_tf`` — its two
    #     ANCHORS L₂(0), L₂(pos_anchor_Y) mapped, then re-interpolated in Y², so
    #     λ2_Y stays LINEAR in Y² and the effective form is still a real tanh_2;
    #   * B via ``btgrid_tf.pos_floor_tf`` — the exact wall condition, as an
    #     identity. λ4 stays a free, sign-free fit parameter underneath it.
    # No λ∞ condition is needed (λ·tanh(c/λ) is even in λ), hence no λ∞ floor and
    # no boundary anywhere. Reduces to plain tanh_2 wherever the tune is already
    # physical, so tanh_2 card closure / λ-response validation carries over.
    # ``pos_floor`` and ``pos_anchor_Y`` are SHAPE CONSTANTS: freeze them.
    # SCETlib cannot produce this form, so it is an evaluation (numerator) form
    # only — never a card/denominator form.
    "tanh_2_pos": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
        # NOT 0: pos_floor_tf at floor 0 degenerates to relu, whose gradient is
        # EXACTLY 0 on the whole negative side — the dead-gradient failure this
        # form exists to avoid. 3e-3 against coefficients of O(0.1–0.5) costs
        # ≤ 0.1% in the identity region (error is O(floor²/x)).
        "pos_floor": {"value": 3e-3, "sigma": None},
        # The |Y| out to which λ2_Y ≥ 0 is guaranteed; matches np_damping_wall's
        # Y_MAX (the |Y| acceptance, and the range the external tune mapping was
        # fitted on). Beyond it λ2_Y is HELD, not extrapolated.
        "pos_anchor_Y": {"value": 2.5, "sigma": None},
    },
    "tanh_4": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
    "tanh_6": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
        "lambda6": {"value": 0.0, "sigma": 0.10},
    },
    # tanh_6 times a large-b_T turn-off S(b_T) = 1/(1+exp((b_T-bT_cutoff)/width)),
    # driving F_eff to 0 above b_T ≈ bT_cutoff whatever the λ do — an unphysical
    # λ4 < 0 tune otherwise makes the bare form turn around and blow up there.
    # bT_cutoff -> ∞ at FIXED width reproduces tanh_6 bit-exactly (both -> ∞ does
    # not: only their ratio matters). The two shape params are constants, not
    # physics λ: freeze them in a fit. Implementation: btgrid_tf.sigmoid_cutoff_tf.
    # SCETlib cannot produce this form, so it is an evaluation (numerator) form
    # only — never a card/denominator form.
    "tanh_6_sigmoid": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
        "lambda6": {"value": 0.0, "sigma": 0.10},
        "bT_cutoff": {"value": 2.0, "sigma": None},
        "bT_cutoff_width": {"value": 0.2, "sigma": None},
    },
    # tanh_6 with the DAMPING FOLD: the damping exponent 2·λ∞·b_T·tanh(a) is
    # replaced by |exponent + m| − m, so
    #     F_eff ≤ 1 + ``abs_margin``   for ANY λ, uniformly in b_T
    # — the ``F_eff ≤ 1.2``-style magnitude cap, as an identity rather than a
    # penalty. ``abs_margin`` is the allowed FRACTIONAL EXCURSION (0.2 ⇒ 1.2); it
    # enters the exponent-space fold as ln(1+margin). No λ sign anti-damps beyond
    # it, so NO damping wall is needed (not even a λ∞ floor: λ∞ < 0 folds too). At
    # the default 0 the fold is the plain absolute value and F_eff ≤ 1, equivalently
    # tanh|a| — b_T, λ∞ > 0 make the two identical there.
    #   Reduction: wherever the cap already holds pointwise — i.e. on the whole
    # physical (damping) region — the fold is the IDENTITY, so this form is tanh_6
    # EXACTLY, and a card closure / λ-response validated under tanh_6 carries
    # over unchanged (unlike tanh_6_sigmoid, which only reduces as bT_cutoff→∞).
    #   ``abs_margin`` is a SHAPE CONSTANT, not a physics λ (:data:`CONST_PARAMS`):
    # growing it relaxes the fold (m → ∞ is unfolded tanh_6), which the fit will
    # happily do, so the param model REFUSES to run with it floating — freeze it.
    # SCETlib cannot produce this form, so it is an evaluation (numerator) form
    # only — never a card/denominator form. Implementation: btgrid_tf.abs_fold_tf.
    "tanh_6_abs": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
        "lambda6": {"value": 0.0, "sigma": 0.10},
        "abs_margin": {"value": 0.0, "sigma": None},
    },
    "frac_2": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
    "frac_4": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
    "exp_2": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
    "exp_4": {
        "lambda2": {"value": 0.0, "sigma": 0.50},
        "lambda4": {"value": 0.0, "sigma": 0.50},
        "delta_lambda2": {"value": 0.0, "sigma": 0.20},
        "lambda_inf": {"value": 0.0, "sigma": None},
    },
}
GNU_MODEL_PARAMS = {
    "tanh_1": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
    "tanh_2": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
    # Structurally-damping CS form, the γ_ν partner of the tanh_2_pos F_eff form.
    # On THIS side the wall's tanh_2 condition really is per-λ positivity
    # (λ2_ν ≥ 0, λ4_ν ≥ 0), so running both through ``btgrid_tf.pos_floor_tf``
    # reproduces its feasible set EXACTLY — γ_ν ≤ 0 and monotone become
    # identities, the CS wall block is unnecessary (this form is in
    # ``np_damping_wall.UNCONDITIONAL_FORMS``), and — unlike tanh_6_pos — there is
    # NO precondition on any other λ, because tanh_2 has no b⁶ term to reopen the
    # question. λ∞_ν is unconditioned too (λ·tanh(P/λ) is even in λ).
    # ``pos_floor_nu`` is a SHAPE CONSTANT: freeze it. Numerator-form only.
    "tanh_2_pos": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
        # NOT 0 — see the tanh_2_pos ``pos_floor`` note above (floor 0 = relu =
        # dead gradient). NB tanh_6_pos's registry default IS 0.0 and predates
        # this; that entry is left alone deliberately, but a tanh_6_pos run should
        # set pos_floor_nu > 0 too.
        "pos_floor_nu": {"value": 3e-3, "sigma": None},
    },
    "tanh_6": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
        "lambda6_nu": {"value": 0.0, "sigma": 0.10},
    },
    # CS-side damping fold, the γ_ν analogue of the tanh_6_abs F_eff form above:
    # γ_ν^NP = −|λ∞_ν·tanh(a) + m_ν| + m_ν, so γ_ν ≤ ``abs_margin_nu`` (≤ 0 at the
    # default 0) for ANY λ — the ``γ_ν ≤ 0.2``-style cap. Additive, so unlike
    # F_eff's 1 + margin the value IS the fold offset. tanh_6 EXACTLY wherever the
    # cap already holds. Shape constant, must be frozen; numerator-form only.
    "tanh_6_abs": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
        "lambda6_nu": {"value": 0.0, "sigma": 0.10},
        "abs_margin_nu": {"value": 0.0, "sigma": None},
    },
    # Structurally-damping CS form: tanh_6 with lambda2_nu passed through a
    # smooth, MONOTONE, strictly-positive map (``pos_floor_tf``). Intended to be
    # run with lambda4_nu FROZEN AT 0 and lambda6_nu frozen >= 0 -- then every
    # coefficient of the tanh argument is non-negative, so gamma_nu <= 0 AND
    # monotone in b_T are IDENTITIES for any lambda2_nu, and the wall's CS block
    # is unnecessary (this form is in np_damping_wall.UNCONDITIONAL_FORMS).
    #
    # Unlike ``tanh_6_abs`` the map is monotone, not even, so there is NO mirror
    # degeneracy -- which is what made the *_abs fold unusable at margin 0.
    # ``pos_floor_nu`` is the shape constant setting how close to 0 the effective
    # coefficient can get; it must be frozen (see CONST_PARAMS).
    "tanh_6_pos": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
        "lambda6_nu": {"value": 0.0, "sigma": 0.10},
        "pos_floor_nu": {"value": 0.0, "sigma": None},
    },
    "frac_1": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
    "frac_2": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
    "exp_1": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
    "exp_2": {
        "lambda2_nu": {"value": 0.0, "sigma": 0.10},
        "lambda4_nu": {"value": 0.0, "sigma": 0.50},
        "lambda_inf_nu": {"value": 0.0, "sigma": None},
    },
}

# ---- λ vocabulary, DERIVED from the registries above ------------------------
# The registry is the single source of truth for which λ exist: a model that uses
# a new λ makes that λ appear in EFF_PARAMS / ALL_PARAMS automatically, with no
# second list to keep in sync (and no way for the two to disagree).
#
# The hints below fix only the ORDER, which is not free: ALL_PARAMS drives the
# fit's parameter ordering (``param_model._param_order`` -> the ``parms`` axis and
# the covariance layout) and every printed λ table / plot inset. A λ the registry
# has and a hint does not is appended (sorted), so forgetting to extend a hint is
# cosmetic, never a bug.
_GNU_ORDER = (
    "lambda2_nu",
    "lambda4_nu",
    "lambda6_nu",
    "lambda_inf_nu",
    "abs_margin_nu",
    "pos_floor_nu",
)
_EFF_ORDER = (
    "lambda2",
    "lambda4",
    "lambda6",
    "delta_lambda2",
    "lambda_inf",
    "abs_margin",
    "pos_floor",
    "pos_anchor_Y",
    "bT_cutoff",
    "bT_cutoff_width",
)


def _derive_params(registry, order_hint):
    """λ names the ``registry`` models actually use, ordered by ``order_hint``."""
    used = {p for model in registry.values() for p in model}
    named = [p for p in order_hint if p in used]
    return tuple(named) + tuple(sorted(used.difference(named)))


GNU_PARAMS = _derive_params(GNU_MODEL_PARAMS, _GNU_ORDER)
EFF_PARAMS = _derive_params(EFF_MODEL_PARAMS, _EFF_ORDER)
ALL_PARAMS = GNU_PARAMS + EFF_PARAMS

# Valid model names (canonical + aliases). The single validation set for both
# btgrid_tf (``np_model not in EFF_MODELS``) and the param model. Re-exported by
# btgrid_tf so it need not keep its own copy.
EFF_MODELS = frozenset(EFF_MODEL_PARAMS) | frozenset(_EFF_MODEL_ALIASES)
GNU_MODELS = frozenset(GNU_MODEL_PARAMS) | frozenset(_GNU_MODEL_ALIASES)

# Registry entries that are SHAPE CONSTANTS of a functional form, not physics λ:
# they say WHICH function the form is, so a fit that floats them changes the model
# it is fitting rather than measuring anything. Every one of them relaxes its form
# monotonically (abs_margin → ∞ unfolds tanh_6_abs; bT_cutoff → ∞ removes the
# tanh_6_sigmoid turn-off), so a minimiser given the freedom takes it. They stay in
# the registry — they ARE arguments of the form factor, and travel with the form in
# an :class:`NPTune` like any other λ — but the param model requires them frozen
# (see ``SCETlibNPParamModel._check_shape_constants``).
# NOTE the *_pos floors/anchor are here for a slightly different reason than the
# rest: they do NOT relax their form monotonically (a bigger pos_floor is MORE
# damping, a bigger pos_anchor_Y guarantees positivity further out), so a
# minimiser would not simply run them away. They must still be frozen — pos_floor
# is exactly degenerate with the coefficient it maps in the region where the map
# is the identity, and floating either one means fitting a different model than
# the one named.
CONST_PARAMS = frozenset(
    {
        "abs_margin",
        "abs_margin_nu",
        "bT_cutoff",
        "bT_cutoff_width",
        "pos_floor_nu",
        "pos_floor",
        "pos_anchor_Y",
    }
)


def param_defaults(np_model=None, np_model_nu=None):
    """``{name: {"value", "sigma"}}`` for the λ the given model(s) use.

    Union of the F_eff (``np_model``) and γ_ν (``np_model_nu``) registry rows,
    aliases resolved. Raises ``KeyError`` on an unknown model name."""
    out = {}
    if np_model is not None:
        out.update(EFF_MODEL_PARAMS[_EFF_MODEL_ALIASES.get(np_model, np_model)])
    if np_model_nu is not None:
        out.update(GNU_MODEL_PARAMS[_GNU_MODEL_ALIASES.get(np_model_nu, np_model_nu)])
    return out


def active_params(np_model=None, np_model_nu=None):
    """Names of the λ the chosen NP model(s) actually use (registry keys).

    Pass ``np_model`` for the F_eff (TMD) set, ``np_model_nu`` for the γ_ν (CS)
    set, or both for the union. A λ outside the returned set is inert for that
    model and is NOT a fit parameter. Source of truth: :data:`EFF_MODEL_PARAMS` /
    :data:`GNU_MODEL_PARAMS`, audited against the ``btgrid_tf`` form branches.

    Includes the form's shape constants (:data:`CONST_PARAMS`): they are genuine
    arguments of the form factor, so every EVALUATION needs them; it is only the
    FIT that must hold them fixed (:func:`const_params` names those)."""
    return set(param_defaults(np_model=np_model, np_model_nu=np_model_nu))


def const_params(np_model=None, np_model_nu=None):
    """The chosen model(s)' SHAPE CONSTANTS — :data:`CONST_PARAMS` ∩ active.

    The subset of :func:`active_params` a fit must FREEZE rather than float; see
    :data:`CONST_PARAMS` for why."""
    return active_params(np_model=np_model, np_model_nu=np_model_nu) & CONST_PARAMS


# ---- NPTune: a form pair and its λ, as one inseparable object ----------------


@dataclass(frozen=True)
class NPTune:
    """A complete, self-consistent NP tune: both form names + every λ they use.

    The bundle the scetlib_np CLIs pass around, so a form string and its λ can
    never travel separately. That separation was a real bug: ``--np-model`` used
    to stamp a form name onto a λ dict sourced elsewhere, producing a "tune"
    declaring ``tanh_6_sigmoid`` over λ with no ``bT_cutoff`` — which only
    surfaced as a ``KeyError`` deep inside the model constructor.

    Build with :meth:`create` (never the bare constructor): it validates both
    form names, fills λ the forms need but the caller did not supply from the
    registry defaults, and REJECTS λ the forms ignore rather than silently
    dropping them. ``values`` then holds exactly
    ``active_params(np_model, np_model_nu)`` — no more, no less — so
    :attr:`eff` / :attr:`gnu` are always valid ``btgrid_tf`` kwargs.

    ``filled`` records which λ came from the registry rather than the caller, so
    a CLI can say so out loud instead of inventing values quietly.
    """

    np_model: str
    np_model_nu: str
    values: Mapping
    filled: tuple = field(default=(), compare=False)

    @classmethod
    def create(
        cls,
        np_model,
        np_model_nu,
        values=None,
        *,
        base=None,
        what="tune",
        strict=True,
    ):
        """Validate and complete a tune.

        λ precedence, lowest first: registry default -> ``base`` -> ``values``.
        ``base`` is a fuller λ mapping to inherit from (e.g. a card λ_central);
        only its active-λ entries are used. ``what`` names the source in error
        messages (e.g. ``"--den-lambdas"``). With ``strict`` (the default) a λ in
        ``values`` that the chosen forms ignore is a hard error — setting a λ the
        form factor never reads is a silent no-op otherwise.
        """
        if np_model not in EFF_MODELS:
            raise ValueError(
                f"{what}: unknown np_model {np_model!r}; known: "
                f"{', '.join(sorted(EFF_MODELS))}"
            )
        if np_model_nu not in GNU_MODELS:
            raise ValueError(
                f"{what}: unknown np_model_nu {np_model_nu!r}; known: "
                f"{', '.join(sorted(GNU_MODELS))}"
            )
        pdefs = param_defaults(np_model=np_model, np_model_nu=np_model_nu)
        active = set(pdefs)
        values = dict(values or {})
        if strict:
            inert = [k for k in values if k not in active]
            if inert:
                raise ValueError(
                    f"{what}: "
                    + ", ".join(sorted(inert))
                    + f" not used by np_model={np_model} / np_model_nu="
                    + f"{np_model_nu} (active: "
                    + ", ".join(k for k in ALL_PARAMS if k in active)
                    + ")"
                )
        base = {k: v for k, v in dict(base or {}).items() if k in active}
        out, filled = {}, []
        for p in (k for k in ALL_PARAMS if k in active):
            if p in values:
                out[p] = float(values[p])
            elif p in base:
                out[p] = float(base[p])
            else:
                out[p] = float(pdefs[p]["value"])
                filled.append(p)
        return cls(
            np_model=np_model,
            np_model_nu=np_model_nu,
            values=out,
            filled=tuple(filled),
        )

    @classmethod
    def from_flat(cls, values, np_model, np_model_nu, what="tune"):
        """Build from a flat ``{λ name: value}`` mapping (the param-model names).

        The mapping is treated as a BASE, not a strict override set: λ it carries
        that the chosen forms ignore are DROPPED rather than rejected. That is the
        right reading for a stored λ set — a fitresults holds the λ of the form it
        was fitted with, which need not match the form being evaluated now. λ the
        forms need but the mapping lacks come from the registry defaults and are
        listed in :attr:`filled` for the caller to report."""
        return cls.create(np_model, np_model_nu, base=values, what=what)

    @property
    def eff(self):
        """F_eff kwargs: the TMD λ plus the ``np_model`` selector."""
        d = {k: self.values[k] for k in EFF_PARAMS if k in self.values}
        d[EFF_MODEL_KEY] = self.np_model
        return d

    @property
    def gnu(self):
        """γ_ν^NP kwargs: the CS λ plus the ``np_model_nu`` selector."""
        d = {k: self.values[k] for k in GNU_PARAMS if k in self.values}
        d[GNU_MODEL_KEY] = self.np_model_nu
        return d

    def replace(self, np_model=None, np_model_nu=None, what="tune", **values):
        """A new tune with some forms and/or λ changed, re-validated.

        Changing a form re-derives the active λ set: λ the new forms drop are
        discarded, λ they add are filled from the registry (and recorded in
        :attr:`filled`). Existing λ carry over via ``base``."""
        return type(self).create(
            np_model or self.np_model,
            np_model_nu or self.np_model_nu,
            values,
            base=self.values,
            what=what,
        )

    def as_lambda_central(self):
        """The ``{"eff_params", "gnu_params"}`` mapping the TF core takes as its
        ``lambda_central`` construction argument."""
        return {"eff_params": self.eff, "gnu_params": self.gnu}

    def describe(self):
        """One-line ``np_model/np_model_nu: λ=…`` summary for CLI logging."""
        lam = ", ".join(f"{k}={v:g}" for k, v in self.values.items())
        return f"{self.np_model}/{self.np_model_nu}: {lam}"


def parse_lambda_overrides(spec, what="--lambdas"):
    """Parse a ``"name=val,name=val"`` λ-override string into ``{name: float}``.

    Hard-errors (``ValueError``) on a malformed token, an unknown parameter name
    (not in :data:`ALL_PARAMS`), or a non-float value. An empty / ``None`` spec
    yields ``{}``. Model-awareness (whether a *known* λ is used by the chosen
    model) is the caller's job via :func:`active_params`. ``what`` names the flag
    the spec came from, so a caller with several λ flags (e.g. ``--lambdas`` and
    ``--ratio-lambdas``) points the error at the right one."""
    out = {}
    for tok in (spec or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" not in tok:
            raise ValueError(f"{what}: expected 'name=value', got {tok!r}")
        k, v = tok.split("=", 1)
        k = k.strip()
        if k not in ALL_PARAMS:
            raise ValueError(
                f"{what}: unknown NP parameter {k!r} "
                f"(known: {', '.join(ALL_PARAMS)})"
            )
        try:
            out[k] = float(v)
        except ValueError:
            raise ValueError(f"{what}: {k}={v.strip()!r} is not a float")
    return out


def check_lambda_overrides(
    overrides, np_model=None, np_model_nu=None, what="--lambdas"
):
    """``ValueError`` if any λ in ``overrides`` is inert for the chosen model(s).

    Companion to :func:`parse_lambda_overrides`, which checks the NAMES but not
    whether the model uses them: setting a λ the form factor ignores is a silent
    no-op, so a caller that knows its models makes it a hard error instead."""
    active = active_params(np_model=np_model, np_model_nu=np_model_nu)
    inert = [k for k in overrides if k not in active]
    if inert:
        raise ValueError(
            f"{what}: "
            + ", ".join(inert)
            + f" not used by np_model={np_model} / np_model_nu={np_model_nu} "
            + "(active: "
            + ", ".join(k for k in (*EFF_PARAMS, *GNU_PARAMS) if k in active)
            + ")"
        )


def split_eff_gnu(values):
    """Split a ``{name: value}`` mapping into ``(eff_params, gnu_params)`` dicts
    by membership in EFF_PARAMS / GNU_PARAMS (values floated; names in neither,
    e.g. the model-name keys, are dropped)."""
    eff = {k: float(values[k]) for k in EFF_PARAMS if k in values}
    gnu = {k: float(values[k]) for k in GNU_PARAMS if k in values}
    return eff, gnu


def bin_sum_matrix(src_centers, target_edges, tol=1e-6):
    """(N_target, N_src) 0/1 matrix summing source bins whose centre falls in
    each target bin. Source bins outside every target bin get 0, truncating to
    the target range (e.g. qT > ptVGen_max, |Y| > absY_max)."""
    src = np.asarray(src_centers, dtype=np.float64)
    edges = np.asarray(target_edges, dtype=np.float64)
    W = np.zeros((edges.size - 1, src.size), dtype=np.float64)
    for i in range(edges.size - 1):
        m = (src >= edges[i] - tol) & (src <= edges[i + 1] + tol)
        W[i, m] = 1.0
    return W
