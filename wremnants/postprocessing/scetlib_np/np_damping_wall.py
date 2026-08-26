"""Physical-damping wall for the continuous-λ SCETlib NP param model.

Companion to ``wremnants/postprocessing/np_monotonicity.py``: same rabbit
``Regularizer`` mechanism and hinge-loss (relu²) wall structure, for the
:class:`~wremnants.postprocessing.scetlib_np.param_model.SCETlibNPParamModel`,
whose λ are direct fit parameters (the model-param block of ``x``), not discrete
template nuisances. No ``PARAM_MAP`` / θ-interpolation: the regularizer reads the
physical λ from ``x[:nparams]`` in the model's canonical parameter order.

The model FORMS and λ ORDER are DERIVED from the model, not re-declared here.
``SCETlibNPParamModel.__init__`` publishes itself on the shared ``indata``
(``indata.scetlib_np_param_model``); this regularizer — built afterward with the
same ``indata`` — reads ``model.fit_forms`` (the F_eff / γ_ν NUMERATOR forms the
fit integrates) and ``model._param_order`` (POI-reordering and all). So you pass
the form ONCE, to the model; nothing on the ``-r`` line repeats it, and a
``poi_params`` reorder is tracked automatically. The wall RAISES if no model has
registered, or if a side's form is one it has no walls for (see dispatch below).

Why a wall (see ``param_model.py``): a wrong-sign λ anti-damps the NP form
factors, so the bT-integral diverges / the differential σ(qT) oscillates negative
— a genuine unphysical region. The fit's qT→ptVGen rebin launders it by averaging
the negative differential σ away, so σ_gen / σ_reco / the NLL stay finite and the
minimizer cannot see the unphysical-ness. These λ must lie in the physical
(damping) region; this regularizer adds the fit-time penalty keeping them there.

Walls, per side and per FORM, read off the TF forms the fit integrates
(``btgrid_tf.gamma_nu_NP_tf`` / ``btgrid_tf.F_eff_tf``). The damping criterion is
"the tanh argument is ≥ 0 ∀b" (γ_ν ≤ 0 / F_eff decays), i.e. a polynomial in
u ≡ b² is non-negative on u ≥ 0:

  CS side   γ_ν^NP(b) = −λ∞_ν · tanh( P(u)/λ∞_ν ) ,
            P(u) = λ2_ν·u + λ4_ν·u² + λ6_ν·u³      (λ6_ν ≡ 0 for tanh_2)
      tanh_2:  λ∞_ν > 0 ,  λ2_ν ≥ 0 ,  λ4_ν ≥ 0
      tanh_6:  λ∞_ν > 0 ,  λ2_ν ≥ 0 ,  λ6_ν ≥ 0 (leading) ,
               and λ4_ν ≥ 0  OR  λ4_ν² ≤ 4·λ2_ν·λ6_ν   (interior; λ4_ν may dip
               negative while the b⁶ term keeps P ≥ 0)

  TMD side  F_eff(Y,b) = exp(−2·λ∞·b·tanh(a)) ,
            a·λ∞ = b·Q(u) ,  Q(u) = λ2_Y + B·u + λ6·u² ,
            B = λ4 + λ2_Y³/(3λ∞²) ,  λ2_Y = λ2 + δλ2·Y²    (λ6 ≡ 0 for tanh_2)
      tanh_2:  λ∞ > 0 ,  λ2_Y ≥ 0 ,  B ≥ 0  [≡ 3·λ∞²·λ4 + λ2_Y³ ≥ 0]
      tanh_6:  λ∞ > 0 ,  λ2_Y ≥ 0 ,  λ6 ≥ 0 (leading) ,
               and B ≥ 0  OR  B² ≤ 4·λ2_Y·λ6                (interior)
            evaluated at Y=0 and Y=Y_MAX (covers δλ2 of either sign — λ2_Y is
            monotonic in Y², so the binding |Y| is one of the two extremes).

      tanh_6_sigmoid:  F_eff is the tanh_6 form TIMES the turn-off
            S(b) = 1/(1+exp((b−bT_cutoff)/w)), so F_eff ≤ 1 reads, with
            L(b) ≡ ln(1+exp((b−bT_cutoff)/w)) ≥ 0,
                2·λ∞·b·tanh(a) ≥ −L(b)
            — the threshold on the tanh argument is no longer the constant 0 but a
            b-dependent NEGATIVE floor (the sigmoid buys the form factor the right
            to anti-damp, by more and more as b grows). Minimising that over b has
            no closed form (tanh × sech² × logistic), so instead:

            F_eff > 1 requires tanh a < 0, i.e. Q(u) < 0, which for a QUADRATIC Q
            is one interval (u₋, u₊) — the quadratic formula. Outside it there is
            nothing to check. On it, bound each side at its worst:
              demand ≤ 2·λ∞·√u₊·|tanh a_min|   (largest b × largest |tanh a|)
              supply ≥ L(√u₋)                   (L increasing → least at left edge)
            and a_min is analytic too: d/db[b·Q(b²)] = 0 is λ2_Y+3B·b²+5λ6·b⁴ = 0,
            again a quadratic in b². Hence the single closed-form condition

                2·λ∞·√u₊·|tanh a_min|  ≤  L(√u₋)     (margin: a_min ≥ margin)

            replacing tanh_6's interior discriminant. λ2_Y ≥ margin (small-b turn-on)
            and λ6 ≥ margin (keeps the window bounded) are UNCHANGED — the turn-off
            does not help at small b, and λ6 is frozen positive in practice anyway.

            Self-gating exactly like the term it replaces: B ≥ 0 or B² ≤ 4·λ2_Y·λ6
            → no window → ``relu(−tanh(a_min−margin))`` is exactly 0 → the term
            vanishes (no spurious offset).

            Its FEASIBLE SET reduces to tanh_6's as bT_cutoff → ∞ (supply → 0 ⇒
            tanh a_min ≥ 0 ⇒ Q ≥ 0), mirroring the model's own reduction property —
            verified over 3000 random tunes at margin=0, zero verdict
            disagreements. The penalty VALUES differ (different functions, only the
            accept/reject boundary coincides).

            Sound but slightly pessimistic: worst demand (right edge) is paired with
            least supply (left edge) though they occur at different b. Measured on a
            20k-tune scan: 0 unsound passes, 0.15 % conservative rejections.

            λ∞ may float: every term above reads the live λ∞, nothing is
            precomputed. (The asymptotic worry — at 2·λ∞·w ≥ 1 the turn-off cannot
            beat an e^{+2λ∞b} tail — needs λ6 < 0 to arise, and λ6 ≥ margin is
            still enforced, so the window is always bounded and the question never
            comes up.)

      tanh_6_abs:  NO WALL AT ALL — the form is unconditionally damping. Its
            damping exponent is folded through ``btgrid_tf.abs_fold_tf``
            (|exponent + margin| − margin), so F_eff ≤ e^margin and γ_ν ≤ margin
            for EVERY λ, of either sign, λ∞ included. There is no infeasible λ to
            wall off, so both this side's penalty terms AND its λ∞ floor are
            dropped and the contribution is EXACTLY 0. The form is supported here
            (rather than raising) only so a ``-r`` line inherited across the
            fit/hessian/saturated steps of ``fitterSCETlibNP.py`` stays harmless;
            the constructor says so out loud. NOTE what the fold does and does not
            buy: the PREDICTION is physical for any λ, the λ themselves are not, so
            a folded tune is a damping-but-kinked form factor that plain tanh_6
            cannot produce — judge it with the postfit fold-activity diagnostic in
            ``param_model_diagnostics``, not with this wall's silence.

The tanh_2 walls are EXACT in λ∞ (read off the actual TF forms, not the AN
λ∞-normalised parametrisation). They are the λ6 = 0 reduction of the tanh_6 walls
(then the leading coefficient is λ4_ν / B, hence the simpler ≥0 limit), and of the
``np_monotonicity.py`` monotonicity walls. NOTE this file uses the DAMPING
criterion (P ≥ 0), which lets the interior coefficient dip negative within the
discriminant bound; ``np_monotonicity.py`` uses the stricter MONOTONICITY
criterion (√3 in place of √4). The interior relu²(relu(−λ4_ν)² − 4·λ2_ν·λ6_ν)
penalty form is self-gating: it vanishes for λ4_ν ≥ 0 (no division by λ6).

THE MARGIN, AND WHAT IT DOES NOT MEAN. Every condition above is enforced as
``coeff ≥ margin`` rather than ``coeff ≥ 0``, so a POSITIVE margin is STRICTER than
physical and a NEGATIVE one buys slack. (That sign convention is the opposite of the
``tanh_6_abs`` fold's ``abs_margin``, where positive = more room; the two knobs are
not interchangeable.) Three things to keep in mind:

  * The margin is NOT commensurate across terms. One number is compared against
    quantities of different DIMENSION and scale: λ2_Y [GeV, ~0.25], λ6 [GeV⁵, 0.016],
    λ6_ν [GeV⁶, 7e-4], and (sigmoid) a dimensionless tanh argument. margin = −0.005
    simultaneously means "0.5 % of λ2" and "31 % of λ6" and "0.5 % of the tanh
    argument". There is no reading of one absolute number that is natural for all of
    them; the unified minP/minQ terms at least put their cushion in λ2-like
    coefficient units, which is why they are written that way.
  * The margin does NOT bound the form factor. It constrains COEFFICIENTS, and how
    much F_eff excursion a given coefficient dip produces depends on where in b_T the
    dip sits. Measured (λ∞=1, λ6 frozen +0.01, random-λ scan, worst sampled
    max_b F_eff over the accepted set): margin 0 → 1.000, −0.001 → 1.016,
    −0.005 → 1.09, −0.01 → 1.24, −0.02 → 1.55, −0.05 → 3.2. So "F_eff ≲ 1.1 of
    headroom" is margin ≈ −0.005, NOT −0.1 — and that is a sampled correspondence,
    not a bound (with λ6 FLOATING, a negative margin admits λ6 < 0, where minQ is not
    the true minimum, and the excursion is unbounded). For a guaranteed cap on the
    function itself, use the ``tanh_6_abs`` fold, whose margin IS the excursion.
  * A condition on a FROZEN coefficient is either satisfied or an irreducible
    constant — the wall still prices it. With λ6 = λ6_ν = 0.01 frozen and margin
    +0.1, the three λ6 terms contribute a fixed 2·(0.09)² that no fit can pay down,
    × exp(2·strength): ≈ 535 of pure offset at strength 5, present in the tanh_6 and
    sigmoid branches alike. Harmless at margin ≤ λ6, and negative margins never
    trigger it, but it makes positive-margin runs' reported χ² incomparable.

Wall hardness is set at fit time by rabbit's ``--regularizationStrength`` (the
penalty × ``exp(2·tau)`` in ``fitter.py``; ``tau`` is a fixed multiplier, NOT a
minimised parameter). A large strength makes this a BARRIER: ≈0 inside the
physical region, steeply rising outside. Free (small strength) vs walled (large
strength) Δχ² is the data–model tension diagnostic: railing against a wall with
large Δχ² is genuine tension, not a masked pathology.

The small-b turn-on walls (λ2_ν ≥ 0, λ2_Y ≥ 0) are a stronger condition than the
large-b limit (they forbid an anti-damping bump near b→0, not just the wrong
asymptote). They can be switched off with the mapping flag ``smallb=0`` — then
ONLY the limiting/interior behaviour and the λ∞ floors are enforced, and the
leading b² coefficient floats either sign (use the postfit σ(qT)≥0 check in
``param_model_diagnostics`` as the real guard then).

Invoke (the model form is taken from ``--paramModel``; nothing repeats it here):

    rabbit_fit.py ... \\
      --regularizationStrength 3 \\
      -r wremnants.postprocessing.scetlib_np.np_damping_wall.NPDampingWall \\
         wremnants.postprocessing.scetlib_np.np_damping_wall.NPDampingMapping \\
         [smallb=0]

(``Y_MAX`` = 2.5 — the |y| ACCEPTANCE, since 2026-08-05; it was 5.0 before, which
is past the kinematic limit ln(√s/Q)=4.96 and past the range our Y-parametrization
was fitted on — see the constant's comment. The λ∞ floor and the damping margin are
likewise fixed module constants, not -r options. ``ymax=`` IS accepted on the -r
line and should be passed EXPLICITLY so a fit's meta records the wall it used.)

References:
  AN-25-085 theory.tex Eqs. eq:npgamma, eq:npf;
  param_model.py / sigma_gen.py docstrings (the σ_gen pipeline & the binning-
  launders-the-pathology discussion).
"""

# rabbit / TF imports deferred to the lazy class factories so the module stays
# importable without rabbit/TF (mirrors np_monotonicity.py). The λ registry
# (params) is numpy-only, so it is safe to import at module level.
from wremnants.postprocessing.scetlib_np.params import active_params

# Forms this wall has damping conditions for. Anything else (frac_*, exp_*,
# tanh_1, tanh_4, signed_lambda, identity, …) raises rather than silently
# applying the wrong tanh_2/tanh_6 walls.
SUPPORTED_FORMS = ("tanh_2", "tanh_6", "tanh_6_sigmoid", "tanh_6_abs")

# Forms that need NO wall: the *_abs damping fold makes them physical for every λ
# (see the module docstring). Accepted here so an inherited -r line is a no-op
# rather than a crash; the per-side penalty is exactly 0.
# Forms that are damping BY CONSTRUCTION, so this regularizer has nothing to add:
# its whole block for that side (λ∞ floor included) is skipped.
#   tanh_6_abs  — the damping FOLD: γ_ν ≤ abs_margin_nu / F_eff ≤ 1+abs_margin for
#                 any λ of either sign.
#   tanh_6_pos  — lambda2_nu mapped onto (0, ∞) by ``btgrid_tf.pos_floor_tf``.
#                 UNCONDITIONAL ONLY IF lambda4_nu is frozen at 0 and lambda6_nu
#                 frozen ≥ 0, which is how the form is meant to be run: all the
#                 argument's coefficients are then non-negative. Floating
#                 lambda4_nu negative would break the guarantee AND silence the
#                 wall — ``_warn_if_conditional_pos`` catches that at init.
#   tanh_2_pos  — the tanh_2 analogue, on BOTH sides, and UNCONDITIONAL with no
#                 strings attached: CS maps lambda2_nu and lambda4_nu (= exactly
#                 this wall's tanh_2 CS conditions), TMD maps lambda2_Y by its
#                 anchors and the combination B = λ4 + λ2_Y³/(3λ∞²) (= exactly
#                 this wall's tanh_2 TMD conditions). Nothing needs freezing
#                 beyond the shape constants, and no λ∞ floor is needed either
#                 (λ·tanh(c/λ) is even in λ). Whichever side runs it, that side's
#                 penalty is 0 because its feasible set is the whole space.
UNCONDITIONAL_FORMS = ("tanh_6_abs", "tanh_6_pos", "tanh_2_pos")

# btgrid_tf form aliases that resolve to a SUPPORTED_FORMS entry (mirrors the
# alias maps in gamma_nu_NP_tf / F_eff_tf). Other aliases ("linear"->frac_1,
# "square_root"->frac_2) resolve to unsupported forms and so fall through to the
# raise. Applied to both sides; only "hyp_tangent" reaches a supported form.
_FORM_ALIASES = {"hyp_tangent": "tanh_2"}

# Fixed knobs (NOT CLI options — edit here to test). Only ``smallb`` is exposed
# on the -r line; these are deliberately constants to keep that line minimal.
Y_MAX = 2.5  # binding |y| for the F_eff Y-evaluation: the |y| ACCEPTANCE.
#
# CHANGED 5.0 -> 2.5 on 2026-08-05 (Luca). The old default demanded physical damping
# at |Y| = 5, which fails three independent tests:
#   1. NOT the kinematic ceiling, just past it. For Z at 13 TeV,
#      Y_max = ln(√s/Q) = ln(13000/91.1876) = 4.96, and at Y = 5 the momentum
#      fraction x1 = (Q/√s)e^Y = 1.04 > 1 — the configuration does not exist.
#   2. Outside the acceptance. The Z analysis measures |Y| <= 2.5; everything
#      beyond is extrapolation the fit has no data to constrain.
#   3. Outside the range our own Y-parametrization was FITTED on. lambda2_Y =
#      lambda2 + delta_lambda2*Y^2 is a mapping of the external tune onto our form,
#      and studies/np-wall-local-minima/scripts/map_to_ours.py fits it on
#      Ys = [0, 0.5, 1, 1.5, 2, 2.5] only. A downward parabola extrapolated past its
#      fit range must eventually cross zero — for the MAP22/lattice external tune
#      (lambda2=0.09754, delta_lambda2=-0.00785) that happens at |Y| = 3.52.
# So the EXTERNAL INPUTS THEMSELVES fail the Y=5 test, and it made no sense to hold
# the fit to a condition its own reference values do not satisfy. That was Luca's
# argument for the change.
#
# ⚠️ REPRODUCIBILITY. Fits run BEFORE 2026-08-05 that recorded no explicit `ymax=`
# token in their meta ran at 5.0, and re-running them now reproduces a DIFFERENT
# wall. Their meta cannot distinguish the two. To keep new fits self-describing
# regardless of any future change to this constant, PASS `ymax=` EXPLICITLY on the
# -r line (hessian/saturated then inherit it through meta). The `*_ymax2p5` runs do.
#
# ⚠️ NOT propagated to the summary tools. fit_summary_table.py / fit_summary_plotly.py
# still evaluate their "max F_eff / max g_nu" physicality columns at |Y| = 5, so a fit
# walled to 2.5 can still be reported UNPHYSICAL there on large-|Y| extrapolation.
# The two are answering different questions; do not read a disagreement as a bug.
LAMBDA_INF_FLOOR = 1e-3  # positive floor on the λ∞ saturation scales
NP_DAMPING_MARGIN = 0.0  # enforce each damping coefficient ≥ 0, the physical
#              condition itself, with no invented cushion (Luca 2026-08-10:
#              "0.005 is too arbitrary"). Set > 0 per run with the `margin=`
#              token on the -r line if you want the cushion back.
#
#              WHAT THE CUSHION WAS FOR, and what dropping it costs. The wall is
#              a hinge, relu2(margin - coeff): its gradient VANISHES at the knee.
#              So a coefficient the data pushes downward settles where the data's
#              pull balances a restoring force that is going to zero — i.e. a hair
#              INSIDE the knee, never exactly on it. Measured on the 2D fits at
#              margin = 5e-3, the leak is consistently ~4e-4:
#                  MSHT20 nominal  λ2_ν      = 0.00456  (margin - 0.00044)
#                  CT18Z  nominal  λ2_ν      = 0.00460  (margin - 0.00040)
#                  MSHT20+lattice  λ2_Y@Y2.5 = 0.00468  (margin - 0.00032)
#              At margin = 0 that same leak puts those coefficients at about
#              -4e-4: very slightly ANTI-damping rather than exactly 0. That is
#              ~500x smaller than the |λ2_ν| ≈ 0.23 excursions the wall exists to
#              prevent, so it is cosmetic for the prediction — but "λ ≥ 0" is now
#              satisfied only to O(1e-4), not exactly, and a postfit λ of -4e-4
#              should be read as "on the wall", not as a physical result.
#
#              A cushion is not the only cure: a penalty whose gradient does not
#              vanish at the knee (linear hinge, or a barrier) would sit exactly
#              at 0. Both were rejected here — the linear hinge puts a gradient
#              discontinuity right where the minimizer spends its time, which
#              this trust-region setup is already fragile about.
#
#              NB fits made before 2026-08-10 used 5e-3 and their postfit λ are
#              therefore offset by that much; they are not directly comparable
#              with margin = 0 runs. AND THE OLD ONES CANNOT TELL YOU WHICH THEY
#              USED: a run that relies on the default stores no `margin=` token,
#              so its meta records only `ymax=...` and the margin in force is
#              unrecoverable from the file — it is whatever this constant happened
#              to be on the day. Re-running such a fit after this change silently
#              gives it a different wall. The construction banner below now prints
#              the effective margin so at least the log carries it; pass `margin=`
#              explicitly on the -r line for anything whose provenance matters.


def _make_mapping_class():
    from rabbit.mappings.mapping import BaseMapping

    class NPDampingMapping(BaseMapping):
        """Vestigial BaseMapping carrying the wall's options to the regularizer.

        Options (``key=value`` tokens, optional):
            smallb=<0|1>   enforce the small-b turn-on walls λ2_ν≥0 and λ2_Y≥0
                           (default 1). smallb=0 drops them, keeping ONLY the
                           large-b limit/interior walls and the λ∞ floors — i.e.
                           constrain the limiting behaviour but let the leading
                           small-b coefficient float either sign.
            margin=<float> signed damping cushion; every condition is enforced
                           ``≥ margin`` (default ``NP_DAMPING_MARGIN``). margin>0 is
                           STRICTER than physical, margin=0 is bare positivity,
                           margin<0 permits a controlled amount of anti-damping —
                           note this is the OPPOSITE sign convention to the
                           tanh_6_abs fold's ``abs_margin``. Applied to ALL terms
                           (turn-on, interior, the sigmoid window gate), which are
                           NOT in commensurate units, and it does NOT bound the form
                           factor: "F_eff ≲ 1.1" is margin ≈ −0.005, not −0.1. It
                           also prices conditions on FROZEN λ. Read the module
                           docstring's margin section before scanning it.

        The model forms and λ order are derived from the registered
        ``SCETlibNPParamModel`` (see module docstring). The binding |Y| (``Y_MAX``)
        and the λ∞ floor (``LAMBDA_INF_FLOOR``) are FIXED module constants, not CLI
        options — edit them in this file to test, kept off the -r line on purpose.
        """

        def __init__(
            self, indata, key, smallb=True, margin=NP_DAMPING_MARGIN, ymax=Y_MAX
        ):
            super().__init__(indata, key)
            self.indata = indata
            self.smallb = bool(smallb)
            self.margin = float(margin)
            self.ymax = float(ymax)

        @classmethod
        def parse_args(cls, indata, *args):
            smallb = True
            margin = NP_DAMPING_MARGIN
            ymax = Y_MAX
            for a in args:
                if "=" not in a:
                    raise ValueError(
                        f"NPDampingMapping: arg must be 'key=value', got '{a}'"
                    )
                k, v = a.split("=", 1)
                if k == "smallb":
                    smallb = v.strip().lower() not in ("0", "false", "no", "off")
                elif k == "margin":
                    # signed damping cushion: coeff >= margin. margin>0 is stricter
                    # than physical, margin=0 bare positivity, margin<0 allows a
                    # controlled amount of anti-damping. Applied to all terms (small-b
                    # turn-on and interior discriminants); the frozen-λ6 large-b limit
                    # is physical by construction regardless.
                    margin = float(v)
                elif k == "ymax":
                    # |Y| at which the TMD-side conditions are evaluated, via
                    # lambda2_Y = lambda2 + delta_lambda2 * Y^2. The default is 2.5 =
                    # the |y| acceptance edge (since 2026-08-05; it was 5.0, which is
                    # past the kinematic limit 4.96 AND past the fit range of the
                    # Y-parametrization). Because the lever is Y^2 this is by far the
                    # strongest constraint on delta_lambda2 -- at Y=0 the term vanishes
                    # identically, so Y=ymax is the ONLY place delta_lambda2 is walled
                    # at all, and 5 -> 2.5 weakens it by exactly 4x (25 -> 6.25).
                    # ymax=0 evaluates ONLY Y=0, which leaves
                    # delta_lambda2 completely UNCONSTRAINED by this wall -- use the
                    # postfit sigma(qT)>=0 / F_eff checks in param_model_diagnostics as
                    # the guard then. Requested by Luca 2026-08-03.
                    ymax = float(v)
                    if ymax < 0:
                        raise ValueError(
                            f"NPDampingMapping: ymax must be >= 0, got {ymax}"
                        )
                else:
                    raise ValueError(
                        f"NPDampingMapping: unknown key '{k}'; supported: "
                        f"'smallb=<0|1>', 'margin=<float>', 'ymax=<float>' "
                        f"(eps is a fixed module constant)."
                    )
            return cls(
                indata,
                f"{cls.__name__} smallb={int(smallb)} margin={margin:g}",
                smallb=smallb,
                margin=margin,
            )

    return NPDampingMapping


def _make_regularizer_class():
    import tensorflow as tf

    from rabbit.regularization.regularizer import Regularizer

    class NPDampingWall(Regularizer):
        """Hinge-loss penalty enforcing NP damping, per-side and per-form
        (tanh_2 / tanh_6); see the module docstring."""

        def __init__(self, mapping, dtype):
            super().__init__(mapping, dtype)
            self.dtype = dtype
            self.mapping = mapping
            self.indata = mapping.indata
            # ymax / eps are fixed module constants (not CLI options); smallb and
            # margin come from the mapping (margin default = NP_DAMPING_MARGIN).
            self.ymax = float(getattr(mapping, "ymax", Y_MAX))
            self.eps = LAMBDA_INF_FLOOR
            self.margin = float(getattr(mapping, "margin", NP_DAMPING_MARGIN))
            self.enforce_small_b = bool(getattr(mapping, "smallb", True))

            # Forms + λ order are DERIVED from the SCETlibNPParamModel that
            # published itself on the shared indata (built before this regularizer
            # in rabbit_fit; see param_model.py). No -r-line repetition.
            model = getattr(self.indata, "scetlib_np_param_model", None)
            if model is None:
                raise ValueError(
                    "NPDampingWall: no SCETlibNPParamModel registered on indata. "
                    "This wall derives the NP form and λ order from the param "
                    "model, so the fit must use "
                    "--paramModel ...scetlib_np.param_model.SCETlibNPParamModel "
                    "(which publishes indata.scetlib_np_param_model)."
                )

            # NAMES only. Positions are resolved in set_expectations against the
            # vector actually being fit, because the saturated goodness-of-fit path
            # reorders and resizes it (see there). Caching positions here is what
            # made the wall silently read histogram bins as λ.
            self._order = tuple(model._param_order)
            self._pidx = None

            # FIT (numerator) forms — the ones the fit integrates, which the wall
            # must constrain (NOT the card/denominator form). Resolve aliases and
            # fail on any form we have no walls for.
            forms = model.fit_forms
            self._np_model_nu = self._resolve_form(
                forms["np_model_nu"], side="CS (γ_ν, np_model_nu)"
            )
            self._np_model = self._resolve_form(
                forms["np_model"], side="TMD (F_eff, np_model)"
            )
            # *_abs sides need no wall (module docstring): say so, loudly, since a
            # silent zero penalty otherwise reads as "I ran walled".
            for side, form in (
                ("CS (γ_ν)", self._np_model_nu),
                ("TMD (F_eff)", self._np_model),
            ):
                if form in UNCONDITIONAL_FORMS:
                    print(
                        f"[NPDampingWall] {side} form {form!r} is unconditionally "
                        "damping by construction; this side contributes EXACTLY 0 "
                        "to the penalty — the run is effectively UNWALLED there.",
                        flush=True,
                    )

            # tanh_6_pos is unconditional ONLY with lambda4_nu frozen at 0 and
            # lambda6_nu frozen ≥ 0. That precondition is enforced in
            # SCETlibNPParamModel._check_pos_form_preconditions, which is where the
            # rabbit freeze list lives; by the time the wall is built it has already
            # raised, so there is nothing left to check here.

            # Required λ come from the CENTRAL REGISTRY, keyed on the resolved
            # forms — the SAME source (active_params) the model uses to build
            # _param_order. Each model has its own λ vocabulary (e.g. tanh_2 has
            # no λ6*), so the wall can never require a λ the chosen model omits.
            required = active_params(
                np_model=self._np_model, np_model_nu=self._np_model_nu
            )
            missing = [p for p in sorted(required) if p not in self._order]
            if missing:
                raise ValueError(
                    f"NPDampingWall: model param order {self._order} is missing "
                    f"λ {missing} required by the fit forms "
                    f"(np_model={self._np_model!r}, np_model_nu={self._np_model_nu!r})."
                )

            # Announce the settings actually in force. A run that relies on the
            # defaults records NO margin=/ymax= token in its meta, so without this
            # the wall a result was produced under is unrecoverable once the
            # constants here change -- and they have (Y_MAX 5.0 -> 2.5 on
            # 2026-08-05, NP_DAMPING_MARGIN 5e-3 -> 0 on 2026-08-10). The log is
            # not as good as provenance in the file, but it is not nothing.
            print(
                f"[NPDampingWall] margin={self.margin:g} ymax={self.ymax:g} "
                f"smallb={int(self.enforce_small_b)} "
                f"forms: CS={self._np_model_nu!r} TMD={self._np_model!r} "
                f"λ order={self._order}",
                flush=True,
            )

            self._cast = lambda v: tf.constant(v, dtype=self.dtype)
            # Model-param block is x[:nparams]; nparams resolved at set_expectations.
            self._nparams = None

        @staticmethod
        def _resolve_form(form, side):
            resolved = _FORM_ALIASES.get(form, form)
            # UNCONDITIONAL_FORMS are accepted too — they need no wall, so "no
            # walls for it" is the right answer, not an error. Without this the
            # *_pos forms raised here on any -r line (the *_abs forms happen to be
            # in SUPPORTED_FORMS as well, which is why it never showed up).
            if resolved not in SUPPORTED_FORMS + UNCONDITIONAL_FORMS:
                raise NotImplementedError(
                    f"NPDampingWall: no damping walls for {side} form {form!r}"
                    + (f" (resolves to {resolved!r})" if resolved != form else "")
                    + "; supported: "
                    + f"{sorted(set(SUPPORTED_FORMS + UNCONDITIONAL_FORMS))}"
                    + ". Add walls for it or run that side with a supported form."
                )
            return resolved

        def set_expectations(self, initial_params, initial_observables, parms=None):
            # Resolve every λ BY NAME against the vector we are about to be
            # evaluated on. Positions are NOT stable across a fit session: the
            # saturated goodness-of-fit path swaps the model for
            # CompositeParamModel([original, SaturatedProjectModel]), whose
            # [POIs | POUs] layout inserts one POI per projected bin ahead of the
            # original model's block and displaces every λ by that count.
            #
            # This used to cache _pidx from model._param_order at construction and
            # index the full vector with it. In a 2D ptll x yll fit that meant the
            # saturated pass read saturated_ch0_ptll0..6 (all ≈0.94, far above the
            # 0.005 margin) as the 7 λ: every condition trivially satisfied, penalty
            # identically zero, the reference model effectively UNWALLED while the
            # nominal one it is compared against was walled. It ran to
            # lambda2_nu = -0.23, a point that costs 1225 nll units when the wall is
            # really applied. Silent since the walls were introduced.
            if parms is not None and hasattr(self, "resolve_indices"):
                self._pidx = self.resolve_indices(
                    parms, self._order, who="NPDampingWall"
                )
                self._nparams = len(self._order)
                return

            # COMPATIBILITY with a rabbit that does not pass parameter names yet
            # (anything before the arm_regularizers change). Fall back to the
            # historical positional mapping so this module keeps working, but say
            # so — with that rabbit the wall CANNOT verify what it is reading, and
            # on the saturated path it is silently wrong in exactly the way
            # described above. The width guard below still catches the composite
            # layout wherever set_expectations is actually reached.
            print(
                "[NPDampingWall] WARNING: this rabbit does not pass parameter "
                "names to set_expectations, so lambda positions cannot be "
                "verified. Falling back to positional indexing. The saturated "
                "goodness-of-fit path is UNRELIABLE on this combination — use a "
                "rabbit with arm_regularizers() for any walled GoF result.",
                flush=True,
            )
            self._pidx = {name: i for i, name in enumerate(self._order)}
            nsyst = len(self.indata.systs)
            self._nparams = int(initial_params.shape[0]) - nsyst
            if self._nparams != len(self._order):
                raise ValueError(
                    f"NPDampingWall: the fit's model-param block is {self._nparams} "
                    f"wide but the model param order has {len(self._order)} entries "
                    f"{self._order}. A wrapping/composite param model (e.g. the "
                    "saturated goodness-of-fit path) reorders/resizes the block in "
                    "a way this wall's flat indexing cannot follow."
                )

        def _lam(self, params, name):
            # λ stored directly in the model-param block (allowNegativeParam=True),
            # so x[index] IS the physical λ — no theta interpolation. The index is
            # resolved by name in set_expectations against the CURRENT layout; it
            # must never be cached across layouts (see there).
            return params[self._pidx[name]]

        def constraint_spec(self, params, observables):
            """The damping conditions as HARD constraints ``values >= lower``.

            Same feasible set the penalty hinges on, handed to a constrained
            minimizer instead of being added to the loss. Implemented for
            tanh_2/tanh_2, which is the configuration the alpha_s fits use;
            other vocabularies raise rather than silently enforcing a subset,
            because a dropped condition looks exactly like a satisfied one.

            Conditions, for each side that is not an UNCONDITIONAL_FORM:
              CS   lambda_inf_nu >= eps,  lambda4_nu >= m,  [lambda2_nu >= m]
              TMD  lambda_inf    >= eps,  and per |Y| in {0, ymax}:
                   [lambda2_Y >= m]  and  3*lambda_inf^2*lambda4 + lambda2_Y^3 >= m
            The bracketed small-b turn-on terms are present only when smallb=1.
            Both per-Y families are MONOTONIC in Y^2 (d/dY^2 of lambda2_Y is
            delta_lambda2, and of the cubic is 3*lambda2_Y^2*delta_lambda2, whose
            sign cannot flip), so the two endpoints bound the whole acceptance
            and no interior-in-Y condition exists to miss.
            """
            for side, form in (
                ("np_model_nu", self._np_model_nu),
                ("np_model", self._np_model),
            ):
                if form not in UNCONDITIONAL_FORMS and form != "tanh_2":
                    raise NotImplementedError(
                        f"[NPDampingWall] constraint_spec supports tanh_2 (or an "
                        f"unconditional form) per side; {side}={form!r}. Use the "
                        f"penalty path (a non-constrained minimizer) for that form."
                    )

            eps = self._cast(self.eps)
            m = self._cast(self.margin)
            vals, lows = [], []

            def add(expr, lower):
                vals.append(tf.reshape(tf.cast(expr, self.dtype), [1]))
                lows.append(tf.reshape(tf.cast(lower, self.dtype), [1]))

            if self._np_model_nu not in UNCONDITIONAL_FORMS:
                add(self._lam(params, "lambda_inf_nu"), eps)
                add(self._lam(params, "lambda4_nu"), m)
                if self.enforce_small_b:
                    add(self._lam(params, "lambda2_nu"), m)

            if self._np_model not in UNCONDITIONAL_FORMS:
                l2 = self._lam(params, "lambda2")
                l4 = self._lam(params, "lambda4")
                dl2 = self._lam(params, "delta_lambda2")
                linf = self._lam(params, "lambda_inf")
                add(linf, eps)
                linf2 = linf * linf
                y_sqs = (0.0,) if self.ymax <= 0 else (0.0, self.ymax * self.ymax)
                for y_sq in y_sqs:
                    l2Y = l2 + dl2 * self._cast(y_sq)
                    if self.enforce_small_b:
                        add(l2Y, m)
                    add(self._cast(3.0) * linf2 * l4 + l2Y**3, m)

            if not vals:
                return None
            return tf.concat(vals, axis=0), tf.concat(lows, axis=0)

        def compute_nll_penalty(self, params, observables):
            zero = self._cast(0.0)
            eps = self._cast(self.eps)
            m = self._cast(self.margin)  # positive cushion: enforce coeff ≥ margin
            two = self._cast(2.0)
            three = self._cast(3.0)
            four = self._cast(4.0)
            nine = self._cast(9.0)
            ten = self._cast(10.0)
            twenty = self._cast(20.0)
            thirtysix = self._cast(36.0)
            # sqrt(0) has an infinite derivative; floor every radicand.
            tiny = self._cast(1e-30)

            def relu2(x):  # hinge: 0 if x ≤ 0 else x²
                return tf.square(tf.maximum(zero, x))

            # Each damping condition "coeff ≥ 0" is enforced as "coeff ≥ margin"
            # (margin=0 → bare ≥0): relu2(margin - coeff). The soft wall's gradient
            # vanishes at the knee, so a weakly-constrained coeff settles a hair
            # past it; the margin moves the knee so the equilibrium stays damping.
            def wall(coeff):  # coeff ≥ margin
                return relu2(m - coeff)

            pens = []

            # ---- CS-side γ_ν^NP damping: P(u)=λ2_ν·u+λ4_ν·u²+λ6_ν·u³ ≥ 0 ∀u≥0.
            # An UNCONDITIONAL_FORMS side (the *_abs fold) is skipped ENTIRELY, λ∞_ν
            # floor included: γ_ν ≤ margin holds for every λ, of either sign.
            if self._np_model_nu not in UNCONDITIONAL_FORMS:
                l2nu = self._lam(params, "lambda2_nu")
                l4nu = self._lam(params, "lambda4_nu")
                linfnu = self._lam(params, "lambda_inf_nu")
                pens.append(relu2(eps - linfnu))  # λ∞_ν > 0 (saturation-scale regime)
                if self._np_model_nu == "tanh_2":
                    pens.append(wall(l4nu))  # λ4_ν ≥ margin  (large-b leading)
                    if self.enforce_small_b:
                        pens.append(wall(l2nu))  # λ2_ν ≥ margin  (small-b turn-on)
                else:  # tanh_6 — λ6_ν only exists in the tanh_6 vocabulary
                    l6nu = self._lam(params, "lambda6_nu")
                    pens.append(wall(l6nu))  # λ6_ν ≥ margin  (large-b leading)
                    if self.enforce_small_b:
                        # small-b turn-on AND interior dip, unified: the minimum over
                        # u≥0 of the argument-per-u P(u)/u = λ2_ν + λ4_ν·u + λ6_ν·u² is
                        #   minP = λ2_ν − relu(−λ4_ν)²/(4·λ6_ν)   (interior vertex if
                        # λ4_ν<0, else the u=0 endpoint λ2_ν). This is in λ2_ν
                        # (coefficient) units, so wall(minP) applies ONE natural-unit
                        # cushion to whichever of turn-on/interior binds. At margin=0 ≡
                        # the old (division-free) small-b + discriminant pair (which was
                        # this × 4λ6_ν). C¹ in λ4_ν.
                        # FLOOR the denominator (same guard as l6_safe in the sigmoid
                        # branch): λ6_ν = 0 — the tanh_2 limit, e.g. a fit with λ6_ν
                        # FROZEN at 0 — otherwise gives 0/0 = NaN whenever λ4_ν ≥ 0
                        # (relu2(-λ4_ν) = 0), poisoning the whole penalty. With the floor
                        # it reduces to the u=0 endpoint minP = λ2_ν, which is the right
                        # answer there: P(u) = λ2_ν·u + λ4_ν·u² ≥ 0 ∀u≥0 already.
                        minP = l2nu - relu2(-l4nu) / (four * tf.maximum(l6nu, eps))
                        pens.append(wall(minP))
                    else:
                        # smallb=0: leave the u=0 turn-on free (λ2_ν floats either
                        # sign); keep only the interior-discriminant condition
                        # (scaled form + m).
                        pens.append(relu2(relu2(-l4nu) - four * l2nu * l6nu + m))

            # ---- TMD-side F_eff damping: Q(u)=λ2_Y+B·u+λ6·u² ≥ 0 ∀u≥0, evaluated
            # at the binding |Y| extremes. cubic ≡ 3·λ∞²·B = 3·λ∞²·λ4 + λ2_Y³, so
            # all conditions stay division-free (multiply through by 3·λ∞² > 0).
            # Skipped entirely (λ∞ floor included) for an UNCONDITIONAL_FORMS side.
            if self._np_model not in UNCONDITIONAL_FORMS:
                l2 = self._lam(params, "lambda2")
                l4 = self._lam(params, "lambda4")
                dl2 = self._lam(params, "delta_lambda2")
                linf = self._lam(params, "lambda_inf")
                sigmoid = self._np_model == "tanh_6_sigmoid"
                # λ6 only exists in the tanh_6 vocabularies; read it only when used.
                l6 = (
                    self._lam(params, "lambda6")
                    if self._np_model in ("tanh_6", "tanh_6_sigmoid")
                    else None
                )
                if sigmoid:
                    b_cut = self._lam(params, "bT_cutoff")
                    w_cut = self._lam(params, "bT_cutoff_width")
                pens.append(relu2(eps - linf))  # λ∞ > 0
                # λ6 ≥ margin (large-b leading coefficient; for the sigmoid it also
                # keeps the anti-damping window BOUNDED, so the asymptotic 2λ∞w<1
                # question never arises). OUTSIDE the |Y| loop: λ6 carries no Y
                # dependence, so charging it once per Y evaluation double-counted the
                # same condition and made this wall pull twice as hard as the CS-side
                # λ6_ν one. Same condition, same feasible set, correct weight.
                if l6 is not None:
                    pens.append(wall(l6))
                linf2 = linf * linf
                # ymax=0 -> Y=0 only. Listing 0.0 twice would double-charge the
                # identical condition and make this wall pull twice as hard.
                y_sqs = (0.0,) if self.ymax <= 0 else (0.0, self.ymax * self.ymax)
                for y_sq in y_sqs:
                    l2Y = l2 + dl2 * self._cast(y_sq)
                    if self._np_model == "tanh_2":
                        if self.enforce_small_b:
                            pens.append(wall(l2Y))  # λ2_Y ≥ margin (small-b turn-on)
                        cubic = three * linf2 * l4 + l2Y**3  # 3·λ∞²·B
                        pens.append(wall(cubic))  # B ≥ margin  (large-b leading)
                    elif sigmoid:  # tanh_6_sigmoid — see the module docstring
                        # λ2_Y ≥ margin is the small-b turn-on, where the turn-off is
                        # ≈1 and helps nothing. (λ6 ≥ margin is charged once, above
                        # the loop.)
                        if self.enforce_small_b:
                            pens.append(wall(l2Y))
                        B = l4 + l2Y**3 / (three * linf2)
                        # window {Q<0} = (u₋, u₊): roots of the quadratic Q. Divides by
                        # λ6 directly, as the tanh_6 branch's minQ does: λ6 is the
                        # leading coefficient of a form that needs it positive (and the
                        # wall enforces λ6 ≥ margin), so λ6 → 0 is a non-form, not a
                        # case to paper over. (This used to floor λ6 at
                        # LAMBDA_INF_FLOOR — an unrelated λ∞ constant, and 1e-3 is 10%
                        # of the λ6 = 0.01 in use, so it silently moved the roots.)
                        disc = B * B - four * l2Y * l6
                        sq = tf.sqrt(tf.maximum(disc, tiny))
                        u_hi = tf.maximum((-B + sq) / (two * l6), zero)
                        u_lo = tf.maximum((-B - sq) / (two * l6), zero)
                        # worst point of the tanh argument: d/db[b·Q(b²)]=0 is a
                        # quadratic in b² too.
                        d2 = tf.sqrt(tf.maximum(nine * B * B - twenty * l2Y * l6, tiny))
                        bs2 = tf.maximum((-three * B + d2) / (ten * l6), zero)
                        bs = tf.sqrt(tf.maximum(bs2, tiny))
                        a_min = bs * (l2Y + B * bs2 + l6 * bs2 * bs2) / linf
                        # anti = max over the window of (−tanh a) = −tanh(a_min), i.e.
                        # the PEAK anti-damping strength, gated to exactly 0 when
                        # a_min ≥ margin (no window worth charging for → this whole
                        # term vanishes, no spurious offset). Sign convention matches
                        # ``wall()``: +margin STRICTER, −margin slack. (It used to be
                        # a_min + m, i.e. looser with +m — the opposite of every other
                        # term, so one margin value pulled two ways.) The margin here
                        # is in tanh-argument units; see the module docstring on why
                        # that is not commensurate with the coefficient walls.
                        anti = tf.maximum(zero, -tf.tanh(a_min - m))
                        demand = two * linf * tf.sqrt(tf.maximum(u_hi, tiny)) * anti
                        supply = tf.math.softplus(
                            (tf.sqrt(tf.maximum(u_lo, tiny)) - b_cut) / w_cut
                        )
                        pens.append(relu2(demand - supply))
                    else:  # tanh_6 (λ6 ≥ margin is charged once, above the loop)
                        if self.enforce_small_b:
                            # small-b turn-on AND interior dip, unified: the minimum
                            # over u≥0 of Q(u) = λ2_Y + B·u + λ6·u²,
                            # B = λ4 + λ2_Y³/(3λ∞²), is
                            #   minQ = λ2_Y − relu(−B)²/(4·λ6)  (interior vertex if
                            # B<0, else the u=0 endpoint λ2_Y). In λ2_Y (coefficient)
                            # units, so wall(minQ) applies ONE natural-unit cushion to
                            # whichever of turn-on/interior binds. At margin=0 ≡ the
                            # old small-b + discriminant pair (which was this × 4λ6).
                            # C¹ in B.
                            B = l4 + l2Y**3 / (three * linf2)
                            # floor the denominator, same reason as minP on the CS side:
                            # λ6 = 0 (the tanh_2 limit / a frozen-at-0 λ6) would give
                            # 0/0 = NaN whenever B ≥ 0.
                            minQ = l2Y - relu2(-B) / (four * tf.maximum(l6, eps))
                            pens.append(wall(minQ))
                        else:
                            # smallb=0: leave the u=0 turn-on free (λ2_Y floats); keep
                            # only the interior discriminant (scaled, division-free
                            # form + m).
                            cubic = three * linf2 * l4 + l2Y**3
                            bound = thirtysix * linf2 * linf2 * l2Y * l6
                            pens.append(relu2(relu2(-cubic) - bound + m))

            # Both sides unconditional (e.g. tanh_6_abs / tanh_6_abs) ⇒ no terms at
            # all; tf.add_n([]) would raise, so return an explicit zero.
            if not pens:
                return zero
            return tf.add_n(pens)

    return NPDampingWall


# PEP-562 lazy class resolution: rabbit's loader does
#     module = importlib.import_module(...); cls = getattr(module, class_name)
# so the classes are synthesised on first attribute access, keeping the module
# importable without rabbit / TF (matches np_monotonicity.py).
def __getattr__(name):
    if name == "NPDampingMapping":
        cls = _make_mapping_class()
        globals()["NPDampingMapping"] = cls
        return cls
    if name == "NPDampingWall":
        cls = _make_regularizer_class()
        globals()["NPDampingWall"] = cls
        return cls
    raise AttributeError(name)
