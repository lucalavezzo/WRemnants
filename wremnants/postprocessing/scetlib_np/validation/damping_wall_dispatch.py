"""Smoke test for the per-side, per-form dispatch in ``np_damping_wall``.

Builds the real :class:`NPDampingWall` regularizer against a mock indata+model
(no fit) and asserts the relu² penalty is ~0 inside the physical (damping) region
and > 0 outside, for both the tanh_2 and tanh_6 forms on each side, plus the
construction-time guards. Covers:

  * tanh_2 (unchanged): λ4_ν ≥ 0 (CS leading), 3λ∞²λ4+λ2_Y³ ≥ 0 (TMD), λ2 turn-on;
    λ6 / λ6_ν inert (sign ignored).
  * tanh_6: λ6(_nu) ≥ 0 (leading), and the INTERIOR discriminant bound that lets
    λ4_ν / B dip negative while the b⁶ term keeps the form damping
    (λ4_ν² ≤ 4λ2_νλ6_ν ; cubic² ≤ 36λ∞⁴λ2_Yλ6). The contrast case shows a λ4_ν that
    tanh_2 forbids but tanh_6 allows.
  * mixed sides (CS tanh_6 / TMD tanh_2), ``smallb=0`` dropping the turn-on walls,
    ``hyp_tangent`` alias → tanh_2, unknown form → NotImplementedError, and the
    missing-model → ValueError guard.
  * the margin: a cushion lifts a coeff==0 tune off zero; the tanh_6_sigmoid window
    gate moves the SAME way as every other term (+margin stricter, −margin slack);
    and the λ6 wall is charged ONCE, not once per |Y| evaluation.
  * tanh_6_abs (the damping fold): that side contributes EXACTLY 0 for any λ, without
    silencing a walled partner side.

Standalone (mocks rabbit's BaseMapping interface), so it needs TF + rabbit but no
input HDF5:  python -m wremnants.postprocessing.scetlib_np.validation.damping_wall_dispatch
"""

import tensorflow as tf

from wremnants.postprocessing.scetlib_np import np_damping_wall as W
from wremnants.postprocessing.scetlib_np.params import ALL_PARAMS

ORDER = list(ALL_PARAMS)  # l2nu,l4nu,l6nu,linfnu, l2,l4,l6,dl2,linf


class _MockModel:
    def __init__(self, np_model, np_model_nu):
        self._param_order = tuple(ORDER)
        self.fit_forms = {"np_model": np_model, "np_model_nu": np_model_nu}


class _MockIndata:
    def __init__(self, model=None):
        self.systs = []
        self.channel_info = {}  # BaseMapping iterates this (empty = no-op)
        self.procs = []
        if model is not None:
            self.scetlib_np_param_model = model


def _make_wall(np_model, np_model_nu, smallb=True, margin=0.0):
    """A wall on a mock model. ``margin=0`` (NOT the module's UNDER-TEST
    ``NP_DAMPING_MARGIN``) by default: these checks are about the per-side,
    per-form DISPATCH, and the "physical -> 0" assertions only hold at bare
    positivity — a cushion m > 0 deliberately lifts a coeff==0 tune off zero
    (``wall(0) = m²``), which the dedicated margin check below covers."""
    indata = _MockIndata(_MockModel(np_model, np_model_nu))
    mapping = W.NPDampingMapping(indata, "key", smallb=smallb, margin=margin)
    wall = W.NPDampingWall(mapping, tf.float64)
    wall.set_expectations(tf.zeros([len(ORDER)], tf.float64), None)
    return wall


def _pen(wall, **overrides):
    base = dict(
        lambda2_nu=0.15,
        lambda4_nu=0.0,
        lambda6_nu=0.0,
        lambda_inf_nu=2.0,
        lambda2=0.4,
        lambda4=0.4,
        lambda6=0.0,
        delta_lambda2=0.0,
        lambda_inf=1.0,
        # shape constants (params.CONST_PARAMS): present because ORDER is
        # ALL_PARAMS, inert for every form walled here
        abs_margin=0.0,
        abs_margin_nu=0.0,
        bT_cutoff=2.0,
        bT_cutoff_width=0.2,
    )
    base.update(overrides)
    v = tf.constant([base[n] for n in ORDER], dtype=tf.float64)
    return float(wall.compute_nll_penalty(v, None).numpy())


def main():
    tol = 1e-12
    fails = []

    def chk(name, cond):
        print(("PASS" if cond else "FAIL"), name)
        if not cond:
            fails.append(name)

    # ---- tanh_2 (both sides) — current behaviour preserved
    w2 = _make_wall("tanh_2", "tanh_2")
    chk("tanh2 physical -> 0", _pen(w2) < tol)
    chk("tanh2 CS l4nu<0 -> penalty", _pen(w2, lambda4_nu=-0.1) > tol)
    chk("tanh2 CS l2nu<0 -> penalty", _pen(w2, lambda2_nu=-0.1) > tol)
    chk("tanh2 CS linfnu<0 -> penalty", _pen(w2, lambda_inf_nu=-0.1) > tol)
    chk("tanh2 TMD cubic<0 (l4=-1) -> penalty", _pen(w2, lambda4=-1.0) > tol)
    chk("tanh2 TMD l2Y<0 -> penalty", _pen(w2, lambda2=-0.1) > tol)
    chk("tanh2 l6nu<0 inert -> 0", _pen(w2, lambda6_nu=-5.0) < tol)
    chk("tanh2 l6<0 inert -> 0", _pen(w2, lambda6=-5.0) < tol)

    # ---- tanh_6 (both sides)
    w6 = _make_wall("tanh_6", "tanh_6")
    chk("tanh6 physical -> 0", _pen(w6, lambda6_nu=0.001, lambda6=0.05) < tol)
    chk("tanh6 CS l6nu<0 -> penalty", _pen(w6, lambda6_nu=-0.001, lambda6=0.05) > tol)
    # interior within discriminant: 4*0.15*0.01=6e-3 > (-0.05)^2=2.5e-3 -> OK
    chk(
        "tanh6 CS l4nu<0 within bound -> 0",
        _pen(w6, lambda4_nu=-0.05, lambda2_nu=0.15, lambda6_nu=0.01, lambda6=0.05)
        < tol,
    )
    chk(
        "tanh6 CS l4nu<0 beyond bound -> penalty",
        _pen(w6, lambda4_nu=-0.2, lambda2_nu=0.15, lambda6_nu=0.01, lambda6=0.05) > tol,
    )
    chk("contrast: l4nu=-0.05 forbidden under tanh2", _pen(w2, lambda4_nu=-0.05) > tol)
    chk("tanh6 TMD l6<0 -> penalty", _pen(w6, lambda6_nu=0.001, lambda6=-0.05) > tol)
    # interior: cubic(l4=-0.5)=-1.436 ; bound=36*l2Y*l6 ; l6=0.5 -> 7.2>1.436² OK
    chk(
        "tanh6 TMD l4<0 within bound (l6=0.5) -> 0",
        _pen(w6, lambda6_nu=0.001, lambda4=-0.5, lambda6=0.5) < tol,
    )
    chk(
        "tanh6 TMD l4<0 beyond bound (l6=0.05) -> penalty",
        _pen(w6, lambda6_nu=0.001, lambda4=-0.5, lambda6=0.05) > tol,
    )

    # ---- mixed sides
    wm = _make_wall("tanh_2", "tanh_6")  # TMD tanh_2, CS tanh_6
    chk("mixed CS=tanh6 uses l6nu", _pen(wm, lambda6_nu=-0.001, lambda6=0.05) > tol)
    chk("mixed TMD=tanh2 uses cubic", _pen(wm, lambda6_nu=0.001, lambda4=-1.0) > tol)

    # ---- smallb=0 drops the turn-on walls
    w2s = _make_wall("tanh_2", "tanh_2", smallb=False)
    chk("smallb=0: l2nu<0 not penalised (l4nu=0)", _pen(w2s, lambda2_nu=-0.1) < tol)
    chk("smallb=1: l2nu<0 penalised", _pen(w2, lambda2_nu=-0.1) > tol)

    # ---- alias resolution
    wa = _make_wall("tanh_2", "hyp_tangent")
    chk("alias hyp_tangent -> tanh_2", wa._np_model_nu == "tanh_2")

    # ---- the margin cushion: a coeff sitting AT 0 is no longer free, by m² per
    # term (this is why the checks above pin margin=0, and why the module default
    # NP_DAMPING_MARGIN > 0 is not a "physical -> 0" regime).
    w2m = _make_wall("tanh_2", "tanh_2", margin=5e-3)
    chk("margin>0 lifts a coeff==0 tune off zero", _pen(w2m) > tol)
    chk(
        "margin>0 leaves a comfortably-damping tune at 0",
        _pen(w2m, lambda4_nu=0.2, lambda2_nu=0.2) < tol,
    )

    # ---- tanh_6_sigmoid: the window gate must move the SAME way as every other
    # term (+margin stricter, −margin slack). It used to be inverted (a_min + m),
    # so one margin value pulled the coefficient walls and the window gate in
    # opposite directions. Probe a tune that HAS an anti-damping window (λ4 < 0
    # deep enough that B² > 4λ2_Y λ6) so the gate is the live term.
    win = dict(lambda6_nu=0.001, lambda2=0.05, lambda4=-0.5, lambda6=0.05)
    p_strict = _pen(_make_wall("tanh_6_sigmoid", "tanh_6", margin=+0.05), **win)
    p_bare = _pen(_make_wall("tanh_6_sigmoid", "tanh_6", margin=0.0), **win)
    p_slack = _pen(_make_wall("tanh_6_sigmoid", "tanh_6", margin=-0.05), **win)
    print(
        f"    sigmoid window penalty: margin−0.05 {p_slack:.6g} | 0 {p_bare:.6g} | "
        f"+0.05 {p_strict:.6g}"
    )
    chk("sigmoid window gate: +margin stricter than 0", p_strict > p_bare)
    chk("sigmoid window gate: −margin looser than 0", p_slack < p_bare)

    # ---- λ6 is charged ONCE, not once per |Y| evaluation (it has no Y dependence).
    # Compare against tanh_2, whose branch has no λ6 term at all: the tanh_6 penalty
    # at a tune where only the λ6 wall is active must equal ONE relu2, not two.
    # λ6 = 0.01 (below the margin so its wall is live) but NOT 0: both tanh_6 and
    # the sigmoid divide by λ6 directly, so λ6 = 0 is a NaN non-form, not a test case.
    w6m = _make_wall("tanh_6", "tanh_6", margin=0.05)
    only_l6 = dict(
        lambda2_nu=0.5,
        lambda4_nu=0.5,
        lambda6_nu=0.5,
        lambda2=0.5,
        lambda4=0.5,
        lambda6=0.01,
        delta_lambda2=0.0,
    )
    expect = (0.05 - 0.01) ** 2
    print(f"    λ6-only penalty: {_pen(w6m, **only_l6):.6g} (one relu2 = {expect:.6g})")
    chk(
        "λ6 wall charged once (not once per |Y|)",
        abs(_pen(w6m, **only_l6) - expect) < 1e-12,
    )

    # ---- tanh_6_abs: unconditionally damping ⇒ that side contributes EXACTLY 0
    # (its λ∞ floor included). The prediction is bounded by the fold for any λ, so
    # there is no infeasible region left to wall; see the module docstring.
    wabs = _make_wall("tanh_6_abs", "tanh_6_abs")
    wild = dict(
        lambda2=-3.0,
        lambda4=-0.61,
        lambda6=-0.5,
        delta_lambda2=-3.35,
        lambda_inf=-1.0,
        lambda2_nu=-0.5,
        lambda4_nu=-0.5,
        lambda6_nu=-0.5,
        lambda_inf_nu=-2.0,
    )
    chk("abs both sides: physical λ -> 0", _pen(wabs) == 0.0)
    chk("abs both sides: wildly unphysical λ -> 0", _pen(wabs, **wild) == 0.0)
    wmix = _make_wall("tanh_6_abs", "tanh_6")  # TMD folded, CS still walled
    chk(
        "abs TMD side inert (l4<0, l6<0 change nothing)",
        _pen(wmix, lambda6_nu=0.001, lambda4=-0.5, lambda6=-0.05)
        == _pen(wmix, lambda6_nu=0.001, lambda4=0.4, lambda6=0.016),
    )
    chk(
        "abs TMD side does not silence the CS wall",
        _pen(wmix, lambda6_nu=-0.001) > tol,
    )

    # ---- unknown form raises
    try:
        _make_wall("frac_2", "tanh_2")
        chk("unknown TMD form raises", False)
    except NotImplementedError:
        chk("unknown TMD form raises", True)

    # ---- missing model raises
    try:
        m = W.NPDampingMapping(_MockIndata(model=None), "key")
        W.NPDampingWall(m, tf.float64)
        chk("missing model raises", False)
    except ValueError:
        chk("missing model raises", True)

    print()
    if fails:
        print(f"{len(fails)} FAILURES:", fails)
        raise SystemExit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
