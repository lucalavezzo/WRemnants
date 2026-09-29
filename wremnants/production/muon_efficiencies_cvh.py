import narf.clingutils

# Declares wrem::cvh_muon_in_bad_module / wrem::cvh_muon_sf / wrem::cvh_veto_leak,
# the two possible treatments of the CVH-refit efficiency holes of the badly
# aligned modules (TIB-L2 detId 369141860, fixed in WMass/cmssw eb96caef, and
# TID-2 detId 402666798, still uncorrected). See the header for the measurement.
narf.clingutils.Declare('#include "muon_efficiencies_cvh.hpp"')


def apply_bad_module_veto(
    df, ptCut=15.0, etaCut=2.4, ptMargin=1.0, filter_name="badCvhModuleVeto"
):
    """Reject events with a muon crossing one of the pathological modules.

    Applied to DATA AND MC alike: the affected (eta, phi') rectangles are simply
    removed from the measurement, which makes the data-only CVH refit
    inefficiency irrelevant instead of correcting for it.

    The decision uses the uncorrected muon kinematics (Muon_pt/eta/phi/charge)
    rather than the corrected ones, because Muon_corrected* is an alias of the
    CVH refit and is undefined (charge = -99) exactly for the muons that fall in
    the holes -- in data those muons are absent from 'vetoMuons' altogether,
    so a veto built on them would fire in MC only.

    For the same reason the muon selection here is deliberately looser than the
    veto muon definition (only Muon_looseId and the kinematic acceptance, no
    global/tracker quality since that too is evaluated on corrected quantities):
    it has to be a superset of 'vetoMuons', so that every event whose selection
    could be changed by a failed refit is removed. ptMargin lowers the pt
    threshold to absorb the (sub-percent) difference between corrected and
    uncorrected pt at the veto threshold.
    """
    df = df.Define(
        "Muon_inBadCvhModule",
        "wrem::cvh_muons_in_bad_module(Muon_eta, Muon_phi, Muon_charge, Muon_pt)",
    )
    df = df.Filter(
        f"!ROOT::VecOps::Any(Muon_inBadCvhModule && Muon_looseId"
        f" && Muon_pt > {ptCut - ptMargin} && abs(Muon_eta) < {etaCut})",
        filter_name,
    )
    return df


def define_cvh_weight(df, muons, out_col="weight_cvhSF"):
    """Define the per-event CVH efficiency weight = product of the per-muon SF
    over all reconstructed muons in the final state (1 for W, 2 for Z).

    Alternative to apply_bad_module_veto, not to be combined with it.

    muons: iterable of (eta_expr, phi_expr, charge_expr, pt_expr) tuples, each a
    column name or an inline RDF expression, e.g.
    ("trigMuons_eta0", "trigMuons_phi0", "trigMuons_charge0", "trigMuons_pt0").
    charge and pt are needed to undo the track bending, i.e. to evaluate the SF
    map in the module azimuth phi' = phi - q*C/pt (see muon_efficiencies_cvh.hpp).
    Call on MC only; multiply the result into the nominal weight.
    """
    expr = "*".join(
        f"wrem::cvh_muon_sf({eta}, {phi}, {charge}, {pt})"
        for eta, phi, charge, pt in muons
    )
    df = df.Define(out_col, expr)
    return df, out_col


def define_cvh_veto_leak(df, veto_col="vetoMuons", prefix="cvhVetoLeak"):
    """Simulate the dimuon events that leak into the W selection in data because
    the second muon's CVH refit failed.

    Companion of define_cvh_weight: that one corrects the muons that have to be
    reconstructed, this one the muon that has to be missed. An MC event with
    exactly two veto muons, one of them in a hotspot, is split into the branch
    where the muon survives (dropped, as it keeps two veto muons) and the branch
    where it is lost, which is kept with weight 1 - SF_cvh and with that muon
    removed from the veto collection -- exactly the state a failed refit leaves
    the event in, in data. See the header for why this is done instead of
    up-weighting the events whose second muon was never reconstructed.

    Call on MC only, between select_veto_muons(..., nMuons=-1), i.e. with its
    'exactly one veto muon' filter disabled, and select_good_muons, which builds
    on the veto collection; then require Sum(vetoMuons) == 1 as usual. Multiply
    the returned weight column into the nominal weight.

    Returns (df, weight column name).
    """
    df = df.Define(
        prefix,
        f"wrem::cvh_veto_leak({veto_col}, Muon_correctedEta, Muon_correctedPhi,"
        " Muon_correctedCharge, Muon_correctedPt)",
    )
    weight_col = f"weight_{prefix}"
    df = df.Define(weight_col, f"{prefix}.weight")
    # index of the muon declared lost, -1 for the events that are left alone,
    # and the flag for the (negligible) two-hotspot case, kept as columns so
    # both can be histogrammed
    df = df.Define(f"{prefix}_index", f"{prefix}.index")
    df = df.Define(f"{prefix}_ambiguous", f"{prefix}.ambiguous")
    # not a self-referencing Redefine, which would be ambiguous to read
    df = df.Define(
        f"{prefix}_{veto_col}", f"wrem::cvh_drop_muon({veto_col}, {prefix}_index)"
    )
    df = df.Redefine(veto_col, f"{prefix}_{veto_col}")
    return df, weight_col
