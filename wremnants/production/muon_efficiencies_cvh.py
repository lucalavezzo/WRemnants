import narf.clingutils
from wums import logging

logger = logging.child_logger(__name__)

# Declares wrem::cvh_muon_sf / wrem::cvh_dimuon_sf, the hard-coded data/MC
# efficiency SF for the CVH-refit glued-module bug (TIB-L2 module 369141860,
# fixed in WMass/cmssw eb96caef). See the header for the measurement. MC only.
narf.clingutils.Declare('#include "muon_efficiencies_cvh.hpp"')


def define_cvh_weight(df, muons, out_col="weight_cvhSF"):
    """Define the per-event CVH efficiency weight = product of the per-muon SF
    over all reconstructed muons in the final state (1 for W, 2 for Z).

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
