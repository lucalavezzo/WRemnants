#ifndef WREMNANTS_MUON_EFFICIENCIES_CVH_H
#define WREMNANTS_MUON_EFFICIENCIES_CVH_H

// Data/MC muon efficiency scale factor for the CVH-refit efficiency hole
// introduced by the bug fixed in WMass/cmssw PR#46, commit eb96caef.
//
// The bug: a garbage-aligned glued (double-sided) module in the real data
// tracker alignment -- TIB layer 2, detId 369141860, whose mono face is
// rotated ~0.75 rad from ideal -- corrupts the composite frame used to
// anchor the CVH fit, killing tracks that cross it. It is DATA-ONLY (MC is
// reconstructed with ideal geometry, which has no such pathology), so MC
// over-counts efficiency in a localized (eta,phi) cell.
//
// Measured by an A/B refit (with vs without the fix) on 2016G SingleMuon,
// ~3.4M events/arm. The loss is:
//   * confined to the module's phi stripe  0.80 < phi < 1.00 rad,
//   * a smooth dip in eta, core SF ~0.585 at eta ~ 0.35-0.40,
//   * charge-independent and pT-independent (checked 25-65 GeV).
// So the correction is a 1D function of eta, gated on phi. Applied as a
// per-muon multiplicative weight on MC.
//
// Values hard-coded from the Delta-eta = 0.05 measurement (phi 0.80-1.00,
// pT 25-65 GeV). Bins outside the measured region return 1.

namespace wrem {

// Per-muon SF: eps_data(buggy) / eps_MC, = 1 outside the affected cell.
inline double cvh_muon_sf(float eta, float phi) {
  // phi stripe of the pathological module (radians, phi in [-pi,pi])
  if (phi < 0.80f || phi >= 1.00f)
    return 1.0;
  // eta dip, 0.05-wide bins from 0.05 to 0.70; else no effect
  if (eta < 0.05f || eta >= 0.70f)
    return 1.0;
  const int ibin = static_cast<int>((eta - 0.05f) / 0.05f); // 0..12
  // SF per bin, index 0 -> [0.05,0.10) ... index 12 -> [0.65,0.70)
  static const double sf[13] = {
      0.982, // [0.05,0.10)
      0.933, // [0.10,0.15)
      0.910, // [0.15,0.20)
      0.823, // [0.20,0.25)
      0.718, // [0.25,0.30)
      0.645, // [0.30,0.35)
      0.585, // [0.35,0.40)  core
      0.606, // [0.40,0.45)
      0.653, // [0.45,0.50)
      0.777, // [0.50,0.55)
      0.847, // [0.55,0.60)
      0.955, // [0.60,0.65)
      0.980  // [0.65,0.70)
  };
  return sf[ibin];
}

// Event-level weight for a dimuon (Z) final state: both reconstructed muons
// must survive, so the per-muon SFs multiply.
inline double cvh_dimuon_sf(float eta1, float phi1, float eta2, float phi2) {
  return cvh_muon_sf(eta1, phi1) * cvh_muon_sf(eta2, phi2);
}

} // namespace wrem

#endif
