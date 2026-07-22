#ifndef WREMNANTS_MUON_EFFICIENCIES_CVH_H
#define WREMNANTS_MUON_EFFICIENCIES_CVH_H

#include <cmath>

// Data/MC muon efficiency scale factor for the CVH-refit efficiency holes.
//
// Some badly aligned modules in the real data tracker alignment corrupt the
// composite frame used to anchor the CVH fit and kill the tracks that cross
// them. The known case is the glued (double-sided) module in TIB layer 2,
// detId 369141860, whose mono face is rotated ~0.75 rad from ideal (fixed in
// WMass/cmssw PR#46, commit eb96caef). The effect is DATA-ONLY -- MC is
// reconstructed with ideal geometry and has no such pathology -- so MC
// over-counts the efficiency in a localized cell and has to be downweighted.
//
// The correction is a map of eps_data/eps_MC in (eta, phi'), measured directly
// from the single-muon CVH refit efficiency in Z->mumu events with
//   mz_dilepton.py --cvhEfficiencyHists,
// and written out by scripts/analysisTools/make_cvh_efficiency_sf.py. Two
// regions of the (eta,phi') plane are affected; everywhere else the data/MC
// agreement is at the 1e-5 level and the SF is exactly 1.
//
// phi' is the azimuth at the module rather than at the vertex: a muon of
// charge q and transverse momentum pt is bent by q*C/pt on its way out to the
// module radius r, with C = 0.3*B*r/2. Undoing it (see cvh_module_phi) makes
// the map charge- and pt-independent, and is what keeps the sharp edges of the
// module from being smeared over the two charges.

namespace wrem {

// Bending to the module radius, C = 0.3*B*r/2 in rad*GeV. Fitted from the
// charge splitting of the hole position, C = 0.214 +- 0.014, which corresponds
// to r = 38 cm, i.e. TIB layer 2 (layers at r = 25.5, 33.9, 41.9, 49.8 cm).
inline constexpr double cvh_bending_C = 0.214;

// Azimuth of the muon where it crosses the module, in [0, 2*pi).
// Accepts phi either in [-pi,pi] (nanoAOD) or in [0,2*pi).
inline float cvh_module_phi(float phi, int charge, float pt,
                            double bendC = cvh_bending_C) {
  constexpr float twopi = 2.f * static_cast<float>(M_PI);
  float p = phi - static_cast<float>(charge * bendC) / pt;
  while (p < 0.f)
    p += twopi;
  while (p >= twopi)
    p -= twopi;
  return p;
}

// A rectangular region of the (eta, phi') plane holding the measured SF map,
// row-major in [ieta][iphi].
struct CvhSFRegion {
  double etaLow, etaBinWidth;
  int nEta;
  double phiLow, phiBinWidth;
  int nPhi;
  const double *sf;
};

// clang-format off
// BEGIN GENERATED -- do not edit by hand, see make_cvh_efficiency_sf.py

// hotspot1 (TIB-L2, detId 369141860)
//   eta in [0.00, 0.75), phi' in [0.698, 1.047)
inline const double cvh_sf_region0[15 * 4] = {
    1.000, 0.991, 0.987, 1.000, // eta ~ +0.025
    1.000, 0.973, 0.975, 1.000, // eta ~ +0.075
    1.000, 0.934, 0.940, 1.000, // eta ~ +0.125
    1.000, 0.892, 0.886, 0.995, // eta ~ +0.175
    0.994, 0.800, 0.771, 0.989, // eta ~ +0.225
    0.993, 0.744, 0.636, 0.977, // eta ~ +0.275
    0.992, 0.690, 0.504, 0.959, // eta ~ +0.325
    0.992, 0.717, 0.421, 0.941, // eta ~ +0.375
    0.995, 0.770, 0.396, 0.927, // eta ~ +0.425
    0.997, 0.859, 0.492, 0.926, // eta ~ +0.475
    1.000, 0.925, 0.636, 0.930, // eta ~ +0.525
    1.000, 0.971, 0.783, 0.949, // eta ~ +0.575
    1.000, 0.990, 0.912, 0.967, // eta ~ +0.625
    1.000, 1.000, 0.967, 0.989, // eta ~ +0.675
    1.000, 1.000, 0.993, 0.994, // eta ~ +0.725
};

// hotspot2 (uncorrected)
//   eta in [-1.80, -1.45), phi' in [5.061, 5.411)
inline const double cvh_sf_region1[7 * 4] = {
    1.000, 1.000, 1.000, 0.983, // eta ~ -1.775
    0.992, 0.990, 0.990, 0.922, // eta ~ -1.725
    0.986, 0.968, 0.964, 0.872, // eta ~ -1.675
    0.958, 0.953, 0.939, 0.887, // eta ~ -1.625
    0.931, 0.923, 0.885, 0.902, // eta ~ -1.575
    0.965, 0.962, 0.943, 0.956, // eta ~ -1.525
    1.000, 1.000, 0.991, 0.994, // eta ~ -1.475
};

inline const CvhSFRegion cvh_sf_regions[] = {
    {0.0000, 0.0500, 15, 0.6981, 0.0873, 4, cvh_sf_region0},
    {-1.8000, 0.0500, 7, 5.0615, 0.0873, 4, cvh_sf_region1},
};

// END GENERATED
// clang-format on

// Per-muon SF: eps_data/eps_MC, = 1 outside the affected regions. Apply to MC
// only, as a multiplicative weight, once per reconstructed muon.
inline double cvh_muon_sf(float eta, float phi, int charge, float pt) {
  const float phiModule = cvh_module_phi(phi, charge, pt);
  for (const auto &r : cvh_sf_regions) {
    // floor, not truncation: the regions can sit at negative eta
    const int ieta =
        static_cast<int>(std::floor((eta - r.etaLow) / r.etaBinWidth));
    if (ieta < 0 || ieta >= r.nEta)
      continue;
    const int iphi =
        static_cast<int>(std::floor((phiModule - r.phiLow) / r.phiBinWidth));
    if (iphi < 0 || iphi >= r.nPhi)
      continue;
    return r.sf[ieta * r.nPhi + iphi];
  }
  return 1.0;
}

// Event-level weight for a dimuon (Z) final state: both reconstructed muons
// must survive, so the per-muon SFs multiply.
inline double cvh_dimuon_sf(float eta1, float phi1, int charge1, float pt1,
                            float eta2, float phi2, int charge2, float pt2) {
  return cvh_muon_sf(eta1, phi1, charge1, pt1) *
         cvh_muon_sf(eta2, phi2, charge2, pt2);
}

} // namespace wrem

#endif
