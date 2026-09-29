#ifndef WREMNANTS_MUON_EFFICIENCIES_CVH_H
#define WREMNANTS_MUON_EFFICIENCIES_CVH_H

#include <ROOT/RVec.hxx>
#include <algorithm>
#include <cmath>
#include <cstddef>

// Treatment of the CVH-refit efficiency holes.
//
// Two badly aligned modules in the real data tracker alignment kill the tracks
// that cross them, each in its own way:
//
//   * the glued (double-sided) module in TIB layer 2, detId 369141860, whose
//     mono face is rotated ~0.75 rad from ideal, corrupting the composite frame
//     that anchors the CVH fit (fixed in WMass/cmssw PR#46, commit eb96caef);
//
//   * the r-phi sensor on TID disk -2, detId 402666798, sitting 5.3 mm out of
//     its disk plane, which makes the CVH propagation abort mid-track (still
//     uncorrected).
//
// The effect is DATA-ONLY -- MC is reconstructed with ideal geometry and has no
// such pathology.
//
// Both affected regions are localized rectangles in (eta, phi'), measured
// directly from the single-muon CVH refit efficiency in Z->mumu events with
//   mz_dilepton.py --cvhEfficiencyHists,
// and written out by
// scripts/analysisTools/w_mass_13TeV/make_cvh_efficiency_sf.py. Everywhere
// outside them the data/MC agreement is at the 1e-5 level.
//
// phi' is the azimuth at the module rather than at the vertex: a muon of
// charge q and transverse momentum pt is bent by q*C/pt on its way out to the
// module radius r, with C = 0.3*B*r/2. Undoing it (see cvh_module_phi) makes
// the regions charge- and pt-independent, and is what keeps the sharp edges of
// the module from being smeared over the two charges.
//
// Two treatments are provided:
//
//  * cvh_muon_sf (--cvhBadModules sf, the histmaker default): the measured
//    SF = eps_data/eps_MC, applied as a per-muon weight on MC only. It has to
//    be paired with cvh_veto_leak below, which carries the mirror effect on the
//    muon that is meant *not* to be reconstructed -- currently done in
//    mw_with_mu_eta_pt.py only, the dilepton selections having no veto to leak
//    through.
//
//  * cvh_muon_in_bad_module (--cvhBadModules veto): a geometric veto. Events
//    with any muon crossing one of the regions are dropped from DATA AND MC
//    alike, so the affected phase space simply leaves the measurement and no
//    efficiency correction is needed at all. Because the decision is made from
//    the uncorrected (standard tracker fit) kinematics it is blind to whether
//    the CVH refit succeeded, hence identical in data and MC by construction.
//    Costs ~1.3% of the muons, and is the cross-check on the SF treatment.

namespace wrem {

// Bending to the module radius, C = 0.3*B*r/2 in rad*GeV, i.e. C = 0.57*r with
// r in metres for B = 3.8 T. It belongs to the module, not to the analysis, so
// each region below carries its own value (bendC), taken from the radius at
// which its module sits: rho = 34.96 cm for the TIB one (C = 0.199) and
// 36.97 cm for the TID one (C = 0.211). make_cvh_efficiency_sf.py --fitBending
// measures the same quantity from the charge splitting of each hole and agrees
// with both within ~1 sigma, but only usefully constrains the deep TIB one.
//
// This constant is just the fallback for callers that ask for a module phi
// without naming a region; the regions themselves never use it.
inline constexpr double cvh_bending_C = 0.199;

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
// row-major in [ieta][iphi]. phi' is built with this region's own bendC.
struct CvhSFRegion {
  double etaLow, etaBinWidth;
  int nEta;
  double phiLow, phiBinWidth;
  int nPhi;
  double bendC;
  const double *sf;
};

// The same regions as plain (eta, phi') rectangles, for the geometric veto.
struct CvhModuleRegion {
  double etaLow, etaHigh, phiLow, phiHigh, bendC;
};

// clang-format off
// BEGIN GENERATED -- do not edit by hand, see make_cvh_efficiency_sf.py

// hotspot1 (TIB-L2, detId 369141860)
//   eta in [0.00, 0.75), phi' in [0.698, 1.047), C = 0.199 (r = 35 cm)
inline const double cvh_sf_region0[15 * 4] = {
    1.000, 0.991, 0.987, 1.000, // eta ~ +0.025
    1.000, 0.973, 0.975, 1.000, // eta ~ +0.075
    1.000, 0.933, 0.940, 1.000, // eta ~ +0.125
    1.000, 0.892, 0.886, 0.995, // eta ~ +0.175
    0.995, 0.800, 0.770, 0.990, // eta ~ +0.225
    0.993, 0.744, 0.635, 0.978, // eta ~ +0.275
    0.993, 0.690, 0.503, 0.960, // eta ~ +0.325
    0.993, 0.717, 0.419, 0.943, // eta ~ +0.375
    0.995, 0.771, 0.394, 0.928, // eta ~ +0.425
    0.997, 0.860, 0.490, 0.926, // eta ~ +0.475
    1.000, 0.925, 0.635, 0.931, // eta ~ +0.525
    1.000, 0.971, 0.782, 0.950, // eta ~ +0.575
    1.000, 0.991, 0.912, 0.967, // eta ~ +0.625
    1.000, 1.000, 0.966, 0.989, // eta ~ +0.675
    1.000, 1.000, 0.993, 0.994, // eta ~ +0.725
};

// hotspot2 (TID-2, detId 402666798)
//   eta in [-1.80, -1.45), phi' in [5.061, 5.411), C = 0.211 (r = 37 cm)
inline const double cvh_sf_region1[7 * 4] = {
    1.000, 1.000, 1.000, 0.983, // eta ~ -1.775
    0.992, 0.990, 0.990, 0.922, // eta ~ -1.725
    0.986, 0.968, 0.964, 0.872, // eta ~ -1.675
    0.958, 0.953, 0.939, 0.887, // eta ~ -1.625
    0.931, 0.923, 0.885, 0.902, // eta ~ -1.575
    0.965, 0.962, 0.943, 0.956, // eta ~ -1.525
    1.000, 1.000, 0.991, 0.994, // eta ~ -1.475
};

// last column of each row is bendC = 0.3*B*r/2, the bending to that module
inline const CvhSFRegion cvh_sf_regions[] = {
    {0.0000, 0.0500, 15, 0.6981, 0.0873, 4, 0.199, cvh_sf_region0},
    {-1.8000, 0.0500, 7, 5.0615, 0.0873, 4, 0.211, cvh_sf_region1},
};

// Bounding rectangles of the same regions, used by the geometric veto.
inline const CvhModuleRegion cvh_bad_modules[] = {
    {0.0000, 0.7500, 0.6981, 1.0472, 0.199}, // hotspot1 (TIB-L2, detId 369141860)
    {-1.8000, -1.4500, 5.0615, 5.4105, 0.211}, // hotspot2 (TID-2, detId 402666798)
};

// END GENERATED
// clang-format on

// True if the muon crosses one of the pathological modules. Meant to be fed
// the UNCORRECTED (standard tracker fit) kinematics, which exist whether or
// not the CVH refit succeeded; that is what makes the decision identical in
// data and MC. phi may be given in [-pi,pi] or [0,2*pi).
inline bool cvh_muon_in_bad_module(float eta, float phi, int charge, float pt) {
  for (const auto &r : cvh_bad_modules) {
    if (eta < r.etaLow || eta >= r.etaHigh)
      continue; // cheap test first, phi' costs a division per region
    const float phiModule = cvh_module_phi(phi, charge, pt, r.bendC);
    if (phiModule >= r.phiLow && phiModule < r.phiHigh)
      return true;
  }
  return false;
}

// Per-muon mask over a muon collection, for use in an event-level veto.
// Templated on the charge container only, whose element type varies with the
// nanoAOD version (Muon_charge is not necessarily an RVec<int>).
template <typename Vcharge>
inline ROOT::VecOps::RVec<bool> cvh_muons_in_bad_module(
    const ROOT::VecOps::RVec<float> &eta, const ROOT::VecOps::RVec<float> &phi,
    const Vcharge &charge, const ROOT::VecOps::RVec<float> &pt) {
  ROOT::VecOps::RVec<bool> out(eta.size(), false);
  for (std::size_t i = 0; i < eta.size(); ++i)
    out[i] = cvh_muon_in_bad_module(eta[i], phi[i], static_cast<int>(charge[i]),
                                    pt[i]);
  return out;
}

// Per-muon SF: eps_data/eps_MC, = 1 outside the affected regions. Apply to MC
// only, as a multiplicative weight, once per reconstructed muon. Alternative
// to the veto above, not to be combined with it.
inline double cvh_muon_sf(float eta, float phi, int charge, float pt) {
  for (const auto &r : cvh_sf_regions) {
    // floor, not truncation: the regions can sit at negative eta
    const int ieta =
        static_cast<int>(std::floor((eta - r.etaLow) / r.etaBinWidth));
    if (ieta < 0 || ieta >= r.nEta)
      continue;
    // each module is unbent to its own radius, so phi' is per region
    const float phiModule = cvh_module_phi(phi, charge, pt, r.bendC);
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

// ---------------------------------------------------------------------------
// The same holes seen from the other side: leakage into the W selection.
//
// A dimuon event enters the single-muon (W) selection when the second muon
// escapes the veto. In data that also happens when its CVH refit fails -- a
// muon with charge = -99 is not in 'vetoMuons' -- while in MC, with ideal
// geometry, it cannot. So DY, and any other two-muon process, leaks into the W
// selection in data only.
//
// The migration is simulated rather than reweighted. An MC event with exactly
// two veto muons, one of them inside a hotspot, is split into its two outcomes:
//
//   * the muon survives, probability s = cvh_muon_sf. Two veto muons remain,
//     the event stays out of the W selection and is dropped, as before.
//   * the muon is lost, probability 1 - s. It is removed from 'vetoMuons' and
//     the event enters the W selection carrying weight 1 - s.
//
// Only the second branch has to be materialised, the first one being rejected
// by the 'exactly one veto muon' requirement anyway. This is the deterministic
// (Rao-Blackwellised) form of throwing a random number per event: same
// expectation, no Monte Carlo noise on top, and the response to a variation of
// s stays smooth and exact instead of moving discrete events across a
// threshold.
//
// Three properties are worth spelling out.
//
//  * The parent population is the right one. These are ordinary, well measured
//    muons that happen to cross one module, not the soft or out-of-acceptance
//    muons that fail the veto for reconstruction reasons, and they carry the
//    corresponding kinematics, isolation and recoil.
//
//  * Removing the muon from 'vetoMuons' is precisely what a failed refit does
//    in data, so everything built on that collection -- the good muon
//    selection, jet cleaning -- follows on its own. The muon stays in the
//    PF/DeepMET, again as in data, since the refit does not touch the PF
//    reconstruction: these events therefore have Z-like MET and sit low in mT,
//    unlike the genuinely unreconstructed ones.
//
//  * The weight is bounded by 1 - s <= 0.61, the deepest cell of the map,
//    instead of the O(10)-O(100) up-weighting needed to express the same
//    migration as a multiple of the ordinary anti-veto leakage; and no MC-truth
//    veto efficiency map enters any more, so the veto SF variations stay
//    faithful.
//
// The survival factor of the muon that is kept is deliberately not included
// here: it is already applied to the good muon by cvh_muon_sf, via
// weight_cvhSF.
//
// No uncertainty is assigned to the refit SF map itself, here or in
// cvh_muon_sf.
struct CvhVetoLeak {
  int index = -1;         // muon removed from the veto collection, -1 if none
  double weight = 1.0;    // probability of the branch that is kept
  bool ambiguous = false; // both veto muons in a hotspot, see below
};

// Returns the muon to drop from the veto collection and the weight of that
// branch. Events with a number of veto muons other than two are left untouched
// (weight 1): with one there is nothing to split, and with none or more than
// two, losing a single muon cannot bring the event into the W selection, so the
// downstream 'exactly one veto muon' filter discards them either way.
//
// An event with both veto muons in a hotspot would need two branches, one per
// choice of surviving muon; only the one with the larger loss probability is
// kept. Two hotspot muons are down by ~1e-4 with respect to one, so the
// neglected branch is a per-cent correction to a per-mille background. It is
// flagged in 'ambiguous' so the rate can be checked directly.
//
// Feed it the corrected kinematics: 'vetoMuons' already requires
// Muon_correctedCharge != -99, so they are defined for every muon it selects,
// and they are what cvh_muon_sf is applied to elsewhere. charge and pt only
// serve to undo the bending out to the module radius.
template <typename Vcharge>
inline CvhVetoLeak cvh_veto_leak(const ROOT::VecOps::RVec<bool> &vetoMuons,
                                 const ROOT::VecOps::RVec<float> &eta,
                                 const ROOT::VecOps::RVec<float> &phi,
                                 const Vcharge &charge,
                                 const ROOT::VecOps::RVec<float> &pt) {
  CvhVetoLeak out;

  int n = 0, i0 = -1, i1 = -1;
  for (std::size_t i = 0; i < vetoMuons.size(); ++i) {
    if (!vetoMuons[i])
      continue;
    if (++n == 1)
      i0 = static_cast<int>(i);
    else if (n == 2)
      i1 = static_cast<int>(i);
    else
      return out; // more than two, losing one is not enough
  }
  if (n != 2)
    return out;

  const double p0 =
      1.0 - cvh_muon_sf(eta[i0], phi[i0], static_cast<int>(charge[i0]), pt[i0]);
  const double p1 =
      1.0 - cvh_muon_sf(eta[i1], phi[i1], static_cast<int>(charge[i1]), pt[i1]);
  if (p0 <= 0.0 && p1 <= 0.0)
    return out; // neither muon in a hotspot, stays a dimuon event

  out.ambiguous = (p0 > 0.0 && p1 > 0.0);
  out.index = (p0 >= p1) ? i0 : i1;
  out.weight = std::max(p0, p1);
  return out;
}

// The veto collection with the lost muon taken out. index < 0 leaves it alone.
inline ROOT::VecOps::RVec<bool>
cvh_drop_muon(const ROOT::VecOps::RVec<bool> &mask, int index) {
  if (index < 0)
    return mask;
  ROOT::VecOps::RVec<bool> out(mask);
  out[index] = false;
  return out;
}

} // namespace wrem

#endif
