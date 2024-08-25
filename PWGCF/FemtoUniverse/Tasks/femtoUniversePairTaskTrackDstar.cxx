#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"Dstar", "Track"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

/// Returns deltaPhi value within the range [-pi/2, 3/2*pi]
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

/// Returns deltaPhi value within the range [0, pi]
///
double wrapDeltaPhi0PI(double phiD, double phiDbar)
{
  double deltaPhi = 0.0;
  deltaPhi = RecoDecay::constrainAngle(phiDbar - phiD, 0.0);
  if (deltaPhi < 0.) {
    deltaPhi = deltaPhi + o2::constants::math::TwoPI;
  }
  if (deltaPhi > o2::constants::math::TwoPI) {
    deltaPhi = o2::constants::math::TwoPI - deltaPhi;
  }
  return deltaPhi;
}

struct femtoUniversePairTaskTrackDstar {

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaCombinedProton{"ConfNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCProton{"ConfNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};

    Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
    Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
  } ConfBothTracks;

  /// Particle 1 --- IDENTIFIED TRACK
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
    Configurable<int> ConfPDGCodeTrack{"ConfPDGCodeTrack", 2212, "Particle 2 - PDG code"};
    Configurable<int> ConfPIDTrack{"ConfPIDTrack", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int8_t> ConfTrackSign{"ConfTrackSign", 1, "Track sign"};
    Configurable<bool> ConfIsTrackIdentified{"ConfIsTrackIdentified", true, "Enable PID for the track"};
  } ConfTrack;

  /// Particle 2 --- Dstar/AntiDstar 
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodeDstar{"ConfPDGCodeDstar", 413, "Dstar - PDG code"};
    Configurable<int> ConfPDGCodeAntiDstar{"ConfPDGCodeAntiDstar", -413, "Anti Dstar - PDG code"};
    Configurable<float> ConfMinPtDstarAntiDstar{"ConfMinPtDstarAntiDstar", 3.0, "Dstar/AntiDstar sel. - min. pT"}; //need to change
    Configurable<float> ConfMaxPtDstarAntiDstar{"ConfMaxPtDstarAntiDstar", 8.0, "Dstar/AntiDstar sel. - max. pT"}; //need to change
    Configurable<float> ConfMinInvMassDstarAntiDstar{"ConfMinInvMassDstarAntiDstar", 0.13, "Dstar/AntiDstar sel. - min. invMass"}; // need to change
    Configurable<float> ConfMaxInvMassDstarAntiDstar{"ConfMaxInvMassDstarAntiDstar", 0.15, "Dstar/AntiDstar sel. - max. invMass"}; //need to change
  } ConfDstar;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfSignalRegionMin{"ConfSignalRegionMin", 0.144, "Min. inv. mass for Dstar/AntiDstar in the signal region"};
    Configurable<float> ConfSignalRegionMax{"ConfSignalRegionMax", 0.146, "Max. inv. mass for Dstar/AntiDstar in the signal region"};
    // Configurable<float> ConfMinInvMassLeftSB{"ConfMinInvMassLeftSB", 1.642, "Min. inv. mass for Dstar/AntiDstar in the left sideband region"};
    // Configurable<float> ConfMaxInvMassLeftSB{"ConfMaxInvMassLeftSB", 1.754, "Max. inv. mass for Dstar/AntiDstar in the left sideband region"};
    Configurable<float> ConfMinInvMassRightSB{"ConfMinInvMassRightSB", 0.147, "Min. inv. mass for Dstar/AntiDstar in the right sideband region"};
    Configurable<float> ConfMaxInvMassRightSB{"ConfMaxInvMassRightSB", 0.154, "Max. inv. mass for Dstar/AntiDstar in the right sideband region"};
  } ConfDstarAntiDstarSideBand; //need to change 

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits"}; 
  Configurable<bool> ConfUseAllDstar{"ConfUseAllDstar", false, "Include cand. which are both Dstar and Anti Dstar cand."};
  Configurable<bool> ConfUsePtCutForDstarAntiDstar{"ConfUsePtCutForDstarAntiDstar", false, "Include pT cut for Dstar/Anti Dstar in same and mixed processes."};
  Configurable<bool> ConfUseMassCutForDstarAntiDstar{"ConfUseMassCutForDstarAntiDstar", false, "Switch to save Dstar/Anti Dstar within declared inv. mass range"};


  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.ConfTrackSign));

  /// Partitions for particle 2
  Partition<FemtoFullParticles> partsAllDstar = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kDstar));
  //===============================================================================================================================================================================
  //________________________________________________Doubt____________________________________________________________________________________________________________________________
  Partition<FemtoFullParticles> partsOnlyDstarAntiDstar = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kDstar)) && ((aod::femtouniverseparticle::mLambda > 0.142f && aod::femtouniverseparticle::mLambda < 0.150f) || (aod::femtouniverseparticle::mAntiLambda > 0.142f && aod::femtouniverseparticle::mAntiLambda < 0.150f));

  /// Partition for Dstardaughters
  Partition<FemtoFullParticles> partsDstarChildren = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kDstarChild);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTrack;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kDstar, 0> trackHistoPartDstar;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDTrack;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarInvMassBins{"ConfDTempFitVarInvMassBins", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis ConfMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis ConfPtBins{"ConfPtBins", {360, 0., 36.}, "binning pT"};
  ConfigurableAxis ConfInvMassBins{"ConfInvMassBins", {500, 0., 5.0}, "binning inv. mass"};
  ConfigurableAxis ConfInvMassFinerBins{"ConfInvMassFinerBins", {120, 1.5848, 2.1848}, "finer binning of inv. mass"};

  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventFemtoCont;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventFemtoCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::same, femtoUniverseAngularContainer::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::mixed, femtoUniverseAngularContainer::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kDstar> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kDstar> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry registry{"registry",
                             {{"hInvMassDstar", ";#it{M}(#pi^{-}K^{+}#pi^{+}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {ConfInvMassBins}}},
                              {"hInvMassAntiDstar", ";#it{M}(#pi^{+}K^{-}#pi^{-}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {ConfInvMassBins}}},
                              {"hPtDstarCand", "3-prong candidates;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtDstar", "D^{s} cand.;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtAntiDstar", "#bar{D^{s}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtDstarAntiDstar", "#bar{D^{s}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPhiDstarCand", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
                              {"hPhiDstar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
                              {"hPhiAntiDstar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
                              {"hEtaDstarCand", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaDstar", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaAntiDstar", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hDecayLengthDstar", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hDecayLengthAntiDstar", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hPtDaughters", ";#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{300, 0., 12.}}}},
                              {"hSignDaughters", ";sign ;counts", {HistType::kTH1F, {{10, -2.5, 2.5}}}},
                              {"hbetaDaughters", "; p (GeV/#it{c}); TOF #beta", {HistType::kTH2F, {{300, 0., 15.}, {200, 0., 2.}}}},
                              {"hdEdxDaughters", "; p (GeV/#it{c}); TPC dE/dx (KeV/cm)", {HistType::kTH2F, {{300, 0., 15.}, {500, 0., 500.}}}},
                              {"hDCAxyDaughters", "; #it{DCA}_{xy} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}},
                              {"hDCAzDaughters", "; #it{DCA}_{z} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}}}};

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCProton -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedProton -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5

    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
        // if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.3) { // 0.0-0.3
      if (TMath::Abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < 3.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < 2.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCPion -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
          // if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.ConfPDGCodeTrack) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }


  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartDstar.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarInvMassBins, ConfBothTracks.ConfIsMC, ConfDstar.ConfPDGCodeDstar);
    if (!ConfTrack.ConfIsSame) {
      trackHistoPartTrack.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfBothTracks.ConfIsMC, ConfTrack.ConfPDGCodeTrack);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    sameEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);



    sameEventFemtoCont.setPDGCodes(ConfDstar.ConfPDGCodeDstar, ConfTrack.ConfPDGCodeTrack);
    sameEventAngularCont.setPDGCodes(ConfDstar.ConfPDGCodeDstar, ConfTrack.ConfPDGCodeTrack);
    mixedEventFemtoCont.setPDGCodes(ConfDstar.ConfPDGCodeDstar, ConfTrack.ConfPDGCodeTrack);
    mixedEventAngularCont.setPDGCodes(ConfDstar.ConfPDGCodeDstar, ConfTrack.ConfPDGCodeTrack);

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiCutMin.value, ConfCPRdeltaPhiCutMax.value, ConfCPRdeltaEtaCutMin.value, ConfCPRdeltaEtaCutMax.value, ConfCPRChosenRadii.value, ConfCPRPlotPerRadii.value);
    }

    vPIDTrack = ConfTrack.ConfPIDTrack.value;
    kNsigma = ConfBothTracks.ConfTrkPIDnSigmaMax.value;

    // Dstar/AntiDstar histograms
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassVsPt", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassVsPtFinerBinning", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassFinerBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassVsPtOnlyDstarAntiDstar", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaPhiSigSig", "SxS correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiDstarBgAntiDstarSig", "B(Dstar)x S(AntiDstar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiDstarSigAntiDstarBg", "S(Dstar)x B(AntiDstar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiBgBg", "BxB correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hPtProng1VsPtProng2VsPtProng3", "3-prong candidates;#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c};#it{p}_{T} (GeV/#it{c})", {HistType::kTH3F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}}); //need to check
    registry.add("hDeltaEtaDeltaPhi", "3-prong candidates;#Delta #eta;#Delta #varphi (rad)", {HistType::kTH2F, {{29, -2., 2.}, {29, 0.0, o2::constants::math::PI}}});
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }
  
  void processDstar(o2::aod::FDCollision& col, FemtoFullParticles&)
  {
    auto groupPartsOnlyDstarAntiDstar = partsOnlyDstarAntiDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // LOGF(info, " Size of only Dstar: %d", groupPartsOnlyDstarAntiDstar.size());
    auto groupPartsAllDstar = partsAllDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // LOGF(info,"Size of AllDstar: %d", groupPartsAllDstar.size() );
    auto groupPartsDstarChildren = partsDstarChildren->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // LOGF(info,"Size of Dstar Childs: %d", groupPartsDstarChildren.size() );
    

    // loop over all D mesons
    for (auto const& dstar : groupPartsAllDstar) {

      if (dstar.mLambda() > 0.0f) {
        registry.fill(HIST("hMassVsPt"), dstar.mLambda(), dstar.pt());
        registry.fill(HIST("hMassVsPtFinerBinning"), dstar.mLambda(), dstar.pt());
      }

      if (dstar.mAntiLambda() > 0.0f) {
        registry.fill(HIST("hMassVsPt"), dstar.mAntiLambda(), dstar.pt());
        registry.fill(HIST("hMassVsPtFinerBinning"), dstar.mAntiLambda(), dstar.pt());
      }

      registry.fill(HIST("hPtDstarCand"), dstar.pt());
      registry.fill(HIST("hPhiDstarCand"), dstar.phi());
      registry.fill(HIST("hEtaDstarCand"), dstar.eta());
    }

    // loop over Dstar/AntiDstar mesons (ONLY)
    for (auto const& dsAntids : groupPartsOnlyDstarAntiDstar) {

      registry.fill(HIST("hPtDstar"), dsAntids.pt());
      
      float massDstarParticle = dsAntids.mLambda();
      float massDstarAntiparticle = dsAntids.mAntiLambda(); 

      LOGF(info,"dsAntids mass D* (Piar Task): %f, massAntiD*: %f",massDstarParticle, massDstarAntiparticle);

      if (massDstarParticle > 0.0 && massDstarAntiparticle < 0.0) {
        registry.fill(HIST("hInvMassVsPtOnlyDstarAntiDstar"), massDstarParticle, dsAntids.pt());
        if (massDstarParticle > ConfDstar.ConfMinInvMassDstarAntiDstar && massDstarParticle < ConfDstar.ConfMaxInvMassDstarAntiDstar) {
          registry.fill(HIST("hInvMassDstar"), dsAntids.mLambda());
        }
        registry.fill(HIST("hPtDstar"), dsAntids.pt());
        registry.fill(HIST("hPhiDstar"), dsAntids.phi());
        registry.fill(HIST("hEtaDstar"), dsAntids.eta());
      }
      if (massDstarParticle < 0.0 && massDstarAntiparticle > 0.0) {
        registry.fill(HIST("hInvMassVsPtOnlyDstarAntiDstar"), massDstarAntiparticle, dsAntids.pt());
        if (massDstarAntiparticle > ConfDstar.ConfMinInvMassDstarAntiDstar && massDstarAntiparticle < ConfDstar.ConfMaxInvMassDstarAntiDstar) {
          registry.fill(HIST("hInvMassAntiDstar"), massDstarAntiparticle);
        }
        registry.fill(HIST("hPtAntiDstar"), dsAntids.pt());
        registry.fill(HIST("hPhiAntiDstar"), dsAntids.phi());
        registry.fill(HIST("hEtaAntiDstar"), dsAntids.eta());
      }
    }

    // loop over D mesons childen
    for (auto const& daughDstar : groupPartsDstarChildren) {
      registry.fill(HIST("hPtDaughters"), daughDstar.pt());
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackDstar, processDstar, "Enable processing Dstar", true);


template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsTrack, PartitionType groupPartsDstar, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto& dstarcandidate : groupPartsDstar) {
      trackHistoPartDstar.fillQA<isMC, false>(dstarcandidate);
    }

    if (!ConfTrack.ConfIsSame) {
      for (auto& track : groupPartsTrack) {
        if (ConfTrack.ConfIsTrackIdentified) {
          if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
            continue;
          }
        }
        trackHistoPartTrack.fillQA<isMC, false>(track);
      }
    }
    /// Now build the combinations
    for (auto& [track, dstarcandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsDstar))) {
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Set pT cut for Dstar/AntiDstar candidates
      if (ConfUsePtCutForDstarAntiDstar) {
        if (dstarcandidate.pt() < ConfDstar.ConfMinPtDstarAntiDstar && dstarcandidate.pt() > ConfDstar.ConfMaxPtDstarAntiDstar) {
          continue;
        }
      }
      // // Set inv. mass cut for Dstar/AntiDstar candidates
      if (ConfUseMassCutForDstarAntiDstar) {
        if ((dstarcandidate.mLambda() < ConfDstarAntiDstarSideBand.ConfSignalRegionMin && dstarcandidate.mLambda() > ConfDstarAntiDstarSideBand.ConfSignalRegionMax) || (dstarcandidate.mAntiLambda() < ConfDstarAntiDstarSideBand.ConfSignalRegionMin && dstarcandidate.mAntiLambda() > ConfDstarAntiDstarSideBand.ConfSignalRegionMax)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, dstarcandidate, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, dstarcandidate, parts)) {
        continue;
      }
      sameEventFemtoCont.setPair<isMC>(track, dstarcandidate, multCol, ConfBothTracks.ConfUse3D);
      sameEventAngularCont.setPair<isMC>(track, dstarcandidate, multCol, ConfBothTracks.ConfUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        FemtoFullParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsAllDstarAntiDstar = partsAllDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsOnlyDstarAntiDstar = partsOnlyDstarAntiDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    if (ConfUseAllDstar) {
      doSameEvent<false>(thegroupPartsTrack, thegroupPartsAllDstarAntiDstar, parts, col.magField(), col.multNtr());
    } else {
      doSameEvent<false>(thegroupPartsTrack, thegroupPartsOnlyDstarAntiDstar, parts, col.magField(), col.multNtr());
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackDstar, processSameEvent, "Enable processing same event", true);

  
/// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsTrack partition for the identified passed by the process function
  /// \param groupPartsDstar partition for Dstar meson passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsTrack, PartitionType groupPartsDstar, PartType parts, float magFieldTesla, int multCol)
  {

    for (auto& [track, dstarcandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsDstar))) {
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Set pT cut for Dstar/AntiDstar candidates
      if (ConfUsePtCutForDstarAntiDstar) {
        if (dstarcandidate.pt() < ConfDstar.ConfMinPtDstarAntiDstar && dstarcandidate.pt() > ConfDstar.ConfMaxPtDstarAntiDstar) {
          continue;
        }
      }
      // // Set inv. mass cut for Dstar/AntiDstar candidates
      if (ConfUseMassCutForDstarAntiDstar) {
        if ((dstarcandidate.mLambda() < ConfDstarAntiDstarSideBand.ConfSignalRegionMin && dstarcandidate.mLambda() > ConfDstarAntiDstarSideBand.ConfSignalRegionMax) || (dstarcandidate.mAntiLambda() < ConfDstarAntiDstarSideBand.ConfSignalRegionMin && dstarcandidate.mAntiLambda() > ConfDstarAntiDstarSideBand.ConfSignalRegionMax)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, dstarcandidate, parts, magFieldTesla, femtoUniverseContainer::EventType::mixed)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(track, dstarcandidate, multCol, ConfBothTracks.ConfUse3D);
      mixedEventAngularCont.setPair<isMC>(track, dstarcandidate, multCol, ConfBothTracks.ConfUse3D);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FDCollisions& cols,
                         FemtoFullParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsAllDstar = partsAllDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsOnlyDstarAntiDstar = partsOnlyDstarAntiDstar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsDstar.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      if (ConfUseAllDstar) {
        doMixedEvent<false>(groupPartsTrack, groupPartsAllDstar, parts, magFieldTesla1, multiplicityCol);
      } else {
        doMixedEvent<false>(groupPartsTrack, groupPartsOnlyDstarAntiDstar, parts, magFieldTesla1, multiplicityCol);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackDstar, processMixedEvent, "Enable processing mixed events", false);
};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackDstar>(cfgc),
  };
  return workflow;
}