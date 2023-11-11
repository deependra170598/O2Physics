#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"
#include "Framework/HistogramRegistry.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis;

struct PtDependentSpectraDstar {

  Configurable<int> numberOfPtBins{"numberOfPtBins", 10, "number of pT bins"};
  Configurable<std::vector<double>> binsPtForDstar{"binsPtForDstar", std::vector<double>{hf_cuts_dstar_to_pi_d0::vecBinsPt}, "pT bin limits for Dstar"};

  HistogramRegistry pTDependentRegistry{
    "pTDependentRegistry",
    {{"hDeltaInvMassDstarIntegrated", "#Delta _{M} D*;#Delta #it{M} _{inv} D* (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0.13, 0.16}}}},
     {"hInvMassD0Integrated", "#it{M}_{inv}D^{0};#it{M}_{inv}D^{0} (GeV/#it{c});entries ", {HistType::kTH1F, {{500, 0., 5.0}}}},
     {"hOverFlowPtBin", "#it{M}_{inv}D^{0};#it{M}_{inv}D^{0} (GeV/#it{c});entries ", {HistType::kTH1F, {{100, 0.13, 0.16}}}}}};

  using SelectedDstar = soa::Join<aod::HfDstarCand, aod::HfSelDstarToD0Pi>;

  double massPi, massK, massD0;
  static constexpr std::string_view HistSuffixs[25] = {"pTBin1", "pTBin2", "pTBin3", "pTBin4", "pTBin5", "pTBin6", "pTBin7", "pTBin8", "pTBin9", "pTBin10", "pTBin11", "pTBin12", "pTBin13", "pTBin14", "pTBin15", "pTBin16", "pTBin17", "pTBin18", "pTBin19", "pTBin20", "pTBin21"};

  void init(InitContext& InitContext)
  {
    massPi = o2::analysis::pdg::MassPiPlus;
    massK = o2::analysis::pdg::MassKPlus;
    massD0 = o2::analysis::pdg::MassD0;
    for (int i = 0; i < 50; i++) {
      if (i >= numberOfPtBins) {
        continue;
      }
      std::string HistSuffix = static_cast<std::string>(HistSuffixs[i]);
      std::string HistName = "hDeltaInvMassDstar" + HistSuffix;
      pTDependentRegistry.add(HistName.c_str(), "#Delta _{M} D*;#Delta #it{M} _{inv} D* (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0.13, 0.16}}}, true);
    }
  }

  void process(aod::Tracks const&, aod::Hf2Prongs const&, aod::HfD0fromDstar const& rowsD0cand, SelectedDstar const& rowsDstarCand)
  {
    for (auto& rowDstarCand : rowsDstarCand) {
      auto statusDstar = rowDstarCand.isSelDstarToD0Pi();
      auto candD0 = rowsD0cand.iteratorAt(rowDstarCand.globalIndex());
      if (statusDstar != 1) {
        continue;
      }

      double invDstarMass, invD0Mass, invDeltaMDstar, invAntiDstarMass, invD0BarMass;

      auto dstarpT = rowDstarCand.pt();
      auto pTBin = findBin(binsPtForDstar, dstarpT);
      LOGF(info, "pTBin: %d", pTBin);
      if (rowDstarCand.prongPi().sign() > 0.) {
        invDstarMass = rowDstarCand.dstarInvMass(std::array{massPi, massPi, massK});
        invD0Mass = candD0.d0m(std::array{massPi, massK});
        invDeltaMDstar = invDstarMass - invD0Mass;

      } else if (rowDstarCand.prongPi().sign() < 0.) {
        invAntiDstarMass = rowDstarCand.dstarInvMass(std::array{massPi, massK, massPi});
        invD0BarMass = candD0.d0m(std::array{massK, massPi});
        invDeltaMDstar = invAntiDstarMass - invD0BarMass;
      }
      switch (pTBin) {
        case 0:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin1"), invDeltaMDstar);
          break;
        case 1:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin2"), invDeltaMDstar);
          break;
        case 2:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin3"), invDeltaMDstar);
          break;
        case 3:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin4"), invDeltaMDstar);
          break;
        case 4:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin5"), invDeltaMDstar);
          break;
        case 5:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin6"), invDeltaMDstar);
          break;
        case 6:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin7"), invDeltaMDstar);
          break;
        case 7:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin8"), invDeltaMDstar);
          break;
        case 8:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin9"), invDeltaMDstar);
          break;
        case 9:
          pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin10"), invDeltaMDstar);
          break;
        /*case 10:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin11"),invDeltaMDstar);
            break;
        case 11:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin12"),invDeltaMDstar);
            break;
        case 12:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin13"),invDeltaMDstar);
            break;
        case 13:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin14"),invDeltaMDstar);
            break;
        case 14:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin15"),invDeltaMDstar);
            break;
        case 15:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin16"),invDeltaMDstar);
            break;
        case 16:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin17"),invDeltaMDstar);
            break;
        case 17:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin18"),invDeltaMDstar);
            break;
        case 18:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin19"),invDeltaMDstar);
            break;
        case 19:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin20"),invDeltaMDstar);
            break;
        case 20:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin21"),invDeltaMDstar);
            break;
        case 21:
            pTDependentRegistry.fill(HIST("hDeltaInvMassDstarpTBin22"),invDeltaMDstar);
            break;*/
        default:
          pTDependentRegistry.fill(HIST("hOverFlowPtBin"), invDeltaMDstar);
          break;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PtDependentSpectraDstar>(cfgc)};
}