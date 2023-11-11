#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

struct SpectraDstar {
  // Before selection
  OutputObj<TH1F> hBeforeInvDeltaMdstar{TH1F("hBeforeInvDeltaMdstar", "#Delta _{M} D*;#Delta _{M} D* (GeV/#it{c});entries", 100, 0.13, 0.16)};
  OutputObj<TH1F> hBeforeInvD0{TH1F("hBeforeInvD0", "#it{M} _{inv} D^{0};#it{M}_{inv} D^{0} (GeV/#it{c});entries", 500, 0., 5.)};

  // After selection
  OutputObj<TH1F> hInvDeltaMdstar{TH1F("hInvDeltaMdstar", "#Delta _{M} D*;#Delta _{M} D* (GeV/#it{c});entries", 100, 0.13, 0.16)};
  OutputObj<TH1F> hInvD0{TH1F("hInvD0", "#it{M}_{inv} D^{0};M _{inv} D^{0} (GeV/#it{c});entries", 500, 0., 5.)};

  OutputObj<TH1F> hPtDStar{TH1F("hPtDStar", "D* candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 6.)};
  OutputObj<TH1F> hPtD0{TH1F("hPtD0", "D^{0} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 6.)};

  OutputObj<TH1F> hEtaDstar{TH1F("hEtaDstar", "D* candidates;candidate #it{#eta}; entries", 100, -1.0, 1.0)};
  OutputObj<TH1F> hEtaD0{TH1F("hEtaD0", "D0 candidates;candidate #it{#eta}; entries", 100, -1.0, 1.0)};

  OutputObj<TH1F> hphiDstar{TH1F("hphiDstar", "D* candidates;candidate #it{#phi} Radian; entries", 100, 0., 2.0 * M_PI)};
  OutputObj<TH1F> hphiD0{TH1F("hphiD0", "D0 candidates;candidate #it{#phi} Radian; entries", 100, 0., 2.0 * M_PI)};

  OutputObj<TH1F> hPtPi{TH1F("hPtPi", "#pi_{s} prong;#it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hEtaPi{TH1F("hEtaPi", "#pi_{s} prong; #it{#eta}; entries ", 100, -1.0, 1.0)};
  OutputObj<TH1F> hphiPi{TH1F("hphiPi", "#pi_{s} prong; #it{#phi} Radian; entries", 100, 0., 2.0 * M_PI)};

  OutputObj<TH1F> hPtD0Prong0{TH1F("hPtD0Prong0", "D^{0} prong0;prong0 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0Prong1{TH1F("hPtD0Prong1", "D^{0} prong1;prong1 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hEtaD0Prong0{TH1F("hEtaD0Prong0", "D^{0} prong0;prong0 #it{#eta}", 100, -1.0, 1.0)};
  OutputObj<TH1F> hEtaD0Prong1{TH1F("hEtaD0Prong1", "D^{0} prong1;prong1 #it{#eta}", 100, -1.0, 1.0)};
  OutputObj<TH1F> hphiD0Prong0{TH1F("hphiD0Prong0", "D^{0} prong0;prong0 #it{#phi}", 100, 0., 2.0 * M_PI)};
  OutputObj<TH1F> hphiD0Prong1{TH1F("hphiD0Prong1", "D^{0} prong1;prong1 #it{#phi}", 100, 0., 2.0 * M_PI)};

  OutputObj<TH2F> hDCAxyD0Candidate{TH2F("hDCAxyD0Candidate", "D0 candidates;candidates #it{p}_{T} (GeV/#it{c};candidate DCA_{xy} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAzD0Candidate{TH2F("hDCAzD0Candidate", "D0 candidates;candidates #it{p}_{T} (GeV/#it{c};candidate DCA_{z} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAxySoftPi{TH2F("hDCAxySoftPi", "#pi_{s} prong;#it{p}_{T} (GeV/#it{c};DCA_{xy} (#mu m) ", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAzSoftPi{TH2F("hDCAzSoftPi", "#pi_{s} prong;#it{p}_{T} (GeV/#it{c};DCA_{z} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAxyD0Prong0{TH2F("hDCAxyD0Prong0", "D^{0} prong0;prong0 #it{p}_{T} (GeV/#it{c};DCA_{xy} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAzD0Prong0{TH2F("hDCAzD0Prong0", "D^{0} prong0;prong0 #it{p}_{T} (GeV/#it{c};DCA_{z} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAxyD0prong1{TH2F("hDCAxyD0prong1", "D^{0} prong1;prong1 #it{p}_{T} (GeV/#it{c};DCA_{xy} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};
  OutputObj<TH2F> hDCAzD0Prong1{TH2F("hDCAzD0Prong1", "D^{0} prong1;prong1 #it{p}_{T} (GeV/#it{c};DCA_{z} (#mu m)", 100, 0., 20., 200, -500.0, 500.0)};

  using SelectedDstar = soa::Join<aod::HfDstarCand, aod::HfSelDstarToD0Pi>;

  double massPi, massK, massD0;

  void init(InitContext& InitContext)
  {
    massPi = o2::analysis::pdg::MassPiPlus;
    massK = o2::analysis::pdg::MassKPlus;
    massD0 = o2::analysis::pdg::MassD0;
    hInvDeltaMdstar->Sumw2();
  }

  void process(aod::Tracks const&, aod::Hf2Prongs const&, aod::HfD0fromDstar const& rowsD0cand, SelectedDstar const& rowsDstarCand)
  {
    for (auto& rowDstarCand : rowsDstarCand) {
      auto statusDstar = rowDstarCand.isSelDstarToD0Pi();
      auto candD0 = rowsD0cand.iteratorAt(rowDstarCand.globalIndex());

      if (rowDstarCand.prongPi().sign() > 0.) {
        auto invDstarMass = rowDstarCand.dstarInvMass(std::array{massPi, massPi, massK});
        auto invD0Mass = candD0.d0m(std::array{massPi, massK});
        auto invDeltaMDstar = invDstarMass - invD0Mass;
        hBeforeInvDeltaMdstar->Fill(invDeltaMDstar);
        hBeforeInvD0->Fill(invD0Mass);
        if (statusDstar == 1) {
          // invMass
          hInvDeltaMdstar->Fill(invDeltaMDstar);
          hInvD0->Fill(invD0Mass);
          // pt Spectra
          hPtDStar->Fill(rowDstarCand.pt());
          hPtD0->Fill(candD0.d0pt());
          hPtPi->Fill(rowDstarCand.prongPi().pt());
          hPtD0Prong0->Fill(candD0.ptProng0());
          hPtD0Prong1->Fill(candD0.ptProng1());
          // eta
          hEtaDstar->Fill(rowDstarCand.eta());
          hEtaD0->Fill(candD0.d0eta());
          hEtaPi->Fill(rowDstarCand.prongPi().eta());
          hEtaD0Prong0->Fill(candD0.prong0().eta());
          hEtaD0Prong1->Fill(candD0.prong1().eta());
          // phi
          hphiDstar->Fill(rowDstarCand.phi());
          hphiD0->Fill(candD0.d0phi());
          hphiPi->Fill(rowDstarCand.prongPi().phi());
          hphiD0Prong0->Fill(candD0.prong0().phi());
          hphiD0Prong1->Fill(candD0.prong1().phi());
          // DCAxy
          hDCAxyD0Candidate->Fill(candD0.d0pt(), candD0.d0impactParameterXY());
          hDCAxyD0Prong0->Fill(candD0.ptProng0(), candD0.impactParameter0());
          hDCAxyD0prong1->Fill(candD0.ptProng1(), candD0.impactParameter1());
          hDCAxySoftPi->Fill(rowDstarCand.ptSoftpiProng(), rowDstarCand.impParamSoftPiProng());
          // DCAz ????
        }

      } else if (rowDstarCand.prongPi().sign() < 0.) {
        auto invAntiDstarMass = rowDstarCand.dstarInvMass(std::array{massPi, massK, massPi});
        auto invD0BarMass = candD0.d0m(std::array{massK, massPi});
        auto invDeltaMDstar = invAntiDstarMass - invD0BarMass;
        hBeforeInvDeltaMdstar->Fill(invDeltaMDstar);
        hBeforeInvD0->Fill(invD0BarMass);
        if (statusDstar == 1) {
          // invMass
          hInvDeltaMdstar->Fill(invDeltaMDstar);
          hInvD0->Fill(invD0BarMass);
          // pt Spectra
          hPtDStar->Fill(rowDstarCand.pt());
          hPtD0->Fill(candD0.d0pt());
          hPtPi->Fill(rowDstarCand.prongPi().pt());
          hPtD0Prong0->Fill(candD0.ptProng0());
          hPtD0Prong1->Fill(candD0.ptProng1());
          // eta
          hEtaDstar->Fill(rowDstarCand.eta());
          hEtaD0->Fill(candD0.d0eta());
          hEtaPi->Fill(rowDstarCand.prongPi().eta());
          hEtaD0Prong0->Fill(candD0.prong0().eta());
          hEtaD0Prong1->Fill(candD0.prong1().eta());
          // phi
          hphiDstar->Fill(rowDstarCand.phi());
          hphiD0->Fill(candD0.d0phi());
          hphiPi->Fill(rowDstarCand.prongPi().phi());
          hphiD0Prong0->Fill(candD0.prong0().phi());
          hphiD0Prong1->Fill(candD0.prong1().phi());
          // DCAxy
          hDCAxyD0Candidate->Fill(candD0.d0pt(), candD0.d0impactParameterXY());
          hDCAxyD0Prong0->Fill(candD0.ptProng0(), candD0.impactParameter0());
          hDCAxyD0prong1->Fill(candD0.ptProng1(), candD0.impactParameter1());
          hDCAxySoftPi->Fill(rowDstarCand.ptSoftpiProng(), rowDstarCand.impParamSoftPiProng());
          // DCAz ????
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SpectraDstar>(cfgc)};
}