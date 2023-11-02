#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;


struct SpectraDstar{
    // Before selection
    OutputObj<TH1F> hBeforeInvDeltaMdstar{TH1F("hBeforeInvDeltaMdstar","#Delta _{M} D*;#Delta _{M} D* (GeV/#it{c});entries",500,0.,5.)};
    OutputObj<TH1F> hBeforeInvD0{TH1F("hBeforeInvD0","#Delta _{M} D^{0};#Delta _{M} D^{0} (GeV/#it{c});entries",500,0.,5.)};

    // After selection
    OutputObj<TH1F> hPtDStar{TH1F("hPtDStar", "D* candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 6.)};
    OutputObj<TH1F> hPtD0{TH1F("hPtD0", "D^{0} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 6.)};

    OutputObj<TH1F> hEtaDstar{TH1F("hEtaDstar","D* candidates;candidate #it{#eta}; entries",100,-1.0,1.0)};
    OutputObj<TH1F> hEtaD0{TH1F("hEtaD0","D0 candidates;candidate #it{#eta}; entries",100,-1.0,1.0)};

    OutputObj<TH1F> hphiDstar{TH1F("hphiDstar","D* candidates;candidate #it{#phi}; entries",100,-1.0*M_PI,1.0*M_PI)};
    OutputObj<TH1F> hphiD0{TH1F("hphiD0","D0 candidates;candidate #it{#phi}; entries",100,-1.0*M_PI,1.0*M_PI)};
    // OutputObj<TH1F> hPtPi{TH1F("hPtPi", "#pi candidates;#it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
    // OutputObj<TH1F> hPtD0Prong0{TH1F("hPtD0Prong0", "D^{0} candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
    // OutputObj<TH1F> hPtD0Prong1{TH1F("hPtD0Prong1", "D^{0} candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
    OutputObj<TH1F> hInvDeltaMdstar{TH1F("hInvDeltaMdstar","#Delta _{M} D*;#Delta _{M} D* (GeV/#it{c});entries",500,0.,5.)};
    OutputObj<TH1F> hInvD0{TH1F("hInvD0","#Delta _{M} D^{0};#Delta _{M} D^{0} (GeV/#it{c});entries",500,0.,5.)};

    
    using SelectedDstar= soa::Join<aod::HfDstarCand,aod::HfSelDstarToD0Pi>;

    double massPi, massK, massD0;

    void init (InitContext & InitContext){
        massPi = o2::analysis::pdg::MassPiPlus;
        massK = o2::analysis::pdg::MassKPlus;
        massD0 = o2::analysis::pdg::MassD0;
    }

    void process(aod::Tracks const&,aod::Hf2Prongs const&, aod::HfD0fromDstar const& rowsD0cand,SelectedDstar const& rowsDstarCand){
        for(auto & rowDstarCand:rowsDstarCand){
            auto statusDstar = rowDstarCand.isSelDstarToD0Pi();
            auto candD0 = rowsD0cand.iteratorAt(rowDstarCand.globalIndex());

            
            if(rowDstarCand.prongPi().sign()> 0.){
                auto invDstarMass = rowDstarCand.dstarInvMass(std::array{massPi,massPi,massK});
                auto invD0Mass = candD0.d0m(std::array{massPi,massK});
                auto invDeltaMDstar = invDstarMass - invD0Mass;
                hBeforeInvDeltaMdstar->Fill(invDeltaMDstar);
                hBeforeInvD0->Fill(invD0Mass);
                if(statusDstar == 1){
                    //invMass
                    hInvDeltaMdstar->Fill(invDeltaMDstar);
                    hInvD0->Fill(invD0Mass);
                    //pt Spectra
                    hPtDStar->Fill(rowDstarCand.pt());
                    hPtD0->Fill(candD0.d0pt());
                    //eta
                    hEtaDstar->Fill(rowDstarCand.eta());
                    hEtaD0->Fill(candD0.d0eta());
                    //phi
                    hphiDstar->Fill(rowDstarCand.phi());
                    hphiD0->Fill(candD0.d0phi());
                }

            }else if(rowDstarCand.prongPi().sign()< 0.){
                auto invAntiDstarMass = rowDstarCand.dstarInvMass(std::array{massPi,massK,massPi});
                auto invD0BarMass = candD0.d0m(std::array{massK,massPi});
                auto invDeltaMDstar = invAntiDstarMass - invD0BarMass;
                hBeforeInvDeltaMdstar->Fill(invDeltaMDstar);
                hBeforeInvD0->Fill(invD0BarMass);
                if(statusDstar == 1){
                    //invMass
                    hInvDeltaMdstar->Fill(invDeltaMDstar);
                    hInvD0->Fill(invD0BarMass);
                    //pt Spectra
                    hPtDStar->Fill(rowDstarCand.pt());
                    hPtD0->Fill(candD0.d0pt());
                    // eta
                    hEtaDstar->Fill(rowDstarCand.eta());
                    hEtaD0->Fill(candD0.d0eta());
                    //phi
                    hphiDstar->Fill(rowDstarCand.phi());
                    hphiD0->Fill(candD0.d0phi());
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