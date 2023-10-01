#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace  o2::framework;

// Struct to applying Dstar selection cuts
struct HfCandidateSelctorDstar{
    // Produces<aod::HfSelDstar> hfSelDstar;
    // Produces<aod::HfSelD0> hfSelD0Candidate;
    HistogramRegistry registry{"DstarSelectorRegistry"};

    void init(InitContext& initContext){
        registry.add("D0InvMass",";#it{InvMass} (GeV/#it{c}^{2});counts",{HistType::kTH1F,{{500,0.0,5.0}}});
        registry.add("DstarInvMass",";#it{InvMass} (GeV/#it{c}^{2});counts",{HistType::kTH1F,{{500,0.0,5.0}}});

    }

    void processQA(aod::HfCandDStarBase const & DstarCands,aod::HfCand2Prong const &){
        LOG(info)<<"selector: processQA called";
        int i=0;
        for(auto & iDstar:DstarCands){
            i++;
            if(i<2)LOG(info)<<"selector: loop start";
            auto iD0 = iDstar.prondD0();
            LOGF(info,"D0p = %f",iD0.p());
            // registry.fill(HIST("D0InvMass"),iD0.m(std::array{0.14,0.49}));
            // registry.fill(HIST("DstarInvMass"),iDstar.dstarmass(std::array{0.14,1.86}));
            if(i>5) continue;
        }
    }
    PROCESS_SWITCH(HfCandidateSelctorDstar,processQA,"Process QA",true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
    return WorkflowSpec{
        adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}
