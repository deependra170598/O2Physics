#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;


struct DstarTest{

    void process(aod::Hf2Prong const &, aod::HfD0fromDstar const & rowsD0cand, aod::HfCandDStarBase const & rowsDstarCand){

        for(auto & rowDstarCand:rowsDstarCand){
            
            LOGF(info,"DstarTest: collisionId from HfCandDStarBase Table = %d",rowDstarCand.collisionId());
            auto iterator = rowDstarCand.globalIndex();
            auto rowD0Cand = rowsD0cand.iteratorAt(iterator);
            LOGF(info,"DstarTest: collisionId from HfD0fromDstar Table = %d",rowD0Cand.collisionId());
            auto D0prong = rowDstarCand.prongD0();
            LOGF(info,"DstarTest: collisionId from HfProng2 Table = %d",D0prong.collisionId() );
        }

    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
    return WorkflowSpec
    {
        adaptAnalysisTask<DstarTest>(cfgc)
    };
}