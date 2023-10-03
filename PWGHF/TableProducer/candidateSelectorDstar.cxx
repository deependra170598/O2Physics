#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;


// Struct to applying Dstar selection cuts
struct HfCandidateSelctorDstar {

  void process(aod::Hf2Prongs const&, aod::HfCand2Prong const &, aod::HfCandDStarBase const & rowsDStarCands)
  {
    LOG(info) << "selector: processQA called";
    int i = 0;
    for (const auto& iDstar : rowsDStarCands) {
      ++i;
      if (i < 2){LOG(info) << "selector: loop start";}
      
      LOGF(info,"selector : Collision id of Dstar: %d",iDstar.collisionId());
      auto iD0 = iDstar.prongD0();
      auto iD0Cand = iDstar.prondD0Cand();
      
      LOGF(info,"selector : Collision id of D0: %d",iD0.collisionId());
      LOGF(info,"selector : Collision id of D0Cand: %d",iD0Cand.collisionId());
      if(i>100){continue;}
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}
