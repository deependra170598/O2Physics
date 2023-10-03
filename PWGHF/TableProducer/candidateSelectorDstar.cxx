#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

using DStarCands = soa::Join<aod::HfCand2Prong, aod::HfCandDStarBase>;

// Struct to applying Dstar selection cuts
struct HfCandidateSelctorDstar {
  // Produces<aod::HfSelDstar> hfSelDstar;
  // Produces<aod::HfSelD0> hfSelD0Candidate;
//   HistogramRegistry registry{"DstarSelectorRegistry"};

//   void init(InitContext& initContext)
//   {
//     registry.add("D0InvMass", ";#it{InvMass} (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{500, 0.0, 5.0}}});
//     registry.add("DstarInvMass", ";#it{InvMass} (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{500, 0.0, 5.0}}});
//   }

  void process(aod::Collision const & rowCollision, aod::Hf2Prongs const&, DStarCands const & rowsDStarCands)
  {
    LOG(info) << "selector: processQA called";
    
    LOGF(info,"selector: collision global Index: %d",rowCollision.globalIndex());
    int i = 0;
    for (const auto& iDstar : rowsDStarCands) {
      i++;
      if (i < 2){LOG(info) << "selector: loop start";}
      
      LOGF(info,"selector : Collision id of Dstar: %d",iDstar.collisionId());
      auto iD0 = iDstar.prongD0();
      LOGF(info,"selector : Collision id of D0: %d",iD0.collisionId());
    //   LOGF(info, "D0p = %f", iD0.p());
      // registry.fill(HIST("D0InvMass"),iD0.m(std::array{0.14,0.49}));
      // registry.fill(HIST("DstarInvMass"),iDstar.dstarmass(std::array{0.14,1.86}));
      if (i > 5)
        continue;
    }
  }

  // void processD0Only(aod::Collision const &, aod::HfCand2Prong const & D0Cands){
  //   LOG(info)<<"processD0Only callled";
  //   for(const auto & iD0: D0Cands){
  //     LOGF(info,"selector: collision Id of D0 = %d", iD0.collisionId());
  //   }
  // }
  // PROCESS_SWITCH(HfCandidateSelctorDstar, processD0Only, "Process D0 only", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}
