#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

struct TestingPartitionInLoop{
    
    Preslice<aod::Tracks> perCol = aod::track::collisionId;
    HistogramRegistry registry{"registry",{},OutputObjHandlingPolicy::AnalysisObject,true,true};
    void init(InitContext &){
        //registry.add("one","one",{kTH1F,{{20,-10,10}}});
        //registry.add("two","two",{kTH1F,{{20,-10,10}}});

    }
    void process(aod::Collisions const & cols, aod::Tracks const & tracks){
        
            Partition<aod::Collisions> colsleft =aod::collision::posZ< static_cast<float>(0);
            colsleft.bindTable(cols);
            
            

            for(auto & col: colsleft){
                auto Slicetracks=tracks.sliceBy(perCol,col.globalIndex());
                LOGF(info,"Left Track sixe %d",Slicetracks.size());

            }

            
        
    }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec{
        adaptAnalysisTask<TestingPartitionInLoop>(cfgc)
    };
}


