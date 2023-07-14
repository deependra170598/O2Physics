#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"

#include "Framework/AnalysisDataModel.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Framework/Expressions.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_cuts_DStar_to_D0_Pi;

namespace o2::aod{
    using D0Table= soa::Join<HfCand2Prong,HfSelD0>;
    using ContestingPionTrack= soa::Join<aod::Tracks,TracksDCA>;
}
using SelectedPionTrack= soa::Filtered<aod::ContestingPionTrack>;
using SelectedD0=soa::Filtered<aod::D0Table>;

#define Prong0pT nsqrt((aod::hf_cand::pxProng0)*(aod::hf_cand::pxProng0) + (aod::hf_cand::pyProng0)*(aod::hf_cand::pyProng0))
#define Prong1pT nsqrt((aod::hf_cand::pxProng1)*(aod::hf_cand::pxProng1) + (aod::hf_cand::pyProng1)*(aod::hf_cand::pyProng1))


struct HfCandidateCreatorDstar{
    // topological cuts
    Configurable<std::vector<double>> binsPt{"binsPt",std::vector<double>{hf_cuts_DStar_to_D0_Pi::vecBinsPt},"pT bin limits"};
    Configurable<LabeledArray<double>> ConfCuts{"ConfCuts",{hf_cuts_DStar_to_D0_Pi::cuts[0],nBinsPt,nCutVars,labelsPt,labelsCutVar},"DStar candidate selection per pT bin"};
    
    // PDG mass for two prong
    double massPi = RecoDecay::getMassPDG(kPiPlus);
    double masska = RecoDecay::getMassPDG(kKMinus);

    HistogramRegistry DStarRegistry{"DStarRegistry",{},OutputObjHandlingPolicy::AnalysisObject,true,true};


    void init(InitContext &){
        DStarRegistry.add("DStarInvMass","DStarInvMass",{kTH1D,{{100,0.0,5.0}}});
        DStarRegistry.add("AllTrackPt","AllTrackPt",{kTH1D,{{100,0.0,5.0}}});
        DStarRegistry.add("D0Pt","D0Pt",{kTH1D,{{100,0.0,5.0}}});
        DStarRegistry.add("DeltaMStar","DeltaMStar",{kTH1D,{{200,0.0,1.0}}});
        for(int i=0;i<nBinsPt;i++){
            std::string HistName="DifferentPTBins/InvMassDeltMDStar_"+static_cast<std::string>(ptBins[i]);
            std::string HistTitle="InvMassDeltMDStar_#"+std::to_string(i);
            DStarRegistry.add(HistName.c_str(),HistTitle.c_str(),{kTH1D,{{200,0.0,1.0}}});
        }
    }
    static constexpr std::string_view ptBins[nBinsPt] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10"};
    
    // Data Selection Filter
    Filter SoftTrackFilter= (aod::track::dcaZ<0.2f) && (aod::track::dcaXY<0.2f);
    Filter D0Filter= (aod::hf_sel_candidate_d0::isSelD0>=1) && (aod::hf_sel_candidate_d0::isSelD0bar>=1);

    void process(aod::Collision const & col, SelectedPionTrack const & tracks, SelectedD0 const & d0Table){

        // For loop over different pT Bins
        // for(int i=0;i<nBinsPt;i++){
            // =========================================================================================================================
            Partition<SelectedPionTrack> SoftTrack_0= (aod::track::pt> static_cast<float> (ConfCuts->get(uint32_t(0),"min_pT_soft_pi"))) 
                                                &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint32_t(0),"max_pT_soft_pi")));
            SoftTrack_0.bindTable(tracks);
            
            Partition<SelectedD0> SelectedD0_0=Prong0pT> static_cast<float>(ConfCuts->get(uint32_t(0),"min_pT_prong0")) 
                                        && Prong1pT> static_cast<float> (ConfCuts->get(uint32_t(0),"min_pT_prong1"));
            SelectedD0_0.bindTable(d0Table);
            // =========================================================================================================================
            Partition<SelectedPionTrack> SoftTrack_1= (aod::track::pt> static_cast<float> (ConfCuts->get(uint32_t(1),"min_pT_soft_pi"))) 
                                                &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint32_t(1),"max_pT_soft_pi")));
            SoftTrack_1.bindTable(tracks);
            
            Partition<SelectedD0> SelectedD0_1=Prong0pT> static_cast<float>(ConfCuts->get(uint32_t(1),"min_pT_prong0")) 
                                        && Prong1pT> static_cast<float> (ConfCuts->get(uint32_t(1),"min_pT_prong1"));
            SelectedD0_1.bindTable(d0Table);
            // =========================================================================================================================
            Partition<SelectedPionTrack> SoftTrack_2= (aod::track::pt> static_cast<float> (ConfCuts->get(uint32_t(2),"min_pT_soft_pi"))) 
                                                &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint32_t(2),"max_pT_soft_pi")));
            SoftTrack_2.bindTable(tracks);
            
            Partition<SelectedD0> SelectedD0_2=Prong0pT> static_cast<float>(ConfCuts->get(uint32_t(2),"min_pT_prong0")) 
                                        && Prong1pT> static_cast<float> (ConfCuts->get(uint32_t(2),"min_pT_prong1"));
            SelectedD0_2.bindTable(d0Table);
            // =========================================================================================================================
            // =========================================================================================================================
            // =========================================================================================================================
            // =========================================================================================================================
            // =========================================================================================================================
            // =========================================================================================================================

            // Testing
            std::vector<Partition<SelectedD0>*>PartitionVec1 {&SelectedD0_0,&SelectedD0_1,&SelectedD0_2};
            std::vector<Partition<SelectedPionTrack>*>PartitionVec2 {&SoftTrack_0,&SoftTrack_1,&SoftTrack_2};
        for(int i=0;i<3;i++){
            for(auto & [t,d0]: combinations(CombinationsFullIndexPolicy(*PartitionVec2[i],*PartitionVec1[i]))){    
                // safe gaurd
                if((t.collisionId()!= d0.collisionId()) && (t.collisionId() != col.globalIndex())) continue;
                if(t.globalIndex()==d0.prong0Id() || t.globalIndex()==d0.prong1Id()) continue;

                // LOGF(info,"2. Tracks did not contribute to D0");
                std::array<float, 3> pVect = {t.px(), t.py(), t.pz()};
                std::array<float, 3> pVecD0=d0.pVector();
                double massD0= d0.m(array{massPi,masska});
                double massDStar = RecoDecay::m(std::array{pVect, pVecD0}, std::array{massPi, massD0});
                auto DeltaMDStar = RecoDecay::m(std::array{pVect, pVecD0}, std::array{massPi, massD0}) - massD0;
                
                // Filling Hist
                DStarRegistry.fill(HIST("AllTrackPt"),t.pt());
                DStarRegistry.fill(HIST("D0Pt"),d0.pt());
                DStarRegistry.fill(HIST("DStarInvMass"),massDStar);
                DStarRegistry.fill(HIST("DeltaMStar"),DeltaMDStar);


                auto Dstar_px=t.px()+d0.px();
                auto Dstar_py=t.py()+d0.py();
                auto DStar_pT=sqrt(Dstar_px*Dstar_px + Dstar_py*Dstar_py);
                
                // Filling pT bin wise
                if(i==0 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_0"),DeltaMDStar);
                }else if(i==1 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_1"),DeltaMDStar);
                }else if(i==2 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_2"),DeltaMDStar);
                }else if(i==3 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_3"),DeltaMDStar);
                }else if(i==4 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_4"),DeltaMDStar);
                }else if(i==5 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_5"),DeltaMDStar);
                }else if(i==6 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_6"),DeltaMDStar);
                }else if(i==7 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_7"),DeltaMDStar);
                }else if(i==8 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_8"),DeltaMDStar);
                }else if(i==9 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_9"),DeltaMDStar);
                }else if(i==10 && DStar_pT>binsPt->at(i) && DStar_pT<binsPt->at(i+1)){
                    DStarRegistry.fill(HIST("DifferentPTBins/InvMassDeltMDStar_10"),DeltaMDStar);
                }

            }
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc)};
}


