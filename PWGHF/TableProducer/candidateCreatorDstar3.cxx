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
    ConfigurableAxis DeltaInvMassAxis{"DeltaInvMassAxis",{100,0.13,0.16},"Delta_M_Star invariant mass axis"};
    Configurable<float> DCAXY_Threshold_of_SoftPi{"DCAXYThresholdOfSoftPi",0.2f,"Maximum DCAXY value for soft pi"};
    Configurable<float> DCAZ_Threshold_of_SoftPi{"DCAZThresholdOfSoftPi",0.2f,"Maximum DCAZ value for soft pi"};

    // PDG mass for two prong
    double massPi = RecoDecay::getMassPDG(kPiPlus);
    double masska = RecoDecay::getMassPDG(kKMinus);

    HistogramRegistry DStarRegistry{"DStarRegistry",{},OutputObjHandlingPolicy::AnalysisObject,true,true};


    void init(InitContext &){
        AxisSpec DeltaInvMassAxisSpec={DeltaInvMassAxis,"#Delta M (GeV/c^2)"};
        DStarRegistry.add("DStarInvMass","DStarInvMass",{kTH1D,{{500,0.0,5.0}}});
        DStarRegistry.add("AllTrackPt","AllTrackPt",{kTH1D,{{500,0.0,5.0}}});
        DStarRegistry.add("D0Pt","D0Pt",{kTH1D,{{100,0.0,5.0}}});
        DStarRegistry.add("DeltaMStar","DeltaMStar",{kTH1D,{DeltaInvMassAxis}});
        for(int i=0;i<nBinsPt;i++){
            std::string HistName="DifferentPTBins/InvMassDeltMDStar_"+static_cast<std::string>(ptBins[i]);
            std::string HistTitle="InvMassDeltMDStar_#"+std::to_string(i);
            DStarRegistry.add(HistName.c_str(),HistTitle.c_str(),{kTH1D,{DeltaInvMassAxis}});
        }
    }
    static constexpr std::string_view ptBins[nBinsPt] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10"};
    
    // Data Selection Filter
    Filter SoftTrackFilter= (aod::track::dcaZ<DCAZ_Threshold_of_SoftPi) && (aod::track::dcaXY<DCAXY_Threshold_of_SoftPi);
    Filter D0Filter= (aod::hf_sel_candidate_d0::isSelD0>=1) || (aod::hf_sel_candidate_d0::isSelD0bar>=1);

    // table slicer according to collisionId
    Preslice<aod::D0Table> D0Slice=aod::hf_cand::collisionId;
    Preslice<aod::ContestingPionTrack> Trackslice= aod::track::collisionId;

    SliceCache cache;

    void process(aod::Collisions const & cols, SelectedPionTrack const & tracks, SelectedD0 const & d0Table){
        // Defining Filters for each pT bins
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_0= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(0),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(0),"max_pT_soft_pi")));
        SoftTrack_0.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_0=Prong0pT> static_cast<float>(ConfCuts->get(uint(0),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(0),"min_pT_prong1"));
        SelectedD0_0.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_1= (aod::track::pt> static_cast<float> (ConfCuts->get(uint(1),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(1),"max_pT_soft_pi")));
        SoftTrack_1.bindTable(tracks);
        
        Partition<SelectedD0> SelectedD0_1=Prong0pT> static_cast<float>(ConfCuts->get(uint(1),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(1),"min_pT_prong1"));
        SelectedD0_1.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_2= (aod::track::pt> static_cast<float> (ConfCuts->get(uint(2),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(2),"max_pT_soft_pi")));
        SoftTrack_2.bindTable(tracks);
        
        Partition<SelectedD0> SelectedD0_2=Prong0pT> static_cast<float>(ConfCuts->get(uint(2),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(2),"min_pT_prong1"));
        SelectedD0_2.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_3= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(3),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(3),"max_pT_soft_pi")));
        SoftTrack_3.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_3=Prong0pT> static_cast<float>(ConfCuts->get(uint(3),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(3),"min_pT_prong1"));
        SelectedD0_3.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_4= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(4),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(4),"max_pT_soft_pi")));
        SoftTrack_4.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_4=Prong0pT> static_cast<float>(ConfCuts->get(uint(4),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(4),"min_pT_prong1"));
        SelectedD0_4.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_5= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(5),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(5),"max_pT_soft_pi")));
        SoftTrack_5.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_5=Prong0pT> static_cast<float>(ConfCuts->get(uint(5),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(5),"min_pT_prong1"));
        SelectedD0_5.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_6= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(6),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(6),"max_pT_soft_pi")));
        SoftTrack_6.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_6=Prong0pT> static_cast<float>(ConfCuts->get(uint(6),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(6),"min_pT_prong1"));
        SelectedD0_6.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_7= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(7),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(7),"max_pT_soft_pi")));
        SoftTrack_7.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_7=Prong0pT> static_cast<float>(ConfCuts->get(uint(7),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(7),"min_pT_prong1"));
        SelectedD0_7.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_8= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(8),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(8),"max_pT_soft_pi")));
        SoftTrack_8.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_8=Prong0pT> static_cast<float>(ConfCuts->get(uint(8),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(8),"min_pT_prong1"));
        SelectedD0_8.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_9= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(9),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(9),"max_pT_soft_pi")));
        SoftTrack_9.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_9=Prong0pT> static_cast<float>(ConfCuts->get(uint(9),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(9),"min_pT_prong1"));
        SelectedD0_9.bindTable(d0Table);
        // =========================================================================================================================
        Partition<SelectedPionTrack> SoftTrack_10= (aod::track::pt > static_cast<float> (ConfCuts->get(uint(10),"min_pT_soft_pi"))) 
                                            &&(aod::track::pt<= static_cast<float>( ConfCuts->get(uint(10),"max_pT_soft_pi")));
        SoftTrack_10.bindTable(tracks);

        Partition<SelectedD0> SelectedD0_10=Prong0pT> static_cast<float>(ConfCuts->get(uint(10),"min_pT_prong0")) 
                                    && Prong1pT> static_cast<float> (ConfCuts->get(uint(10),"min_pT_prong1"));
        SelectedD0_10.bindTable(d0Table);
        // =========================================================================================================================

        // Testing
        
        // LOGF(info,"collision loop starts");
        for(auto & col: cols){
            auto D0_0=SelectedD0_0->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_1=SelectedD0_1->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_2=SelectedD0_2->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_3=SelectedD0_3->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_4=SelectedD0_4->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_5=SelectedD0_5->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_6=SelectedD0_6->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_7=SelectedD0_7->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_8=SelectedD0_8->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_9=SelectedD0_9->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);
            auto D0_10=SelectedD0_10->sliceByCached(aod::hf_cand::collisionId ,col.globalIndex(),cache);

            auto trks_0=SoftTrack_0->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_1=SoftTrack_1->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_2=SoftTrack_2->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_3=SoftTrack_3->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_4=SoftTrack_4->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_5=SoftTrack_5->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_6=SoftTrack_6->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_7=SoftTrack_7->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_8=SoftTrack_8->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_9=SoftTrack_9->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);
            auto trks_10=SoftTrack_10->sliceByCached(aod::track::collisionId ,col.globalIndex(),cache);

            // LOGF(info,"slicing done");
            std::vector<o2::soa::Filtered<SelectedD0>> Vector_D0_Partition{D0_0,D0_1,D0_2,D0_3,D0_4,D0_5,D0_6,D0_7,D0_8,D0_9,D0_10,};
            std::vector<o2::soa::Filtered<SelectedPionTrack>> Vector_SoftPi_Partition{trks_0,trks_1,trks_2,trks_3,trks_4,trks_5,trks_6,trks_7,trks_8,trks_9,trks_10};

            for(int i=0;i<nBinsPt;i++){
                for(auto & [t,d0]: combinations(CombinationsFullIndexPolicy(Vector_SoftPi_Partition[i],Vector_D0_Partition[i]))){    
                    // safe gaurd
                    if((t.collisionId()!= d0.collisionId()) && (t.collisionId() != col.globalIndex())) continue;
                    if(t.globalIndex()==d0.prong0Id() || t.globalIndex()==d0.prong1Id()) continue;

                    // LOGF(info,"2. Tracks did not contribute to D0");
                    std::array<float, 3> pVect = {t.px(), t.py(), t.pz()};
                    std::array<float, 3> pVecD0=d0.pVector();
//============================================================================================
                    // Which mass you are finding D0 or D0Bar? Correct it.                    //
                    double massD0=-999.;                                                      //
                    if(d0.isSelD0()>=1 && d0.isSelD0bar()<1 && t.sign()>0){                   //
                        massD0= d0.m(array{massPi,masska});                                   //
                    }else if(d0.isSelD0bar()>=1 && d0.isSelD0()<1 && t.sign()<0){             //
                        massD0= d0.m(array{masska,massPi});                                   //
                    }else if(d0.isSelD0()>=1 && d0.isSelD0bar()>=1){                          //
                        //need to change it.  It is not coorect                               //
                        massD0= d0.m(array{massPi,masska});                                   //
                    }                                                                         //
//=============================================================================================

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

                }// Dstar reconstruction loop
            }// pT bin loop
        } // collision loop
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc)};
}


