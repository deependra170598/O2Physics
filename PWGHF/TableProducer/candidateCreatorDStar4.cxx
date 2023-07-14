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
    using D0Candidate=D0Table::iterator;
    using ContestingPionTrack= soa::Join<aod::Tracks,TracksDCA>;
}
using SelectedPionTrack= soa::Filtered<aod::ContestingPionTrack>;
using SelectedD0=soa::Filtered<aod::D0Table>;

#define Prong0pT nsqrt((aod::hf_cand::pxProng0)*(aod::hf_cand::pxProng0) + (aod::hf_cand::pyProng0)*(aod::hf_cand::pyProng0))
#define Prong1pT nsqrt((aod::hf_cand::pxProng1)*(aod::hf_cand::pxProng1) + (aod::hf_cand::pyProng1)*(aod::hf_cand::pyProng1))


struct HfCandidateCreatorDstar{
    // topological cuts
    Configurable<int>NumnberOfpTbinsForD0{"NumnberOfpTbinsForD0",11," no of pt bins of D0"};
    Configurable<std::vector<double>> binsPtD0{"binsPtD0",std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt},"pT bin limits for D0"};

    Configurable<int>NumnberOfpTbinsForDStar{"NumnberOfpTbinsForDStar",11,"no of pt bins of DStar"};
    Configurable<std::vector<double>> binsPt{"binsPt",std::vector<double>{hf_cuts_DStar_to_D0_Pi::vecBinsPt},"pT bin limits for D*"};
    Configurable<LabeledArray<double>> ConfCuts{"ConfCuts",{hf_cuts_DStar_to_D0_Pi::cuts[0],nBinsPt,nCutVars,labelsPt,labelsCutVar},"DStar candidate selection per pT bin"};

    ConfigurableAxis DeltaInvMassAxis{"DeltaInvMassAxis",{500,0.13,0.16},"Delta_M_Star invariant mass axis"};
    ConfigurableAxis DStarInvMassAxis{"DStarInvMassAxis",{500,0.0,5.0},"DStar Invarient mass axis"};
    ConfigurableAxis DStarEtaAxis{"DStarEtaAxis",{500,-0.8,0.8},"DStar Eta Axis"};
    ConfigurableAxis DStarPhiAxis{"DStarPhiAxis",{500,-PI,PI},"DStar Phi Axis"};
    ConfigurableAxis DStarPtAxis{"DStarPtAxis",{500,0.0,5.0},"DStar pT Axis"};

    ConfigurableAxis D0InvMassAxis{"D0InvMassAxis",{500,0.0,5.0},"D0 Invarient mass axis"};
    ConfigurableAxis D0EtaAxis{"D0EtaAxis",{500,-0.8,0.8},"D0 Eta Axis"};
    ConfigurableAxis D0PhiAxis{"D0PhiAxis",{500,0,2*PI},"D0 Phi Axis"};
    ConfigurableAxis D0PtAxis{"D0PtAxis",{500,0.0,5.0},"D0 pT Axis"};

    Configurable<float> DCAXY_Threshold_of_SoftPi{"DCAXYThresholdOfSoftPi",0.2f,"Maximum DCAXY value for soft pi"};
    Configurable<float> DCAZ_Threshold_of_SoftPi{"DCAZThresholdOfSoftPi",0.2f,"Maximum DCAZ value for soft pi"};

    // PDG mass for two prong
    double massPi = RecoDecay::getMassPDG(kPiPlus);
    double masska = RecoDecay::getMassPDG(kKMinus);

    HistogramRegistry DStarRegistry{"DStarRegistry",{},OutputObjHandlingPolicy::AnalysisObject,true,true};
    HistogramRegistry D0Registry{"D0Registry",{},OutputObjHandlingPolicy::QAObject,true,true};


    void init(InitContext &){
        
        AxisSpec DeltaInvMassAxisSpec={DeltaInvMassAxis,"#Delta M_{inv} (GeV/c^2)"};
        AxisSpec DStarInvMassAxisSpec={DStarInvMassAxis,"#it{M}_{inv} (GeV/c^2) "};
        AxisSpec DStarEtaAxisSpec={DStarEtaAxis,"#eta"};
        AxisSpec DStarPhiAxisSpec={DStarPhiAxis,"#phi (rad)"};
        AxisSpec DStarPtAxisSpec={DStarPtAxis,"#it{p}_{T} (GeV/#it{c})"};

        AxisSpec D0InvMassAxisSpec={D0InvMassAxis,"#it{M}_{inv} (GeV/c^2) "};
        AxisSpec D0EtaAxisSpec={D0EtaAxis,"#eta"};
        AxisSpec D0PhiAxisSpec={D0PhiAxis,"#phi (rad)"};
        AxisSpec D0PtAxisSpec={D0PtAxis,"#it{p}_{T} (GeV/#it{c})"};

        // DStar histograms
        DStarRegistry.add("DStarInvMass","#it{M}_{inv} of D*",{kTH1D,{DStarInvMassAxisSpec}},true);
        DStarRegistry.add("DStar_Eta", "#eta distribution of reconstructed D*", {kTH1D, {DStarEtaAxisSpec}},true);
        DStarRegistry.add("DStar_phi", "#phi distribution of reconstructed D*", {kTH1D, {DStarPhiAxisSpec}},true);
        DStarRegistry.add("DStarPt", "#it{p}_{T} distribution of reconstructed D*", {kTH1D, {DStarPtAxisSpec}},true);
        DStarRegistry.add("DeltaMStar","#Delta M_{inv} distribution of reconstructed D* for inclusive #it{p}_{T}",{kTH1D,{DeltaInvMassAxisSpec}},true);
        
        // D0 Histograms
        D0Registry.add("D0InvMass","M_{inv} of two prong candidates with D0Sel or/or/& D0BarSel>=1",{kTH1D,{D0InvMassAxisSpec}},true);
        D0Registry.add("D0Pt","#it{p}_{T} distribution of two prong candidates with D0Sel or/& D0BarSel>=1",{kTH1D,{D0PtAxisSpec}},true);
        D0Registry.add("D0Eta","#eta distribution of two prong candidates with D0Sel or/& D0BarSel>=1",{kTH1D,{D0EtaAxisSpec}},true);
        D0Registry.add("D0Phi","#phi distribution of two prong candidates with D0Sel or/& D0BarSel>=1",{kTH1D,{D0PhiAxisSpec}},true);

        NpTBinsD0=NumnberOfpTbinsForD0;
        NpTBinsDStar=NumnberOfpTbinsForDStar;
        for(int i=0;i<NpTBinsDStar;i++){
            std::string HistNameDStar="DifferentPTBins/InvMassDeltMDStar_"+static_cast<std::string>(ptBins[i]);
            std::string HistTitleDStar="#Delta #it{M}_{inv} of D* for #it{p}_{T} bin #"+std::to_string(i);
            DStarRegistry.add(HistNameDStar.c_str(),HistTitleDStar.c_str(),{kTH1D,{DeltaInvMassAxisSpec}},true);
        }
        for(int i=0;i<NpTBinsD0;i++){
            std::string HistNameD0="DifferentPTBins/InvMassMD0_"+static_cast<std::string>(ptBins[i]);
            std::string HistTitleD0="#it{M}_{inv} of D0 for #it{p}_{T} bin #"+std::to_string(i);
            D0Registry.add(HistNameD0.c_str(),HistTitleD0.c_str(),{kTH1D,{D0InvMassAxisSpec}},true);
        }
    }
    int NpTBinsD0;
    int NpTBinsDStar;
    static constexpr std::string_view ptBins[25] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"};
    

    // Data Selection Filter
    Filter SoftTrackFilter= (aod::track::dcaZ<DCAZ_Threshold_of_SoftPi) && (aod::track::dcaXY<DCAXY_Threshold_of_SoftPi);
    Filter D0Filter= (aod::hf_sel_candidate_d0::isSelD0>=1) || (aod::hf_sel_candidate_d0::isSelD0bar>=1);

    
    // template<const int Index,const int MaxNoOfptBins,typename T2>
    // constexpr void FillHistD0( T2 candidate){
    //     // const int MaxNoOfptBins=25; // Max number of allowed pt bins
    //     if constexpr(Index<MaxNoOfptBins){
    //         if(candidate.pt()>binsPtD0->at(Index) && candidate.pt()<=binsPtD0->at(Index+1))
    //         D0Registry.fill(HIST(ptBins[Index]),candidate.pt());
    //         FillHistD0<Index+1,MaxNoOfptBins,T2>(candidate);
    //     }else{return;}
    // }

    void ProcessD0Only(SelectedD0 const & d0Table){
        for( auto & d0: d0Table)
        {
            D0Registry.fill(HIST("D0Pt"),d0.pt());
            D0Registry.fill(HIST("D0Eta"),d0.eta());
            D0Registry.fill(HIST("D0Phi"),d0.phi());
            // FillHistD0<0,11,SelectedD0::iterator>(d0);

//============================================================================================
                    // Which mass you are finding D0 or D0Bar? Correct it.                    //
                    double massD0=-999.;                                                      //
                    if(d0.isSelD0()>=1 && d0.isSelD0bar()<1){                   //
                        massD0= d0.m(array{massPi,masska});                                   //
                    }else if(d0.isSelD0bar()>=1 && d0.isSelD0()<1){             //
                        massD0= d0.m(array{masska,massPi});                                   //
                    }else if(d0.isSelD0()>=1 && d0.isSelD0bar()>=1){                          //
                    // For time being we are considering it as D0 but need to change it. It is not coorect. What is solution?                                //
                        massD0= d0.m(array{massPi,masska});                                   //
                    }                                                                         //
//=============================================================================================
            D0Registry.fill(HIST("D0InvMass"),massD0);

            int ptBin=analysis::findBin(binsPtD0,d0.pt());
            switch(ptBin) {
                case 0:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_0"),massD0);
                    break;
                case 1:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_1"),massD0);
                    break;
                case 2:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_2"),massD0);
                    break;
                case 3:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_3"),massD0);
                    break;
                case 4:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_4"),massD0);
                    break;
                case 5:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_5"),massD0);
                    break;
                case 6:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_6"),massD0);
                    break;
                case 7:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_7"),massD0);
                    break;
                case 8:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_8"),massD0);
                    break;
                case 9:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_9"),massD0);
                    break;
                case 10:
                    D0Registry.fill(HIST("DifferentPTBins/InvMassMD0_10"),massD0);
                    break;
                }
        }
    }

    PROCESS_SWITCH(HfCandidateCreatorDstar,ProcessD0Only,"generate D0 QA",true);
    
    // table slicer according to collisionId
    Preslice<aod::D0Table> D0Slice=aod::hf_cand::collisionId;
    Preslice<aod::ContestingPionTrack> Trackslice= aod::track::collisionId;

    SliceCache cache;

    void process(aod::Collisions const & cols, SelectedPionTrack const & tracks, SelectedD0 const & d0Table){

        // Defining Partition for each pT bins
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

            // DStar Reconstruction loop
            for(int i=0;i<NpTBinsDStar;i++){
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
                    // For time being we are considering it as D0 but need to change it. It is not coorect. What is solution?                                //
                        massD0= d0.m(array{massPi,masska});                                   //
                    }                                                                         //
//=============================================================================================

                    double massDStar = RecoDecay::m(std::array{pVect, pVecD0}, std::array{massPi, massD0});
                    auto DeltaMDStar = RecoDecay::m(std::array{pVect, pVecD0}, std::array{massPi, massD0}) - massD0;
                    
                    // Filling InvMass Hist
                    DStarRegistry.fill(HIST("DStarInvMass"), massDStar);
                    DStarRegistry.fill(HIST("DeltaMStar"), DeltaMDStar);

                    // Find momentum of DStar
                    auto Dstar_px = t.px() + d0.px();
                    auto Dstar_py = t.py() + d0.py();
                    auto Dstar_pz = t.pz()+ d0.pz();

                    auto DStar_pT = sqrt(Dstar_px * Dstar_px + Dstar_py * Dstar_py);
                    auto DStar_P = sqrt(Dstar_px * Dstar_px + Dstar_py * Dstar_py + Dstar_pz*Dstar_pz);
                    // auto DStar_E = sqrt(pow(DStar_P,2) + pow(massDStar,2));
                    auto DStar_eta = 0.5*log((DStar_P+Dstar_pz)/(DStar_P-Dstar_pz));
                    auto DStar_phi = atan2(Dstar_py,Dstar_px);

                    DStarRegistry.fill(HIST("DStarPt"), DStar_pT);
                    DStarRegistry.fill(HIST("DStar_Eta"), DStar_eta);
                    DStarRegistry.fill(HIST("DStar_phi"), DStar_phi);
                    
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


