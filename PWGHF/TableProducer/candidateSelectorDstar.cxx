#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
// using namespace o2::aod::HFCandDStarProng;
// using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod;
// using namespace o2::analysis::hf_cuts_d0_to_pi_k;
// using namespace o2::analysis::hf_cuts_dstar_to_pi_d0;
using namespace o2::analysis;


// Struct to applying Dstar selection cuts
struct HfCandidateSelctorDstar {
  Produces<aod::HfSelD0> hfSelD0Candidate;
  Produces<aod::HfSelDstar> hfSelDstarCandidate;

  struct D0Configurables: ConfigurableGroup{
    Configurable <double> ptD0CandMin{"ptD0CandMin",0.,"Lower bound of D0 candidate pT"};
    Configurable <double> ptD0CandMax{"ptD0CandMax",50.,"Upper bound of D0 candidate pT"};

    // topological cuts
    Configurable<std::vector<double>> binsPtForD0{"binsPtForD0",std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for D0"};
    Configurable<LabeledArray<double>> cutsForD0{"cutsForD0", {hf_cuts_d0_to_pi_k::cuts[0], 
                                                               hf_cuts_d0_to_pi_k::nBinsPt, 
                                                               hf_cuts_d0_to_pi_k::nCutVars, 
                                                               hf_cuts_d0_to_pi_k::labelsPt,
                                                               hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};

  }d0configurable;

  struct DstarConfigurable: ConfigurableGroup{
    Configurable<std::vector<double>> binsPtForDstar{"binsPtForDstar",std::vector<double>{hf_cuts_dstar_to_pi_d0::vecBinsPt},"pT bin limits for Dstar"};
    Configurable<LabeledArray<double>> cutsForDstar{"cutsForDstar",{hf_cuts_dstar_to_pi_d0::cuts[0],
                                                                    hf_cuts_dstar_to_pi_d0::nBinsPt,
                                                                    hf_cuts_dstar_to_pi_d0::nCutVars,
                                                                    hf_cuts_dstar_to_pi_d0::labelsPt,
                                                                    hf_cuts_dstar_to_pi_d0::labelsCutVar},"Dstar candidate selection per pT bin"};

  }dstarconfigurable;

  struct CommonConfigurable: ConfigurableGroup{
    // TPC PID
    Configurable <double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
    Configurable <double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
    Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
    Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};

    // TOF PID
    Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
    Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
    Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
    Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};

    // AND logic for TOF+TPC PID (as in Run2)
    Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Use AND logic for TPC and TOF PID"};

    // CCDB configuration
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  }commonconfigurable;

  // o2::ccdb::CcdbApi ccdbApi;

  // TrackSelectorPi selectorPion;
  // TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidKa>;


  void init (InitContext& initContext){

  }

  
  bool selectionTopolD0(aod::HfCand2Prong::iterator const & candidate){
    auto candpT = candidate.pt();
    auto pTBin = findBin(d0configurable.binsPtForD0,candpT);
    if(pTBin == -1){
      return false;
    }

    // check that the candidate pT is within the analysis range
    if(candpT < d0configurable.ptD0CandMin || candpT >= d0configurable.ptD0CandMax){
      return false;
    }
    // product of daughter impact parameters
    if(candidate.impactParameterProduct() > d0configurable.cutsForD0->get(pTBin,"d0d0")){
      return false;
    }
    // // cosine of pointing angle
    if(candidate.cpaXY() < d0configurable.cutsForD0->get(pTBin,"cos pointing angle xy")){
      return false;
    }
    // normalised decay length in XY plane
    if(candidate.decayLengthXYNormalised() < d0configurable.cutsForD0->get(pTBin,"normalized decay length XY")){
      return false;
    }
    // candidate DCA
    // if(candidate.chi2PCA() > d0configurable.cutsForD0->get(pTBin,"DCA")){
    //   return false;
    // }

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)

    if(std::abs(candidate.impactParameterNormalised0())< 0.5 || std::abs(candidate.impactParameterNormalised1())<0.5){
      return false;
    }
    double decayLengthCut = std::min((candidate.p()*0.0066)+0.01,d0configurable.cutsForD0->get(pTBin,"minimum decay length"));
    if(candidate.decayLength()*candidate.decayLength() < decayLengthCut * decayLengthCut){
      return false;
    }
    if(candidate.decayLength()> d0configurable.cutsForD0->get(pTBin,"decay length")){
      return false;
    }
    if(candidate.decayLengthXY() > d0configurable.cutsForD0->get(pTBin,"decay length XY")){
      return false;
    }

    //.............Why is this if condition commented?
    if(candidate.decayLengthNormalised()* candidate.decayLengthNormalised() < 1.0){
      // return false; // add back when getter fixed
    }
    return true;
  }


  
  


  


  void process(aod::Hf2Prongs const&, aod::HfCand2Prong const &, aod::HfCandDStarBase const & rowsDStarCands)
  {
    LOG(info) << "selector: processQA called";
  }

  
  

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}
