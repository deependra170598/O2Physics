#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "Framework/O2DatabasePDGPlugin.h"

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
  Service<o2::framework::O2DatabasePDG> o2ServicePDG;
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
    Configurable <double> ptDstarCandMin{"ptDstarCandMin",0.,"Lower bound of Dstar candidate pT"};
    Configurable <double> ptDstarCandmax{"ptDstarCandmax",50.,"upper bound of Dstar candidate pT"};

    // topological cuts
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

  double massPi, massK, massD0;
  void init (InitContext& initContext){
    massPi = o2ServicePDG->Mass(211);
    massK = o2ServicePDG->Mass(311);
    massD0 = o2ServicePDG->Mass(421);
  }


  /// Conjugate-independent topological cuts on D0
  /// @brief Topological selection on D0 candidate from Dstar
  /// @tparam T table iterator type of the candidate
  /// @param candidate candidate instance(object)
  /// @return true or false depending on selected or not
  template<typename T>
  bool selectionTopolD0(const T & candidate){
    auto candpT = candidate.d0pt();
    auto pTBin = findBin(d0configurable.binsPtForD0,candpT);
    if(pTBin == -1){
      return false;
    }

    // check that the candidate pT is within the analysis range
    if(candpT < d0configurable.ptD0CandMin || candpT >= d0configurable.ptD0CandMax){
      return false;
    }
    // product of daughter impact parameters
    if(candidate.d0impactParameterProduct() > d0configurable.cutsForD0->get(pTBin,"d0d0")){
      return false;
    }
    // // cosine of pointing angle
    if(candidate.d0cpaXY() < d0configurable.cutsForD0->get(pTBin,"cos pointing angle xy")){
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
    double decayLengthCut = std::min((candidate.d0p()*0.0066)+0.01,d0configurable.cutsForD0->get(pTBin,"minimum decay length"));
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


  /// Conjugate-independent topological cuts on Dstar
  /// @brief selection on Dstar candidate
  /// @tparam T1 table iterator type of the D0 candidate
  /// @tparam T2 table iterator type of the Dstar candidate
  /// @param d0Candidate D0 candidate object
  /// @param dstarCandidate Dstar candidate object
  /// @return true or false depending on selected or not
  template <typename T1, typename T2>
  bool selectionOnDstar(const T1 & d0Candidate, const T2 & dstarCandidate ){
    auto dstarpT = dstarCandidate.candDStarpt();
    auto pTBin = findBin(dstarconfigurable.binsPtForDstar,dstarpT);
    if(pTBin = -1){
      return false;
    }
    // check that the candidate pT is within the analysis range
    if(dstarpT < dstarconfigurable.ptDstarCandMin || dstarpT >= dstarconfigurable.ptDstarCandmax){
      return false;
    }
    // selction on DCA of softpi
    if(std::abs(dstarCandidate.impParamSoftPiProng()) > dstarconfigurable.cutsForDstar->get(pTBin,"d0Softpi")){
      return false;
    }
    if(std::abs(dstarCandidate.normalisedImpParamSoftPiProng()) > dstarconfigurable.cutsForDstar->get(pTBin,"d0SoftPiNormalised")){
      return false;
    }
    // selection on pT of soft Pi
    if((dstarCandidate.ptSoftpiProng() < dstarconfigurable.cutsForDstar->get(pTBin,"minpTSoftPi"))||(dstarCandidate.ptSoftpiProng() > dstarconfigurable.cutsForDstar->get(pTBin,"maxpTSoftPi"))){
      return false;
    }
    // selection on pT of D0Prong1,2
    auto d0prong0 = d0Candidate.prong0();
    if(d0prong0.pt() < dstarconfigurable.cutsForDstar->get(pTBin,"minpTD0_Prong0")){
      return false;
    }
    auto d0prong1 = d0Candidate.prong1();
    if(d0prong1.pt() < dstarconfigurable.cutsForDstar->get(pTBin,"minpTD0_Prong1")){
      return false;
    }
    // selection on D0Candidate
    if(!selectionTopolD0(d0Candidate)){
      return false;
    }
    return true;
  }
  

  


  bool selectionTopolConjugate(aod::HfD0fromDstar::iterator const & d0candidate ,aod::HfCandDStarBase::iterator const & dstarCandidate){
    auto dstarpT = dstarCandidate.candDStarpt();
    auto pTBin = findBin(dstarconfigurable.binsPtForDstar,dstarpT);
    if(pTBin = -1){
      return false;
    }
    auto SoftPi = dstarCandidate.prongPi();
    auto ProngD0 = dstarCandidate.prongD0();
    auto PosProng0 = ProngD0.prong0();
    auto NegProng1 = ProngD0.prong1();
    // Selection of Dstar+
    if(SoftPi.sign() > 0.){
      auto Mpipik = std::array{massPi,massPi,massK};
      auto InvDstar = dstarCandidate.dstarInvMass(Mpipik);
      auto Mpik = std::array{massPi,massK};
      auto InvD0 = d0candidate.d0m(Mpik);
      if(std::abs(InvDstar-InvD0)> dstarconfigurable.cutsForDstar->get(pTBin,"DeltaMDstar")){
        return false;
      }
    }else if(SoftPi.sign()<0){
      auto Mpikpi = std::array{massPi,massK,massPi};
      auto InvDstar = dstarCandidate.dstarInvMass(Mpikpi);
      auto Mkpi = std::array{massK,massPi};
      auto InvD0 = d0candidate.d0m(Mkpi);
      if(std::abs(InvDstar-InvD0)>dstarconfigurable.cutsForDstar->get(pTBin,"DeltaMDstar")){
        return false;
      }
    }
    

    return true;
  }






  void process(aod::Hf2Prong const &, aod::HfD0fromDstar const & rowsD0cand, aod::HfCandDStarBase const & rowsDstarCand)
  {
    LOG(info) << "selector: processQA called";
  }

  
  

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}
