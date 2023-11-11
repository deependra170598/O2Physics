// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file candidateSelectiorDstar.cxx
/// \brief Selection on D* decay candidates
///
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"

#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::analysis;

// Struct to applying Dstar selection cuts
struct HfCandidateSelctorDstar {
  Produces<aod::HfSelDstarToD0Pi> hfSelDstarCandidate;

  // Configurable specific to D0
  Configurable<double> ptD0CandMin{"ptD0CandMin", 0., "Lower bound of D0 candidate pT"};
  Configurable<double> ptD0CandMax{"ptD0CandMax", 50., "Upper bound of D0 candidate pT"};
  Configurable<std::vector<double>> binsPtForD0{"binsPtForD0", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for D0"};
  Configurable<LabeledArray<double>> cutsForD0{"cutsForD0", {hf_cuts_d0_to_pi_k::cuts[0], hf_cuts_d0_to_pi_k::nBinsPt, hf_cuts_d0_to_pi_k::nCutVars, hf_cuts_d0_to_pi_k::labelsPt, hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};

  // Configurable specific to Dstar
  Configurable<double> ptDstarCandMin{"ptDstarCandMin", 0., "Lower bound of Dstar candidate pT"};
  Configurable<double> ptDstarCandmax{"ptDstarCandmax", 50., "upper bound of Dstar candidate pT"};
  Configurable<std::vector<double>> binsPtForDstar{"binsPtForDstar", std::vector<double>{hf_cuts_dstar_to_pi_d0::vecBinsPt}, "pT bin limits for Dstar"};
  Configurable<LabeledArray<double>> cutsForDstar{"cutsForDstar", {hf_cuts_dstar_to_pi_d0::cuts[0], hf_cuts_dstar_to_pi_d0::nBinsPt, hf_cuts_dstar_to_pi_d0::nCutVars, hf_cuts_dstar_to_pi_d0::labelsPt, hf_cuts_dstar_to_pi_d0::labelsCutVar}, "Dstar candidate selection per pT bin"};

  // common Configurable
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};

  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};

  // AND logic for TOF+TPC PID
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Use AND logic for TPC and TOF PID"};

  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  HfHelper::D0FromDstar hfHelper;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidKa>;
  // using TracksSel = soa::Join<aod::Tracks, aod::TracksPidPi, aod::TracksPidKa>;

  double massPi, massK, massD0;
  void init(InitContext& initContext)
  {
    massPi = o2::analysis::pdg::MassPiPlus;
    massK = o2::analysis::pdg::MassKPlus;
    massD0 = o2::analysis::pdg::MassD0;

    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorKaon = selectorPion;
  }

  /// Conjugate-independent topological cuts on D0
  /// @brief Topological selection on D0 candidate from Dstar
  /// @tparam T table iterator type of the candidate
  /// @param candidate candidate instance(object)
  /// @return true or false depending on selected or not
  template <typename T>
  bool selectionTopolD0(const T& candidate)
  {
    auto candpT = candidate.d0pt();
    auto pTBin = findBin(binsPtForD0, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptD0CandMin || candpT >= ptD0CandMax) {
      return false;
    }
    // product of daughter impact parameters
    if (candidate.d0impactParameterProduct() > cutsForD0->get(pTBin, "d0d0")) {
      return false;
    }
    // cosine of pointing angle
    if(candidate.d0cpa() < cutsForD0->get(pTBin,"cos pointing angle")){
      return false;
    }
    // cosine of pointing angle XY
    if (candidate.d0cpaXY() < cutsForD0->get(pTBin, "cos pointing angle xy")) {
      return false;
    }
    // normalised decay length in XY plane
    if (candidate.decayLengthXYNormalised() < cutsForD0->get(pTBin, "normalized decay length XY")) {
      return false;
    }

    // Note: follwoing two cuts are not defined in namespace: hf_cuts_d0_to_pi_k of  SelectionCuts.h, while are defined in namespace: hf_cuts_dstar_to_pi_d0
    // Chi2PCA of secondary vertex reconstruction
    if (candidate.chi2PCA() > cutsForDstar->get(pTBin, "Chi2PCA")) {
      return false;
    }
    if(candidate.d0impactParameterXY() > cutsForD0->get(pTBin,"DCA")){
      return false;
    }
    // d0NormalisedProng0,1
    if (std::abs(candidate.impactParameterNormalised0()) < cutsForDstar->get(pTBin, "d0NormalisedProng0") || std::abs(candidate.impactParameterNormalised1()) < cutsForDstar->get(pTBin, "d0NormalisedProng1")) {
      return false;
    }

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)

    double decayLengthCut = std::min((candidate.d0p() * 0.0066) + 0.01, cutsForD0->get(pTBin, "minimum decay length"));
    if (candidate.decayLength() * candidate.decayLength() < decayLengthCut * decayLengthCut) {
      return false;
    }
    if (candidate.decayLength() > cutsForD0->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cutsForD0->get(pTBin, "decay length XY")) {
      return false;
    }

    //.............Why is this if condition commented?
    if (candidate.decayLengthNormalised() * candidate.decayLengthNormalised() < 1.0) {
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
  bool selectionOnDstar(const T1& d0Candidate, const T2& dstarCandidate)
  {
    auto dstarpT = dstarCandidate.pt();
    auto pTBin = findBin(binsPtForDstar, dstarpT);
    if (pTBin == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (dstarpT < ptDstarCandMin || dstarpT >= ptDstarCandmax) {
      return false;
    }
    // selction on DCA of softpi
    if (std::abs(dstarCandidate.impParamSoftPiProng()) > cutsForDstar->get(pTBin, "d0Softpi")) {
      return false;
    }
    if (std::abs(dstarCandidate.normalisedImpParamSoftPiProng()) > cutsForDstar->get(pTBin, "d0SoftPiNormalised")) {
      return false;
    }
    // selection on pT of soft Pi
    if ((dstarCandidate.ptSoftpiProng() < cutsForDstar->get(pTBin, "minpTSoftPi")) || (dstarCandidate.ptSoftpiProng() > cutsForDstar->get(pTBin, "maxpTSoftPi"))) {
      return false;
    }
    //--------------------------------------------------------------------
    // selection on pT of D0Prong1,2
    // auto d0prong0 = d0Candidate.template prong0_as<TracksSel>();
    // auto d0prong1 = d0Candidate.template prong1_as<TracksSel>();
    // if (d0prong0.pt() < cutsForDstar->get(pTBin, "minpTD0_Prong0")) {
    //   return false;
    // }
    // if (d0prong1.pt() < cutsForDstar->get(pTBin, "minpTD0_Prong1")) {
    //   return false;
    // }
    //---------------------------------------------------------------------
    // selection on D0Candidate
    if (!selectionTopolD0(d0Candidate)) {
      return false;
    }
    return true;
  }

  /// @brief Conjugate-dependent topological cuts
  /// @tparam T1 Table iterator type of D0 Candidate
  /// @tparam T2 Table iterator type of Dstar Candidate
  /// @param d0candidate D0 candidate object
  /// @param dstarCandidate Dstar candidate object
  /// @return true/false depending on success of selection
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& d0candidate, const T2& dstarCandidate)
  {
    auto dstarpT = dstarCandidate.pt();
    auto pTBin = findBin(binsPtForDstar, dstarpT);
    if (pTBin == -1) {
      return false;
    }
    auto SoftPi = dstarCandidate.template prongPi_as<TracksSel>();

    if (SoftPi.sign() > 0.) { // Selection of D*+
      // auto Mpipik = std::array{massPi,massPi,massK};
      // auto InvDstar = dstarCandidate.dstarInvMass(Mpipik);
      auto InvDstar = hfHelper.invMassDstar(dstarCandidate);

      // auto Mpik = std::array{massPi,massK};
      // auto InvD0 = d0candidate.d0m(Mpik);
      auto InvD0 = hfHelper.invMassD0ToPiK(d0candidate);
      if(std::abs(InvD0 - o2::analysis::pdg::MassD0) > cutsForD0->get(pTBin,"m")){
        return false;
      }
      // cut on daughter pT
      auto d0prong0 = d0candidate.template prong0_as<TracksSel>();
      auto d0prong1 = d0candidate.template prong1_as<TracksSel>();
      if(d0prong0.pt() < cutsForD0->get(pTBin, "pT Pi") || d0prong1.pt() < cutsForD0->get(pTBin, "pT K")){
        return false;
      }
      // LOGF(info,"+ %d,%f,%f",SoftPi.sign(),InvDstar,InvD0);
      if (std::abs(InvDstar - InvD0) > cutsForDstar->get(pTBin, "DeltaMDstar")) {
        return false;
      }
      // cut on D0 daughter DCA - need to add secondary vertex constraint here
      if (std::abs(d0candidate.impactParameter0()) > cutsForD0->get(pTBin, "d0pi") || std::abs(d0candidate.impactParameter1()) > cutsForD0->get(pTBin, "d0K")) {
        return false;
      }
      // cut on cos(theta*)
      if (std::abs(hfHelper.cosThetaStarD0(d0candidate)) > cutsForD0->get(pTBin, "cos theta*")) {
        return false;
      }

    } else if (SoftPi.sign() < 0) { // Selection of D*-
      // auto Mpikpi = std::array{massPi,massK,massPi};
      // auto InvDstar = dstarCandidate.dstarInvMass(Mpikpi);
      auto InvDstar = hfHelper.invMassAntiDstar(dstarCandidate);

      // auto Mkpi = std::array{massK,massPi};
      // auto InvD0 = d0candidate.d0m(Mkpi);
      auto InvD0 = hfHelper.invMassD0barToKPi(d0candidate);
      if(std::abs(InvD0 - o2::analysis::pdg::MassD0) > cutsForD0->get(pTBin,"m")){
        return false;
      }
      // LOGF(info,"- %d,%f,%f",SoftPi.sign(),InvDstar,InvD0);
      if (std::abs(InvDstar - InvD0) > cutsForDstar->get(pTBin, "DeltaMDstar")) {
        return false;
      }
      // cut on D0 daughter DCA - need to add secondary vertex constraint here
      if (std::abs(d0candidate.impactParameter0()) > cutsForD0->get(pTBin, "d0K") || std::abs(d0candidate.impactParameter1()) > cutsForD0->get(pTBin, "d0pi")) {
        return false;
      }
      // cut on cos(theta*)
      if (std::abs(hfHelper.cosThetaStarD0bar(d0candidate)) > cutsForD0->get(pTBin, "cos theta*")) {
        return false;
      }
    }
    return true;
  }

  void process(TracksSel const&, aod::Hf2Prongs const&, aod::HfD0fromDstar const& rowsD0cand, aod::HfDstarCand const& rowsDstarCand)
  {
    // LOG(info) << "selector: processQA called";
    for (auto& DstarCand : rowsDstarCand) {
      auto dstarIndex = DstarCand.globalIndex();
      auto D0Cand = rowsD0cand.iteratorAt(dstarIndex);

      // final selection flag: 0 - rejected, 1 - accepted
      int statusDstar = 0, statusD0Flag = 0, statusTopol = 0, statusCand = 0, statusPID = 0;

      if (!(D0Cand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusD0Flag = 1;
      if (!selectionOnDstar(D0Cand, DstarCand)) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusTopol = 1;
      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for Dstar
      bool topoDstar = selectionTopolConjugate(D0Cand, DstarCand);
      if (!topoDstar) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusCand = 1;

      // track-level PID selection
      int pidTrackPosKaon = -1;
      int pidTrackPosPion = -1;
      int pidTrackNegKaon = -1;
      int pidTrackNegPion = -1;

      if (usePidTpcAndTof) {
        pidTrackPosKaon = selectorKaon.statusTpcAndTof(D0Cand.prong0_as<TracksSel>());
        pidTrackPosPion = selectorPion.statusTpcAndTof(D0Cand.prong0_as<TracksSel>());
        pidTrackNegKaon = selectorKaon.statusTpcAndTof(D0Cand.prong1_as<TracksSel>());
        pidTrackNegPion = selectorPion.statusTpcAndTof(D0Cand.prong1_as<TracksSel>());
      } else {
        pidTrackPosKaon = selectorKaon.statusTpcOrTof(D0Cand.prong0_as<TracksSel>());
        pidTrackPosPion = selectorPion.statusTpcOrTof(D0Cand.prong0_as<TracksSel>());
        pidTrackNegKaon = selectorKaon.statusTpcOrTof(D0Cand.prong1_as<TracksSel>());
        pidTrackNegPion = selectorPion.statusTpcOrTof(D0Cand.prong1_as<TracksSel>());
      }

      int pidDstar = -1;
      if (DstarCand.prongPi_as<TracksSel>().sign() > 0.) {
        if (pidTrackPosPion == TrackSelectorPID::Accepted && pidTrackNegKaon == TrackSelectorPID::Accepted) {
          pidDstar = 1; // accept D*+
        } else if (pidTrackPosPion == TrackSelectorPID::Rejected && pidTrackNegKaon == TrackSelectorPID::Rejected) {
          pidDstar = 0; // reject D*+
        }
      } else if (DstarCand.prongPi_as<TracksSel>().sign() < 0.) {
        if (pidTrackNegPion == TrackSelectorPID::Accepted && pidTrackPosKaon == TrackSelectorPID::Accepted) {
          pidDstar = 1; // Accept D*-
        } else if (pidTrackNegPion == TrackSelectorPID::Rejected && pidTrackPosKaon == TrackSelectorPID::Rejected) {
          pidDstar = 0; // reject D*-
        }
      }

      if (pidDstar == 0) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        continue;
      }

      if ((pidDstar == -1 || pidDstar == 1) && topoDstar) {
        statusDstar = 1; // identified as dstar
      }

      statusPID = 1;
      hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelctorDstar>(cfgc)};
}