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

/// \file CorrelationTables.h
/// \brief Correlation table definitions.
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_
#define PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D-Dbar correlation pairs
namespace hf_correlation_d_dbar
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtD, ptD, float);
DECLARE_SOA_COLUMN(PtDbar, ptDbar, float);
DECLARE_SOA_COLUMN(MD, mD, float);
DECLARE_SOA_COLUMN(MDbar, mDbar, float);
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int);
} // namespace hf_correlation_d_dbar

DECLARE_SOA_TABLE(DDbarPair, "AOD", "DDBARPAIR",
                  aod::hf_correlation_d_dbar::DeltaPhi,
                  aod::hf_correlation_d_dbar::DeltaEta,
                  aod::hf_correlation_d_dbar::PtD,
                  aod::hf_correlation_d_dbar::PtDbar);

DECLARE_SOA_TABLE(DDbarRecoInfo, "AOD", "DDBARRECOINFO",
                  aod::hf_correlation_d_dbar::MD,
                  aod::hf_correlation_d_dbar::MDbar,
                  aod::hf_correlation_d_dbar::SignalStatus);

// definition of columns and tables for D0-Hadron correlation pairs
namespace hf_correlation_d0_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);       //! DeltaPhi between D0 and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);       //! DeltaEta between D0 and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);                 //! Transverse momentum of D0
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);       //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);                   //! Invariant mass of D0
DECLARE_SOA_COLUMN(MDbar, mDbar, float);             //! Invariant mass of D0bar
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int); //! Tag for D0,D0bar
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);           //! Pool Bin for the MixedEvent

enum ParticleTypeData {
  D0Only = 1,        // Identified as D0
  D0barOnly,         // Identified as D0bar
  D0D0barBoth,       // Identified as both D0 and D0bar
  D0OnlySoftPi = 11, // Identified as D0 with soft pion
  D0barOnlySoftPi,   // Identified as D0bar with soft pion
  D0D0barBothSoftPi  // Identified as both D0 and D0bar with soft pion
};

enum ParticleTypeMcRec {
  D0Sig = 0, // D0 signal
  D0Ref,     // D0 reflection
  D0Bg,      // D0 background
  D0barSig,  // D0bar signal
  D0barRef,  // D0bar reflection
  D0barBg,   // D0bar background
  SoftPi     // pairs including soft pion
};
} // namespace hf_correlation_d0_hadron

DECLARE_SOA_TABLE(DHadronPair, "AOD", "DHADRONPAIR", //! D0-Hadrons pairs Informations
                  aod::hf_correlation_d0_hadron::DeltaPhi,
                  aod::hf_correlation_d0_hadron::DeltaEta,
                  aod::hf_correlation_d0_hadron::PtD,
                  aod::hf_correlation_d0_hadron::PtHadron,
                  aod::hf_correlation_d0_hadron::PoolBin);

DECLARE_SOA_TABLE(DHadronRecoInfo, "AOD", "DHADRONRECOINFO", //! D0-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_d0_hadron::MD,
                  aod::hf_correlation_d0_hadron::MDbar,
                  aod::hf_correlation_d0_hadron::SignalStatus);

// Note: definition of columns and tables for Lc-Hadron correlation pairs
namespace hf_correlation_lc_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);       //! DeltaPhi between Lc and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);       //! DeltaEta between Lc and Hadrons
DECLARE_SOA_COLUMN(PtLc, ptLc, float);               //! Transverse momentum of Lc
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);       //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MLc, mLc, float);                 //! Invariant mass of Lc
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int); //! Tag for LcToPKPi/LcToPiKP
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);           //! Pool Bin for the MixedEvent
} // namespace hf_correlation_lc_hadron

DECLARE_SOA_TABLE(LcHadronPair, "AOD", "LCHPAIR", //! Lc-Hadrons pairs Informations
                  aod::hf_correlation_lc_hadron::DeltaPhi,
                  aod::hf_correlation_lc_hadron::DeltaEta,
                  aod::hf_correlation_lc_hadron::PtLc,
                  aod::hf_correlation_lc_hadron::PtHadron,
                  aod::hf_correlation_lc_hadron::PoolBin);

DECLARE_SOA_TABLE(LcHadronRecoInfo, "AOD", "LCHRECOINFO", //! Lc-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_lc_hadron::MLc,
                  aod::hf_correlation_lc_hadron::SignalStatus);

// definition of columns and tables for Ds-Hadron correlation pairs
namespace hf_correlation_ds_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float); //! DeltaPhi between Ds and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float); //! DeltaEta between Ds and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);           //! Transverse momentum of Ds
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float); //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);             //! Invariant mass of Ds
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);     //! Pool Bin for the MixedEvent
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);  //! Used in MC-Rec, Ds Signal
DECLARE_SOA_COLUMN(IsPrompt, isPrompt, bool);  //! Used in MC-Rec, Ds Prompt or Non-Prompt
} // namespace hf_correlation_ds_hadron

DECLARE_SOA_TABLE(DsHadronPair, "AOD", "DSHPAIR", //! Ds-Hadrons pairs Informations
                  aod::hf_correlation_ds_hadron::DeltaPhi,
                  aod::hf_correlation_ds_hadron::DeltaEta,
                  aod::hf_correlation_ds_hadron::PtD,
                  aod::hf_correlation_ds_hadron::PtHadron,
                  aod::hf_correlation_ds_hadron::PoolBin);

DECLARE_SOA_TABLE(DsHadronRecoInfo, "AOD", "DSHRECOINFO", //! Ds-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_ds_hadron::MD,
                  aod::hf_correlation_ds_hadron::IsSignal);

DECLARE_SOA_TABLE(DsHadronGenInfo, "AOD", "DSHGENINFO", //! Ds-Hadrons pairs Generated Informations
                  aod::hf_correlation_ds_hadron::IsPrompt);

// definition of columns and tables for Dplus properties
namespace hf_dplus_meson
{
DECLARE_SOA_COLUMN(Phi, phi, float);               //! Phi of D+
DECLARE_SOA_COLUMN(Eta, eta, float);               //! Eta of D+
DECLARE_SOA_COLUMN(PtD, ptD, float);               //! Transverse momentum of D+
DECLARE_SOA_COLUMN(MD, mD, float);                 //! Invariant mass of D+
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);         //! Pool Bin of event defined using zvtx and multiplicity
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);     //! Global index for the collision
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t); //! Timestamp for the collision
} // namespace hf_dplus_meson

DECLARE_SOA_TABLE(Dplus, "AOD", "DPLUS", //! D+-meson properties
                  aod::hf_dplus_meson::Phi,
                  aod::hf_dplus_meson::Eta,
                  aod::hf_dplus_meson::PtD,
                  aod::hf_dplus_meson::MD,
                  aod::hf_dplus_meson::PoolBin,
                  aod::hf_dplus_meson::GIndexCol,
                  aod::hf_dplus_meson::TimeStamp);

// definition of columns and tables for associated hadron properties
namespace hf_assoc_tracks
{
DECLARE_SOA_COLUMN(Phi, phi, float);               //! Phi of hadron
DECLARE_SOA_COLUMN(Eta, eta, float);               //! Eta of hadron
DECLARE_SOA_COLUMN(PtH, ptH, float);               //! Transverse momentum of hadron
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);         //! Pool Bin of event defined using zvtx and multiplicity
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);     //! Global index for the collision
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t); //! Timestamp for the collision
} // namespace hf_assoc_tracks

DECLARE_SOA_TABLE(Hadron, "AOD", "HADRON", //! Associated hadron properties
                  aod::hf_assoc_tracks::Phi,
                  aod::hf_assoc_tracks::Eta,
                  aod::hf_assoc_tracks::PtH,
                  aod::hf_assoc_tracks::PoolBin,
                  aod::hf_assoc_tracks::GIndexCol,
                  aod::hf_assoc_tracks::TimeStamp);

// definition of columns and tables for Dplus-Hadron correlation pairs
namespace hf_correlation_dplus_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);        //! DeltaPhi between D+ and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);        //! DeltaEta between D+ and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);                  //! Transverse momentum of D+
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);        //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);                    //! Invariant mass of D+
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, bool); //! Used in MC-Rec, D+ Signal
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);            //! Pool Bin of event defined using zvtx and multiplicity
} // namespace hf_correlation_dplus_hadron

DECLARE_SOA_TABLE(DplusHadronPair, "AOD", "DPLUSHPAIR", //! D+-Hadrons pairs Informations
                  aod::hf_correlation_dplus_hadron::DeltaPhi,
                  aod::hf_correlation_dplus_hadron::DeltaEta,
                  aod::hf_correlation_dplus_hadron::PtD,
                  aod::hf_correlation_dplus_hadron::PtHadron,
                  aod::hf_correlation_dplus_hadron::PoolBin);

DECLARE_SOA_TABLE(DplusHadronRecoInfo, "AOD", "DPLUSHRECOINFO", //! D+-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_dplus_hadron::MD,
                  aod::hf_correlation_dplus_hadron::SignalStatus);


// definition of columns and tables for Dplus properties
namespace hf_dstar_meson
{
  DECLARE_SOA_COLUMN(PhiDstar, phiDstar, float);
  DECLARE_SOA_COLUMN(EtaDstar, etaDstar, float);
  DECLARE_SOA_COLUMN(PtDstar, ptDstar, float);
  DECLARE_SOA_COLUMN(MDstar, mDstar, float);
  DECLARE_SOA_COLUMN(PoolBin, poolBin, int);
  DECLARE_SOA_COLUMN(TimeStamp,timeStamp, int64_t);
  DECLARE_SOA_COLUMN(IsPrompt,isPrompt,bool);
  DECLARE_SOA_INDEX_COLUMN(Collision, collision);
}// hf_dstar_meson

DECLARE_SOA_TABLE(Dstar, "AOD","DSTAR", //! D* meson properties (This table is not compatible with HfCandDstar & HfSelDstarToD0Pi)
                  aod::hf_dstar_meson::PhiDstar,
                  aod::hf_dstar_meson::EtaDstar,
                  aod::hf_dstar_meson::PtDstar,
                  aod::hf_dstar_meson::MDstar,
                  aod::hf_dstar_meson::PoolBin,
                  aod::hf_dstar_meson::TimeStamp,
                  aod::hf_dstar_meson::CollisionId);

// definition of columns and tables for Dstar-Hadron correlation pair
namespace hf_correlation_dstar_hadron
{
  DECLARE_SOA_COLUMN(DeltaPhiDstar, deltaPhiDstar, float);
  DECLARE_SOA_COLUMN(DeltaEta, deltaEtaDstar, float);
  DECLARE_SOA_COLUMN(MatchingStatus, matchingStatus, bool);
}// hf_correlation_dstar_hadron

DECLARE_SOA_TABLE(DstarHadronPair, "AOD","DSTARHPAIR", // D* Hadrons pairs Informations
                  aod::hf_correlation_dstar_hadron::DeltaPhiDstar,
                  aod::hf_correlation_dstar_hadron::DeltaEta,
                  aod::hf_dstar_meson::PtDstar,
                  aod::hf_assoc_tracks::PtH,
                  aod::hf_dstar_meson::MDstar,
                  aod::hf_dstar_meson::PoolBin);

// Why should we use these two tables below when 
// we already have "HfCandDstarMcRec" & "HfCandDstarMcGen" in "CanidateReconstructionTables.h"
DECLARE_SOA_TABLE(DstarHadronRecoInfo,"AOD","DSTARHRECOINFO", // D* Hadrons pairs Reconstructed Informations
                aod::hf_dstar_meson::IsPrompt,
                aod::hf_correlation_dstar_hadron::MatchingStatus);

DECLARE_SOA_TABLE(DstarHadronGenInfo, "AOD","DSTARHGENINFO", // D* Hadrons pairs Generator level Informations
                aod::hf_dstar_meson::IsPrompt,
                aod::hf_correlation_dstar_hadron::MatchingStatus);

// Note: Table for selection of Lc in a collision
namespace hf_selection_lc_collision
{
DECLARE_SOA_COLUMN(LcSel, lcSel, bool); //! Selection flag for Lc in a collision
} // namespace hf_selection_lc_collision

DECLARE_SOA_TABLE(LcSelection, "AOD", "LCINCOLL", // Selection of Lc in collisions
                  aod::hf_selection_lc_collision::LcSel);

// Table for selection of Dmeson in a collision
namespace hf_selection_dmeson_collision
{
DECLARE_SOA_COLUMN(DmesonSel, dmesonSel, bool); //! Selection flag for D meson in a collision
} // namespace hf_selection_dmeson_collision

DECLARE_SOA_TABLE(DmesonSelection, "AOD", "DINCOLL", // Selection of D meson in collisions
                  aod::hf_selection_dmeson_collision::DmesonSel);
} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_CORRELATIONTABLES_H_
