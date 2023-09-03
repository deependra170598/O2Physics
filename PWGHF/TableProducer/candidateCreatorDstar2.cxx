#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/trackUtilities.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;


namespace o2::aod{
  using HfDstarsWStatus= soa::Join<aod::HfDstars, aod::HfCutStatusDstar>;
  using HfDstarsWStatus_nd_PvRefitInfo =soa::Join<aod::HfDstars, aod::HfCutStatusDstar, aod::HfPvRefitDstar>;
  using Hf2ProngsWStatus = soa:: Join<aod::Hf2Prongs, aod::HfCutStatus2Prong>;

  namespace HFCandDStarProng{
    DECLARE_SOA_COLUMN(CandP,candp,float);
    DECLARE_SOA_COLUMN(CandPt,candpt,float);
    // DECLARE_SOA_COLUMN(DecayL,decayl,float);
    // DECLARE_SOA_COLUMN(DecayLXY,decaylXY,float);
    DECLARE_SOA_COLUMN(DecayLError,decaylerror,float);
    DECLARE_SOA_COLUMN(DecayLXYError,decaylXYerror,float);
    DECLARE_SOA_COLUMN(PCAx,pcax,float);
    DECLARE_SOA_COLUMN(PCAy,pcay,float);
    DECLARE_SOA_COLUMN(PCAz,pcaz,float);
    DECLARE_SOA_COLUMN(D0Mass,d0mass,float);
    DECLARE_SOA_COLUMN(DStarMass,dstarmass,float);
    // DECLARE_SOA_COLUMN(DeltaM,deltam,float);
    DECLARE_SOA_DYNAMIC_COLUMN(DeltaM,deltam,[](float x,float y){return x-y;});
    DECLARE_SOA_COLUMN(Chi2PCA,ch2pca,float);
    DECLARE_SOA_COLUMN(D0P,d0p,float);
    DECLARE_SOA_COLUMN(D0Pt,d0pt,float);
  }
  DECLARE_SOA_TABLE(HfCandDStarProngs,"AOD","HFDSTARCAND",
                    o2::soa::Index<>,
                    hf_track_index::CollisionId,
                    hf_track_index::Prong0Id,
                    hf_track_index::ProngD0Id,
                    hf_track_index::FlagDstarToD0Pi,
                    HFCandDStarProng::CandP,
                    HFCandDStarProng::CandPt,
                    // HFCandDStarProng::DecayL,
                    // HFCandDStarProng::DecayLXY,
                    //Primary vertex
                    collision::PosX,collision::PosY,collision::PosZ,
                    //Secondary vertex
                    HFCandDStarProng::PCAx,HFCandDStarProng::PCAy,HFCandDStarProng::PCAz,
                    //Deacy lenght
                    HFCandDStarProng::DecayLError,HFCandDStarProng::DecayLXYError,
                    //, HFCandDStarProng::D0Mass, HFCandDStarProng::DStarMass, HFCandDStarProng::DeltaM<HFCandDStarProng::DStarMass,HFCandDStarProng::D0Mass>
                    HFCandDStarProng::Chi2PCA,
                    HFCandDStarProng::D0P, HFCandDStarProng::D0Pt
                    );
}

/// Reconstruction of D* decay candidates
struct HfCandidateCreatorDstar{
  Produces<aod::HfCandDStarProngs> DStarCandtables;
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};

  // magnetic field setting from CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  // vertexing
  Configurable<bool> DoRefit{"DoPvRefit",false,"Do Primary Vertex Refit. Optional"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"}; // ........... what is unit of this?
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"}; // ..........What it DZ?
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  

  Service<o2::ccdb::BasicCCDBManager> ccdb; // From utilsBfieldCCDB.h
  o2::base::MatLayerCylSet* lut; // From DCAFitterN.h
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double bz;
  float toMicrometers = 10000.; // from cm to µm

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massK = RecoDecay::getMassPDG(kKPlus);
  double massPiK{0.};
  double massKPi{0.};
  double massD0 = RecoDecay::getMassPDG(pdg::Code::kD0);

  // Refit Container
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVYY{TH1F("hCovPVYY", "2-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVYY{TH1F("hCovSVYY", "2-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVXZ{TH1F("hCovPVXZ", "2-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, -1.e-4, 1.e-4)};
  OutputObj<TH1F> hCovSVXZ{TH1F("hCovSVXZ", "2-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, -1.e-4, 0.2)};
  OutputObj<TH1F> hCovPVZZ{TH1F("hCovPVZZ", "2-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVZZ{TH1F("hCovSVZZ", "2-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  
  OutputObj<TH2F> hDcaXYProngsD0{TH2F("hDcaXYProngsD0", "DCAxy of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngsD0{TH2F("hDcaZProngsD0", "DCAz of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDCAXYPi{TH2F("hDCAXYPi","DCAxy of Soft Pi;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDCAZPi{TH2F("hDCAZPi","DCAz of Soft Pi;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};


  OutputObj<TH1F> hMassD0{TH1F("hMassD0", "D0 candidates;inv. mass (#pi D^{0}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hDeltaMassDStar{TH1F("hDeltaMassDStar", "D* candidates;inv. mass (#pi D^{0}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtPi{TH1F("hPtPi", "#pi candidates;#it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0Prong0{TH1F("hPtD0Prong0", "D^{0} candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0Prong1{TH1F("hPtD0Prong1", "D^{0} candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD0{TH1F("hPtD0", "D^{0} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtDStar{TH1F("hPtDStar","D* candidates;candidate #it{p}_{T} (GeV/#it{c});entries",500, 0., 5.)};
  /// @brief This function initialize the ccdb setting and rin function MatLayerCylSet::rectifyPtrFromFile(..args..)
  /// @param  
  void init (InitContext const&){
    // if(doprocessPvRefit && processNoPvRefit){ //............Warning! remove this if any of this function is removed
    //   LOGP(fatal, "Only one process function between processPvRefit and processNoPvRefit can be enabled at a time.");
    // }
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true); //.............Doubt What does it do?
    ccdb->setLocalObjectValidityChecking();  // set the flag to check object validity before CCDB query
    auto rectification = ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut); // retrieve an object of type T from CCDB as stored under path; will use the timestamp member
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(rectification);
    runNumber =0;
    bz = 0;
  }

  void process(aod::Collisions const&,
   aod::HfDstarsWStatus_nd_PvRefitInfo const& rowsTrackIndexDstar,
   aod::Hf2ProngsWStatus const&,
   aod::TracksWCov const& tracks,
   aod::BCsWithTimestamps const& bcWithTimeStamps){

    // D0-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df;
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over suspected DStar Candidates
    for(const auto& rowTrackIndexDstar : rowsTrackIndexDstar){
      auto trackPi = rowTrackIndexDstar.prong0_as<aod::TracksWCov>();
      auto prongD0 = rowTrackIndexDstar.prongD0_as<aod::Hf2ProngsWStatus>();
      auto trackD0Prong0 = prongD0.prong0_as<aod::TracksWCov>();
      auto trackD0Prong1 = prongD0.prong1_as<aod::TracksWCov>();
      auto collision = rowTrackIndexDstar.collision();

      // Extracts primary vertex position and covariance matrix from a collision
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      // Extracts track parameters and covariance matrix from a track
      auto trackPi_ParVar = getTrackParCov(trackPi);
      // These will be used in DCA Fitter to reconstruct secondary vertex
      auto trackD0Prong0_ParVarPos1 = getTrackParCov(trackD0Prong0); // from trackUtilities.h  
      auto trackD0Prong1_ParVarNeg1 = getTrackParCov(trackD0Prong1);
      
      
      // auto collisionPiId = trackPi.collisionId();
      // auto collisionD0Id = trackD0Prong0.collisionId();
      // LOGF(info, "Pi collision %ld, D0 collision %ld", collisionPiId, collisionD0Id);
      //..................................................Doubt: Should I apply a condition of (collisionPiId == collisionD0Id)

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,  .....................Doubt: Which propagator? Is it propagator to PCA? (Line 762 in traclindexSkimCreator)
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if(runNumber != bc.runNumber()){  // ...............Doubt: In this way, we will always have same Number "0"!
      LOG(info)<< ">>>>>>>>>>>> Current run number: " << runNumber; // Doubt: Still we are using previous value of runNumber=0, Should we update current runNumber?
      initCCDB(bc,runNumber,ccdb,isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);  //Sets up the grp object for magnetic field (w/o matCorr for propagation)
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df.setBz(bz);

      // reconstruct the 2-prong secondary vertex
      if(df.process(trackD0Prong0_ParVarPos1,trackD0Prong1_ParVarNeg1)==0){
        continue; // ................Doubt: what does this process function check?
      }
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      hCovSVYY->Fill(covMatrixPCA[2]);
      hCovSVXZ->Fill(covMatrixPCA[3]);
      hCovSVZZ->Fill(covMatrixPCA[5]);

      // Doubt:................Below, track object are at secondary vertex!
      // < track param propagated to V0 candidate (no check for the candidate validity). propagateTracksToVertex must be called in advance
      auto trackD0Prong0_ParVar0 = df.getTrack(0); // Doubt:..............According to above comment, whare is df.propagateTracksToVertex() called?
      auto trackD0Prong1_ParVar1 = df.getTrack(1); //Doubt:............... why not df.getTrackParamAtPCA()?

      // get track momenta Doubt:................ These momenta are at secondary vertex!
      std::array<float, 3> pVecD0Prong0;
      std::array<float, 3> pVecD0Prong1;
      trackD0Prong0_ParVar0.getPxPyPzGlo(pVecD0Prong0);
      trackD0Prong1_ParVar1.getPxPyPzGlo(pVecD0Prong1);

      
      // This modifies track momenta! ..............Doubt: which point is this momentum at?
      if(DoRefit){
        /// use PV refit
        /// Using it in the *HfCand3ProngBase/HfCand2ProngBase* all dynamic columns shall take it into account
        // coordinates 
        primaryVertex.setX(rowTrackIndexDstar.pvRefitX());
        primaryVertex.setY(rowTrackIndexDstar.pvRefitY());
        primaryVertex.setZ(rowTrackIndexDstar.pvRefitZ());
        // covariance matrix 
        primaryVertex.setSigmaX2(rowTrackIndexDstar.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexDstar.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexDstar.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexDstar.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexDstar.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexDstar.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov(); /// Here covMatrixPV Updated!
      }
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);

      // get track impact parameters
      o2::dataformats::DCA impactParameter0; // GPUROOTCartesianFwd.h
      o2::dataformats::DCA impactParameter1; 
      // Propagating D0 prongs to DCA
      trackD0Prong0_ParVar0.propagateToDCA(primaryVertex,bz,&impactParameter0);
      trackD0Prong1_ParVar1.propagateToDCA(primaryVertex,bz,&impactParameter1);
      
      // What the below function will do ? Can we use this function to get D0Prongs' param at Pv
      // So that momenta of these can be added to soft pi momentum? 
      // trackD0Prong1_ParVar1.propagateParamToDCA()

      // Propagating Soft Pi to DCA
      o2::dataformats::DCA impactParameterPi;
      trackPi_ParVar.propagateToDCA(primaryVertex,bz,&impactParameterPi);
      //trackPi_ParVar.propagateParamToDCA(primaryVertex,bz,&impactParameterPi); //Doubt............ How this line is different from above?
      hDcaXYProngsD0->Fill(trackD0Prong0.pt(),impactParameter0.getY()*toMicrometers);
      hDcaXYProngsD0->Fill(trackD0Prong1.pt(),impactParameter1.getY()*toMicrometers);
      hDcaZProngsD0->Fill(trackD0Prong0.pt(),impactParameter0.getZ()*toMicrometers);
      hDcaZProngsD0->Fill(trackD0Prong1.pt(),impactParameter1.getZ()*toMicrometers);

      hDCAXYPi->Fill(trackPi.pt(),impactParameterPi.getY()*toMicrometers);
      hDCAZPi->Fill(trackPi.pt(),impactParameterPi.getZ()*toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      //Calculates the XX element of a XYZ covariance matrix after rotation of the coordinate system by phi around the z-axis and by minus theta around the new y-axis.
      // ...............Doubt: Why frame was rotated by -theta ?
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      //Calculation Inv mass of D0 with tracks momenta at secondary vertex
      auto arrayMomenta = std::array{pVecD0Prong0,pVecD0Prong1};
      massPiK = RecoDecay::m(arrayMomenta,std::array{massPi,massK});
      massKPi = RecoDecay::m(arrayMomenta,std::array{massK,massPi});

      // Calculation of kinematics for inv mass
      std::array<float,3> pVecPi={trackPi.px(),trackPi.py(),trackPi.pz()};
      auto pVecD0 = RecoDecay::pVec(pVecD0Prong0,pVecD0Prong1);
      auto pD0 = RecoDecay::p(pVecD0);
      // D0 pt vector
      auto pxD0 = pVecD0.at(0);
      auto pyD0 = pVecD0.at(1);
      auto ptVecD0 = std::array{pxD0,pyD0};
      // D0 pt magnitude
      auto ptD0 = RecoDecay::pt(ptVecD0);

      auto InvMassDStar = RecoDecay::m(std::array{pVecPi,pVecD0},std::array{massPi,massD0});

      // Soft pi momentum vector
      std::array<float,3> SoftPipVec;
      trackPi_ParVar.getPxPyPzGlo(SoftPipVec);

      // Dstar momentum vector
      auto pVecDStar = RecoDecay::pVec(pVecD0,SoftPipVec);
      // Magnitude Dstar momentum 
      auto pDStar = RecoDecay::p(pVecDStar);
      auto pxDStar = pVecDStar.at(0);
      auto pyDStar = pVecDStar.at(1);
      // Dstar ptVec
      std::array<float,2> ptVecDStar{pxDStar,pyDStar};
      // DStar pt magnitude
      auto ptDStar = RecoDecay::pt(ptVecDStar);
      
        // LOGF(info,"Table filling start");
      // Fill candidate Table here for DStar
      DStarCandtables(collision.globalIndex(),trackD0Prong0.globalIndex(),trackD0Prong1.globalIndex(),rowTrackIndexDstar.flagDstarToD0Pi(),pDStar,ptDStar,
                      primaryVertex.getX(),primaryVertex.getY(),primaryVertex.getZ(),
                      secondaryVertex.At(0),secondaryVertex.At(1),secondaryVertex.At(2),
                      errorDecayLength,errorDecayLengthXY,chi2PCA,pD0,ptD0);

      if(fillHistograms){
        hDeltaMassDStar->Fill(InvMassDStar-massPiK);
        hDeltaMassDStar->Fill(InvMassDStar-massKPi);
        hMassD0->Fill(massPiK);
        hMassD0->Fill(massKPi);
        hPtD0->Fill(ptD0);
        hPtPi->Fill(RecoDecay::pt(std::array{trackPi.px(),trackPi.py()})); 
        hPtD0Prong0->Fill(RecoDecay::pt(std::array{pVecD0Prong0.at(0),pVecD0Prong0.at(1)}));
        hPtD0Prong1->Fill(RecoDecay::pt(std::array{pVecD0Prong1.at(0),pVecD0Prong1.at(1)}));
        hPtDStar->Fill(ptDStar);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc)};
}