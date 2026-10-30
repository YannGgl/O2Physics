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
///
/// \brief Example how to lopps over V0s ...
/// \author
/// \since

/* Usage :
   o2-analysis-tracks-extra-v002-converter --configuration json://myconfig.json | \
   o2-analysistutorial-task-yann --configuration json://myconfig.json -b --aod-file AO2D.root
*/

#include <cmath>
#include "Common/Core/PID/PIDTOF.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CommonConstants/PhysicsConstants.h>

#include <TPDGCode.h>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include "Common/DataModel/TrackSelectionTables.h"


using namespace o2;
using namespace o2::framework;

using myTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
using myTrack = myTracks::iterator;

using myMcTracks = soa::Join<myTracks, o2::aod::McTrackLabels>;

using mySelectedEvents = soa::Join<aod::Collisions, aod::EvSels>;
using mySelectedEvent = mySelectedEvents::iterator;

using myPIDTracks = soa::Join<myTracks, aod::pidTPCFullPr, aod::pidTOFFullPr>;
//using myPIDTracks = soa::Join<myTracks, aod::pidTPCPr, aod::pidTOFPr>;
using myPIDTrack = myPIDTracks::iterator;
// using myPIDMCTracks = soa::Join<myPIDTracks, o2::aod::McTrackLabels>;
using myPIDMCTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>;
//using myPIDMCTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPr, aod::pidTOFPr, aod::McTrackLabels>;
using myPIDMCTrack = myPIDMCTracks::iterator;

using myTracksDCA = soa::Join<myTracks, aod::TracksDCA>;
using myMcTracksDCA = soa::Join<myMcTracks, aod::TracksDCA>;
using myPIDTracksDCA = soa::Join<myPIDTracks, aod::TracksDCA>;
using myPIDMCTracksDCA = soa::Join<myPIDMCTracks, aod::TracksDCA>;


// Code actually in analysis mode : eta_C (Minv C [1 ; 3.8]GeV/c^2 w/ 2800 bins), correct mass hypothesis, 200 MeV/c pt cut + NO 0.5 GeV/c p cut for all tracks, 
// Efficiency for : Proton (proton = PID 1, anti-proton = PID 2) 
// PID (in conf.json file) : ON
// Dmesons : REMOVED (UNcommented 'if(!Dmeson)' condition)
// beta 0.4 cut = OFF

// Used Event selection : sel8, kNoSameBunchPileup, |zPV| < 10cm
// ES mode (new/old) : NEW
// PID mode (new/old) : NEW
// DCA mode : ON




struct taskYann {

  Configurable<float> cfgZmax{"cfgZmax", 10.f, "Restriction on the PV position |PVz|<zMax (cm)"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.9f, "Restriction in pseudorapidity |eta|<EtaMax"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Restriction in transverse momentum pt>ptMin (GeV/c)"};
  Configurable<bool> cfgPIDMode{"cfgPIDMode", true, "Default mode (true) : PID. Switch to false to deactivate PID."};
  Configurable<float> cfgTPCnSigmaMax{"cfgTPCnSigmaMax", 3.0f, "TPC proton +nSigma"};
  Configurable<float> cfgTPCnSigmaMin{"cfgTPCnSigmaMin", -3.0f, "TPC proton -nSigma"};
  Configurable<float> cfgTOFnSigma{"cfgTOFnSigma", 3.0f, "TOF proton nSigma"};
  Configurable<int> cfgitsNCls{"cfgitsNCls", 7, "Nb of clusters of ITS (7 max)"};

  OutputObj<TH1F> hVtxBeforeES{
    TH1F("hVtxBeforeES", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};
  OutputObj<TH1F> hVtxAfterES{
    TH1F("hVtxAfterES", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};
  //OutputObj<TF1> hExpSignalPrTPC{
    //TF1("hExpSignalPrTPC", ";",)};
  

  /*
  OutputObj<TH2F> hTpc_initial{
    TH2F("hTpc_initial", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};
  OutputObj<TH2F> hTof_initial{
    TH2F("hTof_initial", "TOF; Momentum (GeV); beta", 500, -5., 5., 220, 0., 1.1)};

  OutputObj<TH2F> hTpc1{
    TH2F("hTpc1", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};
  OutputObj<TH2F> hTof1{
    TH2F("hTof1", "TOF; Momentum (GeV); beta", 500, -5., 5., 220, 0., 1.1)};

  
  OutputObj<TH2F> hTpc2{
    TH2F("hTpc2", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};
  OutputObj<TH2F> hTof2{
    TH2F("hTof2", "TOF; Momentum (GeV); beta", 500, -5., 5., 220, 0., 1.1)};
  */

  ConfigurableAxis pTBins{"pT-bins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.4f, 2.8f, 3.2f, 3.6f, 4.0f, 4.8f, 5.6f, 6.5f, 7.5f, 9.0f, 11.0f, 13.0f, 15.0f, 19.0f, 23.0f, 30.0f, 40.0f, 50.0f}, "pt axis for analysis"};

  void init(o2::framework::InitContext& /*ic*/)
  {
    AxisSpec pTAxis = {pTBins, "#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c})"};

    registryMC.add("hEtaCToppBarPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", {HistType::kTH1F, {pTAxis}});
    registryMC.add("hEtaCAsMotherOfpPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });
    registryMC.add("hEtaCAsMotherOfpBarPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });
    registryMC.add("hpWithEtaCMotherPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });
    registryMC.add("hpBarWithEtaCMotherPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });

    registryMC_Gen.add("hEtaCToppBarPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", {HistType::kTH1F, {pTAxis}});
    registryMC_Gen.add("hEtaCAsMotherOfpPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", {HistType::kTH1F, {pTAxis}});
    registryMC_Gen.add("hEtaCAsMotherOfpBarPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", {HistType::kTH1F, {pTAxis}});
    registryMC_Gen.add("hpWithEtaCMotherPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });
    registryMC_Gen.add("hpBarWithEtaCMotherPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon}", { HistType::kTH1F, { pTAxis } });
  }


  HistogramRegistry registry{
    "registry",
    {
      // {"hMinv_Total", "Invariant mass; m (GeV); Entries", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      {"hMinv_PID", "; #it{M_{candidate}  }[GeV/#it{c^{2}}]; #it{Entries  per  1 }MeV/#it{c^{2}}", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      {"hTPC_Initial", "TPC without PID; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4. }, { 200, 0., 1000 } } } },
      {"hTOF_Initial", "TOF without PID; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      {"hTOFInitialBetterRes", "TOF without PID; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 1000, -5., 5. }, { 440, 0., 1.1 } } } },
      {"hTPC", "TPC; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4. }, { 200, 0., 1000 } } } },
      {"hTOF", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      {"hTOFBetterRes", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 1000, -5., 5. }, { 440, 0., 1.1 } } } },
      //{"hTOF_over04", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      //{"hTOF_under04", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      {"hMinvPT", "; #it{M_{p#bar{p}}  }[GeV/#it{c^{2}}]; #it{p_{T}  }[Gev/#it{c}]#it{ of candidate}", { HistType::kTH2F, { { /*500*/ /*2800*/200, 1., /*1.5 2.3*/ 3.8 }, { 50, 0., 10. } } } },
      {"hPTCandidate", "; #it{p_{T}  }[Gev/#it{c}]; #it{Entries  per  100 }MeV/#it{c}", { HistType::kTH1F, { { 330, 0., 33. } } } },
      //{"hTPCnSigma", "; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4. }, { 200, 0., 1000 } } } },
      {"hTPCExpSignalPr", "; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 800, -4., 4. }, { 400, 0., 1000 } } } },
      {"hdcaXY_Before", "; DCA xy [cm]; Entries per 10#mum", { HistType::kTH1F, { {4000, -1000., 1000.} } } },
      {"hdcaXY_After", "; DCA xy [cm]; Entries per 10#mum", { HistType::kTH1F, { {4000, -1000., 1000.} } } },
      {"hdcaZ_Before", "; DCA z [cm]; Entries per 10#mum", { HistType::kTH1F, { {4000, -1000., 1000.} } } },
      {"hdcaZ_After", "; DCA z [cm]; Entries per 10#mum", { HistType::kTH1F, { {4000, -1000., 1000.} } } },
      {"hDCAXYPTDepCut", "; p_{T}; DCA xy [cm]", { HistType::kTH2F, { { 330, 0., 33. }, {1000, 0., 1.} } } },
    }
  };

  HistogramRegistry registryMC{
    "registryMC",
    {
      {"hTOFMC1", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5.}, { 220, 0., 1.1 } } } },
      {"hTPCMC1", "TPC; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4.}, { 200, 0., 1000 } } } },
      //{"hTofMC2", "TOF; Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5.}, { 220, 0., 1.1 } } } },
      //{"hTpcMC2", "TPC; Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4.}, { 200, 0., 1000 } } } },
      {"hMinvMC", "Invariant mass; m  [GeV /c^{2}]; Entries  per  1 MeV/c^{2}", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      {"hMinvMC_MatchedEtaC", "Invariant mass; m [GeV /c^{2}]; Entries  per  1 MeV/c^{2}", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      {"hMinvMC_MatchedJPsi", "Invariant mass; m [GeV /c^{2}]; Entries  per  1 MeV/c^{2}", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      {"hMinvMC_MatchedPsi2S", "Invariant mass; m [GeV /c^{2}]; Entries  per  1 MeV/c^{2}", { HistType::kTH1F, { { /*500*/ 2800, 1., /*1.5 2.3*/ 3.8 } } } },
      //{"hPt_Matched_p", "Proton transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hpPt_MatchedWrong", "Particles transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      //{"hpBarPt_MatchedWrong", "Particles transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hpPt_MatchedRight", "Particles transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hEtaCAsMotherPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon,  per 200 }MeV/#it{c}", { HistType::kTH1F, { { 165, 0., 33. } } } },
      //{"hEtaCToppBarPt_Matched", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon,  per 200 }MeV/#it{c}", { HistType::kTH1F, { { 165, 0., 33. } } } },
      {"hMinvPTMC", "; #it{M_{p#bar{p}}  }[GeV/#it{c^{2}}]; #it{p_{T}  }[Gev/#it{c}]#it{ of candidate}", { HistType::kTH2F, { { /*500*/ /*2800*/200, 1., /*1.5 2.3*/ 3.8 }, { 50, 0., 10. } } } },
      {"hPTCandidateMC", "; #it{p_{T}  }[Gev/#it{c}]; #it{Entries  per  100 }MeV/#it{c}", { HistType::kTH1F, { { 330, 0., 33. } } } },
    }
  };

  HistogramRegistry registryTime{
    "registryTime",
    {
      {"hCollTimeRes", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      {"hTrackTimeRes", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      //{"hCollTimeRes_beta_over04", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      //{"hTrackTimeRes_beta_over04", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      //{"hCollTimeRes_beta_under04", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      //{"hTrackTimeRes_beta_under04", ";Time resolution [ns];Entries", { HistType::kTH1F, { { 5000, 0., 1000 } } } },
      //{"hTTR_beta", ";Time resolution [ns]; #beta", { HistType::kTH2F, { { 5000, 0., 1000 }, { 220, 0., 1.1 } } } },
    }
  };
  

  

  // Collision selector
  // True when collision(=PV) within Zmax(=10cm) distance of geometrical center
  // False otherwise
  template <typename CollisionTypeContainer>
  bool isCollisionAccepted(CollisionTypeContainer const& collision)
  {
    auto z = collision.posZ();
    if (std::abs(z) > cfgZmax)
      return false;

    // Some other selections

    //if (collision.collisionTimeRes() > 3.)
      //return false;

    return true;
  }



  // Tracks DCA
  template <typename TrackType>
  bool isPrimaryTrack(TrackType const& track) 
  {
    //float dcaXYPtCut = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
    float dcaPtCut = 0.004f + 0.013f / track.pt();

    if (std::abs(track.dcaXY()) < dcaPtCut && std::abs(track.dcaZ()) < dcaPtCut)
      return true;
    return false;

  }


  // beta = v*c (Lorentz coeff) where v = d/t (particle velocity)
  template <typename TrackType>
  float tofBeta(TrackType const& track)
  {
    auto length = track.length();
    auto tofSignal = o2::pid::tof::TOFSignal<TrackType>::GetTOFSignal(track);
    return length / tofSignal * o2::constants::physics::invLightSpeedCm2PS;
  }


  // Track-quality selector
  // False if one of the conditions for the track of the particle is not satisfied :
  // True if all conditions are satisfied
  template <typename TrackType>
  bool isTrackAccepted(TrackType const& track)
  {
    //if (!track.isWithinBeamPipe())
    if (!track.hasITS()) // ITS detectors = 7 cylinders (tracking) 
      return false;

    if (!track.hasTPC()) // TPC detectors
      return false;

    if (std::abs(track.eta()) > cfgEtaMax) // rapidity (eta/y) 
      return false;

    if (track.pt() < cfgPtMin) // transverse momentum pt
      return false;

    if (track.itsNCls() < cfgitsNCls ) // nb of clusters of ITS (7 max)
      return false;                   // "we require full track of primary particle"

    // Some other selections

    //if (track.trackTimeRes() < 0.8)
      //return false;


    // Implementer "Coups dans 1 ou 3 premieres couches ITS + 5,6 coups global"

    return true;
  }

  // Function for selecting tracks under straight lines f(p) = ap + b 
  // combining the use of TPC and TOF detectors
  template <typename TrackType>
  bool isTrackProton(TrackType const& track)
  {
    auto mom = track.tpcInnerParam();
    // auto sign = track.sign();

    if (!track.hasTOF() || mom<0.5) { // TPC
      auto dedx = track.tpcSignal();
      double a[4] = {-5000, -208, -7.5, -108};
      double b[4] = {2000, 250, 70, 180};
      if ( dedx >= a[0]*mom + b[0]
      && dedx >= a[1]*mom + b[1]
      && dedx >= a[2]*mom + b[2]
      && dedx >= a[3]*mom + b[3] ) {
        // registry.get<TH2>(HIST("hTPC"))->Fill(sign*mom, dedx);
        return true;
      }
      return false;
    }
    else { // TOF
      auto beta = tofBeta(track);
      double a[4] = {0.2, 0.01, 1., 0.066};
      double b[4] = {0.6, 0.95, 0., 0.79};
      //double b[5] = {0.6, 0.95, 0., 0.79, 0.4}; // removes beta < 0.4 
      //double b[6] = {0.6, 0.95, 0., 0.79, 0.425, 0.330}; // To isolate first horizontal band 
      
      if ( beta <= a[0]*mom + b[0]
      && beta <= a[1]*mom + b[1]
      && beta <= a[2]*mom + b[2]
      && beta <= a[3]*mom + b[3]
      //&& beta >= b[4] 
      //&& beta <= b[4] // To isolate first horizontal band 
      //&& beta >= b[5] // To isolate first horizontal band 
      ) {
        // registry.get<TH2>(HIST("hTOF"))->Fill(sign*mom, beta);
        return true;
      }
      return false;

      /*if (beta <= 0.425 && beta >= 0.105) {
        if ( (beta <= 0.335 && beta >= 0.262)
        || (beta <= 0.198 && beta >= 0.190)
        || (beta <= 0.135 && beta >= 0.128)
        ) {
          return false;
        }
        return true;
      }
      return false;*/
    }
  }

  // Function for PID
  template <typename TrackType>
  bool isTrackProtonNew(TrackType const& track)
  {
    auto mom = track.tpcInnerParam();

    if (!track.hasTOF() || mom<0.5) { // TPC
      if (track.tpcNSigmaPr() < cfgTPCnSigmaMax && track.tpcNSigmaPr() > cfgTPCnSigmaMin)
        return true;
      return false;
    }
    else { // TOF
      if (std::abs(track.tofNSigmaPr()) < cfgTOFnSigma)
        return true;
      return false;
    }

  }
  
  // Function to select MC Truth for antiproton
  template <typename McParticleType>
  bool isTrueProton(McParticleType const& mcParticle, const PDG_t& type)
  {
    // Type : kProton or kProtonBar
    if (mcParticle.pdgCode() == type)
      return true;

    return false;
  }
  

  // Function for invariant mass calculation, returns m^2
  template <typename TrackType>
  double InvariantMass(TrackType const& track1, TrackType const& track2, double const& m1, double const& m2)
  {
    auto p1 = track1.p();
    auto p2 = track2.p();
    auto E1 = sqrt(m1*m1 + p1*p1);
    auto E2 = sqrt(m2*m2 + p2*p2);

    return (E1+E2)*(E1+E2) - ( (track1.px() + track2.px())*(track1.px() + track2.px()) + 
            (track1.py() + track2.py())*(track1.py() + track2.py()) + 
            (track1.pz() + track2.pz())*(track1.pz() + track2.pz()) );
  }

  // Function for transverse momentum calculation for candidate (i.e. from a pair of reconstructed tracks), returns pT
  template <typename TrackType>
  double pTcandidate(TrackType const& track1, TrackType const& track2) 
  {
    return sqrt( (track1.px() + track2.px())*(track1.px() + track2.px()) +
                 (track1.py() + track2.py())*(track1.py() + track2.py()) );
  }





  // Main method of the task
  // Subscribing to the table Collision, iterating over all tracks of a given collision
  template <typename CollisionTypeContainer, typename TrackTypeContainer>
  void process(CollisionTypeContainer const& collision, TrackTypeContainer const& tracks)
  {

    // "Loading bar" (print smthg each X collisions)
    static int ncol = 0;
    if (ncol % 10000 == 0)
      LOG(info) << "Collision: " << ncol;
    ncol++;

    /*if constexpr (requires {collision.has_sel8();}) {
      if (!collision.sel8()) {
        return;
      }
      LOG(info) << "Sel8() OK?";
    }*/

    if (!isCollisionAccepted(collision))
      return;

    // Collision counter...
    hVtxAfterES->Fill(collision.posZ());

    // Collisions time resolution distribution
    //registryTime.get<TH1>(HIST("hCollTimeRes"))->Fill(collision.collisionTimeRes());
  

    // First (outer) loop over all tracks in a given collision
    for (auto const& track1 : tracks) {
      
      /*
      if constexpr (requires {track1.has_mcParticle();}) {
        if (track1.has_mcParticle() /*&& isTrackProtonNew(track1)*//*) {
          auto const& p = track1.mcParticle();

          if ((std::abs(p.eta()) < cfgEtaMax) && (p.pt() > cfgPtMin) && p.has_mothers()) {
            auto const& mother = p.template mothers_first_as<aod::McParticles>();

            // W/out PID, just particles with EtaC as mother
            if (mother.pdgCode() == 441) {
              registryMC.get<TH1>(HIST("hEtaCAsMotherPt_Matched"))->Fill(mother.pt());
              
              // PID proton
              if (track1.sign() > 0 && isTrackProtonNew(track1)) 
                registryMC.get<TH1>(HIST("hEtaCAsMotherOfpPt_Matched"))->Fill(mother.pt());
              
              // PID anti-proton
              if (track1.sign() < 0 && isTrackProtonNew(track1)) 
                registryMC.get<TH1>(HIST("hEtaCAsMotherOfpBarPt_Matched"))->Fill(mother.pt());
          
            }
          }
        }
      }
      */

      
      registry.get<TH1>(HIST("hdcaXY_Before"))->Fill(track1.dcaXY());
      registry.get<TH1>(HIST("hdcaZ_Before"))->Fill(track1.dcaZ());

      float dcaXYPtCut = 0.0105f + 0.0350f / pow(track1.pt(), 1.1f);
      registry.get<TH2>(HIST("hDCAXYPTDepCut"))->Fill(track1.pt(), dcaXYPtCut);
      

      // Track quality
      if (!isTrackAccepted(track1))
        continue;
      

      // DCA mode
      if (!isPrimaryTrack(track1))
        continue;

      registry.get<TH1>(HIST("hdcaXY_After"))->Fill(track1.dcaXY());
      registry.get<TH1>(HIST("hdcaZ_After"))->Fill(track1.dcaZ());
      

      

      if constexpr (requires {track1.has_mcParticle();}) {
        if (track1.has_mcParticle()) {
          auto const& p = track1.mcParticle();

          if ((std::abs(p.eta()) < cfgEtaMax) && (p.pt() > cfgPtMin) && p.has_mothers()) {
            auto const& mother = p.template mothers_first_as<aod::McParticles>();

            // W/out PID, just particles with EtaC as mother
            if (mother.pdgCode() == 441 && std::abs(mother.y()) < 0.5) {
              registryMC.get<TH1>(HIST("hEtaCAsMotherPt_Matched"))->Fill(mother.pt());
              
              // PID proton
              // PID mode
              if (track1.sign() > 0 && isTrackProtonNew(track1)/*isTrackProton(track1)*/ && isTrueProton(p, kProton)) {
                registryMC.get<TH1>(HIST("hEtaCAsMotherOfpPt_Matched"))->Fill(mother.pt());
                registryMC.get<TH1>(HIST("hpWithEtaCMotherPt_Matched"))->Fill(p.pt());
              }
              
              // PID anti-proton
              if (track1.sign() < 0 && isTrackProtonNew(track1)/*isTrackProton(track1)*/ && isTrueProton(p, kProtonBar)) {
                registryMC.get<TH1>(HIST("hEtaCAsMotherOfpBarPt_Matched"))->Fill(mother.pt());
                registryMC.get<TH1>(HIST("hpBarWithEtaCMotherPt_Matched"))->Fill(p.pt());
              }
            }
          }
        }
      }




      // Tracks time resolution distribution
      //registryTime.get<TH1>(HIST("hTrackTimeRes"))->Fill(track1.trackTimeRes());

      auto mom1 = track1.tpcInnerParam();
      auto dedx1 = track1.tpcSignal();
      auto beta1 = tofBeta(track1);
      auto sign1 = track1.sign();

      /*
      // Using nSigma for TPC
      if (std::abs(track1.tpcNSigmaPr()) > 3.0) {
        // registry.get<TH2>(HIST("hTPCnSigma"))->Fill(sign1*mom1, dedx1);
        continue;
      }
      */
      
      // TPC, TOF histograms after "primary track quality", but before PID
      if (track1.hasTOF() && mom1 > 0.5) {
        registry.get<TH2>(HIST("hTOF_Initial"))->Fill(sign1*mom1, beta1);
        registry.get<TH2>(HIST("hTOFInitialBetterRes"))->Fill(sign1*mom1, beta1);
      }
      registry.get<TH2>(HIST("hTPC_Initial"))->Fill(sign1*mom1, dedx1);

      // Sign condition to select particles/antiparticles
      // Select : sign > 0 = Antiparticles ; sign < 0 = Particles
      if (sign1 < 0)
        continue;



     

      // MC Matching of generated (anti-)protons, counts nb of generated particles
      if constexpr (requires {track1.has_mcParticle();}) {
        if (track1.has_mcParticle()) {
          auto const& p1 = track1.mcParticle();

          if (p1.pt() > cfgPtMin) {
            // Select : kProtonBar = Anti-proton efficiency ; kProton = Proton efficiency
            if (p1.pdgCode() == kProton) 
              registryMC.get<TH1>(HIST("hpPt_MatchedWrong"))->Fill(p1.pt());
            //if (p1.pdgCode() == kProtonBar) 
              //registryMC.get<TH1>(HIST("hpBarPt_MatchedWrong"))->Fill(p1.pt());
          }
        }
      }

      // PID
      if (cfgPIDMode) {
        // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
        //if (!isTrackProton(track1))
        if (!isTrackProtonNew(track1))
          continue;
      }



      // Efficiency of proton par moved to processAntonin()..
      // MC Matching of generated (anti-)protons, counts nb of generated particles
      if constexpr (requires {track1.has_mcParticle();}) {
        if (track1.has_mcParticle()) {
          auto const& p1 = track1.mcParticle();

          if (p1.pt() > cfgPtMin) {
            // Select : kProtonBar = Anti-proton efficiency ; kProton = Proton efficiency
            if (p1.pdgCode() == kProton) 
              registryMC.get<TH1>(HIST("hpPt_MatchedRight"))->Fill(p1.pt());
          }
        }
      }



      // TPC and TOF histograms after PID 1
      if (!track1.hasTOF() || mom1<0.5) { // TPC
        registry.get<TH2>(HIST("hTPC"))->Fill(sign1*mom1, dedx1);
        double tpcExpSignalPr = track1.tpcSignal() - (track1.tpcNSigmaPr() * track1.tpcExpSigmaPr());
        registry.get<TH2>(HIST("hTPCExpSignalPr"))->Fill(sign1*mom1, tpcExpSignalPr);
      }
      else { // TOF
        registry.get<TH2>(HIST("hTOF"))->Fill(sign1*mom1, beta1);
        registry.get<TH2>(HIST("hTOFBetterRes"))->Fill(sign1*mom1, beta1);
        registryTime.get<TH1>(HIST("hCollTimeRes"))->Fill(collision.collisionTimeRes());
        //registryTime.get<TH2>(HIST("hTTR_beta"))->Fill(track1.trackTimeRes(), beta1);
        /*if (beta1 < 0.4) {
          registryTime.get<TH1>(HIST("hCollTimeRes_beta_under04"))->Fill(collision.collisionTimeRes());
          registryTime.get<TH1>(HIST("hTrackTimeRes_beta_under04"))->Fill(track1.trackTimeRes());
          registry.get<TH2>(HIST("hTOF_under04"))->Fill(sign1*mom1, beta1);
        }
        else {
          registryTime.get<TH1>(HIST("hCollTimeRes_beta_over04"))->Fill(collision.collisionTimeRes());
          registryTime.get<TH1>(HIST("hTrackTimeRes_beta_over04"))->Fill(track1.trackTimeRes());
          registry.get<TH2>(HIST("hTOF_over04"))->Fill(sign1*mom1, beta1);
        }*/
      }

      // Tracks time resolution distribution
      registryTime.get<TH1>(HIST("hTrackTimeRes"))->Fill(track1.trackTimeRes());

      // MC TPC and TOF after PID 1
      if constexpr (requires {track1.has_mcParticle();}) {
        if (track1.has_mcParticle()) {
          auto const& p1 = track1.mcParticle();
          if (isTrueProton(p1, kProton)) {
            if (track1.hasTOF() && mom1 > 0.5)
              registryMC.get<TH2>(HIST("hTOFMC1"))->Fill(sign1*mom1, beta1);
            registryMC.get<TH2>(HIST("hTPCMC1"))->Fill(sign1*mom1, dedx1);
          }
        }
      }

      /*
      static int ntracks = 0;
      if (ntracks % 10000 == 0)
        // LOG(info) << "Tracks : " << ntracks;
      ntracks++;
      */

      // Second (inner) loop over all tracks in a given collision
      for (auto const& track2 : tracks) {
        
        if (!isTrackAccepted(track2)) 
          continue;

        // DCA mode
        if (!isPrimaryTrack(track2))
          continue;

        auto sign2 = track2.sign();
        if (sign2 == sign1)
          continue;

        double m = 0;
        m = sqrt(InvariantMass(track1, track2, 0.9383, /*0.1396*/ 0.9383));

        // PID for track2, only for etaC meson analysis
        // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
        //if (!isTrackProton(track2))
        if (!isTrackProtonNew(track2))
          continue;
        
        // Dmeson to discard D+/- mesons from invariant mass for Enriched MC data
        bool Dmeson = false;
        

        // MC Truth and Matched block
        // All conditions on track2 only for etaC meson
        // int idx=-1;
        if constexpr (requires {track1.has_mcParticle();}) {
          int idxEtaC=-1;
          int idxJpsi=-1;
          int idxPsi2S=-1;
          if (track1.has_mcParticle() && track2.has_mcParticle()) {
            auto const& p1 = track1.mcParticle();

            if (p1.has_mothers()) {
              auto const& mother = p1.template mothers_first_as<aod::McParticles>();

              if (std::abs(mother.pdgCode()) == 421)
                  Dmeson = true;
            }
            auto const& p2 = track2.mcParticle();

            // For argument "type" of isTrueProton(track, type)
            // Select : kProton = proton PID ; kProtonBar = antiproton PID
            if (isTrueProton(p1, kProton) && isTrueProton(p2, kProtonBar)) {
              registryMC.get<TH1>(HIST("hMinvMC"))->Fill(m);
              registryMC.get<TH1>(HIST("hPTCandidateMC"))->Fill(pTcandidate(track1, track2));
              registryMC.get<TH2>(HIST("hMinvPTMC"))->Fill(m, pTcandidate(track1, track2));
            }
            if (isTrueProton(p1, kProton) && p1.has_mothers() && (std::abs(p1.eta()) < cfgEtaMax) && (p1.pt() > cfgPtMin)) {
              auto const& mother = p1.template mothers_first_as<aod::McParticles>();

              // Select : kLambda0Bar = antiproton PID ; kLambda0 = proton PID
              //if (mother.pdgCode() == kLambda0)
              if (mother.pdgCode() == 441 && std::abs(mother.y()) < 0.5 && mother.pt() > cfgPtMin) // PDG code of EtaC meson
                idxEtaC = mother.globalIndex();
              if (mother.pdgCode() == 443) // PDG code of J/Psi meson
                idxJpsi = mother.globalIndex();
              if (mother.pdgCode() == 100443) // PDG code of Psi2S meson
                idxPsi2S = mother.globalIndex();
            }
          }
          if (track2.has_mcParticle()) {
            auto const& p2 = track2.mcParticle();

            if (isTrueProton(p2, kProtonBar) && p2.has_mothers() && (std::abs(p2.eta()) < cfgEtaMax) && (p2.pt() > cfgPtMin)) {
              auto const& mother = p2.template mothers_first_as<aod::McParticles>();

              if (mother.globalIndex() == idxEtaC && std::abs(mother.y()) < 0.5 && mother.pt() > cfgPtMin) {
                registryMC.get<TH1>(HIST("hMinvMC_MatchedEtaC"))->Fill(m);
                registryMC.get<TH1>(HIST("hEtaCToppBarPt_Matched"))->Fill(mother.pt());
                //LOG(info) << mother.pdgCode();
                //LOG(info) << Form("mother.pdgCode() : %d", mother.pdgCode());
              }
              if (mother.globalIndex() == idxJpsi)
                registryMC.get<TH1>(HIST("hMinvMC_MatchedJPsi"))->Fill(m);
              if (mother.globalIndex() == idxPsi2S)
                registryMC.get<TH1>(HIST("hMinvMC_MatchedPsi2S"))->Fill(m);
                // LOG(info) << "idx: " << idx;
            }
          }
        }

        // DMesons removal for MC data (enriched in D Mesons)
        if (!Dmeson) {
          registry.get<TH1>(HIST("hMinv_PID"))->Fill(m);
          registry.get<TH1>(HIST("hPTCandidate"))->Fill(pTcandidate(track1, track2));
          registry.get<TH2>(HIST("hMinvPT"))->Fill(m, pTcandidate(track1, track2));
        }

        
      }
    }
  }


  // Antonin's analysis

  HistogramRegistry registryAntonin{
    "registryAntonin",
    {
      {"hpTPC", "p(TPC); Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4. }, { 200, 0., 1000 } } } },
      {"hpbarTPC", "#bar{p}(TPC); Momentum p  [GeV / c]; dE/dx", { HistType::kTH2F, { { 400, -4., 4. }, { 200, 0., 1000 } } } },
      {"hpTOF", "p(TOF); Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      {"hpbarTOF", "#bar{p}(TOF); Momentum p  [GeV / c]; #beta", { HistType::kTH2F, { { 500, -5., 5. }, { 220, 0., 1.1 } } } },
      {"hMinv_TPCTPC", "p(TPC) + #bar{p}(TPC); M_{inv} ( p#bar{p} )  [GeV /c^{2}]; Entries  per  1 MeV/c", { HistType::kTH1F, { { 2800, 1., 3.8 } } } },
      {"hMinv_TPCTOF", "p(TPC) + #bar{p}(TOF); M_{inv} ( p#bar{p} )  [GeV /c^{2}]; Entries  per  1 MeV/c", { HistType::kTH1F, { { 2800, 1., 3.8 } } } },
      {"hMinv_TOFTPC", "p(TOF) + #bar{p}(TPC); M_{inv} ( p#bar{p} )  [GeV /c^{2}]; Entries  per  1 MeV/c", { HistType::kTH1F, { { 2800, 1., 3.8 } } } },
      {"hMinv_TOFTOF", "p(TOF) + #bar{p}(TOF); M_{inv} ( p#bar{p} )  [GeV /c^{2}]; Entries  per  1 MeV/c", { HistType::kTH1F, { { 2800, 1., 3.8 } } } },
      {"hpPt_Matched", "Particles transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hpBarPt_Matched", "Particles transverse momentum; p_{T} [GeV / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
    }
  };
  
  template <typename CollisionTypeContainer, typename TrackTypeContainer>
  void processAntonin(CollisionTypeContainer const& collision, TrackTypeContainer const& tracks)
  {

    if (!isCollisionAccepted(collision))
      return;

    
    for (auto const& track1 : tracks) {
      
      // Track quality
      if (!isTrackAccepted(track1))
        continue;

      // DCA mode
      if (!isPrimaryTrack(track1))
        continue;

      auto mom1 = track1.tpcInnerParam();
      auto dedx1 = track1.tpcSignal();
      auto beta1 = tofBeta(track1);
      auto sign1 = track1.sign();

      // p < 0.5 GeV/c cut
      //if (mom1 < 0.5)
        //continue;

      // Histograms TPC, TOF
      if (sign1 > 0) { // Particles
        // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
        if (isTrackProtonNew(track1)/*isTrackProton(track1)*/) { // Proton
          
          // MC Matching of generated (anti-)protons, counts nb of generated particles
          if constexpr (requires {track1.has_mcParticle();}) {
            if (track1.has_mcParticle()) {
              auto const& p1 = track1.mcParticle();

              if (p1.pt() > cfgPtMin && std::abs(p1.eta()) < cfgEtaMax) {
                // Select : kProtonBar = Anti-proton efficiency ; kProton = Proton efficiency
                if (p1.pdgCode() == kProton) 
                  registryAntonin.get<TH1>(HIST("hpPt_Matched"))->Fill(p1.pt());
              }
            }
          }

          if (!track1.hasTOF() || mom1<0.5) { // TPC
            registryAntonin.get<TH2>(HIST("hpTPC"))->Fill(sign1*mom1, dedx1);
          }
          else { // TOF
            registryAntonin.get<TH2>(HIST("hpTOF"))->Fill(sign1*mom1, beta1);
          }
        }
      }
      else { // Anti-particles
        // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
        if (isTrackProtonNew(track1)/*isTrackProton(track1)*/) { // Anti-proton
          
          // MC Matching of generated (anti-)protons, counts nb of generated particles
          if constexpr (requires {track1.has_mcParticle();}) {
            if (track1.has_mcParticle()) {
              auto const& p1 = track1.mcParticle();

              if (p1.pt() > cfgPtMin && std::abs(p1.eta()) < cfgEtaMax) {
                // Select : kProtonBar = Anti-proton efficiency ; kProton = Proton efficiency
                if (p1.pdgCode() == kProtonBar) 
                  registryAntonin.get<TH1>(HIST("hpBarPt_Matched"))->Fill(p1.pt());
              }
            }
          }
          
          if (!track1.hasTOF() || mom1<0.5) { // TPC
            registryAntonin.get<TH2>(HIST("hpbarTPC"))->Fill(sign1*mom1, dedx1);
          }
          else { // TOF
            registryAntonin.get<TH2>(HIST("hpbarTOF"))->Fill(sign1*mom1, beta1);
          }
        }        
      }


      // Select : sign > 0 = Antiparticles ; sign < 0 = Particles
      if (sign1 < 0)
        continue;

      // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
      if (!isTrackProtonNew(track1))
      //if (!isTrackProton(track1))
        continue;

      for (auto const& track2 : tracks) {

        // Track quality
        if (!isTrackAccepted(track2))
          continue;

        // DCA mode
        if (!isPrimaryTrack(track2))
          continue;

        auto sign2 = track2.sign();
        auto mom2 = track2.tpcInnerParam();

        if (sign2 == sign1)
          continue;

        // p < 0.5 GeV/c cut
        //if (mom2 < 0.5)
          //continue;


        // PID for track2, only for etaC meson analysis
        // PID mode : isTrackProton = OLD ; isTrackProtonNew = NEW
        if (!isTrackProtonNew(track2))
        //if (!isTrackProton(track2))
          continue;

        // Dmeson to discard D+/- mesons from invariant mass for Enriched MC data
        bool Dmeson = false;

        double m = 0;
        m = sqrt(InvariantMass(track1, track2, 0.9383, /*0.1396*/ 0.9383));


        if constexpr (requires {track1.has_mcParticle();}) {
          if (track1.has_mcParticle()) {
            auto const& p1 = track1.mcParticle();
            if (p1.has_mothers()) {
              auto const& mother = p1.template mothers_first_as<aod::McParticles>();
              if (std::abs(mother.pdgCode()) == 421)
                Dmeson = true;
            }
          }
        }


        // DMesons removal for MC data (enriched in D Mesons)
        if (!Dmeson) {
          if (!track1.hasTOF() || mom1<0.5) { // TPC 1
            if(!track2.hasTOF() || mom2<0.5)  // TPC 2
              registryAntonin.get<TH1>(HIST("hMinv_TPCTPC"))->Fill(m);
            else // TOF 2
              registryAntonin.get<TH1>(HIST("hMinv_TPCTOF"))->Fill(m);
          }
          else { // TOF 1
            if(!track2.hasTOF() || mom2<0.5)  // TPC 2
              registryAntonin.get<TH1>(HIST("hMinv_TOFTPC"))->Fill(m);
            else // TOF 2
              registryAntonin.get<TH1>(HIST("hMinv_TOFTOF"))->Fill(m);
          }
        }
        

      }
    }
  }


  // ES mode = OLD
  // PID mode : myTracks = OLD ; myPIDTracks = NEW
  // DCA mode : myTracks = OFF ; myTracksDCA = ON, CHECK !!!
  void processData(aod::Collision const& collision, myPIDTracksDCA const& tracks)
  //void processData(aod::Collision const& collision, myTracksDCA const& tracks)
  //void processData(aod::Collision const& collision, myPIDTracks const& tracks)
  //void processData(aod::Collision const& collision, myTracks const& tracks)
  {
    process(collision, tracks);
    processAntonin(collision, tracks);
  }
  PROCESS_SWITCH(taskYann, processData, "Process data", true);

  
  // ES mode = OLD
  // PID mode : myMcTracks = OLD ; myPIDMCTracks = NEW
  // DCA mode : myMcTracks = OFF ; myMcTracksDCA = ON, CHECK !!!
  void processMcRec(aod::Collision const& collision, myPIDMCTracksDCA const& tracks, aod::McParticles&)
  //void processMcRec(aod::Collision const& collision, myMcTracksDCA const& tracks, aod::McParticles&)
  //void processMcRec(aod::Collision const& collision, myPIDMCTracks const& tracks, aod::McParticles&)
  //void processMcRec(aod::Collision const& collision, myMcTracks const& tracks, aod::McParticles&)
  {
    process(collision, tracks);
    processAntonin(collision, tracks);
  }
  PROCESS_SWITCH(taskYann, processMcRec, "Process MC at the reconstruction level", false);


  HistogramRegistry registryMC_Gen{
    "registryMC_Gen",
    {
      {"hZ", "Z position of generated collision; Z (cm); Entries", { HistType::kTH1F, { { 100, -20., 20.} } } },
      // {"hPtMCLbar", "Transverse momentum; p_{T}  [Gev / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      // {"hPtMCL", "Transverse momentum; p_{T}  [Gev / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hpPt_Generated", "; p_{T}  [Gev / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hpBarPt_Generated", "; p_{T}  [Gev / c]; Entries  per  100 MeV/c", { HistType::kTH1F, { { 330, 0., 33. } } } },
      {"hEtaCPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon,  per 200 }MeV/#it{c}", { HistType::kTH1F, { { 165, 0., 33. } } } },
      //{"hEtaCToppBarPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon,  per 200 }MeV/#it{c}", { HistType::kTH1F, { { 165, 0., 33. } } } },
      {"hEtaCAsMotherPt_Generated", ";#it{p}_{T}^{#it{MC Truth}} (GeV/#it{c});#it{A#times#varepsilon,  per 200 }MeV/#it{c}", { HistType::kTH1F, { { 165, 0., 33. } } } },
    }
  };

  void processMcGen(aod::McCollision const& collision, aod::McParticles& particles)
  {

    // Generated collision counter
    if (std::abs(collision.posZ()) > cfgZmax)
      return;
    registryMC_Gen.get<TH1>(HIST("hZ"))->Fill(collision.posZ());


    // Loop directly over all MC particles in this collision
    for (auto const& particle : particles) {
      // pT > 0.5 GeV/c cut
      if (particle.pt() > cfgPtMin) {

        // Directly EtaC particle (pdgCode : 441) w/ rapidity cut |y| < 0.5
        if (particle.pdgCode() == 441 && std::abs(particle.y()) < 0.5) 
          registryMC_Gen.get<TH1>(HIST("hEtaCPt_Generated"))->Fill(particle.pt());

        // EtaC as mother of a particle
        if ((std::abs(particle.eta()) < cfgEtaMax) && particle.has_mothers()) {
          auto const& mother = particle.template mothers_first_as<aod::McParticles>();

          if (mother.pdgCode() == 441 && std::abs(mother.y()) < 0.5 && mother.pt() > cfgPtMin)
            registryMC_Gen.get<TH1>(HIST("hEtaCAsMotherPt_Generated"))->Fill(mother.pt());

          // EtaC as a mother of ppBar pair
          for (auto const& p2 : particles) {
            if ((std::abs(p2.eta()) < cfgEtaMax) && (p2.pt() > cfgPtMin)) {
              int idxEtaC = -1;
              if (isTrueProton(particle, kProton) && particle.has_mothers()) {
                auto const& mother1 = particle.template mothers_first_as<aod::McParticles>();

                if (mother1.pdgCode() == 441 && std::abs(mother1.y()) < 0.5) // TODO? : Add pT cut on etaC mother
                  idxEtaC = mother1.globalIndex();
              }
              if (isTrueProton(p2, kProtonBar) && p2.has_mothers()) {
                auto const& mother2 = p2.template mothers_first_as<aod::McParticles>();

                if (mother2.globalIndex() == idxEtaC && std::abs(mother2.y()) < 0.5) // TODO? : Add pT cut on etaC mother
                  registryMC_Gen.get<TH1>(HIST("hEtaCToppBarPt_Generated"))->Fill(mother2.pt());
              }
            }
          }
        }

      }
     

      if (!particle.isPhysicalPrimary()) // Primary = From collision vertex ?? or from PV ?? Remove ??
        continue;

      // Efficiency part
      if (std::abs(particle.eta()) < cfgEtaMax && particle.pt() > cfgPtMin) { // Produced with the good acceptance and a good enough transverse momentum
        // Efficiency for anti-proton or proton
        // Select : kProtonBar = Anti-proton ; kProton = Proton
        if (particle.pdgCode() == kProton) // Proton/Anti-proton
          registryMC_Gen.get<TH1>(HIST("hpPt_Generated"))->Fill(particle.pt());
        if (particle.pdgCode() == kProtonBar) // Proton/Anti-proton
          registryMC_Gen.get<TH1>(HIST("hpBarPt_Generated"))->Fill(particle.pt());
        
        if (particle.has_mothers()) {
          auto const& mother = particle.template mothers_first_as<aod::McParticles>();

          if (mother.pdgCode() == 441 && std::abs(mother.y()) < 0.5) { // TODO? : Add pT cut on etaC mother
            if (particle.pdgCode() == kProton) { // EtaC from proton
              registryMC_Gen.get<TH1>(HIST("hEtaCAsMotherOfpPt_Generated"))->Fill(mother.pt());
              registryMC_Gen.get<TH1>(HIST("hpWithEtaCMotherPt_Generated"))->Fill(particle.pt());
            }
            if (particle.pdgCode() == kProtonBar) { // EtaC from anti-proton
              registryMC_Gen.get<TH1>(HIST("hEtaCAsMotherOfpBarPt_Generated"))->Fill(mother.pt());
              registryMC_Gen.get<TH1>(HIST("hpBarWithEtaCMotherPt_Generated"))->Fill(particle.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(taskYann, processMcGen, "Process MC at the generator level", false);


  // ES mode = NEW (switch in confMerge.json)
  // PID mode : myMcTracks = OLD ; myPIDMCTracks = NEW
  // DCA mode : myPIDMCTracks = OFF ; myPIDMCTracksDCA = ON
  void processSelectedData(mySelectedEvent const& collision, myPIDMCTracksDCA const& tracks, aod::McParticles&)
  //void processSelectedData(mySelectedEvent const& collision, myPIDMCTracks const& tracks, aod::McParticles&)
  //void processSelectedData(mySelectedEvent const& collision, myMcTracks const& tracks, aod::McParticles&)
  //void processSelectedData(mySelectedEvent const& collision, myTracks const& tracks)
  {
    // Collision counter before ES...
    if (isCollisionAccepted(collision))
      hVtxBeforeES->Fill(collision.posZ());

    if (!collision.sel8())
      return;

    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup))
    //if (!collision.selection_bit(aod::evsel::kNoCollInTimeRangeNarrow))
    //if (!collision.selection_bit(aod::evsel::kNoCollInTimeRangeStrict))
      return;

    process(collision, tracks);
    processAntonin(collision, tracks);
    //LOG(info) << "Sel8() OK";
  }
  PROCESS_SWITCH(taskYann, processSelectedData, "Process data w/ Selected Events", false);




  HistogramRegistry registryBC{
    "registryBC",
    {
      {"hCollDistri", ";#it{N_{Coll}}; Entries", { HistType::kTH1F, { { 20, 0., 10.} } } },
    }
  };


  void processBC(aod::BC const&, aod::Collisions const& collisions)
  {
    int ncol = 0;
    for (auto const& coll : collisions) {
      ncol++;
    }
    registryBC.get<TH1>(HIST("hCollDistri"))->Fill(ncol);
  }
  PROCESS_SWITCH(taskYann, processBC, "Process Bunch Crossings (BC)", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskYann>(cfgc),
  };

}
