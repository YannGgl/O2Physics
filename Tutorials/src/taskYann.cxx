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

#include "Common/Core/PID/PIDTOF.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CommonConstants/PhysicsConstants.h>

#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;

using myTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
using myTrack = myTracks::iterator;

using myMcTracks = soa::Join<myTracks, o2::aod::McTrackLabels>;

struct taskYann {

  Configurable<float> cfgZmax{"cfgZmax", 10., "Restriction on the PV position |PVz|<zMax (cm)"};

  OutputObj<TH1F> hVtx{
    TH1F("hVtx", "Primary vertex position after selection; Z (cm)", 100, -20., 20.)};

  OutputObj<TH2F> hTpc{
    TH2F("hTpc", "TPC; Momentum (GeV); dE/dx", 400, -4., 4., 200, 0., 1000)};
  OutputObj<TH2F> hTof{
    TH2F("hTof", "TOF; Momentum (GeV); beta", 500, -5., 5., 220, 0., 1.1)};
  
  HistogramRegistry registry{
    "registry",
    {
      {"hMinvTOF", "Invariant mass; m (GeV); Counts", { HistType::kTH1F, { { 200, 1., 1.5 } } } },
      {"hMinvTOF_PID", "Invariant mass; m (GeV); Counts", { HistType::kTH1F, { { 200, 1., 1.5 } } } },
      {"hMinvTPC", "Invariant mass; m (GeV); Counts", { HistType::kTH1F, { { 200, 1., 1.5 } } } },
      {"hMinvTPC_PID", "Invariant mass; m (GeV); Counts", { HistType::kTH1F, { { 200, 1., 1.5 } } } }
    }
  };
  /*
  OutputObj<TH1F> hMinv{
    TH1F("hMinv", "Invariant mass; m (GeV); Counts", 400, 0., 2.5)};
  OutputObj<TH1F> hMinvPID{
    TH1F("hMinvPID", "Invariant mass; m (GeV); Counts", 400, 0., 2.5)};
  */

  void init(o2::framework::InitContext& ic)
  {
  }

  // Collision selector
  // True when collision(=PV) within Zmax(=10cm) distance of geometrical center
  // False otherwise
  bool isCollisionAccepted(aod::Collision const& collision)
  {
    auto z = collision.posZ();
    if (abs(z) > cfgZmax)
      return false;

    // Some other selections

    return true;
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
    if (!track.hasITS()) // PID detectors
      return false;

    if (!track.hasTPC()) // TPC detectors
      return false;

    if (abs(track.tgl()) > 0.9) // eta ?
      return false;

    if (track.itsNCls() < 7) // nb of clusters before TOF ?
      return false;

    // Some other selections

    return true;
  }


  // Function for selecting tracks under straight lines f(x) = ax + b for TOF
  template <typename TrackType>
  bool isTrackProtonTOF(TrackType const& track)
  {
    double a[4] = {0.2, 0.01, 1., 0.066};
    double b[4] = {0.6, 0.95, 0., 0.79};
    auto mom = track.tpcInnerParam();
    auto beta = tofBeta(track);
    if (beta <= a[0]*mom + b[0]
    && beta <= a[1]*mom + b[1]
    && beta <= a[2]*mom + b[2]
    && beta <= a[3]*mom + b[3]
    )
      return true;

    return false;
  }

  // Function for selecting tracks under straight lines f(x) = ax + b for TPC
  template <typename TrackType>
  bool isTrackProtonTPC(TrackType const& track)
  {
    double a[4] = {-5000, -208, -7.5, -108};
    double b[4] = {2000, 250, 70, 180};
    auto dedx = track.tpcSignal();
    auto mom = track.tpcInnerParam();
    if (dedx >= a[0]*mom + b[0]
    && dedx >= a[1]*mom + b[1]
    && dedx >= a[2]*mom + b[2]
    && dedx >= a[3]*mom + b[3]
    )
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


  // Main method of the task
  // Subscribing to the table Collision, iterating over all tracks of a given collision
  template <typename TrackTypeContainer>
  void process(aod::Collision const& collision, TrackTypeContainer const& tracks)
  {

    // "Loading bar" (print smthg each 1000 collisions)
    static int ncol = 0;
    if (ncol % 1000 == 0)
      LOG(info) << "Collision: " << ncol;
    ncol++;


    if (!isCollisionAccepted(collision))
      return;

    // Collision counter...
    hVtx->Fill(collision.posZ());


    // Loop over all track in a given collision
    for (auto& track : tracks) {
      if (!isTrackAccepted(track))
        continue;

      // Choose particle (sign > 0) or antiparticle (sign < 0) for 1st track
      auto sign = track.sign();
      if (sign > 0) 
        continue;

      // TPC part
      auto mom = track.tpcInnerParam();
      auto dedx = track.tpcSignal();
      for (auto track2 = track + 1 ; track2 != tracks.end() ; ++track2) {
            if (isTrackAccepted(track2) && track2.sign() != sign)
              registry.get<TH1>(HIST("hMinvTPC"))->Fill(sqrt(InvariantMass(track, track2, 0.94, 0.14))); // Invariant mass histogram w/out PID
      }
      if (isTrackProtonTPC(track)) { // PID
        hTpc->Fill(sign * mom, dedx); 

        for (auto track2 = track + 1 ; track2 != tracks.end() ; ++track2) {
          if (track2.sign() == sign)
            continue;
          if (isTrackAccepted(track2))
            registry.get<TH1>(HIST("hMinvTPC_PID"))->Fill(sqrt(InvariantMass(track, track2, 0.94, 0.14))); // Invariant mass histogram w/ PID
        }
      }

      // TOF part
      if (track.hasTOF()) {
        for (auto track2 = track + 1 ; track2 != tracks.end() ; ++track2) {
            if (isTrackAccepted(track2) && track2.sign() != sign)
              registry.get<TH1>(HIST("hMinvTOF"))->Fill(sqrt(InvariantMass(track, track2, 0.94, 0.14))); // Invariant mass histogram w/out PID
        }
        if (isTrackProtonTOF(track)) { // PID
          auto beta = tofBeta(track);
          hTof->Fill(sign * mom, beta); // TOF histogram w/ PID

          for (auto track2 = track + 1 ; track2 != tracks.end() ; ++track2) {
            if (track2.sign() == sign)
              continue;
            if (isTrackAccepted(track2))
              registry.get<TH1>(HIST("hMinvTOF_PID"))->Fill(sqrt(InvariantMass(track, track2, 0.94, 0.14))); // Invariant mass histogram w/ PID
          }
        }
      }
    }
  }

  void processData(aod::Collision const& collision, myTracks const& tracks)
  {
    process(collision, tracks);
  }
  PROCESS_SWITCH(taskYann, processData, "Process data", true);

  void processMcRec(aod::Collision const& collision, myMcTracks const& tracks, aod::McParticles&)
  {
    process(collision, tracks);
  }
  PROCESS_SWITCH(taskYann, processMcRec, "Process MC at the reconstruction level", false);

  void processMcGen(aod::McParticles& particles)
  {
  }
  PROCESS_SWITCH(taskYann, processMcGen, "Process MC at the generator level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskYann>(cfgc),
  };

}
