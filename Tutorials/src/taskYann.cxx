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
      {"hMinv_Total", "Invariant mass; m (GeV); Entries", { HistType::kTH1F, { { 500, 1., 1.5 } } } },
      {"hMinv_PID", "Invariant mass; m (GeV); Entries", { HistType::kTH1F, { { 500, 1., 1.5 } } } }
    }
  };

  HistogramRegistry registryMC{
    "registryMC",
    {
      {"hTofMC", "TOF; Momentum (GeV); beta", { HistType::kTH2F, { { 500, -5., 5.}, { 220, 0., 1.1 } } } },
      {"hTpcMC", "TPC; Momentum (GeV); dE/dx", { HistType::kTH2F, { { 400, -4., 4.}, { 200, 0., 1000 } } } },
      {"hMinvMC", "Invariant mass; m (GeV); Entries", { HistType::kTH1F, { { 500, 1., 1.5 } } } },
      {"hMinvMC_Matched", "Invariant mass; m (GeV); Entries", { HistType::kTH1F, { { 500, 1., 1.5 } } }}
    }
  };
  

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
    if (!track.hasITS()) // ITS (tracking?) detectors
      return false;

    if (!track.hasTPC()) // TPC detectors
      return false;

    if (abs(track.tgl()) > 0.9) // eta ?
      return false;

    if (track.itsNCls() < 7) // nb of clusters of ITS ?
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

  
  // Function to select MC Truth for antiproton
  template <typename McParticleType>
  bool isTrueProton(McParticleType const& mcParticle)
  {
    if (mcParticle.pdgCode() == kProtonBar)
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

    // "Loading bar" (print smthg each X collisions)
    static int ncol = 0;
    if (ncol % 10000 == 0)
      LOG(info) << "Collision: " << ncol;
    ncol++;


    if (!isCollisionAccepted(collision))
      return;

    // Collision counter...
    hVtx->Fill(collision.posZ());


    // First (outer) loop over all tracks in a given collision
    for (auto const& track1 : tracks) {
      
      // Track quality + sign condition to select particles/antiparticles
      // sign > 0 = Antiparticles ; sign < 0 = Particles
      auto sign = track1.sign();
      if (!isTrackAccepted(track1) || sign > 0)
        continue;

      auto mom = track1.tpcInnerParam();
      auto dedx = track1.tpcSignal();
      auto beta = tofBeta(track1);
      double m = 0;

      // Second (inner) loop over all tracks in a given collision
      for (auto const& track2 : tracks) {
        if (!isTrackAccepted(track2) || track2.sign() == sign) 
          continue;

        m = sqrt(InvariantMass(track1, track2, 0.94, 0.14));

        // Without PID
        registry.get<TH1>(HIST("hMinv_Total"))->Fill(m);

        // MC Truth block
        if constexpr (requires {track1.has_mcParticle();}) {
          int idx=-1;
          if (track1.has_mcParticle()) {
            auto const& particle = track1.mcParticle();
            if (isTrueProton(particle)) {
              registryMC.get<TH1>(HIST("hMinvMC"))->Fill(m);
              registryMC.get<TH2>(HIST("hTpcMC"))->Fill(sign * mom, dedx);
              if (particle.has_mothers()) {
                auto const& mother = particle.template mothers_first_as<aod::McParticles>();
                if (mother.pdgCode() == kLambda0Bar)
                  idx = mother.globalIndex();
              }
            }
            if (track1.hasTOF())
              registryMC.get<TH2>(HIST("hTofMC"))->Fill(sign * mom, beta);
          }
          if (track2.has_mcParticle()) {
            auto const& particle = track2.mcParticle();
            if (particle.has_mothers()) {
              auto const& mother = particle.template mothers_first_as<aod::McParticles>();
              if (mother.globalIndex() == idx)
                registryMC.get<TH1>(HIST("hMinvMC_Matched"))->Fill(m);
            }
          }
        }

        // With PID
        if (!track1.hasTOF() || mom<0.5) { // This block is for PID with TPC
          if (!isTrackProtonTPC(track1) || mom>1.5)
            continue;
          registry.get<TH1>(HIST("hMinv_PID"))->Fill(m);
          hTpc->Fill(sign * mom, dedx);
        }
        else { // This block is for PID with TOF
          if (!isTrackProtonTOF(track1))
            continue;
          registry.get<TH1>(HIST("hMinv_PID"))->Fill(m);
          hTof->Fill(sign * mom, beta);
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
