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
//
// \file JetD0AngSubstructure.cxx
//
// \brief Analysis task for the reconstruction and study of charged jets
//        containing D_0 mesons in pp collisions.
// \inherited from D0 fragmentation and Ds
// \P. Dhankher

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "TVector3.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Definition of a custom AOD table to store jet–D0 quantities
namespace o2::aod
{
namespace exp
{
// Jet-related quantities
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(JetAng, jetAng, float);
// D0 candidate quantities
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
// ML scores
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);
} // namespace exp
namespace mcp
{
// MCP quantities (Particle Level)
// Jets
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, float);
DECLARE_SOA_COLUMN(JetAng, jetAng, float);
// D0 candidates (Heavy Flavour)
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
// DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfPrompt, hfPrompt, bool);
DECLARE_SOA_COLUMN(HfMatch, hfMatch, bool);

} // namespace mcp
namespace mcd
{
// MCD quantities (Detector Level)
// Jets
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, float);
DECLARE_SOA_COLUMN(JetAng, jetAng, float);
// D0 candidates (Heavy Flavour)
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfPrompt, hfPrompt, bool);
DECLARE_SOA_COLUMN(HfMatch, hfMatch, bool);
// Other
DECLARE_SOA_COLUMN(HfMatchedFrom, hfMatchedFrom, int);
DECLARE_SOA_COLUMN(HfSelectedAs, hfSelectedAs, int);
// ML scores
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);

} // namespace mcd

// AOD table definition
DECLARE_SOA_TABLE(EXPJetObjTable, "AOD", "EXPJETOBJTABLE",
                  exp::JetHfDist,
                  exp::JetPt,
                  exp::JetEta,
                  exp::JetPhi,
                  exp::JetNConst,
                  exp::JetAng,
                  exp::HfPt,
                  exp::HfEta,
                  exp::HfPhi,
                  exp::HfMass,
                  exp::HfY,
                  exp::HfMlScore0,
                  exp::HfMlScore1,
                  exp::HfMlScore2);

DECLARE_SOA_TABLE(MCPJetObjTable, "AOD", "MCPJETOBJTABLE",
                  mcp::JetHfDist,
                  mcp::JetPt,
                  mcp::JetEta,
                  mcp::JetPhi,
                  mcp::JetNConst,
                  mcp::JetAng,
                  mcp::HfPt,
                  mcp::HfEta,
                  mcp::HfPhi,
                  // mcp::HfMass,
                  mcp::HfY,
                  mcp::HfPrompt,
                  mcp::HfMatch);

DECLARE_SOA_TABLE(MCDJetObjTable, "AOD", "MCDJETOBJTABLE",
                  mcd::JetHfDist,
                  mcd::JetPt,
                  mcd::JetEta,
                  mcd::JetPhi,
                  mcd::JetNConst,
                  mcd::JetAng,
                  mcd::HfPt,
                  mcd::HfEta,
                  mcd::HfPhi,
                  mcd::HfMass,
                  mcd::HfY,
                  mcd::HfPrompt,
                  mcd::HfMatch,
                  mcd::HfMlScore0,
                  mcd::HfMlScore1,
                  mcd::HfMlScore2,
                  mcd::HfMatchedFrom,
                  mcd::HfSelectedAs);

DECLARE_SOA_TABLE(MatchJetDistanceTable, "AOD", "MATCHTABLE",
                  mcp::JetHfDist,
                  mcp::JetPt,
                  mcp::JetEta,
                  mcp::JetPhi,
                  mcp::JetNConst,
                  mcp::JetAng,
                  mcp::HfPt,
                  mcp::HfEta,
                  mcp::HfPhi,
                  mcp::HfY,
                  mcp::HfPrompt,
                  mcd::JetHfDist,
                  mcd::JetPt,
                  mcd::JetEta,
                  mcd::JetPhi,
                  mcd::JetNConst,
                  mcd::JetAng,
                  mcd::HfPt,
                  mcd::HfEta,
                  mcd::HfPhi,
                  mcd::HfMass,
                  mcd::HfY,
                  mcd::HfPrompt,
                  mcd::HfMlScore0,
                  mcd::HfMlScore1,
                  mcd::HfMlScore2,
                  mcd::HfMatchedFrom,
                  mcd::HfSelectedAs);

} // namespace o2::aod

// Helps to avoid typos in histogram names when using the "HIST" macro
namespace histnames
{
#define HNAME(name) constexpr const char* name = #name;
HNAME(h_exp_collision_counter);
HNAME(h_exp_d0_jet_counter);
HNAME(h_exp_d0_jet_projection);
HNAME(h_exp_d0_jet_distance_vs_projection);
HNAME(h_exp_d0_jet_distance);
HNAME(h_exp_d0_jet_pt);
HNAME(h_exp_d0_jet_eta);
HNAME(h_exp_d0_jet_phi);
HNAME(h_exp_d0_jet_ang);
HNAME(h_exp_d0_mass);
HNAME(h_exp_d0_eta);
HNAME(h_exp_d0_phi);
HNAME(h_mc_collision_counter);
HNAME(h_mc_d0_jet_counter);
#undef HNAME
} // namespace histnames

struct JetD0AngSubstructure {

  // Output table producer
  Produces<aod::EXPJetObjTable>        ObjJetTable;
  Produces<aod::MCDJetObjTable>        mcdJetTable;
  Produces<aod::MCPJetObjTable>        mcpJetTable;
  Produces<aod::MatchJetDistanceTable> matchJetTable;

  // MC Matching Tables
  using JetD0MCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetD0MCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;
  // Type Aliases for the constituent jets in the table (to be used for angularity calculation)
  using JetD0MCDTableConstituent = JetD0MCDTable::iterator;
  using JetD0MCPTableConstituent = JetD0MCPTable::iterator;
  // Slices for access to proper HF MCD jet collision that is associated to MCCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetD0MCDTable>                 d0MCDJetsPerCollisionPreslice    = aod::jet::collisionId;
  Preslice<JetD0MCPTable>                 d0MCPJetsPerMCCollisionPreslice  = aod::jet::mcCollisionId;

  HistogramRegistry registry{"registry", {}};

  // Configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float>       kappa{"kappa", 1.0, "angularity kappa"}; // to do: configurable from json
  Configurable<float>       alpha{"alpha", 1.0, "angularity alpha"};

  std::vector<int> eventSelectionBits;
  int              trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection     = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    addHistograms();
  }

  /*
   * Histogram registry
   *
   * Contains:
   *  - Event and track histograms
   *  - Jet kinematic distributions
   *  - D0–jet substructure observables
   */
  void addHistograms()
  {
    // Data Histograms: Events
    registry.add(histnames::h_exp_collision_counter, "# of collisions;", {HistType::kTH1F, {{2, 0., 2.}}});
    auto expCollisionCounter = registry.get<TH1>(HIST(histnames::h_exp_collision_counter));
    expCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    expCollisionCounter->GetXaxis()->SetBinLabel(2, "sel8 + zcut");
    // Data Histograms: Jets
    registry.add(histnames::h_exp_d0_jet_counter, ";# of D^{0} jets;", {HistType::kTH1F, {{6, 0., 3.0}}});
    auto expJetCounter = registry.get<TH1>(HIST(histnames::h_exp_d0_jet_counter));
    expJetCounter->GetXaxis()->SetBinLabel(1, "Charged jets with D0");

    registry.add(histnames::h_exp_d0_jet_projection, ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::h_exp_d0_jet_distance_vs_projection, ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add(histnames::h_exp_d0_jet_distance, ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::h_exp_d0_jet_pt, ";p_{T,D^{0} jet};dN/dp_{T,D^{0} jet}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add(histnames::h_exp_d0_jet_eta, ";#eta_{T,D^{0} jet};dN/d#eta_{D^{0} jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::h_exp_d0_jet_phi, ";#phi_{T,D^{0} jet};dN/d#phi_{D^{0} jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add(histnames::h_exp_d0_jet_ang, ";angularity_{D^{0} jet};dN/d(angularity)", {HistType::kTH1F, {{1000, 0., 10.}}});
    // Data Histograms: D0 Candidates
    registry.add(histnames::h_exp_d0_mass, ";m_{D^{0}});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::h_exp_d0_eta, ";#eta_{D^{0}});dN/d#eta_{D_{}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::h_exp_d0_phi, ";#phi_{D^{0}});dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
    // MC Histograms: Events
    registry.add(histnames::h_mc_collision_counter, "# of collisions;", {HistType::kTH1F, {{4, 0., 2.}}});
    auto mcCollisionCounter = registry.get<TH1>(HIST(histnames::h_mc_collision_counter));
    mcCollisionCounter->GetXaxis()->SetBinLabel(1, "mccollisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(2, "z_cut");
    mcCollisionCounter->GetXaxis()->SetBinLabel(3, "collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(4, "sel8");
    // MC Histograms: Jets
    registry.add(histnames::h_mc_d0_jet_counter, ";# of D^{0} jets;", {HistType::kTH1F, {{6, 0., 3.0}}});
    auto jetCounter = registry.get<TH1>(HIST(histnames::h_mc_d0_jet_counter));
    jetCounter->GetXaxis()->SetBinLabel(1, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(2, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(3, "particle matched jets");
    jetCounter->GetXaxis()->SetBinLabel(4, "detector matched jets");
    jetCounter->GetXaxis()->SetBinLabel(5, "mcd matched to mcp loop");
    jetCounter->GetXaxis()->SetBinLabel(6, "mcp matched to mcd loop");
    // MC Histograms: D0 candidates
  };

  template <typename T, typename U>
  float jetCalculateAngularityEXP(T const& jet, U const& /*tracks*/)
  {
    float tAngularity = 0.0;
    for (auto& constituent : jet.template tracks_as<U>()) {
      tAngularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent) / (jet.r() / 100.f), alpha);
    }
    tAngularity /= std::pow(jet.pt(), kappa);
    return tAngularity;
  }

  // template <typename T>
  float jetCalculateAngularityMCD(
    JetD0MCDTableConstituent const& jet,
    aod::JetTracks const&           tracks)
  {
    float a = 0.f;
    for (auto id : jet.tracksIds()) {
      auto trk = tracks.iteratorAt(id);
      a += std::pow(trk.pt(), kappa) *
           std::pow(jetutilities::deltaR(jet, trk) / (jet.r() / 100.f), alpha);
    }
    return a / std::pow(jet.pt(), kappa);
  }

  float jetCalculateAngularityMCP(
    JetD0MCPTableConstituent const& jet,
    aod::JetParticles const&        particles)
  {
    float a = 0.f;
    for (auto id : jet.tracksIds()) {
      auto p = particles.iteratorAt(id);
      a += std::pow(p.pt(), kappa) *
           std::pow(jetutilities::deltaR(jet, p) / (jet.r() / 100.f), alpha);
    }
    return a / std::pow(jet.pt(), kappa);
  }

  void processDataChargedSubstructure(aod::JetCollision const&                                            collision,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      aod::CandidatesD0Data const&, aod::JetTracks const& tracks)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST(histnames::h_exp_collision_counter), 0.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST(histnames::h_exp_collision_counter), 1.5);

    // Loop over jets containing D0 candidates
    for (const auto& jet : jets) {
      // number of charged jets with D0
      registry.fill(HIST(histnames::h_exp_d0_jet_counter), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Loop over D0 candidates associated to the jet
      for (const auto& d0Candidate : jet.candidates_as<aod::CandidatesD0Data>()) {
        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());

        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double zParallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, d0Candidate);

        float angularity = jetCalculateAngularityEXP(jet, tracks);

        // filling histograms
        registry.fill(HIST(histnames::h_exp_d0_jet_projection), zParallel);
        registry.fill(HIST(histnames::h_exp_d0_jet_distance_vs_projection), axisDistance, zParallel);
        registry.fill(HIST(histnames::h_exp_d0_jet_distance), axisDistance);
        registry.fill(HIST(histnames::h_exp_d0_jet_pt), jet.pt());
        registry.fill(HIST(histnames::h_exp_d0_jet_eta), jet.eta());
        registry.fill(HIST(histnames::h_exp_d0_jet_phi), jet.phi());
        registry.fill(HIST(histnames::h_exp_d0_jet_ang), angularity);
        registry.fill(HIST(histnames::h_exp_d0_mass), d0Candidate.m());
        registry.fill(HIST(histnames::h_exp_d0_eta), d0Candidate.eta());
        registry.fill(HIST(histnames::h_exp_d0_phi), d0Candidate.phi()); // add more axis

        // filling table
        ObjJetTable(axisDistance,
                    jet.pt(),
                    jet.eta(),
                    jet.phi(),
                    jet.tracks_as<aod::JetTracks>().size(), angularity,
                    d0Candidate.pt(),
                    d0Candidate.eta(),
                    d0Candidate.phi(),
                    d0Candidate.m(),
                    d0Candidate.y(),
                    d0Candidate.mlScores()[0],
                    d0Candidate.mlScores()[1],
                    d0Candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found in jet
      } // end of D0 candidates loop

    } // end of jets loop

  } // end of process function

  PROCESS_SWITCH(JetD0AngSubstructure, processDataChargedSubstructure, "charged HF jet substructure", false);

  // void processMonteCarloEfficiency(aod::JetMcCollisions const&  mccollisions,
  //                                  aod::JetCollisionsMCD const& collisions,
  //                                  JetD0MCDTable const&         mcdjets,
  //                                  JetD0MCPTable const&         mcpjets,
  //                                  aod::CandidatesD0MCD const& /*mcdCandidates*/,
  //                                  aod::CandidatesD0MCP const& /*mcpCandidates*/,
  //                                  aod::JetTracks const& tracks,
  //                                  aod::JetParticles const& /*particles*/)
  // {
  //   for (const auto& mccollision : mccollisions) {

  //     registry.fill(HIST(histnames::h_mc_collision_counter), 0.0);
  //     // skip collisions outside of |z| < vertexZCut
  //     if (std::abs(mccollision.posZ()) > vertexZCut) {
  //       continue;
  //     }
  //     registry.fill(HIST(histnames::h_mc_collision_counter), 1.0);

  //     // reconstructed collisions associated to same mccollision
  //     const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
  //     for (const auto& collision : collisionsPerMCCollision) {

  //       registry.fill(HIST(histnames::h_mc_collision_counter), 2.0);
  //       if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
  //         continue;
  //       }
  //       registry.fill(HIST(histnames::h_mc_collision_counter), 3.0);

  //       // d0 detector level jets associated to the current same collision
  //       const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
  //       for (const auto& mcdjet : d0mcdJetsPerCollision) {

  //         registry.fill(HIST(histnames::h_mc_d0_jet_counter), 0.5);

  //         // obtain leading HF candidate in jet
  //         auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

  //         if (mcdjet.has_matchedJetCand()) {
  //           registry.fill(HIST(histnames::h_mc_d0_jet_counter), 1.5);
  //         }

  //         // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
  //         int matchedFrom  = 0;
  //         int decayChannel = o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
  //         int selectedAs   = 0;

  //         if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
  //           matchedFrom = 1;
  //         } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
  //           matchedFrom = -1;
  //         }
  //         // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
  //         if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
  //           selectedAs = 1;
  //         } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
  //           selectedAs = -1;
  //         }

  //         float angularity = jetCalculateAngularity(mcdjet, tracks);

  //         mcdJetTable(jetutilities::deltaR(mcdjet, mcdd0cand),
  //                     mcdjet.pt(),
  //                     mcdjet.eta(),
  //                     mcdjet.phi(),
  //                     mcdjet.tracks_as<aod::JetTracks>().size(), // detector level jet
  //                     angularity,
  //                     mcdd0cand.pt(),
  //                     mcdd0cand.eta(),
  //                     mcdd0cand.phi(),
  //                     mcdd0cand.m(),
  //                     mcdd0cand.y(),
  //                     (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0 candidate
  //                     mcdjet.has_matchedJetCand(),
  //                     mcdd0cand.mlScores()[0],
  //                     mcdd0cand.mlScores()[1],
  //                     mcdd0cand.mlScores()[2],  // // Machine Learning PID scores: background, prompt, non-prompt
  //                     matchedFrom, selectedAs); // D0 = +1, D0bar = -1, neither = 0
  //       }
  //     } // end of reconstructed collisions loop (detector level mc collisions matching with real collisions)

  //     // d0 particle level jets associated to same mccollision
  //     const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
  //     for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

  //       registry.fill(HIST(histnames::h_mc_d0_jet_counter), 0.0);

  //       // obtain leading HF particle in jet
  //       auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

  //       if (mcpjet.has_matchedJetCand()) {
  //         registry.fill(HIST(histnames::h_mc_d0_jet_counter), 1.0);
  //       }
  //       float angularity = jetCalculateAngularity(mcpjet, tracks);
  //       // store data in MC detector level table (calculate angular distance in eta-phi plane on the fly)
  //       mcpJetTable(jetutilities::deltaR(mcpjet, mcpd0cand),
  //                   mcpjet.pt(),
  //                   mcpjet.eta(),
  //                   mcpjet.phi(),
  //                   mcpjet.tracks_as<aod::JetParticles>().size(), // particle level jet
  //                   angularity,
  //                   mcpd0cand.pt(),
  //                   mcpd0cand.eta(),
  //                   mcpd0cand.phi(),
  //                   // mcpd0cand.m(),
  //                   mcpd0cand.y(),
  //                   (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level D0
  //                   mcpjet.has_matchedJetCand());
  //     } // End of particle level jets loop, related to detector level collisions and jets.
  //   }
  // }

  // PROCESS_SWITCH(JetD0AngSubstructure, processMonteCarloEfficiency, "non-matched and matched MC HF and jets", false);

  // template <typename TMCPJetsPerMCCollisionPreslice, typename TJetsMCD, typename TJetsMCP, typename TCandidatesMCD, typename TCandidatesMCP>
  void processMonteCarlo(aod::JetMcCollisions const&  mccollisions,
                         aod::JetCollisionsMCD const& collisions,
                         JetD0MCDTable const& /*mcdjets*/,
                         JetD0MCPTable const& mcpjets,
                         aod::CandidatesD0MCD const& /*mcdCandidates*/,
                         aod::CandidatesD0MCP const& /*mcpCandidates*/,
                         aod::JetTracks const&    jettracks,
                         aod::JetParticles const& jetparticles)
  {
    for (const auto& mccollision : mccollisions) {
      registry.fill(HIST(histnames::h_mc_collision_counter), 0.0);
      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST(histnames::h_mc_collision_counter), 1.0);

      // hf particle level jets associated to same mccollision
      const JetD0MCPTable mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : mcpJetsPerMCCollision) {

        registry.fill(HIST(histnames::h_mc_d0_jet_counter), 0.0);

        // obtain leading HF particle in jet
        auto mcpcand = mcpjet.template candidates_first_as<aod::CandidatesD0MCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST(histnames::h_mc_d0_jet_counter), 1.0);

          // loop over detector level matched to current particle level
          for (const auto& mcdjet : mcpjet.template matchedJetCand_as<JetD0MCDTable>()) {
            registry.fill(HIST(histnames::h_mc_d0_jet_counter), 2.0);

            // apply collision sel8 selection on detector level jet's collision
            const auto& collision = collisions.iteratorAt(mcdjet.collisionId());
            registry.fill(HIST(histnames::h_mc_collision_counter), 2.0);
            if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
              continue;
            }
            registry.fill(HIST(histnames::h_mc_collision_counter), 3.0);

            // obtain leading HF candidate in jet
            auto mcdcand = mcdjet.template candidates_first_as<aod::CandidatesD0MCD>();

            // reflection information for storage: HF = +1, HFbar = -1, neither = 0
            int matchedFrom  = 0;
            int decayChannel = 0;
            if (jethfutilities::isD0Table<aod::CandidatesD0MCD>()) {
              decayChannel = o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
            } else if (jethfutilities::isLcTable<aod::CandidatesD0MCD>()) {
              decayChannel = o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi;
            }
            int selectedAs = 0;

            if (mcdcand.flagMcMatchRec() == decayChannel) { // matched to HF on truth level
              matchedFrom = 1;
            } else if (mcdcand.flagMcMatchRec() == -decayChannel) { // matched to HFbar on truth level
              matchedFrom = -1;
            }
            // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
            if (mcdcand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as HF
              selectedAs = 1;
            } else if (mcdcand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as HFbar
              selectedAs = -1;
            }

            float mcpAngularity = jetCalculateAngularityMCP(mcpjet, jetparticles);
            float mcdAngularity = jetCalculateAngularityMCD(mcdjet, jettracks);
            // float mcpAngularity = 0.;
            // float mcdAngularity = 0.;

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpcand),
                          mcpjet.pt(),
                          mcpjet.eta(),
                          mcpjet.phi(),
                          mcpjet.template tracks_as<aod::JetParticles>().size(), // particle level jet
                          mcpAngularity,
                          mcpcand.pt(),
                          mcpcand.eta(),
                          mcpcand.phi(),
                          mcpcand.y(),
                          (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level HF
                          jetutilities::deltaR(mcdjet, mcdcand),
                          mcdjet.pt(),
                          mcdjet.eta(),
                          mcdjet.phi(),
                          mcdjet.template tracks_as<aod::JetTracks>().size(), // detector level jet
                          mcdAngularity,
                          mcdcand.pt(),
                          mcdcand.eta(),
                          mcdcand.phi(),
                          mcdcand.m(),
                          mcdcand.y(),
                          (mcdcand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level HF
                          mcdcand.mlScores()[0],
                          mcdcand.mlScores()[1],
                          mcdcand.mlScores()[2], // Machine Learning PID scores: background, prompt, non-prompt
                          matchedFrom,
                          selectedAs); // HF = +1, HFbar = -1, neither = 0
          }
        } else {
          // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
          float mcpAngularity = jetCalculateAngularityMCP(mcpjet, jetparticles);
          // float mcpAngularity = 0.;
          matchJetTable(jetutilities::deltaR(mcpjet, mcpcand),
                        mcpjet.pt(),
                        mcpjet.eta(),
                        mcpjet.phi(),
                        mcpjet.template tracks_as<aod::JetParticles>().size(), // particle level jet
                        mcpAngularity,

                        mcpcand.pt(),
                        mcpcand.eta(),
                        mcpcand.phi(),
                        mcpcand.y(),
                        (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level HF
                        -2, -2, -2, -2, -2,                                       // detector level jet
                        -2, -2, -2, -2, -2, -2, -2,                               // detector level HF
                        -2, -2, -2,                                               // Machine Learning PID scores: background, prompt, non-prompt
                        -2, -2);                                                  // HF = +1, HFbar = -1, neither = 0
        }
      } // end of mcpjets loop
    } // end of mccollisions loop
  };

  PROCESS_SWITCH(JetD0AngSubstructure, processMonteCarlo, "Store all simulated D0 jets information with matched candidate (if any found)", false);
};
// Workflow definition
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetD0AngSubstructure>(
      cfgc,
      TaskName{"jet-d0-ang-substructure"})};
}
