#include <array>
#include <random>
#include <vector>

#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"

int main() {
  Particle::addParticleType("pion+", 0.13957, 1);
  Particle::addParticleType("pion-", 0.13957, -1);
  Particle::addParticleType("kaon+", 0.49367, 1);
  Particle::addParticleType("kaon-", 0.49367, -1);
  Particle::addParticleType("proton+", 0.93827, 1);
  Particle::addParticleType("proton-", 0.93827, -1);
  Particle::addParticleType("k*", 0.89166, 0, 0.050);

  std::random_device rd;
  gRandom->SetSeed(rd());

  std::array<TH1*, 12> histo_array;

  // particle histograms
  TH1I* particle_types_histogram =
      new TH1I("particle_types_histogram", "Particle types", 7, 0, 7);
  histo_array[0] = particle_types_histogram;

  TH1F* azimutal_angles_histogram = new TH1F(
      "azimutal_angles_histogram", "Azimutal angles", 1e5, 0, TMath::Pi());
  histo_array[1] = azimutal_angles_histogram;

  TH1F* polar_angles_histogram = new TH1F(
      "polar_angles_histogram", "Polar angles", 1e5, 0, TMath::Pi() * 2.);
  histo_array[2] = polar_angles_histogram;

  TH1F* momentum_histogram =
      new TH1F("momentum_histogram", "Momentum", 1e4, 0, 10);
  histo_array[3] = momentum_histogram;

  TH1F* momentum_xy_istogram =
      new TH1F("momentum_xy_histogram", "Momentum xy", 1e4, 0, 10);
  histo_array[4] = momentum_xy_istogram;

  TH1F* energy_histogram = new TH1F("energy_histogram", "Energy", 1e5, 0, 4);
  histo_array[5] = energy_histogram;

  // invariant mass histograms
  TH1F* invm_all_h =
      new TH1F("invm_all_h", "Invariant mass, all particles", 1e5, 0, 1000);
  invm_all_h->Sumw2();
  histo_array[6] = invm_all_h;

  TH1F* invm_opposite_charge_h =
      new TH1F("invm_opposite_charge_h", "Invariant mass, opposite charge", 1e5,
               0, 1000);
  invm_opposite_charge_h->Sumw2();
  histo_array[7] = invm_opposite_charge_h;

  TH1F* invm_same_charge_h = new TH1F(
      "invm_same_charge_h", "Invariant mass, same charge", 1e5, 0, 1000);
  invm_same_charge_h->Sumw2();
  histo_array[8] = invm_same_charge_h;

  TH1F* invm_pion_kaon_opposite_h = new TH1F(
      "invm_pionp_kaonn_h",
      "Invariant mass, pion+ and kaon- or pion- and kaon+", 1e5, 0, 1000);
  invm_pion_kaon_opposite_h->Sumw2();
  histo_array[9] = invm_pion_kaon_opposite_h;

  TH1F* invm_pion_kaon_same_h = new TH1F(
      "invm_pion_kaon_same_h",
      "Invariant mass, pion+ and kaon+ or pion- and kaon-", 1e5, 0, 1000);
  invm_pion_kaon_same_h->Sumw2();
  histo_array[10] = invm_pion_kaon_same_h;

  TH1F* invm_decayed_h =
      new TH1F("invm_decayed_h", "Invariant mass, decayed particles from K*",
               1e5, 0, 1000);
  invm_decayed_h->Sumw2();
  histo_array[11] = invm_decayed_h;

  TH1F* check_phi =
      new TH1F("check_phi", "Check phi", 1e5, 0, TMath::Pi() * 2.);

  for (int i{}; i < 1e5; ++i) {
    std::vector<Particle> event_particles(100);

    event_particles.reserve(150);

    for (int j{}; j < 1e2; ++j) {
      auto phi = gRandom->Uniform(0, TMath::Pi() * 2.);
      auto theta = gRandom->Uniform(0, TMath::Pi());
      auto momentum = gRandom->Exp(1);  // GeV

      check_phi->Fill(phi);

      // convert polar to cartesian coordinates
      event_particles[j].setMomentum(
          Momentum{PolarVector{momentum, theta, phi}});

      auto x = gRandom->Rndm();
      if (x <= 0.4) {
        event_particles[j].setIndex("pion+");
      } else if (x <= 0.8) {
        event_particles[j].setIndex("pion-");
      } else if (x <= 0.85) {
        event_particles[j].setIndex("kaon+");
      } else if (x <= 0.9) {
        event_particles[j].setIndex("kaon-");
      } else if (x <= 0.945) {
        event_particles[j].setIndex("proton+");
      } else if (x <= 0.99) {
        event_particles[j].setIndex("proton-");
      } else {
        event_particles[j].setIndex("k*");

        auto decay_into = gRandom->Rndm();

        if (decay_into <= 0.5) {
          Particle pion_positive{"pion+"};
          Particle kaon_negative{"kaon-"};

          event_particles[j].decayToBody(pion_positive, kaon_negative);

          event_particles.push_back(pion_positive);
          event_particles.push_back(kaon_negative);
        } else {
          Particle pion_negative{"pion-"};
          Particle kaon_positive{"kaon+"};

          event_particles[j].decayToBody(pion_negative, kaon_positive);

          event_particles.push_back(pion_negative);
          event_particles.push_back(kaon_positive);
        }
      }
    }

    auto event_particles_begin = event_particles.begin();
    auto event_particles_end = event_particles.end();

    for (auto it = event_particles_begin; it < event_particles_end; ++it) {
      if (it < event_particles_begin + 100) {
        // type
        particle_types_histogram->Fill(it->getIndex().value());

        auto momentum = it->getMomentum();
        auto polar_momentum = momentum.getPolar();

        // azimutal angle
        azimutal_angles_histogram->Fill(polar_momentum.theta);

        // polar angle
        polar_angles_histogram->Fill(polar_momentum.phi);

        // momentum
        momentum_histogram->Fill(std::sqrt(momentum * momentum));

        // momentum xy
        momentum_xy_istogram->Fill(
            std::sqrt(momentum.x * momentum.x + momentum.y + momentum.y));

        // energy
        energy_histogram->Fill(it->getEnergy());
      }

      // invariant mass
      if (it->getName() != "k*") {
        for (auto invm_it = it + 1; invm_it < event_particles_end; ++invm_it) {
          if (invm_it->getName() == "k*") {
            continue;
          }

          auto invariant_mass = it->getInvariantMass(*invm_it);
          // invariant mass with all particles
          invm_all_h->Fill(invariant_mass);

          // invariant mass with opposite charge particles
          if (invm_it->getCharge() * it->getCharge() < 0) {
            invm_opposite_charge_h->Fill(invariant_mass);
          }

          // invariant mass with same charge particles
          if (invm_it->getCharge() * it->getCharge() > 0) {
            invm_same_charge_h->Fill(invariant_mass);
          }

          // invariant mass with pion+ and kaon- or pion- and kaon+
          if ((it->getName() == "pion+" && invm_it->getName() == "kaon-") ||
              (it->getName() == "kaon-" && invm_it->getName() == "pion+") ||
              (it->getName() == "pion-" && invm_it->getName() == "kaon+") ||
              (it->getName() == "kaon+" && invm_it->getName() == "pion-")) {
            invm_pion_kaon_opposite_h->Fill(invariant_mass);
          }

          // invariant mass with pion+ and kaon+ or piaon- and kaon-
          if ((it->getName() == "pion+" && invm_it->getName() == "kaon+") ||
              (it->getName() == "kaon+" && invm_it->getName() == "pion+") ||
              (it->getName() == "pion-" && invm_it->getName() == "kaon-") ||
              (it->getName() == "kaon-" && invm_it->getName() == "pion-")) {
            invm_pion_kaon_opposite_h->Fill(invariant_mass);
          }

          // invariant mass with decayed particles
          if (it > event_particles_begin + 99 && invm_it == it + 1) {
            invm_decayed_h->Fill(invariant_mass);
          }
        }
      }
    }
  }

  TFile* file = new TFile("file.root", "RECREATE");
  for (auto const& p : histo_array) {
    p->Write();
  }
  check_phi->Write();
  file->Close();
}