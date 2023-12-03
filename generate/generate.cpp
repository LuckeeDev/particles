#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "TBenchmark.h"
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom.h"

void generate(int n_gen, const char* file_name) {
  gBenchmark->Start("Benchmark");

  R__LOAD_LIBRARY(ParticleType_cpp.so)
  R__LOAD_LIBRARY(ResonanceType_cpp.so)
  R__LOAD_LIBRARY(Particle_cpp.so)

  Particle::addParticleType("pion+", 0.13957, 1);
  Particle::addParticleType("pion-", 0.13957, -1);
  Particle::addParticleType("kaon+", 0.49367, 1);
  Particle::addParticleType("kaon-", 0.49367, -1);
  Particle::addParticleType("proton+", 0.93827, 1);
  Particle::addParticleType("proton-", 0.93827, -1);
  Particle::addParticleType("k*", 0.89166, 0, 0.050);

  gRandom->SetSeed();

  TList* histo_list = new TList();

  // particle histograms
  TH1I* particle_types_histogram =
      new TH1I("particle_types_histogram", "Particle types", 7, 0, 7);
  histo_list->Add(particle_types_histogram);  // 0

  TH1F* azimutal_angles_histogram = new TH1F(
      "azimutal_angles_histogram", "Azimutal angles", 1e3, 0, TMath::Pi());
  histo_list->Add(azimutal_angles_histogram);  // 1

  TH1F* polar_angles_histogram = new TH1F(
      "polar_angles_histogram", "Polar angles", 1e3, 0, TMath::Pi() * 2.);
  histo_list->Add(polar_angles_histogram);  // 2

  TH1F* momentum_histogram =
      new TH1F("momentum_histogram", "Momentum", 1e3, 0, 9);
  histo_list->Add(momentum_histogram);  // 3

  TH1F* momentum_xy_histogram =
      new TH1F("momentum_xy_histogram", "Momentum xy", 1e3, 0, 9);
  histo_list->Add(momentum_xy_histogram);  // 4

  TH1F* energy_histogram = new TH1F("energy_histogram", "Energy", 1e4, 0, 4);
  histo_list->Add(energy_histogram);  // 5

  // invariant mass histograms
  TH1F* invm_all_h =
      new TH1F("invm_all_h", "Invariant mass, all particles", 1e4, 0, 7);
  invm_all_h->Sumw2();
  histo_list->Add(invm_all_h);  // 6

  TH1F* invm_opposite_charge_h = new TH1F(
      "invm_opposite_charge_h", "Invariant mass, opposite charge", 1e4, 0, 7);
  invm_opposite_charge_h->Sumw2();
  histo_list->Add(invm_opposite_charge_h);  // 7

  TH1F* invm_same_charge_h =
      new TH1F("invm_same_charge_h", "Invariant mass, same charge", 1e4, 0, 7);
  invm_same_charge_h->Sumw2();
  histo_list->Add(invm_same_charge_h);  // 8

  TH1F* invm_pion_kaon_opposite_h =
      new TH1F("invm_pion_kaon_opposite_h",
               "Invariant mass, pion+ and kaon- or pion- and kaon+", 1e4, 0, 7);
  invm_pion_kaon_opposite_h->Sumw2();
  histo_list->Add(invm_pion_kaon_opposite_h);  // 9

  TH1F* invm_pion_kaon_same_h =
      new TH1F("invm_pion_kaon_same_h",
               "Invariant mass, pion+ and kaon+ or pion- and kaon-", 1e4, 0, 7);
  invm_pion_kaon_same_h->Sumw2();
  histo_list->Add(invm_pion_kaon_same_h);  // 10

  TH1F* invm_decayed_h =
      new TH1F("invm_decayed_h", "Invariant mass, decayed particles from K*",
               1e3, 0.6, 1.2);
  invm_decayed_h->Sumw2();
  histo_list->Add(invm_decayed_h);  // 11

  for (int i{}; i < n_gen; ++i) {
    std::vector<Particle> event_particles(100);

    event_particles.reserve(150);

    for (int j{}; j < 1e2; ++j) {
      auto phi = gRandom->Uniform(0, TMath::Pi() * 2.);
      auto theta = gRandom->Uniform(0, TMath::Pi());
      auto rho = gRandom->Exp(1);  // GeV

      // convert polar to cartesian coordinates
      event_particles[j].setMomentum(Momentum{PolarVector{rho, theta, phi}});

      auto x = gRandom->Uniform(0, 1);

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

        auto decay_into = gRandom->Uniform(0, 1);

        Particle decay_product_1{};
        Particle decay_product_2{};

        // TODO: fill invariant mass histograms with data from decayed particles
        if (decay_into <= 0.5) {
          decay_product_1.setIndex("pion+");
          decay_product_2.setIndex("kaon-");
        } else {
          decay_product_1.setIndex("pion-");
          decay_product_2.setIndex("kaon+");
        }

        event_particles[j].decayToBody(decay_product_1, decay_product_2);

        auto invariant_mass_products =
            decay_product_1.getInvariantMass(decay_product_2);

        invm_decayed_h->Fill(invariant_mass_products);

        event_particles.push_back(decay_product_1);
        event_particles.push_back(decay_product_2);
      }

      auto const& new_particle = event_particles[j];

      // type
      particle_types_histogram->Fill(new_particle.getIndex().value());

      auto momentum = new_particle.getMomentum();
      auto polar_momentum = momentum.getPolar();

      // azimutal angle
      azimutal_angles_histogram->Fill(polar_momentum.theta);

      // polar angle
      polar_angles_histogram->Fill(polar_momentum.phi);

      // momentum
      momentum_histogram->Fill(std::sqrt(momentum * momentum));

      // momentum xy
      momentum_xy_histogram->Fill(
          std::sqrt(momentum.x * momentum.x + momentum.y + momentum.y));

      // energy
      energy_histogram->Fill(new_particle.getEnergy());

      // invariant mass
      if (new_particle.getName() != "k*") {
        for (auto invm_i = j - 1; invm_i >= 0; --invm_i) {
          auto const& invm_particle = event_particles[invm_i];

          if (invm_particle.getName() == "k*") {
            continue;
          }

          auto invariant_mass = new_particle.getInvariantMass(invm_particle);

          // invariant mass with all particles
          invm_all_h->Fill(invariant_mass);

          // invariant mass with opposite charge particles
          if (invm_particle.getCharge() * new_particle.getCharge() < 0) {
            invm_opposite_charge_h->Fill(invariant_mass);
          }

          // invariant mass with same charge particles
          if (invm_particle.getCharge() * new_particle.getCharge() > 0) {
            invm_same_charge_h->Fill(invariant_mass);
          }

          // invariant mass with pion+ and kaon- or pion- and kaon+
          if ((new_particle.getName() == "pion+" &&
               invm_particle.getName() == "kaon-") ||
              (new_particle.getName() == "kaon-" &&
               invm_particle.getName() == "pion+") ||
              (new_particle.getName() == "pion-" &&
               invm_particle.getName() == "kaon+") ||
              (new_particle.getName() == "kaon+" &&
               invm_particle.getName() == "pion-")) {
            invm_pion_kaon_opposite_h->Fill(invariant_mass);
          }

          // invariant mass with pion+ and kaon+ or piaon- and kaon-
          if ((new_particle.getName() == "pion+" &&
               invm_particle.getName() == "kaon+") ||
              (new_particle.getName() == "kaon+" &&
               invm_particle.getName() == "pion+") ||
              (new_particle.getName() == "pion-" &&
               invm_particle.getName() == "kaon-") ||
              (new_particle.getName() == "kaon-" &&
               invm_particle.getName() == "pion-")) {
            invm_pion_kaon_same_h->Fill(invariant_mass);
          }
        }
      }
    }

    // decayed particles loop
    auto event_particles_begin = event_particles.begin();
    auto event_particles_end = event_particles.end();

    for (auto it = event_particles_begin + 100; it < event_particles_end;
         ++it) {
      auto const& decayed_particle = *it;

      for (auto invm_it = event_particles_begin; invm_it < it; ++invm_it) {
        auto const& invm_particle = *invm_it;

        if (invm_particle.getName() == "k*") {
          continue;
        }

        auto invariant_mass = decayed_particle.getInvariantMass(invm_particle);

        // invariant mass with all particles
        invm_all_h->Fill(invariant_mass);

        // invariant mass with opposite charge particles
        if (invm_particle.getCharge() * decayed_particle.getCharge() < 0) {
          invm_opposite_charge_h->Fill(invariant_mass);
        }

        // invariant mass with same charge particles
        if (invm_particle.getCharge() * decayed_particle.getCharge() > 0) {
          invm_same_charge_h->Fill(invariant_mass);
        }

        // invariant mass with pion+ and kaon- or pion- and kaon+
        if ((decayed_particle.getName() == "pion+" &&
             invm_particle.getName() == "kaon-") ||
            (decayed_particle.getName() == "kaon-" &&
             invm_particle.getName() == "pion+") ||
            (decayed_particle.getName() == "pion-" &&
             invm_particle.getName() == "kaon+") ||
            (decayed_particle.getName() == "kaon+" &&
             invm_particle.getName() == "pion-")) {
          invm_pion_kaon_opposite_h->Fill(invariant_mass);
        }

        // invariant mass with pion+ and kaon+ or piaon- and kaon-
        if ((decayed_particle.getName() == "pion+" &&
             invm_particle.getName() == "kaon+") ||
            (decayed_particle.getName() == "kaon+" &&
             invm_particle.getName() == "pion+") ||
            (decayed_particle.getName() == "pion-" &&
             invm_particle.getName() == "kaon-") ||
            (decayed_particle.getName() == "kaon-" &&
             invm_particle.getName() == "pion-")) {
          invm_pion_kaon_same_h->Fill(invariant_mass);
        }
      }
    }
  }

  TFile* file = new TFile(file_name, "RECREATE");

  histo_list->Write();

  file->Close();

  gBenchmark->Show("Benchmark");
}