#include <array>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"

// Define repeated constants
const char* HEIGHT_LABEL = "Height (p0)";
const char* MEAN_LABEL = "Mean (p1)";
const char* STDDEV_LABEL = "Std. dev. (p2)";
const char* PARTICLE_NAMES[7]{"pion+",   "pion-",   "kaon+", "kaon-",
                              "proton+", "proton-", "K*"};

// Define fit functions
Double_t gauss(Double_t* xx, Double_t* par) {
  Double_t x = xx[0];
  Double_t val =
      par[0] * TMath::Exp(-(x - par[1]) * (x - par[1]) / (2 * par[2] * par[2]));
  return val;
}

Double_t exp(Double_t* xx, Double_t* par) {
  Double_t x = xx[0];
  Double_t val = par[0] * TMath::Exp(-x / par[1]);
  return val;
}

Double_t uniform(Double_t* xx, Double_t* par) { return par[0]; }

void analyse(const char* file_name) {
  // Set histogram options
  gStyle->SetOptFit(001);
  gStyle->SetOptStat("e");

  TFile* file = new TFile(file_name, "READ");

  std::array<TH1*, 12> histo_array;

  auto key_list = file->GetListOfKeys();

  for (int i{}; i < 12; ++i) {
    auto key = (TKey*)key_list->At(i);
    auto histo = (TH1*)file->Get(key->GetName());
    histo_array[i] = histo;
  }

  std::array<double, 12> expected_entries{
      1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 5e8, 2.5e8, 2.5e8, 4.46e7, 4.45e7, 1e5};

  // Print expected and real entries
  std::cout << "ENTRIES" << '\n';

  for (int i{}; i < 12; ++i) {
    auto title = histo_array[i]->GetTitle();
    auto reached_entries_i = histo_array[i]->GetEntries();

    std::cout << title << ": expected " << expected_entries[i] << ", got "
              << reached_entries_i << '\n';
  }

  auto particle_types_histogram = histo_array[0];

  // Print generated percentages
  std::cout << "\nPARTICLE OCCURRENCES" << '\n';

  for (int i{}; i < 7; ++i) {
    auto occurrences = particle_types_histogram->GetBinContent(i + 1);

    std::cout << "Particle " << PARTICLE_NAMES[i] << ": " << occurrences << " ("
              << occurrences / particle_types_histogram->GetEntries() * 100
              << "%)" << '\n';
  }

  // Create first canvas, with particle types, angles and momentum
  TCanvas* particles_canvas = new TCanvas();
  particles_canvas->Divide(2, 2);

  particles_canvas->cd(1);
  auto x_axis = particle_types_histogram->GetXaxis();
  x_axis->CenterLabels();
  for (int i{0}; i < 7; ++i) {
    x_axis->SetBinLabel(i + 1, PARTICLE_NAMES[i]);
  }
  particle_types_histogram->SetXTitle("Particle type");
  particle_types_histogram->SetYTitle("Occurrences");
  particle_types_histogram->Draw();

  TF1* azimutal_fit = new TF1("azimutal_fit", uniform, 0., TMath::Pi(), 1);
  azimutal_fit->SetParameter(0, 10000);
  azimutal_fit->SetParName(0, HEIGHT_LABEL);

  particles_canvas->cd(2);
  auto azimutal_angles_histogram = histo_array[1];
  azimutal_angles_histogram->SetXTitle("Azimutal angle (rad)");
  azimutal_angles_histogram->SetYTitle("Occurrences");
  azimutal_angles_histogram->Fit(azimutal_fit, "Q");

  std::cout << "\nAZIMUTAL FIT" << '\n'
            << HEIGHT_LABEL << ": " << azimutal_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << azimutal_fit->GetChisquare() / azimutal_fit->GetNDF() << '\n'
            << "Probability: " << azimutal_fit->GetProb() << '\n';

  TF1* polar_fit = new TF1("polar_fit", uniform, 0., TMath::TwoPi(), 1);
  polar_fit->SetParameter(0, 10000);
  polar_fit->SetParName(0, HEIGHT_LABEL);

  particles_canvas->cd(3);
  auto polar_angles_histogram = histo_array[2];
  polar_angles_histogram->SetXTitle("Polar angle (rad)");
  polar_angles_histogram->SetYTitle("Occurrences");
  polar_angles_histogram->Fit(polar_fit, "Q");

  std::cout << "\nPOLAR FIT" << '\n'
            << HEIGHT_LABEL << ": " << polar_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << polar_fit->GetChisquare() / polar_fit->GetNDF() << '\n'
            << "Probability: " << polar_fit->GetProb() << '\n';

  TF1* momentum_fit = new TF1("momentum fit", exp, 0., 9., 2);
  momentum_fit->SetParameter(0, 90000);
  momentum_fit->SetParName(0, HEIGHT_LABEL);
  momentum_fit->SetParameter(1, 1.);
  momentum_fit->SetParName(1, MEAN_LABEL);

  particles_canvas->cd(4);
  auto momentum_histogram = histo_array[3];
  momentum_histogram->SetXTitle("Momentum (GeV)");
  momentum_histogram->SetYTitle("Occurrences");
  momentum_histogram->Fit(momentum_fit, "Q");

  std::cout << "\nMOMENTUM FIT" << '\n'
            << HEIGHT_LABEL << ": " << momentum_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << momentum_fit->GetChisquare() / momentum_fit->GetNDF() << '\n'
            << "Probability: " << momentum_fit->GetProb() << '\n';

  // Create second canvas, with invariant mass
  TCanvas* inv_mass_canvas = new TCanvas();
  inv_mass_canvas->Divide(2, 2);
  inv_mass_canvas->cd(1);
  TF1* k_star_fit = new TF1("k_star_fit", gauss, 0.6, 1.2, 3);
  k_star_fit->SetParameter(0, 500);
  k_star_fit->SetParName(0, HEIGHT_LABEL);
  k_star_fit->SetParameter(1, 0.9);
  k_star_fit->SetParName(1, MEAN_LABEL);
  k_star_fit->SetParameter(2, 0.05);
  k_star_fit->SetParName(2, STDDEV_LABEL);
  auto invm_decayed_h = histo_array[11];
  invm_decayed_h->SetXTitle("Invariant mass (GeV)");
  invm_decayed_h->SetYTitle("Occurrences");
  invm_decayed_h->Fit(k_star_fit, "Q");

  auto invm_opposite_charge_h = histo_array[7];
  auto invm_same_charge_h = histo_array[8];
  TH1F* invm_subtraction_12 = new TH1F(*(TH1F*)invm_opposite_charge_h);
  invm_subtraction_12->SetTitle(
      "Invariant mass, all (opposite charge - same charge)");
  invm_subtraction_12->SetName("invm_subtraction_12");
  invm_subtraction_12->SetXTitle("Invariant mass (GeV)");
  invm_subtraction_12->SetYTitle("Occurrences");
  invm_subtraction_12->Add(invm_opposite_charge_h, invm_same_charge_h, 1, -1);
  invm_subtraction_12->SetEntries(invm_opposite_charge_h->GetEntries());
  invm_subtraction_12->SetAxisRange(0.6, 1.2);

  TF1* invm_12_fit = new TF1("invm fit 12", gauss, 0., 7., 3);
  invm_12_fit->SetParameter(0, 7.996);
  invm_12_fit->SetParName(0, HEIGHT_LABEL);
  invm_12_fit->SetParameter(1, 0.8919);
  invm_12_fit->SetParName(1, MEAN_LABEL);
  invm_12_fit->SetParameter(2, 0.04989);
  invm_12_fit->SetParName(2, STDDEV_LABEL);

  inv_mass_canvas->cd(2);
  invm_subtraction_12->Fit(invm_12_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN ALL PARTICLES (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << HEIGHT_LABEL << ": " << invm_12_fit->GetParameter(0) << '\n'
            << MEAN_LABEL << ": " << invm_12_fit->GetParameter(1) << '\n'
            << STDDEV_LABEL << ": " << invm_12_fit->GetParameter(2) << '\n'
            << "Chi square/NDF :"
            << invm_12_fit->GetChisquare() / invm_12_fit->GetNDF() << '\n'
            << "Probability: " << invm_12_fit->GetProb() << '\n'
            << "K* mass: " << invm_12_fit->GetParameter(1) << '\n'
            << "K* width: " << invm_12_fit->GetParameter(2) << '\n';

  auto invm_pion_kaon_opposite_h = histo_array[9];
  auto invm_pion_kaon_same_h = histo_array[10];
  TH1F* invm_subtraction_34 = new TH1F(*(TH1F*)invm_pion_kaon_opposite_h);
  invm_subtraction_34->SetTitle(
      "Invariant mass, kaon & pion (opposite charge - same charge)");
  invm_subtraction_34->SetName("invm_subtraction_34");
  invm_subtraction_34->SetXTitle("Invariant mass (GeV)");
  invm_subtraction_34->SetYTitle("Occurrences");
  invm_subtraction_34->Add(invm_pion_kaon_opposite_h, invm_pion_kaon_same_h, 1,
                           -1);
  invm_subtraction_34->SetEntries(invm_pion_kaon_opposite_h->GetEntries());
  invm_subtraction_34->SetAxisRange(0.6, 1.2);

  TF1* invm_34_fit = new TF1("invm fit 34", gauss, 0., 7., 3);
  invm_34_fit->SetParameter(0, 7.996);
  invm_34_fit->SetParName(0, HEIGHT_LABEL);
  invm_34_fit->SetParameter(1, 0.8919);
  invm_34_fit->SetParName(1, MEAN_LABEL);
  invm_34_fit->SetParameter(2, 0.04989);
  invm_34_fit->SetParName(2, STDDEV_LABEL);

  inv_mass_canvas->cd(3);
  invm_subtraction_34->Fit(invm_34_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN KAON AND PION (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << HEIGHT_LABEL << ": " << invm_34_fit->GetParameter(0) << '\n'
            << MEAN_LABEL << ": " << invm_34_fit->GetParameter(1) << '\n'
            << STDDEV_LABEL << ": " << invm_34_fit->GetParameter(2) << '\n'
            << "Chi square/NDF :"
            << invm_34_fit->GetChisquare() / invm_34_fit->GetNDF() << '\n'
            << "Probability: " << invm_34_fit->GetProb() << '\n'
            << "K* mass: " << invm_34_fit->GetParameter(1) << '\n'
            << "K* width: " << invm_34_fit->GetParameter(2) << '\n';
}