#include <array>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TMath.h"
#include "TROOT.h"

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

  std::cout << "ENTRIES" << '\n';

  for (int i{}; i < 12; ++i) {
    auto title = histo_array[i]->GetTitle();
    auto reached_entries_i = histo_array[i]->GetEntries();

    std::cout << title << ": expected " << expected_entries[i] << ", got "
              << reached_entries_i << '\n';
  }

  auto particle_types_histo = histo_array[0];

  std::cout << "\nPERCENTAGES" << '\n';

  for (int i{}; i < 7; ++i) {
    std::cout << "Particle " << i << ": "
              << particle_types_histo->GetBinContent(i) /
                     particle_types_histo->GetEntries()
              << '\n';
  }

  TF1* azimutal_fit = new TF1("azimutal_fit", uniform, 0., TMath::Pi(), 1);
  azimutal_fit->SetParameter(0, 10000);

  histo_array[1]->Fit(azimutal_fit, "Q");

  std::cout << "\nAZIMUTAL FIT" << '\n'
            << "Parameter 0: " << azimutal_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << azimutal_fit->GetChisquare() / azimutal_fit->GetNDF() << '\n'
            << "Probability: " << azimutal_fit->GetProb() << '\n';

  TF1* polar_fit = new TF1("polar_fit", uniform, 0., TMath::TwoPi(), 1);
  polar_fit->SetParameter(0, 10000);

  histo_array[2]->Fit(polar_fit, "Q");

  std::cout << "\nPOLAR FIT" << '\n'
            << "Parameter 0: " << polar_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << polar_fit->GetChisquare() / polar_fit->GetNDF() << '\n'
            << "Probability: " << polar_fit->GetProb() << '\n';

  TF1* momentum_fit = new TF1("momentum fit", exp, 0., 9., 2);
  momentum_fit->SetParameter(0, 90000);
  momentum_fit->SetParameter(1, 1.);

  histo_array[3]->Fit(momentum_fit, "Q");

  std::cout << "\nMOMENTUM FIT" << '\n'
            << "Parameter 0: " << momentum_fit->GetParameter(0) << '\n'
            << "Chi square/NDF: "
            << momentum_fit->GetChisquare() / momentum_fit->GetNDF() << '\n'
            << "Probability" << momentum_fit->GetProb() << '\n';

  TH1F* invm_subtraction_12 = new TH1F(*(TH1F*)histo_array[7]);
  invm_subtraction_12->SetTitle(
      "Invariant mass, all (opposite charge - same charge)");
  invm_subtraction_12->SetName("invm_subtraction_12");
  invm_subtraction_12->Add(histo_array[7], histo_array[8], 1, -1);

  TH1F* invm_subtraction_34 = new TH1F(*(TH1F*)histo_array[9]);
  invm_subtraction_34->SetTitle(
      "Invariant mass, kaon & pion (opposite charge - same charge)");
  invm_subtraction_34->SetName("invm_subtraction_34");
  invm_subtraction_34->Add(histo_array[9], histo_array[10], 1, -1);

  TF1* invm_12_fit = new TF1("invm fit 12", gauss, 0., 7., 3);
  invm_12_fit->SetParameter(0, 7.996);
  invm_12_fit->SetParameter(1, 0.8919);
  invm_12_fit->SetParameter(2, 0.04989);

  invm_subtraction_12->Fit(invm_12_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN ALL PARTICLES (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << "Parameter 0: " << invm_12_fit->GetParameter(0) << '\n'
            << "Parameter 1: " << invm_12_fit->GetParameter(1) << '\n'
            << "Parameter 2: " << invm_12_fit->GetParameter(2) << '\n'
            << "Chi square/NDF :"
            << invm_12_fit->GetChisquare() / invm_12_fit->GetNDF() << '\n'
            << "Probability" << invm_12_fit->GetProb() << '\n'
            << "K* mass: " << invm_12_fit->GetParameter(1) << '\n'
            << "K* width: " << invm_12_fit->GetParameter(2) << '\n';

  TF1* invm_34_fit = new TF1("invm fit 12", gauss, 0., 7., 3);
  invm_34_fit->SetParameter(0, 7.996);
  invm_34_fit->SetParameter(1, 0.8919);
  invm_34_fit->SetParameter(2, 0.04989);

  invm_subtraction_34->Fit(invm_34_fit, "Q");

  std::cout << "\nINVARIANT MASS BETWEEN KAON AND PION (OPPOSITE CHARGE - SAME "
               "CHARGE) FIT"
            << '\n'
            << "Parameter 0: " << invm_34_fit->GetParameter(0) << '\n'
            << "Parameter 1: " << invm_34_fit->GetParameter(1) << '\n'
            << "Parameter 2: " << invm_34_fit->GetParameter(2) << '\n'
            << "Chi square/NDF :"
            << invm_34_fit->GetChisquare() / invm_34_fit->GetNDF() << '\n'
            << "Probability" << invm_34_fit->GetProb() << '\n'
            << "K* mass: " << invm_34_fit->GetParameter(1) << '\n'
            << "K* width: " << invm_34_fit->GetParameter(2) << '\n';

  TCanvas* particles_canvas = new TCanvas();
  particles_canvas->Divide(2, 2);

  particles_canvas->cd(1);
  histo_array[0]->Draw();

  particles_canvas->cd(2);
  histo_array[3]->Draw();

  particles_canvas->cd(3);
  histo_array[1]->Draw();

  particles_canvas->cd(4);
  histo_array[2]->Draw();
}