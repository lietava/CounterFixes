#include <random>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include <vector>
#include <iostream>

double sigma_pt(double pt) {
  return pt * sqrt((0.02 / pt) * (0.02 / pt) + 0.01 * 0.01);
}

void test4()
{
  //gSystem->AddIncludePath("-I/Users/romanlietava/anal/RooUnfold/build/");
  //gSystem->Load("/Users/romanlietava/anal/RooUnfold/build/libRooUnfold.dylib");
  TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", 0, 10);
  fpt->SetParameter(0, 0.5);

  int nEvents = 10000;
  double mu = 50.;
  int r = 20;
  double p = r / (r + mu);

  std::mt19937 gen(12345);
  std::negative_binomial_distribution<int> nbd(r, p);

  TRandom3 rand(0);
  constexpr int NPT = 2;
  constexpr int NMULT = 2;
  TH2D* hTrue = new TH2D("hTrue","",NPT,0,10,NMULT,0,100);
  TH2D* hReco = new TH2D("hReco","",NPT,0,10,NMULT,0,100);

  RooUnfoldResponse response(hReco, hTrue, "response", "response");

  double eff = 0.8;
  double fakeMean = 0.5;

  struct TrackInfo {
    double ptTrue;
    double ptReco;
  bool isReco;
  };

  for (int ev = 0; ev < nEvents; ++ev) {
    int multTrue = nbd(gen);

    std::vector<TrackInfo> tracks;
    tracks.reserve(multTrue);

    int multReco = 0;

    for (int i = 0; i < multTrue; ++i) {
      double ptTrue = fpt->GetRandom();
      bool isReco = (rand.Rndm() < eff);
      double ptReco = -1.0;

      if (isReco) {
        ptReco = rand.Gaus(ptTrue, sigma_pt(ptTrue));
        if (ptReco > 0.0 && ptReco < 10.0) {
          ++multReco;
        } else {
          isReco = false;
          ptReco = -1.0;
        }
      }

      tracks.push_back({ptTrue, ptReco, isReco});
    }

    int nFake = rand.Poisson(fakeMean);
    multReco += nFake;

    for (const auto& trk : tracks) {
      hTrue->Fill(trk.ptTrue, multTrue);

      if (trk.isReco) {
        hReco->Fill(trk.ptReco, multReco);
        response.Fill(trk.ptReco, multReco, trk.ptTrue, multTrue);
      } else {
        response.Miss(trk.ptTrue, multTrue);
      }
    }

    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();
      //hReco->Fill(ptFake, multReco);
      //response.Fake(ptFake, multReco);
    }
    //
  }

  RooUnfoldBayes unfold(&response, hReco, 1000);
  TH2* hUnfold = dynamic_cast<TH2*>(unfold.Hreco());
  if (!hUnfold) {
    std::cout << "Error: unfolded histogram is not TH2" << std::endl;
    return;
  }

  hUnfold->SetName("hUnfold");
  hUnfold->SetTitle("Unfolded; p_{T}; mult");

  TFile* fout = new TFile("output.root", "RECREATE");
  hTrue->Write();
  hReco->Write();
  response.Write("response2D");
  hUnfold->Write();
  fout->Close();
  hTrue->Print("all");
  hUnfold->Print("all");
  //hMatrix->Draw("COLZ");
}