#include <random>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooUnfoldResponse.h"

double sigma_pt(double pt) {
  return pt * sqrt((0.02 / pt) * (0.02 / pt) + 0.01 * 0.01);
}

void test3()
{
  //gSystem->AddIncludePath("-I/Users/romanlietava/anal/RooUnfold/build/");
  //gSystem->Load("/Users/romanlietava/anal/RooUnfold/build/libRooUnfold.dylib");
  TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", 0, 10);
  fpt->SetParameter(0, 0.5);

  double mu = 50.;
  int r = 20;
  double p = r / (r + mu);

  std::mt19937 gen(12345);
  std::negative_binomial_distribution<int> nbd(r, p);

  TRandom3 rand(0);
  constexpr int NDIM = 10;

  TH1F* hPtTrue   = new TH1F("hPtTrue",   "True p_{T};p_{T};Counts", NDIM, 0, 10);
  TH1F* hPtReco   = new TH1F("hPtReco",   "Reco p_{T};p_{T};Counts", NDIM, 0, 10);
  TH1F* hPtRes    = new TH1F("hPtRes",    "p_{T}-p_{T}^{rec};#Deltap_{T};Counts", 1000, -1, 1);
  TH1F* hMult     = new TH1F("hMult",     "Multiplicity Gen;N;Counts", 100, 0, 100);
  TH1F* hMultRec  = new TH1F("hMultRec",  "Multiplicity Rec;N;Counts", 100, 0, 100);
  TH2F* hMatrix   = new TH2F("hMatrix",   "Response matrix;Reco p_{T};True p_{T}", NDIM, 0, 10, NDIM, 0, 10);

  RooUnfoldResponse response(hPtReco, hPtTrue, "response", "response");

  double eff = 0.8;
  double fakeMean = 0.5;

  for (int ev = 0; ev < 10000; ++ev) {
    int mult = nbd(gen);
    hMult->Fill(mult);

    int multRec = 0;

    // true particles
    for (int i = 0; i < mult; ++i) {
      double ptTrue = fpt->GetRandom();
      hPtTrue->Fill(ptTrue);

      if (rand.Rndm() < eff) {
        double ptReco = rand.Gaus(ptTrue, sigma_pt(ptTrue));
        if (ptReco > 0.0 && ptReco < 10.0) {
          hPtReco->Fill(ptReco);
          hPtRes->Fill(ptTrue - ptReco);
          hMatrix->Fill(ptReco, ptTrue);
          response.Fill(ptReco, ptTrue);
          ++multRec;
        } else {
          // reconstructed outside analysis range
          response.Miss(ptTrue);
        }
      } else {
        // inefficiency: true track not reconstructed
        response.Miss(ptTrue);
      }
    }

    // fake reco tracks
    int nFake = rand.Poisson(fakeMean);
    multRec += nFake;

    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();   // toy fake spectrum
      if (ptFake > 0.0 && ptFake < 10.0) {
        //hPtReco->Fill(ptFake);
        //response.Fake(ptFake);
      }
    }

    hMultRec->Fill(multRec);
  }
  RooUnfoldBayes unfold(&response, hPtReco, 100);
  unfold.SetVerbose(0);
  TH1* hUnfold = dynamic_cast<TH1*>(unfold.Hunfold());
  if (!hUnfold) {
    std::cout << "Error: unfolded histogram is not TH1" << std::endl;
    return;
  }

  hPtTrue->Print("all");
  hUnfold->Print("all");
  TFile* fout = new TFile("output.root", "RECREATE");
  hPtTrue->Write();
  hPtReco->Write();
  hPtRes->Write();
  hMult->Write();
  hMultRec->Write();
  hMatrix->Write();
  response.Write("response");
  fout->Close();

  //hMatrix->Draw("COLZ");
}