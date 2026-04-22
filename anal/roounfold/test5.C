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

constexpr int NPT = 2;
constexpr int NMULT = 2;

int flattenTrue(int iPtTrue, int iMultTrue)
{
  return iPtTrue + NPT * iMultTrue;
}

int flattenReco(int iPtReco, int iMultReco)
{
  return iPtReco + NPT * iMultReco;
}

double sigma_pt(double pt) {
  return pt * sqrt((0.02 / pt) * (0.02 / pt) + 0.01 * 0.01);
}

void test5()
{
  double R[NPT][NMULT][NPT][NMULT] = {};
  double t[NPT][NMULT] = {};
  double t_new[NPT][NMULT] = {};
  double d[NPT][NMULT] = {};

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
 
  TH2D* hTrue = new TH2D("hTrue","",NPT,0,10,NMULT,0,100);
  TH2D* hReco = new TH2D("hReco","",NPT,0,10,NMULT,0,100);
  //
  int nRecoBins = NPT*NMULT;
  int nTrueBins = NPT*NMULT; 
  TH1D* hTrueF = new TH1D("hTrueF","", nTrueBins, 0.5, nTrueBins + 0.5);
  TH1D* hRecoF = new TH1D("hRecoF","", nRecoBins, 0.5, nRecoBins + 0.5);
  TH2D* hRespFlat = new TH2D("hRespFlat",
                           "response;reco;true",
                           nRecoBins, 0.5, nRecoBins + 0.5,
                           nTrueBins, 0.5, nTrueBins + 0.5);
  //
  const Int_t dim = 4;
 // binning
  Int_t bins[dim]    = {20, 10, 20, 10};   // ptR, multR, ptT, multT
  Double_t xmin[dim] = {0,  0,  0,  0};
  Double_t xmax[dim] = {10, 100, 10, 100};

  THnSparseD* hResp4D = new THnSparseD("hResp4D",
                                      "4D response;ptR;multR;ptT;multT",
                                      dim, bins, xmin, xmax);
  //
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
        hResp4D->Fill(trk.ptReco, multReco, trk.ptTrue, multTrue);
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
  //
  //
  RooUnfoldBayes unfold(&response, hReco, 1000);
  unfold.SetVerbose(0);
  TH2* hUnfold = dynamic_cast<TH2*>(unfold.Hunfold());
  if (!hUnfold) {
    std::cout << "Error: unfolded histogram is not TH2" << std::endl;
    return;
  }

  hUnfold->SetName("hUnfold");
  hUnfold->SetTitle("Unfolded; p_{T}; mult");
  // 
  // Flattened
  for(int i = 0; i < nTrueBins; i ++) {
    int iPt   = i % NPT;
    int iMult = i / NPT;
    double cont = hTrue->GetBinContent(iPt + 1, iMult + 1);
    hTrueF->SetBinContent(i + 1,cont);
  }
  for(int i = 0; i < nRecoBins; i ++) {
    int iPt   = i % NPT;
    int iMult = i / NPT;
    double cont = hReco->GetBinContent(iPt + 1, iMult + 1);
    hRecoF->SetBinContent(i + 1,cont);
  }
  RooUnfoldResponse responseF(hRecoF, hTrueF, "response", "response");
  for(int i = 0; i < nRecoBins; i++){
    int ipt   = i % NPT;
    int imult = i / NPT;
    for(int j = 0; j < nTrueBins; j++){
      int jpt   = j % NPT;
      int jmult = j / NPT;
      //double cont = hResp4D(ipt, imult, jpt, jmult);
      int coord[4];
      coord[0] = ipt + 1;
      coord[1] = imult + 1;
      coord[2] = jpt + 1;
      coord[3] = jmult + 1;

      Long64_t globalBin = hResp4D->GetBin(coord);
      double cont = hResp4D->GetBinContent(globalBin);
      R[ipt][imult][jpt][jmult] = cont;
      //hRespFlat->SetBinContent(i + 1, j + 1, cont);
      responseF.Fill(i + 0.5, j + 0.5, cont);
    }
  }
  RooUnfoldBayes unfoldF(&responseF, hRecoF, 1000);
  unfoldF.SetVerbose(0);
  TH1* hUnfoldF = dynamic_cast<TH1*>(unfold.Hunfold());
  if (!hUnfold) {
    std::cout << "Error: unfolded histogram is not TH2" << std::endl;
    return;
  }
  hUnfoldF->SetName("hUnfoldF");
  hUnfoldF->SetTitle("Unfolded F; p_{T}; mult");
  hUnfoldF->Print("all");
  //return;
  //
  TFile* fout = new TFile("output.root", "RECREATE");
  hTrue->Write();
  hReco->Write();
  response.Write("response2D");
  hUnfold->Write();
  hResp4D->Write();
  fout->Close();
  hTrue->Print("all");
  hUnfold->Print("all");
  //hMatrix->Draw("COLZ");
  //
  for (int ix_t = 0; ix_t < NPT; ++ix_t) {
  for (int iy_t = 0; iy_t < NMULT; ++iy_t) {

    double sum = 0.0;

    for (int ix_r = 0; ix_r < NPT; ++ix_r) {
      for (int iy_r = 0; iy_r < NMULT; ++iy_r) {

        double denom = 0.0;

        for (int kx = 0; kx < NPT; ++kx) {
          for (int ky = 0; ky < NMULT; ++ky) {
            denom += R[ix_r][iy_r][kx][ky] * t[kx][ky];
          }
        }

        if (denom > 0) {
          sum += R[ix_r][iy_r][ix_t][iy_t] * d[ix_r][iy_r] / denom;
        }
      }
    }

    t_new[ix_t][iy_t] = t[ix_t][iy_t] * sum;
  }
}
}