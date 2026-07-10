#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"

void plotPileupFractionVsFT0M(TString inputAO2D = "AO2Ddata.root",
                              TString outputFile = "PileupFractionVsFT0M.root",
                              int nBins = 200,
                              double ft0mMin = 0.0,
                              double ft0mMax = -1.0,
                              int nCentBins = 100,
                              double centMin = 0.0,
                              double centMax = 100.0)
{
  TFile input(inputAO2D);
  if (input.IsZombie() || !input.GetListOfKeys()) {
    std::cerr << "Cannot open AO2D file: " << inputAO2D << std::endl;
    return;
  }

  if (ft0mMax <= ft0mMin) {
    float multFT0M = 0.0f;
    double maxSeen = ft0mMin;
    for (auto key : *input.GetListOfKeys()) {
      TString keyName = key->GetName();
      if (!keyName.Contains("DF_")) {
        continue;
      }
      auto* collTree = dynamic_cast<TTree*>(input.Get((keyName + "/O2npcollisiontabl").Data()));
      if (!collTree || !collTree->GetBranch("fMultFT0M")) {
        continue;
      }
      collTree->SetBranchAddress("fMultFT0M", &multFT0M);
      for (Long64_t i = 0; i < collTree->GetEntries(); ++i) {
        collTree->GetEntry(i);
        if (multFT0M > maxSeen) {
          maxSeen = multFT0M;
        }
      }
    }
    ft0mMax = maxSeen > ft0mMin ? 1.05 * maxSeen : ft0mMin + 1.0;
  }

  auto* hAll = new TH1D("hFT0MAll", "all collisions;FT0M;collisions", nBins, ft0mMin, ft0mMax);
  auto* hPileup = new TH1D("hFT0MPileup", "same-bunch pileup collisions;FT0M;collisions", nBins, ft0mMin, ft0mMax);
  auto* hCentAll = new TH1D("hCentFT0MAll", "all collisions;FT0M centrality (%);collisions", nCentBins, centMin, centMax);
  auto* hCentPileup = new TH1D("hCentFT0MPileup", "same-bunch pileup collisions;FT0M centrality (%);collisions", nCentBins, centMin, centMax);

  Long64_t nDF = 0;
  Long64_t nRows = 0;
  Long64_t nPileup = 0;

  for (auto key : *input.GetListOfKeys()) {
    TString keyName = key->GetName();
    if (!keyName.Contains("DF_")) {
      continue;
    }
    ++nDF;

    auto* collTree = dynamic_cast<TTree*>(input.Get((keyName + "/O2npcollisiontabl").Data()));
    if (!collTree) {
      continue;
    }
    if (!collTree->GetBranch("fMultFT0M") ||
        !collTree->GetBranch("fCentFT0M") ||
        !collTree->GetBranch("fNoSameBunchPileup")) {
      std::cerr << "NPCollisionTABLE in " << keyName
                << " is missing fMultFT0M, fCentFT0M and/or fNoSameBunchPileup" << std::endl;
      return;
    }

    float multFT0M = 0.0f;
    float centFT0M = 0.0f;
    bool noSameBunchPileup = false;
    collTree->SetBranchAddress("fMultFT0M", &multFT0M);
    collTree->SetBranchAddress("fCentFT0M", &centFT0M);
    collTree->SetBranchAddress("fNoSameBunchPileup", &noSameBunchPileup);

    for (Long64_t i = 0; i < collTree->GetEntries(); ++i) {
      collTree->GetEntry(i);
      hAll->Fill(multFT0M);
      hCentAll->Fill(centFT0M);
      if (!noSameBunchPileup) {
        hPileup->Fill(multFT0M);
        hCentPileup->Fill(centFT0M);
        ++nPileup;
      }
      ++nRows;
    }
  }

  auto* hRatio = (TH1D*)hPileup->Clone("hRatioPileupAll");
  hRatio->SetTitle("same-bunch pileup / all;FT0M;pileup / all");
  hRatio->Divide(hAll);
  hRatio->SetMinimum(0.0);
  hRatio->SetMaximum(1.05);

  auto* hCentRatio = (TH1D*)hCentPileup->Clone("hRatioPileupAllCentFT0M");
  hCentRatio->SetTitle("same-bunch pileup / all;FT0M centrality (%);pileup / all");
  hCentRatio->Divide(hCentAll);
  hCentRatio->SetMinimum(0.0);
  hCentRatio->SetMaximum(1.05);

  auto* cRatio = new TCanvas("cPileupFractionVsFT0M", "pileup fraction vs FT0M", 900, 700);
  hRatio->Draw("E");
  auto* cCentRatio = new TCanvas("cPileupFractionVsFT0MCentrality", "pileup fraction vs FT0M centrality", 900, 700);
  hCentRatio->Draw("E");

  TFile output(outputFile, "RECREATE");
  hAll->Write();
  hPileup->Write();
  hRatio->Write();
  hCentAll->Write();
  hCentPileup->Write();
  hCentRatio->Write();
  cRatio->Write();
  cCentRatio->Write();
  output.Close();

  TString pngName = outputFile;
  pngName.ReplaceAll(".root", ".png");
  cRatio->SaveAs(pngName);
  TString centPngName = outputFile;
  centPngName.ReplaceAll(".root", "_CentFT0M.png");
  cCentRatio->SaveAs(centPngName);

  std::cout << "Wrote " << outputFile << ", " << pngName
            << " and " << centPngName << std::endl;
  std::cout << "  DF directories: " << nDF << std::endl;
  std::cout << "  collisions read: " << nRows << std::endl;
  std::cout << "  same-bunch pileup collisions: " << nPileup << std::endl;
  std::cout << "  FT0M axis: [" << ft0mMin << ", " << ft0mMax << "]" << std::endl;
  std::cout << "  FT0M centrality axis: [" << centMin << ", " << centMax << "]" << std::endl;
}
