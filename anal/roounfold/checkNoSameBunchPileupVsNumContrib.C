#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

void checkNoSameBunchPileupVsNumContrib(TString inputAO2D = "AO2Ddata.root",
                                        TString outputFile = "NoSameBunchPileupVsNumContrib.root")
{
  TFile input(inputAO2D);
  if (input.IsZombie() || !input.GetListOfKeys()) {
    std::cerr << "Cannot open AO2D file: " << inputAO2D << std::endl;
    return;
  }

  auto* hNoSameVsNumContrib = new TH2D("hNoSameVsNumContrib",
                                       "NoSameBunchPileup vs NumContrib;NumContrib;NoSameBunchPileup",
                                       501, -0.5, 500.5, 2, -0.5, 1.5);
  auto* hNumContribAll = new TH1D("hNumContribAll",
                                  "NumContrib all collisions;NumContrib;collisions",
                                  501, -0.5, 500.5);
  auto* hNumContribNoSame = new TH1D("hNumContribNoSame",
                                     "NumContrib, NoSameBunchPileup true;NumContrib;collisions",
                                     501, -0.5, 500.5);
  auto* hNumContribSame = new TH1D("hNumContribSame",
                                   "NumContrib, NoSameBunchPileup false;NumContrib;collisions",
                                   501, -0.5, 500.5);

  Long64_t nDF = 0;
  Long64_t nMissingTables = 0;
  Long64_t nRows = 0;
  Long64_t nNoSame = 0;
  Long64_t nSame = 0;

  for (auto key : *input.GetListOfKeys()) {
    TString keyName = key->GetName();
    if (!keyName.Contains("DF_")) {
      continue;
    }
    ++nDF;

    auto* collTree = dynamic_cast<TTree*>(input.Get((keyName + "/O2npcollisiontabl").Data()));
    if (!collTree) {
      ++nMissingTables;
      continue;
    }

    if (!collTree->GetBranch("fNumContrib") || !collTree->GetBranch("fNoSameBunchPileup")) {
      std::cerr << "NPCollisionTABLE in " << keyName
                << " is missing fNumContrib and/or fNoSameBunchPileup" << std::endl;
      return;
    }

    UShort_t numContrib = 0;
    bool noSameBunchPileup = false;
    collTree->SetBranchAddress("fNumContrib", &numContrib);
    collTree->SetBranchAddress("fNoSameBunchPileup", &noSameBunchPileup);

    for (Long64_t i = 0; i < collTree->GetEntries(); ++i) {
      collTree->GetEntry(i);
      const int noSame = noSameBunchPileup ? 1 : 0;

      hNoSameVsNumContrib->Fill(numContrib, noSame);
      hNumContribAll->Fill(numContrib);
      if (noSameBunchPileup) {
        hNumContribNoSame->Fill(numContrib);
        ++nNoSame;
      } else {
        hNumContribSame->Fill(numContrib);
        ++nSame;
      }
      ++nRows;
    }
  }

  if (nDF <= 0) {
    std::cerr << "No DF_* directories found in " << inputAO2D << std::endl;
    return;
  }
  if (nRows <= 0) {
    std::cerr << "No rows read from NPCollisionTABLE in " << inputAO2D << std::endl;
    return;
  }

  TFile output(outputFile, "RECREATE");
  hNoSameVsNumContrib->Write();
  hNumContribAll->Write();
  hNumContribNoSame->Write();
  hNumContribSame->Write();
  output.Close();

  std::cout << "Wrote " << outputFile << std::endl;
  std::cout << "  DF directories: " << nDF
            << " missing NPCollisionTABLE: " << nMissingTables << std::endl;
  std::cout << "  collisions read: " << nRows << std::endl;
  std::cout << "  fNoSameBunchPileup true: " << nNoSame << std::endl;
  std::cout << "  fNoSameBunchPileup false: " << nSame << std::endl;
}
