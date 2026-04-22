#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TClass.h>

#include <iostream>
#include <memory>
#include <string>

void readNPtablesRDF(const char* fname = "AO2D.root")
{
  std::unique_ptr<TFile> f(TFile::Open(fname));
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file " << fname << std::endl;
    return;
  }

  TChain collChain("O2npcollisiontable");
  TChain recoChain("O2nprecochargedcand");

  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;

  int nColl = 0;
  int nReco = 0;

  while ((key = (TKey*)next())) {
    TString name = key->GetName();
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl || !cl->InheritsFrom(TDirectory::Class())) {
      continue;
    }
    if (!name.BeginsWith("DF_")) {
      continue;
    }

    TString collPath = Form("%s/%s/O2npcollisiontabl", fname, name.Data());
    TString recoPath = Form("%s/%s/O2nprecochargedca", fname, name.Data());

    if (collChain.Add(collPath.Data(), 0) > 0) {
      ++nColl;
      std::cout << "Added " << collPath << "\n";
    }
    if (recoChain.Add(recoPath.Data(), 0) > 0) {
      ++nReco;
      std::cout << "Added " << recoPath << "\n";
    }
  }

  std::cout << "collision trees added: " << nColl << "\n";
  std::cout << "reco trees added:      " << nReco << "\n";

  if (nColl == 0) {
    std::cerr << "No collision trees found\n";
    return;
  }
  if (nReco == 0) {
    std::cerr << "No reco trees found\n";
    return;
  }

  ROOT::RDataFrame dfColl(collChain);
  ROOT::RDataFrame dfReco(recoChain);

  std::cout << "Collision entries = " << *dfColl.Count() << "\n";
  std::cout << "Reco entries      = " << *dfReco.Count() << "\n";

  std::cout << "\nCollision columns:\n";
  for (auto const& c : dfColl.GetColumnNames()) {
    std::cout << c << "\n";
  }

  std::cout << "\nReco columns:\n";
  for (auto const& c : dfReco.GetColumnNames()) {
    std::cout << c << "\n";
  }

  std::cout << "\nCollision preview:\n";
  dfColl.Display("",10)->Print();

  std::cout << "\nReco preview:\n";
  dfReco.Display("",10)->Print();
}