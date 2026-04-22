// readNPtablesFragments.C
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

constexpr int NREAD = 100;

bool startsWith(const std::string& s, const std::string& prefix)
{
  return s.rfind(prefix, 0) == 0;
}

std::vector<std::string> getDFdirs(TFile* f)
{
  std::vector<std::string> out;

  TIter next(f->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    std::string name = key->GetName();
    std::string cls  = key->GetClassName();
    if (cls == "TDirectoryFile" && startsWith(name, "DF_")) {
      out.push_back(name);
    }
  }
  return out;
}

void readNPCollsFromTree(TTree* t)
{
  if (!t) {
    return;
  }

  std::cout << "\nReading tree: " << t->GetName()
            << " entries=" << t->GetEntries() << std::endl;

  TTreeReader r(t);

  TTreeReaderValue<int> runNumber(r, "fRunNumber");
  TTreeReaderValue<unsigned long long> globalBC(r, "fGlobalBC");
  TTreeReaderValue<unsigned short> numContrib(r, "fNumContrib");
  TTreeReaderValue<int> multNTracksGlobal(r, "fMultNTracksGlobal");
  TTreeReaderValue<float> centFT0M(r, "fCentFT0M");
  TTreeReaderValue<float> multFT0M(r, "fMultFT0M");

  Long64_t i = 0;
  while (r.Next()) {
    std::cout
      << "row=" << i
      << " run=" << *runNumber
      << " globalBC=" << *globalBC
      << " numContrib=" << *numContrib
      << " multNTracksGlobal=" << *multNTracksGlobal
      << " centFT0M=" << *centFT0M
      << " multFT0M=" << *multFT0M
      << std::endl;
    ++i;
    if(i > NREAD) {
      //break;
    }
  }
}

void readNPRecoFromTree(TTree* t)
{
  if (!t) {
    return;
  }

  std::cout << "\nReading tree: " << t->GetName()
            << " entries=" << t->GetEntries() << std::endl;

  TTreeReader r(t);

  TTreeReaderValue<int> collIdx(r, "fCollIdx");
  TTreeReaderValue<float> pt(r, "fPt");

  Long64_t i = 0;
  while (r.Next()) {
    std::cout
      << "row=" << i
      << " collIdx=" << *collIdx
      << " pt=" << *pt
      << std::endl;
    ++i;
    if(i > NREAD) {
      //break;
    }
  }
}

//void readNPtables(const char* filename = "data/659826/AO2D.root",
void readNPtables(const char* filename = "AO2D.root",
                           const char* collTreeName = "O2npcollisiontabl",
                           const char* recoTreeName = "O2mprecochargedca")
{
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::cout << "Top-level content:\n";
  f->ls();

  auto dfdirs = getDFdirs(f);

  std::cout << "\nFound " << dfdirs.size() << " DF directories\n";
  for (auto const& dname : dfdirs) {
    std::cout << "  " << dname << std::endl;
  }

  int nread = 0;
  for (auto const& dname : dfdirs) {
    TDirectory* dir = dynamic_cast<TDirectory*>(f->Get(dname.c_str()));
    if (!dir) {
      continue;
    }

    std::cout << "\n==============================" << std::endl;
    std::cout << "Directory: " << dname << std::endl;
    dir->ls();

    TTree* tColl = dynamic_cast<TTree*>(dir->Get(collTreeName));
    if (tColl) {
      std::cout << "\nFound " << collTreeName << " in " << dname << std::endl;
      readNPCollsFromTree(tColl);
    }

    TTree* tReco = dynamic_cast<TTree*>(dir->Get(recoTreeName));
    if (tReco) {
      std::cout << "\nFound " << recoTreeName << " in " << dname << std::endl;
      readNPRecoFromTree(tReco);
    }
  }

  f->Close();
}