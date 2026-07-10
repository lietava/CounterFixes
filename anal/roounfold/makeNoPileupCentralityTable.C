#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TString.h"
#include "TTree.h"

struct NpCentBCKey {
  int runNumber;
  ULong64_t globalBC;

  bool operator==(const NpCentBCKey& other) const
  {
    return runNumber == other.runNumber && globalBC == other.globalBC;
  }
};

struct NpCentBCKeyHash {
  std::size_t operator()(const NpCentBCKey& key) const
  {
    const std::size_t h1 = std::hash<int>{}(key.runNumber);
    const std::size_t h2 = std::hash<ULong64_t>{}(key.globalBC);
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
  }
};

struct NpCentCollision {
  int dfIndex = -1;
  TString dfName;
  Long64_t collisionIndex = -1;
  int runNumber = 0;
  ULong64_t globalBC = 0;
  UShort_t numContrib = 0;
  int multNTracksGlobal = 0;
  int multNTracksNP = 0;
  float multFT0M = 0.0f;
  float centFT0M = 0.0f;
  float centFT0MNoPileup = 101.0f;
  bool noSameBunchPileup = false;
};

void copyDirectoryWithEnrichedCollisionTable(TDirectory* inputDir,
                                             TDirectory* outputDir,
                                             const std::vector<NpCentCollision>& collisions)
{
  for (auto keyObject : *inputDir->GetListOfKeys()) {
    auto* key = dynamic_cast<TKey*>(keyObject);
    if (!key) {
      continue;
    }

    TObject* object = key->ReadObj();
    if (!object) {
      continue;
    }

    const TString objectName = object->GetName();
    const TString inputPath = inputDir->GetPath();
    const bool isDirectory = object->InheritsFrom(TDirectory::Class());
    const bool isNpcollisionTable =
      object->InheritsFrom(TTree::Class()) &&
      objectName == "O2npcollisiontabl" &&
      inputPath.Contains("/DF_");

    if (isDirectory) {
      auto* inputSubDir = dynamic_cast<TDirectory*>(object);
      auto* outputSubDir = outputDir->mkdir(objectName);
      if (!outputSubDir) {
        std::cerr << "Cannot create output directory " << objectName << std::endl;
        delete object;
        return;
      }
      copyDirectoryWithEnrichedCollisionTable(inputSubDir, outputSubDir, collisions);
      delete object;
      continue;
    }

    outputDir->cd();

    if (isNpcollisionTable) {
      auto* inputTree = dynamic_cast<TTree*>(object);
      TString dfName = inputDir->GetName();

      std::vector<const NpCentCollision*> dfCollisions;
      for (const auto& coll : collisions) {
        if (coll.dfName == dfName) {
          dfCollisions.push_back(&coll);
        }
      }

      if (static_cast<Long64_t>(dfCollisions.size()) != inputTree->GetEntries()) {
        std::cerr << "Entry mismatch for " << dfName << "/O2npcollisiontabl: input has "
                  << inputTree->GetEntries() << " rows, computed centrality has "
                  << dfCollisions.size() << " rows" << std::endl;
        delete object;
        return;
      }

      const bool hadNoSameBunchPileup = inputTree->GetBranch("fNoSameBunchPileup");
      const bool hadCentFT0MNoPileup = inputTree->GetBranch("fCentFT0MNoPileup");
      inputTree->ResetBranchAddresses();
      inputTree->SetBranchStatus("*", 1);
      if (hadNoSameBunchPileup) {
        inputTree->SetBranchStatus("fNoSameBunchPileup", 0);
      }
      if (hadCentFT0MNoPileup) {
        inputTree->SetBranchStatus("fCentFT0MNoPileup", 0);
      }

      auto* outputTree = inputTree->CloneTree(0);
      outputTree->SetAutoSave(0);

      bool outNoSameBunchPileup = false;
      float outCentFT0MNoPileup = 101.0f;
      outputTree->Branch("fNoSameBunchPileup", &outNoSameBunchPileup, "fNoSameBunchPileup/O");
      outputTree->Branch("fCentFT0MNoPileup", &outCentFT0MNoPileup, "fCentFT0MNoPileup/F");

      for (Long64_t i = 0; i < inputTree->GetEntries(); ++i) {
        inputTree->GetEntry(i);
        outNoSameBunchPileup = dfCollisions[static_cast<std::size_t>(i)]->noSameBunchPileup;
        outCentFT0MNoPileup = dfCollisions[static_cast<std::size_t>(i)]->centFT0MNoPileup;
        outputTree->Fill();
      }
      outputTree->Write("", TObject::kOverwrite);
      inputTree->SetBranchStatus("*", 1);
      delete outputTree;
      delete object;
      continue;
    }

    if (object->InheritsFrom(TTree::Class())) {
      auto* inputTree = dynamic_cast<TTree*>(object);
      inputTree->ResetBranchAddresses();
      inputTree->SetBranchStatus("*", 1);
      auto* outputTree = inputTree->CloneTree(-1, "fast");
      if (!outputTree) {
        std::cerr << "Cannot copy tree " << objectName << " in " << inputPath << std::endl;
        delete object;
        return;
      }
      outputTree->Write(key->GetName(), TObject::kOverwrite);
      delete outputTree;
      delete object;
      continue;
    }

    object->Write(key->GetName(), TObject::kOverwrite);
    delete object;
  }
}

void makeNoPileupCentralityTable(TString inputAO2D = "AO2Ddata.root",
                                 TString outputFile = "NoPileupCentrality.root")
{
  if (inputAO2D == outputFile) {
    std::cerr << "Input and output file names are identical. Use a different output file." << std::endl;
    return;
  }

  TFile input(inputAO2D);
  if (input.IsZombie() || !input.GetListOfKeys()) {
    std::cerr << "Cannot open AO2D file: " << inputAO2D << std::endl;
    return;
  }

  std::vector<NpCentCollision> collisions;
  std::unordered_map<NpCentBCKey, int, NpCentBCKeyHash> collisionsPerBC;

  int dfIndex = 0;
  Long64_t nDF = 0;
  Long64_t nMissingTables = 0;
  bool useStoredNoSameBunchPileup = false;
  bool sawStoredNoSameBunchPileup = false;
  bool sawMissingNoSameBunchPileup = false;

  for (auto key : *input.GetListOfKeys()) {
    TString keyName = key->GetName();
    if (!keyName.Contains("DF_")) {
      continue;
    }
    ++nDF;

    auto* collTree = dynamic_cast<TTree*>(input.Get((keyName + "/O2npcollisiontabl").Data()));
    if (!collTree) {
      ++nMissingTables;
      ++dfIndex;
      continue;
    }

    if (!collTree->GetBranch("fRunNumber") ||
        !collTree->GetBranch("fGlobalBC") ||
        !collTree->GetBranch("fNumContrib") ||
        !collTree->GetBranch("fMultNTracksGlobal") ||
        !collTree->GetBranch("fMultNTracksNP") ||
        !collTree->GetBranch("fMultFT0M") ||
        !collTree->GetBranch("fCentFT0M")) {
      std::cerr << "NPCollisionTABLE in " << keyName
                << " is missing one of required branches: fRunNumber, fGlobalBC, "
                << "fNumContrib, fMultNTracksGlobal, fMultNTracksNP, "
                << "fMultFT0M, fCentFT0M" << std::endl;
      return;
    }

    int runNumber = 0;
    ULong64_t globalBC = 0;
    UShort_t numContrib = 0;
    int multNTracksGlobal = 0;
    int multNTracksNP = 0;
    float multFT0M = 0.0f;
    float centFT0M = 0.0f;
    bool noSameBunchPileup = false;

    collTree->SetBranchAddress("fRunNumber", &runNumber);
    collTree->SetBranchAddress("fGlobalBC", &globalBC);
    collTree->SetBranchAddress("fNumContrib", &numContrib);
    collTree->SetBranchAddress("fMultNTracksGlobal", &multNTracksGlobal);
    collTree->SetBranchAddress("fMultNTracksNP", &multNTracksNP);
    collTree->SetBranchAddress("fMultFT0M", &multFT0M);
    collTree->SetBranchAddress("fCentFT0M", &centFT0M);
    const bool hasStoredNoSameBunchPileup = collTree->GetBranch("fNoSameBunchPileup");
    if (hasStoredNoSameBunchPileup) {
      sawStoredNoSameBunchPileup = true;
      collTree->SetBranchAddress("fNoSameBunchPileup", &noSameBunchPileup);
    } else {
      sawMissingNoSameBunchPileup = true;
    }

    const Long64_t nEntries = collTree->GetEntries();
    collisions.reserve(collisions.size() + static_cast<std::size_t>(nEntries));

    for (Long64_t i = 0; i < nEntries; ++i) {
      collTree->GetEntry(i);

      NpCentCollision out;
      out.dfIndex = dfIndex;
      out.dfName = keyName;
      out.collisionIndex = i;
      out.runNumber = runNumber;
      out.globalBC = globalBC;
      out.numContrib = numContrib;
      out.multNTracksGlobal = multNTracksGlobal;
      out.multNTracksNP = multNTracksNP;
      out.multFT0M = multFT0M;
      out.centFT0M = centFT0M;
      out.noSameBunchPileup = noSameBunchPileup;
      collisions.push_back(out);

      ++collisionsPerBC[NpCentBCKey{runNumber, globalBC}];
    }

    ++dfIndex;
  }

  if (nDF <= 0) {
    std::cerr << "No DF_* directories found in " << inputAO2D << std::endl;
    return;
  }
  if (collisions.empty()) {
    std::cerr << "No NPCollisionTABLE rows found in " << inputAO2D << std::endl;
    return;
  }
  useStoredNoSameBunchPileup = sawStoredNoSameBunchPileup && !sawMissingNoSameBunchPileup;

  std::vector<std::pair<float, std::size_t>> noPileupMult;
  noPileupMult.reserve(collisions.size());
  Long64_t nRejectedSameBC = 0;

  for (std::size_t i = 0; i < collisions.size(); ++i) {
    auto& coll = collisions[i];
    if (!useStoredNoSameBunchPileup) {
      const auto countIt = collisionsPerBC.find(NpCentBCKey{coll.runNumber, coll.globalBC});
      coll.noSameBunchPileup = (countIt != collisionsPerBC.end() && countIt->second == 1);
    }
    if (coll.noSameBunchPileup) {
      noPileupMult.push_back(std::make_pair(coll.multFT0M, i));
    } else {
      ++nRejectedSameBC;
    }
  }

  std::sort(noPileupMult.begin(), noPileupMult.end(),
            [](const std::pair<float, std::size_t>& lhs,
               const std::pair<float, std::size_t>& rhs) {
              return lhs.first > rhs.first;
            });

  const double nNoPileup = static_cast<double>(noPileupMult.size());
  if (nNoPileup > 0.0) {
    std::size_t first = 0;
    while (first < noPileupMult.size()) {
      std::size_t last = first;
      while (last + 1 < noPileupMult.size() &&
             noPileupMult[last + 1].first == noPileupMult[first].first) {
        ++last;
      }

      // Use the upper percentile edge for tied FT0M values, so the most peripheral
      // tie group reaches 100 instead of sitting at its average rank.
      const float percentile = static_cast<float>(
        (static_cast<double>(last) + 1.0) * 100.0 / nNoPileup);
      for (std::size_t rank = first; rank <= last; ++rank) {
        collisions[noPileupMult[rank].second].centFT0MNoPileup = percentile;
      }
      first = last + 1;
    }
  }

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Cannot create output file: " << outputFile << std::endl;
    return;
  }

  copyDirectoryWithEnrichedCollisionTable(&input, &output, collisions);
  output.Close();

  std::cout << "Wrote " << outputFile << std::endl;
  std::cout << "  DF directories: " << nDF
            << " missing NPCollisionTABLE: " << nMissingTables << std::endl;
  std::cout << "  collisions read: " << collisions.size() << std::endl;
  std::cout << "  no-same-bunch-pileup collisions: " << noPileupMult.size() << std::endl;
  std::cout << "  same-bunch-pileup rejected collisions: " << nRejectedSameBC << std::endl;
  std::cout << "  pileup source: "
            << (useStoredNoSameBunchPileup
                  ? "stored fNoSameBunchPileup"
                  : "computed from raw fRunNumber/fGlobalBC")
            << std::endl;
  if (sawStoredNoSameBunchPileup && sawMissingNoSameBunchPileup) {
    std::cout << "  warning: fNoSameBunchPileup was present only in some DF trees; "
              << "used raw fRunNumber/fGlobalBC for consistency" << std::endl;
  }
  std::cout << "  output layout: full input AO2D copy with enriched DF_*/O2npcollisiontabl" << std::endl;
}
