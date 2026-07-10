#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TTree.h"

constexpr int NPT   = 10;
constexpr int NMULT = 80;

struct Matrix2D_t {
  double data[NPT][NMULT];
  Matrix2D_t() { std::fill(&data[0][0], &data[0][0] + NPT * NMULT, 0.0); }

  Matrix2D_t& operator=(const Matrix2D_t& other) {
    if (this != &other) {
      std::copy(&other.data[0][0],
                &other.data[0][0] + NPT * NMULT,
                &data[0][0]);
    }
    return *this;
  }

  double& operator()(int i, int j) { return data[i][j]; }
  const double& operator()(int i, int j) const { return data[i][j]; }

  double* operator[](int i) { return data[i]; }
  const double* operator[](int i) const { return data[i]; }
};

struct TrackInfo {
  double ptTrue;
  double ptReco;
  bool isReco;
};

using Matrix4D_t = double[NPT][NMULT][NPT][NMULT];

void printMatrix2D(const Matrix2D_t& m, const std::string& name = "")
{
  std::cout << "===> " << name << std::endl;
  for (int i = 0; i < NPT; ++i) {
    for (int j = 0; j < NMULT; ++j) {
      std::cout << m[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
Matrix2D_t poissonRandomizeMatrix(const Matrix2D_t& input, TRandom3& rand)
{
  Matrix2D_t randomized;
  for (int ix = 0; ix < NPT; ++ix) {
    for (int iy = 0; iy < NMULT; ++iy) {
      double mean = input[ix][iy];
      randomized[ix][iy] = mean > 0.0 ? rand.PoissonD(mean) : 0.0;
    }
  }
  return randomized;
}

Matrix2D_t makeFlatPrior(double value = 1.0)
{
  Matrix2D_t prior;
  for (int ix = 0; ix < NPT; ++ix) {
    for (int iy = 0; iy < NMULT; ++iy) {
      prior[ix][iy] = value;
    }
  }
  return prior;
}
void printMatrix4D(Matrix4D_t& m, const std::string& name = "")
{
  std::cout << "===> " << name << std::endl;
  for (int jt = 0; jt < NMULT; ++jt) {
    for (int it = 0; it < NPT; ++it) {
      std::cout << "ptTrue:" << it << " multTrue:" << jt << std::endl;
      for (int jr = 0; jr < NMULT; ++jr) {
        for (int ir = 0; ir < NPT; ++ir) {
          std::cout << m[ir][jr][it][jt] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  std::cout << std::endl;
}

double sigma_pt(double pt)
{
  return pt * std::sqrt((0.02 / pt) * (0.02 / pt) + 0.01 * 0.01);
}

void copyTH2ToMatrix2D(const TH2D* h, Matrix2D_t& m)
{
  for (int ix = 0; ix < NPT; ++ix) {
    for (int iy = 0; iy < NMULT; ++iy) {
      m[ix][iy] = h->GetBinContent(ix + 1, iy + 1);
    }
  }
}

double sumMatrix2D(const Matrix2D_t& m)
{
  double sum = 0.0;
  for (int ix = 0; ix < NPT; ++ix) {
    for (int iy = 0; iy < NMULT; ++iy) {
      sum += m[ix][iy];
    }
  }
  return sum;
}
//=====
void fillChainFromAO2D(TChain &chain, const TString &fileName)
{
  TFile file(fileName);
  if (file.IsZombie() || !file.GetListOfKeys()) {
    std::cerr << "Cannot open AO2D file: " << fileName << std::endl;
    return;
  }

  for (auto key : *file.GetListOfKeys())
  {
    TString keyName = key->GetName();
    if (keyName.Contains("DF_"))
    {
      chain.Add((fileName + "/" + keyName + "/" + chain.GetName()).Data());
    }
  }
}

void generateTrainingSample(int nTrainEvents,
                            TF1* fpt,
                            std::negative_binomial_distribution<int>& nbd,
                            std::mt19937& genTrain,
                            TRandom3& randTrain,
                            double eff,
                            double fakeMean,
                            double ptMin,
                            double ptMax,
                            double multMin,
                            double multMax,
                            TH2D* hTruthTrainTmp,
                            TH2D* hRecoTrainTmp,
                            TH2D* hMissTmp,
                            TH1D* hMultTrueTrain,
                            TH1D* hMultRecoTrain,
                            TH1D* hPtTrueTrain,
                            TH1D* hPtRecoTrain,
                            Matrix4D_t& Rcount,
                            Matrix2D_t& FakeTrain,
                            Matrix2D_t& RecoAllTrain)
{
  for (int ev = 0; ev < nTrainEvents; ++ev) {
    int multTrue = nbd(genTrain);
    if (multTrue < multMin || multTrue >= multMax) {
      continue;
    }

    hMultTrueTrain->Fill(multTrue);
    int jt = hMultTrueTrain->GetXaxis()->FindBin(multTrue) - 1;

    std::vector<TrackInfo> tracks;
    tracks.reserve(multTrue);

    int multReco = 0;

    for (int i = 0; i < multTrue; ++i) {
      //double ptTrue = randTrain.Rndm() * 10.0;
      double ptTrue = fpt->GetRandom();
      bool isReco   = (randTrain.Rndm() < eff);
      double ptReco = -1.0;

      if (ptTrue >= ptMin && ptTrue < ptMax) {
        hPtTrueTrain->Fill(ptTrue);
      }

      if (isReco) {
        ptReco = randTrain.Gaus(ptTrue, sigma_pt(ptTrue));
        //ptReco = ptTrue; // debug
        if (ptReco > ptMin && ptReco < ptMax) {
          ++multReco;
        } else {
          isReco = false;
          ptReco = -1.0;
        }
      }

      tracks.push_back({ptTrue, ptReco, isReco});
    }

    int nFake = randTrain.Poisson(fakeMean);
    multReco += nFake;

    if (multReco < multMin || multReco >= multMax) {
      continue;
    }

    hMultRecoTrain->Fill(multReco);
    int jr = hMultRecoTrain->GetXaxis()->FindBin(multReco) - 1;

    for (const auto& trk : tracks) {
      if (trk.ptTrue < ptMin || trk.ptTrue >= ptMax) {
        continue;
      }

      hTruthTrainTmp->Fill(trk.ptTrue, multTrue);
      int it = hTruthTrainTmp->GetXaxis()->FindBin(trk.ptTrue) - 1;

      if (trk.isReco) {
        if (trk.ptReco > ptMin && trk.ptReco < ptMax) {
          hRecoTrainTmp->Fill(trk.ptReco, multReco);
          hPtRecoTrain->Fill(trk.ptReco);

          int ir = hRecoTrainTmp->GetXaxis()->FindBin(trk.ptReco) - 1;
          Rcount[ir][jr][it][jt] += 1.0;
          RecoAllTrain[ir][jr] += 1.0;
        } else {
          hMissTmp->Fill(trk.ptTrue, multTrue);
        }
      } else {
        hMissTmp->Fill(trk.ptTrue, multTrue);
      }
    }

    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();
      if (ptFake >= ptMin && ptFake < ptMax) {
        hRecoTrainTmp->Fill(ptFake, multReco);
        hPtRecoTrain->Fill(ptFake);

        int ir = hRecoTrainTmp->GetXaxis()->FindBin(ptFake) - 1;
        FakeTrain[ir][jr] += 1.0;
        RecoAllTrain[ir][jr] += 1.0;
      }
    }
  }
}

bool readNonPromptCascadeTrainingFromAO2D(const TString& trainAO2D,
                                           double ptMin,
                                           double ptMax,
                                           double multMin,
                                           double multMax,
                                           TH2D* hTruthTrainTmp,
                                           TH2D* hRecoTrainTmp,
                                           TH2D* hMissTmp,
                                           TH1D* hMultTrueTrain,
                                           TH1D* hMultRecoTrain,
                                           TH1D* hPtTrueTrain,
                                           TH1D* hPtRecoTrain,
                                           Matrix4D_t& Rcount,
                                           Matrix2D_t& FakeTrain,
                                           Matrix2D_t& FeedInTrain,
                                           Matrix2D_t& RecoAllTrain)
{
  TChain trainChain("O2npmcchargedtabl");
  fillChainFromAO2D(trainChain, trainAO2D);

  if (trainChain.GetEntries() <= 0) {
    std::cerr << "No nonPromptCascade NPMCChargedTABLE entries found in "
              << trainAO2D << std::endl;
    return false;
  }

  if (!trainChain.GetBranch("fPtGen") ||
      !trainChain.GetBranch("fPtRec") ||
      !trainChain.GetBranch("fMultGen") ||
      !trainChain.GetBranch("fMultNTracksNP")) {
    std::cerr << "NPMCChargedTABLE is missing one of the required branches: "
              << "fPtGen, fPtRec, fMultGen, fMultNTracksNP" << std::endl;
    return false;
  }

  float ptTrue = 0.0f;
  float ptReco = 0.0f;
  int multTrue = 0;
  int multReco = 0;

  trainChain.SetBranchAddress("fPtGen", &ptTrue);
  trainChain.SetBranchAddress("fPtRec", &ptReco);
  trainChain.SetBranchAddress("fMultGen", &multTrue);
  trainChain.SetBranchAddress("fMultNTracksNP", &multReco);

  Long64_t nMatched = 0;
  Long64_t nMissed = 0;
  Long64_t nFake = 0;
  Long64_t nFeedIn = 0;
  Long64_t nIgnored = 0;

  for (Long64_t i = 0; i < trainChain.GetEntries(); ++i) {
    trainChain.GetEntry(i);

    // Sentinel convention from nonPromptCascade NPMCChargedTABLE:
    //   matched: ptGen>=0, ptRec>=0; fake: ptGen=-1/-2/-3; feed-in: ptGen=-4; missed: ptRec=-1/-2.
    const bool isFeedIn = (ptTrue <= -3.5f);
    const bool hasTruth = (ptTrue >= ptMin && ptTrue < ptMax &&
                           multTrue >= multMin && multTrue < multMax);
    const bool hasReco = (ptReco > ptMin && ptReco < ptMax &&
                          multReco >= multMin && multReco < multMax);

    if (hasTruth) {
      hTruthTrainTmp->Fill(ptTrue, multTrue);
      hPtTrueTrain->Fill(ptTrue);
      hMultTrueTrain->Fill(multTrue);
    }

    if (hasReco) {
      hRecoTrainTmp->Fill(ptReco, multReco);
      hPtRecoTrain->Fill(ptReco);
      hMultRecoTrain->Fill(multReco);
    }

    if (hasTruth && hasReco) {
      // Matched row: fills the detector response R(reco pt,mult | truth pt,mult).
      int it = hTruthTrainTmp->GetXaxis()->FindBin(ptTrue) - 1;
      int jt = hTruthTrainTmp->GetYaxis()->FindBin(multTrue) - 1;
      int ir = hRecoTrainTmp->GetXaxis()->FindBin(ptReco) - 1;
      int jr = hRecoTrainTmp->GetYaxis()->FindBin(multReco) - 1;
      Rcount[ir][jr][it][jt] += 1.0;
      RecoAllTrain[ir][jr] += 1.0;
      ++nMatched;
    } else if (hasTruth) {
      // Missed row: truth particle has no usable reco entry; enters response normalization.
      hMissTmp->Fill(ptTrue, multTrue);
      ++nMissed;
    } else if (hasReco) {
      int ir = hRecoTrainTmp->GetXaxis()->FindBin(ptReco) - 1;
      int jr = hRecoTrainTmp->GetYaxis()->FindBin(multReco) - 1;
      if (isFeedIn) {
        // Feed-in row: reco track from outside truth fiducial volume, ptGen=-4; subtracted separately.
        FeedInTrain[ir][jr] += 1.0;
        ++nFeedIn;
      } else {
        // Fake row: reco track without a usable matched truth particle, ptGen=-1/-2/-3.
        FakeTrain[ir][jr] += 1.0;
        ++nFake;
      }
      // All reconstructed rows provide the denominator for non-signal subtraction.
      RecoAllTrain[ir][jr] += 1.0;
    } else {
      ++nIgnored;
    }
  }

  std::cout << "Read nonPromptCascade NPMCChargedTABLE from " << trainAO2D << std::endl;
  std::cout << "  matched: " << nMatched
            << " missed: " << nMissed
            << " fake: " << nFake
            << " feed-in: " << nFeedIn
            << " ignored: " << nIgnored << std::endl;

  return true;
}

void generatePseudoData(int nDataEvents,
                        TF1* fpt,
                        std::negative_binomial_distribution<int>& nbd,
                        std::mt19937& genData,
                        TRandom3& randData,
                        double eff,
                        double fakeMean,
                        double ptMin,
                        double ptMax,
                        double multMin,
                        double multMax,
                        TH2D* hTruthDataTmp,
                        TH2D* hRecoDataTmp,
                        TH1D* hMultTrueData,
                        TH1D* hMultRecoData,
                        TH1D* hPtTrueData,
                        TH1D* hPtRecoData)
{
  for (int ev = 0; ev < nDataEvents; ++ev) {
    int multTrue = nbd(genData);
    if (multTrue < multMin || multTrue >= multMax) {
      continue;
    }

    hMultTrueData->Fill(multTrue);

    std::vector<TrackInfo> tracks;
    tracks.reserve(multTrue);

    int multReco = 0;

    for (int i = 0; i < multTrue; ++i) {
      double ptTrue = fpt->GetRandom();
      bool isReco   = (randData.Rndm() < eff);
      double ptReco = -1.0;

      if (ptTrue >= ptMin && ptTrue < ptMax) {
        hPtTrueData->Fill(ptTrue);
      }

      if (isReco) {
        ptReco = randData.Gaus(ptTrue, sigma_pt(ptTrue));
        if (ptReco > ptMin && ptReco < ptMax) {
          ++multReco;
        } else {
          isReco = false;
          ptReco = -1.0;
        }
      }

      tracks.push_back({ptTrue, ptReco, isReco});
    }

    int nFake = randData.Poisson(fakeMean);
    multReco += nFake;

    if (multReco < multMin || multReco >= multMax) {
      continue;
    }

    hMultRecoData->Fill(multReco);

    for (const auto& trk : tracks) {
      if (trk.ptTrue >= ptMin && trk.ptTrue < ptMax) {
        hTruthDataTmp->Fill(trk.ptTrue, multTrue);
      }

      if (trk.isReco && trk.ptReco > ptMin && trk.ptReco < ptMax) {
        hRecoDataTmp->Fill(trk.ptReco, multReco);
        hPtRecoData->Fill(trk.ptReco);
      }
    }

    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();
      if (ptFake >= ptMin && ptFake < ptMax) {
        hRecoDataTmp->Fill(ptFake, multReco);
        hPtRecoData->Fill(ptFake);
      }
    }
  }
}

bool readRecoDataFromAO2D(const TString& dataAO2D,
                          double ptMin,
                          double ptMax,
                          double multMin,
                          double multMax,
                          double centMin,
                          double centMax,
                          int sameBCPolicy,
                          TH2D* hRecoDataTmp,
                          TH1D* hMultRecoData,
                          TH1D* hMultRecoDataCand,
                          TH1D* hMultRecoDataCollWithCand,
                          TH1D* hPtRecoData)
{
  const auto totalTimerStart = std::chrono::steady_clock::now();
  auto elapsedSeconds = [](const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  };

  TFile file(dataAO2D);
  if (file.IsZombie() || !file.GetListOfKeys()) {
    std::cerr << "Cannot open AO2D data file: " << dataAO2D << std::endl;
    return false;
  }

  Long64_t nCollisions = 0;
  Long64_t nFilled = 0;
  Long64_t nIgnored = 0;
  Long64_t nMissingTables = 0;
  Long64_t nDF = 0;
  Long64_t nCentralityRejected = 0;
  Long64_t nSameBCRejected = 0;
  double timeCollisionRead = 0.0;
  double timeCollisionSelect = 0.0;
  double timeCandidateRead = 0.0;
  double timeCollWithCand = 0.0;

  for (auto key : *file.GetListOfKeys()) {
    TString keyName = key->GetName();
    if (!keyName.Contains("DF_")) {
      continue;
    }
    ++nDF;

    auto* collTree = dynamic_cast<TTree*>(file.Get((keyName + "/O2npcollisiontabl").Data()));
    auto* candTree = dynamic_cast<TTree*>(file.Get((keyName + "/O2nprecochargedca").Data()));
    if (!collTree || !candTree) {
      ++nMissingTables;
      continue;
    }

    if (!collTree->GetBranch("fMultNTracksNP") || !collTree->GetBranch("fCentFT0M")) {
      std::cerr << "NPCollisionTABLE in " << keyName
                << " is missing required branch fMultNTracksNP and/or fCentFT0M" << std::endl;
      return false;
    }
    if (sameBCPolicy != 0 &&
        (!collTree->GetBranch("fNoSameBunchPileup") ||
         !collTree->GetBranch("fCentFT0MNoPileup"))) {
      std::cerr << "NPCollisionTABLE in " << keyName
                << " is missing required branch fNoSameBunchPileup and/or fCentFT0MNoPileup for sameBCPolicy != 0" << std::endl;
      return false;
    }
    if (!candTree->GetBranch("fPtRec") || !candTree->GetBranch("fIndexNPCollisionTable")) {
      std::cerr << "NPRecoChargedCand in " << keyName
                << " is missing required branches fPtRec and/or fIndexNPCollisionTable" << std::endl;
      return false;
    }

    int multReco = 0;
    float centFT0M = 0.0f;
    float centFT0MNoPileup = 101.0f;
    bool noSameBunchPileup = true;
    const Long64_t nCollEntries = collTree->GetEntries();
    std::vector<int> recoMultByCollision(nCollEntries, -1);
    std::vector<float> centFT0MByCollision(nCollEntries, 0.0f);
    std::vector<float> centFT0MNoPileupByCollision;
    std::vector<char> noSameBunchPileupByCollision;
    std::vector<char> acceptedCollision(nCollEntries, 0);
    std::vector<char> hasSelectedCandidate(nCollEntries, 0);
    collTree->SetBranchAddress("fMultNTracksNP", &multReco);
    collTree->SetBranchAddress("fCentFT0M", &centFT0M);

    if (sameBCPolicy != 0) {
      centFT0MNoPileupByCollision.resize(nCollEntries, 101.0f);
      noSameBunchPileupByCollision.resize(nCollEntries, 0);
      collTree->SetBranchAddress("fCentFT0MNoPileup", &centFT0MNoPileup);
      collTree->SetBranchAddress("fNoSameBunchPileup", &noSameBunchPileup);
    }

    auto timerStart = std::chrono::steady_clock::now();
    for (Long64_t i = 0; i < nCollEntries; ++i) {
      collTree->GetEntry(i);
      recoMultByCollision[i] = multReco;
      centFT0MByCollision[i] = centFT0M;

      if (sameBCPolicy != 0) {
        centFT0MNoPileupByCollision[i] = centFT0MNoPileup;
        noSameBunchPileupByCollision[i] = noSameBunchPileup ? 1 : 0;
      }
    }
    timeCollisionRead += elapsedSeconds(timerStart);

    timerStart = std::chrono::steady_clock::now();
    for (Long64_t i = 0; i < nCollEntries; ++i) {
      multReco = recoMultByCollision[i];
      const bool acceptedBC = (sameBCPolicy == 0) ? true : noSameBunchPileupByCollision[i];
      const float centralityForCut = (sameBCPolicy == 0)
        ? centFT0MByCollision[i]
        : centFT0MNoPileupByCollision[i];
      const bool inCentrality = (centralityForCut >= centMin &&
                                 centralityForCut < centMax);
      acceptedCollision[i] = (inCentrality && acceptedBC) ? 1 : 0;
      if (!acceptedBC) {
        ++nSameBCRejected;
        continue;
      }
      if (!inCentrality) {
        ++nCentralityRejected;
        continue;
      }
      if (multReco >= multMin && multReco < multMax) {
        hMultRecoData->Fill(multReco);
      }
      ++nCollisions;
    }
    timeCollisionSelect += elapsedSeconds(timerStart);

    float ptReco = 0.0f;
    int npCollisionId = -1;
    candTree->SetBranchAddress("fPtRec", &ptReco);
    candTree->SetBranchAddress("fIndexNPCollisionTable", &npCollisionId);

    timerStart = std::chrono::steady_clock::now();
    for (Long64_t i = 0; i < candTree->GetEntries(); ++i) {
      candTree->GetEntry(i);
      if (npCollisionId < 0 || npCollisionId >= static_cast<int>(recoMultByCollision.size())) {
        ++nIgnored;
        continue;
      }
      if (!acceptedCollision[npCollisionId]) {
        ++nIgnored;
        continue;
      }
      multReco = recoMultByCollision[npCollisionId];
      if (ptReco > ptMin && ptReco < ptMax &&
          multReco >= multMin && multReco < multMax) {
        hRecoDataTmp->Fill(ptReco, multReco);
        hMultRecoDataCand->Fill(multReco);
        hasSelectedCandidate[npCollisionId] = 1;
        hPtRecoData->Fill(ptReco);
        ++nFilled;
      } else {
        ++nIgnored;
      }
    }
    timeCandidateRead += elapsedSeconds(timerStart);

    timerStart = std::chrono::steady_clock::now();
    for (Long64_t i = 0; i < nCollEntries; ++i) {
      if (!acceptedCollision[i]) {
        continue;
      }
      multReco = recoMultByCollision[i];
      if (hasSelectedCandidate[i] && multReco >= multMin && multReco < multMax) {
        hMultRecoDataCollWithCand->Fill(multReco);
      }
    }
    timeCollWithCand += elapsedSeconds(timerStart);
  }

  if (nDF <= 0) {
    std::cerr << "No DF_* directories found in " << dataAO2D << std::endl;
    return false;
  }
  if (nCollisions <= 0) {
    std::cerr << "No nonPromptCascade NPCollisionTABLE entries found in "
              << dataAO2D << std::endl;
    return false;
  }
  if (nFilled <= 0) {
    std::cerr << "No nonPromptCascade NPRecoChargedCand entries filled from "
              << dataAO2D << std::endl;
    return false;
  }

  std::cout << "Read nonPromptCascade reco data from " << dataAO2D << std::endl;
  std::cout << "  DF directories: " << nDF
            << " missing tables: " << nMissingTables
            << " centrality-accepted collisions read: " << nCollisions << std::endl;
  std::cout << "  centrality FT0M accepted range: [" << centMin << ", " << centMax
            << ") rejected collisions: " << nCentralityRejected << std::endl;
  std::cout << "  pileup policy: " << sameBCPolicy
            << (sameBCPolicy == 0 ? " (none)" : " (stored fNoSameBunchPileup)")
            << " rejected collisions: " << nSameBCRejected << std::endl;
  std::cout << "  reco candidates filled: " << nFilled
            << " ignored: " << nIgnored << std::endl;
  std::cout << "  timing seconds: collisionRead=" << timeCollisionRead
            << " collisionSelect=" << timeCollisionSelect
            << " candidateRead=" << timeCandidateRead
            << " collWithCandidate=" << timeCollWithCand
            << " total=" << elapsedSeconds(totalTimerStart) << std::endl;
  return true;
}

bool readTruthDataFromMCAO2D(const TString& truthAO2D,
                             double ptMin,
                             double ptMax,
                             double multMin,
                             double multMax,
                             TH2D* hTruthDataTmp,
                             TH1D* hMultTrueData,
                             TH1D* hPtTrueData)
{
  TChain truthChain("O2npmcchargedtabl");
  fillChainFromAO2D(truthChain, truthAO2D);

  if (truthChain.GetEntries() <= 0) {
    std::cerr << "No nonPromptCascade NPMCChargedTABLE entries found in "
              << truthAO2D << std::endl;
    return false;
  }

  if (!truthChain.GetBranch("fPtGen") || !truthChain.GetBranch("fMultGen")) {
    std::cerr << "NPMCChargedTABLE is missing required branches fPtGen and/or fMultGen" << std::endl;
    return false;
  }

  float ptTrue = 0.0f;
  int multTrue = 0;
  Long64_t nTruth = 0;
  Long64_t nIgnored = 0;

  truthChain.SetBranchAddress("fPtGen", &ptTrue);
  truthChain.SetBranchAddress("fMultGen", &multTrue);

  for (Long64_t i = 0; i < truthChain.GetEntries(); ++i) {
    truthChain.GetEntry(i);
    if (ptTrue >= ptMin && ptTrue < ptMax &&
        multTrue >= multMin && multTrue < multMax) {
      hTruthDataTmp->Fill(ptTrue, multTrue);
      hPtTrueData->Fill(ptTrue);
      hMultTrueData->Fill(multTrue);
      ++nTruth;
    } else {
      ++nIgnored;
    }
  }

  std::cout << "Read closure truth from " << truthAO2D << std::endl;
  std::cout << "  truth particles filled: " << nTruth
            << " ignored: " << nIgnored << std::endl;
  return true;
}


bool readRecoAndTruthDataFromMCAO2D(const TString& dataAO2D,
                                    double ptMin,
                                    double ptMax,
                                    double multMin,
                                    double multMax,
                                    TH2D* hTruthDataTmp,
                                    TH2D* hRecoDataTmp,
                                    TH1D* hMultTrueData,
                                    TH1D* hMultRecoDataCand,
                                    TH1D* hPtTrueData,
                                    TH1D* hPtRecoData)
{
  TChain dataChain("O2npmcchargedtabl");
  fillChainFromAO2D(dataChain, dataAO2D);

  if (dataChain.GetEntries() <= 0) {
    std::cerr << "No nonPromptCascade NPMCChargedTABLE entries found in "
              << dataAO2D << std::endl;
    return false;
  }

  if (!dataChain.GetBranch("fPtGen") ||
      !dataChain.GetBranch("fPtRec") ||
      !dataChain.GetBranch("fMultGen") ||
      !dataChain.GetBranch("fMultNTracksNP")) {
    std::cerr << "NPMCChargedTABLE is missing one of the required branches: "
              << "fPtGen, fPtRec, fMultGen, fMultNTracksNP" << std::endl;
    return false;
  }

  float ptTrue = 0.0f;
  float ptReco = 0.0f;
  int multTrue = 0;
  int multReco = 0;
  Long64_t nTruth = 0;
  Long64_t nReco = 0;
  Long64_t nIgnored = 0;

  dataChain.SetBranchAddress("fPtGen", &ptTrue);
  dataChain.SetBranchAddress("fPtRec", &ptReco);
  dataChain.SetBranchAddress("fMultGen", &multTrue);
  dataChain.SetBranchAddress("fMultNTracksNP", &multReco);

  for (Long64_t i = 0; i < dataChain.GetEntries(); ++i) {
    dataChain.GetEntry(i);

    const bool hasTruth = (ptTrue >= ptMin && ptTrue < ptMax &&
                           multTrue >= multMin && multTrue < multMax);
    const bool hasReco = (ptReco > ptMin && ptReco < ptMax &&
                          multReco >= multMin && multReco < multMax);

    if (hasTruth) {
      // Truth closure target from the same MC sample.
      hTruthDataTmp->Fill(ptTrue, multTrue);
      hPtTrueData->Fill(ptTrue);
      hMultTrueData->Fill(multTrue);
      ++nTruth;
    }

    if (hasReco) {
      // Reco pseudo-data from the MC table. This includes matched reco rows plus fake/feed-in rows.
      hRecoDataTmp->Fill(ptReco, multReco);
      hPtRecoData->Fill(ptReco);
      hMultRecoDataCand->Fill(multReco);
      ++nReco;
    }

    if (!hasTruth && !hasReco) {
      ++nIgnored;
    }
  }

  std::cout << "Read reco and truth data from MC NPMCChargedTABLE in " << dataAO2D << std::endl;
  std::cout << "  truth particles filled: " << nTruth
            << " reco rows filled: " << nReco
            << " ignored rows: " << nIgnored << std::endl;
  std::cout << "  Note: hMultRecoData/hMultRecoDataCollWithCand are not filled in this mode; "
            << "O2npmcchargedtabl is row-based, not a collision table." << std::endl;
  return true;
}

// useAO2DTraining: true = train response from MC AO2D NPMCChargedTABLE; false = generate toy training.
// useAO2Data: 1 = read real reco data from dataAO2D; 0 = generate independent toy data;
//             2 = use RecoTrain and TruthTrain closure input;
//             3 = read reco from dataAO2D and truth from truthAO2D;
//             4 = read reco and truth from dataAO2D O2npmcchargedtabl.
// priorMode: 0 = flat positive prior; 1 = TruthTrain fixed-point check;
//            2 = RecoTrain-shaped prior; 3 = Poisson-randomized TruthTrain prior;
//            4 = Poisson-randomized flat positive prior.
// subtractFakeFeed: true = subtract fake/feed-in from reco input; false = bypass them.
// centMin/centMax: FT0M centrality percentile range for data-like reco collision tables.
// sameBCPolicy: 0 = use fCentFT0M and no pileup cut;
//               1 = require fNoSameBunchPileup and use fCentFT0MNoPileup.
void matrixH(TString trainAO2D = "AO2D.root", bool useAO2DTraining = false,
             TString dataAO2D = "AO2Ddata.root", TString truthAO2D = "AO2Dtruth.root",
             int useAO2Data = 2, int priorMode = 0, bool subtractFakeFeed = true,
             double centMin = 0.0, double centMax = 100.0, int sameBCPolicy = 0)
{
  if (sameBCPolicy < 0 || sameBCPolicy > 1) {
    std::cerr << "Invalid sameBCPolicy = " << sameBCPolicy
              << ". Use 0 = fCentFT0M without pileup cut, "
              << "or 1 = fNoSameBunchPileup with fCentFT0MNoPileup." << std::endl;
    return;
  }

  const double ptMin   = 0.0;
  const double ptMax   = 10.0;
  const double multMin = 0.0;
  const double multMax = static_cast<double>(NMULT); // keep integer multiplicity bins unit-width

  const int nTrainEvents = 1000000;
  const int nDataEvents  = 1000000;
  const int nIter        = 10;

  const double eff      = 0.8;
  const double fakeMean = 3.;

  TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", ptMin, ptMax);
  fpt->SetParameter(0, 0.5);

  double mu = 12.0;
  int r = 20;
  double p = r / (r + mu);

  std::mt19937 genTrain(12345);
  std::mt19937 genData(67890);
  std::negative_binomial_distribution<int> nbd(r, p);

  TRandom3 randTrain(1);
  TRandom3 randData(2);

  // ------------------------------------------------------------------
  // Histograms used as the primary containers during generation
  // ------------------------------------------------------------------
  TH2D* hTruthTrainTmp = new TH2D("hTruthTrainTmp", "train truth", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoTrainTmp  = new TH2D("hRecoTrainTmp",  "train reco",  NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hTruthDataTmp  = new TH2D("hTruthDataTmp",  "data truth",  NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoDataTmp   = new TH2D("hRecoDataTmp",   "data reco",   NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hMissTmp       = new TH2D("hMissTmp",       "miss",        NPT, ptMin, ptMax, NMULT, multMin, multMax);

  TH1D* hMultTrueTrain = new TH1D("hMultTrueTrain", "train mult true", NMULT, multMin, multMax);
  TH1D* hMultRecoTrain = new TH1D("hMultRecoTrain", "train mult reco", NMULT, multMin, multMax);
  TH1D* hMultTrueData  = new TH1D("hMultTrueData",  "data mult true",  NMULT, multMin, multMax);
  TH1D* hMultRecoData  = new TH1D("hMultRecoData",  "data mult reco all collisions",  NMULT, multMin, multMax);
  TH1D* hMultRecoDataCand = new TH1D("hMultRecoDataCand", "data mult reco selected candidates", NMULT, multMin, multMax);
  TH1D* hMultRecoDataCollWithCand = new TH1D("hMultRecoDataCollWithCand", "data mult reco collisions with selected candidates", NMULT, multMin, multMax);

  TH1D* hPtTrueTrain = new TH1D("hPtTrueTrain", "train pt true", NPT, ptMin, ptMax);
  TH1D* hPtRecoTrain = new TH1D("hPtRecoTrain", "train pt reco", NPT, ptMin, ptMax);
  TH1D* hPtTrueData  = new TH1D("hPtTrueData",  "data pt true",  NPT, ptMin, ptMax);
  TH1D* hPtRecoData  = new TH1D("hPtRecoData",  "data pt reco",  NPT, ptMin, ptMax);

  // Response tensor counts. Static storage keeps the large 4D array off the stack.
  static Matrix4D_t Rcount;
  std::fill(&Rcount[0][0][0][0], &Rcount[0][0][0][0] + NPT * NMULT * NPT * NMULT, 0.0);
  Matrix2D_t FakeTrain;
  Matrix2D_t FeedInTrain;
  Matrix2D_t RecoAllTrain;

  // ------------------------------------------------------------------
  // 1. TRAINING SAMPLE
  // ------------------------------------------------------------------
  if (useAO2DTraining) {
    std::cout << "Reading training sample from " << trainAO2D << std::endl;
    if (trainAO2D.IsNull()) {
      std::cerr << "AO2D training requested, but trainAO2D is empty" << std::endl;
      return;
    }
    if (!readNonPromptCascadeTrainingFromAO2D(trainAO2D, ptMin, ptMax, multMin, multMax,
                                    hTruthTrainTmp, hRecoTrainTmp, hMissTmp,
                                    hMultTrueTrain, hMultRecoTrain,
                                    hPtTrueTrain, hPtRecoTrain, Rcount,
                                    FakeTrain, FeedInTrain, RecoAllTrain)) {
      return;
    }
  } else {
    std::cout << "Generating training sample with " << nTrainEvents << " events" << std::endl;
    generateTrainingSample(nTrainEvents, fpt, nbd, genTrain, randTrain, eff,
                           fakeMean, ptMin, ptMax, multMin, multMax,
                           hTruthTrainTmp, hRecoTrainTmp, hMissTmp,
                           hMultTrueTrain, hMultRecoTrain, hPtTrueTrain,
                           hPtRecoTrain, Rcount, FakeTrain, RecoAllTrain);
  }

  // ------------------------------------------------------------------
  // 2. COPY TRAINING HISTOGRAMS TO MATRICES
  // ------------------------------------------------------------------
  Matrix2D_t TruthTrain;
  Matrix2D_t RecoTrain;
  Matrix2D_t Miss;

  copyTH2ToMatrix2D(hTruthTrainTmp, TruthTrain);
  copyTH2ToMatrix2D(hRecoTrainTmp,  RecoTrain);
  copyTH2ToMatrix2D(hMissTmp,       Miss);

  printMatrix2D(TruthTrain, "TruthTrain");
  printMatrix2D(RecoTrain,  "RecoTrain");
  printMatrix2D(Miss,       "Miss");

  // ------------------------------------------------------------------
  // 3. NORMALIZE RESPONSE
  // ------------------------------------------------------------------
  static Matrix4D_t Rprob;
  std::fill(&Rprob[0][0][0][0], &Rprob[0][0][0][0] + NPT * NMULT * NPT * NMULT, 0.0);
  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      double norm = Miss[it][jt];
      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          norm += Rcount[ir][jr][it][jt];
        }
      }

      if (norm <= 0.0) {
        std::cout << "norm <= 0 for true bin " << it << ", " << jt << std::endl;
        continue;
      }

      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          Rprob[ir][jr][it][jt] = Rcount[ir][jr][it][jt] / norm;
        }
      }
    }
  }

  //printMatrix4D(Rprob, "Rprob");

  // ------------------------------------------------------------------
  // 4. DATA SAMPLE
  // useAO2Data = 1: read real reco data from dataAO2D
  //             0: generate independent pseudo-data
  //             2: use RecoTrain and TruthTrain for closure/debug
  //             3: read reco from dataAO2D and closure truth from truthAO2D
  //             4: read reco and truth from dataAO2D O2npmcchargedtabl
  // ------------------------------------------------------------------
  if (useAO2Data == 1 || useAO2Data == 3) {
    std::cout << "Reading reco data sample from " << dataAO2D << std::endl;
    if (dataAO2D.IsNull()) {
      std::cerr << "AO2D data requested, but dataAO2D is empty" << std::endl;
      return;
    }
    if (!readRecoDataFromAO2D(dataAO2D, ptMin, ptMax, multMin, multMax,
                              centMin, centMax, sameBCPolicy,
                              hRecoDataTmp, hMultRecoData,
                              hMultRecoDataCand, hMultRecoDataCollWithCand,
                              hPtRecoData)) {
      return;
    }
    if (useAO2Data == 3) {
      std::cout << "Reading closure truth sample from " << truthAO2D << std::endl;
      if (truthAO2D.IsNull()) {
        std::cerr << "Closure truth requested, but truthAO2D is empty" << std::endl;
        return;
      }
      if (!readTruthDataFromMCAO2D(truthAO2D, ptMin, ptMax, multMin, multMax,
                                   hTruthDataTmp, hMultTrueData, hPtTrueData)) {
        return;
      }
    }
  } else if (useAO2Data == 4) {
    if (centMin > 0.0 || centMax < 100.0) {
      std::cout << "Centrality cut requested, but useAO2Data=4 reads O2npmcchargedtabl, "
                << "which has no centrality column; centrality cut is not applied." << std::endl;
    }
    if (sameBCPolicy != 0) {
      std::cout << "sameBCPolicy requested, but useAO2Data=4 reads O2npmcchargedtabl, "
                << "which has no fNoSameBunchPileup/fCentFT0MNoPileup columns; pileup cut is not applied." << std::endl;
    }
    std::cout << "Reading reco and truth data sample from MC table in " << dataAO2D << std::endl;
    if (dataAO2D.IsNull()) {
      std::cerr << "MC AO2D data requested, but dataAO2D is empty" << std::endl;
      return;
    }
    if (!readRecoAndTruthDataFromMCAO2D(dataAO2D, ptMin, ptMax, multMin, multMax,
                                        hTruthDataTmp, hRecoDataTmp,
                                        hMultTrueData, hMultRecoDataCand,
                                        hPtTrueData, hPtRecoData)) {
      return;
    }
  } else if (useAO2Data == 0) {
    std::cout << "Generating data sample with " << nDataEvents << " events" << std::endl;
    generatePseudoData(nDataEvents, fpt, nbd, genData, randData, eff, fakeMean,
                       ptMin, ptMax, multMin, multMax, hTruthDataTmp,
                       hRecoDataTmp, hMultTrueData, hMultRecoData, hPtTrueData,
                       hPtRecoData);
  } else if (useAO2Data == 2) {
    std::cout << "Using RecoTrain and TruthTrain as the data sample for closure/debug" << std::endl;
  } else {
    std::cerr << "Unknown useAO2Data option " << useAO2Data
              << ": use 1 for reco data, 0 for generated data, 2 for RecoTrain, 3 for MC2 reco+truth closure, 4 for MC-table reco+truth closure" << std::endl;
    return;
  }

  Matrix2D_t TruthData;
  Matrix2D_t RecoData;

  copyTH2ToMatrix2D(hTruthDataTmp, TruthData);
  copyTH2ToMatrix2D(hRecoDataTmp,  RecoData);

  if (useAO2Data == 2) {
    TruthData = TruthTrain;
    RecoData = RecoTrain;
  }

  // Reco input to Bayes can be corrected for non-signal reco rows measured in training.
  // FakeTrain and FeedInTrain are kept separate in output. When subtraction is disabled
  // for option 2, use the matched reco projection of Rcount for a strict response closure.
  Matrix2D_t RecoDataUnfold = RecoData;
  if (!subtractFakeFeed && useAO2Data == 2) {
    std::cout << "Bypassing fake/feed-in subtraction; using matched Rcount projection" << std::endl;
    for (int ir = 0; ir < NPT; ++ir) {
      for (int jr = 0; jr < NMULT; ++jr) {
        RecoDataUnfold[ir][jr] = 0.0;
        for (int it = 0; it < NPT; ++it) {
          for (int jt = 0; jt < NMULT; ++jt) {
            RecoDataUnfold[ir][jr] += Rcount[ir][jr][it][jt];
          }
        }
      }
    }
  } else if (subtractFakeFeed) {
    std::cout << "Subtracting fake/feed-in contribution from reco input" << std::endl;
    for (int ir = 0; ir < NPT; ++ir) {
      for (int jr = 0; jr < NMULT; ++jr) {
        double nonSignalFraction = 0.0;
        if (RecoAllTrain[ir][jr] > 0.0) {
          nonSignalFraction = (FakeTrain[ir][jr] + FeedInTrain[ir][jr]) / RecoAllTrain[ir][jr];
        }
        nonSignalFraction = std::clamp(nonSignalFraction, 0.0, 1.0);
        RecoDataUnfold[ir][jr] = RecoData[ir][jr] * (1.0 - nonSignalFraction);
      }
    }
  } else {
    std::cout << "Fake/feed-in subtraction disabled; using reco input as provided" << std::endl;
  }

  // ------------------------------------------------------------------
  // 5. MANUAL ITERATIVE BAYES
  // ------------------------------------------------------------------
  TRandom3 randPoisson(123);
  Matrix2D_t Tunfold;
  bool randomizeFlatPrior = false;
  if (priorMode == 0) {
    std::cout << "Using flat positive prior" << std::endl;
    Tunfold = makeFlatPrior();
  } else if (priorMode == 1) {
    std::cout << "Using TruthTrain fixed-point prior" << std::endl;
    Tunfold = TruthTrain;
  } else if (priorMode == 2) {
    std::cout << "Using RecoTrain-shaped prior" << std::endl;
    Tunfold = RecoTrain;
  } else if (priorMode == 3) {
    std::cout << "Using Poisson-randomized TruthTrain prior" << std::endl;
    Tunfold = poissonRandomizeMatrix(TruthTrain, randPoisson);
  } else if (priorMode == 4) {
    std::cout << "Using Poisson-randomized flat positive prior" << std::endl;
    Tunfold = makeFlatPrior();
    randomizeFlatPrior = true;
  } else {
    std::cerr << "Unknown priorMode " << priorMode
              << ": use 0 flat, 1 TruthTrain, 2 RecoTrain, 3 Poisson TruthTrain, 4 Poisson flat" << std::endl;
    return;
  }
  double sumReco = 0.0;
  double sumPrior = 0.0;

  for (int ir = 0; ir < NPT; ++ir) {
    for (int jr = 0; jr < NMULT; ++jr) {
      sumReco += RecoDataUnfold[ir][jr];
    }
  }

  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      sumPrior += Tunfold[it][jt];
    }
  }

  if (sumPrior > 0.0) {
    double scale = sumReco / sumPrior;
    for (int it = 0; it < NPT; ++it) {
      for (int jt = 0; jt < NMULT; ++jt) {
        Tunfold[it][jt] *= scale;
      }
    }
  }
  if (randomizeFlatPrior) {
    Tunfold = poissonRandomizeMatrix(Tunfold, randPoisson);
  }
  printMatrix2D(Tunfold, "Tunfold 0th");

  Matrix2D_t effs;
  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      double sum = 0.0;
      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          sum += Rprob[ir][jr][it][jt];
        }
      }
      effs[it][jt] = sum;
      std::cout << "eff[" << it << "][" << jt << "] = " << sum << std::endl;
    }
  }

  double sumTruthData = 0.0;
  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      sumTruthData += TruthData[it][jt];
    }
  }

  TH1D* hTruthL1RelIter = new TH1D("hTruthL1RelIter", "Truth L1 relative distance;iteration;#Sigma|unfold-truth| / #Sigmatruth", nIter, 0, nIter);
  TH1D* hRecoL1RelIter = new TH1D("hRecoL1RelIter", "Reco L1 relative distance;iteration;#Sigma|R#timesunfold-reco| / #Sigmareco", nIter, 0, nIter);

  for (int iter = 0; iter < nIter; ++iter) {
    Matrix2D_t Tnew;
    Matrix2D_t Denom;

    for (int ir = 0; ir < NPT; ++ir) {
      for (int jr = 0; jr < NMULT; ++jr) {
        double sum = 0.0;
        for (int it = 0; it < NPT; ++it) {
          for (int jt = 0; jt < NMULT; ++jt) {
            sum += Rprob[ir][jr][it][jt] * Tunfold[it][jt];
          }
        }
        Denom[ir][jr] = sum;
      }
    }

    for (int it = 0; it < NPT; ++it) {
      for (int jt = 0; jt < NMULT; ++jt) {
        double corr = 0.0;

        for (int ir = 0; ir < NPT; ++ir) {
          for (int jr = 0; jr < NMULT; ++jr) {
            if (Denom[ir][jr] > 0.0) {
              corr += Rprob[ir][jr][it][jt] * RecoDataUnfold[ir][jr] / Denom[ir][jr];
            }
          }
        }

        if (effs[it][jt] > 0.0) {
          Tnew[it][jt] = Tunfold[it][jt] * corr / effs[it][jt];
        } else {
          Tnew[it][jt] = 0.0;
        }
      }
    }

    double truthAbsDiff = 0.0;
    double recoAbsDiff = 0.0;
    double sumRecoFromTnew = 0.0;
    for (int it = 0; it < NPT; ++it) {
      for (int jt = 0; jt < NMULT; ++jt) {
        truthAbsDiff += std::abs(Tnew[it][jt] - TruthData[it][jt]);
      }
    }
    for (int ir = 0; ir < NPT; ++ir) {
      for (int jr = 0; jr < NMULT; ++jr) {
        double recoFromTnew = 0.0;
        for (int it = 0; it < NPT; ++it) {
          for (int jt = 0; jt < NMULT; ++jt) {
            recoFromTnew += Rprob[ir][jr][it][jt] * Tnew[it][jt];
          }
        }
        sumRecoFromTnew += recoFromTnew;
        recoAbsDiff += std::abs(recoFromTnew - RecoDataUnfold[ir][jr]);
      }
    }

    double truthRelDiff = (sumTruthData > 0.0) ? truthAbsDiff / sumTruthData : 0.0;
    double recoRelDiff = (sumReco > 0.0) ? recoAbsDiff / sumReco : 0.0;
    hTruthL1RelIter->SetBinContent(iter + 1, truthRelDiff);
    hRecoL1RelIter->SetBinContent(iter + 1, recoRelDiff);
    std::cout << "iter " << iter
              << " truthL1Rel=" << truthRelDiff
              << " recoL1Rel=" << recoRelDiff
              << " sumT=" << sumMatrix2D(Tnew)
              << " sumRecoFromT=" << sumRecoFromTnew
              << " sumRecoData=" << sumReco << std::endl;

    printMatrix2D(Tnew, "Tnew " + std::to_string(iter));
    Tunfold = Tnew;
  }

  // ------------------------------------------------------------------
  // 6. FINAL HISTOGRAMS FROM MATRICES
  // ------------------------------------------------------------------
  TH2D* hTruthTrain = new TH2D("hTruthTrain", "Training truth;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoTrain  = new TH2D("hRecoTrain",  "Training reco;p_{T};mult",  NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hTruthData  = new TH2D("hTruthData",  "Pseudo-data truth;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoData   = new TH2D("hRecoData",   "Pseudo-data reco;p_{T};mult",  NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hUnfold     = new TH2D("hUnfold",     "Manual Bayes unfolded;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hMiss       = new TH2D("hMiss",       "Missed;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hFakeTrain  = new TH2D("hFakeTrain",  "Training fakes;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hFeedInTrain = new TH2D("hFeedInTrain", "Training feed-in;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoAllTrain = new TH2D("hRecoAllTrain", "Training reco all;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);
  TH2D* hRecoDataUnfold = new TH2D("hRecoDataUnfold", "Non-signal-corrected reco data;p_{T};mult", NPT, ptMin, ptMax, NMULT, multMin, multMax);

  TH2D* hResponseFlat = new TH2D("hResponseFlat", "Flattened response;reco flat;truth flat",
                                 NPT * NMULT, 0, NPT * NMULT, NPT * NMULT, 0, NPT * NMULT);

  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      hTruthTrain->SetBinContent(it + 1, jt + 1, TruthTrain[it][jt]);
      hTruthData->SetBinContent(it + 1, jt + 1, TruthData[it][jt]);
      hUnfold->SetBinContent(it + 1, jt + 1, Tunfold[it][jt]);
      hMiss->SetBinContent(it + 1, jt + 1, Miss[it][jt]);
    }
  }

  for (int ir = 0; ir < NPT; ++ir) {
    for (int jr = 0; jr < NMULT; ++jr) {
      hRecoTrain->SetBinContent(ir + 1, jr + 1, RecoTrain[ir][jr]);
      hRecoData->SetBinContent(ir + 1, jr + 1, RecoData[ir][jr]);
      hFakeTrain->SetBinContent(ir + 1, jr + 1, FakeTrain[ir][jr]);
      hFeedInTrain->SetBinContent(ir + 1, jr + 1, FeedInTrain[ir][jr]);
      hRecoAllTrain->SetBinContent(ir + 1, jr + 1, RecoAllTrain[ir][jr]);
      hRecoDataUnfold->SetBinContent(ir + 1, jr + 1, RecoDataUnfold[ir][jr]);
    }
  }

  for (int ir = 0; ir < NPT; ++ir) {
    for (int jr = 0; jr < NMULT; ++jr) {
      int iRecoFlat = jr * NPT + ir + 1;
      for (int it = 0; it < NPT; ++it) {
        for (int jt = 0; jt < NMULT; ++jt) {
          int iTruthFlat = jt * NPT + it + 1;
          hResponseFlat->SetBinContent(iRecoFlat, iTruthFlat, Rprob[ir][jr][it][jt]);
        }
      }
    }
  }

  TH1D* hTruthTrainPtProj = hTruthTrain->ProjectionX("hTruthTrainPtProj");
  TH1D* hRecoTrainPtProj  = hRecoTrain->ProjectionX("hRecoTrainPtProj");
  TH1D* hTruthDataPtProj  = hTruthData->ProjectionX("hTruthDataPtProj");
  TH1D* hRecoDataPtProj   = hRecoData->ProjectionX("hRecoDataPtProj");
  TH1D* hUnfoldPtProj     = hUnfold->ProjectionX("hUnfoldPtProj");
  TH1D* hMissPtProj       = hMiss->ProjectionX("hMissPtProj");
  TH1D* hFakeTrainPtProj  = hFakeTrain->ProjectionX("hFakeTrainPtProj");
  TH1D* hFeedInTrainPtProj = hFeedInTrain->ProjectionX("hFeedInTrainPtProj");
  TH1D* hRecoAllTrainPtProj = hRecoAllTrain->ProjectionX("hRecoAllTrainPtProj");
  TH1D* hRecoDataUnfoldPtProj = hRecoDataUnfold->ProjectionX("hRecoDataUnfoldPtProj");
  TH1D* hTruthTrainMultProj = hTruthTrain->ProjectionY("hTruthTrainMultProj");
  TH1D* hRecoTrainMultProj  = hRecoTrain->ProjectionY("hRecoTrainMultProj");
  TH1D* hTruthDataMultProj  = hTruthData->ProjectionY("hTruthDataMultProj");
  TH1D* hRecoDataMultProj   = hRecoData->ProjectionY("hRecoDataMultProj");
  TH1D* hUnfoldMultProj     = hUnfold->ProjectionY("hUnfoldMultProj");

  auto setPoissonErrors = [](TH1D* h) {
    if (!h->GetSumw2N()) {
      h->Sumw2();
    }
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      const double content = h->GetBinContent(ib);
      h->SetBinError(ib, content > 0.0 ? std::sqrt(content) : 0.0);
    }
  };

  setPoissonErrors(hTruthTrainPtProj);
  setPoissonErrors(hTruthDataPtProj);
  setPoissonErrors(hUnfoldPtProj);
  setPoissonErrors(hTruthTrainMultProj);
  setPoissonErrors(hTruthDataMultProj);
  setPoissonErrors(hUnfoldMultProj);

  auto makeMultEventProj = [](TH1D* source, const char* name, const char* title) {
    TH1D* h = (TH1D*)source->Clone(name);
    h->SetTitle(title);
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      const double mult = h->GetBinCenter(ib);
      if (mult > 0.0) {
        h->SetBinContent(ib, h->GetBinContent(ib) / mult);
        h->SetBinError(ib, h->GetBinError(ib) / mult);
      } else {
        h->SetBinContent(ib, 0.0);
        h->SetBinError(ib, 0.0);
      }
    }
    return h;
  };

  TH1D* hTruthTrainMultEventProj = makeMultEventProj(hTruthTrainMultProj, "hTruthTrainMultEventProj", "TruthTrain mult projection divided by multiplicity");
  TH1D* hTruthDataMultEventProj  = makeMultEventProj(hTruthDataMultProj,  "hTruthDataMultEventProj",  "TruthData mult projection divided by multiplicity");
  TH1D* hUnfoldMultEventProj     = makeMultEventProj(hUnfoldMultProj,     "hUnfoldMultEventProj",     "Unfold mult projection divided by multiplicity");
  TH1D* hMissMultProj       = hMiss->ProjectionY("hMissMultProj");
  TH1D* hFakeTrainMultProj  = hFakeTrain->ProjectionY("hFakeTrainMultProj");
  TH1D* hFeedInTrainMultProj = hFeedInTrain->ProjectionY("hFeedInTrainMultProj");
  TH1D* hRecoAllTrainMultProj = hRecoAllTrain->ProjectionY("hRecoAllTrainMultProj");
  TH1D* hRecoDataUnfoldMultProj = hRecoDataUnfold->ProjectionY("hRecoDataUnfoldMultProj");

  auto makeSum = [](TH1D* first, TH1D* second, const char* name, const char* title) {
    TH1D* h = (TH1D*)first->Clone(name);
    h->SetTitle(title);
    if (!h->GetSumw2N()) {
      h->Sumw2();
    }
    h->Add(second);
    return h;
  };

  auto makeRatio = [](TH1D* numerator, TH1D* denominator, const char* name, const char* title) {
    TH1D* h = (TH1D*)numerator->Clone(name);
    h->SetTitle(title);
    if (!h->GetSumw2N()) {
      h->Sumw2();
    }
    h->Divide(denominator);
    h->SetOption("E1");
    return h;
  };

  auto makeSignalDenom = [](TH1D* all, TH1D* fake, TH1D* feed, const char* name, const char* title) {
    TH1D* h = (TH1D*)all->Clone(name);
    h->SetTitle(title);
    if (!h->GetSumw2N()) {
      h->Sumw2();
    }
    h->Add(fake, -1.0);
    h->Add(feed, -1.0);
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      if (h->GetBinContent(ib) < 0.0) {
        h->SetBinContent(ib, 0.0);
        h->SetBinError(ib, 0.0);
      }
    }
    return h;
  };

  TH1D* hFakeFeedTrainPtProj = makeSum(hFakeTrainPtProj, hFeedInTrainPtProj,
                                       "hFakeFeedTrainPtProj",
                                       "Training fake+feed-in;p_{T};entries");
  TH1D* hFakeFeedTrainMultProj = makeSum(hFakeTrainMultProj, hFeedInTrainMultProj,
                                         "hFakeFeedTrainMultProj",
                                         "Training fake+feed-in;mult;entries");
  TH1D* hSignalTrainPtProj = makeSignalDenom(hRecoAllTrainPtProj, hFakeTrainPtProj, hFeedInTrainPtProj,
                                             "hSignalTrainPtProj",
                                             "Training matched signal reco;p_{T};entries");
  TH1D* hSignalTrainMultProj = makeSignalDenom(hRecoAllTrainMultProj, hFakeTrainMultProj, hFeedInTrainMultProj,
                                               "hSignalTrainMultProj",
                                               "Training matched signal reco;mult;entries");

  TH1D* hFakeFractionPtProj = makeRatio(hFakeTrainPtProj, hRecoAllTrainPtProj,
                                        "hFakeFractionPtProj",
                                        "fake / all reco training;p_{T};fake / all");
  TH1D* hFeedInFractionPtProj = makeRatio(hFeedInTrainPtProj, hRecoAllTrainPtProj,
                                          "hFeedInFractionPtProj",
                                          "feed-in / all reco training;p_{T};feed-in / all");
  TH1D* hFakeFeedFractionPtProj = makeRatio(hFakeFeedTrainPtProj, hRecoAllTrainPtProj,
                                            "hFakeFeedFractionPtProj",
                                            "(fake+feed-in) / all reco training;p_{T};(fake+feed-in) / all");
  TH1D* hFakeFeedOverSignalPtProj = makeRatio(hFakeFeedTrainPtProj, hSignalTrainPtProj,
                                              "hFakeFeedOverSignalPtProj",
                                              "(fake+feed-in) / matched signal training;p_{T};(fake+feed-in) / signal");
  TH1D* hRecoDataPurityPtProj = makeRatio(hRecoDataUnfoldPtProj, hRecoDataPtProj,
                                          "hRecoDataPurityPtProj",
                                          "corrected reco data / raw reco data;p_{T};purity applied to data");

  TH1D* hFakeFractionMultProj = makeRatio(hFakeTrainMultProj, hRecoAllTrainMultProj,
                                          "hFakeFractionMultProj",
                                          "fake / all reco training;mult;fake / all");
  TH1D* hFeedInFractionMultProj = makeRatio(hFeedInTrainMultProj, hRecoAllTrainMultProj,
                                            "hFeedInFractionMultProj",
                                            "feed-in / all reco training;mult;feed-in / all");
  TH1D* hFakeFeedFractionMultProj = makeRatio(hFakeFeedTrainMultProj, hRecoAllTrainMultProj,
                                              "hFakeFeedFractionMultProj",
                                              "(fake+feed-in) / all reco training;mult;(fake+feed-in) / all");
  TH1D* hFakeFeedOverSignalMultProj = makeRatio(hFakeFeedTrainMultProj, hSignalTrainMultProj,
                                                "hFakeFeedOverSignalMultProj",
                                                "(fake+feed-in) / matched signal training;mult;(fake+feed-in) / signal");
  TH1D* hRecoDataPurityMultProj = makeRatio(hRecoDataUnfoldMultProj, hRecoDataMultProj,
                                            "hRecoDataPurityMultProj",
                                            "corrected reco data / raw reco data;mult;purity applied to data");

  TH1D* hRatioPt = (TH1D*)hUnfoldPtProj->Clone("hRatioPt");
  if (!hRatioPt->GetSumw2N()) {
    hRatioPt->Sumw2();
  }
  hRatioPt->SetTitle("Unfold / Truth p_{T}");
  hRatioPt->SetOption("E1");
  hRatioPt->Divide(hTruthDataPtProj);

  TH1D* hRatioMult = (TH1D*)hUnfoldMultProj->Clone("hRatioMult");
  if (!hRatioMult->GetSumw2N()) {
    hRatioMult->Sumw2();
  }
  hRatioMult->SetTitle("Unfold / Truth mult");
  hRatioMult->SetOption("E1");
  hRatioMult->Divide(hTruthDataMultProj);

  TH1D* hRatioMultEvent = (TH1D*)hUnfoldMultEventProj->Clone("hRatioMultEvent");
  if (!hRatioMultEvent->GetSumw2N()) {
    hRatioMultEvent->Sumw2();
  }
  hRatioMultEvent->SetTitle("Unfold / Truth event-style mult");
  hRatioMultEvent->SetOption("E1");
  hRatioMultEvent->Divide(hTruthDataMultEventProj);

  TString centMinLabel = TString::Format("%g", centMin);
  TString centMaxLabel = TString::Format("%g", centMax);
  centMinLabel.ReplaceAll(".", "p");
  centMaxLabel.ReplaceAll(".", "p");
  centMinLabel.ReplaceAll("-", "m");
  centMaxLabel.ReplaceAll("-", "m");
  TString outputName = TString::Format("manualBayes4D_cent%s_%s_pileup%d.root",
                                       centMinLabel.Data(), centMaxLabel.Data(),
                                       sameBCPolicy);

  TFile fout(outputName, "RECREATE");
  hTruthTrain->Write();
  hRecoTrain->Write();
  hTruthData->Write();
  hRecoData->Write();
  hUnfold->Write();
  hMiss->Write();
  hFakeTrain->Write();
  hFeedInTrain->Write();
  hRecoAllTrain->Write();
  hRecoDataUnfold->Write();
  hResponseFlat->Write();
  hTruthL1RelIter->Write();
  hRecoL1RelIter->Write();

  hPtTrueTrain->Write();
  hPtRecoTrain->Write();
  hPtTrueData->Write();
  hPtRecoData->Write();
  hMultTrueTrain->Write();
  hMultRecoTrain->Write();
  hMultTrueData->Write();
  hMultRecoData->Write();
  hMultRecoDataCand->Write();
  hMultRecoDataCollWithCand->Write();

  hTruthTrainPtProj->Write();
  hRecoTrainPtProj->Write();
  hTruthDataPtProj->Write();
  hRecoDataPtProj->Write();
  hUnfoldPtProj->Write();
  hMissPtProj->Write();
  hFakeTrainPtProj->Write();
  hFeedInTrainPtProj->Write();
  hRecoAllTrainPtProj->Write();
  hRecoDataUnfoldPtProj->Write();
  hTruthTrainMultProj->Write();
  hRecoTrainMultProj->Write();
  hTruthDataMultProj->Write();
  hRecoDataMultProj->Write();
  hUnfoldMultProj->Write();
  hTruthTrainMultEventProj->Write();
  hTruthDataMultEventProj->Write();
  hUnfoldMultEventProj->Write();
  hMissMultProj->Write();
  hFakeTrainMultProj->Write();
  hFeedInTrainMultProj->Write();
  hRecoAllTrainMultProj->Write();
  hRecoDataUnfoldMultProj->Write();
  hRatioPt->Write();
  hRatioMult->Write();
  hRatioMultEvent->Write();

  auto* qaDir = fout.mkdir("QA");
  if (qaDir) {
    qaDir->cd();
    hFakeTrainPtProj->Write();
    hFeedInTrainPtProj->Write();
    hFakeFeedTrainPtProj->Write();
    hSignalTrainPtProj->Write();
    hRecoAllTrainPtProj->Write();
    hRecoDataPtProj->Write();
    hRecoDataUnfoldPtProj->Write();
    hFakeFractionPtProj->Write();
    hFeedInFractionPtProj->Write();
    hFakeFeedFractionPtProj->Write();
    hFakeFeedOverSignalPtProj->Write();
    hRecoDataPurityPtProj->Write();

    hFakeTrainMultProj->Write();
    hFeedInTrainMultProj->Write();
    hFakeFeedTrainMultProj->Write();
    hSignalTrainMultProj->Write();
    hRecoAllTrainMultProj->Write();
    hRecoDataMultProj->Write();
    hRecoDataUnfoldMultProj->Write();
    hFakeFractionMultProj->Write();
    hFeedInFractionMultProj->Write();
    hFakeFeedFractionMultProj->Write();
    hFakeFeedOverSignalMultProj->Write();
    hRecoDataPurityMultProj->Write();

    hRatioPt->Write();
    hRatioMult->Write();
    hRatioMultEvent->Write();
    hTruthL1RelIter->Write();
    hRecoL1RelIter->Write();
    fout.cd();
  }

  fout.Close();


  std::cout << "Done. Output written to " << outputName << "\n";
}
