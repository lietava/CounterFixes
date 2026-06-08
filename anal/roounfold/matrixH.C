#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "TF1.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

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
                          TH2D* hRecoDataTmp,
                          TH1D* hMultRecoData,
                          TH1D* hPtRecoData)
{
  TChain collChain("O2npcollisiontabl");
  TChain candChain("O2nprecochargedca");
  fillChainFromAO2D(collChain, dataAO2D);
  fillChainFromAO2D(candChain, dataAO2D);

  if (collChain.GetEntries() <= 0) {
    std::cerr << "No nonPromptCascade NPCollisionTABLE entries found in "
              << dataAO2D << std::endl;
    return false;
  }
  if (candChain.GetEntries() <= 0) {
    std::cerr << "No nonPromptCascade NPRecoChargedCand entries found in "
              << dataAO2D << std::endl;
    return false;
  }
  if (!collChain.GetBranch("fMultNTracksNP")) {
    std::cerr << "NPCollisionTABLE is missing required branch fMultNTracksNP" << std::endl;
    return false;
  }
  if (!candChain.GetBranch("fPtRec") || !candChain.GetBranch("fIndexNPCollisionTable")) {
    std::cerr << "NPRecoChargedCand is missing required branches fPtRec and/or fIndexNPCollisionTable" << std::endl;
    return false;
  }

  int multReco = 0;
  std::vector<int> recoMultByCollision(collChain.GetEntries(), -1);
  collChain.SetBranchAddress("fMultNTracksNP", &multReco);
  for (Long64_t i = 0; i < collChain.GetEntries(); ++i) {
    collChain.GetEntry(i);
    recoMultByCollision[i] = multReco;
    if (multReco >= multMin && multReco < multMax) {
      hMultRecoData->Fill(multReco);
    }
  }

  float ptReco = 0.0f;
  int npCollisionId = -1;
  Long64_t nFilled = 0;
  Long64_t nIgnored = 0;
  candChain.SetBranchAddress("fPtRec", &ptReco);
  candChain.SetBranchAddress("fIndexNPCollisionTable", &npCollisionId);

  for (Long64_t i = 0; i < candChain.GetEntries(); ++i) {
    candChain.GetEntry(i);
    if (npCollisionId < 0 || npCollisionId >= static_cast<int>(recoMultByCollision.size())) {
      ++nIgnored;
      continue;
    }
    multReco = recoMultByCollision[npCollisionId];
    if (ptReco > ptMin && ptReco < ptMax &&
        multReco >= multMin && multReco < multMax) {
      hRecoDataTmp->Fill(ptReco, multReco);
      hPtRecoData->Fill(ptReco);
      ++nFilled;
    } else {
      ++nIgnored;
    }
  }

  std::cout << "Read nonPromptCascade real data from " << dataAO2D << std::endl;
  std::cout << "  reco candidates filled: " << nFilled
            << " ignored: " << nIgnored << std::endl;
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

// useAO2DTraining: true = train response from MC AO2D NPMCChargedTABLE; false = generate toy training.
// useAO2Data: 1 = read real reco data from dataAO2D; 0 = generate independent toy data;
//             2 = use RecoTrain closure input; 3 = read reco from dataAO2D and truth from truthAO2D.
void matrixH(TString trainAO2D = "AO2D.root", bool useAO2DTraining = true,
             TString dataAO2D = "AO2Ddata.root", TString truthAO2D = "AO2Dtruth.root",
             int useAO2Data = 1)
{
  const double ptMin   = 0.0;
  const double ptMax   = 10.0;
  const double multMin = 0.0;
  const double multMax = 80.0;

  const int nTrainEvents = 1000000;
  const int nDataEvents  = 1000000;
  const int nIter        = 100;

  const double eff      = 0.8;
  const double fakeMean = 1.5;

  TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", ptMin, ptMax);
  fpt->SetParameter(0, 0.5);

  double mu = 50.0;
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
  TH1D* hMultRecoData  = new TH1D("hMultRecoData",  "data mult reco",  NMULT, multMin, multMax);

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
  //             2: use RecoTrain for closure/debug
  //             3: read reco from dataAO2D and closure truth from truthAO2D
  // ------------------------------------------------------------------
  if (useAO2Data == 1 || useAO2Data == 3) {
    std::cout << "Reading reco data sample from " << dataAO2D << std::endl;
    if (dataAO2D.IsNull()) {
      std::cerr << "AO2D data requested, but dataAO2D is empty" << std::endl;
      return;
    }
    if (!readRecoDataFromAO2D(dataAO2D, ptMin, ptMax, multMin, multMax,
                              hRecoDataTmp, hMultRecoData, hPtRecoData)) {
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
  } else if (useAO2Data == 0 || useAO2Data == 2) {
    std::cout << "Generating data sample with " << nDataEvents << " events" << std::endl;
    generatePseudoData(nDataEvents, fpt, nbd, genData, randData, eff, fakeMean,
                       ptMin, ptMax, multMin, multMax, hTruthDataTmp,
                       hRecoDataTmp, hMultTrueData, hMultRecoData, hPtTrueData,
                       hPtRecoData);
  } else {
    std::cerr << "Unknown useAO2Data option " << useAO2Data
              << ": use 1 for reco data, 0 for generated data, 2 for RecoTrain, 3 for MC2 reco+truth closure" << std::endl;
    return;
  }

  Matrix2D_t TruthData;
  Matrix2D_t RecoData;

  copyTH2ToMatrix2D(hTruthDataTmp, TruthData);
  copyTH2ToMatrix2D(hRecoDataTmp,  RecoData);

  if (useAO2Data == 2) {
    RecoData = RecoTrain;
  }

  // Reco input to Bayes is corrected for non-signal reco rows measured in training.
  // FakeTrain and FeedInTrain are kept separate in output, but both are removed here.
  Matrix2D_t RecoDataUnfold = RecoData;
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

  // ------------------------------------------------------------------
  // 5. MANUAL ITERATIVE BAYES
  // ------------------------------------------------------------------
  TRandom3 randPoisson(123);
  //Matrix2D_t Tunfold = TruthTrain;
  //Matrix2D_t Tunfold =  poissonRandomizeMatrix(TruthTrain,randPoisson);
  //Matrix2D_t Tunfold = RecoTrain;
  Matrix2D_t Tunfold = makeFlatPrior();
  printMatrix2D(Tunfold, "Tunfold 0th");

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

  TH1D* hTruthPt     = hTruthData->ProjectionX("hTruthPt");
  TH1D* hRecoPt      = hRecoData->ProjectionX("hRecoPt");
  TH1D* hUnfoldPt    = hUnfold->ProjectionX("hUnfoldPt");
  TH1D* hTruthMult   = hTruthData->ProjectionY("hTruthMult");
  TH1D* hRecoMult    = hRecoData->ProjectionY("hRecoMult");
  TH1D* hUnfoldMult  = hUnfold->ProjectionY("hUnfoldMult");

  TH1D* hRatioPt = (TH1D*)hUnfoldPt->Clone("hRatioPt");
  hRatioPt->SetTitle("Unfold / Truth p_{T}");
  hRatioPt->Divide(hTruthPt);

  TH1D* hRatioMult = (TH1D*)hUnfoldMult->Clone("hRatioMult");
  hRatioMult->SetTitle("Unfold / Truth mult");
  hRatioMult->Divide(hTruthMult);

  TFile fout("manualBayes4D.root", "RECREATE");
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

  hPtTrueTrain->Write();
  hPtRecoTrain->Write();
  hPtTrueData->Write();
  hPtRecoData->Write();
  hMultTrueTrain->Write();
  hMultRecoTrain->Write();
  hMultTrueData->Write();
  hMultRecoData->Write();

  hTruthPt->Write();
  hRecoPt->Write();
  hUnfoldPt->Write();
  hTruthMult->Write();
  hRecoMult->Write();
  hUnfoldMult->Write();
  hRatioPt->Write();
  hRatioMult->Write();

  fout.Close();


  std::cout << "Done. Output written to manualBayes4D.root\n";
}
