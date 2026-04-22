
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"

constexpr int NPT = 2;
constexpr int NMULT = 2;
// using Matrix2D = std::vector<std::vector<double>>;
// typedef double Matrix2D_t[NPT][NMULT];
struct Matrix2D_t {
  double data[NPT][NMULT];
  Matrix2D_t() { std::fill(&data[0][0], &data[0][0] + NPT * NMULT, 0.0); }
  // assignment operator
  Matrix2D_t &operator=(const Matrix2D_t &other) {
    if (this != &other) {
      std::copy(&other.data[0][0], &other.data[0][0] + NPT * NMULT,
                &data[0][0]);
    }
    return *this;
  }

  // round bracket access
  double &operator()(int i, int j) { return data[i][j]; }
  const double &operator()(int i, int j) const { return data[i][j]; }

  // square bracket access
  double *operator[](int i) { return data[i]; }
  const double *operator[](int i) const { return data[i]; }
};
struct TrackInfo {
  double ptTrue;
  double ptReco;
  bool isReco;
};
// using Matrix2D = std::vector<std::vector<double>>;
// using Matrix2D_t = double[NPT][NMULT];
//                          reco         true (gen)
using Matrix4D_t = double[NPT][NMULT][NPT][NMULT];
void printMatrix2D(Matrix2D_t &m, std::string name = "") {
  std::cout << "===> " << name << std::endl;
  for (int i = 0; i < NPT; i++) {
    for (int j = 0; j < NMULT; j++) {
      std::cout << m[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
void printMatrix4D(Matrix4D_t &m, std::string name = "") {
  std::cout << "===> " << name << std::endl;
  for (int jt = 0; jt < NMULT; jt++) {
    for (int it = 0; it < NPT; it++) {
      std::cout << "ptTrue:" << it << " jMULTTrue:" << jt << std::endl;
      for (int jr = 0; jr < NMULT; jr++) {
        for (int ir = 0; ir < NPT; ir++) {
          // std::cout << ir << " " << jr << " " << it << " " << jt <<
          // std::endl;
          std::cout << m[ir][jr][it][jt] << " "; // << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }
  std::cout << std::endl;
}
//
double sigma_pt(double pt) {
  return pt * std::sqrt((0.02 / pt) * (0.02 / pt) + 0.01 * 0.01);
}
//
void matrix() {
  const double ptMin = 0.0;
  const double ptMax = 10.0;
  const double multMin = 0.0;
  const double multMax = 100.0;

  const int nTrainEvents = 1000;
  const int nDataEvents = 1000;
  const int nIter = 10;

  const double eff = 0.5;
  const double fakeMean = 1.5;

  TF1 *fpt = new TF1("fpt", "x*exp(-x/[0])", ptMin, ptMax);
  fpt->SetParameter(0, 0.5);

  double mu = 50.0;
  int r = 20;
  double p = r / (r + mu);

  std::mt19937 genTrain(12345);
  std::mt19937 genData(67890);
  std::negative_binomial_distribution<int> nbd(r, p);

  TRandom3 randTrain(1);
  TRandom3 randData(2);

  // Response tensor counts:
  // Rcount[ir][jr][it][jt]
  Matrix4D_t Rcount = {};
  // Miss counts per truth bin
  Matrix2D_t Miss;
  // Priors / truth / reco / unfolded
  Matrix2D_t TruthTrain;
  Matrix2D_t RecoTrain;
  Matrix2D_t TruthData;
  Matrix2D_t RecoData;
  // Helper bin finders
  auto findPtBin = [&](double pt) -> int {
    if (pt < ptMin || pt >= ptMax)
      return -1;
    int b = int((pt - ptMin) / (ptMax - ptMin) * NPT);
    if (b < 0 || b >= NPT)
      return -1;
    return b;
  };

  auto findMultBin = [&](double mult) -> int {
    if (mult < multMin || mult >= multMax)
      return -1;
    int b = int((mult - multMin) / (multMax - multMin) * NMULT);
    if (b < 0 || b >= NMULT)
      return -1;
    return b;
  };

  // ------------------------------------------------------------------
  // 1. TRAINING SAMPLE: build response tensor Rcount and training truth
  // ------------------------------------------------------------------
  for (int ev = 0; ev < nTrainEvents; ++ev) {
    int multTrue = nbd(genTrain);
    int jt = findMultBin(multTrue);
    if (jt < 0)
      continue;

    std::vector<TrackInfo> tracks;
    tracks.reserve(multTrue);

    int multReco = 0;

    for (int i = 0; i < multTrue; ++i) {
      double ptTrue = fpt->GetRandom();
      ptTrue = randTrain.Rndm() * 10.;
      bool isReco = (randTrain.Rndm() < eff);
      double ptReco = -1.0;

      if (isReco) {
        ptReco = randTrain.Gaus(ptTrue, sigma_pt(ptTrue));
        ptReco = ptTrue;  // debug
        if (ptReco > ptMin && ptReco < ptMax) {
          ++multReco;
        } else {
          isReco = false;
          ptReco = -1.0;
        }
      }

      tracks.push_back({ptTrue, ptReco, isReco});
    }

    // int nFake = randTrain.Poisson(fakeMean);
    // multReco += nFake;

    int jr = findMultBin(multReco);
    if (jr < 0) {
      std::cout << " jr not found for multReco:" << multReco << " Exiting."
                << std::endl;
      exit(1);
    }

    for (const auto &trk : tracks) {
      int it = findPtBin(trk.ptTrue);
      if (it < 0) {
        std::cout << " it not found. Exiting." << std::endl;
        exit(1);
        // continue;
      }

      TruthTrain[it][jt] += 1.0;

      if (trk.isReco) {
        int ir = findPtBin(trk.ptReco);
        if (ir >= 0) {
          Rcount[ir][jr][it][jt] += 1.0;
          RecoTrain[ir][jr] += 1.0;
        } else {
          Miss[it][jt] += 1.0;
          std::cout << "ir not found for:" << trk.ptReco << std::endl;
          exit(1);
        }
      } else {
        Miss[it][jt] += 1.0;
      }
    }

    // fake tracks are not associated with truth bins in the tensor
    // if you want explicit fake handling in manual Bayes, that is an extra term
    // for this first implementation we keep reco pseudo-data containing fakes,
    // but the response tensor only models matched + misses.
  }
  printMatrix2D(TruthTrain, "TruthTrain");
  printMatrix2D(Miss, "Miss");
  printMatrix2D(RecoTrain,"RecoTrain");
  // ------------------------------------------------------------------
  // 2. NORMALIZE RESPONSE:
  //    Rprob[ir][jr][it][jt] = P(reco bin ir,jr | true bin it,jt)
  // ------------------------------------------------------------------
  Matrix4D_t Rprob = {};
  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      double norm = Miss[it][jt];
      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          norm += Rcount[ir][jr][it][jt];
        }
      }
      // std::cout << it << " " << jt << " norm:" << norm << std::endl;
      if (norm <= 0.0) {
        std::cout << "norm < 0. Exiting." << std::endl;
        exit(1);
        // continue;
      }

      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          Rprob[ir][jr][it][jt] = Rcount[ir][jr][it][jt] / norm;
        }
      }
    }
  }
  printMatrix4D(Rprob, "Rprob");
  // ------------------------------------------------------------------
  // 3. PSEUDO-DATA SAMPLE: build truth + reco to unfold
  // ------------------------------------------------------------------
  for (int ev = 0; ev < nDataEvents; ++ev) {
    int multTrue = nbd(genData);
    int jt = findMultBin(multTrue);
    if (jt < 0)
      continue;

    std::vector<TrackInfo> tracks;
    tracks.reserve(multTrue);

    int multReco = 0;

    for (int i = 0; i < multTrue; ++i) {
      double ptTrue = fpt->GetRandom();
      bool isReco = (randData.Rndm() < eff);
      double ptReco = -1.0;

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

    int jr = findMultBin(multReco);
    if (jr < 0)
      continue;

    for (const auto &trk : tracks) {
      int it = findPtBin(trk.ptTrue);
      if (it >= 0) {
        TruthData[it][jt] += 1.0;
      }

      if (trk.isReco) {
        int ir = findPtBin(trk.ptReco);
        if (ir >= 0) {
          RecoData[ir][jr] += 1.0;
        }
      }
    }

    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();
      int ir = findPtBin(ptFake);
      if (ir >= 0) {
        RecoData[ir][jr] += 1.0;
      }
    }
  }
  RecoData = RecoTrain;
  //RecoData = TruthTrain;
  // ------------------------------------------------------------------
  // 4. MANUAL ITERATIVE BAYES
  //
  // t_new[it][jt] = t[it][jt] * sum_{ir,jr} Rprob[ir][jr][it][jt]
  //                                * d[ir][jr] / sum_{kt,lt}
  //                                Rprob[ir][jr][kt][lt] * t[kt][lt]
  // ------------------------------------------------------------------

  // Start prior from training truth
  //Matrix2D_t Tunfold = RecoTrain;
  Matrix2D_t Tunfold = TruthTrain;
  printMatrix2D(Tunfold, "Tunfold 0th");
  // Normalize prior to pseudo-data total truth-like scale if wanted
  double sumReco = 0.0;
  double sumPrior = 0.0;
  for (int ir = 0; ir < NPT; ++ir)
    for (int jr = 0; jr < NMULT; ++jr)
      sumReco += RecoData[ir][jr];

  for (int it = 0; it < NPT; ++it)
    for (int jt = 0; jt < NMULT; ++jt)
      sumPrior += Tunfold[it][jt];

  if (sumPrior > 0.0) {
    double scale = sumReco / sumPrior;
    for (int it = 0; it < NPT; ++it)
      for (int jt = 0; jt < NMULT; ++jt)
        Tunfold[it][jt] *= scale;
  }
  // effs
  Matrix2D_t effs;
  for (int it = 0; it < NPT; ++it) {
    for (int jt = 0; jt < NMULT; ++jt) {
      double sum = 0;
      for (int ir = 0; ir < NPT; ++ir) {
        for (int jr = 0; jr < NMULT; ++jr) {
          sum += Rprob[ir][jr][it][jt];
        }
      }
      std::cout << "eff:" << sum << std::endl;
      effs[it][jt] = sum;
    }
  }
  //
  for (int iter = 0; iter < nIter; ++iter) {
    Matrix2D_t Tnew;
    // Predicted reco: d_pred[ir][jr] = sum_{it,jt} Rprob * Tunfold
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
              corr += Rprob[ir][jr][it][jt] * RecoData[ir][jr] / Denom[ir][jr];
            } else {
              std::cout << "Warning: denom <= 0:" << Denom[ir][jr]
                        << " pTr:" << ir << " mult rec:" << jr << std::endl;
              std::cout << Rprob[ir][jr][it][jt] << " " << RecoData[ir][jr]
                        << std::endl;
            }
          }
        }
        if (effs[it][jt] > 0.) {
          Tnew[it][jt] = Tunfold[it][jt] * corr / effs[it][jt];
        } else {
          Tnew[it][jt] = 0.0;
        }
      }
    }
    printMatrix2D(Tnew, "Tnew");
    Tunfold = Tnew;
  }

  // ------------------------------------------------------------------
  // 5. Convert arrays to ROOT histograms
  // ------------------------------------------------------------------
  TH2D *hTruthTrain = new TH2D("hTruthTrain", "Training truth;p_{T};mult", NPT,
                               ptMin, ptMax, NMULT, multMin, multMax);

  TH2D *hRecoTrain = new TH2D("hRecoTrain", "Training truth;p_{T};mult", NPT,
                              ptMin, ptMax, NMULT, multMin, multMax);

  TH2D *hTruthData = new TH2D("hTruthData", "Pseudo-data truth;p_{T};mult", NPT,
                              ptMin, ptMax, NMULT, multMin, multMax);

  TH2D *hRecoData = new TH2D("hRecoData", "Pseudo-data reco;p_{T};mult", NPT,
                             ptMin, ptMax, NMULT, multMin, multMax);

  TH2D *hUnfold = new TH2D("hUnfold", "Manual Bayes unfolded;p_{T};mult", NPT,
                           ptMin, ptMax, NMULT, multMin, multMax);
  TH2D *hMiss = new TH2D("hMiss", "Missed ;p_{T};mult", NPT, ptMin, ptMax,
                         NMULT, multMin, multMax);
  // flattened response matrix for debugging
  TH2D *hResponseFlat =
      new TH2D("hResponseFlat", "Flattened response;reco flat;truth flat",
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
      hRecoData->SetBinContent(ir + 1, jr + 1, RecoData[ir][jr]);
      hRecoTrain->SetBinContent(ir + 1, jr + 1, RecoTrain[ir][jr]);
    }
  }

  for (int ir = 0; ir < NPT; ++ir) {
    for (int jr = 0; jr < NMULT; ++jr) {
      int iRecoFlat = jr * NPT + ir + 1;
      for (int it = 0; it < NPT; ++it) {
        for (int jt = 0; jt < NMULT; ++jt) {
          int iTruthFlat = jt * NPT + it + 1;
          hResponseFlat->SetBinContent(iRecoFlat, iTruthFlat,
                                       Rprob[ir][jr][it][jt]);
        }
      }
    }
  }

  // ------------------------------------------------------------------
  // 6. Some simple diagnostics
  // ------------------------------------------------------------------
  TH1D *hTruthPt = hTruthData->ProjectionX("hTruthPt");
  TH1D *hRecoPt = hRecoData->ProjectionX("hRecoPt");
  TH1D *hUnfoldPt = hUnfold->ProjectionX("hUnfoldPt");

  TH1D *hTruthMult = hTruthData->ProjectionY("hTruthMult");
  TH1D *hRecoMult = hRecoData->ProjectionY("hRecoMult");
  TH1D *hUnfoldMult = hUnfold->ProjectionY("hUnfoldMult");

  // ratios
  TH1D *hRatioPt = (TH1D *)hUnfoldPt->Clone("hRatioPt");
  hRatioPt->SetTitle("Unfold / Truth p_{T}");
  hRatioPt->Divide(hTruthPt);

  TH1D *hRatioMult = (TH1D *)hUnfoldMult->Clone("hRatioMult");
  hRatioMult->SetTitle("Unfold / Truth mult");
  hRatioMult->Divide(hTruthMult);

  // ------------------------------------------------------------------
  // 7. Write output
  // ------------------------------------------------------------------
  TFile fout("manualBayes4D.root", "RECREATE");
  hTruthTrain->Write();
  hRecoTrain->Write();
  hTruthData->Write();
  hRecoData->Write();
  hUnfold->Write();
  hResponseFlat->Write();
  hTruthPt->Write();
  hRecoPt->Write();
  hUnfoldPt->Write();
  hTruthMult->Write();
  hRecoMult->Write();
  hUnfoldMult->Write();
  hRatioPt->Write();
  hRatioMult->Write();
  hMiss->Write();
  fout.Close();
  // hTruthTrain->Print("all");
  // hUnfold->Print("all");

  std::cout << "Done. Output written to manualBayes4D.root\n";
}