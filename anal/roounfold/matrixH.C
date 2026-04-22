#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

constexpr int NPT   = 5;
constexpr int NMULT = 10;

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

void matrixH()
{
  const double ptMin   = 0.0;
  const double ptMax   = 10.0;
  const double multMin = 0.0;
  const double multMax = 100.0;

  const int nTrainEvents = 100000;
  const int nDataEvents  = 1000;
  const int nIter        = 10;

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

  // Response tensor counts
  Matrix4D_t Rcount = {};

  // ------------------------------------------------------------------
  // 1. TRAINING SAMPLE
  // ------------------------------------------------------------------
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
        } else {
          hMissTmp->Fill(trk.ptTrue, multTrue);
        }
      } else {
        hMissTmp->Fill(trk.ptTrue, multTrue);
      }
    }
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
  Matrix4D_t Rprob = {};
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

  printMatrix4D(Rprob, "Rprob");

  // ------------------------------------------------------------------
  // 4. PSEUDO-DATA SAMPLE
  // ------------------------------------------------------------------
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
    //multReco += nFake;

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

    // for (int i = 0; i < nFake; ++i) {
    //   double ptFake = fpt->GetRandom();
    //   if (ptFake >= ptMin && ptFake < ptMax) {
    //     hRecoDataTmp->Fill(ptFake, multReco);
    //     hPtRecoData->Fill(ptFake);
    //   }
    // }
  }

  Matrix2D_t TruthData;
  Matrix2D_t RecoData;

  copyTH2ToMatrix2D(hTruthDataTmp, TruthData);
  copyTH2ToMatrix2D(hRecoDataTmp,  RecoData);

  // optional debugging line from your original code:
  RecoData = RecoTrain;

  // ------------------------------------------------------------------
  // 5. MANUAL ITERATIVE BAYES
  // ------------------------------------------------------------------
  Matrix2D_t Tunfold = TruthTrain;
  printMatrix2D(Tunfold, "Tunfold 0th");

  double sumReco = 0.0;
  double sumPrior = 0.0;

  for (int ir = 0; ir < NPT; ++ir) {
    for (int jr = 0; jr < NMULT; ++jr) {
      sumReco += RecoData[ir][jr];
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
              corr += Rprob[ir][jr][it][jt] * RecoData[ir][jr] / Denom[ir][jr];
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

    printMatrix2D(Tnew, "Tnew");
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