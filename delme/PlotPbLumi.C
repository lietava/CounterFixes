// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TestCTPScalers.C
/// \brief create CTP scalers, test it and add to database
/// \author Roman Lietava
// root -b -q "GetScalers.C(\"519499\", 1656286373953)"
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsCTP/Configuration.h"
#include <DataFormatsParameters/GRPLHCIFData.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"
#include <string>
#include <map>
#include <iostream>
#endif
using namespace o2::ctp;
using namespace std::chrono;
//
// if fillN = 0: pileup correction not done
// QCDB =1 : use for ongoing run
//
void PlotPbLumi(int runNumber, int fillN = 0, int QCDB = 0, std::string ccdbHost = "http://alice-ccdb.cern.ch")
{ // "http://ccdb-test.cern.ch:8080"
  std::string mCCDBPathCTPScalers = "CTP/Calib/Scalers";
  std::string mQCDBPathCTPScalers = "qc/CTP/Scalers";
  std::string mCCDBPathCTPConfig = "CTP/Config/Config";
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  ccdbMgr.setURL(ccdbHost);
  // Timestamp
  auto soreor = ccdbMgr.getRunDuration(runNumber);
  uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  std::cout << "Timestamp:" << timeStamp << std::endl;
  // Filling
  std::map<string, string> metadata;
  int nbc = 0;
  if (fillN) {
    std::string sfill = std::to_string(fillN);
    metadata["fillNumber"] = sfill;
    auto lhcifdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    auto bfilling = lhcifdata->getBunchFilling();
    std::vector<int> bcs = bfilling.getFilledBCs();
    nbc = bcs.size();
    std::cout << "Number of interacting bc:" << nbc << std::endl;
  }
  if (QCDB) { // use this option for ongoing run
    mCCDBPathCTPScalers = mQCDBPathCTPScalers;
    ccdbMgr.setURL("http://ccdb-test.cern.ch:8080");
    timeStamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    std::cout << "For scalers using Current time:" << timeStamp << std::endl;
  }
  // Scalers
  std::string srun = std::to_string(runNumber);
  metadata.clear(); // can be empty
  metadata["runNumber"] = srun;
  auto scl = ccdbMgr.getSpecific<CTPRunScalers>(mCCDBPathCTPScalers, timeStamp, metadata);
  if (scl == nullptr) {
    LOG(info) << "CTPRunScalers not in database, timestamp:" << timeStamp;
    return;
  }
  scl->convertRawToO2();
  std::vector<CTPScalerRecordO2> recs = scl->getScalerRecordO2();
  //
  // CTPConfiguration ctpcfg;
  auto ctpcfg = ccdbMgr.getSpecific<CTPConfiguration>(mCCDBPathCTPConfig, timeStamp, metadata);
  if (ctpcfg == nullptr) {
    LOG(info) << "CTPRunConfig not in database, timestamp:" << timeStamp;
    return;
  }
  std::vector<int> clslist = ctpcfg->getTriggerClassList();
  // std::vector<uint32_t> clslist = scl->getClassIndexes();
  std::map<int, int> clsIndexToScaler;
  std::cout << "Classes:";
  int i = 0;
  for (auto const& cls : clslist) {
    std::cout << cls << " ";
    clsIndexToScaler[cls] = i;
    i++;
  }
  std::cout << std::endl;
  std::vector<CTPClass> ctpcls = ctpcfg->getCTPClasses();
  int tsc = 255;
  int tce = 255;
  int vch = 255;
  int zncc = 255;
  for (auto const& cls : ctpcls) {
    if (cls.name.find("CMTVXTSC-B-NOPF") != std::string::npos && tsc == 255) {
      int itsc = cls.getIndex();
      tsc = clsIndexToScaler[itsc];
      // tsc = scl->getScalerIndexForClass(itsc);
      std::cout << cls.name << ":" << tsc << ":" << itsc << std::endl;
    }
    if (cls.name.find("CMTVXTCE-B-NOPF-CRU") != std::string::npos) {
      int itce = cls.getIndex();
      tce = clsIndexToScaler[itce];
      // tce = scl->getScalerIndexForClass(itce);
      std::cout << cls.name << ":" << tce << ":" << itce << std::endl;
    }
    if (cls.name.find("CMTVXVCH-B-NOPF-CRU") != std::string::npos) {
      int ivch = cls.getIndex();
      vch = clsIndexToScaler[ivch];
      // vch = scl->getScalerIndexForClass(ivch);
      std::cout << cls.name << ":" << vch << ":" << ivch << std::endl;
    }
    if (cls.name.find("C1ZNC-B") != std::string::npos) {
      int izncc = cls.getIndex();
      zncc = clsIndexToScaler[izncc];
    }
  }
  if (tsc == 255 || tce == 255 || vch == 255 || zncc == 255) {
    std::cout << " One of dcalers not available, check config to find alternative)" << std::endl;
    return;
  }
  //
  // Times
  double_t frev = 11245;
  double_t time0 = recs[0].epochTime;
  double_t timeL = recs[recs.size() - 1].epochTime;
  double_t Trun = timeL - time0;
  double_t orbit0 = recs[0].intRecord.orbit;
  double_t lumiZNCuncor = recs[recs.size() - 1].scalers[zncc].l1Before - recs[0].scalers[zncc].l1Before;
  double_t lumiTCEuncor = recs[recs.size() - 1].scalers[tce].lmBefore - recs[0].scalers[tce].lmBefore;

  int n = recs.size() - 1;

  std::cout << " Run duration:" << Trun << " Scalers size:" << n + 1 << " ZNC lumi (no pileupcorr):" << lumiZNCuncor/28 << " TCE lumi:" << lumiTCEuncor << std::endl;
  Double_t x[n], znc[n], zncpp[n], ex[n], eznc[n];
  Double_t tcetsctoznc[n], tcetoznc[n], vchtoznc[n];
  Double_t etcetsctoznc[n], etcetoznc[n], evchtoznc[n];

  double_t lumiZNC = 0;
  double_t lumiTCE = 0;
  for (int i = 0; i < n; i++) {
    ex[i] = 0;
    x[i] = (double_t)(recs[i + 1].intRecord.orbit + recs[i].intRecord.orbit) / 2. - orbit0;
    x[i] *= 88e-6;
    // x[i] = (double_t)(recs[i+1].epochTime + recs[i].epochTime)/2.;
    double_t tt = (double_t)(recs[i + 1].intRecord.orbit - recs[i].intRecord.orbit);
    tt = tt * 88e-6;
    //
    double_t znci = (double_t)(recs[i + 1].scalersInps[25] - recs[i].scalersInps[25]);
    double_t zncipp = znci;
    double_t mu = 0;
    if (fillN) {
      mu = -TMath::Log(1. - znci / tt / nbc / frev);
      zncipp = mu * nbc * frev * tt;
    }
    zncpp[i] = zncipp / 28.;
    lumiZNC += zncipp / 28.;
    znc[i] = zncipp / 28. / tt;
    eznc[i] = TMath::Sqrt(zncipp) / 28. / tt;
    //
    auto hadtsc = recs[i + 1].scalers[tsc].lmBefore - recs[i].scalers[tsc].lmBefore;
    auto hadtce = recs[i + 1].scalers[tce].lmBefore - recs[i].scalers[tce].lmBefore;
    auto had = hadtsc + hadtce;
    auto hadtcecorr = hadtce;
     if (fillN) {
      mu = -TMath::Log(1. - hadtce / tt / nbc / frev);
      hadtcecorr = mu * nbc * frev * tt;
    }
    lumiTCE += hadtcecorr;
    tcetsctoznc[i] = (double_t)(had) / zncpp[i];
    etcetsctoznc[i] = TMath::Sqrt(tcetsctoznc[i] * (1 - tcetsctoznc[i]) / zncpp[i]);
    had = recs[i + 1].scalers[tce].lmBefore - recs[i].scalers[tce].lmBefore;
    double_t tcec = had;
    tcetoznc[i] = (double_t)(had) / zncpp[i];
    had = recs[i + 1].scalers[vch].lmBefore - recs[i].scalers[vch].lmBefore;
    vchtoznc[i] = (double_t)(had) / zncpp[i];
    // std::cout << "mu:" << mu << " zncpp corr:" <<  zncipp << "  zncraw:" << znci << " tce:" << tcec << " tce/had" << tcec/zncpp[i] << std::endl;
  }
  std::cout << " lumi ZNCcorr:" << lumiZNC << " TCEcor:" << lumiTCE << std::endl;
  return;
  //
  TFile myfile("file.root", "RECREATE");
  gStyle->SetMarkerSize(0.5);
  // TGraph* gr1 = new TGraph(n, x, znc);
  TGraphErrors* gr1 = new TGraphErrors(n, x, znc, ex, eznc);
  TGraph* gr2 = new TGraph(n, x, tcetsctoznc);
  // TGraphErrors* gr2 = new TGraphErrors(n, x, tcetsctoznc, ex, etcetsctoznc); // nom and denom are strongly correlated
  TGraph* gr3 = new TGraph(n, x, tcetoznc);
  TGraph* gr4 = new TGraph(n, x, vchtoznc);
  gr1->SetMarkerStyle(20);
  gr2->SetMarkerStyle(21);
  gr3->SetMarkerStyle(23);
  gr4->SetMarkerStyle(23);
  gr1->SetTitle("R=ZNC/28 rate [Hz]; time[sec]; R");
  gr2->SetTitle("R=(TSC+TCE)*TVTX*B*28/ZNC; time[sec]; R");
  // gr2->SetTitle("R=TSC*TVTX*B*28/ZNC; time[sec]; R");
  //  gr2->GetHistogram()->SetMaximum(1.1);
  //  gr2->GetHistogram()->SetMinimum(0.9);
  gr3->SetTitle("R=(TCE)*TVTX*B*28/ZNC; time[sec]; R");
  // gr3->GetHistogram()->SetMaximum(0.6);
  // gr3->GetHistogram()->SetMinimum(0.4);
  gr4->SetTitle("R=(VCH)*TVTX*B*28/ZNC; time[sec]; R");
  // gr4->GetHistogram()->SetMaximum(0.6);
  // gr4->GetHistogram()->SetMinimum(0.4);
  TCanvas* c1 = new TCanvas("c1", srun.c_str(), 200, 10, 800, 500);
  TF1* fun = new TF1("poly0", "[0]+x*[1]");
  c1->Divide(2, 2);
  c1->cd(1);
  gr1->Draw("AP");
  c1->cd(2);
  fun->SetParameter(0, 1);
  gr2->Fit("poly0", "FM");
  gr2->Draw("AP");
  c1->cd(3);
  fun->SetParameter(0, 0.5);
  gr3->Fit("poly0");
  gr3->Draw("AP");
  c1->cd(4);
  gr4->Fit("poly0");
  gr4->Draw("AP");
  // getRate test:
  double tt = timeStamp / 1000.;
  std::pair<double, double> r1 = scl->getRateGivenT(tt, 25, 7);
  std::cout << "ZDC input getRateGivetT:" << r1.first / 28. << " " << r1.second / 28. << std::endl;
  std::pair<double, double> r2 = scl->getRateGivenT(tt, tce, 1);
  std::cout << "LM before TCE class getRateGivetT:" << r2.first << " " << r2.second << std::endl;
  gr1->Write();
  gr2->Write();
  gr3->Write();
  gr4->Write();
}
