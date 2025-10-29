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
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <string>
#include <map>
#include <iostream>
#endif
using namespace o2::ctp;
void PlotLumi(int runNumber = 564400, std::string inplumi = "MTVX")
{ //
  // what = 1: znc rate
  // what = 2: (TCE+TSC)/ZNC
  // what = 3: TCE/ZNC
  std::string mCCDBPathCTPScalers = "CTP/Calib/Scalers";
  std::string mCCDBPathCTPConfig = "CTP/Config/Config";
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  // Timestamp
  auto soreor = ccdbMgr.getRunDuration(runNumber);
  uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  std::cout << "Timestamp:" << timeStamp << std::endl;
  // Filling
  std::map<string, string> metadata;
  auto lhcifdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
  auto fillN = lhcifdata->getFillNumber();
  std::string sfill = std::to_string(fillN);
  metadata["fillNumber"] = sfill;
  auto bfilling = lhcifdata->getBunchFilling();
  std::vector<int> bcs = bfilling.getFilledBCs();
  int nbc = bcs.size();
  std::cout << "Number of interacting bc:" << nbc << std::endl;
  // Scalers
  std::string srun = std::to_string(runNumber);
  metadata.clear(); // can be empty
  metadata["runNumber"] = srun;
  //ccdbMgr.setURL("http://ccdb-test.cern.ch:8080");
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
    //std::cout << cls << " ";
    clsIndexToScaler[cls] = i;
    i++;
  }
  //std::cout << std::endl;
  std::vector<CTPClass> ctpcls = ctpcfg->getCTPClasses();
  int tvx = 255;
  for (auto const& cls : ctpcls) {
    if (cls.name.find(inplumi) != std::string::npos) {
      int itvx = cls.getIndex();
      tvx = clsIndexToScaler[itvx];
      // tsc = scl->getScalerIndexForClass(itsc);
      std::cout << cls.name << ":" << tvx << ":" << itvx << std::endl;
    }
  }
  if(tvx == 255) {
    std::cout << " One of dcalers not available, check config to find alternative)" << std::endl;
    return;
  }
  //
  // Anal
  //
  // Times
  double_t frev = 11245;
  double_t time0 = recs[0].epochTime;
  double_t timeL = recs[recs.size() - 1].epochTime;
  double_t Trun = timeL - time0;
  double_t orbit0 = recs[0].intRecord.orbit;
  int n = recs.size() - 1;
  std::cout << " Run duration:" << Trun << " Scalers size:" << n + 1 << std::endl;
  Double_t x[n];
  Double_t tvxa[n];
  for (int i = 0; i < n; i++) {
    x[i] = (double_t)(recs[i + 1].intRecord.orbit + recs[i].intRecord.orbit) / 2. - orbit0;
    x[i] *= 88e-6;
    // x[i] = (double_t)(recs[i+1].epochTime + recs[i].epochTime)/2.;
    double_t tt = (double_t)(recs[i + 1].intRecord.orbit - recs[i].intRecord.orbit);
    tt = tt * 88e-6;
    //
    //
    auto had = recs[i + 1].scalers[tvx].lmBefore - recs[i].scalers[tvx].lmBefore;
    // std::cout << recs[i+1].scalers[tce].lmBefore << std::endl;
    tvxa[i] = (double_t)(had) / tt;
  }
  //
  gStyle->SetMarkerSize(0.5);
  TGraph* gr1 = new TGraph(n, x, tvxa);
  gr1->SetMarkerStyle(20);
  gr1->SetTitle("Rate; time[sec]; R");
  TCanvas* c1 = new TCanvas("c1", srun.c_str(), 200, 10, 800, 500);
  gr1->Draw("AP");
  //std::pair<double, double> r2 = scl->getRateGivenT(tt, tce, 1);
  //std::cout << "LM before TCE class getRateGivetT:" << r2.first << " " << r2.second << std::endl;
}
