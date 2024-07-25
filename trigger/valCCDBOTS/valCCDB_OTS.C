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
// O2 includes


#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <memory>
#include <bitset>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Configuration.h"
#include "DataFormatsCTP/Scalers.h"

//
//
const int Ndim = 128;
ofstream* mylog;
//
//================================
//
//
// main structure bcInfo alias ZorroHelper
//
struct ZorroHelper {
  uint64_t bcAOD, bcEvSel, trigMask[2], selMask[2];
  void print() const;
  //ClassDefNV(ZorroHelper, 1);
};
void ZorroHelper::print() const
{
  *mylog << bcEvSel << " " << bcAOD << " TM:" << std::hex << trigMask[1] << " " << trigMask[0]  << " SM:" << selMask[1] << " " << selMask[0] << std::dec;
  *mylog << std::endl;
};
#pragma link C++ class ZorroHelper + ;
#pragma link C++ class std::vector < ZorroHelper> + ;

//
// bcInfos
//
typedef ZorroHelper bcInfo;
struct bcInfos
{
  int Ndimused = 64;
  std::vector<std::string> labels;
  const std::string ccdbTest = "http://ccdb-test.cern.ch:8080";
  const std::string ccdbProd = "http://alice-ccdb.cern.ch";
  const std::string mCCDBPathOTS= "Users/m/mpuccio/EventFiltering/OTS";
  const std::string mCCDBPathOTSrl= "Users/r/rlietava/EventFiltering/OTS";
  const std::string mCCDBPathCTPConfig = "CTP/Config/Config";
  const std::string mCCDBPathCTPScalers = "CTP/Calib/Scalers";
  int runNumber;
  //const std::string mCCDBPathTrigScalersSkimmed = "Users/m/mpuccio/EventFiltering/OTS/SelectionCounters";
  bcInfos() = default;
  std::vector<bcInfo> bcs;
  std::array<int, Ndim> selectionCounters{0};  // counting bits from disk or CCDB file
  std::array<int, Ndim> triggerCounters{0};    // counting bits from disk or CCDB file
  std::array<double_t,Ndim> downscaleFactors{0};
  int totalSelected = 0;
  int totalTriggered = 0;
  TH1* mFilterCounters = nullptr;       // from CCDB
  TH1* mSelectionCounters = nullptr;    // from CCDB
  TH1* mInspectedTVX = nullptr;         // from CCDB
  std::array<double_t,Ndim> eff{0};
  std::array<double_t,Ndim> err{0};
  uint64_t TVXCTP{0},TVXCCDB{0};
  //
  void getDataCCDB(int run, int Nmax = 0, bool print = 0);
  double_t getDuration();
  void selectionEfficiency();
  void selectionEfficiencyBiased() const;
  void extractLabels();
  void checkConsistency();
  //
  void printbcInfoVect();
  void printbcInfoVect(std::vector<bcInfo>& bcs);
};
//
void bcInfos::getDataCCDB(int runNumber, int Nmax, bool print)
{
  this->runNumber = runNumber;
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  //ccdbMgr.setURL(ccdbProd);
  //
  auto soreor = ccdbMgr.getRunDuration(runNumber);
  ULong64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  //timeStamp = 1688349854935;
  //
  std::map<std::string,std::string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  *mylog << "run:" << metadata["runNumber"] << " Timestamp:" << timeStamp << std::endl;
  //
  auto scl = ccdbMgr.getSpecific<o2::ctp::CTPRunScalers>(mCCDBPathCTPScalers, timeStamp, metadata);
  scl->convertRawToO2();
  std::vector<o2::ctp::CTPScalerRecordO2> recs = scl->getScalerRecordO2();
  auto ctpcfg = ccdbMgr.getSpecific<o2::ctp::CTPConfiguration>(mCCDBPathCTPConfig, timeStamp, metadata);
  std::vector<o2::ctp::CTPClass> ctpcls = ctpcfg->getCTPClasses();
  std::string className = "MTVX";
  if(runNumber < 534202) {
    className = "minbias_TVX_L0";
  }
  int tvx = 255;
  for (auto const& cls : ctpcls) {
    if (cls.name.find(className) != std::string::npos && tvx == 255) {
      int itvx = cls.getIndex();
      tvx = scl->getScalerIndexForClass(itvx);
      std::cout << cls.name << ":" << tvx << ":" << itvx << std::endl;
    }
  }
  if( tvx != 255) {
    uint64_t scal0;
    uint64_t scalL;
    if(className.find("0") == std::string::npos) {
      scal0 = recs[0].scalers[tvx].lmBefore ;
      scalL = recs[recs.size() - 1].scalers[tvx].lmBefore;
    } else {
      scal0 = recs[0].scalers[tvx].l0Before ;
      scalL = recs[recs.size() - 1].scalers[tvx].l0Before;
    }
    std::cout << "CTP TVX:" << scalL-scal0 << std::endl;
    TVXCTP = scalL-scal0;
  } else {
    std::cout << "Class not found:" << className << std::endl;
  }
  //
  std::string path = mCCDBPathOTS+"/InspectedTVX";
  mInspectedTVX = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if(mInspectedTVX == nullptr) {
    *mylog << "InspectedTVX not found" << std::endl;
    return;
  }
  std::cout << "TVX CCDB:" << mInspectedTVX->GetEntries() << std::endl;
  TVXCCDB = mInspectedTVX->GetEntries();
  //
  path = mCCDBPathOTS+"/FilterCounters";
  mFilterCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mFilterCounters == nullptr ) {
    *mylog << "FilterScalers not found" << std::endl;
    return;
  }
  int Nfiltered = mFilterCounters->GetNbinsX();
  *mylog << "FilterCounters N bins:" << mFilterCounters->GetNbinsX() << std::endl;
  path = mCCDBPathOTS+"/SelectionCounters";
  mSelectionCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mSelectionCounters == nullptr ) {
    *mylog << "SelectionCounters not found" << std::endl;
    return;
  }
  int Nselected = mSelectionCounters->GetNbinsX();
  *mylog << "SelectionCounters N bins:" << mSelectionCounters->GetNbinsX() << std::endl;
  if(Nfiltered != Nselected) {
    *mylog << "ERROR Nfiltered != Nselected" << std::endl;
    exit(1);
  }
  Ndimused = Nselected - 2;
  // bcinfo
  path = mCCDBPathOTS+"/ZorroHelpers";
  std::vector<bcInfo>* bcinfos = ccdbMgr.getSpecific<std::vector<bcInfo>>(path, timeStamp, metadata);
  *mylog << " bcInfos size:" << (*bcinfos).size() << std::endl;
  int nread = Nmax;
  if (Nmax == 0) {
    nread = (*bcinfos).size();
  }
  for(int i = 0; i < nread; i++) {
    bcInfo bci = (*bcinfos)[i];
    bcs.push_back(bci);
    for(int m = 0; m < Ndimused; m++) {
      int index = m/64;
      int mask = m % 64;
      if(bci.trigMask[index] & (1ull << mask)) {
        triggerCounters[m]++;
      }
      if(bci.selMask[index] & (1ull << mask)) {
        selectionCounters[m]++;
      }
    }
  }
  std::sort(bcs.begin(), bcs.end(), [](const bcInfo& a, const bcInfo& b) { return a.bcEvSel < b.bcEvSel; });
  if(print) {
    extractLabels();
    for(int i = 0; i < Ndimused; i++) {
      std::string label = "-";
      if(i < (int)labels.size()) {
        label = labels[i];
      }
      *mylog << i << " " << label << " " << triggerCounters[i] << " " << mFilterCounters->GetBinContent(i+2) << "||" << selectionCounters[i] << " " << mSelectionCounters->GetBinContent(i+2) << std::endl;
    }
  }
}
double_t bcInfos::getDuration()
{
  ULong64_t tt =   bcs.back().bcEvSel - bcs.front().bcEvSel;
  double_t t = tt*25.e-9;
  return t;
}
void bcInfos::printbcInfoVect()
{
  printbcInfoVect(bcs);
};
//
// utils
//
void bcInfos::extractLabels()
{
  labels.clear();
  TH1* hist1 = mFilterCounters;
  *mylog << "Labels Nbins:" <<  hist1->GetNbinsX() << std::endl;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      //std::string labelc = label.substr(1,6);
      labels.push_back(label);
    }
  }
  //*mylog << "Label size CCDB:" << labels.size() << std::endl;
  //for(size_t i = 0; i < labels.size(); i++) {
  //  *mylog << i << " "  << labels[i] << std::endl;
  //}
}
void bcInfos::selectionEfficiency()
{
  *mylog << "===> Selection efficiency - downscaling" << std::endl;
  //std::array<double_t,Ndim> eff{0};
  for(int i =0; i < Ndimused; i++) {
    double_t before = mFilterCounters->GetBinContent(i+2);
    double_t after = selectionCounters[i];
    eff[i] = 0;
    err[i] = 0;
    if( before != 0.) {
      eff[i] = after/before;
      err[i] = sqrt(eff[i]*(1-eff[i])/before);
    } else {
      //*mylog << "after:" << after << " before:" << before << std::endl;
    }
    std::string label = "?";
    if(i < (int) labels.size()) {
      label = labels[i];
    }
    downscaleFactors[i] = eff[i];
    *mylog << "BIT:" << std::setw(2) << i << "  " << std::setw(25) << label << "   b:" << std::setw(10) << before << " a:" << std::setw(10) << after << " eff:" << eff[i] << std::endl;
  }
  //printSelectionEfficiency(eff);
}
void bcInfos::selectionEfficiencyBiased() const
{
  *mylog << "===> Selection efficiency - downscaling 2" << std::endl;
  double_t before = 0;
  double_t after = 0;
  for(int i =0; i < Ndimused; i++) {
    before = triggerCounters[i];
    after = selectionCounters[i];
    std::string label = "-";
    if(i < (int)labels.size()) {
      label = labels[i];
    }
    *mylog << "BIT:" << std::setw(2) << i << "  " << std::setw(25) << label << "   b:" << std::setw(10) << before << " a:" << std::setw(10) << after << " eff:" << after/before << std::endl;
  }  
  //printSelectionEfficiency(eff);
}
void bcInfos::checkConsistency()
{
  //*mylog << i << " " << label << " " << triggerCounters[i] << " " << mFilterCounters->GetBinContent(i+2) << "||" << selectionCounters[i] << " " << mSelectionCounters->GetBinContent(i+2) << std::endl;
  *mylog << "===> RUN:" << runNumber << " TVXCTP: " << TVXCTP << " TVXCCDB:" << TVXCCDB << " CCDB/TVX:" << (double_t) TVXCCDB/TVXCTP << std::endl;
  for(int i = 0; i < Ndimused; i++) {
    if(selectionCounters[i] != mSelectionCounters->GetBinContent(i+2)){
        std::string label = "-";
        if(i < (int)labels.size()) {
          label = labels[i];
        }
        *mylog << "ERROR:";
        *mylog << std::setw(3) << i << std::setw(25) << label << " TrigBitCnt/FilterHisto:" << std::setw(15) << triggerCounters[i] << "/" << std::setw(15) << mFilterCounters->GetBinContent(i+2);
        *mylog  << " SelBitCnt/SelHisto:" << std::setw(15) << selectionCounters[i] << "/" << mSelectionCounters->GetBinContent(i+2) << std::endl;
    }
    if(selectionCounters[i] == 0) {
      *mylog << "WARN " << i << " " << 0 << std::endl;
    }
  }
}
//
// main
//
void valCCDB_CTP()
{
  ofstream ff;
  ff.open("file.txt");
  mylog = &ff;
  //
  int runNum = 529691;
  // loop over files
  std::string inputList = "list/inputList.txt";
  std::ifstream files(inputList.data());
  std::string file;
  int nfiles = 0;
  int nperiods = 0;
  while (std::getline(files, file) && nperiods == 0) {
    std::cout << "Doing:" << file << std::endl;
    *mylog << "Doing:" << file << std::endl;
    std::ifstream fileruns(file.data());
    std::string runline;
    while (std::getline(fileruns, runline) && nfiles < 100000) {
      std::cout << runline << std::endl;
      *mylog<< runline << std::endl; 
      std::vector<std::string> tokens = o2::utils::Str::tokenize(runline, '/');
      size_t ntokens = tokens.size();
      //std::cout << runline << " tokens:" << ntokens << " " << tokens[4] << std::endl;
      nfiles++;
      runNum = std::stoi(tokens[4]);
      std::cout << runNum << std::endl;
      bcInfos bcis;
      bcis.getDataCCDB(runNum,0,0);
      bcis.extractLabels();
      bcis.checkConsistency();
    }
    nperiods++;
  }
}
