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
  ULong64_t bcAOD, bcEvSel, trigMask[2], selMask[2];
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
  std::array<int, Ndim> selectionCounters{0};  // counting bits from CCDB file
  std::array<int, Ndim> triggerCounters{0};    // counting bits from CCDB file
    std::array<int, Ndim> selectionCountersD{0};  // counting bits from disk file
  std::array<int, Ndim> triggerCountersD{0};    // counting bits from disk file
  std::array<double_t,Ndim> downscaleFactors{0};
  int totalSelected = 0;
  int totalTriggered = 0;
  TH1* mFilterCounters = nullptr;       // from CCDB
  TH1* mSelectionCounters = nullptr;    // from CCDB
  TH1* mFilterCountersD = nullptr;       // from CCDB
  TH1* mSelectionCountersD = nullptr;    // from CCDB
  TH1* mInspectedTVX = nullptr;         // from CCDB
  std::array<double_t,Ndim> eff{0};
  std::array<double_t,Ndim> err{0};
  uint64_t TVXCTP{0},TVXCCDB{0};
  //
  void getData(TFile& inputFile, int Nmax = 0);
  void getDataCCDB(int run, int Nmax = 0, bool print = 0);
  void getHistosFRomAnal();
  double_t getDuration();
  void selectionEfficiency();
  void selectionEfficiencyBiased() const;
  void extractLabels();
  void checkConsistency();
  void checkConsistencyD();
  //
  void printbcInfoVect();
  void printbcInfoVect(std::vector<bcInfo>& bcs);
};
//
// read data from root file
//
void bcInfos::getData(TFile& inputFile, int Nmax)
{
  int i = 0;
  for (auto key : *(inputFile.GetListOfKeys())) {
    TTree* cefpTree = (TTree*)inputFile.Get(Form("%s/selectedBC", key->GetName()));
    if (!cefpTree)
      continue;
    bcInfo bcAO2D;
    cefpTree->SetBranchAddress("bcAO2D", &bcAO2D.bcAOD);
    cefpTree->SetBranchAddress("bcEvSel", &bcAO2D.bcEvSel);
    if(cefpTree->GetBranch("selMask") && cefpTree->GetBranch("triMask")) {
      cefpTree->SetBranchAddress("selMask", &bcAO2D.selMask[0]);
      cefpTree->SetBranchAddress("triMask", &bcAO2D.trigMask[0]);
    } else {
      cefpTree->SetBranchAddress("selMask0", &bcAO2D.selMask[0]);
      cefpTree->SetBranchAddress("triMask0", &bcAO2D.trigMask[0]);
      cefpTree->SetBranchAddress("selMask1", &bcAO2D.selMask[1]);
      cefpTree->SetBranchAddress("triMask1", &bcAO2D.trigMask[1]);
    }
    for (int i = 0; i < cefpTree->GetEntries(); i++) {
      if((i < Nmax) || (Nmax == 0)) {
        cefpTree->GetEntry(i);
        bcs.push_back(bcAO2D);
        // Check consistency
        //if(~bcAO2D.trigMask & bcAO2D.selMask) {
        //  *mylog << "ERROR selMask is not subset of trigMask:";
        //  bcAO2D.print();
        //}
        // Counters
        for (int j = 0; j < Ndimused; j++) {
          int index = j/64;
          int mask = j % 64;
          if (bcAO2D.selMask[index] & (1ull << mask))
            selectionCountersD[j]++;
          if (bcAO2D.trigMask[index] & (1ull << mask))
            triggerCountersD[j]++;
        }
      }
    }
  }
  std::sort(bcs.begin(), bcs.end(), [](const bcInfo& a, const bcInfo& b) { return a.bcEvSel < b.bcEvSel; });
  //
  *mylog << "bcs Size:" << bcs.size() << std::endl;
  for (int i = 0; i < Ndimused; i++) {
    totalSelected += selectionCountersD[i];
    totalTriggered += triggerCountersD[i];
    *mylog << i << " Trig:" << triggerCountersD[i] << " Sel:" << selectionCountersD[i] << std::endl;
  }
  *mylog << "Total original triggers:" << totalTriggered << " selected triggers:" << totalSelected ;
  *mylog << " Duration:" << getDuration() << std::endl;
}
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
  mInspectedTVX->SetDirectory(0);
  //
  path = mCCDBPathOTS+"/FilterCounters";
  mFilterCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mFilterCounters == nullptr ) {
    *mylog << "FilterScalers not found" << std::endl;
    return;
  }
  mFilterCounters->SetDirectory(0);
  int Nfiltered = mFilterCounters->GetNbinsX();
  *mylog << "FilterCounters N bins:" << mFilterCounters->GetNbinsX() << std::endl;
  path = mCCDBPathOTS+"/SelectionCounters";
  mSelectionCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mSelectionCounters == nullptr ) {
    *mylog << "SelectionCounters not found" << std::endl;
    return;
  }
  mSelectionCounters->SetDirectory(0);
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
void bcInfos::getHistosFRomAnal()
{ 
  std::string filename = "AnalysisResults_fullrun.root";
  *mylog << "Extracting labels from file:" << filename << std::endl;
  TFile file(filename.c_str());
  if(!file.IsOpen()) {
    *mylog << "File AnalysisResults.root can not be opned" << std::endl;
    return;
  }
  // Extract the histograms
  TH1* hist1 = dynamic_cast<TH1*>(file.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  if(hist1 == nullptr) {
    *mylog << " Can not find labels" << std::endl;
    return;
  }
  //*mylog << "ptr: " << originalBCs.mFilterCounters << std::endl;
  TH1* hist2 = dynamic_cast<TH1*>(file.Get("central-event-filter-task/scalers/mScalers;1")); // Replace with the correct path
  if(hist2 == nullptr) {
    *mylog << " Can not find labels" << std::endl;
    return;
  }
  mFilterCountersD = hist2;
  mSelectionCountersD = hist1;
  mFilterCountersD->SetDirectory(0);
  mSelectionCountersD->SetDirectory(0);
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
void bcInfos::checkConsistencyD()
{
  for(int i = 0; i < Ndimused; i++) {
    if(selectionCounters[i] != selectionCountersD[i]) {
      *mylog << "sel counters from disk vs ccdb diff" << std::endl;
    }
    if(triggerCounters[i] != triggerCountersD[i]) {
      *mylog << "trig counters from disk vs ccdb diff" << std::endl;
    }
    if(mSelectionCounters->GetBinContent(i+2) != mSelectionCountersD->GetBinContent(i+2)) {
      *mylog << "sel histos from disk vs ccdb diff" << std::endl;
    }
    if(mFilterCounters->GetBinContent(i+2) != mFilterCountersD->GetBinContent(i+2)) {
      *mylog << "trig histos from disk vs ccdb diff" << std::endl;
    }
  }
  *mylog << "checkConsistencyD finished" << std::endl;
}
//
// main
//
void valCCDB_OTS_disk()
{
  ofstream ff;
  ff.open("file.txt");
  mylog = &ff;
  //
  int runNum = 552205;
  //file = "/home/rl/CounterFixes/trigger/valCCDBOTS/list/list_LHC24aj.txt";
  //std::cout << runline << " tokens:" << ntokens << " " << tokens[4] << std::endl;
  std::cout << runNum << std::endl;
  bcInfos bcis;
  bcis.getDataCCDB(runNum,0,0);
  std::string oribcs = "bcRanges_fullrun.root";
  TFile orifile(oribcs.data());
  bcis.getData(orifile);
  bcis.getHistosFRomAnal();
  bcis.extractLabels();
  bcis.checkConsistency();
  bcis.checkConsistencyD();
}
