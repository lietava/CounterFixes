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
#include <vector>
#include <map>
#include "TGrid.h"
#include "TSystem.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TTree.h"
#include "CommonUtils/StringUtils.h"
const int Ndim = 128;
struct bcInfo {
  ULong64_t bcAOD, bcEvSel, trigMask[2], selMask[2];
};
struct bcInfos {
  bcInfo bci;
  TH1* scalers = nullptr;
  TH1* filters = nullptr;
  std::array<int, Ndim> selectionCounters{0};  // counting bits from disk or CCDB file
  std::array<int, Ndim> triggerCounters{0};    // counting bits from disk or CCDB file
  int Ndimused = 0;
  std::vector<std::string> labels;
  void readbcInfos(std::string arfile, std::string bcrfile, bool useAlien = 1);
};
struct bcInfos2 {
  bcInfos unskimmed;
  bcInfos skimmed;
};
int getRunToFileMap(std::string list, std::map<int,std::string>& run2file)
{
  std::ifstream files(list.data());
  if(!files.is_open()) {
    std:;cout << "Can not open file:" << list << std::endl;
    return 1;
  }
  std::string file;
  while (std::getline(files, file)) {
    std::vector<std::string> tokens = o2::utils::Str::tokenize(file, '/');
    size_t ntokens = tokens.size();
    //std::cout << file << " tokens:" << ntokens << " " << std::endl;
    if(ntokens < 5) {
      std::cout << "skipping:" << file << std::endl;
      continue;
    }
    int runNum = std::stoi(tokens[4]);
    //std::cout << runNum << std::endl;
    if(run2file.count(runNum)) {
      std::cout << "ERROR:" << runNum << " twice in " << list << std::endl; 
    } else {
      run2file[runNum] = file;
    }
  }
  return 0;
}
TH1F* addHistWithNames(std::string const& name, std::vector<std::string>& labels, int Ndimh)
{
  TH1F *h = new TH1F(name.c_str(), name.c_str(), Ndimh, 0, Ndimh);
  for (int i = 0; i < h->GetXaxis()->GetNbins(); i++) {
    h->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  h->GetXaxis()->SetLabelSize(0.025);
  h->GetXaxis()->LabelsOption("v");
  return h;
}
//
void bcInfos::readbcInfos(std::string arfile, std::string bcrfile, bool useAlien)
{
    if (useAlien) {
    TGrid::Connect("alien://");
  }
  std::string file;
  file = useAlien ? "alien://" + arfile : arfile;
  std::unique_ptr<TFile> scalersFile{TFile::Open((file).data(), "READ")};
  file = useAlien ? "alien://" + bcrfile : bcrfile;
  std::unique_ptr<TFile> bcrFile{TFile::Open((file).data(), "READ")};
  scalers = (TH1*)scalersFile->Get("central-event-filter-task/scalers/mScalers");
  filters = (TH1*)scalersFile->Get("central-event-filter-task/scalers/mFiltered");
  for (int i = 1; i <= scalers->GetNbinsX(); i++) {
    std::string label = scalers->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
    }
    }
  //
  Ndimused = scalers->GetNbinsX() - 2;
  std::vector<bcInfo> bcs;
    for (auto key : *(bcrFile->GetListOfKeys())) {
    TTree* cefpTree = (TTree*)bcrFile->Get(Form("%s/selectedBC", key->GetName()));
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
      cefpTree->GetEntry(i);
      bcs.push_back(bcAO2D);
      // Check consistency
      if(~bcAO2D.trigMask[0] & bcAO2D.selMask[0]) {
          std::cout << "ERROR selMask is not subset of trigMask [0]";
      }
      if(~bcAO2D.trigMask[1] & bcAO2D.selMask[1]) {
          std::cout << "ERROR selMask is not subset of trigMask [1]";
      }
      // Counters
      for (int j = 0; j < Ndimused; j++) {
        int index = j/64;
        int mask = j % 64;
        if (bcAO2D.selMask[index] & (1ull << mask))
          selectionCounters[j]++;
        if (bcAO2D.trigMask[index] & (1ull << mask))
          triggerCounters[j]++;
      }
    }
  }
  bcrFile->Close();
  scalersFile->Close();
}
//
/*
void efficiency(TFile* outputFile = nullptr, std::vector<std::string> labels) 
{
  gStyle->SetOptStat(0);
  // downscales
  std::string hname = std::to_string(run) + " DS factors";
  TH1F* hds = addHistWithNames(hname, labels, Ndimused);
  //TH1F hds(hname.c_str(), hname.c_str(), Ndimused, 0, Ndimused);
  std::string hnamese = std::to_string(run) + " Skimming eff";
  TH1F* hse = addHistWithNames(hnamese, labels, Ndimused);
  for(int i = 0; i < Ndimused; i++) {
    double_t ds = 0;
    if(scalers->GetBinContent(i+2)) {
      ds = filters->GetBinContent(i+2)/scalers->GetBinContent(i+2);
    }
    double_t ef = 0;
    if(filters->GetBinContent(i+2)) {
      ef = triggerCounters[i]/filters->GetBinContent(i+2);
    }
    hds->SetBinContent(i+1,ds);
    hse->SetBinContent(i+1,ef);
  }
  //hds.Write();
  outputFile->WriteObject(hds,hname.c_str());
  outputFile->WriteObject(hse,hnamese.c_str());
  // eff fact = bcRanges count/selected - biased unless ds=1 for all skimmedEfficiency
}
*/
//
// unskimmed bcRange, skimmed bcRange
// 
//
void calcSkimEfficiency(std::string periodName = "LHC24aj")
{
  std::map<int,std::string> armap_us,armap_sk;
  std::map<int,std::string> bcrmap_us,bcrmap_sk;
  int year = 2000 + std::stoi(periodName.substr(3, 2));
  std::cout << "year:" << year << " period:" << periodName << std::endl;
  // unskimmed bcRanges
  std::string bcrlist_us = "bcrlist_unskim_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ bcRanges_fullrun.root | sed '/skim/d'> %s", year, periodName.data(), bcrlist_us.data()));
  // skimmed bcRanges
  std::string bcrlist_sk = "bcrlist_skim_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ ctf_skim_full/bcRanges_fullrun.root > %s", year, periodName.data(), bcrlist_sk.data()));
  // unskimmed AR
  std::string arlist_us = "arlist_unskim_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ AnalysisResults_fullrun.root | sed '/skim/d' | sed '/Calib/d' > %s", year, periodName.data(), arlist_us.data()));
  // skimmed AR
  std::string arlist_sk = "arlist_skim_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ AnalysisResults_fullrun.root | sed '/skim/d' | sed '/Calib/d' > %s", year, periodName.data(), arlist_sk.data()));
  //
  getRunToFileMap(bcrlist_us,bcrmap_us);
  getRunToFileMap(bcrlist_sk,bcrmap_sk);
  getRunToFileMap(arlist_us,armap_us);
  getRunToFileMap(arlist_sk,armap_sk);
  //
  if(bcrmap_us.size() != bcrmap_us.size()) {
    std::cout << "bcr unskimmed run list != bcr skimmed run list" << std::endl;
  }
  if(armap_us.size() != armap_sk.size()) {
    std::cout << "ar unskimmed run list != ar skimmed run list" << std::endl;
  }
  if(bcrmap_us.size() != armap_us.size()) {
    std::cout << "bcr unskimmed run list != ar unskimmed run list" << std::endl;
  }
  //
  int NFILES = 2;
  int nf = 0;
  for(auto const& run: armap_us) {
    bcInfos2 data;
    std::cout << "Run:" << run.first << std::endl;
    if( bcrmap_us.count(run.first)) {
      data.unskimmed.readbcInfos(run.second,bcrmap_us[run.first]);
    } else {
      std::cout << "Run:" << run.first << " not in bcRanges us" << std::endl;
    }
    if( bcrmap_sk.count(run.first) && armap_sk.count(run.first)) {
      data.skimmed.readbcInfos(armap_sk[run.first],bcrmap_sk[run.first]);
    } else {
      std::cout << "Run:" << run.first << " not in ar or bcRanges sk" << std::endl;
    }
    nf++;
    if( nf == NFILES ) {
      break;
    }
  }
  return;
  //
  //
  TFile* outputFile = new TFile("output.root", "RECREATE");
  outputFile->Close();
}
