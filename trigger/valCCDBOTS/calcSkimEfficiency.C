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
#include <map>
#include "TGrid.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "CommonUtils/StringUtils.h"
const int Ndim = 128;
struct bcInfo {
  ULong64_t bcAOD, bcEvSel, trigMask[2], selMask[2];
};
int getRunToFileMap(std::string list, std::map<int,std::string>& run2file)
{
  std::ifstream files(list.data());
  if(!files.is_open()) {
    std:;cout << "Can not open file:" << list << std::endl;
    return 1;
  }
  int nf = 0;
  std::string file;
  while (std::getline(files, file) && nf < 3) {
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
    nf++;
  }
  return 0;
}
void efficiency(int run, std::string arfile, std::string bcrfile, bool useAlien = 1) 
{
  if (useAlien) {
    TGrid::Connect("alien://");
  }
  std::string file;
  file = useAlien ? "alien://" + arfile : arfile;
  std::unique_ptr<TFile> scalersFile{TFile::Open((file).data(), "READ")};
  file = useAlien ? "alien://" + bcrfile : bcrfile;
  std::unique_ptr<TFile> bcrFile{TFile::Open((file).data(), "READ")};
  TH1* scalers = (TH1*)scalersFile->Get("central-event-filter-task/scalers/mScalers");
  TH1* filters = (TH1*)scalersFile->Get("central-event-filter-task/scalers/mFiltered");
  //
  //int Ndimused = 64;
  int Ndimused = scalers->GetNbinsX() - 2;
  std::array<int, Ndim> selectionCounters{0};  // counting bits from disk or CCDB file
  std::array<int, Ndim> triggerCounters{0};    // counting bits from disk or CCDB file
  std::array<double_t,Ndim> downscaleFactors{0};
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
      //if(~bcAO2D.trigMask & bcAO2D.selMask) {
      //  *mylog << "ERROR selMask is not subset of trigMask:";
      //  bcAO2D.print();
      //}
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
}
void calcSkimEfficiency(std::string ARList = "ARList.txt", std::string bcRList = "bcRList.txt", std::string periodName = "LHC24aj")
{
  std::map<int,std::string> armap;
  std::map<int,std::string> bcrmap;
  int year = 2000 + std::stoi(periodName.substr(3, 2));
  std::cout << "year:" << year << " period:" << periodName << std::endl;
  std::string bcrlist = "bcrlist_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ ctf_skim_full/bcRanges_fullrun.root > %s", year, periodName.data(), bcrlist.data()));
  std::string arlist = "arlist_"+periodName+".txt";
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ AnalysisResults_fullrun.root | sed '/skim/d' | sed '/Calib/d' > %s", year, periodName.data(), arlist.data()));
  getRunToFileMap(bcrlist,bcrmap);
  getRunToFileMap(arlist,armap);
  int NFILES = 1;
  int nf = 0;
  for(auto const& run: armap) {
    if( bcrmap.count(run.first)) {
      efficiency(run.first,run.second,bcrmap[run.first]);
    } else {
      std::cout << "Run:" << run.first << " not in bcRanges" << std::endl;
    }
    nf++;
    if( nf == NFILES ) {
      break;
    }
  }
}
