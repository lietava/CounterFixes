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
//
//
//================================
//

// bcInfos
//
struct bcInfos
{
  bcInfos() = default;
  TH1* mFilterCounters = nullptr;       // from CCDB
  TH1* mSelectionCounters = nullptr;    // from CCDB
  //
};

//
//===========================
//
struct effUtils
{
  effUtils() = default;
  bcInfos originalBCs, skimmedBCs;
  //
  void extractLabelsAnal(std::vector<std::string>& labels);
  void printTrigsAndSels();
};

void effUtils::extractLabelsAnal(std::vector<std::string>& labels)
{
  //TFile file("AnalysisResults.root");
  std::cout << "Extracting labels from file" << std::endl;
  TFile file("AnalysisResults_fullrun.root");
  if(!file.IsOpen()) {
    std::cout << "File AnalysisResults.root can not be opned" << std::endl;
    return;
  }
  // Extract the histograms
  TH1* h = dynamic_cast<TH1*>(file.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  if(h == nullptr) {
    std::cout << " Can not find labels" << std::endl;
    return;
  }
  std::cout << "size:" << sizeof(*h) << std::endl;
  originalBCs.mFilterCounters = (TH1*) h->Clone();
  originalBCs.mFilterCounters->SetDirectory(0);
  std::cout << "size:" << sizeof(*originalBCs.mFilterCounters) << std::endl;

  std::cout << "===> " << originalBCs.mFilterCounters->GetBinContent(2) << std::endl;
  std::cout << "ptr: " << originalBCs.mFilterCounters << std::endl;
  //return;

  std::cout << "Exctracting labels done" << std::endl;
}
void effUtils::printTrigsAndSels()
{

  for(int  i = 0; i < 64; i++) {
    std::cout << "ptr: " << originalBCs.mFilterCounters << std::endl;
    TH1F* h = (TH1F*) originalBCs.mFilterCounters;
    if(originalBCs.mFilterCounters != nullptr) {
      std::cout << i <<  " OFC:" << h->GetBinContent(i); ;
    }
    if(originalBCs.mSelectionCounters) {
     //std::cout << " OSC:" << originalBCs.mSelectionCounters->GetBinContent(i+2);
    }
    if(skimmedBCs.mFilterCounters) {
      //std::cout << " SFC:" << skimmedBCs.mFilterCounters->GetBinContent(i+2);
    }
    if(skimmedBCs.mSelectionCounters) {
      //std::cout << "SSC:" << skimmedBCs.mSelectionCounters->GetBinContent(i+2);
    }
    std::cout << std::endl;
  }
}
//
// main
//
void eff3H(std::string original = "bcRanges_fullrun.root", std::string skimmed = "bcRanges_fullrun-skimmed.root")
{
  effUtils eff;

  std::vector<std::string> labs;
  eff.extractLabelsAnal(labs);
/*
  TFile file("AnalysisResults_fullrun.root");
  if(!file.IsOpen()) {
    std::cout << "File AnalysisResults.root can not be opned" << std::endl;
    return;
  }
  // Extract the histograms
  eff.originalBCs.mFilterCounters = dynamic_cast<TH1*>(file.Get("central-event-filter-task/scalers/mFiltered;1")); // Replace with the correct path
  if(eff.originalBCs.mFilterCounters == nullptr) {
    std::cout << " Can not find labels" << std::endl;
    return;
  }
  std::cout << "===> " << eff.originalBCs.mFilterCounters->GetBinContent(2) << std::endl;
  std::cout << "ptr: " << eff.originalBCs.mFilterCounters << std::endl;
*/
  eff.printTrigsAndSels();
  return;
}
