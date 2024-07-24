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
// O2 includes// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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

#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <regex>

void checkConsistency(std::string anFileName, std::string bcRangesFile, TFile* outputFile = nullptr)
{
  gStyle->SetOptStat(0);
  std::string runNumber = "";
  std::regex re("/5[0-9]*");
  std::smatch match;
  if (std::regex_search(anFileName, match, re)) {
    // Remove the leading '/'
    runNumber = match.str().substr(1);
  }

  // Load the root files
  TFile file1(anFileName.c_str());
  TFile rangeFile(bcRangesFile.c_str());

  // Extract the histograms
  TH1* hist0 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mScalers;1"));
  TH1* hist1 = dynamic_cast<TH1*>(file1.Get("central-event-filter-task/scalers/mFiltered;1"));

  if (!hist0 || !hist1) {
    std::cerr << "Error: Failed to extract histograms from the root files." << std::endl;
    return;
  }

  std::vector<std::string> labels;
  std::vector<int> bits;
  std::vector<double> selected_bins1;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    // if (label != "Total number of events" && label != "Filtered events" && hist0->GetBinContent(i) == hist1->GetBinContent(i) && hist0->GetBinContent(i) > 0) {
      labels.push_back(label);
      selected_bins1.push_back(hist1->GetBinContent(i));
      bits.push_back(i - 2);
    // }
  }
  TH1D hOriginal("hOriginal", "AnalysisResults;;Number of events", labels.size(), 0, labels.size());       // Histogram for the original values
  TH1D hFromBCRange("hFromBCRange", "bcRanges;;Number of events", labels.size(), 0, labels.size()); // Histogram for the skimmed values
  TH1D hRatio("hRatio", (runNumber + ";;bcRanges / AnalysisResults").data(), labels.size(), 0, labels.size());                        // Histogram for the ratio of the two

  // Fill the histograms
  for (size_t i = 0; i < labels.size(); i++) {
    hOriginal.SetBinContent(i + 1, selected_bins1[i]);
    hOriginal.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hFromBCRange.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hRatio.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  int nEvents = std::round(hist1->GetBinContent(hist1->GetNbinsX()));
  int numSelectedBCs = 0;
  double testPrecision = 0.;

  ULong64_t bcStart{0ull}, bcEnd{0ull}, bcAO2D{0ull}, bcEvSel{0ull}, triMask[2]{0ull}, selMask[2]{0ull};

  for (auto key : *rangeFile.GetListOfKeys()) {
    auto dir = dynamic_cast<TDirectory*>(rangeFile.Get(key->GetName()));
    if (!dir) {
      continue;
    }
    TTree* selectedBC = dynamic_cast<TTree*>(dir->Get("selectedBC"));
    selectedBC->SetBranchAddress("bcAO2D", &bcAO2D);
    selectedBC->SetBranchAddress("bcEvSel", &bcEvSel);
    if (selectedBC->GetBranch("triMask")) {
      selectedBC->SetBranchAddress("triMask", &triMask[0]);
      selectedBC->SetBranchAddress("selMask", &selMask[0]);
    } else {
      selectedBC->SetBranchAddress("triMask0", &triMask[0]);
      selectedBC->SetBranchAddress("triMask1", &triMask[1]);
      selectedBC->SetBranchAddress("selMask0", &selMask[0]);
      selectedBC->SetBranchAddress("selMask1", &selMask[1]);
    }
    numSelectedBCs += selectedBC->GetEntries();
    for (int i = 0; i < selectedBC->GetEntries(); i++) {
      selectedBC->GetEntry(i);
      for (uint32_t ibits{0}; ibits < bits.size(); ++ibits) {
        int maskId = bits[ibits] / 64;
        int maskBit = bits[ibits] % 64;
        if (selMask[maskId] & (1ull << maskBit)) {
          hFromBCRange.Fill(ibits + 1.e-6);
          hRatio.Fill(ibits + 1.e-6);
        }
      }
      testPrecision += bool(selMask[0] | selMask[1]);
    }
  }
  hRatio.Divide(&hOriginal);

  if (nEvents != int(std::round(testPrecision))) {
    std::cerr << "[FATAL]: merging failed, number of events seen by CEFF (" << nEvents << ") and stored selected BCs (" << int(std::round(testPrecision)) << ") do not match! Run " << runNumber << std::endl;
  }

  bool localFile = !outputFile;
  if (!outputFile) {
    outputFile = new TFile("output.root", "RECREATE");
  }
  if (!runNumber.empty()) {
    outputFile->mkdir(runNumber.data());
  }
  outputFile->cd(runNumber.data());

  hOriginal.Write();
  hFromBCRange.Write();
  hRatio.Write();

  if (localFile) {
    outputFile->Close();
  }
}

void checkConsistency2(std::string listName = "period.txt") {
  std::string periodName = listName.substr(0, listName.find_first_of('.'));
  std::ifstream file(listName);
  std::string line;
  TFile* outputFile = new TFile((periodName + ".root").data(), "RECREATE");
  TCanvas c1("c1", "c1", 800, 600);
  c1.SetGridx();
  c1.SetGridy();
  int lineCounter = 0;
  while (std::getline(file, line)) {
    lineCounter++;
  }
  file.clear();
  file.seekg(0, std::ios::beg);
  int counter = 0;
  while (std::getline(file, line)) {
    size_t pos = line.find(",");
    std::string anFileName = line;
    std::string path = line.substr(0, line.find_last_of('/') + 1);
    std::string bcRangesFile = path + "bcRanges_fullrun.root";
    checkConsistency(anFileName, bcRangesFile, outputFile);
    TH1* hRatio = (TH1*)gDirectory->Get("hRatio");
    hRatio->Draw();
    std::string suffix = counter == 0 ? "(" : (counter == lineCounter - 1 ? ")" : "");
    c1.Print((periodName + ".pdf" + suffix).data());
    counter++;
  }
  outputFile->Close();
}
void checkConsistency()
{
  std::string periodName = "535069";
  TFile* outputFile = new TFile((periodName + ".root").data(), "RECREATE");
  std::string anFileName = "AnalysisResults_fullrun.root";
  std::string bcRangesFile = "bcRanges_fullrun-skimmed.root";
  checkConsistency(anFileName, bcRangesFile, outputFile);
}
