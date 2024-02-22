#include <iostream>
#include <TFile.h>
#include <TTree.h>

#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"

using o2::InteractionRecord;
using o2::dataformats::IRFrame;
void getData(std::map<uint64_t,uint64_t>& gbc2mask,const char* filename = "bcRangesUnskimmed.root")
{
  TFile inputFile(filename, "READ");
  TTree* treeRanges = nullptr;
  TTree* treeSelectedBC = nullptr;
  int bcRanges = 0;
  int selectedBC = 0;
  std::cout << std::hex;
  for (auto directoryKey : *inputFile.GetListOfKeys()) {
    std::cout << "Processing directory " << directoryKey->GetName() << std::endl;
    auto dirName = directoryKey->GetName();
    treeRanges = dynamic_cast<TTree*>(inputFile.Get(Form("%s/O2bcranges", dirName)));
    treeSelectedBC = dynamic_cast<TTree*>(inputFile.Get(Form("%s/selectedBC", dirName)));
    if (!treeRanges || !treeSelectedBC) {
      std::cerr << "Error: could not find the required trees in directory " << dirName << std::endl;
      continue;
    }  
    //
    // Get the branches we need from the trees
    ULong64_t fBCstart, fBCend;
    treeRanges->SetBranchAddress("fBCstart", &fBCstart);
    treeRanges->SetBranchAddress("fBCend", &fBCend);
    ULong64_t bcAO2D,bcEvSel,triMask,selMask;
    treeSelectedBC->SetBranchAddress("bcAO2D", &bcAO2D);
    treeSelectedBC->SetBranchAddress("bcEvSel", &bcEvSel);
    treeSelectedBC->SetBranchAddress("triMask", &triMask);
    treeSelectedBC->SetBranchAddress("selMask", &selMask);
    // Loop over the entries in the ranges tree and check if the BC range is valid
    int nEntriesRanges = treeRanges->GetEntries();
    bcRanges += nEntriesRanges;
    for (int iEntryRanges = 0; iEntryRanges < nEntriesRanges; ++iEntryRanges) {
      treeRanges->GetEntry(iEntryRanges);
      //std::cout << "start:0x" << fBCstart << " end:0x" << fBCend << " " << std::dec << fBCend - fBCstart << std::hex << std::endl; 
    }
    int nEntriesSelectedBC = treeSelectedBC->GetEntries();
    selectedBC += nEntriesSelectedBC;
    for (int iEntrySelectedBC = 0; iEntrySelectedBC < nEntriesSelectedBC; ++iEntrySelectedBC) {
      treeSelectedBC->GetEntry(iEntrySelectedBC);
      //std::cout << "bcAO2D:0x" << bcAO2D << " bcEvSel:0x" << bcEvSel << " " << std::dec << (int64_t)bcAO2D - (int64_t) bcEvSel << std::hex;
      //std::cout << " 0x" << triMask << " 0x" << selMask << std::endl;
      gbc2mask[bcEvSel] = selMask; 
    }
  }
}
void eff()
{
  std::map<uint64_t,uint64_t> trigsUnskimmed;	
  getData(trigsUnskimmed);
  std::map<uint64_t,uint64_t> trigsSkimmed;	
  getData(trigsSkimmed,"bcRangesSkimmed.root");
  //
  std::cout << std::dec;
  std::cout << "US:" << trigsUnskimmed.size() << std::endl;
  std::cout << "SK:" << trigsSkimmed.size() << std::endl;
  std::cout << "====================" << std::endl;
  std::map<uint64_t,uint64_t>::iterator it;
  for(int i = 0; i < 64; i++) {  
  uint64_t mask = 1ull<<i;
  std::map<uint64_t,uint64_t>::iterator bg = trigsUnskimmed.begin();	
  std::vector<int> bc1;
  std::map<uint64_t,int> gbc1;
  for(it = bg; it != trigsUnskimmed.end(); it++) {
    if(it->second & mask) {
      InteractionRecord ir;
      ir.setFromLong(it->first);
      bc1.push_back(ir.bc);
      gbc1[it->first] = ir.bc;
    }
  }
  std::vector<int> bc2;
  std::map<uint64_t,int> gbc2;
  bg = trigsSkimmed.begin();	
  for(it = bg; it != trigsSkimmed.end(); it++) {
    if(it->second & mask) {
      InteractionRecord ir;
      ir.setFromLong(it->first);
      bc2.push_back(ir.bc);
      gbc2[it->first] = ir.bc;
    }
  }
  std::vector<int> bcfound;
  std::vector<int> bcfound1;
  std::vector<int> bcfoundm1;
  for(auto const& bc: gbc1) {
    if(gbc2.find(bc.first) != gbc2.end()){
      bcfound.push_back(bc.second);
    }
    for(int del = 1; del < 0; del++) {
      if(gbc2.find(bc.first-del) != gbc2.end()){
        bcfoundm1.push_back(bc.second);
      }
      if(gbc2.find(bc.first+del) != gbc2.end()){
        bcfound1.push_back(bc.second);
      }
    }
  }
  //
  std::cout << "=====>" << i << "<== US:" << bc1.size() << " SK:" << bc2.size() << std::endl;
  std::cout << "US: "; 
  for(auto const& bc: bc1) std::cout << bc << " ";
  std::cout << std::endl;
  std::cout << "SK: "; 
  for(auto const& bc: bc2) std::cout << bc << " ";
  std::cout << std::endl;
  std::cout << "Found: "; 
  for(auto const& bc: bcfound) std::cout << bc << " ";
  std::cout << std::endl;
  if(bcfoundm1.size()+bcfound1.size()) {
    std::cout << "!!! ";
    std::cout << "EFF:  " <<bcfoundm1.size() << " " << bcfound.size() << " " << bcfound1.size() << " / " << bc1.size() << std::endl;
  } else {
    std::cout << "EFF bit:" << i << " " << bcfound.size() << " / " << bc1.size() << " = " << (float)bcfound.size()/bc1.size() << std::endl;
  }
  }
}
