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
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include "CCDB/BasicCCDBManager.h"
//
#pragma link C++ class std::vector < std::array < uint64_t, 2>> + ;
//
std::vector<std::string> labels = {"fPHOSnbar", "fPHOSPair", "fPHOSElectron", "fPHOSPhoton", "fOmegaLargeRadius", "fSingleXiYN", "fQuadrupleXi", "fhadronXi", "fTripleXi", "fGammaHighPtDCAL", "fGammaHighPtEMCAL", "fJetFullHighPt", "fLD", "fPD", "fLLL", "fPLL", "fPPL", "fPPP", "fHighFt0cFv0Flat", "fHighFt0cFv0Mult", "fHighFt0Flat", "fHighFt0Mult", "fHighTrackMult", "fHfDoubleCharmMix", "fHfDoubleCharm3P", "fHfSoftGamma3P", "fHfFemto2P", "fHfBeauty4P", "fHfFemto3P", "fHfBeauty3P", "fHfSoftGamma2P", "fHfDoubleCharm2P", "fDiMuon", "fDiElectron", "fUDdiff", "fHe"};
//
int runNUMBER = 2240;
const int Ndim = 64;
const int Nfreq = 10;
const int Ndim_used = 64;
ofstream* mylog;
//
template <class T, std::size_t N>
ostream& operator<<(ostream& o, const array<T, N>& arr)
{
  copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, " "));
  return o;
}
template <class T>
ostream& operator<<(ostream& o, const std::vector<T>&arr)
{
  copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, " "));
  return o;
}
void printMap(std::map<uint64_t,int>& m)
{
  for(auto const& i : m) {
    *mylog << i.first << " " << i.second << std::endl;
  }
};
//
// main structure
//
struct bcInfo {
  bcInfo() = default;
  ULong64_t bcAOD, bcEvSel, trigMask, selMask;
  void print() const;
};
void bcInfo::print() const
{
  *mylog << bcEvSel << " " << bcAOD << " " << std::hex << trigMask << " " << selMask << " " << std::dec;
  for(size_t i = 0; i < Ndim; i++) {
    ULong64_t mask = 1ull << i;
    if(mask & trigMask) {
      if(i < labels.size()) {
        *mylog << labels[i] << " ";
      } else {
        *mylog << "XXX ";
      }
    }
  }
  *mylog << std::endl;
};
struct bcInfos
{
  const std::string ccdbTest = "http://ccdb-test.cern.ch:8080";
  const std::string ccdbProd = "http://alice-ccdb.cern.ch";
  const std::string mCCDBPathOTSmp= "Users/m/mpuccio/EventFiltering/OTS";
  const std::string mCCDBPathOTS= "Users/r/rlietava/EventFiltering/OTS";
  //const std::string mCCDBPathTrigScalersSkimmed = "Users/m/mpuccio/EventFiltering/OTS/SelectionCounters";
  bcInfos() = default;
  std::vector<bcInfo> bcs;
  std::vector<bcInfo> bcs_cleaned;
  std::array<int, Ndim> selectionCounters{0}, triggerCounters{0};
  int totalSelected = 0;
  int totalTriggered = 0;
  std::array<bool,Ndim> zeros{0};  // Take care of empty or not used bits
  std::map<uint64_t,std::vector<bcInfo>> evsel2aod;
  TH1* mFilterCounters = nullptr;
  TH1* mSelectionCounters = nullptr;
  bool skimmed = 0;
  std::array<double_t,Ndim> eff{0};
  //
  void getTrigScalers(int runNumber = runNUMBER);
  uint64_t dataSize(uint64_t window);
  void getData(TFile& inputFile, int Nmax = 0);
  void getDataCCDB(int run);
  double_t getDuration();
  void frequencyBC() const;
  void frequencyBCSorted(std::array<int,Nfreq>& freq);
  void frequencyCol() const;
  void evSel2AOD();
  //void AOD2evSel();
  void evSel2AODInverse(std::vector<bcInfo>& bcs);
  void cleanBC(std::vector<bcInfo>& bcs_cleaned);
  void cleanBC();
  void selectionEfficiency( std::array<double,Ndim>& arr,std::array<double,Ndim>& err ) const;
  void selectionEfficiencyBiased() const;
  void corMatrix();
  void getArrayForBit(std::vector<uint64_t>& bctrigs, std::vector<uint64_t>& bcselec, int bit);
  //
  void printbcInfoVect();
  void printbcInfoVect(std::vector<bcInfo>& bcs);
  void printEvSel2AOD() const;
  void printSelectionEfficiency(std::array<double_t, Ndim>& eff) const;
  void printCorMatrix(float_t (*mat)[Ndim][Ndim]) const;
  void printCorMatrixSimple(float_t (*mat)[Ndim][Ndim]) const;
};
//
void bcInfos::getTrigScalers(int runNumber)
{
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  //ccdbMgr.setURL(ccdbProd);
  //
  auto soreor = ccdbMgr.getRunDuration(runNumber);
  uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  //timeStamp = 1688349854935;
  //
  std::map<std::string,std::string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  *mylog << "run:" << metadata["runNumber"] << " Timestamp:" << timeStamp << std::endl;
  //mFilterCounters = ccdbMgr.getSpecific<TH1>(mCCDBPathTrigScalers, timeStamp, metadata);
  if( mFilterCounters == nullptr ) {
    *mylog << "FilterCounters not found" << std::endl;
    return;
  }
  *mylog << "FilterCounters N bins:" << mFilterCounters->GetNbinsX() << std::endl;
  //trigScalers->Draw();
}
//
// estimate data size in unit of bcs
//
uint64_t bcInfos::dataSize(uint64_t window)
{
  uint64_t start = bcs.front().bcEvSel;
  uint64_t maxsize = bcs.back().bcEvSel - start;
  *mylog << "Max size:" <<  maxsize << std::endl;
  uint64_t size = 0;
  auto it = begin(bcs);
  ++it;
  while(it != end(bcs)){
    if( (it->bcEvSel - start) > window) {
      size += window;
    } else {
      size += it->bcEvSel - start;
    }
    start = it->bcEvSel;
    ++it;
  }
  uint64_t lastsize = bcs.back().bcEvSel - start;
  if(lastsize > window) {
    lastsize = window;
  }
  size += lastsize;
  *mylog << "Selection window:" << window << " data ssize in BC:" << size << " fraction of max:" << (double_t)size/(double_t) maxsize <<  std::endl;
  return size;
}
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
    cefpTree->SetBranchAddress("selMask", &bcAO2D.selMask);
    cefpTree->SetBranchAddress("triMask", &bcAO2D.trigMask);
    for (int i = 0; i < cefpTree->GetEntries(); i++) {
      if((i < Nmax) || (Nmax == 0)) {
        cefpTree->GetEntry(i);
        bcs.push_back(bcAO2D);
        // Check consistency
        if(~bcAO2D.trigMask & bcAO2D.selMask) {
          *mylog << "ERROR selMask is not subset of trigMask:";
          bcAO2D.print();
        }
        // Counters
        for (ULong64_t j = 0; j < Ndim; j++) {
          if (bcAO2D.selMask & (1ull << j))
            selectionCounters[j]++;
          if (bcAO2D.trigMask & (1ull << j))
            triggerCounters[j]++;
        }
      }
    }
  }
  std::sort(bcs.begin(), bcs.end(), [](const bcInfo& a, const bcInfo& b) { return a.bcEvSel < b.bcEvSel; });
  //
  *mylog << "bcs Size:" << bcs.size() << std::endl;
  for (int i = 0; i < Ndim; i++) {
    totalSelected += selectionCounters[i];
    totalTriggered += triggerCounters[i];
  }
  *mylog << "Total original triggers:" << totalTriggered << " selected triggers:" << totalSelected ;
  *mylog << " Duration:" << getDuration() << std::endl;
}
void bcInfos::getDataCCDB(int runNumber)
{
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  //ccdbMgr.setURL(ccdbProd);
  //
  auto soreor = ccdbMgr.getRunDuration(536663);
  uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  //timeStamp = 1688349854935;
  //
  std::map<std::string,std::string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  *mylog << "run:" << metadata["runNumber"] << " Timestamp:" << timeStamp << std::endl;
  std::string path = mCCDBPathOTS+"/FilterCounters";
  mFilterCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mFilterCounters == nullptr ) {
    *mylog << "FilterScalers not found" << std::endl;
    return;
  }
  *mylog << "FilterCounters N bins:" << mFilterCounters->GetNbinsX() << std::endl;
  path = mCCDBPathOTS+"/SelectionCounters";
  mSelectionCounters = ccdbMgr.getSpecific<TH1>(path, timeStamp, metadata);
  if( mSelectionCounters == nullptr ) {
    *mylog << "SelectionCounters not found" << std::endl;
    return;
  }
  *mylog << "SelectionCounters N bins:" << mSelectionCounters->GetNbinsX() << std::endl;
  // bcinfo
  path = mCCDBPathOTS+"/FilterBitMasks";
  std::vector<array<uint64_t,2>>* fbm = ccdbMgr.getSpecific<std::vector<array<uint64_t,2>>>(path, timeStamp, metadata);
  *mylog << " FBM size:" << (*fbm).size() << std::endl;
  path = mCCDBPathOTS+"/SelectionBitMasks";
  std::vector<array<uint64_t,2>>* sbm = ccdbMgr.getSpecific<std::vector<array<uint64_t,2>>>(path, timeStamp, metadata);
  *mylog << " SBM size:" << (*sbm).size() << std::endl;
  path = mCCDBPathOTS+"/SelectedBCs";
  std::vector<array<uint64_t,2>>* sbc = ccdbMgr.getSpecific<std::vector<array<uint64_t,2>>>(path, timeStamp, metadata);
  *mylog << " SBC size:" << (*sbc).size() << std::endl;
  if((*fbm).size() != (*sbm).size()) {
    *mylog << "ERROR: FBM size != SBM size" << std::endl;
    return;
  }
  if((*fbm).size() != (*sbc).size()) {
    *mylog << "ERROR: FBM size != SBC size" << std::endl;
    return;
  }
  for(size_t i = 0; i < (*fbm).size(); i++) {
    bcInfo bci;
    bci.bcAOD = (*sbc)[i][0];
    bci.bcEvSel = (*sbc)[i][1];
    bci.trigMask = (*fbm)[i][0];
    bci.selMask = (*sbm)[i][0];
    bcs.push_back(bci);
  }
}
double_t bcInfos::getDuration()
{
  uint64_t tt =   bcs.back().bcEvSel - bcs.front().bcEvSel;
  double_t t = tt*25.e-9;
  return t;
}
void bcInfos::printbcInfoVect(std::vector<bcInfo>& bcs)
{
  for(auto const& i: bcs) {
    i.print();
  }
};
void bcInfos::printbcInfoVect()
{
  printbcInfoVect(bcs);
};
//
// utils
//
void bcInfos::printEvSel2AOD() const
{
  for(auto const& i: evsel2aod) {
    *mylog << "===>evSel BC:" << i.second.size() << " " << i.first << ":" << std::endl;
    for(auto const& j: i.second){
      j.print();
    }
  }
};
//
// Frequency: # of Cols with 1,2,3,... collision TVXs
//
void bcInfos::frequencyCol() const
{
  std::map<uint64_t, int> bcFreq;
  for(auto const& bc: bcs) {
    if( bcFreq.find(bc.bcAOD) != bcFreq.end()) {
      bcFreq[bc.bcAOD] += 1;
    } else {
      bcFreq[bc.bcAOD] = 0;
    }
  }
  //printMap(bcFreq);
  //
  std::array<int,Nfreq> freq{0};
  for(auto const& bc: bcFreq) {
    if(bc.second < 15) {
      freq[bc.second]++;
    }  else  {
      freq[Nfreq-1]++;
    }
  }
  float sum = 0;
  float i = 1;
  
  *mylog << "Freq of bcAOD:" << std::endl;
  for(auto const& n: freq) {
    *mylog << n << " ";
    sum += n*i;
    i += 1;
  }
  *mylog << std::endl;
  *mylog << "Fraction of singles:" << (float)freq[0]/bcs.size() << " sum:" << sum << " size:"  << bcs.size() << std::endl;
}

//
// Frequency: # of TVX with 1,2,3,... collision BCs
//
void bcInfos::frequencyBC() const
{
  std::map<uint64_t, int> bcFreq;
  for(auto const& bc: bcs) {
    if( bcFreq.find(bc.bcEvSel) != bcFreq.end()) {
      bcFreq[bc.bcEvSel] += 1;
    } else {
      bcFreq[bc.bcEvSel] = 0;
    }
  }
  //printMap(bcFreq);
  //
  std::array<int,Nfreq> freq{0};
  for(auto const& bc: bcFreq) {
    if(bc.second < 15) {
      freq[bc.second]++;
    }  else  {
      freq[Nfreq-1]++;
    }
  }
  float sum = 0;
  float i = 1;
  for(auto const& n: freq) {
    *mylog << n << " ";
    sum += n*i;
    i += 1;
  }
  *mylog << std::endl;
  *mylog << "Fraction of singles:" << (float)freq[0]/bcs.size() << " sum:" << sum << " size:"  << bcs.size() << std::endl;
}
//
// Frequency: # of TVX with 1,2,3,... collision BCs
// Assuming sorted in Vollision BC (bcEvSel)
void bcInfos::frequencyBCSorted(std::array<int,Nfreq>& freq)
{
  std::map<uint64_t, int> bcFreq;
  uint64_t prev = 0;
  int count = 0;
  for(auto const& bc: bcs) {
    //*mylog << prev << " " << bc.bcEvSel << std::endl;
    if(bc.bcEvSel != prev) {
      if(prev) {
        bcFreq[prev] = count;
        //*mylog << "saving "  << prev << " " << count << std::endl;
      }
      if( count > 5) {
        //*mylog << "count:" << count << " " << prev << std::endl;
      }
      count = 0;
      prev = bc.bcEvSel;
      //*mylog << "in " << std::endl;
    } else {
      count++;
    }
  }
  bcFreq[prev] = count;
  //*mylog << "saving "  << prev << " " << count << std::endl;
  //printMap(bcFreq);
  //
  for(auto const& bc: bcFreq) {
    if(bc.second < Nfreq) {
      //*mylog << " *** " << bc.first << " " << bc.second << std::endl;
      //*mylog << freq << std::endl;
      freq[bc.second]++;
    }  else  {
      freq[Nfreq-1]++;
    }
  }
  (*mylog) << "Frequency of evSel (TVX), i.e. number of singles, number of doubles. ..." << std::endl;
  float sum = 0;
  int npairs = 0;
  float i = 1;
  for(auto const& n: freq) {
    *mylog << n << " ";
    sum += n*i;
    if(i > 1) {
      npairs += n*TMath::Binomial(i,2);
    }
    i += 1.;
  }
  *mylog << std::endl;
  *mylog << "Fraction of singles:" << (float)freq[0]/bcs.size() << " sum:" << sum << " pairs"  << npairs << " size:"  << bcs.size() << std::endl;
}
//
// Selection efficiency - should be raun on cleaned bcs ?
//
void bcInfos::printSelectionEfficiency(std::array<double_t, Ndim>& eff) const
{
  /*
  for(size_t i = 0; i < Ndim; i++ ) {
    if( i < labels.size() ){
      *mylog << labels[i] << " bit "  << i << " eff:" << eff[i] << std::endl;
    } else {
      *mylog << "XXX bit " << ":"  << i << " eff:" << eff[i] << std::endl;
    }
  }
  */
  TH1F * h = new TH1F("Dowscaling","Downscaling", Ndim_used,0,Ndim_used);
  for(size_t i = 0 ; i < Ndim_used; i ++) {
    h->Fill(i,eff[i]);
  }
}
void bcInfos::selectionEfficiency( std::array<double_t,Ndim>& eff,std::array<double_t,Ndim>& err ) const
{
  *mylog << "===> Selection efficiency - downscaling" << std::endl;
  //std::array<double_t,Ndim> eff{0};
  for(int i =0; i < Ndim; i++) {
    double_t before = mFilterCounters->GetBinContent(i+2);
    double_t after = 0;
    uint64_t mask = (1ull << i);
    for(auto const& bc: bcs) {
      if(bc.selMask & mask) {
        after++;
      }
    }
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
    *mylog << "BIT:" << std::setw(2) << i << "  " << std::setw(25) << label << "   b:" << std::setw(10) << before << " a:" << std::setw(10) << after << " eff:" << eff[i] << std::endl;
  }
  //printSelectionEfficiency(eff);
}
void bcInfos::selectionEfficiencyBiased() const
{
  *mylog << "===> Selection efficiency - downscaling 2" << std::endl;
  std::array<double_t,Ndim> eff{0};
  for(int i =0; i < Ndim; i++) {
    double_t before = 0;
    double_t after = 0;
    uint64_t mask = (1ull << i);
    for(auto const& bc: bcs) {
      if(bc.trigMask & mask) {
        before++;
      }
      if(bc.selMask & mask) {
        after++;
      }
    }
    eff[i] = after/before;
    std::string label = "?";
    if(i < (int) labels.size()) {
      label = labels[i];
    }
    *mylog << "BIT:" << std::setw(2) << i << "  " << std::setw(25) << label << "   b:" << std::setw(10) << before << " a:" << std::setw(10) << after << " eff:" << after/before << std::endl;
  }
  //printSelectionEfficiency(eff);
}
//
// Correlation matrix
//
void bcInfos::printCorMatrix(float_t (*mat)[Ndim][Ndim]) const
{
  *mylog << std::dec << "Matrix " << Ndim << std::endl;
  for(int i = 0; i < Ndim; i++) {
    if(!zeros[i]) {
      for(int j = 0; j < Ndim; j++) {
        if(!zeros[j]) {
          *mylog <<   std::setw(10) << std::setprecision(2) << (*mat)[i][j] << " ";
        }
      }
      *mylog << std::endl;
    }
  }
}
void bcInfos::printCorMatrixSimple(float_t (*mat)[Ndim][Ndim]) const
{
  *mylog << std::dec << "Matrix " << Ndim << std::endl;
  for(int i = 0; i < Ndim; i++) {
    double_t sum = 0;
    if(!zeros[i]) {
      for(int j = 0; j < Ndim; j++) {
        if(!zeros[j]) {
          if(i != j) {
            //*mylog <<   std::setw(10) << std::setprecision(2) << (*mat)[i][j] << " ";
            sum += ((*mat)[i][j]);
          }
        }
      }
    }
    *mylog << "covMatrix:" << i << " sum off diag:" << sum << std::endl;
  }
}
void bcInfos::corMatrix()
{
  int NTotTrigs = totalTriggered;
  float_t cm[Ndim][Ndim];
  //std::array<std::array<int,Ndim>,Ndim> cmm;
  float_t sigma[Ndim], mean[Ndim];
  for(int i = 0; i < Ndim; i++) {
    uint64_t maski = 1ull << i;
    int sum = 0;
    for(auto const& bc: bcs) {
      sum += (bc.trigMask & maski) != 0;
    }
    if (sum) {
      float_t mu = (float_t)sum / (float_t) NTotTrigs;
      mean[i] = mu;
      sigma[i] = sqrt(mu *(1.- mu));
      cm[i][i] = 1.;
      //*mylog << i << " " << i << " MAT:" << cm[i][i] << " mean:" << mu << std::endl;
    } else {
      zeros[i] = 1;
    }
  }
  //*mylog << "zeros: " << zeros << std::endl;
  //
  for(int i = 0; i < Ndim; i++) {
    uint64_t maski = 1ull << i;
    if(!zeros[i]) {
    for(int j = i + 1; j < Ndim; j++) {
      if (!zeros[j]) {
      //if(1) {
      uint64_t maskj = 1ull << j;
      int sumij = 0;
      for(auto const& bc: bcs) {
        int it = (bc.trigMask & maski) != 0;
        int jt = (bc.trigMask & maskj) != 0;
        sumij += it*jt;
      }
      cm[i][j] = ((float_t)(sumij)/NTotTrigs - mean[i]*mean[j])/sigma[i]/sigma[j];
      //*mylog << i << " " << j << " MAT:" << cm[i][j] << " " << sumij << std::endl;
      cm[j][i] = cm[i][j];
    }
  }
  }
  }
  printCorMatrix(&cm);
  printCorMatrixSimple(&cm);
}
//
// Map collisions to ao2d bcs
//
void bcInfos::evSel2AOD()
{
  //std::map<uint64_t,std::vector<bcInfo>> evsel2aod;
  std::array<int,8> counters{0};
  uint64_t prev = 0;
  for(auto const& bc: bcs) {
    //*mylog << " bc.bcEvSel:" << bc.bcEvSel << std::endl;
    //continue;
    if(bc.bcEvSel != prev) {
      evsel2aod[bc.bcEvSel].push_back(bc);
      //*mylog << "saving "  << prev << " " << bc.bcEvSel << std::endl;
      prev = bc.bcEvSel;
    } else {
      //*mylog << "saving else "  << prev << " " << bc.bcEvSel << std::endl;
      //*mylog << bc.bcEvSel << " bcs " << bc.bcAOD << std::endl;
        evsel2aod[prev].push_back(bc);
    }
  }
  //printEvSel2AOD();
  //*mylog << "LUT Counters:" << counters << std::endl;
}
void bcInfos::evSel2AODInverse(std::vector<bcInfo>& bcs)
{
  for(auto const& bc: evsel2aod) {
    for(auto const& bcinfo: bc.second) {
      bcs.push_back(bcinfo);
    }
  }
}
void bcInfos::cleanBC(std::vector<bcInfo>& bcs_cleaned)
{
  // use set instead of vector ?
  std::array<int,8> counters{0};
  int tot = 0;
  if(evsel2aod.size() == 0) {
    *mylog << "Creating evsel2aod" << std::endl;
    evSel2AOD();
  }
  for(auto const& bc: evsel2aod) {
    const int nbc = bc.second.size();
    std::vector<int> removed(nbc,0);
    for(int i = 0; i < nbc; i++) {
      if( !removed[i]) {
        for(int j = i+1; j < nbc; j++) {
          bool bceq = (bc.second[i].bcAOD == bc.second[j].bcAOD);
          bool treq = (bc.second[i].trigMask == bc.second[j].trigMask);
          bool seeq = (bc.second[i].selMask == bc.second[j].selMask);
          uint16_t lut = bceq + 0x2*treq + 0x4*seeq;
          counters[lut] += 1;
          tot++;
          switch(lut)
          {
            //
            case 0:
            // TVX has two different bcAOD with two different triggers and selections (sel) - splitted vtx ?
              break;
            case 1:
            // bcAOD is equal but two different triggers with same sel ?
              break;
            case 2:
            // different bcAOD , same trigger, different sel - split vertex ?
              break;
            case 3:
            // same bcAOD, same trigger , different sel ?
              break;
            case 4:
            // different bcAOD and trigger but same sel - split vertex ?
              break;
            case 5:
            // same bcAOD, different trigger , same sel ?
              break;
            case 6:
            // different bcAOD, same trigger and sel - split vertex ?
              break;
            case 7:
            // exactly same record, why is it here ?
              removed[j] = 1;
              break;
          }
        }
        bcs_cleaned.push_back(bc.second[i]);
      }
    }
  }
  *mylog << "LUT counters:" << counters << " tot:" << tot << std::endl;
  *mylog << "bcs cleaned size:" << bcs_cleaned.size() << std::endl;
}
void bcInfos::cleanBC()
{
  cleanBC(bcs_cleaned);
}
void bcInfos::getArrayForBit(std::vector<uint64_t>& bctrigs, std::vector<uint64_t>& bcselec, int bit)
{
  // create arrays
  evSel2AOD();
  std::vector<bcInfo> bcs_cleaned;
  cleanBC(bcs_cleaned);
  uint64_t mask = 1ull << bit;
  for(auto const&bcinfo: bcs_cleaned) {
    if(bcinfo.trigMask & mask) {
      bctrigs.push_back(bcinfo.bcEvSel);
    }
    if(bcinfo.selMask & mask) {
      bcselec.push_back(bcinfo.bcEvSel);
    }
  }
  *mylog << "BIT:" << bit << " # trigs:" << bctrigs.size() << " # selec:" << bcselec.size() << std::endl;
  //
}
//
//================================
//
struct Hists
{
  Hists() = default;
  TObjArray hists; 
  void addHist(std::string const&name,int n, int n1, int n2)
  {
    TH1F *h = new TH1F(name.c_str(),name.c_str(),n,n1,n2); 
    h->Sumw2();
    hists.Add(h);
  }
  void addHist(TH1* h)
  {
    hists.Add(h);
  }
  void writeHists()
  {
    TFile *file;
    file = new TFile("Histos.root","RECREATE","FILE");
    hists.Write(); 
  }
};

//
//===========================
//
struct effUtils : Hists
{
  effUtils() = default;
  bcInfos originalBCs, skimmedBCs;
  std::array<double_t,Ndim> effSkimmed{0};
  std::array<double_t,Ndim> effSkimmedC{0};
  std::array<double_t,Ndim> effSkimmedError{0};
  std::array<double_t,Ndim> downscaleFactors{0};
  std::array<double_t,Ndim> downscaleFactorsError{0};
  //
  void extractLabelsAnal(std::vector<std::string>& labels);
  void extractLabels(std::vector<std::string>& labels);
  void readFiles(TFile& originalFile, TFile& skimmedFile);
  void readFiles(TFile& file);
  void skimmedEfficiency();
  void skimmedEfficiency_cleaned();
  void getDownscaleFactors() {originalBCs.selectionEfficiency(downscaleFactors,downscaleFactorsError);};
  void correlatev0(int i, int j);
  void correlate(int i, int j);
  void correlateAll(int delta, int ds);
  void printCorrelations(int bit,std::vector<int>& cor, int dist, int Nprint = 10);
  //
  void printDownscaleFactors();
  void fillFreq(TH1* histo, std::array<int,Nfreq>& arr);
  void fillSkimmedEff(TH1 *histo);
  int runNumber = runNUMBER;
};
void effUtils::extractLabelsAnal(std::vector<std::string>& labels)
{
  //TFile file("AnalysisResults.root");
  TFile file("AnalysisResults_fullrun.root");
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
  labels.clear();
  *mylog << "Labels Nbins:" <<  hist1->GetNbinsX() << std::endl;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
    }
  }
  *mylog << "Label size:" << labels.size() << std::endl;
}
void effUtils::extractLabels(std::vector<std::string>& labels)
{
  labels.clear();
  TH1* hist1 = originalBCs.mFilterCounters;
  *mylog << "Labels Nbins:" <<  hist1->GetNbinsX() << std::endl;
  for (int i = 1; i <= hist1->GetNbinsX(); i++) {
    std::string label = hist1->GetXaxis()->GetBinLabel(i);
    if (label != "Total number of events" && label != "Filtered events") {
      labels.push_back(label);
    }
  }
  *mylog << "Label size:" << labels.size() << std::endl;
}
void effUtils::readFiles(TFile& originalFile, TFile& skimmedFile)
{
  *mylog << "=== Unskimmed:" << std::endl;
  originalBCs.getData(originalFile);
  //originalBCs.mCCDBPathTrigScalers = mCCDBPathTrigScalersOrigina;
  *mylog << "=== Skimmed:" << std::endl;
  skimmedBCs.getData(skimmedFile);
  //skimmedBCs.mCCDBPathTrigScalers = mCCDBPathTrigScalersSkimmed;
}
void effUtils::readFiles(TFile& file)
{
  skimmedBCs.getData(file);
  originalBCs.getDataCCDB(runNumber);
}
void effUtils::skimmedEfficiency()
{
  *mylog << "===> Skimming efficiency" << std::endl;
  std::array<double_t,Ndim> before{0},after{0};
  for(auto const& bc: originalBCs.bcs) {
    for(int i =0; i < Ndim; i++) {
      uint64_t mask = (1ull << i);
      if(bc.selMask & mask) {
        before[i]++;
      }
    }
  }
  for(auto const& bc: skimmedBCs.bcs) {
    for(int i =0; i < Ndim; i++) {
      uint64_t mask = (1ull << i);
      if(bc.trigMask & mask) {
        after[i]++;
      }
    }
  }
  for(int i = 0; i < Ndim; i++) {
    //after[i] = skimmedBCs.mFilterCounters->GetBinContent(i+2);
    double_t e = 0;
    if(before[i] != 0) {
      double_t e = after[i]/before[i];
      effSkimmed[i] = e;
      effSkimmedError[i] = sqrt(e*(1-e)/before[i]);
    }
    std::string label = "?";
    if(i < (int) labels.size()) {
      label = labels[i];
    }
    //*mylog << "BIT:" << i << " b:" << before[i] << " a:" << after[i] << " eff:" << effSkimmed[i] << "+-" << effSkimmedError[i] << std::endl;
    *mylog << "BIT:" << std::setw(2) << i << "  " << std::setw(25) << label  << "   b:" << std::setw(10) << before[i] << " a:" << std::setw(10) << after[i] << " eff:" << effSkimmed[i] << std::endl;

  }
}
void effUtils::skimmedEfficiency_cleaned()
{
  *mylog << "===> Skimming efficiency of cleaned" << std::endl;
  for(int i =0; i < Ndim; i++) {
    double_t before = 0;
    double_t after = 0;
    uint64_t mask = (1ull << i);
    for(auto const& bc: originalBCs.bcs_cleaned) {
      if(bc.selMask & mask) {
        before++;
      }
    }
    for(auto const& bc: skimmedBCs.bcs_cleaned) {
      if(bc.trigMask & mask) {
        after++;
      }
    }
    effSkimmedC[i] = after/before;
    *mylog << "BIT:" << i << " b:" << before << " a:" << after << " eff:" << effSkimmedC[i] << std::endl;
  }
}
void effUtils::correlatev0(int ibit, int jbit)
{
  int delta = 3;
  const int ndimcor = 2*delta + 1;
  std::vector<int> corr(ndimcor,0);
  uint64_t maski = 1ull << ibit;
  uint64_t maskj = 1ull << jbit;
  std::vector<bcInfo>& v1 = originalBCs.bcs;
  std::vector<bcInfo>& v2 = originalBCs.bcs;
  int ndatai = v1.size();
  int ndataj = v2.size();

  *mylog << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
  //
  for(int i = 0; i < ndatai; i++) {
    int posi = v1[i].bcEvSel;
    bool trigi = v1[i].selMask & maski;
    if(trigi) {
      for(int j = 0; j < ndataj; j++) {
        int posj = v2[j].bcEvSel;
        if( posj+delta < posi) {
          continue;
        }
        if(posi+delta < posj) {
          continue;
        }
        int dist = posi - posj;
        //*mylog << dist << std::endl;
        bool trigj = v2[j].trigMask & maskj;
        if(trigi && trigj) {
          corr[dist + delta] += 1;
        }
      }
    }
  }
  *mylog << "Correlations:" << corr << std::endl;
}
//
// array instead of vectors - it is 10/13.8 faster
//
void effUtils::correlate(int ibit, int jbit)
{
  int delta = 10;
  const int ndimcor = 2*delta + 1;
  std::vector<int> corr(ndimcor,0);
  uint64_t maski = 1ull << ibit;
  uint64_t maskj = 1ull << jbit;
  //std::vector<bcInfo>& v1 = originalBCs.bcs_cleaned;
  //std::vector<bcInfo>& v2 = skimmedBCs.bcs_cleaned;
  std::vector<bcInfo>& v1 = originalBCs.bcs;
  std::vector<bcInfo>& v2 = skimmedBCs.bcs;
  int ndatai = v1.size();
  int ndataj = v2.size();
  *mylog << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
  //
  uint64_t *arrbc1 = new uint64_t[ndatai];
  uint64_t *arrtr1 = new uint64_t[ndatai];
  int i = 0;
  for(auto const& bcinfo: v1) {
    arrbc1[i] = bcinfo.bcEvSel;
    arrtr1[i] = bcinfo.selMask;
    i++;
  }
  uint64_t *arrbc2 = new uint64_t[ndataj];
  uint64_t *arrtr2 = new uint64_t[ndataj];
  int j = 0;
  for(auto const& bcinfo: v2) {
    arrbc2[j] = bcinfo.bcEvSel;
    arrtr2[j] = bcinfo.trigMask;
    j++;
  }
  //
  for(int i = 0; i < ndatai; i++) {
    int posi = arrbc1[i];
    bool trigi = arrtr1[i] & maski;
    if(trigi) {
      for(int j = 0; j < ndataj; j++) {
        int posj = arrbc2[j];
        if( posj+delta < posi) {
          continue;
        }
        if(posi+delta < posj) {
          continue;
        }
        int dist = posi - posj;
        //*mylog << dist << std::endl;
        bool trigj = arrtr2[j] & maskj;
        if(trigi && trigj) {
          //*mylog << dist+delta << std::endl;
          corr[dist + delta] += 1;
        }
      }
    }
  }
  *mylog << "Correlations:" << corr << std::endl;
}
void effUtils::correlateAll(int delta, int ds)
{
  int Nprint = 10;
  *mylog << "Starting correlation all. Correlation window +/-" << delta << std::endl;
  //*mylog << "Printing only +/-" << Nprint << std::endl;
  const int ndimcor = 2*delta + 1;
  std::vector<bcInfo>& v1 = originalBCs.bcs;
  std::vector<bcInfo>& v2 = skimmedBCs.bcs;
  int ndatai = v1.size();
  int ndataj = v2.size();
  *mylog << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
  //
  uint64_t *arrbc1 = new uint64_t[ndatai];
  uint64_t *arrtr1 = new uint64_t[ndatai];
  int idata = 0, i = 0;
  for(auto const& bcinfo: v1) {
    if( (i%ds) == 0 ) {
      arrbc1[idata] = bcinfo.bcEvSel;
      arrtr1[idata] = bcinfo.selMask;
      idata++;
    }
    i++;
  }
  uint64_t *arrbc2 = new uint64_t[ndataj];
  uint64_t *arrtr2 = new uint64_t[ndataj];
  int jdata = 0, j = 0;
  for(auto const& bcinfo: v2) {
    if( (j%ds) == 0 && ds) {	  
      arrbc2[jdata] = bcinfo.bcEvSel;
      arrtr2[jdata] = bcinfo.trigMask;
      jdata++;
    }
    j++;
  }
  *mylog << "idata:" << idata << " jdata:" << jdata << std::endl;
  //
  for(int l = 0; l < Ndim; l++) {
    std::vector<int> corr(ndimcor,0);
    uint64_t maskl = 1ull << l;
    for(int i = 0; i < idata; i++) {
      int posi = arrbc1[i];
      bool trigi = arrtr1[i] & maskl;
      if(trigi) {
        for(int j = 0; j < jdata; j++) {
          int posj = arrbc2[j];
          if( posj+delta < posi) {
            continue;
          }
          if(posi+delta < posj) {
            continue;
          }
          int dist = posi - posj;
          //*mylog << dist << std::endl;
          bool trigj = arrtr2[j] & maskl;
          if(trigi && trigj) {
            //*mylog << dist+delta << std::endl;
            corr[dist + delta] += 1;
          }
        }
      }
    }
    *mylog << "BIT:" << l ; //<< " Correlations:" << corr << std::endl;
    printCorrelations(l,corr,delta,delta);
  }
}
void effUtils::printCorrelations(int bit, std::vector<int>& corr, int dist, int Nprint)
{
  if( (int)corr.size() < dist+Nprint) {
    *mylog << "printCorrelations par not compatible" << std::endl;
    return;
  }
  for(int i = 0; i < 2*Nprint + 1; i++) {
    *mylog << " " << corr[dist - Nprint + i];
  }
  *mylog << std::endl;
  //
  std::string name = "["+ std::to_string(bit) +"] ";
  if( bit < (int)labels.size() ) {
    name += labels[bit];
  } else {
  }
  TH1F* h = new TH1F(name.c_str(), name.c_str(), 2*dist+1, -dist, dist);
  for(int i = 0; i < 2*dist + 1; i++) {
    h->Fill(i-dist,corr[i]);
  }
  addHist(h);
}
void effUtils::fillFreq(TH1* h, std::array<int,Nfreq>& arr)
{
  int Nfreq = 10;
  float_t sum = 0;
  for(int i = 0;i < Nfreq; i ++) {
    h->Fill(i,arr[i]);
    //*mylog << "fU:" << freqU[i] << std::endl;;
    sum += (i+1)*arr[i];
  }
  *mylog << "h->INtegral:" << h->Integral() << " wsum:" << sum << std::endl;
  h->Scale(1./sum);
  addHist(h);
}
void effUtils::fillSkimmedEff(TH1* h)
{
  h->SetAxisRange(0,1.2,"Y");
  for(int i = 0; i <Ndim_used; i++)
  {
    //h->Fill(i,effSkimmed[i]);
    h->SetBinContent(i+1,effSkimmed[i]);
    h->SetBinError(i+1,effSkimmedError[i]);
  }
  addHist(h);
}
void effUtils::printDownscaleFactors()
{
  getDownscaleFactors();
  TH1F * h = new TH1F("Downscale Factors","Downscale Factors", Ndim_used,0,Ndim_used);
  h->SetAxisRange(0,1.1,"Y");
  for(int i = 0; i < Ndim_used; i++) {
    //h->Fill(i,downscaleFactors[i]);
    h->SetBinContent(i+1,downscaleFactors[i]);
    h->SetBinError(i+1,downscaleFactorsError[i]);
  }
  addHist(h);
}
//
// main
//
void eff3H(std::string original = "bcRanges.root", std::string skimmed = "bcRanges_fullrun-skimmed.root")
{
  ofstream ff;
  ff.open("file.txt");
  mylog = &ff;
  int deltaCor = 1000;
  TH1F * hfreqU = new TH1F("Freq UnSkimmed","Freq UnSkimmed", Nfreq,0,Nfreq);
  TH1F * hfreqS = new TH1F("Freq Skimmed","Freq Skimmed", Nfreq,0,Nfreq);
  TH1F * hEffSkimmed = new TH1F("Skimming eff","Skimming eff", Ndim_used,0,Ndim_used);
  //
  TFile originalFile(original.data());
  TFile skimmedFile(skimmed.data());
  effUtils eff;
  //eff.addHist("Freq Unskimm",Nfreq,0,Nfreq);
  //
  //eff.readFiles(originalFile,skimmedFile);
  eff.readFiles(skimmedFile);
  eff.extractLabels(labels);
  //eff.originalBCs.dataSize(1);
  //return;
  //return;
  *mylog << "=== Unskimmed" << std::endl;
  //eff.originalBCs.frequencyBC();
  eff.originalBCs.evSel2AOD();
  //eff.originalBCs.cleanBC();
  std::array<int,Nfreq> freqU{0};
  eff.originalBCs.frequencyBCSorted(freqU);
  eff.fillFreq(hfreqU, freqU);
  eff.originalBCs.frequencyCol();
  //eff.originalBCs.printbcInfoVect(bcs_cleaned);
  //eff.originalBCs.selectionEfficiency();
  eff.printDownscaleFactors();
  //return ;
  //
  *mylog << "=== Skimmed" << std::endl;
  eff.skimmedBCs.evSel2AOD();
  //eff.skimmedBCs.cleanBC();
  std::array<int,Nfreq> freqS{0};
  eff.skimmedBCs.frequencyBCSorted(freqS);
  eff.fillFreq(hfreqS,freqS);
  eff.skimmedEfficiency();
  eff.fillSkimmedEff(hEffSkimmed);
  //
  // correlations
  //
  //eff.originalBCs.corMatrix();
  //
  clock_t start = clock();
  eff.correlateAll(deltaCor,1);
  clock_t stop = clock();
  double_t time = double_t(stop - start)/ (double_t)CLOCKS_PER_SEC;
  *mylog << "Time:" << time << std::endl;
  eff.writeHists();
}
