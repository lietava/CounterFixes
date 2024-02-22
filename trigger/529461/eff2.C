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

#include <TFile.h>
#include <TTree.h>
std::vector<std::string> labels = {"fPHOSnbar", "fPHOSPair", "fPHOSElectron", "fPHOSPhoton", "fOmegaLargeRadius", "fSingleXiYN", "fQuadrupleXi", "fhadronXi", "fTripleXi", "fGammaHighPtDCAL", "fGammaHighPtEMCAL", "fJetFullHighPt", "fLD", "fPD", "fLLL", "fPLL", "fPPL", "fPPP", "fHighFt0cFv0Flat", "fHighFt0cFv0Mult", "fHighFt0Flat", "fHighFt0Mult", "fHighTrackMult", "fHfDoubleCharmMix", "fHfDoubleCharm3P", "fHfSoftGamma3P", "fHfFemto2P", "fHfBeauty4P", "fHfFemto3P", "fHfBeauty3P", "fHfSoftGamma2P", "fHfDoubleCharm2P", "fDiMuon", "fDiElectron", "fUDdiff", "fHe"};
//
// main structure
//
struct bcInfo {
  ULong64_t bcAOD, bcEvSel, trigMask, selMask;
  void print() const {
    std::cout << bcEvSel << " " << bcAOD << " " << std::hex << trigMask << " " << selMask << " " << std::dec;
    for(int i = 0; i < 64; i++) {
      ULong64_t mask = 1ull << i;
      if(mask & trigMask) {
        if(i < labels.size()) {
          std::cout << labels[i] << " ";
        } else {
          std::cout << "XXX ";
        }
      }
    }
    std::cout << std::endl;
  }
};
//
// utils
//
void printEvSel2AOD(std::map<uint64_t,std::vector<bcInfo>>& m)
{
  for(auto const& i: m) {
    std::cout << "===>evSel BC:" << i.second.size() << " " << i.first << ":" << std::endl;
    for(auto const& j: i.second){
      j.print();
    }
  }
};
void printbcInfoVect(std::vector<bcInfo>& v)
{
  for(auto const& i: v) {
    i.print();
  }
};
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
    std::cout << i.first << " " << i.second << std::endl;
  }
};
//
// read data from root file
//
void getData(TFile& inputFile, std::vector<bcInfo>& bcs, std::array<int, 64>& selectionCounters, std::array<int, 64>& triggerCounters, int Nmax = 0)
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
        for (ULong64_t j = 0; j < 64; j++) {
          if (bcAO2D.selMask & (1ull << j))
            selectionCounters[j]++;
          if (bcAO2D.trigMask & (1ull << j))
            triggerCounters[j]++;
        }
      }
    }
  }
}
//
// Frequency: # of TVX with 1,2,3,... collision BCs
//
void frequencyBC(std::vector<bcInfo>& bcs)
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
  std::array<int,16> freq{0};
  for(auto const& bc: bcFreq) {
    if(bc.second < 15) {
      freq[bc.second]++;
    }  else  {
      freq[15]++;
    }
  }
  float sum = 0;
  float i = 1;
  for(auto const& n: freq) {
    std::cout << n << " ";
    sum += n*i;
    i += 1;
  }
  std::cout << std::endl;
  std::cout << "Fraction of singles:" << (float)freq[0]/bcs.size() << " sum:" << sum << " size:"  << bcs.size() << std::endl;
}
//
// Frequency: # of TVX with 1,2,3,... collision BCs
// Assuming sorted in Vollision BC (bcEvSel)
void frequencyBCSorted(std::vector<bcInfo>& bcs)
{
  std::map<uint64_t, int> bcFreq;
  uint64_t prev = 0;
  int count = 0;
  for(auto const& bc: bcs) {
    //std::cout << prev << " " << bc.bcEvSel << std::endl;
    if(bc.bcEvSel != prev) {
      if(prev) {
        bcFreq[prev] = count;
        //std::cout << "saving "  << prev << " " << count << std::endl;
      }
      count = 0;
      prev = bc.bcEvSel;
      //std::cout << "in " << std::endl;
    } else {
      count++;
    }
  }
  bcFreq[prev] = count;
  //std::cout << "saving "  << prev << " " << count << std::endl;
  //printMap(bcFreq);
  //
  std::array<int,16> freq{0};
  for(auto const& bc: bcFreq) {
    if(bc.second < 15) {
      //std::cout << " *** " << bc.first << " " << bc.second << std::endl;
      //std::cout << freq << std::endl;
      freq[bc.second]++;
    }  else  {
      freq[15]++;
    }
  }
  float sum = 0;
  float i = 1;
  for(auto const& n: freq) {
    std::cout << n << " ";
    sum += n*i;
    i += 1;
  }
  std::cout << std::endl;
  std::cout << "Fraction of singles:" << (float)freq[0]/bcs.size() << " sum:" << sum << " size:"  << bcs.size() << std::endl;
}
//
const int Ndim = 64;
//
// Selection efficiency - should be raun on cleaned bcs ?
//
void printSelectionEfficiency(std::array<double_t, Ndim>& eff)
{
  for(int i = 0; i < Ndim; i++ ) {
    if( i < labels.size() ){
      std::cout << labels[i] << " bit "  << i << " eff:" << eff[i] << std::endl;
    } else {
      std::cout << "XXX bit " << ":"  << i << " eff:" << eff[i] << std::endl;
    }
  }
}
void selectionEfficiency(std::vector<bcInfo>& bcs)
{
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
  }
  printSelectionEfficiency(eff);
}
//
// Correlation matrix
//
std::array<bool,Ndim> zeros{0};  // Take care of empty or not used bits
void printCorMatrix(float_t (*mat)[Ndim][Ndim])
{
  std::cout << std::dec << "Matrix " << Ndim << std::endl;
  for(int i = 0; i < Ndim; i++) {
    if(!zeros[i]) {
      for(int j = 0; j < Ndim; j++) {
        if(!zeros[j]) {
          std::cout <<   std::setw(10) << std::setprecision(2) << (*mat)[i][j] << " ";
        }
      }
      std::cout << std::endl;
    }
  }
}
void corMatrix(std::vector<bcInfo>& bcs, int NTotTrigs)
{
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
      std::cout << i << " " << i << " MAT:" << cm[i][i] << " mean:" << mu << std::endl;
    } else {
      zeros[i] = 1;
    }
  }
  std::cout << "zeros: " << zeros << std::endl;
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
      std::cout << i << " " << j << " MAT:" << cm[i][j] << " " << sumij << std::endl;
      cm[j][i] = cm[i][j];
    }
  }
  }
  }
  printCorMatrix(&cm);
}
//
// Map collisions to ao2d bcs
//
void evSel2AOD(std::vector<bcInfo>& bcs, int skim = 0)
{
  std::cout << "Skimming bcs:" << skim << std::endl;
  std::map<uint64_t,std::vector<bcInfo>> evsel2aod;
  std::array<int,8> counters{0};
  uint64_t prev = 0;
  for(auto const& bc: bcs) {
    std::cout << " bc.bcEvSel:" << bc.bcEvSel << std::endl;
    //continue;
    if(bc.bcEvSel != prev) {
      evsel2aod[bc.bcEvSel].push_back(bc);
      //std::cout << "saving "  << prev << " " << bc.bcEvSel << std::endl;
      prev = bc.bcEvSel;
    } else {
      //std::cout << "saving else "  << prev << " " << bc.bcEvSel << std::endl;
      //std::cout << bc.bcEvSel << " bcs " << bc.bcAOD << std::endl;
      if(skim) {
        for(auto const& bcmap: evsel2aod[bc.bcEvSel]) {
          bool bceq = (bcmap.bcAOD == bc.bcAOD);
          bool treq = (bcmap.trigMask == bc.trigMask);
          bool seeq = (bcmap.selMask == bc.selMask);
          //uint16_t lupt = bceq + 0x2*treq + 0x4*seeq;
          int lupt = bceq + 0x2*treq + 0x4*seeq;
	  counters[lupt] += 1;
          switch (lupt)
          {
            case 0x7: // all equal - remove, i.e. do not save
              break;
            default:
              std::cout  << " case:" << lupt << " " << bc.bcAOD << std::endl;
              evsel2aod[prev].push_back(bc);
          }
        }
      } else {
        evsel2aod[prev].push_back(bc);
      }
    }
  }
  //printEvSel2AOD(evsel2aod);
  std::cout << "Counters:" << counters << std::endl;
}
//
// main
//
void eff2(std::string original = "bcRanges_fullrun.root", std::string skimmed = "bcRanges_fullrun-skimmed.root")
{
  TFile originalFile(original.data());
  TFile skimmedFile(skimmed.data());
  std::vector<bcInfo> originalBCs, skimmedBCs;
  std::array<int, 64> originalSelected{0}, originalTriggered{0}, skimmedSelected{0}, skimmedTriggered{0};
  int originalTotal = 0, skimmedTotal = 0;
  getData(originalFile, originalBCs, originalSelected, originalTriggered,00);
  getData(skimmedFile, skimmedBCs, skimmedSelected, skimmedTriggered,00);

  std::cout << "Original BCs: " << originalBCs.size() << std::endl;
  for (int i = 0; i < 64; i++) {
    // std::cout << "* BIT " << i << " (" << labels[i] << ") selected: " << originalSelected[i] << " triggered: " << originalTriggered[i] << std::endl;
    originalTotal += originalSelected[i];
  }
  std::cout << "* Total original triggers: " << originalTotal << std::endl;

  std::cout << "Skimmed BCs: " << skimmedBCs.size() << std::endl;
  for (int i = 0; i < 64; i++) {
    // std::cout << "* BIT " << i << " (" << labels[i] << ") selected: " << originalSelected[i] << " triggered: " << originalTriggered[i] << std::endl;
    skimmedTotal += skimmedSelected[i];
  }
  std::cout << "* Total skimmed triggers: " << skimmedTotal << std::endl;
  //
  std::sort(originalBCs.begin(), originalBCs.end(), [](const bcInfo& a, const bcInfo& b) { return a.bcEvSel < b.bcEvSel; });
  std::sort(skimmedBCs.begin(), skimmedBCs.end(), [](const bcInfo& a, const bcInfo& b) { return a.bcEvSel < b.bcEvSel; });
  //printbcInfoVect(originalBCs);
  //frequencyBC(originalBCs);
  frequencyBCSorted(originalBCs);
  evSel2AOD(originalBCs,1);
  //selectionEfficiency(originalBCs);
  //evsel2AOD(skimmedBCs);
  //corMatrix(originalBCs,originalTotal);
  // for (int i = 0; i < originalBCs.size() - 1; i++) {
  //   if (originalBCs[i].bcEvSel == originalBCs[i + 1].bcEvSel) {
  //     std::cout << "Duplicate BC in original: " << originalBCs[i].bcEvSel << "\t" << originalBCs[i].bcAOD << " - " << originalBCs[i + 1].bcAOD << std::endl;
  //   }
  // }

  // TFile outputFile("diffSkimmingBCs.root", "recreate");

}
