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
#include <iomanip>
#include <TFile.h>
#include <TTree.h>
std::vector<std::string> labels = {"fPHOSnbar", "fPHOSPair", "fPHOSElectron", "fPHOSPhoton", "fOmegaLargeRadius", "fSingleXiYN", "fQuadrupleXi", "fhadronXi", "fTripleXi", "fGammaHighPtDCAL", "fGammaHighPtEMCAL", "fJetFullHighPt", "fLD", "fPD", "fLLL", "fPLL", "fPPL", "fPPP", "fHighFt0cFv0Flat", "fHighFt0cFv0Mult", "fHighFt0Flat", "fHighFt0Mult", "fHighTrackMult", "fHfDoubleCharmMix", "fHfDoubleCharm3P", "fHfSoftGamma3P", "fHfFemto2P", "fHfBeauty4P", "fHfFemto3P", "fHfBeauty3P", "fHfSoftGamma2P", "fHfDoubleCharm2P", "fDiMuon", "fDiElectron", "fUDdiff", "fHe"};
//
const int Ndim = 64;
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
    std::cout << i.first << " " << i.second << std::endl;
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
  std::cout << bcEvSel << " " << bcAOD << " " << std::hex << trigMask << " " << selMask << " " << std::dec;
  for(size_t i = 0; i < Ndim; i++) {
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
};
struct bcInfos
{
  bcInfos() = default;
  std::vector<bcInfo> bcs;
  std::vector<bcInfo> bcs_cleaned;
  std::array<int, Ndim> selectionCounters{0}, triggerCounters{0};
  int totalSelected = 0;
  int totalTriggered = 0;
  std::array<bool,Ndim> zeros{0};  // Take care of empty or not used bits
  std::map<uint64_t,std::vector<bcInfo>> evsel2aod;
  //
  void getData(TFile& inputFile, int Nmax = 0);
  void frequencyBC() const;
  void frequencyBCSorted();
  void frequencyBCSorted(std::vector<bcInfo>& bcsa);
  void evSel2AOD();
  //void AOD2evSel();
  void evSel2AODInverse(std::vector<bcInfo>& bcs);
  void cleanBC(std::vector<bcInfo>& bcs_cleaned);
  void cleanBC();
  void selectionEfficiency() const;
  void corMatrix();
  void getArrayForBit(std::vector<uint64_t>& bctrigs, std::vector<uint64_t>& bcselec, int bit);
  //
  void printbcInfoVect();
  void printbcInfoVect(std::vector<bcInfo>& bcs);
  void printEvSel2AOD() const;
  void printSelectionEfficiency(std::array<double_t, Ndim>& eff) const;
  void printCorMatrix(float_t (*mat)[Ndim][Ndim]) const;
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
    cefpTree->SetBranchAddress("selMask", &bcAO2D.selMask);
    cefpTree->SetBranchAddress("triMask", &bcAO2D.trigMask);
    for (int i = 0; i < cefpTree->GetEntries(); i++) {
      if((i < Nmax) || (Nmax == 0)) {
        cefpTree->GetEntry(i);
        bcs.push_back(bcAO2D);
        // Check consistrmcy
        if(~bcAO2D.trigMask & bcAO2D.selMask) {
          std::cout << "ERROR selMask is not subset of trigMask:";
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
  std::cout << "bcs Size:" << bcs.size() << std::endl;
  for (int i = 0; i < Ndim; i++) {
    totalSelected += selectionCounters[i];
  }
  for (int i = 0; i < Ndim; i++) {
    totalTriggered += triggerCounters[i];
  }
  std::cout << "Total original triggers:" << totalTriggered << " selected triggers:" << totalSelected << std::endl;
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
    std::cout << "===>evSel BC:" << i.second.size() << " " << i.first << ":" << std::endl;
    for(auto const& j: i.second){
      j.print();
    }
  }
};

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
void bcInfos::frequencyBCSorted(std::vector<bcInfo>& bcsa)
{
  std::map<uint64_t, int> bcFreq;
  uint64_t prev = 0;
  int count = 0;
  for(auto const& bc: bcsa) {
    //std::cout << prev << " " << bc.bcEvSel << std::endl;
    if(bc.bcEvSel != prev) {
      if(prev) {
        bcFreq[prev] = count;
        //std::cout << "saving "  << prev << " " << count << std::endl;
      }
      if( count > 5) {
        //std::cout << "count:" << count << " " << prev << std::endl;
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
  std::cout << "Frequncy of evSel (TVX), i.e. number of singles, number of doubles. ..." << std::endl;
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
void bcInfos::frequencyBCSorted()
{
  frequencyBCSorted(bcs);
}
//
// Selection efficiency - should be raun on cleaned bcs ?
//
void bcInfos::printSelectionEfficiency(std::array<double_t, Ndim>& eff) const
{
  for(size_t i = 0; i < Ndim; i++ ) {
    if( i < labels.size() ){
      std::cout << labels[i] << " bit "  << i << " eff:" << eff[i] << std::endl;
    } else {
      std::cout << "XXX bit " << ":"  << i << " eff:" << eff[i] << std::endl;
    }
  }
}
void bcInfos::selectionEfficiency() const
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
    std::cout << "BIT:" << i << " b:" << before << " a:" << after << " eff:" << eff[i] << std::endl;
  }
  //printSelectionEfficiency(eff);
}
//
// Correlation matrix
//
void bcInfos::printCorMatrix(float_t (*mat)[Ndim][Ndim]) const
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
void bcInfos::evSel2AOD()
{
  //std::map<uint64_t,std::vector<bcInfo>> evsel2aod;
  std::array<int,8> counters{0};
  uint64_t prev = 0;
  for(auto const& bc: bcs) {
    //std::cout << " bc.bcEvSel:" << bc.bcEvSel << std::endl;
    //continue;
    if(bc.bcEvSel != prev) {
      evsel2aod[bc.bcEvSel].push_back(bc);
      //std::cout << "saving "  << prev << " " << bc.bcEvSel << std::endl;
      prev = bc.bcEvSel;
    } else {
      //std::cout << "saving else "  << prev << " " << bc.bcEvSel << std::endl;
      //std::cout << bc.bcEvSel << " bcs " << bc.bcAOD << std::endl;
        evsel2aod[prev].push_back(bc);
    }
  }
  //printEvSel2AOD();
  //std::cout << "LUT Counters:" << counters << std::endl;
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
    std::cout << "Creating evsel2aod" << std::endl;
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
            // TVX has two different bcAOD with two different triggers and selections (sel)
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
  std::cout << "LUT counters:" << counters << " tot:" << tot << std::endl;
  std::cout << "bcs cleaned size:" << bcs_cleaned.size() << std::endl;
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
  std::cout << "BIT:" << bit << " # trigs:" << bctrigs.size() << " # selec:" << bcselec.size() << std::endl;
  //
}
//
//===========================
//
struct effUtils
{
  effUtils() = default;
  bcInfos originalBCs, skimmedBCs;
  void readFiles(TFile& originalFile, TFile& skimmedFile);
  void correlateFFT();
  void correlatev0(int i, int j);
  void correlate(int i, int j);
  void correlateAll(int delta);
  void printCorrelations(std::vector<int>& cor, int dist, int Nprint = 10);
};
void effUtils::readFiles(TFile& originalFile, TFile& skimmedFile)
{
  std::cout << "=== Unskimmed:" << std::endl;
  originalBCs.getData(originalFile, 0);
  std::cout << "=== Skimmed:" << std::endl;
  skimmedBCs.getData(skimmedFile,0);
}
void effUtils::correlateFFT()
{
  std::vector<uint64_t> bctrigsunsk;
  std::vector<uint64_t> bctrigsskim;
  std::vector<uint64_t> bcselecunsk;
  std::vector<uint64_t> bcselecskim;
  originalBCs.getArrayForBit(bctrigsunsk, bcselecunsk, 3);
  skimmedBCs.getArrayForBit(bctrigsskim, bcselecskim, 3);
  std::cout << "first unsk " << bcselecunsk[0] << std::endl;
  std::cout << "first skim " << bctrigsskim[0] << std::endl;
  int nunsk = bcselecunsk.size();
  int nskim = bctrigsskim.size();
  uint64_t first = bcselecunsk[0];
  if( first > bctrigsskim[0]) {
    first = bctrigsskim[0];
  }
  uint64_t last = bcselecunsk[nunsk-1];
  std::cout << " time span:" << last - first << std::endl;
  //double_t *unsk = &bcselecunsk[0];
  //double_t *skim = &bctrigsskim[0];
  double_t *unsk = new double_t(bcselecunsk.size());
  int i = 0;
  for(auto const& bc: bcselecunsk) {
    unsk[i] = bcselecunsk[i] - first;
  }
  double_t *skim = new double_t(bctrigsskim.size());
  i = 0;
  for(auto const& bc: bctrigsskim) {
    skim[i] = bctrigsskim[i] - first;
  }
  // not finished - need version which can work with zero supressed data
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

  std::cout << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
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
        //std::cout << dist << std::endl;
        bool trigj = v2[j].trigMask & maskj;
        if(trigi && trigj) {
          corr[dist + delta] += 1;
        }
      }
    }
  }
  std::cout << "Correlations:" << corr << std::endl;
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
  std::vector<bcInfo>& v1 = originalBCs.bcs_cleaned;
  std::vector<bcInfo>& v2 = skimmedBCs.bcs_cleaned;
  int ndatai = v1.size();
  int ndataj = v2.size();
  std::cout << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
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
        //std::cout << dist << std::endl;
        bool trigj = arrtr2[j] & maskj;
        if(trigi && trigj) {
          //std::cout << dist+delta << std::endl;
          corr[dist + delta] += 1;
        }
      }
    }
  }
  std::cout << "Correlations:" << corr << std::endl;
}
void effUtils::correlateAll(int delta)
{
  int Nprint = 10;
  std::cout << "Starting correlation all. Correlation window +/-" << delta << std::endl;
  //std::cout << "Printing only +/-" << Nprint << std::endl;
  const int ndimcor = 2*delta + 1;
  std::vector<bcInfo>& v1 = originalBCs.bcs_cleaned;
  std::vector<bcInfo>& v2 = skimmedBCs.bcs_cleaned;
  int ndatai = v1.size();
  int ndataj = v2.size();
  std::cout << "ndatai:" << ndatai << " ndataj:" << ndataj << std::endl;
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
  for(int l = 0; l < Ndim; l++) {
    std::vector<int> corr(ndimcor,0);
    uint64_t maskl = 1ull << l;
    for(int i = 0; i < ndatai; i++) {
      int posi = arrbc1[i];
      bool trigi = arrtr1[i] & maskl;
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
          //std::cout << dist << std::endl;
          bool trigj = arrtr2[j] & maskl;
          if(trigi && trigj) {
            //std::cout << dist+delta << std::endl;
            corr[dist + delta] += 1;
          }
        }
      }
    }
    std::cout << "BIT:" << l ; //<< " Correlations:" << corr << std::endl;
    printCorrelations(corr,delta);
  }
}
void effUtils::printCorrelations( std::vector<int>& corr, int dist, int Nprint)
{
  if( (int)corr.size() < dist+Nprint) {
    std::cout << "printCorrelations par not compatible" << std::endl;
    return;
  }
  for(int i = 0; i < 2*Nprint + 1; i++) {
    std::cout << " " << corr[dist - Nprint + i];
  }
  std::cout << std::endl;
}
//
// main
//
void eff3(std::string original = "bcRanges_fullrun.root", std::string skimmed = "bcRanges_fullrun-skimmed.root")
{
  TFile originalFile(original.data());
  TFile skimmedFile(skimmed.data());
  effUtils eff;
  eff.readFiles(originalFile,skimmedFile);
  //eff.originalBCs.printbcInfoVect();
  eff.originalBCs.cleanBC();
  eff.skimmedBCs.cleanBC();
  clock_t start = clock();
  //eff.correlate(0,0);
  eff.correlateAll(5000);
  clock_t stop = clock();
  double_t time = double_t(stop - start)/ (double_t)CLOCKS_PER_SEC;
  std::cout << "Time:" << time << std::endl;
  return;
  std::cout << "=== Unskimmed" << std::endl;
  //eff.originalBCs.frequencyBC();
  eff.originalBCs.frequencyBCSorted();
  eff.originalBCs.evSel2AOD();
  std::vector<bcInfo> bcs_cleaned;
  eff.originalBCs.cleanBC(bcs_cleaned);
  eff.originalBCs.frequencyBCSorted(bcs_cleaned);
  //eff.originalBCs.printbcInfoVect(bcs_cleaned);
  //eff.originalBCs.selectionEfficiency();
  return ;
  //
  std::cout << "=== Skimmed" << std::endl;
  eff.skimmedBCs.frequencyBCSorted();
  eff.skimmedBCs.evSel2AOD();
  std::vector<bcInfo> bcs_cleaned2;
  eff.skimmedBCs.cleanBC(bcs_cleaned2);
  eff.skimmedBCs.frequencyBCSorted(bcs_cleaned2);
  //
  //corMatrix(originalBCs,originalTotal);

}
