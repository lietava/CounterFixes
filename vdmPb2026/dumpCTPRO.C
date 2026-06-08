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

/// \file checkCTPDigits.C
/// \brief create CTP config, test it and add to database
/// \author Roman Lietava

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <vector>
#include "TKey.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "DataFormatsCTP/Configuration.h"
#endif
// Check if trigger class mask has corresponding input class mask
// Tp be generalised to use CTP config
using namespace o2::ctp;
int dumpCTPRO(bool files = 1)
{
  if (files == 0) {
    return 0;
  }
  const char* ctpinputs[48] = { " T0A", " T0C", " TVX", " TSC", " TCE", " VBA", " VOR", " VIR", " VNC", " VCH", "11", "12", " UCE", "DMC", " USC", " UVX", " U0C", " U0A", "COS", "LAS", "EMC", " PH0", "23", "24", "ZED", "ZNC", "PHL", "PHH", "PHM", "30", "31", "32", "33", "34", "35", "36", "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DG1", "DJ2", "DG2", "45", "46", "47", "48"};
  
  // CTP digits
  TFile* fileDigits = TFile::Open("ctpdigits.root");
  //
  fileDigits->ls();
  o2::ctp::CTPDigit* dig = new o2::ctp::CTPDigit;
  //
  TTreeReader reader("o2sim", fileDigits);
  TTreeReaderArray<o2::ctp::CTPDigit> ctpdigs(reader, "CTPDigits");
  bool firstE = true;
  //
  int ORBIT = 3564;
  int ninps = 0;
  int nClass = 0;
  uint32_t orbitmax = 0;
  uint32_t orbitmin = 0xffffffff;
  while (reader.Next()) {
    if (ctpdigs.GetSetupStatus() < 0) {
      std::cout << "Error:" << std::dec << ctpdigs.GetSetupStatus() << " for:" << ctpdigs.GetBranchName() << std::endl;
      return 1;
    }
    //std::cout << "size:" << std::dec << ctpdigs.GetSize() << std::endl;
    size_t i;
    for (i = 0; i < ctpdigs.GetSize(); i++) {
      o2::ctp::CTPDigit* dig = &ctpdigs[i];
      uint64_t inpMask = dig->CTPInputMask.to_ullong();
      uint64_t trgMask = dig->CTPClassMask.to_ullong();
      auto ir2 = dig->intRecord;
      //if(trgMask && !inpMask) std::cout << "===>"; 
      if(1) {
	std::string inpnames = "";
	for(int j = 0; j < 48; j++) {
	   uint64_t mask = (1ull << j);
	   if(inpMask & mask) {
	     //std::cout << j << " " << ctpinputs[j] << std::endl;		   
             inpnames += ctpinputs[j];
	     inpnames += " ";
	   }
	}
        if(trgMask) std::cout << "TCR:";	
      	std::cout << std::hex << ir2.orbit << " " << ir2.bc << " " << std::hex << inpMask << " " << trgMask << std::dec << " (" << ir2.orbit << " " << ir2.bc << ") ";
        std::cout << inpnames;
	std::cout << std::endl;
      }
    }
  }
  return 0;
}
