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

/// \file TestCTPScalers.C
/// \brief create CTP scalers, test it and add to database
/// \author Roman Lietava
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsCTP/Configuration.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <map>
#include <iostream>
#endif
using namespace o2::ctp;
//void GetAndSave(std::string ccdbHost = "http://ccdb-test.cern.ch:8080")
void GetAndSave(std::string ccdbHost = ": http://alice-ccdb.cern.ch/")
{
  std::string CCDBPathCTPScalers = "CTP/Calib/Scalers";
  std::string CCDBPathCTPConfig = "CTP/Config/Config";
  std::vector<string> runs = {""};
  std::vector<long> timestamps = {1659363659341}; // scalers
  int i = 0;
  bool doscalers = 0;
  bool doconfig = 1;
  CTPRunManager mng;
  // mng.setCCDBHost(ccdbHost);
  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL(ccdbHost);
  for (auto const& run : runs) {
    CTPConfiguration ctpcfg;
    CTPRunScalers scl;
    map<string, string> metadata; // can be empty
    metadata["runNumber"] = run;
    if(doscalers) {
      CTPRunScalers* ctpscalers = mgr.getSpecific<CTPRunScalers>(CCDBPathCTPScalers, timestamps[i], metadata);
      if (ctpscalers == nullptr) {
        std::cout << run << " CTPRunScalers not in database, timestamp:" << timestamps[i] << std::endl;
      } else {
        ctpscalers->printStream(std::cout);
        //ctpscalers->convertRawToO2();
        std::string name = run + ".root";
        TFile* myFile = TFile::Open(name.c_str(), "RECREATE");
        myFile->WriteObject(ctpscalers, "CTPRunScalers");
        // myFile->Write();
        std::cout << run << " ok" << std::endl;
      }
   }
   if(doconfig) {
     CTPConfiguration* ctpcfg = mgr.getSpecific<CTPConfiguration>(CCDBPathCTPConfig,timestamps[i], metadata);
     if (ctpcfg == nullptr) {
       std::cout << "CTP Config not found" << std::endl;
     } else {
       std::string name = run + "_cfg.root";
       TFile* myFile = TFile::Open(name.c_str(), "RECREATE");
       myFile->WriteObject(ctpcfg, "CTPConfig");
     }
   }
  i++;
  }
}
