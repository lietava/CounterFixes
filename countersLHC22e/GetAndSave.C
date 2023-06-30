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
int prepareUpload(std::vector<int>& runs, std::vector<long>& startTS, std::vector<long>& stopTS)
{
   using namespace std::chrono_literals;
   std::chrono::seconds min5 = 300s;
   long time5min = std::chrono::duration_cast<std::chrono::milliseconds>(min5).count();
   FILE *fptr = fopen("command.txt", "w");
   // 
   int i = 0;
   for(auto const& run: runs) {
    long start = startTS[i] - time5min;
    long end = stopTS[i] + time5min;    
    printf("o2-ccdb-upload -f %d.root --starttimestamp %ld --endtimestamp %ld  -k \"CTPRunScalers\" --path CTP/Calib/Scalers --host alice-ccdb.cern.ch -m \"JIRA=O2-3684;runNumber=%d\"\n", run, start, end, run);
    fprintf(fptr, "o2-ccdb-upload -f %d.root --starttimestamp %ld --endtimestamp %ld  -k \"CTPRunScalers\" --path CTP/Calib/Scalers --host alice-ccdb.cern.ch -m \"JIRA=O2-3684;runNumber=%d\"\n", run, start, end, run);
    i++;
   }
   fclose(fptr);
   return 0;
}
void GetAndSave(std::string ccdbHost = "http://ccdb-test.cern.ch:8080")
{
  std::string ccdbHostProd = "http://alice-ccdb.cern.ch/";
  std::string CCDBPathCTPScalers = "CTP/Calib/Scalers";
  std::string CCDBPathCTPConfig = "CTP/Config/Config";
  std::vector<string> runs = {"519041","519043", "519045","519497","519498","519499","519502","519503","519504","519506","519507"};
  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL(ccdbHost);
  o2::ccdb::CcdbApi api;
  api.init(ccdbHostProd);
  //auto hd = api.retrieveHeaders("RCT/Info/RunInformation", std::map<std::string,std::string>(), 519041);
  //return;
  std::vector<long> startTS;
  std::vector<long> stopTS;
  std::vector<int> runis;
  for (auto const& run : runs) {
    int runnumber = std::stoi(run);
    auto hd = api.retrieveHeaders("RCT/Info/RunInformation", std::map<std::string,std::string>(), runnumber);
    long timestampS = 0;
    long timestampE = 0;
    if(hd.size()) {
     timestampS = std::stol(hd["SOR"]);
     timestampE = std::stol(hd["EOR"]);
     std::cout << "Found sor:" << timestampS << " eor:" << timestampE << std::endl;
     // for(auto const& i: hd) {
     //  std::cout << i.first << " " << i.second << std::endl;
     // }
    } else {
      std::cout << " something wrong with " << run << std::endl;
      std::cout << "size:" << hd.size() << std::endl;
      for(auto const& i: hd) {
       std::cout << i.first << " " << i.second << std::endl;
      }
      return;	    
    }
    CTPConfiguration ctpcfg;
    CTPRunScalers scl;
    map<string, string> metadata; // can be empty
    metadata["runNumber"] = run;
    CTPRunScalers* ctpscalers = mgr.getSpecific<CTPRunScalers>(CCDBPathCTPScalers, timestampS, metadata);
    if (ctpscalers == nullptr) {
       std::cout << run << " CTPRunScalers not in database, timestamp sor:" << timestampS << std::endl;
    } else {
        //ctpscalers->printStream(std::cout);
        //ctpscalers->convertRawToO2();
        std::string name = run + ".root";
        TFile* myFile = TFile::Open(name.c_str(), "RECREATE");
        myFile->WriteObject(ctpscalers, "CTPRunScalers");
        // myFile->Write();
        std::cout << run << " ok" << std::endl;
	//
	int runi = std::stoi(run);
	runis.push_back(runi);
	startTS.push_back(timestampS);
	stopTS.push_back(timestampE);
   }
  }
  prepareUpload(runis, startTS, stopTS);
}
