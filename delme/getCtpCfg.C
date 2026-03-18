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

/// \file getCtpCfg.C
/// \brief get ctpcfg from CCDB
/// \author Roman Lietava
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Configuration.h"
#include <string>
#include <map>
#include <iostream>
#endif
using namespace o2::ctp;
void getCtpCfg(int runNumber = 566888, uint64_t timeStamp = 1760205014264+10)
{ 
  std::string ccdbTest = "http://ccdb-test.cern.ch:8080";
  std::string ccdbProd = "http://alice-ccdb.cern.ch";
  std::string mCCDBPathCtpCfg = "CTP/Config/CtpCfg";
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  ccdbMgr.setURL(ccdbTest);
  if(timeStamp == 0) {
    // Timestamp
    auto soreor = ccdbMgr.getRunDuration(runNumber);
    uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  }
  std::cout << "Timestamp:" << timeStamp << std::endl;
  if(timeStamp == 0) {
    return;
  }
  // CtpCfg
  std::map<std::string,std::string> metadata;
  std::string srun = std::to_string(runNumber);
  metadata["runNumber"] = srun;
  //ccdbMgr.setURL(ccdbTest);
  auto ctpcfg = ccdbMgr.getSpecific<CtpCfg>(mCCDBPathCtpCfg, timeStamp, metadata);
  if (ctpcfg == nullptr) {
    LOG(info) << "CtpCfg not in database, timestamp:" << timeStamp;
    return;
  }
  std::vector<int> inputList = ctpcfg->listOfUsedInputs();
  for(auto const& inp: inputList) {
     std::cout << inp << " ";
  }
  std::cout << std::endl;
}
