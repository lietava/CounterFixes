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
// root -b -q "GetScalers.C(\"519499\", 1656286373953)"
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsCTP/Configuration.h"
#include <string>
#include <map>
#include <iostream>
#endif
using namespace o2::ctp;
void getLumi(int runNumber, uint64_t timeStamp = 0, std::string className= "")
{ //
  std::string ccdbTest = "http://ccdb-test.cern.ch:8080";
  std::string ccdbProd = "http://alice-ccdb.cern.ch";
  std::string mCCDBPathCTPScalers = "CTP/Calib/Scalers";
  std::string mCCDBPathCTPConfig = "CTP/Config/Config";
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
  ccdbMgr.setURL(ccdbProd);
  if(timeStamp == 0) {
    // Timestamp
    auto soreor = ccdbMgr.getRunDuration(runNumber);
    uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
  }
  std::cout << "Timestamp:" << timeStamp << std::endl;
  // Scalers
  std::map<std::string,std::string> metadata;
  std::string srun = std::to_string(runNumber);
  metadata["runNumber"] = srun;
  ccdbMgr.setURL(ccdbTest);
  auto scl = ccdbMgr.getSpecific<CTPRunScalers>(mCCDBPathCTPScalers, timeStamp, metadata);
  std::cout << "Counters size:" << scl->getScalerRecordRaw().size() << std::endl;
  if (scl == nullptr) {
    LOG(info) << "CTPRunScalers not in database, timestamp:" << timeStamp;
    return;
  }
  //scl->printClasses(std::cout);
  //return;
  scl->convertRawToO2();
  std::vector<CTPScalerRecordO2> recs = scl->getScalerRecordO2();
  //scl->printFromZero(std::cout);
  //
  // CTPConfiguration ctpcfg;
  auto ctpcfg = ccdbMgr.getSpecific<CTPConfiguration>(mCCDBPathCTPConfig, timeStamp, metadata);
  if (ctpcfg == nullptr) {
    LOG(info) << "CTPRunConfig not in database, timestamp:" << timeStamp;
    return;
  }
  ctpcfg->printStream(std::cout);
  //ctpcfg->printClasses(std::cout);
  //return;
  std::vector<CTPClass> ctpcls = ctpcfg->getCTPClasses();
  int tvx = 255;
  for (auto const& cls : ctpcls) {
    if (cls.name.find(className) != std::string::npos && tvx == 255) {
      int itvx = cls.getIndex();
      tvx = scl->getScalerIndexForClass(itvx);
      std::cout << cls.name << ":" << tvx << ":" << itvx << std::endl;
    }
  }
  if (tvx == 255) {
    std::cout << "Scalers  for:" << className << " not available, check config to find alternative)" << std::endl;
    return;
  }
  //
  auto orbit0 = recs[0].intRecord.orbit;
  double_t time0 = recs[0].epochTime;
  double_t timeL = recs[recs.size() - 1].epochTime;
  double_t orbitL = recs[recs.size() - 1].intRecord.orbit;
  double_t scal0 = 0;
  double_t scalL = 0;

  if(className.find("0") == std::string::npos) {
    scal0 = recs[0].scalers[tvx].lmBefore ;
    scalL = recs[recs.size() - 1].scalers[tvx].lmBefore;
  } else {
    scal0 = recs[0].scalers[tvx].l0Before ;
    scalL = recs[recs.size() - 1].scalers[tvx].l0Before;
  }
  double_t timeTot = (timeL - time0);
  double_t lumiTot = scalL - scal0;
  double_t arate = lumiTot/timeTot;
  std::cout << "Orbit0:" << orbit0 << " time0:" << (uint64_t)recs[0].epochTime << std::endl;
  std::cout << "OrbitLast:" << orbitL << " timeLast:" << (uint64_t)recs[recs.size()-1].epochTime << std::endl;
  std::cout << "run:" << runNumber << " duration:" << timeTot << " secs Total TVX:" << lumiTot << " Average rate:" << arate << std::endl;
}
