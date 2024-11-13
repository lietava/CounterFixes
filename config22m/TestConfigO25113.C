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

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <string>
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsCTP/Configuration.h>
#include "CTPWorkflowScalers/ctpCCDBManager.h"
#endif
using namespace o2::ctp;

void TestConfigO25113()
{
  const std::string ccdbConfig = "CTP/Config/Config";
  o2::ctp::ctpCCDBManager::setCCDBHost("https://alice-ccdb.cern.ch");
  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");
  o2::ccdb::CcdbApi apitest;
  apitest.init("http://ccdb-test.cern.ch:8080");
  //uint64_t timestamp = 1660276134898+10000;
  //std::string run = "523186";
  std::vector<int> runs = {523142, 523148, 523182, 523186, 523298, 523306, 523308, 523309, 523669};
  for(auto runNumber: runs) {
    std::string run = std::to_string(runNumber);
    std::cout << "Doing run:" << run << std::endl;
    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    auto soreor = ccdbMgr.getRunDuration(runNumber);
    uint64_t timestamp = (soreor.second - soreor.first) / 2 + soreor.first;
    bool ok;
    auto ctpcfg = o2::ctp::ctpCCDBManager::getConfigFromCCDB(timestamp, run, ok);
    if (ok == 0) {
     std::cout << "Can not get config for run:" << run << std::endl;
    }
    CTPConfiguration ctpconfig;
    ctpconfig.loadConfigurationRun3(ctpcfg.getConfigString());
    //ctpconfig.printStream(std::cout);
    auto& triggerclasses = ctpconfig.getCTPClasses();
    std::cout << "Found " << triggerclasses.size() << " trigger classes" << std::endl;
    for (const auto& trg : triggerclasses) {
      if (trg.cluster->maskCluster[o2::detectors::DetID::EMC]) {
        // Class triggering EMCAL cluster
        LOG(info) << "Found trigger class for EMCAL cluster: " << trg.name << " with input mask " << std::bitset<64>(trg.descriptor->getInputsMask());
        //trg.descriptor->getInputsMask();
      }
    }
    //
    // validity of original files
    //
    std::map<std::string, std::string> metadata;
    metadata["runNumber"] = run;
    auto headersS = api.retrieveHeaders(ccdbConfig, metadata, timestamp);
    const auto valSF = headersS.find("Valid-From");
    const auto valSU = headersS.find("Valid-Until");
    long valF = std::stol(valSF->second);
    long valU = std::stol(valSU->second);
    std::cout << "scalers: valS:" << valF << " " << valU << std::endl;
    //
    // write to ccdb
    std::cout << "Writing to CCDB startValGRP:" << valF << " endValGRP:" << valU << std::endl;
    //apitest.storeAsTFileAny(&ctpconfig, ccdbConfig, metadata, valF, valU);
    api.storeAsTFileAny(&ctpconfig, ccdbConfig, metadata, valF, valU);
  }
  return;
}
