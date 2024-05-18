
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
//

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <array>

#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"

#include "CCDB/BasicCCDBManager.h"

const std::string kBaseCCDBPath = "Users/r/rlietava/EventFiltering/OTS/";

#pragma link C++ class std::vector < std::array < uint64_t, 2>> + ;

void test()
{

  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");
  std::vector<std::array<uint64_t, 2>> testtest;
  testtest.push_back({1ull,2ull});
  std::cout << "vector size:" << testtest.size() << std::endl;
  std::pair<int64_t, int64_t> duration = {10,100};
  std::map<std::string, std::string> metadata;
  metadata["runNumber"] = "333";
  if(0) {
    api.storeAsTFileAny(&testtest, kBaseCCDBPath + "test", metadata, duration.first, duration.second);
  }
  uint64_t timeStamp = duration.first + 10;
  if(0) {
    std::vector<std::array<uint64_t, 2>>* fbm = api.retrieveFromTFileAny<std::vector<std::array<uint64_t, 2>>>(kBaseCCDBPath + "test", metadata, timeStamp);
    std::cout << "vector size:" << (*fbm).size() << std::endl;
    std::cout << (*fbm)[0][0] << " " << (*fbm)[0][1] << std::endl;
  }
  if(1) {
    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    std::vector<std::array<uint64_t, 2>>* fbm = ccdbMgr.getSpecific<std::vector<array<uint64_t,2>>>(kBaseCCDBPath + "test", timeStamp, metadata);
    std::cout << "vector size:" << (*fbm).size() << std::endl;
    std::cout << (*fbm)[0][0] << " " << (*fbm)[0][1] << std::endl;
  }
}
