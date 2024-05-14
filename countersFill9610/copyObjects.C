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
void copyObjects(bool cfg = 0, bool scl = 0, bool writeToCCDB = false) {
  std::vector<int> runNbVect = {551394,551398};
  const std::string ccdbScalers = "CTP/Calib/Scalers";
  const std::string ccdbConfig = "CTP/Config/Config";
  // write
  o2::ccdb::CcdbApi apiProd;
  apiProd.init("http://alice-ccdb.cern.ch");
  // read
  o2::ccdb::CcdbApi apiTest;
  apiTest.init("http://ccdb-test.cern.ch:8080");
  //
  o2::ccdb::BasicCCDBManager& apiCCDB = o2::ccdb::BasicCCDBManager::instance();
  int isclt = 0;
  int isclw = 0;
  int icfgt = 0;
  int icfgw = 0;
  for (int runNb : runNbVect) {
    std::cout << "run:" << runNb << std::endl;
    apiCCDB.setURL("http://alice-ccdb.cern.ch");
    std::pair<int64_t, int64_t> pp = apiCCDB.getRunDuration(runNb);
    std::map<std::string, std::string> metadata;
    std::string runNbStr = std::to_string(runNb);
    metadata["runNumber"] = runNbStr;
    long tstart = pp.first + 60 * 1000;
    std::cout << "Run start:" << tstart << std::endl;
    apiCCDB.setURL("http://ccdb-test.cern.ch:8080");
    if (scl) {
      auto scl = apiCCDB.getSpecific<o2::ctp::CTPRunScalers>(ccdbScalers,
                                                             tstart, metadata);
      // validity of original files
      auto headersS =
          apiTest.retrieveHeaders(ccdbScalers, metadata, tstart + 2000);
      const auto valSF = headersS.find("Valid-From");
      const auto valSU = headersS.find("Valid-Until");
      long valF = std::stol(valSF->second);
      long valU = std::stol(valSU->second);
      // write to ccdb
      if (writeToCCDB) {
        std::cout << "Scalers writing to CCDB startValGRP:" << valF
                  << " endValGRP:" << valU << std::endl;
        apiProd.storeAsTFileAny(scl, ccdbScalers, metadata, valF, valU);
        isclw++;
      } else {
        std::cout << "Scalers test (not writing) startValGRP:" << valF
                  << " endValGRP:" << valU << std::endl;
        isclt++;
      }
    }
    //
    if (cfg) {
      auto cfg = apiCCDB.getSpecific<o2::ctp::CTPConfiguration>(
          ccdbConfig, tstart, metadata);
      // validity of original files
      auto headersS =
          apiTest.retrieveHeaders(ccdbConfig, metadata, tstart + 2000);
      const auto valSF = headersS.find("Valid-From");
      const auto valSU = headersS.find("Valid-Until");
      long valF = std::stol(valSF->second);
      long valU = std::stol(valSU->second);
      // write to ccdb
      if (writeToCCDB) {
        std::cout << "Config writing to CCDB startValGRP:" << valF
                  << " endValGRP:" << valU << std::endl;
        apiProd.storeAsTFileAny(cfg, ccdbConfig, metadata, valF, valU);
        icfgw++;
      } else {
        std::cout << "Config test (not writing) startValGRP:" << valF
                  << " endValGRP:" << valU << std::endl;
        icfgt++;
      }
    }
  }
  std::cout << "Number of runs in the list:" << runNbVect.size() << std::endl;
  std::cout << "Number of scalers only tested (NOT written to ccdb):" << isclt
            << " written to CCDB:" << isclw << std::endl;
  std::cout << "Number of configs only tested (NOT written to ccdb):" << icfgt
            << " written to CCDB:" << icfgw << std::endl;
  return;
}
