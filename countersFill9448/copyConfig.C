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
void copyConfig() {
  std::vector<int> runNbVect = {549286, 549287, 549288, 549289};
  constexpr int LHCMaxBunches = 3564;                              // max N bunches
  constexpr double LHCRFFreq = 400.789e6;                          // LHC RF frequency in Hz
  constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;      // bunch spacing in ns (10 RFbuckets)
  constexpr double LHCOrbitNS = LHCMaxBunches * LHCBunchSpacingNS; // orbit duration in ns
  constexpr double LHCOrbitMUS = LHCOrbitNS * 1e-3;                // orbit duration in \mus
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
  int ic = 0;
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
    auto cfg = apiCCDB.getSpecific<o2::ctp::CTPConfiguration>(ccdbConfig, tstart, metadata);
    // validity of original files
    auto headersS = apiTest.retrieveHeaders(ccdbConfig, metadata, tstart + 2000);
    const auto valSF = headersS.find("Valid-From");
    const auto valSU = headersS.find("Valid-Until");
    long valF = std::stol(valSF->second);
    long valU = std::stol(valSU->second);
    // write to ccdb
    std::cout << "Writing to CCDB startValGRP:" << valF << " endValGRP:" << valU << std::endl;
    //apiProd.storeAsTFileAny(scl, ccdbConfig, metadata, valF, valU);
    ic++;
  }
  std::cout << "Number of runs processed:" << ic << std::endl;
  return;
}


