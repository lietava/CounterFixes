#include <iostream>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"

//
void checkFirstOrbit() {
  //std::vector<int> runNbVect = {554523,554524,554526,554537,554538};
  std::vector<int> runNbVect = {554524,554526,554538};

  constexpr int LHCMaxBunches = 3564;                              // max N bunches
  constexpr double LHCRFFreq = 400.789e6;                          // LHC RF frequency in Hz
  constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;      // bunch spacing in ns (10 RFbuckets)
  constexpr double LHCOrbitNS = LHCMaxBunches * LHCBunchSpacingNS; // orbit duration in ns
  constexpr double LHCOrbitMUS = LHCOrbitNS * 1e-3;                // orbit duration in \mus
  o2::ccdb::BasicCCDBManager& ccdb = o2::ccdb::BasicCCDBManager::instance();
  for (uint runNb : runNbVect) {
    std::cout << "run: " << runNb << std::endl;
    std::pair<int64_t, int64_t> pp = ccdb.getRunDuration(runNb);
    long ts = pp.first + 60 * 1000;
    std::cout << "Run start:" << pp.first << std::endl;
    //
    std::map<std::string, std::string> metadata;
    std::string runNbStr = std::to_string(runNb);
    metadata["runNumber"] = runNbStr;
    o2::ctp::CTPRunScalers *scl = nullptr; 
    scl = ccdb.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", ts, metadata);
    if(scl == nullptr) {
       std::cout << " ===> ERROR: scalers not read" << std::endl; 
       continue;
    }
    if(scl->getRunNumber() != runNb) {
       std::cout << "===> ERROR scaler run number inconsistent" << std::endl;
       continue;
    }
    //scl->convertRawToO2();
    std::pair<unsigned long, unsigned long> ol = scl->getOrbitLimitFromRaw();
    unsigned long orbitFirst = ol.first;
    //
    std::vector<Long64_t>* oreset = nullptr;
    oreset = ccdb.getSpecific<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
    if(oreset == nullptr) {
       std::cout << " ===> ERROR: orbit reset not read" << std::endl; 
       continue;
    }
    //
    // First orbit
    auto or1 = ccdb.getSpecific<std::vector<long>>("CTP/Calib/FirstRunOrbit",ts);
    std::cout  << "CTP/Calib/FirstRunOrbit record:"<<  (*or1)[0] << " " << (*or1)[1] << " " << (*or1)[2] << std::endl;
    //
    Long64_t oresetTime = (*oreset)[0] / 1000;
    std::cout << "Orbit reset: " << oresetTime << std::endl;
    Long64_t startTimeRun_accordingToOrbit = oresetTime + orbitFirst * LHCOrbitMUS / 1000;
    std::cout << " First orbit scalers: " << orbitFirst << std::endl;
    std::cout << " First orbit CCDB: " << (*or1)[2] << " diff:" << (*or1)[2] - orbitFirst << std::endl;
    std::cout << " startTime run                 = " << pp.first << std::endl;
    std::cout << " startTimeRun_accordingToOrbit = " << startTimeRun_accordingToOrbit << std::endl;
    if (std::abs(startTimeRun_accordingToOrbit - pp.first) > 240 * 1000) { // 4 minutes
      std::cout << " FATAL : run " << runNb << " has start time from run and start time from orbit that differ too much (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    } else {
      std::cout << " INFO : run " << runNb << " has start time from run and start time from orbit that differ (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    }
  }
  
  return;
}
  
  
