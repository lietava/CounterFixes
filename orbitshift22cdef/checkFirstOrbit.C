void checkFirstOrbit() {
  //std::vector<int> runNbVect = {523148};
  std::vector<int> runNbVect = {
    520259, 520294, 520471, 520472, 520473,
    517619, 517620, 517623, 517677, 517678, 517679, 517685, 517690, 517693, 517737, 517748, 517751, 517753, 517758, 517767,
    518541, 518543, 518546, 518547,
    519041, 519043, 519045, 519497, 519498, 519499, 519502, 519503, 519504, 519506, 519507, 519903, 519904, 519905, 519906
    };

  constexpr int LHCMaxBunches = 3564;                              // max N bunches
  constexpr double LHCRFFreq = 400.789e6;                          // LHC RF frequency in Hz
  constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;      // bunch spacing in ns (10 RFbuckets)
  constexpr double LHCOrbitNS = LHCMaxBunches * LHCBunchSpacingNS; // orbit duration in ns
  constexpr double LHCOrbitMUS = LHCOrbitNS * 1e-3;                // orbit duration in \mus
  o2::ccdb::BasicCCDBManager& ccdb = o2::ccdb::BasicCCDBManager::instance();
  for (int runNb : runNbVect) {
   if(runNb <= 526802)
   {
    std::cout << "run: " << runNb << std::endl;
    std::pair<int64_t, int64_t> pp = ccdb.getRunDuration(runNb);
    std::map<std::string, std::string> metadata;
    std::string runNbStr = std::to_string(runNb);
    metadata["runNumber"] = runNbStr;
    long ts = pp.first + 60 * 1000;
    if(runNb <= 526802) {
      ccdb.setURL("http://ccdb-test.cern.ch:8080");
    }
    std::cout << "URL scalers:" << ccdb.getURL() << std::endl;
    o2::ctp::CTPRunScalers *scl = nullptr;
    scl = ccdb.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", ts, metadata);
    if(scl == nullptr) {
       std::cout << " ===> ERROR: scalers not read" << std::endl;
       continue;
    }
    if(scl->getRunNUmber() != runNb) {
       std::cout << "===> ERROR scaler run number inconsistent" << std::endl;
       continue;
    }
    //scl->convertRawToO2();
    std::pair<unsigned long, unsigned long> ol = scl->getOrbitLimitFromRaw();
    unsigned long orbitFirst = ol.first;
    std::map<std::string, std::string> metadataOR;
    ccdb.setURL("http://alice-ccdb.cern.ch");
    std::cout << "URL orbit reset:" << ccdb.getURL() << std::endl;
    std::vector<Long64_t>* oreset = nullptr;
    oreset = ccdb.getSpecific<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts, metadataOR);
    if(oreset == nullptr) {
       std::cout << " ===> ERROR: orbit reset not read" << std::endl;
       continue;
    }
    Long64_t oresetTime = (*oreset)[0] / 1000;
    std::cout << "Orbit reset: " << oresetTime << std::endl;
    Long64_t startTimeRun_accordingToOrbit = oresetTime + orbitFirst * LHCOrbitMUS / 1000;
    std::cout << runNb << " First orbit scalers: " << orbitFirst << std::endl;
    std::cout << " startTime run                 = " << pp.first << std::endl;
    std::cout << " startTimeRun_accordingToOrbit = " << startTimeRun_accordingToOrbit << std::endl;
    if (std::abs(startTimeRun_accordingToOrbit - pp.first) > 240 * 1000) { // 4 minutes
      std::cout << " FATAL : run " << runNb << " has start time from run and start time from orbit that differ too much (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    } else {
      std::cout << " INFO : run " << runNb << " has start time from run and start time from orbit that differ (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    }
   }
  }

  return;
}


