void fixOrbits22CTPcdef() {
  //std::vector<int> runNbVect = {517619};
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
  const std::string ccdbOrbitReset = "CTP/Calib/OrbitReset";
  const std::string ccdbScalers = "CTP/Calib/Scalers";
  // write
  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");
  o2::ccdb::CcdbApi apitest;
  apitest.init("http://ccdb-test.cern.ch:8080");
  int ic = 0;
  for (int runNb : runNbVect) {
    if (runNb > 526802) continue;
    o2::ccdb::BasicCCDBManager& ccdb = o2::ccdb::BasicCCDBManager::instance();
    std::cout << "run:" << runNb << std::endl;
    ccdb.getRunDuration(runNb);
    std::pair<int64_t, int64_t> pp = ccdb.getRunDuration(runNb);
    std::map<std::string, std::string> metadata;
    std::string runNbStr = std::to_string(runNb);
    metadata["runNumber"] = runNbStr;
    long tstart = pp.first + 60 * 1000;
    auto scl = ccdb.getSpecific<o2::ctp::CTPRunScalers>(ccdbScalers, tstart, metadata);
    //scl->convertRawToO2();
    //std::pair<unsigned long, unsigned long> ol = scl->getOrbitLimit();
    //std::pair<unsigned long, unsigned long> tl = scl->getTimeLimit();
    if(runNb >= 517619 && runNb <= 518547) {
      //auto sclraw = scl->getScalerRecordRaw();
      //std::cout << "Scaler size:" << sclraw.size() << std::endl;
      auto ttlim = scl->getTimeLimitFromRaw();
      double_t t1 = (double_t) (ttlim.first/1000000);
      double_t t2 = (double_t) (ttlim.second/1000000);
      scl->setEpochTime(t1,0);
      scl->setEpochTime(t2,1);
    }
    std::pair<unsigned long, unsigned long> olr = scl->getOrbitLimitFromRaw();
    std::pair<unsigned long, unsigned long> tlr = scl->getTimeLimitFromRaw();
    std::map<std::string, std::string> metadataOR;
    std::vector<Long64_t>* oreset = ccdb.getSpecific<std::vector<Long64_t>>(ccdbOrbitReset, tstart, metadataOR);
    Long64_t oresetTime = (*oreset)[0] / 1000;
    std::cout << " Orbit reset " << oresetTime << " Time first scalers " << tlr.first << " ccdbi sor: " << pp.first << std::endl;
    unsigned long orbitscalers = olr.first;
    unsigned long orbitfromreset = 0;
    if(runNb >= 517619 && runNb <= 518547) {
      //tlr.first = tlr.first/1000;
    }
    if (tlr.first >= oresetTime) {
      orbitfromreset = (double) (tlr.first - oresetTime)/LHCOrbitMUS*1000.;
    } else {
      std::cout << "error: scaler time < orbit reset time, skipping" << std::endl;
      continue;
    }
    std::cout << "Orbit from reset:" << orbitfromreset << " orbit from scalers:" << orbitscalers << std::endl;
    uint32_t offset = 0;
    offset = orbitscalers - orbitfromreset;
    std::cout << "orbit offset:" << offset << " " << offset*88e-6 << " secs" << std::endl;
    unsigned long orbitscalers_c = orbitscalers - offset;
    std::cout << "First and last before:" << olr.first << " " << olr.second << std::endl;
    scl->addOrbitOffset(offset);
    scl->getTimeLimitFromRaw();
    std::pair<unsigned long, unsigned long> olrn = scl->getOrbitLimitFromRaw();
    std::cout << "First and last after:" << olrn.first << " " << olrn.second << std::endl;
    // check
    Long64_t startTimeRun_accordingToOrbit = oresetTime + olrn.first * LHCOrbitMUS / 1000;
    if (std::abs(startTimeRun_accordingToOrbit - pp.first) > 240 * 1000) { // 4 minutes
      std::cout << " FATAL : run " << runNb << " has start time from run and start time from orbit that differ too much (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    } else {
      std::cout << " INFO : run " << runNb << " has start time from run and start time from orbit that differ (by " << std::abs(startTimeRun_accordingToOrbit - pp.first) / 1000 << " secs)" << std::endl;
    }
    //
    // validity of original files
    //
    auto headersS = api.retrieveHeaders(ccdbScalers, metadata, tstart + 2000);
    const auto valSF = headersS.find("Valid-From");
    const auto valSU = headersS.find("Valid-Until");
    long valF = std::stol(valSF->second);
    long valU = std::stol(valSU->second);
    std::cout << "scalers: valS:" << valF << " " << valU << std::endl;
    std::cout << "Comparison between scalers validity and scalers tl first:" << (int) (tlr.first - valF) << " last:" << (int) (valU-tlr.second )<< std::endl;
    //
    // write to ccdb
    std::cout << "Writing to CCDB startValGRP:" << valF << " endValGRP:" << valU << std::endl;
    //apitest.storeAsTFileAny(scl, ccdbScalers, metadata, valF, valU);
    //api.storeAsTFileAny(scl, ccdbScalers, metadata, valF, valU);
    ic++;
  }
  std::cout << "Number of runs processed:" << ic << std::endl;
  return;
}


