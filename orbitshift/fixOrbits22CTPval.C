void fixOrbits22CTPval() {
  //std::vector<int> runNbVect = {523148};
  std::vector<int> runNbVect = { 523141,523142,523148,523182,523186,523298,523306,523308,523309,523397,523399,523401,523441,523541,523559,523669,523671,523677,523728,523731,523779,523783,523786,523788,523789,523792,523797,523821,526463,526465,526466,526467,526468,526486,526505,526512,526525,526526,526528,526559,526596,526606,526612,526638,526639,526641,526643,526647,526649,526713,526714,526715,526716,526719,526720,526776,526860,526865,526886,526938,526963,526964,526966,526967,526968,527015,527016,527028,527031,527033,527034,527038,527039,527041,527057,527076,527108,527109,527228,527237,527240,527259,527260,527261,527262,527345,527347,527349,527446,527518,527523,527690,527694,527731,527734,527736,527821,527825,527826,527828,527848,527850,527852,527863,527864,527865,527869,527871,527895,527898,527899,527902,527963,527976,527978,527979,528021,528026,528036,528093,528094,528097,528105,528107,528109,528110,528231,528232,528233,528263,528266,528292,528294,528316,528319,528328,528329,528330,528332,528336,528347,528359,528379,528381,528386,528448,528451,528461,528463,528529,528530,528531,528534,528537,528543,528602,528604,528617,528781,528782,528783,528784,528798,528801,529077,529078,529084,529088,529115,529116,529117,529128,529129,529208,529209,529210,529211,529235,529237,529242,529248,529252,529270,529306,529310,529317,529320,529324,529338,529341,529450,529452,529454,529458,529460,529461,529462,529542,529552,529554,529662,529663,529664,529674,529675,529690,529691};
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
  //o2::ccdb::CcdbApi apitest;
  //apitest.init("http://ccdb-test.cern.ch:8080");

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
    std::pair<unsigned long, unsigned long> olr = scl->getOrbitLimitFromRaw();
    std::pair<unsigned long, unsigned long> tlr = scl->getTimeLimitFromRaw();
    std::map<std::string, std::string> metadataOR;
    std::vector<Long64_t>* oreset = ccdb.getSpecific<std::vector<Long64_t>>(ccdbOrbitReset, tstart, metadataOR);
    Long64_t oresetTime = (*oreset)[0] / 1000;
    std::cout << " Orbit reset " << oresetTime << " Time first scalers " << tlr.first << " ccdb: " << pp.first << std::endl;
    unsigned long orbitscalers = olr.first;
    unsigned long orbitfromreset = 0;
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
    // validity of original files
    auto headersS = api.retrieveHeaders(ccdbScalers, metadata, tstart + 2000);
    const auto valSF = headersS.find("Valid-From");
    const auto valSU = headersS.find("Valid-Until");
    long valF = std::stol(valSF->second);
    long valU = std::stol(valSU->second);
    std::cout << "scalers: valS:" << valF << " " << valU << std::endl;
    std::cout << "Comparison between scalers validity and scalers tl first:" << (int) (tlr.first - valF) << " last:" << (int) (valU-tlr.second )<< std::endl;
    // write to ccdb
    std::cout << "Writing to CCDB startValGRP:" << valF << " endValGRP:" << valU << std::endl;
    //apitest.storeAsTFileAny(scl, ccdbScalers, metadata, valF, valU);
    api.storeAsTFileAny(scl, ccdbScalers, metadata, valF, valU);
    ic++;
  }
  std::cout << "Number of runs processed:" << ic << std::endl;
  return;
}


