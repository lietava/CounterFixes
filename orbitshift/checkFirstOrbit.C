void checkFirstOrbit() {
  //std::vector<int> runNbVect = {523148};
  std::vector<int> runNbVect = { 523141,523142,523148,523182,523186,523298,523306,523308,523309,523397,523399,523401,523441,523541,523559,523669,523671,523677,523728,523731,523779,523783,523786,523788,523789,523792,523797,523821,526463,526465,526466,526467,526468,526486,526505,526512,526525,526526,526528,526559,526596,526606,526612,526638,526639,526641,526643,526647,526649,526713,526714,526715,526716,526719,526720,526776,526860,526865,526886,526938,526963,526964,526966,526967,526968,527015,527016,527028,527031,527033,527034,527038,527039,527041,527057,527076,527108,527109,527228,527237,527240,527259,527260,527261,527262,527345,527347,527349,527446,527518,527523,527690,527694,527731,527734,527736,527821,527825,527826,527828,527848,527850,527852,527863,527864,527865,527869,527871,527895,527898,527899,527902,527963,527976,527978,527979,528021,528026,528036,528093,528094,528097,528105,528107,528109,528110,528231,528232,528233,528263,528266,528292,528294,528316,528319,528328,528329,528330,528332,528336,528347,528359,528379,528381,528386,528448,528451,528461,528463,528529,528530,528531,528534,528537,528543,528602,528604,528617,528781,528782,528783,528784,528798,528801,529077,529078,529084,529088,529115,529116,529117,529128,529129,529208,529209,529210,529211,529235,529237,529242,529248,529252,529270,529306,529310,529317,529320,529324,529338,529341,529450,529452,529454,529458,529460,529461,529462,529542,529552,529554,529662,529663,529664,529674,529675,529690,529691};		

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
      //ccdb.setURL("http://ccdb-test.cern.ch:8080");
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
    std::cout << " First orbit scalers: " << orbitFirst << std::endl;
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
  
  
