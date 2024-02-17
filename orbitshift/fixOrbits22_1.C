void fixOrbits22() {
 
  std::vector<int> runNbVect = {523148};
  //std::vector<int> runNbVect = {523148, 523186, 523298, 523306, 523308, 523309, 523397, 523399, 523401, 523441, 523541, 523559, 523669, 523671, 523677, 523728, 523731, 523779, 523783, 523786, 523788, 523789, 523792, 523797, 523821, 523897, 526465, 526466, 526467, 526468, 526486, 526505, 526508, 526510, 526512, 526525, 526526, 526528, 526534, 526559, 526596, 526606, 526612, 526638, 526639, 526641, 526643, 526647, 526649, 526689, 526712, 526713, 526714, 526715, 526716, 526719, 526720, 526776};
  constexpr int LHCMaxBunches = 3564;                              // max N bunches
  constexpr double LHCRFFreq = 400.789e6;                          // LHC RF frequency in Hz
  constexpr double LHCBunchSpacingNS = 10 * 1.e9 / LHCRFFreq;      // bunch spacing in ns (10 RFbuckets)
  constexpr double LHCOrbitNS = LHCMaxBunches * LHCBunchSpacingNS; // orbit duration in ns
  constexpr double LHCOrbitMUS = LHCOrbitNS * 1e-3;                // orbit duration in \mus
  const std::string ccdbOrbitReset = "CTP/Calib/OrbitReset";
  const std::string ccdbScalers = "CTP/Calib/Scalers";
  // read
  o2::ccdb::BasicCCDBManager& ccdb = o2::ccdb::BasicCCDBManager::instance();
  // write
  o2::ccdb::CcdbApi api;
  api.init("http://ccdb-test.cern.ch:8080");
  int ic = 0;
  for (int runNb : runNbVect) {
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
    std::cout << " Orbit reset " << oresetTime << " Time first " << tlr.first << std::endl;
    unsigned long orbitscalers = olr.first;
    unsigned long orbitfromreset = 0;
    if (tlr.first >= oresetTime) {
      orbitfromreset = (tlr.first - oresetTime)/LHCOrbitMUS;
    } else {
      std::cout << " scalter time > orbit reset time, skipping" << std::endl;
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
    std::map<std::string, std::string> metadataGRP;
    metadataGRP["runNumber"] = runNbStr;
    //o2::parameters::GRPECSObject* grp = ccdb.retrieveFromTFileAny<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", metadata, tstart + 20000);
    o2::parameters::GRPECSObject* grp = ccdb.getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", tstart + 20000, metadata);
    //headersGRP = ccdb.retrieveHeaders("GLO/Config/GRPECS", metadataGRP, tstart + 2000);
    auto headersGRP = ccdb.getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", tstart + 2000, metadataGRP);
    //const auto validFromGRP = headersGRP->find("Valid-From");  // start validity
    const auto validFromGRP = headersGRP->getTimeStart();  // start validity
    //const auto validUntilGRP = headersGRP.find("Valid-Until"); // end validity
    const auto validUntilGRP = headersGRP->getTimeEnd();  // start validity
    long startValGRP = validFromGRP;
    std::cout << "startValGRP:" << startValGRP << std::endl;
    //long endValGRP = std::stol(validUntilGRP->second);
    long endValGRP = validUntilGRP;
    std::cout << "endValGRP:" << endValGRP << std::endl;
    // write to ccdb
    //api.storeAsTFileAny(scl, ccdbScalers, metadata, validFromGRP, validUntilGRP);
    ic++;
  }
  std::cout << "Number of runs processed:" << ic << std::endl;
  return;
}
  
  
