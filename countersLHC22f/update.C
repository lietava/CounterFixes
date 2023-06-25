const double orbitDuration = 88.924596234; // us

void update(int period = 0, bool debug = false) {

  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch"); // alice-ccdb.cern.ch
  auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    
  int runNumber = -1;
  std::vector<int> runs = {517619, 517620, 517623, 517677, 517678, 517679, 517685, 517690, 517693, 517737, 517748, 517751, 517753, 517758, 517767, 518541, 518543, 518546, 518547, 519041, 519043, 519045, 519497, 519498, 519499, 519502, 519503, 519504, 519506, 519507, 519903, 519904, 519905, 519906, 520259, 520294, 520471, 520472, 520473, 523141, 523142, 523148, 523182, 523186, 523298, 523306, 523308, 523309, 523397, 523399, 523401, 523441, 523541, 523559, 523669, 523671, 523677, 523728, 523731, 523779, 523783, 523786, 523788, 523789, 523792, 523797, 523821, 523897, 526463, 526465, 526466, 526467, 526468, 526486, 526505, 526512, 526525, 526526, 526528, 526559, 526596, 526606, 526612, 526638, 526639, 526641, 526643, 526647, 526649, 526689, 526713, 526714, 526715, 526716, 526719, 526720, 526776, 526860, 526865, 526886, 526938, 526963, 526964, 526966, 526967, 526968, 527015, 527016, 527028, 527031, 527033, 527034, 527038, 527039, 527041, 527057, 527076, 527108, 527109, 527228, 527237, 527240, 527259, 527260, 527261, 527262, 527345, 527347, 527349, 527446, 527518, 527523, 527671, 527690, 527694, 527731, 527734, 527736, 527777, 527799, 527821, 527825, 527826, 527828, 527848, 527850, 527852, 527863, 527864, 527865, 527869, 527871, 527895, 527898, 527899, 527902, 527976, 527978, 527979, 528021, 528026, 528036, 528093, 528094, 528097, 528105, 528107, 528109, 528110, 528231, 528232, 528233, 528263, 528266, 528292, 528294, 528316, 528319, 528328, 528329, 528330, 528332, 528336, 528347, 528359, 528379, 528381, 528386, 528448, 528451, 528461, 528463, 528529, 528530, 528531, 528534, 528537, 528543, 528602, 528604, 528617, 528781, 528782, 528783, 528784, 528798, 528801, 528991, 528997, 529003, 529005, 529006, 529009, 529015, 529035, 529037, 529038, 529039, 529043, 529066, 529067, 529077, 529078, 529084, 529088, 529115, 529116, 529117, 529128, 529129, 529208, 529209, 529210, 529211, 529235, 529237, 529242, 529248, 529252, 529270, 529306, 529317, 529320, 529324, 529337, 529341, 529397, 529399, 529403, 529414, 529418, 529450, 529452, 529454, 529458, 529460, 529461, 529462, 529542, 529552, 529554, 529610, 529662, 529663, 529664, 529674, 529675, 529690, 529691};

  for (int i = 0; i < runs.size(); ++i) {
    runNumber = runs[i];
    if (runNumber != 520294) {
      continue;
    }
    if (period > 1) {
      if (runNumber < 527799) {
	continue;
      }
    }
    else {
      if (runNumber >= 527799) {
	continue;
      }
    }

    std::cout << std::endl << "processing run " << runNumber << std::endl;

    uint64_t startRunCTP = 0;
    uint64_t endRunCTP = 0;
    TString period = "ppp";
    
    if (runNumber >= 517616  && runNumber <= 517767) {
      period = "LHC22c";
    }
    else if (runNumber >= 518541 && runNumber <= 518547) {
      period = "LHC22d";
    }
    else if (runNumber >= 519041  && runNumber <= 520099) {
      period = "LHC22e";
    }
    else if (runNumber >= 520143 && runNumber <= 520473) {
      period = "LHC22f";
    }
    else if (runNumber >= 523141 && runNumber <= 523898) {
      period = "LHC22m";
    }
    else if (runNumber >= 526383  && runNumber <= 527523) {
      period = "LHC22o";
    }
    else if (runNumber >= 527671  && runNumber <= 527777) {
      period = "LHC22o-test";
    }
    else if (runNumber >= 527799  && runNumber <= 528543) {
      period = "LHC22o";
    }
    else if (runNumber >= 528563 && runNumber <= 528801) {
      period = "LHC22p";
    }
    else if (runNumber >= 528991 && runNumber <= 529043) {
      period = "LHC22q";
    }
    else if (runNumber >= 529066 && runNumber <= 529341) {
      period = "LHC22r";
    }
    else if (runNumber >= 529397 && runNumber <= 529418) {
      period = "LHC22s";
    }
    else if (runNumber >= 529450 && runNumber <= 529691) {
      period = "LHC22t";
    }  
    
    TGrid::Connect("alien");
    int tfIdx = 1;
    if (runNumber == 519503 || runNumber == 527039|| runNumber == 527865) {
      tfIdx = 2;
    }
    TGridResult* gr = gGrid->Command(Form("find /alice/data/2022/%s/%d/raw o2_ctf_*tf000000000%d*root", period.Data(), runNumber, tfIdx));
    TGridResult* gr1; 
    TString path = gr->GetKey(0, "turl");
    TString path1;
    
    // periods 22cde have some CTFs in "JUN" as period
    if (period == "LHC22c" || period == "LHC22d" || period == "LHC22e") {
      gr1 = gGrid->Command(Form("find /alice/data/2022/JUN/%d/raw o2_ctf_*tf000000000%d*root", runNumber, tfIdx));
      path1 = gr1->GetKey(0, "turl");
    }
    if (path.Length() == 0 && path1.Length() == 0) {
      std::cout << "Impossible to find first CTF file, cannot check this run" << std::endl;
      return;
    }
    else {
      std::cout << "path = " << path << std::endl;
      std::cout << "path1 = " << path1 << std::endl;
    }
    int firstOrbitCTFfile = -1;
    int firstOrbitCTFfile1 = -1;
    if (period != "LHC22c") {
      TPRegexp("alien:///alice/data/2022/LHC22\[a-z]/\[0-9]+/raw/\[0-9]+/o2_ctf_run\[0-9]+_orbit(\[0-9]+)_tf000000000\[0-9]_epn\[0-9]+.root").Substitute(path, "$1");
      firstOrbitCTFfile = path.Atoi();  

      if (path1.Length() != 0) {
	TPRegexp("alien:///alice/data/2022/JUN/\[0-9]+/raw/\[0-9]+/o2_ctf_run\[0-9]+_orbit(\[0-9]+)_tf000000000\[0-9]_epn\[0-9]+.root").Substitute(path1, "$1");
	firstOrbitCTFfile1 = path1.Atoi();
      }
    }
    else {
      TPRegexp("alien:///alice/data/2022/LHC22\[a-z]/\[0-9]+/raw/\[0-9]+/o2_ctf_run\[0-9]+_orbit(\[0-9]+)_tf000000000\[0-9].root").Substitute(path, "$1");
      firstOrbitCTFfile = path.Atoi();  

      if (path1.Length() != 0) {
	TPRegexp("alien:///alice/data/2022/JUN/\[0-9]+/raw/\[0-9]+/o2_ctf_run\[0-9]+_orbit(\[0-9]+)_tf000000000\[0-9].root").Substitute(path1, "$1");
	firstOrbitCTFfile1 = path1.Atoi();
      }
    }
      
    std::cout << "firstOrbitCTFfile = " << firstOrbitCTFfile << " firstOrbitCTFfile1 = " << firstOrbitCTFfile1 << std::endl;
    // if part of the CTFs were stored in "JUN", we take the smallest orbit as a start
    if (firstOrbitCTFfile1 > 0 && firstOrbitCTFfile1 < firstOrbitCTFfile) {
      firstOrbitCTFfile = firstOrbitCTFfile1;
    }

    std::cout << "Final firstOrbitCTFfile = " << firstOrbitCTFfile << std::endl;
    //    continue;

    auto soreor = ccdbMgr.getRunDuration(runNumber);
    uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
    std::cout << "Processing run " << runNumber << " which started at " << soreor.first << " and ended at " << soreor.second << std::endl;
    std::cout << " query will be done for timeStamp " << timeStamp << std::endl;
    o2::parameters::GRPECSObject* glo = ccdbMgr.getForTimeStamp<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", timeStamp);

    if (runNumber >= 519041 && runNumber <= 519507) { 
      ccdbMgr.setURL("http://ccdb-test.cern.ch:8080");
    }
    //auto scl = ccdbMgr.getForTimeStamp<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", timeStamp);
    std::string srun = std::to_string(runNumber);
    map<string, string> metadata; // can be empty
    metadata["runNumber"] = srun;
    auto scl = ccdbMgr.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", timeStamp, metadata);
    std::cout << " scl:" << scl << std::endl;
    if (runNumber >= 519041 && runNumber <= 519507) { 
      ccdbMgr.setURL("http://alice-ccdb.cern.ch");
    }
    
    auto* orbitReset = ccdbMgr.getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", timeStamp);
    int64_t orbitResetMUS = (*orbitReset)[0];
    std::cout << "orbitReset [mus] = " << orbitResetMUS << std::endl;
    scl->convertRawToO2();
    std::vector<o2::ctp::CTPScalerRecordO2> mScalerRecordO2 = scl->getScalerRecordO2();
    int n = mScalerRecordO2.size();
    std::cout << " mScalerRecordO2 size:" << n << std::endl;
    std::vector<int64_t> vOrbit;
    std::vector<int64_t> vScaler;
    if (n != 0) {
      std::int64_t totScalers = 0;
      int i = 0;
      for (auto& record : mScalerRecordO2){
	//if (debug) {
	  //record.printStream(std::cout);
	//}
	std::vector<o2::ctp::CTPScalerO2>& scalers = record.scalers;
	o2::InteractionRecord& intRecord = record.intRecord;
	vOrbit.push_back(intRecord.orbit);
	std::cout << " i: " << i << " " << scalers.size() << std::endl;
	//scalers.at(i).printStream(std::cout);
	//if (debug) {
	//std::cout << i << " orbit = " << intRecord.orbit << " scalers = " << scalers[0].lmBefore << std::endl;
	//}
	//vScaler.push_back(scalers[0].lmBefore); // use scalers for class 0 (usually TVX). TODO: extract info on class id from trigger config
	//totScalers += scalers[0].lmBefore;
	//std::cout << " qqq " << std::endl;
	++i;
      }
    }
    //    std::cout << " rrr " << std::endl;
    int64_t startOrbit = vOrbit.front();
    int64_t endOrbit = vOrbit.back();
    //std::cout << " sss " << std::endl;
    int64_t durationInOrbits = endOrbit - startOrbit;
    std::cout << "startOrbit = " << startOrbit << " endOrbit = " << endOrbit << std::endl;
    std::cout << "firstOrbitCTFfile = " << firstOrbitCTFfile << std::endl;
    if (firstOrbitCTFfile != startOrbit) {
      std::cout << "startOrbit in Scalers seems wrong, we will anyway use the one from CTF files" << std::endl;
    }
    startRunCTP = (uint64_t)firstOrbitCTFfile * orbitDuration * 1e-3 + (uint64_t)orbitResetMUS * 1e-3;
    endRunCTP = (uint64_t)(firstOrbitCTFfile + durationInOrbits) * orbitDuration * 1e-3 + (uint64_t)orbitResetMUS * 1e-3;
    
    std::cout << "startRunCTP = " << startRunCTP << " endRunCTP = " << endRunCTP << std::endl;
    std::cout << "Adding these to GRPECS" << std::endl;    
    glo->setTimeStartCTP(startRunCTP);
    glo->setTimeEndCTP(endRunCTP);

    uint64_t startTimeGLO = glo->getTimeStart();
    uint64_t endTimeGLO = glo->getTimeEnd();

    if (std::abs(int(startTimeGLO - startRunCTP)) > 60000 * 5) {
      std::cout << "Large difference between GLO start (" << startTimeGLO << ") and what we get from CTP (" << startRunCTP << " --> " << std::abs(int(startTimeGLO - startRunCTP)) / 1000 /60 << " min, PLEASE CHECK" << std::endl;
    }

    if (std::abs(int(endTimeGLO - endRunCTP)) > 60000 * 5) {
      std::cout << "Large difference between GLO end (" << endTimeGLO << ") and what we get from CTP (" << endRunCTP << " --> " << std::abs(int(endTimeGLO - endRunCTP)) / 1000 /60 << " min, PLEASE CHECK" << std::endl;
    }
	
    std::map<std::string, std::string> md1;
    std::string runString = std::to_string(runNumber);
    md1["runNumber"] = runString;
    std::map<std::string, std::string> hd = api.retrieveHeaders("GLO/Config/GRPECS", md1, timeStamp);
    const auto validFromGRP = hd.find("Valid-From");
    const auto validUntilGRP = hd.find("Valid-Until");
    long startValGRP = std::stol(validFromGRP->second);
    long endValGRP = std::stol(validUntilGRP->second);
    std::cout << "GRPECS will be valid from " << startValGRP << " to " << endValGRP << std::endl;
    
    // uploading new GLO
    std::map<std::string, std::string> md;
    md["JIRA"] = "O2-3826";
    md["runNumber"] = runString;  
    //api.storeAsTFileAny(glo, "GLO/Config/GRPECS", md, startValGRP, endValGRP);
    
    std::map<std::string, std::string> hdRCT = api.retrieveHeaders("RCT/Info/RunInformation", std::map<std::string, std::string>(), runNumber);
    const auto startRCT = hdRCT.find("SOR");
    const auto endRCT = hdRCT.find("EOR");
    std::map<std::string, std::string> mdRCT;
    mdRCT["JIRA"] = "O2-3826";
    mdRCT["SOR"] = startRCT->second;
    mdRCT["EOR"] = endRCT->second;
    mdRCT["SOX"] = std::to_string(startRunCTP);
    mdRCT["EOX"] = std::to_string(endRunCTP);
    std::cout << "Adding SOX = " << startRunCTP << " and EOX = " << endRunCTP << " to RCT" << std::endl; 
    //api.storeAsTFileAny(glo, "RCT/Info/RunInformation", mdRCT, runNumber, runNumber + 1);
  }
}
