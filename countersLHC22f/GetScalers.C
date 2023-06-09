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

/// \file TestCTPScalers.C
/// \brief create CTP scalers, test it and add to database
/// \author Roman Lietava
// root -b -q "GetScalers.C(\"519499\", 1656286373953)"
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
using namespace o2::ctp;
void GetScalers(std::string srun, long time, std::string ccdbHost = "http://ccdb-test.cern.ch:8080")
{
  std::map<std::string, std::string> metadata;
  metadata["runNumber"] = srun;
  // auto hd = cdb.retrieveHeaders("RCT/Info/RunInformation", {}, runNumber);
  // auto hd = cdb.retrieveHeaders("RCT/Info/RunInformation", metadata);
  // std::cout << stol(hd["SOR"]) << "\n";
  CTPConfiguration ctpcfg;
  CTPRunScalers scl;
  CTPRunManager mng;
  mng.setCCDBHost(ccdbHost);
  bool ok;
  // ctpcfg = mng.getConfigFromCCDB(time, srun);
  // ctpcfg.printStream(std::cout);
  // return;
  scl = mng.getScalersFromCCDB(time, srun, ok);
  if (ok == 1) {
    scl.convertRawToO2();
    //scl.printO2(std::cout);
    //scl.printFromZero(std::cout);
    scl.printIntegrals();
    //scl.printRates();
    std::vector<o2::ctp::CTPScalerRecordO2> mScalerRecordO2 = scl.getScalerRecordO2();
    int n = mScalerRecordO2.size();
    std::vector<int64_t> vOrbit;
    std::vector<int64_t> vScaler;
    if (n != 0) {
	std::int64_t totScalers = 0;
	int i = 0;
	for (auto& record : mScalerRecordO2){
	   std::vector<o2::ctp::CTPScalerO2>& scalers = record.scalers;
	   o2::InteractionRecord& intRecord = record.intRecord;
	   vOrbit.push_back(intRecord.orbit);
	   if(scalers.size()) {
	   	vScaler.push_back(scalers[0].lmBefore);
		std::cout << i << " ok" << std::endl;
	   	++i;
	   } else {
		std::cout << "I:" << i << " scalers size 0" << std::endl;
	   }
	}
}
  } else {
    std::cout << "Can not find run, please, check parameters" << std::endl;
  }
}
