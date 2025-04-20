#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DataFormatsCTP/Scalers.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsCTP/Configuration.h"
#endif
#include <TFile.h>

const double orbitDuration = 88.924596234; 

void check(TString filename, int run, bool debug = false) {

  TFile* f = new TFile(Form("%s", filename.Data()));
  o2::ctp::CTPRunScalers* scl = (o2::ctp::CTPRunScalers*)f->Get("CTPRunScalers");
  scl->convertRawToO2();
  std::vector<o2::ctp::CTPScalerRecordO2> mScalerRecordO2 = scl->getScalerRecordO2();

  float ir = 0.f;
  long duration = 0;
  
  int n = mScalerRecordO2.size();
  if (n != 0) {
    std::int64_t totScalers = 0;
    std::vector<int64_t> vOrbit;
    std::vector<int64_t> vScaler;
    int i = 0;
    for (auto& record : mScalerRecordO2) {
      if (debug) {
	record.printStream(std::cout);
      }
      std::vector<o2::ctp::CTPScalerO2>& scalers = record.scalers;
      o2::InteractionRecord& intRecord = record.intRecord;
      vOrbit.push_back(intRecord.orbit);
      if (debug) {
	std::cout << i << " orbit = " << intRecord.orbit << " scalers = " << scalers[0].lmBefore << std::endl;
      }
      vScaler.push_back(scalers[0].lmBefore); // use scalers for class 0 (usually TVX). TODO: extract info on class id from trigger config
      totScalers += scalers[0].lmBefore;
      ++i;
    }
    duration = std::round((vOrbit.back() - vOrbit.front()) * orbitDuration * 1e-6); // s
    ir = float(vScaler.back() - vScaler.front()) / duration;
std::cout << "run " << run << " orbit.front = " << vOrbit.front() << " orbit.back = " << vOrbit.back() << "  duration = " << duration << " s scalers = " << vScaler.back() - vScaler.front() << " IR = " << ir << " Hz" << std::endl;
  }

}
