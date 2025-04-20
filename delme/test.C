#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsCTP/Configuration.h>
#endif
using namespace o2::ctp;

void test(int runNumber = 559713, uint64_t startInterval = 0, uint64_t endInterval = 0) {
    const double orbitDuration = 88.924596234; // mus
    int debug = 0;
    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    o2::ccdb::CcdbApi ccdb;
    ccdb.init("http://alice-ccdb.cern.ch"); // alice-ccdb.cern.ch
    auto soreor = ccdbMgr.getRunDuration(runNumber);
    uint64_t timeStamp = (soreor.second - soreor.first) / 2 + soreor.first;
    o2::ctp::CTPRunScalers* ctpscalers = nullptr;
    std::map<string, string> metadataCTP;
    std::string runString = std::to_string(runNumber);
    metadataCTP["runNumber"] = runString;
    ctpscalers = ccdbMgr.getSpecific<o2::ctp::CTPRunScalers>("CTP/Calib/Scalers", timeStamp, metadataCTP);
    ctpscalers->convertRawToO2();
    std::pair<unsigned long, unsigned long> pp = ctpscalers->getOrbitLimit();
    std::cout << "orbit limits:" << pp.first << " " << pp.second << std::endl;
    uint64_t orb = 22660704;
    std::cout << ctpscalers->getRate(orb,25,7).second << std::endl;
    uint64_t ts = 1731095698;
    std::cout << ctpscalers->getRateGivenT(ts,25,7).second << std::endl;

}
