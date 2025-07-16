#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include <vector>
#include <iostream>
int test()
{
    auto & cc = o2::ccdb::BasicCCDBManager::instance();
    std::map<std::string, std::string> metadata;
    metadata["runNumber"]="535645";
    auto forbit = cc.getSpecific<std::vector<long>>("CTP/Calib/FirstRunOrbit",1683020656031,metadata);
    auto orbitr = cc.getSpecific<std::vector<long>>("CTP/Calib/OrbitReset",1683020656031);
    for(int i = 0;i < 3; i++) {
        std::cout << "first orbit:" << (*forbit)[i] << " orbit reset:" << (*orbitr)[i] << std::endl;
    }
    return 0;
}
