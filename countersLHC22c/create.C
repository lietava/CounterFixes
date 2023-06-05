
#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fairlogger/Logger.h>
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsCTP/Configuration.h"
#include <string>
#include <map>
#include <iostream>
#include "TFile.h"
#endif
using namespace o2::ctp;
int readFile(std::vector<int>& runs, std::vector<long>& startTS, std::vector<long>& stopTS, std::vector<uint32_t>& count)
{
   std::string name = "runs.txt";
   fstream in;
   in.open(name);
   std::string line;
   while(getline(in,line)) {
      //std::cout << line << std::endl;
      TString tline(line);
      TObjArray *items = tline.Tokenize("\tab");
      int size = items->GetEntriesFast();
      //std::cout << "size:" << size << std::endl;
      for(int i = 0; i < size; i++) {
	 TString ss = ((TObjString*)items->At(i))->String();     
	 //std::cout << i << " " << ss << std::endl;
      }
      if(size) {
        TString ss = ((TObjString*)items->At(0))->String();
        runs.push_back(ss.Atoi());
        ss = ((TObjString*)items->At(2))->String();
        startTS.push_back(ss.Atoi());
        ss = ((TObjString*)items->At(3))->String();
        stopTS.push_back(ss.Atoi());
        ss = ((TObjString*)items->At(4))->String();
	//std::cout << ss << " " << ss.Atoi() << std::endl;
        count.push_back(ss.Atoi());
      }     
   }
   return 0;
}
void create() 
{
   using namespace std::chrono_literals;
   std::chrono::seconds min5 = 300s;
   long time5min = std::chrono::duration_cast<std::chrono::milliseconds>(min5).count();
   //	
   std::vector<int> runs;
   std::vector<uint32_t> count;
   std::vector<long> startTS;
   std::vector<long> stopTS;
   readFile(runs,startTS,stopTS,count);
   std::cout << "count size:" << count.size() << std::endl;
   int i = 0;
   o2::InteractionRecord ir0 = {0,0};
   CTPScalerRaw sc0 = {0,0,0,0,0,0,0};
   for(auto const c : count) {
      std::cout << c << " " << startTS[i] << std::endl;
      //
      CTPScalerRaw sc = {0,c,c,c,c,c,c};
      // 1st scaler record
      CTPScalerRecordRaw screc0;
      screc0.intRecord = ir0;
      screc0.epochTime = (double_t) (startTS[i] - time5min);
      screc0.scalers.push_back(sc0);
      //  2nd = last scaler record
      CTPScalerRecordRaw screc;
      uint32_t orbit = stopTS[i]/88e-6;
      screc.intRecord = {0,orbit};
      screc.epochTime = (double_t) (stopTS[i] + time5min);
      screc.scalers.push_back(sc);
      // Full scalers consisting of two records
      CTPRunScalers runsc;
      runsc.setRunNumber(runs[i]);
      runsc.setClassMask(0x1);
      runsc.addScalerRacordRaw(screc0);
      runsc.addScalerRacordRaw(screc);
      std::string name = std::to_string(runs[i]) + ".root";
      TFile* myFile = TFile::Open(name.c_str(), "RECREATE");
      myFile->WriteObject(&runsc, "CTPRunScalers");
      i++;
      //if(i == 1) break;
   }
}
