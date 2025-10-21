void fillChainFromAO2D(TChain &chain, const TString &fileName)
{
  TFile file(fileName);
  for (auto key : *file.GetListOfKeys())
  {
    TString keyName = key->GetName();
    if (keyName.Contains("DF_"))
    {
      chain.Add((fileName + "/" + keyName + "/" + chain.GetName()).Data());
    }
  }
}
bool isNumeric(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}
void analMult(TString normalisationFilename = "dataMB/AnalysisResults.root")
{
  //ROOT::EnableImplicitMT();
  //
  std::set<int> runs = {556160,556182,556741,556767,557415,557876,556152,556164,556734,557374,557425,557862,557897,557913,557926};
  int mRunMuber = 0;
  std::map<int,TH1*> centPerRunHist;
  //
  TFile *f = TFile::Open(normalisationFilename);
  TDirectory *dir = (TDirectory*)f->Get("non-prompt-cascade-task/mult");  // your directory name

  TIter next(dir->GetListOfKeys());
  int i = 0;
  int Nloops = 0;
  TKey *key;
  while ((key = (TKey*)next())) {
      bool del = 1;
      TObject *obj = key->ReadObj();
      if (obj->InheritsFrom("TH2")) {
          TH2 *h = (TH2*)obj;
          h->SetDirectory(0);
          //h->Draw();   // or do any analysis
          //gPad->WaitPrimitive();  // wait for click before moving on
          //std::cout << h->GetName() << " " << h->GetEntries() <<  std::endl;
          std::string name = h->GetName();
          if(name.find("MultVsCentZoom") != std::string::npos) {
            int runNumber = 0;
            std::string runS = name.substr(name.size() - 6);
            //std::cout << name << " " << h->GetEntries() << " run:" << runS << std::endl;
            if(isNumeric(runS)) {
              runNumber = std::stoi(runS.c_str());
            }
            if(runs.count(runNumber)) {
              std::cout << runNumber << " in small skimmed" << std::endl;
              centPerRunHist[runNumber] = h->ProjectionX();
              int biny = h->GetYaxis()->FindBin(3100.);
              TH1* hx = h->ProjectionX("hx",biny,biny);
              int maxBin = hx->GetMaximumBin();
              double centAtMax = hx->GetBinCenter(maxBin);
              std::cout << "cent:" << centAtMax << std::endl;
              del = 0;
            }
          }
      }
      if(del) {
        //delete h;
        delete obj;
      }
      if(Nloops !=0 && i > Nloops) break;
      i++;
  }
  TFile output("outputMult.root", "recreate");
  for(auto const &h: centPerRunHist) {
    h.second->Write();
  }

}
