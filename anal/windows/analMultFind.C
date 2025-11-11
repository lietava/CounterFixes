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
void analMultFind(TString normalisationFilename = "dataMB/AnalysisResults.root")
{
  //ROOT::EnableImplicitMT();
  //
  //std::vector<int> runs = {556160,556182,556741,556767,557415,557876,556152,556164,556734,557374,557425,557862,557897,557913,557926};
  std::vector<int> runs = {556160,556182,556741,556767,557415,557876,556152,556164,556734,557374,557425,557862,557897,557926};
  constexpr float CENT = 0.8;
  std::map<float,float> centVsRun; 
  std::map<float,float> tracksVsRun;
  int mRunMuber = 0;
  std::map<int,TH1*> centPerRunHist;
  std::map<int,TH1*> tracksPerRunHist;
  // Track mults for cent bins
  constexpr int NCENTBINS = 11;
  double centbins[NCENTBINS + 1]{0, 0.8,10,20,30,40,50,60,70,80,90,100};
  //double centbins[NCENTBINS + 1]{ 0, 0.01,10,20,30,40,50,60,70,80,90,100 };


  //
  TFile *f = TFile::Open(normalisationFilename);
  std::string dir = "non-prompt-cascade-task/mult/";  // your directory name

  std::map<int, std::array<TH1*, NCENTBINS>> centBinsTracks;
  TH1* test;
  bool onlyone = 1;
  for(auto const& runNumber: runs) {
    std::string histnameM = dir + "hMultVsCentZoom_run" + std::to_string(runNumber);
    TH2 *h = (TH2F*) f->Get(histnameM.c_str());
    if(h == nullptr){
      std::cout << histnameM << " not found. Exiting." << std::endl;
      exit(1);
    }
    std::string histnameB = dir + "hNTracksVsCent_run" + std::to_string(runNumber);
    TH2* hTvsC = (TH2F*)f->Get(histnameB.c_str());
    if (hTvsC == nullptr) {
        std::cout << histnameB << " not found. Exiting." << std::endl;
        exit(1);
    }
    std::string histnameT = dir + "hNTracksVsCentZoom_run" + std::to_string(runNumber);
    TH2 *hTvsCZ = (TH2F*) f->Get(histnameT.c_str());
    if(hTvsCZ == nullptr){
      std::cout << histnameT << " not found. Exiting." << std::endl;
      exit(1);
    }
    if (onlyone) {
        std::cout << " Zoomed x range:" << hTvsCZ->GetXaxis()->GetXmax() << std::endl;
        onlyone = 0;
    }

    //h->Draw();   // or do any analysis
    //gPad->WaitPrimitive();  // wait for click before moving on
    //centPerRunHist[runNumber] = h->ProjectionX();
    std::string namex = "hcent_T0M3100_" + std::to_string(runNumber);
    int biny = h->GetYaxis()->FindBin(3100.);
    TH1* hx = h->ProjectionX(namex.c_str(), biny, biny);
    int maxBin = hx->GetMaximumBin();
    double centAtMax = hx->GetBinCenter(maxBin);
    centVsRun[runNumber] = centAtMax;
    centPerRunHist[runNumber] = hx;
    //
    //tracksPerRunHist[runNumber] = hTvsCZ->ProjectionY();
    int binx = hTvsCZ->GetXaxis()->FindBin(centAtMax);
    TH1* hy = hTvsCZ->ProjectionY("hy",binx,binx);
    tracksPerRunHist[runNumber] = hTvsCZ->ProjectionY();
    float meanTracks  = hy->GetMean();
    tracksVsRun[runNumber] = meanTracks;
    std::cout << runNumber << " cent:" << centAtMax << " <Tracks>:" << meanTracks << std::endl;
    //
    //  <tracks> per Cent
    //
    for (int ic = 0; ic < NCENTBINS; ic++) {
        TH2* h = hTvsC;
        if (centbins[ic + 1] < hTvsCZ->GetXaxis()->GetXmax()) {
            h = hTvsCZ;
        }
        int binx1 = h->GetXaxis()->FindBin(centbins[ic]);
        int binx2 = h->GetXaxis()->FindBin(centbins[ic+1]);
        std::string name = "tracksCent_" + to_string(ic);
        test = h->ProjectionY(name.c_str(), binx1, binx2);
        test->SetName(name.c_str());
        centBinsTracks[runNumber][ic] = test;
        //std::cout << "Tracks:" << test->GetMean() << "[" << centbins[ic] << "-" << centbins[ic+1] << "]" << std::endl;
    }
  }
  //
  // Cent/Tracks versus Run
  //
  TGraph* gCentVsRun = new TGraph(centVsRun.size());
  gCentVsRun->SetTitle("CentVsRun");
  gCentVsRun->SetName("CentVsRun");
  int i = 0;
  for(const auto&[key, val]: centVsRun) {
    gCentVsRun->SetPoint(i++,key,val);
  }
  //
  TGraph* gTracksVsRun = new TGraph(tracksVsRun.size());
  gTracksVsRun->SetTitle("TracksVsRun");
  gTracksVsRun->SetName("TracksVsRun");
  i = 0;
  for(const auto&[key, val]: tracksVsRun) {
    gTracksVsRun->SetPoint(i++,key,val);
  }
  //
  // Output
  //
  TFile output("outputMult.root", "recreate");
  //test->Write();
  gCentVsRun->Write();
  gTracksVsRun->Write();
  for(auto const &h: centPerRunHist) {
    h.second->Write();
  }
  for(auto const &h: tracksPerRunHist) {
    h.second->Write();
  }
  std::array<TH1*, NCENTBINS> total;
  for (auto& h : total) {
      h = nullptr;
  }
  for (auto &[runNumber, a ]: centBinsTracks) {
      //std::cout << runNumber << std::endl;
      for (int i = 0; i < NCENTBINS; i++) {
          TH1* h = a[i];
          if (!total[i]) {
              std::string name = "hTracks_[" + std::to_string(centbins[i]) + "-" + std::to_string(centbins[i+1]) + "]";
              //std::cout << "creating:" << i << std::endl;
              total[i] = (TH1*)h->Clone(name.c_str());
          }
          else {
            total[i]->Add(h);
            //std::cout << "adding:" << i << " " << runNumber << std::endl;
          }
      }
  }
  //
  //percentilesNTracks(total[0], 10);
  //
  std::vector<float> x(NCENTBINS), y(NCENTBINS),ex(NCENTBINS),ey(NCENTBINS);
  for (int i = 0; i < NCENTBINS; i++) {
      if (total[i]) {
        total[i]->Write();
        std::cout << i <<  " <Tracks>:" << total[i]->GetMean() << " rms:" << total[i]->GetRMS() << " emean:" << total[i]->GetMeanError()  << " stat:" << total[i]->GetEntries() << std::endl;
        x[i] = i;
        y[i] = total[i]->GetMean();
        ex[i] = 0;
        ey[i] = total[i]->GetRMS();
      }
  }
  TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetName("NTracks");
  g->SetTitle("NTracks for cent bins");
  g->Write();
}
