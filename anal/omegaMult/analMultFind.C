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

  //
  TFile *f = TFile::Open(normalisationFilename);
  std::string dir = "non-prompt-cascade-task/mult/";  // your directory name

  for(auto const& runNumber: runs) {
    std::string histnameM = dir + "hMultVsCentZoom_run" + std::to_string(runNumber);
    TH2 *h = (TH2F*) f->Get(histnameM.c_str());
    if(h == nullptr){
      std::cout << histnameM << " not found. Exiting." << std::endl;
      exit(1);
    }
    std::string histnameT = dir + "hNTracksVsCentZoom_run" + std::to_string(runNumber);
    TH2 *hTvsC = (TH2F*) f->Get(histnameT.c_str());
    if(hTvsC == nullptr){
      std::cout << histnameT << " not found. Exiting." << std::endl;
      exit(1);
    }
    //h->Draw();   // or do any analysis
    //gPad->WaitPrimitive();  // wait for click before moving on
    //centPerRunHist[runNumber] = h->ProjectionX();
    int biny = h->GetYaxis()->FindBin(3100.);
    TH1* hx = h->ProjectionX("hx",biny,biny);
    int maxBin = hx->GetMaximumBin();
    double centAtMax = hx->GetBinCenter(maxBin);
    centVsRun[runNumber] = centAtMax;
    centPerRunHist[runNumber] = hx;
    //
    //tracksPerRunHist[runNumber] = hTvsC->ProjectionY();
    int binx = hTvsC->GetXaxis()->FindBin(centAtMax);
    TH1* hy = hTvsC->ProjectionY("hy",binx,binx);
    tracksPerRunHist[runNumber] = hTvsC->ProjectionY();
    float meanTracks  = hy->GetMean();
    tracksVsRun[runNumber] = meanTracks;
    std::cout << runNumber << " cent:" << centAtMax << " <Tracks>:" << meanTracks << std::endl;
    //
    int binx1 = hTvsC->GetXaxis()->FindBin(centbins[0]);
    int binx2 = hTvsC->GetXaxis()->FindBin(centbins[1]);
    float cent0 = hTvsC->ProjectionX("hp",binx1,binx2)->Integral("width");
    std::cout << "Tracks:" << cent0 << std::endl;
  }
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
  TFile output("outputMult.root", "recreate");
  gCentVsRun->Write();
  gTracksVsRun->Write();
  for(auto const &h: centPerRunHist) {
    h.second->Write();
  }
  for(auto const &h: tracksPerRunHist) {
    h.second->Write();
  }
}
