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
void percentilesNTracks(TH1* tracks, int n = 10, float firstbin = 0)
{
    int Ntot = tracks->GetEntries();
    std::cout << tracks->GetName() << ":" << Ntot << std::endl;
    int lastbin = tracks->GetNbinsX();
    float Ntarget = tracks->Integral(firstbin,lastbin) / n;
    int j = 0;
    float cum = 0;
    float mean = 0;
    float sumw = 0;
    for (int i = firstbin; i <= lastbin; i++) {
        cum += tracks->GetBinContent(i);
        sumw += tracks->GetBinContent(i);
        //std::cout << tracks->GetBinContent(i) << std::endl;
        mean += tracks->GetXaxis()->GetBinCenter(i) * tracks->GetBinContent(i);
        if (cum >= Ntarget) {
            mean = mean / sumw;
            std::cout << j << " " << Ntarget << " bin:" << i << " ntracks:" << tracks->GetXaxis()->GetBinCenter(i) << " mean:" << mean << std::endl;
            j++;
            Ntarget = Ntarget + tracks->Integral(firstbin, lastbin) / n;
            sumw = 0;
            mean = 0;
        }
    }
    mean = mean / sumw;
    std::cout << j << " " << Ntarget << " bin:" << lastbin << " ntracks:" << tracks->GetXaxis()->GetBinCenter(lastbin) << " mean:" << mean << std::endl;

}
void doNTracksPercentiles()
{
    TFile* f = TFile::Open("dataMB/AnalysisResults.root");
    std::string hName = "non-prompt-cascade-task/hNTracksVsCent";
    TH2* h2 = (TH2F*)f->Get(hName.c_str());
    auto h = h2->ProjectionY();
    if (h == nullptr) {
        std::cout << hName << " not found. Exiting." << std::endl;
        exit(1);
    }
    percentilesNTracks(h, 10, 30);
}
//
// Track multiplicity estiamtor
//
void analMultFind(TString normalisationFilename = "dataMB/AnalysisResults.root")
{
  doNTracksPercentiles();
  //return;
  //ROOT::EnableImplicitMT();
  //
  //std::vector<int> runs = {556160,556182,556741,556767,557415,557876,556152,556164,556734,557374,557425,557862,557897,557913,557926};
  std::vector<int> runs = {556160,556182,556741,556767,557415,557876,556152,556164,556734,557374,557425,557862,557897,557926};
  constexpr double CENT = 0.8;
  std::map<double,std::array<double,2>> centVsRun; 
  std::map<double,std::array<double,2>> centMinMaxVsRun; 
  std::map<double, std::array<double,2>> tracksVsRun;
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
  double cmin = 100;
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
    //TH1* hx = h->ProjectionX(namex.c_str(), biny, h->GetYaxis()->GetNbins());
    int maxBin = hx->GetMaximumBin();
    double centAtMax = hx->GetBinCenter(maxBin);
    double centMean = hx->GetMean();
    double centRMS = hx->GetRMS();
    double xmin = 0;
    double xmax = 100.;
    for(int i = 1; i < hx->GetNbinsX(); i++) {
      if(hx->GetBinContent(i) != 0.) {
        if(xmin == 0) {
          xmin = hx->GetBinCenter(i);
          break;
        }
      }
    }
    for(int i = hx->GetNbinsX(); i > 0; i--) {
      if(hx->GetBinContent(i) != 0.) {
        if(xmax == 100.) {
          xmax =  hx->GetBinCenter(i);
          break;
        }
      }
    }
    if(cmin > xmin) {
      cmin = xmin;
    }
    centMinMaxVsRun[runNumber][0] = xmin - centMean;
    centMinMaxVsRun[runNumber][1] = xmax - centMean;
    std::cout << runNumber << " cmean:" << centMean << " cmin:" << xmin << " cmax:" << xmax << std::endl; 
    centVsRun[runNumber][0] = centMean;
    centVsRun[runNumber][1] = centRMS;
    centPerRunHist[runNumber] = hx;
    //
    //
    //tracksPerRunHist[runNumber] = hTvsCZ->ProjectionY();
    int binx = hTvsCZ->GetXaxis()->FindBin(centAtMax);
    TH1* hy = hTvsCZ->ProjectionY("hy",binx,binx);
    tracksPerRunHist[runNumber] = hTvsCZ->ProjectionY();
    double meanTracks  = hy->GetMean();
    tracksVsRun[runNumber][0] = meanTracks;
    tracksVsRun[runNumber][1] = hy->GetMeanError();

    //std::cout << runNumber << " cent:" << centAtMax << " <Tracks>:" << meanTracks << std::endl;
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
  std::cout << "===> cmin:" << cmin << std::endl;
  //
  // Cent/Tracks versus Run
  //
  TGraphErrors* gCentVsRun = new TGraphErrors(centVsRun.size());
  TGraphAsymmErrors* gCentVsRunMinMax = new TGraphAsymmErrors(centVsRun.size());
  gCentVsRun->SetTitle("CentVsRun RMS Errors");
  gCentVsRunMinMax->SetName("CentVsRun MinMax errors");
  int i = 0;
  for(const auto&[key, val]: centVsRun) {
    gCentVsRun->SetPoint(i,key,val[0]);
    gCentVsRun->SetPointError(i,0,val[1]);
    gCentVsRunMinMax->SetPoint(i,key,val[0]);
    gCentVsRunMinMax->SetPointError(i, 0, 0, centMinMaxVsRun[key][0], centMinMaxVsRun[key][1]);
    i++;
  }
  //
  TGraphErrors* gTracksVsRun = new TGraphErrors(tracksVsRun.size());
  gTracksVsRun->SetTitle("TracksVsRun");
  gTracksVsRun->SetName("TracksVsRun");
  i = 0;
  for(const auto&[key, val]: tracksVsRun) {
    gTracksVsRun->SetPoint(i,key,val[0]);
    gTracksVsRun->SetPointError(i++,0,val[1]);
  }
  //
  if(0) {
    TGaxis::SetMaxDigits(6);
    gCentVsRun->SetMarkerStyle(20);
    gCentVsRun->SetMarkerSize(0.0);
    gCentVsRunMinMax->SetMarkerStyle(20);
    gCentVsRunMinMax->SetMarkerSize(0.7);
    gCentVsRunMinMax->SetMarkerColor(kBlue);
    gCentVsRunMinMax->SetLineColor(kBlue);
    gCentVsRunMinMax->SetTitle("Centrality vs Run");
    gCentVsRunMinMax->GetXaxis()->SetTitle("Run number");
    gCentVsRunMinMax->GetYaxis()->SetTitle("Centrality[%]");
    gCentVsRunMinMax->Draw("AP");
    gCentVsRun->Draw("P");
    TLatex tex;
    tex.SetNDC();               // use normalized (0–1) coordinates
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.50, 0.85, "Blue: Min - Max range");
    tex.DrawLatex(0.50, 0.80, "Black: RMS");
  }
  TGaxis::SetMaxDigits(6);
  gTracksVsRun->SetMarkerStyle(20);
  gTracksVsRun->SetMarkerSize(0.7);
  gTracksVsRun->GetHistogram()->SetMinimum(17);
  gTracksVsRun->Draw("AP");

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
  std::vector<double> x(NCENTBINS), y(NCENTBINS),ex(NCENTBINS),ey(NCENTBINS);
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
