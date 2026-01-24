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
TH1* doNTracksPercentiles()
{
    TFile* f = TFile::Open("dataMB/AnalysisResults.root");
    std::string hName = "non-prompt-cascade-task/hNTracksVsCent";
    TH2* h2 = (TH2F*)f->Get(hName.c_str());
    auto h = h2->ProjectionY();
    if (h == nullptr) {
        std::cout << hName << " not found. Exiting." << std::endl;
        exit(1);
    }
    return h;
}
void analNTrackFind(TString normalisationFilename = "dataMB/AnalysisResults.root")
{
  //auto h =  doNTracksPercentiles();
  //percentilesNTracks(h, 5, 0);
  //return;
  //ROOT::EnableImplicitMT();
  //
  // Track mults for cent bins
  constexpr int NCENTBINS = 5; // - 1
  // NTrack bins
  double centbins[NCENTBINS + 1]{0, 2,4,7,12,100};
  //double centbins[NCENTBINS + 1]{ 0, 10,11,12,13,14, 16, 18, 20,25,50,100 };
  //
  TFile *f = TFile::Open(normalisationFilename);
  std::string dir = "non-prompt-cascade-task/";  // your directory name
  std::string histnameTracks = dir + "hNTracksVsCent";
  TH2* h2tracks = (TH2F*)f->Get(histnameTracks.c_str());
  if (h2tracks == nullptr) {
      std::cout << histnameTracks << " not found. Exiting." << std::endl;
      exit(1);
  }
  //
  std::array<TH1*, NCENTBINS> binsTracks;
  std::array<float, NCENTBINS> norm;
  std::vector<float> xntracks(NCENTBINS);
  auto ntracks = h2tracks->ProjectionY();
  float ntot = ntracks->GetEntries();
  std::cout << "Total:" << ntot << std::endl;
  for (int ic = 0; ic < NCENTBINS; ic++) {
      int biny1 = ntracks->GetXaxis()->FindBin(centbins[ic]);
      int biny2 = ntracks->GetXaxis()->FindBin(centbins[ic + 1]);
      std::string name = "tracksNTracks_" + to_string(ic);
      //std::cout << biny1 << " " << biny2 << std::endl;
      binsTracks[ic] = new TH1F(name.data(), name.data(), biny2 - biny1, biny1, biny2);
      for (int i = biny1; i < biny2; i++) {
          //std::cout << i << " " << ntracks->GetBinContent(i) << std::endl;
          binsTracks[ic]->SetBinContent(i - biny1 + 1, ntracks->GetBinContent(i));
          binsTracks[ic]->SetBinError(i -  biny1 + 1, ntracks->GetBinError(i));
      }
      std::cout << "Mean:" << binsTracks[ic]->GetMean() << " Inte:" << binsTracks[ic]->Integral()/ntot << std::endl;
      norm[ic] = binsTracks[ic]->Integral() / ntot;
      xntracks[ic] = binsTracks[ic]->GetMean();
  }
  TFile *finput = TFile::Open("outputMBNtracks.root","RECREATE");
  std::vector<float> yields3(NCENTBINS), eyields(NCENTBINS), ex(NCENTBINS);
  for (int i = 0; i < NCENTBINS; i++) {
      std::string hname = "cent" + std::to_string(i) + "/Omegant/spectrumnt";
      //std::cout << hname << std::endl;
      TH1F* h = (TH1F*)finput->Get(hname.c_str());
      if (!h) {
          std::cout << "Can not get:" << hname << std::endl;
          exit(1);
      }
      //std::cout << i << std::endl;
      //yields3[i] = 1. * h->Integral();
      double err = 0;
      //std::cout << h->IntegralAndError(1, h->GetNbinsX(), err) << std::endl;
      yields3[i] = h->IntegralAndError(1, h->GetNbinsX(), err) / norm[i];
      eyields[i] = err;
      ex[i] = 0;
      std::cout << xntracks[i] << ":" << yields3[i] << " e:" << eyields[i] << std::endl;
  }
  yields3[NCENTBINS - 1] = 0.;
  TGraphErrors* grun3 = new TGraphErrors(NCENTBINS - 1, xntracks.data(), yields3.data(), ex.data(), eyields.data());
  grun3->SetName("omVsNTracks");
  grun3->SetTitle("Omegas vs NTracks");
  //
  // Output
  //
  TFile output("outputNTracks.root", "recreate");
  ntracks->Write();
  grun3->Write();
}
