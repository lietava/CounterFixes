double sigma_pt(double pt) {
  return pt * sqrt((0.02/pt)*(0.02/pt) + 0.01*0.01);
}
int smearMultiplicity(TRandom3& rand, int nTrue, double eff = 0.8, double fakeMean = 0.5) {
  int nRec = 0;
  for (int i = 0; i < nTrue; ++i) {
    if (rand.Rndm() < eff) {
      ++nRec;
    }
  }
  nRec += rand.Poisson(fakeMean);
  return nRec;
}
void test2()
{
  TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", 0, 10);
  fpt->SetParameter(0, 0.5);
  double eff = 0.8;
  double fakeMean = 0.;
  //
  double mu = 50.;
  int r = 20;
  double p = r/(r + mu);
  std::mt19937 gen(12345);  // RNG
  std::negative_binomial_distribution<int> nbd(r, p);
  TRandom3 rand(0);
  //
  TH1F* hPt = new TH1F("hPt", "p_{T}", 100, 0, 10);
  TH1F* hPtReco = new TH1F("hPtReco", "Reco p_{T}", 100, 0, 10);
  TH1F* hPtRes = new TH1F("hPtRes","hPtRes",1000,-1,1);
  TH1F* hMult = new TH1F("hMult", "Multiplicity Gen", 100, 0, 100);
  TH1F* hMultRec = new TH1F("hMultRec", "Multiplicity Rec", 100, 0, 100);


  for (int ev = 0; ev < 10000; ++ev) {

    // multiplicity
    int mult = nbd(gen);  // number of particles in coll 
    hMult->Fill(mult);
    // generate tracks
    int multRec  = 0;
    for (int i = 0; i < mult; ++i) {
      double pt = fpt->GetRandom();
      hPt->Fill(pt);
      if(rand.Rndm() < eff ) {
        double ptRec = rand.Gaus(pt, sigma_pt(pt));
        hPtReco->Fill(ptRec);  // reco
        hPtRes->Fill(pt-ptRec);
        multRec++;
      }
    } 
    // fake reconstructed tracks: contribute only to reco multiplicity
    int nFake = rand.Poisson(fakeMean);
    multRec += nFake;

    // optional: give fake tracks a reco pT spectrum too
    for (int i = 0; i < nFake; ++i) {
      double ptFake = fpt->GetRandom();   // toy choice
      hPtReco->Fill(ptFake);
    }
    hMultRec->Fill(multRec);
  }
  TFile* fout = new TFile("output.root","RECREATE");
  hPt->Write();
  hPtReco->Write();
  hMult->Write();
  hPtRes->Write();
  hMultRec->Write();
  hPtRes->Draw();
}