{
  TFile *f=new TFile("Histos.root");
  TH1F* h1 = (TH1F*)f->Get("Skimming efficiency");
}
