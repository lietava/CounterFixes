int plot()
{
  TFile *f = new TFile("eff.root");
  TH1F* hfU = (TH1F*) f->Get("Freq Unskimm");
  TH1F* hfS = (TH1F*) f->Get("Freq Skimmed");
  TH1F* hEffSkimmed = (TH1F*) f->Get("Skimming eff");

 //Double_t scale = hfU->GetXaxis()->GetBinWidth(1)/(hfU->Integral());
 //hfU->Draw("HIST");
 //hfS->Draw("HIST");
 hEffSkimmed->Draw("HIST");
 return 0;
}
