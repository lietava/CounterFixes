const Double_t cutdummy= -99999.0;
double sigma_pt(double pt) {
  double a = 0.01;   // low-pt term
  double b = 0.01;   // constant
  double c = 0.001;  // high-pt term
  return pt * sqrt((a/pt)*(a/pt) + b*b + (c*pt)*(c*pt));
}
Double_t smear (Double_t pt){
  Double_t xeff= sqrt(pt)/sqrt(10.)/3.;  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  double sigma = sigma_pt(pt);
  Double_t xsmear= gRandom->Gaus(0.,sigma);     // bias and smear
  return pt+xsmear;
}
int generateMultiplicity() {
  double mean = 30;   // <N>
  double k    = 2.0;  // controls width (smaller = broader)
  return gRandom->NegativeBinomial(k, mean/(mean + k));
}
void test()
{
    auto *f0 = new TH1F("f0","f0",100,0,10);
    auto *g0 = new TH1F("g0","g0",100,-0,10);
    RooUnfoldResponse response (40, -10.0, 10.0);
    //
    TF1* fpt = new TF1("fpt", "x*exp(-x/[0])", 0, 10);
    fpt->SetParameter(0, 0.5); // "temperature" ~ 0.5 GeV
    //
    for (Int_t i= 0; i<100000; i++) {
        Double_t pt= fpt->GetRandom();
        f0->Fill(pt);
        Double_t ptm= smear (pt);
        if (ptm!=cutdummy){
        g0->Fill(ptm);
        response.Fill (ptm, pt);
    }
    else{
        response.Miss (pt);
    }
    }
    auto* c = new TCanvas();
    f0->SetStats(0);
    f0->SetTitle("");
    f0->SetFillColor(7);
    f0->Draw();
    g0->SetFillColor(42);
    g0->Draw("same");
    auto* leg = new TLegend(.55,0.7,.9,.9);
    leg->AddEntry(f0,"True Distribution");
    leg->AddEntry(g0,"Predicted Measured");
    leg->Draw();
    c->Draw();
    c->SaveAs("true-response.png");
    //return;
    //
    auto* R = response.HresponseNoOverflow();
    auto* c1 = new TCanvas();
    R->SetStats(0);
    R->Draw("colz");
    c1->Draw();
    c1->SaveAs("response.png");

}