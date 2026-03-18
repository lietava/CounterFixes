//void thr(const char* filename = "AnalysisResults_567576_config2026.root")
void thr(const char* filename = "AnalysisResults_616860train.root")

{
    //const char* histname = "lf-strangeness-filter/EventsvsMultiplicity/AllEventsvsMultiplicityFT0MNorm";
    //const char* histnameO = "lf-strangeness-filter/EventsvsMultiplicity/AllEventsvsMultiplicityFT0MwOmegaNorm";
    const char* histname = "lf-strangeness-filter/EventsvsMultiplicity/AllEventsvsMultiplicityTracksGlob";
    const char* histnameO = "lf-strangeness-filter/EventsvsMultiplicity/AllEventsvsMultiplicityTracksGlobwOmega";


  // Open file
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
      printf("❌ Cannot open file: %s\n", filename);
      return;
  }

  // Get histogram
  TH1* h = nullptr;
  f->GetObject(histname, h);
  std::cout << h << std::endl;
  if (!h) {
      printf("❌ Histogram %s not found in %s\n", histname, filename);
      return;
  }
  // Get histogram
  TH1* hO = nullptr;
  f->GetObject(histnameO, hO);
  std::cout << hO << std::endl;
  if (!hO) {
      printf("❌ Histogram %s not found in %s\n", histnameO, filename);
      return;
  }
  // --- 1) NORMALISE TO TOTAL PROBABILITY = 1 ---
  double integral = h->Integral();
  if (integral == 0) {
    printf("⚠️ Histogram integral is zero, cannot normalize.\n");
  } else {
    h->Scale(1.0 / integral);
  }
  // --- 2) BUILD CUMULATIVE PROBABILITY (CDF) ---
  int nbins = h->GetNbinsX();
  TH1* hCdf = (TH1*)h->Clone("hCdf");
  hCdf->Reset("ICES");
  hCdf->SetTitle("CDF (cumulative probability);x;F(x)");

  double cum = 0.0;
  double cum2 = 0;
  for (int i = 1; i <= nbins; ++i) {
    cum += h->GetBinContent(nbins-i+1);
    cum2 += hO->GetBinContent(i);
    hCdf->SetBinContent(nbins-i+1, cum);
    std::cout << hO->GetBinCenter(i) << " " << hO->GetEntries() - cum2 << std::endl; 
  }
  // After this, if normalization worked, last bin should be ~1

  // --- 3) DRAW BOTH ---
  if(1) {
     // Optional cosmetics
    h->SetTitle("PDF (normalised histogram);x;P(x)");
    h->SetLineWidth(2);
    TCanvas* c = new TCanvas("c", "PDF and CDF", 500, 1000);
    c->Divide(1,2);

    c->cd(1);
    h->Draw("HIST");

    c->cd(2);
    hCdf->SetLineWidth(2);
    hCdf->Draw("HIST");

    printf("✅ Normalized histogram and CDF created.\n");
    printf("   PDF hist name : %s\n", h->GetName());
    printf("   CDF hist name : %s\n", hCdf->GetName());
    int nbinsCdf = hCdf->GetNbinsX();

    printf("Bin   x_center          content\n");
    printf("-------------------------------------\n");

    for (int i = 1; i <= nbinsCdf; ++i) {
        double x = hCdf->GetBinCenter(i);
        double y = hCdf->GetBinContent(i);
        //printf("%4d  %12.6f   %12.6f\n", i, x, y);
        std::cout << x << " " << y << std::endl;
    }
  }
  if(0) {
     // Optional cosmetics
    //h->SetTitle("PDF (normalised histogram);x;P(x)");
    h->SetLineWidth(2);
    TCanvas* c = new TCanvas("c", "MB and HM Omega", 500, 1000);
    c->Divide(1,2);

    c->cd(1);
    hO->Draw("HIST");

    c->cd(2);
    h->SetLineWidth(2);
    h->Draw("HIST");
  }
}