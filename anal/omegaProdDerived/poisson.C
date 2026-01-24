void poisson()
{
    float mu = 0.0246;
    //
    const int nf = 7;
    float xf[] =  {1,2,3,4,5,6,7};
    float exf[] = {0,0,0,0,0,0,0};
    float yf[] = {725868, 11360, 748, 276, 109, 24, 7};
    //
    const int n = 2;
    float x[] =  {1,2};
    float ex[] = {0,0};
    float y[] = {725868, 11360};
    //
    // const int n = 3;
    // float x[] =  {1,2,3};
    // float ex[] = {0,0,0};
    // float y[] = {704338, 3253, 12};
    float ey[n],eyf[nf];
    for(int i = 0; i < n; i++) {
        ey[i] = sqrt(y[i]);
        std::cout << x[i] << " " << y[i] << " " << ey[i] << std::endl;
    }
    for(int i = 0; i < nf; i++){
        eyf[i] = sqrt(yf[i]);
    }
    TGraphErrors* gf = new TGraphErrors(nf, xf, yf, exf, eyf);
    //g->SetTitle("Coll ditribution with kNoSameBunchPileup");
    gf->GetXaxis()->SetTitle("numContrib");
    gf->Draw("AP"); 
    gf->SetTitle("Coll ditribution run 551760");
    gf->SetMarkerSize(1);
    gf->SetMarkerStyle(20);
    //
    TGraphErrors* g = new TGraphErrors(n, x, y, ex, ey);
    float xmin = 0.9;
    float xmax = nf + 0.1;
    TF1* fPois = new TF1("fPois", "[0]*TMath::Poisson(x,[1])", xmin, xmax);
    fPois->SetParNames("Norm", "Mean");
    fPois->SetParameters(100., 3.); // initial guesses: N≈total counts, μ≈where it peaks
    g->Fit("fPois", "R");   // "R" = respect the [xmin,xmax] range of the function
    double norm = fPois->GetParameter(0);
    fPois->Draw("same");
    //
    TGraphErrors* g2 = new TGraphErrors(n, x, y, ex, ey);
    TF1* fPoisAlt = new TF1("fPoisAlt","[0]*TMath::Poisson(x,0.0246)", xmin, xmax);
    fPoisAlt->SetParameters(norm);  // Norm, mu
    g2->Fit("fPoisAlt", "R");   // "R" = respect the [xmin,xmax] range of the function
    fPoisAlt->SetLineColor(kBlue);
    // //
    fPoisAlt->Draw("SAME");
}
