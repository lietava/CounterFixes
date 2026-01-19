void poisson()
{
    const int n = 7;
    float x[] =  {1,2,3,4,5,6,7};
    float ex[] = {0,0,0,0,0,0,0};
    float y[] = {725868, 11360, 748, 276, 109, 24, 7};
    // const int n = 3;
    // float x[] =  {1,2,3};
    // float ex[] = {0,0,0};
    // float y[] = {704338, 3253, 12};
    float ey[n];
    for(int i = 0; i < n; i++) {
        ey[i] = sqrt(y[i]);
        std::cout << x[i] << " " << y[i] << " " << ey[i] << std::endl;
    }
    TGraphErrors* g = new TGraphErrors(n, x, y, ex, ey);
    //g->SetTitle("Coll ditribution with kNoSameBunchPileup");
    g->SetTitle("Coll ditribution run 551760");

    g->SetMarkerSize(1);
    g->SetMarkerStyle(20);

    float xmin = 0.9;
    float xmax = n + 0.1;
    TF1* fPois = new TF1("fPois", "[0]*TMath::Poisson(x,[1])", xmin, xmax);
    fPois->SetParNames("Norm", "Mean");
    fPois->SetParameters(100., 3.); // initial guesses: N≈total counts, μ≈where it peaks
    g->Fit("fPois", "R");   // "R" = respect the [xmin,xmax] range of the function
    g->GetXaxis()->SetTitle("numContrib");
    g->Draw("AP"); 

}
