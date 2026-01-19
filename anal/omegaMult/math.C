void math()
{
    int Nbins = 30;
    TH1D hMult4pi("mult","mult",200,0.,200.);
    TH1D hMultF("multF","multF",Nbins,0.,Nbins);
    TH1D hMultC("multC","multC",Nbins,0.,Nbins);
    TH2D hCent("CentVsNTrk","CentVsNTrk",10000,0.,100.,Nbins,0,Nbins);
    std::map<int,array<int,3>> events;
    //std::map<int,int> multF;
    std::map<int,int> multC;
    std::vector<std::pair<int,int>> multF;


    int N = 10000000;
     TRandom3 rng(0);
    double mu = 100;
    //
    int r = 20;
    double p = r/(r + mu);
    std::mt19937 gen(12345);  // RNG
    std::negative_binomial_distribution<int> nbd(r, p);

    TStopwatch sw;
    sw.Start();
   
    for (int i = 0; i < N; ++i) {
        int np = nbd(gen);   
        hMult4pi.Fill(np);
        int nForw = 0;
        int nCent = 0;
        //std::cout << "===> " << np << std::endl;
        for(int j = 0; j < np; j++) {
            double y = rng.Uniform(-10.0, 10.0);
            if(abs(y) > 9.) {
                nForw++;
            }
            if(abs(y) < 1.0) {
                nCent++;
            }
            //std::cout << nForw << " " << nCent << std::endl;
        }
        if(nForw != 0.) {
            hMultF.Fill(nForw);
        }
        hMultC.Fill(nCent);
        events[i][0] = np;
        events[i][1] = nCent;
        events[i][2] = nForw;
        //multF[i] = nForw;
        multF.push_back({i,nForw});
        multC[i] = nCent;
    }
    //
    //
    std::sort(multF.begin(), multF.end(),
    [](const auto &a, const auto &b) {
        return a.second < b.second;
    });
    double i = 0.;
    for(const auto& m: multF) {
        //std::cout << i * 100./(double) N << " " <<  multC[m.first] << " " << m.first << " " << m.second << std::endl;
        hCent.Fill(100. - i * 100./(double) N, multC[m.first]);
        i+=1.;
    }
    TProfile* prf = hCent.ProfileX("Mean track vs CentFT0M");
    TGraphErrors* gp = new TGraphErrors(prf);
    TGraph *g = new TGraph();
    int j = 0;
    for (int i = 0; i < gp->GetN(); i++) {
        double x, y;
        gp->GetPoint(i, x, y);
        if (y == 0) continue;   // skip zero points 
        g->SetPoint(j, x, y);
        j++;
    }
    g->GetXaxis()->SetTitle("Centrality[%]");
    g->GetYaxis()->SetTitle("<dN/dy>");
    sw.Stop();
    std::cout << "Real time:" << sw.RealTime() << " CPU:" << sw.CpuTime() << std::endl;
    g->SetMarkerSize(0.7);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    //
    //
    TFile output("outputMath.root","RECREATE");
    hMult4pi.Write();
    hMultC.Write();
    hMultF.Write();
    hCent.Write();
    // TCanvas can("can","can",1000,400);
    // can.Divide(2,2);
    // can.cd(1);
    // hMult4pi.Draw();
    // can.Draw();   
}
