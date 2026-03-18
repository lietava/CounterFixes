void toy()
{
    int Nbins = 30;
    TH1D hMult4pi("mult","mult",200,0.,200.);
    TH1D hMultF("multF","multF",Nbins,0.,Nbins);
    TH1D hMultC("multC","multC",Nbins,0.,Nbins);
    TH2D hCent("CentVsNTrk","CentVsNTrk",10000,0.,100.,Nbins,0,Nbins);
    TH2D hMult("MultVsNTrk","Forward Trk Vs Cent Trk",5000,0.,5000.,300,0,300.); // mult =  ft0
    TH1D hRap("hRap","hRap", 20,-10,10);

    std::map<int,array<int,3>> events;
    //std::map<int,int> multF;
    std::map<int,int> multC;
    std::vector<std::pair<int,int>> multF;


    int N = 100000000; // number of collisions
    int sigma = 3;
     TRandom3 rng(0);
    double mu = 50; // avarega mult in pp
    TF1 fEta("fEta","1/(1+exp((abs(x)-[0])/[1]))",-10,10);
    fEta.SetParameters(2.5,0.6); // plateau width and edge steepness
    //
    int r = 20;
    double p = r/(r + mu);
    std::mt19937 gen(12345);  // RNG
    std::negative_binomial_distribution<int> nbd(r, p);

    TStopwatch sw;
    sw.Start();
   
    int ncol = 0;
    for (int i = 0; i < N; ++i) {
        int np = nbd(gen);  // number of particles in coll 
        hMult4pi.Fill(np);
        int nForw = 0;
        int nCent = 0;
        //std::cout << "===> " << np << std::endl;
        for(int j = 0; j < np; j++) {
            //double y = rng.Uniform(-10.0, 10.0);
            //double y = rng.Gaus(0.,sigma);
            double y = fEta.GetRandom();
            hRap.Fill(y);
            if((y > 3.5)&& (y < 4.9)) {
                nForw += 50;
            }
            if((y > -3.3)&& (y < - 2.1)) {
                nForw += 50;
            }
            if(abs(y) < 0.8) {
                if(rng.Uniform(0,1) < 0.5)
                {
                    nCent += 1;
                }
            }
            //std::cout << nForw << " " << nCent << std::endl;
        }
        hMultF.Fill(nForw);
        hMult.Fill(nForw, nCent);
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
    sw.Stop();

    TProfile* prfCent = hCent.ProfileX("Mean track vs CentFT0M");
    TProfile* prfMult = hMult.ProfileX("Mean track central vs mean track forward");
    prfMult->SetMarkerSize(0.7);
    prfMult->SetMarkerStyle(20);
    prfMult->Draw("P");
    std::cout << "Real time:" << sw.RealTime() << " CPU:" << sw.CpuTime() << std::endl;
    //
    //
    TFile output("outputToy.root","RECREATE");

    hMult4pi.Write();
    hMultC.Write();
    hMultF.Write();
    hCent.Write();
    hMult.Write();
    prfCent->Write();
    prfMult->Write();
    hRap.Scale(1./N);
    hRap.Write();
    // TCanvas can("can","can",1000,400);
    // can.Divide(2,2);
    // can.cd(1);
    // hMult4pi.Draw();
    // can.Draw();   
}
