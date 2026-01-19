void mult()
{
    TFile fileMB("dataMB/AnalysisResults.root");
    //TFile fileMB("../omegaProdDerived/AnalysisResults.root");

    TH2D* hTrkVsCnt = (TH2D*) fileMB.Get("non-prompt-cascade-task/hNTracksVsCentZoom");
    if(hTrkVsCnt == nullptr) {
        std::cout << "Can not get hist" << std::endl;
        exit(1);
    }
    TFile fileOut("outputmult.root","recreate");
    TH1* hyp = hTrkVsCnt->ProjectionX("hyp");
    hyp->Write();
    TProfile* p = hTrkVsCnt->ProfileX("Mean track vs CentFT0M");
    TGraphErrors* gp = new TGraphErrors(p);
    TGraph *g = new TGraph();

    int j = 0;
    for (int i = 0; i < gp->GetN(); i++) {
        double x, y;
        gp->GetPoint(i, x, y);
        if (y == 0) continue;   // skip zero points 
        g->SetPoint(j, x, y);
        j++;
    }
    g->SetMarkerSize(0.7);
    g->SetMarkerStyle(20);
    //TF1 *f = new TF1("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0.1, 1);
    //TF1 *f = new TF1("f", "[0]*exp(-x*x*[1])", 0., 2);
    //f->SetParameters(21,-2,0.2,0);
    //    TF1 *f = new TF1("f", "[0]*exp(-x*x*[1])", 0., 2);
    //g ->Fit(f,"V");
    g->GetXaxis()->SetTitle("Centrality[%]");
    g->GetYaxis()->SetTitle("<dN/dy>");
    g->Draw("AP");
    p->Write();
    int b1 = hTrkVsCnt->GetXaxis()->FindBin(0.1);
    TH1* hY_1 = hTrkVsCnt->ProjectionY("hCent_1",1,b1);
    std::cout << "0.1 bin:" << b1 << ":" << hY_1->GetMean() << std::endl;

    hY_1->Write();
    for(int i = 1; i < hTrkVsCnt->GetNbinsX(); i ++) {
        std::string cnt = "cent_"+to_string(i);
        TH1* hY = hTrkVsCnt->ProjectionY(cnt.data(),i,i);
        if (hY->GetEntries()) {
            //std:;cout << hTrkVsCnt->GetXaxis()->GetBinCenter(i) << " " << hY->GetMean(1) << std::endl;
            hY->Write();
        }
    }
}
