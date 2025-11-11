void theplot()
{
    TFile* f = TFile::Open("HEPData-ins1471838-v1-Table_45.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file!" << std::endl;
        return;
    }
    TH1F* hrun2 = (TH1F*)f->Get("Table 45/Hist1D_y1");
    if (!hrun2) {
        std::cerr << "Can not get published omega" << std::endl;
        return;
    }
    TGraphErrors* g = (TGraphErrors*)f->Get("Table 45/Graph1D_y1");
    if (!g) {
        std::cerr << "Graph not found!" << std::endl;
        return;
    }

    int Npoints = g->GetN();
    std::vector<float> xrun2(Npoints), yrun2(Npoints), exrun2(Npoints), eyrun2(Npoints);
    std::cout << "Graph has " << Npoints << " points." << std::endl;
    for (int i = 0; i < Npoints; i++) {
        double x, y, ex, ey;
        g->GetPoint(i, x, y);
        xrun2[i] = x;
        yrun2[i] = y;
        exrun2[i] = g->GetErrorX(i);
        eyrun2[i] = g->GetErrorY(i);
        std::cout << x << " " << y << " " << exrun2[i] << " " << eyrun2[i] << std::endl;
    }
    TFile* ftracks = TFile::Open("outputMult.root");
    if (!ftracks || ftracks->IsZombie()) {
        std::cerr << "Error: cannot open file!" << std::endl;
        return;
    }
    TGraphErrors* gdatatracks = (TGraphErrors*)ftracks->Get("NTracks");
    if (!gdatatracks) {
        std::cerr << "Graph not found!" << std::endl;
        return;
    }
    int Npoints3 = gdatatracks->GetN();
    std::vector<float> xrun3(Npoints3), yrun3(Npoints3);
    for (int i = 0; i < gdatatracks->GetN(); i++) {
        double x, y, ex, ey;
        gdatatracks->GetPoint(i, x, y);
        xrun3[i] = x;
        yrun3[i] = y;
        //ey = g->GetErrorY(i);
        //std::cout << x << " " << y << " " << ex << " " << ey << std::endl;
    }
    TFile* fomegas = TFile::Open("outputMB.root");
    if (!fomegas || fomegas->IsZombie()) {
        std::cerr << "Error: cannot open file!" << std::endl;
        return;
    }
    int NCENTBINS = 11;
    std::vector<float> yields3(NCENTBINS), eyields(NCENTBINS), ex(NCENTBINS);
    for (int i = 0; i < NCENTBINS; i++) {
        std::string hname = "cent" + std::to_string(i) + "/Omegant/spectrumnt";
        //std::cout << hname << std::endl;
        TH1F* h = (TH1F*)fomegas->Get(hname.c_str());
        if (!h) {
            std::cout << "Can not get:" << hname << std::endl;
            exit(1);
        }
        //std::cout << i << std::endl;
        //yields3[i] = 1. * h->Integral();
        double err = 0;
        yields3[i] = h->IntegralAndError(1, h->GetNbinsX(), err);
        eyields[i] = err;
        ex[i] = 0;
        std::cout << yrun3[i] << ":" << yields3[i] << " e:" << eyields[i] << std::endl;
    }
    TGraphErrors* grun3 = new TGraphErrors(NCENTBINS - 1, yrun3.data(), yields3.data(), ex.data(), eyields.data());
    //
    // NTracks
    // 
    //
    TFile* ftracksfile = TFile::Open("outputNTracks.root");
    if (!ftracksfile || ftracksfile->IsZombie()) {
        std::cerr << "Error: cannot open file!" << std::endl;
        return;
    }
    TH1F* omvstrk = (TH1F*) ftracksfile->Get("omVsNTracks");
    if (!omvstrk) {
        std::cout << "Can not get omvstrk" << std::endl;
        exit(1);
    }
    //
    TFile* fout = new TFile("outputCompRun2.root", "recreate");
    grun3->Write();
    //
    float dx[2] = {0, 40};
    float dy[2] = {0, 0.01};
    TGraph* dummy = new TGraph(2, dx, dy);
    dummy->Draw("AP");
    g->SetMarkerSize(0.7);
    g->SetMarkerStyle(21);
    g->Draw("P");
    grun3->SetMarkerSize(0.7);
    grun3->SetMarkerStyle(20);
    grun3->SetMarkerColor(kBlue);
    grun3->Draw("P");
    omvstrk->SetMarkerSize(1.2);
    omvstrk->SetMarkerStyle(22);
    omvstrk->SetMarkerColor(kRed);
    omvstrk->Draw("P");
    dummy->SetTitle("Omega yield vs NTracksGlobal");
}; 