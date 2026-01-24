void comp()
{
    TFile fileTRG("outputTRG.root");
    TH1F* hTRG = (TH1F*) fileTRG.Get("Omegant/spectrumnt")->Clone("hTRG");
    hTRG->SetDirectory(nullptr); 
    hTRG->SetName("trg");
    hTRG->SetTitle("trg");
    TFile fileMB("outputMB.root");
    TH1F* hMB = (TH1F*) fileMB.Get("cent0/Omegant/spectrumnt")->Clone("hMB");
    hMB->SetDirectory(nullptr); 
    hMB->SetName("mb");
    hMB->SetTitle("mb");

    TFile fileOut("outputComp.root","recreate");
    TH1F* hRatio = (TH1F*) hTRG->Clone("Ratio");
    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle("Ratio");
    hRatio->Divide(hMB);
    for(int i  = 0; i < hRatio->GetNbinsX(); i ++){
        std::cout << i << " " << hRatio->GetBinContent(i) << " " << hRatio->GetBinError(i) << std::endl;
    }
    //
    hRatio->GetYaxis()->SetTitle("");
    float ratio = 1600.0 / 2560.0;
    int W = 1000;
    int H = W * ratio;  
    TCanvas* cComp = new TCanvas("trgovercebt","TRG vs Cent",W,H);
    cComp->Divide(2,1);
    cComp->cd(1);
    //gPad->SetLeftMargin(0.15);
    hMB->SetLineColor(kBlue);
    hMB->SetLineWidth(2);
    hMB->Draw();
    hTRG->SetLineWidth(2);
    hTRG->Draw("same");
    cComp->cd(2);
    //gPad->SetLeftMargin(0.15);
    hRatio->SetLineWidth(2);
    hRatio->Draw();
    cComp->Draw();
    //
    hTRG->Write();
    hMB->Write();
    hRatio->Write();
    cComp->Write();
}
