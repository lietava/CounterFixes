void comp()
{
    TFile fileTRG("outputTRG.root");
    TH1F* hTRG = (TH1F*) fileTRG.Get("Omegant/spectrumnt");
    hTRG->SetName("trg");
    hTRG->SetTitle("trg");
    TFile fileMB("outputMB.root");
    TH1F* hMB = (TH1F*) fileMB.Get("cent0/Omegant/spectrumnt");
    hMB->SetName("mb");
    hMB->SetTitle("mb");

    TFile fileOut("outputComp.root","recreate");
    TH1F* hRatio = (TH1F*) hTRG->Clone("Ratio TRG/MB");
    hRatio->Divide(hMB);
    std::cout << hRatio->GetBinContent(2) << " " << hRatio->GetBinError(2) << std::endl;
    hTRG->Write();
    hMB->Write();
    hRatio->Write();
}
