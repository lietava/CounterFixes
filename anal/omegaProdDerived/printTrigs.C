void printTrigs()
{
    TFile file("results/AnalysisResults.root");
    TH1* h = (TH1*)file.Get("non-prompt-cascade-task/Zorro/556767/Selections");
    if(h == nullptr) {
        std::cout <<  "Histo not found" << std::endl;
        exit(1);
    }  
    TAxis* ax = h->GetXaxis();
    int nb = ax->GetNbins();

    for (int i = 1; i <= nb; i++) {
        std::cout << "Bin " << i 
              << " label = " << ax->GetBinLabel(i) 
              << std::endl;
    }  
}

