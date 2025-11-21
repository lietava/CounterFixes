int utils( TString normalisationFilename = "results/AnalysisResults.root", int Z = 1, std::string comment = "", int ntoisman = 1)
{
    std::cout << comment << std::endl;
    std::cout << "File:" << normalisationFilename.Data() << std::endl;
    TFile normalisationFile(normalisationFilename.Data());
    //
    ZorroSummary *zorroSummary;
    if(Z ==1 ) {
        zorroSummary = (ZorroSummary *)normalisationFile.Get("non-prompt-cascade-task/ZorroSummary");
    } else {
        zorroSummary = (ZorroSummary *)normalisationFile.Get("non-prompt-cascade-task/zorroSummary");
    }
    if(zorroSummary == nullptr) {
        std::cout << "Can not get zorro" << std::endl;
        exit(1);
    } else {
        std::string text = zorroSummary->getTOInames();
        std::cout << "TOI:" <<  text << std::endl;
        int count;
        if(text == "") {
            count = ntoisman;
        } else {
            count = std::count(text.begin(), text.end(), ',') + 1;
        }
        std::cout << "Number of TOIs:" << count << std::endl;
        for(int i = 0; i < count; i++) {
            double nf = zorroSummary->getNormalisationFactor(i);
            std::cout << i << " Norm factor " << nf <<  comment << std::endl;
        }
    }
    return 0;
}
