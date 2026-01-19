std::vector<std::string> names = {"fOmegaHighMult","fLeadingPtTrack","fHighFt0cFv0Flat","fHighFt0cFv0Mult","fHighFt0Flat","fHighFt0Mult","fHighMultFv0","fHighTrackMult"};
int getNormFromFile( TString normalisationFilename = "AnalysisResults.root", int Z = 1, std::string comment = "", int ntoisman = 1)
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
            std::cout << i << " " << names[i] << " Norm factor " << nf <<  comment << std::endl;
        }
    }
    return 0;
}
int getNorm()
{
    getNormFromFile("AnalysisResults.root");   
    return 0;
    if(1) {
        getNormFromFile("../omegaProd0104/AnalysisResults.root");
        getNormFromFile("../omegaProd0605/AnalysisResults.root");
        getNormFromFile("../omegaProd0506/AnalysisResults.root");
        getNormFromFile("../omegaProd0406/AnalysisResults.root");
        getNormFromFile("../omegaProd1505/AnalysisResults.root");
        getNormFromFile("../omegaProd0805/AnalysisResults0.root");
        getNormFromFile("../omegaProd0805/AnalysisResults1.root");
        getNormFromFile("../omegaProd0805/AnalysisResults.root");

    } else if(0) {
        // merged ((run1) + (run2)) versus (run1+run2)
        std::cout << "merged ((run1) + (run2)) versus (run1+run2)" << std::endl;
        getNormFromFile( "results550889/AnalysisResults.root");
        getNormFromFile( "results556152/AnalysisResults.root");
        getNormFromFile( "resultsBoth2/AnalysisResults.root");
        getNormFromFile( "resultsMerged2runs/AnalysisResults.root"); 
    } else if(0) {
        // (toi1) + (toi2) versus (toi1+toi2)
        std::cout << "(toi1) + (toi2) versus (toi1+toi2) ; 2 runs in data" << std::endl;
        getNormFromFile( "resultsBoth2/AnalysisResults.root");
        getNormFromFile( "resultsOmega2/AnalysisResults.root");
        getNormFromFile( "resultsTrk2/AnalysisResults.root");
        getNormFromFile( "resultsMerged2tois/AnalysisResults.root");
    } else {
        // file produces with old tag hrun 474352
        getNormFromFile("/Users/rl/anal/TEST/omegaAnalMax/dataTRG/474352/AnalysisResults.root", 0, "===> H run: 474352 old tag, run recently");
        // file with latest code
        getNormFromFile("/Users/rl/anal/TEST/omegaHMTrig/dataTRG/455519/AnalysisResults2107.root", 1, "===> H run: 455519 tag with extended tables",2);
        // old tag old file
        getNormFromFile("/Users/rl/anal/TEST/omegaAnalMax/dataTRG/348293/AnalysisResults.root", 0, "===> H run: 348293 old tag,old data");
        //
        getNormFromFile("/Users/rl/anal/TEST/omegaAnalMax/dataTRG/472810/AnalysisResults.root", 0, "===> H run: 472810 old tag,old data");
        // latest
        getNormFromFile("/Users/rl/anal/TEST/omegaHMTrig/dataTRG/478154/AnalysisResults.root", 1, "===> H run: 478154 last tag with extended tables, only HM",2);
        // latest
        getNormFromFile("/Users/rl/anal/TEST/omegaHMTrig/dataTRG/479741/AnalysisResults.root", 1, "===> H run: 479741 last tag with extended tables",2);
    }

    return 0;
}
