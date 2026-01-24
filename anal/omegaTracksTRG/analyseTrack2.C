/// fBachKaonTPCNSigma>-3&&fBachKaonNClusTPC>80&&std::abs(fMassXi-1.32171)>0.008
constexpr double pt[8]{1., 1.4, 1.8, 2.3, 2.8, 3.3, 3.8, 4.8};

void fillChainFromAO2D(TChain &chain, const TString &fileName)
{
  TFile file(fileName);
  for (auto key : *file.GetListOfKeys())
  {
    TString keyName = key->GetName();
    if (keyName.Contains("DF_"))
    {
      chain.Add((fileName + "/" + keyName + "/" + chain.GetName()).Data());
    }
  }
}

void analyseTrack2(TString dataFilename = "dataTRG/513390/AO2D.root", TString mcFilename = "mcMB/AO2D.root", TString normalisationFilename = "dataTRG/513390/AnalysisResults.root")
{
  TString selectionString{"std::abs(fProtonEta) < 0.9&&std::abs(fPionEta)<0.9&&fCascCosPA > 0.9995&&std::abs(fBachEta)<0.9&&fBachNClusTPC>70&&std::abs(fMassXi-1.32171)>0.008&& std::abs(fPvZ)<10"};

  ROOT::EnableImplicitMT();
  std::cout << "Starting" << std::endl;
  //constexpr int NCENTBINS = 11;
  //double centbins[NCENTBINS + 1]{0,10,20,30,40,50,60,70,80,90,100,150};
      
  constexpr int NCENTBINS = 1;
  double centbins[NCENTBINS + 1] {60,100};
  
  //
  constexpr int NPTBINS = 7;
  double val[NPTBINS]{0.0009488760957, 0.0006952320908, 0.0004759373703, 0.0002954698723, 0.0001662988815, 8.73277117e-05, 3.635092251e-05};
  double stat[NPTBINS]{4.808611514e-05, 2.719983815e-05, 1.605100658e-05, 1.091321782e-05, 7.543390105e-06, 5.044361557e-06, 2.234297454e-06};
  double syst[NPTBINS]{5.649579873e-05, 2.975781163e-05, 2.84376873e-05, 9.258206825e-06, 3.867393999e-06, 4.169908233e-06, 2.43093823e-06};
  TH1D published_stat("published_stat", ";#it{p}_{T} (GeV/c);d^{2}#it{N}/d#it{p}_{T}d#it{y} (c/GeV)", NPTBINS, pt);
  TH1D published_syst("published_syst", ";#it{p}_{T} (GeV/c);d^{2}#it{N}/d#it{p}_{T}d#it{y} (c/GeV)", NPTBINS, pt);
  for (int i{1}; i <= NPTBINS; ++i)
  {
    published_stat.SetBinContent(i, val[i - 1]);
    published_stat.SetBinError(i, stat[i - 1]);
    published_syst.SetBinContent(i, val[i - 1]);
    published_syst.SetBinError(i, syst[i - 1]);
  }
  published_stat.SetMarkerStyle(20);
  published_stat.SetMarkerSize(0.8);
  published_stat.SetMarkerColor(kRed);
  published_stat.SetLineColor(kRed);
  published_syst.SetMarkerStyle(20);
  published_syst.SetMarkerSize(0.8);
  published_syst.SetMarkerColor(kRed);
  published_syst.SetLineColor(kRed);
  published_syst.SetFillStyle(0);
  std::cout << "init finished" << std::endl;
  //
  std::string names[2]{"", "nt"};
  TFile normalisationFile(normalisationFilename.Data());
  // ZorroSummary *zorroSummary = (ZorroSummary *)normalisationFile.Get("non-prompt-cascade-task/zorroSummary");
  TH1 *hCounterTVX = static_cast<TH1 *>(normalisationFile.Get("lumi-task/hCounterTVX"));
  if(hCounterTVX == nullptr) {
    hCounterTVX = static_cast<TH1 *>(normalisationFile.Get("eventselection-run3/luminosity/hCounterTVX"));
    if(hCounterTVX == nullptr) {
      std::cout << "TVX counter not available" << std::endl;
      exit(1);
    }
  }
  double normalisations[2]{0.7558 / hCounterTVX->GetEntries(), 0.7558 / hCounterTVX->GetEntries()};
  std::cout << "Normalisation finished:" << hCounterTVX->GetEntries() << std::endl;
  std::array<float, NCENTBINS > normCent;
  for(int i = 0; i < NCENTBINS; i++){
    normCent[i] = 1;
  }
  //
  TChain genChain("O2npcasctablegen");
  fillChainFromAO2D(genChain, mcFilename);
  ROOT::RDataFrame genDF(genChain);
  auto genDF2 = genDF.Define("pz", "fgPt*std::sinh(fgEta)").Define("eOmega", "std::sqrt(fgPt*fgPt + pz*pz + 1.67245 * 1.67245)").Define("yOmega", "0.5*std::log((eOmega + pz)/(eOmega - pz))");
  auto genFilteredDF = genDF2.Filter("std::abs(fPDGcode)==3334 && std::abs(yOmega) < 0.5");
  // loop over cent bins
  auto genPtHist = genFilteredDF.Histo1D({"genPtHist", ";#it{p}_{T} (GeV/c)", NPTBINS, pt}, "fgPt");
  ROOT::RDF::RResultPtr<TH1D> mcPtHist[2][NCENTBINS];
  TH1D* dataPtHist[2][NCENTBINS];
  
  TFile output("outputMB.root", "recreate");
  constexpr int maxIter{2};
 
  //
  TChain dataChainT(Form("O2npcasctable%s", names[0].c_str()));
  fillChainFromAO2D(dataChainT, dataFilename);
  TChain mcChainT(Form("O2npcasctablemc%s", names[0].c_str()));
  fillChainFromAO2D(mcChainT, mcFilename);
  TChain dataChainN(Form("O2npcasctable%s", names[1].c_str()));
  fillChainFromAO2D(dataChainN, dataFilename);
  TChain mcChainN(Form("O2npcasctablemc%s", names[1].c_str()));
  fillChainFromAO2D(mcChainN, mcFilename);
  std::vector<std::reference_wrapper<TChain>> dataChains = {dataChainT, dataChainN};
  std::vector<std::reference_wrapper<TChain>> mcChains = {mcChainT, mcChainN};
  for(int imult{0}; imult < NCENTBINS; imult++) 
  {
    std::cout << "Doing multiplicity " << imult << ":" << centbins[imult] << "-" << centbins[imult+1] << std::endl;
    auto multdir = output.mkdir(Form("cent%i",imult));
    multdir->cd();
    for (int iNt{0}; iNt < maxIter; ++iNt)
    {
      auto dir = multdir->mkdir(Form("Omega%s", names[iNt].c_str()));
      dir->cd();
      TChain& dataChain = dataChains[iNt];
      TChain& mcChain = mcChains[iNt];
      ROOT::RDataFrame dataDF(dataChain);
      double centMin = centbins[imult];
      double centMax = centbins[imult+1];
      //auto dataFilteredCentDF = dataDF.Filter([=](float cent){return (cent > centMin) && (cent < centMax);},{"fCentFT0M"});
      auto dataFilteredCentDF = dataDF.Filter(Form("(fMultNTracksGlobal > %f) && (fMultNTracksGlobal <%f) && fNoSameBunchPileup && fSel8", centMin, centMax));
      auto dataFilteredDF = dataFilteredCentDF.Filter(selectionString.Data());
      
      ROOT::RDataFrame mcDFNt(mcChain);
      auto mcFilteredDFNt = mcDFNt.Filter((selectionString + "&&std::abs(fPDGcode)==3334").Data());
      auto dataPtMassHist = dataFilteredDF.Histo2D({Form("dataPtMassHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Invariant mass (GeV/c^{2})", NPTBINS, pt, 50, 1.65, 1.71}, "fCascPt", "fMassOmega");
      auto mcPtMassHist = mcFilteredDFNt.Histo2D({Form("mcPtMassHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Invariant mass (GeV/c^{2})", NPTBINS, pt, 60, 1.65, 1.71}, "fCascPt", "fMassOmega");
      mcPtHist[iNt][imult] = mcFilteredDFNt.Histo1D({Form("mcPtHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c)", NPTBINS, pt}, "fCascPt");
      dataPtHist[iNt][imult] = new TH1D(Form("dataPtHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);d#it{N}_{raw}/d#it{p}_{T}", NPTBINS, pt);

      TH1D* purity = new TH1D(Form("purity%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Purity", NPTBINS, pt);
      dir->mkdir("fits")->cd();
      for (int i{0}; i < NPTBINS; ++i)
      {
        auto dataPtMassHistSlice = dataPtMassHist->ProjectionY(Form("dataPtMassHistSlice%s", names[iNt].data()), i + 1, i + 1);
        auto mcPtMassHistSlice = mcPtMassHist->ProjectionY(Form("mcPtMassHistSlice%s", names[iNt].data()), i + 1, i + 1);
        RooRealVar m("m", "m_{#Omega}", 1.65, 1.71, "GeV/#it{c}^{2}");
        RooDataHist data("data", "data", m, dataPtMassHistSlice);
        RooDataHist mcnt("mcnt", "mcnt", m, mcPtMassHistSlice);
        RooRealVar mean("mean", "#mu", 1.67245, 1.66, 1.69);
        RooRealVar sigma("sigma", "#sigma", 0.002, 0.0001, 0.01);
        RooRealVar n0("n0", "n_{0}", 0.1, 0., 20.);
        RooRealVar alpha0("alpha0", "#alpha_{0}", 1., 0., 10.);
        RooRealVar n1("n1", "n_{1}", 0.1, 0., 20.);
        RooRealVar alpha1("alpha1", "#alpha_{1}", 1., 0., 10.);
        RooCrystalBall cb("cb", "cb", m, mean, sigma, alpha0, n0, alpha1, n1);
        RooRealVar fsig("fsig", "f_{sig}", 0.5, 0., 1.);
        RooRealVar tau("tau", "#tau", -2, -200, 4);
        RooExponential exp("exp", "exp", m, tau);
        RooAddPdf model("model", "model", RooArgList(cb, exp), RooArgList(fsig));
        cb.fitTo(mcnt);
        RooPlot *mcframe = m.frame();
        mcframe->SetTitle(Form("MC %1.1f#leq #it{p}_{T} < %1.1f GeV/#it{c}", pt[i], pt[i + 1]));
        mcnt.plotOn(mcframe);
        cb.plotOn(mcframe);
        cb.paramOn(mcframe, RooFit::Format("TEU", RooFit::AutoPrecision(1)));
        mcframe->Write(Form("mcfit_%d", i));
        n0.setConstant();
        n1.setConstant();
        alpha0.setConstant();
        alpha1.setConstant();
        model.fitTo(data);
        RooPlot *dataframe = m.frame();
        dataframe->SetTitle(Form("Data %1.1f#leq #it{p}_{T} < %1.1f GeV/#it{c}", pt[i], pt[i + 1]));
        data.plotOn(dataframe);
        model.plotOn(dataframe);
        model.plotOn(dataframe, RooFit::Components("exp"), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue));
        model.paramOn(dataframe, RooFit::Format("TEU", RooFit::AutoPrecision(1)));
        m.setRange("fitRange", 1.66, 1.685);
        auto modelInt = model.createIntegral(RooArgSet(m), RooFit::NormSet(RooArgSet(m)), RooFit::Range("fitRange"));
        auto expInt = exp.createIntegral(RooArgSet(m), RooFit::NormSet(RooArgSet(m)), RooFit::Range("fitRange"));
        double totalIntegral = modelInt->getVal();
        double backgroundIntegral = expInt->getVal() * (1. - fsig.getVal());
        purity->SetBinContent(i + 1, 1. - backgroundIntegral / totalIntegral);
        dataframe->Write(Form("datafit_%d", i));
        dataPtHist[iNt][imult]->SetBinContent(i + 1, fsig.getVal() * data.sumEntries());
        dataPtHist[iNt][imult]->SetBinError(i + 1, fsig.getError() * data.sumEntries());
      }
      dataPtHist[iNt][imult]->Scale(1., "width");

      dir->cd();
      purity->Write();
      TH1D efficiency(Form("efficiency%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Efficiency", NPTBINS, pt);
      efficiency.Divide(mcPtHist[iNt][imult].GetPtr(), genPtHist.GetPtr(), 1., 1., "b");
      efficiency.Write();

      TCanvas cComp(Form("spectra_comparison%s", names[iNt].c_str()));
      TH1 *spectrum = (TH1 *)dataPtHist[iNt][imult]->Clone(Form("spectrum%s", names[iNt].data()));
      //
      double aeff = mcPtHist[iNt][imult].GetPtr()->GetEntries()/genPtHist.GetPtr()->GetEntries();
      std::cout << iNt << " aver eff:" << aeff << " 3 Raw:" << spectrum->GetEntries() << std::endl;
      int Nbins = spectrum->GetNbinsX();
      int NbinsE = efficiency.GetNbinsX();
    
      if(Nbins != NbinsE) {
        std::cout << "FATAL: spectrum N bins:" << Nbins << " eff Nbins:" << NbinsE << std::endl;
        exit(1);
      }
      for(int i = 0; i < Nbins; i++){
        std::cout << i << " Raw:" << spectrum->GetBinContent(i) << " eff:" << efficiency.GetBinContent(i) << std::endl;
      }
      spectrum->Divide(&efficiency);
      spectrum->Scale(normalisations[iNt] / normCent[imult]);
      published_stat.Draw("e x0");
      published_syst.Draw("e2 same");
      spectrum->Draw("same");
      dir->cd();
      cComp.Write();
      spectrum->Write();

      TH1D *ratio_stat_nt = (TH1D *)spectrum->Clone(Form("ratio_stat%s", names[iNt].data()));
      ratio_stat_nt->SetTitle(";#it{p}_{T} (GeV/c);Data/Published");
      ratio_stat_nt->Divide(&published_stat);
      ratio_stat_nt->Write();
    }

    multdir->cd();
    if (maxIter > 1) {
      TH1D* str_eff_mc = new TH1D("str_eff_mc", ";#it{p}_{T} (GeV/c);Strangeness tracking efficiency MC", NPTBINS, pt);
      str_eff_mc->Divide(mcPtHist[0][imult].GetPtr(), mcPtHist[1][imult].GetPtr(), 1., 1., "b");
      str_eff_mc->Write();
      TH1D* str_eff_data = new TH1D("str_eff_data", ";#it{p}_{T} (GeV/c);Strangeness tracking efficiency data", NPTBINS, pt);
      str_eff_data->Divide(dataPtHist[0][imult], dataPtHist[1][imult], 1., 1., "b");
      str_eff_data->Write();
      TH1D* eff_ratio = new TH1D("eff_ratio", ";#it{p}_{T} (GeV/c);ST efficiency ratio MC / data", NPTBINS, pt);
      eff_ratio->Divide(str_eff_data, str_eff_mc);
      eff_ratio->Fit("pol0");
      eff_ratio->Write();
    }
  }
  output.cd();
  TH1D* rawYieldvsCent[2];
  rawYieldvsCent[0] = new TH1D("hOmegaTracked","hOmegaTracked", NCENTBINS , centbins);
  rawYieldvsCent[1] = new TH1D("hOmegaNT","hOmegaNT", NCENTBINS , centbins);
  for(int i = 0; i < maxIter; i++){
    for(int imult = 0; imult < NCENTBINS ; imult++) {
      double error = 0;
      double integral = dataPtHist[i][imult]->IntegralAndError(0,NPTBINS,error);
      rawYieldvsCent[i]->SetBinContent(imult + 1, integral);
      rawYieldvsCent[i]->SetBinError(imult + 1, error);
    }
    rawYieldvsCent[i]->Write();
  }
}