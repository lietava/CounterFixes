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

void analyseTRG(TString dataFilename = "dataTRG/513390/AO2D.root", TString mcFilename = "mcMB/AO2D.root", TString normalisationFilename = "dataTRG/513390/AnalysisResults.root")
{
  TString selectionString{"std::abs(fProtonEta) < 0.9&&std::abs(fPionEta)<0.9&&fCascCosPA > 0.9995&&std::abs(fBachEta)<0.9&&fBachNClusTPC>70&&std::abs(fMassXi-1.32171)>0.008&& std::abs(fPvZ)<10"};

  ROOT::EnableImplicitMT();

  double val[7]{0.0009488760957, 0.0006952320908, 0.0004759373703, 0.0002954698723, 0.0001662988815, 8.73277117e-05, 3.635092251e-05};
  double stat[7]{4.808611514e-05, 2.719983815e-05, 1.605100658e-05, 1.091321782e-05, 7.543390105e-06, 5.044361557e-06, 2.234297454e-06};
  double syst[7]{5.649579873e-05, 2.975781163e-05, 2.84376873e-05, 9.258206825e-06, 3.867393999e-06, 4.169908233e-06, 2.43093823e-06};
  TH1D published_stat("published_stat", ";#it{p}_{T} (GeV/c);d^{2}#it{N}/d#it{p}_{T}d#it{y} (c/GeV)", 7, pt);
  TH1D published_syst("published_syst", ";#it{p}_{T} (GeV/c);d^{2}#it{N}/d#it{p}_{T}d#it{y} (c/GeV)", 7, pt);
  for (int i{1}; i <= 7; ++i)
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

  TFile normalisationFile(normalisationFilename.Data());
  TH1 *hCounterTVX = static_cast<TH1 *>(normalisationFile.Get("bc-selection-task/hCounterTVX"));
  if(hCounterTVX == nullptr) {
    hCounterTVX = static_cast<TH1 *>(normalisationFile.Get("eventselection-run3/luminosity/hCounterTVX"));
    if(hCounterTVX == nullptr) {
      std::cout << "TVX counter not available" << std::endl;
      exit(1);
    }
  }
  ZorroSummary *zorroSummary = (ZorroSummary *)normalisationFile.Get("non-prompt-cascade-task/ZorroSummary");
  if(zorroSummary == nullptr) {
    std::cout << " No zorroSummary" << std::endl;
    exit(1);
  }
  std::cout << "TOIs:" << zorroSummary->getTOInames() << std::endl;
  // find triggers
  uint32_t toiHM = 1;
  uint32_t toiTO = 2;
  float normHM = 0.8/100.;
  //double normalisations[2]{0.7558 / zorroSummary->getNormalisationFactor(0) / 0.857, 0.7558 /zorroSummary->getNormalisationFactor(0)};
  //                                       Tracked                                            Nontracked
  double normalisations[2]{0.7558 / zorroSummary->getNormalisationFactor(1) / 0.857, 0.7558 /zorroSummary->getNormalisationFactor(0)/ normHM };
  std::cout << "Zorro normalisation 0 " << zorroSummary->getNormalisationFactor(0) << " 1:" <<  zorroSummary->getNormalisationFactor(1) << std::endl;
  //double normalisations[2]{0.7558 / hCounterTVX->GetEntries(), 0.7558 / hCounterTVX->GetEntries()};

  TChain genChain("O2npcasctablegen");
  fillChainFromAO2D(genChain, mcFilename);
  ROOT::RDataFrame genDF(genChain);
  auto genDF2 = genDF.Define("pz", "fgPt*std::sinh(fgEta)").Define("eOmega", "std::sqrt(fgPt*fgPt + pz*pz + 1.67245 * 1.67245)").Define("yOmega", "0.5*std::log((eOmega + pz)/(eOmega - pz))");
  auto genFilteredDF = genDF2.Filter("std::abs(fPDGcode)==3334 && std::abs(yOmega) < 0.5");
  auto genPtHist = genFilteredDF.Histo1D({"genPtHist", ";#it{p}_{T} (GeV/c)", 7, pt}, "fgPt");
  ROOT::RDF::RResultPtr<TH1D> mcPtHist[2];

  TH1D* dataPtHist[2];
  ROOT::RDF::RResultPtr<TH1D> dataCentHist[2];
  ROOT::RDF::RResultPtr<TH1D> dataFT0MHist[2];

  TFile output("outputTRG.root", "recreate");
  std::string names[2]{"", "nt"};
  constexpr int maxIter{2};
  for (int iNt{0}; iNt < maxIter; ++iNt)
  {
    auto dir = output.mkdir(Form("Omega%s", names[iNt].c_str()));
    dir->cd();
    TChain dataChain(Form("O2npcasctable%s", names[iNt].c_str()));
    fillChainFromAO2D(dataChain, dataFilename);
    TChain mcChain(Form("O2npcasctablemc%s", names[iNt].c_str()));
    fillChainFromAO2D(mcChain, mcFilename);

    ROOT::RDataFrame dataDF(dataChain);
    auto dataFilteredDF = dataDF.Filter(selectionString.Data());
    uint32_t toitrg = 0;
    float centCut = 0;
    if(iNt == 0) {
      toitrg = toiTO;
      centCut = 101;
    } else {
      toitrg = toiHM;
      centCut = 0.8;
    }
    dataFilteredDF = dataFilteredDF.Filter( Form("fToiMask & %d && fCentFT0M < %f && fNoSameBunchPileup",toitrg,centCut));
    dataCentHist[iNt] =  dataFilteredDF.Histo1D({"Cent","Cent",1000,0.,10.},"fCentFT0M");
    dataFT0MHist[iNt] =  dataFilteredDF.Histo1D({"FT0M","FT0M",7000,0.,7000.},"fMultFT0M");

    ROOT::RDataFrame mcDFNt(mcChain);
    auto mcFilteredDFNt = mcDFNt.Filter((selectionString + "&&std::abs(fPDGcode)==3334").Data());

    auto dataPtMassHist = dataFilteredDF.Histo2D({Form("dataPtMassHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Invariant mass (GeV/c^{2})", 7, pt, 50, 1.65, 1.71}, "fCascPt", "fMassOmega");
    auto mcPtMassHist = mcFilteredDFNt.Histo2D({Form("mcPtMassHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Invariant mass (GeV/c^{2})", 7, pt, 60, 1.65, 1.71}, "fCascPt", "fMassOmega");
    mcPtHist[iNt] = mcFilteredDFNt.Histo1D({Form("mcPtHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c)", 7, pt}, "fCascPt");
    dataPtHist[iNt] = new TH1D(Form("dataPtHist%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);d#it{N}_{raw}/d#it{p}_{T}", 7, pt);

    TH1D* purity = new TH1D(Form("purity%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Purity", 7, pt);
    dir->mkdir("fits")->cd();
    for (int i{0}; i < 7; ++i)
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
      dataPtHist[iNt]->SetBinContent(i + 1, fsig.getVal() * data.sumEntries());
      dataPtHist[iNt]->SetBinError(i + 1, fsig.getError() * data.sumEntries());
    }
    dataPtHist[iNt]->Scale(1., "width");

    dir->cd();
    purity->Write();
    TH1D efficiency(Form("efficiency%s", names[iNt].data()), ";#it{p}_{T} (GeV/c);Efficiency", 7, pt);
    efficiency.Divide(mcPtHist[iNt].GetPtr(), genPtHist.GetPtr(), 1., 1., "b");
    efficiency.Write();

    TCanvas cComp(Form("spectra_comparison%s", names[iNt].c_str()));
    TH1 *spectrum = (TH1 *)dataPtHist[iNt]->Clone(Form("spectrum%s", names[iNt].data()));
    spectrum->Divide(&efficiency);
    spectrum->Scale(normalisations[iNt]);
    published_stat.Draw("e x0");
    published_syst.Draw("e2 same");
    spectrum->Draw("same");
    dir->cd();
    cComp.Write();
    spectrum->Write();
    dataCentHist[iNt]->Write();
    dataFT0MHist[iNt]->Write();

    TH1D *ratio_stat_nt = (TH1D *)spectrum->Clone(Form("ratio_stat%s", names[iNt].data()));
    ratio_stat_nt->SetTitle(";#it{p}_{T} (GeV/c);Data/Published");
    ratio_stat_nt->Divide(&published_stat);
    ratio_stat_nt->Write();
  }

  output.cd();
  if (maxIter > 1) {
    TH1D* str_eff_mc = new TH1D("str_eff_mc", ";#it{p}_{T} (GeV/c);Strangeness tracking efficiency MC", 7, pt);
    str_eff_mc->Divide(mcPtHist[0].GetPtr(), mcPtHist[1].GetPtr(), 1., 1., "b");
    str_eff_mc->Write();
    TH1D* str_eff_data = new TH1D("str_eff_data", ";#it{p}_{T} (GeV/c);Strangeness tracking efficiency data", 7, pt);
    str_eff_data->Divide(dataPtHist[0], dataPtHist[1], 1., 1., "b");
    str_eff_data->Write();
    TH1D* eff_ratio = new TH1D("eff_ratio", ";#it{p}_{T} (GeV/c);ST efficiency ratio MC / data", 7, pt);
    eff_ratio->Divide(str_eff_data, str_eff_mc);
    eff_ratio->Fit("pol0");
    eff_ratio->Write();
  }
}