#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"

struct DnDetaCentBin {
  double centMin;
  double centMax;
  TString fileName;
};

double integralWithError(TH1D* hist, double& error)
{
  if (!hist) {
    error = 0.0;
    return 0.0;
  }
  return hist->IntegralAndError(1, hist->GetNbinsX(), error);
}

void plotDnDetaVsCentrality(TString inputDir = ".",
                            TString outputFile = "DnDetaVsCentrality.root",
                            double etaWidth = 1.6)
{
  std::vector<DnDetaCentBin> bins = {
    {0.0, 1.0, "manualBayes4D_cent0_1_pileup1.root"},
    {1.0, 5.0, "manualBayes4D_cent1_5_pileup1.root"},
    {5.0, 10.0, "manualBayes4D_cent5_10_pileup1.root"},
    {10.0, 15.0, "manualBayes4D_cent10_15_pileup1.root"},
    {15.0, 20.0, "manualBayes4D_cent15_20_pileup1.root"},
    {20.0, 30.0, "manualBayes4D_cent20_30_pileup1.root"},
    {30.0, 40.0, "manualBayes4D_cent30_40_pileup1.root"},
    {40.0, 50.0, "manualBayes4D_cent40_50_pileup1.root"},
    {50.0, 70.0, "manualBayes4D_cent50_70_pileup1.root"},
    {70.0, 100.0, "manualBayes4D_cent70_100_pileup1.root"},
  };

  if (etaWidth <= 0.0) {
    std::cerr << "etaWidth must be positive" << std::endl;
    return;
  }

  if (!inputDir.IsNull() && !inputDir.EndsWith("/")) {
    inputDir += "/";
  }

  std::vector<double> centEdges;
  centEdges.reserve(bins.size() + 1);
  centEdges.push_back(bins.front().centMin);
  for (const auto& bin : bins) {
    centEdges.push_back(bin.centMax);
  }

  auto* hDnDeta = new TH1D("hDnDetaVsCentrality",
                           ";FT0M centrality (%);dN_{ch}/d#eta",
                           static_cast<int>(bins.size()), centEdges.data());
  auto* gDnDeta = new TGraphErrors(static_cast<int>(bins.size()));
  gDnDeta->SetName("gDnDetaVsCentrality");
  gDnDeta->SetTitle(";FT0M centrality (%);dN_{ch}/d#eta");

  // Table 6.1 average over the four periods LHC22o, LHC24al, LHC24ao, LHC24af.
  const std::vector<double> table61Average = {
    21.6925, 18.3925, 15.7000, 13.8625, 12.4775,
    10.8500, 9.0925, 7.6400, 5.8975, 3.7500
  };
  const std::vector<double> table61AverageErr = {
    0.1613, 0.1313, 0.1139, 0.1027, 0.0925,
    0.0815, 0.0702, 0.0602, 0.0487, 0.0314
  };

  for (std::size_t i = 0; i < bins.size(); ++i) {
    const auto& bin = bins[i];
    const TString filePath = inputDir + bin.fileName;
    TFile input(filePath);
    if (input.IsZombie()) {
      std::cerr << "Cannot open " << filePath << std::endl;
      continue;
    }

    auto* hYield = dynamic_cast<TH1D*>(input.Get("hUnfoldMultProj"));
    auto* hEvents = dynamic_cast<TH1D*>(input.Get("hUnfoldMultEventProj"));
    if (!hYield || !hEvents) {
      std::cerr << filePath
                << " is missing hUnfoldMultProj and/or hUnfoldMultEventProj"
                << std::endl;
      continue;
    }

    double yieldErr = 0.0;
    double eventsErr = 0.0;
    const double yield = integralWithError(hYield, yieldErr);
    const double events = integralWithError(hEvents, eventsErr);
    if (events <= 0.0) {
      std::cerr << filePath << " has zero event integral" << std::endl;
      continue;
    }

    const double dNdeta = yield / events / etaWidth;
    const double relYieldErr = (yield > 0.0) ? yieldErr / yield : 0.0;
    const double relEventsErr = eventsErr / events;
    const double dNdetaErr = dNdeta * std::sqrt(relYieldErr * relYieldErr +
                                                relEventsErr * relEventsErr);

    const double cent = 0.5 * (bin.centMin + bin.centMax);
    const double centErr = 0.5 * (bin.centMax - bin.centMin);
    hDnDeta->SetBinContent(static_cast<int>(i) + 1, dNdeta);
    hDnDeta->SetBinError(static_cast<int>(i) + 1, dNdetaErr);
    gDnDeta->SetPoint(static_cast<int>(i), cent, dNdeta);
    gDnDeta->SetPointError(static_cast<int>(i), centErr, dNdetaErr);

    std::cout << Form("%5.1f-%5.1f%%  dN/deta = %.6g +/- %.6g  events = %.6g  yield = %.6g",
                      bin.centMin, bin.centMax, dNdeta, dNdetaErr, events, yield)
              << std::endl;
  }

  hDnDeta->SetMarkerStyle(20);
  hDnDeta->SetMarkerSize(1.35);
  hDnDeta->SetLineWidth(2);
  gDnDeta->SetMarkerStyle(20);
  gDnDeta->SetMarkerSize(1.35);
  gDnDeta->SetLineWidth(2);
  gDnDeta->SetMarkerColor(kBlack);
  gDnDeta->SetLineColor(kBlack);

  auto* gTableAverage = new TGraphErrors(static_cast<int>(bins.size()));
  gTableAverage->SetName("gTable6p1Average");
  gTableAverage->SetTitle(";FT0M centrality (%);dN_{ch}/d#eta");
  gTableAverage->SetMarkerStyle(24);
  gTableAverage->SetMarkerSize(1.1);
  gTableAverage->SetMarkerColor(kRed + 1);
  gTableAverage->SetLineColor(kRed + 1);
  gTableAverage->SetLineWidth(2);

  auto* hRatioToTableAverage = new TH1D("hRatioToTable6p1Average",
                                        ";FT0M centrality (%);this work / Table 6.1 avg.",
                                        static_cast<int>(bins.size()), centEdges.data());
  auto* gRatioToTableAverage = new TGraphErrors(static_cast<int>(bins.size()));
  gRatioToTableAverage->SetName("gRatioToTable6p1Average");
  gRatioToTableAverage->SetTitle(";FT0M centrality (%);this work / Table 6.1 avg.");
  gRatioToTableAverage->SetMarkerStyle(20);
  gRatioToTableAverage->SetMarkerSize(1.1);
  gRatioToTableAverage->SetLineWidth(2);

  for (std::size_t i = 0; i < bins.size(); ++i) {
    const double cent = 0.5 * (bins[i].centMin + bins[i].centMax);
    const double centErr = 0.5 * (bins[i].centMax - bins[i].centMin);
    const double refCent = cent + 0.08 * (bins[i].centMax - bins[i].centMin);
    const double ref = table61Average[i];
    const double refErr = table61AverageErr[i];
    const double y = hDnDeta->GetBinContent(static_cast<int>(i) + 1);
    const double yErr = hDnDeta->GetBinError(static_cast<int>(i) + 1);
    const double ratio = (ref > 0.0) ? y / ref : 0.0;
    const double relYErr = (y > 0.0) ? yErr / y : 0.0;
    const double relRefErr = (ref > 0.0) ? refErr / ref : 0.0;
    const double ratioErr = ratio * std::sqrt(relYErr * relYErr + relRefErr * relRefErr);

    gTableAverage->SetPoint(static_cast<int>(i), refCent, ref);
    gTableAverage->SetPointError(static_cast<int>(i), 0.0, refErr);
    hRatioToTableAverage->SetBinContent(static_cast<int>(i) + 1, ratio);
    hRatioToTableAverage->SetBinError(static_cast<int>(i) + 1, ratioErr);
    gRatioToTableAverage->SetPoint(static_cast<int>(i), cent, ratio);
    gRatioToTableAverage->SetPointError(static_cast<int>(i), centErr, ratioErr);
  }
  hRatioToTableAverage->SetMinimum(0.7);
  hRatioToTableAverage->SetMaximum(1.3);
  hRatioToTableAverage->SetMarkerStyle(20);
  hRatioToTableAverage->SetMarkerSize(1.1);
  hRatioToTableAverage->SetLineWidth(2);

  auto* cDnDeta = new TCanvas("cDnDetaVsCentrality", "dN/deta vs centrality", 900, 850);
  auto* padTop = new TPad("padDnDeta", "dN/deta", 0.0, 0.32, 1.0, 1.0);
  auto* padBottom = new TPad("padRatio", "ratio", 0.0, 0.0, 1.0, 0.32);
  padTop->SetBottomMargin(0.02);
  padBottom->SetTopMargin(0.04);
  padBottom->SetBottomMargin(0.28);
  padTop->Draw();
  padBottom->Draw();

  padTop->cd();
  auto* hTopFrame = padTop->DrawFrame(-1.0, 0.0, 100.0, 23.0,
                                      ";FT0M centrality (%);dN_{ch}/d#eta");
  hTopFrame->GetXaxis()->SetLabelSize(0.0);
  gTableAverage->Draw("P E1 SAME");
  hDnDeta->Draw("E1 SAME");
  gDnDeta->Draw("P E1 SAME");
  auto* legend = new TLegend(0.55, 0.72, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(gDnDeta, "this work", "pe");
  legend->AddEntry(gTableAverage, "Table 6.1 period avg.", "pe");
  legend->Draw();

  padBottom->cd();
  hRatioToTableAverage->GetYaxis()->SetNdivisions(505);
  hRatioToTableAverage->GetYaxis()->SetTitleSize(0.08);
  hRatioToTableAverage->GetYaxis()->SetLabelSize(0.075);
  hRatioToTableAverage->GetXaxis()->SetTitleSize(0.09);
  hRatioToTableAverage->GetXaxis()->SetLabelSize(0.08);
  hRatioToTableAverage->Draw("E1");
  cDnDeta->Modified();
  cDnDeta->Update();

  TFile output(outputFile, "RECREATE");
  hDnDeta->Write();
  gDnDeta->Write();
  gTableAverage->Write();
  hRatioToTableAverage->Write();
  gRatioToTableAverage->Write();
  cDnDeta->Write();
  output.Close();

  TString pngName = outputFile;
  pngName.ReplaceAll(".root", ".png");
  cDnDeta->SaveAs(pngName);

  std::cout << "Wrote " << outputFile << " and " << pngName << std::endl;
  std::cout << "Used etaWidth = " << etaWidth << std::endl;
}
