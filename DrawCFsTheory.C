#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "DreamPlot.h"
#include "TFile.h"
#include "TDatabasePDG.h"

int main(int argc, char *argv[]) {
  if (!argv[1]) {
    std::cout << "SysFileMissing RadFile from latest task missing\n";
    return -1;
  }
  // if (!argv[2]) {
  //  std::cout << "pp  RadFile  from Bernie missing\n";
  //  return -1;
  // }
  if (!argv[2]) {
    std::cout << "Source Name\n";
    return -1;
  }

  const char *ppFile= argv[1];
  const char *sourceName= argv[2];
  const char *pdWF200= argv[3];
  const char *pdWF400= argv[4];
  const char *pdWF800= argv[5];
  const char *CoulombOnlyPd= argv[6];
  TFile *pdWF2002= TFile::Open(pdWF200, "read");
  TFile *pdWF4002= TFile::Open(pdWF400, "read");
  TFile *pdWF8002= TFile::Open(pdWF800, "read");
  TFile *CoulombOnlyPd2= TFile::Open(CoulombOnlyPd, "read");

  TGraph *g1= (TGraph*) pdWF2002->Get("gKitty");
  TGraph *g2= (TGraph*) pdWF4002->Get("gKitty");
  TGraph *g3= (TGraph*) pdWF8002->Get("gKitty");
  TGraph *g4= (TGraph*) CoulombOnlyPd2->Get("Ck_pd_Coulumb_only");
  //const char *pLNLOFile = argv[2];

  double yMin1= 1234567;
  double yMax1= 0;
  double x1, y1;
  double yMin2= 1234567;
  double yMax2= 0;
  double x2, y2;
  double yMin3= 1234567;
  double yMax3= 0;
  double x3, y3;
   double x4, y4;
/*
  for (int iBin= 0; iBin < g1->GetN(); iBin++) {
    g1->GetPoint(iBin, x1, y1);

    g1->SetPoint(iBin,x1, 0.75 * y1);
  }
  for (int iBin= 0; iBin < g2->GetN(); iBin++) {
      g2->GetPoint(iBin, x2, y2);

      g2->SetPoint(iBin,x2, 0.75 * y2);
    }

  for (int iBin= 0; iBin < g3->GetN(); iBin++) {
      g3->GetPoint(iBin, x3, y3);
      g3->SetPoint(iBin, x3, 0.75 * y3);
    }
      for (int iBin= 0; iBin < g4->GetN(); iBin++) {
      g4->GetPoint(iBin, x4, y4);
      g4->SetPoint(iBin, x4, 0.75 * y4);
    }

*/
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile *ppHMFile= TFile::Open(ppFile, "read");

  TGraphErrors *mTppHMSys= (TGraphErrors*) ppHMFile->Get("grError");
  TH1F *mTppHMStat= (TH1F*) ppHMFile->Get("histDefault");

  double yMin= 1234567;
  double yMax= 0;
  double x, y;
  for (int iBin= 0; iBin < mTppHMSys->GetN(); iBin++) {
    mTppHMSys->GetPoint(iBin, x, y);
    if (y < yMin) {
      yMin= y;
    }
    if (yMax < y) {
      yMax= y;
    }
    mTppHMSys->SetPointError(iBin, 0.5 * mTppHMSys->GetErrorX(iBin),
                             mTppHMSys->GetErrorY(iBin));
  }

  TFile *out= TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4= new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  //TLegend *leg= new TLegend(0.48, 0.35, 0.75, 0.45);
  TLegend *leg= new TLegend(0.6, 0.45, 0.70, 0.75);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
//  leg->SetNColumns(2);
  leg->SetTextSizePixels(40);

  mTppHMSys->SetTitle("; #it{k}*  (MeV/#it{c}); #it{C}(#it{k}*)");
  mTppHMSys->GetXaxis()->SetTitleSize(40);
  mTppHMSys->GetYaxis()->SetTitleSize(40);
  mTppHMSys->GetXaxis()->SetTitleOffset(1.35);
  mTppHMSys->GetYaxis()->SetTitleOffset(1.4);
  mTppHMSys->GetXaxis()->SetLabelSize(40);
  mTppHMSys->GetYaxis()->SetLabelSize(40);
  mTppHMSys->GetXaxis()->SetLabelOffset(.02);
  mTppHMSys->GetYaxis()->SetLabelOffset(.02);
  mTppHMSys->SetFillColorAlpha(kOrange - 3, 0.6);
  mTppHMSys->SetFillStyle(3102);
  mTppHMSys->SetLineColor(kOrange - 3);
  //mTppHMSys->GetXaxis()->SetRangeUser(0, 420);
  mTppHMSys->GetXaxis()->SetLimits(0, 310);
  mTppHMSys->GetYaxis()->SetRangeUser(0.0, 9.0);

  mTppHMStat->SetMarkerColor(kGray + 2);
  mTppHMStat->SetLineColor(kGray + 2);
  mTppHMStat->SetFillColorAlpha(kOrange - 3, 0.0);
  mTppHMStat->SetLineWidth(1);
  mTppHMStat->SetFillStyle(3102);
  mTppHMStat->SetMarkerStyle(20);
  mTppHMStat->SetMarkerSize(1.3);

  //mTppHMStat->Draw("2AP");
  //mTppHMSys->Draw("pez same");
  g1->SetMarkerColor(kBlack);
  g1->SetLineColor(kBlack + 2);
  g2->SetMarkerColor(kBlue + 2);
  g2->SetLineColor(kBlue + 2);
  g3->SetMarkerColor(kRed + 2);
  g3->SetLineColor(kRed + 2);
  g4->SetMarkerColor(kMagenta+ 2);
  g4->SetLineColor(kMagenta + 2);
  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);
  g4->SetLineWidth(2);


  mTppHMSys->Draw("2AP");
  mTppHMStat->Draw("pez same");
  g1->Draw("same");
  g2->Draw("same");
 // g3->Draw("same");
  g4->Draw("same");

  leg->AddEntry(mTppHMStat,
                "p#minus#kern[0.4]{d} #oplus #bar{p}#minus#kern[0.4]{#bard}",
                "pef");
  leg->AddEntry(g1,"p#minus#kern[0.4]{d} D wave ","pef");
  leg->AddEntry(g2,"p#minus#kern[0.4]{d} Q wave ","pef");
 // leg->AddEntry(g3,"p#minus#kern[0.4]{d} #Lambda = 800 MeV","pef");
  leg->AddEntry(g4,"p#minus#kern[0.4]{d} Pure Coulomb","pef");
  //leg->AddEntry(mTppHMSys, "Sys. uncertainties", "pef");

  TLatex BeamText;
  BeamText.SetTextFont(43);
  BeamText.SetTextSize(40);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.48, 0.86,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(
      0.48,
      0.79,
      "High-mult. (0#kern[-0.65]{ }#minus#kern[-0.65]{ } 0.17#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
  //BeamText.DrawLatex(0.48, 0.79, "High Mult");
  //BeamText.DrawLatex(0.48, 0.72, TString::Format("%s", sourceName).Data());
// BeamText.DrawLatex(0.48, 0.65, "Avg Ref. Mult_{#eta < 0.8} ~23");

  leg->Draw("same");
  c4->SaveAs(Form("%s/FinalCFwithStyematicErrrorforWFModelWithCoulombDwave-QWave.png", gSystem->pwd()));

  mTppHMSys->Write();
  mTppHMStat->Write();
  out->Write();
  c4->Write();
  out->Close();
  ppHMFile->Close();
}