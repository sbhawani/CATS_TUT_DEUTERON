#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"

void Plots() {

  const char *pdFile = "/home/sbhawani/GentleFemto/install/DreamFunctions/NanoAODLatest40MeV.root";
  //const char *pdFile2 = argv[2];
  //const char *outputFileName= argv[3];

  TFile *CF_pdFile1 = TFile::Open(pdFile, "read");
  TList*PairDist = (TList*)CF_pdFile1->FindObjectAny("AntiPairDist");
  TList*PairReweighted = (TList*)PairDist->FindObject("PairReweighted");
  //TFile *CFpdFile2  = TFile::Open(pdFile2, "read");
  //TGraphErrors *GraphCF1 = (TGraphErrors*) CF_pdFile1->Get("Graph_from_hCk_Reweighted_1MeV");
  //TCanvas*Canv = (TCanvas*) CFpdFile2->FindObjectAny("Canvas_1");
  //TGraphErrors *GraphCF2 = (TGraphErrors*) Canv->GetListOfPrimitives()->FindObject("Graph_from_hCk_Reweighted_0MeV");
  TH1F *GraphCF1 = (TH1F*) CF_pdFile1->Get("hCk_ReweightedMeV_1");
  //TH1F *GraphCF1 = (TH1F*) PairReweighted->FindObject("CFDist_Particle1_Particle3_clone_Shifted_FixShifted_Rebinned_20_Reweighted");
 // TH1F *GraphCF2 = (TH1F*) Canv->GetListOfPrimitives()->FindObject("hCk_Rebinned_7");
  //TH1F *GraphCF2 = (TH1F*) Canv->GetListOfPrimitives()->FindObject("CFDist_Particle1_Particle3_clone");
  //TH1F  *GraphCF1 = (TH1F*)GraphCF->Clone("GraphCF1");
  //GraphCF1->Sumw2();
  if (!GraphCF1) {
    std::cout << "hCk_RebinnedMeV_0 from first file1 missing\n";
    return -1;
  }
  /*if (!GraphCF2) {
    std::cout << "hCk_RebinnedMeV_0 from first file2 missing\n";
    return -1;
  }*/

  /*TGraphErrors *GraphCF2 = new TGraphErrors();
  for(int uBin=0;uBin<GraphCF22->GetNbinsX();uBin++){
     x = GraphCF22->GetBinCenter(uBin);
     y = GraphCF22->GetBinContent(uBin);
    CF_xErr = GraphCF22->GetXaxis()->GetBinWidth(uBin);
    CF_yErr = GraphCF22->GetBinError(uBin);
    GraphCF2->SetPoint(uBin,x,y);
    GraphCF2->SetPointError(uBin,CF_xErr,CF_yErr);
  }*/


  std::cout<<"Reached-0\n";
  //TFile *out = TFile::Open(Form("%s.root", outputFileName), "recreate");
 // out->cd();
  TCanvas *c4 = new TCanvas("c8", "c8", 1200, 800);
  //TCanvas *c4 = new TCanvas();
  TLegend *leg = new TLegend(0.6, 0.3, 0.73, 0.4);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
//  leg->SetNColumns(2);
  //leg->SetTextSizePixels(40);
  std::cout<<"Reached-1\n";
  GraphCF1->SetTitle("; #it{k}*  (MeV); #it{C}(k*)");
  GraphCF1->GetXaxis()->SetTitleSize(40);
  GraphCF1->GetYaxis()->SetTitleSize(40);
  GraphCF1->GetXaxis()->SetTitleOffset(1.35);
  GraphCF1->GetYaxis()->SetTitleOffset(1.4);
  GraphCF1->GetXaxis()->SetLabelSize(40);
  GraphCF1->GetYaxis()->SetLabelSize(40);
  GraphCF1->GetXaxis()->SetLabelOffset(.02);
  GraphCF1->GetYaxis()->SetLabelOffset(.02);
  GraphCF1->SetMarkerColor(kRed + 1);
  GraphCF1->SetLineColor(kRed + 1);
  GraphCF1->SetMarkerStyle(8);
  GraphCF1->SetMarkerSize(1.0);
  GraphCF1->GetXaxis()->SetRangeUser(0.0, 400);

  //GraphCF1->GetXaxis()->SetLimits(0.0, 0.400);
  GraphCF1->GetYaxis()->SetRangeUser(0.0, 1.5);

 // GraphCF2->SetMarkerColor(kBlue + 1);
 // GraphCF2->SetLineColor(kBlue + 1);
  //GraphCF1->SetFillColorAlpha(kBlue + 1, 0.7);

 // GraphCF2->Rebin(5);
 // GraphCF2->Scale(0.2);
  //GraphCF2->SetMarkerStyle(8);
  //GraphCF2->SetMarkerSize(1.0);
  std::cout<<"Reached-2\n";
  GraphCF1->Draw("ep");
  // GraphCF1->Draw("AZ");
  //GraphCF2->Draw("2PSAME");

  std::cout<<"Reached-3\n";

 leg->AddEntry(GraphCF1, "p#minus#kern[0.4]{d} #oplus #bar{p}#minus#kern[0.4]{#bar{d}}");
  //leg->AddEntry(GraphCF2,
          //       "p#minus#kern[0.4]{d} #oplus #bar{p}#minus#kern[0.4]{#bar{d}} Michael");

  //leg->AddEntry(GraphCF1,
                 //  "#bar{p}#minus#kern[0.4]{#bar{d}} Bhawani");
   // leg->AddEntry(GraphCF2,
                 //  "#bar{p}#minus#kern[0.4]{#bar{d}} Michael");
  //leg->AddEntry(GraphCF1,
               //  "p#minus#kern[0.4]{d} #it{C}(k*)");
  // leg->AddEntry(GraphCF1,
           //      "p#minus#kern[0.4]{d} Avg Ref. Mult_{#eta < 0.8} ~30.1", "pef");
  // leg->AddEntry(GraphCF2, "p#minus#kern[0.4]{p} Avg Ref. Mult_{#eta < 0.8} ~23",
            //     "pef");

   leg->Draw("same");
   //c4->SaveAs(Form("%s/NanopDApAdZoomForPAG.pdf", gSystem->pwd()));
   c4->SaveAs("NanopDApAdZoomForPAG.pdf");

  // GraphCF1->Write();
   //GraphCF2->Write();
  // c4->Write();
   //out->Write();
   //out->Close();
   CF_pdFile1->Close();
  // CFpdFile2->Close();
}
