#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

void LabelCFHisto(TF1 *h1, TString Title, EColor color, double size) {
    h1->SetMarkerColor(color);
    h1->SetLineColor(color);
    h1->SetMarkerSize(size);
    h1->SetMarkerStyle(20);
    h1->SetTitle(Title.Data());
}
void LabelCFHisto(TH1F *h1, TString Title, EColor color, double size) {
    h1->SetMarkerColor(color);
    h1->SetLineColor(color);
    h1->SetMarkerSize(size);
    h1->SetMarkerStyle(20);
    h1->SetTitle(Title.Data());
}
void MakeCFPlot1(TString Title, TString BLType, TH1F*h, TF1 *h1, TF1 *g1, TF1 *g2, double a, double b, double c, double aE, double bE, double cE, double chi) {

    LabelCFHisto(g1, "", kBlack, 0);
    LabelCFHisto(g2, "", kRed, 0);
    LabelCFHisto(h1, "", kGreen, 0);
    LabelCFHisto(h, "data", kBlue, 1);

    TCanvas *CFPlot1 = new TCanvas("CFPlot1", "CFplot1", 0, 0, 800, 600);
    //TString pdfCFName = "FromCATS.png";
    TString pdfCFName = "FitResultReso_Pol0" + Title + ".png";
    //TString pdfCFName = "FitResults_NLO_LowR.png";
    //TString pdfCFName = "FitResults_Usmani_LowR.png";
    TLegend *leg = new TLegend( 0.45, 0.47, 0.75, 0.58);
    leg->SetTextSize(0.035);
    leg->AddEntry(g1, "BaseLine pol0");
    leg->AddEntry(g2, "Modeled Ck");
    leg->AddEntry(h1, "Fit");
    leg->AddEntry(h, "Data(pd #oplus #bar{p}-#bar{d})");
    leg->SetLineColor(0);

    TString textgUNO1    = "BL : ";
    TString textgUNO2    = "a = ";
    TString textgUNO3    = "b = ";
    TString textgUNO4    = "c = ";
    TString textgUNO5    = "#chi^{2}/NDF = ";
    TLegend *leggUNO = new TLegend(0.35, 0.25, 0.8, 0.45);
    leggUNO->SetTextSize(0.03);
    TString textgUNO    = "Source size (Reso) = 0.92 fm";
    textgUNO1    = textgUNO1 + BLType;
    textgUNO2    = textgUNO2 + Form(" %.4f #pm", a) + Form("%.6f", aE);
    textgUNO3    = textgUNO3 + Form(" %.3e #pm", b) + Form(" %.3e MeV^{-1}", bE);
    textgUNO4    = textgUNO4 + Form(" %.3e #pm", c) + Form(" %.3e MeV^{-2}", cE);
    textgUNO5    = textgUNO5 + Form(" %.4f", chi);
    leggUNO->AddEntry((TObject*)0, textgUNO, "");
    leggUNO->AddEntry((TObject*)0, textgUNO1, "");
    leggUNO->AddEntry((TObject*)0, textgUNO2, "");
    leggUNO->AddEntry((TObject*)0, textgUNO3, "");
    leggUNO->AddEntry((TObject*)0, textgUNO4, "");
    leggUNO->AddEntry((TObject*)0, textgUNO5, "");
    leggUNO->SetLineColor(0);

    /*g->Add(g1);
    //g->Add(g3);
    g->Add(g2);
    */
    h->SetTitle(Title);
    h->GetXaxis()->SetTitle("#it{k*} (MeV)");
    h->GetYaxis()->SetTitle("#it{C}(k*)");
    h-> GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetXaxis()->SetRangeUser(0, 432);
    h->GetYaxis()->SetRangeUser(-0.2, 1.4);

    h->Draw("Ep");
    g1->Draw("same");
    g2->Draw("same");
    h1->Draw("same");
    leg->Draw("same");

    leggUNO->Draw("same");
    CFPlot1->Print(pdfCFName);
    delete CFPlot1;
    return;
}

void MakeCFPlot1(TString Title, TString BLType, TH1F*h, TF1 *h1, TF1 *g1, TF1 *g2, double a, double b, double aE, double bE, double chi) {

    LabelCFHisto(g1, "", kBlack, 0);
    LabelCFHisto(g2, "", kRed, 0);
    LabelCFHisto(h1, "", kGreen, 0);
    LabelCFHisto(h, "data", kBlue, 1);

    TCanvas *CFPlot1 = new TCanvas("CFPlot1", "CFplot1", 0, 0, 800, 600);
    //TString pdfCFName = "FromCATS.png";
    TString pdfCFName = "FitResult_Pol1" + Title + ".png";
    //TString pdfCFName = "FitResults_NLO_LowR.png";
    //TString pdfCFName = "FitResults_Usmani_LowR.png";
    TLegend *leg = new TLegend( 0.45, 0.47, 0.75, 0.58);
    leg->SetTextSize(0.035);
    leg->AddEntry(g1, "BaseLine pol1");
    leg->AddEntry(g2, "Modeled Ck");
    leg->AddEntry(h1, "Fit");
    leg->AddEntry(h, "pd #oplus #bar{p}-#bar{d}");
    leg->SetLineColor(0);

    TString textgUNO1    = "BL : ";
    TString textgUNO2    = "a = ";
    TString textgUNO3    = "b = ";
    TString textgUNO4    = "#chi^{2}/NDF = ";
    TLegend *leggUNO = new TLegend(0.35, 0.25, 0.8, 0.45);
    leggUNO->SetTextSize(0.03);
    TString textgUNO    = "Source size (Core +Reso) = 1.0 fm";
    textgUNO1    = textgUNO1 + BLType;
    textgUNO2    = textgUNO2 + Form(" %.4f #pm", a) + Form("%.6f", aE);
    textgUNO3    = textgUNO3 + Form(" %.3e #pm", b) + Form(" %.3e MeV^{-1}", bE);
    textgUNO4    = textgUNO4 + Form(" %.4f", chi);
    leggUNO->AddEntry((TObject*)0, textgUNO, "");
    leggUNO->AddEntry((TObject*)0, textgUNO1, "");
    leggUNO->AddEntry((TObject*)0, textgUNO2, "");
    leggUNO->AddEntry((TObject*)0, textgUNO3, "");
    leggUNO->AddEntry((TObject*)0, textgUNO4, "");
    leggUNO->SetLineColor(0);

    /*g->Add(g1);
    //g->Add(g3);
    g->Add(g2);
    */
    h->SetTitle(Title);
    h->GetXaxis()->SetTitle("#it{k*} (MeV)");
    h->GetYaxis()->SetTitle("#it{C}(k*)");
    h-> GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetXaxis()->SetRangeUser(0.0, 300);
    h->GetYaxis()->SetRangeUser(-0.1, 1.3);

    h->Draw("Ep");
    g1->Draw("same");
    g2->Draw("same");
    h1->Draw("same");
    leg->Draw("same");

    leggUNO->Draw("same");
    CFPlot1->Print(pdfCFName);
    delete CFPlot1;
    return;
}