
#include "ProtonDeuteronCF_fromWFModel.h"
#include "ExtendedCk.h"
#include "Basics.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>

void ProtonDeuteronCF_fromWFModelNLO() {
  const double SourceSize = 1.073;
  const double kMin = 2.5;
  const double kMax = 252.5;
  const unsigned NumMomBins = 50;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);
  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
//home/sbhawani/cernbox/ProtonDeuteron/WFs/PionlessEFT/Pd_CorrctNorm/MathematicaCalculation/
  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  ExternalWF = Init_pd_SebastianNewNLO("/home/sbhawani/cernbox/ProtonDeuteron/WFs/PionlessEFT/LO_NLO/", Kitty, 0, 3200);
  for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  }
  Kitty.SetChannelWeight(0, 1.0 / 3.0);
  Kitty.SetChannelWeight(1, 2.0 / 3.0);
  Kitty.SetQ1Q2(1);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/cernbox/ProtonDeuteron/Outputs/CATSOutput/PionlessEFTCFs/LOandNLO/NLOLO-L800-l0p5.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.1; RAD < 99; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}
void ProtonDeuteronCF_fromWFModel() {
  const double SourceSize = 1.073;
  const double kMin = 2.5;
  const double kMax = 132.5;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);
  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
//home/sbhawani/cernbox/ProtonDeuteron/WFs/PionlessEFT/Pd_CorrctNorm/MathematicaCalculation/
  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  ExternalWF = Init_pd_SebastianNew("/home/sbhawani/cernbox/ProtonDeuteron/WFs/PionlessEFT/Pd_CorrctNorm2/", Kitty, 1, 800);
  for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  }
  Kitty.SetChannelWeight(0, 1.0 / 3.0);
  Kitty.SetChannelWeight(1, 2.0 / 3.0);
  Kitty.SetQ1Q2(1);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/cernbox/ProtonDeuteron/Outputs/CATSOutput/PionlessEFTCFs/PdCFs2/Quartet_CorrNormNew-L800-l0p5.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.1; RAD < 99; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, real(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}

void NeuteronDeuteronCF_fromWFModel() {
  const double SourceSize = 0.94;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);
  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
// Kitty.SetQ1Q2(1);
  //Kitty.SetGamow(true);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  ExternalWF = Init_nd_Sebastian("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/WFndFromSebastian/", Kitty, 400);
  for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  }
  Kitty.SetChannelWeight(0, 0.333);
  Kitty.SetChannelWeight(1, 0.666);
  Kitty.SetQ1Q2(0);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Output/QA_nd/LatestLatestLatestnd_pdCFs/TotalQD_nDCF_QDLatestLatestLatest400.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.1; RAD < 99; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}


void ProtonDeuteronCF_Coulomb_fromWFModel() {
  const double SourceSize = 0.94;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);
  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  //Kitty.SetGamow(true);
  printf("reached in coulomb -1\n");
  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  ExternalWF = Init_pd_Coulomb_Sebastian("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/WFFromSebastian/", Kitty, 1, 200);
  printf("reached in coulomb -2 \n");
  for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  }
  printf("WF ready\n");
  Kitty.SetChannelWeight(0, 0.33);
  Kitty.SetChannelWeight(1, 0.66);
  Kitty.SetQ1Q2(1);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Output/QA_pd/QA_Pd_Coulomb_Q_L200_0p2.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.1; RAD < 200; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}


void LabelDCAHisto(TH1F *Histo, TString Title, EColor color) {

  Histo->SetMarkerColor(color);
  Histo->SetLineColor(color);
  Histo->SetMarkerSize(0.8);
  Histo->SetMarkerStyle(20);
// Histo->SetTitle(Title.Data());
  Histo->GetXaxis()->SetTitle("k*(MeV/c)");
  Histo->GetYaxis()->SetTitle("#delta(k*)");
  Histo->GetXaxis()->SetTitleOffset(1.0);

}

void MakePlotforPhase(TString Name, TH1F *h1, TH1F *h2, TH1F *h3, int Cutoff, double photonmass, double xRangeLow, double xRangeHigh, double yRangeLow, double yRangeHigh ) {

  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName = Form("%s__PhaseShift.png", Name.Data()); //plots are output as .pdf If you prefer other formats simply change the ending
  //gPad->SetLogy();
  LabelDCAHisto(h1, "#Rgothic [#delta]" , kBlack);
  LabelDCAHisto(h2, "#Jgothic [#delta]", kRed);
  LabelDCAHisto(h3, "|#delta|", kGreen);

  //h1->Draw("EPZ");
  h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  h1->GetXaxis()->SetTitle("k* (MeV/c)");
  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetXaxis()->SetTitleOffset(0.75);
  h1->GetYaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleOffset(0.95);
  h1->GetYaxis()->SetTitle("#delta(k*)");

  h1->GetXaxis()->SetRangeUser(xRangeLow, xRangeHigh);
  h1->GetYaxis()->SetRangeUser(yRangeLow, yRangeHigh);

  TLegend *leg1 = new TLegend(0.15, 0.8, 0.7, 0.85);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.06);
  leg1->SetLineColor(0);
  leg1->SetNColumns(3);
  leg1->AddEntry(h1, " #Rgothic [#delta]");
  leg1->AddEntry(h2, " #Jgothic [#delta]");
  leg1->AddEntry(h3, " |#delta|");
  h1->Draw("PE");
  leg1->Draw("PEZSAME");

  h2->DrawCopy("PEZSAME");
  h3->DrawCopy("PEZSAME");

  Plot->Print(pdfName);
  delete Plot;
}

void ProtonDeuteron_PhaseShift_fromWFModel() {
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
//  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  DLM_Histo<float> OutputRelPhase_D;
  DLM_Histo<float> OutputImPhase_D;
  DLM_Histo<float> OutputAbsPhase_D;
  DLM_Histo<float> OutputRelPhase_Q;
  DLM_Histo<float> OutputImPhase_Q;
  DLM_Histo<float> OutputAbsPhase_Q;
  printf("going to load WF \n");
  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  ExternalWF = Init_pd_PhaseShift_Sebastian("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/PhaseShiftFromSebastian/", OutputRelPhase_D,
               OutputImPhase_D, OutputAbsPhase_D, OutputRelPhase_Q,
               OutputImPhase_Q, OutputAbsPhase_Q, 0, 200);
  printf(" loaded WF \n");
  //for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
  //Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  //}
  printf("WF ready\n");

  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Output/QA_pd/QA_Pd_PhaseShift_Q_L200_0p2.root"), "recreate");
  printf("File Created\n");

//Histos for Signal extraction
  TGraph gWF_Re_D;
  TGraph gWF_Im_D;
  TGraph gWF_Abs_D;

  TH1F *hWF_Re_D;
  TH1F *hWF_Im_D;
  TH1F *hWF_Abs_D;


  TH1F *hWF_Re_Q;
  TH1F *hWF_Im_Q;
  TH1F *hWF_Abs_Q;

  float binsx_k[27] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0};

  hWF_Re_D = new TH1F("hWF_Re_D", "hWF_Re_D", 26, binsx_k);
  hWF_Im_D = new TH1F("hWF_Im_D", "hWF_Im_D", 26, binsx_k);
  hWF_Abs_D = new TH1F("hWF_Abs_D", "hWF_Abs_D", 26, binsx_k);


  hWF_Re_Q = new TH1F("hWF_Re_Q", "hWF_Re_Q", 26, binsx_k);
  hWF_Im_Q = new TH1F("hWF_Im_Q", "hWF_Im_Q", 26, binsx_k);
  hWF_Abs_Q = new TH1F("hWF_Abs_Q", "hWF_Abs_Q", 26, binsx_k);


  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];

  unsigned COUNTER = 0;
  gWF_Re_D.SetName(TString::Format("gWF_Re_D"));
  gWF_Im_D.SetName(TString::Format("gWF_Im_D"));
  gWF_Abs_D.SetName(TString::Format("gWF_Abs_D"));
  double kval = 0.0;
  int hkval = 0;
  for (double RAD = 1; RAD <= 26; RAD ++) {
    kval = RAD * 5.0;
    hkval = RAD;
    printf("RelDelta_D(%.2f) = %.2f ImDelta_D(%.2f) = %.2f AbsDelta_D(%.2f) = %.2f\n", kval, OutputRelPhase_D.Eval(&kval), kval, OutputImPhase_D.Eval(&kval), kval, OutputAbsPhase_D.Eval(&kval));
    gWF_Re_D.SetPoint(COUNTER, kval, OutputRelPhase_D.Eval(&kval));
    gWF_Im_D.SetPoint(COUNTER, kval, OutputImPhase_D.Eval(&kval));
    gWF_Abs_D.SetPoint(COUNTER, kval, OutputAbsPhase_D.Eval(&kval));

    hWF_Re_D->SetBinContent(hkval, OutputRelPhase_D.Eval(&kval));
    hWF_Re_D->SetBinError(hkval, 0.01);
    hWF_Im_D->SetBinContent(hkval, OutputImPhase_D.Eval(&kval));
    hWF_Im_D->SetBinError(hkval, 0.01);
    hWF_Abs_D->SetBinContent(hkval, OutputAbsPhase_D.Eval(&kval));
    hWF_Abs_D->SetBinError(hkval, 0.01);

    hWF_Re_Q->SetBinContent(hkval, OutputRelPhase_Q.Eval(&kval));
    hWF_Re_Q->SetBinError(hkval, 0.01);
    hWF_Im_Q->SetBinContent(hkval, OutputImPhase_Q.Eval(&kval));
    hWF_Im_Q->SetBinError(hkval, 0.01);
    hWF_Abs_Q->SetBinContent(hkval, OutputAbsPhase_Q.Eval(&kval));
    hWF_Abs_Q->SetBinError(hkval, 0.01);

    COUNTER++;
  }


  gWF_Re_D.Write();
  gWF_Im_D.Write();
  gWF_Abs_D.Write();

  hWF_Re_D->Write();
  hWF_Im_D->Write();
  hWF_Abs_D->Write();

  hWF_Re_Q->Write();
  hWF_Im_Q->Write();
  hWF_Abs_Q->Write();
  MakePlotforPhase("Doublet", hWF_Re_D, hWF_Im_D, hWF_Abs_D, 200, 0.5, 0.0, 130.0, -150, 150 );
  MakePlotforPhase("Quartet", hWF_Re_Q, hWF_Im_Q, hWF_Abs_Q, 200, 0.5, 0.0, 130.0, -150, 150 );
  delete OutputFile;

  //delete hWF_Re_D;
  /*  delete hWF_Im_D;
    delete hWF_Abs_D;

    delete hWF_Re_Q;
    delete hWF_Im_Q;
    delete hWF_Abs_Q;*/
}

void ProtonProtonCF_Sebastian_fromWFModel() {
  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);

  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;
  // ExternalWF = Init_pp_Sebastian("/home/sbhawani/Desktop/PdCFTheory/ppCalculationFromSebastian/ppPionlessEFT13092020/", Kitty, 0, 800);
  ExternalWF = Init_ppPionLess_Sebastian("/home/sbhawani/Desktop/PdCFTheory/ppCalculationFromSebastian/ppPionlessEFT13092020/", Kitty, 0, 800);


  //for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
  Kitty.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
  // }
  printf("WF ready\n");

  Kitty.SetSpin(0, 0);
  Kitty.SetSpin(1, 1);
  Kitty.SetChannelWeight(0, 1.0 / 4.0);
  Kitty.SetChannelWeight(1, 3.0 / 4.0);

  Kitty.SetQ1Q2(0);
  Kitty.SetPdgId(2212, 2212);

  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/Output/QA_pp_Sebastian_pionless_17092020_Coulomb.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.001; RAD < 300; RAD += 0.0005) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}

void ProtonProtonAV18CF_Sebastian_fromWFModel() {
  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);

  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;

  ExternalWF = Init_ppAV18_Sebastian("/home/sbhawani/Desktop/PdCFTheory/NormCorrectedWF/ppAV18Sebastian/MathematicaCalculation/", Kitty, 0, 800);

  //for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
  Kitty.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
// }
  printf("WF ready\n");

  Kitty.SetQ1Q2(1);
  Kitty.SetPdgId(2212, 2212);
  Kitty.SetChannelWeight(0, 1.0 / 4.0);
  Kitty.SetChannelWeight(1, 3.0 / 4.0);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/OutputCorrNorm/ppAV18CorrectedNorm17112020.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.001; RAD < 300; RAD += 0.0005) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}

void ProtonProtonPionLessCF_Sebastian_fromWFModel() {
  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 52;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);

  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;

  //ExternalWF = Init_ppPionLessNew_Sebastian("/home/sbhawani/Desktop/PdCFTheory/NormCorrectedWF/ppPionLessSebastian/MathematicaCalculation/", Kitty, 1, 800);
  ExternalWF = Init_ppPionLessMomSpace_Sebastian("/home/sbhawani/Desktop/PdCFTheory/NormCorrectedWF/ppPionLessSebastianMomSpace/MathematicaCalculation/", Kitty, 1, 800);


  //for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
  Kitty.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
// }
  printf("WF ready\n");

  Kitty.SetQ1Q2(1);
  Kitty.SetPdgId(2212, 2212);
  Kitty.SetChannelWeight(0, 1.0 / 4.0);
  Kitty.SetChannelWeight(1, 3.0 / 4.0);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/OutputCorrNorm/ppPionless1p0l0p25MomSpaceUpdated.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.001; RAD < 300; RAD += 0.0005) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}

void NeuteronProtonAV18CF_Sebastian_fromWFModel() {
  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);

  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;

  ExternalWF = Init_npAV18_Sebastian("/home/sbhawani/Desktop/PdCFTheory/ppCalculationFromSebastian/AV18FromPionlessEFT/", Kitty, 0, 800);

  for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
  }
  printf("WF ready\n");

  Kitty.SetQ1Q2(0);
  Kitty.SetPdgId(2212, 2212);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Output/QA_pd/QA_npAV18_Sebastian_CorrectedWithoutFirstBin.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.01; RAD < 100; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}

void ProtonDeuteronCF_Coulomb_CATS() {
  const double SourceSize = 1.0;
  const double kMin = 0;
  const double kMax = 400;
  const unsigned NumMomBins = 100;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);
  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
  Kitty.SetPdgId(2212, 1000010020);
  Kitty.SetRedMass( Mass_p * Mass_d / (Mass_p + Mass_d));
  // Kitty.SetRedMass( Mass_p/2.0);
  Kitty.SetQ1Q2(1);
  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0, 0);
  Kitty.SetSpin(0, 0);
  Kitty.SetChannelWeight(0, 1.);
  Kitty.KillTheCat();
  DLM_Ck Ck_pd(Kitty.GetNumSourcePars(), 0, Kitty);
  Ck_pd.Update();


  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/Systematics/CoulombCheckLaura/Ck_pdCoulomb_1p0.root"), "recreate");
  printf("File Created\n");
  TGraphErrors gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPointError(uBin, 20.0, 0.0);
  }
  gKitty.Write();
  delete OutputFile;
  delete cPars;
}

void ProtonDeuteronCF_CoulombvsRadii_CATS() {
  const double SourceSize[20] = {1.073,1.2,1.5,1.8,2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0,12.0,15.0,20.0,40.0,70.0};
  const double kMin = 0;
  const double kMax = 400;
  const unsigned NumMomBins = 100;

  TGraph *gR1 = new TGraph();
   TGraph *gR2 = new TGraph();
   TGraph *gR3 = new TGraph();
   TGraph *gR4 = new TGraph();
   TGraph *gR5 = new TGraph();
   TGraph *gR6 = new TGraph();
   TGraph *gR7 = new TGraph();
   TGraph *gR8 = new TGraph();
   TGraph *gR9 = new TGraph();
   TGraph *gR10 = new TGraph();
   TGraph *gR11 = new TGraph();
   TGraph *gR12 = new TGraph();
   TGraph *gR13 = new TGraph();
   TGraph *gR14 = new TGraph();
   TGraph *gR15 = new TGraph();
   TGraph *gR16 = new TGraph();
   TGraph *gR17 = new TGraph();

   TH1F *hR1 = new TH1F("hR1", "hR1", NumMomBins, 0.0,
                        kMax);
   TH1F *hR2 = new TH1F("hR2", "hR2", NumMomBins, 0.0,
                        kMax);
   TH1F *hR3 = new TH1F("hR3", "hR3", NumMomBins, 0.0,
                        kMax);
   TH1F *hR4 = new TH1F("hR4", "hR4", NumMomBins, 0.0,
                        kMax);
   TH1F *hR5 = new TH1F("hR5", "hR5", NumMomBins, 0.0,
                        kMax);
   TH1F *hR6 = new TH1F("hR6", "hR6", NumMomBins, 0.0,
                        kMax);
   TH1F *hR7 = new TH1F("hR7", "hR7", NumMomBins, 0.0,
                        kMax);
   TH1F *hR8 = new TH1F("hR8", "hR8", NumMomBins, 0.0,
                        kMax);
   TH1F *hR9 = new TH1F("hR9", "hR9", NumMomBins, 0.0,
                        kMax);
   TH1F *hR10 = new TH1F("hR10", "hR10", NumMomBins, 0.0,
                         kMax);
   TH1F *hR11 = new TH1F("hR11", "hR11", NumMomBins, 0.0,
                         kMax);
   TH1F *hR12 = new TH1F("hR12", "hR12", NumMomBins, 0.0,
                         kMax);

   TH1F *hR13 = new TH1F("hR13", "hR13", NumMomBins, 0.0,
                         kMax);
   TH1F *hR14 = new TH1F("hR14", "hR14", NumMomBins, 0.0,
                         kMax);
   TH1F *hR15 = new TH1F("hR15", "hR15", NumMomBins, 0.0,
                         kMax);
   TH1F *hR16 = new TH1F("hR16", "hR16", NumMomBins, 0.0,
                         kMax);
   TH1F *hR17 = new TH1F("hR17", "hR17", NumMomBins, 0.0,
                         kMax);
   gR1->SetName("gR1");
   gR2->SetName("gR2");
   gR3->SetName("gR3");
   gR4->SetName("gR4");
   gR5->SetName("gR5");
   gR6->SetName("gR6");
   gR7->SetName("gR7");
   gR8->SetName("gR8");
   gR9->SetName("gR9");
   gR10->SetName("gR10");
   gR11->SetName("gR11");
   gR12->SetName("gR12");
   gR13->SetName("gR13");
   gR14->SetName("gR14");
   gR15->SetName("gR15");
   gR16->SetName("gR16");
   gR17->SetName("gR17");

   float CkR1 = 0.0;
   float CkR2= 0.0;
   float CkR3 = 0.0;
   float CkR4 = 0.0;
   float CkR5 = 0.0;
   float CkR6 = 0.0;
   float CkR7 = 0.0;
   float CkR8 = 0.0;
   float CkR9 = 0.0;
   float CkR10 = 0.0;
   float CkR11 = 0.0;
   float CkR12 = 0.0;
   float CkR13 = 0.0;
   float CkR14 = 0.0;
   float CkR15 = 0.0;
   float CkR16 = 0.0;
   float CkR17 = 0.0;

   float R1 = 1.0;
   float R2= 1.2;
   float R3 = 1.5;
   float R4 = 1.8;
   float R5 = 2.0;
   float R6 = 2.5;
   float R7 = 3.0;
   float R8 = 4.0;
   float R9 = 5.0;
   float R10 = 6.0;
   float R11 = 8.0;
   float R12 = 10.0;
   float R13 = 12.0;
   float R14 = 15.0;
   float R15 = 20.0;
   float R16 = 40.0;
   float R17 = 70.0;

  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  for(int i = 1;i<18;i++){
    cPars->SetParameter(0, SourceSize[i-1]);
    CATS Kitty;
    Kitty.SetAnaSource(GaussSource, *cPars);
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    Kitty.SetAnaSource(GaussSource, *cPars);
    Kitty.SetAnaSource(0, SourceSize[i-1]);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetPdgId(2212, 1000010020);
    Kitty.SetRedMass( Mass_p * Mass_d / (Mass_p + Mass_d));
    // Kitty.SetRedMass( Mass_p/2.0);
    Kitty.SetQ1Q2(1);
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);
    Kitty.KillTheCat();
    DLM_Ck Ck_pd(Kitty.GetNumSourcePars(), 0, Kitty);
    Ck_pd.Update();


    TFile* OutputFile = new TFile(
      TString::Format("/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/Systematics/CoulombCheckLaura/Quartet_Ck_pdCoulomb_1p0.root"), "recreate");
    printf("File Created\n");
    for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
      if(i==1){
        gR1->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        hR1->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==2){
      gR2->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR2->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==3){
      gR3->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR3->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==4){
      gR4->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR4->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==5){
      gR5->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR5->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==6){
      gR6->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR6->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==7){
      gR7->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR7->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i ==8){
      gR8->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR8->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==9){
      gR9->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR9->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==10){
      gR10->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR10->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==11){
      gR11->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR11->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==12){
      gR12->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR12->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==13){
      gR13->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR13->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i==14){
      gR14->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR14->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i == 15){
      gR15->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR15->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i== 16){
      gR16->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR16->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }else if(i== 17){
      gR17->SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
      hR17->SetBinContent(uBin+1, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    }
    }
  }

   TString OutputFileName = "";
       OutputFileName =
           OutputFileName
               + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/CoulombVsRadPlots.root";

   TFile *OutputFile = new TFile(OutputFileName, "recreate");
   printf("File Created\n");
   TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
   TString pdfName = "";
   pdfName =pdfName+ Form("%s CoulombVsRadPlots.pdf","/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
    //plots are output as .pdf If you prefer other formats simply change the ending
   printf("reached -2!");
   //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
     gR1->GetXaxis()->SetTitle("k (MeV/c)");
     gR1->GetXaxis()->SetTitleSize(0.045);
     gR1->GetXaxis()->SetTitleOffset(0.5);
     gR1->GetYaxis()->SetTitleSize(0.045);
     gR1->GetYaxis()->SetTitleOffset(0.65);
     gR1->GetYaxis()->SetTitle("#it{C}(k)");

     gR1->GetXaxis()->SetRangeUser(0, 400);
     gR1->GetYaxis()->SetRangeUser(0, 50);


     hR1->GetXaxis()->SetTitle("k (MeV/c)");
     hR1->GetXaxis()->SetTitleSize(0.045);
     hR1->GetXaxis()->SetTitleOffset(0.5);
     hR1->GetYaxis()->SetTitleSize(0.045);
     hR1->GetYaxis()->SetTitleOffset(0.65);
     hR1->GetYaxis()->SetTitle("#it{C}(k)");

     hR1->GetXaxis()->SetRangeUser(0, 400);
     hR1->GetYaxis()->SetRangeUser(0, 50);

     gR1->Draw("ACP");
     gR1->Write();
     gR2->Write();
     gR3->Write();
     gR4->Write();
     gR5->Write();
     gR6->Write();
     gR7->Write();
     gR8->Write();
     gR9->Write();
     gR10->Write();
     gR11->Write();
     gR12->Write();
     gR13->Write();
     gR14->Write();
     gR15->Write();
     gR16->Write();
     gR17->Write();


     hR1->Write();
     hR2->Write();
     hR3->Write();
     hR4->Write();
     hR5->Write();;
     hR6->Write();
     hR7->Write();
     hR8->Write();
     hR9->Write();
     hR10->Write();
     hR11->Write();
     hR12->Write();
     hR13->Write();
     hR14->Write();
     hR15->Write();
     hR16->Write();
     hR17->Write();
     printf("reached -3!");
     Plot->Print(pdfName);

     delete gR1;
     delete gR2;
     delete gR3;
     delete gR4;
     delete gR5;
     delete gR6;
     delete gR7;
     delete gR8;
     delete gR9;
     delete gR10;
     delete gR11;
     delete gR12;
     delete gR13;
     delete gR14;
     delete gR15;
     delete gR16;
     delete gR17;

     delete hR1;
     delete hR2;
     delete hR3;
     delete hR4;
     delete hR5;
     delete hR6;
     delete hR7;
     delete hR8;
     delete hR9;
     delete hR10;
     delete hR11;
     delete hR12;
     delete hR13;
     delete hR14;
     delete hR15;
     delete hR16;
     delete hR17;
  delete OutputFile;
  delete cPars;
}

void CATS_GaussSource1(CATS & Kitty, const double & SourceSize) {
  //object containing the source or potential parameters. Arguments:
  //(tSource or tPotential, NumberOfParameters, ReadyForMulti-threading)
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  //set up the parameters
  cPars.SetParameter(0, SourceSize);
  //set the source to a Gaussian function, using the parameters above. These parameters are copied into CATS,
  //thus it is okay if cPars is deleted afterwards.
  Kitty.SetAnaSource(GaussSource, cPars);
  //this is important, since if its `false` it is assumed that the source will be sampled from a transport model, and the Gaussian function will not be used
  Kitty.SetUseAnalyticSource(true);
  //if true, the source is automatically renormalized in the range 0-64 fm. Nice to dummy proof the source, but problematic for sources with large tails
  //for the Gaussian example above, both options should be completely identical
  Kitty.SetAutoNormSource(false);
}

void CATS_pp_Basic1(CATS & Kitty, const unsigned & NumMomBins, const double & kMin, const double & kMax) {
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetExcludeFailedBins(false);
  //the charge of the pair given in q1 * q2, where q1 and q2 are the charges of the individual particles
  Kitty.SetQ1Q2(1);
  //pid of the particle. Important only to set properly the quantum statistics (on for identical guys)
  //one could also force the QS by: Kitty.SetQuantumStatistics(true/false);
  Kitty.SetPdgId(2212, 2212);
  //reduced mass of the pair
  Kitty.SetRedMass( 0.5 * Mass_p );
}
//initialize the interaction for pp, using the AV18 potential.
//by default we have only s-waves included, but one can include the p and d waves if desired
void CATS_pp_AV181(CATS & Kitty, const bool & pwaves, const bool & dwaves) {
  //the 4 channels for pp are:
  //s=0: 1S0 + 3D1
  //s=1: 3P0
  //s=1: 3P1
  //s=1: 3P2
  //note that for s=0 the p-waves are Pauli blocked, for s=1 these are the s and d waves
  if (pwaves) {
    Kitty.SetNumChannels(4);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 2);
    Kitty.SetNumPW(2, 2);
    Kitty.SetNumPW(3, 2);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetSpin(2, 1);
    Kitty.SetSpin(3, 1);
    Kitty.SetChannelWeight(0, 3. / 12.);
    Kitty.SetChannelWeight(1, 1. / 12.);
    Kitty.SetChannelWeight(2, 3. / 12.);
    Kitty.SetChannelWeight(3, 5. / 12.);
  } else {
    //important: even with the p-waves switched off, physics wise the spin 1 state still exists and
    //the p-waves are there, just in there asymptotic state (free wave). To include this in the computation,
    //CATS still needs a second channel, even if it is `empty`!
    Kitty.SetNumChannels(2);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1.0);
    Kitty.SetChannelWeight(1, 0.0);
  }

  //to set up the strong interaction, one can use the predefined functions available in DLM_Potentials.h
  //the main idea is to always pass the function fDlmPot but with different input parameters, based on which the interaction is set up automatically
  //this works not only for pp, but many other systems are included.
  //The fDlmPot is only used as an interface to easily use different potentials, that are hard coded in DLM_Potentials.h
  //Feel free to expand the data base of this file if you think others will benefit from it. For further details contact Dimi
  //The input parameters are by default 9 and are defined as follows:
  //0: potential flag (defines which potential to use, see the enumerators in DLM_Potentials.h for more info)
  //1: a second flag, that can be used if needed (depending on the definition of the potential)
  //2: total isospin
  //3: 2 x isospin of particle 1
  //4: 2 x isospin of particle 2
  //5: total spin
  //6: l quantum number
  //7: j quantum number
  CATSparameters cPars_pp_1S0(CATSparameters::tPotential, 8, true);
  cPars_pp_1S0.SetParameter(0, NN_ReidV8); //choose the AV18
  cPars_pp_1S0.SetParameter(1, v18_Coupled3P2); //default option, which takes the 3P2 channel from a coupled-channel computation, but in CATS only the first diagonal potential elements is used
  cPars_pp_1S0.SetParameter(2, 1);
  cPars_pp_1S0.SetParameter(3, 1);
  cPars_pp_1S0.SetParameter(4, 1);
  cPars_pp_1S0.SetParameter(5, 0);
  cPars_pp_1S0.SetParameter(6, 0);
  cPars_pp_1S0.SetParameter(7, 0);

  //copy all settings from cPars_pp_1S0, and just change quantum numbers s,l,j
  CATSparameters cPars_pp_3P0(cPars_pp_1S0);
  cPars_pp_3P0.SetParameter(5, 1);
  cPars_pp_3P0.SetParameter(6, 1);
  cPars_pp_3P0.SetParameter(7, 0);

  CATSparameters cPars_pp_3P1(cPars_pp_1S0);
  cPars_pp_3P1.SetParameter(5, 1);
  cPars_pp_3P1.SetParameter(6, 1);
  cPars_pp_3P1.SetParameter(7, 1);

  CATSparameters cPars_pp_3P2(cPars_pp_1S0);
  cPars_pp_3P2.SetParameter(5, 1);
  cPars_pp_3P2.SetParameter(6, 1);
  cPars_pp_3P2.SetParameter(7, 2);

  CATSparameters cPars_pp_1D2(cPars_pp_1S0);
  cPars_pp_1D2.SetParameter(5, 0);
  cPars_pp_1D2.SetParameter(6, 2);
  cPars_pp_1D2.SetParameter(7, 2);

  //plug in the strong potential for each channel and partial wave
  //the arguments are: #WhichChannel,#WhichPartialWave,#PotentialFunction,#PotentialParameters
  Kitty.SetShortRangePotential(0, 0, fDlmPot, cPars_pp_1S0);
  if (pwaves) {
    Kitty.SetShortRangePotential(1, 1, fDlmPot, cPars_pp_3P0);
    Kitty.SetShortRangePotential(2, 1, fDlmPot, cPars_pp_3P1);
    Kitty.SetShortRangePotential(3, 1, fDlmPot, cPars_pp_3P2);
  }
  if (dwaves) {
    Kitty.SetShortRangePotential(0, 2, fDlmPot, cPars_pp_1D2);
  }
  //if later on you would like to switch some contribution off, this can be done with:
  //Kitty.RemoveShortRangePotential(#WhichChannel,#WhichPartialWave);
}


void CATS_pp_AV181_SwaveOnly(CATS & Kitty, const bool & pwaves, const bool & dwaves) {
  //the 4 channels for pp are:
  //s=0: 1S0 + 3D1
  //s=1: 3P0
  //s=1: 3P1
  //s=1: 3P2
  //note that for s=0 the p-waves are Pauli blocked, for s=1 these are the s and d waves
  if (pwaves) {
    Kitty.SetNumChannels(4);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 2);
    Kitty.SetNumPW(2, 2);
    Kitty.SetNumPW(3, 2);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetSpin(2, 1);
    Kitty.SetSpin(3, 1);
    Kitty.SetChannelWeight(0, 3. / 12.);
    Kitty.SetChannelWeight(1, 1. / 12.);
    Kitty.SetChannelWeight(2, 3. / 12.);
    Kitty.SetChannelWeight(3, 5. / 12.);
  } else {
    //important: even with the p-waves switched off, physics wise the spin 1 state still exists and
    //the p-waves are there, just in there asymptotic state (free wave). To include this in the computation,
    //CATS still needs a second channel, even if it is `empty`!
    Kitty.SetNumChannels(2);
    if (dwaves) Kitty.SetNumPW(0, 3);
    else Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1.0 / 4.0);
    Kitty.SetChannelWeight(1, 3.0 / 4.0);
  }

  //to set up the strong interaction, one can use the predefined functions available in DLM_Potentials.h
  //the main idea is to always pass the function fDlmPot but with different input parameters, based on which the interaction is set up automatically
  //this works not only for pp, but many other systems are included.
  //The fDlmPot is only used as an interface to easily use different potentials, that are hard coded in DLM_Potentials.h
  //Feel free to expand the data base of this file if you think others will benefit from it. For further details contact Dimi
  //The input parameters are by default 9 and are defined as follows:
  //0: potential flag (defines which potential to use, see the enumerators in DLM_Potentials.h for more info)
  //1: a second flag, that can be used if needed (depending on the definition of the potential)
  //2: total isospin
  //3: 2 x isospin of particle 1
  //4: 2 x isospin of particle 2
  //5: total spin
  //6: l quantum number
  //7: j quantum number
  CATSparameters cPars_pp_1S0(CATSparameters::tPotential, 8, true);
  cPars_pp_1S0.SetParameter(0, NN_ReidV8); //choose the AV18
  cPars_pp_1S0.SetParameter(1, v18_Coupled3P2); //default option, which takes the 3P2 channel from a coupled-channel computation, but in CATS only the first diagonal potential elements is used
  cPars_pp_1S0.SetParameter(2, 1);
  cPars_pp_1S0.SetParameter(3, 1);
  cPars_pp_1S0.SetParameter(4, 1);
  cPars_pp_1S0.SetParameter(5, 0);
  cPars_pp_1S0.SetParameter(6, 0);
  cPars_pp_1S0.SetParameter(7, 0);

  //copy all settings from cPars_pp_1S0, and just change quantum numbers s,l,j
  CATSparameters cPars_pp_3P0(cPars_pp_1S0);
  cPars_pp_3P0.SetParameter(5, 1);
  cPars_pp_3P0.SetParameter(6, 1);
  cPars_pp_3P0.SetParameter(7, 0);

  CATSparameters cPars_pp_3P1(cPars_pp_1S0);
  cPars_pp_3P1.SetParameter(5, 1);
  cPars_pp_3P1.SetParameter(6, 1);
  cPars_pp_3P1.SetParameter(7, 1);

  CATSparameters cPars_pp_3P2(cPars_pp_1S0);
  cPars_pp_3P2.SetParameter(5, 1);
  cPars_pp_3P2.SetParameter(6, 1);
  cPars_pp_3P2.SetParameter(7, 2);

  CATSparameters cPars_pp_1D2(cPars_pp_1S0);
  cPars_pp_1D2.SetParameter(5, 0);
  cPars_pp_1D2.SetParameter(6, 2);
  cPars_pp_1D2.SetParameter(7, 2);

  //plug in the strong potential for each channel and partial wave
  //the arguments are: #WhichChannel,#WhichPartialWave,#PotentialFunction,#PotentialParameters
  Kitty.SetShortRangePotential(0, 0, fDlmPot, cPars_pp_1S0);
  if (pwaves) {
    Kitty.SetShortRangePotential(1, 1, fDlmPot, cPars_pp_3P0);
    Kitty.SetShortRangePotential(2, 1, fDlmPot, cPars_pp_3P1);
    Kitty.SetShortRangePotential(3, 1, fDlmPot, cPars_pp_3P2);
  }
  if (dwaves) {
    Kitty.SetShortRangePotential(0, 2, fDlmPot, cPars_pp_1D2);
  }
  //if later on you would like to switch some contribution off, this can be done with:
  //Kitty.RemoveShortRangePotential(#WhichChannel,#WhichPartialWave);
}

void Ck_pp_CATS_CoulombOnly() {

  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/Output/CF_PP_CoulombOnly.root"), "recreate");
  printf("File Created\n");

  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;


  CATS Kitty_pp;

  Kitty_pp.SetQ1Q2(1);
  Kitty_pp.SetPdgId(2212, 2212);
// Kitty_pp.SetChannelWeight(0, 1.0 / 4.0);
  //Kitty_pp.SetChannelWeight(1, 3.0 / 4.0);
  //initialize a Gaussian source with a source size of 1.2
  CATS_GaussSource1(Kitty_pp, SourceSize);
  //you can change the parameters of the source at any time using the function:
  //Kitty_pp.SetAnaSource(#WhichParameter,#Value);
  //set up the CATS object for pp (Coulomb and QS) with 100 bins in the range 0-400 MeV

  //Kitty_pp.SetAnaSource(GaussSource, *cPars);
  Kitty_pp.SetMomBins(NumMomBins, kMin, kMax);
  //Kitty_pp.SetAnaSource(GaussSource, *cPars);
  Kitty_pp.SetAnaSource(0, SourceSize);
  Kitty_pp.SetUseAnalyticSource(true);
  Kitty_pp.SetMomentumDependentSource(false);
  Kitty_pp.SetThetaDependentSource(false);
  Kitty_pp.SetExcludeFailedBins(false);
  Kitty_pp.SetPdgId(2212, 1000010020);
  Kitty_pp.SetRedMass( Mass_p * Mass_d / (Mass_p + Mass_d));
  // Kitty.SetRedMass( Mass_p/2.0);
  Kitty_pp.SetQ1Q2(1);
  Kitty_pp.SetNumChannels(1);
  Kitty_pp.SetNumPW(0, 0);
  Kitty_pp.SetSpin(0, 0);
  Kitty_pp.SetChannelWeight(0, 1.);
  //set up the cats interaction including only s-waves
  //CATS_pp_AV181_SwaveOnly(Kitty_pp, false, false);
  //compute the correlation function
  Kitty_pp.KillTheCat();
  //set up a `histogram` for the pp interaction, based on the above CATS object
  //the arguments are the number of source/potential parameters to be controlled by DLM_Ck
  DLM_Ck Ck_pp(Kitty_pp.GetNumSourcePars(), 0, Kitty_pp);
  //btw, now you can change the source parameters also by using:
  //Ck_pp.SetSourcePar(#WhichParameter,#Value);
  //this function reads the CATS objects and fills the bins of the DLM_Ck object
  Ck_pp.Update();

  //save the correlation (theoretical at the moment) in a file
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
  }
  gKitty.Write();
  delete OutputFile;
}

void Ck_pp_CATS() {

  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/Output/CF_AV18_pp_CATSQ1Q2_1.root"), "recreate");
  printf("File Created\n");

  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;


  CATS Kitty_pp;

  Kitty_pp.SetQ1Q2(1);
  Kitty_pp.SetPdgId(2212, 2212);
// Kitty_pp.SetChannelWeight(0, 1.0 / 4.0);
  //Kitty_pp.SetChannelWeight(1, 3.0 / 4.0);
  //initialize a Gaussian source with a source size of 1.2
  CATS_GaussSource1(Kitty_pp, SourceSize);
  //you can change the parameters of the source at any time using the function:
  //Kitty_pp.SetAnaSource(#WhichParameter,#Value);
  //set up the CATS object for pp (Coulomb and QS) with 100 bins in the range 0-400 MeV
  CATS_pp_Basic1(Kitty_pp, NumMomBins, kMin, kMax);
  //set up the cats interaction including only s-waves
  CATS_pp_AV181_SwaveOnly(Kitty_pp, false, false);
  //compute the correlation function
  Kitty_pp.KillTheCat();
  //set up a `histogram` for the pp interaction, based on the above CATS object
  //the arguments are the number of source/potential parameters to be controlled by DLM_Ck
  DLM_Ck Ck_pp(Kitty_pp.GetNumSourcePars(), 0, Kitty_pp);
  //btw, now you can change the source parameters also by using:
  //Ck_pp.SetSourcePar(#WhichParameter,#Value);
  //this function reads the CATS objects and fills the bins of the DLM_Ck object
  Ck_pp.Update();

  //save the correlation (theoretical at the moment) in a file
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
  }
  gKitty.Write();

  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  //TGraph* gWF_Im_D = new TGraph [NumMomBins];
  //TGraph* gWF_Im_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  double k = 0;
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    //let's take output for sebastian
    string k_str = to_string(Kitty_pp.GetMomentum(uBin));
    k_str.erase ( k_str.find_last_not_of('0') + 1, std::string::npos );
    string location = "/home/sbhawani/Desktop/AV18_CATS_WF/";
    size_t file_size = 1024;
    ofstream out_file(location + "WF-pp-k" + k_str + "-AV18-CATS.dat");
    out_file.precision (10);
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty_pp.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty_pp.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty_pp.GetMomentum(uBin)));

    for (double RAD = 0.01; RAD < 300; RAD += 0.5) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty_pp.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty_pp.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty_pp.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      out_file << RAD << "       " << real(Kitty_pp.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)) << "       " << imag(Kitty_pp.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)) << endl;
      //out_file.close();
      COUNTER++;
    }
    out_file.close();
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  /*
      //the same for s and p-waves
      CATS Kitty_pp_s;
      CATS_GaussSource1(Kitty_pp_s, SourceSize);
      CATS_pp_Basic1(Kitty_pp_s, 26, 0, 130);
      CATS_pp_AV181(Kitty_pp_s, false, false);
      Kitty_pp_s.KillTheCat();
      DLM_Ck Ck_pp_s(Kitty_pp_s.GetNumSourcePars(), 0, Kitty_pp_s);
      Ck_pp_s.Update();

     //the same for s and p-waves
      CATS Kitty_pp_sp;
      CATS_GaussSource1(Kitty_pp_sp, SourceSize);
      CATS_pp_Basic1(Kitty_pp_sp, 26, 0, 130);
      CATS_pp_AV181(Kitty_pp_sp, true, false);
      Kitty_pp_sp.KillTheCat();
      DLM_Ck Ck_pp_sp(Kitty_pp_sp.GetNumSourcePars(), 0, Kitty_pp_sp);
      Ck_pp_sp.Update();


      //and now the same including d-waves
      CATS Kitty_pp_spd;
      CATS_GaussSource1(Kitty_pp_spd, SourceSize);
      CATS_pp_Basic1(Kitty_pp_spd, 26, 0, 130);
      CATS_pp_AV181(Kitty_pp_spd, true, true);

      Kitty_pp_spd.KillTheCat();
      DLM_Ck Ck_pp_spd(Kitty_pp_spd.GetNumSourcePars(), 0, Kitty_pp_spd);

      Ck_pp_spd.Update();
  */

  delete OutputFile;
}

void ProtonProtonAV18_CATS_WF() {
  const double SourceSize = 1.2;
  const double kMin = 0;
  const double kMax = 130;
  const unsigned NumMomBins = 26;
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, SourceSize);

  CATS Kitty;
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetMomBins(NumMomBins, kMin, kMax);
  Kitty.SetAnaSource(GaussSource, *cPars);
  Kitty.SetAnaSource(0, SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetExcludeFailedBins(false);

  DLM_Histo<complex<double>>*** ExternalWF = NULL;

  ExternalWF = Init_AV18WF_CATS_Sebastian("/home/sbhawani/Desktop/", Kitty, 0, 800);

  //for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
  Kitty.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
// }
  printf("WF ready\n");

  Kitty.SetQ1Q2(1);
  Kitty.SetPdgId(2212, 2212);
  Kitty.SetChannelWeight(0, 1.0 / 4.0);
  Kitty.SetChannelWeight(1, 3.0 / 4.0);
  Kitty.KillTheCat();
  TFile* OutputFile = new TFile(
    TString::Format("/home/sbhawani/Desktop/AV18_CATS_WF/Output/CF_USING_WF_FROM_CATS_AV18_pp.root"), "recreate");
  printf("File Created\n");
  TGraph gKitty;
  gKitty.SetName(TString::Format("gKitty"));
  gKitty.Set(NumMomBins);
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    printf("C(%.2f) = %.2f\n", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gKitty.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
  }
  gKitty.Write();
  TGraph* gWF_Re_D = new TGraph [NumMomBins];
  TGraph* gWF_Re_Q = new TGraph [NumMomBins];
  TGraph* gWF_Re_R = new TGraph [NumMomBins];
  for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
    unsigned COUNTER = 0;
    gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f", Kitty.GetMomentum(uBin)));
    gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f", Kitty.GetMomentum(uBin)));
    for (double RAD = 0.001; RAD < 300; RAD += 0.05) {
      gWF_Re_D[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 0, 0, RAD, true)));
      gWF_Re_Q[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalRadialWaveFunction(uBin, 1, 0, RAD, true)));
      gWF_Re_R[uBin].SetPoint(COUNTER, RAD, abs(Kitty.EvalReferenceRadialWF(uBin, 0, RAD, true)));
      COUNTER++;
    }
    gWF_Re_D[uBin].Write();
    gWF_Re_Q[uBin].Write();
    gWF_Re_R[uBin].Write();
  }
  CleanUpWfHisto(Kitty, ExternalWF);
  delete [] gWF_Re_D;
  delete [] gWF_Re_Q;
  delete [] gWF_Re_R;
  delete OutputFile;
  delete cPars;
}
