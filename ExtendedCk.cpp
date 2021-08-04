
#include "Basics.h"
#include "ExtendedCk.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TSystem.h"

#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "Plots.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"
#include "DLM_Random.h"
#include "DLM_CppTools.h"

#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

///An example how to compute the pLambda correlation function using the Usmani potential and the Lednicky model
void Ck_pL_Ledni_Usmani() {
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 1.2);

    //the CATS object to model pLambda (using the Usmani potential)
    CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins, kMin, kMax);
    Kitty_pL.SetAnaSource(GaussSource, SOURCE_PARS);
    Kitty_pL.SetAnaSource(0, SOURCE_PARS.GetParameter(0));
    Kitty_pL.SetUseAnalyticSource(true);
    Kitty_pL.SetMomentumDependentSource(false);
    Kitty_pL.SetThetaDependentSource(false);
    Kitty_pL.SetExcludeFailedBins(false);
    Kitty_pL.SetQ1Q2(0);
    Kitty_pL.SetPdgId(2212, 3122);
    Kitty_pL.SetRedMass( (Mass_p * Mass_L) / (Mass_p + Mass_L) );
    //two spin channels (0 and 1)
    Kitty_pL.SetNumChannels(2);
    Kitty_pL.SetNumPW(0, 1);
    Kitty_pL.SetNumPW(1, 1);
    Kitty_pL.SetSpin(0, 0);
    Kitty_pL.SetSpin(1, 1);
    //the weights are based on the singlet/triplet configuration
    Kitty_pL.SetChannelWeight(0, 1. / 4.);
    Kitty_pL.SetChannelWeight(1, 3. / 4.);

    //the standard input to the predefined potentials using fDlmPot
    //this is the configuration for Usmani, where apart from the flag,
    //the the only relevant parameter is the spin (s)
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8] = {pL_UsmaniOli, 0, 0, 0, 0, 0, 0, 0};
    double PotPars3S1[8] = {pL_UsmaniOli, 0, 0, 0, 0, 1, 0, 0};
    //define a set of parameters for each channel
    CATSparameters cPotPars1S0(CATSparameters::tPotential, 8, true);
    cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1(CATSparameters::tPotential, 8, true);
    cPotPars3S1.SetParameters(PotPars3S1);

    //give to CATS the relevant information for the interaction in each channel
    Kitty_pL.SetShortRangePotential(0, 0, fDlmPot, cPotPars1S0);
    Kitty_pL.SetShortRangePotential(1, 0, fDlmPot, cPotPars3S1);
    //Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.KillTheCat();

    //The DLM_Ck object is in essence a histogram, that can be used to read a CATS object, or some function
    //Here we define the DLM_Ck for the Usmani potential, from the CATS object
    //number of source/potential pars to be "accessed"
    DLM_Ck Ck_Usmani(1, 0, Kitty_pL);
    //Ck_Usmani.SetSourcePar(0,2);
    Ck_Usmani.Update();

    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 4, NumMomBins, kMin, kMax, Lednicky_SingletTriplet);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 2.88);
    Ck_Lednicky.SetPotPar(1, 2.92);
    Ck_Lednicky.SetPotPar(2, 1.66);
    Ck_Lednicky.SetPotPar(3, 3.78);
    Ck_Lednicky.Update();

    RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_Usmani.root", "gCk_pL_Usmani", &Ck_Usmani);
    RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_Usmani.root", "gCk_pL_Lednicky", &Ck_Lednicky);
}

void Ck_pL_Ledni_QA1() {

    TString B2effplot = "TestForLD.pdf";
    TCanvas *binplot2 = new TCanvas("binplot2", "Graph Draw Options", 1000, 1000);
    binplot2->SetLeftMargin(0.1505);
    binplot2->SetRightMargin(0.035);


    TString OutputFileName = "../OutputFiles/Ck_pL_Ledni_QA1.root";
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 1.2);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 4, NumMomBins, kMin, kMax, Lednicky_SingletTriplet);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 2.88);
    Ck_Lednicky.SetPotPar(1, 2.92);
    Ck_Lednicky.SetPotPar(2, 1.66);
    Ck_Lednicky.SetPotPar(3, 3.78);
    Ck_Lednicky.Update();




    DLM_Ck Ck_Lednicky1(1, 6, NumMomBins, kMin, kMax, Lednicky_2channel);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky1.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky1.SetPotPar(0, 2.88);
    Ck_Lednicky1.SetPotPar(1, 2.92);
    Ck_Lednicky1.SetPotPar(2, 1.66);
    Ck_Lednicky1.SetPotPar(3, 3.78);
    Ck_Lednicky1.SetPotPar(4, 0.25);
    Ck_Lednicky1.SetPotPar(5, 0.75);
    Ck_Lednicky1.Update();


    const unsigned NumBins = Ck_Lednicky.GetNbins();
    TGraph graph;
    TGraph graph1;
    graph.SetName("GeneralLednicky");
    graph1.SetName("TwoChannelForLD");

    graph.Set(Ck_Lednicky.GetNbins());
    graph1.Set(Ck_Lednicky1.GetNbins());
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        graph.SetPoint(uBin, Ck_Lednicky.GetBinCenter(0, uBin), Ck_Lednicky.GetBinContent(uBin));
        graph1.SetPoint(uBin, Ck_Lednicky1.GetBinCenter(0, uBin), Ck_Lednicky1.GetBinContent(uBin));
    }
    graph.SetTitle("; #it{k}* (MeV); #it{C}(k*)");
    graph.GetXaxis()->SetLimits(0.0, 160);//our concern is R from 0.8 to 5 fm as we know that this range is most useful to study coalescence (available source size falls within this range :P)
    graph.GetYaxis()->SetRangeUser(0, 40);
    graph.SetLineColor(kBlack);
    graph.SetLineWidth(2);
    graph1.SetLineColor(kRed + 2);
    graph1.SetLineWidth(2);

    TLegend *legend = new TLegend(0.6, 0.65, 0.88, 0.85);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->SetLineColor(0);
    legend->AddEntry(&graph, "GeneralLednicky");
    legend->AddEntry(&graph1, "Lednicky_2channel" );

    graph.Draw("ACP");
    graph1.Draw("same");
    legend->Draw("same");
    binplot2->SaveAs(B2effplot);
    binplot2->Close();
    RootFile_DlmCk(OutputFileName, "Ck_pL_QA1", &Ck_Lednicky);
    // RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_QA1.root", "gCk_pL_Lednicky", &Ck_Lednicky);
}

void Ck_pL_Ledni_QA2() {
    TString OutputFileName = "../OutputFiles/Ck_pL_Ledni_QA2.root";
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 1.2);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 6, NumMomBins, kMin, kMax, Lednicky_2channel);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 2.88);
    Ck_Lednicky.SetPotPar(1, 2.92);
    Ck_Lednicky.SetPotPar(2, 1.66);
    Ck_Lednicky.SetPotPar(3, 3.78);
    Ck_Lednicky.SetPotPar(4, 0.25);
    Ck_Lednicky.SetPotPar(5, 0.75);
    Ck_Lednicky.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pL_QA2", &Ck_Lednicky);
    //RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_QA2.root", "gCk_pL_Lednicky", &Ck_Lednicky);
}
void Ck_np_Ledni_Usmani() {
    TString OutputFileName = "/home/sbhawani/Desktop/PdCFTheory/WaveFunctions/Output/QA_pd/Ck_np_Ledni_QA2_ExperimentEcattering.root";
    const double SourceSize = 1.2;
    const double kMin = 0;
    const double kMax = 130;
    const unsigned NumMomBins = 26;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, SourceSize);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 6, NumMomBins, kMin, kMax, Lednicky_2channel);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, -23.749);
    Ck_Lednicky.SetPotPar(1, -2.81);
    Ck_Lednicky.SetPotPar(2, 5.419);
    Ck_Lednicky.SetPotPar(3, 1.753);
    Ck_Lednicky.SetPotPar(4, 0.25);
    Ck_Lednicky.SetPotPar(5, 0.75);
    Ck_Lednicky.Update();
    RootFile_DlmCk(OutputFileName, "Ck_np_QA2_ExperimentScattering", &Ck_Lednicky);
    //RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_QA2.root", "gCk_pL_Lednicky", &Ck_Lednicky);
}
//Here we reproduce Haidenabuar resutls using LD scattering parameters
void Ck_LD_Ledni_QA2() {
    TString OutputFileName = "../OutputFiles/Ck_LD_LedniR_0p9.root";
    const unsigned NumMomBins = 160;
    const double kMin = 0;
    const double kMax = 160;

    TString B2effplot = "LambdaDeuteronForOtonR1p2.png";
    TCanvas *binplot2 = new TCanvas("binplot2", "Graph Draw Options", 1000, 1000);
    binplot2->SetLeftMargin(0.1505);
    binplot2->SetRightMargin(0.035);
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 1.2);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky1(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky2(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky3(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 16.8);
    Ck_Lednicky.SetPotPar(1, 3.2);
    Ck_Lednicky.SetPotPar(2, 0.0);
    Ck_Lednicky.SetPotPar(3, 0.0);
    Ck_Lednicky.SetPotPar(4, 0.333);
    Ck_Lednicky.SetPotPar(5, 0.0);
    Ck_Lednicky.Update();

///χEFT(NLO)
    Ck_Lednicky1.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky1.SetPotPar(0, 16.8);
    Ck_Lednicky1.SetPotPar(1, 3.2);
    Ck_Lednicky1.SetPotPar(2, -7.6);
    Ck_Lednicky1.SetPotPar(3, 3.6);
    Ck_Lednicky1.SetPotPar(4, 0.333);
    Ck_Lednicky1.SetPotPar(5, 0.666);
    Ck_Lednicky1.Update();

//NSC97f
    Ck_Lednicky2.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky2.SetPotPar(0, 16.8);
    Ck_Lednicky2.SetPotPar(1, 3.2);
    Ck_Lednicky2.SetPotPar(2, -10.8);
    Ck_Lednicky2.SetPotPar(3, 3.8);
    Ck_Lednicky2.SetPotPar(4, 0.333);
    Ck_Lednicky2.SetPotPar(5, 0.666);
    Ck_Lednicky2.Update();
/// Alexander
    Ck_Lednicky3.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky3.SetPotPar(0, 16.8);
    Ck_Lednicky3.SetPotPar(1, 3.2);
    Ck_Lednicky3.SetPotPar(2, -17.3);
    Ck_Lednicky3.SetPotPar(3, 3.6);
    Ck_Lednicky3.SetPotPar(4, 0.333);
    Ck_Lednicky3.SetPotPar(5, 0.666);
    Ck_Lednicky3.Update();
    RootFile_DlmCk(OutputFileName, "2S ", &Ck_Lednicky);
    RootFile_DlmCk(OutputFileName, "2S + 4SE", &Ck_Lednicky1);
    RootFile_DlmCk(OutputFileName, "2S + 4Sf", &Ck_Lednicky2);
    RootFile_DlmCk(OutputFileName, "2S + 4SA", &Ck_Lednicky3);

    //void RootFile_DlmCk(const TString & RootFileName, const TString & GraphName, DLM_Ck * CkToPlot) {

    const unsigned NumBins = Ck_Lednicky.GetNbins();
    TGraph graph;
    TGraph graph1;
    TGraph graph2;
    TGraph graph3;
    graph.SetName("2S ");
    graph1.SetName("2S + 4SE");
    graph2.SetName("2S + 4Sf");
    graph3.SetName("2S + 4SA");

    graph.Set(Ck_Lednicky.GetNbins());
    graph1.Set(Ck_Lednicky1.GetNbins());
    graph2.Set(Ck_Lednicky2.GetNbins());
    graph3.Set(Ck_Lednicky3.GetNbins());
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        graph.SetPoint(uBin, Ck_Lednicky.GetBinCenter(0, uBin), Ck_Lednicky.GetBinContent(uBin));
        graph1.SetPoint(uBin, Ck_Lednicky1.GetBinCenter(0, uBin), Ck_Lednicky1.GetBinContent(uBin));
        graph2.SetPoint(uBin, Ck_Lednicky2.GetBinCenter(0, uBin), Ck_Lednicky2.GetBinContent(uBin));
        graph3.SetPoint(uBin, Ck_Lednicky3.GetBinCenter(0, uBin), Ck_Lednicky3.GetBinContent(uBin));
    }
    graph.SetTitle("; #it{k}* (MeV); #it{C}(k*)");
    graph.GetXaxis()->SetLimits(0.0, 160);//our concern is R from 0.8 to 5 fm as we know that this range is most useful to study coalescence (available source size falls within this range :P)
    graph.GetYaxis()->SetRangeUser(0, 40);
    graph.SetLineColor(kBlack);
    graph.SetLineWidth(2);
    graph1.SetLineColor(kRed);
    graph1.SetLineWidth(2);
    graph2.SetLineColor(kGreen);
    graph2.SetLineWidth(2);
    graph3.SetLineColor(kBlue);
    graph3.SetLineWidth(2);

    TLegend *legend = new TLegend(0.6, 0.65, 0.88, 0.85);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->SetLineColor(0);
    legend->AddEntry(&graph, "2S ");
    legend->AddEntry(&graph1, "2S + 4SE" );
    legend->AddEntry(&graph2, "2S + 4Sf" );
    legend->AddEntry(&graph3, "2S + 4SA");

    graph.Draw("ACP");
    graph1.Draw("same");
    graph2.Draw("same");
    graph3.Draw("same");
    legend->Draw("same");
    binplot2->SaveAs(B2effplot);
    binplot2->Close();
}
void Ck_pd_Lenicky_fromSebastian() {

    TString OutputFileName = "../OutputFiles/Ck_pD_LedniOnlyforPd.root";
    const unsigned NumMomBins = 50;
    const double kMin = 0;
    const double kMax = 450;

    TString B2effplot = "Pd_LenickyUsingSebastianSlidesFromLendniOnly.png";

    TCanvas *binplot2 = new TCanvas("binplot2", "Graph Draw Options", 1000, 1000);
    binplot2->SetLeftMargin(0.1505);
    binplot2->SetRightMargin(0.035);
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 0.9);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 6, NumMomBins, kMin, kMax,  LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky1(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky2(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky3(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    DLM_Ck Ck_Lednicky4(1, 6, NumMomBins, kMin, kMax, LednickyOnlyForPdLD);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.

//van Oers, Brockman (1967)
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 11.4);
    Ck_Lednicky.SetPotPar(1, 0);
    Ck_Lednicky.SetPotPar(2, 1.2);
    Ck_Lednicky.SetPotPar(3, 0.0);
    Ck_Lednicky.SetPotPar(4, 0.33);
    Ck_Lednicky.SetPotPar(5, 0.666);
    Ck_Lednicky.Update();

///Arvieux (1973)
    Ck_Lednicky1.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky1.SetPotPar(0, 11.8);
    Ck_Lednicky1.SetPotPar(1, 0);
    Ck_Lednicky1.SetPotPar(2, 2.73);
    Ck_Lednicky1.SetPotPar(3, 0);
    Ck_Lednicky1.SetPotPar(4, 0.333);
    Ck_Lednicky1.SetPotPar(5, 0.666);
    Ck_Lednicky1.Update();

//Huttelet al.(1983)
    Ck_Lednicky2.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky2.SetPotPar(0, 11.1);
    Ck_Lednicky2.SetPotPar(1, 0);
    Ck_Lednicky2.SetPotPar(2, 4.0);
    Ck_Lednicky2.SetPotPar(3, 0);
    Ck_Lednicky2.SetPotPar(4, 0.333);
    Ck_Lednicky2.SetPotPar(5, 0.666);
    Ck_Lednicky2.Update();
/// Kievskyet al.(1997)
    Ck_Lednicky3.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky3.SetPotPar(0, 13.8);
    Ck_Lednicky3.SetPotPar(1, 0);
    Ck_Lednicky3.SetPotPar(2, 0.0);
    Ck_Lednicky3.SetPotPar(3, 0);
    Ck_Lednicky3.SetPotPar(4, 0.333);
    Ck_Lednicky3.SetPotPar(5, 0.666);
    Ck_Lednicky3.Update();

/// Black et al.(1999)
    Ck_Lednicky4.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky4.SetPotPar(0, 14.7);
    Ck_Lednicky4.SetPotPar(1, 0);
    Ck_Lednicky4.SetPotPar(2, -0.13);
    Ck_Lednicky4.SetPotPar(3, 0);
    Ck_Lednicky4.SetPotPar(4, 0.333);
    Ck_Lednicky4.SetPotPar(5, 0.666);
    Ck_Lednicky4.Update();

    RootFile_DlmCk(OutputFileName, "Brockman ", &Ck_Lednicky);
    RootFile_DlmCk(OutputFileName, "Arvieux", &Ck_Lednicky1);
    RootFile_DlmCk(OutputFileName, "Huttel", &Ck_Lednicky2);
    RootFile_DlmCk(OutputFileName, "Kievsky", &Ck_Lednicky3);
    RootFile_DlmCk(OutputFileName, "Black", &Ck_Lednicky4);

    //void RootFile_DlmCk(const TString & RootFileName, const TString & GraphName, DLM_Ck * CkToPlot) {

    const unsigned NumBins = Ck_Lednicky.GetNbins();
    TGraph graph;
    TGraph graph1;
    TGraph graph2;
    TGraph graph3;
    TGraph graph4;

    graph.SetName("Brockman ");
    graph1.SetName("Arvieux");
    graph2.SetName("Huttel");
    graph3.SetName("Kievsky");
    graph4.SetName("Black");

    graph.Set(Ck_Lednicky.GetNbins());
    graph1.Set(Ck_Lednicky1.GetNbins());
    graph2.Set(Ck_Lednicky2.GetNbins());
    graph3.Set(Ck_Lednicky3.GetNbins());
    graph4.Set(Ck_Lednicky4.GetNbins());
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        printf("k*= %0.1f Ck1 = %0.1f Ck2 = %0.1f Ck2 = %0.1fCk4 = %0.1f Ck5 = %0.1f \n", Ck_Lednicky.GetBinCenter(0, uBin), Ck_Lednicky.GetBinContent(uBin), Ck_Lednicky1.GetBinContent(uBin), Ck_Lednicky2.GetBinContent(uBin), Ck_Lednicky3.GetBinContent(uBin), Ck_Lednicky4.GetBinContent(uBin));
        graph.SetPoint(uBin, Ck_Lednicky.GetBinCenter(0, uBin), Ck_Lednicky.GetBinContent(uBin));
        graph1.SetPoint(uBin, Ck_Lednicky1.GetBinCenter(0, uBin), Ck_Lednicky1.GetBinContent(uBin));
        graph2.SetPoint(uBin, Ck_Lednicky2.GetBinCenter(0, uBin), Ck_Lednicky2.GetBinContent(uBin));
        graph3.SetPoint(uBin, Ck_Lednicky3.GetBinCenter(0, uBin), Ck_Lednicky3.GetBinContent(uBin));
        graph4.SetPoint(uBin, Ck_Lednicky4.GetBinCenter(0, uBin), Ck_Lednicky4.GetBinContent(uBin));
    }
    graph.SetTitle("; #it{k}* (MeV); #it{C}(k*)");
    graph.GetXaxis()->SetLimits(0.0, 400);//our concern is R from 0.8 to 5 fm as we know that this range is most useful to study coalescence (available source size falls within this range :P)
    graph.GetYaxis()->SetRangeUser(0, 60);
    graph.SetLineColor(kBlack);
    graph.SetLineWidth(2);
    graph1.SetLineColor(kRed + 2);
    graph1.SetLineWidth(2);
    graph2.SetLineColor(kGreen + 2);
    graph2.SetLineWidth(2);
    graph3.SetLineColor(kBlue + 2);
    graph3.SetLineWidth(2);
    graph4.SetLineColor(kMagenta + 2);
    graph4.SetLineWidth(2);

    TLegend *legend = new TLegend(0.6, 0.65, 0.88, 0.85);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->SetLineColor(0);
    legend->AddEntry(&graph, "Brockman ");
    legend->AddEntry(&graph1, "Arvieux" );
    legend->AddEntry(&graph2, "Huttel" );
    legend->AddEntry(&graph3, "Kievsky");
    legend->AddEntry(&graph4, "Black");
    graph.Draw("ACP");
    graph1.Draw("same");
    graph2.Draw("same");
    graph3.Draw("same");
    graph4.Draw("same");
    legend->Draw("same");
    binplot2->SaveAs(B2effplot);
    binplot2->Close();
}
void Ck_PD_Ledni_QA2() {
    TString OutputFileName = "../OutputFiles/Ck_PD_LedniQDWaves.root";
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 150;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 0.94);



    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 6, NumMomBins, kMin, kMax, Lednicky_2channel);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, -0.06);
    Ck_Lednicky.SetPotPar(1, 0.0);
    Ck_Lednicky.SetPotPar(2, 10.9);
    Ck_Lednicky.SetPotPar(3, 0);
    Ck_Lednicky.SetPotPar(4, 0.333);
    Ck_Lednicky.SetPotPar(5, 0.666);
    Ck_Lednicky.Update();
    RootFile_DlmCk(OutputFileName, "2S + 4S PD", &Ck_Lednicky);
    //RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_QA2.root", "gCk_pL_Lednicky", &Ck_Lednicky);
}
void PlotCF(const TString& GraphName, DLM_Ck* CkToPlot) {

    TString B2effplot = "LambdaDeuteronGenuine.pdf";
    TCanvas *binplot2 = new TCanvas("binplot2", "Graph Draw Options", 1000, 1000);
    binplot2->SetLeftMargin(0.1505);
    binplot2->SetRightMargin(0.035);

    const unsigned NumBins = CkToPlot->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        graph.SetPoint(uBin, CkToPlot->GetBinCenter(0, uBin), CkToPlot->GetBinContent(uBin));
    }
    graph.SetTitle("; #it{k}* (MeV); #it{C}(k*)");
    graph.GetXaxis()->SetLimits(0.0, 405);//our concern is R from 0.8 to 5 fm as we know that this range is most useful to study coalescence (available source size falls within this range :P)
    graph.GetYaxis()->SetRangeUser(0, 13);
    graph.SetLineColor(kBlue + 2);
    graph.SetLineWidth(2);
    graph.Draw("ALP");
    binplot2->SaveAs(B2effplot);
    binplot2->Close();

}
//save a DLM_Ck to a file in the form of a TGraph
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot) {
    TFile* RootFile = new TFile(RootFileName, "update");
    if (!RootFile) RootFile = new TFile(RootFileName, "recreate");
    const unsigned NumBins = CkToPlot->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        graph.SetPoint(uBin, CkToPlot->GetBinCenter(0, uBin), CkToPlot->GetBinContent(uBin));
    }
    graph.Write("", TObject::kOverwrite);
    delete RootFile;
}
//save a DLM_CkDecomposition to a file in the form of a TGraph
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_CkDecomposition* CkToPlot, bool PlotIndividualContributions) {
    TFile* RootFile = new TFile(RootFileName, "update");
    if (!RootFile) RootFile = new TFile(RootFileName, "recreate");
    const unsigned NumBins = CkToPlot->GetCk()->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
        graph.SetPoint(uBin, CkToPlot->GetCk()->GetBinCenter(0, uBin), CkToPlot->EvalCk(CkToPlot->GetCk()->GetBinCenter(0, uBin)));
    }
    graph.Write("", TObject::kOverwrite);

    //this bit is to get the individual contributions (associated to lambda pars) to the total correlation function
    if (PlotIndividualContributions) {
        TGraph graphSignal;
        graphSignal.SetName(GraphName + "_signal");
        graphSignal.Set(NumBins);
        for (unsigned uBin = 0; uBin < NumBins; uBin++) {
            double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0, uBin);
            //the signal is defined as C(k)-1
            graphSignal.SetPoint(uBin, MOMENTUM, CkToPlot->EvalSignalSmeared(MOMENTUM));
        }
        graphSignal.Write("", TObject::kOverwrite);

        TGraph graphMain;
        graphMain.SetName(GraphName + "_main");
        graphMain.Set(NumBins);
        for (unsigned uBin = 0; uBin < NumBins; uBin++) {
            double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0, uBin);
            //this is the contribution of the main channel only
            graphMain.SetPoint(uBin, MOMENTUM, CkToPlot->EvalSignalSmearedMain(MOMENTUM));
        }
        graphMain.Write("", TObject::kOverwrite);

        for (unsigned uChild = 0; uChild < CkToPlot->GetNumChildren(); uChild++) {
            TGraph graphChild;
            graphChild.SetName(GraphName + "_child" + uChild);
            graphChild.Set(NumBins);
            for (unsigned uBin = 0; uBin < NumBins; uBin++) {
                double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0, uBin);
                //this is the contribution of each secondary channel, e.g. feed-down or misid
                graphChild.SetPoint(uBin, MOMENTUM, CkToPlot->EvalSignalSmearedChild(uChild, MOMENTUM));
            }
            graphChild.Write("", TObject::kOverwrite);
        }
    }
    delete RootFile;
}

//initialize the Gauss source for the CATS object
void CATS_GaussSource(CATS& Kitty, const double& SourceSize) {
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


DLM_CleverMcLevyResoTM* CATS_ResoSource_pd(CATS& Kitty, const double& SourceSize, const double& Alpha, const bool& FixAlpha,
        const double& fracreso, const double& taureso, const double& massreso, const double& cutoff) {
//CATS_ResoSource_pd(Kitty_pd, SourceSize, 2.0, true, 0.6422, 1.65, 1362, 200)
//Kitty_pd, SourceSize, 2.0, true, 0.6422, 1.65, 1362, 200
    DLM_CleverMcLevyResoTM* CleverMcLevyResoTM_pd = new DLM_CleverMcLevyResoTM();;

    if (FixAlpha) CleverMcLevyResoTM_pd->InitStability(1, Alpha - 1e-6, Alpha + 1e-6); //only a single grid point for the Alpha
    else CleverMcLevyResoTM_pd->InitStability(21, 1, 2);
    //this sets the possible source sizes between 0.4 and 3 fm. Typically should provide good cover of all ranges, but take care you do not go outside
    CleverMcLevyResoTM_pd->InitScale(52, 0.4, 3.0);
    CleverMcLevyResoTM_pd->InitRad(400, 0, 64);
    CleverMcLevyResoTM_pd->InitType(2);
    CleverMcLevyResoTM_pd->SetUpReso(0, fracreso);

    Float_t k_D, fP1, fP2, fM1, fM2, Tau1, Tau2, AngleRcP1, AngleRcP2, AngleP1P2;
    DLM_Random RanGen(11);
    double RanVal1, RanVal2;

    // I mimic the pd source using pOmega (Well we don't have deuteron in EPOS so omega is best suiatble in place of deuteron)
    TFile* F_EposDisto_p_Omega_Reso = new TFile("/home/sbhawani/Desktop/CATS_CF_Pd/RootFiles/eposdisto-preso-omega.root");
    TNtuple* T_EposDisto_d_pReso = (TNtuple*)F_EposDisto_p_Omega_Reso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_d_pReso = T_EposDisto_d_pReso->GetEntries();
    T_EposDisto_d_pReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_d_pReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_d_pReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_d_pReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_d_pReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_d_pReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_d_pReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_d_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_d_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_d_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < N_EposDisto_d_pReso; uEntry++) {
        T_EposDisto_d_pReso->GetEntry(uEntry);
        Tau1 = 1.65;
        //the values were provided from the statistical hadronization model, we rewrite the EPOS output
        Tau2 = 0;
        fM1 = 1362;
        //reject daughters with too high relative momenta
        if (k_D > cutoff) continue;
        //sample randomly the fly length of the resonance
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        //plug in the result for the primordial-resonance case
        CleverMcLevyResoTM_pd->AddBGT_PR(RanVal1, cos(AngleRcP2));
        //plug in the result for the primordial-resonance case. For identical particle it is symmetric, apart the sign of the cosine.vid
    }
    delete F_EposDisto_p_Omega_Reso;
    CleverMcLevyResoTM_pd->InitNumMcIter(1000000);

    Kitty.SetAnaSource(CatsSourceForwarder, CleverMcLevyResoTM_pd, 2);
    //and this, as aways, is the way to communicate with CATS and set the source size
    Kitty.SetAnaSource(0, SourceSize);
    Kitty.SetAnaSource(1, Alpha);

    //this is important, since if its `false` it is assumed that the source will be sampled from a transport model, and the Gaussian function will not be used
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAutoNormSource(false);

    return CleverMcLevyResoTM_pd;
}

//Basic initialization of a CATS object for pp, WITHOUT any setup of the source or interaction
void CATS_pp_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax) {
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
//Basic initialization of a CATS object for pp, WITHOUT any setup of the source or interaction
void CATS_pd_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax) {
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 1000010020);
    Kitty.SetRedMass( Mass_p * Mass_d / (Mass_p + Mass_d));
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);
}

//pi^+-proton feed down for pd
void CATS_PiPronton_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax) {
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 211);
    Kitty.SetRedMass( Mass_p * Mass_PiPlus / (Mass_p + Mass_PiPlus));
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);
}

//initialize the interaction for pp, using the AV18 potential.
//by default we have only s-waves included, but one can include the p and d waves if desired
void CATS_pp_AV18(CATS& Kitty, const bool& pwaves, const bool& dwaves) {
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
        Kitty.SetChannelWeight(0, 1. / 4.);
        Kitty.SetChannelWeight(1, 3. / 4.);
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
void CATS_pL_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax) {
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass( (Mass_p * Mass_L) / (Mass_p + Mass_L) );
}
void CATS_pL_Usmani(CATS& Kitty) {
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1. / 4.);
    Kitty.SetChannelWeight(1, 3. / 4.);

    CATSparameters cPars_pL_Usmani1S0(CATSparameters::tPotential, 8, true);
    cPars_pL_Usmani1S0.SetParameter(0, pL_UsmaniOli);
    cPars_pL_Usmani1S0.SetParameter(5, 0);
    cPars_pL_Usmani1S0.SetParameter(7, 0);

    CATSparameters cPars_pL_Usmani3S1(CATSparameters::tPotential, 8, true);
    cPars_pL_Usmani3S1.SetParameter(0, pL_UsmaniOli);
    cPars_pL_Usmani3S1.SetParameter(5, 1);
    cPars_pL_Usmani3S1.SetParameter(7, 1);

    Kitty.SetShortRangePotential(0, 0, fDlmPot, cPars_pL_Usmani1S0);
    Kitty.SetShortRangePotential(1, 0, fDlmPot, cPars_pL_Usmani3S1);
}

void Ck_dL_Ledni() {

    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 300;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 0.9);

    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 2, NumMomBins, kMin, kMax, Lednicky_Singlet);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 19.1);
    Ck_Lednicky.SetPotPar(1, 3.2);
    Ck_Lednicky.Update();
    PlotCF("gCk_dL_Lednicky", &Ck_Lednicky);
}
void CATS_pL_NLO(CATS& Kitty) {
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1. / 4.);
    Kitty.SetChannelWeight(1, 3. / 4.);

    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    DLM_Histo<complex<double>>*** ExternalWF = nullptr;
    ExternalWF = Init_pL_Haidenbauer(
                     TString::Format("%s/Desktop/RadiiFit/WaveFunctions/Haidenbauer/pLambdaNLO/",
                                     HomeDir.Data()),
                     Kitty, 10, 600);
    for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
        Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                      ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(Kitty.GetNumChannels(), ExternalWF);
    std::cout << "aree Bhai, pL Potential is NLO ChiEFT\n";
}

void CATS_pL_LO(CATS& Kitty) {
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1. / 4.);
    Kitty.SetChannelWeight(1, 3. / 4.);
    TString HomeDir = gSystem->GetHomeDirectory().c_str();
    DLM_Histo<complex<double>>*** ExternalWF = nullptr;
    ExternalWF = Init_pL_Haidenbauer(
                     TString::Format("%s/Desktop/RadiiFit/WaveFunctions/Haidenbauer/pLambdaLO_600/",
                                     HomeDir.Data()),
                     Kitty, 0, 600);
    for (unsigned uCh = 0; uCh < Kitty.GetNumChannels(); uCh++) {
        Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0],
                                      ExternalWF[1][uCh][0]);
    }
    CleanUpWfHisto(Kitty.GetNumChannels(), ExternalWF);
    std::cout << "pL Potential is LO ChiEFT\n";
}

void CATS_pXim_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax) {
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    Kitty.SetRedMass( (Mass_p * Mass_Xim) / (Mass_p + Mass_Xim) );
}
void CATS_pXim_Hal(CATS& Kitty) {
    //0: S=0 I=0
    //1: S=1 I=0
    //2: S=0 I=1
    //3: S=1 I=1
    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetNumPW(2, 1);
    Kitty.SetNumPW(3, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetSpin(2, 0);
    Kitty.SetSpin(3, 1);
    Kitty.SetChannelWeight(0, 1. / 8.);
    Kitty.SetChannelWeight(1, 3. / 8.);
    Kitty.SetChannelWeight(2, 1. / 8.);
    Kitty.SetChannelWeight(3, 3. / 8.);

    CATSparameters cPars_pXim_HalI0S0(CATSparameters::tPotential, 8, true);
    cPars_pXim_HalI0S0.SetParameter(0, pXim_HALQCD1);
    cPars_pXim_HalI0S0.SetParameter(1, 12);
    cPars_pXim_HalI0S0.SetParameter(2, 0);
    cPars_pXim_HalI0S0.SetParameter(3, -1);
    cPars_pXim_HalI0S0.SetParameter(4, 1);
    cPars_pXim_HalI0S0.SetParameter(5, 0);
    cPars_pXim_HalI0S0.SetParameter(6, 0);
    cPars_pXim_HalI0S0.SetParameter(7, 0);

    CATSparameters cPars_pXim_HalI0S1(CATSparameters::tPotential, 8, true);
    cPars_pXim_HalI0S1.SetParameter(0, pXim_HALQCD1);
    cPars_pXim_HalI0S1.SetParameter(1, 12);
    cPars_pXim_HalI0S1.SetParameter(2, 0);
    cPars_pXim_HalI0S1.SetParameter(3, -1);
    cPars_pXim_HalI0S1.SetParameter(4, 1);
    cPars_pXim_HalI0S1.SetParameter(5, 1);
    cPars_pXim_HalI0S1.SetParameter(6, 0);
    cPars_pXim_HalI0S1.SetParameter(7, 1);

    CATSparameters cPars_pXim_HalI1S0(CATSparameters::tPotential, 8, true);
    cPars_pXim_HalI1S0.SetParameter(0, pXim_HALQCD1);
    cPars_pXim_HalI1S0.SetParameter(1, 12);
    cPars_pXim_HalI1S0.SetParameter(2, 1);
    cPars_pXim_HalI1S0.SetParameter(3, 1);
    cPars_pXim_HalI1S0.SetParameter(4, 1);
    cPars_pXim_HalI1S0.SetParameter(5, 0);
    cPars_pXim_HalI1S0.SetParameter(6, 0);
    cPars_pXim_HalI1S0.SetParameter(7, 0);

    CATSparameters cPars_pXim_HalI1S1(CATSparameters::tPotential, 8, true);
    cPars_pXim_HalI1S1.SetParameter(0, pXim_HALQCD1);
    cPars_pXim_HalI1S1.SetParameter(1, 12);
    cPars_pXim_HalI1S1.SetParameter(2, 1);
    cPars_pXim_HalI1S1.SetParameter(3, 1);
    cPars_pXim_HalI1S1.SetParameter(4, 1);
    cPars_pXim_HalI1S1.SetParameter(5, 1);
    cPars_pXim_HalI1S1.SetParameter(6, 0);
    cPars_pXim_HalI1S1.SetParameter(7, 1);

    Kitty.SetShortRangePotential(0, 0, fDlmPot, cPars_pXim_HalI0S0);
    Kitty.SetShortRangePotential(1, 0, fDlmPot, cPars_pXim_HalI0S1);
    Kitty.SetShortRangePotential(2, 0, fDlmPot, cPars_pXim_HalI1S0);
    Kitty.SetShortRangePotential(3, 0, fDlmPot, cPars_pXim_HalI1S1);
}


TH2F* GetSmearMatrix(TString filename, TString histoname) {
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if (!histo) {printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n", histoname.Data(), filename.Data()); return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}
TH1F* GetExpCorrelation(TString filename, TString histoname) {
    TFile* FileROOT = new TFile(filename, "read");
    if (!FileROOT) {printf("\033[1;31mERROR:\033[0m The file '%s' does not exist\n", filename.Data()); return NULL;}
    TH1F* histo = (TH1F*)FileROOT->Get(histoname);
    if (!histo) {printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n", histoname.Data(), filename.Data()); return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT; FileROOT = NULL;
    histoCopy->SetName(Name);
    return histoCopy;
}

DLM_CkDecomposition* Ck_Fitter_pd;

void GetCATS_PL_Usmani_CF(const TString& SourceType) {
    TString OutputFileName = "../OutputFiles/Ck_pL_needed_forPd_LowR" + SourceType + ".root";
    const double SourceSize_pL = 0.92;
    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;

    const double FitMin = 0;
    const double FitMax = 350;

    CATS Kitty_pL_Usmani1;
    CATS_GaussSource(Kitty_pL_Usmani1, SourceSize_pL);
    CATS_pL_Basic(Kitty_pL_Usmani1, 100, 0, 400);
    CATS_pL_Usmani(Kitty_pL_Usmani1);
    printf(" Executing KillTheCat for Kitty_pL_Usmani\n");
    Kitty_pL_Usmani1.SetNotifications(CATS::nWarning);
    Kitty_pL_Usmani1.KillTheCat();
    DLM_Ck Ck_pL_Usmani1(Kitty_pL_Usmani1.GetNumSourcePars(), 0, Kitty_pL_Usmani1);
    Ck_pL_Usmani1.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pL_Usmani", &Ck_pL_Usmani1);
}
void GetCATS_Ld_Lednicky_CF(const TString& SourceType) {
    TString OutputFileName = "../OutputFiles/Ck_Ld_LednickyRightHandHaidenabur_needed_forPd_LowR" + SourceType + ".root";
    const double SourceSize_pL = 0.9;
    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;

    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    SOURCE_PARS.SetParameter(0, 0.9);
    DLM_Ck Ck_Lednicky(1, 2, NumMomBins, kMin, kMax, Lednicky_Singlet);
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 10.01);
    Ck_Lednicky.SetPotPar(1, 3.045);
    Ck_Lednicky.Update();

    printf(" Executing KillTheCat for Kitty_Ld_Lednicky\n");
    RootFile_DlmCk(OutputFileName, "Ck_Ld_lednicky", &Ck_Lednicky);
}
void GetCATS_PL_LO_CF(const TString& SourceType) {
    TString OutputFileName = "../OutputFiles/Ck_pL_LO_needed_forPd_LowR" + SourceType + ".root";
    const double SourceSize_pL = 0.9;
    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;

    const double FitMin = 0;
    const double FitMax = 350;

    CATS Kitty_pL_LO;
    CATS_GaussSource(Kitty_pL_LO, SourceSize_pL);
    CATS_pL_Basic(Kitty_pL_LO, 100, 0, 400);
    CATS_pL_LO(Kitty_pL_LO);
    printf(" Executing KillTheCat for Kitty_pL_LO\n");
    Kitty_pL_LO.SetNotifications(CATS::nWarning);
    Kitty_pL_LO.KillTheCat();
    DLM_Ck Ck_pL_LO(Kitty_pL_LO.GetNumSourcePars(), 0, Kitty_pL_LO);
    Ck_pL_LO.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pL_LO", &Ck_pL_LO);
}
void GetCATS_PL_NLO_CF(const TString& SourceType) {
    TString OutputFileName = "../OutputFiles/Ck_pL_NLO_needed_forPd_LowR" + SourceType + ".root";
    const double SourceSize_pL = 0.9;
    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;

    const double FitMin = 0;
    const double FitMax = 350;

    CATS Kitty_pL_NLO;
    CATS_GaussSource(Kitty_pL_NLO, SourceSize_pL);
    CATS_pL_Basic(Kitty_pL_NLO, 100, 0, 400);
    CATS_pL_NLO(Kitty_pL_NLO);
    printf(" Executing KillTheCat for Kitty_pL_NLO\n");
    Kitty_pL_NLO.SetNotifications(CATS::nWarning);
    Kitty_pL_NLO.KillTheCat();
    DLM_Ck Ck_pL_NLO(Kitty_pL_NLO.GetNumSourcePars(), 0, Kitty_pL_NLO);
    Ck_pL_NLO.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pL_NLO", &Ck_pL_NLO);
}

double BaseLine(double *x, double *par) {
    double& MOM = *x;
    double Baseline = (par[0] + par[1] * MOM);// + par[2] * MOM * MOM);
    //double Baseline = (par[0]); //+ par[1] * MOM + par[2] * MOM * MOM);
    return Baseline;
}
/*double BaseLine(double *x, double *par) {
    double& MOM = *x;
    double Baseline = (par[0] + par[1] * MOM);
    return Baseline;
}*/
double CkSignal(double *x, double *par) {
    double& MOM = *x;
    double Correlation = par[0] * Ck_Fitter_pd->EvalCk(MOM);
    return Correlation;
}
double Fitter_pd(double* x, double* par) {
    double& MOM = *x;
    return BaseLine(x, &par[1]) + CkSignal(x, &par[0]);
}

void Ck_pd_Decomposition(const TString& SourceType) {
    //timer
    TCanvas *CFplot = new TCanvas ("CFplot", "CFPlots", 1000, 1000);
    CFplot->SetLeftMargin(0.1505);
    CFplot->SetRightMargin(0.035);

    std::cout << "HI its here\n";
    DLM_Timer TIMER;
    printf("\033[1;37mExecuting Ck_pp_Decomposition for SourceType=='%s'\033[0m\n\n", SourceType.Data());
    if (SourceType != "Gauss" && SourceType != "CoreReso") {printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SourceType.Data()); return;}
    TString OutputFileName = "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/LowpT1p4/FITCATS/LatestCkFITPd_R_1p0Correct.root";
    //set the source size slightly smaller for the Gauss+Reso scenario.
    const double SourceSize = SourceType == "Gauss" ? 0.92 : 0.93;
    const double SourceSize_pL = 1.0;
    const double SourceSize_pXim = 1.2;

    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 300;
    const double FitMin = 0;
    const double FitMax = 300;
    const TString BL = "Pol1";

    CATS Kitty_pd;
    CATS Kitty_PiprotonFeed;
    CATS_GaussSource(Kitty_PiprotonFeed, 1.0);//fm
    if (SourceType == "Gauss") CATS_GaussSource(Kitty_pd, 0.944);
    //initialize a Gaussian core + reso source
    DLM_CleverMcLevyResoTM* CRS_pd = NULL;
    if (SourceType == "CoreReso") CRS_pd = CATS_ResoSource_pd(Kitty_pd, SourceSize, 2.0, true, 0.6422, 1.65, 1362, 200);


    CATS_pd_Basic(Kitty_pd, NumMomBins, kMin, kMax);
    printf(" Executing KillTheCat for Kitty_pd\n");
    Kitty_pd.KillTheCat();
    DLM_Ck Ck_pd(Kitty_pd.GetNumSourcePars(), 0, Kitty_pd);
    Ck_pd.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pd_Coulumb_only", &Ck_pd);

    CATS_PiPronton_Basic(Kitty_PiprotonFeed, NumMomBins, kMin, kMax);
    printf(" Executing KillTheCat for Kitty_PiprotonFeed\n");
    Kitty_PiprotonFeed.KillTheCat();

    DLM_Ck Ck_PiProtonFeed(Kitty_PiprotonFeed.GetNumSourcePars(), 0, Kitty_PiprotonFeed);
    Ck_PiProtonFeed.Update();
    RootFile_DlmCk("../OutputFiles/PiprotonFeedToPdR_0p1.root", "Ck_Piprotonfeed_Coulumb_only", &Ck_PiProtonFeed);

    //we will include feed-downs from Lambda and Xi, thus we need to compute pL and pXi correlations as well
    /* CATS Kitty_pL_Usmani;
     CATS_GaussSource(Kitty_pL_Usmani, SourceSize_pL);
     CATS_pL_Basic(Kitty_pL_Usmani, 100, 0, 400);
     CATS_pL_Usmani(Kitty_pL_Usmani);
     printf(" Executing KillTheCat for Kitty_pL_Usmani\n");
     //CATS_pL_NLO(Kitty_pL_Usmani);
     //printf(" Executing KillTheCat for Kitty_pL_NLO\n");
     //CATS_pL_LO(Kitty_pL_Usmani);
     //printf(" Executing KillTheCat for Kitty_pL_LO\n");

     Kitty_pL_Usmani.SetNotifications(CATS::nWarning);
     Kitty_pL_Usmani.KillTheCat();
     DLM_Ck Ck_pL_Usmani(Kitty_pL_Usmani.GetNumSourcePars(), 0, Kitty_pL_Usmani);
     Ck_pL_Usmani.Update();
     RootFile_DlmCk(OutputFileName, "Ck_pL_Usmani", &Ck_pL_Usmani);
    */
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    SOURCE_PARS.SetParameter(0, 0.9);
    DLM_Ck Ck_Lednicky(1, 2, NumMomBins, kMin, kMax, Lednicky_Singlet);
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 10.1);
    Ck_Lednicky.SetPotPar(1, 3.045);
    Ck_Lednicky.Update();

    PlotCF("LD", &Ck_Lednicky);
    //some inputs snearing matrices
    TString FileName_MomSmear = "../Files/pp13TeV_MB.root";
    TString HistoName_MomSmear = "hSigmaMeV_Proton_Proton";
    TH2F* hResolution_pd = GetSmearMatrix(FileName_MomSmear, HistoName_MomSmear);

    // decay metrices
    TString FileName_Feed_pL_pp = "../Files/Decay_matrix_dp_dL.root";
    TString HistoName_Feed_pL_pp = "hRes_dp_dL";
    TH2F* hFeedDown_pL_pp = GetSmearMatrix(FileName_Feed_pL_pp, HistoName_Feed_pL_pp);


    // wrong k* calculation matrix
    TString FileName_Feed_pd_ppi = "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/CATSOutput/DecayMatrices/DeltaDecayCorrect1.root";
    TString HistoName_Feed_pd_ppi = "hRes_dp_ppi";
    TH2F* hFeedDown_pd_ppi = GetSmearMatrix(FileName_Feed_pd_ppi, HistoName_Feed_pd_ppi);

    TString FileName_Feed_pXim_pL = "../Files/DecayMatrices_Oli.root";
    TString HistoName_Feed_pXim_pL = "hRes_pL_pXim";
    TH2F* hFeedDown_pXim_pL = GetSmearMatrix(FileName_Feed_pXim_pL, HistoName_Feed_pXim_pL);

    //the object is initialized by: a name of your choice, number of contributions to C(k) APART from the primary, DLM_Ck, resolution matrix
    DLM_CkDecomposition CorrCkDec_pd("pd", 3, Ck_pd, hResolution_pd);
    DLM_CkDecomposition CkDec_pL("dLambda", 2, Ck_Lednicky, NULL);
    DLM_CkDecomposition CkDec_PiProton("PiProton", 0, Ck_PiProtonFeed, hResolution_pd);
    //DLM_CkDecomposition CkDec_pXim("pXim", 2, Ck_pXim_Hal, NULL);
    //example lambda pars from Phys. Lett. B 797 (2019) 134822
    /*  double lambda_pL_d = 0.086;
      double lambda_pd_flat = 0.03;
      double lambda_pd_misid = 0.12;*/
    double lambda_pL_d = 0.094;
    double lambda_pd_flat = 0.06;
    double lambda_pd_misid = 0.01;
    //CorrCkDec_pd.AddContribution(0, lambda_pd_flat, DLM_CkDecomposition::cFeedDown);
   // CorrCkDec_pd.AddContribution(1, lambda_pd_misid, DLM_CkDecomposition::cFake);
    CorrCkDec_pd.AddContribution(0, lambda_pL_d, DLM_CkDecomposition::cFeedDown, &CkDec_pL, hFeedDown_pL_pp);
    CorrCkDec_pd.AddContribution(1, lambda_pd_flat, DLM_CkDecomposition::cFeedDown);
    CorrCkDec_pd.AddContribution(2, lambda_pd_misid, DLM_CkDecomposition::cFake);
    
    double lambda_pL_flat =  0.541;
    double lambda_pL_misid = 0.042;
    //CkDec_pL.AddContribution(0, lambda_pXim_pL, DLM_CkDecomposition::cFeedDown, &CkDec_pXim, hFeedDown_pXim_pL);
    CkDec_pL.AddContribution(0, lambda_pL_flat, DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1, lambda_pL_misid, DLM_CkDecomposition::cFake);


    CorrCkDec_pd.Update(true, true);
    CkDec_pL.Update(true, true);
    //CkDec_pXim.Update(true, true);

    RootFile_DlmCk(OutputFileName, "CorrCkDec_pd", &CorrCkDec_pd, true);

    std::cout << " Warning!!! Hey Buddy look at the CF you took from experiment,, this should be pdCF in principle otherwise you are doing wrong" << std::endl;

    TString DataFileName = "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/Systematics/systematicslow/CF_pd_Var0.root";
    TString DataHistoName = "hCk_ReweightedpdVar0MeV_1";
    TH1F* hExpCk = GetExpCorrelation(DataFileName, DataHistoName);

    Kitty_pd.SetNotifications(CATS::nWarning);

    //set the pointer to the object that we want to use in ExampleFitter_pp
    Ck_Fitter_pd = &CorrCkDec_pd;
    printf(" Fitting begins now!\n");

    //pointer to the DLM_CkDecomposition of the p-p to perform the fit, which has 3 fit parameters
    TF1* fit_pd = new TF1("fit_pd", Fitter_pd, FitMin, FitMax, 2);
    //source size
    fit_pd->FixParameter(0, 1.0);
    fit_pd->SetParameter(1, 0.001);
    fit_pd->SetParLimits(1, 0.005, 2.0);

    fit_pd->SetParameter(2, 0.001);
    fit_pd->SetParLimits(2, -1e-3, 2.0);

    //fit_pd->SetParameter(3, 0.001);
    //fit_pd->SetParLimits(3, -1e-3, 2);
    //in this example we fit with norm only
    hExpCk->Fit(fit_pd, "S,R, M");

    TF1* Bline = new TF1("Bline", BaseLine, FitMin, FitMax, 2);
    Bline->SetLineColor(kGreen);
    Bline->SetParameters(fit_pd->GetParameter(1), fit_pd->GetParameter(2));


    TF1* Signal = new TF1("Signal", CkSignal, FitMin, FitMax, 1);
    Signal->SetParameter(0, fit_pd->GetParameter(0));
    Signal->SetLineColor(kBlack);

    //hExpCk->GetYaxis()->SetRangeUser(0.0, 1.9);
    //hExpCk->GetXaxis()->SetRangeUser(0, 400);

    //printf(" The extracted source radius is: %.3f +/- %.3f fm\n", fit_pd->GetParameter(0), fit_pd->GetParError(0));
    printf("a: %.3f +/- %.3f\n", fit_pd->GetParameter(1), fit_pd->GetParError(1));
    printf("b: %.3e +/- %.3e 1/MeV\n", fit_pd->GetParameter(2), fit_pd->GetParError(2));
    printf("c: %.3e +/- %.3e\n", fit_pd->GetParameter(3), fit_pd->GetParError(3));
    printf(" The extracted BL slope is: %.3e +/- %.3e 1/MeV\n", fit_pd->GetParameter(2), fit_pd->GetParError(2));
    printf(" chi2/ndf = %.2f/%i = %.2f\n", fit_pd->GetChisquare(), fit_pd->GetNDF(), fit_pd->GetChisquare() / double(fit_pd->GetNDF()));
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double aE = 0.0;
    double bE = 0.0;
    double cE = 0.0;
    double chi = 0.0;

    a = fit_pd->GetParameter(1);
    aE = fit_pd->GetParError(1);

    b = fit_pd->GetParameter(2);
    bE = fit_pd->GetParError(2);

    c = fit_pd->GetParameter(3);
    cE = fit_pd->GetParError(3);

    chi = fit_pd->GetChisquare() / (fit_pd->GetNDF());

    TFile* fOutput = new TFile(OutputFileName, "update");

    hResolution_pd->Write("", TObject::kOverwrite);
    hFeedDown_pL_pp->Write("", TObject::kOverwrite);
    hFeedDown_pXim_pL->Write("", TObject::kOverwrite);

    hExpCk->Write("", TObject::kOverwrite);
    fit_pd->Write("", TObject::kOverwrite);
    Signal->Write("", TObject::kOverwrite);
    Bline->Write("", TObject::kOverwrite);

  //  GetCATS_PL_Usmani_CF("Gauss");
  //  GetCATS_PL_NLO_CF("Gauss");
  //  GetCATS_PL_LO_CF("Gauss");
  //  GetCATS_Ld_Lednicky_CF("Gauss");

    if (BL == "Pol2") {
        MakeCFPlot1("Usmani", "Pol0", hExpCk, fit_pd, Bline, Signal, a,  b,  c, aE,  bE,  cE,  chi);
    }
    if (BL == "Pol1") {
        MakeCFPlot1("High Purity Sample", "Pol1", hExpCk, fit_pd, Bline, Signal, a,  b, aE,  bE, chi);
    }
    delete hResolution_pd;
    delete hFeedDown_pL_pp;
    delete hFeedDown_pXim_pL;
    delete hExpCk;
    delete fit_pd;
    delete Bline;
    delete Signal;
    delete fOutput;

    //don't forget to delete the source ones its not needed
    delete CRS_pd;
    //print out the timer information
    long long ExeTime = TIMER.Stop() / 1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime, strtime, 0, true, 6);
    printf("\nExecution time of Ck_pp_Decomposition for SourceType='%s': %s\n", SourceType.Data(), strtime);
    delete [] strtime;

}
