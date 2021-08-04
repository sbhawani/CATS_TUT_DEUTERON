
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"


TMultiGraph *g0 = new TMultiGraph();
TGraph* gUsmani = new TGraph();
TGraph* gLO = new TGraph();
TGraph* gLDGenuine = new TGraph();
TGraph* gNLO = new TGraph();
TGraph* gLD = new TGraph();
TH1F *hCkNLO_Smear = new TH1F("hCkNLO_Smear", "hCkNLO", 100, 0, 400);
TH1F *hCkUS_Smeared = new TH1F("hCkUS_Smeared", "hCkUS_Tranformed", 100, 0, 400);
TH1F *hCkLO_Smeared = new TH1F("hCkLO_Smeared", "hCkLO_Tranformed", 100, 0, 400);
TH1F *hCkLO_Genuine = new TH1F("hCkLO_Genuine", "hCkLDGenuine", 100, 0, 400);
TH1F *hCkNLO_Smeared = new TH1F("hCkNLO_Smeared", "hCkNLO_Tranformed", 100, 0, 400);
TH1F *hCkLD_Smeared = new TH1F("hCkLD_Smeared", "hCkLD_Tranformed", 100, 0, 400);
void SaveRootFile(TString Name, TGraph *Histo) {
    TFile *myfile = TFile::Open(Name.Data(), "RECREATE");
    Histo->Write();
}

void Convert_Histo(TGraph* graph, TH1F* h, double Lambda, TGraph* graph1) {
	auto nPoints = graph->GetN(); // number of points in your TGraph
	for (int i = 0; i < nPoints; ++i) {
		double x, y;
		graph->GetPoint(i, x, y);
		graph1->SetPoint(i, x, y /Lambda + 1);
		h->Fill(x, y/Lambda + 1); // ?
	}
	//return h;
}

void LabelCFgraph(TGraph *Histo, TString Title, EColor color) {
	//Histo->SetMarkerColor(color);
	Histo->SetLineColor(color);
	Histo->SetLineWidth(2);
	//Histo->SetMarkerStyle(20);
	//Histo->SetTitle(Title.Data());
	//Histo->GetXaxis()->SetTitle("k* (MeV/c)");
	//Histo->GetYaxis()->SetTitle("C(k*)");
	//Histo->GetXaxis()->SetTitleOffset(0.9);
	//Histo->GetYaxis()->SetTitleOffset(0.8);
}
void LabelCFHisto(TH1F *h1, TString Title, EColor color, double size) {
	//h1->SetMarkerColor(color);
	h1->SetLineColor(color);
	h1->SetLineWidth(2);
	//h1->SetMarkerSize(size);
	//h1->SetMarkerStyle(20);
	h1->SetTitle(Title.Data());
	//h1->GetXaxis()->SetRangeUser(0, 400);
	//h1->GetXaxis()->SetTitleOffset(0.9);
	//h1->GetYaxis()->SetTitleOffset(0.8);
}
void MakeCFPlot(TString Title, TGraph *g1, TGraph *g2, TGraph *g3, TGraph *g4, TMultiGraph*g) {

	LabelCFgraph(g1, "Usmani", kBlack);
	LabelCFgraph(g2, "LO", kRed);
	LabelCFgraph(g3, "NLO", kGreen);
	LabelCFgraph(g4, "From Haidenbuar #lambda D", kBlue);

	TCanvas *CFPlot = new TCanvas("CFPlot", "CFplot", 0, 0, 800, 600);
	//TString pdfCFName = "RightHandHaidenBaurFromCATS.png";
	TString pdfCFName = "RightHandHaidenBaurTransformed.png";


	TLegend *leg = new TLegend( 0.55, 0.50, 0.75, 0.69);
	leg->SetTextSize(0.040);
	leg->AddEntry(g1, "p#Lambda Usmani");
	leg->AddEntry(g2, "p#Lambda LO #chi EFT ");
	leg->AddEntry(g3, "p#Lambda NLO #chi EFT ");
	leg->AddEntry(g4, "Haidenbuar #Lambdad (RHS)");
	leg->SetLineColor(0);

	TLegend *leggUNO = new TLegend(0.50, 0.7, 0.75, 0.81);
	leggUNO->SetTextSize(0.04);
	const char * textgUNO    = " Source size = 1.2 fm";
	leggUNO->AddEntry((TObject*)0, textgUNO, "");
	leggUNO->SetLineColor(0);

	g->Add(g1);
	g->Add(g2);
	g->Add(g3);
	g->Add(g4);
	g->SetTitle(Title);
	g->GetXaxis()->SetTitle("#it{k*} (GeV/c)");
	g->GetYaxis()->SetTitle("#it{C}(k*)");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->SetTitleOffset(0.8);
	g->GetXaxis()->SetRangeUser(0, 400);
	g->GetYaxis()->SetRangeUser(0.8, 1.7);
	//g->GetYaxis()->SetRangeUser(0.0, 40);
	g->SetTitle("Tranformed CF");
	g->Draw("ACP");
	leg->Draw("same");
	leggUNO->Draw("same");
	CFPlot->Print(pdfCFName);
	delete CFPlot;
	return;
}
void MakeCFPlot0(TString Title, TGraph *g) {



	LabelCFgraph(g, "Genuine #lambda D", kBlue);

	TCanvas *CFPlot = new TCanvas("CFPlot", "CFplot", 0, 0, 800, 600);
	//TString pdfCFName = "RightHandHaidenBaurFromCATS.png";
	TString pdfCFName = "AfterFitGenuineLD_CF.pdf";


	TLegend *leg = new TLegend( 0.55, 0.50, 0.75, 0.69);
	leg->SetTextSize(0.040);


	leg->AddEntry(g, "Genuine #Lambda d");
	leg->SetLineColor(0);

	//TLegend *leggUNO = new TLegend(0.50, 0.7, 0.75, 0.81);
	//leggUNO->SetTextSize(0.04);
	//const char * textgUNO    = " Source size = 0.9 fm";
	//leggUNO->AddEntry((TObject*)0, textgUNO, "");
	//leggUNO->SetLineColor(0);

	g->SetTitle(Title);
	g->GetXaxis()->SetTitle("#it{k*} (GeV/c)");
	g->GetYaxis()->SetTitle("#it{C}(k*)");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->SetTitleOffset(0.8);
	g->GetXaxis()->SetRangeUser(0, 400);
	g->GetYaxis()->SetRangeUser(0.8, 40.01);
	//g->GetYaxis()->SetRangeUser(0.0, 40);
	g->SetTitle("Genuine #Lambda d ");
	g->Draw("ACP");
	leg->Draw("same");
	//leggUNO->Draw("same");
	CFPlot->Print(pdfCFName);
	SaveRootFile("Genuine4MeV.root", g);
	delete CFPlot;
	return;
}
void MakeCFPlot(TString Title, TGraph *g2, TGraph *g4, TMultiGraph *g) {


	LabelCFgraph(g2, "LO", kRed);
	LabelCFgraph(g4, "From Haidenbuar #lambda D", kBlue);

	TCanvas *CFPlot = new TCanvas("CFPlot", "CFplot", 0, 0, 800, 600);
	//TString pdfCFName = "RightHandHaidenBaurFromCATS.png";
	TString pdfCFName = "RightHandHaidenBaurTransformedComparision.png";


	TLegend *leg = new TLegend( 0.55, 0.50, 0.75, 0.69);
	leg->SetTextSize(0.040);

	leg->AddEntry(g2, "Haidenbuaer #Lambdad (LHS)");
	leg->AddEntry(g4, "#Lambdad From FIT");
	leg->SetLineColor(0);

	TLegend *leggUNO = new TLegend(0.50, 0.7, 0.75, 0.81);
	leggUNO->SetTextSize(0.04);
	const char * textgUNO    = " Source size = 1.2 fm";
	leggUNO->AddEntry((TObject*)0, textgUNO, "");
	leggUNO->SetLineColor(0);

	g->Add(g2);
	g->Add(g4);
	g->SetTitle(Title);
	g->GetXaxis()->SetTitle("#it{k*} (GeV/c)");
	g->GetYaxis()->SetTitle("#it{C}(k*)");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->SetTitleOffset(0.8);
	g->GetXaxis()->SetRangeUser(0, 400);
	g->GetYaxis()->SetRangeUser(0.8, 1.7);
	//g->GetYaxis()->SetRangeUser(0.0, 40);
	g->SetTitle("Tranformed ComparisionCF");
	g->Draw("ACP");
	leg->Draw("same");
	//leggUNO->Draw("same");
	CFPlot->Print(pdfCFName);
	delete CFPlot;
	return;
}
void MakeCFPlot1(TString Title, TH1F*h, TH1F *h1, TH1F *g1, TH1F *g2) {

	LabelCFHisto(g1, "", kBlack, 0);
	LabelCFHisto(g2, "", kRed, 0);
	LabelCFHisto(h1, "", kGreen, 0);
	LabelCFHisto(h, "data", kBlue, 1);

	TCanvas *CFPlot1 = new TCanvas("CFPlot1", "CFplot1", 0, 0, 800, 600);
	//TString pdfCFName = "FromCATS.png";
	TString pdfCFName = "FitResultRightHandsideHaidenbaur" + Title + ".png";
	//TString pdfCFName = "FitResults_NLO_LowR.png";
	//TString pdfCFName = "FitResults_Usmani_LowR.png";
	TLegend *leg = new TLegend( 0.45, 0.4, 0.75, 0.55);
	leg->SetTextSize(0.035);
	leg->AddEntry(g1, "BaseLine pol1");
	leg->AddEntry(g2, "Modeled Ck");
	leg->AddEntry(h1, "Fit");
	leg->AddEntry(h, "Data");
	leg->SetLineColor(0);

	TLegend *leggUNO = new TLegend(0.35, 0.28, 0.75, 0.35);
	leggUNO->SetTextSize(0.03);
	const char * textgUNO    = "Source size (Gauss) = 0.92 fm";
	//const char * textgUNO1    = "Feed-down Usa";
	const char * textgUNO2    = "BL norm is: 1.021 +/- 0.031";
	const char * textgUNO3    = "BL slope is: 7.991e-05 +/- 1.084e-04 1/MeV";
	const char * textgUNO4    = "#chi^{2}/ndf 80.08/15 = 5.34";
	leggUNO->AddEntry((TObject*)0, textgUNO, "");
	//leggUNO->AddEntry((TObject*)0, textgUNO1, "");
	//leggUNO->AddEntry((TObject*)0, textgUNO2, "");
	//leggUNO->AddEntry((TObject*)0, textgUNO3, "");
	leggUNO->AddEntry((TObject*)0, textgUNO4, "");
	leggUNO->SetLineColor(0);

	/*g->Add(g1);
	//g->Add(g3);
	g->Add(g2);
	*/
	h->SetTitle(Title);
	h->GetXaxis()->SetTitle("#it{k*} (GeV/c)");
	h->GetYaxis()->SetTitle("#it{C}(k*)");
	h-> GetXaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(0.9);
	h->GetYaxis()->SetTitleOffset(0.8);
	h->GetXaxis()->SetRangeUser(0, 432);
	h->GetYaxis()->SetRangeUser(-0.1, 1.4);

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

void TransformCk() {

	TCanvas *binplot = new TCanvas ("binplot", "CF", 1000, 700);
	binplot->SetLeftMargin(0.1505);
	binplot->SetRightMargin(0.035);
	gStyle->SetOptStat(0);

	TString FileNameDataGenuineLD;
	TGraph * Ck_LD_Genuine1;
	FileNameDataGenuineLD = "/home/sbhawani/Desktop/CATS_Tutorial/OutputFiles/FinralFitResultsCoulombOnlyGauss.root";
	TFile* FileGenuineLD = TFile::Open(FileNameDataGenuineLD.Data(), "READ"); // number 1
	Ck_LD_Genuine1 = (TGraph*)FileGenuineLD->FindObjectAny("CorrCkDec_pd_child2");
	TGraph* Ck_LD_Genuine = (TGraph*)Ck_LD_Genuine1->Clone("Ck_pPion_Genuine");

	Convert_Histo(Ck_LD_Genuine, hCkLO_Smeared, 0.6, gLDGenuine);
	MakeCFPlot0("Genuine P_{#Lambda}d", gLDGenuine);



	bool UseFromCATS = false;
	bool PlotsPriorFit = true;
	bool PlotFitResutls = false;

	if (PlotsPriorFit) {

		TString FileNameData;
		TString FileNameData1;
		TString FileNameData2;
		TString FileNameData3;
		TString Ck[4] = {"Ck_pL_Usmani", "Ck_pL_LO", "Ck_pL_NLO", "Ck_Ld_lednicky"};
		TString Ck1[3] = {"CorrCkDec_pd_child0", "CorrCkDec_pd_child0", "CorrCkDec_pd_child0"};

		TGraph * Ck_US1;
		TGraph * Ck_LO1;
		TGraph * Ck_NLO1;
		TGraph * Ck_LD1;
		if (UseFromCATS) {
			FileNameData = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pL_needed_forPd_LowRGauss.root";
			FileNameData1 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pL_LO_needed_forPd_LowRGauss.root";
			FileNameData2 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pL_NLO_needed_forPd_LowRGauss.root";
			//FileNameData3 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_Ld_Lednicky_needed_forPd_LowRGauss.root";
			FileNameData3 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_Ld_LednickyRightHandHaidenabur_needed_forPd_LowRGauss.root";
		} else {
			FileNameData = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pd_DecompositionUsmani_LowRGauss.root";
			//FileNameData1 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pd_DecompositionLO_LowRGauss.root";
			FileNameData2 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pd_DecompositionNLO_LowRGauss.root";
			FileNameData1 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/LambdaDeuteronCkUsedGauss.root";
			FileNameData3 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/LambdaDeuteronRightHandSideHaidenbaurCkUsedGauss.root";
		}

		TFile* FileData = TFile::Open(FileNameData.Data(), "READ"); // number 1
		TFile* FileData1 = TFile::Open(FileNameData1.Data(), "READ"); // number 2
		TFile* FileData2 = TFile::Open(FileNameData2.Data(), "READ"); // number 3
		TFile* FileData3 = TFile::Open(FileNameData3.Data(), "READ"); // number 4

		if (UseFromCATS) {
			Ck_US1 = (TGraph*)FileData->FindObjectAny(Ck[0].Data());
			Ck_LO1 = (TGraph*)FileData1->FindObjectAny(Ck[1].Data());
			Ck_NLO1 = (TGraph*)FileData2->FindObjectAny(Ck[2].Data());
			Ck_LD1 = (TGraph*)FileData3->FindObjectAny(Ck[3].Data());
		} else {
			Ck_US1 = (TGraph*)FileData->FindObjectAny(Ck1[0].Data());
			Ck_LO1 = (TGraph*)FileData1->FindObjectAny(Ck1[1].Data());
			Ck_NLO1 = (TGraph*)FileData2->FindObjectAny(Ck1[2].Data());
			Ck_LD1 = (TGraph*)FileData3->FindObjectAny(Ck1[2].Data());
		}

		TGraph* Ck_US = (TGraph*)Ck_US1->Clone("Ck_US");
		TGraph* Ck_LO = (TGraph*)Ck_LO1->Clone("Ck_LO");
		TGraph* Ck_NLO = (TGraph*)Ck_NLO1->Clone("Ck_NLO");
		TGraph* Ck_LD = (TGraph*)Ck_LD1->Clone("Ck_NLO");

		if (UseFromCATS) {
			MakeCFPlot(" From CATS ", Ck_US, Ck_LO, Ck_NLO, Ck_LD, g0);
		} else {
			//Convert_Histo(Ck_US, hCkUS_Smeared, 0.084, gUsmani);
			Convert_Histo(Ck_LO, hCkLO_Smeared, 0.084, gLO);
			//Convert_Histo(Ck_NLO, hCkNLO_Smeared, 0.084, gNLO);
			Convert_Histo(Ck_LD, hCkLD_Smeared, 0.084, gLD);
			//MakeCFPlot("Transformed into P_{#Lambda}d", gUsmani, gLO, gNLO , gLD , g0);
			MakeCFPlot("Transformed into P_{#Lambda}d", gLO, gLD , g0);
		}
		//MakeCFPlot(gUsmani, gLO, gNLO , g0);
		//hCkLO_Smeared->Draw("HIST");
		//hCkNLO_Smeared->Draw("HIST");
		//hCkUS_Smeared->Draw("HIST");
		//Ck_NLO->Draw(" ");
		//delete FileData;
		//delete FileData1;
		//delete FileData2;
		//delete Ck_US;
		//delete Ck_LO;
		//delete Ck_NLO;

	}
///
	if (PlotFitResutls) {

		bool NLO = false;
		bool LO = false;
		bool Gauss = true;
		bool LD = true;
		TString FileNameData2;
		TString FileNameData21;
		TString FileNameData22;
		TString FileNameData23;
		TString Ck[4] = {"Bline", "Signal", "hCk_ReweightedMeV_0", "fit_pd"};


		TH1F*  Ck_Coulmb1;// = new TGraph();
		TH1F* Ck_CorrCk1;// = new TGraph();
		TH1F* fit_pp1;// = new TGraph();
		TH1F * hCk1;

		if (Gauss) {
			FileNameData2 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Usmani_Gauss.root";
			FileNameData21 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/LO_Gauss.root";
			FileNameData22 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/NLO_Gauss.root";
			FileNameData23 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/LambdaDeuteronCkUsedGauss.root";
		} else {
			FileNameData2 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Usmani_CoreReso.root";
			FileNameData21 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/LO_CoreReso.root";
			FileNameData22 = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/NLO_CoreReso.root";
		}

		TFile* FileData2 = TFile::Open(FileNameData2.Data(), "READ"); // number 1
		TFile* FileData21 = TFile::Open(FileNameData21.Data(), "READ"); // number 2
		TFile* FileData22 = TFile::Open(FileNameData22.Data(), "READ"); // number 3
		TFile* FileData23 = TFile::Open(FileNameData23.Data(), "READ"); // number 3


		if (LO) {
			Ck_Coulmb1 = (TH1F*)FileData21->FindObjectAny(Ck[0].Data());
			Ck_CorrCk1 = (TH1F*)FileData21->FindObjectAny(Ck[1].Data());
			hCk1 = (TH1F*)FileData21->FindObjectAny(Ck[2].Data());
			fit_pp1 = (TH1F*)FileData21->FindObjectAny(Ck[3].Data());
			std::cout << "Plotting for LO\n";
		}  else if (NLO) {
			Ck_Coulmb1 = (TH1F*)FileData22->FindObjectAny(Ck[0].Data());
			Ck_CorrCk1 = (TH1F*)FileData22->FindObjectAny(Ck[1].Data());
			hCk1 = (TH1F*)FileData22->FindObjectAny(Ck[2].Data());
			fit_pp1 = (TH1F*)FileData22->FindObjectAny(Ck[3].Data());
			std::cout << "Plotting for NLO\n";
		} else if (LD) {
			Ck_Coulmb1 = (TH1F*)FileData23->FindObjectAny(Ck[0].Data());
			Ck_CorrCk1 = (TH1F*)FileData23->FindObjectAny(Ck[1].Data());
			hCk1 = (TH1F*)FileData23->FindObjectAny(Ck[2].Data());
			fit_pp1 = (TH1F*)FileData23->FindObjectAny(Ck[3].Data());
			std::cout << "Plotting for LD HaidenBaur\n";
		} else {

			Ck_Coulmb1 = (TH1F*)FileData2->FindObjectAny(Ck[0].Data());
			Ck_CorrCk1 = (TH1F*)FileData2->FindObjectAny(Ck[1].Data());
			hCk1 = (TH1F*)FileData2->FindObjectAny(Ck[2].Data());
			fit_pp1 = (TH1F*)FileData2->FindObjectAny(Ck[3].Data());
			std::cout << "Plotting for Usmani\n";
		}

		TH1F* Ck_Coulmb = (TH1F*)Ck_Coulmb1->Clone("Ck_Coulmb");
		TH1F* Ck_CorrCk = (TH1F*)Ck_CorrCk1->Clone("Ck_CorrCk");
		TH1F* fit_pp = (TH1F*)fit_pp1->Clone("fit_pp");
		TH1F* hCk = (TH1F*)hCk1->Clone("hCk");

		//MakeCFPlot("Fit Results", Ck_Coulmb, fit_pp, Ck_CorrCk, g0);
		if (LO) {
			MakeCFPlot1("LO_chiEFT", hCk, fit_pp, Ck_Coulmb, Ck_CorrCk);
		} else if (NLO) {
			MakeCFPlot1("NLO_chiEFT", hCk, fit_pp, Ck_Coulmb, Ck_CorrCk);
		}  else if (LD) {
			MakeCFPlot1("Lambda_d From Haidenbuar", hCk, fit_pp, Ck_Coulmb, Ck_CorrCk);
			//Convert_Histo(Ck_LD, hCkLD_Smeared, 0.084, gLD);
		} else {
			MakeCFPlot1("#Lambda_d From Haidenbuar", hCk, fit_pp, Ck_Coulmb, Ck_CorrCk);
		}

		delete Ck_Coulmb;
		delete Ck_CorrCk;
		delete fit_pp;
		delete hCk;

	}

	delete binplot;
}