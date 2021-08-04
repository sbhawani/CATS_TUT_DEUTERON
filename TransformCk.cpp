
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"

#include "Basics.h"
#include "ExtendedCk.h"
#include "TransformCk.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TSystem.h"

#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"
#include "DLM_Random.h"
#include "DLM_CppTools.h"

#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

TH1F *hCkNLO_Smear = new TH1F("hCkNLO_Smear", "hCkNLO", 100, 0, 400);
TH1F *hCkNLO_Smeared = new TH1F("hCkNLO_Smeared", "hCkNLO_Tranformed", 100, 0, 400);

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


void Convert_Histo(TGraph* graph, TH1F* h) {
	auto nPoints = graph->GetN(); // number of points in your TGraph
	for (int i = 0; i < nPoints; ++i) {
		double x, y;
		graph->GetPoint(i, x, y);
		h->Fill(x, y); // ?
	}
	//return h;
}


/*void Smear(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared) {
	if (!SmearMatrix) {
		CkSmeared[0] = CkToSmear[0];
	} else {
		for (unsigned uBinSmear = 0; uBinSmear < CkToSmear->GetNbins(); uBinSmear++) {
			CkSmeared->SetBinContent(uBinSmear, 0);
			unsigned FirstBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisFirst];
			unsigned LastBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisLast];
			if (LastBin >= CkToSmear->GetNbins()) LastBin = CkToSmear->GetNbins() - 1;
			for (unsigned uBinTrue = FirstBin; uBinTrue <= LastBin; uBinTrue++) {
				//as the response matrix is normalized to the size of the bin, during the integration we multiply for it
				CkSmeared->Add(uBinSmear,   SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue)*
				               CkToSmear->GetBinSize(uBinTrue)*CkSmeared->GetBinSize(uBinSmear));
			}
		}
	}
}*/

void TRansform(const TH1F* CkToSmear,  const DLM_ResponseMatrix* SmearMatrix, TH1F* CkSmeared) {
	if (!SmearMatrix) {
		CkSmeared[0] = CkToSmear[0];
	} else {
		for (unsigned uBinSmear = 0; uBinSmear < CkToSmear->GetNbinsX(); uBinSmear++) {
			CkSmeared->SetBinContent(uBinSmear, 0);
			unsigned FirstBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisFirst];
			unsigned LastBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisLast];
			//unsigned FirstBin = SmearMatrix->GetXaxis()->();//first
			//unsigned LastBin = SmearMatrix->GetXaxis()->();//last;
			if (LastBin >= CkToSmear->GetNbinsX()) LastBin = CkToSmear->GetNbinsX() - 1;
			for (unsigned uBinTrue = FirstBin; uBinTrue <= LastBin; uBinTrue++) {
				//as the response matrix is normalized to the size of the bin, during the integration we multiply for it
				CkSmeared->Add(uBinSmear,  SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue)*
				               CkToSmear->GetBinWidth(uBinTrue)*CkSmeared->GetBinWidth(uBinSmear));
			}
		}
	}
}
/*void Transform(const TH1F* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, TH1F* CkSmeared) {//,const TH2F* SmearMatrix
	if (!SmearMatrix) {
		CkSmeared[0] = CkToSmear[0];
	} else {
		for (unsigned uBinSmear = 0; uBinSmear < CkToSmear->GetNbinsX(); uBinSmear++) {
			CkSmeared->SetBinContent(uBinSmear, 0);
			unsigned FirstBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisFirst];
			unsigned LastBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisLast];
			if (LastBin >= CkToSmear->GetNbinsX()) LastBin = CkToSmear->GetNbinsX() - 1;
			for (unsigned uBinTrue = FirstBin; uBinTrue <= LastBin; uBinTrue++) {
				//as the response matrix is normalized to the size of the bin, during the integration we multiply for it
				CkSmeared->Add(uBinSmear,   SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue)*
				               CkToSmear->GetSize(uBinTrue)*CkSmeared->GetSize(uBinSmear));
			}
		}
	}
}*/
void TransformCk() {

// tranforamation matrix

	TString FileName_Feed_pL_pp = "../Files/Decay_matrix_dp_dL.root";
	TString HistoName_Feed_pL_pp = "hRes_dp_dL";
	TH2F* hTMatrix = GetSmearMatrix(FileName_Feed_pL_pp, HistoName_Feed_pL_pp);

	TString FileNameData = "/home/sbhawani/Desktop/CATS_TUT_Deuteron/OutputFiles/Ck_pL_NLO_needed_forPd_Gauss.root";
	TFile* FileData = TFile::Open(FileNameData.Data(), "READ"); // number 3

	TGraph * Ck_NLO1;
	TString Ck[3] = {"Ck_pL_NLO", "Ck_pL_LO", "Ck_pL_Usmani"};
	Ck_NLO1 = (TGraph*)FileData->FindObjectAny(Ck[0].Data());
	TGraph* Ck_NLO = (TGraph*)Ck_NLO1->Clone("Ck_NLO");
	//Ck_NLO->Draw();
	DLM_Histo<double>* CkSmeared1;

	Convert_Histo(Ck_NLO, hCkNLO_Smear);// ->Draw();
	hCkNLO_Smear->Draw("HIST");
	TRansform(hCkNLO_Smear, hTMatrix, CkSmeared1);
	CkSmeared1->Draw("HIST");
}

