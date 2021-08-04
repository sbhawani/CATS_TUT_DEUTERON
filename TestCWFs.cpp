#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include "TestCWFs.h"

#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"

using namespace std;

#define SIGN(a) (((a) < 0) ? (-1) : (1))

const double precision = 1E-10, sqrt_precision = 1E-5;

#include "complex_functions.H"
#include "cwfcomp.cpp"
#include "test_rec_rel.cpp"

TGraph* gWF_C_F = new TGraph();

void TestCWFs() {
	TFile* OutputFile = new TFile(
	    TString::Format("/home/sbhawani/Desktop/CATS_TUT_Deuteron/Output_cwfcomp/QA_cwfcomp.root"), "recreate");
	printf("File Created\n");
	cout << "You're welcome here\n";

	string is_it_normalized_str;
	//cin >> is_it_normalized_str;

	const bool is_it_normalized = true ;//(is_it_normalized_str == "true") ? (true) : (false);

	complex<double> l = (1, 0.1);
	complex<double> eta = (50, 50);


	double R = 100.156;

	int Nz, Nl;
	Nz = 100;
	Nl = 3;

	const double step = 2.0 / static_cast<double> (Nz);

	class Coulomb_wave_functions cwf(is_it_normalized, l, eta), cwf_lp1(is_it_normalized, l + 1, eta);

	ofstream out_file("testComplexCoulomb.output");
	out_file.precision (10);
	cout << "Reached-1\n";
	if (is_it_normalized) out_file << "is_it_normalized:true l:" << l << " eta:" << eta << " R:" << R << "  Nz:" << Nz << "  Nl:" << Nl << endl << endl;
	if (!is_it_normalized) out_file << "is_it_normalized:false l:" << l << " eta:" << eta << " R:" << R << " Nz:" << Nz << "  Nl:" << Nl << endl << endl;
	cout << "Reached-2\n";
	for (int i = 0 ; i < Nz ; i++) {
		unsigned COUNTER = 0;
		const double x = i * step;
		const complex<double> z = R * exp (complex<double> (0, x * M_PI));
		complex<double> F, dF, Hp, dHp, Hm, dHm, G, dG;
		cwf.F_dF (z, F, dF);
		cwf.G_dG (z, G, dG);
		cwf.H_dH (1, z, Hp, dHp);
		cwf.H_dH (-1, z, Hm, dHm);

		out_file << "z:" << z << endl;
		out_file << "F:" << F << " F':" << dF << endl;
		out_file << "G:" << G << " G':" << dG << endl;
		out_file << "H+:" << Hp << " H+':" << dHp << endl;
		out_file << "H-:" << Hm << " H-':" << dHm << endl;
		//out_file << "Wronskian test: " << Wronskian_test (z, cwf, cwf_lp1) << endl << endl;
		cout << "Reached-3\n";

		//gWF_C_F->SetPoint(COUNTER, z, F);
	}
	gWF_C_F->Write();
	delete  gWF_C_F;
	delete OutputFile;

	/* TGraph* gWF_Re_D = new TGraph [NumMomBins];
	 TGraph* gWF_Re_Q = new TGraph [NumMomBins];
	 TGraph* gWF_Re_R = new TGraph [NumMomBins];
	 for(unsigned uBin=0; uBin<NumMomBins; uBin++){
	     unsigned COUNTER=0;
	     gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f",Kitty.GetMomentum(uBin)));
	     gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f",Kitty.GetMomentum(uBin)));
	     gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f",Kitty.GetMomentum(uBin)));
	     for(double RAD=0.1; RAD<95; RAD+=0.05){
	         gWF_Re_D[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,0,0,RAD,true)));
	         gWF_Re_Q[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,1,0,RAD,true)));
	         gWF_Re_R[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalReferenceRadialWF(uBin,0,RAD,true)));
	         COUNTER++;
	     }
	     gWF_Re_D[uBin].Write();
	     gWF_Re_Q[uBin].Write();
	     gWF_Re_R[uBin].Write();
	 }
	 CleanUpWfHisto(Kitty,ExternalWF);
	 delete [] gWF_Re_D;
	 delete [] gWF_Re_Q;
	 delete [] gWF_Re_R;
	 delete OutputFile;
	 delete cPars;*/
}


