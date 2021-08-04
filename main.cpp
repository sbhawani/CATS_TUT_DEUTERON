#include<iostream>
#include "Basics.h"
#include "ProtonDeuteronCF_fromWFModel.h"
#include "examples.h"
#include "ExtendedCk.h"
#include "TransformCk.h"
#include "examples.h"
#include "TestCWFs.h"

#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
void ComparePionPion() {
    TGraph* grCATS = Basics_PiPiCATS(1, 1);
    TGraph* grTheory = Basics_PiPiTheory(1, 1);
    TFile* OutputFile = new TFile("ComparePionPion.root", "recreate");
    grCATS->Write();
    grTheory->Write();

    delete grCATS;
    delete grTheory;
    delete OutputFile;
}

int main(int argc, char *argv[]) {

    printf("\nHello to this lovely CATS tutorial!\n");
    printf(" To find a bug: continue with this tutorial\n");
    printf("-------------------------------------------\n");
  ///  ProtonDeuteronCF_CoulombvsRadii_CATS();
    //ComparePionPion();
    //Basics_ProtonLambda();
    ProtonDeuteronCF_fromWFModelNLO();
    ProtonDeuteronCF_fromWFModel();
    //Ck_pL_Ledni_Usmani();
    //Ck_pd_Decomposition("Gauss");
  // ProtonDeuteronCF_fromWFModel();
   // NeuteronDeuteronCF_fromWFModel();
    //ProtonDeuteronCF_Coulomb_fromWFModel();
   // ProtonDeuteronCF_Coulomb_CATS();
    //ProtonDeuteron_PhaseShift_fromWFModel();
    //Ck_pd_Lenicky_fromSebastian();
   // Ck_np_Ledni_Usmani();
   // ProtonProtonCF_Sebastian_fromWFModel();
   // ProtonProtonAV18CF_Sebastian_fromWFModel();
   //ProtonProtonPionLessCF_Sebastian_fromWFModel();
   // NeuteronProtonAV18CF_Sebastian_fromWFModel();
    //Ck_pp_CATS();//Ck_pp_CATS
   // Ck_pp_CATS_CoulombOnly();
   //// ProtonProtonAV18_CATS_WF();// using WF from CATS
    //ProtonDeuteronCF_Coulomb_CATS();
    //Ck_pL_Ledni_QA1();
    //Ck_pL_Ledni_QA2();
   // Ck_LD_Ledni_QA2();
   // Ck_PD_Ledni_QA2();
    //Ck_pd_Decomposition("CoreReso");
    //TransformCk();
   // Test Complex Coulomb Wave Function"
    //TestCWFs();
    return 0;
}
