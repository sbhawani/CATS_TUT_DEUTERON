
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"

//classes to handle Ck in general
#include "DLM_CkDecomposition.h"
//pre-defined potential functions
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
void ProtonDeuteronCF_fromWFModelNLO();
void ProtonDeuteronCF_fromWFModel();
void NeuteronDeuteronCF_fromWFModel();
void ProtonDeuteronCF_Coulomb_fromWFModel();
void ProtonDeuteron_PhaseShift_fromWFModel();
void ProtonDeuteronCF_Coulomb_CATS();
void ProtonProtonCF_Sebastian_fromWFModel();
void ProtonProtonAV18CF_Sebastian_fromWFModel();
void NeuteronProtonAV18CF_Sebastian_fromWFModel();
void Ck_pp_CATS();
void CATS_GaussSource1(CATS& Kitty, const double& SourceSize);
void CATS_pp_AV181(CATS& Kitty, const bool& pwaves, const bool& dwaves);
void CATS_pp_AV181_SwaveOnly(CATS& Kitty, const bool& pwaves, const bool& dwaves);
void CATS_pp_Basic1(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax);
void LabelDCAHisto(TH1F *Histo, TString Title, EColor color);
void MakePlotforPhase(TString Name, TH1F *h1, TH1F *h2, TH1F *h3, int Cutoff, double photonmass, double xRangeLow, double xRangeHigh, double yRangeLow, double yRangeHigh);
void ProtonProtonAV18_CATS_WF();
void ProtonProtonPionLessCF_Sebastian_fromWFModel();
void Ck_pp_CATS_CoulombOnly();
void ProtonDeuteronCF_CoulombvsRadii_CATS();
