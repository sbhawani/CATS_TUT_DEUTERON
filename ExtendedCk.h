
class TString;
class DLM_Ck;
class DLM_CkDecomposition;

void Ck_pL_Ledni_Usmani();
//save in a file the correlation function from DLM_Ck
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot);
//save in a file the correlation function from DLM_CkDecomposition
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_CkDecomposition* CkToPlot, bool PlotIndividualContributions=false);

void Ck_pd_Decomposition(const TString& SourceType);
void PlotCF(const TString& GraphName, DLM_Ck* CkToPlot);
void Ck_dL_Ledni();
void Ck_pL_Ledni_QA1();
void Ck_pL_Ledni_QA2();
void Ck_LD_Ledni_QA2(); 
void Ck_PD_Ledni_QA2();

void Ck_np_Ledni_Usmani();
void Ck_pd_Lenicky_fromSebastian();