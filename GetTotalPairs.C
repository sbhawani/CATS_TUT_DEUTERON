

void GetTotalPairs() {

    TString FolderName[4] = {"HMDeuteronDCA0", "HMAntiDeuteronDCA0"};
    TString FolderNameMC[4] = {"HMDeuteronMC0", "HMAntiDeuteronMC0"};
    TString Particle[4] = {"Deuteron", "Antideuteron"};

    TString FileNameData;
    TString FileNameDataMC;

    //FileNameData = "~/cernbox/ProtonDeuteron/AnalysisResultsRootfiles/AODs/ResultsFromAODTrain/Data/FinalCuts/FullpTRange/AnalysisResults.root";
      FileNameData = "/home/sbhawani/cernbox/ProtonDeuteron/AnalysisResultsRootfiles/AODs/ResultsFromAODTrain/Data/FinalCuts/LowPTRange/AnalysisResults.root";

    TFile *FileData = TFile::Open(FileNameData.Data(), "READ");

    TDirectoryFile *MyTaskDirectory = (TDirectoryFile*)(FileData->FindObjectAny("HMResults0"));
    TList* MyTaskDirectory2 = (TList*)(MyTaskDirectory->FindObjectAny("HMResults0"));

    TList* MyTaskDirectory3 = (TList*)(MyTaskDirectory2->FindObject("Particle0_Particle2"));
    TList *MyTaskDirectory4 = (TList*)(MyTaskDirectory2->FindObject("Particle1_Particle3"));

    TH1F  * H1 = (TH1F*)MyTaskDirectory3->FindObject("SEDist_Particle0_Particle2");
    TH1F  * H2 = (TH1F*)MyTaskDirectory4->FindObject("SEDist_Particle1_Particle3");

    TH1F  *h1 = (TH1F*)H1->Clone("h1");
    TH1F  *h2 = (TH1F*)H2->Clone("h2");
    h1->Draw();
    int fbin1 = h1->FindBin(0.001);
    int lbin1 = h1->FindBin(0.20);
    int fbin2 = h2->FindBin(0.001);
    int lbin2 = h2->FindBin(0.20);
    float histoIn1 = h1->Integral(fbin1, lbin1);
    float histoIn2 = h2->Integral(fbin2, lbin2);


    std::cout << "fbin1 = " << fbin1 << std::endl;
    std::cout << "lbin2 = " << lbin2 << std::endl;
    std::cout << "pdPairs = " << histoIn1 << std::endl;
    std::cout << "ApAdPairs = " << histoIn2 << std::endl;
    std::cout << "Total = " << histoIn1 + histoIn2 << std::endl;
    //std::cout << "pd - pair for 3He = " << 0.333*(histoIn1) << std::endl;
   // std::cout << "apad - pair for 3He = " << 0.333*( histoIn2) << std::endl;
  //  std::cout << "pd apad - pair for 3He = " << 0.333*(histoIn1 + histoIn2) << std::endl;

    std::cout << "percentage Loss = " << (histoIn1 + histoIn2-8432)/8432 * 100 << std::endl;
    return;
}
