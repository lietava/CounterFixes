#include <iostream>

#include "TFile.h"
#include "TObjString.h"
#include "TString.h"

void printMatrixHOptions(const char* fileName = "manualBayes4D_useAO2Data3_cent0_100_pileup0.root")
{
  TFile file(fileName, "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    return;
  }

  auto* options = dynamic_cast<TObjString*>(file.Get("matrixH_options"));
  if (!options) {
    std::cerr << "Object matrixH_options not found in " << fileName << std::endl;
    return;
  }

  std::cout << options->String().Data() << std::endl;
}
