
bool isNumeric(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}
void findRuns(TString file = "/home/rl/TEST/skimmed2024/results/AnalysisResults.root")
{
  //ROOT::EnableImplicitMT();
  //
  int mRunMuber = 0;
  //
  TFile *f = TFile::Open(file);
  TDirectory *dir = (TDirectory*)f->Get("non-prompt-cascade-task/Zorro");  // your directory name

  TIter next(dir->GetListOfKeys());
  int i = 0;
  int Nloops = 100;
  TKey *key;
  while ((key = (TKey*)next())) {
      TObject *obj = key->ReadObj();
      std::string name =  obj->GetName();
      if(isNumeric(name)) {
        std::cout << name << ",";  //std::endl;
      }
      delete obj;
      if(Nloops !=0 && i > Nloops) break;
      i++;
  }
  std::cout << std::endl;
}
