//#include "FFT.h"
int run()
{
  //Compile Macro
  Int_t workedornot = gSystem->CompileMacro("eff3.C","-kfo");
  cout<<"--------------------------------------------------"<<endl;
  cout<<endl;
  if( workedornot == 0 ){
    cout<<"********************************"<<endl;
    cout<<" eff3.C compilation failed! "<<endl;
    cout<<"********************************"<<endl;
    return 1;
  }
  //Load Class
  int ret=gSystem->Load("eff3_C.so");
  cout << "Load ret: " << ret << endl;
  //FFT *v0 = new FFT;
  //return;
  //v0->FFT();
}
