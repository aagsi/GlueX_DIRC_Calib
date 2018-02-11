// macro to process analysis TTree with DSelector
#include <iostream> 

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"


////////////////////
// proto types//////
////////////////////
// check file existance
bool exists_test (const std::string& name);

void runDselector(){


TString sample = Form("/data.local/dirc/halld/analysis/justin_1/GlueX_DIRC_Calib/justin_%d.root", 1);
cout<<"sample data path= " <<sample<<endl;
string path_sample = (string)sample;
cout<<"exist sample)" <<exists_test(path_sample)<<endl;
if (!exists_test(path_sample)) cout<<"sample not found "<<endl;;

	// Load DSelector library

	int proof_Nthreads = 4;
	DPROOFLiteManager::Process_Tree("sample", "pimkpks__B3_M16_Tree", "/data.local/dirc/halld/analysis/justin_1/GlueX_DIRC_Calib/DSelector_justin_1_analyzer.C+", 4);

}



//////////////////////////
// check file existance //
//////////////////////////
bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

