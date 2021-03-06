// macro to process analysis TTree with DSelector
// root Load_DSelector.C runDSelector.C
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

void runDSelector(bool proof = 1){
    string path ="/data.local/dirc/halld/analysis/justin_1/GlueX_DIRC_Calib/";
    string SampleName = path;
    string DSelectorName = path;
    string TreeName ="pimkpks__B3_M16_Tree";
    string ConfigName = path;
    
    SampleName += Form("data_%d.root", 1);
    DSelectorName += Form("DSelector_justin_%d_analyzer.C+", 1);
    ConfigName += Form("config_%d.in", 1);
    
    cout<<"########### Sample used= " <<SampleName<<endl;
    cout<<"########### DSelector used= " <<DSelectorName<<endl;
    cout<<"########### Config file used= " <<ConfigName<<endl;
    cout<<"exist Sample)" <<exists_test(SampleName)<<endl;
    cout<<"exist DSelector)" <<exists_test(DSelectorName)<<endl;
    cout<<"exist ConfigName)" <<exists_test(ConfigName)<<endl;
    if (!exists_test(SampleName)) cout<<"Sample not found "<<endl;
    if (!exists_test(DSelectorName)) cout<<"DSelector not found "<<endl;
    if (!exists_test(ConfigName)) cout<<"Config file not found "<<endl;
    
    
    int proof_Nthreads = 4;
    string outputHistFileName = "hist_data_ks.root";
    string outputTreeFileName = "tree_data_ks.root";
    DPROOFLiteManager::Process_Tree( SampleName, TreeName, DSelectorName, proof_Nthreads, outputHistFileName,  outputTreeFileName, ConfigName);
    
    
    //    if(false) { // add TTree to chain and use PROOFLiteManager
    //        TChain *chain = new TChain("pimkpks__B3_M16_Tree");
    //        chain->Add("/data.local/dirc/halld/analysis/justin_1/GlueX_DIRC_Calib/justin_1.root");
    //        string outputHistFileName = "hist_ks.root";
    //        string outputTreeFileName = "tree_ks.root";
    //        DPROOFLiteManager::Process_Chain(chain, DSelectorName, outputHistFileName, outputTreeFileName, SampleName, proof_Nthreads);
    //    }
    

}

//////////////////////////
// check file existance //
//////////////////////////
bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

