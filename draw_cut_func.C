

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRotation.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <TLegend.h>
//#include "/u/aali/dirc/prttools/prttools.C"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <cmath>

void draw_cut_func(){
    
    // EXAMPLE CUT PARAMETERS:
    TF1 *fFunc_dEdxCut_SelectHeavy = new TF1("fFunc_dEdxCut_SelectHeavy", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.); // dFunc_dEdxCut_SelectHeavy
    fFunc_dEdxCut_SelectHeavy->SetParameters(4.0, 2.5, 1.25);
    TF1 *fFunc_dEdxCut_SelectLight = new TF1("fFunc_dEdxCut_SelectLight", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.);// dFunc_dEdxCut_SelectLight
    fFunc_dEdxCut_SelectLight->SetParameters(4.0, 2.0, 2.5);
    
    
    TFile *file_histo;
    TCanvas *can_histo;
    TH1F *histo;
    TF1 *funtion;
    
    TString can_name = "Canvas_1";
    TString file_path = "/Users/ahmed/dirc_calib/histo.root";
    file_histo  = new TFile(file_path,"read");
    can_histo=(TCanvas*)file_histo->Get(can_name);
    histo =(TH1F*)can_histo->GetPrimitive("Proton_dEdx_P");
    
    TLegend *legend_histo = new TLegend(0.107769, 0.713904, 0.358396, 0.882353);
    legend_histo->SetHeader("CDC dEdx Proton","C"); // option "C" allows to center the header
    
    std::cout<<"############"<< "no problem 5 " <<std::endl;
    gStyle->SetOptStat(0);
//    histo->SetLineColor(kPink);
//    histo->SetLineStyle(1);
//    histo->GetXaxis()->SetTitle("P_{p}");
//    histo->GetYaxis()->SetTitle("dEdx CDC [KeV/cm]");
//    histo->GetXaxis()->SetTitleSize(0.05);
//    histo->GetYaxis()->SetTitleSize(0.05);
//    histo->GetXaxis()->SetTitleOffset(0.9);
//    histo->GetYaxis()->SetTitleOffset(1.0);
//    histo->SetFillColor(kPink);
//    histo->SetFillStyle(3001);
    std::cout<<"############"<< "no problem 6 " <<std::endl;
    TCanvas * c= new TCanvas("c","c", 400, 500);
    histo->Draw("colz");
    fFunc_dEdxCut_SelectHeavy->Draw("same");
    fFunc_dEdxCut_SelectLight->Draw("same");
    //legend_histo->Draw();
    
}
