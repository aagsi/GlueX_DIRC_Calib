#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRotation.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TLegend.h>
#include "/u/aali/dirc/prttools/prttools.C"
//#include "/Users/ahmed/dirc/prttools/prttools.C"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TList.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#define PI 3.14159265
//pimkpks_20_x2cut2.root
////////////////////
// proto types//////
////////////////////
// fitting functions
Double_t* FitHisto(TH1I *KsMass_KinFit);
// file existance
bool exists_test (const std::string& name);
// sideband method
Double_t fline(Double_t *x, Double_t *par);
Bool_t reject, reject_sd;


//.x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C
//DPROOFLiteManager::Process_Tree("/cache/halld/RunPeriod-2017-01/analysis/ver08/tree_pimkpks__B3_M16/merged/tree_pimkpks__B3_M16_03*", "pimkpks__B3_M16_Tree", "DSelector_justin_1_analyzer.C+", 8, "pimkpks_20_x2cut2_test.root")


void fitks() {
    TFile *ffile_data;
    TH1I *KsMass_KinFit;
    Int_t path_length_cut = 20;
    TString nid = Form("_%2.0d", path_length_cut);
    TString data_path = Form("/u/aali/ks/pimkpks_%d_x2cut2_test.root", path_length_cut);
    ffile_data  = new TFile(data_path, "READ");
    KsMass_KinFit=(TH1I*)ffile_data->Get("KsMass_KinFit");
    KsMass_KinFit->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}]");
    
    TLegend* legend = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    prt_canvasAdd("r_sideband"+nid,800,400);
    gStyle->SetOptStat(0);


    Double_t binwidth = KsMass_KinFit->GetBinWidth(1);
    Double_t ks_mass_fit=0;
    Double_t ks_sigma_fit=0;
    //TSpectrum *fSpect= new TSpectrum(10);
    //Int_t nfound = fSpect->Search(KsMass_KinFit,1,"",0.9);
    //if(nfound>0) ks_mass_fit = fSpect->GetPositionX()[0];
    TF1 *fFit = new TF1("fFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    ks_mass_fit =  KsMass_KinFit->GetXaxis()->GetBinCenter(KsMass_KinFit->GetMaximumBin());
    fFit->SetParameters(500,ks_mass_fit,0.001);
    fFit->SetParNames("p0","KsMass","#sigma","p3","p4");
    //fFit->SetParLimits(0,100,100000);
    //fFit->SetParLimits(1,ks_mass_fit-0.05,ks_mass_fit+0.05);
    //fFit->SetParLimits(2,0.0005,0.8);
    KsMass_KinFit->Fit("fFit","M0","",ks_mass_fit-0.06,ks_mass_fit+0.06);
    //KsMass_KinFit->Fit("fFit","R");
    Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
    ks_mass_fit = fFit->GetParameter(1);
    ks_sigma_fit = fFit->GetParameter(2);
    Double_t KsMass_minus_5_sgma = ks_mass_fit-5*ks_sigma_fit;
    Double_t KsMass_plus_5_sgma = ks_mass_fit+5*ks_sigma_fit;
    Double_t KsMass_minus_3_sgma = ks_mass_fit-3*ks_sigma_fit;
    Double_t KsMass_plus_3_sgma = ks_mass_fit+3*ks_sigma_fit;
    Double_t r_min = ks_mass_fit-8*ks_sigma_fit;
    Double_t r_max = ks_mass_fit+8*ks_sigma_fit;
    Double_t sumundercurve = fFit->Integral(KsMass_minus_3_sgma,KsMass_plus_3_sgma)/binwidth;
    /*
    Double_t fitmin=KsMass_minus_5_sgma, fitmax=KsMass_plus_5_sgma;
    //Double_t fitmin=0.3, fitmax=0.7;
    TF1 *test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
    test_voigt->SetParameters(100,0.5,0.001,0.02,800,-1000,0.);
    test_voigt->SetParNames("signalyield","KsMass","#Gamma","#sigma","BGconstant","BGslope","curve");
    KsMass_KinFit->Fit("test_voigt","M","",fitmin,fitmax);
    //KsMass_KinFit->Fit("TMath::Voigt","M","",ks_mass_fit-0.06,ks_mass_fit+0.06);
    // now extract signal and BG yields
    Double_t binwidth = KsMass_KinFit->GetBinWidth(1);
    Double_t sigbg=test_voigt->Integral(fitmin,fitmax)/binwidth;
    // create pure bkg function
    TF1 *ff2=new TF1("ff2","pol2");
    ff2->SetLineColor(kRed+1);
    ff2->SetLineStyle(2);
    ff2->SetNpx(500);  // some style setting
    ff2->SetRange(fitmin,fitmax);
    // copy parameters from full fcn
    for (int i=0; i<3; ++i) ff2->SetParameter(i,test_voigt->GetParameter(4+i));
    // add it to the histogram (so that it is shown when drawing the histogram)
    KsMass_KinFit->GetListOfFunctions()->Add(ff2);
    // compute bkg intergral
    Double_t bkg = ff2->Integral(fitmin,fitmax)/binwidth;
    // #signals = #total - #bkg
    if (bkg>0) sigbg -= bkg;
    cout << "Signal yield: " << sigbg << "   Background: " << bkg <<endl;
     */
    ////////////////
    // SideBand  ///
    ////////////////
    TF1 *fl = new TF1("fl",fline,0.6,1,6);
    fl->FixParameter (2,KsMass_minus_5_sgma);  // del peak mine
    fl->FixParameter (3,KsMass_plus_5_sgma); // del peak max
    fl->FixParameter (4,r_min); // del smaller than 0.73
    fl->FixParameter (5,r_max); // del  biger than 0.9

    //fl->SetParameters(0,p3);
    //fl->SetParameters(1,p4);
    //fl->SetParNames("plinear0","plinear1");
    //fl->SetParLimits(1,p3-5,p3+5);
    //fl->SetParLimits(1,p4-0,p4+0);
    ////////////////////////////////
    //fit only the linear background excluding the signal area and boarders
    reject = kTRUE;
    KsMass_KinFit->Fit(fl,"0+");
    KsMass_KinFit->Fit(fl,"0+");
    reject = kFALSE;
    //store 2 separate functions for visualization
    TF1 *fleft = new TF1("fleft",fline, r_min , KsMass_minus_5_sgma ,6);
    fleft->SetParameters(fl->GetParameters());
    KsMass_KinFit->GetListOfFunctions()->Add(fleft);
    gROOT->GetListOfFunctions()->Remove(fleft);
    TF1 *fright = new TF1("fright",fline, KsMass_plus_5_sgma ,r_max,6);
    fright->SetParameters(fl->GetParameters());
    KsMass_KinFit->GetListOfFunctions()->Add(fright);
    gROOT->GetListOfFunctions()->Remove(fright);

    TF1 *fcenter = new TF1("fcenter",fline, KsMass_minus_3_sgma ,KsMass_plus_3_sgma,6);
    fcenter->SetParameters(fl->GetParameters());
    KsMass_KinFit->GetListOfFunctions()->Add(fcenter);
    gROOT->GetListOfFunctions()->Remove(fcenter);

    fleft->SetLineColor(kBlue);
    fright->SetLineColor(kBlue);
    fcenter->SetLineColor(kMagenta);
    Double_t B_center = fl->Integral(KsMass_minus_3_sgma,KsMass_plus_3_sgma)/binwidth;
    Double_t B_lift = fl->Integral(r_min,KsMass_minus_5_sgma)/binwidth;
    Double_t B_right = fl->Integral(KsMass_plus_5_sgma,r_max)/binwidth;

    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SIDBAND @@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
    std::cout<<"############  BG center integration = "<< B_center <<std::endl;
    std::cout<<"############  BG lift integration = "<< B_lift <<std::endl;
    std::cout<<"############  BG right integration= "<< B_right <<std::endl;
    std::cout<<"############  BG right + BG lift = "<< B_right + B_lift <<std::endl;
    std::cout<<"############  signal/ (signal + Bkg) (gaus+pol1 estimation)= "<< (sumundercurve-B_center) *100/ (sumundercurve )<<std::endl;
    std::cout<<"############  signal/ Bkg (gaus+pol1 estimation)= "<< (sumundercurve-B_center) *100/ B_center<<std::endl;
    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;




    Double_t fitmin=KsMass_minus_3_sgma, fitmax=KsMass_plus_3_sgma;
    TF1 *test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
    test_voigt->SetParameters(100,0.5,0.001,0.02,800,-1000,0.);
    test_voigt->SetParNames("signalyield","KsMass","#Gamma","#sigma","BGconstant","BGslope","curve");
    test_voigt->SetParameters(53,0.49, 0.003, 0.001, 141, -1266, 1585);
    KsMass_KinFit->Fit("test_voigt","M+","",fitmin,fitmax);
    // now extract signal and BG yields
    Double_t sigbg=test_voigt->Integral(fitmin,fitmax)/binwidth;
    Double_t fraction= (sigbg-B_center) *100/ sigbg ;
    
    std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;


    ////////////////
    // side band ///
    ////////////////
    if (true) {
        
        legend->AddEntry(KsMass_KinFit,"data","f");
        KsMass_KinFit->SetTitle(Form("select path lenght >  %d mm", path_length_cut));
        KsMass_KinFit->Draw();
        prt_canvasGet("r_sideband"+nid)->Update();
        /*
        TLine *line = new TLine(0,0,0,1000);
        line->SetX1(ks_mass_fit);
        line->SetX2(ks_mass_fit);
        line->SetY1(gPad->GetUymin());
        line->SetY2(gPad->GetUymax());
        line->SetLineColor(kRed);
        line->Draw();
        */

        //legend->Draw();
        //prt_waitPrimitive("r_sideband"+nid);// wait

        TLine *line_sigm_r = new TLine(0,0,0,1000);
        line_sigm_r->SetX1(KsMass_plus_3_sgma);
        line_sigm_r->SetX2(KsMass_plus_3_sgma);
        line_sigm_r->SetY1(gPad->GetUymin());
        line_sigm_r->SetY2(gPad->GetUymax());
        line_sigm_r->SetLineColor(46);
        line_sigm_r->Draw();

        TLine *line_sigm_r2 = new TLine(0,0,0,1000);
        line_sigm_r2->SetX1(KsMass_plus_5_sgma);
        line_sigm_r2->SetX2(KsMass_plus_5_sgma);
        line_sigm_r2->SetY1(gPad->GetUymin());
        line_sigm_r2->SetY2(gPad->GetUymax());
        line_sigm_r2->SetLineColor(46);
        line_sigm_r2->Draw();

        // TLine *line_sigm_r3 = new TLine(0,0,0,1000);
        // line_sigm_r3->SetX1(out_array[9]);
        // line_sigm_r3->SetX2(out_array[9]);
        // line_sigm_r3->SetY1(gPad->GetUymin());
        // line_sigm_r3->SetY2(gPad->GetUymax());
        // line_sigm_r3->SetLineColor(kBlue);
        // line_sigm_r3->Draw();

        TLine *line_sigm_l = new TLine(0,0,0,1000);
        line_sigm_l->SetX1(KsMass_minus_5_sgma);
        line_sigm_l->SetX2(KsMass_minus_5_sgma);
        line_sigm_l->SetY1(gPad->GetUymin());
        line_sigm_l->SetY2(gPad->GetUymax());
        line_sigm_l->SetLineColor(46);
        line_sigm_l->Draw();

        TLine *line_sigm_l2 = new TLine(0,0,0,1000);
        line_sigm_l2->SetX1(KsMass_minus_3_sgma);
        line_sigm_l2->SetX2(KsMass_minus_3_sgma);
        line_sigm_l2->SetY1(gPad->GetUymin());
        line_sigm_l2->SetY2(gPad->GetUymax());
        line_sigm_l2->SetLineColor(46);
        line_sigm_l2->Draw();

        // TLine *line_sigm_l3 = new TLine(0,0,0,1000);
        // line_sigm_l3->SetX1(out_array[8]);
        // line_sigm_l3->SetX2(out_array[8]);
        // line_sigm_l3->SetY1(gPad->GetUymin());
        // line_sigm_l3->SetY2(gPad->GetUymax());
        // line_sigm_l3->SetLineColor(kBlue);
        // line_sigm_l3->Draw();

        TBox *cutObj1 = new TBox(KsMass_minus_5_sgma, gPad->GetUymin(), KsMass_minus_3_sgma, gPad->GetUymax());
        cutObj1->SetFillColor(46-4);
        cutObj1->SetFillStyle(3004);
        cutObj1->SetLineWidth(0);
        cutObj1->Draw("lsames");

        TBox *cutObj2 = new TBox(KsMass_plus_3_sgma, gPad->GetUymin(), KsMass_plus_5_sgma, gPad->GetUymax());
        cutObj2->SetFillColor(46-4);
        cutObj2->SetFillStyle(3004);
        cutObj2->SetLineWidth(0);
        cutObj2->Draw("lsames");

        prt_canvasGet("r_sideband"+nid)->Update();
        
   
   
      TPaveText *pt = new TPaveText(0.5384226,1029.375,0.6323293,2089.833);

   //pt->AddLine(0.572238,745.977,0.648483,1788.15);

   
   pt->AddText(Form("#frac{Signal}{Signal+Bkg}=  %2.2f %%", fraction));
   
   
   
   
   
   pt->Draw();
   
   
      prt_canvasGet("r_sideband"+nid)->Update();
        
        
        
        
        
    }
}

//////////////////////////
// check file existance //
//////////////////////////

bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

////////////////////////
// side calculation ////
////////////////////////
Double_t fline(Double_t *x, Double_t *par)
{
    if ((reject && (x[0] > par[2] && x[0] < par[3])) || (x[0]< par[4] || x[0]> par[5])) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}







