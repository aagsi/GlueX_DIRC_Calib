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
//pimkpKs_20_x2cut2.root
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
//DPROOFLiteManager::Process_Tree("/cache/halld/RunPeriod-2017-01/analysis/ver08/tree_pimkpKs__B3_M16/merged/tree_pimkpKs__B3_M16_03*", "pimkpKs__B3_M16_Tree", "DSelector_justin_1_analyzer.C+", 8, "pimkpKs_20_x2cut2_test.root")

bool bool_ks = true;

void fitks() {
    prt_savepath="ks";
    TString data_path = "/u/aali/ks/pimkpks_path_2_chi2.root" ;

    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TFile *ffile_data;
    TH1I *KsMass_KinFit;
    TH1I *Ks_KinFit_ChiSq;
    TH1I *Ks_Path_Length;
    TH2I *PiMinus_p_theta;
    //TString data_path = Form("/u/aali/Ks/pimkpKs_%d_x2cut2_test_PID.root", path_length_cut);

    ffile_data  = new TFile(data_path, "READ");
    KsMass_KinFit=(TH1I*)ffile_data->Get("KsMass_KinFit");
    Ks_KinFit_ChiSq=(TH1I*)ffile_data->Get("KinFitChiSq");
    Ks_Path_Length=(TH1I*)ffile_data->Get("dHist_DetachedPathLength");
    PiMinus_p_theta=(TH2I*)ffile_data->Get("cartizian_theta_mom");
    KsMass_KinFit->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}]");

    Ks_KinFit_ChiSq->GetXaxis()->SetTitle("Kinematic Fit #chi^{2}/NDF");
    Ks_KinFit_ChiSq->GetYaxis()->SetTitle("Counts[#]");
    Ks_KinFit_ChiSq->SetFillColor(kMagenta);
    Ks_KinFit_ChiSq->SetFillStyle(3001);

    Ks_Path_Length->GetXaxis()->SetTitle("Ks Path Length [cm]");
    Ks_Path_Length->GetYaxis()->SetTitle("Counts[#]");
    Ks_Path_Length->SetFillColor(kGreen);
    Ks_Path_Length->SetFillStyle(3001);

    PiMinus_p_theta->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Polar angle X charge [deg]");
    PiMinus_p_theta->GetYaxis()->SetTitle("P [GeV/c]");




    TLegend* legend = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    prt_canvasAdd("r_Ks",800,400);
    gStyle->SetOptStat(0);


    Double_t binwidth = KsMass_KinFit->GetBinWidth(1);
    Double_t Ks_mass_fit=0;
    Double_t Ks_sigma_fit=0;
    //TSpectrum *fSpect= new TSpectrum(10);
    //Int_t nfound = fSpect->Search(KsMass_KinFit,1,"",0.9);
    //if(nfound>0) Ks_mass_fit = fSpect->GetPositionX()[0];
    TF1 *fFit = new TF1("fFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    Ks_mass_fit =  KsMass_KinFit->GetXaxis()->GetBinCenter(KsMass_KinFit->GetMaximumBin());
    fFit->SetParameters(500,Ks_mass_fit,0.001);
    fFit->SetParNames("p0","KsMass","#sigma","p3","p4");
    //fFit->SetParLimits(0,100,100000);
    //fFit->SetParLimits(1,Ks_mass_fit-0.05,Ks_mass_fit+0.05);
    //fFit->SetParLimits(2,0.0005,0.8);
    KsMass_KinFit->Fit("fFit","M0","",Ks_mass_fit-0.06,Ks_mass_fit+0.06);
    //KsMass_KinFit->Fit("fFit","R");
    Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
    Ks_mass_fit = fFit->GetParameter(1);
    Ks_sigma_fit = fFit->GetParameter(2);
    Double_t KsMass_minus_5_sgma = Ks_mass_fit-5*Ks_sigma_fit;
    Double_t KsMass_plus_5_sgma = Ks_mass_fit+5*Ks_sigma_fit;
    Double_t KsMass_minus_3_sgma = Ks_mass_fit-3*Ks_sigma_fit;
    Double_t KsMass_plus_3_sgma = Ks_mass_fit+3*Ks_sigma_fit;
    Double_t r_min = Ks_mass_fit-8*Ks_sigma_fit;
    Double_t r_max = Ks_mass_fit+8*Ks_sigma_fit;
    Double_t sumundercurve = fFit->Integral(KsMass_minus_3_sgma,KsMass_plus_3_sgma)/binwidth;
    /*
    Double_t fitmin=KsMass_minus_5_sgma, fitmax=KsMass_plus_5_sgma;
    //Double_t fitmin=0.3, fitmax=0.7;
    TF1 *test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
    test_voigt->SetParameters(100,0.5,0.001,0.02,800,-1000,0.);
    test_voigt->SetParNames("signalyield","KsMass","#Gamma","#sigma","BGconstant","BGslope","curve");
    KsMass_KinFit->Fit("test_voigt","M","",fitmin,fitmax);
    //KsMass_KinFit->Fit("TMath::Voigt","M","",Ks_mass_fit-0.06,Ks_mass_fit+0.06);
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

    if (false) {
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
    }
    
    Double_t fraction=-1;
    if (bool_ks== true) {
        //Double_t fitmin=KsMass_minus_3_sgma, fitmax=KsMass_plus_3_sgma;
        Double_t fitmin=0.36, fitmax=0.63;
        TF1 *test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
        test_voigt->SetParNames("signalyield","KsMass","#Gamma","#sigma","BGconstant","BGslope","curve");
        test_voigt->SetParameters(53,0.49, 0.003, 0.001, 141, -1266, 1585);
        test_voigt->SetNpx(10000);
        KsMass_KinFit->Fit("test_voigt","M+","",fitmin,fitmax);
        // now extract signal and BG yields
        //Double_t sigbg=test_voigt->Integral(fitmin,fitmax)/binwidth;
        //Double_t fraction= (sigbg-B_center) *100/ sigbg ;
        //std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;

        // now extract signal and BG yields
        Double_t sigbg=test_voigt->Integral(KsMass_minus_3_sgma,KsMass_plus_3_sgma)/binwidth;
        // create pure bkg function
        TF1 *ff2=new TF1("ff2","pol2");
        ff2->SetLineColor(kBlue+1);
        ff2->SetLineStyle(2);
        ff2->SetNpx(1000);  // some style setting
        ff2->SetRange(fitmin,fitmax);
        // copy parameters from full fcn
        for (int i=0; i<3; ++i) ff2->SetParameter(i,test_voigt->GetParameter(4+i));
        // add it to the histogram (so that it is shown when drawing the histogram)
        KsMass_KinFit->GetListOfFunctions()->Add(ff2);
        // compute bkg intergral
        Double_t bkg = ff2->Integral(KsMass_minus_3_sgma,KsMass_plus_3_sgma)/binwidth;
        // #signals = #total - #bkg
        fraction= (sigbg-bkg) *100/ sigbg ;
        std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;
    }

    ////////////////
    // side band ///
    ////////////////
    if (true) {

        legend->AddEntry(KsMass_KinFit,"data","f");
        //KsMass_KinFit->SetTitle(Form("#Lambda Invariant Mass, select path lenght >  %d cm", path_length_cut));
        KsMass_KinFit->SetTitle("Ks Invariant Mass");
        KsMass_KinFit->Draw();
        prt_canvasGet("r_Ks")->Update();
        /*
        TLine *line = new TLine(0,0,0,1000);
        line->SetX1(Ks_mass_fit);
        line->SetX2(Ks_mass_fit);
        line->SetY1(gPad->GetUymin());
        line->SetY2(gPad->GetUymax());
        line->SetLineColor(kRed);
        line->Draw();
        */

        //legend->Draw();
        //prt_waitPrimitive("r_Ks"+nid);// wait

        TLine *line_ChiSq = new TLine(0,0,0,1000);
        line_ChiSq->SetX1(KsMass_plus_3_sgma);
        line_ChiSq->SetX2(KsMass_plus_3_sgma);
        line_ChiSq->SetY1(gPad->GetUymin());
        line_ChiSq->SetY2(gPad->GetUymax());
        line_ChiSq->SetLineColor(46);
        line_ChiSq->Draw();

        //~ TLine *line_ChiSq2 = new TLine(0,0,0,1000);
        //~ line_ChiSq2->SetX1(KsMass_plus_5_sgma);
        //~ line_ChiSq2->SetX2(KsMass_plus_5_sgma);
        //~ line_ChiSq2->SetY1(gPad->GetUymin());
        //~ line_ChiSq2->SetY2(gPad->GetUymax());
        //~ line_ChiSq2->SetLineColor(46);
        //~ line_ChiSq2->Draw();

        // TLine *line_ChiSq3 = new TLine(0,0,0,1000);
        // line_ChiSq3->SetX1(out_array[9]);
        // line_ChiSq3->SetX2(out_array[9]);
        // line_ChiSq3->SetY1(gPad->GetUymin());
        // line_ChiSq3->SetY2(gPad->GetUymax());
        // line_ChiSq3->SetLineColor(kBlue);
        // line_ChiSq3->Draw();

        //~ TLine *line_sigm_l = new TLine(0,0,0,1000);
        //~ line_sigm_l->SetX1(KsMass_minus_5_sgma);
        //~ line_sigm_l->SetX2(KsMass_minus_5_sgma);
        //~ line_sigm_l->SetY1(gPad->GetUymin());
        //~ line_sigm_l->SetY2(gPad->GetUymax());
        //~ line_sigm_l->SetLineColor(46);
        //~ line_sigm_l->Draw();

        TLine *line_path_length = new TLine(0,0,0,1000);
        line_path_length->SetX1(KsMass_minus_3_sgma);
        line_path_length->SetX2(KsMass_minus_3_sgma);
        line_path_length->SetY1(gPad->GetUymin());
        line_path_length->SetY2(gPad->GetUymax());
        line_path_length->SetLineColor(46);
        line_path_length->Draw();

        // TLine *line_sigm_l3 = new TLine(0,0,0,1000);
        // line_sigm_l3->SetX1(out_array[8]);
        // line_sigm_l3->SetX2(out_array[8]);
        // line_sigm_l3->SetY1(gPad->GetUymin());
        // line_sigm_l3->SetY2(gPad->GetUymax());
        // line_sigm_l3->SetLineColor(kBlue);
        // line_sigm_l3->Draw();

        //~ TBox *cutObj1 = new TBox(KsMass_minus_5_sgma, gPad->GetUymin(), KsMass_minus_3_sgma, gPad->GetUymax());
        //~ cutObj1->SetFillColor(46-4);
        //~ cutObj1->SetFillStyle(3004);
        //~ cutObj1->SetLineWidth(0);
        //~ cutObj1->Draw("lsames");
        //~
        //~ TBox *cutObj2 = new TBox(KsMass_plus_3_sgma, gPad->GetUymin(), KsMass_plus_5_sgma, gPad->GetUymax());
        //~ cutObj2->SetFillColor(46-4);
        //~ cutObj2->SetFillStyle(3004);
        //~ cutObj2->SetLineWidth(0);
        //~ cutObj2->Draw("lsames");

        // prt_canvasGet("r_Ks"+nid)->Update();
        Int_t Entries = KsMass_KinFit->GetEntries();
        TPaveText *pt = new TPaveText(0.5384226, 1550.24,0.6323293,3600.01);
        //pt->AddLine(0.572238,745.977,0.648483,1788.15);
        pt->AddText(Form("#frac{Signal}{Signal+Bkg}=  %2.2f %%", fraction));
        pt->AddText(Form("Entries =  %d", Entries));
        pt->Draw();


        prt_canvasGet("r_Ks")->Update();

    }
    prt_canvasAdd("r_PiMinus_p_theta",800,400);
    PiMinus_p_theta->Draw("colz");
    prt_canvasGet("r_PiMinus_p_theta")->Update();
    Int_t Entries2 = PiMinus_p_theta->GetEntries();
    TPaveText *pt2 = new TPaveText(1, 9,11.0,11);
    pt2->AddText(Form("number tracks hitting the TOF=  %d", Entries2));
    pt2->Draw();
    prt_canvasGet("r_PiMinus_p_theta")->Update();


    prt_canvasAdd("r_Ks_KinFit_ChiSq",800,400);
    Ks_KinFit_ChiSq->Draw();
    prt_canvasGet("r_Ks_KinFit_ChiSq")->Update();
    TLine *line_ChiSq = new TLine(0,0,0,1000);
    line_ChiSq->SetX1(2.0);
    line_ChiSq->SetX2(2.0);
    line_ChiSq->SetY1(gPad->GetUymin());
    line_ChiSq->SetY2(gPad->GetUymax());
    line_ChiSq->SetLineColor(46);
    line_ChiSq->Draw();
    prt_canvasGet("r_Ks_KinFit_ChiSq")->Update();




    prt_canvasAdd("r_Ks_Path_Length",800,400);
    Ks_Path_Length->Draw();
    prt_canvasGet("r_Ks_Path_Length")->Update();
    TLine *line_path_length = new TLine(0,0,0,1000);
    line_path_length->SetX1(2.0);
    line_path_length->SetX2(2.0);
    line_path_length->SetY1(gPad->GetUymin());
    line_path_length->SetY2(gPad->GetUymax());
    line_path_length->SetLineColor(46);
    line_path_length->Draw();
    prt_canvasGet("r_Ks_Path_Length")->Update();


    prt_canvasSave(2,0);
    prt_canvasDel("*");
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

















