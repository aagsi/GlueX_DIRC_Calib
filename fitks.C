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

bool bool_ks = false;
bool bool_phi = false;
bool bool_lambda = true;

void fitks() {

    TString data_path,mass_name;
    if (bool_ks) {
        data_path = "/u/aali/ks/ks_final.root" ;
        mass_name= "Ks";
    }
    if(bool_lambda) {
        data_path = "/u/aali/ks/lampda_final.root" ;
        mass_name= "Lampda";
    }

    if(bool_phi) {
        data_path = "/u/aali/ks/phi_final.root" ;
        mass_name= "Phi";
    }

    prt_savepath=mass_name;
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TF1 *fFunc_dEdxCut_SelectLight = new TF1("fFunc_dEdxCut_SelectLight", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.);// dFunc_dEdxCut_SelectLight
    fFunc_dEdxCut_SelectLight->SetParameters(4.0, 2.0, 2.5);



    TFile *ffile_data;
    TH1I *particleMass_KinFit;
    TH1I *particle_KinFit_ChiSq;
    TH1I *particle_Path_Length;

    TH1I *RF_selection, *RF_selection_cut;
    TH1I *MissingMassSquared;

    TH2I *final_stat_p_theta;
    TH2I *dEdx_CDC_Proton;

    ffile_data  = new TFile(data_path, "READ");
    particleMass_KinFit=(TH1I*)ffile_data->Get(mass_name+"Mass_KinFit");
    particle_KinFit_ChiSq=(TH1I*)ffile_data->Get("KinFitChiSq");
    particle_Path_Length=(TH1I*)ffile_data->Get("dHist_DetachedPathLength");
    final_stat_p_theta=(TH2I*)ffile_data->Get("cartizian_theta_mom");
    RF_selection=(TH1I*)ffile_data->Get("dHist_RF");
    RF_selection_cut=(TH1I*)ffile_data->Get("dHist_RF_cut");
    MissingMassSquared=(TH1I*)ffile_data->Get("MissingMassSquared");

    // Histo style
    RF_selection->GetXaxis()->SetTitle("#Deltat_{Beam#gamma - RF}");
    RF_selection->GetYaxis()->SetTitle("Counts[#]");
    RF_selection->SetTitle("#Deltat_{Beam#gamma - RF} Cut");
    RF_selection->SetFillColor(kMagenta);
    RF_selection->SetFillStyle(3001);

    RF_selection_cut->SetFillColor(kGreen);
    RF_selection_cut->SetFillStyle(3001);

    MissingMassSquared->GetXaxis()->SetTitle("Missing Mass Squared (GeV/c^{2})^{2}");
    MissingMassSquared->GetYaxis()->SetTitle("Counts[#]");
    MissingMassSquared->SetFillColor(kBlue);
    MissingMassSquared->SetFillStyle(3001);
    MissingMassSquared->SetTitle("Missing Mass Squared Cut");
    //particleMass_KinFit->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}]");

    particle_KinFit_ChiSq->GetXaxis()->SetTitle("Kinematic Fit #chi^{2}/NDF");
    particle_KinFit_ChiSq->GetYaxis()->SetTitle("Counts[#]");
    particle_KinFit_ChiSq->SetFillColor(kMagenta);
    particle_KinFit_ChiSq->SetFillStyle(3001);


    particle_Path_Length->GetYaxis()->SetTitle("Counts[#]");
    particle_Path_Length->SetFillColor(kGreen);
    particle_Path_Length->SetFillStyle(3001);

    //final_stat_p_theta->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Polar angle X charge [deg]");
    final_stat_p_theta->GetYaxis()->SetTitle("P [GeV/c]");

    prt_canvasAdd("r_particle",800,400);
    gStyle->SetOptStat(0);

    Double_t binwidth = particleMass_KinFit->GetBinWidth(1);
    Double_t particle_mass_fit=0;
    Double_t particle_sigma_fit=0;

    TF1 *fFit = new TF1("fFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    particle_mass_fit =  particleMass_KinFit->GetXaxis()->GetBinCenter(particleMass_KinFit->GetMaximumBin());
    fFit->SetParameters(500,particle_mass_fit,0.001);
    fFit->SetParNames("p0","particleMass","#sigma","p3","p4");

    particleMass_KinFit->Fit("fFit","M0","",particle_mass_fit-0.06,particle_mass_fit+0.059);
    //particleMass_KinFit->Fit("fFit","R");
    Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
    particle_mass_fit = fFit->GetParameter(1);
    particle_sigma_fit = fFit->GetParameter(2);
    Double_t particleMass_minus_5_sgma = particle_mass_fit-5*particle_sigma_fit;
    Double_t particleMass_plus_5_sgma = particle_mass_fit+5*particle_sigma_fit;
    Double_t particleMass_minus_3_sgma = particle_mass_fit-3*particle_sigma_fit;
    Double_t particleMass_plus_3_sgma = particle_mass_fit+3*particle_sigma_fit;
    Double_t r_min = particle_mass_fit-8*particle_sigma_fit;
    Double_t r_max = particle_mass_fit+8*particle_sigma_fit;
    Double_t sumundercurve = fFit->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;

    /*
    Double_t fitmin=particleMass_minus_5_sgma, fitmax=particleMass_plus_5_sgma;
    //Double_t fitmin=0.3, fitmax=0.7;
    TF1 *test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
    test_voigt->SetParameters(100,0.5,0.001,0.02,800,-1000,0.);
    test_voigt->SetParNames("signalyield","particleMass","#Gamma","#sigma","BGconstant","BGslope","curve");
    particleMass_KinFit->Fit("test_voigt","M","",fitmin,fitmax);
    //particleMass_KinFit->Fit("TMath::Voigt","M","",particle_mass_fit-0.06,particle_mass_fit+0.06);
    // now extract signal and BG yields
    Double_t binwidth = particleMass_KinFit->GetBinWidth(1);
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
    particleMass_KinFit->GetListOfFunctions()->Add(ff2);
    // compute bkg intergral
    Double_t bkg = ff2->Integral(fitmin,fitmax)/binwidth;
    // #signals = #total - #bkg
    if (bkg>0) sigbg -= bkg;
    cout << "Signal yield: " << sigbg << "   Background: " << bkg <<endl;
     */

    TLegend* legend = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    legend->AddEntry(particleMass_KinFit,"data","f");
    //particleMass_KinFit->SetTitle(Form("#Lambda Invariant Mass, select path lenght >  %d cm", path_length_cut));
    particleMass_KinFit->Draw();
    prt_canvasGet("r_particle")->Update();
    TLine *line_mass_plus_3_sgma = new TLine(0,0,0,1000);
    line_mass_plus_3_sgma->SetX1(particleMass_plus_3_sgma);
    line_mass_plus_3_sgma->SetX2(particleMass_plus_3_sgma);
    line_mass_plus_3_sgma->SetY1(gPad->GetUymin());
    line_mass_plus_3_sgma->SetY2(gPad->GetUymax());
    line_mass_plus_3_sgma->SetLineColor(46);
    line_mass_plus_3_sgma->Draw();
    prt_canvasGet("r_particle")->Update();

    TLine *line_mass_minus_3_sgma = new TLine(0,0,0,1000);
    line_mass_minus_3_sgma->SetX1(particleMass_minus_3_sgma);
    line_mass_minus_3_sgma->SetX2(particleMass_minus_3_sgma);
    line_mass_minus_3_sgma->SetY1(gPad->GetUymin());
    line_mass_minus_3_sgma->SetY2(gPad->GetUymax());
    line_mass_minus_3_sgma->SetLineColor(46);
    line_mass_minus_3_sgma->Draw();
    prt_canvasGet("r_particle")->Update();

    TF1 *test_voigt;
    TF1 *ff2;
    Double_t sigbg;
    Double_t bkg;
    Int_t Entries ;
    TPaveText *pt, *pt2, *pt3;


    Double_t fitmin=0, fitmax=0;
    Double_t fraction=-1;
    if (bool_ks== true) {
        //Double_t fitmin=particleMass_minus_3_sgma, fitmax=particleMass_plus_3_sgma;
        fitmin=0.36, fitmax=0.63;
        test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
        test_voigt->SetParNames("signalyield","particleMass","#Gamma","#sigma","BGconstant","BGslope","curve");
        test_voigt->SetParameters(53,0.49, 0.003, 0.001, 141, -1266, 1585);
        test_voigt->SetNpx(10000);
        particleMass_KinFit->Fit("test_voigt","M+","",fitmin,fitmax);
        // now extract signal and BG yields
        //Double_t sigbg=test_voigt->Integral(fitmin,fitmax)/binwidth;
        //Double_t fraction= (sigbg-B_center) *100/ sigbg ;
        //std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;

        // now extract signal and BG yields
        sigbg=test_voigt->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;
        // create pure bkg function
        ff2=new TF1("ff2","pol2");
        ff2->SetLineColor(kBlue+1);
        ff2->SetLineStyle(2);
        ff2->SetNpx(1000);  // some style setting
        ff2->SetRange(fitmin,fitmax);
        // copy parameters from full fcn
        for (int i=0; i<3; ++i) ff2->SetParameter(i,test_voigt->GetParameter(4+i));
        // add it to the histogram (so that it is shown when drawing the histogram)
        particleMass_KinFit->GetListOfFunctions()->Add(ff2);
        // compute bkg intergral
        bkg = ff2->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;
        // #signals = #total - #bkg
        fraction= (sigbg-bkg) *100/ sigbg ;
        std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;
        particleMass_KinFit->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}]");
        final_stat_p_theta->GetXaxis()->SetTitle("#pi^{#plus}#pi^{#minus} Polar angle X charge [deg]");
        particleMass_KinFit->SetTitle("Ks Invariant Mass");

        Entries = particleMass_KinFit->GetEntries();
        pt = new TPaveText(0.5384226, 1550.24,0.6323293,3600.01);
        pt->AddText(Form("#frac{Signal}{Signal+Bkg}=  %2.2f %%", fraction));
        //pt->AddText(Form("Entries =  %d", Entries));
        pt->Draw();
        prt_canvasGet("r_particle")->Update();

        dEdx_CDC_Proton=(TH2I*)ffile_data->Get("Hist_ParticleID_pid_precut/Step0__Photon_Proton_->_Pi-_K+_KShort_Proton/Proton/dEdxVsP_CDC");
        dEdx_CDC_Proton->SetTitle("CDC dE/dx");
        particle_Path_Length->GetXaxis()->SetTitle("Ks Path Length [cm]");
    }

    if (bool_lambda== true) {
        //Double_t fitmin=LampdaMass_minus_3_sgma, fitmax=LampdaMass_plus_3_sgma;
        fitmin=1.08, fitmax=1.18;
        test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+pol2(4)",fitmin,fitmax);
        test_voigt->SetParNames("signalyield","LampdaMass","#Gamma","#sigma","BGconstant","BGslope","curve");
        test_voigt->SetParameters(3.74132e+01, 1.11587e+00, 8.64176e-06, 4.19434e-03, 1.41593e+03, -2.72943e+01, -1.18555e+03 );
        test_voigt->SetNpx(500);
        particleMass_KinFit->Fit("test_voigt","M+","",fitmin,fitmax);
        // now extract signal and BG yields
        sigbg=test_voigt->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;
        // create pure bkg function
        ff2=new TF1("ff2","pol2");
        ff2->SetLineColor(kBlue+1);
        ff2->SetLineStyle(2);
        ff2->SetNpx(1000);  // some style setting
        ff2->SetRange(fitmin,fitmax);
        // copy parameters from full fcn
        for (int i=0; i<3; ++i) ff2->SetParameter(i,test_voigt->GetParameter(4+i));
        // add it to the histogram (so that it is shown when drawing the histogram)
        particleMass_KinFit->GetListOfFunctions()->Add(ff2);
        // compute bkg intergral
        bkg = ff2->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;
        // #signals = #total - #bkg
        fraction= (sigbg-bkg) *100/ sigbg ;
        std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;
        particleMass_KinFit->GetXaxis()->SetTitle("P #pi^{#minus} Invariant Mass [GeV/c^{2}]");
        final_stat_p_theta->GetXaxis()->SetTitle("P #pi^{#minus} Polar angle X charge [deg]");
        particleMass_KinFit->SetTitle("#Lambda Invariant Mass");
        //prt_canvasGet("r_particle")->Update();
        Entries = particleMass_KinFit->GetEntries();
        pt = new TPaveText(1.13048,1918.68,1.17568,3595.92);
        pt->AddText(Form("#frac{Signal}{Signal+Bkg}=  %2.2f %%", fraction));
        //pt->AddText(Form("Entries =  %d", Entries));
        pt->Draw();
        prt_canvasGet("r_particle")->Update();
        particle_Path_Length->GetXaxis()->SetTitle("#Lambda Path Length [cm]");
    }


    if (bool_phi == true) {
        fitmin=particleMass_minus_3_sgma;
        fitmax=1.1;
        test_voigt = new TF1("test_voigt","[0]*TMath::Voigt(x-[1],[2],[3],4.)+cheb3(4)",fitmin,fitmax); //cheb3
        test_voigt->SetParameters( 2.37102e+02, 1.01973e+00, 1.33778e-03, 5.93537e-03);
        test_voigt->SetNpx(1000);
        particleMass_KinFit->Fit("test_voigt","M+","",fitmin,fitmax);
        // now extract signal and BG yields
        sigbg=test_voigt->Integral(particleMass_minus_3_sgma,particleMass_plus_3_sgma)/binwidth;
        // create pure bkg function
        ff2=new TF1("ff2","cheb3");
        ff2->SetLineColor(kBlue+1);
        ff2->SetLineStyle(2);
        ff2->SetNpx(1000);  // some style setting
        ff2->SetRange(fitmin,fitmax);
        // copy parameters from full fcn
        //for (int i=0; i<4; ++i) ff2->SetParameter(i,test_voigt->GetParameter(4+i));
        ff2->SetParameters(-1.06302e+05 , 3.86587e+03, 1.72174e+05, -6.98597e+04);
        // add it to the histogram (so that it is shown when drawing the histogram)
        particleMass_KinFit->GetListOfFunctions()->Add(ff2);
        // compute bkg intergral
        bkg = ff2->Integral(particleMass_minus_3_sgma, particleMass_plus_3_sgma)/binwidth;
        // #signals = #total - #bkg
        fraction= (sigbg-bkg) *100/ sigbg ;
        std::cout<<"############  signal/ (signal + Bkg)= "<< fraction <<std::endl;
        std::cout<<"############  bkg= "<< bkg <<std::endl;
        particleMass_KinFit->GetXaxis()->SetTitle("K^{#plus} K^{#minus} Invariant Mass [GeV/c^{2}]");
        final_stat_p_theta->GetXaxis()->SetTitle("K^{#plus} K^{#minus} Polar angle X charge [deg]");
        particleMass_KinFit->SetTitle("#phi Invariant Mass");

        prt_canvasGet("r_particle")->Update();
        Entries = particleMass_KinFit->GetEntries();
        pt = new TPaveText(0.908835,28298.4,0.996867,59771.9);
        pt->AddText(Form("#frac{Signal}{Signal+Bkg}=  %2.2f %%", fraction));
        //pt->AddText(Form("Entries =  %d", Entries));
        pt->Draw();
        prt_canvasGet("r_particle")->Update();
        dEdx_CDC_Proton=(TH2I*)ffile_data->Get("Hist_ParticleID_pid_precut/Step0__Photon_Proton_->_K+_K-_Proton/Proton/dEdxVsP_CDC");
        dEdx_CDC_Proton->SetTitle("CDC dE/dx");

    }

    if (true) {
        Double_t pt1_x1, pt1_x2, pt1_x3, pt1_x4;
        Double_t pt2_x1, pt2_x2, pt2_x3, pt2_x4;
        Double_t pt3_x1, pt3_x2, pt3_x3, pt3_x4;
        pt1_x1=-2.5;
        pt1_x2=9.14439;
        pt1_x3=11.8045;
        pt1_x4=10.3476;
        if (bool_ks) {
            pt2_x1=12.6175;
            pt2_x2= 14479.2;
            pt2_x3= 24.9138;
            pt2_x4=15971.9;

            pt3_x1=6.5;
            pt3_x2=1.12264e+06;
            pt3_x3=14.9;
            pt3_x4=1.25172e+06;
        }

        if (bool_phi) {
            pt2_x1=12.6175;
            pt2_x2= 409417;
            pt2_x3= 24.8747;
            pt2_x4=445070;
        }
        
        if (bool_lambda) {
            pt2_x1=12.6175;
            pt2_x2= 3361.76;
            pt2_x3= 24.9138;
            pt2_x4= 3649.79;

            pt3_x1=6.5;
            pt3_x2=142756;
            pt3_x3=14.9;
            pt3_x4=154089;
        }

        prt_canvasAdd("r_final_stat_p_theta",800,400);
        final_stat_p_theta->SetTitle("Final state Momentum Vs Polar angle distribution");
        final_stat_p_theta->Draw("colz");
        prt_canvasGet("r_final_stat_p_theta")->Update();
        Int_t Entries_F = (final_stat_p_theta->GetEntries()/2)*100/Entries;
        Int_t Entries_hist =   final_stat_p_theta->GetEntries()/2.0;
        TPaveText *pt1 = new TPaveText(pt1_x1, pt1_x2, pt1_x3, pt1_x4);
        pt1->AddText(Form("Number of reconstructed final state tracks/2 in TOF ~ %d %% ", Entries_F));
        //pt1->AddText(Form("Entries =  %d", Entries_hist));
        pt1->Draw();
        prt_canvasGet("r_final_stat_p_theta")->Update();

        prt_canvasAdd("r_particle_KinFit_ChiSq",800,400);
        particle_KinFit_ChiSq->SetTitle("Kinematic Fit #chi^{2}/NDF");
        particle_KinFit_ChiSq->Draw();
        prt_canvasGet("r_particle_KinFit_ChiSq")->Update();
        Double_t chi_line= 2.0 ;
        if (bool_lambda)chi_line =5.0;
        TLine *line_ChiSq = new TLine(0,0,0,1000);
        line_ChiSq->SetX1(chi_line);
        line_ChiSq->SetX2(chi_line);
        line_ChiSq->SetY1(gPad->GetUymin());
        line_ChiSq->SetY2(gPad->GetUymax());
        line_ChiSq->SetLineColor(46);
        line_ChiSq->Draw();
        prt_canvasGet("r_particle_KinFit_ChiSq")->Update();




        pt2 = new TPaveText(pt2_x1, pt2_x2, pt2_x3, pt2_x4);
        pt2->AddText(Form("select ChiSq/NDF <  %2.2f", chi_line));
        pt2->Draw();
        prt_canvasGet("r_particle_KinFit_ChiSq")->Update();

        if (bool_lambda||bool_ks) {
            prt_canvasAdd("r_particle_Path_Length",800,400);
            particle_Path_Length->SetTitle("Path Length");
            particle_Path_Length->Draw();
            prt_canvasGet("r_particle_Path_Length")->Update();
            TLine *line_path_length = new TLine(0,0,0,1000);
            Double_t length_line= 2.0 ;
            if (bool_lambda)length_line =0.5;

            line_path_length->SetX1(length_line);
            line_path_length->SetX2(length_line);
            line_path_length->SetY1(gPad->GetUymin());
            line_path_length->SetY2(gPad->GetUymax());
            line_path_length->SetLineColor(46);
            line_path_length->Draw();
            prt_canvasGet("r_particle_Path_Length")->Update();
            pt3 = new TPaveText(pt3_x1, pt3_x2, pt3_x3, pt3_x4);
            pt3->AddText(Form("select path length >  %2.2f [cm]", length_line));
            pt3->Draw();
            prt_canvasGet("r_particle_Path_Length")->Update();



        }
        if (bool_phi||bool_ks) {
            prt_canvasAdd("r_dEdx_CDC_Proton",800,400);
            dEdx_CDC_Proton->Draw("colz");
            fFunc_dEdxCut_SelectLight->Draw("same");
        }

        prt_canvasAdd("r_RF",800,400);
        RF_selection->Draw();
        RF_selection_cut->Draw("same");
        prt_canvasGet("r_RF")->Update();

        if (bool_phi||bool_ks) {
            TLine *line_RF_lift = new TLine(0,0,0,1000);
            line_RF_lift->SetX1(-2.0);
            line_RF_lift->SetX2(-2.0);
            line_RF_lift->SetY1(gPad->GetUymin());
            line_RF_lift->SetY2(gPad->GetUymax());
            line_RF_lift->SetLineColor(46);
            line_RF_lift->Draw();

            TLine *line_RF_right = new TLine(0,0,0,1000);
            line_RF_right->SetX1(2.0);
            line_RF_right->SetX2(2.0);
            line_RF_right->SetY1(gPad->GetUymin());
            line_RF_right->SetY2(gPad->GetUymax());
            line_RF_right->SetLineColor(46);
            line_RF_right->Draw();
            prt_canvasGet("r_RF")->Update();
        } else {
            TLine *line_RF_lambda_lift = new TLine(0,0,0,1000);
            line_RF_lambda_lift->SetX1(-6.0);
            line_RF_lambda_lift->SetX2(-6.0);
            line_RF_lambda_lift->SetY1(gPad->GetUymin());
            line_RF_lambda_lift->SetY2(gPad->GetUymax());
            line_RF_lambda_lift->SetLineColor(46);
            line_RF_lambda_lift->Draw();

            TLine *line_RF_lambda_right = new TLine(0,0,0,1000);
            line_RF_lambda_right->SetX1(6.0);
            line_RF_lambda_right->SetX2(6.0);
            line_RF_lambda_right->SetY1(gPad->GetUymin());
            line_RF_lambda_right->SetY2(gPad->GetUymax());
            line_RF_lambda_right->SetLineColor(46);
            line_RF_lambda_right->Draw();
            prt_canvasGet("r_RF")->Update();
        }


        prt_canvasAdd("r_MissingMassSquared",800,400);
        MissingMassSquared->Draw();
        prt_canvasGet("r_MissingMassSquared")->Update();
        TLine *line_MissingMassSquared_lift = new TLine(0,0,0,1000);
        line_MissingMassSquared_lift->SetX1(-0.01);
        line_MissingMassSquared_lift->SetX2(-0.01);
        line_MissingMassSquared_lift->SetY1(gPad->GetUymin());
        line_MissingMassSquared_lift->SetY2(gPad->GetUymax());
        line_MissingMassSquared_lift->SetLineColor(46);
        line_MissingMassSquared_lift->Draw();

        TLine *line_MissingMassSquared_right = new TLine(0,0,0,1000);
        line_MissingMassSquared_right->SetX1(0.01);
        line_MissingMassSquared_right->SetX2(0.01);
        line_MissingMassSquared_right->SetY1(gPad->GetUymin());
        line_MissingMassSquared_right->SetY2(gPad->GetUymax());
        line_MissingMassSquared_right->SetLineColor(46);
        line_MissingMassSquared_right->Draw();
        prt_canvasGet("r_MissingMassSquared")->Update();

    }

    prt_canvasSave(2,0);
    prt_canvasDel("*");
}










