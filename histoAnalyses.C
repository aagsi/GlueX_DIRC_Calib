
#include <TLegend.h>
#include "glxtools.C"
#include "TNtuple.h"
#include "TStopwatch.h"
//#include "TMultiGraph.h"
//#include "TGraph.h"
//#include "TGraphAsymmErrors.h"
//#include "THStack.h"
//#include "TH3D.h"

TStopwatch timer;

void histoAnalyses(){
    timer.Start();
    gStyle->SetPalette(55);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    
    
    TString path ="/Users/ahmed/GlueX_DIRC_Calib/histo_shift.root";
    cout<<"path= " <<path<<endl;
    TFile *f = new TFile(path, "READ");
    
    TH1F*  histo_xy_cR_corr[40][24];
    TH1F*  histo_xy_cR_Notcorr[40][24];
    TH1F*  histo_xy_cD_corr[40][24];
    TH1F*  histo_xy_cD_Notcorr[40][24];
    TH1F*  histo_xy_c_corr[40][24];
    TH1F*  histo_xy_c_Notcorr[40][24];
    
    TH2F * histo_xy_cR_corr_map        = new TH2F( " histo_xy_cR_corr_map"        , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_cR_Notcorr_map     = new TH2F( " histo_xy_cR_Notcorr_map"     , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_cD_corr_map        = new TH2F( " histo_xy_cD_corr_map"        , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_cD_Notcorr_map     = new TH2F( " histo_xy_cD_Notcorr_map"     , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_corr_map         = new TH2F( " histo_xy_c_corr_map"         , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_Notcorr_map      = new TH2F( " histo_xy_c_Notcorr_map"      , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_diff_corr_map    = new TH2F( " histo_xy_c_diff_corr_map"    , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_diff_Notcorr_map = new TH2F( " histo_xy_c_diff_Notcorr_map" , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    TH2F * histo_xy_sprR_corr_map        = new TH2F( " histo_xy_sprR_corr_map"        , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_sprR_Notcorr_map     = new TH2F( " histo_xy_sprR_Notcorr_map"     , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_sprD_corr_map        = new TH2F( " histo_xy_sprD_corr_map"        , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_sprD_Notcorr_map     = new TH2F( " histo_xy_sprD_Notcorr_map"     , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_spr_corr_map         = new TH2F( " histo_xy_spr_corr_map"         , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_spr_Notcorr_map      = new TH2F( " histo_xy_spr_Notcorr_map"      , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_spr_diff_corr_map    = new TH2F( " histo_xy_spr_diff_corr_map"    , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_spr_diff_Notcorr_map = new TH2F( " histo_xy_spr_diff_Notcorr_map" , "; Bar Hit X bin; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    TF1 *fit_cherenkov = new TF1("fit_cherenkov","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",-0.05,0.05);
    fit_cherenkov->SetNpx(1000);
    
    fit_cherenkov->SetLineColor(kBlack);
    fit_cherenkov->SetParameters(100,9,2);
    fit_cherenkov->SetParNames("p0","mean ","spr");
    fit_cherenkov->SetParLimits(0,0.1,1E6);
    fit_cherenkov->SetParLimits(1,-0.05,0.05);
    fit_cherenkov->SetParLimits(2,0.005,0.014);
    

    for(Int_t i=0; i<40; i++){
        for(Int_t j=0; j<24; j++){
            
            histo_xy_cR_corr[i][j] = (TH1F*)f->Get(Form("histo_xy_cR_corr_%d_%d",i,j));
            histo_xy_cR_Notcorr[i][j] = (TH1F*)f->Get(Form("histo_xy_cR_Notcorr_%d_%d",i,j));
            histo_xy_cD_corr[i][j] = (TH1F*)f->Get(Form("histo_xy_cD_corr_%d_%d",i,j));
            histo_xy_cD_Notcorr[i][j] = (TH1F*)f->Get(Form("histo_xy_cD_Notcorr_%d_%d",i,j));
            histo_xy_c_corr[i][j] = (TH1F*)f->Get(Form("histo_xy_c_corr_%d_%d",i,j));
            histo_xy_c_Notcorr[i][j] = (TH1F*)f->Get(Form("histo_xy_c_Notcorr_%d_%d",i,j));
            
            histo_xy_cR_corr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaR_corr=fit_cherenkov->GetParameter(1)*1000;
            double sprR_corr=fit_cherenkov->GetParameter(2)*1000;
            histo_xy_cR_Notcorr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaR_Notcorr=fit_cherenkov->GetParameter(1)*1000;
            double sprR_Notcorr=fit_cherenkov->GetParameter(2)*1000;
            histo_xy_cD_corr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaD_corr=fit_cherenkov->GetParameter(1)*1000;
            double sprD_corr=fit_cherenkov->GetParameter(2)*1000;
            histo_xy_cD_Notcorr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaD_Notcorr=fit_cherenkov->GetParameter(1)*1000;
            double sprD_Notcorr=fit_cherenkov->GetParameter(2)*1000;
            histo_xy_c_corr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaS_corr=fit_cherenkov->GetParameter(1)*1000;
            double sprS_corr=fit_cherenkov->GetParameter(2)*1000;
            histo_xy_c_Notcorr[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
            double thetaS_Notcorr=fit_cherenkov->GetParameter(1)*1000;
            double sprS_Notcorr=fit_cherenkov->GetParameter(2)*1000;
            
            histo_xy_cR_corr_map->Fill(i,j,thetaR_corr);
            histo_xy_cR_Notcorr_map->Fill(i,j,thetaR_Notcorr);
            histo_xy_cD_corr_map->Fill(i,j,thetaD_corr);
            histo_xy_cD_Notcorr_map->Fill(i,j,thetaD_Notcorr);
            histo_xy_c_corr_map->Fill(i,j,thetaS_corr);
            histo_xy_c_Notcorr_map->Fill(i,j,thetaS_Notcorr);
            histo_xy_c_diff_corr_map->Fill(i,j,thetaR_corr-thetaD_corr);
            histo_xy_c_diff_Notcorr_map->Fill(i,j,thetaR_Notcorr-thetaD_Notcorr);
            //////////////////////
            histo_xy_sprR_corr_map->Fill(i,j,sprR_corr);
            histo_xy_sprR_Notcorr_map->Fill(i,j,sprR_Notcorr);
            histo_xy_sprD_corr_map->Fill(i,j,sprD_corr);
            histo_xy_sprD_Notcorr_map->Fill(i,j,sprD_Notcorr);
            histo_xy_spr_corr_map->Fill(i,j,sprS_corr);
            histo_xy_spr_Notcorr_map->Fill(i,j,sprS_Notcorr);
            histo_xy_spr_diff_corr_map->Fill(i,j,sprR_corr-sprD_corr);
            histo_xy_spr_diff_Notcorr_map->Fill(i,j,sprR_Notcorr-sprD_Notcorr);
        }
    }
    
    TString ok="histoAnalyses";
    TFile file(ok,"recreate");
    
    histo_xy_cR_corr_map->Write();
    histo_xy_cR_Notcorr_map->Write();
    histo_xy_cD_corr_map->Write();
    histo_xy_cD_Notcorr_map->Write();
    histo_xy_c_corr_map->Write();
    histo_xy_c_Notcorr_map->Write();
    histo_xy_c_diff_corr_map->Write();
    histo_xy_c_diff_Notcorr_map->Write();
    histo_xy_sprR_corr_map->Write();
    histo_xy_sprR_Notcorr_map->Write();
    histo_xy_sprD_corr_map->Write();
    histo_xy_sprD_Notcorr_map->Write();
    histo_xy_spr_corr_map->Write();
    histo_xy_spr_Notcorr_map->Write();
    histo_xy_spr_diff_corr_map->Write();
    histo_xy_spr_diff_Notcorr_map->Write();
    
    file.Write();
    file.Close();
    
    
//    TCanvas *cc8 = new TCanvas("cc8","cc8",800,500);
//    histo_xy_tdiff_sigmaR->SetMinimum(0);
//    histo_xy_tdiff_sigmaR->SetMaximum(3);
//    histo_xy_tdiff_sigmaR->Draw("COLZ");
//    cc8->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/histo1/histo_xy_tdiff_sigmaR.png");
//    cc8->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/histo1/histo_xy_tdiff_sigmaR.root");
    
    timer.Stop();
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    
    
}
