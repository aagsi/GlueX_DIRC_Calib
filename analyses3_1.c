//#include "TMultiGraph.h"
//#include "TGraph.h"
//#include "TGraphAsymmErrors.h"
#include <TLegend.h>
#include "glxtools.C"
//#include "THStack.h"
//#include "TH3D.h"
#include "TNtuple.h"

//root analyses.C'("out2.root")'

#include "TStopwatch.h"
TStopwatch timer;

// momentum rotation big range 

//int analyses3(TString infile="/Users/ahmed/GlueX_DIRC_Calib/ok/new/all.root"){// outFile_v3.root
int analyses3_1(TString infile="/Users/ahmed/GlueX_DIRC_Calib/rotat.root"){// outFile_v3.root
    
    timer.Start();
    
    glx_savepath="data3";
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    glx_initDigi();
    
    // calculate cherenkov angle
    // Double_t momentum=3.5;
    Int_t pdg[]= {11,13,211,321,2212};
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    //Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.009),range(5*sigma),noise(0.3);
    
    Double_t fAngleK[10]={0};
    Double_t fAnglePi[10]={0};
    
    Double_t momentum[] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
    //Double_t momentum[] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
    for(int f=0;f<10;f++){
        fAngleK[f] = acos(sqrt(momentum[f]*momentum[f] + mass[3]*mass[3])/momentum[f]/1.4738)-0.00;
        fAnglePi[f]= acos(sqrt(momentum[f]*momentum[f] + mass[2]*mass[2])/momentum[f]/1.4738)-0.00;
    }
    
    
    // histograms
    
    
    TH2F * histo_rotation_map_mean[26];
    TH2F * histo_rotation_map_sigma[26];
    TH2F * histo_rotation_map_occu[26];
    
    TH2F * ratio_mean[26];
    
    
    for(Int_t i=0; i<26; i++) {
        histo_rotation_map_mean[i]= new TH2F( Form("histo_rotation_map_mean_%d",i) , Form("Cherenkov Shift Mean @ Bar %d ; X rotation [mrad]; Y roation [mrad]",i), 20, -20, 20, 20, -20, 20);
        histo_rotation_map_sigma[i]= new TH2F( Form("histo_rotation_map_sigma_%d",i) , Form("Cherenkov Shift Sigma @ Bar %d ; X rotation [mrad]; Y roation [mrad]",i), 20, -20, 20, 20, -20, 20);
        histo_rotation_map_occu[i]= new TH2F( Form("histo_rotation_map_occu_%d",i) , Form("Cherenkov Shift Occupancy @ Bar %d ; X rotation [mrad]; Y roation [mrad]",i), 20, -20, 20, 20, -20, 20);
        
        ratio_mean[i]= new TH2F( Form("ratio_mean_mom_%d",i) , Form("Cherenkov Shift Occupancy @ Bar %d ; X rotation [mrad]; Y roation [mrad]",i),  20, -20, 20, 20, -20, 20);
    }
    
    TH1F* histo_rotation_cell_shift[26][42][42];
    TH1F* histo_rotation_cell_spr[26][42][42];
    
    for(Int_t i=0; i<26; i++) {
        for(Int_t j=0; j<42; j++){
            for(Int_t k=0; k<42; k++){
                histo_rotation_cell_shift[i][j][k] = new TH1F(Form("histo_rotation_cell_shift_%d_%d_%d",i,j,k),Form("Bar %d X rotation %d Y rotation %d ; Measured - Expected [mrad]; Entries [#]",i,j,k) ,100,-20,20);
                histo_rotation_cell_spr[i][j][k] = new TH1F(Form("histo_rotation_cell_spr_%d_%d_%d",i,j,k),Form("Bar %d X rotation %d Y rotation %d ; Measured - Expected [mrad]; Entries [#]",i,j,k) ,100,-20,20);
            }
        }
    }
    
    
    
    // variables
    double diff(-1);
    
    double mean_min(0.818), mean_max(0.834); // 0.817,0.8348);
    double spr_min(5.1), spr_max(11.5); //  6 11.5
    double calc_trk_res(-1);
    
    int x_pos_bin(-1),y_pos_bin(-1);
    double content_histo_pos_xy(-1),content_histo_pos_xy_tmp(-1),average_bin(-1);
    double track_resolution(-1),track_resolution_error(-1), yield_BinCenter(-1);
    double track_pos_resolution(-1),track_pos_resolution_error(-1), pos_BinCenter(-1);
    
    double track_spr_bin(-1),track_spr_error(-1);
    double track_mean_bin(-1),track_mean_error(-1);
    double content_histo_pos_xy_reso_tmp(-1);
    double fit_quality(-1);
    
    int mom_bin_flag(-1);
    int spr_bin_flag= 0;
    
    double mass_phi_mini(1.01),mass_phi_max(1.028);
    double chisq_phi_mini(0),chisq_phi_max(20);
    double miss_mass_phi_mini(-0.005),miss_mass_phi_max(0.005);
    
    //double mass_rho_mini(0.66),mass_rho_max(0.82);
    double mass_rho_mini(0.6),mass_rho_max(0.9);
    
    double chisq_rho_mini(0),chisq_rho_max(10);
    double miss_mass_rho_mini(-0.005),miss_mass_rho_max(0.005);
    
    
    
    // read tree
    
    TFile *f = new TFile(infile);
    TTree *tree_variables = (TTree*)f->Get("tree_variables");
    double track_spr(-1),track_mean(-1), track_yield(-1), track_mom(-1), track_xbar(0),track_ybar(0);
    double track_fit_chisqu(-1),track_fit_NDF(-1);
    int track_pid(-1), track_nbar(-1);
    
    double track_inv_mass(-1),track_missing_mass(-1),track_chi_square(-1),track_TofTrackDist(-1);
    double track_xrotate(0),track_yrotate(0);
    
    tree_variables->SetBranchAddress("track_pid",&track_pid);
    tree_variables->SetBranchAddress("track_spr",&track_spr);
    tree_variables->SetBranchAddress("track_mean",&track_mean);
    tree_variables->SetBranchAddress("track_yield",&track_yield);
    tree_variables->SetBranchAddress("track_mom",&track_mom);
    tree_variables->SetBranchAddress("track_xbar",&track_xbar);
    tree_variables->SetBranchAddress("track_ybar",&track_ybar);
    tree_variables->SetBranchAddress("track_nbar",&track_nbar);
    
    
    tree_variables->SetBranchAddress("track_inv_mass",&track_inv_mass);
    tree_variables->SetBranchAddress("track_missing_mass",&track_missing_mass);
    tree_variables->SetBranchAddress("track_chi_square",&track_chi_square);
    
    tree_variables->SetBranchAddress("track_xrotate",&track_xrotate);
    tree_variables->SetBranchAddress("track_yrotate",&track_yrotate);
    
    
    
    Long64_t nentries = tree_variables->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        if(track_pid==2){
            if(track_inv_mass<mass_rho_mini || track_inv_mass> mass_rho_max)continue;
            if(track_missing_mass < miss_mass_rho_mini || track_missing_mass > miss_mass_rho_max)continue;
            if(track_chi_square<chisq_rho_mini || track_chi_square> chisq_rho_max)continue;
        }
        
        if(track_pid==3){
            if(track_inv_mass<mass_phi_mini || track_inv_mass> mass_phi_max)continue;
            if(track_missing_mass < miss_mass_phi_mini || track_missing_mass > miss_mass_phi_max)continue;
            if(track_chi_square<chisq_phi_mini || track_chi_square> chisq_phi_max)continue;
        }
        
        if(track_pid !=2 ) continue; // select pion !=2
        // mean SPR cut
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        
        double ExAnglePi= acos(sqrt(track_mom*track_mom + mass[2]*mass[2])/track_mom/1.4738);
        double ExMeandiff = track_mean - ExAnglePi;
        
        int barnum = track_nbar;
        
        
        //int  xrotat_flag = (track_xrotate + 0.010)*1000 ;
        //int  yrotat_flag = (track_yrotate + 0.010)*1000 ;
        
        int  xrotat_flag = (track_xrotate + 0.020)*1000 ;
        int  yrotat_flag = (track_yrotate + 0.020)*1000 ;
        
        
        //cout<<"##### start analyses "<<barnum<<" "<<xrotat_flag<<" "<<yrotat_flag<<endl;
        
        histo_rotation_cell_shift[barnum][xrotat_flag][yrotat_flag]->Fill(ExMeandiff*1000);
        histo_rotation_cell_spr[barnum][xrotat_flag][yrotat_flag]->Fill(track_spr*1000);
        
        
    }
    
    cout<<"##### start analyses "<<endl;
    TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
    
    if(false){
        TF1 *fit_gause = new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
        fit_gause->SetLineColor(kBlack);
        fit_gause->SetParameters(100,9,2);
        fit_gause->SetParNames("p0","mean ","sigma");
        fit_gause->SetParLimits(0,0.1,1E6);
        fit_gause->SetParLimits(1,7,11);
        fit_gause->SetParLimits(2,1,10);
        
        
        TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
        for(Int_t i=0; i<26; i++) {
            for(Int_t j=0; j<42; j++){
                for(Int_t k=0; k<42; k++){
                    
                    if(histo_rotation_cell_spr[i][j][k]->GetEntries()<100)continue;
                    histo_rotation_cell_spr[i][j][k]->Fit("fit_gause","M","", -20, 20);
                    double trck_mean_fit = fit_gause->GetParameter(1);
                    double trck_sigma_fit = fit_gause->GetParameter(2);
                    
                    //double trck_mean_fit = histo_rotation_cell_shift[i][j][k]->GetMean();
                    //double trck_sigma_fit = histo_rotation_cell_shift[i][j][k]->GetRMS();
                    
                    if(histo_rotation_cell_spr[i][j][k]->GetEntries()<100){
                        trck_mean_fit = -1000;
                        trck_sigma_fit = 0;
                        
                    }
                    //                    if(trck_mean_fit != -1000){
                    //                        histo_rotation_map_mean[i]->Fill(j-10,k-10,trck_mean_fit);
                    //                        histo_rotation_map_occu[i]->Fill(j-10,k-10);
                    //                    }
                    if(trck_mean_fit != -1000){
                        histo_rotation_map_mean[i]->Fill(j-20,k-20,trck_mean_fit);
                        histo_rotation_map_occu[i]->Fill(j-20,k-20);
                    }
                    //                    cc1->cd();
                    //                    cc1->Update();
                    //                    histo_rotation_cell_spr[i][j][k]->Draw();
                    //                    cc1->Update();
                    //                    TLine *lineMeanSHift= new TLine(0,0,0,1000);
                    //                    lineMeanSHift->SetX1(fit_gause->GetParameter(1));
                    //                    lineMeanSHift->SetX2(fit_gause->GetParameter(1));
                    //                    lineMeanSHift->SetY1(gPad->GetUymin());
                    //                    lineMeanSHift->SetY2(gPad->GetUymax());
                    //                    lineMeanSHift->SetLineColor(kRed);
                    //                    lineMeanSHift->Draw();
                    //                    cc1->Update();
                    //                    cc1->WaitPrimitive();
                    
                }
            }
        }
        
        
        int counter_MeanShiftYield =0;
        for(Int_t k=0; k<26; k++){
            cc1->Clear();
            cc1->cd();
            cc1->Update();
            //if(k!=3) continue;
            ratio_mean[k] = (TH2F*)histo_rotation_map_mean[k]->Clone();
            ratio_mean[k]->SetTitle(Form("SPR @ Bar %d", k));
            ratio_mean[k]->Divide(histo_rotation_map_occu[k]);
            //        cc2->cd();
            //        cc2->Update();
            //        ratio_mean[k]->SetMinimum(-5);
            //        ratio_mean[k]->SetMaximum(5);
            //        ratio_mean[k]->Draw("COLZ");
            //        cc2->Update();
            //        cc2->WaitPrimitive();
            
            TString num_string=Form("_%d",counter_MeanShiftYield);
            //glx_canvasAdd("r_SPR"+num_string,800,400);
            ratio_mean[k]->SetMinimum(0);
            ratio_mean[k]->SetMaximum(15);
            ratio_mean[k]->Draw("COLZ");
            //glx_canvasGet("r_SPR"+num_string)->Update();
            
            cc1->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/histo3/spr_rotation_%d.png",counter_MeanShiftYield));
            cc1->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/histo3/spr_rotation_%d.root",counter_MeanShiftYield));
            
            ++counter_MeanShiftYield;
        }
    }
    
    
    
    if(true){
        
        TF1 *fit_gause = new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
        fit_gause->SetLineColor(kBlack);
        fit_gause->SetParameters(100,9,2);
        fit_gause->SetParNames("p0","mean ","sigma");
        fit_gause->SetParLimits(0,0.1,1E6);
        fit_gause->SetParLimits(1,-10,10);
        fit_gause->SetParLimits(2,1,10);
        
        
        TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
        for(Int_t i=0; i<26; i++) {
            for(Int_t j=0; j<42; j++){
                for(Int_t k=0; k<42; k++){

                    //cout<<"##### No Problem 1 "<<endl;
                    //if(histo_rotation_cell_shift[i][j][k]->GetEntries()<100)continue;
                    

                    //histo_rotation_cell_shift[i][j][k]->Fit("fit_gause","M","", -20, 20);
                    //double trck_mean_fit = fit_gause->GetParameter(1);
                    //double trck_sigma_fit = fit_gause->GetParameter(2);
                    
                    double trck_mean_fit = histo_rotation_cell_shift[i][j][k]->GetMean();
                    double trck_sigma_fit = histo_rotation_cell_shift[i][j][k]->GetRMS();

                    if(histo_rotation_cell_shift[i][j][k]->GetEntries()<5){
                        trck_mean_fit = -1000;
                        trck_sigma_fit = 0;

                    }

                    if(trck_mean_fit != -1000){
                        histo_rotation_map_mean[i]->Fill(j-20,k-20,trck_mean_fit);
                        histo_rotation_map_occu[i]->Fill(j-20,k-20);
                    }

                    
                    
                    //                    cc1->cd();
                    //                    cc1->Update();
                    //                    histo_rotation_cell_shift[i][j][k]->Draw();
                    //                    cc1->Update();
                    //                    TLine *lineMeanSHift= new TLine(0,0,0,1000);
                    //                    lineMeanSHift->SetX1(trck_mean_fit);
                    //                    lineMeanSHift->SetX2(trck_mean_fit);
                    //                    lineMeanSHift->SetY1(gPad->GetUymin());
                    //                    lineMeanSHift->SetY2(gPad->GetUymax());
                    //                    lineMeanSHift->SetLineColor(kRed);
                    //                    lineMeanSHift->Draw();
                    //                    cc1->Update();
                    //                    cc1->WaitPrimitive();

                }
            }
        }
        cout<<"##### No Problem 2 "<<endl;
        
        int counter_MeanShiftYield =0;
        for(Int_t k=0; k<26; k++){
            cc1->Clear();
            cc1->cd();
            cc1->Update();
            //if(k!=3) continue;
            ratio_mean[k] = (TH2F*)histo_rotation_map_mean[k]->Clone();
            ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ Bar %d", k));
            ratio_mean[k]->Divide(histo_rotation_map_occu[k]);
            //        cc2->cd();
            //        cc2->Update();
            //        ratio_mean[k]->SetMinimum(-5);
            //        ratio_mean[k]->SetMaximum(5);
            //        ratio_mean[k]->Draw("COLZ");
            //        cc2->Update();
            //        cc2->WaitPrimitive();

            TString num_string=Form("_%d",counter_MeanShiftYield);
            //glx_canvasAdd("r_MeanShiftYield"+num_string,800,400);
            ratio_mean[k]->SetMinimum(-6);
            ratio_mean[k]->SetMaximum(10);
            ratio_mean[k]->Draw("COLZ");
            //glx_canvasGet("r_MeanShiftYield"+num_string)->Update();
            
            cc1->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/histo3/mean_shift_%d.png",counter_MeanShiftYield));
            cc1->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/histo3/mean_shift_%d.root",counter_MeanShiftYield));
            ++counter_MeanShiftYield;
        }
    }
    
    //glx_canvasSave(2,0);
    //glx_canvasDel("*");
    
    
    cout<<"##### No Problem 3 "<<endl;
    cout<<"####### @ 3.5 GeV/c fAngleK "<< fAngleK[1]<<"  fAnglePi "<<fAnglePi[1]<<endl;
    
    timer.Stop();
    
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    
    
    return 0;
}
