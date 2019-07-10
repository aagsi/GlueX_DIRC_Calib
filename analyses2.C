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

int analyses2(TString infile="out2.root"){// outFile_v3.root
    
    timer.Start();
    
    glx_savepath="data";
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
    const int nbar =26;
    const int pos_min(-100), pos_max(100);
    const int nbin_yield =100;
    const int nbin_mom =10;
    
    TH1F * histo_cherenkov = new TH1F("histo_cherenkov","histo_cherenkov", 100,0.6,1);
    TH1F * histo_tdiff= new TH1F("histo_tdiff","histo_tdiff", 500,-10,10);
    
    TH1F *hist_ev_rho_mass = new TH1F("hist_ev_rho_mass","; #pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.3, 1.2);
    TH1F *hist_ev_phi_mass = new TH1F("hist_ev_phi_mass","; k^{#plus}k^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.9, 1.2);
    TH1F *hist_ev_missing_mass_phi = new TH1F("hist_ev_missing_mass_phi",";#phi Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.03, 0.03);
    TH1F *hist_ev_missing_mass_rho = new TH1F("hist_ev_missing_mass_rho",";#rho Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.03, 0.03);
    TH1F *hist_ev_chi_phi = new TH1F("hist_ev_chi_phi","; #phi Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    TH1F *hist_ev_chi_rho = new TH1F("hist_ev_chi_rho","; #rho Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    
    hist_ev_rho_mass->SetTitle("#rho Invariant Mass");
    hist_ev_phi_mass->SetTitle("#phi Invariant Mass");
    
    
    TH1F* histo_tmp_pos = new TH1F("histo_tmp_pos","; X Bar Hit [cm]; entries [#]",40,-100,100);
    TH1F* histo_nbar = new TH1F("histo_nbar","; X Bar Hit [cm]; entries [#]",24,0,24);
    TH1F* histo_mom = new TH1F("histo_mom","; Track Momentum [GeV/c]; entries [#]",20,0,20);
    TH1F* histo_yield = new TH1F("histo_yield","; Photon Yield [cm]; entries [#]",100,0,100);
    
    
    TH1F* histo_delta_shift_mom[42][26][16];
    TH1F* histo_delta_shift_spr[42][26][16];
    TH1F* histo_delta_shift_nph[42][26][60];
    
    for(Int_t i=0; i<42; i++) {
        for(Int_t j=0; j<26; j++){
            for(Int_t k=0; k<16; k++){
                histo_delta_shift_mom[i][j][k] = new TH1F(Form("histo_delta_shift_%d_%d_%d",i,j,k),Form("xbin %d nbar %d mom flag %d ; Measured - Expected [mrad]; Entries [#]",i,j,k) ,100,-20,20);
                histo_delta_shift_spr[i][j][k] = new TH1F(Form("histo_delta_shift_spr_%d_%d_%d",i,j,k),Form("xbin %d nbar %d spr flag %d ; Measured - Expected [mrad]; Entries [#]",i,j,k) ,100,-20,20);
            }
            for(Int_t m=0; m<60; m++){
                histo_delta_shift_nph[i][j][m] = new TH1F(Form("histo_delta_shift_nph_%d_%d_%d",i,j,m),Form("xbin %d nbar %d nph %d ; Measured - Expected [mrad]; Entries [#]",i,j,m) ,100,-20,20);
                
            }
        }
    }
    
    
    //    TH1F* histo_delta_shift_yield[42][26][22][50];
    //    for(Int_t i=0; i<42; i++) { // x segment bin
    //        for(Int_t j=0; j<26; j++){ // bar num
    //            for(Int_t k=0; k<22; k++){ // mom bin
    //                for(Int_t l=0; l<50; l++){// yield
    //                    histo_delta_shift_yield[i][j][k][l] = new TH1F(Form("histo_delta_shift_yield_%d_%d_%d_%d",i,j,k,l),Form("xbin %d nbar %d mom flag %d yield %d; Measured - Expected [mrad]; Entries [#]",i,j,k,l) ,100,-20,20);
    //                }
    //            }
    //        }
    //    }
    
    TH2F * histo_pos_xy_TrkMeanShift_average = new TH2F( "histo_pos_xy_TrkMeanShift_average" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_pos_xy_TrkMeanShift_occu_average = new TH2F( "histo_pos_xy_TrkMeanShift_occu_average" , "; Bar Hit X [cm]; Bar number [#]",40, 0, 40, 24, 0, 24);
    TH2F * histo_pos_xy_TrkSigmaShift_average = new TH2F( "histo_pos_xy_TrkSigmaShift_average" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_pos_xy_TrkSigmaShift_occu_average = new TH2F( "histo_pos_xy_TrkSigmaShift_occu_average" , "; Bar Hit X [cm]; Bar number [#]",40, 0, 40, 24, 0, 24);
    TH2F * ratio_mean_average = new TH2F( "ratio_mean_average" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * ratio_sigma_average = new TH2F( "ratio_sigma_average" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    
    
    
    TH2F * histo_pos_xy_TrkMeanShift[60];
    TH2F * histo_pos_xy_TrkMeanShift_occu[60];
    TH2F * ratio_mean[60];
    for(Int_t i=0; i<60; i++) {
        histo_pos_xy_TrkMeanShift[i]= new TH2F( Form("histo_pos_xy_TrkMeanShift_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
        histo_pos_xy_TrkMeanShift_occu[i]= new TH2F( Form("histo_pos_xy_TrkMeanShift_occu_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
        ratio_mean[i]= new TH2F( Form("ratio_mean_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
    }
    
    
    TH2F * histo_pos_xy_TrkSigmaShift[60];
    TH2F * histo_pos_xy_TrkSigmaShift_occu[60];
    TH2F * ratio_sigma[60];
    for(Int_t i=0; i<60; i++) {
        histo_pos_xy_TrkSigmaShift[i]= new TH2F( Form("histo_pos_xy_TrkSigmaShift_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
        histo_pos_xy_TrkSigmaShift_occu[i]= new TH2F( Form("histo_pos_xy_TrkSigmaShift_occu_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
        ratio_sigma[i]= new TH2F( Form("ratio_sigma_mom_%d",i) , Form("@ Momentum %d ; Bar Hit X [cm]; Bar number [#]",i), 40, 0, 40, 24, 0, 24);
    }
    
    
    
    //
    //    TH1F*  histo_track_yield_bar_mom[26][10];
    //    for(Int_t i=0; i<26; i++) {
    //        for(Int_t j=0; j<10; j++) {
    //            int k =j+2;
    //            int l= j+3;
    //            histo_track_yield_bar_mom[i][j] = new TH1F(Form("histo_track_yield_bar_mom_%d_%d",i,j),Form("bar %d mom %d - %d GeV/c ;Photon Yield ; Photon Yield; Entries [#]",i,k,l) ,100 ,0,100);
    //        }
    //    }
    //
    //    // independant
    //    const int pos_bin_shift(50); // best 200
    //    TH2F * histo_pos_xy_occupancy = new TH2F( "histo_pos_xy_occupancy" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //    TH2F * histo_pos_xy_shift_tmp = new TH2F( "histo_pos_xy_shift_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //    TH2F * histo_pos_xy_shift = new TH2F( "histo_pos_xy_shift" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //    TH2F * histo_pos_xy_shiftEx_tmp = new TH2F( "histo_pos_xy_shiftEx_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //    TH2F * histo_pos_xy_shiftEx = new TH2F( "histo_pos_xy_shiftEx" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //    TH2F * histo_pos_xy_shiftEx_tmp_positive = new TH2F( "histo_pos_xy_shiftEx_tmp_positive" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //    TH2F * histo_pos_xy_shiftEx_positive = new TH2F( "histo_pos_xy_shiftEx_positive" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //    TH2F * histo_pos_xy_shiftEx_tmp_negative  = new TH2F( "histo_pos_xy_shiftEx_tmp_negative " , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //    TH2F * histo_pos_xy_shiftEx_negative  = new TH2F( "histo_pos_xy_shiftEx_negative " , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //    TH2F * histo_pos_xy_occupancy_postiveShift = new TH2F( "histo_pos_xy_occupancy_postiveShift" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //    TH2F * histo_pos_xy_occupancy_negativeShift = new TH2F( "histo_pos_xy_occupancy_negativeShift" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    //
    //
    //    //
    //    TH1F* histo_track_mean = new TH1F("histo_track_mean","; Mean per Track [rad]; entries [#]",250,0.817,0.8348);
    //    TH1F* histo_track_spr = new TH1F("histo_track_spr","; SPR per Track [m rad]; entries [#]",250,5.1,20);
    //
    //    //
    //    TH1F*  histo_track_yield[10];
    //    for(Int_t i=0; i<10; i++){
    //        int kk = i+2;
    //        int ll =i+3;
    //        histo_track_yield[i] = new TH1F(Form("histo_track_yield_%d",i), Form("Photon Yield @ momentum %d - %d GeV/c ; Photon Yield; Entries [#]",kk,ll) ,100 ,0,100);
    //    }
    //
    //
    //    TH1F*  histo_track_mean_bar_mom[26][10];
    //    for(Int_t i=0; i<26; i++){
    //        for(Int_t j=0; j<10; j++){
    //            int k = j+2;
    //            int l =j+3;
    //            histo_track_mean_bar_mom[i][j] = new TH1F(Form("histo_track_mean_bar_mom_%d_%d_%d",i,k,l), Form("Track Mean @ bar %d mom %d - %d ; Mean [rad]; Entries [#]",i,k,l) ,50,0.817,0.8348);
    //        }
    //    }
    
    
    // graphs
    
    
    
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
    /*
     std::vector<int> *vpx = 0;
     std::vector<int> *vpy = 0;
     std::vector<int> *vpz = 0;
     std::vector<double> *vtdiff = 0;
     std::vector<double> *vtangle = 0;
     
     TBranch *bvpx = 0;
     TBranch *bvpy = 0;
     TBranch *bvpz = 0;
     TBranch *bvtdiff = 0;
     TBranch *bvtangle = 0;
     */
    
    tree_variables->SetBranchAddress("track_pid",&track_pid);
    tree_variables->SetBranchAddress("track_spr",&track_spr);
    tree_variables->SetBranchAddress("track_mean",&track_mean);
    tree_variables->SetBranchAddress("track_yield",&track_yield);
    tree_variables->SetBranchAddress("track_mom",&track_mom);
    tree_variables->SetBranchAddress("track_xbar",&track_xbar);
    tree_variables->SetBranchAddress("track_ybar",&track_ybar);
    tree_variables->SetBranchAddress("track_nbar",&track_nbar);
    //tree_variables->SetBranchAddress("track_fit_chisqu",&track_fit_chisqu);
    //tree_variables->SetBranchAddress("track_fit_NDF",&track_fit_NDF);
    
    tree_variables->SetBranchAddress("track_inv_mass",&track_inv_mass);
    tree_variables->SetBranchAddress("track_missing_mass",&track_missing_mass);
    tree_variables->SetBranchAddress("track_chi_square",&track_chi_square);
    //tree_variables->SetBranchAddress("track_TofTrackDist",&track_TofTrackDist);
    
    /*
     tree_variables->SetBranchAddress("vpx",&vpx,&bvpx);
     tree_variables->SetBranchAddress("vpy",&vpy,&bvpy);
     tree_variables->SetBranchAddress("vpz",&vpz,&bvpz);
     
     tree_variables->SetBranchAddress("vtdiff",&vtdiff,&bvtdiff);
     tree_variables->SetBranchAddress("vtangle",&vtangle,&bvtangle);
     */
    
    
    /*
     ///////////////////////////////////////
     ///////////////////////////////////////
     ///////////////////////////////////////
     double pion_counter(1),kaon_counter(1);
     double sum1,sum2,noise = 0.5;
     TSpectrum *spect = new TSpectrum(10);
     double minChangle = 0.6;
     double maxChangle = 0.9;
     TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
     double cut_cangle=0.04;
     double cherenkovreco[5],spr[5];
     
     TGaxis::SetMaxDigits(3);
     double sigma[]={0.01,0.01,0.01,0.010,0.01,0.01};
     double mAngle[5];
     TF1  *fAngle[5];
     TH1F *hAngle[5], *hLnDiff[5];
     double init_mon=3.5;
     for(int i=0; i<5; i++){
     hAngle[i]= new TH1F(Form("hAngle_%d",i),  "cherenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
     //hNph[i] = new TH1F(Form("hNph_%d",i),";detected photons [#];entries [#]", 150,0,150);
     
     mAngle[i] = acos(sqrt(init_mon * init_mon + glx_mass[i]*glx_mass[i])/init_mon/1.473);
     fAngle[i] = new TF1(Form("fAngle_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
     fAngle[i]->SetParameter(0,1);        // const
     fAngle[i]->SetParameter(1,mAngle[i]);// mean
     fAngle[i]->SetParameter(2,sigma[i]); // sigma
     hAngle[i]->SetMarkerStyle(20);
     hAngle[i]->SetMarkerSize(0.8);
     hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",110,-120,120); //,80,-150,150);
     //hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",80,-150,150);
     }
     hAngle[2]->SetLineColor(4);
     hAngle[3]->SetLineColor(2);
     hAngle[2]->SetMarkerColor(kBlue+1);
     hAngle[3]->SetMarkerColor(kRed+1);
     fAngle[2]->SetLineColor(4);
     fAngle[3]->SetLineColor(2);
     hLnDiff[2]->SetLineColor(4);
     hLnDiff[3]->SetLineColor(2);
     ///////////////////////////////////////
     ///////////////////////////////////////
     ///////////////////////////////////////
     */
    Long64_t nentries = tree_variables->GetEntries();
    //Double_t mean_array[]={0.824512,0.825366,0.826363,0.826648,0.826861,0.826576 ,0.826149}; // 1st version
    //Double_t mean_array[]={0.8245,0.826,0.8265,0.82685,0.8269,0.8269,0.8269}; // 2nd version
    
    // 3rd version
    //double diff_array_detailed[26][10]={0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165};
    
    
    
    
    // 4D resolution moentum position
    //double trk_xpos,trk_mom,trk_shift;
    //int trk_nbar;
    //auto f2 = TFile::Open("4d_diff.root","RECREATE");
    //TNtuple ntuple("ntuple","data from ascii file","trk_xpos:trk_nbar:trk_mom:trk_shift");
    
    
    
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        /*
         bvpx->GetEntry(i);
         bvpy->GetEntry(i);
         bvpz->GetEntry(i);
         bvtdiff->GetEntry(i);
         bvtangle->GetEntry(i);
         */
        
        if(track_pid==2){
            hist_ev_rho_mass->Fill(track_inv_mass);
            hist_ev_missing_mass_rho->Fill(track_missing_mass);
            hist_ev_chi_rho->Fill(track_chi_square);
            if(track_inv_mass<mass_rho_mini || track_inv_mass> mass_rho_max)continue;
            if(track_missing_mass < miss_mass_rho_mini || track_missing_mass > miss_mass_rho_max)continue;
            if(track_chi_square<chisq_rho_mini || track_chi_square> chisq_rho_max)continue;
        }
        
        if(track_pid==3){
            hist_ev_phi_mass->Fill(track_inv_mass);
            hist_ev_missing_mass_phi->Fill(track_missing_mass);
            hist_ev_chi_phi->Fill(track_chi_square);
            if(track_inv_mass<mass_phi_mini || track_inv_mass> mass_phi_max)continue;
            if(track_missing_mass < miss_mass_phi_mini || track_missing_mass > miss_mass_phi_max)continue;
            if(track_chi_square<chisq_phi_mini || track_chi_square> chisq_phi_max)continue;
        }
        
        //if(!(track_mom> 4.5   || track_mom<3.5 ))continue;
        //if(!(track_yield==30 ))continue;
        
        //if(track_nbar!=5) continue;
        
        /*
         // wall cut
         if(track_nbar<4 || track_nbar>8 ) continue;
         if(track_xbar>10 || track_xbar < -10 ) continue;
         
         // momentum cut
         if(track_mom<3.5 || track_mom>4.0 )continue;
         
         //////////////////////////////////////////
         //////// calculate separation power //////
         //////////////////////////////////////////
         double percentage = kaon_counter/pion_counter*100.0;
         //cout<<"##### percentage  "<<percentage<<endl;
         if (!(percentage <100 && track_pid==2) || track_pid==3){
         sum1=0;
         sum2=0;
         for (UInt_t j = 0; j < vtangle->size(); ++j){
         double time_diff = vtdiff->at(j);
         if(fabs(time_diff)>3)continue;
         double tangle = vtangle->at(j);
         
         for(int p=0; p<5; p++){
         mAngle[p] = acos(sqrt(track_mom * track_mom + glx_mass[p]*glx_mass[p])/track_mom/1.473);
         fAngle[p]->SetParameter(1,mAngle[p]);
         }
         hAngle[track_pid]->Fill(tangle);
         if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))>cut_cangle)continue;
         
         sum1 += TMath::Log(fAngle[2]->Eval(tangle)+noise);
         sum2 += TMath::Log(fAngle[3]->Eval(tangle)+noise);
         }
         
         double sum = sum1-sum2;
         //cout<<"##### sum  "<<sum<<endl;
         hLnDiff[track_pid]->Fill(sum);
         
         if(track_pid==3) kaon_counter++;
         if(track_pid==2) pion_counter++;
         
         
         }
         
         continue;
         */
        
        //if(track_nbar>5) continue;
        if(track_pid !=2 ) continue; // select pion !=2
        //fit_quality=track_fit_chisqu/track_fit_NDF;
        //histo_chiNDF->Fill(fit_quality);
        
        //if(track_fit_chisqu>100)continue;
        //if(track_yield<10)continue;
        
        
        
        // momentum cut
        //if(track_mom> 4.5   || track_mom<3.5) continue;
        //if(track_mom<3.5 || track_mom>4.0 )continue;
        // mean SPR cut
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        // wall cut
        //if(track_nbar<4 || track_nbar>8 ) continue;
        //if(track_xbar>10 || track_xbar < -10 ) continue;
        
        //if(fit_quality>2) continue;
        //histo_chiNDF_cut->Fill(fit_quality);
        
        // mean tail study
        //if(track_mean< mean_max) continue;
        //if(track_spr*1000< spr_max) continue;
        
        //        ////////////////////////////////
        //        ////// New Implimentation //////
        //        ////////////////////////////////
        //
        //        for (UInt_t j = 0; j < vpx->size(); ++j) {
        //            glx_hdigi[vpx->at(j)]->Fill(vpy->at(j),vpz->at(j));
        //        }
        //        for (UInt_t j = 0; j < vtangle->size(); ++j) {
        //            double time_diff = vtdiff->at(j);
        //            if(fabs(time_diff)>3)continue;
        //            histo_tdiff->Fill(time_diff);
        //            histo_cherenkov->Fill(vtangle->at(j));
        //        }
        
        
        
        //diff = fAnglePi-track_mean; // not used because there are systematic shifts
        diff = track_mean-   0.82608; // default value will be changed
        
        if(track_mom>1 && track_mom<1.5){
            mom_bin_flag=0;
            //diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>1.5 && track_mom<2){
            mom_bin_flag=1;
        }
        else if(track_mom>2 && track_mom<2.5){
            mom_bin_flag=2;
        }
        else if(track_mom>2.5 && track_mom<3){
            mom_bin_flag=3;
        }
        else if(track_mom>3 && track_mom<3.5){
            mom_bin_flag=4;
        }
        else if(track_mom>3.5 && track_mom<4){
            mom_bin_flag=5;
        }
        else if(track_mom>4 && track_mom<4.5){
            mom_bin_flag=6;
        }
        else if(track_mom>4.5 && track_mom<5){
            mom_bin_flag=7;
        }
        else if(track_mom>5 && track_mom<5.5){
            mom_bin_flag=8;
        }
        else if(track_mom>5.5 && track_mom<6){
            mom_bin_flag=9;
        }
        else if(track_mom>6 && track_mom<6.5){
            mom_bin_flag=10;
        }
        else if(track_mom>6.5 && track_mom<7){
            mom_bin_flag=11;
        }
        else if(track_mom>7 && track_mom<7.5){
            mom_bin_flag=12;
        }
        else if(track_mom>7.5 && track_mom<8){
            mom_bin_flag=13;
        }
        else if(track_mom>8 && track_mom<8.5){
            mom_bin_flag=14;
        }
        else{
            continue;
        }
        
        //        // SPR Flag
        if(track_spr>0.006 && track_spr<0.007){
            spr_bin_flag=0;
        }
        else if(track_spr>0.007 && track_spr<0.008){
            spr_bin_flag=1;
        }
        else if(track_spr>0.008 && track_spr<0.009){
            spr_bin_flag=2;
        }
        else if(track_spr>0.009 && track_spr<0.010){
            spr_bin_flag=3;
        }
        else if(track_spr>0.010 && track_spr<0.011){
            spr_bin_flag=4;
        }
        else if(track_spr>0.011 && track_spr<0.012){
            spr_bin_flag=5;
        }
        else if(track_spr>0.012 && track_spr<0.013){
            spr_bin_flag=6;
        }
        else if(track_spr>0.013 && track_spr<0.014){
            spr_bin_flag=7;
        }
        else{
            continue;
        }
        //cout<<spr_bin_flag<<endl;
        
        
        //        histo_track_mean->Fill(track_mean);
        //        histo_track_spr->Fill(track_spr*1000);
        //        histo_track_yield_bar_mom[track_nbar][mom_bin_flag]->Fill(track_yield);
        //        histo_track_mean_bar_mom[track_nbar][mom_bin_flag]->Fill(track_mean);
        
        //        // shift map
        
        //        if(true){
        //            histo_pos_xy_occupancy->Fill(track_xbar,track_ybar);
        //            x_pos_bin = histo_pos_xy_occupancy->GetXaxis()->FindBin(track_xbar);
        //            y_pos_bin = histo_pos_xy_occupancy->GetYaxis()->FindBin(track_ybar);
        //            content_histo_pos_xy=histo_pos_xy_occupancy->GetBinContent(x_pos_bin,y_pos_bin);
        //
        //
        //            histo_pos_xy_shift_tmp->Fill(track_xbar,track_ybar,diff*1000);
        //            content_histo_pos_xy_tmp=histo_pos_xy_shift_tmp->GetBinContent(x_pos_bin,y_pos_bin);
        //            average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        //            histo_pos_xy_shift->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        //
        //
        double ExAnglePi= acos(sqrt(track_mom*track_mom + mass[2]*mass[2])/track_mom/1.4738);
        double ExMeandiff = track_mean - ExAnglePi;
        //            histo_pos_xy_shiftEx_tmp->Fill(track_xbar,track_ybar,ExMeandiff*1000);
        //            content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp->GetBinContent(x_pos_bin,y_pos_bin);
        //            average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        //            histo_pos_xy_shiftEx->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        //
        //            if (ExMeandiff > 0) {
        //
        //                histo_pos_xy_occupancy_postiveShift->Fill(track_xbar,track_ybar);
        //                x_pos_bin = histo_pos_xy_occupancy_postiveShift->GetXaxis()->FindBin(track_xbar);
        //                y_pos_bin = histo_pos_xy_occupancy_postiveShift->GetYaxis()->FindBin(track_ybar);
        //                content_histo_pos_xy=histo_pos_xy_occupancy_postiveShift->GetBinContent(x_pos_bin,y_pos_bin);
        //
        //                histo_pos_xy_shiftEx_tmp_positive->Fill(track_xbar,track_ybar,ExMeandiff*1000);
        //                content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp_positive->GetBinContent(x_pos_bin,y_pos_bin);
        //                average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        //                histo_pos_xy_shiftEx_positive->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        //
        //            }else{
        //                histo_pos_xy_occupancy_negativeShift->Fill(track_xbar,track_ybar);
        //                x_pos_bin = histo_pos_xy_occupancy_negativeShift->GetXaxis()->FindBin(track_xbar);
        //                y_pos_bin = histo_pos_xy_occupancy_negativeShift->GetYaxis()->FindBin(track_ybar);
        //                content_histo_pos_xy=histo_pos_xy_occupancy_negativeShift->GetBinContent(x_pos_bin,y_pos_bin);
        //
        //                histo_pos_xy_shiftEx_tmp_negative ->Fill(track_xbar,track_ybar,ExMeandiff*1000);
        //                content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp_negative ->GetBinContent(x_pos_bin,y_pos_bin);
        //                average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        //                histo_pos_xy_shiftEx_negative ->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        //            }
        //        }
        
        // Fill ntuple
        //double trk_xpos=track_xbar;
        //int trk_nbar=track_nbar;
        //double trk_mom=track_mom;
        //double trk_shift=ExMeandiff;
        //ntuple.Fill(trk_xpos,trk_nbar,trk_mom,trk_shift);
        
        
        
        
        
        
        histo_nbar->Fill(track_nbar);
        histo_mom->Fill(track_mom);
        int x_pos_tmp_bin = histo_tmp_pos->GetXaxis()->FindBin(track_xbar);
        int mom_tmp_bin = histo_mom->GetXaxis()->FindBin(track_mom);
        
        histo_delta_shift_mom[x_pos_tmp_bin][track_nbar][mom_bin_flag]->Fill(ExMeandiff*1000);
        int yield_bin_flag = track_yield;
        histo_delta_shift_nph[x_pos_tmp_bin][track_nbar][yield_bin_flag]->Fill(ExMeandiff*1000);
        //int spr_bin_flag= track_spr*1000;
        histo_delta_shift_spr[x_pos_tmp_bin][track_nbar][spr_bin_flag]->Fill(ExMeandiff*1000);
        
        
        //int yied_int=track_yield;
        //histo_delta_shift_yield[x_pos_tmp_bin][track_nbar][mom_tmp_bin][yied_int]->Fill(ExMeandiff*1000);
    }
    
    cout<<"##### start analyses "<<endl;
    //f2->Write();
    //f2->Close();
    
    //glx_canvasAdd("r_nbar",800,400);
    //histo_nbar->Draw();
    
    //glx_canvasAdd("r_mom",800,400);
    //histo_mom->Draw();
    
    /*
     TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
     TF1 *fit_gause1 = new TF1("fit_gause1","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
     fit_gause1->SetLineColor(kBlack);
     fit_gause1->SetParameters(100,9,2);
     fit_gause1->SetParNames("p0","mean ","sigma");
     fit_gause1->SetParLimits(0,0.1,1E6);
     fit_gause1->SetParLimits(1,-10,10);
     fit_gause1->SetParLimits(2,1,10);
     
     //double shift_array_detailed[42][26][22]={0};
     for(Int_t i=0; i<42; i++) {
     for(Int_t j=0; j<26; j++){
     for(Int_t k=0; k<16; k++){
     
     histo_delta_shift_mom[i][j][k]->Fit("fit_gause1","M","", -20, 20);
     double trck_mean_fit = fit_gause1->GetParameter(1);
     double trck_sigma_fit = fit_gause1->GetParameter(2);
     //shift_array_detailed[i][j[k]=fit_gause1->GetParameter(1);
     if(histo_delta_shift_mom[i][j][k]->GetEntries()<200){
     trck_mean_fit = -1000;
     trck_sigma_fit = 0;
     
     }
     if(trck_mean_fit != -1000){
     histo_pos_xy_TrkMeanShift_average->Fill(i,j,trck_mean_fit);
     histo_pos_xy_TrkMeanShift_occu_average->Fill(i,j);
     }
     
     histo_pos_xy_TrkSigmaShift_average->Fill(i,j,trck_sigma_fit);
     histo_pos_xy_TrkSigmaShift_occu_average->Fill(i,j);
     
     histo_pos_xy_TrkMeanShift[k]->Fill(i,j,trck_mean_fit);
     histo_pos_xy_TrkMeanShift_occu[k]->Fill(i,j);
     
     histo_pos_xy_TrkSigmaShift[k]->Fill(i,j,trck_sigma_fit);
     histo_pos_xy_TrkSigmaShift_occu[k]->Fill(i,j);
     
     
     //                cc1->cd();
     //                cc1->Update();
     //                histo_delta_shift_mom[i][j][k]->Draw();
     //                cc1->Update();
     //                TLine *lineMeanSHift= new TLine(0,0,0,1000);
     //                lineMeanSHift->SetX1(fit_gause1->GetParameter(1));
     //                lineMeanSHift->SetX2(fit_gause1->GetParameter(1));
     //                lineMeanSHift->SetY1(gPad->GetUymin());
     //                lineMeanSHift->SetY2(gPad->GetUymax());
     //                lineMeanSHift->SetLineColor(kRed);
     //                lineMeanSHift->Draw();
     //                cc1->Update();
     //                cc1->WaitPrimitive();
     }
     }
     }
     
     //printing diff array per bar per mom
     //    cout<<"shift_array_detailed[i][j][k]={";
     //    for(Int_t i=0; i<42; i++) {
     //        for(Int_t j=0; j<26; j++){
     //            for(Int_t k=0; k<22; k++){
     //                cout<<","<<fit_gause1->GetParameter(1);
     //            }
     //        }
     //    }
     //    cout<<"};"<<endl;
     
     double v11(0) ,v22(0);
     TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
     for(Int_t k=0; k<16; k++){
     v11 = v11+0.5;
     v22 = v11+0.5;
     ratio_mean[k] = (TH2F*)histo_pos_xy_TrkMeanShift[k]->Clone();
     ratio_mean[k]->GetXaxis()->SetTitle(" ");
     ratio_mean[k]->GetYaxis()->SetTitle(" ");
     ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ Momentun [%1.1f - %1.1f] GeV/c", v11,v22));
     ratio_mean[k]->Divide(histo_pos_xy_TrkMeanShift_occu[k]);
     cc2->cd();
     cc2->Update();
     ratio_mean[k]->SetMinimum(-5);
     ratio_mean[k]->SetMaximum(5);
     ratio_mean[k]->Draw("COLZ");
     cc2->Update();
     cc2->WaitPrimitive();
     }
     
     
     double v1(0) ,v2(0);
     TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
     for(Int_t k=0; k<16; k++){
     v1 = v1+0.5;
     v2 = v1+0.5;
     ratio_sigma[k] = (TH2F*)histo_pos_xy_TrkSigmaShift[k]->Clone();
     ratio_sigma[k]->GetXaxis()->SetTitle(" ");
     ratio_sigma[k]->GetYaxis()->SetTitle(" ");
     ratio_sigma[k]->SetTitle(Form("Sigma of  #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ Momentun [%1.1f - %1.1f] GeV/c",v1, v2));
     ratio_sigma[k]->Divide(histo_pos_xy_TrkSigmaShift_occu[k]);
     cc3->cd();
     cc3->Update();
     ratio_sigma[k]->Draw("COLZ");
     cc3->Update();
     cc3->WaitPrimitive();
     }
     
     glx_canvasAdd("r_ratio_mean",800,400);
     ratio_mean_average = (TH2F*)histo_pos_xy_TrkMeanShift_average->Clone();
     ratio_mean_average->GetXaxis()->SetTitle(" ");
     ratio_mean_average->GetYaxis()->SetTitle(" ");
     ratio_mean_average->SetTitle("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
     ratio_mean_average->Divide(histo_pos_xy_TrkMeanShift_occu_average);
     ratio_mean_average->Draw("COLZ");
     
     //glx_canvasAdd("r_histo_pos_xy_TrkMeanShift",800,400);
     //histo_pos_xy_TrkMeanShift_average->Draw("colz");
     
     glx_canvasAdd("r_ratio_sigma",800,400);
     ratio_sigma_average = (TH2F*)histo_pos_xy_TrkSigmaShift_average->Clone();
     ratio_sigma_average->GetXaxis()->SetTitle(" ");
     ratio_sigma_average->GetYaxis()->SetTitle(" ");
     ratio_sigma_average->SetTitle("Sigma of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
     ratio_sigma_average->Divide(histo_pos_xy_TrkSigmaShift_occu_average);
     ratio_sigma_average->Draw("COLZ");
     
     
     */
    
    /////////////////////////////////////////////////////////////////////////////////////////
    
    /*
     TCanvas *cc4 = new TCanvas("cc4","cc4",800,500);
     TF1 *fit_gause2 = new TF1("fit_gause2","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
     fit_gause2->SetLineColor(kBlack);
     fit_gause2->SetParameters(100,9,2);
     fit_gause2->SetParNames("p0","mean ","sigma");
     fit_gause2->SetParLimits(0,0.1,1E6);
     fit_gause2->SetParLimits(1,-10,10);
     fit_gause2->SetParLimits(2,1,10);
     
     //double shift_array_detailed[42][26][22]={0};
     for(Int_t i=0; i<42; i++) {
     for(Int_t j=0; j<26; j++){
     for(Int_t k=0; k<16; k++){
     
     //if(histo_delta_shift_spr[i][j][k]->GetEntries()<200)continue;
     
     histo_delta_shift_spr[i][j][k]->Fit("fit_gause2","M","", -20, 20);
     double trck_mean_fit = fit_gause2->GetParameter(1);
     double trck_sigma_fit = fit_gause2->GetParameter(2);
     //shift_array_detailed[i][j[k]=fit_gause1->GetParameter(1);
     if(histo_delta_shift_spr[i][j][k]->GetEntries()<200){
     trck_mean_fit = -1000;
     trck_sigma_fit = 0;
     
     }
     if(trck_mean_fit != -1000){
     histo_pos_xy_TrkMeanShift_average->Fill(i,j,trck_mean_fit);
     histo_pos_xy_TrkMeanShift_occu_average->Fill(i,j);
     }
     
     histo_pos_xy_TrkSigmaShift_average->Fill(i,j,trck_sigma_fit);
     histo_pos_xy_TrkSigmaShift_occu_average->Fill(i,j);
     
     histo_pos_xy_TrkMeanShift[k]->Fill(i,j,trck_mean_fit);
     histo_pos_xy_TrkMeanShift_occu[k]->Fill(i,j);
     
     histo_pos_xy_TrkSigmaShift[k]->Fill(i,j,trck_sigma_fit);
     histo_pos_xy_TrkSigmaShift_occu[k]->Fill(i,j);
     
     //                cc4->cd();
     //                cc4->Update();
     //                histo_delta_shift_spr[i][j][k]->Draw();
     //                cc4->Update();
     //                TLine *lineMeanSHift= new TLine(0,0,0,1000);
     //                lineMeanSHift->SetX1(fit_gause2->GetParameter(1));
     //                lineMeanSHift->SetX2(fit_gause2->GetParameter(1));
     //                lineMeanSHift->SetY1(gPad->GetUymin());
     //                lineMeanSHift->SetY2(gPad->GetUymax());
     //                lineMeanSHift->SetLineColor(kRed);
     //                lineMeanSHift->Draw();
     //                cc4->Update();
     //                cc4->WaitPrimitive();
     }
     }
     }
     int v11(5) ,v22(5);
     TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
     for(Int_t k=0; k<16; k++){
     v11 = v11+1;
     v22 = v11+1;
     ratio_mean[k] = (TH2F*)histo_pos_xy_TrkMeanShift[k]->Clone();
     ratio_mean[k]->GetXaxis()->SetTitle(" ");
     ratio_mean[k]->GetYaxis()->SetTitle(" ");
     ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ SPR [%d - %d] mrad", v11,v22));
     ratio_mean[k]->Divide(histo_pos_xy_TrkMeanShift_occu[k]);
     cc2->cd();
     cc2->Update();
     ratio_mean[k]->SetMinimum(-5);
     ratio_mean[k]->SetMaximum(5);
     ratio_mean[k]->Draw("COLZ");
     cc2->Update();
     cc2->WaitPrimitive();
     }
     
     
     int v1(5) ,v2(5);
     TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
     for(Int_t k=0; k<16; k++){
     v1 = v1+1;
     v2 = v1+1;
     ratio_sigma[k] = (TH2F*)histo_pos_xy_TrkSigmaShift[k]->Clone();
     ratio_sigma[k]->GetXaxis()->SetTitle(" ");
     ratio_sigma[k]->GetYaxis()->SetTitle(" ");
     ratio_sigma[k]->SetTitle(Form("Sigma of  #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ SPR [%d - %d] mrad",v1, v2));
     ratio_sigma[k]->Divide(histo_pos_xy_TrkSigmaShift_occu[k]);
     cc3->cd();
     cc3->Update();
     ratio_sigma[k]->Draw("COLZ");
     cc3->Update();
     cc3->WaitPrimitive();
     }
     
     glx_canvasAdd("r_ratio_mean",800,400);
     ratio_mean_average = (TH2F*)histo_pos_xy_TrkMeanShift_average->Clone();
     ratio_mean_average->GetXaxis()->SetTitle(" ");
     ratio_mean_average->GetYaxis()->SetTitle(" ");
     ratio_mean_average->SetTitle("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
     ratio_mean_average->Divide(histo_pos_xy_TrkMeanShift_occu_average);
     ratio_mean_average->Draw("COLZ");
     
     //glx_canvasAdd("r_histo_pos_xy_TrkMeanShift",800,400);
     //histo_pos_xy_TrkMeanShift_average->Draw("colz");
     
     glx_canvasAdd("r_ratio_sigma",800,400);
     ratio_sigma_average = (TH2F*)histo_pos_xy_TrkSigmaShift_average->Clone();
     ratio_sigma_average->GetXaxis()->SetTitle(" ");
     ratio_sigma_average->GetYaxis()->SetTitle(" ");
     ratio_sigma_average->SetTitle("Sigma of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
     ratio_sigma_average->Divide(histo_pos_xy_TrkSigmaShift_occu_average);
     ratio_sigma_average->Draw("COLZ");
     
     ///////////////////////////////////////////////////////////////////
     
     
     */
    
    
    TF1 *fit_gause2 = new TF1("fit_gause2","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
    fit_gause2->SetLineColor(kBlack);
    fit_gause2->SetParameters(100,9,2);
    fit_gause2->SetParNames("p0","mean ","sigma");
    fit_gause2->SetParLimits(0,0.1,1E6);
    fit_gause2->SetParLimits(1,-10,10);
    fit_gause2->SetParLimits(2,1,10);
    
    if(true){
        ////////////////////////////////
        // Cherenkov Shift Yield Study//
        ////////////////////////////////
        
        // Reset Histograms
        histo_pos_xy_TrkSigmaShift_average->Reset();
        histo_pos_xy_TrkSigmaShift_occu_average->Reset();
        ratio_mean_average->Reset();
        ratio_sigma_average->Reset();
        
        for(Int_t k=0; k<60; k++){
            histo_pos_xy_TrkMeanShift[k]->Reset();
            histo_pos_xy_TrkMeanShift_occu[k]->Reset();
            histo_pos_xy_TrkSigmaShift[k]->Reset();
            histo_pos_xy_TrkSigmaShift_occu[k]->Reset();
            ratio_mean[k]->Reset();
            ratio_sigma[k]->Reset();
        }
        
        TCanvas *cc5 = new TCanvas("cc5","cc5",800,500);
        
        //double shift_array_detailed[42][26][22]={0};
        for(Int_t i=0; i<42; i++) {
            for(Int_t j=0; j<26; j++){
                for(Int_t k=0; k<60; k++){
                    if(histo_delta_shift_nph[i][j][k]->GetEntries()<500)continue;
                    histo_delta_shift_nph[i][j][k]->Fit("fit_gause2","M","", -20, 20);
                    double trck_mean_fit = fit_gause2->GetParameter(1);
                    double trck_sigma_fit = fit_gause2->GetParameter(2);
                    //shift_array_detailed[i][j[k]=fit_gause1->GetParameter(1);
                    if(histo_delta_shift_nph[i][j][k]->GetEntries()<200){
                        trck_mean_fit = -1000;
                        trck_sigma_fit = 0;
                        
                    }
                    if(trck_mean_fit != -1000){
                        histo_pos_xy_TrkMeanShift_average->Fill(i,j,trck_mean_fit);
                        histo_pos_xy_TrkMeanShift_occu_average->Fill(i,j);
                    }
                    
                    histo_pos_xy_TrkSigmaShift_average->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu_average->Fill(i,j);
                    
                    histo_pos_xy_TrkMeanShift[k]->Fill(i,j,trck_mean_fit);
                    histo_pos_xy_TrkMeanShift_occu[k]->Fill(i,j);
                    
                    histo_pos_xy_TrkSigmaShift[k]->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu[k]->Fill(i,j);
                    
                                    cc5->cd();
                                    cc5->Update();
                                    histo_delta_shift_nph[i][j][k]->Draw();
                                    cc5->Update();
                                    TLine *lineMeanSHift= new TLine(0,0,0,1000);
                                    lineMeanSHift->SetX1(fit_gause2->GetParameter(1));
                                    lineMeanSHift->SetX2(fit_gause2->GetParameter(1));
                                    lineMeanSHift->SetY1(gPad->GetUymin());
                                    lineMeanSHift->SetY2(gPad->GetUymax());
                                    lineMeanSHift->SetLineColor(kRed);
                                    lineMeanSHift->Draw();
                                    cc5->Update();
                                    cc5->WaitPrimitive();
                    
                    
                }
            }
        }
        
        
        if(true){
            //    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
            int counter_MeanShiftYield =0;
            for(Int_t k=1; k<60; k++){
                ratio_mean[k] = (TH2F*)histo_pos_xy_TrkMeanShift[k]->Clone();
                //ratio_mean[k]->GetXaxis()->SetTitle(" ");
                //ratio_mean[k]->GetYaxis()->SetTitle(" ");
                ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ photon yield %d", k));
                ratio_mean[k]->Divide(histo_pos_xy_TrkMeanShift_occu[k]);
                //        cc2->cd();
                //        cc2->Update();
                //        ratio_mean[k]->SetMinimum(-5);
                //        ratio_mean[k]->SetMaximum(5);
                //        ratio_mean[k]->Draw("COLZ");
                //        cc2->Update();
                //        cc2->WaitPrimitive();
                
                TString num_string=Form("_%d",counter_MeanShiftYield);
                glx_canvasAdd("r_MeanShiftYield"+num_string,800,400);
                ratio_mean[k]->SetMinimum(-5);
                ratio_mean[k]->SetMaximum(5);
                ratio_mean[k]->Draw("COLZ");
                glx_canvasGet("r_MeanShiftYield"+num_string)->Update();
                ++counter_MeanShiftYield;
            }
            
            //    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
            int counter_SigmaShiftYield =0;
            for(Int_t k=1; k<60; k++){
                ratio_sigma[k] = (TH2F*)histo_pos_xy_TrkSigmaShift[k]->Clone();
                //ratio_sigma[k]->GetXaxis()->SetTitle(" ");
                //ratio_sigma[k]->GetYaxis()->SetTitle(" ");
                ratio_sigma[k]->SetTitle(Form("Sigma of  #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ photon yield %d",k));
                ratio_sigma[k]->Divide(histo_pos_xy_TrkSigmaShift_occu[k]);
                //        cc3->cd();
                //        cc3->Update();
                //        ratio_sigma[k]->Draw("COLZ");
                //        cc3->Update();
                //        cc3->WaitPrimitive();
                TString num_string=Form("_%d",counter_SigmaShiftYield);
                glx_canvasAdd("r_SigmaShiftYield"+num_string,800,400);
                ratio_sigma[k]->SetMinimum(0);
                ratio_sigma[k]->SetMaximum(5);
                ratio_sigma[k]->Draw("COLZ");
                glx_canvasGet("r_SigmaShiftYield"+num_string)->Update();
                ++counter_SigmaShiftYield;
            }
            
            glx_canvasSave(2,0);
            glx_canvasDel("*");
        }
    }
    
    
    if(false){
        //////////////////////////////
        // Cherenkov Shift SPR Study//
        //////////////////////////////
        cout<<"No proplem b4 Reset"<<endl;
        // Reset Histograms
        histo_pos_xy_TrkSigmaShift_average->Reset();
        histo_pos_xy_TrkSigmaShift_occu_average->Reset();
        ratio_mean_average->Reset();
        ratio_sigma_average->Reset();
        
        for(Int_t k=0; k<60; k++){
            histo_pos_xy_TrkMeanShift[k]->Reset();
            histo_pos_xy_TrkMeanShift_occu[k]->Reset();
            histo_pos_xy_TrkSigmaShift[k]->Reset();
            histo_pos_xy_TrkSigmaShift_occu[k]->Reset();
            ratio_mean[k]->Reset();
            ratio_sigma[k]->Reset();
        }
        
        cout<<"No proplem after Reset"<<endl;
        
        //TCanvas *cc5 = new TCanvas("cc5","cc5",800,500);
        //    TF1 *fit_gause2 = new TF1("fit_gause2","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-20,20);
        //    fit_gause2->SetLineColor(kBlack);
        //    fit_gause2->SetParameters(100,9,2);
        //    fit_gause2->SetParNames("p0","mean ","sigma");
        //    fit_gause2->SetParLimits(0,0.1,1E6);
        //    fit_gause2->SetParLimits(1,-10,10);
        //    fit_gause2->SetParLimits(2,1,10);
        //double shift_array_detailed[42][26][22]={0};
        for(Int_t i=0; i<42; i++) {
            for(Int_t j=0; j<26; j++){
                for(Int_t k=0; k<16; k++){
                    //if(histo_delta_shift_spr[i][j][k]->GetEntries()<500)continue;
                    histo_delta_shift_spr[i][j][k]->Fit("fit_gause2","M","", -20, 20);
                    double trck_mean_fit = fit_gause2->GetParameter(1);
                    double trck_sigma_fit = fit_gause2->GetParameter(2);
                    //shift_array_detailed[i][j[k]=fit_gause1->GetParameter(1);
                    if(histo_delta_shift_spr[i][j][k]->GetEntries()<200){
                        trck_mean_fit = -1000;
                        trck_sigma_fit = 0;
                        
                    }
                    if(trck_mean_fit != -1000){
                        histo_pos_xy_TrkMeanShift_average->Fill(i,j,trck_mean_fit);
                        histo_pos_xy_TrkMeanShift_occu_average->Fill(i,j);
                    }
                    
                    histo_pos_xy_TrkSigmaShift_average->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu_average->Fill(i,j);
                    
                    histo_pos_xy_TrkMeanShift[k]->Fill(i,j,trck_mean_fit);
                    histo_pos_xy_TrkMeanShift_occu[k]->Fill(i,j);
                    
                    histo_pos_xy_TrkSigmaShift[k]->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu[k]->Fill(i,j);
                    
                    //                cc5->cd();
                    //                cc5->Update();
                    //                histo_delta_shift_spr[i][j][k]->Draw();
                    //                cc5->Update();
                    //                TLine *lineMeanSHift= new TLine(0,0,0,1000);
                    //                lineMeanSHift->SetX1(fit_gause2->GetParameter(1));
                    //                lineMeanSHift->SetX2(fit_gause2->GetParameter(1));
                    //                lineMeanSHift->SetY1(gPad->GetUymin());
                    //                lineMeanSHift->SetY2(gPad->GetUymax());
                    //                lineMeanSHift->SetLineColor(kRed);
                    //                lineMeanSHift->Draw();
                    //                cc5->Update();
                    //                cc5->WaitPrimitive();
                    
                    
                }
            }
        }
        
        
        if(true){
            //    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
            int v11(5) ,v22(5);
            int counter_MeanShiftSPR =0;
            for(Int_t k=1; k<16; k++){
                v11 = v11+1;
                v22 = v11+1;
                ratio_mean[k] = (TH2F*)histo_pos_xy_TrkMeanShift[k]->Clone();
                //ratio_mean[k]->GetXaxis()->SetTitle(" ");
                //ratio_mean[k]->GetYaxis()->SetTitle(" ");
                ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ SPR [%d - %d] mrad", v11,v22));
                ratio_mean[k]->Divide(histo_pos_xy_TrkMeanShift_occu[k]);
                //        cc2->cd();
                //        cc2->Update();
                //        ratio_mean[k]->SetMinimum(-5);
                //        ratio_mean[k]->SetMaximum(5);
                //        ratio_mean[k]->Draw("COLZ");
                //        cc2->Update();
                //        cc2->WaitPrimitive();
                
                TString num_string=Form("_%d",counter_MeanShiftSPR);
                glx_canvasAdd("r_MeanShiftSPR"+num_string,800,400);
                ratio_mean[k]->SetMinimum(-5);
                ratio_mean[k]->SetMaximum(5);
                ratio_mean[k]->Draw("COLZ");
                glx_canvasGet("r_MeanShiftSPR"+num_string)->Update();
                ++counter_MeanShiftSPR;
            }
            
            //    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
            int counter_SigmaShiftSPR =0;
            int v1(5) ,v2(5);
            for(Int_t k=1; k<16; k++){
                v1 = v1+1;
                v2 = v1+1;
                ratio_sigma[k] = (TH2F*)histo_pos_xy_TrkSigmaShift[k]->Clone();
                //ratio_sigma[k]->GetXaxis()->SetTitle(" ");
                //ratio_sigma[k]->GetYaxis()->SetTitle(" ");
                ratio_sigma[k]->SetTitle(Form("Sigma of  #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ SPR [%d - %d] mrad",v1, v2));
                ratio_sigma[k]->Divide(histo_pos_xy_TrkSigmaShift_occu[k]);
                //        cc3->cd();
                //        cc3->Update();
                //        ratio_sigma[k]->Draw("COLZ");
                //        cc3->Update();
                //        cc3->WaitPrimitive();
                TString num_string=Form("_%d",counter_SigmaShiftSPR);
                glx_canvasAdd("r_SigmaShiftSPR"+num_string,800,400);
                ratio_sigma[k]->SetMinimum(0);
                ratio_sigma[k]->SetMaximum(5);
                ratio_sigma[k]->Draw("COLZ");
                glx_canvasGet("r_SigmaShiftSPR"+num_string)->Update();
                ++counter_SigmaShiftSPR;
            }
            
            glx_canvasSave(2,0);
            glx_canvasDel("*");
        }
    }
    
    
    if(false){
        ///////////////////////////////////
        // Cherenkov Shift Momentum Study//
        ///////////////////////////////////
        cout<<"No proplem b4 Reset"<<endl;
        // Reset Histograms
        histo_pos_xy_TrkSigmaShift_average->Reset();
        histo_pos_xy_TrkSigmaShift_occu_average->Reset();
        ratio_mean_average->Reset();
        ratio_sigma_average->Reset();
        
        for(Int_t k=0; k<60; k++){
            histo_pos_xy_TrkMeanShift[k]->Reset();
            histo_pos_xy_TrkMeanShift_occu[k]->Reset();
            histo_pos_xy_TrkSigmaShift[k]->Reset();
            histo_pos_xy_TrkSigmaShift_occu[k]->Reset();
            ratio_mean[k]->Reset();
            ratio_sigma[k]->Reset();
        }
        
        cout<<"No proplem after Reset"<<endl;
        
        //TCanvas *cc5 = new TCanvas("cc5","cc5",800,500);
        for(Int_t i=0; i<42; i++) {
            for(Int_t j=0; j<26; j++){
                for(Int_t k=0; k<16; k++){
                    //if(histo_delta_shift_mom[i][j][k]->GetEntries()<500)continue;
                    histo_delta_shift_mom[i][j][k]->Fit("fit_gause2","M","", -20, 20);
                    double trck_mean_fit = fit_gause2->GetParameter(1);
                    double trck_sigma_fit = fit_gause2->GetParameter(2);
                    //shift_array_detailed[i][j[k]=fit_gause1->GetParameter(1);
                    if(histo_delta_shift_mom[i][j][k]->GetEntries()<200){
                        trck_mean_fit = -1000;
                        trck_sigma_fit = 0;
                        
                    }
                    if(trck_mean_fit != -1000){
                        histo_pos_xy_TrkMeanShift_average->Fill(i,j,trck_mean_fit);
                        histo_pos_xy_TrkMeanShift_occu_average->Fill(i,j);
                    }
                    
                    histo_pos_xy_TrkSigmaShift_average->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu_average->Fill(i,j);
                    
                    histo_pos_xy_TrkMeanShift[k]->Fill(i,j,trck_mean_fit);
                    histo_pos_xy_TrkMeanShift_occu[k]->Fill(i,j);
                    
                    histo_pos_xy_TrkSigmaShift[k]->Fill(i,j,trck_sigma_fit);
                    histo_pos_xy_TrkSigmaShift_occu[k]->Fill(i,j);
                    
                    //                cc5->cd();
                    //                cc5->Update();
                    //                histo_delta_shift_mom[i][j][k]->Draw();
                    //                cc5->Update();
                    //                TLine *lineMeanSHift= new TLine(0,0,0,1000);
                    //                lineMeanSHift->SetX1(fit_gause2->GetParameter(1));
                    //                lineMeanSHift->SetX2(fit_gause2->GetParameter(1));
                    //                lineMeanSHift->SetY1(gPad->GetUymin());
                    //                lineMeanSHift->SetY2(gPad->GetUymax());
                    //                lineMeanSHift->SetLineColor(kRed);
                    //                lineMeanSHift->Draw();
                    //                cc5->Update();
                    //                cc5->WaitPrimitive();
                    
                    
                }
            }
        }
        
        
        if(true){
            //    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
            double v11(0) ,v22(0);
            int counter_MeanShiftMom =0;
            for(Int_t k=0; k<16; k++){
                v11 = v11+0.5;
                v22 = v11+0.5;
                ratio_mean[k] = (TH2F*)histo_pos_xy_TrkMeanShift[k]->Clone();
                //ratio_mean[k]->GetXaxis()->SetTitle(" ");
                //ratio_mean[k]->GetYaxis()->SetTitle(" ");
                ratio_mean[k]->SetTitle(Form("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ Momentum [%1.1f - %1.1f] GeV/c", v11,v22));
                ratio_mean[k]->Divide(histo_pos_xy_TrkMeanShift_occu[k]);
                //        cc2->cd();
                //        cc2->Update();
                //        ratio_mean[k]->SetMinimum(-5);
                //        ratio_mean[k]->SetMaximum(5);
                //        ratio_mean[k]->Draw("COLZ");
                //        cc2->Update();
                //        cc2->WaitPrimitive();
                
                TString num_string=Form("_%d",counter_MeanShiftMom);
                glx_canvasAdd("r_MeanShiftMom"+num_string,800,400);
                ratio_mean[k]->SetMinimum(-5);
                ratio_mean[k]->SetMaximum(5);
                ratio_mean[k]->Draw("COLZ");
                glx_canvasGet("r_MeanShiftMom"+num_string)->Update();
                ++counter_MeanShiftMom;
            }
            
            //    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
            int counter_SigmaShiftMom =0;
            double v1(0) ,v2(0);
            for(Int_t k=0; k<16; k++){
                v1 = v1+0.5;
                v2 = v1+0.5;
                ratio_sigma[k] = (TH2F*)histo_pos_xy_TrkSigmaShift[k]->Clone();
                //ratio_sigma[k]->GetXaxis()->SetTitle(" ");
                //ratio_sigma[k]->GetYaxis()->SetTitle(" ");
                ratio_sigma[k]->SetTitle(Form("Sigma of  #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution @ Momentum [%1.1f - %1.1f] GeV/c",v1, v2));
                ratio_sigma[k]->Divide(histo_pos_xy_TrkSigmaShift_occu[k]);
                //        cc3->cd();
                //        cc3->Update();
                //        ratio_sigma[k]->Draw("COLZ");
                //        cc3->Update();
                //        cc3->WaitPrimitive();
                TString num_string=Form("_%d",counter_SigmaShiftMom);
                glx_canvasAdd("r_SigmaShiftMom"+num_string,800,400);
                ratio_sigma[k]->SetMinimum(0);
                ratio_sigma[k]->SetMaximum(5);
                ratio_sigma[k]->Draw("COLZ");
                glx_canvasGet("r_SigmaShiftMom"+num_string)->Update();
                ++counter_SigmaShiftMom;
            }
            
            glx_canvasSave(2,0);
            glx_canvasDel("*");
        }
    }
    
    
    
    return 0;
    
    glx_canvasAdd("r_ratio_mean",800,400);
    ratio_mean_average = (TH2F*)histo_pos_xy_TrkMeanShift_average->Clone();
    ratio_mean_average->GetXaxis()->SetTitle(" ");
    ratio_mean_average->GetYaxis()->SetTitle(" ");
    ratio_mean_average->SetTitle("Mean Value of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
    ratio_mean_average->Divide(histo_pos_xy_TrkMeanShift_occu_average);
    ratio_mean_average->SetMinimum(-5);
    ratio_mean_average->SetMaximum(5);
    ratio_mean_average->Draw("COLZ");
    
    //glx_canvasAdd("r_histo_pos_xy_TrkMeanShift",800,400);
    //histo_pos_xy_TrkMeanShift_average->Draw("colz");
    
    glx_canvasAdd("r_ratio_sigma",800,400);
    ratio_sigma_average = (TH2F*)histo_pos_xy_TrkSigmaShift_average->Clone();
    ratio_sigma_average->GetXaxis()->SetTitle(" ");
    ratio_sigma_average->GetYaxis()->SetTitle(" ");
    ratio_sigma_average->SetTitle("Sigma of #theta_{c}^{Measured} - #theta_{c}^{Expected} distribution ");
    ratio_sigma_average->Divide(histo_pos_xy_TrkSigmaShift_occu_average);
    ratio_sigma_average->SetMinimum(0);
    ratio_sigma_average->SetMaximum(5);
    ratio_sigma_average->Draw("COLZ");
    
    ///////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    //    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
    //    for(Int_t l=0; l<50; l++){
    //        ratio_yield[l] = (TH2F*)histo_pos_xy_TrkMeanShift_yield[l]->Clone();
    //        ratio_yield[l]->GetXaxis()->SetTitle(" ");
    //        ratio_yield[l]->GetYaxis()->SetTitle(" ");
    //        //ratio_yield[l]->SetTitle("histo_pos_xy_TrkMeanShift_yield/histo_pos_xy_TrkMeanShift_occu_average");
    //        ratio_yield[l]->Divide(histo_pos_xy_TrkMeanShift_occu_yield[l]);
    //        ratio_yield[l]->Draw("COLZ");
    //
    //    }
    
    
    
    cout<<"##### commint 5 "<<endl;
    
    /////////////////////////////////////
    
    //    if(true){
    //
    //        // yield per bar per mom
    //        TCanvas *canvas4 = new TCanvas("canvas4","canvas4",800,500);
    //        for(int i=0;i<nbar;i++){
    //            for(int j=0;j<nbin_mom;j++){
    //                //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
    //                if(histo_track_yield_bar_mom[i][j]->GetEntries()<1)continue;
    //                canvas4->cd();
    //                canvas4->Update();
    //                histo_track_yield_bar_mom[i][j]->Draw();
    //                canvas4->Update();
    //                canvas4->WaitPrimitive();
    //            }
    //
    //        }
    //        delete canvas4;
    //
    //
    //        // diff array per bar per mom  used
    //        TCanvas *canvas1 = new TCanvas("canvas1","canvas1",800,500);
    //        TF1 *fit_gause2 = new TF1("fit_gause2","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    //        fit_gause2->SetLineColor(kBlack);
    //        fit_gause2->SetParameters(100,9,2);
    //        fit_gause2->SetParNames("p0","mean ","sigma");
    //        fit_gause2->SetParLimits(0,0.1,1E6);
    //        fit_gause2->SetParLimits(1,0.818,0.834);
    //        fit_gause2->SetParLimits(2,0.0001,0.005);
    //
    //        //double diff_array_detailed[26][10]={0.82608};
    //        for(int i=0;i<26;i++){
    //            for(int j=0;j<10;j++){
    //                //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
    //                if(histo_track_mean_bar_mom[i][j]->GetEntries()<100)continue;
    //                histo_track_mean_bar_mom[i][j]->Fit("fit_gause2","M","", 0, 30);
    //                diff_array_detailed[i][j]=fit_gause2->GetParameter(1);
    //                canvas1->cd();
    //                canvas1->Update();
    //                histo_track_mean_bar_mom[i][j]->Draw();
    //                canvas1->Update();
    //                TLine *lineEXm1= new TLine(0,0,0,1000);
    //                lineEXm1->SetX1(diff_array_detailed[i][j]);
    //                lineEXm1->SetX2(diff_array_detailed[i][j]);
    //                lineEXm1->SetY1(gPad->GetUymin());
    //                lineEXm1->SetY2(gPad->GetUymax());
    //                lineEXm1->SetLineColor(kRed);
    //                lineEXm1->Draw();
    //                canvas1->Update();
    //                canvas1->WaitPrimitive();
    //                delete histo_track_mean_bar_mom[i][j];
    //            }
    //        }
    //        delete canvas1;
    //
    //        // printing diff array per bar per mom
    //        //    cout<<"diff_array_detailed[i][j]={";
    //        //    for(int i=0;i<26;i++){
    //        //        for(int j=0;j<10;j++){
    //        //            cout<<","<<fit_gause2->GetParameter(1);
    //        //        }
    //        //    }
    //        //    cout<<"};"<<endl;
    //
    //    }
    
    //    if(true){
    //        cout<<"##### Histograms "<<endl;
    //        //        glx_drawDigi("m,p,v\n",0);
    //        //        glx_canvasAdd("r_cherenkov",800,400);
    //        //        histo_cherenkov->Draw();
    //        //
    //        //        glx_canvasAdd("r_tdiff",800,400);
    //        //        histo_tdiff->Draw();
    //
    //        glx_canvasAdd("r_pos_shift",800,400);
    //        //histo_pos_xy_shift->SetMinimum(5);
    //        //histo_pos_xy_shift->SetMaximum(10);
    //        histo_pos_xy_shift->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_shift->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_shiftEx",800,400);
    //        //histo_pos_xy_shiftEx->SetMinimum(5);
    //        //histo_pos_xy_shiftEx->SetMaximum(10);
    //        histo_pos_xy_shiftEx->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_shiftEx->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_shiftEx_positive",800,400);
    //        //histo_pos_xy_shiftEx_positive->SetMinimum(5);
    //        //histo_pos_xy_shiftEx_positive->SetMaximum(10);
    //        histo_pos_xy_shiftEx_positive->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_shiftEx_positive->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_shiftEx_negative",800,400);
    //        //histo_pos_xy_shiftEx_negative->SetMinimum(5);
    //        //histo_pos_xy_shiftEx_negative->SetMaximum(10);
    //        histo_pos_xy_shiftEx_negative->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_shiftEx_negative->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_shiftEx_occu_positive",800,400);
    //        //histo_pos_xy_occupancy_postiveShift->SetMinimum(5);
    //        //histo_pos_xy_occupancy_postiveShift->SetMaximum(10);
    //        histo_pos_xy_occupancy_postiveShift->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_occupancy_postiveShift->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_shiftEx_occu_negative",800,400);
    //        //histo_pos_xy_occupancy_negativeShift->SetMinimum(5);
    //        //histo_pos_xy_occupancy_negativeShift->SetMaximum(10);
    //        histo_pos_xy_occupancy_negativeShift->GetYaxis()->SetRangeUser(-100,0);
    //        histo_pos_xy_occupancy_negativeShift->Draw("colz");
    //
    //
    //        glx_canvasAdd("r_rho_mass",800,400);
    //        hist_ev_rho_mass->Draw();
    //        glx_canvasGet("r_rho_mass")->Update();
    //        TLine *lin_rho_mass_max= new TLine(0,0,0,1000);
    //        lin_rho_mass_max->SetX1(mass_rho_max);
    //        lin_rho_mass_max->SetX2(mass_rho_max);
    //        lin_rho_mass_max->SetY1(gPad->GetUymin());
    //        lin_rho_mass_max->SetY2(gPad->GetUymax());
    //        lin_rho_mass_max->SetLineColor(kBlack);
    //        lin_rho_mass_max->Draw();
    //        TLine *lin_rho_mass_min= new TLine(0,0,0,1000);
    //        lin_rho_mass_min->SetX1(mass_rho_mini);
    //        lin_rho_mass_min->SetX2(mass_rho_mini);
    //        lin_rho_mass_min->SetY1(gPad->GetUymin());
    //        lin_rho_mass_min->SetY2(gPad->GetUymax());
    //        lin_rho_mass_min->SetLineColor(kBlack);
    //        lin_rho_mass_min->Draw();
    //        glx_canvasGet("r_rho_mass")->Update();
    //
    //        glx_canvasAdd("r_phi_mass",800,400);
    //        hist_ev_phi_mass->Draw();
    //        glx_canvasGet("r_phi_mass")->Update();
    //        TLine *lin_phi_mass_max= new TLine(0,0,0,1000);
    //        lin_phi_mass_max->SetX1(mass_phi_max);
    //        lin_phi_mass_max->SetX2(mass_phi_max);
    //        lin_phi_mass_max->SetY1(gPad->GetUymin());
    //        lin_phi_mass_max->SetY2(gPad->GetUymax());
    //        lin_phi_mass_max->SetLineColor(kBlack);
    //        lin_phi_mass_max->Draw();
    //        TLine *lin_phi_mass_min= new TLine(0,0,0,1000);
    //        lin_phi_mass_min->SetX1(mass_phi_mini);
    //        lin_phi_mass_min->SetX2(mass_phi_mini);
    //        lin_phi_mass_min->SetY1(gPad->GetUymin());
    //        lin_phi_mass_min->SetY2(gPad->GetUymax());
    //        lin_phi_mass_min->SetLineColor(kBlack);
    //        lin_phi_mass_min->Draw();
    //        glx_canvasGet("r_phi_mass")->Update();
    //
    //        glx_canvasAdd("r_missing_mass_phi",800,400);
    //        hist_ev_missing_mass_phi->Draw();
    //        glx_canvasGet("r_missing_mass_phi")->Update();
    //        TLine *lin_phi_miss_mass_max= new TLine(0,0,0,1000);
    //        lin_phi_miss_mass_max->SetX1(miss_mass_phi_max);
    //        lin_phi_miss_mass_max->SetX2(miss_mass_phi_max);
    //        lin_phi_miss_mass_max->SetY1(gPad->GetUymin());
    //        lin_phi_miss_mass_max->SetY2(gPad->GetUymax());
    //        lin_phi_miss_mass_max->SetLineColor(kBlack);
    //        lin_phi_miss_mass_max->Draw();
    //        TLine *lin_phi_miss_mass_min= new TLine(0,0,0,1000);
    //        lin_phi_miss_mass_min->SetX1(miss_mass_phi_mini);
    //        lin_phi_miss_mass_min->SetX2(miss_mass_phi_mini);
    //        lin_phi_miss_mass_min->SetY1(gPad->GetUymin());
    //        lin_phi_miss_mass_min->SetY2(gPad->GetUymax());
    //        lin_phi_miss_mass_min->SetLineColor(kBlack);
    //        lin_phi_miss_mass_min->Draw();
    //        glx_canvasGet("r_missing_mass_phi")->Update();
    //
    //        glx_canvasAdd("r_missing_mass_rho",800,400);
    //        hist_ev_missing_mass_rho->Draw();
    //        glx_canvasGet("r_missing_mass_rho")->Update();
    //        TLine *lin_rho_miss_mass_max= new TLine(0,0,0,1000);
    //        lin_rho_miss_mass_max->SetX1(miss_mass_rho_max);
    //        lin_rho_miss_mass_max->SetX2(miss_mass_rho_max);
    //        lin_rho_miss_mass_max->SetY1(gPad->GetUymin());
    //        lin_rho_miss_mass_max->SetY2(gPad->GetUymax());
    //        lin_rho_miss_mass_max->SetLineColor(kBlack);
    //        lin_rho_miss_mass_max->Draw();
    //        TLine *lin_rho_miss_mass_min= new TLine(0,0,0,1000);
    //        lin_rho_miss_mass_min->SetX1(miss_mass_rho_mini);
    //        lin_rho_miss_mass_min->SetX2(miss_mass_rho_mini);
    //        lin_rho_miss_mass_min->SetY1(gPad->GetUymin());
    //        lin_rho_miss_mass_min->SetY2(gPad->GetUymax());
    //        lin_rho_miss_mass_min->SetLineColor(kBlack);
    //        lin_rho_miss_mass_min->Draw();
    //        glx_canvasGet("r_missing_mass_rho")->Update();
    //
    //        glx_canvasAdd("r_chi_phi",800,400);
    //        hist_ev_chi_phi->Draw();
    //        glx_canvasGet("r_chi_phi")->Update();
    //        TLine *lin_chi_phi_max= new TLine(0,0,0,1000);
    //        lin_chi_phi_max->SetX1(chisq_phi_max);
    //        lin_chi_phi_max->SetX2(chisq_phi_max);
    //        lin_chi_phi_max->SetY1(gPad->GetUymin());
    //        lin_chi_phi_max->SetY2(gPad->GetUymax());
    //        lin_chi_phi_max->SetLineColor(kBlack);
    //        lin_chi_phi_max->Draw();
    //        TLine *lin_chi_phi_min= new TLine(0,0,0,1000);
    //        lin_chi_phi_min->SetX1(chisq_phi_mini);
    //        lin_chi_phi_min->SetX2(chisq_phi_mini);
    //        lin_chi_phi_min->SetY1(gPad->GetUymin());
    //        lin_chi_phi_min->SetY2(gPad->GetUymax());
    //        lin_chi_phi_min->SetLineColor(kBlack);
    //        lin_chi_phi_min->Draw();
    //        glx_canvasGet("r_chi_phi")->Update();
    //
    //        glx_canvasAdd("r_chi_rho",800,400);
    //        hist_ev_chi_rho->Draw();
    //        glx_canvasGet("r_chi_rho")->Update();
    //        TLine *lin_chi_rho_max= new TLine(0,0,0,1000);
    //        lin_chi_rho_max->SetX1(chisq_rho_max);
    //        lin_chi_rho_max->SetX2(chisq_rho_max);
    //        lin_chi_rho_max->SetY1(gPad->GetUymin());
    //        lin_chi_rho_max->SetY2(gPad->GetUymax());
    //        lin_chi_rho_max->SetLineColor(kBlack);
    //        lin_chi_rho_max->Draw();
    //        TLine *lin_chi_rho_min= new TLine(0,0,0,1000);
    //        lin_chi_rho_min->SetX1(chisq_rho_mini);
    //        lin_chi_rho_min->SetX2(chisq_rho_mini);
    //        lin_chi_rho_min->SetY1(gPad->GetUymin());
    //        lin_chi_rho_min->SetY2(gPad->GetUymax());
    //        lin_chi_rho_min->SetLineColor(kBlack);
    //        lin_chi_rho_min->Draw();
    //        glx_canvasGet("r_chi_rho")->Update();
    //    }
    
    cout<<"##### commint 6 "<<endl;
    /*
     ////////////////////////////////
     //////// Cherenkove angle //////
     ////////////////////////////////
     
     
     glx_canvasAdd("hAngle",800,400);
     
     //scal
     if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
     if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
     
     for(int i=0; i<5; i++){
     if(hAngle[i]->GetEntries()<20) continue;
     
     int nfound = spect->Search(hAngle[i],1,"goff",0.9);
     if(nfound>0) cherenkovreco[i] = spect->GetPositionX()[0];
     else cherenkovreco[i] =  hAngle[i]->GetXaxis()->GetBinCenter(hAngle[i]->GetMaximumBin());
     if(cherenkovreco[i]>0.85) cherenkovreco[i]=0.82;
     
     if(i==2)  fit->SetLineColor(kBlue);
     if(i==3)  fit->SetLineColor(kRed);
     fit->SetParameters(100,cherenkovreco[i],0.010,10);
     fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
     fit->SetParLimits(0,0.1,1E6);
     fit->SetParLimits(1,cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);
     fit->SetParLimits(2,0.005,0.030); // width
     hAngle[i]->Fit("fgaus","I","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
     hAngle[i]->Fit("fgaus","M","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
     
     cherenkovreco[i] = fit->GetParameter(1);
     spr[i] = fit->GetParameter(2);
     }
     
     gStyle->SetOptTitle(0);
     gStyle->SetOptStat(0);
     gStyle->SetOptFit(0);
     
     hAngle[2]->GetXaxis()->SetRangeUser(0.7,0.9);
     hAngle[2]->GetYaxis()->SetRangeUser(0,1.2);
     hAngle[2]->Draw();
     hAngle[3]->Draw("same");
     // fAngle[3]->Draw("same");
     // fAngle[2]->Draw("same");
     
     
     TLine *line = new TLine(0,0,0,1000);
     line->SetX1(mAngle[3]);
     line->SetX2(mAngle[3]);
     line->SetY1(0);
     line->SetY2(1.2);
     line->SetLineColor(kRed);
     line->Draw();
     
     TLine *line2 = new TLine(0,0,0,1000);
     line2->SetX1(mAngle[2]);
     line2->SetX2(mAngle[2]);
     line2->SetY1(0);
     line2->SetY2(1.2);
     line2->SetLineColor(kBlue);
     line2->Draw();
     
     TLine *line3 = new TLine(0,0,0,1000);
     line3->SetLineStyle(2);
     line3->SetX1(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
     line3->SetX2(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
     line3->SetY1(0);
     line3->SetY2(1.2);
     line3->SetLineColor(1);
     line3->Draw();
     
     TLine *line4 = new TLine(0,0,0,1000);
     line4->SetLineStyle(2);
     line4->SetX1(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
     line4->SetX2(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
     line4->SetY1(0);
     line4->SetY2(1.2);
     line4->SetLineColor(1);
     line4->Draw();
     
     TLegend *leg = new TLegend(0.1,0.5,0.4,0.85);
     leg->SetFillColor(0);
     leg->SetFillStyle(0);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->AddEntry(hAngle[2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
     leg->AddEntry(hAngle[3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
     leg->AddEntry(hAngle[2],Form("#sigma_{c}^{#pi} = %2.1f mrad",spr[2]*1000),"");
     leg->AddEntry(hAngle[3],Form("#sigma_{c}^{K} = %2.1f mrad",spr[3]*1000),"");
     leg->Draw();
     
     TLegend *lnpa = new TLegend(0.7,0.67,0.9,0.85);
     lnpa->SetFillColor(0);
     lnpa->SetFillStyle(0);
     lnpa->SetBorderSize(0);
     lnpa->SetFillStyle(0);
     lnpa->AddEntry(hAngle[2],"pions","lp");
     lnpa->AddEntry(hAngle[3],"kaons","lp");
     lnpa->Draw();
     
     gStyle->SetOptTitle(0);
     gStyle->SetOptStat(0);
     gStyle->SetOptFit(0);
     
     //////////////////////////////////////////
     //////// calculate separation power //////
     //////////////////////////////////////////
     
     glx_canvasAdd("r_separation",800,400);
     
     
     TF1 *ff;
     double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0;
     if(hLnDiff[3]->GetEntries()>100){
     hLnDiff[3]->Fit("gaus","S");
     ff = hLnDiff[3]->GetFunction("gaus");
     ff->SetLineColor(1);
     m1=ff->GetParameter(1);
     s1=ff->GetParameter(2);
     }
     if(hLnDiff[2]->GetEntries()>100){
     hLnDiff[2]->Fit("gaus","S");
     ff = hLnDiff[2]->GetFunction("gaus");
     ff->SetLineColor(1);
     m2=ff->GetParameter(1);
     s2=ff->GetParameter(2);
     }
     if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
     
     hLnDiff[2]->SetTitle(Form("sep = %2.2f s.d.",sep));
     hLnDiff[3]->SetTitle(Form("sep = %2.2f s.d.",sep));
     hLnDiff[2]->Draw();
     hLnDiff[3]->Draw("same");
     
     TLegend *lnpl = new TLegend(0.7,0.67,0.9,0.85);
     lnpl->SetFillColor(0);
     lnpl->SetFillStyle(0);
     lnpl->SetBorderSize(0);
     lnpl->SetFillStyle(0);
     lnpl->AddEntry(hLnDiff[2],"pions","lp");
     lnpl->AddEntry(hLnDiff[3],"kaons","lp");
     lnpl->Draw();
     
     */
    
    glx_canvasSave(2,0);
    //glx_canvasDel("*");
    
    
    
    cout<<"####### @ 3.5 GeV/c fAngleK "<< fAngleK[1]<<"  fAnglePi "<<fAnglePi[1]<<endl;
    
    timer.Stop();
    
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    
    
    return 0;
}
