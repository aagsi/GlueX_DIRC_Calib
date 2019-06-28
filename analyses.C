#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <TLegend.h>
#include "glxtools.C"
#include "THStack.h"
//#include "TH3D.h"
#include "TNtuple.h"
#include "TStopwatch.h"


TStopwatch timer;

//root analyses.C'("out2.root")'
int analyses(TString infile="outFile_v3.root"){
    
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
    
    TH1F*  histo_track_yield_bar_allMom[26];
    for(Int_t i=0; i<26; i++) {
        histo_track_yield_bar_allMom[i] = new TH1F(Form("histo_track_yield_bar_allMom_%d",i),Form(";Photon Yield @ bar %d ; Photon Yield; Entries [#]",i) ,100 ,0,100);
    }
    
    TH1F*  histo_track_resolution_bar_allMom[26][102];
    for(Int_t i=0; i<26; i++) {
        for(Int_t j=0; j<102; j++) {
            histo_track_resolution_bar_allMom[i][j] = new TH1F(Form("histo_track_resolution_bar_allMom_%d_%d",i,j), Form("Cherenkov track resolution @ bar %d @ yield bin %d all mom;Expected #theta_{c}- Measured #theta_{c} [m rad]; Entries [#]",i,j) , 100, -50, 50 );
        }
    }
    //////////////////
    
    TH1F*  histo_track_yield_bar_mom[26][10];
    for(Int_t i=0; i<26; i++) {
        for(Int_t j=0; j<10; j++) {
            int k =j+2;
            int l= j+3;
            histo_track_yield_bar_mom[i][j] = new TH1F(Form("histo_track_yield_bar_mom_%d_%d",i,j),Form("bar %d mom %d - %d GeV/c ;Photon Yield ; Photon Yield; Entries [#]",i,k,l) ,100 ,0,100);
        }
    }
    
    TH1F*  histo_track_resolution_bar_mom[26][10][102];
    for(Int_t i=0; i<26; i++) {
        for(Int_t j=0; j<10; j++) {
            for(Int_t k=0; k<102; k++) {
                int m =j+2;
                int l= j+3;
                
                histo_track_resolution_bar_mom[i][j][k] = new TH1F(Form("histo_track_resolution_bar_allMom_%d_%d_%d",i,j,k), Form("Cherenkov track resolution @ bar %d mom %d - %d GeV/c @ yield bin %d ;Expected #theta_{c}- Measured #theta_{c} [m rad]; Entries [#]",i,m,l,k) , 100, -50, 50 );
            }
        }
    }
    // deleted
    //////////////////
    
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
    
    //TH1F*  histo_chiNDF = new TH1F("histo_chiNDF",";ChiSquare/NDF; entries [#]",200 ,0,10);
    //TH1F*  histo_chiNDF_cut = new TH1F("histo_chiNDF_cut",";ChiSquare/NDF; entries [#]",200 ,0,10);
    
    const int pos_bin(100), pos_min(-100), pos_max(100);
    
    TH2F * histo_pos_xy = new TH2F( "histo_pos_xy" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield_tmp = new TH2F( "histo_pos_xy_yield_tmp" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield = new TH2F( "histo_pos_xy_yield" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_spr_tmp = new TH2F( "histo_pos_xy_spr_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_spr = new TH2F( "histo_pos_xy_spr" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_reso_tmp = new TH2F( "histo_pos_xy_reso_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_reso = new TH2F( "histo_pos_xy_reso" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    // independant
    const int pos_bin_shift(50); // best 200
    TH2F * histo_pos_xy_occupancy = new TH2F( "histo_pos_xy_occupancy" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    TH2F * histo_pos_xy_shift_tmp = new TH2F( "histo_pos_xy_shift_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    TH2F * histo_pos_xy_shift = new TH2F( "histo_pos_xy_shift" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    TH2F * histo_pos_xy_shiftEx_tmp = new TH2F( "histo_pos_xy_shiftEx_tmp" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    TH2F * histo_pos_xy_shiftEx = new TH2F( "histo_pos_xy_shiftEx" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    TH2F * histo_pos_xy_shiftEx_tmp_positive = new TH2F( "histo_pos_xy_shiftEx_tmp_positive" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    TH2F * histo_pos_xy_shiftEx_positive = new TH2F( "histo_pos_xy_shiftEx_positive" , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    TH2F * histo_pos_xy_shiftEx_tmp_negative  = new TH2F( "histo_pos_xy_shiftEx_tmp_negative " , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    TH2F * histo_pos_xy_shiftEx_negative  = new TH2F( "histo_pos_xy_shiftEx_negative " , "; Bar Hit X [cm]; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    TH2F * histo_pos_xy_occupancy_postiveShift = new TH2F( "histo_pos_xy_occupancy_postiveShift" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    TH2F * histo_pos_xy_occupancy_negativeShift = new TH2F( "histo_pos_xy_occupancy_negativeShift" , "; Bar Hit X [cm] ; Bar Hit Y [cm]", pos_bin_shift, pos_min, pos_max, pos_bin_shift, pos_min, pos_max);
    
    
    TH1F* histo_track_mean = new TH1F("histo_track_mean","; Mean per Track [rad]; entries [#]",250,0.817,0.8348);
    TH1F* histo_track_spr = new TH1F("histo_track_spr","; SPR per Track [m rad]; entries [#]",250,5.1,20);
    
    const int nbin_yield =100;
    const int nbin_mom =10;
    
    TH1F*  histo_track_yield[10];
    TH1F*  histo_track_mean_mom[10];
    TH1F*  histo_track_spr_mom[10];
    
    TH1F*  histo_track_resolution_bin[10][100];
    TH1F*  histo_track_spr_bin[10][100];
    TH1F*  histo_track_mean_bin[10][100];
    
    for(Int_t i=0; i<10; i++){
        int kk = i+2;
        int ll =i+3;
        histo_track_yield[i] = new TH1F(Form("histo_track_yield_%d",i), Form("Photon Yield @ momentum %d - %d GeV/c ; Photon Yield; Entries [#]",kk,ll) ,100 ,0,100);
        histo_track_mean_mom[i] = new TH1F(Form("histo_track_mean_mom_%d",i), Form("Track Mean @ momentum %d - %d GeV/c ; Mean [rad]; Entries [#]",kk,ll) ,250,0.817,0.8348);
        histo_track_spr_mom[i] = new TH1F(Form("histo_track_spr_mom_%d",i), Form("Track SPR @ momentum %d - %d GeV/c ; SPR [mrad]; Entries [#]",kk,ll) ,250,5.1,20);
        
        for(Int_t j=0; j<100; j++) {
            int k = i+2;
            int l =i+3;
            histo_track_resolution_bin[i][j] = new TH1F(Form("histo_track_resolution_%d_mom_%d_%d",j,k,l ), Form("Cherenkov track resolution @ yield bin %d momentum %d - %d GeV/c;Expected #theta_{c}- Measured #theta_{c} [m rad]; Entries [#]",j,k,l) , 100, -50, 50 );
            histo_track_spr_bin[i][j] = new TH1F(Form("histo_spr_resolution_%d_mom_%d",j,i), Form("SPR @ yield bin %d momentum %d - %d GeV/c ;SPR [m rad];  [#]",j,k,l) ,250,5.1,20);
            histo_track_mean_bin[i][j] = new TH1F(Form("histo_mean_resolution_%d_mom_%d",j,i), Form("#theta_{c}^{tr} @ yield bin %d momentum %d - %d GeV/c ;#theta_{c}^{tr}  [rad]; Entries [#]",j,k,l), 100,0.817,0.8348);
        }
    }
    
    TH1F*  histo_track_mean_bar[26];
    for(Int_t i=0; i<26; i++){
        histo_track_mean_bar[i] = new TH1F(Form("histo_track_mean_bar_%d",i), Form("Track Mean @ bar %d ; Mean [rad]; Entries [#]",i) ,150,0.817,0.8348);
    }
    
    
    TH1F*  histo_track_mean_bar_mom[26][10];
    for(Int_t i=0; i<26; i++){
        for(Int_t j=0; j<10; j++){
            int k = j+2;
            int l =j+3;
            histo_track_mean_bar_mom[i][j] = new TH1F(Form("histo_track_mean_bar_mom_%d_%d_%d",i,k,l), Form("Track Mean @ bar %d mom %d - %d ; Mean [rad]; Entries [#]",i,k,l) ,50,0.817,0.8348);
        }
    }
    
    TH1F * histo_track_pos_resolution_bin[pos_bin][pos_bin];
    TH1F * histo_track_pos_mom_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_spr_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_mean_bin[pos_bin][pos_bin];
    
    for(Int_t x=0; x<pos_bin; x++) {
        for(Int_t y=0; y<pos_bin; y++) {
            
            histo_track_pos_resolution_bin[x][y] = new TH1F(Form("histo_track_resolution_xbin_%d_ybin_%d",x,y), Form("Cherenkov track resolution @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 100, -50, 50 );
            histo_track_pos_mom_bin[x][y] = new TH1F(Form("histo_track_pos_mom_bin_%d_ybin_%d",x,y), Form("Cherenkov track mom @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 100, 0, 10 );
            
            //histo_track_pos_spr_bin[x][y] = new TH1F(Form("histo_track_pos_spr_xbin_%d_ybin_%d",x,y), Form("SPR @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 250,5.1,20);
            //histo_track_pos_mean_bin[x][y] = new TH1F(Form("histo_track_pos_mean_xbin_%d_ybin_%d",x,y), Form("#theta_{c}^{tr} @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) ,  100,0.817,0.8348);
        }
    }
    
    // graphs
    
    TGraphAsymmErrors *g_yield_mom = new TGraphAsymmErrors();
    g_yield_mom->SetMarkerColor(kBlue);
    g_yield_mom->SetMarkerStyle(21);
    g_yield_mom->SetLineColor(kBlack);
    
    TGraph *g_pi[nbin_mom];
    TGraph *g_k[nbin_mom];
    
    TGraphAsymmErrors *graph_reso[nbin_mom];
    TGraphAsymmErrors *graph_spr[nbin_mom];
    TGraphAsymmErrors *graph_mean[nbin_mom];
    
    for(Int_t j=0; j<nbin_mom; j++){
        
        g_pi[j] = new TGraph();
        g_pi[j]->SetMarkerColor(kBlue);
        g_pi[j]->SetMarkerStyle(0);
        g_pi[j]->SetLineColor(kBlue);
        g_pi[j]->SetLineWidth(3);
        g_pi[j]->SetLineColor(j+1);
        
        g_k[j] = new TGraph();
        g_k[j]->SetMarkerColor(kBlue);
        g_k[j]->SetMarkerStyle(0);
        g_k[j]->SetLineColor(kRed);
        g_k[j]->SetLineWidth(3);
        g_k[j]->SetLineColor(j+1);
        g_k[j]->SetLineStyle(2);
        
        graph_reso[j] = new TGraphAsymmErrors();
        graph_reso[j]->SetTitle("Cherenkov track resolution vs Photon yield");
        graph_reso[j]->SetMarkerColor(j+1);
        graph_reso[j]->SetMarkerStyle(21);
        graph_reso[j]->SetLineColor(j+1);
        
        graph_spr[j] = new TGraphAsymmErrors();
        graph_spr[j]->SetTitle("Cherenkov Track SPR vs Photon yield");
        graph_spr[j]->SetMarkerColor(j+1);
        graph_spr[j]->SetMarkerStyle(21);
        graph_spr[j]->SetLineColor(j+1);
        
        graph_mean[j] = new TGraphAsymmErrors();
        graph_mean[j]->SetTitle("Cherenkov Track Mean vs Photon yield");
        graph_mean[j]->SetMarkerColor(j+1);
        graph_mean[j]->SetMarkerStyle(21);
        graph_mean[j]->SetLineColor(j+1);
        
    }
    
    TGraphAsymmErrors *graph_pos_reso = new TGraphAsymmErrors();
    graph_pos_reso->SetTitle("Cherenkov track resolution vs Photon yield");
    graph_pos_reso->SetMarkerColor(4);
    graph_pos_reso->SetMarkerStyle(21);
    graph_pos_reso->SetLineStyle(0);
    
    TGraphAsymmErrors *graph_reso_allMom[nbar];
    for(Int_t j=0; j<nbar; j++){
        graph_reso_allMom[j] = new TGraphAsymmErrors();
        graph_reso_allMom[j]->SetTitle("Cherenkov track resolution vs Photon yield");
        graph_reso_allMom[j]->SetMarkerColor(j+10);
        graph_reso_allMom[j]->SetMarkerStyle(21);
        graph_reso_allMom[j]->SetLineColor(1);
    }
    
    
    TGraphAsymmErrors *graph_reso_mom[nbar][nbin_mom];
    for(Int_t i=0; i<nbar; i++){
        for(Int_t j=0; j<nbin_mom; j++){
            graph_reso_mom[i][j] = new TGraphAsymmErrors();
            graph_reso_mom[i][j]->SetTitle("Cherenkov track resolution vs Photon yield");
            graph_reso_mom[i][j]->SetMarkerColor(j+10);
            graph_reso_mom[i][j]->SetMarkerStyle(21);
            graph_reso_mom[i][j]->SetLineColor(1);
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
    
    Long64_t nentries = tree_variables->GetEntries();
    Double_t mean_array[]={0.824512,0.825366,0.826363,0.826648,0.826861,0.826576 ,0.826149}; // 1st version
    //Double_t mean_array[]={0.8245,0.826,0.8265,0.82685,0.8269,0.8269,0.8269}; // 2nd version
    // 3rd version
    double diff_array_detailed[26][10]={0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165,0.827165};
    
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
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
        
        //if(track_nbar!=5) continue;
        
        if(track_pid !=2 ) continue; // select pion !=2
        
        //if(track_yield<10)continue;
        
        histo_track_mean->Fill(track_mean);
        histo_track_spr->Fill(track_spr*1000);
        
        
        // mean SPR cut
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        // wall cut
        //if(track_nbar<4 || track_nbar>8 ) continue;
        //if(track_xbar>10 || track_xbar < -10 ) continue;
        
        // momentum cut
        //if(track_mom<3.5 || track_mom>4.0 )continue;
        
        
        //diff = fAnglePi-track_mean; // not used because there are systematic shifts
        diff = track_mean-   0.82608; // default value will be changed
        
        if(track_mom>1 && track_mom<2){
            mom_bin_flag=0;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>2 && track_mom<3){
            mom_bin_flag=1;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>3 && track_mom<4){
            mom_bin_flag=2;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>4 && track_mom<5){
            mom_bin_flag=3;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>5 && track_mom<6){
            mom_bin_flag=4;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>6 && track_mom<7){
            mom_bin_flag=5;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>7 && track_mom<8){
            mom_bin_flag=6;
            diff = track_mean - diff_array_detailed[track_nbar][mom_bin_flag];
            //diff = track_mean - mean_array[mom_bin_flag];
        }
        else{
            continue;
            
        }
        // histogram trk mean per momentum per bar
        histo_track_mean_bar_mom[track_nbar][mom_bin_flag]->Fill(track_mean);
        
        // histogram trk mean per bar
        histo_track_mean_bar[track_nbar]->Fill(track_mean);
        
        // resolution per bar
        histo_track_yield_bar_allMom[track_nbar]->Fill(track_yield);
        int xbin_yield_allMom = histo_track_yield_bar_allMom[track_nbar]->GetXaxis()->FindBin(track_yield);
        histo_track_resolution_bar_allMom[track_nbar][xbin_yield_allMom]->Fill(diff*1000);
        
        // resolution per mom bar bar
        histo_track_yield_bar_mom[track_nbar][mom_bin_flag]->Fill(track_yield);
        int xbin_yield_mom= histo_track_yield_bar_mom[track_nbar][mom_bin_flag]->GetXaxis()->FindBin(track_yield);
        histo_track_resolution_bar_mom[track_nbar][mom_bin_flag][xbin_yield_mom]->Fill(diff*1000);
        
        
        // histogram yield and trk mean
        histo_track_yield[mom_bin_flag]->Fill(track_yield);
        histo_track_mean_mom[mom_bin_flag]->Fill(track_mean);
        histo_track_spr_mom[mom_bin_flag]->Fill(track_spr*1000);
        
        int xbin_yield = histo_track_yield[mom_bin_flag]->GetXaxis()->FindBin(track_yield);
        //cout<<xbin_yield<<endl;
        histo_track_resolution_bin[mom_bin_flag][xbin_yield]->Fill(diff*1000);
        histo_track_spr_bin[mom_bin_flag][xbin_yield]->Fill(track_spr*1000);
        histo_track_mean_bin[mom_bin_flag][xbin_yield]->Fill(track_mean);
        
        //Yield map
        histo_pos_xy->Fill(track_xbar,track_ybar);
        x_pos_bin = histo_pos_xy->GetXaxis()->FindBin(track_xbar);
        y_pos_bin = histo_pos_xy->GetYaxis()->FindBin(track_ybar);
        content_histo_pos_xy=histo_pos_xy->GetBinContent(x_pos_bin,y_pos_bin);
        
        histo_pos_xy_yield_tmp->Fill(track_xbar,track_ybar,track_yield);
        content_histo_pos_xy_tmp=histo_pos_xy_yield_tmp->GetBinContent(x_pos_bin,y_pos_bin);
        average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        histo_pos_xy_yield->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        
        //SPR map
        histo_pos_xy_spr_tmp->Fill(track_xbar,track_ybar,track_spr*1000);
        content_histo_pos_xy_tmp=histo_pos_xy_spr_tmp->GetBinContent(x_pos_bin,y_pos_bin);
        average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
        histo_pos_xy_spr->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
        
        //Resolution map
        int xbin_pos = histo_pos_xy->GetXaxis()->FindBin(track_xbar);
        int ybin_pos = histo_pos_xy->GetYaxis()->FindBin(track_ybar);
        //cout<<"#######"<< xbin_pos <<"   "<<ybin_pos<<endl;
        histo_track_pos_resolution_bin[xbin_pos][ybin_pos]->Fill(diff*1000);
        histo_track_pos_mom_bin[xbin_pos][ybin_pos]->Fill(track_mom);
        //histo_track_pos_spr_bin[xbin_pos][ybin_pos]->Fill(track_spr*1000);
        //histo_track_pos_mean_bin[xbin_pos][ybin_pos]->Fill(track_mean);
        
        // shift map
        
        if(true){
            histo_pos_xy_occupancy->Fill(track_xbar,track_ybar);
            x_pos_bin = histo_pos_xy_occupancy->GetXaxis()->FindBin(track_xbar);
            y_pos_bin = histo_pos_xy_occupancy->GetYaxis()->FindBin(track_ybar);
            content_histo_pos_xy=histo_pos_xy_occupancy->GetBinContent(x_pos_bin,y_pos_bin);
            
            
            histo_pos_xy_shift_tmp->Fill(track_xbar,track_ybar,diff*1000);
            content_histo_pos_xy_tmp=histo_pos_xy_shift_tmp->GetBinContent(x_pos_bin,y_pos_bin);
            average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
            histo_pos_xy_shift->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
            
            
            double ExAnglePi= acos(sqrt(track_mom*track_mom + mass[2]*mass[2])/track_mom/1.4738);
            double ExMeandiff = track_mean - ExAnglePi;
            histo_pos_xy_shiftEx_tmp->Fill(track_xbar,track_ybar,ExMeandiff*1000);
            content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp->GetBinContent(x_pos_bin,y_pos_bin);
            average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
            histo_pos_xy_shiftEx->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
            
            if (ExMeandiff > 0) {
                
                histo_pos_xy_occupancy_postiveShift->Fill(track_xbar,track_ybar);
                x_pos_bin = histo_pos_xy_occupancy_postiveShift->GetXaxis()->FindBin(track_xbar);
                y_pos_bin = histo_pos_xy_occupancy_postiveShift->GetYaxis()->FindBin(track_ybar);
                content_histo_pos_xy=histo_pos_xy_occupancy_postiveShift->GetBinContent(x_pos_bin,y_pos_bin);
                
                histo_pos_xy_shiftEx_tmp_positive->Fill(track_xbar,track_ybar,ExMeandiff*1000);
                content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp_positive->GetBinContent(x_pos_bin,y_pos_bin);
                average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
                histo_pos_xy_shiftEx_positive->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
                
            }else{
                histo_pos_xy_occupancy_negativeShift->Fill(track_xbar,track_ybar);
                x_pos_bin = histo_pos_xy_occupancy_negativeShift->GetXaxis()->FindBin(track_xbar);
                y_pos_bin = histo_pos_xy_occupancy_negativeShift->GetYaxis()->FindBin(track_ybar);
                content_histo_pos_xy=histo_pos_xy_occupancy_negativeShift->GetBinContent(x_pos_bin,y_pos_bin);
                
                histo_pos_xy_shiftEx_tmp_negative ->Fill(track_xbar,track_ybar,ExMeandiff*1000);
                content_histo_pos_xy_tmp=histo_pos_xy_shiftEx_tmp_negative ->GetBinContent(x_pos_bin,y_pos_bin);
                average_bin= content_histo_pos_xy_tmp/content_histo_pos_xy;
                histo_pos_xy_shiftEx_negative ->SetBinContent(x_pos_bin,y_pos_bin,average_bin);
            }
        }
        
        
    }
    
    cout<<"##### start analyses "<<endl;
    if(true){
        
        // yield per bar per mom
        TCanvas *canvas4 = new TCanvas("canvas4","canvas4",800,500);
        for(int i=0;i<nbar;i++){
            for(int j=0;j<nbin_mom;j++){
                //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
                if(histo_track_yield_bar_mom[i][j]->GetEntries()<1)continue;
                canvas4->cd();
                canvas4->Update();
                histo_track_yield_bar_mom[i][j]->Draw();
                canvas4->Update();
                canvas4->WaitPrimitive();
            }
            
        }
        delete canvas4;
        
        
        // diff array per bar per mom  used
        TCanvas *canvas1 = new TCanvas("canvas1","canvas1",800,500);
        TF1 *fit_gause2 = new TF1("fit_gause2","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
        fit_gause2->SetLineColor(kBlack);
        fit_gause2->SetParameters(100,9,2);
        fit_gause2->SetParNames("p0","mean ","sigma");
        fit_gause2->SetParLimits(0,0.1,1E6);
        fit_gause2->SetParLimits(1,0.818,0.834);
        fit_gause2->SetParLimits(2,0.0001,0.005);
        
        //double diff_array_detailed[26][10]={0.82608};
        for(int i=0;i<26;i++){
            for(int j=0;j<10;j++){
                //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
                if(histo_track_mean_bar_mom[i][j]->GetEntries()<100)continue;
                histo_track_mean_bar_mom[i][j]->Fit("fit_gause2","M","", 0, 30);
                diff_array_detailed[i][j]=fit_gause2->GetParameter(1);
                canvas1->cd();
                canvas1->Update();
                histo_track_mean_bar_mom[i][j]->Draw();
                canvas1->Update();
                TLine *lineEXm1= new TLine(0,0,0,1000);
                lineEXm1->SetX1(diff_array_detailed[i][j]);
                lineEXm1->SetX2(diff_array_detailed[i][j]);
                lineEXm1->SetY1(gPad->GetUymin());
                lineEXm1->SetY2(gPad->GetUymax());
                lineEXm1->SetLineColor(kRed);
                lineEXm1->Draw();
                canvas1->Update();
                canvas1->WaitPrimitive();
                delete histo_track_mean_bar_mom[i][j];
            }
        }
        delete canvas1;
        
        // printing diff array per bar per mom
        //    cout<<"diff_array_detailed[i][j]={";
        //    for(int i=0;i<26;i++){
        //        for(int j=0;j<10;j++){
        //            cout<<","<<fit_gause2->GetParameter(1);
        //        }
        //    }
        //    cout<<"};"<<endl;
        
    }
    if(false){
        // diff array per mom  Not used
        TF1 *fit_gause = new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
        fit_gause->SetLineColor(kBlack);
        fit_gause->SetParameters(100,9,2);
        fit_gause->SetParNames("p0","mean ","sigma");
        fit_gause->SetParLimits(0,0.1,1E6);
        
        TCanvas *canvas2 = new TCanvas("canvas2","canvas2",800,500);
        for(int i=0;i<nbin_mom;i++){
            if(histo_track_mean_mom[i]->GetEntries()<1)continue;
            
            fit_gause->SetParLimits(1,0.80,0.84);
            fit_gause->SetParLimits(2,0.0001,500);
            histo_track_mean_mom[i]->Fit("fit_gause","MQ0","", 0, 30) ;
            
            canvas2->cd();
            canvas2->Update();
            histo_track_mean_mom[i]->Draw();
            canvas2->Update();
            TLine *lineEX= new TLine(0,0,0,1000);
            lineEX->SetX1(mean_array[i]);
            lineEX->SetX2(mean_array[i]);
            lineEX->SetY1(gPad->GetUymin());
            lineEX->SetY2(gPad->GetUymax());
            lineEX->SetLineColor(kBlack);
            lineEX->Draw();
            canvas2->Update();
            
            canvas2->WaitPrimitive();
            
            int binmax = histo_track_mean_mom[i]->GetMaximumBin();
            double x = histo_track_mean_mom[i]->GetXaxis()->GetBinCenter(binmax);
            //cout<<"#############  "<< fit_gause->GetParameter(1)<<endl;
            delete histo_track_mean_mom[i];
        }
        delete canvas2;
        
        // diff array per bar Not used not used
        Double_t meanbar_array[26] {0.828,0.827,0.8277,0.826,0.8269,0.8254,0.8239,0.827,0.8275,0.8255,0.824,0.8285,0.829,0.8266,0.828,0.827,0.827,0.829,0.826,0.828,0.829,0.827,0.826,0.826};
        TCanvas *canvas3 = new TCanvas("canvas3","canvas3",800,500);
        for(int i=0;i<26;i++){
            //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
            if(histo_track_mean_bar[i]->GetEntries()<1)continue;
            canvas3->cd();
            canvas3->Update();
            histo_track_mean_bar[i]->Draw();
            canvas3->Update();
            TLine *lineEXm= new TLine(0,0,0,1000);
            lineEXm->SetX1(meanbar_array[i]);
            lineEXm->SetX2(meanbar_array[i]);
            lineEXm->SetY1(gPad->GetUymin());
            lineEXm->SetY2(gPad->GetUymax());
            lineEXm->SetLineColor(kBlack);
            lineEXm->Draw();
            canvas3->Update();
            canvas3->WaitPrimitive();
            delete histo_track_mean_bar[i];
        }
        delete canvas3;
        
        // yield per mom
        TCanvas *canvas5 = new TCanvas("canvas5","canvas5",800,500);
        for(int i=0;i<nbin_mom;i++){
            //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
            if(histo_track_yield[i]->GetEntries()<1)continue;
            canvas5->cd();
            canvas5->Update();
            histo_track_yield[i]->Draw();
            canvas5->Update();
            canvas5->WaitPrimitive();
        }
        delete canvas5;
    }
    
    Double_t spr_array[10]={0};
    TF1 *fit_gause_spr = new TF1("fit_gause_spr","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    fit_gause_spr->SetLineColor(kBlack);
    fit_gause_spr->SetParameters(100,9,2);
    fit_gause_spr->SetParNames("p0","mean ","sigma");
    fit_gause_spr->SetParLimits(0,0.1,1E6);
    
    TCanvas *canvas6 = new TCanvas("canvas6","canvas6",800,500);
    for(int i=0;i<nbin_mom;i++){
        if(histo_track_spr_mom[i]->GetEntries()<1)continue;
        fit_gause_spr->SetParLimits(1,5,12);
        fit_gause_spr->SetParLimits(2,0.0001,500);
        histo_track_spr_mom[i]->Fit("fit_gause_spr","MQ)","", 0, 30) ;
        spr_array[i]=fit_gause_spr->GetParameter(1);
        canvas6->cd();
        canvas6->Update();
        histo_track_spr_mom[i]->Draw();
        canvas6->Update();
        canvas6->WaitPrimitive();
        delete histo_track_spr_mom[i];
    }
    delete canvas6;
    
    cout<<"##### End of diff array printing, yield histograming, track mean histograming "<<endl;
    TCanvas *cc = new TCanvas("cc","cc",800,500);
    
    // fitting functions
    TF1 *fit_track_resolution = new TF1("fit_track_resolution","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-50,50);
    fit_track_resolution->SetLineColor(kBlack);
    fit_track_resolution->SetParameters(100,0,2);
    fit_track_resolution->SetParNames("p0","mean","resolution");
    fit_track_resolution->SetParLimits(0,0.1,1E6);
    fit_track_resolution->SetParLimits(1,-1,1);
    fit_track_resolution->SetParLimits(2,1,20); //5
    fit_track_resolution->SetNpx(1000);
    
    TF1 *fit_track_spr = new TF1("fit_track_spr","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    fit_track_spr->SetLineColor(kBlack);
    fit_track_spr->SetParameters(100,9,2);
    fit_track_spr->SetParNames("p0","mean of sigma","sigma of sigma");
    fit_track_spr->SetParLimits(0,0.1,1E6);
    fit_track_spr->SetParLimits(1,2,12); //(1,8,11);
    fit_track_spr->SetParLimits(2,1,5);
    fit_track_spr->SetNpx(1000);
    
    TF1 *fit_track_mean = new TF1("fit_track_mean","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    fit_track_mean->SetLineColor(kBlack);
    fit_track_mean->SetParameters(100,9,2);
    fit_track_mean->SetParNames("p0","mean of mean","sigma of mean");
    fit_track_mean->SetParLimits(0,0.1,1E6);
    fit_track_mean->SetParLimits(1,0.80,0.84);
    fit_track_mean->SetParLimits(2,0.001,500);
    fit_track_mean->SetNpx(1000);
    
    TF1 *fit_trk_reso = new TF1("fit_trk_reso","[0]*sqrt(([1]/sqrt(x))^2+[2])",7,48);
    fit_trk_reso->SetLineColor(kRed);
    fit_trk_reso->SetParameters(100,9,2);
    fit_trk_reso->SetParNames("p0","SPR","Reso");
    fit_trk_reso->SetParLimits(0,0.1,1E6);
    fit_trk_reso->SetParLimits(1,8,12);
    fit_trk_reso->SetParLimits(2,1,5);
    fit_trk_reso->SetNpx(1000);
    
    int couter[nbin_mom]={0};
    for(int f=0;f<nbin_mom;f++){
        fit_track_resolution->SetLineColor(f+1);
        fit_track_spr->SetLineColor(f+1);
        fit_track_mean->SetLineColor(f+1);
        
        for (int i=0;i<nbin_yield;i++){
            if (histo_track_resolution_bin[f][i]->GetEntries() <1500)continue; //400 //175 // 200 // 500
            histo_track_resolution_bin[f][i]->Fit("fit_track_resolution","M","", -50, 50) ;
            histo_track_spr_bin[f][i]->Fit("fit_track_spr","MQ0","", 0, 30) ;
            histo_track_mean_bin[f][i]->Fit("fit_track_mean","MQ0","", 0, 30) ;
            if(false){
                cc->cd();
                cc->Update();
                histo_track_resolution_bin[f][i]->Draw();
                cc->Update();
                cc->WaitPrimitive();
            }
            if(false){
                cc->cd();
                cc->Update();
                histo_track_spr_bin[f][i]->Draw();
                cc->Update();
                cc->WaitPrimitive();
            }
            
            if(false){
                cc->cd();
                cc->Update();
                histo_track_mean_bin[f][i]->Draw();
                cc->Update();
                cc->WaitPrimitive();
            }
            yield_BinCenter = histo_track_yield[f]->GetXaxis()->GetBinCenter(i);
            track_resolution= fit_track_resolution->GetParameter(2);
            track_resolution_error= fit_track_resolution->GetParError(2);
            
            track_mean_bin= fit_track_mean->GetParameter(1);
            //track_mean_error= fit_track_mean->GetParError(1);
            track_mean_error= fit_track_mean->GetParameter(2);
            
            //spr
            //track_spr_bin= fit_track_spr->GetParameter(2);
            //track_spr_error= fit_track_spr->GetParError(2);
            
            //track_spr_error= fit_track_spr->GetParError(1);
            
            //track_spr_bin= histo_track_spr_bin[f][i]->GetStdDev();
            //track_spr_error= histo_track_spr_bin[f][i]->GetStdDevError();
            
            track_spr_bin= histo_track_spr_bin[f][i]->GetMean();
            //track_spr_error= histo_track_spr_bin[f][i]->GetMeanError();
            track_spr_error= histo_track_spr_bin[f][i]->GetStdDev();
            
            //track_spr_bin= fit_track_spr->GetParameter(1);
            //track_spr_error= fit_track_spr->GetParameter(2);
            
            graph_reso[f]->SetPoint(couter[f], yield_BinCenter, track_resolution);
            graph_reso[f]->SetPointError(couter[f], 1/2, 1/2,track_resolution_error/2,track_resolution_error/2);
            
            graph_spr[f]->SetPoint(couter[f], yield_BinCenter, track_spr_bin);
            graph_spr[f]->SetPointError(couter[f], 1/2, 1/2,track_spr_error/2,track_spr_error/2);
            
            graph_mean[f]->SetPoint(couter[f], yield_BinCenter, track_mean_bin);
            graph_mean[f]->SetPointError(couter[f], 1/2, 1/2,track_mean_error/2,track_mean_error/2);
            
            g_pi[f]->SetPoint(couter[f], yield_BinCenter, fAnglePi[f]);
            g_k[f]->SetPoint(couter[f], yield_BinCenter, fAngleK[f]);
            
            ++couter[f];
        }
    }
    cout<<"##### End resolution per Mom "<<endl;
    
    
    // per bar
    int couter_all[nbar]={0};
    
    for (int i=0;i<nbar;i++){
        fit_track_resolution->SetLineColor(i+1);
        fit_track_spr->SetLineColor(i+1);
        fit_track_mean->SetLineColor(i+1);
        
        for (int j=0;j<nbin_yield;j++){
            if (histo_track_resolution_bar_allMom[i][j]->GetEntries() <4000)continue; // 1500
            histo_track_resolution_bar_allMom[i][j]->Fit("fit_track_resolution","M","", -50, 50) ;
            if(false){
                cc->cd();
                cc->Update();
                histo_track_resolution_bar_allMom[i][j]->Draw();
                cc->Update();
                cc->WaitPrimitive();
            }
            yield_BinCenter = histo_track_yield_bar_allMom[i]->GetXaxis()->GetBinCenter(j);
            track_resolution= fit_track_resolution->GetParameter(2);
            track_resolution_error= fit_track_resolution->GetParError(2);
            graph_reso_allMom[i]->SetPoint(couter_all[i], yield_BinCenter, track_resolution);
            graph_reso_allMom[i]->SetPointError(couter_all[i], 1/2, 1/2,track_resolution_error/2,track_resolution_error/2);
            ++couter_all[i];
        }
    }
    cout<<"##### End resolution per Bar "<<endl;
    
    // per bar per mom
    int couter_bar_mom[nbar][nbin_mom]={0};
    
    for (int i=0;i<nbar;i++){
        for(int j=0;j<nbin_mom;j++){
            
            fit_track_resolution->SetLineColor(i+1);
            fit_track_spr->SetLineColor(i+1);
            fit_track_mean->SetLineColor(i+1);
            
            for (int k=0;k<nbin_yield;k++){
                if (histo_track_resolution_bar_mom[i][j][k]->GetEntries() <1000)continue;
                histo_track_resolution_bar_mom[i][j][k]->Fit("fit_track_resolution","MQ0","", -50, 50) ;
                if(false){
                    cc->cd();
                    cc->Update();
                    histo_track_resolution_bar_mom[i][j][k]->Draw();
                    cc->Update();
                    cc->WaitPrimitive();
                }
                yield_BinCenter = histo_track_yield_bar_mom[i][j]->GetXaxis()->GetBinCenter(k);
                track_resolution= fit_track_resolution->GetParameter(2);
                track_resolution_error= fit_track_resolution->GetParError(2);
                //cout<<"#####"<<yield_BinCenter<<"  "<<track_resolution<<endl;
                graph_reso_mom[i][j]->SetPoint(couter_bar_mom[i][j], yield_BinCenter, track_resolution);
                graph_reso_mom[i][j]->SetPointError(couter_bar_mom[i][j], 1/2, 1/2,track_resolution_error/2,track_resolution_error/2);
                ++couter_bar_mom[i][j];
            }
        }
    }
    
    
    glx_canvasAdd("r_resolution_bin_allMom",800,400);
    TLegend * legend_reso_bar= new TLegend(0.630326, 0.466667,0.889724,0.872);
    TMultiGraph *mg_allMom = new TMultiGraph();
    
    for(int i=0;i<nbar;i++){
        if( graph_reso_allMom[i]->GetN()>22){
            mg_allMom->Add(graph_reso_allMom[i]);
            legend_reso_bar->AddEntry(graph_reso_allMom[i],Form("Bar %d",i) ,"P");
        }
    }
    mg_allMom->SetTitle(" Cherenkov Resolution per Track all Momenta; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
    mg_allMom->Draw("APL");
    legend_reso_bar->Draw();
    //////////////////////
    
    
    glx_canvasAdd("r_resolution_bin_bar_mom",800,400);
    TLegend * legend_reso_bar_mom= new TLegend(0.630326, 0.466667,0.889724,0.872);
    TMultiGraph *mg_bar_mom = new TMultiGraph();
    for(int i=0;i<nbar;i++){
        for(int j=0;j<nbin_mom;j++){
            int k = j+2;
            int l = j+3;
            
            //cout<<"#####  N points ="<<graph_reso_mom[i][j]->GetN()<<endl;
            if( graph_reso_mom[i][j]->GetN()>3){//22
                mg_bar_mom->Add(graph_reso_mom[i][j]);
                legend_reso_bar_mom->AddEntry(graph_reso_mom[i][j],Form("Bar %d Mom %d - %d GeV/c",i,k,l) ,"P");
            }
        }
    }
    mg_bar_mom->SetTitle(" Cherenkov Resolution per track per Bar per Momentum; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
    mg_bar_mom->Draw("APL");
    legend_reso_bar_mom->Draw();
    //////////////////////
    
    
    
    glx_canvasAdd("r_resolution_bin",800,400);
    TMultiGraph *mg = new TMultiGraph();
    TGraphAsymmErrors *graph_tracker_reso = new TGraphAsymmErrors();
    graph_tracker_reso->SetMarkerColor(kBlack);
    graph_tracker_reso->SetMarkerStyle(21);
    TLegend * legend_reso= new TLegend(0.630326, 0.466667,0.889724,0.872);
    //legend_reso->SetHeader("Pions","C");
    int counter2=0;
    double tracker_reso(-1),tracker_reso_error(-1);
    for(int i=0;i<nbin_mom;i++){
        fit_trk_reso->SetLineColor(i+1);
        if(graph_reso[i]->GetN()<1)continue;
        legend_reso->AddEntry(graph_reso[i],Form("%d : %d GeV/c",i+2,i+3) ,"l");
        cout<<"########## i= "<<i<<" spr_array[i]= "<< spr_array[i]<<endl;
        fit_trk_reso->SetParameter(1,spr_array[i]); // 9
        graph_reso[i]->Fit("fit_trk_reso","M","", 7, 48) ;
        tracker_reso=fit_trk_reso->GetParameter(2);
        tracker_reso_error=fit_trk_reso->GetParError(2);
        
        graph_tracker_reso->SetPoint(counter2,i+2+0.5, tracker_reso);
        graph_tracker_reso->SetPointError(counter2, 1/2, 1/2,tracker_reso_error/2,tracker_reso_error/2);
        mg->Add(graph_reso[i]);
        ++counter2;
    }
    mg->SetTitle(" Cherenkov Resolution per Track ; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
    mg->Draw("APL");
    legend_reso->Draw();
    //////////////////////
    
    // trck shift
    glx_canvasAdd("r_pos_shift",800,400);
    //histo_pos_xy_shift->SetMinimum(5);
    //histo_pos_xy_shift->SetMaximum(10);
    histo_pos_xy_shift->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_shift->Draw("colz");
    
    glx_canvasAdd("r_pos_shiftEx",800,400);
    //histo_pos_xy_shiftEx->SetMinimum(5);
    //histo_pos_xy_shiftEx->SetMaximum(10);
    histo_pos_xy_shiftEx->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_shiftEx->Draw("colz");
    
    glx_canvasAdd("r_pos_shiftEx_positive",800,400);
    //histo_pos_xy_shiftEx_positive->SetMinimum(5);
    //histo_pos_xy_shiftEx_positive->SetMaximum(10);
    histo_pos_xy_shiftEx_positive->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_shiftEx_positive->Draw("colz");
    
    glx_canvasAdd("r_pos_shiftEx_negative",800,400);
    //histo_pos_xy_shiftEx_negative->SetMinimum(5);
    //histo_pos_xy_shiftEx_negative->SetMaximum(10);
    histo_pos_xy_shiftEx_negative->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_shiftEx_negative->Draw("colz");
    
    glx_canvasAdd("r_pos_shiftEx_occu_positive",800,400);
    //histo_pos_xy_occupancy_postiveShift->SetMinimum(5);
    //histo_pos_xy_occupancy_postiveShift->SetMaximum(10);
    histo_pos_xy_occupancy_postiveShift->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_occupancy_postiveShift->Draw("colz");
    
    glx_canvasAdd("r_pos_shiftEx_occu_negative",800,400);
    //histo_pos_xy_occupancy_negativeShift->SetMinimum(5);
    //histo_pos_xy_occupancy_negativeShift->SetMaximum(10);
    histo_pos_xy_occupancy_negativeShift->GetYaxis()->SetRangeUser(-100,0);
    histo_pos_xy_occupancy_negativeShift->Draw("colz");
    cout<<"##### End trk shift histograms per bar "<<endl;
    ////////////
    
    if(true){
        ////////////
        glx_canvasAdd("r_resolution_bin_fit",800,400);
        TMultiGraph *mg_reo_fit = new TMultiGraph();
        mg_reo_fit->Add(graph_tracker_reso);
        mg_reo_fit->SetTitle(" Tracker Resolution ; Pion Momentum [GeV/c]; #sigma_{tracker} [m rad]");
        mg_reo_fit->Draw("APL");
        mg_reo_fit->GetHistogram()->GetYaxis()->SetRangeUser(0,5);
        mg_reo_fit->GetHistogram()->GetXaxis()->SetRangeUser(0,10);
        glx_canvasGet("r_resolution_bin_fit")->Update();
        TLine *test2= new TLine(0,0,0,1000);
        test2->Draw();
        
        ////////////
        glx_canvasAdd("r_spr_bin",800,400);
        TMultiGraph *mg2 = new TMultiGraph();
        for(int i=0;i<nbin_mom;i++){
            fit_track_spr->SetLineColor(i+1);
            if(graph_spr[i]->GetN()<1)continue;
            mg2->Add(graph_spr[i]);
        }
        mg2->SetTitle(" SPR per Track ; Photon Yield [#]; SPR [m rad]");
        mg2->Draw("APL");
        mg2->GetHistogram()->GetYaxis()->SetRangeUser(0,15);
        glx_canvasGet("r_spr_bin")->Update();
        legend_reso->Draw();
        
        
        /////////////
        glx_canvasAdd("r_mean_bin",800,400);
        TMultiGraph *mg3 = new TMultiGraph();
        for(int i=0;i<nbin_mom;i++){
            fit_track_mean->SetLineColor(i+1);
            if(graph_mean[i]->GetN()<1)continue;
            mg3->Add(graph_mean[i]);
            if (i==1)mg3->Add(g_pi[i]);
            if (i==1)mg3->Add(g_k[i]);
        }
        mg3->SetTitle(" #theta_{c}^{tr} per Track ; Photon Yield [#]; #theta_{c}^{tr} [rad]");
        mg3->Draw("APL");
        mg3->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.84);
        legend_reso->Draw();
        
        /////////////
        glx_canvasAdd("r_yield",800,400);
        THStack *hs = new THStack("hs","Stacked 1D histograms");
        hs->SetTitle("Photon Yield");
        int counter3=0;
        for(int i=0;i<nbin_mom ;i++){  //  mom_bin_flag
            if( histo_track_yield[i]->GetEntries()<1)continue;
            histo_track_yield[i]->SetLineColor(i+1);
            g_yield_mom->SetPoint(counter3, i+2 , histo_track_yield[i]->GetMean());
            g_yield_mom->SetPointError(counter3 ,1/2,1/2, histo_track_yield[i]->GetStdDev()/2,histo_track_yield[i]->GetStdDev()/2);
            hs->Add(histo_track_yield[i]);
            ++counter3;
        }
        hs->Draw("nostack");
        hs->GetYaxis()->SetTitle("Entries [#]");
        hs->GetXaxis()->SetTitle("Photon Yield [#]");
        legend_reso->Draw();
        
        glx_canvasAdd("r_yield_graph",800,400);
        TMultiGraph *mg4 = new TMultiGraph();
        mg4->Add(g_yield_mom);
        mg4->SetTitle(" Photon Yield; Pion Momentun [GeV/c]; Photon Yield [#] ");
        mg4->Draw("APL");
        mg4->GetHistogram()->GetYaxis()->SetRangeUser(0,45);
        glx_canvasGet("r_yield_graph")->Update();
        TLine *tes= new TLine(0,0,0,0);
        tes->Draw();
        
        glx_canvasAdd("r_pos",800,400);
        histo_pos_xy->GetYaxis()->SetRangeUser(-100,0);
        histo_pos_xy->Draw("colz");
        
        glx_canvasAdd("r_pos_spr",800,400);
        histo_pos_xy_spr->SetMinimum(5);
        histo_pos_xy_spr->SetMaximum(10);
        histo_pos_xy_spr->GetYaxis()->SetRangeUser(-100,0);
        histo_pos_xy_spr->Draw("colz");
        
        glx_canvasAdd("r_pos_yield",800,400);
        histo_pos_xy_yield->SetMinimum(5);
        histo_pos_xy_yield->SetMaximum(37);
        histo_pos_xy_yield->GetYaxis()->SetRangeUser(-100,0);
        histo_pos_xy_yield->Draw("colz");
        
        /////////////
        glx_canvasAdd("r_histo_track_mean",800,400);
        histo_track_mean->Draw();
        glx_canvasGet("r_histo_track_mean")->Update();
        TLine *lin_mean_max= new TLine(0,0,0,1000);
        lin_mean_max->SetX1(mean_max);
        lin_mean_max->SetX2(mean_max);
        lin_mean_max->SetY1(gPad->GetUymin());
        lin_mean_max->SetY2(gPad->GetUymax());
        lin_mean_max->SetLineColor(kBlack);
        lin_mean_max->Draw();
        TLine *lin_mean_min= new TLine(0,0,0,1000);
        lin_mean_min->SetX1(mean_min);
        lin_mean_min->SetX2(mean_min);
        lin_mean_min->SetY1(gPad->GetUymin());
        lin_mean_min->SetY2(gPad->GetUymax());
        lin_mean_min->SetLineColor(kBlack);
        lin_mean_min->Draw();
        glx_canvasGet("r_histo_track_mean")->Update();
        
        /////////////
        glx_canvasAdd("r_spr",800,400);
        histo_track_spr->Draw();
        glx_canvasGet("r_spr")->Update();
        TLine *lin_spr_max= new TLine(0,0,0,1000);
        lin_spr_max->SetX1(spr_max);
        lin_spr_max->SetX2(spr_max);
        lin_spr_max->SetY1(gPad->GetUymin());
        lin_spr_max->SetY2(gPad->GetUymax());
        lin_spr_max->SetLineColor(kBlack);
        lin_spr_max->Draw();
        TLine *lin_spr_min= new TLine(0,0,0,1000);
        lin_spr_min->SetX1(spr_min);
        lin_spr_min->SetX2(spr_min);
        lin_spr_min->SetY1(gPad->GetUymin());
        lin_spr_min->SetY2(gPad->GetUymax());
        lin_spr_min->SetLineColor(kBlack);
        lin_spr_min->Draw();
        glx_canvasGet("r_spr")->Update();
        
        // histograms
        //glx_canvasAdd("r_chiNDF",800,400);
        //histo_chiNDF->Draw();
        //histo_chiNDF_cut->Draw("same");
    }
    
    // resolution map
    cout<<"##### commint 4 "<<endl;
    if(true){
        // 4D resolution moentum position
        Float_t xpos,ypos,mom_pos,reso_pos;
        auto f2 = TFile::Open("reso_pos_mom.root","RECREATE");
        TNtuple ntuple("ntuple","data from ascii file","xpos:ypos:mom_pos:reso_pos");
        
        double mom_avr(-1);
        TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
        int couter2(0);
        for (int x=0;x<pos_bin;x++){
            for (int y=0;y<pos_bin;y++){
                double hentry =histo_track_pos_resolution_bin[x][y]->GetEntries();
                if (hentry <100)continue; //100
                histo_track_pos_resolution_bin[x][y]->Fit("fit_track_resolution","MQ0","", -50, 50) ;
                mom_avr = histo_track_pos_mom_bin[x][y]->GetMean();
                if(false){
                    cc2->cd();
                    cc2->Update();
                    histo_track_pos_resolution_bin[x][y]->Draw();
                    cc2->Update();
                    cc2->WaitPrimitive();
                }
                pos_BinCenter = histo_pos_xy->GetBinContent(x,y);
                
                track_pos_resolution= fit_track_resolution->GetParameter(2);
                track_pos_resolution_error= fit_track_resolution->GetParError(2);
                
                histo_pos_xy_reso->SetBinContent(x,y,track_pos_resolution);
                //histo_pos_xy_reso_4d->SetBinContent(x,y,mom_avr,track_pos_resolution);
                
                double  xpos = histo_pos_xy_reso->GetXaxis()->GetBinCenter(x);
                double  ypos = histo_pos_xy_reso->GetYaxis()->GetBinCenter(y);
                double  mom_pos = mom_avr;
                double  reso_pos = track_pos_resolution;
                if (track_pos_resolution<5 && track_pos_resolution>1 )ntuple.Fill(xpos,ypos,mom_pos,reso_pos);
                
                graph_pos_reso->SetPoint(couter2, mom_pos, track_pos_resolution); // pos_BinCenter
                //graph_pos_reso->SetPointError(couter2, 1/2, 1/2,track_pos_resolution_error/2,track_pos_resolution_error/2);
                
                ++couter2;
            }
        }
        f2->Write();
        //f2->Close();
        
        glx_canvasAdd("r_pos_resolution_map",800,400);
        histo_pos_xy_reso->SetMinimum(1);
        histo_pos_xy_reso->SetMaximum(5);
        histo_pos_xy_reso->GetYaxis()->SetRangeUser(-100,0);
        histo_pos_xy_reso->Draw("colz");
        
        glx_canvasAdd("r_pos_resolution_map_4d",800,400);
        ntuple.SetMarkerStyle(20);
        ntuple.SetMarkerSize(2);
        ntuple.Draw("ypos:xpos:mom_pos:reso_pos","","COLZ");
        // to run with relavent axis title use root d3_reso.C
        
        glx_canvasAdd("r_graph_pos_reso",800,400);
        TMultiGraph *mg_4d = new TMultiGraph();
        mg_4d->Add(graph_pos_reso);
        mg_4d->SetTitle(" Averaged Cherenkov Resolution per Track all DIRC wall all Photon yield; Momentum [GeV/c]; #sigma( #theta_{c}^{tr} ) [m rad]");
        mg_4d->Draw("AP");
        mg_4d->GetHistogram()->GetYaxis()->SetRangeUser(1,6);
        glx_canvasGet("r_graph_pos_reso")->Update();
        TLine *tet= new TLine(0,0,0,0);
        tet->Draw();
    }
    
    // kinematics
    if(true){
        cout<<"##### kinematics histogramming "<<endl;
        glx_canvasAdd("r_rho_mass",800,400);
        hist_ev_rho_mass->Draw();
        glx_canvasGet("r_rho_mass")->Update();
        TLine *lin_rho_mass_max= new TLine(0,0,0,1000);
        lin_rho_mass_max->SetX1(mass_rho_max);
        lin_rho_mass_max->SetX2(mass_rho_max);
        lin_rho_mass_max->SetY1(gPad->GetUymin());
        lin_rho_mass_max->SetY2(gPad->GetUymax());
        lin_rho_mass_max->SetLineColor(kBlack);
        lin_rho_mass_max->Draw();
        TLine *lin_rho_mass_min= new TLine(0,0,0,1000);
        lin_rho_mass_min->SetX1(mass_rho_mini);
        lin_rho_mass_min->SetX2(mass_rho_mini);
        lin_rho_mass_min->SetY1(gPad->GetUymin());
        lin_rho_mass_min->SetY2(gPad->GetUymax());
        lin_rho_mass_min->SetLineColor(kBlack);
        lin_rho_mass_min->Draw();
        glx_canvasGet("r_rho_mass")->Update();
        
        glx_canvasAdd("r_phi_mass",800,400);
        hist_ev_phi_mass->Draw();
        glx_canvasGet("r_phi_mass")->Update();
        TLine *lin_phi_mass_max= new TLine(0,0,0,1000);
        lin_phi_mass_max->SetX1(mass_phi_max);
        lin_phi_mass_max->SetX2(mass_phi_max);
        lin_phi_mass_max->SetY1(gPad->GetUymin());
        lin_phi_mass_max->SetY2(gPad->GetUymax());
        lin_phi_mass_max->SetLineColor(kBlack);
        lin_phi_mass_max->Draw();
        TLine *lin_phi_mass_min= new TLine(0,0,0,1000);
        lin_phi_mass_min->SetX1(mass_phi_mini);
        lin_phi_mass_min->SetX2(mass_phi_mini);
        lin_phi_mass_min->SetY1(gPad->GetUymin());
        lin_phi_mass_min->SetY2(gPad->GetUymax());
        lin_phi_mass_min->SetLineColor(kBlack);
        lin_phi_mass_min->Draw();
        glx_canvasGet("r_phi_mass")->Update();
        
        glx_canvasAdd("r_missing_mass_phi",800,400);
        hist_ev_missing_mass_phi->Draw();
        glx_canvasGet("r_missing_mass_phi")->Update();
        TLine *lin_phi_miss_mass_max= new TLine(0,0,0,1000);
        lin_phi_miss_mass_max->SetX1(miss_mass_phi_max);
        lin_phi_miss_mass_max->SetX2(miss_mass_phi_max);
        lin_phi_miss_mass_max->SetY1(gPad->GetUymin());
        lin_phi_miss_mass_max->SetY2(gPad->GetUymax());
        lin_phi_miss_mass_max->SetLineColor(kBlack);
        lin_phi_miss_mass_max->Draw();
        TLine *lin_phi_miss_mass_min= new TLine(0,0,0,1000);
        lin_phi_miss_mass_min->SetX1(miss_mass_phi_mini);
        lin_phi_miss_mass_min->SetX2(miss_mass_phi_mini);
        lin_phi_miss_mass_min->SetY1(gPad->GetUymin());
        lin_phi_miss_mass_min->SetY2(gPad->GetUymax());
        lin_phi_miss_mass_min->SetLineColor(kBlack);
        lin_phi_miss_mass_min->Draw();
        glx_canvasGet("r_missing_mass_phi")->Update();
        
        glx_canvasAdd("r_missing_mass_rho",800,400);
        hist_ev_missing_mass_rho->Draw();
        glx_canvasGet("r_missing_mass_rho")->Update();
        TLine *lin_rho_miss_mass_max= new TLine(0,0,0,1000);
        lin_rho_miss_mass_max->SetX1(miss_mass_rho_max);
        lin_rho_miss_mass_max->SetX2(miss_mass_rho_max);
        lin_rho_miss_mass_max->SetY1(gPad->GetUymin());
        lin_rho_miss_mass_max->SetY2(gPad->GetUymax());
        lin_rho_miss_mass_max->SetLineColor(kBlack);
        lin_rho_miss_mass_max->Draw();
        TLine *lin_rho_miss_mass_min= new TLine(0,0,0,1000);
        lin_rho_miss_mass_min->SetX1(miss_mass_rho_mini);
        lin_rho_miss_mass_min->SetX2(miss_mass_rho_mini);
        lin_rho_miss_mass_min->SetY1(gPad->GetUymin());
        lin_rho_miss_mass_min->SetY2(gPad->GetUymax());
        lin_rho_miss_mass_min->SetLineColor(kBlack);
        lin_rho_miss_mass_min->Draw();
        glx_canvasGet("r_missing_mass_rho")->Update();
        
        glx_canvasAdd("r_chi_phi",800,400);
        hist_ev_chi_phi->Draw();
        glx_canvasGet("r_chi_phi")->Update();
        TLine *lin_chi_phi_max= new TLine(0,0,0,1000);
        lin_chi_phi_max->SetX1(chisq_phi_max);
        lin_chi_phi_max->SetX2(chisq_phi_max);
        lin_chi_phi_max->SetY1(gPad->GetUymin());
        lin_chi_phi_max->SetY2(gPad->GetUymax());
        lin_chi_phi_max->SetLineColor(kBlack);
        lin_chi_phi_max->Draw();
        TLine *lin_chi_phi_min= new TLine(0,0,0,1000);
        lin_chi_phi_min->SetX1(chisq_phi_mini);
        lin_chi_phi_min->SetX2(chisq_phi_mini);
        lin_chi_phi_min->SetY1(gPad->GetUymin());
        lin_chi_phi_min->SetY2(gPad->GetUymax());
        lin_chi_phi_min->SetLineColor(kBlack);
        lin_chi_phi_min->Draw();
        glx_canvasGet("r_chi_phi")->Update();
        
        glx_canvasAdd("r_chi_rho",800,400);
        hist_ev_chi_rho->Draw();
        glx_canvasGet("r_chi_rho")->Update();
        TLine *lin_chi_rho_max= new TLine(0,0,0,1000);
        lin_chi_rho_max->SetX1(chisq_rho_max);
        lin_chi_rho_max->SetX2(chisq_rho_max);
        lin_chi_rho_max->SetY1(gPad->GetUymin());
        lin_chi_rho_max->SetY2(gPad->GetUymax());
        lin_chi_rho_max->SetLineColor(kBlack);
        lin_chi_rho_max->Draw();
        TLine *lin_chi_rho_min= new TLine(0,0,0,1000);
        lin_chi_rho_min->SetX1(chisq_rho_mini);
        lin_chi_rho_min->SetX2(chisq_rho_mini);
        lin_chi_rho_min->SetY1(gPad->GetUymin());
        lin_chi_rho_min->SetY2(gPad->GetUymax());
        lin_chi_rho_min->SetLineColor(kBlack);
        lin_chi_rho_min->Draw();
        glx_canvasGet("r_chi_rho")->Update();
    }
    cout<<"##### closing  "<<endl;
    
    glx_canvasSave(2,0);
    //glx_canvasDel("*");
    
    
    cout<<"##### histograms saved"<<endl;
    if(false){
        //delete graphs
        delete g_yield_mom;
        delete graph_pos_reso;
        for(Int_t j=0; j<nbin_mom; j++){
            delete g_pi[j];
            delete g_k[j];
            delete graph_reso[j];
            delete graph_spr[j];
            delete graph_mean[j];
        }
        for(Int_t j=0; j<nbar; j++){
            delete graph_reso_allMom[j];
        }
        for(Int_t i=0; i<nbar; i++){
            for(Int_t j=0; j<nbin_mom; j++){
                delete graph_reso_mom[i][j];
            }
        }
        cout<<"##### graphs deleted "<<endl;
        //delete histograms
        for(Int_t i=0; i<26; i++) {
            delete histo_track_yield_bar_allMom[i];
        }
        cout<<"##### no problem 1"<<endl;
        
        for(Int_t i=0; i<26; i++) {
            for(Int_t j=0; j<102; j++) {
                delete histo_track_resolution_bar_allMom[i][j];
            }
        }
        
        cout<<"##### no problem 2"<<endl;
        for(Int_t i=0; i<26; i++) {
            for(Int_t j=0; j<10; j++) {
                delete histo_track_yield_bar_mom[i][j];
            }
        }
        cout<<"##### no problem 3"<<endl;
        for(Int_t i=0; i<26; i++) {
            for(Int_t j=0; j<10; j++) {
                for(Int_t k=0; k<102; k++) {
                    
                    delete histo_track_resolution_bar_mom[i][j][k];
                }
            }
        }
        cout<<"##### no problem 4"<<endl;
        for(Int_t i=0; i<10; i++){
            delete histo_track_yield[i];
            for(Int_t j=0; j<100; j++) {
                delete histo_track_resolution_bin[i][j];
                delete histo_track_spr_bin[i][j];
                delete histo_track_mean_bin[i][j];
            }
        }
        cout<<"##### no problem 5"<<endl;
        for(Int_t x=0; x<pos_bin; x++) {
            for(Int_t y=0; y<pos_bin; y++) {
                delete histo_track_pos_resolution_bin[x][y];
                delete histo_track_pos_mom_bin[x][y];
            }
        }
        delete histo_cherenkov;
        delete histo_tdiff;
        delete hist_ev_rho_mass;
        delete hist_ev_phi_mass;
        delete hist_ev_missing_mass_phi;
        delete hist_ev_missing_mass_rho;
        delete hist_ev_chi_phi;
        delete hist_ev_chi_rho;
        delete histo_pos_xy;
        delete histo_pos_xy_yield_tmp;
        delete histo_pos_xy_yield;
        delete histo_pos_xy_spr_tmp;
        delete histo_pos_xy_spr;
        delete histo_pos_xy_reso_tmp;
        delete histo_pos_xy_reso;
        delete histo_pos_xy_occupancy;
        delete histo_pos_xy_shift_tmp;
        delete histo_pos_xy_shift;
        delete histo_pos_xy_shiftEx_tmp;
        delete histo_pos_xy_shiftEx;
        delete histo_pos_xy_shiftEx_tmp_positive;
        delete histo_pos_xy_shiftEx_positive;
        delete histo_pos_xy_shiftEx_tmp_negative;
        delete histo_pos_xy_shiftEx_negative;
        delete histo_pos_xy_occupancy_postiveShift;
        delete histo_pos_xy_occupancy_negativeShift;
        delete histo_track_mean;
        delete histo_track_spr;
        cout<<"##### histograms deleted"<<endl;
    }
    
    cout<<"####### @ 3.5 GeV/c fAngleK "<< fAngleK[1]<<"  fAnglePi "<<fAnglePi[1]<<endl;
    
    timer.Stop();
    
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    return 0;
}
