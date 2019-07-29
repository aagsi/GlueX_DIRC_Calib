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

int analyses4(TString infile="outFile_v3.root"){// outFile_v3.root   out2.root
    
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
    
    
    
    // variables
    double diff(-1);
    
    int mom_bin_flag(-1);
    int spr_bin_flag(-1);
    
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
    
    tree_variables->SetBranchAddress("vpx",&vpx,&bvpx);
    tree_variables->SetBranchAddress("vpy",&vpy,&bvpy);
    tree_variables->SetBranchAddress("vpz",&vpz,&bvpz);
    
    tree_variables->SetBranchAddress("vtdiff",&vtdiff,&bvtdiff);
    tree_variables->SetBranchAddress("vtangle",&vtangle,&bvtangle);
    
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
    
    Long64_t nentries = tree_variables->GetEntries();
    
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        bvpx->GetEntry(i);
        bvpy->GetEntry(i);
        bvpz->GetEntry(i);
        bvtdiff->GetEntry(i);
        bvtangle->GetEntry(i);
        
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
        
        if(track_mom>1 && track_mom<1.5){
            mom_bin_flag=0;
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
        
        
        // SPR Flag
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
        //if(track_nbar>5) continue;
        if(track_pid !=2 ) continue; // select pion !=2
        
        ////////////////////////////////
        ////// New Implimentation //////
        ////////////////////////////////
        
        for (UInt_t j = 0; j < vpx->size(); ++j) {
            glx_hdigi[vpx->at(j)]->Fill(vpy->at(j),vpz->at(j));
        }
        for (UInt_t j = 0; j < vtangle->size(); ++j) {
            double time_diff = vtdiff->at(j);
            if(fabs(time_diff)>3)continue;
            histo_tdiff->Fill(time_diff);
            histo_cherenkov->Fill(vtangle->at(j));
        }
        
    }
    
    cout<<"##### start analyses "<<endl;
    
    if(true){
        cout<<"##### Histograms "<<endl;
        glx_drawDigi("m,p,v\n",0);
        glx_canvasAdd("r_cherenkov",800,400);
        histo_cherenkov->Draw();
        
        glx_canvasAdd("r_tdiff",800,400);
        histo_tdiff->Draw();
        
        
        
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
    
    
    ////////////////////////////////
    //////// Cherenkove angle //////
    ////////////////////////////////
    
    
    glx_canvasAdd("r_angle",800,400);
    
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
    line->SetX1(cherenkovreco[3]); // mAngle[3]
    line->SetX2(cherenkovreco[3]); // mAngle[3]
    line->SetY1(0);
    line->SetY2(1.2);
    line->SetLineColor(kRed);
    line->Draw();
    
    TLine *line2 = new TLine(0,0,0,1000);
    line2->SetX1(cherenkovreco[2]); // mAngle[2]
    line2->SetX2(cherenkovreco[2]); // mAngle[2]
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
    
    /////////////////////////////////////////
    /////// calculate separation power //////
    /////////////////////////////////////////
    
    glx_canvasAdd("r_separation",800,400);
    
    
    TF1 *ff;
    double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0;
    if(hLnDiff[3]->GetEntries()>10){
        hLnDiff[3]->Fit("gaus","S");
        ff = hLnDiff[3]->GetFunction("gaus");
        ff->SetLineColor(1);
        m1=ff->GetParameter(1);
        s1=ff->GetParameter(2);
    }
    if(hLnDiff[2]->GetEntries()>10){
        hLnDiff[2]->Fit("gaus","S");
        ff = hLnDiff[2]->GetFunction("gaus");
        ff->SetLineColor(1);
        m2=ff->GetParameter(1);
        s2=ff->GetParameter(2);
    }
    if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
    
    cout<<"#######  sep= "<<sep<<endl;
    
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
    
    TLegend *leg_sep = new TLegend(0.14787, 0.570667, 0.348371, 0.749333);
    leg_sep->SetFillColor(0);
    leg_sep->SetFillStyle(0);
    leg_sep->SetBorderSize(0);
    leg_sep->SetFillStyle(0);
    leg_sep->AddEntry(hLnDiff[2],Form("Separation = %1.2f",sep),"");
    
    leg_sep->Draw();
    
    
    glx_canvasSave(2,0);
    //glx_canvasDel("*");
    
    
    
    cout<<"####### @ 3.5 GeV/c fAngleK "<< fAngleK[1]<<"  fAnglePi "<<fAnglePi[1]<<endl;
    
    timer.Stop();
    
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    
    
    return 0;
}
