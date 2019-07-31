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
    
    Double_t fit_angleK[10]={0};
    Double_t fit_anglePi[10]={0};
    
    Double_t momentum[] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
    //Double_t momentum[] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
    for(int f=0;f<10;f++){
        fit_angleK[f] = acos(sqrt(momentum[f]*momentum[f] + mass[3]*mass[3])/momentum[f]/1.4738)-0.00;
        fit_anglePi[f]= acos(sqrt(momentum[f]*momentum[f] + mass[2]*mass[2])/momentum[f]/1.4738)-0.00;
    }
    
    
    // histograms
    const int nbar =26;
    const int pos_min(-100), pos_max(100);
    const int nbin_yield =100;
    const int nbin_mom =10;
    

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
    TH1F* histo_tmp_mom = new TH1F("histo_mom","; Track Momentum [GeV/c]; entries [#]",20,0,10);
    
    
    // variables
    double diff(-1);
    
    int mom_bin_flag(-1);
    int spr_bin_flag(-1);
    
    double mean_min(0.818), mean_max(0.834); // 0.817,0.8348);
    double spr_min(5.1), spr_max(11.5); //  6 11.5
    double calc_trk_res(-1);
    
    
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
    
    
    double noise = 0.5;
    double pion_counter[42][22]={1};
    double kaon_counter[42][22]={1};
    double sum1[42][22],sum2[42][22];
    TSpectrum *spect = new TSpectrum(10);
    double minChisto_Cherenkov = 0.6;
    double maxChisto_Cherenkov = 0.9;
    TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChisto_Cherenkov,maxChisto_Cherenkov);
    double cut_cangle=0.04;
    double cherenkovreco[5],spr[5];
    
    TGaxis::SetMaxDigits(3);
    double sigma[]={0.01,0.01,0.01,0.010,0.01,0.01};

    TH1F *histo_Cherenkov[42][22][5], *histo_LnDiff[42][22][5];

    double theory_angle[5];
    TF1  *fit_angle[5];
    
    double init_mon=3.5;
    for(int i=0; i<5; i++){
        theory_angle[i] = acos(sqrt(init_mon * init_mon + glx_mass[i]*glx_mass[i])/init_mon/1.473);
        fit_angle[i] = new TF1(Form("fit_angle_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
        fit_angle[i]->SetParameter(0,1);        // const
        fit_angle[i]->SetParameter(1,theory_angle[i]);// mean
        fit_angle[i]->SetParameter(2,sigma[i]); // sigma
    }
    
    fit_angle[2]->SetLineColor(4);
    fit_angle[3]->SetLineColor(2);


    for(int i=0; i<42; i++){
        for(int j=0; j<22; j++){
            for(int k=2; k<4; k++){
                histo_LnDiff[i][j][k] = new TH1F(Form("histo_LnDiff_%d_%d_%d",i,j,k),";ln L(#pi) - ln L(K);entries [#]",110,-120,120); //,80,-150,150);
                //histo_LnDiff[i] = new TH1F(Form("histo_LnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",80,-150,150);

                histo_Cherenkov[i][j][k] = new TH1F(Form("histo_Cherenkov_%d_%d_%d",i,j,k),  "cherenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
                histo_Cherenkov[i][j][k]->SetMarkerStyle(20);
                histo_Cherenkov[i][j][k]->SetMarkerSize(0.8);
            }
            histo_LnDiff[i][j][2]->SetLineColor(4);
            histo_LnDiff[i][j][3]->SetLineColor(2);

            histo_Cherenkov[i][j][2]->SetLineColor(4);
            histo_Cherenkov[i][j][3]->SetLineColor(2);
            histo_Cherenkov[i][j][2]->SetMarkerColor(kBlue+1);
            histo_Cherenkov[i][j][3]->SetMarkerColor(kRed+1);
        }
    }

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

        int x_pos_tmp_bin = histo_tmp_pos->GetXaxis()->FindBin(track_xbar);
        int mom_tmp_bin = histo_tmp_mom->GetXaxis()->FindBin(track_mom);

        //////////////////////////////////////////
        //////// calculate separation power //////
        //////////////////////////////////////////

        double percentage = kaon_counter[x_pos_tmp_bin][mom_tmp_bin]/pion_counter[x_pos_tmp_bin][mom_tmp_bin]*100.0;

        //cout<<"##### percentage  "<<percentage<<endl;
        if (!(percentage <100 && track_pid==2) || track_pid==3){
            sum1[x_pos_tmp_bin][mom_tmp_bin]=0;
            sum2[x_pos_tmp_bin][mom_tmp_bin]=0;
            for (UInt_t j = 0; j < vtangle->size(); ++j){
                double time_diff = vtdiff->at(j);
                histo_tdiff->Fill(time_diff);
                if(fabs(time_diff)>3)continue;
                double tangle = vtangle->at(j);

                for(int p=0; p<5; p++){
                    theory_angle[p] = acos(sqrt(track_mom * track_mom + glx_mass[p]*glx_mass[p])/track_mom/1.473);
                    fit_angle[p]->SetParameter(1,theory_angle[p]);
                }
                histo_Cherenkov[x_pos_tmp_bin][mom_tmp_bin][track_pid]->Fill(tangle);
                if(fabs(tangle-0.5*(theory_angle[2]+theory_angle[3]))>cut_cangle)continue;

                sum1[x_pos_tmp_bin][mom_tmp_bin] += TMath::Log(fit_angle[2]->Eval(tangle)+noise);
                sum2[x_pos_tmp_bin][mom_tmp_bin] += TMath::Log(fit_angle[3]->Eval(tangle)+noise);
            }

            double sum = sum1[x_pos_tmp_bin][mom_tmp_bin]-sum2[x_pos_tmp_bin][mom_tmp_bin];
            //cout<<"##### sum  "<<sum<<endl;
            histo_LnDiff[x_pos_tmp_bin][mom_tmp_bin][track_pid]->Fill(sum);

            if(track_pid==3) kaon_counter[x_pos_tmp_bin][mom_tmp_bin]++;
            if(track_pid==2) pion_counter[x_pos_tmp_bin][mom_tmp_bin]++;


        }
        //if(track_nbar>5) continue;
        //if(track_pid !=2 ) continue; // select pion !=2

        /////////////////////
        ////// Hit Map //////
        /////////////////////

        for (UInt_t j = 0; j < vpx->size(); ++j) {
            glx_hdigi[vpx->at(j)]->Fill(vpy->at(j),vpz->at(j));
        }



    }

    TFile file("speration.root","recreate");

    for(int i=0; i<42; i++){
        for(int j=0; j<22; j++){
            for(int k=2; k<4; k++){
                histo_Cherenkov[i][j][k]->Write();
            }
        }
    }

    for(int i=0; i<42; i++){
        for(int j=0; j<22; j++){
            for(int k=2; k<4; k++){
                histo_LnDiff[i][j][k]->Write();
            }
        }
    }

    file.Write();
    file.Close();

    cout<<"##### start analyses "<<endl;
    cout<<"####### @ 3.5 GeV/c fit_angleK "<< fit_angleK[1]<<"  fit_anglePi "<<fit_anglePi[1]<<endl;
    timer.Stop();
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
    
    
    return 0;
    
}
