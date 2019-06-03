#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <TLegend.h>
#include "glxtools.C"
void analyses(){
    
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    
    // calculate cherenkov angle
    Double_t momentum=3.5;
    Int_t pdg[]= {11,13,211,321,2212};
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.009),range(5*sigma),noise(0.3);
    Double_t fAngleK = acos(sqrt(momentum*momentum+ mass[3]*mass[3])/momentum/1.4738)-0.00;
    Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00;
    
    
    
    // histograms
    
    //TH1F*  histo_chiNDF = new TH1F("histo_chiNDF",";ChiSquare/NDF; entries [#]",200 ,0,10);
    //TH1F*  histo_chiNDF_cut = new TH1F("histo_chiNDF_cut",";ChiSquare/NDF; entries [#]",200 ,0,10);
    
    const int pos_bin(100), pos_min(-100), pos_max(100);
    
    TH2F * histo_pos_xy = new TH2F( "histo_pos_xy" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield_tmp = new TH2F( "histo_pos_xy_yield_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield = new TH2F( "histo_pos_xy_yield" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_spr_tmp = new TH2F( "histo_pos_xy_spr_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_spr = new TH2F( "histo_pos_xy_spr" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_reso_tmp = new TH2F( "histo_pos_xy_reso_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_reso = new TH2F( "histo_pos_xy_reso" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH1F* histo_track_mean = new TH1F("histo_track_mean",";Track Mean [rad]; entries [#]",250,0.817,0.8348);
    TH1F* histo_track_spr = new TH1F("histo_track_spr",";Track SPR [m rad]; entries [#]",250,5.1,20);
    
    
    
    const int nbin_yield =100;
    const int nbin_mom =10;
    
    
    TH1F*  histo_track_yield[10];
    TH1F*  histo_track_mean_mom[10];
    TH1F*  histo_track_spr_mom[10];
    
    TH1F*  histo_track_resolution_bin[10][100];
    TH1F*  histo_track_spr_bin[10][100];
    TH1F*  histo_track_mean_bin[10][100];
    
    for(Int_t i=0; i<10; i++){
        histo_track_yield[i] = new TH1F(Form("histo_track_yield_%d",i), Form("Photon Yield @ mom bin %d ; Photon Yield; Entries [#]",i) ,100 ,0,100);
        histo_track_mean_mom[i] = new TH1F(Form("histo_track_mean_mom_%d",i), Form("Track Mean @ mom bin %d ; Mean [rad]; Entries [#]",i) ,250,0.817,0.8348);
        histo_track_spr_mom[i] = new TH1F(Form("histo_track_spr_mom_%d",i), Form("Track SPR @ mom bin %d ; SPR [mrad]; Entries [#]",i) ,250,5.1,20);
        
        
        
        for(Int_t j=0; j<100; j++) {
            histo_track_resolution_bin[i][j] = new TH1F(Form("histo_track_resolution_%d_mom_%d",j,i), Form("Cherenkov track resolution @ yield bin %d mom flag %d;Expected - measured [m rad]; Entries [#]",j,i) , 100, -50, 50 );
            histo_track_spr_bin[i][j] = new TH1F(Form("histo_spr_resolution_%d_mom_%d",j,i), Form("SPR @ yield bin %d mom flag %d ;SPR [m rad];  [#]",j,i) ,250,5.1,20);
            histo_track_mean_bin[i][j] = new TH1F(Form("histo_mean_resolution_%d_mom_%d",j,i), Form("#theta_{c}^{tr} @ yield bin %d mom flag %d ;#theta_{c}^{tr}  [rad]; Entries [#]",j,i), 100,0.817,0.8348);
        }
    }
    
    
    TH1F * histo_track_pos_resolution_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_spr_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_mean_bin[pos_bin][pos_bin];
    
    for(Int_t x=0; x<pos_bin; x++) {
        for(Int_t y=0; y<pos_bin; y++) {
            
            histo_track_pos_resolution_bin[x][y] = new TH1F(Form("histo_track_resolution_xbin_%d_ybin_%d",x,y), Form("Cherenkov track resolution @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 100, -50, 50 );
            //histo_track_pos_spr_bin[x][y] = new TH1F(Form("histo_track_pos_spr_xbin_%d_ybin_%d",x,y), Form("SPR @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 250,5.1,20);
            //histo_track_pos_mean_bin[x][y] = new TH1F(Form("histo_track_pos_mean_xbin_%d_ybin_%d",x,y), Form("#theta_{c}^{tr} @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) ,  100,0.817,0.8348);
        }
    }
    
    // graphs
    
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
        g_pi[j]->SetLineWidth(1);
        
        g_k[j] = new TGraph();
        g_k[j]->SetMarkerColor(kBlue);
        g_k[j]->SetMarkerStyle(0);
        g_k[j]->SetLineColor(kRed);
        g_k[j]->SetLineWidth(1);
        
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
    
    //    TGraph *g_calc_trk_res = new TGraph();
    //    g_calc_trk_res->SetMarkerColor(kBlue);
    //    g_calc_trk_res->SetMarkerStyle(0);
    //    g_calc_trk_res->SetLineColor(kRed);
    //    g_calc_trk_res->SetLineWidth(1);
    
    
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
    
    
    // read tree
    TFile *f = new TFile("outFile.root");
    TTree *tree_variables = (TTree*)f->Get("tree_variables");
    double track_spr(-1),track_mean(-1), track_yield(-1), track_mom(-1), track_xbar(0),track_ybar(0);
    double track_fit_chisqu(-1),track_fit_NDF(-1);
    int track_pid(-1), track_nbar(-1);
    
    tree_variables->SetBranchAddress("track_pid",&track_pid);
    tree_variables->SetBranchAddress("track_spr",&track_spr);
    tree_variables->SetBranchAddress("track_mean",&track_mean);
    tree_variables->SetBranchAddress("track_yield",&track_yield);
    tree_variables->SetBranchAddress("track_mom",&track_mom);
    tree_variables->SetBranchAddress("track_xbar",&track_xbar);
    tree_variables->SetBranchAddress("track_ybar",&track_ybar);
    tree_variables->SetBranchAddress("track_nbar",&track_nbar);
    tree_variables->SetBranchAddress("track_fit_chisqu",&track_fit_chisqu);
    tree_variables->SetBranchAddress("track_fit_NDF",&track_fit_NDF);
    
    Long64_t nentries = tree_variables->GetEntries();
    Double_t mean_array[]={0.824512,0.825366,0.826363,0.826648,0.826861,0.826576 ,0.826149};
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        if(track_pid !=2 ) continue; // select pion !=2
        //fit_quality=track_fit_chisqu/track_fit_NDF;
        //histo_chiNDF->Fill(fit_quality);
        
        //if(track_fit_chisqu>100)continue;
        //if(track_yield<10)continue;
        
        histo_track_mean->Fill(track_mean);
        histo_track_spr->Fill(track_spr*1000);
        // momentum cut
        //if(track_mom> 4.5   || track_mom<3.5) continue;
        // mean SPR cut
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        // wall cut
        //if(track_nbar<4 || track_nbar>8 ) continue;
        //if(track_xbar>10 || track_xbar < -10 ) continue;
        
        //if(fit_quality>2) continue;
        //histo_chiNDF_cut->Fill(fit_quality);
        
        //diff = fAnglePi-track_mean;
        diff = track_mean-   0.82608;
        
        
        if(track_mom>1 && track_mom<2){
            mom_bin_flag=0;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>2 && track_mom<3){
            mom_bin_flag=1;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>3 && track_mom<4){
            mom_bin_flag=2;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>4 && track_mom<5){
            mom_bin_flag=3;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>5 && track_mom<6){
            mom_bin_flag=4;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>6 && track_mom<7){
            mom_bin_flag=5;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else if(track_mom>7 && track_mom<8){
            mom_bin_flag=6;
            diff = track_mean - mean_array[mom_bin_flag];
        }
        else{
            continue;
            
        }
        
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
        //histo_track_pos_spr_bin[xbin_pos][ybin_pos]->Fill(track_spr*1000);
        //histo_track_pos_mean_bin[xbin_pos][ybin_pos]->Fill(track_mean);
    }
    if(false){
        TCanvas *cctest = new TCanvas("cctest","cctest",800,500);
        for(int i=0;i<nbin_mom;i++){
            //cout<<"####### i= "<<i<<"  "<<histo_track_yield[i]->GetEntries()<<endl;
            if(histo_track_yield[i]->GetEntries()<1)continue;
            cctest->cd();
            cctest->Update();
            histo_track_yield[i]->Draw();
            cctest->Update();
            cctest->WaitPrimitive();
        }
        
        TCanvas *cctest2 = new TCanvas("cctest2","cctest2",800,500);
        for(int i=0;i<nbin_mom;i++){
            if(histo_track_mean_mom[i]->GetEntries()<1)continue;
            cctest2->cd();
            cctest2->Update();
            histo_track_mean_mom[i]->Draw();
            cctest2->Update();
            cctest2->WaitPrimitive();
            
            int binmax = histo_track_mean_mom[i]->GetMaximumBin();
            double x = histo_track_mean_mom[i]->GetXaxis()->GetBinCenter(binmax);
            
            cout<<"###### i  "<<i<< " x  "<<x<<endl;
        }
        
        TCanvas *cctest3 = new TCanvas("cctest3","cctest3",800,500);
        for(int i=0;i<nbin_mom;i++){
            if(histo_track_spr_mom[i]->GetEntries()<1)continue;
            cctest3->cd();
            cctest3->Update();
            histo_track_spr_mom[i]->Draw();
            cctest3->Update();
            cctest3->WaitPrimitive();
        }
    }
    
    
    TCanvas *cc = new TCanvas("cc","cc",800,500);
    
    // fitting functions
    TF1 *fit_track_resolution = new TF1("fit_track_resolution","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-50,50);
    fit_track_resolution->SetLineColor(kBlack);
    fit_track_resolution->SetParameters(100,0,2);
    fit_track_resolution->SetParNames("p0","mean","resolution");
    fit_track_resolution->SetParLimits(0,0.1,1E6);
    fit_track_resolution->SetParLimits(1,-1,1);
    fit_track_resolution->SetParLimits(2,1,20); //5
    
    TF1 *fit_track_spr = new TF1("fit_track_spr","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    fit_track_spr->SetLineColor(kBlack);
    fit_track_spr->SetParameters(100,9,2);
    fit_track_spr->SetParNames("p0","mean of sigma","sigma of sigma");
    fit_track_spr->SetParLimits(0,0.1,1E6);
    fit_track_spr->SetParLimits(1,8,11);
    fit_track_spr->SetParLimits(2,1,5);
    
    TF1 *fit_track_mean = new TF1("fit_track_mean","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,30);
    fit_track_mean->SetLineColor(kBlack);
    fit_track_mean->SetParameters(100,9,2);
    fit_track_mean->SetParNames("p0","mean of mean","sigma of mean");
    fit_track_mean->SetParLimits(0,0.1,1E6);
    fit_track_mean->SetParLimits(1,0.80,0.84);
    fit_track_mean->SetParLimits(2,0.001,500);
    
    TF1 *fit_trk_reso = new TF1("fit_trk_reso","[0]*sqrt(([1]/sqrt(x))^2+[2])",7,100);
    fit_trk_reso->SetLineColor(kRed);
    fit_trk_reso->SetParameters(100,9,2);
    fit_trk_reso->SetParNames("p0","SPR","Reso");
    fit_trk_reso->SetParLimits(0,0.1,1E6);
    fit_trk_reso->SetParLimits(1,8,12);
    fit_trk_reso->SetParLimits(2,1,5);
    
    int couter[nbin_mom]={0};
    for(int f=0;f<nbin_mom;f++){
        fit_track_resolution->SetLineColor(f+1);
        for (int i=0;i<nbin_yield;i++){
            if (histo_track_resolution_bin[f][i]->GetEntries() <500)continue; //400 //175 // 200
            histo_track_resolution_bin[f][i]->Fit("fit_track_resolution","MQ0","", -50, 50) ;
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
            
            yield_BinCenter = histo_track_yield[mom_bin_flag]->GetXaxis()->GetBinCenter(i);
            
            track_resolution= fit_track_resolution->GetParameter(2);
            track_resolution_error= fit_track_resolution->GetParError(2);
            
            track_mean_bin= fit_track_mean->GetParameter(1);
            //track_mean_error= fit_track_mean->GetParError(1);
            track_mean_error= fit_track_mean->GetParameter(2);
            
            //spr
            //track_spr_bin= fit_track_spr->GetParameter(2);
            //track_spr_error= fit_track_spr->GetParError(2);
            
            track_spr_bin= fit_track_spr->GetParameter(1);
            //track_spr_error= fit_track_spr->GetParError(1);
            
            //track_spr_bin= histo_track_spr_bin[i]->GetStdDev();
            //track_spr_error= histo_track_spr_bin[i]->GetStdDevError();
            
            //track_spr_bin= histo_track_spr_bin[i]->GetMean();
            //track_spr_error= histo_track_spr_bin[i]->GetMeanError();
            //track_spr_error= histo_track_spr_bin[i]->GetStdDev();
            track_spr_error= fit_track_spr->GetParameter(2);
            
            graph_reso[f]->SetPoint(couter[f], yield_BinCenter, track_resolution);
            graph_reso[f]->SetPointError(couter[f], 1/2, 1/2,track_resolution_error/2,track_resolution_error/2);
            
            graph_spr[f]->SetPoint(couter[f], yield_BinCenter, track_spr_bin);
            graph_spr[f]->SetPointError(couter[f], 1/2, 1/2,track_spr_error/2,track_spr_error/2);
            
            graph_mean[f]->SetPoint(couter[f], yield_BinCenter, track_mean_bin);
            graph_mean[f]->SetPointError(couter[f], 1/2, 1/2,track_mean_error/2,track_mean_error/2);
            
            /////////
            // warning
            g_pi[f]->SetPoint(couter[f], yield_BinCenter, fAnglePi);
            g_k[f]->SetPoint(couter[f], yield_BinCenter, fAngleK);
            /////////
            
            //calc_trk_res=sqrt(track_resolution*track_resolution - (track_spr_bin/sqrt(yield_BinCenter))* (track_spr_bin/sqrt(yield_BinCenter)));
            //g_calc_trk_res->SetPoint(couter, yield_BinCenter, calc_trk_res);
            
            ++couter[f];
        }
    }
    
    
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
        
        legend_reso->AddEntry(graph_reso[i],Form("%d GeV/c",i+2) ,"l");
        
        graph_reso[i]->Fit("fit_trk_reso","M","", 7, 100) ;
        tracker_reso=fit_trk_reso->GetParameter(2);
        tracker_reso_error=fit_trk_reso->GetParError(2);
        
        graph_tracker_reso->SetPoint(counter2,i+2, tracker_reso);
        graph_tracker_reso->SetPointError(counter2, 1/2, 1/2,tracker_reso_error/2,tracker_reso_error/2);
        mg->Add(graph_reso[i]);
        ++counter2;
    }
    mg->SetTitle(" Cherenkov Resolution per Track ; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
    mg->Draw("APL");
    TLine *test= new TLine(0,0,0,1000);
    legend_reso->Draw();
    test->Draw();
    
    
    
    glx_canvasAdd("r_resolution_bin_fit",800,400);
    TMultiGraph *mg_reo_fit = new TMultiGraph();
    mg_reo_fit->Add(graph_tracker_reso);
    mg_reo_fit->SetTitle(" Tracker Resolution ; Pion Momentum [GeV/c]; #sigma_{tracker} [m rad]");
    mg_reo_fit->Draw("APL");
    mg_reo_fit->GetHistogram()->GetYaxis()->SetRangeUser(0,5);
    glx_canvasGet("r_resolution_bin_fit")->Update();
    TLine *test2= new TLine(0,0,0,1000);
    test2->Draw();

    //
    //    ///////////
    //    glx_canvasAdd("r_spr_bin",800,400);
    //    TMultiGraph *mg2 = new TMultiGraph();
    //    mg2->Add(graph_spr);
    //    mg2->SetTitle(" SPR per Track ; Photon Yield [#]; SPR [m rad]");
    //    mg2->Draw("APL");
    //    mg2->GetHistogram()->GetYaxis()->SetRangeUser(0,15);
    //    line_k->Draw();// to fix bug on canvas update
    //    glx_canvasGet("r_spr_bin")->Update();
    //
    //    /////////////
    //    glx_canvasAdd("r_mean_bin",800,400);
    //    TMultiGraph *mg3 = new TMultiGraph();
    //    mg3->Add(graph_mean);
    //    mg3->Add(g_pi);
    //    mg3->Add(g_k);
    //    mg3->SetTitle(" #theta_{c}^{tr} per Track ; Photon Yield [#]; #theta_{c}^{tr} [rad]");
    //    mg3->Draw("APL");
    //    mg3->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.84);
    //    glx_canvasGet("r_mean_bin")->Update();
    //    TLine *line_pi= new TLine(0,0,0,1000);
    //    line_pi->SetY1(fAnglePi);
    //    line_pi->SetY2(fAnglePi);
    //    line_pi->SetX1(gPad->GetUymin());
    //    line_pi->SetX2(gPad->GetUymax());
    //    line_pi->SetLineColor(kBlack);
    //    line_pi->Draw();
    //    TLine *line_k= new TLine(0,0,0,1000);
    //    line_k->SetY1(fAngleK);
    //    line_k->SetY2(fAngleK);
    //    line_k->SetX1(gPad->GetUymin());
    //    line_k->SetX2(gPad->GetUymax());
    //    line_k->SetLineColor(kBlack);
    //    line_k->Draw();
    //    glx_canvasGet("r_mean_bin")->Update();
    
    //    if(false){
    //        // histograms
    //        //glx_canvasAdd("r_chiNDF",800,400);
    //        //histo_chiNDF->Draw();
    //        //histo_chiNDF_cut->Draw("same");
    //
    //
    //
    //        glx_canvasAdd("r_pos",800,400);
    //        histo_pos_xy->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_yield",800,400);
    //        histo_pos_xy_yield->Draw("colz");
    //
    //        glx_canvasAdd("r_pos_spr",800,400);
    //        histo_pos_xy_spr->Draw("colz");
    //
    //
    //
    //        glx_canvasAdd("r_yield",800,400);
    //        histo_track_yield[mom_bin_flag]->Draw();
    //        /////////////
    //
    //        glx_canvasAdd("r_histo_track_mean",800,400);
    //        histo_track_mean->Draw();
    //        glx_canvasGet("r_histo_track_mean")->Update();
    //        TLine *lin_mean_max= new TLine(0,0,0,1000);
    //        lin_mean_max->SetX1(mean_max);
    //        lin_mean_max->SetX2(mean_max);
    //        lin_mean_max->SetY1(gPad->GetUymin());
    //        lin_mean_max->SetY2(gPad->GetUymax());
    //        lin_mean_max->SetLineColor(kBlack);
    //        lin_mean_max->Draw();
    //        TLine *lin_mean_min= new TLine(0,0,0,1000);
    //        lin_mean_min->SetX1(mean_min);
    //        lin_mean_min->SetX2(mean_min);
    //        lin_mean_min->SetY1(gPad->GetUymin());
    //        lin_mean_min->SetY2(gPad->GetUymax());
    //        lin_mean_min->SetLineColor(kBlack);
    //        lin_mean_min->Draw();
    //        glx_canvasGet("r_histo_track_mean")->Update();
    //        /////////////
    //
    //        glx_canvasAdd("r_spr",800,400);
    //        histo_track_spr->Draw();
    //        glx_canvasGet("r_spr")->Update();
    //        TLine *lin_spr_max= new TLine(0,0,0,1000);
    //        lin_spr_max->SetX1(spr_max);
    //        lin_spr_max->SetX2(spr_max);
    //        lin_spr_max->SetY1(gPad->GetUymin());
    //        lin_spr_max->SetY2(gPad->GetUymax());
    //        lin_spr_max->SetLineColor(kBlack);
    //        lin_spr_max->Draw();
    //        TLine *lin_spr_min= new TLine(0,0,0,1000);
    //        lin_spr_min->SetX1(spr_min);
    //        lin_spr_min->SetX2(spr_min);
    //        lin_spr_min->SetY1(gPad->GetUymin());
    //        lin_spr_min->SetY2(gPad->GetUymax());
    //        lin_spr_min->SetLineColor(kBlack);
    //        lin_spr_min->Draw();
    //        glx_canvasGet("r_spr")->Update();
    //
    //
    //    }
    
    
    // resolution map
    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
    int couter2(0);
    for (int x=0;x<pos_bin;x++){
        for (int y=0;y<pos_bin;y++){
            double hentry =histo_track_pos_resolution_bin[x][y]->GetEntries();
            if (hentry <100)continue;
            histo_track_pos_resolution_bin[x][y]->Fit("fit_track_resolution","MQ0","", -50, 50) ;
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
            
            graph_pos_reso->SetPoint(couter2, pos_BinCenter, track_pos_resolution);
            graph_pos_reso->SetPointError(couter2, 1/2, 1/2,track_pos_resolution_error/2,track_pos_resolution_error/2);
            
            ++couter2;
        }
    }
    
    glx_canvasAdd("r_pos_resolution_map",800,400);
    
    histo_pos_xy_reso->SetMinimum(1);
    histo_pos_xy_reso->SetMaximum(5);
    histo_pos_xy_reso->Draw("colz");
    
    
    //delete histograms
    for (int x=0;x<pos_bin;x++){
        for (int y=0;y<pos_bin;y++){
            delete histo_track_pos_resolution_bin[x][y];
        }
    }
    
    //
    //    for(Int_t i=0; i<nbin_mom; i++){
    //        delete histo_track_yield[i];
    //        for(Int_t j=0; j<=nbin_yield+1; j++) {
    //            delete histo_track_resolution_bin[i][j];
    //            delete histo_track_spr_bin[i][j];
    //            delete histo_track_mean_bin[i][j];
    //
    //        }
    //    }
    cout<<"####### fAngleK "<< fAngleK<<"  fAnglePi "<<fAnglePi<<endl;
    
    
    //delete [] mean_array;
}
