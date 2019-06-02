#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
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
    
    // graphs
    TGraph *g_pi = new TGraph();
    g_pi->SetMarkerColor(kBlue);
    g_pi->SetMarkerStyle(0);
    g_pi->SetLineColor(kBlue);
    g_pi->SetLineWidth(1);
    
    TGraph *g_k = new TGraph();
    g_k->SetMarkerColor(kBlue);
    g_k->SetMarkerStyle(0);
    g_k->SetLineColor(kRed);
    g_k->SetLineWidth(1);
    
    TGraph *g_calc_trk_res = new TGraph();
    g_calc_trk_res->SetMarkerColor(kBlue);
    g_calc_trk_res->SetMarkerStyle(0);
    g_calc_trk_res->SetLineColor(kRed);
    g_calc_trk_res->SetLineWidth(1);
    
    TGraphAsymmErrors *graph_reso = new TGraphAsymmErrors();
    graph_reso->SetTitle("Cherenkov track resolution vs Photon yield");
    graph_reso->SetMarkerColor(4);
    graph_reso->SetMarkerStyle(21);
    
    TGraphAsymmErrors *graph_spr = new TGraphAsymmErrors();
    graph_spr->SetTitle("Cherenkov Track SPR vs Photon yield");
    graph_spr->SetMarkerColor(4);
    graph_spr->SetMarkerStyle(21);
    
    TGraphAsymmErrors *graph_mean = new TGraphAsymmErrors();
    graph_mean->SetTitle("Cherenkov Track Mean vs Photon yield");
    graph_mean->SetMarkerColor(4);
    graph_mean->SetMarkerStyle(21);
    
    TGraphAsymmErrors *graph_pos_reso = new TGraphAsymmErrors();
    graph_pos_reso->SetTitle("Cherenkov track resolution vs Photon yield");
    graph_pos_reso->SetMarkerColor(4);
    graph_pos_reso->SetMarkerStyle(21);
    
    // histograms
    int pos_bin(100), pos_min(-100), pos_max(100);
    
    TH2F * histo_pos_xy = new TH2F( "histo_pos_xy" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield_tmp = new TH2F( "histo_pos_xy_yield_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_yield = new TH2F( "histo_pos_xy_yield" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_spr_tmp = new TH2F( "histo_pos_xy_spr_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_spr = new TH2F( "histo_pos_xy_spr" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    TH2F * histo_pos_xy_reso_tmp = new TH2F( "histo_pos_xy_reso_tmp" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    TH2F * histo_pos_xy_reso = new TH2F( "histo_pos_xy_reso" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max);
    
    
    
    
    TH1F* histo_track_mean = new TH1F("histo_track_mean",";Track Mean [rad]; entries [#]",250,0.817,0.8348);
    TH1F* histo_track_spr = new TH1F("histo_track_spr",";Track SPR [m rad]; entries [#]",250,5.1,20);
    
    const int nbin_yield =100; // 100
    TH1F*  histo_track_yield = new TH1F("histo_track_yield",";Photon Yield; entries [#]",nbin_yield ,0,100);
    
    TH1F*  histo_track_resolution_bin[nbin_yield];
    TH1F*  histo_track_spr_bin[nbin_yield];
    TH1F*  histo_track_mean_bin[nbin_yield];
    
    for(Int_t bin=0; bin<=nbin_yield+1; bin++) {
        histo_track_resolution_bin[bin] = new TH1F(Form("histo_track_resolution_%d",bin), Form("Cherenkov track resolution @ yield bin %d ;Expected - measured [m rad]; Entries [#]",bin) , 100, -50, 50 );
        histo_track_spr_bin[bin] = new TH1F(Form("histo_spr_resolution_%d",bin), Form("SPR @ yield bin %d ;SPR [m rad];  [#]",bin) ,250,5.1,20);
        histo_track_mean_bin[bin] = new TH1F(Form("histo_mean_resolution_%d",bin), Form("#theta_{c}^{tr} @ yield bin %d ;#theta_{c}^{tr}  [rad]; Entries [#]",bin), 100,0.817,0.8348);
        
    }
    
    TH1F * histo_track_pos_resolution_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_spr_bin[pos_bin][pos_bin];
    //TH1F * histo_track_pos_mean_bin[pos_bin][pos_bin];
    
    for(Int_t x=0; x<=pos_bin+1; x++) {
        for(Int_t y=0; y<=pos_bin+1; y++) {
            
            histo_track_pos_resolution_bin[x][y] = new TH1F(Form("histo_track_resolution_xbin_%d_ybin_%d",x,y), Form("Cherenkov track resolution @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 100, -50, 50 );
            //histo_track_pos_spr_bin[x][y] = new TH1F(Form("histo_track_pos_spr_xbin_%d_ybin_%d",x,y), Form("SPR @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) , 250,5.1,20);
            //histo_track_pos_mean_bin[x][y] = new TH1F(Form("histo_track_pos_mean_xbin_%d_ybin_%d",x,y), Form("#theta_{c}^{tr} @ xbin %d ybin %d; X [cm]; Y [cm]",x,y) ,  100,0.817,0.8348);
        }
    }
    
    
    // variables
    double diff(-1);
    
    double mean_max(0.834), mean_min (0.818) ; // 0.817,0.8348);
    double spr_max(11.5), spr_min(6); // 11.5 6
    double calc_trk_res(-1);
    
    int x_pos_bin(-1),y_pos_bin(-1);
    double content_histo_pos_xy(-1),content_histo_pos_xy_tmp(-1),average_bin(-1);
    double track_resolution(-1),track_resolution_error(-1), yield_BinCenter(-1);
    double track_pos_resolution(-1),track_pos_resolution_error(-1), pos_BinCenter(-1);
    
    double track_spr_bin(-1),track_spr_error(-1);
    double track_mean_bin(-1),track_mean_error(-1);
    double content_histo_pos_xy_reso_tmp(-1);
    
    
    
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
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        if(track_pid !=2 ) continue; // select pion !=2
        
        //if(track_fit_chisqu>100)continue;
        //if(track_yield<10)continue;
        
        histo_track_mean->Fill(track_mean);
        histo_track_spr->Fill(track_spr*1000);
        
        if(track_mom> 4.5   || track_mom<3.5) continue;
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        
        //if(track_nbar<4 || track_nbar>8 ) continue;
        //if(track_xbar>10 || track_xbar < -10 ) continue;
        
        
        //diff = fAnglePi-track_mean;
        diff = track_mean-   0.82608;
        
        //std::cout<<"couter "<<couter<<"   "<<"track_yield "<<track_yield<<"   "<<"diff"<<diff*1000<<std::endl;
        
        histo_track_yield->Fill(track_yield);
        int xbin_yield = histo_track_yield->GetXaxis()->FindBin(track_yield);
        //cout<<xbin_yield<<endl;
        histo_track_resolution_bin[xbin_yield]->Fill(diff*1000);
        histo_track_spr_bin[xbin_yield]->Fill(track_spr*1000);
        histo_track_mean_bin[xbin_yield]->Fill(track_mean);
        
        
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
        
        //Resolution, spr, mean maps
        int xbin_pos = histo_pos_xy->GetXaxis()->FindBin(track_xbar);
        int ybin_pos = histo_pos_xy->GetYaxis()->FindBin(track_ybar);
        
        histo_track_pos_resolution_bin[xbin_pos][ybin_pos]->Fill(diff*1000);
        //histo_track_pos_spr_bin[xbin_pos][ybin_pos]->Fill(track_spr*1000);
        //histo_track_pos_mean_bin[xbin_pos][ybin_pos]->Fill(track_mean);
    }
    
    
    
    TCanvas *cc = new TCanvas("cc","cc",800,500);
    
    // fitting functions
    TF1 *fit_track_resolution = new TF1("fit_track_resolution","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-50,50);
    fit_track_resolution->SetLineColor(kBlack);
    fit_track_resolution->SetParameters(100,0,2);
    fit_track_resolution->SetParNames("p0","mean","resolution");
    fit_track_resolution->SetParLimits(0,0.1,1E6);
    fit_track_resolution->SetParLimits(1,-1,1);
    fit_track_resolution->SetParLimits(2,1,5);
    
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
    
    TF1 *fit_trk_reso = new TF1("fit_trk_reso","[0]*sqrt(([1]/sqrt(x))^2+[2])",0,100);
    fit_trk_reso->SetLineColor(kRed);
    fit_trk_reso->SetParameters(100,9,2);
    fit_trk_reso->SetParNames("p0","SPR","Reso");
    fit_trk_reso->SetParLimits(0,0.1,1E6);
    fit_trk_reso->SetParLimits(1,8,12);
    fit_trk_reso->SetParLimits(2,1,5);
    
    
    int couter(0);
    for (int i=0;i<nbin_yield;i++){
        if (histo_track_resolution_bin[i]->GetEntries() <200)continue; //400 //175 // 200
        histo_track_resolution_bin[i]->Fit("fit_track_resolution","M","", -50, 50) ;
        histo_track_spr_bin[i]->Fit("fit_track_spr","M","", 0, 30) ;
        histo_track_mean_bin[i]->Fit("fit_track_mean","M","", 0, 30) ;
        if(false){
            cc->cd();
            cc->Update();
            histo_track_resolution_bin[i]->Draw();
            cc->Update();
            cc->WaitPrimitive();
        }
        
        if(false){
            cc->cd();
            cc->Update();
            histo_track_spr_bin[i]->Draw();
            cc->Update();
            cc->WaitPrimitive();
        }
        
        if(false){
            cc->cd();
            cc->Update();
            histo_track_mean_bin[i]->Draw();
            cc->Update();
            cc->WaitPrimitive();
        }
        
        yield_BinCenter = histo_track_yield->GetXaxis()->GetBinCenter(i);
        
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
        
        graph_reso->SetPoint(couter, yield_BinCenter, track_resolution);
        graph_reso->SetPointError(couter, 1/2, 1/2,track_resolution_error/2,track_resolution_error/2);
        
        graph_spr->SetPoint(couter, yield_BinCenter, track_spr_bin);
        graph_spr->SetPointError(couter, 1/2, 1/2,track_spr_error/2,track_spr_error/2);
        
        graph_mean->SetPoint(couter, yield_BinCenter, track_mean_bin);
        graph_mean->SetPointError(couter, 1/2, 1/2,track_mean_error/2,track_mean_error/2);
        
        /////////
        g_pi->SetPoint(couter, yield_BinCenter, fAnglePi);
        g_k->SetPoint(couter, yield_BinCenter, fAngleK);
        /////////
        
        calc_trk_res=sqrt(track_resolution*track_resolution - (track_spr_bin/sqrt(yield_BinCenter))* (track_spr_bin/sqrt(yield_BinCenter)));
        g_calc_trk_res->SetPoint(couter, yield_BinCenter, calc_trk_res);
        
        ++couter;
    }
    
    
    
    if(true){
        // histograms
        glx_canvasAdd("r_resolution_bin",800,400);
        TMultiGraph *mg = new TMultiGraph();
        graph_reso->Fit("fit_trk_reso","M","", 0, 100) ;
        //graph_reso->Fit("fit_trk_reso","M","", 015, 50) ;
        mg->Add(graph_reso);
        //mg->Add(g_calc_trk_res);
        mg->SetTitle(" Cherenkov Resolution per Track ; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
        mg->Draw("APL");
        
        glx_canvasAdd("r_pos",800,400);
        histo_pos_xy->Draw("colz");
        
        glx_canvasAdd("r_pos_yield",800,400);
        histo_pos_xy_yield->Draw("colz");
        
        glx_canvasAdd("r_pos_spr",800,400);
        histo_pos_xy_spr->Draw("colz");
        
        /////////
        glx_canvasAdd("r_mean_bin",800,400);
        TMultiGraph *mg3 = new TMultiGraph();
        mg3->Add(graph_mean);
        mg3->Add(g_pi);
        mg3->Add(g_k);
        mg3->SetTitle(" #theta_{c}^{tr} per Track ; Photon Yield [#]; #theta_{c}^{tr} [rad]");
        mg3->Draw("APL");
        mg3->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.84);
        glx_canvasGet("r_mean_bin")->Update();
        TLine *line_pi= new TLine(0,0,0,1000);
        line_pi->SetY1(fAnglePi);
        line_pi->SetY2(fAnglePi);
        line_pi->SetX1(gPad->GetUymin());
        line_pi->SetX2(gPad->GetUymax());
        line_pi->SetLineColor(kBlack);
        line_pi->Draw();
        TLine *line_k= new TLine(0,0,0,1000);
        line_k->SetY1(fAngleK);
        line_k->SetY2(fAngleK);
        line_k->SetX1(gPad->GetUymin());
        line_k->SetX2(gPad->GetUymax());
        line_k->SetLineColor(kBlack);
        line_k->Draw();
        glx_canvasGet("r_mean_bin")->Update();
        /////////////
        glx_canvasAdd("r_spr_bin",800,400);
        TMultiGraph *mg2 = new TMultiGraph();
        mg2->Add(graph_spr);
        mg2->SetTitle(" SPR per Track ; Photon Yield [#]; SPR [m rad]");
        mg2->Draw("APL");
        mg2->GetHistogram()->GetYaxis()->SetRangeUser(0,15);
        line_k->Draw();// to fix bug on canvas update
        glx_canvasGet("r_spr_bin")->Update();
        /////////////
        
        glx_canvasAdd("r_yield",800,400);
        histo_track_yield->Draw();
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
        
        
    }
    

    // resolution map
    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
    int couter2(0);
    for (int x=0;x<pos_bin;x++){
        for (int y=0;y<pos_bin;y++){
            double hentry =histo_track_pos_resolution_bin[x][y]->GetEntries();
            if (hentry <100)continue;
            histo_track_pos_resolution_bin[x][y]->Fit("fit_track_resolution","M","", -50, 50) ;
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
    histo_pos_xy_reso->Draw("colz");
    
    
    
    // delete histograms
    for (int x=0;x<pos_bin;x++){
        for (int y=0;y<pos_bin;y++){
            delete histo_track_pos_resolution_bin[x][y];
        }
    }

    
    cout<<"####### fAngleK "<< fAngleK<<"  fAnglePi "<<fAnglePi<<endl;
    
}
