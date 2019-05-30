#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "glxtools.C"
void analyses(){

    // calculate cherenkov angle
    Double_t momentum=3.5;
    Int_t pdg[]= {11,13,211,321,2212};
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.009),range(5*sigma),noise(0.3);
    Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
    Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00;
    
    // graphs
    TGraph *g_pi = new TGraph();
    g_pi->SetMarkerColor(kBlue);
    g_pi->SetMarkerStyle(20);
    g_pi->SetLineColor(kBlue);
    g_pi->SetLineWidth(1);
    
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
    
    // histograms
    TH1F* histo_track_mean = new TH1F("histo_track_mean",";track mean [rad]; entries [#]",250,0.817,0.8348);
    TH1F* histo_track_spr = new TH1F("histo_track_spr",";track SPR [m rad]; entries [#]",250,5.1,20);
    
    const int nbin_yield =100; // 100
    TH1F*  histo_track_yield = new TH1F("histo_track_yield",";phton yield; entries [#]",nbin_yield ,0,100);
    
    TH1F*  histo_track_resolution_bin[nbin_yield];
    TH1F*  histo_track_spr_bin[nbin_yield];
    TH1F*  histo_track_mean_bin[nbin_yield];
    
    for(Int_t bin=0; bin<=nbin_yield+1; bin++) {
        histo_track_resolution_bin[bin] = new TH1F(Form("histo_track_resolution_%d",bin), Form("Cherenkov track resolution @ yield bin %d ;Expected - measured [m rad]; Entries [#]",bin) , 100, -50, 50 );
        histo_track_spr_bin[bin] = new TH1F(Form("histo_spr_resolution_%d",bin), Form("SPR @ yield bin %d ;SPR [m rad];  [#]",bin) ,250,5.1,20);
        histo_track_mean_bin[bin] = new TH1F(Form("histo_mean_resolution_%d",bin), Form("#theta_{c}^{tr} @ yield bin %d ;#theta_{c}^{tr}  [rad]; Entries [#]",bin), 100,0.817,0.8348);

    }
    
    // variables
    double diff(-1);
    
    double mean_max(0.834), mean_min (0.818) ; // 0.817,0.8348);
    double spr_max(11.5), spr_min(6);


    
    // read tree
    TFile *f = new TFile("outFile.root");
    TTree *tree_variables = (TTree*)f->Get("tree_variables");
    double track_spr(-1),track_mean(-1), track_yield(-1), track_mom(-1), track_xbar(0),track_ybar(0);
    int track_pid(-1), track_nbar(-1);
    
    tree_variables->SetBranchAddress("track_pid",&track_pid);
    tree_variables->SetBranchAddress("track_spr",&track_spr);
    tree_variables->SetBranchAddress("track_mean",&track_mean);
    tree_variables->SetBranchAddress("track_yield",&track_yield);
    tree_variables->SetBranchAddress("track_mom",&track_mom);
    tree_variables->SetBranchAddress("track_xbar",&track_xbar);
    tree_variables->SetBranchAddress("track_ybar",&track_ybar);
    tree_variables->SetBranchAddress("track_nbar",&track_nbar);
    
    Long64_t nentries = tree_variables->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        if(track_pid !=2 ) continue; // select pion !=2
        
        histo_track_mean->Fill(track_mean);
        histo_track_spr->Fill(track_spr*1000);
        
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        //if(track_mom> 4.5   || track_mom<3.5) continue;
        
        
        
        //diff = fAnglePi-track_mean;
        diff = track_mean-   0.82608;
        
        //std::cout<<"couter "<<couter<<"   "<<"track_yield "<<track_yield<<"   "<<"diff"<<diff*1000<<std::endl;
        
        histo_track_yield->Fill(track_yield);
        int xbin_yield = histo_track_yield->GetXaxis()->FindBin(track_yield);
        //cout<<xbin_yield<<endl;
        histo_track_resolution_bin[xbin_yield]->Fill(diff*1000);
        histo_track_spr_bin[xbin_yield]->Fill(track_spr*1000);
        histo_track_mean_bin[xbin_yield]->Fill(track_mean);

    }
    

    
    TCanvas *cc = new TCanvas("cc","cc",800,500);
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
    
    double track_resolution(-1),track_resolution_error(-1), yield_BinCenter(-1);
    double track_spr_bin(-1),track_spr_error(-1);
    double track_mean_bin(-1),track_mean_error(-1);
    
    double width =histo_track_yield->GetBinWidth(1);
    
    int couter(0);
    for (int i=0;i<nbin_yield;i++){
        if (histo_track_resolution_bin[i]->GetEntries() <200)continue; //400 //175
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
        track_mean_error= fit_track_mean->GetParError(1);

        //track_spr_bin= fit_track_spr->GetParameter(2);
        //track_spr_error= fit_track_spr->GetParError(2);
        
        //track_spr_bin= fit_track_spr->GetParameter(1);
        //track_spr_error= fit_track_spr->GetParError(1);
        
        //track_spr_bin= histo_track_spr_bin[i]->GetStdDev();
        //track_spr_error= histo_track_spr_bin[i]->GetStdDevError();
        
        track_spr_bin= histo_track_spr_bin[i]->GetMean();
        track_spr_error= histo_track_spr_bin[i]->GetMeanError();

        graph_reso->SetPoint(couter, yield_BinCenter, track_resolution);
        graph_reso->SetPointError(couter, width/2, width/2,track_resolution_error/2,track_resolution_error/2);
        
        graph_spr->SetPoint(couter, yield_BinCenter, track_spr_bin);
        graph_spr->SetPointError(couter, width/2, width/2,track_spr_error/2,track_spr_error/2);
        
        graph_mean->SetPoint(couter, yield_BinCenter, track_mean_bin);
        graph_mean->SetPointError(couter, width/2, width/2,track_mean_error/2,track_mean_error/2);
        
        ++couter;
    }
    
    glx_canvasAdd("r_resolution_bin",800,400);
    TMultiGraph *mg = new TMultiGraph();
    //mg->Add(g_pi);
    mg->Add(graph_reso);
    mg->SetTitle(" Cherenkov Resolution per Track ; Photon Yield [#]; #sigma( #theta_{c}^{tr} ) [m rad]");
    mg->Draw("APL");
    //mg->GetHistogram()->GetYaxis()->SetRangeUser(6800,7050);
    
    glx_canvasAdd("r_spr_bin",800,400);
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(graph_spr);
    mg2->SetTitle(" SPR per Track ; Photon Yield [#]; SPR [m rad]");
    mg2->Draw("APL");

    glx_canvasAdd("r_mean_bin",800,400);
    TMultiGraph *mg3 = new TMultiGraph();
    mg3->Add(graph_mean);
    mg3->SetTitle(" #theta_{c}^{tr} per Track ; Photon Yield [#]; #theta_{c}^{tr} [rad]");
    mg3->Draw("APL");
    
    

    
    if(false){
        glx_canvasAdd("r_yield",800,400);
        histo_track_yield->Draw();
        
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
    
}
