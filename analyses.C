#include "TMultiGraph.h"
#include "TGraph.h"
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
    
    // histograms
    
    TH1F* histo_track_mean = new TH1F("histo_track_mean",";track mean [rad]; entries [#]",250,0.80,0.84);
    TH1F* histo_track_spr = new TH1F("histo_track_spr",";track SPR [m rad]; entries [#]",250,5.1,20);
    
    int nbin_yield =100; // 100
    TH1F*  histo_track_yield = new TH1F("histo_track_yield",";phton yield; entries [#]",nbin_yield ,0,100);
    
    TH1F*  histo_track_resolution[nbin_yield];
    for(Int_t bin=0; bin<=nbin_yield+1; bin++) {
        histo_track_resolution[bin] = new TH1F(Form("histo_track_resolution_%d",bin), Form("Cherenkov track resolution @ yield bin %d ;Expected - measured [m rad]; count [#]",bin) , 100, -50, 50 );
    }
    
    
    // variables
    Double_t track_resolution(-1);
    
    double mean_max(0.834), mean_min (0.818) ;
    double spr_max(11), spr_min(6);
    
    
    
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
    
    int couter(0);
    Long64_t nentries = tree_variables->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
        tree_variables->GetEntry(i);
        
        if(track_pid != 2 ) continue; // select pion
        
        histo_track_mean->Fill(track_mean);
        histo_track_spr->Fill(track_spr*1000);
        
        if(track_mean> mean_max   || track_mean<mean_min) continue;
        if(track_spr*1000> spr_max || track_spr*1000<spr_min ) continue;
        
        
        //track_resolution = fAnglePi-track_mean;
        track_resolution = track_mean-   0.82608;
        
        //std::cout<<"couter "<<couter<<"   "<<"track_yield "<<track_yield<<"   "<<"track_resolution"<<track_resolution*1000<<std::endl;
        
        histo_track_yield->Fill(track_yield);
        int xbin_yield = histo_track_yield->GetXaxis()->FindBin(track_yield);
        //cout<<xbin_yield<<endl;
        histo_track_resolution[xbin_yield]->Fill(track_resolution*1000);
        
        
        g_pi->SetPoint(couter, track_yield, track_resolution*1000);
        ++couter;
    }
    
    TCanvas *cc = new TCanvas("cc","cc",800,500);
    TF1 *fit_track_resolution = new TF1("fit_track_resolution","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-50,50);
    fit_track_resolution->SetLineColor(kBlack);
    fit_track_resolution->SetParameters(100,0,2);
    fit_track_resolution->SetParNames("p0","mean","cherenkov track resolution");
    fit_track_resolution->SetParLimits(0,0.1,1E6);
    fit_track_resolution->SetParLimits(1,-1,1);
    fit_track_resolution->SetParLimits(2,1,5);
    
    for (int i=0;i<nbin_yield;i++){

        
        
        histo_track_resolution[i]->Fit("fit_track_resolution","M","", -50, 50) ;
        cc->cd();
        cc->Update();
        histo_track_resolution[i]->Draw();
        cc->Update();
        cc->WaitPrimitive();
    }
    
    
    glx_canvasAdd("r_yield",800,400);
    histo_track_yield->Draw();
    
    
    
    if(false){
        
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
        
        
        
        glx_canvasAdd("r_graph1",800,400);
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(g_pi);
        mg->SetTitle(" yield ; Yield [#]; SPR [m rad]");
        mg->Draw("AP");
        //mg->GetHistogram()->GetYaxis()->SetRangeUser(6800,7050);
        g_pi->Draw("AP");
    }
    
}
