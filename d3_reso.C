
{
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);
    
     int pos_bin(100), pos_min(-100), pos_max(100);

    
    //TH3F * pol = new TH3F( "pol" , "; Bar Hit X ; Bar Hit Y (cm)", pos_bin, pos_min, pos_max, pos_bin, pos_min, pos_max,100,0,10);
    
    
    TFile f("reso_pos_mom.root");
    TNtuple *ntuple = (TNtuple*)f.Get("ntuple");
    
    TCanvas *c1 = new TCanvas("c1","c1",800,400);
    ntuple->SetMarkerStyle(20);
    ntuple->SetMarkerSize(1.5);
    
 
    ntuple->Draw("ypos:xpos:mom_pos:reso_pos >> pol","","COLZ");
    

    pol->GetZaxis()->SetTitle("Bar Hit Y [cm]");
    pol->GetYaxis()->SetTitle("Bar Hit X [cm]");
    pol->GetXaxis()->SetTitle("Momentum GeV/c");
    pol->GetXaxis()->SetTitleOffset(1.5);
    pol->GetYaxis()->SetTitleOffset(1.8);
    pol->GetZaxis()->SetTitleOffset(1.5);
    pol->GetXaxis()->CenterTitle(true);
    pol->GetYaxis()->CenterTitle(true);
    pol->GetZaxis()->CenterTitle(true);
    
    pol->SetTitle("Averaged Cherenkov Track Resolution 4D");

    
    c1->Update();
    gPad->Update();
    

    
    
    
}
