#include "TImage.h"
#include "prttools/prttools.C"

bool exists_test (const std::string& name);
void fNorm(TH1I *p_diff_time_sim, TH1I *p_diff_time_data);

void gluex_plots(){
    
     prt_savepath="gx";
    
    gStyle->SetPalette( kRainBow);
    
    Double_t momentum=4.0;
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
    Double_t fAngleK= acos(sqrt(momentum*momentum + mass[3]*mass[3])/momentum/1.4738)-0.00;
    Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00;
    
    
    cout<<"fAngleP "<< fAngleP<< endl;
     cout<<"fAnglePi "<< fAnglePi<< endl;
     cout<<"diff  "<< fAnglePi-fAngleP<< endl;
    
    TFile *ffile_pi,  *ffile_k;
    TH1I *hist_chere_pi, *hist_nph_pi, *hist_loglikelihood_pi;
    TH1I *hist_chere_k , *hist_nph_k , *hist_loglikelihood_k;
    
    Int_t nf= 41;
    Int_t kth_chere(0), kphi_chere(0);
    Int_t counter_mean_loop(0);
    Double_t separation(0);
    Double_t m1,m2,s1,s2;
    
    TCanvas* c = new TCanvas("c","c",0,0,800,1200);
    
    TH2F * hist_cherenkove_map_nph_pi =  new TH2F("hist_cherenkove_map_nph_pi",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    TH2F * hist_cherenkove_map_chere_pi =  new TH2F("hist_cherenkove_map_chere_pi",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    TH2F * hist_cherenkove_map_spr_pi =  new TH2F("hist_cherenkove_map_spr_pi",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    
    TH2F * hist_cherenkove_map_nph_k =  new TH2F("hist_cherenkove_map_nph_k",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    TH2F * hist_cherenkove_map_chere_k =  new TH2F("hist_cherenkove_map_chere_k",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    TH2F * hist_cherenkove_map_spr_k =  new TH2F("hist_cherenkove_map_spr_k",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    
    TH2F * hist_pi_k_separation =  new TH2F("hist_pi_k_separation",";#Theta[Degree]; #Phi[Degree]", nf, 0, 12, nf, -180 , 20);
    
    
    TF1 *ff;
    
    for (int i=0; i<=120; i+=3) {
        for (int j=20; j>=-180; j-=5) {
            
            if ( i > 50 || i < 30) continue;
            
             if (  j != -80) continue;
            
            //////////////////////
            // Loop over files //
            /////////////////////
            
            Double_t iii= (Float_t) i/10;
            
            Double_t ii = (Float_t) i/10.0;
            TString ii_string = Form("th%.1f_", ii);
            if(i%10==0)
            {
                Int_t ii=i/10;
                //std::cout<<"ii "<<ii<<std::endl;
                ii_string = Form("th%.d_", ii);
                //std::cout<<"ii_string "<<ii_string<<std::endl;
            }
            TString jj_string = Form("phi%d", j);
            if(j==0) jj_string = "phi0";
            if(i==0) ii_string = "th0_";
            TString MCPiSamplePath = "/Users/ahmed/Desktop/gluex/MCPiSample/mom4_pdg2_"+ii_string+jj_string+"/root/hd_root_particle_gun_060000_000.root";
            TString MCKSamplePath = "/Users/ahmed/Desktop/gluex/MCKSample/mom4_pdg3_"+ii_string+jj_string+"/root/hd_root_particle_gun_060000_000.root";
            //cout<<"MCPiSamplePath= " <<MCPiSamplePath<<endl;
            string path_data_Pi = (string)MCPiSamplePath;
            string path_data_K = (string)MCKSamplePath;
            
            //////////////////////////
            // check file existance //
            //////////////////////////
            
            //cout<<"exists_test(path_data_Pi)" <<exists_test(path_data_Pi)<<endl;
            if (!exists_test(path_data_Pi)){
                //std::cout<<"jj_string "<<j<<std::endl;
                cout<<"Not found MCPiSamplePath= " <<MCPiSamplePath<<endl;
            }
            
            if (!exists_test(path_data_K)){
                //std::cout<<"jj_string "<<j<<std::endl;
                cout<<"Not found MCKSamplePath= " <<MCKSamplePath<<endl;
            }
            
            if (!exists_test(path_data_Pi)) continue;
            if (!exists_test(path_data_K)) continue;
            
            //continue;
            
            //////////////////////////
            // Access histograms    //
            //////////////////////////
            
            ffile_pi  = new TFile(MCPiSamplePath, "READ");
            ffile_k  = new TFile(MCKSamplePath, "READ");
            
            hist_chere_pi=(TH1I*)ffile_pi->Get("DIRC/Pi+/hThetaC_Pi+");
            hist_nph_pi=(TH1I*)ffile_pi->Get("DIRC/Pi+/hNphC_Pi+");
            hist_loglikelihood_pi=(TH1I*)ffile_pi->Get("DIRC/Pi+/hLikelihoodDiff_Pi+");
            
            hist_chere_k=(TH1I*)ffile_k->Get("DIRC/K+/hThetaC_K+");
            hist_nph_k=(TH1I*)ffile_k->Get("DIRC/K+/hNphC_K+");
            hist_loglikelihood_k=(TH1I*)ffile_k->Get("DIRC/K+/hLikelihoodDiff_K+");
            
            //if (counter_mean_loop >= 5 ) break;
            
            Double_t chere_mean =hist_chere_pi->GetMean() ;
            Double_t mean_nph =hist_nph_pi->GetMean() ;
            Double_t mean_nph_k =hist_nph_k->GetMean() ;
            //cout<<"counter_mean_loop " << counter_mean_loop<<"  "<< chere_mean <<" "<<iii<<" "<<j<<endl;
            if (hist_chere_pi->GetEntries()== 0){
                //cout<<"####### Warnning empty histogram #####"<<endl;
                //cout<<"MCPiSamplePath= " <<MCPiSamplePath<<endl;
                //chere_mean = 1;
            }
            
            kth_chere  = hist_cherenkove_map_chere_pi->GetXaxis()->FindBin(iii);
            kphi_chere = hist_cherenkove_map_chere_pi->GetYaxis()->FindBin(j);
            //cout<<"####### 1st bin "<<kth_chere<<"  2nd bin= "<< kphi_chere<<" val "<<chere_mean<<"         "<< iii<<""<<j<<endl;
            //cout<<"####### 1st bin "<< iii <<"  "<<kth_chere<<endl;
            
            /////////////////////////////////
            // Fit cherenkov histograms    //
            /////////////////////////////////
            
            TF1 *fFit = new TF1("fFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            TSpectrum *fSpect= new TSpectrum(10);
            Double_t cangle=0;
            Double_t spr=0;
            //gROOT->SetBatch(1);
            Int_t nfound = fSpect->Search(hist_chere_pi,1,"",0.9); //0.6
            if(nfound>0) cangle = fSpect->GetPositionX()[0];
            cangle =  hist_chere_pi->GetXaxis()->GetBinCenter(hist_chere_pi->GetMaximumBin());
            if(cangle>0.85) cangle=0.82;
            fFit->SetParameters(100,cangle,0.010);
            fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit->SetParLimits(0,0.1,1E6);
            fFit->SetParLimits(1,cangle-0.04,cangle+0.04);
            fFit->SetParLimits(2,0.005,0.018); // changed 0.014
            hist_chere_pi->Fit("fFit","M+","",cangle-0.06,cangle+0.06);
            //p_cherenkov_data_copy->Fit("fFit","0","",cangle-0.06,cangle+0.06);
            //hist_chere_pi->Fit("fFit","R");
            Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
            cangle = fFit->GetParameter(1);
            spr = fFit->GetParameter(2);
            Double_t cangle_minus_5_sgma = cangle-5*spr;
            Double_t cangle_plus_5_sgma = cangle+5*spr;
            Double_t cangle_minus_3_sgma = cangle-3*spr;
            Double_t cangle_plus_3_sgma = cangle+3*spr;
            Double_t r_min = cangle-8*spr;
            Double_t r_max = cangle+8*spr;
            Double_t sumundercurve = fFit->Integral(cangle_minus_3_sgma,cangle_plus_3_sgma);
            
            Double_t cangle_k=0;
            Double_t spr_k=0;
            hist_chere_k->Fit("fFit","M+","",cangle-0.06,cangle+0.06);
            Double_t chi_k = fFit->GetChisquare()/fFit->GetNDF();
            cangle_k = fFit->GetParameter(1);
            spr_k = fFit->GetParameter(2);
            
            //////////////////
            // Fill Maps    //
            //////////////////
            
            if (hist_nph_pi->GetEntries()>30){
                
                hist_cherenkove_map_nph_pi->Fill(iii,j, mean_nph);
                hist_cherenkove_map_chere_pi->Fill(iii,j, cangle-fAnglePi);
                hist_cherenkove_map_spr_pi->Fill(iii,j, spr*1000);
                
                hist_cherenkove_map_nph_k->Fill(iii,j, mean_nph_k);
                hist_cherenkove_map_chere_k->Fill(iii,j, cangle_k-fAngleK);
                hist_cherenkove_map_spr_k->Fill(iii,j, spr_k*1000);
                
                fNorm(hist_loglikelihood_k,hist_loglikelihood_pi);
                hist_loglikelihood_k->SetLineColor(2);
                
                
                if(hist_loglikelihood_k->GetEntries()>10){
                    hist_loglikelihood_k->Fit("gaus","S");
                    ff = hist_loglikelihood_k->GetFunction("gaus");
                    m1=ff->GetParameter(1);
                    s1=ff->GetParameter(2);
                }
                if(hist_loglikelihood_pi->GetEntries()>10){
                    hist_loglikelihood_pi->Fit("gaus","S");
                    ff = hist_loglikelihood_pi->GetFunction("gaus");
                    m2=ff->GetParameter(1);
                    s2=ff->GetParameter(2);
                }
                separation = (fabs(m2-m1))/(0.5*(s1+s2));
                
                
                

                
               
                
                
//                                c->cd();
//                                hist_loglikelihood_k->SetLineColor(kRed);
//                                hist_loglikelihood_pi->Draw();
//                                hist_loglikelihood_k->Draw("same");
//                                c->Update();
//                                c->WaitPrimitive();
                
                std::cout<<"separation "<< separation <<std::endl;
                hist_pi_k_separation->Fill(iii,j, separation);
                
//                                c->cd();
//                                hist_chere_pi->SetLineColor(kRed);
//                                hist_chere_k->Draw();
//                                hist_chere_pi->Draw("same");
//                                c->Update();
//                                c->WaitPrimitive();
                
                
                
                
//                hist_chere_k->SetLineColor(4);
//                hist_chere_pi->SetLineColor(2);
//                hist_chere_k->SetMarkerColor( kBlue+1 );
//                hist_chere_pi->SetMarkerColor(kRed+1);
//                hist_chere_k->SetLineColor(4);
//                hist_chere_pi->SetLineColor(2);
//
//
//                c->cd();
//                hist_chere_pi->SetLineColor(kRed);
//                hist_chere_k->Draw();
//                //hist_chere_pi->Draw("same");
//
//                c->Update();
//                TLine *lin_ch_p_v = new TLine(0,0,0,1000);
//                lin_ch_p_v->SetX1(fAngleP);
//                lin_ch_p_v->SetX2(fAngleP);
//                lin_ch_p_v->SetY1(gPad->GetUymin());
//                lin_ch_p_v->SetY2(gPad->GetUymax());
//                lin_ch_p_v->SetLineColor(kRed);
//                TLine *lin_ch_pi_v = new TLine(0,0,0,1000);
//                lin_ch_pi_v->SetX1(fAnglePi);
//                lin_ch_pi_v->SetX2(fAnglePi);
//                lin_ch_pi_v->SetY1(gPad->GetUymin());
//                lin_ch_pi_v->SetY2(gPad->GetUymax());
//                lin_ch_pi_v->SetLineColor(kBlue);
//                lin_ch_p_v->Draw();
//                lin_ch_pi_v->Draw();
//                c->Update();
//
//                c->WaitPrimitive();
                
                
                                                c->cd();
                                                hist_nph_k->SetLineColor(kBlue+1);
                                                hist_nph_k->Draw();
                                                //hist_chere_pi->Draw("same");
                                                c->Update();
                                                c->WaitPrimitive();
                
                
            }else{
              
                hist_cherenkove_map_chere_pi->Fill(iii,j, -10);
                hist_cherenkove_map_chere_k->Fill(iii,j, -10);
                
            }
            
            
            /////////////
            // Close   //
            /////////////
            
            ffile_pi->Close();
            delete ffile_pi;
            ++counter_mean_loop;
        }
    }
    //Pion plots
    prt_canvasAdd("r_cherenkov_reco_map_pi ",800,400);
    hist_cherenkove_map_chere_pi->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_chere_pi->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_chere_pi->SetStats(0);
    hist_cherenkove_map_chere_pi-> SetTitle("Reco Cherenkov angle - Expected (GlueX DIRC)");
    hist_cherenkove_map_chere_pi->SetMinimum(- 0.01);
    hist_cherenkove_map_chere_pi->Draw("colz SPH");
    
    prt_canvasAdd("r_cherenkov_spr_map_pi ",800,400);
    hist_cherenkove_map_spr_pi->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_spr_pi->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_spr_pi->SetStats(0);
    hist_cherenkove_map_spr_pi-> SetTitle("Singel photon resolution (GlueX DIRC)");
    hist_cherenkove_map_spr_pi->Draw("colz");
    
    prt_canvasAdd("r_cherenkov_mean_map_pi ",800,400);
    hist_cherenkove_map_nph_pi->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_nph_pi->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_nph_pi->SetStats(0);
    hist_cherenkove_map_nph_pi-> SetTitle("Photon yield (GlueX DIRC)");
    hist_cherenkove_map_nph_pi->Draw("colz");
    
    // Kaon plots
    prt_canvasAdd("r_cherenkov_reco_map_k ",800,400);
    hist_cherenkove_map_chere_k->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_chere_k->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_chere_k->SetStats(0);
    hist_cherenkove_map_chere_k-> SetTitle("Reco Cherenkov angle - Expected (GlueX DIRC)");
    hist_cherenkove_map_chere_k->SetMinimum(- 0.01);
    hist_cherenkove_map_chere_k->Draw("colz SPH");
    
    prt_canvasAdd("r_cherenkov_spr_map_k ",800,400);
    hist_cherenkove_map_spr_k->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_spr_k->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_spr_k->SetStats(0);
    hist_cherenkove_map_spr_k-> SetTitle("Singel photon resolution (GlueX DIRC)");
    hist_cherenkove_map_spr_k->Draw("colz");
    
    prt_canvasAdd("r_cherenkov_mean_map_k ",800,400);
    hist_cherenkove_map_nph_k->GetXaxis()->SetRangeUser(0.,11.5);
    hist_cherenkove_map_nph_k->GetYaxis()->SetRangeUser(-170,10);
    hist_cherenkove_map_nph_k->SetStats(0);
    hist_cherenkove_map_nph_k-> SetTitle("Photon yield (GlueX DIRC)");
    hist_cherenkove_map_nph_k->Draw("colz");
    
    // sepration
    prt_canvasAdd("r_hist_pi_k_separation ",800,400);
    hist_pi_k_separation->GetXaxis()->SetRangeUser(0.,11.5);
    hist_pi_k_separation->GetYaxis()->SetRangeUser(-170,10);
    hist_pi_k_separation->SetStats(0);
    hist_pi_k_separation-> SetTitle("pi k separation (GlueX DIRC)");
    hist_pi_k_separation->Draw("colz");
    
    prt_canvasSave(2,0);
    prt_canvasDel("*");
    
}

//////////////////////////
// check file existance //
//////////////////////////

bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}

//////////////////////////
// Histo Normalization  //
//////////////////////////

void fNorm(TH1I *p_diff_time_sim, TH1I *p_diff_time_data) {
    p_diff_time_data->Scale(p_diff_time_sim->GetMaximum() /p_diff_time_data->GetMaximum());
}
