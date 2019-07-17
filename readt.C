
{
    
    gStyle->SetPalette(55);
    
    gStyle->SetOptStat(0);
    TH1F*  histo_time_bar_pos[24][40];
    TH1F*  histo_tdiffD_bar_pos[24][40];
    TH1F*  histo_tdiffR_bar_pos[24][40];
    
    TH1F*  histo_tD_bar_pos[24][40];
    TH1F*  histo_tR_bar_pos[24][40];
    
    TH1F*  histo_cD_bar_pos[24][40];
    TH1F*  histo_cR_bar_pos[24][40];
    TH1F*  histo_c_bar_pos[24][40];
    
    TString path ="/Users/ahmed/GlueX_DIRC_Calib/histo.root";//histo_with_time_correction.root
    cout<<"path= " <<path<<endl;
    TFile *f = new TFile(path, "READ");
    
    //    TH2F * histo_pos_xy_entries = new TH2F( "histo_pos_xy_entries" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
    //    TH2F * histo_xy_tdiff_meanD = new TH2F( "histo_xy_tdiff_meanD" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
    //    TH2F * histo_xy_tdiff_sigmaD = new TH2F( "histo_xy_tdiff_sigmaD" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
    //
    //    TH2F * histo_xy_tdiff_meanR = new TH2F( "histo_xy_tdiff_meanR" , "; Bar number ; Bar Hit X", 24, 0, 24, 40, 0, 40);
    //    TH2F * histo_xy_tdiff_sigmaR = new TH2F( "histo_xy_tdiff_sigmaR" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
    
    // maps
    TH2F * histo_pos_xy_entries = new TH2F( "histo_pos_xy_entries" , "; Bar Hit X [cm]; Bar number [#]  ", 40, 0, 40 ,24, 0, 24);
    
    TH2F * histo_xy_tdiff_meanD = new TH2F( "histo_xy_tdiff_meanD" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_tdiff_sigmaD = new TH2F( "histo_xy_tdiff_sigmaD" , "; Bar Hit X [cm]; Bar number [#] ",  40, 0, 40,24, 0, 24);
    TH2F * histo_xy_tdiff_meanR = new TH2F( "histo_xy_tdiff_meanR" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_tdiff_sigmaR = new TH2F( "histo_xy_tdiff_sigmaR" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    TH2F * histo_xy_t_meanD = new TH2F( "histo_xy_t_meanD" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_t_sigmaD = new TH2F( "histo_xy_t_sigmaD" , "; Bar Hit X [cm]; Bar number [#] ",  40, 0, 40,24, 0, 24);
    TH2F * histo_xy_t_meanR = new TH2F( "histo_xy_t_meanR" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_t_sigmaR = new TH2F( "histo_xy_t_sigmaR" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    TH2F * histo_xy_c_meanD = new TH2F( "histo_xy_c_meanD" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_sigmaD = new TH2F( "histo_xy_c_sigmaD" , "; Bar Hit X [cm]; Bar number [#] ",  40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_meanR = new TH2F( "histo_xy_c_meanR" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_sigmaR = new TH2F( "histo_xy_c_sigmaR" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    TH2F * histo_xy_c_mean = new TH2F( "histo_xy_c_mean" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_c_sigma = new TH2F( "histo_xy_c_sigma" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    
    //occupancy
    TH2F * histo_xy_c_occu = new TH2F( "histo_xy_c_occu" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_cD_occu = new TH2F( "histo_xy_cD_occu" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    TH2F * histo_xy_cR_occu = new TH2F( "histo_xy_cR_occu" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);

    
    
    histo_xy_tdiff_meanD->SetTitle("Time Difference Shift Direct photons");
    histo_xy_tdiff_meanR->SetTitle("Time Difference Shift Reflected photons");
    histo_xy_tdiff_sigmaD->SetTitle("Time Difference Sigma Direct photons");
    histo_xy_tdiff_sigmaR->SetTitle("Time Difference Sigma Reflected photons");
    
    histo_xy_t_meanD->SetTitle("Time Spectrum Direct Photons [ns]");
    histo_xy_t_meanR->SetTitle("Time Spectrum Reflected Photons [ns]");
    histo_xy_t_sigmaD->SetTitle("Time Spectrum Sigma Direct Photons [ns]");
    histo_xy_t_sigmaR->SetTitle("Time Spectrum Sigma Reflected Photons [ns]");
    
    histo_xy_c_meanD->SetTitle("#theta_{c} Shift Direct Photons [mrad]");
    histo_xy_c_meanR->SetTitle("#theta_{c} Shift Reflected Photons [mrad]");
    histo_xy_c_sigmaD->SetTitle("SPR Direct Photons [mrad]");
    histo_xy_c_sigmaR->SetTitle("SPR Reflected Photons [mrad]");
    
    histo_xy_c_mean->SetTitle("#theta_{c} Shift [mrad]");
    histo_xy_c_sigma->SetTitle("SPR [mrad]");
    
    
    histo_xy_c_occu->SetTitle("Occupancy");
    histo_xy_cD_occu->SetTitle("Occupancy Direct Photons");
    histo_xy_cR_occu->SetTitle("Occupancy Reflected Photons");
    
    
    
    
    // extra
    TH1F* histo_shiftD_diff = new TH1F("histo_shiftD_diff","; Mean per segment [ns]; entries [#]",100,-5,5);
    TH1F* histo_shiftR_diff = new TH1F("histo_shiftR_diff","; Mean per segment [ns]; entries [#]",100,-5,5);
    
    TH1F* histo_shiftD_sgma = new TH1F("histo_shiftD_sgma","; Sigma per segment [ns]; entries [#]",100,0,5);
    TH1F* histo_shiftR_sgma = new TH1F("histo_shiftR_sgma","; Sigma per segment [ns]; entries [#]",100,0,5);
    
    
    // fitting functions
    
    TF1 *fit_gause = new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-10,10);
    //TF1 *fit_gause =  new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-10,10);
    fit_gause->SetLineColor(kBlack);
    fit_gause->SetParameters(100,9,2);
    fit_gause->SetParNames("p0","mean ","sigma");
    fit_gause->SetParLimits(0,0.1,1E6);
    fit_gause->SetParLimits(1,-3,3);
    fit_gause->SetParLimits(2,0.5,3);
    
    
    TF1 *fit_cherenkov = new TF1("fit_cherenkov","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",-0.05,0.05);
    fit_cherenkov->SetNpx(1000);
    
    fit_cherenkov->SetLineColor(kBlack);
    fit_cherenkov->SetParameters(100,9,2);
    fit_cherenkov->SetParNames("p0","mean ","sigma");
    fit_cherenkov->SetParLimits(0,0.1,1E6);
    fit_cherenkov->SetParLimits(1,-0.05,0.05);
    fit_cherenkov->SetParLimits(2,0.005,0.014);
    
    //graph
    TGraph *g_c_meanD_time= new TGraph();
    TGraph *g_c_sigmaD_time= new TGraph();
    TGraph *g_c_meanR_time= new TGraph();
    TGraph *g_c_sigmaR_time= new TGraph();
    
    
    g_c_meanD_time->SetMarkerColor(kBlue);
    g_c_meanD_time->SetMarkerStyle(21);
    g_c_meanD_time->SetLineColor(kBlue);
    g_c_meanD_time->SetLineWidth(3);
    g_c_meanD_time->SetLineColor(kBlue);
    
    g_c_meanR_time->SetMarkerColor(kRed);
    g_c_meanR_time->SetMarkerStyle(21);
    g_c_meanR_time->SetLineColor(kRed);
    g_c_meanR_time->SetLineWidth(3);
    g_c_meanR_time->SetLineColor(kRed);
    
    g_c_sigmaD_time->SetMarkerColor(kBlue);
    g_c_sigmaD_time->SetMarkerStyle(21);
    g_c_sigmaD_time->SetLineColor(kBlue);
    g_c_sigmaD_time->SetLineWidth(3);
    g_c_sigmaD_time->SetLineColor(kBlue);
    
    g_c_sigmaR_time->SetMarkerColor(kRed);
    g_c_sigmaR_time->SetMarkerStyle(21);
    g_c_sigmaR_time->SetLineColor(kRed);
    g_c_sigmaR_time->SetLineWidth(3);
    g_c_sigmaR_time->SetLineColor(kRed);
    
    
    // arrays
    Double_t shiftD_array[24][40]={0};
    Double_t shiftR_array[24][40]={0};
    double cop[]= {62,62,62,62,62,62,62,62,60,58,56,54,53,52,52,51,50,49,47,46,46,45,45,44,44,43,43,42,41,41,40,40,38,37.5,37,37.5,37,36.5,36,36};
    
    // variables
    double average_bin(-1),content_histo_xy_tdiff_meanD_tmp(-1),content_histo_pos_xy(-1);
    int x_pos_bin(-1),y_pos_bin(-1);
    
    int counter=0;
    int gcouter =0;
    
    // canvas
    TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
    //TCanvas *cc9 = new TCanvas("cc9","cc9",800,500);
    //TCanvas *cc10 = new TCanvas("cc10","cc10",800,500);
    
    for(Int_t j=0; j<40; j++){
        for(Int_t i=0; i<24; i++){
            //if(i !=0)continue;
            //if(j !=21)continue;
            
            histo_time_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_time_bar_pos_%d_%d",i,j));
            histo_tdiffD_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tdiffD_bar_pos_%d_%d",i,j));
            histo_tdiffR_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tdiffR_bar_pos_%d_%d",i,j));
            histo_tdiffR_bar_pos[i][j]->SetLineColor(kRed);
            
            histo_tD_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tD_bar_pos_%d_%d",i,j));
            histo_tR_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tR_bar_pos_%d_%d",i,j));
            histo_tR_bar_pos[i][j]->SetLineColor(kRed);
            
            histo_cD_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_cD_bar_pos_%d_%d",i,j));
            histo_cR_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_cR_bar_pos_%d_%d",i,j));
            histo_cR_bar_pos[i][j]->SetLineColor(kRed);
            
            histo_c_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_c_bar_pos_%d_%d",i,j));
            histo_c_bar_pos[i][j]->SetLineColor(kBlack);
            
            // titles
            histo_time_bar_pos[i][j]->GetXaxis()->SetTitle("Time [ns]");
            histo_tdiffD_bar_pos[i][j]->GetXaxis()->SetTitle("Time Difference Direct photon [ns]");
            histo_tdiffR_bar_pos[i][j]->GetXaxis()->SetTitle("Time Difference Refelected photon [ns]");
            
            int entries = histo_time_bar_pos[i][j]->GetEntries();
            
            
            if (entries<300000){ //300000 //8000
                histo_xy_tdiff_meanR->Fill(j,i,-1000);
                histo_xy_tdiff_meanD->Fill(j,i,-1000);
                histo_xy_c_mean->Fill(j,i,-1000);
                histo_xy_c_meanD->Fill(j,i,-1000);
                histo_xy_c_meanR->Fill(j,i,-1000);
                continue;
            }
            
            /*
             cc9->Clear();
             cc9->cd();
             histo_tdiffD_bar_pos[i][j]->Draw();
             cc9->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffD_bar_pos_%d.png",counter));
             cc9->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffD_bar_pos_%d.root",counter));
             
             cc10->Clear();
             cc10->cd();
             histo_tdiffR_bar_pos[i][j]->Draw();
             cc10->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffR_bar_pos_%d.png",counter));
             cc10->SaveAs(Form("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffR_bar_pos_%d.root",counter));
             ++counter;
             */
            
            //if (entries<8000)continue;
            
            //histo_pos_xy_entries->Fill(j,23-i,entries);
            histo_pos_xy_entries->Fill(j,i,entries);
            if(true){
                cc1->cd();
                cc1->Update();
                histo_time_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                histo_time_bar_pos[i][j]->Draw();
                //                cc1->Update();
                //                TLine *lin_mean_max= new TLine(0,0,0,1000);
                //                lin_mean_max->SetX1(cop[j]);
                //                lin_mean_max->SetX2(cop[j]);
                //                lin_mean_max->SetY1(gPad->GetUymin());
                //                lin_mean_max->SetY2(gPad->GetUymax());
                //                lin_mean_max->SetLineColor(kRed);
                //                lin_mean_max->Draw();
                //                cc1->Update();
                //                cc1->WaitPrimitive();
            }
            // time
            if(true){
                cc1->cd();
                cc1->Update();
                histo_tD_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                histo_tR_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                if(histo_tD_bar_pos[i][j]->GetEntries()>histo_tR_bar_pos[i][j]->GetEntries()){
                    histo_tD_bar_pos[i][j]->Draw();
                    histo_tR_bar_pos[i][j]->Draw("same");
                }else{
                    histo_tR_bar_pos[i][j]->Draw();
                    histo_tD_bar_pos[i][j]->Draw("same");
                }
                //cc1->Update();
                //cc1->WaitPrimitive();
                //////////////
                histo_xy_t_meanD->Fill(j,i,histo_tD_bar_pos[i][j]->GetMean());
                histo_xy_t_sigmaD->Fill(j,i,histo_tD_bar_pos[i][j]->GetRMS());
                histo_xy_t_meanR->Fill(j,i,histo_tR_bar_pos[i][j]->GetMean());
                histo_xy_t_sigmaR->Fill(j,i,histo_tR_bar_pos[i][j]->GetRMS());
                
                // Cherenkov
                cc1->cd();
                cc1->Update();
                histo_c_bar_pos[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
                double thetaS=fit_cherenkov->GetParameter(1)*1000;
                double sprS=fit_cherenkov->GetParameter(2)*1000;
                histo_cR_bar_pos[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
                double thetaR=fit_cherenkov->GetParameter(1)*1000;
                double sprR=fit_cherenkov->GetParameter(2)*1000;
                histo_cD_bar_pos[i][j]->Fit("fit_cherenkov","M","", -0.05, 0.05);
                double thetaD=fit_cherenkov->GetParameter(1)*1000;
                double sprD=fit_cherenkov->GetParameter(2)*1000;
                
                histo_c_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                histo_cR_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                histo_cD_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                
                histo_c_bar_pos[i][j]->Draw();
                histo_cD_bar_pos[i][j]->Draw("same");
                histo_cR_bar_pos[i][j]->Draw("same");
                //cc1->Update();
                //cc1->WaitPrimitive();
                
                //////////////
                histo_xy_c_meanD->Fill(j,i,thetaD);
                histo_xy_c_sigmaD->Fill(j,i,sprD);
                histo_xy_c_meanR->Fill(j,i,thetaR);
                histo_xy_c_sigmaR->Fill(j,i,sprR);
                histo_xy_c_mean->Fill(j,i,thetaS);
                histo_xy_c_sigma->Fill(j,i,sprS);
                
                histo_xy_c_occu->Fill(j,i,histo_c_bar_pos[i][j]->GetEntries());
                histo_xy_cD_occu->Fill(j,i,histo_cD_bar_pos[i][j]->GetEntries());
                histo_xy_cR_occu->Fill(j,i,histo_cR_bar_pos[i][j]->GetEntries());
                
                g_c_meanD_time->SetPoint(gcouter, histo_tD_bar_pos[i][j]->GetMean(), thetaD);
                g_c_sigmaD_time->SetPoint(gcouter, histo_tD_bar_pos[i][j]->GetRMS(), sprD);
                g_c_meanR_time->SetPoint(gcouter, histo_tR_bar_pos[i][j]->GetMean(), thetaR);
                g_c_sigmaR_time->SetPoint(gcouter, histo_tR_bar_pos[i][j]->GetRMS(), sprR);
                
                ++gcouter;
                
                
            }
            
            if(true){
                int binmax = histo_tdiffD_bar_pos[i][j]->GetMaximumBin();
                double x = histo_tdiffD_bar_pos[i][j]->GetXaxis()->GetBinCenter(binmax);
                fit_gause->SetParLimits(1,x-1,x+1);
                histo_tdiffD_bar_pos[i][j]->Fit("fit_gause","M","", x-2, x+2);
                double meanD=fit_gause->GetParameter(1);
                double sigmaD=fit_gause->GetParameter(2);
                histo_shiftD_diff->Fill(meanD);
                histo_shiftD_sgma->Fill(sigmaD);
                //if(i==0)meanD=-10;
                //if(i==1)meanD=10;
                shiftD_array[i][j]=meanD;
                ///////////////
                histo_xy_tdiff_meanD->Fill(j,i,meanD);
                histo_xy_tdiff_sigmaD->Fill(j,i,sigmaD);
                
                //histo_xy_tdiff_meanD->Fill(j,23-i,meanD);
                //histo_xy_tdiff_sigmaD->Fill(j,23-i,sigmaD);
                
                //                cc1->cd();
                //                cc1->Update();
                //                histo_tdiffD_bar_pos[i][j]->Draw();
                //                cc1->Update();
                //                cc1->WaitPrimitive();
            }
            
            if(true){
                int binmax = histo_tdiffR_bar_pos[i][j]->GetMaximumBin();
                double x = histo_tdiffR_bar_pos[i][j]->GetXaxis()->GetBinCenter(binmax);
                fit_gause->SetParLimits(1,x-1,x+1);
                histo_tdiffR_bar_pos[i][j]->Fit("fit_gause","M","", x-2, x+2);
                double meanR=fit_gause->GetParameter(1);
                double sigmaR=fit_gause->GetParameter(2);
                histo_shiftR_diff->Fill(meanR);
                histo_shiftR_sgma->Fill(sigmaR);
                //if(i==0)meanR=-10;
                //if(i==1)meanR=10;
                shiftR_array[i][j]=meanR;
                ///////////////
                histo_xy_tdiff_meanR->Fill(j,i,meanR);
                histo_xy_tdiff_sigmaR->Fill(j,i,sigmaR);
                //histo_xy_tdiff_meanR->Fill(j,23-i,meanR);
                //histo_xy_tdiff_sigmaR->Fill(j,23-i,sigmaR);
                
                //                                                cc1->cd();
                //                                                cc1->Update();
                //                                                histo_tdiffR_bar_pos[i][j]->Draw();
                //                                                cc1->Update();
                //                                                cc1->WaitPrimitive();
            }
        }
    }
    
    //    TCanvas *cc11 = new TCanvas("cc11","cc11",800,500);
    //    histo_xy_t_meanD->SetMinimum(0);
    //    histo_xy_t_meanD->SetMaximum(100);
    //    histo_xy_t_meanD->Draw("COLZ");
    //    cc11->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_meanD.png");
    //    cc11->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_meanD.root");
    //
    //    TCanvas *cc12 = new TCanvas("cc12","cc12",800,500);
    //    histo_xy_t_sigmaD->SetMinimum(0);
    //    histo_xy_t_sigmaD->SetMaximum(20);
    //    histo_xy_t_sigmaD->Draw("COLZ");
    //    cc12->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_sigmaD.png");
    //    cc12->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_sigmaD.root");
    //
    //    TCanvas *cc13 = new TCanvas("cc13","cc13",800,500);
    //    histo_xy_t_meanR->SetMinimum(0);
    //    histo_xy_t_meanR->SetMaximum(100);
    //    histo_xy_t_meanR->Draw("COLZ");
    //    cc13->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_meanR.png");
    //    cc13->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_meanR.root");
    //
    //    TCanvas *cc14 = new TCanvas("cc14","cc14",800,500);
    //    histo_xy_t_sigmaR->SetMinimum(0);
    //    histo_xy_t_sigmaR->SetMaximum(20);
    //    histo_xy_t_sigmaR->Draw("COLZ");
    //    cc14->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_sigmaR.png");
    //    cc14->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_t_sigmaR.root");
    //////////////////////
    
    // Cherenkov shift and SPR histo
    TCanvas *cc15 = new TCanvas("cc15","cc15",800,500);
    histo_xy_c_meanD->SetMinimum(-5);
    histo_xy_c_meanD->SetMaximum(5);
    histo_xy_c_meanD->Draw("COLZ");
    cc15->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_meanD.png");
    cc15->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_meanD.root");
    
    TCanvas *cc16 = new TCanvas("cc16","cc16",800,500);
    histo_xy_c_sigmaD->SetMinimum(0);
    histo_xy_c_sigmaD->SetMaximum(14);
    histo_xy_c_sigmaD->Draw("COLZ");
    cc16->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigmaD.png");
    cc16->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigmaD.root");
    
    TCanvas *cc17 = new TCanvas("cc17","cc17",800,500);
    histo_xy_c_meanR->SetMinimum(-5);
    histo_xy_c_meanR->SetMaximum(5);
    histo_xy_c_meanR->Draw("COLZ");
    cc17->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_meanR.png");
    cc17->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_meanR.root");
    
    TCanvas *cc18 = new TCanvas("cc18","cc18",800,500);
    histo_xy_c_sigmaR->SetMinimum(0);
    histo_xy_c_sigmaR->SetMaximum(14);
    histo_xy_c_sigmaR->Draw("COLZ");
    cc18->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigmaR.png");
    cc18->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigmaR.root");
    
    TCanvas *cc19 = new TCanvas("cc19","cc19",800,500);
    histo_xy_c_mean->SetMinimum(-5);
    histo_xy_c_mean->SetMaximum(5);
    histo_xy_c_mean->Draw("COLZ");
    cc19->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_mean.png");
    cc19->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_mean.root");
    
    TCanvas *cc20 = new TCanvas("cc20","cc20",800,500);
    histo_xy_c_sigma->SetMinimum(0);
    histo_xy_c_sigma->SetMaximum(14);
    histo_xy_c_sigma->Draw("COLZ");
    cc20->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigma.png");
    cc20->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_sigma.root");
    

    // graphs
    TCanvas *cc21 = new TCanvas("cc21","cc21",800,500);
    TMultiGraph *mg1= new TMultiGraph();
    mg1->Add(g_c_meanD_time);
    mg1->Add(g_c_meanR_time);
    mg1->SetTitle(" correlation between cherenkove shift and photon time; photon time [ns] ; #theta_{c} shift [mrad]");
    mg1->Draw("APL");
    //mg1->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.84);
    //legend_reso->Draw();
    cc21->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/mg_mean.png");
    cc21->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/mg_mean.root");
    
    
    TCanvas *cc22 = new TCanvas("cc22","cc22",800,500);
    TMultiGraph *mg2= new TMultiGraph();
    mg2->Add(g_c_sigmaD_time);
    mg2->Add(g_c_sigmaR_time);
    mg2->SetTitle(" correlation between cherenkove shift and photon time; photon time [ns] ; #theta_{c} shift [mrad]");
    mg2->Draw("APL");
    //mg2->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.84);
    //legend_reso->Draw();
    cc22->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/mg_mean.png");
    cc22->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/mg_mean.root");
    
    
    // Entries histo
    TCanvas *cc23 = new TCanvas("cc23","cc23",800,500);
    histo_xy_c_occu->Draw("COLZ");
    cc23->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_occu.png");
    cc23->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_c_occu.root");
    
    TCanvas *cc24 = new TCanvas("cc24","cc24",800,500);
    histo_xy_cD_occu->Draw("COLZ");
    cc24->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_cD_occu.png");
    cc24->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_cD_occu.root");
    
    TCanvas *cc25 = new TCanvas("cc25","cc25",800,500);
    histo_xy_cR_occu->Draw("COLZ");
    cc25->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_cR_occu.png");
    cc25->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_cR_occu.root");
    

    
    
    
    
    /*
     //
     //    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
     //    histo_shiftR_diff->SetLineColor(2);
     //    histo_shiftD_diff->Draw();
     //    histo_shiftR_diff->Draw("same");
     //
     //    TCanvas *cc4 = new TCanvas("cc4","cc4",800,500);
     //    histo_shiftR_sgma->SetLineColor(2);
     //    histo_shiftD_sgma->Draw();
     //    histo_shiftR_sgma->Draw("same");
     
     TCanvas *cc5 = new TCanvas("cc5","cc5",800,500);
     histo_xy_tdiff_meanD->SetMinimum(-1.5);
     histo_xy_tdiff_meanD->SetMaximum(1.5);
     histo_xy_tdiff_meanD->Draw("COLZ");
     cc5->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_meanD.png");
     cc5->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_meanD.root");
     
     TCanvas *cc6 = new TCanvas("cc6","cc6",800,500);
     histo_xy_tdiff_sigmaD->SetMinimum(0);
     histo_xy_tdiff_sigmaD->SetMaximum(3);
     histo_xy_tdiff_sigmaD->Draw("COLZ");
     cc6->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_sigmaD.png");
     cc6->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_sigmaD.root");
     //
     TCanvas *cc7 = new TCanvas("cc7","cc7",800,500);
     histo_xy_tdiff_meanR->SetMinimum(-1.5);
     histo_xy_tdiff_meanR->SetMaximum(1.5);
     histo_xy_tdiff_meanR->Draw("COLZ");
     cc7->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_meanR.png");
     cc7->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_meanR.root");
     //
     TCanvas *cc8 = new TCanvas("cc8","cc8",800,500);
     histo_xy_tdiff_sigmaR->SetMinimum(0);
     histo_xy_tdiff_sigmaR->SetMaximum(3);
     histo_xy_tdiff_sigmaR->Draw("COLZ");
     cc8->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_sigmaR.png");
     cc8->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_xy_tdiff_sigmaR.root");
     
     
     //    TCanvas *cc9 = new TCanvas("cc9","cc9",800,500);
     //    histo_tdiffD_bar_pos[6][20]->Draw();
     //    cc9->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffD_bar_pos.png");
     //    cc9->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffD_bar_pos.root");
     //
     //
     //    TCanvas *cc10 = new TCanvas("cc10","cc10",800,500);
     //    histo_tdiffR_bar_pos[6][20]->Draw();
     //    cc10->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffR_bar_pos.png");
     //    cc10->SaveAs("/Users/ahmed/GlueX_DIRC_Calib/withtcorrection/histo_tdiffR_bar_pos.root");
     
     
     
     
     //    cout<<"double shiftD_tdiff[24][40]={"<<endl;
     //    for(Int_t i=0; i<24; i++){
     //        for(Int_t j=0; j<40; j++){
     //            cout<<shiftD_array[i][j] <<",";
     //        }
     //    }
     //    cout<<" } ##### end"<<endl;
     //
     //    cout<<"double shiftR_tdiff[24][40]={"<<endl;
     //    for(Int_t i=0; i<24; i++){
     //        for(Int_t j=0; j<40; j++){
     //            cout<<shiftR_array[i][j] <<",";
     //        }
     //    }
     //    cout<<" } ##### end"<<endl;
     
     
     */
    
    
    //cout<<shiftD_array[0][20] <<endl;
    //cout<<shiftD_array[1][20] <<endl;
    
    }
