
{
    
    gStyle->SetPalette(55);
    
    gStyle->SetOptStat(0);
    TH1F*  histo_time_bar_pos[24][40];
    TH1F*  histo_tdiffD_bar_pos[24][40];
    TH1F*  histo_tdiffR_bar_pos[24][40];
    
    TString path ="/Users/ahmed/GlueX_DIRC_Calib/histo1.root";
    cout<<"path= " <<path<<endl;
    TFile *f = new TFile(path, "READ");
    
//    TH2F * histo_pos_xy_entries = new TH2F( "histo_pos_xy_entries" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
//    TH2F * histo_pos_xy_meanD = new TH2F( "histo_pos_xy_meanD" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
//    TH2F * histo_pos_xy_sigmaD = new TH2F( "histo_pos_xy_sigmaD" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
//
//    TH2F * histo_pos_xy_meanR = new TH2F( "histo_pos_xy_meanR" , "; Bar number ; Bar Hit X", 24, 0, 24, 40, 0, 40);
//    TH2F * histo_pos_xy_sigmaR = new TH2F( "histo_pos_xy_sigmaR" , "; Bar number ; Bar Hit X ", 24, 0, 24, 40, 0, 40);
    
    
    TH2F * histo_pos_xy_entries = new TH2F( "histo_pos_xy_entries" , "; Bar Hit X [cm]; Bar number [#]  ", 40, 0, 40 ,24, 0, 24);
    TH2F * histo_pos_xy_meanD = new TH2F( "histo_pos_xy_meanD" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_pos_xy_sigmaD = new TH2F( "histo_pos_xy_sigmaD" , "; Bar Hit X [cm]; Bar number [#] ",  40, 0, 40,24, 0, 24);
    
    TH2F * histo_pos_xy_meanR = new TH2F( "histo_pos_xy_meanR" , "; Bar Hit X [cm]; Bar number [#]", 40, 0, 40,24, 0, 24);
    TH2F * histo_pos_xy_sigmaR = new TH2F( "histo_pos_xy_sigmaR" , "; Bar Hit X [cm]; Bar number [#] ", 40, 0, 40,24, 0, 24);
    
    
    
    
    
    TH1F* histo_shiftD_diff = new TH1F("histo_shiftD_diff","; Mean per segment [ns]; entries [#]",100,-5,5);
    TH1F* histo_shiftR_diff = new TH1F("histo_shiftR_diff","; Mean per segment [ns]; entries [#]",100,-5,5);
    
    TH1F* histo_shiftD_sgma = new TH1F("histo_shiftD_sgma","; Sigma per segment [ns]; entries [#]",100,0,5);
    TH1F* histo_shiftR_sgma = new TH1F("histo_shiftR_sgma","; Sigma per segment [ns]; entries [#]",100,0,5);
    
    
    histo_pos_xy_meanD->SetTitle("Time Difference Shift Direct photons");
    histo_pos_xy_meanR->SetTitle("Time Difference Shift Reflected photons");
    histo_pos_xy_sigmaD->SetTitle("Time Difference Sigma Direct photons");
    histo_pos_xy_sigmaR->SetTitle("Time Difference Sigma Reflected photons");
    
    TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
    
    TF1 *fit_gause = new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-10,10);
    //TF1 *fit_gause =  new TF1("fit_gause","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",-10,10);
    
    fit_gause->SetLineColor(kBlack);
    fit_gause->SetParameters(100,9,2);
    fit_gause->SetParNames("p0","mean ","sigma");
    fit_gause->SetParLimits(0,0.1,1E6);
    fit_gause->SetParLimits(1,-3,3);
    fit_gause->SetParLimits(2,0.5,3);
    
    Double_t shiftD_array[24][40]={0};
    Double_t shiftR_array[24][40]={0};
    
    
    double average_bin(-1),content_histo_pos_xy_meanD_tmp(-1),content_histo_pos_xy(-1);
    int x_pos_bin(-1),y_pos_bin(-1);
    
    
    
    double cop[]= {62,62,62,62,62,62,62,62,60,58,56,54,53,52,52,51,50,49,47,46,46,45,45,44,44,43,43,42,41,41,40,40,38,37.5,37,37.5,37,36.5,36,36};
    
    
    
    for(Int_t j=0; j<40; j++){
        for(Int_t i=0; i<24; i++){
            //if(i !=0)continue;
            //if(j !=21)continue;
            
            histo_time_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_time_bar_pos_%d_%d",i,j));
            histo_tdiffD_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tdiffD_bar_pos_%d_%d",i,j));
            histo_tdiffR_bar_pos[i][j] = (TH1F*)f->Get(Form("histo_tdiffR_bar_pos_%d_%d",i,j));
            
            
            histo_time_bar_pos[i][j]->GetXaxis()->SetTitle("Time [ns]");
            histo_tdiffD_bar_pos[i][j]->GetXaxis()->SetTitle("Time Difference Direct photon [ns]");
            histo_tdiffR_bar_pos[i][j]->GetXaxis()->SetTitle("Time Difference Refelected photon [ns]");
            
            int entries = histo_time_bar_pos[i][j]->GetEntries();
            
            
            if (entries<300000)continue; //300000 //8000
            //if (entries<8000)continue;
            
            histo_pos_xy_entries->Fill(j,23-i,entries);
            
            //histo_pos_xy_entries->Fill(i,j,entries);
            
            if(true){
                cc1->cd();
                cc1->Update();
                histo_time_bar_pos[i][j]->SetTitle(Form("Bar %d X Bin %d",i,j));
                histo_time_bar_pos[i][j]->Draw();
                cc1->Update();
                TLine *lin_mean_max= new TLine(0,0,0,1000);
                lin_mean_max->SetX1(cop[j]);
                lin_mean_max->SetX2(cop[j]);
                lin_mean_max->SetY1(gPad->GetUymin());
                lin_mean_max->SetY2(gPad->GetUymax());
                lin_mean_max->SetLineColor(kRed);
                lin_mean_max->Draw();
                cc1->Update();
                cc1->WaitPrimitive();
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
                //histo_pos_xy_meanD->Fill(i,j,meanD);
                //histo_pos_xy_sigmaD->Fill(i,j,sigmaD);
                
                histo_pos_xy_meanD->Fill(j,23-i,meanD);
                histo_pos_xy_sigmaD->Fill(j,23-i,sigmaD);
                
                                                                cc1->cd();
                                                                cc1->Update();
                                                                histo_tdiffD_bar_pos[i][j]->Draw();
                                                                cc1->Update();
                                                                cc1->WaitPrimitive();
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
                //histo_pos_xy_meanR->Fill(i,j,meanR);
                //histo_pos_xy_sigmaR->Fill(i,j,sigmaR);
                
                histo_pos_xy_meanR->Fill(j,23-i,meanR);
                histo_pos_xy_sigmaR->Fill(j,23-i,sigmaR);
                
                //                                                cc1->cd();
                //                                                cc1->Update();
                //                                                histo_tdiffR_bar_pos[i][j]->Draw();
                //                                                cc1->Update();
                //                                                cc1->WaitPrimitive();
            }
        }
    }
    
    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
    histo_pos_xy_entries->Draw("COLZ");
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
    histo_pos_xy_meanD->SetMinimum(-1.5);
    histo_pos_xy_meanD->SetMaximum(1.5);
    histo_pos_xy_meanD->Draw("COLZ");
    
        TCanvas *cc6 = new TCanvas("cc6","cc6",800,500);
        histo_pos_xy_sigmaD->SetMinimum(0);
        histo_pos_xy_sigmaD->SetMaximum(3);
        histo_pos_xy_sigmaD->Draw("COLZ");
    //
    TCanvas *cc7 = new TCanvas("cc7","cc7",800,500);
    histo_pos_xy_meanR->SetMinimum(-1.5);
    histo_pos_xy_meanR->SetMaximum(1.5);
    histo_pos_xy_meanR->Draw("COLZ");
    //
        TCanvas *cc8 = new TCanvas("cc8","cc8",800,500);
        histo_pos_xy_sigmaR->SetMinimum(0);
        histo_pos_xy_sigmaR->SetMaximum(3);
        histo_pos_xy_sigmaR->Draw("COLZ");
    
    
    
    
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
    
    

    
    
    //cout<<shiftD_array[0][20] <<endl;
    //cout<<shiftD_array[1][20] <<endl;
    
}
