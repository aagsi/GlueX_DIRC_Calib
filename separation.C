#include <TLegend.h>
#include "glxtools.C"
#include "TNtuple.h"
#include "TStopwatch.h"
TStopwatch timer;

void separation(){
    timer.Start();
    TString path ="histo_f5.root";
    cout<<"path= " <<path<<endl;
    TFile *f = new TFile(path, "READ");
    gStyle->SetPalette(55);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();
    bool save_bool(true);
    TH2F*  histo_separation_mom_map[22];
    TH2F*  histo_separation_mom_map_occ[22];
    TH2F*  histo_separation_mom_map_ratio[22];
    double v11(0) ,v22(0);
    for(Int_t i=0; i<22; i++) {
        v11 = v11+0.5;
        v22 = v11+0.5;
        histo_separation_mom_map[i] = new TH2F(Form("histo_separation_mom_map_%d",i),Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ; X bin number [#]; Bar number [#]",v11,v22) , 40, 0, 40, 24, 0, 24);
        histo_separation_mom_map_occ[i] = new TH2F(Form("histo_separation_mom_map_occ_%d",i),Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ; X bin number [#]; Bar number [#]",v11,v22) , 40, 0, 40, 24, 0, 24);
        histo_separation_mom_map_ratio[i] = new TH2F(Form("histo_separation_mom_map_ratio_%d",i),Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ; X bin number [#]; Bar number [#]",v11,v22) , 40, 0, 40, 24, 0, 24);
    }
    TH2F * histo_cherenkov_shift_map = new TH2F( "histo_cherenkov_shift_map" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_cherenkov_shift_map_occ = new TH2F( "histo_cherenkov_shift_map_occ" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_cherenkov_shift_map_ratio = new TH2F( "histo_cherenkov_shift_map_ratio" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_cherenkov_spr_map = new TH2F( "histo_cherenkov_spr_map" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_cherenkov_spr_map_occ = new TH2F( "histo_cherenkov_spr_map_occ" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);
    TH2F * histo_cherenkov_spr_map_ratio = new TH2F( "histo_cherenkov_spr_map_ratio" , "; X bin number [#]; Bar number [#]", 40, 0, 40, 24, 0, 24);

    /////////////////////////////////////
    // calculate cherenkov angle
    // Double_t momentum=3.5;
    Int_t pdg[]= {11,13,211,321,2212};
    Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    Double_t fit_angleK[10]={0};
    Double_t fit_anglePi[10]={0};
    Double_t momentum[] = {2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
    //Double_t momentum[] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
    for(int f=0;f<10;f++){
        fit_angleK[f] = acos(sqrt(momentum[f]*momentum[f] + mass[3]*mass[3])/momentum[f]/1.4738)-0.00;
        fit_anglePi[f]= acos(sqrt(momentum[f]*momentum[f] + mass[2]*mass[2])/momentum[f]/1.4738)-0.00;
    }

    //////////////////////////////////////
    TSpectrum *spect = new TSpectrum(10);
    double minChisto_Cherenkov = 0.6;
    double maxChisto_Cherenkov = 0.9;
    TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChisto_Cherenkov,maxChisto_Cherenkov);
    double cut_cangle=0.04;
    double cherenkovreco[5],spr[5];
    TGaxis::SetMaxDigits(3);
    double sigma[]={0.01,0.01,0.01,0.010,0.01,0.01};
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
    //////////////////////////////////////

    TH1F *histo_LnDiff[42][26][22][5], *histo_Cherenkov[42][26][22][5];
    for(int i=0; i<42; i++){
        for(int d=0; d<26; d++){
            for(int j=0; j<22; j++){
                for(int k=2; k<4; k++){
                    histo_LnDiff[i][d][j][k] = (TH1F*)f->Get(Form("histo_LnDiff_%d_%d_%d_%d",i,d,j,k));
                    histo_Cherenkov[i][d][j][k] = (TH1F*)f->Get(Form("histo_Cherenkov_%d_%d_%d_%d",i,d,j,k));
                }
            }
        }
    }

    // canvas
    TCanvas *cc1 = new TCanvas("cc1","cc1",800,500);
    int counter=0;
    for(int i=10; i<32; i++){ // xbin       0   42
        for(int d=0; d<14; d++){ // bar     0   26
            for(int j=6; j<12; j++){ // mom 0   22

                /////////////////////////////////////////
                /////// calculate separation power //////
                /////////////////////////////////////////

                cc1->Clear();
                cc1->cd();
                cc1->Update();
                if( (histo_LnDiff[i][d][j][2]->GetEntries()<75 ) || (histo_LnDiff[i][d][j][3]->GetEntries()<75)) continue;
                //if( (histo_LnDiff[i][d][j][2]->GetEntries()<200 ) || (histo_LnDiff[i][d][j][3]->GetEntries()<200)) continue;

                TF1 *ff;
                double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0;
                if(histo_LnDiff[i][d][j][3]->GetEntries()>100){
                    histo_LnDiff[i][d][j][3]->Fit("gaus","S");
                    ff = histo_LnDiff[i][d][j][3]->GetFunction("gaus");
                    ff->SetLineColor(1);
                    m1=ff->GetParameter(1);
                    s1=ff->GetParameter(2);
                }
                if(histo_LnDiff[i][d][j][2]->GetEntries()>100){
                    histo_LnDiff[i][d][j][2]->Fit("gaus","S");
                    ff = histo_LnDiff[i][d][j][2]->GetFunction("gaus");
                    ff->SetLineColor(1);
                    m2=ff->GetParameter(1);
                    s2=ff->GetParameter(2);
                }
                if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
                cout<<"#######  sep= "<<sep<<endl;

                histo_separation_mom_map[j]->Fill(i,d,sep);
                histo_separation_mom_map_occ[j]->Fill(i,d);
                ////////////////////////////////

                histo_LnDiff[i][d][j][2]->SetTitle(Form("sep = %2.2f s.d.",sep));
                histo_LnDiff[i][d][j][3]->SetTitle(Form("sep = %2.2f s.d.",sep));
                histo_LnDiff[i][d][j][2]->Draw();
                histo_LnDiff[i][d][j][3]->Draw("same");

                TLegend *lnpl = new TLegend(0.7,0.67,0.9,0.85);
                lnpl->SetFillColor(0);
                lnpl->SetFillStyle(0);
                lnpl->SetBorderSize(0);
                lnpl->SetFillStyle(0);
                lnpl->AddEntry(histo_LnDiff[i][d][j][2],"pions","lp");
                lnpl->AddEntry(histo_LnDiff[i][d][j][3],"kaons","lp");
                lnpl->Draw();

                TLegend *leg_sep = new TLegend(0.14787, 0.570667, 0.348371, 0.749333);
                leg_sep->SetFillColor(0);
                leg_sep->SetFillStyle(0);
                leg_sep->SetBorderSize(0);
                leg_sep->SetFillStyle(0);
                leg_sep->AddEntry(histo_LnDiff[i][d][j][2],Form("Separation = %1.2f",sep),"");
                leg_sep->Draw();

                if(save_bool){
                    cc1->SaveAs(Form("histo_like_%d.png",counter));
                    cc1->SaveAs(Form("histo_like_%d.root",counter));
                }

                //cc1->Update();
                //cc1->WaitPrimitive();

                ////////////////////////////////
                //////// Cherenkove angle //////
                ////////////////////////////////

                cc1->Clear();
                cc1->cd();
                cc1->Update();

                //scal
                if(histo_Cherenkov[i][d][j][2]->GetMaximum()>0) histo_Cherenkov[i][d][j][2]->Scale(1/histo_Cherenkov[i][d][j][2]->GetMaximum());
                if(histo_Cherenkov[i][d][j][3]->GetMaximum()>0) histo_Cherenkov[i][d][j][3]->Scale(1/histo_Cherenkov[i][d][j][3]->GetMaximum());

                for(int s=2; s<4; s++){
                    if(histo_Cherenkov[i][d][j][s]->GetEntries()<20) continue;

                    int nfound = spect->Search(histo_Cherenkov[i][d][j][s],1,"goff",0.9);
                    if(nfound>0) cherenkovreco[s] = spect->GetPositionX()[0];
                    else cherenkovreco[s] =  histo_Cherenkov[i][d][j][s]->GetXaxis()->GetBinCenter(histo_Cherenkov[i][d][j][s]->GetMaximumBin());
                    if(cherenkovreco[s]>0.85) cherenkovreco[s]=0.82;

                    if(s==2)  fit->SetLineColor(kBlue);
                    if(s==3)  fit->SetLineColor(kRed);
                    fit->SetParameters(100,cherenkovreco[s],0.010,10);
                    fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                    fit->SetParLimits(0,0.1,1E6);
                    fit->SetParLimits(1,cherenkovreco[s]-2*cut_cangle,cherenkovreco[s]+2*cut_cangle);
                    fit->SetParLimits(2,0.005,0.030); // width
                    histo_Cherenkov[i][d][j][s]->Fit("fgaus","I","",cherenkovreco[s]-cut_cangle,cherenkovreco[s]+cut_cangle);
                    histo_Cherenkov[i][d][j][s]->Fit("fgaus","M","",cherenkovreco[s]-cut_cangle,cherenkovreco[s]+cut_cangle);

                    cherenkovreco[s] = fit->GetParameter(1);
                    spr[s] = fit->GetParameter(2);

                    if(s==2){
                        double track_mom = j/2.0;
                        double ExAnglePi= acos(sqrt(track_mom*track_mom + mass[s]*mass[s])/track_mom/1.4738);
                        double ExMeandiff = cherenkovreco[s] - ExAnglePi;

                        histo_cherenkov_shift_map->Fill(i,d, ExMeandiff *1000 ); // for Pions
                        histo_cherenkov_shift_map_occ->Fill(i,d);

                        histo_cherenkov_spr_map->Fill(i,d, spr[s] *1000 ); // for Pions
                        histo_cherenkov_spr_map_occ->Fill(i,d);
                    }
                }

                gStyle->SetOptTitle(0);
                gStyle->SetOptStat(0);
                gStyle->SetOptFit(0);

                histo_Cherenkov[i][d][j][2]->GetXaxis()->SetRangeUser(0.7,0.9);
                histo_Cherenkov[i][d][j][2]->GetYaxis()->SetRangeUser(0,1.2);
                histo_Cherenkov[i][d][j][2]->Draw();
                histo_Cherenkov[i][d][j][3]->Draw("same");

                TLine *line = new TLine(0,0,0,1000);
                line->SetX1(cherenkovreco[3]); // theory_angle[3]
                line->SetX2(cherenkovreco[3]); // theory_angle[3]
                line->SetY1(0);
                line->SetY2(1.2);
                line->SetLineColor(kRed);
                line->Draw();

                TLine *line2 = new TLine(0,0,0,1000);
                line2->SetX1(cherenkovreco[2]); // theory_angle[2]
                line2->SetX2(cherenkovreco[2]); // theory_angle[2]
                line2->SetY1(0);
                line2->SetY2(1.2);
                line2->SetLineColor(kBlue);
                line2->Draw();

                TLine *line3 = new TLine(0,0,0,1000);
                line3->SetLineStyle(2);
                line3->SetX1(0.5*(theory_angle[2]+theory_angle[3])-cut_cangle);
                line3->SetX2(0.5*(theory_angle[2]+theory_angle[3])-cut_cangle);
                line3->SetY1(0);
                line3->SetY2(1.2);
                line3->SetLineColor(1);
                line3->Draw();

                TLine *line4 = new TLine(0,0,0,1000);
                line4->SetLineStyle(2);
                line4->SetX1(0.5*(theory_angle[2]+theory_angle[3])+cut_cangle);
                line4->SetX2(0.5*(theory_angle[2]+theory_angle[3])+cut_cangle);
                line4->SetY1(0);
                line4->SetY2(1.2);
                line4->SetLineColor(1);
                line4->Draw();

                TLegend *leg = new TLegend(0.1,0.5,0.4,0.85);
                leg->SetFillColor(0);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(histo_Cherenkov[i][d][j][2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
                leg->AddEntry(histo_Cherenkov[i][d][j][3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
                leg->AddEntry(histo_Cherenkov[i][d][j][2],Form("#sigma_{c}^{#pi} = %2.1f mrad",spr[2]*1000),"");
                leg->AddEntry(histo_Cherenkov[i][d][j][3],Form("#sigma_{c}^{K} = %2.1f mrad",spr[3]*1000),"");
                leg->Draw();

                TLegend *lnpa = new TLegend(0.7,0.67,0.9,0.85);
                lnpa->SetFillColor(0);
                lnpa->SetFillStyle(0);
                lnpa->SetBorderSize(0);
                lnpa->SetFillStyle(0);
                lnpa->AddEntry(histo_Cherenkov[i][d][j][2],"pions","lp");
                lnpa->AddEntry(histo_Cherenkov[i][d][j][3],"kaons","lp");
                lnpa->Draw();

                gStyle->SetOptTitle(0);
                gStyle->SetOptStat(0);
                gStyle->SetOptFit(0);

                if(save_bool){
                    cc1->SaveAs(Form("histo_cher_%d.png",counter));
                    cc1->SaveAs(Form("histo_cher_%d.root",counter));
                }
                //cc1->Update();
                //cc1->WaitPrimitive();
                ++counter;
            }
        }
    }

    TCanvas *cc2 = new TCanvas("cc2","cc2",800,500);
    double mom_flag1(0) ,mom_flag2(0);
    for(Int_t k=6; k<12; k++){
        mom_flag1 = mom_flag1+0.5;
        mom_flag2 = mom_flag1+0.5;
        histo_separation_mom_map_ratio[k] = (TH2F*)histo_separation_mom_map[k]->Clone();
        histo_separation_mom_map_ratio[k]->Divide(histo_separation_mom_map_occ[k]);
        histo_separation_mom_map_ratio[k]->SetTitle(Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ", mom_flag1,mom_flag2));
        histo_separation_mom_map_occ[k]->SetTitle(Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ", mom_flag1,mom_flag2));
        histo_separation_mom_map[k]->SetTitle(Form("#pi/K Separation @ Momentum [%1.1f - %1.1f] GeV/c ", mom_flag1,mom_flag2));
        histo_separation_mom_map_ratio[k]->SetMinimum(0);
        histo_separation_mom_map_ratio[k]->SetMaximum(3);
        histo_separation_mom_map_ratio[k]->Draw("COLZ");
        cc2->Update();
        cc2->SaveAs(Form("histo_sep_mom_%d.png",k));
        cc2->SaveAs(Form("histo_sep_mom_%d.root",k));
    }


    TCanvas *cc3 = new TCanvas("cc3","cc3",800,500);
    histo_cherenkov_shift_map_ratio = (TH2F*)histo_cherenkov_shift_map->Clone();
    histo_cherenkov_shift_map_ratio->Divide(histo_cherenkov_shift_map_occ);
    histo_cherenkov_shift_map_ratio->SetTitle("#theta_{c} Shift ");
    histo_cherenkov_shift_map_ratio->SetMinimum(-5);
    histo_cherenkov_shift_map_ratio->SetMaximum(5);
    histo_cherenkov_shift_map_ratio->Draw("COLZ");
    cc3->Update();
    cc3->SaveAs("histo_shift_%d.png");
    cc3->SaveAs("histo_shift_%d.root");


/////////////////

    TCanvas *cc4 = new TCanvas("cc4","cc4",800,500);
    histo_cherenkov_spr_map_ratio = (TH2F*)histo_cherenkov_spr_map->Clone();
    histo_cherenkov_spr_map_ratio->Divide(histo_cherenkov_spr_map_occ);
    histo_cherenkov_spr_map_ratio->SetTitle("SPR ");
    histo_cherenkov_spr_map_ratio->SetMinimum(6);
    histo_cherenkov_spr_map_ratio->SetMaximum(13);
    histo_cherenkov_spr_map_ratio->Draw("COLZ");
    cc4->Update();
    cc4->SaveAs("histo_spr_%d.png");
    cc4->SaveAs("histo_spr_%d.root");

    //    if(true){
    //        cout<<"##### Histograms "<<endl;
    //        glx_drawDigi("m,p,v\n",0);
    //        glx_canvasAdd("r_cherenkov",800,400);
    //        histo_cherenkov->Draw();
    //
    //        glx_canvasAdd("r_tdiff",800,400);
    //        histo_tdiff->Draw();
    //
    //
    //
    //        glx_canvasAdd("r_rho_mass",800,400);
    //        hist_ev_rho_mass->Draw();
    //        glx_canvasGet("r_rho_mass")->Update();
    //        TLine *lin_rho_mass_max= new TLine(0,0,0,1000);
    //        lin_rho_mass_max->SetX1(mass_rho_max);
    //        lin_rho_mass_max->SetX2(mass_rho_max);
    //        lin_rho_mass_max->SetY1(gPad->GetUymin());
    //        lin_rho_mass_max->SetY2(gPad->GetUymax());
    //        lin_rho_mass_max->SetLineColor(kBlack);
    //        lin_rho_mass_max->Draw();
    //        TLine *lin_rho_mass_min= new TLine(0,0,0,1000);
    //        lin_rho_mass_min->SetX1(mass_rho_mini);
    //        lin_rho_mass_min->SetX2(mass_rho_mini);
    //        lin_rho_mass_min->SetY1(gPad->GetUymin());
    //        lin_rho_mass_min->SetY2(gPad->GetUymax());
    //        lin_rho_mass_min->SetLineColor(kBlack);
    //        lin_rho_mass_min->Draw();
    //        glx_canvasGet("r_rho_mass")->Update();
    //
    //        glx_canvasAdd("r_phi_mass",800,400);
    //        hist_ev_phi_mass->Draw();
    //        glx_canvasGet("r_phi_mass")->Update();
    //        TLine *lin_phi_mass_max= new TLine(0,0,0,1000);
    //        lin_phi_mass_max->SetX1(mass_phi_max);
    //        lin_phi_mass_max->SetX2(mass_phi_max);
    //        lin_phi_mass_max->SetY1(gPad->GetUymin());
    //        lin_phi_mass_max->SetY2(gPad->GetUymax());
    //        lin_phi_mass_max->SetLineColor(kBlack);
    //        lin_phi_mass_max->Draw();
    //        TLine *lin_phi_mass_min= new TLine(0,0,0,1000);
    //        lin_phi_mass_min->SetX1(mass_phi_mini);
    //        lin_phi_mass_min->SetX2(mass_phi_mini);
    //        lin_phi_mass_min->SetY1(gPad->GetUymin());
    //        lin_phi_mass_min->SetY2(gPad->GetUymax());
    //        lin_phi_mass_min->SetLineColor(kBlack);
    //        lin_phi_mass_min->Draw();
    //        glx_canvasGet("r_phi_mass")->Update();
    //
    //        glx_canvasAdd("r_missing_mass_phi",800,400);
    //        hist_ev_missing_mass_phi->Draw();
    //        glx_canvasGet("r_missing_mass_phi")->Update();
    //        TLine *lin_phi_miss_mass_max= new TLine(0,0,0,1000);
    //        lin_phi_miss_mass_max->SetX1(miss_mass_phi_max);
    //        lin_phi_miss_mass_max->SetX2(miss_mass_phi_max);
    //        lin_phi_miss_mass_max->SetY1(gPad->GetUymin());
    //        lin_phi_miss_mass_max->SetY2(gPad->GetUymax());
    //        lin_phi_miss_mass_max->SetLineColor(kBlack);
    //        lin_phi_miss_mass_max->Draw();
    //        TLine *lin_phi_miss_mass_min= new TLine(0,0,0,1000);
    //        lin_phi_miss_mass_min->SetX1(miss_mass_phi_mini);
    //        lin_phi_miss_mass_min->SetX2(miss_mass_phi_mini);
    //        lin_phi_miss_mass_min->SetY1(gPad->GetUymin());
    //        lin_phi_miss_mass_min->SetY2(gPad->GetUymax());
    //        lin_phi_miss_mass_min->SetLineColor(kBlack);
    //        lin_phi_miss_mass_min->Draw();
    //        glx_canvasGet("r_missing_mass_phi")->Update();
    //
    //        glx_canvasAdd("r_missing_mass_rho",800,400);
    //        hist_ev_missing_mass_rho->Draw();
    //        glx_canvasGet("r_missing_mass_rho")->Update();
    //        TLine *lin_rho_miss_mass_max= new TLine(0,0,0,1000);
    //        lin_rho_miss_mass_max->SetX1(miss_mass_rho_max);
    //        lin_rho_miss_mass_max->SetX2(miss_mass_rho_max);
    //        lin_rho_miss_mass_max->SetY1(gPad->GetUymin());
    //        lin_rho_miss_mass_max->SetY2(gPad->GetUymax());
    //        lin_rho_miss_mass_max->SetLineColor(kBlack);
    //        lin_rho_miss_mass_max->Draw();
    //        TLine *lin_rho_miss_mass_min= new TLine(0,0,0,1000);
    //        lin_rho_miss_mass_min->SetX1(miss_mass_rho_mini);
    //        lin_rho_miss_mass_min->SetX2(miss_mass_rho_mini);
    //        lin_rho_miss_mass_min->SetY1(gPad->GetUymin());
    //        lin_rho_miss_mass_min->SetY2(gPad->GetUymax());
    //        lin_rho_miss_mass_min->SetLineColor(kBlack);
    //        lin_rho_miss_mass_min->Draw();
    //        glx_canvasGet("r_missing_mass_rho")->Update();
    //
    //        glx_canvasAdd("r_chi_phi",800,400);
    //        hist_ev_chi_phi->Draw();
    //        glx_canvasGet("r_chi_phi")->Update();
    //        TLine *lin_chi_phi_max= new TLine(0,0,0,1000);
    //        lin_chi_phi_max->SetX1(chisq_phi_max);
    //        lin_chi_phi_max->SetX2(chisq_phi_max);
    //        lin_chi_phi_max->SetY1(gPad->GetUymin());
    //        lin_chi_phi_max->SetY2(gPad->GetUymax());
    //        lin_chi_phi_max->SetLineColor(kBlack);
    //        lin_chi_phi_max->Draw();
    //        TLine *lin_chi_phi_min= new TLine(0,0,0,1000);
    //        lin_chi_phi_min->SetX1(chisq_phi_mini);
    //        lin_chi_phi_min->SetX2(chisq_phi_mini);
    //        lin_chi_phi_min->SetY1(gPad->GetUymin());
    //        lin_chi_phi_min->SetY2(gPad->GetUymax());
    //        lin_chi_phi_min->SetLineColor(kBlack);
    //        lin_chi_phi_min->Draw();
    //        glx_canvasGet("r_chi_phi")->Update();
    //
    //        glx_canvasAdd("r_chi_rho",800,400);
    //        hist_ev_chi_rho->Draw();
    //        glx_canvasGet("r_chi_rho")->Update();
    //        TLine *lin_chi_rho_max= new TLine(0,0,0,1000);
    //        lin_chi_rho_max->SetX1(chisq_rho_max);
    //        lin_chi_rho_max->SetX2(chisq_rho_max);
    //        lin_chi_rho_max->SetY1(gPad->GetUymin());
    //        lin_chi_rho_max->SetY2(gPad->GetUymax());
    //        lin_chi_rho_max->SetLineColor(kBlack);
    //        lin_chi_rho_max->Draw();
    //        TLine *lin_chi_rho_min= new TLine(0,0,0,1000);
    //        lin_chi_rho_min->SetX1(chisq_rho_mini);
    //        lin_chi_rho_min->SetX2(chisq_rho_mini);
    //        lin_chi_rho_min->SetY1(gPad->GetUymin());
    //        lin_chi_rho_min->SetY2(gPad->GetUymax());
    //        lin_chi_rho_min->SetLineColor(kBlack);
    //        lin_chi_rho_min->Draw();
    //        glx_canvasGet("r_chi_rho")->Update();
    //    }

    cout<<"##### End Analyses "<<endl;
    //cout<<"####### @ 3.5 GeV/c fit_angleK "<< fit_angleK[1]<<"  fit_anglePi "<<fit_anglePi[1]<<endl;
    timer.Stop();
    printf(" RT=%7.3f s, Cpu=%7.3f s",timer.RealTime(),timer.CpuTime());
}
