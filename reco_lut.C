#define glx__sim

#include "/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/classes/DrcEvent.h"
#include "/lustre/nyx/panda/aali/gluex/gluex_top/halld_recon/halld_recon-4.2.0/src/plugins/Analysis/lut_dirc/DrcLutNode.h"

#include "glxtools.C"
#define PI 3.14159265


// gPDF = 1 Creat PDF per pix
// gPDF = 2 Calculate Sepration from PDF

// gCherenkov_Correction = 1 Creat histo per PMT
// gCherenkov_Correction = 2 Apply per PMT correction

void reco_lut(TString infile="vol/tree_060772.root",TString inlut="lut/lut_12/lut_all_avr.root", int gPDF=0, int gCherenkov_Correction=0, int xbar=-1, int ybar=-1, double moms=3.75){
    if(!glx_initc(infile,1,"data/reco_lut_sim")) return;
    const int nodes = glx_maxch;
    const int luts = 24;
    TFile *fLut = new TFile(inlut);
    TTree *tLut=(TTree *) fLut->Get("lut_dirc") ;
    TClonesArray* cLut[luts];
    for(int l=0; l<luts; l++){
        cLut[l] = new TClonesArray("DrcLutNode");
        tLut->SetBranchAddress(Form("LUT_%d",l),&cLut[l]);
    }
    tLut->GetEntry(0);
    DrcLutNode *lutNode[luts][nodes];
    for(int l=0; l<luts; l++){
        for(int i=0; i<nodes; i++) lutNode[l][i] = (DrcLutNode*) cLut[l]->At(i);
    }
    TGaxis::SetMaxDigits(4);
    
    TVector3 fnX1 = TVector3 (1,0,0);
    TVector3 fnY1 = TVector3( 0,1,0);
    TVector3 fnZ1 = TVector3( 0,0,1);
    double radiatorL = 489.712; //4*122.5;
    double barend = -294.022; // 4*1225-1960; -294.022
    
    double minChangle = 0.6;
    double maxChangle = 0.9;
    double sum1,sum2,noise = 0.5;
    double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
    double evtime,luttheta,tangle,lenz;
    int64_t pathid;
    TVector3 posInBar,momInBar,momInBar_unit,dir,dird;
    double cherenkovreco[5],spr[5];
    
    TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    TSpectrum *spect = new TSpectrum(10);
    TH1F *hAngle[5], *hLnDiff[5], *hNph[5], *hNph_p[5], *hNph_n[5] ;
    TF1  *fAngle[5];
    double mAngle[5];
    TH1F *hDiff = new TH1F("hDiff",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
    TH1F *hDiffT = new TH1F("hDiffT",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
    TH1F *hDiffD = new TH1F("hDiffD",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
    TH1F *hDiffR = new TH1F("hDiffR",";t_{calc}-t_{measured} [ns];entries [#]", 400,-40,40);
    TH1F *hTime = new TH1F("hTime",";propagation time [ns];entries [#]",   1000,0,200);
    TH1F *hCalc = new TH1F("hCalc",";calculated time [ns];entries [#]",   1000,0,200);
    TH1F *hNphC = new TH1F("hNphC",";detected photons [#];entries [#]", 150,0,150);
    hDiff->SetMinimum(0);
    TGaxis::SetMaxDigits(3);
    double sigma[]={0.01,0.01,0.01,0.010,0.01,0.01};
    for(int i=0; i<5; i++){
        double momentum=3.75;
        hAngle[i]= new TH1F(Form("hAngle_%d",i),  "cherenkov angle;#theta_{C} [rad];entries/N_{max} [#]", 250,0.6,1);
        hNph[i] = new TH1F(Form("hNph_%d",i),";detected photons [#];entries [#]", 150,0,150);
        //hNph_p[i] = new TH1F(Form("hNph_p_%d",i),";detected photons (+) [#];entries [#]", 150,0,150);
        //hNph_n[i] = new TH1F(Form("hNph_n_%d",i),";detected photons (-) [#];entries [#]", 150,0,150);
        mAngle[i] = acos(sqrt(momentum * momentum + glx_mass[i]*glx_mass[i])/momentum/1.473);  //1.4738
        fAngle[i] = new TF1(Form("fAngle_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
        fAngle[i]->SetParameter(0,1);        // const
        fAngle[i]->SetParameter(1,mAngle[i]);// mean
        fAngle[i]->SetParameter(2,sigma[i]);    // sigma
        hAngle[i]->SetMarkerStyle(20);
        hAngle[i]->SetMarkerSize(0.8);
        hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",120,-120,120); // 120,-60,60
    }
    hAngle[2]->SetLineColor(4);
    hAngle[3]->SetLineColor(2);
    hAngle[2]->SetMarkerColor(kBlue+1);
    hAngle[3]->SetMarkerColor(kRed+1);
    fAngle[2]->SetLineColor(4);
    fAngle[3]->SetLineColor(2);
    hLnDiff[2]->SetLineColor(4);
    hLnDiff[3]->SetLineColor(2);
    int evtcount=0,count[5]={0};
    TCanvas *cc = new TCanvas("cc","cc",800,500);
    TLine *gLine = new TLine();
    // cuts
    double cut_cangle=0.04;
    double cut_tdiffd=3;//3;
    double cut_tdiffr=4;//3.5;
    const int nbins=20;
    //TFile file("outFile.root","recreate");
    double inv_mass, missing_mass, chi_square, TofTrackDist;
    TH1F *hist_ev_rho_mass = new TH1F("hist_ev_rho_mass","; #pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.3, 1.2);
    TH1F *hist_ev_phi_mass = new TH1F("hist_ev_phi_mass","; k^{#plus}k^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.9, 1.2);
    TH1F *hist_ev_missing_mass_phi = new TH1F("hist_ev_missing_mass_phi",";#phi Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    TH1F *hist_ev_missing_mass_rho = new TH1F("hist_ev_missing_mass_rho",";#rho Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    TH1F *hist_ev_chi_phi = new TH1F("hist_ev_chi_phi","; #phi Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    TH1F *hist_ev_chi_rho = new TH1F("hist_ev_chi_rho","; #rho Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    
    TH1F *hist_ev_rho_mass_cut = new TH1F("hist_ev_rho_mass_cut","; #pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.3, 1.2);
    TH1F *hist_ev_phi_mass_cut = new TH1F("hist_ev_phi_mass_cut","; k^{#plus}k^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.9, 1.2);
    TH1F *hist_ev_missing_mass_phi_cut = new TH1F("hist_ev_missing_mass_phi_cut",";#phi Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    TH1F *hist_ev_missing_mass_rho_cut = new TH1F("hist_ev_missing_mass_rho_cut",";#rho Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    TH1F *hist_ev_chi_phi_cut = new TH1F("hist_ev_chi_phi_cut","; #phi Kinematic Fit #chi^{2}/NDF ;entries [#]", 100, 0, 45);
    TH1F *hist_ev_chi_rho_cut = new TH1F("hist_ev_chi_rho_cut","; #rho Kinematic Fit #chi^{2}/NDF ;entries [#]", 100, 0, 45);
    
    hist_ev_rho_mass->SetTitle("#rho Invariant Mass");
    hist_ev_phi_mass->SetTitle("#phi Invariant Mass");
    hist_ev_rho_mass_cut->SetTitle("#rho Invariant Mass cut");
    hist_ev_phi_mass_cut->SetTitle("#phi Invariant Mass cut");
    
    TH2F * hExtrapolatedBarHitXY = new TH2F( "hExtrapolatedBarHitXY" , "; Bar Hit X (cm); Bar Hit Y (cm)", 200, -100, 100, 200, -100, 100);
    TH2F * hExtrapolatedBarHitXY_cut = new TH2F( "hExtrapolatedBarHitXY_cut" , "; Bar Hit X (cm); Bar Hit Y (cm)", 200, -100, 100, 200, -100, 100);
    
    //TH2F * mom_theta = new TH2F( "mom_theta" , "; Momentum [Gev/c]; Polar Angle [degree]", 100, 0, 12, 100, 0, 12);
    //TH2F * mom_theta_cut = new TH2F( "mom_theta_cut" , "; Momentum [GeV/c]; Polar Angle [degree]", 100, 0, 12, 100, 0, 12);
    
    TH2F * mom_theta_phi= new TH2F("mom_theta_phi", " ;kaons #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    TH2F * mom_theta_rho= new TH2F("mom_theta_rho", " ;pions #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    
    TH2F * mom_theta_phi_cut= new TH2F("mom_theta_phi_cut", " ;kaons #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    TH2F * mom_theta_rho_cut= new TH2F("mom_theta_rho_cut", " ;pions #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    
    TH1F * hmom_phi = new TH1F("hmom_phi",";kaon Momentum [GeV/c];entries [#]", 100,0,12);
    TH1F * hmom_rho = new TH1F("hmom_rho",";pions Momentum [GeV/c];entries [#]", 100,0,12);
    
    TH1F*  hdir_x = new TH1F("hdir_x",";dir x component ;entries [#]", 100,-1.0,1.0);
    TH1F*  hdir_y = new TH1F("hdir_y",";dir y component ;entries [#]", 100,-1.0,1.0);
    TH1F*  hdir_z = new TH1F("hdir_z",";dir z component;entries [#]", 100,-1.0,1.0);
    
    
    ///////////////////////
    /// cherenkove PDF  ///
    ///////////////////////
    TH1F*  fHistCh_k[glx_nch], *fHistCh_pi[glx_nch], *fHistCh_read_k[glx_nch], *fHistCh_read_pi[glx_nch];
    TFile *ffile_cherenkov_pdf;
    TString cherenkov_pdf_path;
    
    for(Int_t i=0; i<glx_nch; i++) {
        fHistCh_k[i] = new TH1F(Form("fHistCh_k_%d",i),Form("fHistCh_k_%d;#theta_{C} [rad];entries [#]",i), 2000,0.6,1); //2000
        fHistCh_pi[i] = new TH1F(Form("fHistCh_pi_%d",i),Form("fHistCh_pi_%d;#theta_{C} [rad];entries [#]",i), 2000,0.6,1); //2000
    }
    // read pdf
    if (gPDF==2) {
        //cherenkov_data_k_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/pdf/histo_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        cherenkov_pdf_path ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_cherenkovPDF.root";
        cout<<"cherenkov_pdf_path= " <<cherenkov_pdf_path<<endl;
        ffile_cherenkov_pdf  = new TFile(cherenkov_pdf_path, "READ");
        for(Int_t pix=0; pix<glx_nch; pix++) {
            fHistCh_read_k[pix] = (TH1F*)ffile_cherenkov_pdf->Get(Form("fHistCh_k_%d",pix));
            fHistCh_read_pi[pix] = (TH1F*)ffile_cherenkov_pdf->Get(Form("fHistCh_pi_%d",pix));
        }
    }
    ////////////////////////////
    /// cherenkove correction///
    ////////////////////////////
    int PMT_num=108; //108
    double referance_angle = mAngle[2]; // pi
    //double referance_angle = mAngle[3]; // k
    
    double referance_angle_pi = mAngle[2];
    double referance_angle_k = mAngle[3];
    
    
    TH1F*  fHistMCP_k[PMT_num], *fHistMCP_pi[PMT_num], *fHistMCP_read_k[PMT_num], *fHistMCP_read_pi[PMT_num];
    TFile *ffile_cherenkov_correction;
    TString cherenkov_correction_path;
    // correction array
    double array_correction[PMT_num];
    // creat histograms Cherenkov per PMT
    for(Int_t i=0; i<PMT_num; i++) {
        fHistMCP_k[i] = new TH1F(Form("fHistMCP_k_%d",i),Form("fHistMCP_k_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
        fHistMCP_pi[i] = new TH1F(Form("fHistMCP_pi_%d",i),Form("fHistMCP_pi_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
    }
    // Read Cherenkov per PMT
    TF1 *fit_PMT = new TF1("fit_PMT","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    fit_PMT->SetLineColor(kMagenta);
    TH1F*  hdiff_test = new TH1F("hdiff_test",";mean diff ;entries [#]", 100,0,0.1);
    TH1F*  hsigma_test = new TH1F("hsigma_test",";sigam ;entries [#]", 100,0,20);
    TH1F*  hmean_test = new TH1F("hmean_test",";mean ;entries [#]", 250,0.6,1);
    //gStyle->SetPalette(kLightTemperature);
    gStyle->SetPalette(kThermometer);
    glx_canvasAdd("r_pmt_correction",800,400);
    if (gCherenkov_Correction==2) {
        cherenkov_correction_path ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/cherenkov_correction.root";
        cout<<"cherenkov_correction_path= " <<cherenkov_correction_path<<endl;
        ffile_cherenkov_correction  = new TFile(cherenkov_correction_path, "READ");
        for(Int_t PMT=0; PMT<PMT_num; PMT++) {
            if(PMT<=10 || (PMT>=90 && PMT<=96)) continue; // dummy pmts
            fHistMCP_read_k[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistMCP_k_%d",PMT));
            fHistMCP_read_pi[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistMCP_pi_%d",PMT));
            fit_PMT->SetParameters(100,0.82,0.010,10);
            fit_PMT->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fit_PMT->SetParLimits(0,0.1,1E6);
            fit_PMT->SetParLimits(1,0.82-2*cut_cangle,0.82+2*cut_cangle);
            fit_PMT->SetParLimits(2,0.005,0.030); // width
            fHistMCP_read_pi[PMT]->Fit("fit_PMT","M","",0.82-cut_cangle,0.82+ cut_cangle);
            double histo_cor_entries = fHistMCP_read_pi[PMT]->GetEntries();
            double mean_cherenkov_cor=  fit_PMT->GetParameter(1);
            double sigma_cherenkov_cor= fit_PMT->GetParameter(2);
            double delta_cherenkov_cor= referance_angle - mean_cherenkov_cor;
            double val_1 = (delta_cherenkov_cor)*1000.0 ;
            // cout<<"##########"<< fabs(mean_cherenkov_cor - referance_angle) << "  "<<sigma_cherenkov_cor*1000<<endl;
            hsigma_test->Fill(sigma_cherenkov_cor*1000);
            hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
            hmean_test->Fill(mean_cherenkov_cor);
            // correction condition
            if ( fabs(delta_cherenkov_cor) <0.01 && (sigma_cherenkov_cor*1000<12 && sigma_cherenkov_cor*1000> 7) && histo_cor_entries>0 ) {
                array_correction[PMT]= delta_cherenkov_cor;
                cout<<"##########"<< "shift "<<val_1<<endl;
                // Fill PMT
                for(Int_t m=0; m<64; m++)
                    for(Int_t n=0; n<64; n++){
                        glx_hdigi[PMT]->SetBinContent(m, n,val_1);
                    }
                // default correction
                //array_correction[PMT]= -0.004;

                fHistMCP_read_pi[PMT]->Draw();
                glx_canvasGet("r_pmt_correction")->Update();
                TLine *lin_ref = new TLine(0,0,0,1000);
                lin_ref->SetX1(referance_angle);
                lin_ref->SetX2(referance_angle);
                lin_ref->SetY1(gPad->GetUymin());
                lin_ref->SetY2(gPad->GetUymax());
                lin_ref->SetLineColor(kBlue);
                lin_ref->Draw();
                glx_canvasGet("r_pmt_correction")->Update();
                glx_waitPrimitive("r_pmt_correction");
                //hsigma_test->Fill(sigma_cherenkov_cor*1000);
                //hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
                //hmean_test->Fill(mean_cherenkov_cor);
            }
        }
    }
/*
    glx_canvasAdd("r_diff_test",800,400);
    hdiff_test->Draw();
    glx_canvasAdd("r_hsigma_test",800,400);
    hsigma_test->Draw();
    glx_canvasAdd("r_hmean_test",800,400);
    hmean_test->Draw();
    
    double max_digi(10);//30
    double min_digi(-10);//-30
    glx_drawDigi("m,p,v\n",0, max_digi,min_digi);
    glx_canvasSave(0);
    return;
  */  
    
    //////////////////
    /// Reco Method //
    //////////////////
    TString outFile;
    if(gPDF==1){
        outFile= "created_cherenkovPDF.root";
    } else if(gPDF==2){ outFile= "outFile_separation_PDF.root";
    } else if(gCherenkov_Correction==1){ outFile= "cherenkov_correction.root";
    } else {
        outFile= "outFile.root";
    }
    TFile file(outFile,"recreate");
    int pion_counter =0;
    DrcHit hit;
    for (int e = 0; e < glx_ch->GetEntries(); e++){
        glx_ch->GetEntry(e);
        for (int t = 0; t < glx_events->GetEntriesFast(); t++){
            // cut in Event criteria here
            
            
            glx_nextEventc(e,t,1000);
            posInBar = glx_event->GetPosition();
            momInBar = glx_event->GetMomentum();
            double momentum = momInBar.Mag();
            int pdgId = glx_findPdgId(glx_event->GetPdg());
            int bar = glx_event->GetId();
            //if(count[pdgId]>1000) continue;
            
            //std::cout<<"######### No Problem "<<pdgId<<std::endl;
            inv_mass=  glx_event->GetInvMass();
            missing_mass=  glx_event->GetMissMass();
            chi_square=  glx_event->GetChiSq();
            TofTrackDist=  glx_event->GetTofTrackDist();
            
            ////////////////////////
            //////// selection//////
            ////////////////////////
            
            if (pdgId == 2) pion_counter++ ;
            if (pdgId == 2 && pion_counter%20 !=0) continue;
            // pion =2  ,kaon=3
            if(true){
                if (pdgId == 2 && chi_square> 20) continue;
                if (pdgId == 3 && chi_square> 20)continue;
                if (pdgId == 2 && (inv_mass< 0.66  || inv_mass> 0.82 )  ) continue;
                if (pdgId == 3 && (inv_mass< 1.015 || inv_mass> 1.025) ) continue;
                if(missing_mass< -0.01 || missing_mass> 0.01 )continue;
            }
            
            momInBar_unit=momInBar.Unit();
            double dir_x =momInBar_unit.X();
            double dir_y =momInBar_unit.Y();
            double dir_z =momInBar_unit.Z();
            //hdir_x->Fill(dir_x);
            //hdir_y->Fill(dir_y);
            //hdir_z->Fill(dir_z);
            //cout<<"=========>" << dir_x << "  "<< dir_y<< endl;
            hExtrapolatedBarHitXY->Fill(posInBar.X(), posInBar.Y());
            if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest
            if(momInBar.Mag()<3.5 || momInBar.Mag()>4.0 ) continue;
            //if(momInBar.Mag()<2.8 || momInBar.Mag()>3 ) continue;
            //if(momInBar.Mag()<3.9 || momInBar.Mag()>4.1 ) continue;
            int bin = (100+posInBar.X())/200.*nbins;
            // not used
            // if(bar<0 || bar>=luts || (bar!=ybar && ybar!=-1)) continue;
            // if(bin<0 || bin>nbins || (bin!=xbar && xbar!=-1)) continue;
            // commented
            if(bar<0 || bar>=luts || (bar<4 || bar>8)) continue;
            //if(bin<0 || bin>nbins || (bin<7 || bin>13)) continue;
            //std::cout<<"##################### bar "<<bar<<" "<<"####### bin "<<bin<<std::endl;
            if ( posInBar.X()>10 || posInBar.X() < -10 ) continue;
            //if (bar != 6 )continue ;
            hExtrapolatedBarHitXY_cut->Fill(posInBar.X(), posInBar.Y());
            //if(fabs(dir_x)>0.01 )continue;
            //if(dir_y<-0.05)continue;
            hdir_x->Fill(dir_x);
            hdir_y->Fill(dir_y);
            hdir_z->Fill(dir_z);
            
            if (pdgId == 2){
                hist_ev_rho_mass_cut->Fill(inv_mass);
                hist_ev_missing_mass_rho_cut->Fill(missing_mass);
                hist_ev_chi_rho_cut->Fill(chi_square);
                
                if (glx_event->GetPdg() > 0 ) mom_theta_rho_cut->Fill(momInBar.Theta()*180/PI, momInBar.Mag());
                if (glx_event->GetPdg() < 0 ) mom_theta_rho_cut->Fill(-1.0 * momInBar.Theta()*180/PI, momInBar.Mag());
            }
            
            if (pdgId == 3){
                hist_ev_phi_mass_cut->Fill(inv_mass);
                hist_ev_missing_mass_phi_cut->Fill(missing_mass);
                hist_ev_chi_phi_cut->Fill(chi_square);
                
                if (glx_event->GetPdg() > 0 ) mom_theta_phi_cut->Fill(momInBar.Theta()*180/PI, momInBar.Mag());
                if (glx_event->GetPdg() < 0 ) mom_theta_phi_cut->Fill(-1.0 * momInBar.Theta()*180/PI, momInBar.Mag());
            }
            
            
            if(glx_event->GetParent()>0) continue;
            // if(hLnDiff[pdgId]->GetEntries()>200) continue;
            
            for(int p=0; p<5; p++){
                mAngle[p] = acos(sqrt(momentum * momentum + glx_mass[p]*glx_mass[p])/momentum/1.473);  //1.4738
                fAngle[p]->SetParameter(1,mAngle[p]);// mean
            }
            
            sum1=0;
            sum2=0;
            int nph=0;
            int nph_p=0;
            int nph_n=0;
            int nphc=0;
            //      hNphC->Fill(glx_event->GetHitSize());
            
            bool goodevt=0;
            for(int h = 0; h < glx_event->GetHitSize(); h++){
                hit = glx_event->GetHit(h);
                int ch = hit.GetChannel();
                int pmt = hit.GetPmtId();
                int pix = hit.GetPixelId();
                
                double hitTime = hit.GetLeadTime()-glx_event->GetTime();
                
                if(ch>glx_nch) continue;
                //if(hitTime>40) continue;
                nphc++;
                
                
                /////////////////////////////////////
                // Reflection condition may change //
                /////////////////////////////////////
                bool reflected = hitTime>40;
                
                
                lenz = fabs(barend-posInBar.X());
                double rlenz = 2*radiatorL - lenz;
                double dlenz = lenz;
                if(reflected) lenz = 2*radiatorL - lenz;
                
                
                bool isGood(false);
                
                double p1,p2;
                
                for(int i = 0; i < lutNode[bar][ch]->Entries(); i++){
                    dird   = lutNode[bar][ch]->GetEntry(i);
                    evtime = lutNode[bar][ch]->GetTime(i);
                    pathid = lutNode[bar][ch]->GetPathId(i);
                    bool samepath(false);
                    if(fabs(pathid-hit.GetPathId())<0.0001) samepath=true;
                    p1=hit.GetPathId();
                    if(samepath) p2=pathid;
                    //if(!samepath) continue;
                    
                    for(int r=0; r<2; r++){
                        if(!reflected && r==1) continue;
                        
                        if(r) lenz = rlenz;
                        else lenz = dlenz;
                        
                        for(int u = 0; u < 4; u++){
                            if(u == 0) dir = dird;
                            if(u == 1) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z());
                            if(u == 2) dir.SetXYZ( dird.X(), dird.Y(), -dird.Z());
                            if(u == 3) dir.SetXYZ( dird.X(),-dird.Y(), -dird.Z());
                            if(r) dir.SetXYZ( -dir.X(), dir.Y(), dir.Z());
                            if(dir.Angle(fnY1) < criticalAngle || dir.Angle(fnZ1) < criticalAngle) continue;
                            
                            luttheta = dir.Angle(TVector3(-1,0,0));
                            if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
                            //tangle = momInBar.Angle(dir);//-0.004; //correction
                            
                            ////////////////////
                            // PMT Correction //
                            ////////////////////
                            
                            if(gCherenkov_Correction != 2) tangle = momInBar.Angle(dir);
                            if(gCherenkov_Correction == 2) tangle = momInBar.Angle(dir) + array_correction[pmt];
                            
                            //double bartime = lenz/cos(luttheta)/20.4; //198 //203.767 for 1.47125
                            double bartime = lenz/cos(luttheta)/19.6; //203.767 for 1.47125
                            double totalTime = bartime+evtime;
                            // hTime->Fill(hitTime);
                            // hCalc->Fill(totalTime);
                            
                            if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))<cut_cangle){
                                hDiff->Fill(totalTime-hitTime);
                                //if(samepath)
                                {
                                    hDiffT->Fill(totalTime-hitTime);
                                    if(r) hDiffR->Fill(totalTime-hitTime);
                                    else hDiffD->Fill(totalTime-hitTime);
                                }
                            }
                            
                            ///////////////
                            // Time Cut  //
                            ///////////////
                            if(!r && fabs(totalTime-hitTime)>cut_tdiffd) continue;
                            if(r && fabs(totalTime-hitTime) >cut_tdiffr) continue;
                            
                            // cherenkove PDF per PIX
                            if(gPDF ==1 && pdgId == 3) fHistCh_k[ch]->Fill(tangle); // good after time cut
                            if(gPDF ==1 && pdgId == 2) fHistCh_pi[ch]->Fill(tangle); // good after time cut
                            // cherenkove correction per PMT
                            if(gCherenkov_Correction ==1 && pdgId == 3) fHistMCP_k[pmt]->Fill(tangle);
                            if(gCherenkov_Correction ==1 && pdgId == 2) fHistMCP_pi[pmt]->Fill(tangle);
                            
                            // fill cherenkove histo
                            hAngle[pdgId]->Fill(tangle);
                            
                            ////////////////////
                            // Cherenkov Cut  //
                            ////////////////////
                            //if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))>cut_cangle) continue;
                            
                            if(fabs(tangle-0.5*(referance_angle_k+referance_angle_pi))>cut_cangle) continue;
                            
                            isGood=true;
                            hTime->Fill(hitTime);
                            hCalc->Fill(totalTime);
                            
                            
                            if(gPDF ==0 || gPDF ==1 ){ //continue; //cut_cangle  0.2
                                sum1 += TMath::Log(fAngle[2]->Eval(tangle)+noise);
                                sum2 += TMath::Log(fAngle[3]->Eval(tangle)+noise);
                            }
                            
                            if(gPDF ==2){
                                // use histograms
                                Int_t kk = fHistCh_read_k[ch]->GetXaxis()->FindBin(tangle);
                                Int_t kpi = fHistCh_read_pi[ch]->GetXaxis()->FindBin(tangle);
                                if (fHistCh_read_pi[ch]->GetBinContent(kpi) > 0 )sum1 += TMath::Log(fHistCh_read_pi[ch]->GetBinContent(kpi));
                                if (fHistCh_read_k[ch]->GetBinContent(kk) > 0 )sum2 += TMath::Log(fHistCh_read_k[ch]->GetBinContent(kk));
                                
                                //if (sum1 != 0 || sum2!=0 )std::cout<<"No Problem  separation  " <<kpi<<" "<<kk<<"  sum "<<sum1 <<"  "<< sum2<<std::endl;
                                //std::cout<<"###### No Problem  separation  " << fHistCh_read_k[ch]->GetBinContent(kk) <<"  "<< fHistCh_read_pi[ch]->GetBinContent(kpi)<<std::endl;
                            }
                            
                            
                            if(0){
                                TString x=(sum1>sum2)? " <====== PION" : "";
                                std::cout<<Form("%1.6f  %1.6f | %1.6f  %1.6f        pid %d",TMath::Log(fAngle[2]->Eval(tangle)+noise),TMath::Log(fAngle[3]->Eval(tangle)+noise), sum1, sum2,pdgId)<<"  " <<std::endl;
                                
                                cc->cd();
                                fAngle[2]->Draw("");
                                fAngle[3]->Draw("same");
                                
                                cc->Update();
                                gLine->SetLineWidth(2);
                                gLine->SetX1(tangle);
                                gLine->SetX2(tangle);
                                gLine->SetY1(cc->GetUymin());
                                gLine->SetY2(cc->GetUymax());
                                gLine->Draw();
                                cc->Update();
                                cc->WaitPrimitive();
                            }
                            
                        } // bar ambiguities
                    } // reflection loop
                } // LUT loop
                
                if(isGood){
                    nph++;
                    //if (glx_event->GetPdg() > 0 ) nph_p++;
                    //if (glx_event->GetPdg() < 0 ) nph_n++;
                    if(pmt<108) {
                        glx_hdigi[pmt]->Fill(pix%8, pix/8);
                        goodevt=1;
                    }
                }
            } // hit loop
            
            if(goodevt) evtcount++;
            if(nph<5) continue;
            hNph[pdgId]->Fill(nph);
            
            //hNph_p[pdgId]->Fill(nph_p);
            //hNph_n[pdgId]->Fill(nph_n);
            
            hNphC->Fill(nphc);
            
            double sum = sum1-sum2;
            //cout<<"########### sum  "<<sum <<"  sum1  "<<sum1<<"  sum2  "<<sum2<<endl;
            hLnDiff[pdgId]->Fill(sum);
            
            count[pdgId]++;
            
            if(0 && pdgId==3){
                //	if(!cc)
                TString x=(sum1>sum2)? " <====== Pion" : "";
                // std::cout<<Form("f %1.6f s %1.6f PMT %d pix %d   pid %d",aminf,amins,PMT,pix  ,prt_particle)<<"  "<<x <<std::endl;
                
                std::cout<<"PID "<< glx_event->GetPdg() <<" sum1 "<<sum1<<" sum2 "<<sum2<<" sum "<<sum<<" "<<x<<std::endl;
                
                cc->cd();
                
                if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
                if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
                
                hAngle[2]->Draw("hist");
                hAngle[3]->Draw("hist same");
                fAngle[2]->Draw("same");
                fAngle[3]->Draw("same");
                
                // hAngle[2]->GetYaxis()->SetRangeUser(0,20);
                // hAngle[3]->GetYaxis()->SetRangeUser(0,20);
                
                cc->Update();
                TLine *line = new TLine(0,0,0,1000);
                line->SetX1(mAngle[3]);
                line->SetX2(mAngle[3]);
                line->SetY1(cc->GetUymin());
                line->SetY2(cc->GetUymax());
                line->SetLineColor(kRed);
                line->Draw();
                
                TLine *line2 = new TLine(0,0,0,1000);
                line2->SetX1(mAngle[2]);
                line2->SetX2(mAngle[2]);
                line2->SetY1(cc->GetUymin());
                line2->SetY2(cc->GetUymax());
                line2->SetLineColor(kBlue);
                line2->Draw();
                
                cc->Update();
                cc->WaitPrimitive();
            }
            
            // hAngle[2]->Reset();
            // hAngle[3]->Reset();
            
        }
    } // Event Loop
    
    if(evtcount>0){
        for(int i=0; i<glx_nch; i++){
            int pmt=i/64;
            int pix=i%64;
            double rel = glx_hdigi[pmt]->GetBinContent(pix%8+1,pix/8+1)/(double)evtcount;
            glx_hdigi[pmt]->SetBinContent(pix%8+1, pix/8+1,rel);
        }
    }
    
    //TString nid=Form("_%2.2f_%2.2f",theta,phi);
    TString nid=Form("_%d_%d",xbar,ybar);
    
    glx_drawDigi("m,p,v\n",0);
    glx_cdigi->SetName("hp"+nid);
    glx_canvasAdd(glx_cdigi);
    
    glx_canvasAdd("hAngle"+nid,800,400);
    
    //scal
    if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
    if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
    
    for(int i=0; i<5; i++){
        if(hAngle[i]->GetEntries()<20) continue;
        
        int nfound = spect->Search(hAngle[i],1,"goff",0.9);
        if(nfound>0) cherenkovreco[i] = spect->GetPositionX()[0];
        else cherenkovreco[i] =  hAngle[i]->GetXaxis()->GetBinCenter(hAngle[i]->GetMaximumBin());
        if(cherenkovreco[i]>0.85) cherenkovreco[i]=0.82;
        
        if(i==2)  fit->SetLineColor(kBlue);
        if(i==3)  fit->SetLineColor(kRed);
        fit->SetParameters(100,cherenkovreco[i],0.010,10);
        fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
        fit->SetParLimits(0,0.1,1E6);
        fit->SetParLimits(1,cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);
        fit->SetParLimits(2,0.005,0.030); // width
        hAngle[i]->Fit("fgaus","I","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
        hAngle[i]->Fit("fgaus","M","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
        
        cherenkovreco[i] = fit->GetParameter(1);
        spr[i] = fit->GetParameter(2);
    }
    
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    
    
    
    TF1 *ff;
    double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0;
    if(hLnDiff[3]->GetEntries()>100){
        hLnDiff[3]->Fit("gaus","S");
        ff = hLnDiff[3]->GetFunction("gaus");
        ff->SetLineColor(1);
        m1=ff->GetParameter(1);
        s1=ff->GetParameter(2);
    }
    if(hLnDiff[2]->GetEntries()>100){
        hLnDiff[2]->Fit("gaus","S");
        ff = hLnDiff[2]->GetFunction("gaus");
        ff->SetLineColor(1);
        m2=ff->GetParameter(1);
        s2=ff->GetParameter(2);
    }
    if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
    
    hAngle[2]->GetXaxis()->SetRangeUser(0.7,0.9);
    hAngle[2]->GetYaxis()->SetRangeUser(0,1.2);
    hAngle[2]->Draw();
    hAngle[3]->Draw("same");
    // fAngle[3]->Draw("same");
    // fAngle[2]->Draw("same");
    
    
    TLine *line = new TLine(0,0,0,1000);
    //line->SetX1(mAngle[3]);
    //line->SetX2(mAngle[3]);
    line->SetX1(referance_angle_k);
    line->SetX2(referance_angle_k);
    line->SetY1(0);
    line->SetY2(1.2);
    line->SetLineColor(kRed);
    line->Draw();
    
    TLine *line2 = new TLine(0,0,0,1000);
    //line2->SetX1(mAngle[2]);
    //line2->SetX2(mAngle[2]);
    line2->SetX1(referance_angle_pi);
    line2->SetX2(referance_angle_pi);
    line2->SetY1(0);
    line2->SetY2(1.2);
    line2->SetLineColor(kBlue);
    line2->Draw();
    
    TLine *line3 = new TLine(0,0,0,1000);
    line3->SetLineStyle(2);
    //line3->SetX1(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
    //line3->SetX2(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
    line3->SetX1(0.5*(referance_angle_k+referance_angle_pi)-cut_cangle);
    line3->SetX2(0.5*(referance_angle_k+referance_angle_pi)-cut_cangle);
    line3->SetY1(0);
    line3->SetY2(1.2);
    line3->SetLineColor(1);
    line3->Draw();
    
    TLine *line4 = new TLine(0,0,0,1000);
    line4->SetLineStyle(2);
    //line4->SetX1(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
    //line4->SetX2(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
    line4->SetX1(0.5*(referance_angle_k+referance_angle_pi)+cut_cangle);
    line4->SetX2(0.5*(referance_angle_k+referance_angle_pi)+cut_cangle);
    line4->SetY1(0);
    line4->SetY2(1.2);
    line4->SetLineColor(1);
    line4->Draw();
    
    
    TLegend *leg = new TLegend(0.1,0.5,0.4,0.85);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hAngle[2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
    leg->AddEntry(hAngle[3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
    leg->AddEntry(hAngle[2],Form("#sigma_{c}^{#pi} = %2.1f mrad",spr[2]*1000),"");
    leg->AddEntry(hAngle[3],Form("#sigma_{c}^{K} = %2.1f mrad",spr[3]*1000),"");
    leg->Draw();
    
    TLegend *lnpa = new TLegend(0.7,0.67,0.9,0.85);
    lnpa->SetFillColor(0);
    lnpa->SetFillStyle(0);
    lnpa->SetBorderSize(0);
    lnpa->SetFillStyle(0);
    lnpa->AddEntry(hAngle[2],"pions","lp");
    lnpa->AddEntry(hAngle[3],"kaons","lp");
    lnpa->Draw();
    
    // fAngle[2]->Draw("same");
    // fAngle[3]->Draw("same");
    
    glx_canvasAdd("hTime"+nid,800,400);
    
    hTime->Draw();
    hCalc->SetLineColor(2);
    hCalc->Draw("same");
    TLegend *leg1 = new TLegend(0.5,0.6,0.85,0.80);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(hTime,"measured in geant","lp");
    leg1->AddEntry(hCalc,"calculated","lp");
    leg1->Draw();
    
    glx_canvasAdd("hDiff"+nid,800,400);
    hDiff->SetLineColor(kBlack);
    hDiff->Draw();
    
    // hDiffT->SetLineColor(kRed+1);
    // hDiffT->Draw("same");
    hDiffD->SetLineColor(kGreen+2);
    hDiffD->Draw("same");
    hDiffR->SetLineColor(kBlue+1);
    hDiffR->Draw("same");
    
    double maxTD= hDiffD->GetXaxis()->GetBinCenter(hDiffD->GetMaximumBin());
    double maxTR= hDiffR->GetXaxis()->GetBinCenter(hDiffR->GetMaximumBin());
    double maxTT= hTime->GetXaxis()->GetBinCenter(hTime->GetMaximumBin());
    
    line = new TLine(0,0,0,1000);
    line->SetLineStyle(2);
    line->SetX1(-cut_tdiffd);
    line->SetX2(-cut_tdiffd);
    line->SetY1(0);
    line->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
    line->SetLineColor(1);
    line->Draw();
    
    line2 = new TLine(0,0,0,1000);
    line2->SetLineStyle(2);
    line2->SetX1(cut_tdiffd);
    line2->SetX2(cut_tdiffd);
    line2->SetY1(0);
    line2->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
    line2->SetLineColor(1);
    line2->Draw();
    
    TLegend *leg2 = new TLegend(0.6,0.57,0.9,0.85);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(hDiff,"all","lp");
    // leg2->AddEntry(hDiffT,"MC path in EV","lp");
    // leg2->AddEntry(hDiffD,"MC path in EV for direct photons","lp");
    // leg2->AddEntry(hDiffR,"MC path in EV for reflected photons","lp");
    leg2->AddEntry(hDiffD,"direct photons","lp");
    leg2->AddEntry(hDiffR,"reflected photons","lp");
    
    leg2->Draw();
    
    glx_canvasAdd("hLnDiff"+nid,800,400);
    hLnDiff[2]->SetTitle(Form("sep = %2.2f s.d.",sep));
    hLnDiff[2]->Draw();
    hLnDiff[3]->Draw("same");
    
    TLegend *lnpl = new TLegend(0.7,0.67,0.9,0.85);
    lnpl->SetFillColor(0);
    lnpl->SetFillStyle(0);
    lnpl->SetBorderSize(0);
    lnpl->SetFillStyle(0);
    lnpl->AddEntry(hLnDiff[2],"pions","lp");
    lnpl->AddEntry(hLnDiff[3],"kaons","lp");
    lnpl->Draw();
    
    glx_canvasAdd("hNph"+nid,800,400);
    
    double nph = 0;
    if(hNph[2]->GetEntries()>50){
        nph = glx_fit(hNph[2],40,100,40).X();
        auto rfit = hNph[2]->GetFunction("glx_gaust");
        if(rfit) rfit->SetLineColor(kBlue);
        hNph[2]->SetLineColor(kBlue);
        hNph[2]->Draw();
        //glx_fit(hNph[3],40,100,40).X();
        //hNph[3]->GetFunction("glx_gaust")->SetLineColor(kRed);
        hNph[3]->SetLineColor(kRed);
        hNph[3]->Draw("same");
    }
    
    
    
    
    //////
    
    /*
     
     glx_canvasAdd("hNph_p"+nid,800,400);
     
     double nph_p = 0;
     if(hNph_p[2]->GetEntries()>50){
     nph_p = glx_fit(hNph_p[2],40,100,40).X();
     auto rfit = hNph_p[2]->GetFunction("glx_gaust");
     if(rfit) rfit->SetLineColor(kBlue);
     hNph_p[2]->SetLineColor(kBlue);
     hNph_p[2]->Draw();
     //glx_fit(hNph[3],40,100,40).X();
     //hNph[3]->GetFunction("glx_gaust")->SetLineColor(kRed);
     hNph_p[3]->SetLineColor(kRed);
     hNph_p[3]->Draw("same");
     }
     
     
     glx_canvasAdd("hNph_n"+nid,800,400);
     
     double nph_n = 0;
     if(hNph_n[2]->GetEntries()>50){
     nph_n = glx_fit(hNph_n[2],40,100,40).X();
     auto rfit = hNph_n[2]->GetFunction("glx_gaust");
     if(rfit) rfit->SetLineColor(kBlue);
     hNph_n[2]->SetLineColor(kBlue);
     hNph_n[2]->Draw();
     //glx_fit(hNph[3],40,100,40).X();
     //hNph[3]->GetFunction("glx_gaust")->SetLineColor(kRed);
     hNph_n[3]->SetLineColor(kRed);
     hNph_n[3]->Draw("same");
     }
     
     */
    
    // hNphC->SetLineColor(kBlack);
    // hNphC->Draw("same");
    
    
    TLegend *lnph = new TLegend(0.6,0.65,0.9,0.85);
    lnph->SetFillColor(0);
    lnph->SetFillStyle(0);
    lnph->SetBorderSize(0);
    lnph->SetFillStyle(0);
    // lnph->AddEntry(hNphC,"simulated","lp");
    lnph->AddEntry(hNph[2],"pions","lp");
    lnph->AddEntry(hNph[3],"kaons","lp");
    
    
    //lnph->AddEntry(hNph_p[2],"pions","lp");
    //lnph->AddEntry(hNph_p[3],"kaons","lp");
    
    //lnph->AddEntry(hNph_n[2],"pions","lp");
    //lnph->AddEntry(hNph_n[3],"kaons","lp");
    
    lnph->Draw();
    
    std::cout<<" ###### separation = "<< sep << "  nph = "<<nph <<std::endl;
    std::cout<<"maxTD "<<maxTD<<"  maxTR "<<maxTR<<std::endl;
    
    
    
    
    
    //TFile fc(infile+"_res"+nid+".root","recreate");
    // TFile fc("data/reco_lut_res/res"+nid+".root","recreate");
    // TTree *tc = new TTree("reco","reco");
    // // tc->Branch("theta",&theta,"theta/D");
    // // tc->Branch("phi",&phi,"prt_phi/D");
    // tc->Branch("sep",&sep,"sep/D");
    // tc->Branch("esep",&esep,"esep/D");
    // tc->Branch("moms",&moms,"prt_mom/D");
    // tc->Branch("xbar",&xbar,"xbar/D");
    // tc->Branch("ybar",&ybar,"ybar/D");
    // tc->Branch("nph",&nph,"nph/D");
    // tc->Branch("spr",&spr[3],"spr/D");
    // tc->Branch("maxTD",&maxTD,"maxTD/D");
    // tc->Branch("maxTR",&maxTR,"maxTR/D");
    // tc->Branch("maxTT",&maxTT,"maxTT/D");
    // tc->Fill();
    // tc->Write();
    // fc.Write();
    // fc.Close();
    
    /*
     
     glx_canvasAdd("1",800,400);
     hist_ev_rho_mass_cut->Draw();
     glx_canvasAdd("2",800,400);
     hist_ev_missing_mass_rho_cut->Draw();
     glx_canvasAdd("3",800,400);
     hist_ev_chi_rho_cut->Draw();
     
     glx_canvasAdd("4",800,400);
     hist_ev_phi_mass_cut->Draw();
     glx_canvasAdd("5",800,400);
     hist_ev_missing_mass_phi_cut->Draw();
     glx_canvasAdd("6",800,400);
     hist_ev_chi_phi_cut->Draw();
     
     glx_canvasAdd("7",800,400);
     hist_ev_rho_mass->Draw();
     glx_canvasAdd("8",800,400);
     hist_ev_missing_mass_rho->Draw();
     glx_canvasAdd("9",800,400);
     hist_ev_chi_rho->Draw();
     
     glx_canvasAdd("10",800,400);
     hist_ev_phi_mass->Draw();
     glx_canvasAdd("11",800,400);
     hist_ev_missing_mass_phi->Draw();
     glx_canvasAdd("12",800,400);
     hist_ev_chi_phi->Draw();
     
     
     glx_canvasAdd("13",800,400);
     hExtrapolatedBarHitXY->Draw("colz");
     
     glx_canvasAdd("14",800,400);
     hExtrapolatedBarHitXY_cut->Draw("colz");
     
     
     
     glx_canvasAdd("17",800,400);
     hmom_phi->Draw();
     
     glx_canvasAdd("18",800,400);
     hmom_rho->Draw();
     
     glx_canvasAdd("19",800,400);
     mom_theta_phi->Draw("colz");
     
     glx_canvasAdd("20",800,400);
     mom_theta_rho->Draw("colz");
     
     glx_canvasAdd("21",800,400);
     mom_theta_phi_cut->Draw("colz");
     
     glx_canvasAdd("22",800,400);
     mom_theta_rho_cut->Draw("colz");
     
     glx_canvasAdd("23",800,400);
     hdir_x->Draw();
     glx_canvasAdd("24",800,400);
     hdir_y->Draw();
     glx_canvasAdd("25",800,400);
     hdir_z->Draw();
     
     */
    
    
    
    
    glx_canvasSave(0);
    
    
    if(gPDF ==1) {
        for(Int_t i=0; i<glx_nch; i++) {
            fHistCh_k[i]->Write();
            fHistCh_pi[i]->Write();
        }
    }
    if(gCherenkov_Correction==1) {
        for(Int_t i=0; i<PMT_num; i++) {
            fHistMCP_k[i]->Write();
            fHistMCP_pi[i]->Write();
        }
    }
    
    mom_theta_phi->Write();
    file.Write();
}
