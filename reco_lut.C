#define glx__sim

//#include "/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/classes/DrcEvent.h"

#include "/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/classes_v2/DrcEvent.h"
#include "/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/classes_v2/DrcHit.h"
#include "/lustre/nyx/panda/aali/gluex/gluex_top/halld_recon/halld_recon-4.2.0/src/plugins/Analysis/lut_dirc/DrcLutNode.h"

#include "TMultiGraph.h"
#include "TGraph.h"
#include <vector>

#include "glxtools.C"
#define PI 3.14159265


// gPDF_pix = 1 Create PDF per pix
// gPDF_pix = 2 Calculate Sepration from PDF

// gPDF_pmt = 1 Create PDF per pmt
// gPDF_pmt = 2 Calculate Sepration from pmt PDF

// gCherenkov_Correction = 1 Create histo per PMT
// gCherenkov_Correction = 2 Apply per PMT correction

void reco_lut(TString infile="vol/tree_060772.root", TString inlut="lut/lut_12/lut_all_avr.root", TString justName="NoName.root", int gPDF_pix=0,int gPDF_pmt=0, int gCherenkov_Correction=0, int xbar=-1, int ybar=-1, double moms=3.75){
    if(!glx_initc(infile,1,"data/reco_lut_sim")) return;
    
    
    int num_events=glx_ch->GetEntries();
    cout<<"####### num_events ="<<num_events<<endl;
    
    double momentum=3.75;
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
    
    
    
    TF1 *fit_track = new TF1("fit_track","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    fit_track->SetNpx(1000);
    
    TF1 *fit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    TSpectrum *spect = new TSpectrum(10);
    TH1F *hAngle[5], *hLnDiff[5], *hNph[5], *hNph_p[5], *hNph_n[5] ;
    TF1  *fAngle[5];
    double mAngle[5];
    
    TH1F * histo_cherenkov_track = new TH1F("histo_cherenkov_track","histo_cherenkov_track", 100,0.6,1); //250
    
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
        hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",110,-120,120); //,80,-150,150);//
        //hLnDiff[i] = new TH1F(Form("hLnDiff_%d",i),";ln L(#pi) - ln L(K);entries [#]",80,-150,150);//
        
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
    //TCanvas *cc = new TCanvas("cc","cc",800,500);
    TLine *gLine = new TLine();
    // cuts
    double cut_cangle=0.04;
    double cut_tdiffd=3;//3;
    double cut_tdiffr=3.5;//3.5;
    const int nbins=20;
    
    double xmin(0), xmax(0), ymin(0), ymax(0) ;
    
    //TFile file("outFile.root","recreate");
    
    //    TH1F *hist_ev_rho_mass = new TH1F("hist_ev_rho_mass","; #pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.3, 1.2);
    //    TH1F *hist_ev_phi_mass = new TH1F("hist_ev_phi_mass","; k^{#plus}k^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.9, 1.2);
    //    TH1F *hist_ev_missing_mass_phi = new TH1F("hist_ev_missing_mass_phi",";#phi Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    //    TH1F *hist_ev_missing_mass_rho = new TH1F("hist_ev_missing_mass_rho",";#rho Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    //    TH1F *hist_ev_chi_phi = new TH1F("hist_ev_chi_phi","; #phi Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    //    TH1F *hist_ev_chi_rho = new TH1F("hist_ev_chi_rho","; #rho Kinematic Fit #chi^{2} ;entries [#]", 100, 0, 45);
    //
    //    TH1F *hist_ev_rho_mass_cut = new TH1F("hist_ev_rho_mass_cut","; #pi^{#plus}#pi^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.3, 1.2);
    //    TH1F *hist_ev_phi_mass_cut = new TH1F("hist_ev_phi_mass_cut","; k^{#plus}k^{#minus} Invariant Mass [GeV/c^{2}];entries [#]", 900, 0.9, 1.2);
    //    TH1F *hist_ev_missing_mass_phi_cut = new TH1F("hist_ev_missing_mass_phi_cut",";#phi Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    //    TH1F *hist_ev_missing_mass_rho_cut = new TH1F("hist_ev_missing_mass_rho_cut",";#rho Missing Mass Squared (GeV/c^{2})^{2};entries [#]", 1000, -0.2, 0.2);
    //    TH1F *hist_ev_chi_phi_cut = new TH1F("hist_ev_chi_phi_cut","; #phi Kinematic Fit #chi^{2}/NDF ;entries [#]", 100, 0, 45);
    //    TH1F *hist_ev_chi_rho_cut = new TH1F("hist_ev_chi_rho_cut","; #rho Kinematic Fit #chi^{2}/NDF ;entries [#]", 100, 0, 45);
    //
    //    hist_ev_rho_mass->SetTitle("#rho Invariant Mass");
    //    hist_ev_phi_mass->SetTitle("#phi Invariant Mass");
    //    hist_ev_rho_mass_cut->SetTitle("#rho Invariant Mass cut");
    //    hist_ev_phi_mass_cut->SetTitle("#phi Invariant Mass cut");
    
    
    
    
    //TH2F * mom_theta = new TH2F( "mom_theta" , "; Momentum [Gev/c]; Polar Angle [degree]", 100, 0, 12, 100, 0, 12);
    //TH2F * mom_theta_cut = new TH2F( "mom_theta_cut" , "; Momentum [GeV/c]; Polar Angle [degree]", 100, 0, 12, 100, 0, 12);
    
    //    TH2F * mom_theta_phi= new TH2F("mom_theta_phi", " ;kaons #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    //    TH2F * mom_theta_rho= new TH2F("mom_theta_rho", " ;pions #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    //
    //    TH2F * mom_theta_phi_cut= new TH2F("mom_theta_phi_cut", " ;kaons #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    //    TH2F * mom_theta_rho_cut= new TH2F("mom_theta_rho_cut", " ;pions #theta X charge (deg); #p [GeV/c]", 200, -12, 12, 200, 0, 12);
    //
    //    TH1F * hmom_phi = new TH1F("hmom_phi",";kaon Momentum [GeV/c];entries [#]", 100,0,12);
    //    TH1F * hmom_rho = new TH1F("hmom_rho",";pions Momentum [GeV/c];entries [#]", 100,0,12);
    //
    //    TH1F*  hdir_x = new TH1F("hdir_x",";dir x component ;entries [#]", 100,-1.0,1.0);
    //    TH1F*  hdir_y = new TH1F("hdir_y",";dir y component ;entries [#]", 100,-1.0,1.0);
    //    TH1F*  hdir_z = new TH1F("hdir_z",";dir z component;entries [#]", 100,-1.0,1.0);
    //
    
    
    
    //    TH2F * hExtrapolatedBarHitXY_k = new TH2F( "hExtrapolatedBarHitXY_k" , "; Bar Hit X + (100* charge) (cm); Bar Hit Y (cm)", 200, -200, 200, 200, -200, 200);
    //    TH2F * hExtrapolatedBarHitXY_pi = new TH2F( "hExtrapolatedBarHitXY_pi" , "; Bar Hit X + (100* charge) (cm); Bar Hit Y (cm)", 200, -200, 200, 200, -200, 200);
    //
    //
    //    TH2F * hExtrapolatedBarHitXY_cut = new TH2F( "hExtrapolatedBarHitXY_cut" , "; Bar Hit X + (100* charge) (cm); Bar Hit Y (cm)", 200, -200, 200, 200, -200, 200);
    //
    //    Int_t nf= 121;
    //    TH2F * histo_theta_phi_map_pi =  new TH2F("histo_theta_phi_map_pi",";#Theta multiplied by charge [Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    //    TH2F * histo_theta_phi_map_k =  new TH2F("histo_theta_phi_map_k",";#Theta multiplied by charge [Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    //
    //    TH2F * histo_theta_phi_mom_map_pi =  new TH2F("histo_theta_phi_mom_map_pi ",";#Theta multiplied by charge [Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    //    TH2F * histo_theta_phi_mom_map_k  =  new TH2F("histo_theta_phi_mom_map_k ",";#Theta multiplied by charge [Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    //
    //
    //    TH2F * histo_theta_phi_mom_tmp_map_pi =  new TH2F("histo_theta_phi_mom_tmp_map_pi ",";#Theta[Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    //    TH2F * histo_theta_phi_mom_tmp_map_k  =  new TH2F("histo_theta_phi_mom_tmp_map_k ",";#Theta[Degree]; #Phi[Degree]", nf, -12, 12, nf, -180 , 20);
    
    
    TH1F*  histo_time_bar_pos[24][40];
    TH1F*  histo_tdiffD_bar_pos[24][40];
    TH1F*  histo_tdiffR_bar_pos[24][40];
    
    TH1F* histo_tmp_pos = new TH1F("histo_tmp_pos","; X Bar Hit [cm]; entries [#]",40,-100,100);
    
    for(Int_t j=0; j<40; j++){
      for(Int_t i=0; i<24; i++){
            histo_time_bar_pos[i][j] = new TH1F(Form("histo_time_bar_pos_%d_%d",i,j),Form("histo_time_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 100,0,100);
            histo_tdiffD_bar_pos[i][j] = new TH1F(Form("histo_tdiffD_bar_pos_%d_%d",i,j),Form("histo_tdiffD_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 200,-50,50);
            histo_tdiffR_bar_pos[i][j] = new TH1F(Form("histo_tdiffR_bar_pos_%d_%d",i,j),Form("histo_tdiffR_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 200,-50,50);
        }
    }
    //
    //    //////////////////////////////
    //    /// cherenkove PDF per pix ///
    //    //////////////////////////////
    //    TH1F*  fHistCh_k[glx_nch], *fHistCh_pi[glx_nch], *fHistCh_read_k[glx_nch], *fHistCh_read_pi[glx_nch];
    //    TFile *ffile_cherenkov_pdf_pix;
    //    TString cherenkov_pdf_path_pix;
    //
    //    for(Int_t i=0; i<glx_nch; i++) {
    //        fHistCh_k[i] = new TH1F(Form("fHistCh_k_%d",i),Form("fHistCh_k_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); //2000
    //        fHistCh_pi[i] = new TH1F(Form("fHistCh_pi_%d",i),Form("fHistCh_pi_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); //2000
    //    }
    //    // read pdf per pix
    //    if (gPDF_pix==2) {
    //        //cherenkov_data_k_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/pdf/histo_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
    //        cherenkov_pdf_path_pix ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_cherenkovPDF_pix.root";
    //        cout<<"cherenkov_pdf_path_pix= " <<cherenkov_pdf_path_pix<<endl;
    //        ffile_cherenkov_pdf_pix  = new TFile(cherenkov_pdf_path_pix, "READ");
    //        for(Int_t pix=0; pix<glx_nch; pix++) {
    //            fHistCh_read_k[pix] = (TH1F*)ffile_cherenkov_pdf_pix->Get(Form("fHistCh_k_%d",pix));
    //            fHistCh_read_pi[pix] = (TH1F*)ffile_cherenkov_pdf_pix->Get(Form("fHistCh_pi_%d",pix));
    //        }
    //    }
    //
    //    //////////////////////////////
    //    /// cherenkove PDF per pmt ///
    //    //////////////////////////////
    //    int PMT_num=108; //108
    //    //double norm = 20000;
    //
    //    TH1F*  fHistPMT_PDF_k[PMT_num], *fHistPMT_PDF_pi[PMT_num], *fHistPMT_PDF_read_k[PMT_num], *fHistPMT_PDF_read_pi[PMT_num];
    //    TFile *ffile_cherenkov_pdf_pmt;
    //    TString cherenkov_pdf_path_pmt;
    //
    //    for(Int_t i=0; i<PMT_num; i++) {
    //        fHistPMT_PDF_k[i] = new TH1F(Form("fHistPMT_PDF_k_%d",i),Form("fHistPMT_PDF_k_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1);    // 250
    //        fHistPMT_PDF_pi[i] = new TH1F(Form("fHistPMT_PDF_pi_%d",i),Form("fHistPMT_PDF_pi_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); // 250
    //    }
    //    // read pdf per pmt
    //    if (gPDF_pmt==2) {
    //        //cherenkov_pdf_path_pmt ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_pdf_pmt_corrected.root";
    //        cherenkov_pdf_path_pmt ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_cherenkovPDF_pmt_NotCorrected.root";
    //        cout<<"cherenkov_pdf_path_pmt= " <<cherenkov_pdf_path_pmt<<endl;
    //        ffile_cherenkov_pdf_pmt  = new TFile(cherenkov_pdf_path_pmt, "READ");
    //        for(Int_t PMT=0; PMT<PMT_num; PMT++) {
    //            fHistPMT_PDF_read_k[PMT] = (TH1F*)ffile_cherenkov_pdf_pmt->Get(Form("fHistPMT_PDF_k_%d",PMT));
    //            fHistPMT_PDF_read_pi[PMT] = (TH1F*)ffile_cherenkov_pdf_pmt->Get(Form("fHistPMT_PDF_pi_%d",PMT));
    //
    //            // No success
    //            /*
    //             double scale_k = norm/(fHistPMT_PDF_read_k[PMT]->Integral());
    //             fHistPMT_PDF_read_k[PMT]->Scale(scale_k);
    //             double scale_pi = norm/(fHistPMT_PDF_read_pi[PMT]->Integral());
    //             fHistPMT_PDF_read_k[PMT]->Scale(scale_pi);
    //             */
    //        }
    //    }
    //
    //
    //    ////////////////////////////
    //    /// cherenkove correction///
    //    ////////////////////////////
    //    double referance_angle = mAngle[2]; // pi
    //    //double referance_angle = mAngle[3]; // k
    double referance_angle_pi = mAngle[2];
    double referance_angle_k = mAngle[3];
    //    TGraph *shifted_pi = new TGraph();
    //
    //    shifted_pi->SetMarkerColor(kBlue);
    //    shifted_pi->SetMarkerStyle(20);
    //    shifted_pi->SetLineColor(kBlue);
    //    shifted_pi->SetLineWidth(1);
    //
    //    TH1F*  fHistPMT_k[PMT_num], *fHistPMT_pi[PMT_num], *fHistPMT_read_k[PMT_num], *fHistPMT_read_pi[PMT_num];
    //    TFile *ffile_cherenkov_correction;
    //    TString cherenkov_correction_path;
    //    // correction array
    //    double array_correction[108]={0};
    //    // creat histograms Cherenkov per PMT
    //    for(Int_t i=0; i<PMT_num; i++) {
    //        fHistPMT_k[i] = new TH1F(Form("fHistPMT_k_%d",i),Form("fHistPMT_k_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
    //        fHistPMT_pi[i] = new TH1F(Form("fHistPMT_pi_%d",i),Form("fHistPMT_pi_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
    //    }
    //    // Read Cherenkov per PMT
    //    TF1 *fit_PMT = new TF1("fit_PMT","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    //    fit_PMT->SetLineColor(kBlue);
    //    TH1F*  hdiff_test = new TH1F("hdiff_test",";mean diff ;entries [#]", 100,0,0.1);
    //    TH1F*  hsigma_test = new TH1F("hsigma_test",";sigam ;entries [#]", 100,0,20);
    //    TH1F*  hmean_test = new TH1F("hmean_test",";mean ;entries [#]", 250,0.6,1);
    //    //gStyle->SetPalette(kLightTemperature);
    //    //////////////////////////////////////////////////
    //    //    glx_canvasAdd("r_pmt_correction",800,400);//
    //    /////////////////////////////////////////////////
    //    if (gCherenkov_Correction==2) {
    //        cherenkov_correction_path ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/cherenkov_correction.root";//outFile_separation_PDF_pmt.root //cherenkov_correction.root
    //        cout<<"cherenkov_correction_path= " <<cherenkov_correction_path<<endl;
    //        ffile_cherenkov_correction  = new TFile(cherenkov_correction_path, "READ");
    //        int pmtCounter =0;
    //        for(Int_t PMT=0; PMT<PMT_num; PMT++) {
    //            TString pmt_counter=Form("_%d",pmtCounter);
    //            if(PMT<=10 || (PMT>=90 && PMT<=96)) continue; // dummy pmts
    //            if((PMT ==13 || PMT==14 || PMT==31 || PMT==32 || PMT==33 || PMT==12  || PMT==15 || PMT==34  || PMT==35 || PMT==35 || PMT==16 || PMT==17 || PMT==11 || PMT==102)) continue;
    //            //if(! (PMT==102)) continue; // custmization
    //            ////////////////////////////////////////////////////////////////////////
    //            //            glx_canvasAdd("r_pmt_correction"+pmt_counter,800,400);//
    //            ////////////////////////////////////////////////////////////////////////
    //            fHistPMT_read_k[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistPMT_k_%d",PMT));
    //            fHistPMT_read_pi[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistPMT_pi_%d",PMT));
    //            fit_PMT->SetParameters(100,0.82,0.010,10);
    //            fit_PMT->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    //            fit_PMT->SetParLimits(0,0.1,1E6);
    //            //fit_PMT->SetParLimits(1,0.82-2*cut_cangle,0.82+2*cut_cangle);
    //            fit_PMT->SetParLimits(1,0.809,0.835);
    //            fit_PMT->SetParLimits(2,0.005,0.030); // width
    //
    //            double  rang_min= 0.82-cut_cangle;
    //            double  rang_max= 0.82+cut_cangle;
    //
    //            //////////////////////////////////
    //            /// custumize fitting function ///
    //            //////////////////////////////////
    //
    //            if(PMT==42) rang_min=0.81;
    //            if(PMT==24 || PMT==26 ||PMT==18 ||PMT==11){
    //                rang_min= 0.82-cut_cangle/2;
    //                rang_max= 0.82+cut_cangle/2;
    //            }
    //            if(PMT==107){
    //                rang_min= 0.82-cut_cangle/2;
    //                rang_max= 0.82+cut_cangle;
    //
    //            }
    //            if(PMT==40 || PMT==41){
    //                rang_min= 0.814;
    //                rang_max= 0.84;
    //            }
    //
    //            fHistPMT_read_pi[PMT]->Fit("fit_PMT","M","",rang_min,rang_max);
    //            double histo_cor_entries = fHistPMT_read_pi[PMT]->GetEntries();
    //            double mean_cherenkov_cor=  fit_PMT->GetParameter(1);
    //            double sigma_cherenkov_cor= fit_PMT->GetParameter(2);
    //            double delta_cherenkov_cor= referance_angle - mean_cherenkov_cor;
    //            double val_1 = (delta_cherenkov_cor)*1000.0 ;
    //            shifted_pi->SetPoint(pmtCounter, mean_cherenkov_cor, PMT);
    //
    //            cout<<"##########"<<"PMT= "<<PMT<< "    delta_cherenkov_cor= " << delta_cherenkov_cor<<endl;
    //            //hsigma_test->Fill(sigma_cherenkov_cor*1000);
    //            //hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
    //            //hmean_test->Fill(mean_cherenkov_cor);
    //            // correction condition
    //            // if( ( fabs(delta_cherenkov_cor) <0.016 && (sigma_cherenkov_cor*1000<16 && sigma_cherenkov_cor*1000> 5.5) && histo_cor_entries>0 )) {//0.01  12  7
    //            array_correction[PMT]= delta_cherenkov_cor;
    //            cout<<"##########"<< "shift "<<val_1<<endl;
    //
    //
    //            /*
    //             // Fill PMT
    //             for(Int_t m=0; m<64; m++)
    //             for(Int_t n=0; n<64; n++){
    //             glx_hdigi[PMT]->SetBinContent(m, n,val_1);
    //             }
    //
    //             // default correction
    //             //array_correction[PMT]= -0.004;
    //             fHistPMT_read_pi[PMT]->Draw();
    //             //glx_canvasGet("r_pmt_correction"+pmt_counter)->Update();
    //             glx_canvasGet("r_pmt_correction")->Update();
    //             TLine *lin_ref = new TLine(0,0,0,1000);
    //             lin_ref->SetX1(referance_angle);
    //             lin_ref->SetX2(referance_angle);
    //             lin_ref->SetY1(gPad->GetUymin());
    //             lin_ref->SetY2(gPad->GetUymax());
    //             lin_ref->SetLineColor(kBlue);
    //             lin_ref->Draw();
    //
    //             TLine *lin_ref_k = new TLine(0,0,0,1000);
    //             lin_ref_k->SetX1(referance_angle_k);
    //             lin_ref_k->SetX2(referance_angle_k);
    //             lin_ref_k->SetY1(gPad->GetUymin());
    //             lin_ref_k->SetY2(gPad->GetUymax());
    //             lin_ref_k->SetLineColor(kRed);
    //             lin_ref_k->Draw();
    //
    //             //glx_canvasGet("r_pmt_correction"+pmt_counter)->Update();
    //             glx_canvasGet("r_pmt_correction")->Update();
    //             glx_waitPrimitive("r_pmt_correction");
    //             //hsigma_test->Fill(sigma_cherenkov_cor*1000);
    //             //hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
    //             //hmean_test->Fill(mean_cherenkov_cor);
    //             */
    //            // }
    //            ++pmtCounter;
    //        }
    //        /*
    //         glx_canvasAdd("r_pmt_shift",800,400);
    //         TMultiGraph *mg = new TMultiGraph();
    //         mg->Add(shifted_pi);
    //         mg->SetTitle(" Shift ; Mean [rad]; PMT ID [#]");
    //         mg->Draw("AP");
    //         //mg->GetHistogram()->GetYaxis()->SetRangeUser(6800,7050);
    //
    //         glx_canvasGet("r_pmt_shift")->Update();
    //
    //         TLine *lin_ref_pi = new TLine(0,0,0,1000);
    //         lin_ref_pi->SetX1(referance_angle_pi);
    //         lin_ref_pi->SetX2(referance_angle_pi);
    //         lin_ref_pi->SetY1(gPad->GetUymin());
    //         lin_ref_pi->SetY2(gPad->GetUymax());
    //         lin_ref_pi->SetLineColor(kBlue);
    //         lin_ref_pi->Draw();
    //
    //         glx_canvasGet("r_pmt_shift")->Update();
    //         TLine *lin_ref_k = new TLine(0,0,0,1000);
    //         lin_ref_k->SetX1(referance_angle_k);
    //         lin_ref_k->SetX2(referance_angle_k);
    //         lin_ref_k->SetY1(gPad->GetUymin());
    //         lin_ref_k->SetY2(gPad->GetUymax());
    //         lin_ref_k->SetLineColor(kRed);
    //         lin_ref_k->Draw();
    //         */
    //    }
    
    
    
    /*
     glx_canvasAdd("r_diff_test",800,400);
     hdiff_test->Draw();
     glx_canvasAdd("r_hsigma_test",800,400);
     hsigma_test->Draw();
     glx_canvasAdd("r_hmean_test",800,400);
     hmean_test->Draw();
     //gStyle->SetPalette(kLightTemperature);
     gStyle->SetPalette(kThermometer);
     double max_digi(10);//30
     double min_digi(-10);//-30
     glx_drawDigi("m,p,v\n",0, max_digi,min_digi);
     glx_canvasSave(0);
     
     return;
     */
    //return;
    
    //////////////////
    /// Reco Method //
    //////////////////
    TString outFile;
    if(gPDF_pix==1){outFile= "created_cherenkovPDF_pix.root";
    } else if(gPDF_pmt==1 && gCherenkov_Correction==2){ outFile= "created_pdf_pmt_corrected.root";
    } else if(gPDF_pix==2){ outFile= "outFile_separation_PDF_pix.root";
    } else if(gPDF_pmt==2){ outFile= "outFile_separation_PDF_pmt.root";
    } else if(gPDF_pmt==1){ outFile= "created_cherenkovPDF_pmt_NotCorrected.root";
    } else if(gCherenkov_Correction==1){ outFile= "cherenkov_correction.root";
    } else {
        //outFile= "outFile.root";
        //outFile= "out_"+infile;
        outFile= "out2_"+justName;
        
    }
    cout<<"####"<<outFile<<endl;
    double cop[]= {62,62,62,62,62,62,62,62,60,58,56,54,53,52,52,51,50,49,47,46,46,45,45,44,44,43,43,42,41,41,40,40,38,37.5,37,37.5,37,36.5,36,36};

    
    double shiftD_tdiff[24][40]={ 0,0,0,0,0,0,-0.212497,-0.168497,-0.131315,-0.120108,-0.115487,-0.0974607,-0.121494,-0.0912085,-0.0507834,0.0781914,0.245025,0.544255,0.750692,1.18085,1.71303,1.79146,1.24404,0.877314,0.602609,0.276667,0.0681443,-0.128882,-0.379789,-0.548475,-0.634345,-0.757345,-0.874287,-0.947678,-1.0129,0,0,0,0,0,0,0,0,0,0,0,-0.124542,-0.0886363,-0.0466412,-0.045753,0.0133786,0.0404589,0.0828445,0.0976572,0.103918,0.179398,0.325307,0.60341,0.8109,0.933208,1.18353,1.22501,1.1194,0.920326,0.689237,0.382306,0.237853,0.111665,-0.11591,-0.319721,-0.397971,-0.517027,-0.647768,-0.609641,0,0,0,0,0,0,0,0,0,0,0,0,-0.050176,-0.0435546,-0.0271583,0.0379114,-0.00836784,0.0727276,0.138301,0.124994,0.122776,0.198836,0.318348,0.517275,0.773196,0.809877,0.89285,0.921646,0.875704,0.807228,0.591924,0.372603,0.216679,0.137514,-0.0471617,-0.165383,-0.260127,-0.340788,-0.394445,-0.58131,0,0,0,0,0,0,0,0,0,0,0,0,-0.0353932,-0.0318091,-0.046814,-0.0643999,-0.0876217,-0.022806,0.0741525,0.104005,0.09392,0.142027,0.212933,0.388194,0.592119,0.664085,0.69806,0.702764,0.635391,0.60563,0.42235,0.276258,0.16021,0.0280315,-0.0709147,-0.124245,-0.221172,-0.334673,-0.344158,-0.431136,0,0,0,0,0,0,0,0,0,0,0,0,-0.135951,-0.117562,-0.124928,-0.131315,-0.122391,-0.0783247,-0.00842946,-0.0175169,-0.00531461,-0.0110562,0.0627202,0.236342,0.395865,0.438967,0.478636,0.481237,0.410696,0.38114,0.298232,0.168116,0.0185137,-0.0376193,-0.0811981,-0.1695,-0.262101,-0.405857,-0.452009,-0.39605,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.106184,-0.0930739,-0.0948498,-0.0973551,-0.0751308,-0.0157032,0.0064228,-0.0167574,-0.0401675,0.0488597,0.186276,0.369384,0.396924,0.405641,0.40356,0.348608,0.345339,0.244302,0.137854,0.00994912,-0.029286,-0.0672536,-0.0787963,-0.175779,-0.25277,-0.291674,-0.237686,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.121473,-0.0457286,-0.0574831,-0.059427,-0.0297513,-0.00759922,-0.0174272,-0.0552326,-0.0457817,0.0584857,0.178635,0.318815,0.363596,0.353752,0.336481,0.305996,0.277141,0.188933,0.074431,-0.0205597,-0.0264286,-0.0882313,-0.108045,-0.116541,-0.146015,-0.175727,-0.202849,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0297215,-0.0298446,-0.0192065,-0.00883032,-0.00343533,0.0140206,-0.0445895,-0.089882,-0.0625097,0.0466929,0.190052,0.288683,0.314223,0.28853,0.275963,0.233441,0.20829,0.147475,0.0271816,-0.0247789,-0.00918168,-0.0736942,-0.100536,-0.10446,-0.0480655,-0.0339893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0376528,-0.0085324,-0.0201655,-0.0307386,-0.0181594,-0.0962811,-0.111502,-0.0871813,-0.030972,0.129332,0.224974,0.233403,0.178783,0.174591,0.145225,0.149318,0.0805856,-0.0124332,-0.030286,-0.0552115,-0.102939,-0.121996,-0.0974656,-0.0849368,-0.029955,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0212187,-0.0426328,-0.0196482,-0.00426026,0.0525694,-0.0533956,-0.0645615,-0.0303608,0.0127104,0.153931,0.242974,0.272539,0.233415,0.185906,0.161315,0.158384,0.114722,0.0760059,-0.00179772,0.00177117,-0.055545,-0.103146,-0.0846749,-0.0532756,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0265657,-0.0244398,-0.00667925,-0.0153027,-0.0494108,-0.0599105,-0.0386066,-0.00836632,0.169426,0.23827,0.238785,0.228106,0.205007,0.177146,0.156961,0.140241,0.0602354,0.0611113,0.0471731,0.0238059,-0.0386623,-0.0471459,-0.0861525,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.102169,-0.0759268,-0.106632,-0.128279,-0.100222,-0.0663984,-0.0229217,0.123481,0.204289,0.207343,0.180538,0.172402,0.171432,0.164627,0.121178,0.0371667,0.0049031,0.00457614,0.0515949,-0.0478936,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.220337,-0.175774,-0.148908,-0.0299439,0.0536656,0.0411358,0.026935,-0.0075003,0.00220719,0.0135967,-0.00426324,-0.00879467,-0.0518067,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.160375,-0.104351,-0.0398385,0.112564,0.182259,0.18584,0.0612845,0.0685847,0.0265474,0.07112,0.120471,0.046869,0.0444415,0.0526628,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0115343,0.151956,0.160157,0.170822,0.0302506,0.0151447,0.0848597,0.0263626,0.0996271,0.0571817,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0275436,0.0649307,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ;
    
    double shiftR_tdiff[24][40]={ 0,0,0,0,0,0,-0.581964,-0.664308,-0.652945,-0.534984,-0.510857,-0.302394,-0.416157,-0.28187,0.0186578,0.126143,0.255395,0.555315,0.774699,1.25938,1.80409,1.77441,1.28966,1.02125,0.802201,0.563475,0.319484,0.152015,-0.00706454,-0.150087,-0.223669,-0.25439,-0.318416,-0.45389,-0.482179,0,0,0,0,0,0,0,0,0,0,0,-1.09325,-0.749646,-0.50478,-0.355006,-0.383651,-0.216663,-0.185276,-0.0278354,0.151185,0.304574,0.462995,0.71319,0.921894,1.10713,1.30609,1.35047,1.24878,1.08459,0.844757,0.604089,0.428763,0.328786,0.221729,0.0584101,-0.00166835,-0.125871,-0.221795,-0.270271,0,0,0,0,0,0,0,0,0,0,0,0,-0.351929,-0.657451,-0.280183,-0.162935,-0.157293,-0.174671,-0.0985791,-0.0793121,0.0585537,0.25064,0.386326,0.664939,0.847615,0.896619,1.05586,1.07643,1.0456,0.951191,0.717313,0.564453,0.401,0.372185,0.287688,0.176245,0.0776931,-0.0234355,-0.107355,-0.16243,0,0,0,0,0,0,0,0,0,0,0,0,0.186884,0.346454,0.233309,0.0620987,-0.16476,-0.0556778,0.02804,-0.00272601,0.0978376,0.250333,0.370768,0.622158,0.758156,0.810971,0.906285,0.918509,0.865927,0.826387,0.646789,0.537723,0.418315,0.365248,0.276698,0.245276,0.153009,0.0462721,0.00368508,-0.0769824,0,0,0,0,0,0,0,0,0,0,0,0,0.0243746,0.10191,0.000871238,-0.186645,-0.236495,-0.168105,-0.0836931,-0.0392987,0.0180929,0.110715,0.194889,0.420275,0.588097,0.610185,0.673065,0.700161,0.669483,0.601534,0.431978,0.339815,0.271006,0.285635,0.230478,0.153938,0.0885552,-0.0284336,0.00823109,-0.0853845,0,0,0,0,0,0,0,0,0,0,0,0,0,0.152089,0.0701478,-0.0485647,-0.114893,-0.0384545,0.0565695,0.054406,0.0578533,0.0819917,0.17356,0.388854,0.571162,0.592309,0.637545,0.654218,0.639932,0.618113,0.447969,0.366449,0.325117,0.290085,0.260138,0.210951,0.155837,0.112023,0.100376,-0.00893973,0,0,0,0,0,0,0,0,0,0,0,0,0,0.350998,0.471647,0.406919,0.246531,0.149684,0.161659,0.116077,0.116356,0.128173,0.206709,0.369228,0.53044,0.590443,0.623112,0.62634,0.629609,0.608003,0.505345,0.413999,0.34944,0.320103,0.268243,0.226914,0.188386,0.146957,0.156163,0.107816,0,0,0,0,0,0,0,0,0,0,0,0,0,0.325289,0.142475,0.230344,0.247083,0.0686047,0.0580396,0.0386985,0.0129045,0.0392193,0.0942243,0.280547,0.37676,0.415877,0.439063,0.441126,0.435359,0.408119,0.369388,0.364975,0.30189,0.291175,0.237647,0.191658,0.161091,0.17956,0.14481,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.290689,0.212513,0.170077,0.0875794,-0.00199978,-0.104823,-0.0745239,0.0275235,0.0420727,0.215527,0.31249,0.33236,0.321457,0.348461,0.339956,0.348369,0.301067,0.313453,0.27424,0.26608,0.221358,0.180576,0.165817,0.128184,0.115419,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.418682,0.34355,0.286488,0.170779,0.101755,0.0633065,0.0700362,0.150909,0.196778,0.316952,0.394009,0.436356,0.442196,0.444018,0.439615,0.422846,0.402918,0.404907,0.361663,0.353314,0.302624,0.272173,0.238751,0.200676,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.422634,0.384061,0.33505,0.253009,0.259315,0.197915,0.247598,0.296218,0.447647,0.531396,0.563429,0.589411,0.581974,0.584494,0.577116,0.540153,0.483452,0.456252,0.511952,0.487791,0.396934,0.337398,0.26956,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0481874,0.0766456,0.0909551,0.0958395,0.0320431,0.0533661,0.135747,0.25381,0.357743,0.359198,0.346693,0.337905,0.333841,0.305407,0.299968,0.266037,0.231725,0.258377,0.244533,0.170873,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.013011,0.0363971,0.0675664,0.140749,0.265944,0.237025,0.204179,0.169915,0.148941,0.149919,0.165732,0.14353,0.141517,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.148677,0.197085,0.250668,0.327356,0.40475,0.38439,0.315532,0.301604,0.26037,0.268657,0.306959,0.323411,0.293989,0.293496,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.208064,0.303278,0.378298,0.324918,0.23574,0.186843,0.19051,0.204916,0.280738,0.342376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.293563,0.243974,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ;
    
    /////////////////////////////////
    //////// Creat file and trees ///
    /////////////////////////////////
    
    TFile file(outFile,"recreate");
    //TTree *tree_cut = glx_ch->CloneTree(0);
    
    
    //bool btree= false; // for histo not tree
    bool btree= true; //for tree not histo
    
    TTree tree_variables("tree_variables","tree for cherenkov track resolution");
    double track_spr(-1),track_mean(-1), track_yield(-1), track_mom(-1), track_xbar(0),track_ybar(0),track_fit_chisqu(-1),track_fit_NDF(-1);
    int track_pid(-1), track_nbar(-1), track_x_pos_bin(-1);
    TString track_file="noname";
    std::vector<int> vpx;
    std::vector<int> vpy;
    std::vector<int> vpz;
    std::vector<double> vtdiff;
    std::vector<double> vtime;
    std::vector<double> vtangle;
    std::vector<bool> vreflected;
    
    double track_inv_mass(-1),track_missing_mass(-1),track_chi_square(-1),track_TofTrackDist(-1);
    
    if(btree){
        tree_variables.Branch("track_pid",&track_pid,"track_pid/I");
        tree_variables.Branch("track_spr",&track_spr,"track_spr/D");
        tree_variables.Branch("track_mean",&track_mean,"track_mean/D");
        tree_variables.Branch("track_yield",&track_yield,"track_yield/D");
        tree_variables.Branch("track_mom",&track_mom,"track_mom/D");
        tree_variables.Branch("track_xbar",&track_xbar,"track_xbar/D");
        tree_variables.Branch("track_ybar",&track_ybar,"track_ybar/D");
        tree_variables.Branch("track_nbar",&track_nbar,"track_nbar/I");
        tree_variables.Branch("track_fit_chisqu",&track_fit_chisqu,"track_fit_chisqu/D");
        tree_variables.Branch("track_fit_NDF",&track_fit_NDF,"track_fit_NDF/D");
        tree_variables.Branch("track_file",&track_file,"track_file/C");
        tree_variables.Branch("vpx",&vpx);
        tree_variables.Branch("vpy",&vpy);
        tree_variables.Branch("vpz",&vpz);
        tree_variables.Branch("vtdiff",&vtdiff);
        tree_variables.Branch("vtime",&vtime);
        tree_variables.Branch("vtangle",&vtangle);
        tree_variables.Branch("vreflected",&vreflected);
        
        tree_variables.Branch("track_inv_mass",&track_inv_mass,"track_inv_mass/D");
        tree_variables.Branch("track_missing_mass",&track_missing_mass,"track_missing_mass/D");
        tree_variables.Branch("track_chi_square",&track_chi_square,"track_chi_square/D");
        tree_variables.Branch("track_TofTrackDist",&track_TofTrackDist,"track_TofTrackDist/D");
        tree_variables.Branch("track_x_pos_bin",&track_x_pos_bin,"track_x_pos_bin/I");
    }
    
    
    
    double pion_counter =0;
    double kaon_counter =0;
    DrcHit hit;
    for (int e = 0; e < glx_ch->GetEntries(); e++){
        glx_ch->GetEntry(e);
        //if(e>2000)break;
        for (int t = 0; t < glx_events->GetEntriesFast(); t++){
            // cut in Event criteria here
            
            glx_nextEventc(e,t,1000);
            posInBar = glx_event->GetPosition();
            int x_pos_bin = histo_tmp_pos->GetXaxis()->FindBin(posInBar.X());
            momInBar = glx_event->GetMomentum();
            double momentum = momInBar.Mag();
            int pdgId = glx_findPdgId(glx_event->GetPdg());
            int bar = glx_event->GetId();
            //if(count[pdgId]>1000) continue;
            
            //std::cout<<"######### No Problem "<<pdgId<<std::endl;
            
            double inv_mass=  glx_event->GetInvMass();
            double missing_mass=  glx_event->GetMissMass();
            double chi_square=  glx_event->GetChiSq();
            double TofTrackDist=  glx_event->GetTofTrackDist();
            
            ////////////////////////
            //////// selection//////
            ////////////////////////
            
            if (!(pdgId ==2 || pdgId ==3)) continue;
            
            // pion =2  ,kaon=3
            if(!btree){ 
                if (pdgId == 2 && chi_square> 10) continue;
                if (pdgId == 3 && chi_square> 20)continue;
                if (pdgId == 2 && (inv_mass< 0.66  || inv_mass> 0.82 )  ) continue;
                if (pdgId == 3 && (inv_mass< 1.015 || inv_mass> 1.025) ) continue;
                if(missing_mass< -0.01 || missing_mass> 0.01 )continue;
            }
            
            //            momInBar_unit=momInBar.Unit();
            //            double dir_x =momInBar_unit.X();
            //            double dir_y =momInBar_unit.Y();
            //            double dir_z =momInBar_unit.Z();
            //hdir_x->Fill(dir_x);
            //hdir_y->Fill(dir_y);
            //hdir_z->Fill(dir_z);
            //cout<<"=========>" << dir_x << "  "<< dir_y<< endl;
            
            
            if(glx_event->GetType()!=2) continue; //1-LED 2-beam 0-rest
            if(glx_event->GetParent()>0) continue;
            
            //////////////////////////////
            //////// Momentum Cut //////
            /////////////////////////////
            
            //if(momInBar.Mag()<3.5 || momInBar.Mag()>4.0 ) continue;
            
            /////////////////////////////
            //////// DIRC Wall Cut //////
            /////////////////////////////
            
            //if(momInBar.Mag()<2.8 || momInBar.Mag()>3 ) continue;
            //if(momInBar.Mag()<3.9 || momInBar.Mag()>4.1 ) continue;
            //            int bin = (100+posInBar.X())/200.*nbins;
            // not used
            // if(bar<0 || bar>=luts || (bar!=ybar && ybar!=-1)) continue;
            // if(bin<0 || bin>nbins || (bin!=xbar && xbar!=-1)) continue;
            // commented
            //if(bar<0 || bar>=luts || (bar<4 || bar>8)) continue;
            //if(bin<0 || bin>nbins || (bin<7 || bin>13)) continue;
            //std::cout<<"##################### bar "<<bar<<" "<<"####### bin "<<bin<<std::endl;
            //if ( posInBar.X()>10 || posInBar.X() < -10 ) continue;
            //
            //            if (xmin <  posInBar.X()) xmin=posInBar.X();
            //            if (xmax >  posInBar.X()) xmax=posInBar.X();
            //            if (ymin <  posInBar.Y()) ymin=posInBar.Y();
            //            if (ymax >  posInBar.Y()) ymax=posInBar.Y();
            
            
            
            
            
            
            ///////////////////////
            //////// pos map //////
            ///////////////////////
            
            //            if (pdgId == 2){
            //                if (glx_event->GetPdg() > 0 ) hExtrapolatedBarHitXY_k->Fill(100+posInBar.X(), posInBar.Y());
            //                if (glx_event->GetPdg() < 0 ) hExtrapolatedBarHitXY_k->Fill(-100+ posInBar.X(), posInBar.Y());
            //            }
            //            if (pdgId == 3){
            //                if (glx_event->GetPdg() > 0 ) hExtrapolatedBarHitXY_pi->Fill(100+posInBar.X(), posInBar.Y());
            //                if (glx_event->GetPdg() < 0 ) hExtrapolatedBarHitXY_pi->Fill(-100+ posInBar.X(), posInBar.Y());
            //            }
            
            /////////////////////////////
            //////// theta phi map //////
            /////////////////////////////
            
            //            double theta_mom =  momInBar_unit.Theta()* 180/PI;
            //            double ph_mom =  momInBar_unit.Phi()* 180/PI;
            //            if (glx_event->GetPdg() < 0 ) theta_mom = theta_mom *-1.0 ;
            //
            //            int theta_bin(-1), phi_bin(-1);
            //            double content_histo_theta_phi_map(-1), content_histo_theta_phi_mom_map(-1), average_bin(-1);
            //
            //            if (pdgId == 2){
            //                histo_theta_phi_map_pi->Fill(theta_mom,ph_mom);
            //                histo_theta_phi_mom_tmp_map_pi->Fill(theta_mom,ph_mom,momInBar.Mag());
            //
            //                theta_bin = histo_theta_phi_map_pi->GetXaxis()->FindBin(theta_mom);
            //                phi_bin = histo_theta_phi_map_pi->GetYaxis()->FindBin(ph_mom);
            //                content_histo_theta_phi_map=histo_theta_phi_map_pi->GetBinContent(theta_bin,phi_bin);
            //                content_histo_theta_phi_mom_map=histo_theta_phi_mom_tmp_map_pi->GetBinContent(theta_bin,phi_bin);
            //                average_bin= content_histo_theta_phi_mom_map/content_histo_theta_phi_map;
            //                //cout<< "###### average_bin= "<<average_bin<<" "<<content_histo_theta_phi_map<<" "<<content_histo_theta_phi_mom_map<<" "<<momInBar_unit.Mag()<<endl;
            //                histo_theta_phi_mom_map_pi->SetBinContent(theta_bin,phi_bin,average_bin);
            //            }
            //            if (pdgId == 3){
            //                histo_theta_phi_map_k->Fill(theta_mom,ph_mom);
            //                histo_theta_phi_mom_tmp_map_k->Fill(theta_mom,ph_mom,momInBar.Mag());
            //
            //                theta_bin = histo_theta_phi_map_k->GetXaxis()->FindBin(theta_mom);
            //                phi_bin = histo_theta_phi_map_k->GetYaxis()->FindBin(ph_mom);
            //                content_histo_theta_phi_map=histo_theta_phi_map_k->GetBinContent(theta_bin,phi_bin);
            //                content_histo_theta_phi_mom_map=histo_theta_phi_mom_tmp_map_k->GetBinContent(theta_bin,phi_bin);
            //                average_bin= content_histo_theta_phi_mom_map/content_histo_theta_phi_map;
            //
            //                histo_theta_phi_mom_map_k->SetBinContent(theta_bin,phi_bin,average_bin);
            //            }
            
            
            /////////////////////////////
            //////// DIRC Wall Cut //////
            /////////////////////////////
            /*
             
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
             
             
             if (xmin >  posInBar.X()) xmin=posInBar.X();
             if (xmax <  posInBar.X()) xmax=posInBar.X();
             
             if (ymin >  posInBar.Y()) ymin=posInBar.Y();
             if (ymax <  posInBar.Y()) ymax=posInBar.Y();
             
             */
            
            //tree_cut->Fill(); old
            //glx_event->Clear();
            //cout<< "ID= "<<pdgId<<" theta_mom= "<<theta_mom<<"ph_mom,momentum= "<<ph_mom<<endl;
            
            
            
            //            if (glx_event->GetPdg() > 0 ) hExtrapolatedBarHitXY_cut->Fill(100+posInBar.X(), posInBar.Y());
            //            if (glx_event->GetPdg() < 0 ) hExtrapolatedBarHitXY_cut->Fill(-100+  posInBar.X(), posInBar.Y());
            
            //if(fabs(dir_x)>0.01 )continue;
            //if(dir_y<-0.05)continue;
            //            hdir_x->Fill(dir_x);
            //            hdir_y->Fill(dir_y);
            //            hdir_z->Fill(dir_z);
            //            if (pdgId == 2){
            //                hist_ev_rho_mass_cut->Fill(inv_mass);
            //                hist_ev_missing_mass_rho_cut->Fill(missing_mass);
            //                hist_ev_chi_rho_cut->Fill(chi_square);
            //                if (glx_event->GetPdg() > 0 ) mom_theta_rho_cut->Fill(momInBar.Theta()*180/PI, momInBar.Mag());
            //                if (glx_event->GetPdg() < 0 ) mom_theta_rho_cut->Fill(-1.0 * momInBar.Theta()*180/PI, momInBar.Mag());
            //            }
            //            if (pdgId == 3){
            //                hist_ev_phi_mass_cut->Fill(inv_mass);
            //                hist_ev_missing_mass_phi_cut->Fill(missing_mass);
            //                hist_ev_chi_phi_cut->Fill(chi_square);
            //                if (glx_event->GetPdg() > 0 ) mom_theta_phi_cut->Fill(momInBar.Theta()*180/PI, momInBar.Mag());
            //                if (glx_event->GetPdg() < 0 ) mom_theta_phi_cut->Fill(-1.0 * momInBar.Theta()*180/PI, momInBar.Mag());
            //            }
            //
            
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
            
            // reset variables
            if(btree){
                track_spr=-1; track_mean=-1; track_yield=-1; track_mom=-1; track_xbar=0;track_ybar=0;
                track_pid=-1; track_nbar=-1;
                histo_cherenkov_track->Reset();
                vpx.clear();
                vpy.clear();
                vpz.clear();
                vtdiff.clear();
                vtime.clear();
                vtangle.clear();
                vreflected.clear();
            }
            
            
            for(int h = 0; h < glx_event->GetHitSize(); h++){
                hit = glx_event->GetHit(h);
                int ch = hit.GetChannel();
                int pmt = hit.GetPmtId();
                int pix = hit.GetPixelId();
                double hitTime = hit.GetLeadTime()-glx_event->GetTime();
                
                if(ch>glx_nch) continue;
                //histo_time_bar_pos[bar][x_pos_bin]->Fill(hitTime);
                //if(hitTime>40) continue;
                nphc++;
                
                /////////////////////////////////////
                // Reflection condition may change //
                /////////////////////////////////////
                //bool reflected = hitTime>40;
                bool reflected = true;
                if(hitTime < cop[x_pos_bin]){
                    reflected = false;
                    //if (x_pos_bin==20 && (bar==1 || bar==0)) cout<<"bar= "<<bar<<" seg= "<<x_pos_bin <<" shiftD="<<shiftD_tdiff[bar][x_pos_bin]<<" hitTime= "<<hitTime<<endl;
                    hitTime=hitTime + shiftD_tdiff[bar][x_pos_bin];
                    //cout<<" hitTime= "<<hitTime<<endl;
                    
                }else hitTime=hitTime + shiftR_tdiff[bar][x_pos_bin];
                
                
                
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
                            
                            //if(gCherenkov_Correction == 2) tangle = momInBar.Angle(dir)+ referance_angle - 0.8257;
                            //if(gCherenkov_Correction != 2) tangle = momInBar.Angle(dir);
                            //if(gCherenkov_Correction == 2) tangle = momInBar.Angle(dir) + array_correction[pmt];
                            tangle = momInBar.Angle(dir);
                            
                            //double bartime = lenz/cos(luttheta)/20.4; //198 //203.767 for 1.47125
                            double bartime = lenz/cos(luttheta)/19.6; //203.767 for 1.47125
                            double totalTime = bartime+evtime;
                            // hTime->Fill(hitTime);
                            // hCalc->Fill(totalTime);
                            
                            if(!btree && (fabs(tangle-0.5*(mAngle[2]+mAngle[3]))<0.01 )){ // 10 mrad
                                //hDiff->Fill(totalTime-hitTime);
                                //if(samepath)
                                {
                                    //hDiffT->Fill(totalTime-hitTime);
                                    if(r) {
                                        //hDiffR->Fill(totalTime-hitTime);
                                        histo_tdiffR_bar_pos[bar][x_pos_bin]->Fill(totalTime-hitTime);
                                        
                                    }
                                    else {
                                        //hDiffD->Fill(totalTime-hitTime);
                                        histo_tdiffD_bar_pos[bar][x_pos_bin]->Fill(totalTime-hitTime);
                                        
                                    }
                                    
                                    
                                }
                                
                                if(fabs(totalTime-hitTime)< 5)histo_time_bar_pos[bar][x_pos_bin]->Fill(hitTime);
                            }
                            // skim
                            if(fabs(totalTime-hitTime)> 10) continue;
                            if(tangle > 1.0) continue;
                            
                            if (btree){
                                vtdiff.push_back(totalTime-hitTime);
                                vtime.push_back(hitTime);
                                vtangle.push_back(tangle);
                                vreflected.push_back(reflected);
                            }
                            
                            ///////////////
                            // Time Cut  //
                            ///////////////
                            if(!r && fabs(totalTime-hitTime)>cut_tdiffd) continue; // removed
                            if(r && fabs(totalTime-hitTime) >cut_tdiffr) continue; // removed
                            
                            //////////////////////
                            // Cherenkov track  //
                            //////////////////////
                            
                            if (!btree) histo_cherenkov_track->Fill(tangle);
                            
                            //                            ////////////////
                            //                            // Fill PDF   //
                            //                            ////////////////
                            //                            // cherenkove PDF per PIX
                            //                            if(gPDF_pix ==1 && pdgId == 3) fHistCh_k[ch]->Fill(tangle); // good after time cut
                            //                            if(gPDF_pix ==1 && pdgId == 2) fHistCh_pi[ch]->Fill(tangle); // good after time cut
                            //
                            //                            // cherenkove PDF per PMT
                            //                            if(gPDF_pmt ==1 && pdgId == 3) fHistPMT_PDF_k[pmt]->Fill(tangle); // good after time cut
                            //                            if(gPDF_pmt ==1 && pdgId == 2) fHistPMT_PDF_pi[pmt]->Fill(tangle); // good after time cut
                            //
                            //                            ////////////////////////////
                            //                            // Fill PMT coorrection   //
                            //                            ////////////////////////////
                            //                            // cherenkove correction per PMT
                            //                            if(gCherenkov_Correction ==1 && pdgId == 3) fHistPMT_k[pmt]->Fill(tangle);
                            //                            if(gCherenkov_Correction ==1 && pdgId == 2) fHistPMT_pi[pmt]->Fill(tangle);
                            //
                            //                            // fill cherenkove histo
                            //                            hAngle[pdgId]->Fill(tangle);
                            //
                            ////////////////////
                            // Cherenkov Cut  //
                            ////////////////////
                            if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))>cut_cangle) continue;
                            //if(fabs(tangle-0.5*(referance_angle_k+referance_angle_pi))>cut_cangle) continue; // removed
                            //if(tangle> 0.844 ||tangle < 0.798)  continue;
                            
                            isGood=true;
                            //                            hTime->Fill(hitTime);
                            //                            hCalc->Fill(totalTime);
                            //
                            //
                            //                            if(!(gPDF_pix==2||gPDF_pmt==2 )){ // fixed
                            //                                sum1 += TMath::Log(fAngle[2]->Eval(tangle)+noise);
                            //                                sum2 += TMath::Log(fAngle[3]->Eval(tangle)+noise);
                            //                            }
                            //
                            //                            if(gPDF_pix ==2){
                            //                                // use histograms
                            //                                Int_t kk = fHistCh_read_k[ch]->GetXaxis()->FindBin(tangle);
                            //                                Int_t kpi = fHistCh_read_pi[ch]->GetXaxis()->FindBin(tangle);
                            //                                if (fHistCh_read_pi[ch]->GetBinContent(kpi) > 0 )sum1 += TMath::Log(fHistCh_read_pi[ch]->GetBinContent(kpi));
                            //                                if (fHistCh_read_k[ch]->GetBinContent(kk) > 0 )sum2 += TMath::Log(fHistCh_read_k[ch]->GetBinContent(kk));
                            //
                            //                                //if (sum1 != 0 || sum2!=0 )std::cout<<"No Problem  separation  " <<kpi<<" "<<kk<<"  sum "<<sum1 <<"  "<< sum2<<std::endl;
                            //                                //std::cout<<"###### No Problem  separation  " << fHistCh_read_k[ch]->GetBinContent(kk) <<"  "<< fHistCh_read_pi[ch]->GetBinContent(kpi)<<std::endl;
                            //                            }
                            //
                            //                            if(gPDF_pmt ==2){
                            //                                // use histograms
                            //                                Int_t k_bin = fHistPMT_PDF_read_k[pmt]->GetXaxis()->FindBin(tangle);
                            //                                Int_t pi_bin = fHistPMT_PDF_read_pi[pmt]->GetXaxis()->FindBin(tangle);
                            //                                if (fHistPMT_PDF_read_pi[pmt]->GetBinContent(pi_bin) > 0 )sum1 += TMath::Log(fHistPMT_PDF_read_pi[pmt]->GetBinContent(pi_bin));
                            //                                if (fHistPMT_PDF_read_k[pmt]->GetBinContent(k_bin) > 0 )sum2 += TMath::Log(fHistPMT_PDF_read_k[pmt]->GetBinContent(k_bin));
                            //
                            //
                            //                            }
                            //
                            //
                            //                            if(0){
                            //                                TString x=(sum1>sum2)? " <====== PION" : "";
                            //                                std::cout<<Form("%1.6f  %1.6f | %1.6f  %1.6f        pid %d",TMath::Log(fAngle[2]->Eval(tangle)+noise),TMath::Log(fAngle[3]->Eval(tangle)+noise), sum1, sum2,pdgId)<<"  " <<std::endl;
                            //
                            //                                cc->cd();
                            //                                fAngle[2]->Draw("");
                            //                                fAngle[3]->Draw("same");
                            //
                            //                                cc->Update();
                            //                                gLine->SetLineWidth(2);
                            //                                gLine->SetX1(tangle);
                            //                                gLine->SetX2(tangle);
                            //                                gLine->SetY1(cc->GetUymin());
                            //                                gLine->SetY2(cc->GetUymax());
                            //                                gLine->Draw();
                            //                                cc->Update();
                            //                                cc->WaitPrimitive();
                            //                            }
                            
                        } // bar ambiguities
                    } // reflection loop
                } // LUT loop
                
                if(isGood){
                    nph++;
                    //if (glx_event->GetPdg() > 0 ) nph_p++;
                    //if (glx_event->GetPdg() < 0 ) nph_n++;
                    if(pmt<108) {
                        //glx_hdigi[pmt]->Fill(pix%8, pix/8);
                        if (btree){
                            vpx.push_back(pmt);
                            vpy.push_back(pix%8);
                            vpz.push_back(pix/8);
                        }
                        
                        goodevt=1;
                    }
                }
            } // hit loop
            
            if(goodevt) evtcount++;
            if(nph<5) continue;
            //hNph[pdgId]->Fill(nph);
            
            //hNph_p[pdgId]->Fill(nph_p);
            //hNph_n[pdgId]->Fill(nph_n);
            
            //hNphC->Fill(nphc);
            
            double sum = sum1-sum2;
            //cout<<"########### sum  "<<sum <<"  sum1  "<<sum1<<"  sum2  "<<sum2<<endl;
            //hLnDiff[pdgId]->Fill(sum);
            
            count[pdgId]++;
            
            //            if(0 && pdgId==3){
            //                //    if(!cc)
            //                TString x=(sum1>sum2)? " <====== Pion" : "";
            //                // std::cout<<Form("f %1.6f s %1.6f PMT %d pix %d   pid %d",aminf,amins,PMT,pix  ,prt_particle)<<"  "<<x <<std::endl;
            //
            //                std::cout<<"PID "<< glx_event->GetPdg() <<" sum1 "<<sum1<<" sum2 "<<sum2<<" sum "<<sum<<" "<<x<<std::endl;
            //
            //                cc->cd();
            //
            //                if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
            //                if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
            //
            //                hAngle[2]->Draw("hist");
            //                hAngle[3]->Draw("hist same");
            //                fAngle[2]->Draw("same");
            //                fAngle[3]->Draw("same");
            //
            //                // hAngle[2]->GetYaxis()->SetRangeUser(0,20);
            //                // hAngle[3]->GetYaxis()->SetRangeUser(0,20);
            //
            //                cc->Update();
            //                TLine *line = new TLine(0,0,0,1000);
            //                line->SetX1(mAngle[3]);
            //                line->SetX2(mAngle[3]);
            //                line->SetY1(cc->GetUymin());
            //                line->SetY2(cc->GetUymax());
            //                line->SetLineColor(kRed);
            //                line->Draw();
            //
            //                TLine *line2 = new TLine(0,0,0,1000);
            //                line2->SetX1(mAngle[2]);
            //                line2->SetX2(mAngle[2]);
            //                line2->SetY1(cc->GetUymin());
            //                line2->SetY2(cc->GetUymax());
            //                line2->SetLineColor(kBlue);
            //                line2->Draw();
            //
            //                cc->Update();
            //                cc->WaitPrimitive();
            //            }
            
            // hAngle[2]->Reset();
            // hAngle[3]->Reset();
            
            /////////////////////
            // fill tree here////
            /////////////////////
            if (btree){
                fit_track->SetParameters(100,0.82,0.010,10);
                fit_track->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                fit_track->SetParLimits(0,0.1,1E6);
                fit_track->SetParLimits(1,0.809,0.835);
                fit_track->SetParLimits(2,0.005,0.030);
                if (pdgId==3)fit_track->SetLineColor(kRed);
                else fit_track->SetLineColor(kBlue);
                histo_cherenkov_track->Fit("fit_track","MQ0","", 0.5*(mAngle[2]+mAngle[3])-cut_cangle, 0.5*(mAngle[2]+mAngle[3])-cut_cangle) ;
                
                //cc->cd();
                //histo_cherenkov_track->Draw();
                //cc->Update();
                //cc->WaitPrimitive();
                
                track_mean=  fit_track->GetParameter(1);
                track_spr= fit_track->GetParameter(2);
                track_yield = nph;
                track_mom = momInBar.Mag();
                track_xbar = posInBar.X();
                track_ybar = posInBar.Y();
                track_pid = pdgId;
                track_nbar = bar;
                track_fit_chisqu = fit_track->GetChisquare();
                track_fit_NDF = fit_track->GetNDF();
                track_file= justName;
                
                track_inv_mass= inv_mass;
                track_missing_mass= missing_mass;
                track_chi_square= chi_square;
                track_TofTrackDist= TofTrackDist;
                
                track_x_pos_bin=  x_pos_bin;
                //cout<<"#### mom "<<track_mom<<endl;
                tree_variables.Fill();
            }
            
            
            ///////////////////////////////////
            //////// reduce pions number //////
            ///////////////////////////////////
            //double percentage = kaon_counter/pion_counter*100.0;
            //if (percentage <100 && pdgId == 2 )continue;
            //if (pdgId == 2) pion_counter++;
            //if (pdgId == 3) kaon_counter++;
            
        }
    } // Event Loop
    
    //    if(evtcount>0){
    //        for(int i=0; i<glx_nch; i++){
    //            int pmt=i/64;
    //            int pix=i%64;
    //            double rel = glx_hdigi[pmt]->GetBinContent(pix%8+1,pix/8+1)/(double)evtcount;
    //            glx_hdigi[pmt]->SetBinContent(pix%8+1, pix/8+1,rel);
    //        }
    //    }
    //
    //    //TString nid=Form("_%2.2f_%2.2f",theta,phi);
    //    TString nid=Form("_%d_%d",xbar,ybar);
    //
    //    glx_drawDigi("m,p,v\n",0);
    //    glx_cdigi->SetName("hp"+nid);
    //    glx_canvasAdd(glx_cdigi);
    //
    //    glx_canvasAdd("hAngle"+nid,800,400);
    //
    //    //scal
    //    if(hAngle[2]->GetMaximum()>0) hAngle[2]->Scale(1/hAngle[2]->GetMaximum());
    //    if(hAngle[3]->GetMaximum()>0) hAngle[3]->Scale(1/hAngle[3]->GetMaximum());
    //
    //    for(int i=0; i<5; i++){
    //        if(hAngle[i]->GetEntries()<20) continue;
    //
    //        int nfound = spect->Search(hAngle[i],1,"goff",0.9);
    //        if(nfound>0) cherenkovreco[i] = spect->GetPositionX()[0];
    //        else cherenkovreco[i] =  hAngle[i]->GetXaxis()->GetBinCenter(hAngle[i]->GetMaximumBin());
    //        if(cherenkovreco[i]>0.85) cherenkovreco[i]=0.82;
    //
    //        if(i==2)  fit->SetLineColor(kBlue);
    //        if(i==3)  fit->SetLineColor(kRed);
    //        fit->SetParameters(100,cherenkovreco[i],0.010,10);
    //        fit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    //        fit->SetParLimits(0,0.1,1E6);
    //        fit->SetParLimits(1,cherenkovreco[i]-2*cut_cangle,cherenkovreco[i]+2*cut_cangle);
    //        fit->SetParLimits(2,0.005,0.030); // width
    //        hAngle[i]->Fit("fgaus","I","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
    //        hAngle[i]->Fit("fgaus","M","",cherenkovreco[i]-cut_cangle,cherenkovreco[i]+cut_cangle);
    //
    //        cherenkovreco[i] = fit->GetParameter(1);
    //        spr[i] = fit->GetParameter(2);
    //    }
    //
    //    gStyle->SetOptTitle(0);
    //    gStyle->SetOptStat(0);
    //    gStyle->SetOptFit(0);
    //
    //
    //
    //
    //    TF1 *ff;
    //    double sep=0,esep=0, m1=0,m2=0,s1=0,s2=0;
    //    if(hLnDiff[3]->GetEntries()>100){
    //        hLnDiff[3]->Fit("gaus","S");
    //        ff = hLnDiff[3]->GetFunction("gaus");
    //        ff->SetLineColor(1);
    //        m1=ff->GetParameter(1);
    //        s1=ff->GetParameter(2);
    //    }
    //    if(hLnDiff[2]->GetEntries()>100){
    //        hLnDiff[2]->Fit("gaus","S");
    //        ff = hLnDiff[2]->GetFunction("gaus");
    //        ff->SetLineColor(1);
    //        m2=ff->GetParameter(1);
    //        s2=ff->GetParameter(2);
    //    }
    //    if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
    //
    //    hAngle[2]->GetXaxis()->SetRangeUser(0.7,0.9);
    //    hAngle[2]->GetYaxis()->SetRangeUser(0,1.2);
    //    hAngle[2]->Draw();
    //    hAngle[3]->Draw("same");
    //    // fAngle[3]->Draw("same");
    //    // fAngle[2]->Draw("same");
    //
    //
    //    TLine *line = new TLine(0,0,0,1000);
    //    //line->SetX1(mAngle[3]);
    //    //line->SetX2(mAngle[3]);
    //    line->SetX1(referance_angle_k);
    //    line->SetX2(referance_angle_k);
    //    line->SetY1(0);
    //    line->SetY2(1.2);
    //    line->SetLineColor(kRed);
    //    line->Draw();
    //
    //    TLine *line2 = new TLine(0,0,0,1000);
    //    //line2->SetX1(mAngle[2]);
    //    //line2->SetX2(mAngle[2]);
    //    line2->SetX1(referance_angle_pi);
    //    line2->SetX2(referance_angle_pi);
    //    line2->SetY1(0);
    //    line2->SetY2(1.2);
    //    line2->SetLineColor(kBlue);
    //    line2->Draw();
    //
    //    TLine *line3 = new TLine(0,0,0,1000);
    //    line3->SetLineStyle(2);
    //    //line3->SetX1(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
    //    //line3->SetX2(0.5*(mAngle[2]+mAngle[3])-cut_cangle);
    //    line3->SetX1(0.5*(referance_angle_k+referance_angle_pi)-cut_cangle);
    //    line3->SetX2(0.5*(referance_angle_k+referance_angle_pi)-cut_cangle);
    //    line3->SetY1(0);
    //    line3->SetY2(1.2);
    //    line3->SetLineColor(1);
    //    line3->Draw();
    //
    //    TLine *line4 = new TLine(0,0,0,1000);
    //    line4->SetLineStyle(2);
    //    //line4->SetX1(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
    //    //line4->SetX2(0.5*(mAngle[2]+mAngle[3])+cut_cangle);
    //    line4->SetX1(0.5*(referance_angle_k+referance_angle_pi)+cut_cangle);
    //    line4->SetX2(0.5*(referance_angle_k+referance_angle_pi)+cut_cangle);
    //    line4->SetY1(0);
    //    line4->SetY2(1.2);
    //    line4->SetLineColor(1);
    //    line4->Draw();
    //
    //
    //    TLegend *leg = new TLegend(0.1,0.5,0.4,0.85);
    //    leg->SetFillColor(0);
    //    leg->SetFillStyle(0);
    //    leg->SetBorderSize(0);
    //    leg->SetFillStyle(0);
    //    leg->AddEntry(hAngle[2],Form("#theta_{c}^{#pi} = %2.4f rad",cherenkovreco[2]),"");
    //    leg->AddEntry(hAngle[3],Form("#theta_{c}^{K} = %2.4f rad",cherenkovreco[3]),"");
    //    leg->AddEntry(hAngle[2],Form("#sigma_{c}^{#pi} = %2.1f mrad",spr[2]*1000),"");
    //    leg->AddEntry(hAngle[3],Form("#sigma_{c}^{K} = %2.1f mrad",spr[3]*1000),"");
    //    leg->Draw();
    //
    //    TLegend *lnpa = new TLegend(0.7,0.67,0.9,0.85);
    //    lnpa->SetFillColor(0);
    //    lnpa->SetFillStyle(0);
    //    lnpa->SetBorderSize(0);
    //    lnpa->SetFillStyle(0);
    //    lnpa->AddEntry(hAngle[2],"pions","lp");
    //    lnpa->AddEntry(hAngle[3],"kaons","lp");
    //    lnpa->Draw();
    //
    //    // fAngle[2]->Draw("same");
    //    // fAngle[3]->Draw("same");
    //
    //    glx_canvasAdd("hTime"+nid,800,400);
    //
    //    hTime->Draw();
    //    hCalc->SetLineColor(2);
    //    hCalc->Draw("same");
    //    TLegend *leg1 = new TLegend(0.5,0.6,0.85,0.80);
    //    leg1->SetFillColor(0);
    //    leg1->SetFillStyle(0);
    //    leg1->SetBorderSize(0);
    //    leg1->SetFillStyle(0);
    //    leg1->AddEntry(hTime,"measured in geant","lp");
    //    leg1->AddEntry(hCalc,"calculated","lp");
    //    leg1->Draw();
    //
    //    glx_canvasAdd("hDiff"+nid,800,400);
    //    hDiff->SetLineColor(kBlack);
    //    hDiff->Draw();
    //
    //    // hDiffT->SetLineColor(kRed+1);
    //    // hDiffT->Draw("same");
    //    hDiffD->SetLineColor(kGreen+2);
    //    hDiffD->Draw("same");
    //    hDiffR->SetLineColor(kBlue+1);
    //    hDiffR->Draw("same");
    //
    //    double maxTD= hDiffD->GetXaxis()->GetBinCenter(hDiffD->GetMaximumBin());
    //    double maxTR= hDiffR->GetXaxis()->GetBinCenter(hDiffR->GetMaximumBin());
    //    double maxTT= hTime->GetXaxis()->GetBinCenter(hTime->GetMaximumBin());
    //
    //    line = new TLine(0,0,0,1000);
    //    line->SetLineStyle(2);
    //    line->SetX1(-cut_tdiffd);
    //    line->SetX2(-cut_tdiffd);
    //    line->SetY1(0);
    //    line->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
    //    line->SetLineColor(1);
    //    line->Draw();
    //
    //    line2 = new TLine(0,0,0,1000);
    //    line2->SetLineStyle(2);
    //    line2->SetX1(cut_tdiffd);
    //    line2->SetX2(cut_tdiffd);
    //    line2->SetY1(0);
    //    line2->SetY2(hDiff->GetMaximum()+0.05*hDiff->GetMaximum());
    //    line2->SetLineColor(1);
    //    line2->Draw();
    //
    //    TLegend *leg2 = new TLegend(0.6,0.57,0.9,0.85);
    //    leg2->SetFillColor(0);
    //    leg2->SetFillStyle(0);
    //    leg2->SetBorderSize(0);
    //    leg2->SetFillStyle(0);
    //    leg2->AddEntry(hDiff,"all","lp");
    //    // leg2->AddEntry(hDiffT,"MC path in EV","lp");
    //    // leg2->AddEntry(hDiffD,"MC path in EV for direct photons","lp");
    //    // leg2->AddEntry(hDiffR,"MC path in EV for reflected photons","lp");
    //    leg2->AddEntry(hDiffD,"direct photons","lp");
    //    leg2->AddEntry(hDiffR,"reflected photons","lp");
    //
    //    leg2->Draw();
    //
    //    glx_canvasAdd("hLnDiff"+nid,800,400);
    //    hLnDiff[2]->SetTitle(Form("sep = %2.2f s.d.",sep));
    //    hLnDiff[3]->SetTitle(Form("sep = %2.2f s.d.",sep));
    //    hLnDiff[2]->Draw();
    //    hLnDiff[3]->Draw("same");
    //
    //    TLegend *lnpl = new TLegend(0.7,0.67,0.9,0.85);
    //    lnpl->SetFillColor(0);
    //    lnpl->SetFillStyle(0);
    //    lnpl->SetBorderSize(0);
    //    lnpl->SetFillStyle(0);
    //    lnpl->AddEntry(hLnDiff[2],"pions","lp");
    //    lnpl->AddEntry(hLnDiff[3],"kaons","lp");
    //    lnpl->Draw();
    //
    //    glx_canvasAdd("hNph"+nid,800,400);
    //
    //    double nph = 0;
    //    if(hNph[2]->GetEntries()>50){
    //        nph = glx_fit(hNph[2],40,100,40).X();
    //        auto rfit = hNph[2]->GetFunction("glx_gaust");
    //        if(rfit) rfit->SetLineColor(kBlue);
    //        hNph[2]->SetLineColor(kBlue);
    //        hNph[2]->Draw();
    //        //glx_fit(hNph[3],40,100,40).X();
    //        //hNph[3]->GetFunction("glx_gaust")->SetLineColor(kRed);
    //        hNph[3]->SetLineColor(kRed);
    //        hNph[3]->Draw("same");
    //    }
    
    
    
    
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
    
    
    //    TLegend *lnph = new TLegend(0.6,0.65,0.9,0.85);
    //    lnph->SetFillColor(0);
    //    lnph->SetFillStyle(0);
    //    lnph->SetBorderSize(0);
    //    lnph->SetFillStyle(0);
    //    // lnph->AddEntry(hNphC,"simulated","lp");
    //    lnph->AddEntry(hNph[2],"pions","lp");
    //    lnph->AddEntry(hNph[3],"kaons","lp");
    //
    //
    //    //lnph->AddEntry(hNph_p[2],"pions","lp");
    //    //lnph->AddEntry(hNph_p[3],"kaons","lp");
    //
    //    //lnph->AddEntry(hNph_n[2],"pions","lp");
    //    //lnph->AddEntry(hNph_n[3],"kaons","lp");
    //
    //    lnph->Draw();
    //
    //    std::cout<<" ###### separation = "<< sep << "  nph = "<<nph <<std::endl;
    //    std::cout<<"maxTD "<<maxTD<<"  maxTR "<<maxTR<<std::endl;
    
    
    
    
    
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
    
    /*
     glx_canvasAdd("26",800,400);
     histo_theta_phi_map_pi->Draw("colz");
     glx_canvasAdd("27",800,400);
     histo_theta_phi_map_k->Draw("colz");
     
     glx_canvasAdd("28",800,400);
     histo_theta_phi_mom_map_pi->Draw("colz");
     glx_canvasAdd("29",800,400);
     histo_theta_phi_mom_map_k->Draw("colz");
     
     
     
     glx_canvasAdd("13",800,400);
     hExtrapolatedBarHitXY_k->Draw("colz");
     
     glx_canvasAdd("14",800,400);
     hExtrapolatedBarHitXY_pi->Draw("colz");
     
     
     glx_canvasGet("14")->Update();
     TLine *linex1 = new TLine(0,0,0,1000);
     linex1->SetLineStyle(2);
     linex1->SetX1(xmin);
     linex1->SetX2(xmin);
     linex1->SetY1(gPad->GetUymin());
     linex1->SetY2(gPad->GetUymax());
     linex1->SetLineColor(kRed);
     linex1->Draw();
     
     glx_canvasGet("14")->Update();
     
     TLine *linex2 = new TLine(0,0,0,1000);
     linex2->SetLineStyle(2);
     linex2->SetX1(xmax);
     linex2->SetX2(xmax);
     linex2->SetY1(gPad->GetUymin());
     linex2->SetY2(gPad->GetUymax());
     linex2->SetLineColor(kGreen);
     linex2->Draw();
     
     glx_canvasGet("14")->Update();
     
     TLine *liney1 = new TLine(0,0,0,1000);
     liney1->SetLineStyle(2);
     liney1->SetX1(gPad->GetUxmin());
     liney1->SetX2(gPad->GetUxmax());
     liney1->SetY1(ymin);
     liney1->SetY2(ymin);
     liney1->SetLineColor(kBlack);
     liney1->Draw();
     
     glx_canvasGet("14")->Update();
     TLine *liney2 = new TLine(0,0,0,1000);
     liney2->SetLineStyle(2);
     liney2->SetX1(gPad->GetUxmin());
     liney2->SetX2(gPad->GetUxmax());
     liney2->SetY1(ymax);
     liney2->SetY2(ymax);
     liney2->SetLineColor(kBlue);
     liney2->Draw();
     
     cout<<"###########"<<xmin<<"	"<<xmax<<"	"<<ymin<<"	"<<ymax<<endl;
     */
    
    //glx_canvasGet("14")->Update();
    
    
    /*
     glx_canvasAdd("14",800,400);
     hExtrapolatedBarHitXY_cut->Draw("colz");
     */
    
    //    glx_canvasSave(0);
    //
    //
    //    if(gPDF_pix ==1) {
    //        for(Int_t i=0; i<glx_nch; i++) {
    //            fHistCh_k[i]->Write();
    //            fHistCh_pi[i]->Write();
    //        }
    //    }
    //    if(gCherenkov_Correction==1) {
    //        for(Int_t i=0; i<PMT_num; i++) {
    //            fHistPMT_k[i]->Write();
    //            fHistPMT_pi[i]->Write();
    //        }
    //    }
    //
    //
    //    if(gPDF_pmt==1) {
    //        for(Int_t i=0; i<PMT_num; i++) {
    //            fHistPMT_PDF_k[i]->Write();
    //            fHistPMT_PDF_pi[i]->Write();
    //        }
    //    }
    
    if (!btree){
        for(Int_t i=0; i<24; i++)
            for(Int_t j=0; j<40; j++)  {
                histo_time_bar_pos[i][j]->Write();
            }
        for(Int_t i=0; i<24; i++)
            for(Int_t j=0; j<40; j++)  {
                histo_tdiffD_bar_pos[i][j]->Write();
            }
        for(Int_t i=0; i<24; i++)
            for(Int_t j=0; j<40; j++)  {
                histo_tdiffR_bar_pos[i][j]->Write();
            }
    }
    
    
    //mom_theta_phi->Write();
    
    
    
    file.Write();
    file.Close();
    
    //  tree_cut->Print();
    //  tree_cut->AutoSave();
    
    
    cout << "##########  pion_counter = " << pion_counter <<"   "<<" Kaon_counter = "<<kaon_counter<<endl;
}
