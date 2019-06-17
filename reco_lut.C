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
    TCanvas *cc = new TCanvas("cc","cc",800,500);
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
    
    for(Int_t i=0; i<24; i++)
        for(Int_t j=0; j<40; j++)  {
            histo_time_bar_pos[i][j] = new TH1F(Form("histo_time_bar_pos_%d_%d",i,j),Form("histo_time_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 100,0,100);
            histo_tdiffD_bar_pos[i][j] = new TH1F(Form("histo_tdiffD_bar_pos_%d_%d",i,j),Form("histo_tdiffD_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 200,-50,50);
            histo_tdiffR_bar_pos[i][j] = new TH1F(Form("histo_tdiffR_bar_pos_%d_%d",i,j),Form("histo_tdiffR_bar_pos_%d_%d; # X Bar Hit [cm]; Bar number [#]",i,j), 200,-50,50);
        }
    
    //////////////////////////////
    /// cherenkove PDF per pix ///
    //////////////////////////////
    TH1F*  fHistCh_k[glx_nch], *fHistCh_pi[glx_nch], *fHistCh_read_k[glx_nch], *fHistCh_read_pi[glx_nch];
    TFile *ffile_cherenkov_pdf_pix;
    TString cherenkov_pdf_path_pix;
    
    for(Int_t i=0; i<glx_nch; i++) {
        fHistCh_k[i] = new TH1F(Form("fHistCh_k_%d",i),Form("fHistCh_k_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); //2000
        fHistCh_pi[i] = new TH1F(Form("fHistCh_pi_%d",i),Form("fHistCh_pi_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); //2000
    }
    // read pdf per pix
    if (gPDF_pix==2) {
        //cherenkov_data_k_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/pdf/histo_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        cherenkov_pdf_path_pix ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_cherenkovPDF_pix.root";
        cout<<"cherenkov_pdf_path_pix= " <<cherenkov_pdf_path_pix<<endl;
        ffile_cherenkov_pdf_pix  = new TFile(cherenkov_pdf_path_pix, "READ");
        for(Int_t pix=0; pix<glx_nch; pix++) {
            fHistCh_read_k[pix] = (TH1F*)ffile_cherenkov_pdf_pix->Get(Form("fHistCh_k_%d",pix));
            fHistCh_read_pi[pix] = (TH1F*)ffile_cherenkov_pdf_pix->Get(Form("fHistCh_pi_%d",pix));
        }
    }
    
    //////////////////////////////
    /// cherenkove PDF per pmt ///
    //////////////////////////////
    int PMT_num=108; //108
    //double norm = 20000;
    
    TH1F*  fHistPMT_PDF_k[PMT_num], *fHistPMT_PDF_pi[PMT_num], *fHistPMT_PDF_read_k[PMT_num], *fHistPMT_PDF_read_pi[PMT_num];
    TFile *ffile_cherenkov_pdf_pmt;
    TString cherenkov_pdf_path_pmt;
    
    for(Int_t i=0; i<PMT_num; i++) {
        fHistPMT_PDF_k[i] = new TH1F(Form("fHistPMT_PDF_k_%d",i),Form("fHistPMT_PDF_k_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1);    // 250
        fHistPMT_PDF_pi[i] = new TH1F(Form("fHistPMT_PDF_pi_%d",i),Form("fHistPMT_PDF_pi_%d;#theta_{C} [rad];entries [#]",i), 500,0.6,1); // 250
    }
    // read pdf per pmt
    if (gPDF_pmt==2) {
        //cherenkov_pdf_path_pmt ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_pdf_pmt_corrected.root";
        cherenkov_pdf_path_pmt ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/created_cherenkovPDF_pmt_NotCorrected.root";
        cout<<"cherenkov_pdf_path_pmt= " <<cherenkov_pdf_path_pmt<<endl;
        ffile_cherenkov_pdf_pmt  = new TFile(cherenkov_pdf_path_pmt, "READ");
        for(Int_t PMT=0; PMT<PMT_num; PMT++) {
            fHistPMT_PDF_read_k[PMT] = (TH1F*)ffile_cherenkov_pdf_pmt->Get(Form("fHistPMT_PDF_k_%d",PMT));
            fHistPMT_PDF_read_pi[PMT] = (TH1F*)ffile_cherenkov_pdf_pmt->Get(Form("fHistPMT_PDF_pi_%d",PMT));
            
            // No success
            /*
             double scale_k = norm/(fHistPMT_PDF_read_k[PMT]->Integral());
             fHistPMT_PDF_read_k[PMT]->Scale(scale_k);
             double scale_pi = norm/(fHistPMT_PDF_read_pi[PMT]->Integral());
             fHistPMT_PDF_read_k[PMT]->Scale(scale_pi);
             */
        }
    }
    
    
    ////////////////////////////
    /// cherenkove correction///
    ////////////////////////////
    double referance_angle = mAngle[2]; // pi
    //double referance_angle = mAngle[3]; // k
    double referance_angle_pi = mAngle[2];
    double referance_angle_k = mAngle[3];
    TGraph *shifted_pi = new TGraph();
    
    shifted_pi->SetMarkerColor(kBlue);
    shifted_pi->SetMarkerStyle(20);
    shifted_pi->SetLineColor(kBlue);
    shifted_pi->SetLineWidth(1);
    
    TH1F*  fHistPMT_k[PMT_num], *fHistPMT_pi[PMT_num], *fHistPMT_read_k[PMT_num], *fHistPMT_read_pi[PMT_num];
    TFile *ffile_cherenkov_correction;
    TString cherenkov_correction_path;
    // correction array
    double array_correction[108]={0};
    // creat histograms Cherenkov per PMT
    for(Int_t i=0; i<PMT_num; i++) {
        fHistPMT_k[i] = new TH1F(Form("fHistPMT_k_%d",i),Form("fHistPMT_k_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
        fHistPMT_pi[i] = new TH1F(Form("fHistPMT_pi_%d",i),Form("fHistPMT_pi_%d;#theta_{C} [rad];entries [#]",i), 250,0.6,1);
    }
    // Read Cherenkov per PMT
    TF1 *fit_PMT = new TF1("fit_PMT","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",minChangle,maxChangle);
    fit_PMT->SetLineColor(kBlue);
    TH1F*  hdiff_test = new TH1F("hdiff_test",";mean diff ;entries [#]", 100,0,0.1);
    TH1F*  hsigma_test = new TH1F("hsigma_test",";sigam ;entries [#]", 100,0,20);
    TH1F*  hmean_test = new TH1F("hmean_test",";mean ;entries [#]", 250,0.6,1);
    //gStyle->SetPalette(kLightTemperature);
    //////////////////////////////////////////////////
    //    glx_canvasAdd("r_pmt_correction",800,400);//
    /////////////////////////////////////////////////
    if (gCherenkov_Correction==2) {
        cherenkov_correction_path ="/lustre/nyx/panda/aali/gluex/gluex_top/hdgeant4/hdgeant4-2.1.0/macro/dirc/cherenkov_correction.root";//outFile_separation_PDF_pmt.root //cherenkov_correction.root
        cout<<"cherenkov_correction_path= " <<cherenkov_correction_path<<endl;
        ffile_cherenkov_correction  = new TFile(cherenkov_correction_path, "READ");
        int pmtCounter =0;
        for(Int_t PMT=0; PMT<PMT_num; PMT++) {
            TString pmt_counter=Form("_%d",pmtCounter);
            if(PMT<=10 || (PMT>=90 && PMT<=96)) continue; // dummy pmts
            if((PMT ==13 || PMT==14 || PMT==31 || PMT==32 || PMT==33 || PMT==12  || PMT==15 || PMT==34  || PMT==35 || PMT==35 || PMT==16 || PMT==17 || PMT==11 || PMT==102)) continue;
            //if(! (PMT==102)) continue; // custmization
            ////////////////////////////////////////////////////////////////////////
            //            glx_canvasAdd("r_pmt_correction"+pmt_counter,800,400);//
            ////////////////////////////////////////////////////////////////////////
            fHistPMT_read_k[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistPMT_k_%d",PMT));
            fHistPMT_read_pi[PMT] = (TH1F*)ffile_cherenkov_correction->Get(Form("fHistPMT_pi_%d",PMT));
            fit_PMT->SetParameters(100,0.82,0.010,10);
            fit_PMT->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fit_PMT->SetParLimits(0,0.1,1E6);
            //fit_PMT->SetParLimits(1,0.82-2*cut_cangle,0.82+2*cut_cangle);
            fit_PMT->SetParLimits(1,0.809,0.835);
            fit_PMT->SetParLimits(2,0.005,0.030); // width
            
            double  rang_min= 0.82-cut_cangle;
            double  rang_max= 0.82+cut_cangle;
            
            //////////////////////////////////
            /// custumize fitting function ///
            //////////////////////////////////
            
            if(PMT==42) rang_min=0.81;
            if(PMT==24 || PMT==26 ||PMT==18 ||PMT==11){
                rang_min= 0.82-cut_cangle/2;
                rang_max= 0.82+cut_cangle/2;
            }
            if(PMT==107){
                rang_min= 0.82-cut_cangle/2;
                rang_max= 0.82+cut_cangle;
                
            }
            if(PMT==40 || PMT==41){
                rang_min= 0.814;
                rang_max= 0.84;
            }
            
            fHistPMT_read_pi[PMT]->Fit("fit_PMT","M","",rang_min,rang_max);
            double histo_cor_entries = fHistPMT_read_pi[PMT]->GetEntries();
            double mean_cherenkov_cor=  fit_PMT->GetParameter(1);
            double sigma_cherenkov_cor= fit_PMT->GetParameter(2);
            double delta_cherenkov_cor= referance_angle - mean_cherenkov_cor;
            double val_1 = (delta_cherenkov_cor)*1000.0 ;
            shifted_pi->SetPoint(pmtCounter, mean_cherenkov_cor, PMT);
            
            cout<<"##########"<<"PMT= "<<PMT<< "	delta_cherenkov_cor= " << delta_cherenkov_cor<<endl;
            //hsigma_test->Fill(sigma_cherenkov_cor*1000);
            //hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
            //hmean_test->Fill(mean_cherenkov_cor);
            // correction condition
            // if( ( fabs(delta_cherenkov_cor) <0.016 && (sigma_cherenkov_cor*1000<16 && sigma_cherenkov_cor*1000> 5.5) && histo_cor_entries>0 )) {//0.01  12  7
            array_correction[PMT]= delta_cherenkov_cor;
            cout<<"##########"<< "shift "<<val_1<<endl;
            /*
             // Fill PMT
             for(Int_t m=0; m<64; m++)
             for(Int_t n=0; n<64; n++){
             glx_hdigi[PMT]->SetBinContent(m, n,val_1);
             }
             
             // default correction
             //array_correction[PMT]= -0.004;
             fHistPMT_read_pi[PMT]->Draw();
             //glx_canvasGet("r_pmt_correction"+pmt_counter)->Update();
             glx_canvasGet("r_pmt_correction")->Update();
             TLine *lin_ref = new TLine(0,0,0,1000);
             lin_ref->SetX1(referance_angle);
             lin_ref->SetX2(referance_angle);
             lin_ref->SetY1(gPad->GetUymin());
             lin_ref->SetY2(gPad->GetUymax());
             lin_ref->SetLineColor(kBlue);
             lin_ref->Draw();
             
             TLine *lin_ref_k = new TLine(0,0,0,1000);
             lin_ref_k->SetX1(referance_angle_k);
             lin_ref_k->SetX2(referance_angle_k);
             lin_ref_k->SetY1(gPad->GetUymin());
             lin_ref_k->SetY2(gPad->GetUymax());
             lin_ref_k->SetLineColor(kRed);
             lin_ref_k->Draw();
             
             //glx_canvasGet("r_pmt_correction"+pmt_counter)->Update();
             glx_canvasGet("r_pmt_correction")->Update();
             glx_waitPrimitive("r_pmt_correction");
             //hsigma_test->Fill(sigma_cherenkov_cor*1000);
             //hdiff_test->Fill(fabs(mean_cherenkov_cor-referance_angle));
             //hmean_test->Fill(mean_cherenkov_cor);
             */
            // }
            ++pmtCounter;
        }
        /*
         glx_canvasAdd("r_pmt_shift",800,400);
         TMultiGraph *mg = new TMultiGraph();
         mg->Add(shifted_pi);
         mg->SetTitle(" Shift ; Mean [rad]; PMT ID [#]");
         mg->Draw("AP");
         //mg->GetHistogram()->GetYaxis()->SetRangeUser(6800,7050);
         
         glx_canvasGet("r_pmt_shift")->Update();
         
         TLine *lin_ref_pi = new TLine(0,0,0,1000);
         lin_ref_pi->SetX1(referance_angle_pi);
         lin_ref_pi->SetX2(referance_angle_pi);
         lin_ref_pi->SetY1(gPad->GetUymin());
         lin_ref_pi->SetY2(gPad->GetUymax());
         lin_ref_pi->SetLineColor(kBlue);
         lin_ref_pi->Draw();
         
         glx_canvasGet("r_pmt_shift")->Update();
         TLine *lin_ref_k = new TLine(0,0,0,1000);
         lin_ref_k->SetX1(referance_angle_k);
         lin_ref_k->SetX2(referance_angle_k);
         lin_ref_k->SetY1(gPad->GetUymin());
         lin_ref_k->SetY2(gPad->GetUymax());
         lin_ref_k->SetLineColor(kRed);
         lin_ref_k->Draw();
         */
    }
    
    
    
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
    double shift_tdiff[24][40]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0484555,-0.0354493,0.35056,0.099063,0.106277,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.128436,0.0400649,0.105171,0.0283179,0.136207,0.129174,0.0509148,0.0244918,-0.0016855,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.219853,-0.0737354,0.159231,0.146284,0.123299,0.0353285,0.057264,0.0669868,-0.0431749,0.124274,0.11793,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0443389,-0.010351,0.144584,0.167535,0.129673,0.0894336,0.0149637,0.0334035,0.0312673,0.134101,0.157287,0.163244,0,0,0,0,0,0,0,0,0,0,0,0,-0.0109966,-0.220774,0.169731,0.212142,0.16085,0.126701,0.135627,0.121442,0.0499737,0.122402,0.128811,0.15737,0.11934,0,0,0,0,0,0,0,0,0,0,0,-0.00782742,-0.150989,0.186154,0.21918,0.178342,0.136228,0.123202,0.154193,0.133912,0.163065,0.181085,0.171055,0.0391149,0.0889785,0,0,0,0,0,0,0,0,0,0,0.0835393,0.150198,0.177662,0.183543,0.137075,0.153928,0.181749,0.15866,0.154047,0.19948,0.232563,0.226626,0.0704589,0.086086,0.121518,0.0662947,0,0,0,0,0,0,0,0,0.0841535,0.145702,0.206464,0.169638,0.141065,0.172967,0.174205,0.144694,0.184069,0.199037,0.244997,0.210197,0.0292002,0.0355092,0.0290452,0.0534557,-0.0351165,0,0,0,0,0,0,0,0.085724,0.164682,0.2165,0.167911,0.150495,0.183209,0.169329,0.166109,0.177069,0.227666,0.241485,0.184138,0.0027776,0.0243227,0.00654651,-0.0541308,-0.00404641,0,0,0,0,0,0,0,0.107842,0.19749,0.252776,0.219025,0.185657,0.190596,0.204788,0.17564,0.167063,0.233824,0.253876,0.208454,-0.0317229,0.0565071,-0.00870096,-0.032719,0.00610309,-0.0738033,0,0,0,0,0,0,0.107041,0.252879,0.311105,0.280969,0.215917,0.247506,0.221,0.20466,0.186726,0.234477,0.237015,0.170027,-0.0439488,0.0759515,0.0869067,-0.0135607,0.116568,0.0471306,0.0723177,0,0,0,0,0,0.141757,0.261828,0.319707,0.302958,0.257017,0.256147,0.200369,0.161661,0.125224,0.160671,0.196845,0.14449,-0.0520265,0.11576,0.121756,0.122954,0.163013,0.111975,0.0484454,0.194809,0,0,0,0,0.171518,0.275187,0.307965,0.291526,0.266654,0.230735,0.172033,0.114582,0.0966193,0.137031,0.183111,0.0891136,0.012301,0.127834,0.183852,0.179489,0.204186,0.164188,0.192083,0.106122,0,0,0,0,0.256979,0.337382,0.367115,0.34409,0.244315,0.199646,0.16355,0.116371,0.121487,0.169231,0.175652,0.149268,0.0205011,0.166876,0.222973,0.21423,0.235252,0.180907,0.227778,0.177711,0.292843,0,0,0,0.423643,0.465037,0.479693,0.410611,0.263701,0.242166,0.21457,0.166087,0.16705,0.192644,0.204521,0.191679,0.0869033,0.206789,0.250065,0.272298,0.241632,0.222175,0.268539,0.293759,0.33817,0,0,0,0.725575,0.740799,0.715482,0.621035,0.40466,0.354143,0.324878,0.299437,0.260641,0.279625,0.307112,0.265067,0.163861,0.280992,0.321935,0.355474,0.301971,0.283339,0.373811,0.457592,0.417735,0,0,0,0.950818,0.951673,0.890637,0.745973,0.596449,0.601227,0.465558,0.411163,0.345117,0.368678,0.366612,0.303725,0.196065,0.314516,0.358048,0.329245,0.294593,0.297833,0.345911,0.403051,0.513864,0.38563,0,0,1.36867,1.12123,0.94363,0.839209,0.652887,0.628982,0.540399,0.452434,0.36597,0.411794,0.384593,0.309495,0.148578,0.281512,0.263998,0.229866,0.226089,0.222082,0.343845,0.420341,0.442428,0,0,0,1.92236,1.36569,1.09846,0.887148,0.712891,0.64572,0.545014,0.42746,0.319841,0.378124,0.372523,0.284111,0.128976,0.207101,0.194574,0.16628,0.2005,0.155013,0.259195,0.375682,0.470856,0.407513,0,0,2.01906,1.41681,1.13046,0.896505,0.709332,0.648193,0.517753,0.406326,0.304463,0.345361,0.373136,0.294617,0.0993746,0.188603,0.173755,0.179623,0.195836,0.191101,0.287064,0.400992,0.466,0.361539,0,0,1.4621,1.29693,1.09323,0.831755,0.644312,0.597216,0.490958,0.364789,0.270081,0.315619,0.353508,0.290061,0.11985,0.182116,0.199869,0.160848,0.202368,0.23682,0.322385,0.382495,0.434683,0.405445,0,0,1.14012,1.14301,0.976661,0.79192,0.616518,0.539383,0.474992,0.352494,0.266265,0.309797,0.33003,0.280832,0.145202,0.187874,0.18328,0.196861,0.213956,0.177515,0.273324,0.352289,0.397228,0.37937,0,0,0.827087,0.890285,0.786642,0.664625,0.510495,0.470278,0.404346,0.290087,0.211788,0.268287,0.312761,0.246882,0.126446,0.230824,0.225923,0.179073,0.214367,0.193096,0.261107,0.345096,0.330653,0,0,0,0.515694,0.644616,0.632142,0.542591,0.403364,0.381829,0.294817,0.202702,0.137065,0.221989,0.252789,0.199705,0.105985,0.190114,0.218622,0.212496,0.19969,0.142444,0.220413,0.236997,0.324881,0,0,0,0.312724,0.473409,0.442828,0.412404,0.309941,0.284752,0.220669,0.060527,0.0779171,0.164957,0.215385,0.137745,0.00932066,0.150187,0.191998,0.161405,0.191306,0.0578354,0.224053,0.27969,0.293328,0,0,0,0.177287,0.348568,0.358936,0.307365,0.233995,0.221168,0.204349,0.0965681,0.0239412,0.177434,0.251703,0.125934,0.0330784,0.179473,0.153984,0.167659,0.183589,0.116312,0.237547,0.201601,0,0,0,0,-0.0669073,0.137761,0.171553,0.174695,0.144268,0.159986,0.136592,0.0128444,-0.00955643,0.145326,0.21197,0.149466,-0.0294701,0.101494,0.0940754,0.0808618,0.0282518,0.00325814,0.173595,0.208221,0,0,0,0,-0.231522,-0.179675,-0.0590653,0.0142431,0.0257128,0.0939398,0.105458,-0.0330371,-0.0311701,0.0289467,0.173189,0.0227273,-0.101745,-0.0506557,-0.0609564,-0.0104289,0.0358258,-0.048497,0.189429,0,0,0,0,0,-0.319261,-0.288284,-0.175551,-0.106451,-0.134527,0.00727297,0.0393707,-0.00510716,-0.00323208,0.117724,0.171612,-0.0613118,-0.165411,-0.0928218,-0.140565,-0.110293,-0.0336569,-0.0909058,0,0,0,0,0,0,-0.432025,-0.369398,-0.275893,-0.232665,-0.324574,-0.197646,-0.052544,0.0754222,0.0867731,0.134527,0.0342633,-0.0946818,-0.255951,-0.0927508,-0.0671486,-0.043222,-0.0636675,0,0,0,0,0,0,0,-0.575479,-0.480157,-0.332118,-0.337863,-0.495491,-0.28086,-0.0680482,0.0643614,0.0254752,0.163175,0.109449,-0.196485,-0.257506,-0.0974613,0.00714947,0.0879513,-0.0235813,0,0,0,0,0,0,0,-0.71696,-0.546338,-0.49057,-0.476121,-0.52483,-0.288839,-0.124298,-0.0956534,-0.091014,0.108539,0.183561,-0.225335,-0.330447,-0.00225432,-0.0179166,0,0,0,0,0,0,0,0,0,-0.818417,-0.736128,-0.737093,-0.723545,-0.53314,-0.330092,-0.403355,-0.293202,-0.24218,-0.128041,0.116102,-0.261308,-0.457386,-0.0902638,0,0,0,0,0,0,0,0,0,0,-0.894983,-0.863325,-1.0398,-0.75881,-0.516174,-0.331565,-0.447565,-0.439813,-0.412139,-0.117445,0.148536,-0.251335,0,0,0,0,0,0,0,0,0,0,0,0,-0.71129,-1.09435,-1.02932,-0.753744,-0.490395,-0.364743,-0.44317,-0.60414,-0.586372,-0.197612,0.0720365,-0.288913,0,0,0,0,0,0,0,0,0,0,0,0,-0.184685,-1.10635,-0.780496,-0.584552,-0.534895,-0.210287,-0.255855,-0.509604,-0.509893,-0.41877,0.14744,0,0,0,0,0,0,0,0,0,0,0,0,0,6.54575,0.0786484,-0.380539,-0.503242,-0.477557,-0.0942726,-0.109767,-0.542088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //double shift_tdiff[24][40]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.206375,-0.120977,-0.0632259,-0.0429772,0.00747347,-0.0146532,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.182306,-0.232246,-0.176029,0.0152093,-0.0418729,-0.0245196,-0.0133792,-0.0537889,-0.11011,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.105011,-0.0782876,0.00766007,0.0348279,-0.0191352,0.0217075,-0.044096,-0.0512151,-0.0294671,0.0567046,0.0841607,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0951189,0.0165561,-0.0559123,0.144276,0.00961639,-0.0119877,-0.0032614,-0.00882342,-0.0275106,-0.00472598,-0.00213266,0.073259,0,0,0,0,0,0,0,0,0,0,0,0,-0.070211,-0.0293339,0.129513,0.145992,0.0456247,0.0120037,-0.0105453,0.0260046,0.00662134,0.0349205,0.109562,0.109924,0,0,0,0,0,0,0,0,0,0,0,0,-0.022451,0.0817875,0.113035,0.120863,0.0210432,0.0835056,0.116326,0.0861171,0.0819289,0.124119,0.161047,0.164751,-0.0700882,-0.0365541,0,0,0,0,0,0,0,0,0,0,-0.0233332,0.0862023,0.141587,0.105897,0.0236786,0.108331,0.104791,0.0906122,0.114128,0.127726,0.181853,0.141345,-0.0873594,-0.0414573,-0.0924231,0,0,0,0,0,0,0,0,0,0.0176569,0.105414,0.155845,0.102785,0.0827252,0.108978,0.10364,0.100162,0.118493,0.15658,0.17382,0.113746,-0.104513,-0.0361405,-0.0468388,-0.116422,0,0,0,0,0,0,0,0,0.0470503,0.136206,0.189396,0.15627,0.12275,0.132528,0.140132,0.114471,0.103533,0.17471,0.186681,0.149281,-0.0975165,-0.0622213,-0.0758128,-0.106883,-0.0906574,0,0,0,0,0,0,0,0.0350743,0.185329,0.250495,0.21582,0.155364,0.187518,0.153889,0.139597,0.122184,0.174551,0.173563,0.0997633,-0.114964,-0.028553,-0.0386911,-0.0835532,0.0464993,-0.0853314,0,0,0,0,0,0,0.0682792,0.199881,0.258163,0.239962,0.19224,0.196294,0.145619,0.102063,0.00399847,0.0957003,0.130271,0.0843635,-0.111692,-0.000247798,0.0055788,0.0144015,0.037596,0.003136,0,0,0,0,0,0,0.101523,0.209164,0.248916,0.233186,0.196796,0.169061,0.106062,-0.0123911,-0.0273042,0.077524,0.118173,0.0296335,-0.0570873,0.00240015,0.126288,0.119882,0.137874,0.0978324,0.133246,0,0,0,0,0,0.192402,0.275101,0.30872,0.280945,0.185974,0.140566,0.100579,-0.00243232,-0.000675846,0.10292,0.108936,0.0247199,-0.0420769,0.101744,0.151405,0.146548,0.168609,0.0509927,0.160777,0,0,0,0,0,0.357617,0.40312,0.417512,0.34748,0.201496,0.177813,0.153136,0.105488,0.101254,0.128386,0.136196,0.126364,-0.0250021,0.145194,0.185316,0.137221,0.126271,0.0980243,0.204516,0.217792,0,0,0,0,0.656715,0.680235,0.652814,0.500618,0.339135,0.291788,0.25842,0.238501,0.195156,0.217345,0.246385,0.19604,0.0348437,0.223829,0.26325,0.281999,0.23499,0.210771,0.298814,0.379574,0,0,0,0,0.883868,0.886061,0.821912,0.681943,0.474418,0.468443,0.402174,0.34882,0.284041,0.308701,0.307586,0.24706,0.138672,0.255879,0.294335,0.267102,0.239704,0.231805,0.274282,0.368345,0,0,0,0,1.30366,1.0509,0.878242,0.775212,0.587467,0.510003,0.47769,0.389596,0.305923,0.347708,0.320987,0.250304,0.0825788,0.218644,0.198415,0.168737,0.157952,0.0736943,0.281943,0.372571,0,0,0,0,1.85218,1.2975,0.992593,0.820532,0.645445,0.528325,0.476375,0.366335,0.257308,0.312215,0.307261,0.217868,0.0628457,0.145219,0.125537,0.102688,0.132491,0.0805325,0.192931,0.29885,0,0,0,0,1.95276,1.3382,1.02442,0.829505,0.640167,0.530475,0.454082,0.3424,0.237855,0.282463,0.307611,0.229658,-0.0183328,0.125308,0.109486,0.112934,0.130076,0.0515368,0.227634,0.332633,0.406381,0,0,0,1.38598,1.2341,0.98184,0.767929,0.534755,0.481209,0.425954,0.298727,0.206334,0.252278,0.284897,0.222402,0.0517366,0.12054,0.129048,0.0951738,0.0827893,0.161261,0.255696,0.322637,0.386345,0,0,0,1.04886,1.07806,0.902163,0.725041,0.50761,0.470327,0.407936,0.288444,0.200643,0.242038,0.265568,0.216087,0.028602,0.12573,0.116091,0.121369,0.14353,0.0554943,0.216787,0.2889,0.323208,0,0,0,0.759166,0.824731,0.717299,0.599859,0.438278,0.399089,0.334255,0.226427,0.145044,0.197248,0.244311,0.177303,0.0625023,0.172256,0.159303,0.11263,0.141516,0.053383,0.113653,0.282784,0,0,0,0,0.442804,0.571157,0.519125,0.438746,0.337019,0.314312,0.228918,0.0853099,0.0241698,0.155896,0.189898,0.133067,-0.0067844,0.084858,0.156968,0.149541,0.131152,0.00874533,0.0963787,0.169875,0,0,0,0,0.242857,0.404417,0.370217,0.344157,0.240356,0.217933,0.152823,-0.00255466,-0.0300135,0.0514462,0.152049,0.042377,-0.0555858,0.0918079,0.0837694,0.106264,0.0690981,-0.0122815,0.0840268,0.131297,0,0,0,0,0.0983326,0.277509,0.281437,0.238535,0.167504,0.153848,0.135528,-0.00984144,-0.0434906,0.109989,0.184929,0.0576361,-0.0676725,0.111528,0.0839548,0.105897,0.0729967,-0.0294737,0.169866,0,0,0,0,0,-0.137766,0.067839,0.0993935,0.103282,0.0775381,0.0943456,0.0272432,-0.0575783,-0.0799284,0.0295568,0.156486,0.0393061,-0.0969438,-0.0138838,-0.0203956,-0.0310712,-0.0375451,-0.0529255,0.115933,0,0,0,0,0,-0.293959,-0.241706,-0.132339,-0.0577204,-0.0772528,-0.0114359,-0.00884542,-0.10428,-0.0928766,-0.030287,0.105486,-0.0439717,-0.172406,-0.116128,-0.13064,-0.125332,-0.0374852,-0.113614,0,0,0,0,0,0,-0.382473,-0.362927,-0.234509,-0.167022,-0.199928,-0.10027,-0.0288441,-0.0762747,-0.0624127,0.00566043,0.0470891,-0.113922,-0.240428,-0.159455,-0.197585,-0.187921,-0.1187,0,0,0,0,0,0,0,-0.543142,-0.485574,-0.33891,-0.290941,-0.390789,-0.270084,-0.108127,-0.0296511,-0.0154228,0.0207321,-0.0235342,-0.163968,-0.332387,-0.172064,-0.139598,-0.120121,0,0,0,0,0,0,0,0,-0.63883,-0.551125,-0.392767,-0.386559,-0.614114,-0.358374,-0.141736,-0.00555223,-0.033874,0.0415662,0.0300131,-0.267633,-0.342674,-0.162464,-0.0725263,0,0,0,0,0,0,0,0,0,-0.766841,-0.600636,-0.556217,-0.546904,-0.592107,-0.353903,-0.242189,-0.169991,-0.15218,0.0356908,0.112256,-0.305968,-0.410471,-0.13413,0,0,0,0,0,0,0,0,0,0,-0.864668,-0.805071,-0.784181,-0.797985,-0.590788,-0.408068,-0.498623,-0.4176,-0.318743,-0.171815,0.0551932,-0.344196,0,0,0,0,0,0,0,0,0,0,0,0,-0.911755,-0.929319,-1.13094,-0.809612,-0.575254,-0.383256,-0.505917,-0.497968,-0.479322,-0.228721,0.0510376,-0.310328,0,0,0,0,0,0,0,0,0,0,0,0,-0.574165,-1.10109,-1.11278,-0.791123,-0.569908,-0.416169,-0.509241,-0.679135,-0.647672,-0.255733,-0.0209098,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.22699,-1.19949,-0.853074,-0.611401,-0.672023,-0.275643,-0.337414,-0.575414,-0.58357,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.58213,0.0936823,-0.434353,-0.553283,-0.581149,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    /////////////////////////////////
    //////// Creat file and trees ///
    /////////////////////////////////
    
    TFile file(outFile,"recreate");
    //TTree *tree_cut = glx_ch->CloneTree(0);
    
    TTree tree_variables("tree_variables","tree for cherenkov track resolution");
    double track_spr(-1),track_mean(-1), track_yield(-1), track_mom(-1), track_xbar(0),track_ybar(0),track_fit_chisqu(-1),track_fit_NDF(-1);
    int track_pid(-1), track_nbar(-1);

    TString track_file="noname";

    std::vector<int> vpx;
    std::vector<int> vpy;
    std::vector<int> vpz;

    std::vector<double> vtdiff;
    //    std::vector<double> vtime;
    std::vector<double> vtangle;

    double track_inv_mass(-1),track_missing_mass(-1),track_chi_square(-1),track_TofTrackDist(-1);

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
    //    tree_variables.Branch("vtime",&vtime);
    tree_variables.Branch("vtangle",&vtangle);


    tree_variables.Branch("track_inv_mass",&track_inv_mass,"track_inv_mass/D");
    tree_variables.Branch("track_missing_mass",&track_missing_mass,"track_missing_mass/D");
    tree_variables.Branch("track_chi_square",&track_chi_square,"track_chi_square/D");
    tree_variables.Branch("track_TofTrackDist",&track_TofTrackDist,"track_TofTrackDist/D");

    
    
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
            //            if(true){
            //                if (pdgId == 2 && chi_square> 10) continue;
            //                if (pdgId == 3 && chi_square> 20)continue;
            //                if (pdgId == 2 && (inv_mass< 0.66  || inv_mass> 0.82 )  ) continue;
            //                if (pdgId == 3 && (inv_mass< 1.015 || inv_mass> 1.025) ) continue;
            //                if(missing_mass< -0.01 || missing_mass> 0.01 )continue;
            //            }
            
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
            track_spr=-1; track_mean=-1; track_yield=-1; track_mom=-1; track_xbar=0;track_ybar=0;
            track_pid=-1; track_nbar=-1;
            histo_cherenkov_track->Reset();
            vpx.clear();
            vpy.clear();
            vpz.clear();

            vtdiff.clear();
            //vtime.clear();
            vtangle.clear();
            
            
            for(int h = 0; h < glx_event->GetHitSize(); h++){
                hit = glx_event->GetHit(h);
                int ch = hit.GetChannel();
                int pmt = hit.GetPmtId();
                int pix = hit.GetPixelId();
                double hitTime = hit.GetLeadTime()-glx_event->GetTime();
                hitTime=hitTime+shift_tdiff[bar][x_pos_bin];
                if(ch>glx_nch) continue;
                //histo_time_bar_pos[bar][x_pos_bin]->Fill(hitTime);
                //if(hitTime>40) continue;
                nphc++;
                
                /////////////////////////////////////
                // Reflection condition may change //
                /////////////////////////////////////
                //bool reflected = hitTime>40;
                bool reflected = true;
                if(hitTime < cop[x_pos_bin]) reflected = false;
                
                
                
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
                            if(gCherenkov_Correction != 2) tangle = momInBar.Angle(dir);
                            if(gCherenkov_Correction == 2) tangle = momInBar.Angle(dir) + array_correction[pmt];
                            
                            //double bartime = lenz/cos(luttheta)/20.4; //198 //203.767 for 1.47125
                            double bartime = lenz/cos(luttheta)/19.6; //203.767 for 1.47125
                            double totalTime = bartime+evtime;
                            // hTime->Fill(hitTime);
                            // hCalc->Fill(totalTime);
                            
                            if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))<0.01){
                                hDiff->Fill(totalTime-hitTime);
                                //if(samepath)
                                {
                                    hDiffT->Fill(totalTime-hitTime);
                                    if(r) {
                                        hDiffR->Fill(totalTime-hitTime);
                                        histo_tdiffR_bar_pos[bar][x_pos_bin]->Fill(totalTime-hitTime);
                                        
                                    }
                                    else {
                                        hDiffD->Fill(totalTime-hitTime);
                                        histo_tdiffD_bar_pos[bar][x_pos_bin]->Fill(totalTime-hitTime);
                                        
                                    }
                                    
                                    
                                }
                                
                                if(fabs(totalTime-hitTime)< 5)histo_time_bar_pos[bar][x_pos_bin]->Fill(hitTime);
                            }
                            // skim
                            if(fabs(totalTime-hitTime)> 10) continue;
                            if(tangle > 1.0) continue;
                            if(tangle > 1.0) continue;
                            
                            vtdiff.push_back(totalTime-hitTime);
                            //vtime.push_back(hitTime);
                            vtangle.push_back(tangle);
                            
                            ///////////////
                            // Time Cut  //
                            ///////////////
                            if(!r && fabs(totalTime-hitTime)>cut_tdiffd) continue; // removed
                            if(r && fabs(totalTime-hitTime) >cut_tdiffr) continue; // removed
                            
                            //////////////////////
                            // Cherenkov track  //
                            //////////////////////
                            
                            histo_cherenkov_track->Fill(tangle);
                            
                            ////////////////
                            // Fill PDF   //
                            ////////////////
                            // cherenkove PDF per PIX
                            if(gPDF_pix ==1 && pdgId == 3) fHistCh_k[ch]->Fill(tangle); // good after time cut
                            if(gPDF_pix ==1 && pdgId == 2) fHistCh_pi[ch]->Fill(tangle); // good after time cut
                            
                            // cherenkove PDF per PMT
                            if(gPDF_pmt ==1 && pdgId == 3) fHistPMT_PDF_k[pmt]->Fill(tangle); // good after time cut
                            if(gPDF_pmt ==1 && pdgId == 2) fHistPMT_PDF_pi[pmt]->Fill(tangle); // good after time cut
                            
                            ////////////////////////////
                            // Fill PMT coorrection   //
                            ////////////////////////////
                            // cherenkove correction per PMT
                            if(gCherenkov_Correction ==1 && pdgId == 3) fHistPMT_k[pmt]->Fill(tangle);
                            if(gCherenkov_Correction ==1 && pdgId == 2) fHistPMT_pi[pmt]->Fill(tangle);
                            
                            // fill cherenkove histo
                            hAngle[pdgId]->Fill(tangle);
                            
                            ////////////////////
                            // Cherenkov Cut  //
                            ////////////////////
                            if(fabs(tangle-0.5*(mAngle[2]+mAngle[3]))>cut_cangle) continue;
                            //if(fabs(tangle-0.5*(referance_angle_k+referance_angle_pi))>cut_cangle) continue; // removed
                            //if(tangle> 0.844 ||tangle < 0.798)  continue;
                            
                            isGood=true;
                            hTime->Fill(hitTime);
                            hCalc->Fill(totalTime);
                            
                            
                            if(!(gPDF_pix==2||gPDF_pmt==2 )){ // fixed
                                sum1 += TMath::Log(fAngle[2]->Eval(tangle)+noise);
                                sum2 += TMath::Log(fAngle[3]->Eval(tangle)+noise);
                            }
                            
                            if(gPDF_pix ==2){
                                // use histograms
                                Int_t kk = fHistCh_read_k[ch]->GetXaxis()->FindBin(tangle);
                                Int_t kpi = fHistCh_read_pi[ch]->GetXaxis()->FindBin(tangle);
                                if (fHistCh_read_pi[ch]->GetBinContent(kpi) > 0 )sum1 += TMath::Log(fHistCh_read_pi[ch]->GetBinContent(kpi));
                                if (fHistCh_read_k[ch]->GetBinContent(kk) > 0 )sum2 += TMath::Log(fHistCh_read_k[ch]->GetBinContent(kk));
                                
                                //if (sum1 != 0 || sum2!=0 )std::cout<<"No Problem  separation  " <<kpi<<" "<<kk<<"  sum "<<sum1 <<"  "<< sum2<<std::endl;
                                //std::cout<<"###### No Problem  separation  " << fHistCh_read_k[ch]->GetBinContent(kk) <<"  "<< fHistCh_read_pi[ch]->GetBinContent(kpi)<<std::endl;
                            }
                            
                            if(gPDF_pmt ==2){
                                // use histograms
                                Int_t k_bin = fHistPMT_PDF_read_k[pmt]->GetXaxis()->FindBin(tangle);
                                Int_t pi_bin = fHistPMT_PDF_read_pi[pmt]->GetXaxis()->FindBin(tangle);
                                if (fHistPMT_PDF_read_pi[pmt]->GetBinContent(pi_bin) > 0 )sum1 += TMath::Log(fHistPMT_PDF_read_pi[pmt]->GetBinContent(pi_bin));
                                if (fHistPMT_PDF_read_k[pmt]->GetBinContent(k_bin) > 0 )sum2 += TMath::Log(fHistPMT_PDF_read_k[pmt]->GetBinContent(k_bin));
                                
                                
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
                        
                        vpx.push_back(pmt);
                        vpy.push_back(pix%8);
                        vpz.push_back(pix/8);
                        
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
            
            /////////////////////
            // fill tree here////
            /////////////////////
            
            fit_track->SetParameters(100,0.82,0.010,10);
            fit_track->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fit_track->SetParLimits(0,0.1,1E6);
            fit_track->SetParLimits(1,0.809,0.835);
            fit_track->SetParLimits(2,0.005,0.030);
            if (pdgId==3)fit_track->SetLineColor(kRed);
            else fit_track->SetLineColor(kBlue);
            histo_cherenkov_track->Fit("fit_track","MQ0","", 0.5*(referance_angle_k+referance_angle_pi)-cut_cangle, 0.5*(referance_angle_k+referance_angle_pi)-cut_cangle) ;

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

            tree_variables.Fill();


            ///////////////////////////////////
            //////// reduce pions number //////
            ///////////////////////////////////
            //double percentage = kaon_counter/pion_counter*100.0;
            //if (percentage <100 && pdgId == 2 )continue;
            //if (pdgId == 2) pion_counter++;
            //if (pdgId == 3) kaon_counter++;
            
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
    hLnDiff[3]->SetTitle(Form("sep = %2.2f s.d.",sep));
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
    
    glx_canvasSave(0);
    
    
    if(gPDF_pix ==1) {
        for(Int_t i=0; i<glx_nch; i++) {
            fHistCh_k[i]->Write();
            fHistCh_pi[i]->Write();
        }
    }
    if(gCherenkov_Correction==1) {
        for(Int_t i=0; i<PMT_num; i++) {
            fHistPMT_k[i]->Write();
            fHistPMT_pi[i]->Write();
        }
    }
    
    
    if(gPDF_pmt==1) {
        for(Int_t i=0; i<PMT_num; i++) {
            fHistPMT_PDF_k[i]->Write();
            fHistPMT_PDF_pi[i]->Write();
        }
    }
    
   
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
    
    
    //mom_theta_phi->Write();
    
    
    
    file.Write();
    file.Close();
    
    //  tree_cut->Print();
    //  tree_cut->AutoSave();
    
    
    cout << "##########  pion_counter = " << pion_counter <<"   "<<" Kaon_counter = "<<kaon_counter<<endl;
}
