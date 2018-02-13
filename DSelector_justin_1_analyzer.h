#ifndef DSelector_justin_1_analyzer_h
#define DSelector_justin_1_analyzer_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"
#include "DSelector/DComboTreeHelper.h"

#include "TH1I.h"
#include "TH2I.h"
#include <iostream>   // std::cout
#include <string>     // std::string, std::stod
class DSelector_justin_1_analyzer : public DSelector
{
public:
    
    DSelector_justin_1_analyzer(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_justin_1_analyzer(){}
    
    void Init(TTree *tree);
    Bool_t Process(Long64_t entry);
    
private:
    
    void Get_ComboWrappers(void);
    void Finalize(void);
    
    // BEAM POLARIZATION INFORMATION
    UInt_t dPreviousRunNumber;
    bool dIsPolarizedFlag; //else is AMO
    bool dIsPARAFlag; //else is PERP or AMO
    
    //CREATE REACTION-SPECIFIC PARTICLE ARRAYS
    
    //Step 0
    DParticleComboStep* dStep0Wrapper;
    DBeamParticle* dComboBeamWrapper;
    DChargedTrackHypothesis* dPiMinus1Wrapper;
    DChargedTrackHypothesis* dKPlusWrapper;
    DChargedTrackHypothesis* dProtonWrapper;
    
    //Step 1
    DParticleComboStep* dStep1Wrapper;
    DChargedTrackHypothesis* dPiMinus2Wrapper;
    DChargedTrackHypothesis* dPiPlusWrapper;
    
    // DEFINE YOUR HISTOGRAMS HERE
    // EXAMPLES:
    TH1I* dHist_MissingMassSquared;
    TH1I* dHist_BeamEnergy;
    TH1I* dHist_KsMass_Measured;
    TH1I* dHist_KsMass_KinFit;
    TH1I* dHist_RF, *dHist_RF_cut, *dHist_test;
    // from workshop 2016
    TH1I* dHist_KinFitChiSq, *dHist_KinFitCL;
    TH2I* dHist_Proton_dEdx_P, *dHist_KPlus_dEdx_P, *dHist_PiPlus_dEdx_P, *dHist_PiMinus1_dEdx_P, *dHist_PiMinus2_dEdx_P;
    // DEFINE CUT PARAMETERS HERE
    TF1 *fFunc_dEdxCut_SelectHeavy;
    TF1 *fFunc_dEdxCut_SelectLight;
    double dMinKinFitCL, dMaxKinFitChiSq, dMinBeamEnergy, dMaxBeamEnergy, dMinKsMass, dMaxKsMass;
    Int_t test_val;
    
    // TOOL FOR FLAT TREE OUTPUT
    DComboTreeHelper *dComboTreeHelper;
    
    ClassDef(DSelector_justin_1_analyzer, 0);
};

void DSelector_justin_1_analyzer::Get_ComboWrappers(void)
{
    //Step 0
    dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
    dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
    dPiMinus1Wrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
    dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
    dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
    
    //Step 1
    dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
    dPiMinus2Wrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
    dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_justin_1_analyzer_h
