#ifndef DSelector_pimkpks__B3_M16_analyzer_h
#define DSelector_pimkpks__B3_M16_analyzer_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_pimkpks__B3_M16_analyzer : public DSelector
{
public:
    
    DSelector_pimkpks__B3_M16_analyzer(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_pimkpks__B3_M16_analyzer(){}
    
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
    // from workshop 2016
    TH1I* dHist_KinFitChiSq, *dHist_KinFitCL;
    TH2I* dHist_Proton_dEdx_P, *dHist_KPlus_dEdx_P, *dHist_PiPlus_dEdx_P, *dHist_PiMinus1_dEdx_P, *dHist_PiMinus2_dEdx_P;
    // DEFINE CUT PARAMETERS HERE
    TF1 *fFunc_dEdxCut_SelectHeavy;
    TF1 *fFunc_dEdxCut_SelectLight;
    double dMinKinFitCL, dMaxKinFitChiSq, dMinBeamEnergy, dMaxBeamEnergy, dMinKsMass, dMaxKsMass;
    
    ClassDef(DSelector_pimkpks__B3_M16_analyzer, 0);
};

void DSelector_pimkpks__B3_M16_analyzer::Get_ComboWrappers(void)
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

#endif // DSelector_pimkpks__B3_M16_analyzer_h
