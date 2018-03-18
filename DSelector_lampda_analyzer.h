#ifndef DSelector_lampda_analyzer_h
#define DSelector_lampda_analyzer_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_lampda_analyzer : public DSelector
{
	public:

		DSelector_lampda_analyzer(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_lampda_analyzer(){}

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
		DChargedTrackHypothesis* dKPlusWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;
    
    // DEFINE YOUR HISTOGRAMS HERE
    // EXAMPLES:
    TH1I* dHist_MissingMassSquared;
    TH1I* dHist_BeamEnergy;
    TH1I* dHist_LampdaMass_Measured;
    TH1I* dHist_LampdaMass_KinFit;
    TH1I* dHist_RF, *dHist_RF_cut, *dHist_test, *dHist_StepVertexZ, *dHist_DetachedPathLengthSignificance, *dHist_DetachedLifetime, *dHist_DetachedPathLength;
    // from worLampdahop 2016
    TH1I* dHist_KinFitChiSq, *dHist_KinFitCL;
    TH2I* dHist_Proton_dEdx_P;
    TH2I* dHist_StepVertexYVsX;
    
    TH2I* cartizian_theta_phi;
    TH2I* cartizian_theta_mom;
    // DEFINE CUT PARAMETERS HERE
    TF1 *fFunc_dEdxCut_SelectHeavy;
    TF1 *fFunc_dEdxCut_SelectLight;
    double dMinKinFitCL, dMaxKinFitChiSq, dMinBeamEnergy, dMaxBeamEnergy, dMinLampdaMass, dMaxLampdaMass;
    Int_t test_val;
    TEnv *env;
    
    double beamPhoton_RF_cut ;
    double simple_PathLength_cut;
    double beam_vertex_XYcut;
    double beam_vertex_Z1cut;
    double beam_vertex_Z2cut;
    double ChiSq_NDF_cut;
    double MissingMassSquared_cut;
    

	ClassDef(DSelector_lampda_analyzer, 0);
};

void DSelector_lampda_analyzer::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_lampda_analyzer_h
