#include "DSelector_justin_1_analyzer.h"



void DSelector_justin_1_analyzer::Init(TTree *locTree)
{
    // USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.
    
    // The Init() function is called when the selector needs to initialize a new tree or chain.
    // Typically here the branch addresses and branch pointers of the tree will be set.
    // Init() will be called many times when running on PROOF (once per file to be processed).
    
    //USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
    dOutputFileName = "justin_1_analyzer.root"; //"" for none
    dOutputTreeFileName = ""; //"" for none
    dFlatTreeFileName = "justin_1_analyzer_flat.root"; //output flat tree (one combo per tree entry), "" for none
    dFlatTreeName = ""; //if blank, default name will be chosen
    
    
    //test_val = dOption.Atoi();
    
    env = new TEnv(dOption);
    if (!dOption) return;
    test_val=9;
    test_val = env->GetValue("test_val", test_val);
    
    //    std::string orbits ("-6");
    //    string::size_type sz;     // alias of size_t
    //    test_val = std::stod (orbits,&sz);
    
    
    //Because this function gets called for each TTree in the TChain, we must be careful:
    //We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
    bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
    DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
    //gDirectory now points to the output file with name dOutputFileName (if any)
    if(locInitializedPriorFlag)
        return; //have already created histograms, etc. below: exit
    
    Get_ComboWrappers();
    dPreviousRunNumber = 0;
    
    double dTargetCenterZ;
    dTargetCenterZ = dParticleComboWrapper->Get_TargetCenter().Z();
    
    /*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/
    
    //ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
    //false/true below: use measured/kinfit data
    
    //PID
    //dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
    //below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
    //dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));
    //dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.0, Unknown, SYS_NULL));
    //dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, true, 2.0, Unknown, SYS_NULL));
    
    
    
    dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, "pid_precut"));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, PiPlus, SYS_BCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.5, PiPlus, SYS_FCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, PiPlus, SYS_TOF));
    
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.4, PiMinus, SYS_BCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.5, PiMinus, SYS_FCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.2, PiMinus, SYS_TOF));
    
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.75, KPlus, SYS_BCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.7, KPlus, SYS_FCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.30, KPlus, SYS_TOF));
    
    
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.4, Proton, SYS_BCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.5, Proton, SYS_FCAL));
    dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.2, Proton, SYS_TOF));
    dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, "pid_postcut"));
    
    //void DHistogramAction_ParticleID::Create_Hists(int locStepIndex, Particle_t locPID, string locStepROOTName)
    //void DHistogramAction_ParticleComboKinematics::Create_Hists(int locStepIndex, string locStepROOTName, Particle_t locPID, bool locIsBeamFlag)
    
    dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false ,"test"));
    
    //MASSES
    //dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
    //dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));
    
    //KINFIT RESULTS
    dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));
    
    //CUT MISSING MASS
    dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.04, 0.04));
    
    //BEAM ENERGY
    dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
    //dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));
    
    //KINEMATICS
    dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
    
    //INITIALIZE ACTIONS
    //If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
    Initialize_Actions();
    
    /******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/
    
    //EXAMPLE MANUAL HISTOGRAMS:
    dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
    dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
    dHist_KsMass_Measured = new TH1I("KsMass_Measured", ";#pi^{#plus}#pi^{#minus} Invariant Mass", 50, 0.470, 0.525);
    dHist_KsMass_KinFit = new TH1I("KsMass_KinFit", ";#pi^{#plus}#pi^{#minus} Invariant Mass", 50, 0.470, 0.525);
    
    //added from workshop 2016
    dHist_Proton_dEdx_P = new TH2I("Proton_dEdx_P", " ;p_{proton} GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
    //    dHist_KPlus_dEdx_P = new TH2I("KPlus_dEdx_P", " ;p_K^{#plus} GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
    //    dHist_PiPlus_dEdx_P = new TH2I("PiPlus_dEdx_P", " ;p_#pi^{#plus} GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
    //    dHist_PiMinus1_dEdx_P = new TH2I("PiMinus1_dEdx_P", " ;p_#pi^{#minus} 1st GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
    //    dHist_PiMinus2_dEdx_P = new TH2I("PiMinus2_dEdx_P", " ;p_#pi^{#minus} 2nd GeV/c; dE/dx (keV/cm)", 250, 0.0, 5.0, 250, 0.0, 25.);
    
    dHist_KinFitChiSq = new TH1I("KinFitChiSq", ";Kinematic Fit #chi^{2}/NDF", 250, 0., 25.);
    dHist_KinFitCL = new TH1I("KinFitCL", ";Kinematic Fit Confidence Level", 100, 0., 1.);
    
    dHist_RF=new TH1I("dHist_RF", ";#Deltat_{Beam#gamma - RF}", 1000, -10, 10);
    dHist_RF_cut=new TH1I("dHist_RF_cut", ";#Deltat_{Beam#gamma - RF}", 1000, -10, 10);
    dHist_test=new TH1I("dHist_test", ";dHist_test", 100, -10, 10);
    
    
    // EXAMPLE CUT PARAMETERS:
    fFunc_dEdxCut_SelectHeavy = new TF1("fFunc_dEdxCut_SelectHeavy", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.); // dFunc_dEdxCut_SelectHeavy
    fFunc_dEdxCut_SelectHeavy->SetParameters(4.0, 2.5, 1.25);
    fFunc_dEdxCut_SelectLight = new TF1("fFunc_dEdxCut_SelectLight", "exp(-1.*[0]*x + [1]) + [2]", 0., 10.);// dFunc_dEdxCut_SelectLight
    fFunc_dEdxCut_SelectLight->SetParameters(4.0, 2.0, 2.5);
    dMinKinFitCL = 0.0; //5.73303e-7;
    dMaxKinFitChiSq = 5.0;
    dMinBeamEnergy = 8.4;
    dMaxBeamEnergy = 9.0;
    dMinKsMass = 0.757;
    dMaxKsMass = 0.807;
    
    /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/
    
    //EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
    //The type for the branch must be included in the brackets
    //1st function argument is the name of the branch
    //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
    /*
     dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
     dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
     dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
     dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
     dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
     */
    
    //CREATE HELPER AND INITIALIZE WITH DESIRED COMBINATIONS TO BE STORED
    //dComboTreeHelper = new DComboTreeHelper(dComboWrapper, dFlatTreeInterface, "K+ K- p; K+ K-; K- p");
    
    /************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/
    
    //EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
    //The type for the branch must be included in the brackets
    //1st function argument is the name of the branch
    //2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
    /*
     dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
     dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
     dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
     dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
     */
    
    /************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/
    
    //TO SAVE PROCESSING TIME
    //If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
    //By default, for each event, the data is retrieved for all branches
    //If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
    //Do this by doing something similar to the commented code below
    
    //dTreeInterface->Clear_GetEntryBranches(); //now get none
    //dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_justin_1_analyzer::Process(Long64_t locEntry)
{
    // The Process() function is called for each entry in the tree. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    //
    // This function should contain the "body" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    // Use fStatus to set the return value of TTree::Process().
    // The return value is currently not used.
    
    //CALL THIS FIRST
    DSelector::Process(locEntry); //Gets the data from the tree for the entry
    //cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
    TLorentzVector locProductionX4 = Get_X4_Production();
    
    /******************************************** GET POLARIZATION ORIENTATION ******************************************/
    
    //Only if the run number changes
    //RCDB environment must be setup in order for this to work! (Will return false otherwise)
    UInt_t locRunNumber = Get_RunNumber();
    if(locRunNumber != dPreviousRunNumber)
    {
        dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
        dPreviousRunNumber = locRunNumber;
    }
    
    /********************************************* SETUP UNIQUENESS TRACKING ********************************************/
    
    //ANALYSIS ACTIONS: Reset uniqueness tracking for each action
    //For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
    Reset_Actions_NewEvent();
    
    //PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
    //Sometimes, some content is the exact same between one combo and the next
    //e.g. maybe two combos have different beam particles, but the same data for the final-state
    //When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
    //So, for each quantity you histogram, keep track of what particles you used (for a given combo)
    //Then for each combo, just compare to what you used before, and make sure it's unique
    
    //EXAMPLE 1: Particle-specific info:
    set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search
    set<Int_t> locUsedSoFar_Proton;
    set<Int_t> locUsedSoFar_KPlus;
    set<Int_t> locUsedSoFar_PiPlus;
    set<Int_t> locUsedSoFar_PiMinus1;
    set<Int_t> locUsedSoFar_PiMinus2;
    
    //EXAMPLE 2: Combo-specific info:
    //In general: Could have multiple particles with the same PID: Use a set of Int_t's
    //In general: Multiple PIDs, so multiple sets: Contain within a map
    //Multiple combos: Contain maps within a set (easier, faster to search)
    
    //INSERT USER ANALYSIS UNIQUENESS TRACKING HERE
    set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;
    set<map<Particle_t, set<Int_t> > > locUsedSoFar_KsMass;
    
    /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/
    
    /*
     Int_t locMyInt = 7;
     dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);
     
     TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
     dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);
     
     for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
     dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
     */
    
    /************************************************* LOOP OVER COMBOS *************************************************/
    
    //Loop over combos
    for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
    {
        //Set branch array indices for combo and all combo particles
        dComboWrapper->Set_ComboIndex(loc_i);
        
        // Is used to indicate when combos have been cut
        if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
            continue; // Combo has been cut previously
        
        /********************************************** GET PARTICLE INDICES *********************************************/
        
        //Used for tracking uniqueness when filling histograms, and for determining unused particles
        
        //Step 0
        Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
        Int_t locPiMinus1TrackID = dPiMinus1Wrapper->Get_TrackID();
        Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
        Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
        
        //Step 1
        Int_t locPiMinus2TrackID = dPiMinus2Wrapper->Get_TrackID();
        Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
        
        /*********************************************** GET FOUR-MOMENTUM **********************************************/
        
        // Get P4's: //is kinfit if kinfit performed, else is measured
        //dTargetP4 is target p4
        //Step 0
        TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
        TLorentzVector locPiMinus1P4 = dPiMinus1Wrapper->Get_P4();
        TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
        TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
        //Step 1
        TLorentzVector locPiMinus2P4 = dPiMinus2Wrapper->Get_P4();
        TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
        
        // Get Measured P4's:
        //Step 0
        
        TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
        TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
        TLorentzVector locPiMinus1P4_Measured = dPiMinus1Wrapper->Get_P4_Measured();
        TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
        TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
        //Step 1
        TLorentzVector locPiMinus2P4_Measured = dPiMinus2Wrapper->Get_P4_Measured();
        TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
        
        /********************************************* COMBINE FOUR-MOMENTUM ********************************************/
        
        // DO YOUR STUFF HERE
        TLorentzVector locKsP4_Measured = locPiPlusP4_Measured + locPiMinus2P4_Measured;
        TLorentzVector locKsP4 = locPiPlusP4 + locPiMinus2P4;
        
        
        
        // Combine 4-vectors
        TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
        locMissingP4_Measured -= locPiMinus1P4_Measured + locKPlusP4_Measured + locProtonP4_Measured + locPiMinus2P4_Measured + locPiPlusP4_Measured;
        /******************************************** Test histo *******************************************/
        dHist_test->Fill(test_val);
        /******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/
        
        // Loop through the analysis actions, executing them in order for the active particle combo
        if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
            continue;
        cout<<"######################### dOption=  "<<dOption<<endl;
        //if you manually execute any actions, and it fails a cut, be sure to call:
        //dComboWrapper->Set_IsComboCut(true);
        
        double beamPhoton_RF =locBeamX4_Measured.T() - locProductionX4.T();
        dHist_RF->Fill(beamPhoton_RF);
        
        if (fabs(beamPhoton_RF) > 2) continue;
        dHist_RF_cut->Fill(beamPhoton_RF);
        /**************************************** EXAMPLE: PID dEdx CUT ACTION ************************************************/
        
        // Proton CDC dE/dx histogram and cut
        double locProton_dEdx_CDC = dProtonWrapper->Get_dEdx_CDC()*1e6;
        double locKPlus_dEdx_CDC = dPiPlusWrapper->Get_dEdx_CDC()*1e6;
        double locPiPlus_dEdx_CDC = dPiPlusWrapper->Get_dEdx_CDC()*1e6;
        double locPiMinus1_dEdx_CDC = dPiMinus1Wrapper->Get_dEdx_CDC()*1e6;
        double locPiMinus2_dEdx_CDC = dPiMinus2Wrapper->Get_dEdx_CDC()*1e6;
        
        // remove the compo which dose not pass the dEdx cuts
        
        if(locProton_dEdx_CDC < fFunc_dEdxCut_SelectLight->Eval(locProtonP4_Measured.P())) {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }
        if(locUsedSoFar_Proton.find(locProtonTrackID) == locUsedSoFar_Proton.end())
        {
            dHist_Proton_dEdx_P->Fill(locProtonP4_Measured.P(), locProton_dEdx_CDC);
            locUsedSoFar_Proton.insert(locProtonTrackID);
        }
        
        //        if(locKPlus_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locKPlusP4_Measured.P())) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            //continue;
        //        }
        //        if(locUsedSoFar_KPlus.find(locKPlusTrackID) == locUsedSoFar_KPlus.end())
        //        {
        //            dHist_KPlus_dEdx_P->Fill(locKPlusP4_Measured.P(), locKPlus_dEdx_CDC);
        //            locUsedSoFar_KPlus.insert(locKPlusTrackID);
        //        }
        //
        //        if(locPiPlus_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiPlusP4_Measured.P())) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            //continue;
        //        }
        //        if(locUsedSoFar_PiPlus.find(locPiPlusTrackID) == locUsedSoFar_PiPlus.end())
        //        {
        //            dHist_PiPlus_dEdx_P->Fill(locPiPlusP4_Measured.P(), locPiPlus_dEdx_CDC);
        //            locUsedSoFar_PiPlus.insert(locPiPlusTrackID);
        //        }
        //        if(locPiMinus1_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiMinus1P4_Measured.P())) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            //continue;
        //        }
        //        if(locUsedSoFar_PiMinus1.find(locPiMinus1TrackID) == locUsedSoFar_PiMinus1.end())
        //        {
        //            dHist_PiMinus1_dEdx_P->Fill(locPiMinus1P4_Measured.P(), locPiMinus1_dEdx_CDC);
        //            locUsedSoFar_PiMinus1.insert(locPiMinus1TrackID);
        //        }
        //        if(locPiMinus2_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiMinus2P4_Measured.P())) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            //continue;
        //        }
        //        if(locUsedSoFar_PiMinus2.find(locPiMinus2TrackID) == locUsedSoFar_PiMinus2.end())
        //        {
        //            dHist_PiMinus2_dEdx_P->Fill(locPiMinus2P4_Measured.P(), locPiMinus2_dEdx_CDC);
        //            locUsedSoFar_PiMinus2.insert(locPiMinus2TrackID);
        //        }
        //
        //
        //
        //        if(locProton_dEdx_CDC < fFunc_dEdxCut_SelectLight->Eval(locProtonP4_Measured.P())
        //           || locKPlus_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locKPlusP4_Measured.P())
        //           || locPiPlus_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiPlusP4_Measured.P())
        //           || locPiMinus1_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiMinus1P4_Measured.P())
        //           || locPiMinus2_dEdx_CDC > fFunc_dEdxCut_SelectHeavy->Eval(locPiMinus2P4_Measured.P()) )
        //        {
        //
        //            dComboWrapper->Set_IsComboCut(true);
        //            //continue;
        //        }
        
        /**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/
        
        /*
         TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
         //for arrays below: 2nd argument is value, 3rd is array index
         //NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
         //So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
         dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
         dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
         */
        
        /**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/
        
        //Histogram beam energy (if haven't already)
        if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
        {
            dHist_BeamEnergy->Fill(locBeamP4.E());
            locUsedSoFar_BeamEnergy.insert(locBeamID);
        }
        
        // beam energy cut for SDME
        //        if(locBeamP4.E() < dMinBeamEnergy || locBeamP4.E() > dMaxBeamEnergy) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            continue;
        //        }
        /************************************** HIST, CUT KINFIT CONFIDENCE LEVEL ****************************************/
        
        
        // kinematic fit CL cut
        dHist_KinFitChiSq->Fill(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit());
        dHist_KinFitCL->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        //        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < dMinKinFitCL) {
        //            dComboWrapper->Set_IsComboCut(true);
        //            continue;
        //        }
        
        if(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit()> 4) {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }
        /************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/
        
        //Missing Mass Squared
        double locMissingMassSquared = locMissingP4_Measured.M2();
        
        //Uniqueness tracking: Build the map of particles used for the missing mass
        //For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
        map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
        locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
        locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus1TrackID);
        locUsedThisCombo_MissingMass[KPlus].insert(locKPlusTrackID);
        locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
        locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinus2TrackID);
        locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
        
        //compare to what's been used so far
        if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
        {
            //unique missing mass combo: histogram it, and register this combo of particles
            dHist_MissingMassSquared->Fill(locMissingMassSquared);
            locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
        }
        
        //E.g. Cut
        if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }
        
        
        /**************************************** HISTOGRAM Ks INVARIANT MASS *****************************************/
        
        double locKsMass_Measured = locKsP4_Measured.M();
        double locKsMass_KinFit = locKsP4.M();
        
        
        //Uniqueness tracking:
        map<Particle_t, set<Int_t> > locUsedThisCombo_KsMass;
        locUsedThisCombo_KsMass[PiMinus].insert(locPiMinus2TrackID);
        locUsedThisCombo_KsMass[PiPlus].insert(locPiPlusTrackID);
        
        
        //compare to what's been used so far
        if(locUsedSoFar_KsMass.find(locUsedThisCombo_KsMass) == locUsedSoFar_KsMass.end())
        {
            //unique missing mass combo: histogram it, and register this combo of particles
            dHist_KsMass_Measured->Fill(locKsMass_Measured);
            dHist_KsMass_KinFit->Fill(locKsMass_KinFit);
            locUsedSoFar_KsMass.insert(locUsedThisCombo_KsMass);
        }
        
        /****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/
        
        /*
         //FILL ANY CUSTOM BRANCHES FIRST!!
         Int_t locMyInt_Flat = 7;
         dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);
         
         TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
         dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);
         
         for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
         {
         dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
         TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
         dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
         }
         */
        
        //FILL FLAT TREE
        //dComboTreeHelper->Fill(loc_i); //fills branches for sub-combinations
        //Fill_FlatTree(); //for the active combo
    } // end of combo loop
    
    //FILL HISTOGRAMS: Num combos / events surviving actions
    Fill_NumCombosSurvivedHists();
    
    /******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
    /*
     //Thrown beam: just use directly
     if(dThrownBeam != NULL)
     double locEnergy = dThrownBeam->Get_P4().E();
     
     //Loop over throwns
     for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
     {
     //Set branch array indices corresponding to this particle
     dThrownWrapper->Set_ArrayIndex(loc_i);
     
     //Do stuff with the wrapper here ...
     }
     */
    /****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
    /*
     //Loop over beam particles (note, only those appearing in combos are present)
     for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
     {
     //Set branch array indices corresponding to this particle
     dBeamWrapper->Set_ArrayIndex(loc_i);
     
     //Do stuff with the wrapper here ...
     }
     
     //Loop over charged track hypotheses
     for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
     {
     //Set branch array indices corresponding to this particle
     dChargedHypoWrapper->Set_ArrayIndex(loc_i);
     
     //Do stuff with the wrapper here ...
     }
     
     //Loop over neutral particle hypotheses
     for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
     {
     //Set branch array indices corresponding to this particle
     dNeutralHypoWrapper->Set_ArrayIndex(loc_i);
     
     //Do stuff with the wrapper here ...
     }
     */
    
    /************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
    /*
     Bool_t locIsEventCut = true;
     for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
     //Set branch array indices for combo and all combo particles
     dComboWrapper->Set_ComboIndex(loc_i);
     // Is used to indicate when combos have been cut
     if(dComboWrapper->Get_IsComboCut())
     continue;
     locIsEventCut = false; // At least one combo succeeded
     break;
     }
     if(!locIsEventCut && dOutputTreeFileName != "")
     Fill_OutputTree();
     */
    
    return kTRUE;
}

void DSelector_justin_1_analyzer::Finalize(void)
{
    //Save anything to output here that you do not want to be in the default DSelector output ROOT file.
    
    //Otherwise, don't do anything else (especially if you are using PROOF).
    //If you are using PROOF, this function is called on each thread,
    //so anything you do will not have the combined information from the various threads.
    //Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.
    
    //DO YOUR STUFF HERE
    
    //CALL THIS LAST
    DSelector::Finalize(); //Saves results to the output file
}
