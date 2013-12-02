//
//  Code performs scaling of menu thresholds for
//         a) Muon Isolation
//         b) Muon Pt Assignment
//         c) Data/MC Difference
//
//  Requirements:
//      Menu File: Defines names of trigger algorithms
//      Rate File: Trigger Rates and Thresholds (output by L1Menu2015.C)
//      Rate vs Threshold Root File:  Output of L1Menu2015.C for configuration being considered.
//
//    Next two are used for data mc scaling
//      Rate vs Threshold Root File for 8 TeV Data:  Output of L1Menu2015.C for configuration being considered.
//      Rate vs Threshold Root File for 8 Tev MC:    Output of L1Menu2015.C for configuration being considered.
//
//    L1 to Offline Conversion
//      L12Offline_*.txt  File contains linear fits for conversion of L1 threshold to offline threshold
//
//    Muon Scaling: (Root Files with histograms that are scale factors vs Muon Pt.
//      
//
//  Running:
//      linux>  root init.C  (loads libraries)
//      root[0]: .x runEvalScaling.C
//
//  Output:
//      L1RatesScaled_*.txt  (Scaled Rates and Thresholds and Offline Thresholds
//
runEvalScaling() {

   EvaluateL1Menu *myEval = new EvaluateL1Menu();

// Define Pieces to construct file names used for scaling.
   TString dirPlots  = "./";
	TString dirMenu   = "./";
   TString menuFile  =  dirMenu;
	menuFile         += "L1Menu_14TeV50PU_25nsMC_Fallback_TDRFinal_PS.txt";
	TString sampleKey = "14TeV100PU_50nsMC_"; 
	TString configKey = "Fallback_4p4e34Rates";
	TString menuKey   = "Fallback_4p4e34Rates";
	TString configKey8TeV = "_Fallback_1p6e34Rates";

// Files which have the linear fits for the L1 --> Offline threshold conversion
//	TString L12OffFile= "L12Offline_Current_95-85_LowPtMu_v5.txt";
	TString L12OffFile= "L12Offline_Fallback_95-85_LowPtMu_v5.txt";


/* These contain the scale factor for the muon upgrade
	   MuonScales_NoScale.root (Scale Factors = 1.0; use for Current Trigger)
		MuonScales_Scale_Iso50PU_3pt.root (Isolation and Pt Assignment Scale Factors for upgrade)
*/		
//	TString muonScaleFile = "MuonScales_NoScale.root";
	TString muonScaleFile = "MuonScales_Scale_Iso50PU_3pt.root";

	int nTrigAlg = 21;
	bool scaleDefault = true;
	bool suppressUnusedThr = true;


// Load the default menu
   myEval->LoadL1Menu(menuFile.Data());
	myEval->LoadL1toOffline(L12OffFile.Data());
	
// Load the standard histograms of rate vs threshold
    TString thrFile = dirPlots;
	 thrFile += "L1RateHist_";
	 thrFile += sampleKey;
	 thrFile += configKey;
	 thrFile += "Thr1_rates.root";
    myEval->LoadThresholdPlots(thrFile.Data());	

// Load Rate File   
    TString rateFile = dirMenu;
	 rateFile += "L1Rates_";
	 rateFile += sampleKey;
	 rateFile += menuKey;
	 if(scaleDefault) rateFile += "_Default";
	 rateFile += ".txt";
	 myEval->LoadL1Rates(rateFile.Data(),nTrigAlg);

// Define output, numerator, and denominator for scaling
    TString outThrFile = dirMenu;
	 outThrFile += "L1RateHist_";
	 outThrFile += sampleKey;
	 outThrFile += configKey;
	 outThrFile += configKey8TeV;
	 outThrFile += "Thr1_rates_Scaled.root";
	 
	 TString numeratorFile = dirPlots;
	 numeratorFile += "L1RateHist_8TeV66PUData";
	 numeratorFile += configKey8TeV;
	 numeratorFile += "Thr1_rates.root";

	 TString denominatorFile = dirPlots;
	 denominatorFile += "L1RateHist_8TeV66PU_NoOOTMC";
	 denominatorFile += configKey8TeV;
	 denominatorFile += "Thr1_rates.root";
	  
// Perfom the scaling	  
    myEval->RateScaling(outThrFile.Data(),
	 						   numeratorFile.Data(),
								denominatorFile.Data(),
								muonScaleFile.Data()); 

// Load the scaled histograms of rate vs threshold
    myEval->LoadThresholdPlots(outThrFile.Data());	

// Extract thresholds from scaled histograms
    myEval->FindThresholdsFromScaledRates();
	 
// Write out new Rate values
    TString outRatesFile = dirMenu;
	 outRatesFile += "L1RatesScaled_";
	 outRatesFile += sampleKey;
	 outRatesFile += menuKey;
	 if(scaleDefault) outRatesFile += "_Default";
	 outRatesFile += ".txt";
    myEval->WriteL1Rates(outRatesFile.Data(),suppressUnusedThr);	 
    
}
