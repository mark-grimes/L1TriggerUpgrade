{

	gROOT->ProcessLine(".X ~/trigger_studies/temp/CMSSW_5_3_1/src/UserCode/L1TriggerDPG/macros/initL1Analysis.C");
//	gROOT->ProcessLine(".X ~/trigger_studies/temp/CMSSW_5_3_1/src/UserCode/L1TriggerDPG/macros/upgrade/style_macro.C");
	cout << "Loading Upgrade Analysis Macro" << endl;
	gROOT->ProcessLine(".L ~/trigger_studies/temp/CMSSW_5_3_1/src/UserCode/L1TriggerDPG/macros/upgrade/UpgradeAnalysis_12.C+");
	UpgradeAnalysis_12 a;
	a.OpenWithList("inputFiles/inputFiles_ZBHPF1_UP_2012HPF_66_v3.txt");

   // Some variable definitions
   Int_t nEvents=-1;

	Int_t nLumis = gROOT->ProcessLine("a->FillDistros(nEvents,0,10066, 66, \"lumiStuff/getLumi_out_pixCorrLumi_66PU_stdCorr.txt\", 1)");
	
	gROOT->ProcessLine(".q");

}
