//
//  Script depends on the following libraries being loaded:
//          L1Ntupe_C.so
//          L1Menu2015_C.so
//          EvaluateL1Menu_C.so
//  (The libraries are loaded with init.C)
//
//  example:
//
//   linux> root init.C
//   root [1] .x runL1Menu2015.C
//
//  Inputs (Defined below):
//  -----------------------
//  L1MenuFile                Text File defining the L1 Menu algorithms  			 
//  InputRootFileName         ROOT file containing the L1Tree
//  lsFileName 			      Text File containing the luminosity information on a lumi-section by lumi-section basis.
//  ZeroBiasPrescale 	      ZeroBias prescale used for this sample. For data probably the product of L1 and HLT.  For MC probably 1
//  NumberOfLumiSections      Number of LumiSections in L1Tree. Used for rate normalizations.
//  InstLumi                  Instantaneous luminosity (used if no lsFile)
//  jobTag                    String used for naming of output files.
//  targetLumi                Target Luminosity used for scaling.  Given in units of e32.
//  numEvts                   Number of events to run over. (-1 = all events in file)
//  makeThrPlots              Flag for making rate vs Threshold plots.  0=Skip;  1=make 1-d plots; 2=make 1-d and 2-d plots
//  selectDataInput           Defines which tree to use for information (1=current;  22=upgrade)
//  thresholdPlotsFile        ROOT file containing the rate vs Threshold plots (if not remaking them)
//  targetTotalRate			   Target total L1 Menu rate.
//  targetTolerance           Allowable tolerance on the target rate.  If too small, may not converge.
//  maxIter                   Maximum interation to converge on target rate. 
//
// 
// ==================================================================================================================

// Define the input parameters
	TString L1MenuFile           = "L1Menu_v14m20_std.txt";

// Input file
   std::string InputRootFileName= "/obsidian/users/winer/cms/L1Trigger/Simulations/14TeV/version14/CurrentL1_ZeroBias_8TeV_PU45/L1Tree_All.root";

//	
	TString lsFileName           = "getLumi_out_45_stdCorr_v4.txt"; //""; //
	Float_t ZeroBiasPrescale     = 92.; //1.;
	Float_t NumberOfLumiSections = 143.; //float(1000)/(23.3570304*11246.*1368.); //63.;	
	Float_t InstLumi             = 0.061; //200.;  //units of e32
	TString jobTag               = "8TeV45PUData_1p1E34_Current_v14m20";		
	float targetLumi             = 110.; //units of e32  
	int  numEvts                 = -1; //Number of events to process (-1 == all in file)
	int  makeThrPlots            = 1;
	Int_t selectDataInput        = 1;
	TString thresholdPlotsFile   = "L1RateHist_8TeV45PUData_1p1E34_Current_v14m20Thr1_rates.root"; //Set a default. If makeThrPlot != 0 it will be reset
	Float_t targetTotalRate      = 95.;
	Float_t targetTolerance      = 5.;
	int maxIter                  = 10;
	Float_t totalRate;
	EvaluateL1Menu *myEval;
	L1Menu2015 *myL1Menu;

runL1Menu2015(){
  
   gSystem->Exec("date"); //beginning time stamp
	
// Instantiate the L1 Menu Evaluation Code and open the ROOT file
   myL1Menu = new L1Menu2015(0,targetLumi,NumberOfLumiSections,InstLumi,InputRootFileName,0.,ZeroBiasPrescale,false);
   myL1Menu->Open(InputRootFileName);

// Set flag for whether we use upgraded muons (Dummy...not implemented yet)	
	myL1Menu->UseUpgradeMuons(false);
	
// Turn off Dumping of Events
   myL1Menu->SetDumpEvents(0);	
	
// If threshold plots need to be made, Run the job and create them
   totalRate = 0.;
   if(makeThrPlots>0) {

// Dump Some Events to screen	
	    myL1Menu->SetDumpEvents(0);
		 
       printf("\n ---> Making Rate vs Threshold Plots \n");
       totalRate = myL1Menu->RunL1Ana(L1MenuFile,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);

// Set the file name to correspond to the rates just calculated
       thresholdPlotsFile = "L1RateHist_";
       thresholdPlotsFile += jobTag;
       thresholdPlotsFile += "Thr";
       thresholdPlotsFile += makeThrPlots;
       thresholdPlotsFile += "_rates.root";
       
// Turn off Dumping of Events for Iteration Phase
       myL1Menu->SetDumpEvents(0);

   }
 
// Instantiate the code that will adjust thresholds
    myEval = new EvaluateL1Menu();

// Load the rate vs threshold plots
    myEval->LoadThresholdPlots(thresholdPlotsFile);
	
// Load the L1 Menu	
    myEval->LoadL1Menu(L1MenuFile);

// Determine the thresholds
    myEval->DetermineThresholds();
		
// Write out a temporary L1 Menu file with new thresholds 
    TString tmpMenu = "L1Menu_Tmp_";
	 tmpMenu += jobTag;
	 tmpMenu += ".txt"; 
	 myEval->WriteL1Menu(tmpMenu);	

// Turn off making the rate vs threshold plots...already have them	
    makeThrPlots = 0;

// Get total rate with new thresholds
    totalRate = myL1Menu->RunL1Ana(tmpMenu,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);

// Save a copy of the first pass of the total rates (before any bandwidth scaling)
    TString fName = "L1Rates_";
	 fName += jobTag;
	 TString fNameSave = fName;
	 fName += ".txt";
	 fNameSave  += "_Default.txt";
	 TString cmd = "mv "; cmd += fName; cmd += " "; cmd += fNameSave;
	 gSystem->Exec(cmd.Data());

// Now iterate to fill the bandwidth if necessary
    int numIter = 0;
//	printf("Starting Interation %f %f %f \n",totalRate,targetTotalRate,fabs(totalRate-targetTotalRate));
	
   while(fabs(totalRate-targetTotalRate)>targetTolerance && numIter<maxIter) {

// Scale the bandwidths to get the total rate near the target
      double scaleFactor = targetTotalRate/totalRate;
      myEval->ScaleBandwidth(scaleFactor);

// Determine new thresholds for the new bandwidth
      myEval->DetermineThresholds();
	
// Write out new menu 
      myEval->WriteL1Menu(tmpMenu);

// Get total rate with new menu
      totalRate = myL1Menu->RunL1Ana(tmpMenu,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);
//	   printf("Total Rate Value %f after iteration %i \n",totalRate,numIter);

// Count interations
      numIter++;
   }

// Write out final L1 Menu
   TString finalMenu = "L1Menu_";
   finalMenu += jobTag;
   finalMenu += "_Final.txt";
   myEval->WriteL1Menu(finalMenu);	
		
   printf("Job %s ending...",jobTag.Data());	
   gSystem->Exec("date"); //ending time stamp

   
}
