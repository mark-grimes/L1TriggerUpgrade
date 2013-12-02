#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TText.h"
#include "TH2.h"
#include "TAxis.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>


class EvaluateL1Menu {
	public :

	EvaluateL1Menu(){}

	~EvaluateL1Menu() {}


   void  LoadL1Menu(TString fileName);
	void  LoadL1toOffline(TString fileName);
	void  LoadL1Rates(TString fileName, int nTAlgs=0);
	void  LoadL1RatesAndOfflineThr(TString fileName, int nTAlgs=0);
   void  LoadThresholdPlots(TString fileName);	
	void  WriteL1Menu(TString fileName);
	void  WriteL1Rates(TString fileName, bool suppressThr=true);
   void  DetermineThresholds();
	void  DetermineRates();
	void  FindThresholdsFromScaledRates();
	void  ScaleBandwidth(double scaleFactor);
	void  RateScaling(TString outHistFile, TString numFileName, TString demonFileName, TString muonFileName);
   double FindThreshold(TH1F* byThreshold, double Threshold);
   double FindRate(TString alg, double threshold);
	
	typedef struct {
	   float primTh ;
	   float secTh;
	   float triTh;
	   float quadTh;
	   float etaCut;
	   int minQual;
		float bandwidth;
		int   scalable;
		bool  locked;
	} trigPar;

   typedef struct {
	   float slope1;
		float offset1;
	   float slope2;
		float offset2;
	   float slope3;
		float offset3;
	   float slope4;
		float offset4;						
	} trig2Offline;

	typedef struct {
	   int bitNum;
	   float primTh ;
	   float secTh;
	   float triTh;
	   float quadTh;
	   float rate;
      float rateErr;
	   float cumuRate;
      float cumuRateErr;
	   float pureRate;
      float pureRateErr;	
		float primOffTh;
		float secOffTh;
		float triOffTh;
		float quadOffTh;			
	} trigRate;
	
	std::map<std::string, trigPar> trigParList;
	std::map<std::string, trigRate> trigRateList;
	std::map<std::string, trig2Offline> trig2OfflineList;
	std::map<std::string, int> Prescales;  //keep definitions used in L1Menu2012.
   std::map<std::string, int> BitMapping;

	private :

   TFile* histFile;
   TString RateHeader[2];
	TString RateFooter[5];
	TString orderList[50];
};

void EvaluateL1Menu::LoadL1Menu(TString fileName) {

// Open File
   printf("\n Reading L1 Menu File %s \n",fileName.Data());
   ifstream ifs( fileName );  
	
// Read through Menu
   while(ifs) {
	
	  TString algName;
	  ifs >> algName;
	  
	  ifs >> BitMapping[algName.Data()];	  
	  ifs >> Prescales[algName.Data()];
	  
	  ifs >> trigParList[algName.Data()].primTh;
	  ifs >> trigParList[algName.Data()].secTh;
	  ifs >> trigParList[algName.Data()].triTh;
	  ifs >> trigParList[algName.Data()].quadTh;
	  ifs >> trigParList[algName.Data()].etaCut;
	  ifs >> trigParList[algName.Data()].minQual;
	  ifs >> trigParList[algName.Data()].bandwidth;
	  ifs >> trigParList[algName.Data()].scalable;
	  ifs >> trigParList[algName.Data()].locked;	  
	  
//	  printf("Alg: %20s %2i %2i %6.2f %6.2f %6.2f %6.2f %6.2f %2i %6.2f %2i %2i\n",algName.Data(),BitMapping[algName.Data()],Prescales[algName.Data()],
//	         trigParList[algName.Data()].primTh,trigParList[algName.Data()].secTh,trigParList[algName.Data()].triTh,trigParList[algName.Data()].quadTh,
//				trigParList[algName.Data()].etaCut,trigParList[algName.Data()].minQual,trigParList[algName.Data()].bandwidth,trigParList[algName.Data()].scalable,trigParList[algName.Data()].locked);
	}

  return;
}

void EvaluateL1Menu::LoadL1toOffline(TString fileName) {

// Open File
   printf("\n Reading L1 to Offline Scalling File %s \n",fileName.Data());
   ifstream ifs( fileName );  
	
// Read through Menu
   while(ifs) {
	
	  TString algName;
	  ifs >> algName;	  
	  
	  ifs >> trig2OfflineList[algName.Data()].slope1;
	  ifs >> trig2OfflineList[algName.Data()].offset1;
	  ifs >> trig2OfflineList[algName.Data()].slope2;
	  ifs >> trig2OfflineList[algName.Data()].offset2;
	  ifs >> trig2OfflineList[algName.Data()].slope3;
	  ifs >> trig2OfflineList[algName.Data()].offset3;	  
	  ifs >> trig2OfflineList[algName.Data()].slope4;
	  ifs >> trig2OfflineList[algName.Data()].offset4;
	  
	}

  return;
}


void EvaluateL1Menu::LoadL1Rates(TString fileName, int nTAlgs) {

// Open File
   printf("\n Reading L1 Rates File %s \n",fileName.Data());
   ifstream ifs( fileName );
	char line[150];  


// Read the headers
   for(int i=0; i<2; i++) {
     ifs.getline(line,150);
	  RateHeader[i] = line;
	  printf("%s\n",RateHeader[i].Data());
   }
	
// Read through Menu
   int nRead = 0;
   while(ifs && nRead<nTAlgs) {
	 
	  int bitNumTmp;
	  ifs >> bitNumTmp;
	
	  TString algName;
	  ifs >> algName;
	  orderList[nRead] = algName;
	  
	  trigRateList[algName.Data()].bitNum = bitNumTmp;
	  
	  int tmp;
	  ifs >> tmp; //Prescale Not needed
	  
	  ifs >> trigRateList[algName.Data()].primTh;
	  ifs >> trigRateList[algName.Data()].secTh;
	  ifs >> trigRateList[algName.Data()].triTh;
	  ifs >> trigRateList[algName.Data()].quadTh;
	  ifs >> trigRateList[algName.Data()].rate;
	  ifs >> trigRateList[algName.Data()].rateErr;
	  ifs >> trigRateList[algName.Data()].cumuRate;
	  ifs >> trigRateList[algName.Data()].cumuRateErr;
	  ifs >> trigRateList[algName.Data()].pureRate;
	  ifs >> trigRateList[algName.Data()].pureRateErr;	  
	  
	  nRead++;
	  
	  printf("%2i %20s   %2i   %5.1f %5.1f %5.1f %5.1f  %8.2f  %5.2f   %8.2f  %5.2f   %8.2f  %5.2f \n",trigRateList[algName.Data()].bitNum,algName.Data(),Prescales[algName.Data()],
	         trigRateList[algName.Data()].primTh,trigRateList[algName.Data()].secTh,trigRateList[algName.Data()].triTh,trigRateList[algName.Data()].quadTh,
				trigRateList[algName.Data()].rate,trigRateList[algName.Data()].rateErr,trigRateList[algName.Data()].cumuRate,trigRateList[algName.Data()].cumuRateErr,
				trigRateList[algName.Data()].pureRate,trigRateList[algName.Data()].pureRateErr
				);
	}

// Read the headers
   ifs.getline(line,150); //get end of last line
   for(int i=0; i<5; i++) {
     ifs.getline(line,150);
	  RateFooter[i] = line;
	  printf("%s\n",RateFooter[i].Data());
   }



  return;
}


void EvaluateL1Menu::LoadL1RatesAndOfflineThr(TString fileName, int nTAlgs) {

// Open File
   printf("\n Reading L1 Rates File %s \n",fileName.Data());
   ifstream ifs( fileName );
	char line[150];  


// Read the headers
   for(int i=0; i<2; i++) {
     ifs.getline(line,150);
	  RateHeader[i] = line;
	  printf("%s\n",RateHeader[i].Data());
   }
	
// Read through Menu
   int nRead = 0;
   while(ifs && nRead<nTAlgs) {
	 
	  int bitNumTmp;
	  ifs >> bitNumTmp;
	
	  TString algName;
	  ifs >> algName;
	  orderList[nRead] = algName;
	  
	  trigRateList[algName.Data()].bitNum = bitNumTmp;
	  
	  int tmp;
	  ifs >> tmp; //Prescale Not needed
	  
	  ifs >> trigRateList[algName.Data()].primTh;
	  ifs >> trigRateList[algName.Data()].secTh;
	  ifs >> trigRateList[algName.Data()].triTh;
	  ifs >> trigRateList[algName.Data()].quadTh;

	  ifs >> trigRateList[algName.Data()].primOffTh;
	  ifs >> trigRateList[algName.Data()].secOffTh;
	  ifs >> trigRateList[algName.Data()].triOffTh;
	  ifs >> trigRateList[algName.Data()].quadOffTh;

	  ifs >> trigRateList[algName.Data()].rate;
	  ifs >> trigRateList[algName.Data()].rateErr;
	  ifs >> trigRateList[algName.Data()].cumuRate;
	  ifs >> trigRateList[algName.Data()].cumuRateErr;
	  ifs >> trigRateList[algName.Data()].pureRate;
	  ifs >> trigRateList[algName.Data()].pureRateErr;	  
	  
	  nRead++;
	  
	  printf("%2i %20s   %2i   %5.1f %5.1f %5.1f %5.1f  %5.1f %5.1f %5.1f %5.1f %8.2f  %5.2f   %8.2f  %5.2f   %8.2f  %5.2f \n",trigRateList[algName.Data()].bitNum,algName.Data(),Prescales[algName.Data()],
	         trigRateList[algName.Data()].primTh,trigRateList[algName.Data()].secTh,trigRateList[algName.Data()].triTh,trigRateList[algName.Data()].quadTh,
	         trigRateList[algName.Data()].primOffTh,trigRateList[algName.Data()].secOffTh,trigRateList[algName.Data()].triOffTh,trigRateList[algName.Data()].quadOffTh,
				trigRateList[algName.Data()].rate,trigRateList[algName.Data()].rateErr,trigRateList[algName.Data()].cumuRate,trigRateList[algName.Data()].cumuRateErr,
				trigRateList[algName.Data()].pureRate,trigRateList[algName.Data()].pureRateErr
				);
	}

// Read the headers
   ifs.getline(line,150); //get end of last line
   for(int i=0; i<5; i++) {
     ifs.getline(line,150);
	  RateFooter[i] = line;
	  printf("%s\n",RateFooter[i].Data());
   }



  return;
}



void EvaluateL1Menu::LoadThresholdPlots(TString fileName) {

  histFile = new TFile(fileName);
  printf("\n EvaluateL1Menu:: Loading rate vs Threshold plots from %s\n",fileName.Data());
  
  return;
}

void EvaluateL1Menu::WriteL1Menu(TString fileName) {

  FILE *outMenu = fopen(fileName,"w");
  
  printf("\n EvaluateL1Menu:: Writing New L1 Menu to %s\n \n",fileName.Data());
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {
      TString algName = it->first;
		
      if(algName.Contains("L1")) fprintf(outMenu,"%20s %2i %2i %6.2f %6.2f %6.2f %6.2f %6.2f %2i %6.2f %2i %2i\n",algName.Data(),BitMapping[algName.Data()],Prescales[algName.Data()],
	         trigParList[algName.Data()].primTh,trigParList[algName.Data()].secTh,trigParList[algName.Data()].triTh,trigParList[algName.Data()].quadTh,
				trigParList[algName.Data()].etaCut,trigParList[algName.Data()].minQual,trigParList[algName.Data()].bandwidth,trigParList[algName.Data()].scalable,trigParList[algName.Data()].locked); 
  }

  fclose(outMenu);
  
  return;
}


void EvaluateL1Menu::WriteL1Rates(TString fileName, bool suppressThr) {

  FILE *outRate = fopen(fileName,"w");

// Write the header (copied from unscaled)
   for(int i=0; i<2; i++)  fprintf(outRate,"%s\n",RateHeader[i].Data());
  
  printf("\n EvaluateL1Menu:: Writing Scaled L1 Rates to %s\n \n",fileName.Data());
//  for (std::map<std::string, trigRate>::iterator it=trigRateList.begin(); it != trigRateList.end(); it++) {
   float minThr = -999.;
	if(suppressThr) minThr = 1.;
   for(int i=0; i<23; i++) {
      TString algName = orderList[i];
	  	  
	  if(algName.Contains("L1")) { 
	         fprintf(outRate,"%2i %20s   %2i   %4.0f",trigRateList[algName.Data()].bitNum,algName.Data(),Prescales[algName.Data()],
	         trigRateList[algName.Data()].primTh);
				(trigRateList[algName.Data()].secTh<minThr)  ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].secTh);
				(trigRateList[algName.Data()].triTh<minThr)  ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].triTh);
				(trigRateList[algName.Data()].quadTh<minThr) ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].quadTh);
				
				(trigRateList[algName.Data()].primTh<minThr)  ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].primOffTh);
				(trigRateList[algName.Data()].secTh<minThr)  ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].secOffTh);
				(trigRateList[algName.Data()].triTh<minThr)  ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].triOffTh);
				(trigRateList[algName.Data()].quadTh<minThr) ? fprintf(outRate,"     "): fprintf(outRate," %4.0f",trigRateList[algName.Data()].quadOffTh);
								
				fprintf(outRate,"  %7.1f  %4.1f   %7.1f  %4.1f   %7.1f  %4.1f \n",
				trigRateList[algName.Data()].rate,trigRateList[algName.Data()].rateErr,trigRateList[algName.Data()].cumuRate,trigRateList[algName.Data()].cumuRateErr,
				trigRateList[algName.Data()].pureRate,trigRateList[algName.Data()].pureRateErr
				);
	  }	 		
	}

// Write the footer (copied from unscaled)
   for(int i=0; i<5; i++)  fprintf(outRate,"%s\n",RateFooter[i].Data());
  



  return;
}

void EvaluateL1Menu::DetermineThresholds() {

 
  printf("\n============= EvaluateL1Menu:: Finding new Thresholds ===============================================================\n");
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
        if(!trigParList[algName.Data()].locked  && Prescales[algName.Data()]>0) {
// Get name of threshold plot for this algorithm
          TString hname = "h_";
			 hname += algName(3,algName.Length());
			 hname += "_byThreshold";
          TH1F* tmp = (TH1F*)histFile->Get(hname)->Clone();		
			 
// Get dedicated threshold for this trigger
// If the trigger is prescaled inflate the bandwidth.  The bandwidth allowed is assumed to be after prescale and the rate vs threshold is assumed to be unprescaled.
          double threshold = FindThreshold(tmp,Prescales[algName.Data()]*trigParList[algName.Data()].bandwidth);			 
			 printf("Changing prim. threshold for %20s from %7.3f ---> %7.3f  to achieve target rate of %6.2f kHz\n",algName.Data(),trigParList[algName.Data()].primTh,threshold,trigParList[algName.Data()].bandwidth);
			 
			 float binWidth = tmp->GetBinWidth(1);
			 // use "floor" function to put non-primary thresholds on correct boundaries.
			 if(trigParList[algName.Data()].secTh>0.)  trigParList[algName.Data()].secTh  = binWidth*floor((threshold*(trigParList[algName.Data()].secTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigParList[algName.Data()].triTh>0.)  trigParList[algName.Data()].triTh  = binWidth*floor((threshold*(trigParList[algName.Data()].triTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigParList[algName.Data()].quadTh>0.) trigParList[algName.Data()].quadTh = binWidth*floor((threshold*(trigParList[algName.Data()].quadTh/trigParList[algName.Data()].primTh))/binWidth+0.5);
			 trigParList[algName.Data()].primTh = threshold;		 
		  } else if(trigParList[algName.Data()].locked) {
		    printf("Locked:: prim. threshold for %20s at   %7.3f \n",algName.Data(),trigParList[algName.Data()].primTh);
		  }  
		}  
  }
  printf("=====================================================================================================================\n");
  
  return;
}


void EvaluateL1Menu::DetermineRates() {

 
  double totalRate = 0.; 
  
  printf("\n============= EvaluateL1Menu:: Finding new Rates ===============================================================\n");
  for (std::map<std::string, trigRate>::iterator it=trigRateList.begin(); it != trigRateList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
        if(Prescales[algName.Data()]>0) {
		    printf("For Alg %s ",algName.Data());
// Get name of threshold plot for this algorithm
          TString hname = "h_";
			 hname += algName(3,algName.Length());
			 hname += "_byThreshold";
          TH1F* tmp = (TH1F*)histFile->Get(hname)->Clone();		

// Get Rate at this thresholds
          double threshold = trigRateList[algName.Data()].primTh;
			 printf(" rate for threshold %f ",threshold);
// Get the rate in this bin
          int thBin = tmp->FindBin(threshold);
			 double scaledRate = tmp->GetBinContent(thBin);
			 printf(" Old Rate %f  New Rate %f \n",trigRateList[algName.Data()].rate,scaledRate);

// Scale the pure rates assuming the fraction does not change			 
			 if(trigRateList[algName.Data()].rate>0.) trigRateList[algName.Data()].pureRate *= scaledRate/trigRateList[algName.Data()].rate;
//			 trigRateList[algName.Data()].pureRateErr *=  ?????
			 
// Now store the scaled rate
			 trigRateList[algName.Data()].rate = scaledRate;
			 trigRateList[algName.Data()].rateErr = tmp->GetBinError(thBin);
			 totalRate += scaledRate;
			 
		  } else if(trigParList[algName.Data()].locked) {
		    printf("Disabled for %20s at   %7.3f \n",algName.Data());
		  }  
		}  
  }
  printf("Sum of Individual Rates %f\n",totalRate);
  printf("=====================================================================================================================\n");
  
  return;
}

void EvaluateL1Menu::FindThresholdsFromScaledRates(){

 
  double totalRate = 0.;
  
  printf("\n============= EvaluateL1Menu:: Finding new Thresholds from Scaled Rates  ============================================\n");
  for (std::map<std::string, trigRate>::iterator it=trigRateList.begin(); it != trigRateList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
        if(Prescales[algName.Data()]>0) {
// Get name of threshold plot for this algorithm
          TString hname = "h_";
			 hname += algName(3,algName.Length());
			 hname += "_byThreshold";
          TH1F* tmp = (TH1F*)histFile->Get(hname)->Clone();		
			 
// Get dedicated threshold for this trigger
// If the trigger is prescaled inflate the bandwidth.  The bandwidth allowed is assumed to be after prescale and the rate vs threshold is assumed to be unprescaled.
          double threshold = FindThreshold(tmp,Prescales[algName.Data()]*trigRateList[algName.Data()].rate);			 
			 printf("Changing prim. threshold for %20s from %7.3f ---> %7.3f  to achieve target rate of %6.2f kHz\n",algName.Data(),trigRateList[algName.Data()].primTh,threshold,trigRateList[algName.Data()].rate);
			 
			 float binWidth = tmp->GetBinWidth(1);
			 // use "floor" function to put non-primary thresholds on correct boundaries.
			 if(trigRateList[algName.Data()].secTh>0.)  trigRateList[algName.Data()].secTh  = binWidth*floor((threshold*(trigRateList[algName.Data()].secTh/trigRateList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigRateList[algName.Data()].triTh>0.)  trigRateList[algName.Data()].triTh  = binWidth*floor((threshold*(trigRateList[algName.Data()].triTh/trigRateList[algName.Data()].primTh))/binWidth+0.5);
			 if(trigRateList[algName.Data()].quadTh>0.) trigRateList[algName.Data()].quadTh = binWidth*floor((threshold*(trigRateList[algName.Data()].quadTh/trigRateList[algName.Data()].primTh))/binWidth+0.5);
			 trigRateList[algName.Data()].primTh = threshold;

// Get the rate in this bin
          int thBin = tmp->FindBin(threshold);
			 double scaledRate = tmp->GetBinContent(thBin);

// Scale the pure rates assuming the fraction does not change			 
			 if(trigRateList[algName.Data()].rate>0.) trigRateList[algName.Data()].pureRate *= scaledRate/trigRateList[algName.Data()].rate;
//			 trigRateList[algName.Data()].pureRateErr *=  ?????
			 
// Now store the scaled rate
			 trigRateList[algName.Data()].rate = scaledRate;
			 trigRateList[algName.Data()].rateErr = tmp->GetBinError(thBin);
			 totalRate += scaledRate;

// Calculate the offline Threshold
          trigRateList[algName.Data()].primOffTh = trigRateList[algName.Data()].primTh*trig2OfflineList[algName.Data()].slope1 + trig2OfflineList[algName.Data()].offset1;
			 //printf("Converting primary threshold to Offline.  Slope %f, Offset %f   L1 Thres  %f  Off Thres %f\n",trig2OfflineList[algName.Data()].slope1,trig2OfflineList[algName.Data()].offset1,trigRateList[algName.Data()].primTh,trigRateList[algName.Data()].primOffTh);
			 trigRateList[algName.Data()].secOffTh  = trigRateList[algName.Data()].secTh*trig2OfflineList[algName.Data()].slope2  + trig2OfflineList[algName.Data()].offset2;
			 trigRateList[algName.Data()].triOffTh  = trigRateList[algName.Data()].triTh*trig2OfflineList[algName.Data()].slope3  + trig2OfflineList[algName.Data()].offset3;
			 trigRateList[algName.Data()].quadOffTh = trigRateList[algName.Data()].quadTh*trig2OfflineList[algName.Data()].slope4 + trig2OfflineList[algName.Data()].offset4;			 
			 
		  }  
		}  
  }
  printf("Total Scaled Rate %f\n",totalRate);
  printf("=====================================================================================================================\n");
  
  return;
}
void EvaluateL1Menu::ScaleBandwidth(double scaleFactor) {

 
  printf("\n============== EvaluateL1Menu:: Scaling Bandwidths ===================================\n");
  for (std::map<std::string, trigPar>::iterator it=trigParList.begin(); it != trigParList.end(); it++) {

      TString algName = it->first;		
      if(algName.Contains("L1")) {
		   if(trigParList[algName.Data()].scalable == 1 && !trigParList[algName.Data()].locked && Prescales[algName.Data()]>0) {
		     printf("Changing allocated bandwidth for %20s from %7.3f ---> %7.3f \n",algName.Data(),trigParList[algName.Data()].bandwidth,scaleFactor*trigParList[algName.Data()].bandwidth);
           trigParList[algName.Data()].bandwidth *= scaleFactor; 			  
		   } else if( (trigParList[algName.Data()].scalable == 0 || trigParList[algName.Data()].locked) && Prescales[algName.Data()]>0){
			  printf("Algorithm Rates are locked for   %20s at %7.3f \n",algName.Data(),trigParList[algName.Data()].bandwidth);
			}
		}  
  }
  printf("========================================================================================\n");
  
  return;
}


double EvaluateL1Menu::FindThreshold(TH1F* byThreshold, double targetRate) {

	double threshold = 999.;
   int i=1;
	bool fnd = false;
   while(i<byThreshold->GetNbinsX() && !fnd) {
	   if(byThreshold->GetBinContent(i)<targetRate) {
          threshold = byThreshold->GetBinCenter(i);
			 fnd = true;
			 
// Check whether previous bin was closer to target, if so, use it.
          if((byThreshold->GetBinContent(i-1)-targetRate)<
			    (targetRate-byThreshold->GetBinContent(i))) {
				    threshold = byThreshold->GetBinCenter(i-1);
  			 }		 			  			 
		} else {
		    i++;
		}		   
	}
	if(!fnd) printf("WARNING: did not find an acceptable threshold.\n");

  return threshold;
}


void EvaluateL1Menu::RateScaling(TString outFileName, TString numFileName, TString denomFileName, TString muonFileName){

	TFile* numFile = new TFile(numFileName);
   TFile* denomFile = new TFile(denomFileName);
	TFile* muonFile  = new TFile(muonFileName);

// Get Scale Factors for Muons
    TH1F* muonScalePt   = (TH1F*)muonFile->Get("muonPtScale")->Clone();		   
    TH1F* muonScaleIsol = (TH1F*)muonFile->Get("muonIsoScale")->Clone();		   


// Open the output File
   TFile *outFile = new TFile(outFileName,"RECREATE");
	outFile->cd();
	TH1::SetDefaultSumw2();


// Loop over each histogram to process
   for (std::map<std::string, trigRate>::iterator it=trigRateList.begin(); it != trigRateList.end(); it++) {
		TString trigAlg = it->first;
		trigAlg.Remove(0,3);


//      TFile* defFile = new TFile(defFileName);
		TString hname = "h_"; hname += trigAlg; hname += "_byThreshold";
//		TH1F* defHist = (TH1F*)defFile->Get(hname)->Clone();
		TH1F* defHist = (TH1F*)histFile->Get(hname)->Clone("default");
//      defHist->DrawCopy();

      
		TH1F*  numHist = (TH1F*)numFile->Get(hname)->Clone("num");
//		numHist->SetLineColor(kGreen);
//      numHist->DrawCopy("HISTSAME");


		TH1F*  denomHist = (TH1F*)denomFile->Get(hname)->Clone("denom");
//		denomHist->SetLineColor(kRed);
//		denomHist->DrawCopy("HISTSAME");
		
		defHist->Multiply(numHist);
		defHist->Divide(denomHist);
//		defHist->SetMarkerStyle(1);
//		defHist->Draw("EHISTSAME");


// Scaling for Muon System Upgrades  (Problem: Can't easily scale the EG_Mu trigger because x-axis is EG Et not muon Pt...need to fix)
      if(trigAlg.Contains("Mu") && !trigAlg.Contains("EG_Mu") ) {
         
//  Apply Scale Factor for Muon Pt Assignment
         printf("Applying the Scaling for Muon Pt Assignment for Trigger %s \n",trigAlg.Data());
			defHist->Multiply(muonScalePt);
//	      defHist->Scale(0.69);	   
		}

// Scaling for Muon Isolation  
      if(trigAlg.Contains("IsoMu") || trigAlg.Contains("isoMu")) {
         
//  Apply Scale Factor for Muon Pt Assignment
         printf("Applying the Scaling for Muon Isolation for Trigger %s \n",trigAlg.Data());
			defHist->Multiply(muonScaleIsol);
		   
		}


		
		outFile->cd();
      TH1F* tmp = (TH1F*)outFile->Get(hname);
 		tmp=(TH1F*)defHist->Clone(hname);
//      tmp->Draw();
	}

// Write out the results and close
   outFile->Write();
	outFile->Close();
   
	return;
}
double EvaluateL1Menu::FindRate(TString alg, double threshold) {

// Get name of threshold plot for this algorithm
          TString hname = "h_";
			 hname += alg;
			 hname += "_byThreshold";
          TH1F* tmp = (TH1F*)histFile->Get(hname)->Clone();		


// Get the rate and error for this threshold 
          int thBin      = tmp->FindBin(threshold);
			 double Rate    = tmp->GetBinContent(thBin);
			 double errRate = tmp->GetBinError(thBin);
			 printf("For %s Rate \t \t %5.1f  +/-  %5.1f \n",alg.Data(),Rate,errRate);
  
   return Rate;
}
