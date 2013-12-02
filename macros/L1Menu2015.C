#include "L1UpgradeNtuple.h"
#include "L1AnalysisDataFormat.h"
#include "hist.C"
#include "Style.C"

#include "TLegend.h"
#include "TMath.h"
#include "TText.h"
#include "TH2.h"
#include "TAxis.h"
#include "TString.h"
#include "TRandom2.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <set>

// -- Huge prescale value for seeds "for lower PU"
#define INFTY 10000



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Notes:
  
    -> This needs to be run within the UserCode/L1TriggerDPG package
	 
	 -> In also requires:
	           L1AnalaysisDataFormat.h
				  getLumi_out_pixelCorrLumi_*PU_stdCorr.txt (for HPF Data)
				  
	 -> General running format (from UserCode/L1TriggerDPG/macro/ area)
	 
	    linux>  root initL1Analysis.C      (loading libraries etc) 
		 root>   .L L1Menu2015.C++          (compiles the script) 			  
       root>   RunL1_HFW(Bool_t calcThreshold=false,Bool_t selectDataInput=0, Int_t pMenu2015 = 0, Int_t usedL1MenuThr=0,Int_t whichDataSetToUse=1,Int_ whichFileToUse=0, Float_t targetlumi=200, Int_t pNevts = -1)
 
       Definition of input quantities:
                  calcThreshold = flag for whether rate vs threshold plots are made (uses a lot more CPU)
						selectDataInput = Selects which information to use for the event (see fillDataStructure method for more)
						pMenu2015     = Selects different trigger menus (set in EvalMenu)
						usedL1MenuThr = Selects the thresholds to be used for the algorithms (set in MyInit)
						whichDataSetToUse = Specifies input data set                                                NOTE: This is ugly and should be improved.
						whichFileToUse = Also a file specifier for selecting different versions of the same sample  NOTE:  ditto
						targetlumi    = units of E32
						pNevts        = number of events to run over (-1 ==> run over all events in file)
						
   There are a lot of details that are not yet summarized here....
	
					

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
//HFW
//Cross
TH1F *h_SingleMu_ETM_byThreshold; 
TH2F *h2_SingleMu_ETM_byThreshold;
TH1F *h_SingleMu_CJet_byThreshold;
TH2F *h2_SingleMu_CJet_byThreshold;
TH1F *h_SingleIsoMu_ETM_byThreshold; 
TH2F *h2_SingleIsoMu_ETM_byThreshold;
TH1F *h_SingleIsoMu_CJet_byThreshold;
TH2F *h2_SingleIsoMu_CJet_byThreshold;
TH1F *h_SingleMu_HTM_byThreshold;
TH1F *h_SingleIsoMu_HTM_byThreshold;  
TH2F *h2_SingleMu_HTM_byThreshold;
TH2F *h2_SingleIsoMu_HTM_byThreshold;  


TH1F *h_SingleEG_ETM_byThreshold; 
TH2F *h2_SingleEG_ETM_byThreshold;
TH1F *h_SingleEG_CJet_byThreshold;
TH2F *h2_SingleEG_CJet_byThreshold;
TH1F *h_SingleIsoEG_ETM_byThreshold; 
TH2F *h2_SingleIsoEG_ETM_byThreshold;
TH1F *h_SingleIsoEG_CJet_byThreshold;
TH2F *h2_SingleIsoEG_CJet_byThreshold;
TH1F *h_SingleEG_HTM_byThreshold;
TH1F *h_SingleIsoEG_HTM_byThreshold;  
TH2F *h2_SingleEG_HTM_byThreshold;
TH2F *h2_SingleIsoEG_HTM_byThreshold;  


TH1F *h_Mu_EG_byThreshold;
TH2F *h2_Mu_EG_byThreshold;
TH1F *h_isoMu_EG_byThreshold;
TH2F *h2_isoMu_EG_byThreshold;
TH1F *h_isoMu_isoEG_byThreshold;
TH2F *h2_isoMu_isoEG_byThreshold;

TH1F *h_EG_Mu_byThreshold;
TH1F *h_isoEG_Mu_byThreshold;
TH2F *h2_isoEG_Mu_byThreshold;
TH1F *h_isoEG_isoMu_byThreshold;
TH2F *h2_isoEG_isoMu_byThreshold;

TH1F *h_EG_Tau_byThreshold;
TH2F *h2_EG_Tau_byThreshold;
TH1F *h_isoEG_Tau_byThreshold;
TH2F *h2_isoEG_Tau_byThreshold;
TH1F *h_isoEG_isoTau_byThreshold;
TH2F *h2_isoEG_isoTau_byThreshold;

TH1F *h_Mu_Tau_byThreshold;
TH2F *h2_Mu_Tau_byThreshold;
TH1F *h_isoMu_Tau_byThreshold;
TH2F *h2_isoMu_Tau_byThreshold;
TH1F *h_isoMu_isoTau_byThreshold;
TH2F *h2_isoMu_isoTau_byThreshold;


//Jets
TH1F *h_SingleJet_byThreshold;
TH1F *h_SingleJetC_byThreshold;
TH1F *h_DoubleJet_byThreshold;
TH2F *h2_DoubleJet_byThreshold;
TH1F *h_DoubleFwdJet_byThreshold;
TH2F *h2_DoubleFwdJet_byThreshold;
TH1F *h_QuadJetC_byThreshold;
TH2F *h2A_QuadJetCentral_byThreshold;
TH2F *h2B_QuadJetCentral_byThreshold;
TH1F *h_SixJet_byThreshold;
TH2F *h2A_SixJet_byThreshold;
TH1F *h_SingleCJet_ETM_byThreshold;
TH2F *h2_SingleCJet_ETM_byThreshold;
TH1F *h_DoubleCJet_ETM_byThreshold;
TH2F *h2_DoubleCJet_ETM_byThreshold;

// Taus
TH1F *h_SingleTau_byThreshold;
TH1F *h_SingleIsoTau_byThreshold;
TH1F *h_DoubleTau_byThreshold;
TH2F *h2_DoubleTau_byThreshold;
TH1F *h_isoTau_Tau_byThreshold;
TH2F *h2_isoTau_Tau_byThreshold;
TH1F *h_DoubleIsoTau_byThreshold;
TH2F *h2_DoubleIsoTau_byThreshold;
TH1F *h_SingleTau_ETM_byThreshold;
TH2F *h2_SingleTau_ETM_byThreshold;
TH1F *h_SingleIsoTau_ETM_byThreshold;
TH1F *h_SingleTau_CJet_byThreshold;
TH2F *h2_SingleTau_CJet_byThreshold;
TH1F *h_SingleIsoTau_CJet_byThreshold;
TH1F *h_SingleTau_TwoFJet_byThreshold;
TH2F *h2_SingleTau_TwoFJet_byThreshold;
TH1F *h_SingleTau_HTM_byThreshold;
TH1F *h_SingleIsoTau_HTM_byThreshold;  



//Sums
TH1F *h_HTT_byThreshold;
TH1F *h_ETM_byThreshold;
TH1F *h_HTM_byThreshold;
TH1F *h_HTT_ETM_byThreshold;
TH2F *h2_HTT_ETM_byThreshold;

//EGamma
TH1F *h_SingleEG_byThreshold;
TH1F *h_SingleIsoEG_byThreshold;
TH1F *h_DoubleEG_byThreshold;
TH2F *h2_DoubleEG_byThreshold;
TH1F *h_isoEG_EG_byThreshold;
TH1F *h2_isoEG_EG_byThreshold;
TH1F *h_EG_isoEG_byThreshold;
TH1F *h_DoubleIsoEG_byThreshold;
TH2F *h2_DoubleIsoEG_byThreshold;


//Muons
TH1F *h_SingleMu_byThreshold;
TH1F *h_SingleIsoMu_byThreshold;
TH1F *h_DoubleMu_byThreshold;
TH2F *h2_DoubleMu_byThreshold;
TH1F *h_isoMu_Mu_byThreshold;
TH2F *h2_isoMu_Mu_byThreshold;
TH1F *h_DoubleIsoMu_byThreshold;
TH2F *h2_DoubleIsoMu_byThreshold;



// Plots for Trigger Quantities
TH1F *h_Mu_Nmu,      *h_Mu_Et,     *h_Mu_Eta,     *h_Mu_Phi;
TH1F *h_isoMu_Nmu,      *h_isoMu_Et,     *h_isoMu_Eta,     *h_isoMu_Phi;
TH1F *h_isoEG_Nele,  *h_isoEG_Et,  *h_isoEG_Eta,  *h_isoEG_Phi;
TH1F *h_nIsoEG_Nele, *h_nIsoEG_Et, *h_nIsoEG_Eta, *h_nIsoEG_Phi;
TH1F *h_CJet_Njet, *h_CJet_Et, *h_CJet_Eta, *h_CJet_Phi;
TH1F *h_FJet_Njet, *h_FJet_Et, *h_FJet_Eta, *h_FJet_Phi;
TH1F *h_TJet_Njet, *h_TJet_Et, *h_TJet_Eta, *h_TJet_Phi;
TH1F *h_isoTJet_Njet, *h_isoTJet_Et, *h_isoTJet_Eta, *h_isoTJet_Phi;
TH1F *h_Sum_ETT,   *h_Sum_ETM, *h_Sum_PhiETM;
TH1F *h_Sum_HTT,   *h_Sum_HTM, *h_Sum_PhiHTM;
TH1F *h_puWeight,  *h_Sim_meanInt, *h_lumiSec, *h_bx;      


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Int_t NPAGS = 6;
TH2F *cor_PAGS;
TH1F *h_PAGS_pure;
TH1F *h_PAGS_shared;

const Int_t N128 = 128;			// could be > 128 for "test seeds"
Int_t kOFFSET = 0;
Bool_t TheTriggerBits[N128];	// contains the emulated triggers for each event
TH1F *h_All, *h_Trig;		// one bin for each trigger. Fill bin i if event fires trigger i.
TH1F *h_Pure;		// one bin for each trigger. Fill bin i if event fires trigger i and NO OTHER TRIGGER.
TH2F *h_Corr;
TH1F *h_Cumm;

Int_t Menu2015 = 0;

// Methods to scale L1 jets for new HCAL LUTs and estimate the rate changes 

// correction by 5% overall (from HCAL January 2012)
Double_t CorrectedL1JetPtByFactor(Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr = JetPt;

	if (theL1JetCorrection) {
		JetPtcorr = JetPt*1.05;
	}
	return JetPtcorr;
}

// correction by 8% for forward jets (from HCAL January 2012)
Double_t CorrectedL1FwdJetPtByFactor(Bool_t isFwdJet, Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr = JetPt;

	if (theL1JetCorrection) {
		if (isFwdJet) { JetPtcorr = JetPt*1.08; }
	}
	return JetPtcorr;
}

// correction for HF bins (from HCAL January 2012)
Size_t   JetHFiEtabins   = 13;
Int_t    JetHFiEtabin[]  = {29,30,31,32,33,34,35,36,37,38,39,40,41};
Double_t JetHFiEtacorr[] = {0.982,0.962,0.952, 0.943,0.947,0.939, 0.938,0.935,0.934, 0.935,0.942,0.923,0.914};

Double_t CorrectedL1JetPtByHFtowers(Double_t JetiEta,Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr   = JetPt;

	if (theL1JetCorrection) {
		Int_t    iJetiEtabin = 0;
		for (iJetiEtabin=0; iJetiEtabin<JetHFiEtabins; iJetiEtabin++) {
			if (JetHFiEtabin[iJetiEtabin]==JetiEta) {
				JetPtcorr = JetPt * (1+(1-JetHFiEtacorr[iJetiEtabin]));
			}
		}
	}
	return JetPtcorr;
}

// correction for RCT->GCT bins (from HCAL January 2012)
// HF from 29-41, first 3 HF trigger towers 3 iEtas, last highest eta HF trigger tower 4 iEtas; each trigger tower is 0.5 eta, RCT iEta from 0->21 (left->right)
Double_t JetRCTHFiEtacorr[]  = {0.965,0.943,0.936,0.929}; // from HF iEta=29 to 41 (smaller->higher HF iEta)

Double_t CorrectedL1JetPtByGCTregions(Double_t JetiEta,Double_t JetPt, Bool_t theL1JetCorrection=false) {

	Double_t JetPtcorr   = JetPt;

	if (theL1JetCorrection) {

		if ((JetiEta>=7 && JetiEta<=14)) {
			JetPtcorr = JetPt * 1.05;
		}

		if ((JetiEta>=4 && JetiEta<=6) || (JetiEta>=15 && JetiEta<=17)) {
			JetPtcorr = JetPt * 0.95;
		}

		if (JetiEta==0 || JetiEta==21) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[3]));
		}
		else if (JetiEta==1 || JetiEta==20) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[2]));
		}
		else if (JetiEta==2 || JetiEta==19) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[1]));
		}
		else if (JetiEta==3 || JetiEta==18) {
			JetPtcorr = JetPt * (1+(1-JetRCTHFiEtacorr[0]));
		}
	}

	return JetPtcorr;
}

// methods for the correlation conditions

size_t PHIBINS = 18;
Double_t PHIBIN[] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350};

size_t ETABINS = 23;
Double_t ETABIN[] = {-5.,-4.5,-4.,-3.5,
	-3.,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,
	0,
	0.348,0.696,1.044,1.392,1.74,2.172,3.,
	3.5,4.,4.5,5.};

size_t ETAMUBINS = 65;
Double_t ETAMU[] = { -2.45,-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45 };

Int_t etaMuIdx(Double_t eta) {
	size_t etaIdx = 0.;
	for (size_t idx=0; idx<ETAMUBINS; idx++) {
		if (eta>=ETAMU[idx] and eta<ETAMU[idx+1])
			etaIdx = idx;
	}
	return int(etaIdx);
}

Int_t etaINjetCoord(Double_t eta){
	size_t etaIdx = 0.;
	for (size_t idx=0; idx<ETABINS; idx++) {
		if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
			etaIdx = idx;
	}
	return int(etaIdx);
}

Double_t degree(Double_t radian) {
	if (radian<0)
		return 360.+(radian/TMath::Pi()*180.);
	else
		return radian/TMath::Pi()*180.;
}

Int_t phiINjetCoord(Double_t phi) {
	size_t phiIdx = 0;
	Double_t phidegree = degree(phi);
	for (size_t idx=0; idx<PHIBINS; idx++) {
		if (phidegree>=PHIBIN[idx] and phidegree<PHIBIN[idx+1])
			phiIdx = idx;
		else if (phidegree>=PHIBIN[PHIBINS-1] || phidegree<=PHIBIN[0])
			phiIdx = idx;
	}
	phiIdx = phiIdx + 1;
	if (phiIdx == 18)  phiIdx = 0;
	return int(phiIdx);
}

Bool_t correlateInPhi(Int_t jetphi, Int_t muphi, Int_t delta=1) {

	Bool_t correlationINphi = fabs(muphi-jetphi)<fabs(2 +delta-1) || fabs(muphi-jetphi)>fabs(PHIBINS-2 - (delta-1) );
	return correlationINphi;

}

Bool_t correlateInEta(Int_t mueta, Int_t jeteta, Int_t delta=1) {
	Bool_t correlationINeta = fabs(mueta-jeteta)<2 + delta-1;
	return correlationINeta;
}

// set the errors properly
void CorrectScale(TH1F* h, Float_t scal) {

	Int_t nbins = h -> GetNbinsX();

	for (Int_t i=1; i<= nbins; i++)  {
		Float_t val = h -> GetBinContent(i);
		Float_t er = sqrt(val);
		val = val * scal;
		er = er * scal;
		h -> SetBinContent(i,val);
		h -> SetBinError(i,er);
	}
}

class L1Menu2015 : public L1UpgradeNtuple {
	public :


	L1Menu2015(Int_t aL1Menu,Float_t aTargetLumi, Float_t aNumberOfUserdLumiSections, Float_t aLumiForThisSetOfLumiSections, std::string aL1NtupleFileName,Float_t aAveragePU, Float_t aZeroBiasPrescale,Bool_t aL1JetCorrection) : 
     	theL1Menu(aL1Menu), 
		theTargetLumi(aTargetLumi), 
		theNumberOfUserdLumiSections(aNumberOfUserdLumiSections),
		theLumiForThisSetOfLumiSections(aLumiForThisSetOfLumiSections),
		theL1NtupleFileName(aL1NtupleFileName),
		theAveragePU(aAveragePU),
		theZeroBiasPrescale(aZeroBiasPrescale),
		theL1JetCorrection(aL1JetCorrection),
		usePUWeight(false)
		{}


//	L1Menu2015()
//		{ printf("Here\n"); }

	~L1Menu2015() {}


     Float_t RunL1Ana(
      	TString L1MenuFileName, 					//File name containing L1 Menu Algorithms and thresholds. (See InitL1Menu for format)
      	TString lsFileName,  						//LumiSection Luminosity (used for data)
      	TString jobTag = "Test",					//Job tag for storing output
      	Int_t   procNevts=-1,						//Number of events to run over.  If -1, run over all events in the file
      	Int_t   makeThresholdPlots=0,  			//Flag for whether to calculate the rate vs threshold plots (0=skip all; 1=do 1-D plots; 2=do 1-d and 2-d plots)
      	Int_t  selectInputData=0 					//Selects what input to use for event
                );

	Int_t theL1Menu;

	// The luminosity for which we want the rates, in units 1e32 (this sets also the right prescales for the 2 menus and some descriptions correctly).
	// 70. for 7e33, 50 for 5e33, etc. Use 70.001 for the "emergency columns" to the corresponding target luminosity.
	// For the moment we have pre-scales for 5e33,6e33,7e33 plus emergency pre-scales for 5e33 and 7e33
	Float_t theTargetLumi;

	// the setting below are/will be specific for each L1Ntuple file used
	Float_t theNumberOfUserdLumiSections;
	Float_t theLumiForThisSetOfLumiSections;
	std::string theL1NtupleFileName;
	Float_t theAveragePU;
	Float_t theZeroBiasPrescale;
	Bool_t theL1JetCorrection;


	void InitL1Menu(TString menuFile);
	void FilL1Bits();
	void UsePUWeight(Bool_t flg) {usePUWeight = flg;};
	void UseUpgradeMuons(Bool_t flg) {useUpgradeMuons = flg;};
	void dumpEvent(int dumpCode=0xFFFF);
	int numDumpEvt;
	double calculateHTT();
	double calculateHTM();
	void SetDumpEvents(int num) {numDumpEvt = num;}
	
	L1Analysis::L1AnalysisDataFormat myEvt_;

	std::map<std::string, int> Counts;
	std::map<std::string, int> Prescales;
	std::map<std::string, bool> Biased;

        std::map<std::string, int> BitMapping;
	
	typedef struct {
	   float primTh ;
	   float secTh;
	   float triTh;
	   float quadTh;
	   float etaCut;
	   int minQual;
		float bandwidth;
		int  scalable;
		bool locked;
	} trigPar;
	
	std::map<std::string, trigPar> trigParList;

	std::map<std::string, float> WeightsPAGs;


	void InsertInMenu(std::string L1name, Bool_t value);

	Int_t L1BitNumber(std::string l1name);

	Bool_t EvalMenu(double lumiWeight);
	TH1F   *puWeightsHist;
	Bool_t usePUWeight;
	
	TRandom2 *randGen;
	Bool_t useUpgradeMuons;
	Bool_t UpgradeMuon(Float_t pt,Float_t eta, Float_t ptcut, Float_t etaCut);
	Bool_t UpgradeIsolatedMuon(Float_t pt,Float_t eta, Float_t ptcut, Float_t etaCut);
	double puWeight;
	void   calculatePUWeighting(TString mcHistFile, TString dataHistFile, bool dumpWgts=false);
	double getPUWeight(float mcMeanInt);
	void EvalThresh(Int_t calcThreshold, double lumiWeight);

// -- Cross
	Bool_t Mu_EG(Float_t mucut, Float_t EGcut, Int_t minMuQual = 4 );
	Bool_t Mu_isoEG(Float_t mucut, Float_t EGcut, Int_t minMuQual = 4 );	
	Bool_t EG_Mu(Float_t EGcut, Float_t mucut, Int_t minMuQual = 4 );
	Bool_t isoEG_Mu(Float_t EGcut, Float_t mucut, Int_t minMuQual = 4 );
	Bool_t MuOpen_EG(Float_t mucut, Float_t EGcut );
	Bool_t Mu_JetCentral(Float_t mucut, Float_t jetcut );
	Bool_t Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut );
	Bool_t Mu_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut );
	Bool_t Muer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut = 2.1, Int_t minMuQual=4 );
	Bool_t IsoMuer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut = 2.1, Int_t minMuQual=4 );	
	Bool_t Muer_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut );
	Bool_t Mu_HTT(Float_t mucut, Float_t HTcut );
	Bool_t Muer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut = 2.1, Int_t minMuQual=4 );
	Bool_t Muer_HTM(Float_t mucut, Float_t ETMcut, Float_t etacut = 2.1, Int_t minMuQual=4 );	
	Bool_t IsoMuer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut = 2.1, Int_t minMuQual=4 );
	Bool_t IsoMuer_HTM(Float_t mucut, Float_t ETMcut, Float_t etacut = 2.1, Int_t minMuQual=4 );	
	Bool_t EG_FwdJet(Float_t EGcut, Float_t FWcut ) ;
	Bool_t EG_JetCentral(Float_t EGcut, Float_t jetcut );
	Bool_t IsoEG_JetCentral(Float_t EGcut, Float_t jetcut, Float_t etaCut );	
	Bool_t EG_HT(Float_t EGcut, Float_t HTcut );
	Bool_t EG_DoubleJetCentral(Float_t EGcut, Float_t jetcut );
	Bool_t EG_ETM(Float_t EGcut, Float_t ETMcut);
	Bool_t EG_HTM(Float_t EGcut, Float_t ETMcut);	
	Bool_t IsoEG_ETM(Float_t EGcut, Float_t ETMcut, Float_t etaCut  );
	Bool_t IsoEG_HTM(Float_t EGcut, Float_t ETMcut, Float_t etaCut  );  
	Bool_t DoubleEG_HT(Float_t EGcut, Float_t HTcut );
	Bool_t EGEta2p1_JetCentral(Float_t EGcut, Float_t jetcut);		// delta
	Bool_t EGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut);          // delta
	Bool_t IsoEGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut);          // delta
	Bool_t EGEta2p1_DoubleJetCentral(Float_t EGcut, Float_t jetcut);	// delta
	Bool_t EGEta2p1_DoubleJetCentral_TripleJetCentral(Float_t EGcut, Float_t jetcut2, Float_t jetcut3);   

   Bool_t EG_Tau(Float_t EGcut, Float_t taucut, Float_t etaCut);
   Bool_t isoEG_Tau(Float_t EGcut, Float_t taucut, Float_t etaCut);
   Bool_t isoEG_isoTau(Float_t EGcut, Float_t taucut, Float_t etaCut);	
	Bool_t Mu_Tau(Float_t Mucut, Float_t taucut, Float_t etaCut, Int_t minMuQual);
	Bool_t Mu_isoTau(Float_t Mucut, Float_t taucut, Float_t etaCut, Int_t minMuQual);	
	Bool_t Tau_JetCentral(Float_t taucut, Float_t jetcut, Float_t etaCut);
	Bool_t IsoTau_JetCentral(Float_t taucut, Float_t jetcut, Float_t etaCut);	
	Bool_t Tau_TwoJetForward(Float_t taucut, Float_t jetcut1, Float_t jetcut2);	
	Bool_t Tau_ETM(Float_t taucut, Float_t ETMcut, Float_t etaCut);
	Bool_t IsoTau_ETM(Float_t taucut, Float_t ETMcut, Float_t etaCut);	
	Bool_t Tau_HTM(Float_t taucut, Float_t ETMcut, Float_t etaCut);
	Bool_t IsoTau_HTM(Float_t taucut, Float_t ETMcut, Float_t etaCut);	

	
	Bool_t HTT_HTM(Float_t HTTcut, Float_t HTMcut);
	Bool_t HTT_ETM(Float_t HTTcut, Float_t ETMcut);
	Bool_t JetCentral_ETM(Float_t jetcut, Float_t ETMcut, Float_t etaCut);
	Bool_t DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut, Float_t etaCut);
	Bool_t DoubleMu_EG(Float_t mucut, Float_t EGcut );
	Bool_t Mu_DoubleEG(Float_t mucut, Float_t EGcut);

	Bool_t Muer_TripleJetCentral(Float_t mucut, Float_t jet1, Float_t jet2, Float_t jet3);
	Bool_t Mia(Float_t mucut, Float_t jet1, Float_t jet2);	// delta
	Bool_t Mu_JetCentral_delta(Float_t mucut, Float_t ptcut);	// delta
	Bool_t Mu_JetCentral_deltaOut(Float_t mucut, Float_t ptcut); // delta


// -- Jets 
	Bool_t SingleJet(Float_t cut);	
	Bool_t SingleTauJet(Float_t cut, Float_t etaCut=4.5);
   Bool_t SingleIsoTauJet(Float_t cut, Float_t etaCut=4.5 ); 
	Bool_t SingleJetCentral(Float_t cut);
	Bool_t DoubleJetCentral(Float_t cut1, Float_t cut2);
	Bool_t DoubleJetForward(Float_t cut1, Float_t cut2);	
	Bool_t DoubleJet_Eta1p7_deltaEta4(Float_t cut1, Float_t cut2);
	Bool_t TripleJetCentral(Float_t cut1, Float_t cut2, Float_t cut3);
	Bool_t TripleJet_VBF(Float_t cut1, Float_t cut2, Float_t cut3);

	Bool_t QuadJetCentral(Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4);
	Bool_t MultiJet(Int_t nj, Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4);	
	Bool_t DoubleTauJetEta(Float_t cut1, Float_t cut2, Float_t etaCut=4.5);
	Bool_t isoTau_Tau(Float_t cut1, Float_t cut2, Float_t etaCut=4.5);
	Bool_t DoubleIsoTau(Float_t cut1, Float_t cut2, Float_t etaCut=4.5);


// -- Sums
	Bool_t ETT(Float_t ETTcut);
	Bool_t HTT(Float_t HTTcut);
	Bool_t ETM(Float_t ETMcut);
   Bool_t HTM(Float_t HTMcut);

// -- Egamma
	Bool_t SingleEG(Float_t cut);
	Bool_t SingleEG_Eta(Float_t cut, Float_t etaCut=4.5);
	Bool_t SingleIsoEG_Eta(Float_t cut, Float_t etaCut=4.5);
	Bool_t isoEG_EG(Float_t cut1, Float_t cut2);
	Bool_t EG_isoEG(Float_t cut1, Float_t cut2);
	Bool_t DoubleEG(Float_t cut1, Float_t cut2);
	Bool_t DoubleIsoEG(Float_t cut1, Float_t cut2);
	Bool_t TripleEG(Float_t cut1, Float_t cut2, Float_t cut3);

// -- Muons 
	Bool_t SingleMu(Float_t ptcut, Int_t qualmin=4);
	Bool_t SingleMuEta(Float_t ptcut, Float_t etaCut=2.1, Int_t qualmin=4);
	Bool_t SingleIsoMuEta(Float_t ptcut, Float_t etaCut=2.1, Int_t qualmin=4);
	Bool_t DoubleMu(Float_t cut1, Float_t cut2, Int_t qualmin=4);	// on top of DoubleMu3
//	Bool_t DoubleIsoMu(Float_t cut1, Float_t cut2, Int_t qualmin=4);
	Bool_t DoubleMuHighQEtaCut(Float_t ptcut, Float_t etacut);
	Bool_t TripleMu(Float_t cut1, Float_t cut2, Float_t cut3, Int_t qualmin);	// on top of DoubleMu3
	Bool_t DoubleMuXOpen(Float_t ptcut);	// on top of SingleMu7
	Bool_t Onia(Float_t ptcut1, Float_t ptcut2, Float_t etacut, Int_t delta);   

	void Loop(Int_t calcThreshold, Int_t selectDataInput, TString lsRunFile, TString L1MenuFile, const int n_events_=-1);
	void fillDataStructure(Int_t selectDataInput=0);

	private :

	Bool_t PhysicsBits[128];
	Bool_t first;

	Int_t insert_ibin;
	Bool_t insert_val[100];
	std::string insert_names[100];

	Int_t NBITS_TRIGS;


};


double L1Menu2015::getPUWeight(float mcMeanInt){

    double wgt = 1.0;
	 
	 int bin = puWeightsHist->FindBin(mcMeanInt);
	 wgt = puWeightsHist->GetBinContent(bin);
	 
    return wgt;
}


void L1Menu2015::calculatePUWeighting(TString mcHistFile, TString dataHistFile, bool dumpWgts){

    TFile *mcF   = new TFile(mcHistFile);
	 TH1F  *mcH   = (TH1F*)mcF->Get("pileup")->Clone();
	 mcH->Scale(1./mcH->Integral());
	 
	 TFile *dataF = new TFile(dataHistFile);
	 puWeightsHist = (TH1F*)dataF->Get("pileup")->Clone();
    puWeightsHist->Scale(1./puWeightsHist->Integral());
	 puWeightsHist->Divide(mcH);

	 
//	 mcF->Close();
//	 dataF->Close();
  
    if(dumpWgts) {
	   printf("PU Weights for Reweighting Monte Carlo\n");
		for(int i=1; i<=puWeightsHist->GetNbinsX()+1; i++) printf("PU %f  Wght %f\n",puWeightsHist->GetBinCenter(i),puWeightsHist->GetBinContent(i));
	 }
	 	
    return;
}

double L1Menu2015::calculateHTT(){


   double httValue = 0.;
	
// Calculate our own HT and HTM from the jets that survive the double jet removal.
// Also do not include tau jets because they are already in the jet list (would be double counting)
       for(unsigned int i=0; i<myEvt_.Njet; i++) {
		    if(myEvt_.Bxjet.at(i)==0 && !myEvt_.Taujet.at(i)) {
			    if(myEvt_.Etajet.at(i) > 4 and myEvt_.Etajet.at(i) < 17) {
				    httValue += myEvt_.Etjet.at(i);				 
				 } //in proper eta range
			 } //correct beam crossing
		 } //loop over cleaned jets 	  
	
	return httValue;
}

double L1Menu2015::calculateHTM(){

   double htmValue =0.;
	double htmValueX=0.;
	double htmValueY=0.;

// Calculate our own HT and HTM from the jets that survive the double jet removal.
// Also do not include tau jets because they are already in the jet list (would be double counting)
       for(unsigned int i=0; i<myEvt_.Njet; i++) {
		    if(myEvt_.Bxjet.at(i)==0 && !myEvt_.Taujet.at(i)) {
			    if(myEvt_.Etajet.at(i) > 4 and myEvt_.Etajet.at(i) < 17) {

//  Get the phi angle  towers are 0-17 (this is probably not real mapping but OK for just magnitude of HTM
                float phi = TMath::TwoPi()*(myEvt_.Phijet.at(i)/18.);
					 htmValueX += cos(phi)*myEvt_.Etjet.at(i);
					 htmValueY += sin(phi)*myEvt_.Etjet.at(i);
				 
				 } //in proper eta range
			 } //correct beam crossing
		 } //loop over cleaned jets 	  

   	 htmValue     = sqrt(htmValueX*htmValueX + htmValueY*htmValueY) ; 
				 


	return htmValue;
}

void L1Menu2015::dumpEvent(int dumpCode) {

			
	  // Dump Event Information
   	  printf("\n");
		  printf("---------------------------------------------\n");
		  if( (dumpCode & 0x1) != 0) printf("Run %10i  Lumi Section %i Event %10i \n",myEvt_.Run,myEvt_.LS,myEvt_.Event);


	  // Muon Info
   	  if( (dumpCode & 0x2) != 0) {
   		 printf(" ======== Muons %i =========\n",myEvt_.Nmu);
   		 for(int i=0; i<myEvt_.Nmu; i++) {
	   		 printf(" %1i)",i);
				 printf(" Bx  %3i ",  myEvt_.Bxmu.at(i));
				 printf(" Pt  %6.2f ",myEvt_.Ptmu.at(i));
				 printf(" Eta %6.2f ",myEvt_.Etamu.at(i));
				 printf(" Phi %6.2f ",myEvt_.Phimu.at(i));
				 printf(" Qua %3i ",  myEvt_.Qualmu.at(i));
				 printf(" Iso %6.2f ",myEvt_.Isomu.at(i));
				 printf("\n");
			 } 
		  }


   	  if( (dumpCode & 0x4) != 0) {
   		  printf(" ======== IsoEM =========\n");
   		  for(int i=0; i<myEvt_.Nele; i++) {
				if(myEvt_.Isoel.at(i)) {
	   		  printf(" %1i)",i);
				  printf(" Bx  %2i ",  myEvt_.Bxel.at(i)); 
				  printf(" Et  %5.2f ",myEvt_.Etel.at(i));  	
      		  printf(" Eta %5.2f ",myEvt_.Etael.at(i));
				  printf(" Phi %5.2f ",myEvt_.Phiel.at(i));
				  printf("\n");
				}	
			  }
		  }

   	  if( (dumpCode & 0x8) != 0) {
   		  printf(" ======== NonIsoEM =========\n");
   		  for(int i=0; i<myEvt_.Nele; i++) {
				if(!myEvt_.Isoel.at(i)) {
	   		  printf(" %1i)",i);
				  printf(" Bx  %2i ",  myEvt_.Bxel.at(i)); 
				  printf(" Et  %5.2f ",myEvt_.Etel.at(i));  	
      		  printf(" Eta %5.2f ",myEvt_.Etael.at(i));
				  printf(" Phi %5.2f ",myEvt_.Phiel.at(i));
				  printf("\n");
				}	
			  }
   	  }

   	  if( (dumpCode & 0x40) != 0) {	
   		 printf(" ======== Tau Jets =========\n");
   		 for(int i=0; i<myEvt_.Njet; i++) {
			  if(myEvt_.Taujet.at(i)) {
	   		 printf(" %1i)",i);
				 printf(" Bx  %2i ",  myEvt_.Bxjet.at(i)); 
				 printf(" Et  %5.2f ",myEvt_.Etjet.at(i));		
      		 printf(" Eta %5.2f ",myEvt_.Etajet.at(i));
				 printf(" Phi %5.2f ",myEvt_.Phijet.at(i));
				 (myEvt_.isoTaujet.at(i)) ? printf(" Isol") : printf(" NonIso");
				 printf("\n");
			  }
			 }
   	  }


   	  if( (dumpCode & 0x10) != 0) {
   		 printf(" ======== Central Jets =========\n");
   		 for(int i=0; i<myEvt_.Njet; i++) {
			  if(!myEvt_.Taujet.at(i) && !myEvt_.Fwdjet.at(i) ) {
	   		 printf(" %1i)",i);
				 printf(" Bx  %2i ",  myEvt_.Bxjet.at(i)); 
				 printf(" Et  %5.2f ",myEvt_.Etjet.at(i));		
      		 printf(" Eta %5.2f ",myEvt_.Etajet.at(i));
				 printf(" Phi %5.2f ",myEvt_.Phijet.at(i));
				 printf("\n");
			  }
			 }
		  }  

   	  if( (dumpCode & 0x20) != 0) {
   		 printf(" ======== Forward Jets =========\n");
   		 for(int i=0; i<myEvt_.Njet; i++) {
			  if(myEvt_.Fwdjet.at(i) ) {
	   		 printf(" %1i)",i);
				 printf(" Bx  %2i ",  myEvt_.Bxjet.at(i)); 
				 printf(" Et  %5.2f ",myEvt_.Etjet.at(i));		
      		 printf(" Eta %5.2f ",myEvt_.Etajet.at(i));
				 printf(" Phi %5.2f ",myEvt_.Phijet.at(i));
				 printf("\n");
			  }
			 }
   	  }



	  // Total Energy information
   	  if( (dumpCode & 0x80) != 0) {
   		 printf(" ========= Global Info ===========\n");
			 printf(" EtTot  %6.2f   HTT     %6.2f  \n",  myEvt_.ETT,myEvt_.HTT);   
			 printf(" EtMiss %6.2f   EtMissPhi %4i  \n",myEvt_.ETM,myEvt_.PhiETM);   
			 printf(" HtMiss %6.2f   HtMissPhi %4i  \n",myEvt_.HTM,myEvt_.PhiHTM);   
   		 printf("======================================\n");	
   	  }

  return;

}


void L1Menu2015::fillDataStructure(Int_t selectDataInput) {
   
//	 printf("Entering fillDataStucture %i %i \n",selectDataInput,usePUWeight);
	 myEvt_.Reset();



// Fill a few quantities for this event
	 puWeight = 1.0;
//	 if(usePUWeight) puWeight = event_->puWeight;
//	 if(usePUWeight) puWeight = getPUWeight(sim_->meanInt);	 
	 if(puWeight<0.) printf("WARNING: Bad PU Weight.  Value %f \n",puWeight);
	 h_puWeight->Fill(puWeight);
	 h_lumiSec->Fill(event_->lumi);
	 h_bx->Fill(event_->bx);
//	 if(sim_->meanInt>0.) h_Sim_meanInt->Fill(sim_->meanInt); 


// Grab standard event information
    myEvt_.Run   = event_->run;
	 myEvt_.LS    = event_->lumi;
	 myEvt_.Event = event_->event;

    int nMu_bx0=0, nIsoMu_bx0=0;
    int nIsoEG_bx0=0, nNIsoEG_bx0=0 ;
	 int nCJets_bx0 = 0, nFJets_bx0 = 0;
	 int nTJets_bx0=0, nIsoTJets_bx0=0;
/* =======================================================================================================
/    Select the input source information
/ ---------------------------------------------------------------------------------------------------------
/    case  0: Use Original L1ExtraTree that came with the event
/    case  1: Use reEmulated L1ExtraTree that was produced (note this may not be present or may be identical to the Original tree)
/    case  2: Use Original L1ExtraTree that came with the event except for muons which get from GMT. (For old Ntuple with no quality flag in L1Extra)
/
/	  case 10: Use Original L1Tree (GMT/GT) that came with the event
/    case 11: Use reEmulated L1Tree (GMT/GT) that was produced (note this may not be present or may be identical to the Original tree)
/    case 12: Use Original L1Tree (GMT/GCT) that was produced 
/    case 13: Use reEmulated L1Tree (GMT/GCT) that was produced (note this may not be present or may be identical to the Original tree)
/
/    case 21: Use the L1ExtraUpgradeTree (Assuming Stage 1 Quantities filled in tree)
/    case 22: Use the L1ExtraUpgradeTree (Assuming Stage 2 Quantities filled in tree)
/ ======================================================================================================= */
    switch (selectDataInput) {
	 	 
    case 0:  //Select from default L1ExtraTree

// Grab the iso first  
		 for(unsigned int i=0; i<l1extra_->nIsoEm; i++) {


      	 myEvt_.Bxel.push_back(l1extra_->isoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->isoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->isoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->isoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
		 
		 
		 for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++) {

      	 myEvt_.Bxel.push_back(l1extra_->nonIsoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->nonIsoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->nonIsoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->nonIsoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
			 myEvt_.Isoel.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);	 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


		 for(unsigned int i=0; i< l1extra_->nCenJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->cenJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->cenJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->cenJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->cenJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1extra_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);
		 

		 for(unsigned int i=0; i< l1extra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1extra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

// Add these taus to the central jet list in order for the old system to remaind the same.
//
//  In the upgraded system clusters that appear in the tau list also appear in the the
//  jet list (with perhaps a slightly different energy).  Therefore the trigger algorithms
//  in this macro never treat a tau jet as an average jet (to avoid double counting in the
//   upgrade side of things). However, in the old system the taus appeared only in the jet
//   list and were treated as average jets.  So here we put in another copy of the taus by
//   hand but do not flag them as taus.
//
		 for(unsigned int i=0; i< l1extra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1extra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 


	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1extra_->nMet; i++) {
  		     if(l1extra_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1extra_->et.at(i) ; 
      	    myEvt_.ETM     = l1extra_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1extra_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
   	  
   	 for(unsigned int i=0; i< l1extra_->nMht; i++) {
  		     if(l1extra_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1extra_->ht.at(i) ; 
   	       myEvt_.HTM     = l1extra_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1extra_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra
	
	
// Get the muon information  
		 for(unsigned int i=0; i<l1extra_->nMuons; i++) {

      	 myEvt_.Bxmu.push_back(l1extra_->muonBx.at(i));
      	 myEvt_.Ptmu.push_back(l1extra_->muonEt.at(i));
      	 myEvt_.Phimu.push_back(l1extra_->muonPhi.at(i)); 
      	 myEvt_.Etamu.push_back(l1extra_->muonEta.at(i)); 
			 myEvt_.Qualmu.push_back(l1extra_->muonQuality.at(i));
          myEvt_.Isomu.push_back(l1extra_->muonIso.at(i));

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight);	


       break;

    case 1:  //Select from reEmulated L1ExtraTree

// Grab the iso first  
		 for(unsigned int i=0; i<l1emuextra_->nIsoEm; i++) {


      	 myEvt_.Bxel.push_back(l1emuextra_->isoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1emuextra_->isoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1emuextra_->isoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1emuextra_->isoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
		 
		 
		 for(unsigned int i=0; i<l1emuextra_->nNonIsoEm; i++) {

      	 myEvt_.Bxel.push_back(l1emuextra_->nonIsoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1emuextra_->nonIsoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1emuextra_->nonIsoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1emuextra_->nonIsoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
			 myEvt_.Isoel.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);	 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


		 for(unsigned int i=0; i< l1emuextra_->nCenJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1emuextra_->cenJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1emuextra_->cenJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1emuextra_->cenJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1emuextra_->cenJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1emuextra_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1emuextra_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1emuextra_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1emuextra_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1emuextra_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);
		 

		 for(unsigned int i=0; i< l1emuextra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1emuextra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1emuextra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1emuextra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1emuextra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

// Add these taus to the central jet list in order for the old system to remaind the same.
//
//  In the upgraded system clusters that appear in the tau list also appear in the the
//  jet list (with perhaps a slightly different energy).  Therefore the trigger algorithms
//  in this macro never treat a tau jet as an average jet (to avoid double counting in the
//   upgrade side of things). However, in the old system the taus appeared only in the jet
//   list and were treated as average jets.  So here we put in another copy of the taus by
//   hand but do not flag them as taus.
//
		 for(unsigned int i=0; i< l1emuextra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1emuextra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1emuextra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1emuextra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1emuextra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 



	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1emuextra_->nMet; i++) {
  		     if(l1emuextra_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1emuextra_->et.at(i) ; 
      	    myEvt_.ETM     = l1emuextra_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1emuextra_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
   	  
   	 for(unsigned int i=0; i< l1emuextra_->nMht; i++) {
  		     if(l1emuextra_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1emuextra_->ht.at(i) ; 
   	       myEvt_.HTM     = l1emuextra_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1emuextra_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra
	
	
// Get the muon information  
		 for(unsigned int i=0; i<l1emuextra_->nMuons; i++) {

      	 myEvt_.Bxmu.push_back(l1emuextra_->muonBx.at(i));
      	 myEvt_.Ptmu.push_back(l1emuextra_->muonEt.at(i));
      	 myEvt_.Phimu.push_back(l1emuextra_->muonPhi.at(i)); 
      	 myEvt_.Etamu.push_back(l1emuextra_->muonEta.at(i)); 
			 myEvt_.Qualmu.push_back(l1emuextra_->muonQuality.at(i));
          myEvt_.Isomu.push_back(l1emuextra_->muonIso.at(i));

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight);	


       break;

    case 2:  //Select from default L1ExtraTree except muons (For use with old L1ExtraTrees which were missing Muon Quality)

// Grab the iso first  
		 for(unsigned int i=0; i<l1extra_->nIsoEm; i++) {


      	 myEvt_.Bxel.push_back(l1extra_->isoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->isoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->isoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->isoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
		 
		 
		 for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++) {

      	 myEvt_.Bxel.push_back(l1extra_->nonIsoEmBx.at(i));
      	 myEvt_.Etel.push_back(l1extra_->nonIsoEmEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1extra_->nonIsoEmPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1extra_->nonIsoEmEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
			 myEvt_.Isoel.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);	 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


		 for(unsigned int i=0; i< l1extra_->nCenJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->cenJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->cenJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->cenJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->cenJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1extra_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1extra_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);
		 

		 for(unsigned int i=0; i< l1extra_->nTauJets; i++) {

      	 myEvt_.Bxjet.push_back(l1extra_->tauJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1extra_->tauJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1extra_->tauJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1extra_->tauJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1extra_->nMet; i++) {
  		     if(l1extra_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1extra_->et.at(i) ; 
      	    myEvt_.ETM     = l1extra_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1extra_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
   	  
   	 for(unsigned int i=0; i< l1extra_->nMht; i++) {
  		     if(l1extra_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1extra_->ht.at(i) ; 
   	       myEvt_.HTM     = l1extra_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1extra_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra
	
	
// Get the muon information  (From GMT because old L1ExtraTrees were missing Muon quality)
		 for(int i=0; i<gmt_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmt_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmt_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmt_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmt_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmt_->Qual[i]);
			 myEvt_.Isomu.push_back(gmt_->Isol[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0);	

       break;
		 
		 
		 
// Extract the default quantities from the GT/GMT
// ==============================================		 
	 case 10:  	 


// Get the muon information  
		 for(int i=0; i<gmt_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmt_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmt_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmt_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmt_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmt_->Qual[i]);
			 myEvt_.Isomu.push_back(gmt_->Isol[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0);	

// EG
		 for(int i=0; i< gt_->Nele; i++) {      	 
      	 myEvt_.Bxel.push_back(gt_->Bxel[i]);
      	 myEvt_.Etel.push_back(gt_->Rankel[i]);
      	 myEvt_.Phiel.push_back(gt_->Phiel[i]);
      	 myEvt_.Etael.push_back(gt_->Etael[i]);
      	 myEvt_.Isoel.push_back(gt_->Isoel[i]);
			 
// Histogram Quantities          		 
          if(gt_->Bxel[i]==0) {
			   if(myEvt_.Isoel.at(myEvt_.Nele)) {
			     h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nIsoEG_bx0++;
				} else {
			     h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nNIsoEG_bx0++;
				}  				
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0);
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0);
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(int i=0; i< gt_->Njet; i++) {

      	 myEvt_.Bxjet.push_back(gt_->Bxjet[i]);
      	 myEvt_.Etjet.push_back((gt_->Rankjet[i])*4.);
      	 myEvt_.Phijet.push_back(gt_->Phijet[i]);
      	 myEvt_.Etajet.push_back(gt_->Etajet[i]);
      	 myEvt_.Taujet.push_back(gt_->Taujet[i]);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(gt_->Fwdjet[i]);

          if(gt_->Bxjet[i]==0) {
			    if(myEvt_.Taujet.at(myEvt_.Njet)){
  			      h_TJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nTJets_bx0++;
				 } else if(myEvt_.Fwdjet.at(myEvt_.Njet)){
  			      h_FJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nFJets_bx0++;
				 } else {
  			      h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nCJets_bx0++;				 
				 }	
			 }
			 myEvt_.Njet++;

		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
		 h_FJet_Njet->Fill(nFJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

	   // Fill energy sums
   	 myEvt_.ETT     = (gt_ -> RankETT)/2. ; 
   	 myEvt_.OvETT   = gt_ -> OvETT	; 
   	 myEvt_.HTT     = (gt_ -> RankHTT)/2. ; 
   	 myEvt_.OvHTT   = gt_ -> OvHTT	; 
   	 myEvt_.ETM     = (gt_ -> RankETM)/2. ; 
   	 myEvt_.PhiETM  = gt_ -> PhiETM  ; 
   	 myEvt_.OvETM   = gt_ -> OvETM	; 
   	 myEvt_.HTM     = (gt_ -> RankHTM)*2. ; 
   	 myEvt_.PhiHTM  = gt_ -> PhiHTM  ; 
   	 myEvt_.OvHTM   = gt_ -> OvHTM	; 
		 
// Fill Histograms
       h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
       h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
       h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight); 
       h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
       h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
       h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);  						
   
	    break;

// Extract the quantities from the reEmulated GT/GMT
// ==================================================		 
	 case 11:  	 


// Get the muon information  
		 for(int i=0; i<gmtEmu_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmtEmu_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmtEmu_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmtEmu_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmtEmu_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmtEmu_->Qual[i]);
			 myEvt_.Isomu.push_back(gmtEmu_->Isol[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0);	

// EG
		 for(int i=0; i< gtEmu_->Nele; i++) {      	 
      	 myEvt_.Bxel.push_back(gtEmu_->Bxel[i]);
      	 myEvt_.Etel.push_back(gtEmu_->Rankel[i]);
      	 myEvt_.Phiel.push_back(gtEmu_->Phiel[i]);
      	 myEvt_.Etael.push_back(gtEmu_->Etael[i]);
      	 myEvt_.Isoel.push_back(gtEmu_->Isoel[i]);
			 
// Histogram Quantities          		 
          if(gtEmu_->Bxel[i]==0) {
			   if(myEvt_.Isoel.at(myEvt_.Nele)) {
			     h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nIsoEG_bx0++;
				} else {
			     h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nNIsoEG_bx0++;
				}  				
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0);
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0);
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(int i=0; i< gtEmu_->Njet; i++) {

      	 myEvt_.Bxjet.push_back(gtEmu_->Bxjet[i]);
      	 myEvt_.Etjet.push_back((gtEmu_->Rankjet[i])*4.);
      	 myEvt_.Phijet.push_back(gtEmu_->Phijet[i]);
      	 myEvt_.Etajet.push_back(gtEmu_->Etajet[i]);
      	 myEvt_.Taujet.push_back(gtEmu_->Taujet[i]);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(gtEmu_->Fwdjet[i]);

          if(gtEmu_->Bxjet[i]==0) {
			    if(myEvt_.Taujet.at(myEvt_.Njet)){
  			      h_TJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nTJets_bx0++;
				 } else if(myEvt_.Fwdjet.at(myEvt_.Njet)){
  			      h_FJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nFJets_bx0++;
				 } else {
  			      h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nCJets_bx0++;				 
				 }	
			 }
			 myEvt_.Njet++;

		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
		 h_FJet_Njet->Fill(nFJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

	   // Fill energy sums
   	 myEvt_.ETT     = (gtEmu_ -> RankETT)/2. ; 
   	 myEvt_.OvETT   = gtEmu_ -> OvETT	; 
   	 myEvt_.HTT     = (gtEmu_ -> RankHTT)/2. ; 		 
   	 myEvt_.OvHTT   = gtEmu_ -> OvHTT	; 
   	 myEvt_.ETM     = (gtEmu_ -> RankETM)/2. ; 
   	 myEvt_.PhiETM  = gtEmu_ -> PhiETM  ; 
   	 myEvt_.OvETM   = gtEmu_ -> OvETM	; 
   	 myEvt_.HTM     = (gtEmu_ -> RankHTM)*2. ; 
   	 myEvt_.PhiHTM  = gtEmu_ -> PhiHTM  ; 
   	 myEvt_.OvHTM   = gtEmu_ -> OvHTM	; 
		 
// Fill Histograms
       h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
       h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
       h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight); 
       h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
       h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
       h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);  						
   
	    break;

// Extract the quantities from the standard GCT/GMT
// ==================================================		 
	 case 12:  	 


// Get the muon information  
		 for(int i=0; i<gmt_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmt_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmt_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmt_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmt_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmt_->Qual[i]);
			 myEvt_.Isomu.push_back(gmt_->Isol[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0);	

// EG
		 for(int i=0; i< gct_->IsoEmSize; i++) {      	 
      	 myEvt_.Bxel.push_back(gct_->IsoEmBx[i]);
      	 myEvt_.Etel.push_back(gct_->IsoEmRnk[i]);
      	 myEvt_.Phiel.push_back(gct_->IsoEmPhi[i]);
      	 myEvt_.Etael.push_back(gct_->IsoEmEta[i]);
      	 myEvt_.Isoel.push_back(true);
			 
// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			     h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nIsoEG_bx0++;
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0);

		 
// non isoEG
		 for(int i=0; i< gct_->NonIsoEmSize; i++) {      	 
      	 myEvt_.Bxel.push_back(gct_->NonIsoEmBx[i]);
      	 myEvt_.Etel.push_back(gct_->NonIsoEmRnk[i]);
      	 myEvt_.Phiel.push_back(gct_->NonIsoEmPhi[i]);
      	 myEvt_.Etael.push_back(gct_->NonIsoEmEta[i]);
      	 myEvt_.Isoel.push_back(false);
			 
// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			     h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nNIsoEG_bx0++;  				
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0);		 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(int i=0; i< gct_->CJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gct_->CJetBx[i]);
      	 myEvt_.Etjet.push_back((gct_->CJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gct_->CJetPhi[i]);
      	 myEvt_.Etajet.push_back(gct_->CJetEta[i]);
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

          if(gct_->CJetBx[i]==0) {
  			      h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nCJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);

		 
		 for(int i=0; i< gct_->FJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gct_->FJetBx[i]);
      	 myEvt_.Etjet.push_back((gct_->FJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gct_->FJetPhi[i]);
      	 myEvt_.Etajet.push_back(gct_->FJetEta[i]);
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

          if(gct_->FJetBx[i]==0) {
  			      h_FJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nFJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_FJet_Njet->Fill(nFJets_bx0,puWeight);
		 
		 for(int i=0; i< gct_->TJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gct_->TJetBx[i]);
      	 myEvt_.Etjet.push_back((gct_->TJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gct_->TJetPhi[i]);
      	 myEvt_.Etajet.push_back(gct_->TJetEta[i]);
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

          if(gct_->TJetBx[i]==0) {
  			      h_TJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nTJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
	 		 
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

// Loop over sum Et quantities
       for(int i=0; i<gct_->EtTotSize; i++) {
          if(gct_->EtTotBX[i]==0) myEvt_.ETT = (gct_ -> EtTot[i])/2.;
		 }  
		 myEvt_.OvETT   = 0.	;
		 
// Loop over sum Et quantities
       for(int i=0; i<gct_->EtMissSize; i++) {
          if(gct_->EtMissBX[i]==0) {
			   myEvt_.ETM = (gct_ -> EtMiss[i])/2.;
				myEvt_.PhiETM  = gct_ -> EtMissPhi[i]  ;
		    }
		 }  
		 myEvt_.OvETM   = 0.	;
		 
// Loop over sum Ht quantities
       for(int i=0; i<gct_->EtHadSize; i++) {
          if(gct_->EtHadBX[i]==0) myEvt_.HTT = (gct_ -> EtHad[i])/2.;
		 }  
		 myEvt_.OvHTT   = 0.	;
		 
// Loop over sum Ht miss quantities
       for(int i=0; i<gct_->HtMissSize; i++) {
          if(gct_->HtMissBX[i]==0) {
			   myEvt_.HTM = (gct_ -> HtMiss[i])/2.;
				myEvt_.PhiHTM  = gct_ -> HtMissPhi[i]  ;
		    }
		 }  
       myEvt_.OvHTM   = 0.	;

		 
// Fill Histograms
       h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
       h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
       h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight); 
       h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
       h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
       h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);  						
   
	    break;




// Extract the quantities from the reEmulated GCT/GMT
// ==================================================		 
	 case 13:  	 


// Get the muon information  
		 for(int i=0; i<gmtEmu_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmtEmu_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmtEmu_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmtEmu_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmtEmu_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmtEmu_->Qual[i]);
			 myEvt_.Isomu.push_back(gmtEmu_->Isol[i]);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
			    h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu));
				 h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu));
				 h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu));
				 nMu_bx0++;
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0);	

// EG
		 for(int i=0; i< gctEmu_->IsoEmSize; i++) {      	 
      	 myEvt_.Bxel.push_back(gctEmu_->IsoEmBx[i]);
      	 myEvt_.Etel.push_back(gctEmu_->IsoEmRnk[i]);
      	 myEvt_.Phiel.push_back(gctEmu_->IsoEmPhi[i]);
      	 myEvt_.Etael.push_back(gctEmu_->IsoEmEta[i]);
      	 myEvt_.Isoel.push_back(true);
			 
// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			     h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nIsoEG_bx0++;
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0);

		 
// non isoEG
		 for(int i=0; i< gctEmu_->NonIsoEmSize; i++) {      	 
      	 myEvt_.Bxel.push_back(gctEmu_->NonIsoEmBx[i]);
      	 myEvt_.Etel.push_back(gctEmu_->NonIsoEmRnk[i]);
      	 myEvt_.Phiel.push_back(gctEmu_->NonIsoEmPhi[i]);
      	 myEvt_.Etael.push_back(gctEmu_->NonIsoEmEta[i]);
      	 myEvt_.Isoel.push_back(false);
			 
// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			     h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				  h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				  nNIsoEG_bx0++;  				
			 }	 			 
          
			 myEvt_.Nele++;
		 }
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0);		 
		 // printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);

		 for(int i=0; i< gctEmu_->CJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gctEmu_->CJetBx[i]);
      	 myEvt_.Etjet.push_back((gctEmu_->CJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gctEmu_->CJetPhi[i]);
      	 myEvt_.Etajet.push_back(gctEmu_->CJetEta[i]);
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

          if(gctEmu_->CJetBx[i]==0) {
  			      h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nCJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);

		 
		 for(int i=0; i< gctEmu_->FJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gctEmu_->FJetBx[i]);
      	 myEvt_.Etjet.push_back((gctEmu_->FJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gctEmu_->FJetPhi[i]);
      	 myEvt_.Etajet.push_back(gctEmu_->FJetEta[i]);
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

          if(gctEmu_->FJetBx[i]==0) {
  			      h_FJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nFJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_FJet_Njet->Fill(nFJets_bx0,puWeight);
		 
		 for(int i=0; i< gctEmu_->TJetSize; i++) {

      	 myEvt_.Bxjet.push_back(gctEmu_->TJetBx[i]);
      	 myEvt_.Etjet.push_back((gctEmu_->TJetRnk[i])*4.);
      	 myEvt_.Phijet.push_back(gctEmu_->TJetPhi[i]);
      	 myEvt_.Etajet.push_back(gctEmu_->TJetEta[i]);
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

          if(gctEmu_->TJetBx[i]==0) {
  			      h_TJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				   h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				   h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
					nTJets_bx0++;
			 }
			 myEvt_.Njet++;

		 }
		 h_TJet_Njet->Fill(nTJets_bx0,puWeight);
	 		 
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);

// Loop over sum Et quantities
       for(int i=0; i<gctEmu_->EtTotSize; i++) {
          if(gctEmu_->EtTotBX[i]==0) myEvt_.ETT = (gctEmu_ -> EtTot[i])/2.;
		 }  
		 myEvt_.OvETT   = 0.	;
		 
// Loop over sum Et quantities
       for(int i=0; i<gctEmu_->EtMissSize; i++) {
          if(gctEmu_->EtMissBX[i]==0) {
			   myEvt_.ETM = (gctEmu_ -> EtMiss[i])/2.;
				myEvt_.PhiETM  = gctEmu_ -> EtMissPhi[i]  ;
		    }
		 }  
		 myEvt_.OvETM   = 0.	;
		 
// Loop over sum Ht quantities
       for(int i=0; i<gctEmu_->EtHadSize; i++) {
          if(gctEmu_->EtHadBX[i]==0) myEvt_.HTT = (gctEmu_ -> EtHad[i])/2.;
		 }  
		 myEvt_.OvHTT   = 0.	;
		 
// Loop over sum Ht miss quantities
       for(int i=0; i<gctEmu_->HtMissSize; i++) {
          if(gctEmu_->HtMissBX[i]==0) {
			   myEvt_.HTM = (gctEmu_ -> HtMiss[i])/2.;
				myEvt_.PhiHTM  = gctEmu_ -> HtMissPhi[i]  ;
		    }
		 }  
       myEvt_.OvHTM   = 0.	;

		 
// Fill Histograms
       h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
       h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
       h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight); 
       h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
       h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
       h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);  						
   
	    break;


    case 21:  //Select from L1ExtraUpgradeTree (Stage 1)

// Grab the iso first  
/*		 for(unsigned int i=0; i<l1upgrade_->nIsoEG; i++) {


      	 myEvt_.Bxel.push_back(l1upgrade_->isoEGBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->isoEGEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->isoEGPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->isoEGEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
*/		 
		 
// NOTES:  Stage 1 has EG Relaxed and EG Isolated.  The isolated EG are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.	
	
		 for(unsigned int i=0; i<l1upgrade_->nEG; i++) {

      	 myEvt_.Bxel.push_back(l1upgrade_->egBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->egEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->egPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->egEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord

// Check whether this EG is located in the isolation list
          bool isolated = false;
			 bool fnd = false;
			 int isoEG = 0;
			 while(!fnd && isoEG < l1upgrade_->nIsoEG) {
			    if(l1upgrade_->isoEGPhi.at(isoEG) == l1upgrade_->egPhi.at(i) &&
				    l1upgrade_->isoEGEta.at(isoEG) == l1upgrade_->egEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoEG++;
			 } 
			 myEvt_.Isoel.push_back(isolated);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0 && !isolated) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 else if (myEvt_.Bxel.at(myEvt_.Nele)==0 && isolated ) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	 
		 //printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


// Note:  Taus are in the jet list.  Decide what to do with them. For now
//  leave them the there as jets (not even flagged..)

		 for(unsigned int i=0; i< l1upgrade_->nJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1upgrade_->jetBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->jetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->jetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->jetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1upgrade_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1upgrade_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);


// NOTES:  Stage 1 has Tau Relaxed and TauIsolated.  The isolated Tau are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.			 

		 for(unsigned int i=0; i< l1upgrade_->nTau; i++) {

      	 myEvt_.Bxjet.push_back(l1upgrade_->tauBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->tauEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->tauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->tauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
      	 myEvt_.Fwdjet.push_back(false);

          bool isolated = false;
			 bool fnd = false;
			 int isoTau = 0;
			 while(!fnd && isoTau < l1upgrade_->nIsoTau) {
			    if(l1upgrade_->isoTauPhi.at(isoTau) == l1upgrade_->tauPhi.at(i) &&
				    l1upgrade_->isoTauEta.at(isoTau) == l1upgrade_->tauEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoTau++;
			 } 
			 myEvt_.isoTaujet.push_back(isolated);



// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0 && !isolated) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			 } else if (myEvt_.Bxjet.at(myEvt_.Njet)==0 && isolated) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;			 	 
			 }	 
			 myEvt_.Njet++;

		 }	
		 h_TJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);

/*		 
		 for(unsigned int i=0; i< l1upgrade_->nIsoTau; i++) {

      	 myEvt_.Bxjet.push_back(l1upgrade_->isoTauBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->isoTauEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->isoTauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->isoTauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 		 	 		 
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);
*/
	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1upgrade_->nMet; i++) {
  		     if(l1upgrade_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1upgrade_->et.at(i) ; 
      	    myEvt_.ETM     = l1upgrade_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1upgrade_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
   	  
   	 for(unsigned int i=0; i< l1upgrade_->nMht; i++) {
  		     if(l1upgrade_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1upgrade_->ht.at(i) ; 
   	       myEvt_.HTM     = l1upgrade_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1upgrade_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra
	
/*	
// Get the muon information  (not currently in the L1UpgradeTree)
		 for(unsigned int i=0; i<l1upgrade_->nMuons; i++) {

      	 myEvt_.Bxmu.push_back(l1upgrade_->muonBx.at(i));
      	 myEvt_.Ptmu.push_back(l1upgrade_->muonEt.at(i));
      	 myEvt_.Phimu.push_back(l1upgrade_->muonPhi.at(i)); 
      	 myEvt_.Etamu.push_back(l1upgrade_->muonEta.at(i)); 
			 myEvt_.Qualmu.push_back(l1upgrade_->muonQuality.at(i));
          myEvt_.Isomu.push_back(l1upgrade_->muonIso.at(i));

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	
*/
// Get the muon information  from reEmul GMT
		 for(int i=0; i<gmtEmu_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmtEmu_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmtEmu_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmtEmu_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmtEmu_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmtEmu_->Qual[i]);
			 myEvt_.Isomu.push_back(false);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	



       break;

    case 22:  //Select from L1ExtraUpgradeTree (Stage 2)

// Grab the iso first  
/*		 for(unsigned int i=0; i<l1upgrade_->nIsoEG; i++) {


      	 myEvt_.Bxel.push_back(l1upgrade_->isoEGBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->isoEGEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->isoEGPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->isoEGEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
*/		 
		 
// NOTES:  Stage 1 has EG Relaxed and EG Isolated.  The isolated EG are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.	
	
		 for(unsigned int i=0; i<l1upgrade_->nEG; i++) {

      	 myEvt_.Bxel.push_back(l1upgrade_->egBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->egEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->egPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->egEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord

// Check whether this EG is located in the isolation list
          bool isolated = false;
			 bool fnd = false;
			 int isoEG = 0;
			 while(!fnd && isoEG < l1upgrade_->nIsoEG) {
			    if(l1upgrade_->isoEGPhi.at(isoEG) == l1upgrade_->egPhi.at(i) &&
				    l1upgrade_->isoEGEta.at(isoEG) == l1upgrade_->egEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoEG++;
			 } 
			 myEvt_.Isoel.push_back(isolated);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0 && !isolated) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 else if (myEvt_.Bxel.at(myEvt_.Nele)==0 && isolated ) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }
			 myEvt_.Nele++;
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	 
		 //printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


// Note:  Taus are in the jet list.  Decide what to do with them. For now
//  leave them the there as jets (not even flagged..)
		 for(unsigned int i=0; i< l1upgrade_->nJets; i++) {


// For each jet look for a possible duplicate if so remove it.
         bool duplicate = false;
			for(unsigned int j=0; j<i; j++) {
			   if(  l1upgrade_->jetBx.at(i) == l1upgrade_->jetBx.at(j) &&
				     l1upgrade_->jetEt.at(i) == l1upgrade_->jetEt.at(j) &&
					  l1upgrade_->jetEta.at(i)== l1upgrade_->jetEta.at(j) &&
					  l1upgrade_->jetPhi.at(i)== l1upgrade_->jetPhi.at(j) ) {
				 duplicate = true;
				 //printf("Duplicate jet found and removed \n");
				} 
			}

      	if(!duplicate) { 
      	  myEvt_.Bxjet.push_back(l1upgrade_->jetBx.at(i));
      	  myEvt_.Etjet.push_back(l1upgrade_->jetEt.at(i));
      	  myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->jetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	  myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->jetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	  myEvt_.Taujet.push_back(false);
			  myEvt_.isoTaujet.push_back(false);
//      	  myEvt_.Fwdjet.push_back(false); //COMMENT OUT IF JET ETA FIX
			  
			 
//			  if(fabs(l1upgrade_->jetEta.at(i))>=3.0) printf("Et %f  Eta  %f  iEta  %f Phi %f  iPhi  %f \n",myEvt_.Etjet.at(myEvt_.Njet),l1upgrade_->jetEta.at(i),myEvt_.Etajet.at(myEvt_.Njet),l1upgrade_->jetPhi.at(i),myEvt_.Phijet.at(myEvt_.Njet));
           (fabs(l1upgrade_->jetEta.at(i))>=3.0) ?  myEvt_.Fwdjet.push_back(true) : myEvt_.Fwdjet.push_back(false);
  
// Histogram Quantities          			 
           if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			     h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
			  	  h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				  h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				  nCJets_bx0++;
			  }
			  myEvt_.Njet++;
			}   	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1upgrade_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1upgrade_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

//			 printf("Et %f  Eta  %f  iEta  %f  Phi %f  iPhi  %f \n",myEvt_.Etjet.at(myEvt_.Njet),l1upgrade_->fwdJetEta.at(i),myEvt_.Etajet.at(myEvt_.Njet),l1upgrade_->fwdJetPhi.at(i),myEvt_.Phijet.at(myEvt_.Njet));



// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);


// NOTES:  Stage 1 has Tau Relaxed and TauIsolated.  The isolated Tau are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.			 

		 for(unsigned int i=0; i< l1upgrade_->nTau; i++) {

// remove duplicates
         bool duplicate = false;
			for(unsigned int j=0; j<i; j++) {
			   if(  l1upgrade_->tauBx.at(i) == l1upgrade_->tauBx.at(j) &&
				     l1upgrade_->tauEt.at(i) == l1upgrade_->tauEt.at(j) &&
					  l1upgrade_->tauEta.at(i)== l1upgrade_->tauEta.at(j) &&
					  l1upgrade_->tauPhi.at(i)== l1upgrade_->tauPhi.at(j) ) {
				 duplicate = true;
				 //printf("Duplicate jet found and removed \n");
				} 
			}

         if(!duplicate) { 
      	  myEvt_.Bxjet.push_back(l1upgrade_->tauBx.at(i));
      	  myEvt_.Etjet.push_back(l1upgrade_->tauEt.at(i));
      	  myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->tauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	  myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->tauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	  myEvt_.Taujet.push_back(true);
      	  myEvt_.Fwdjet.push_back(false);

           bool isolated = false;
			  bool fnd = false;
			  int isoTau = 0;
			  while(!fnd && isoTau < l1upgrade_->nIsoTau) {
			    if(l1upgrade_->isoTauPhi.at(isoTau) == l1upgrade_->tauPhi.at(i) &&
				    l1upgrade_->isoTauEta.at(isoTau) == l1upgrade_->tauEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoTau++;
			  } 
			  myEvt_.isoTaujet.push_back(isolated);



// Histogram Quantities          			 
           if(myEvt_.Bxjet.at(myEvt_.Njet)==0 && !isolated) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			  } else if (myEvt_.Bxjet.at(myEvt_.Njet)==0 && isolated) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;			 	 
			  }	 
			  myEvt_.Njet++;
          } // duplicate check 
		 }	
		 h_TJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);

/*		 
		 for(unsigned int i=0; i< l1upgrade_->nIsoTau; i++) {

      	 myEvt_.Bxjet.push_back(l1upgrade_->isoTauBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->isoTauEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->isoTauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->isoTauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 		 	 		 
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);
*/

	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1upgrade_->nMet; i++) {
//  		     if(l1upgrade_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1upgrade_->et.at(i) ; 
      	    myEvt_.ETM     = l1upgrade_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1upgrade_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
//			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
		 

   	 for(unsigned int i=0; i< l1upgrade_->nMht; i++) {
	     if(l1upgrade_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = calculateHTT(); //l1upgrade_->ht.at(i) ; 
   	       myEvt_.HTM     = calculateHTM(); //l1upgrade_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = 0.; //l1upgrade_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra

// Histogram Quantities
       h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
		 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
		 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);	

	
/*	
// Get the muon information  (not currently in the L1UpgradeTree)
		 for(unsigned int i=0; i<l1upgrade_->nMuons; i++) {

      	 myEvt_.Bxmu.push_back(l1upgrade_->muonBx.at(i));
      	 myEvt_.Ptmu.push_back(l1upgrade_->muonEt.at(i));
      	 myEvt_.Phimu.push_back(l1upgrade_->muonPhi.at(i)); 
      	 myEvt_.Etamu.push_back(l1upgrade_->muonEta.at(i)); 
			 myEvt_.Qualmu.push_back(l1upgrade_->muonQuality.at(i));
          myEvt_.Isomu.push_back(l1upgrade_->muonIso.at(i));

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	
*/
// Get the muon information  from reEmul GMT
		 for(int i=0; i<gmtEmu_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmtEmu_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmtEmu_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmtEmu_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmtEmu_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmtEmu_->Qual[i]);
			 myEvt_.Isomu.push_back(false);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	



       break;

    case 23:  //Select from L1ExtraUpgradeTree (Stage 1)

// Grab the iso first  
/*		 for(unsigned int i=0; i<l1upgrade_->nIsoEG; i++) {


      	 myEvt_.Bxel.push_back(l1upgrade_->isoEGBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->isoEGEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->isoEGPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->isoEGEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
       	 myEvt_.Isoel.push_back(true);

// Histogram Quantities          		 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }	 
			 myEvt_.Nele++;

		 }
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	
*/		 
		 
// NOTES:  Stage 1 has EG Relaxed and EG Isolated.  The isolated EG are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.	
	
		 for(unsigned int i=0; i<l1upgrade_->nEG; i++) {
        if(l1upgrade_->egEt.at(i)<64) {  
      	 myEvt_.Bxel.push_back(l1upgrade_->egBx.at(i));
      	 myEvt_.Etel.push_back(l1upgrade_->egEt.at(i));
      	 myEvt_.Phiel.push_back(phiINjetCoord(l1upgrade_->egPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etael.push_back(etaINjetCoord(l1upgrade_->egEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord

// Check whether this EG is located in the isolation list
          bool isolated = false;
			 bool fnd = false;
			 int isoEG = 0;
			 while(!fnd && isoEG < l1upgrade_->nIsoEG) {
			    if(l1upgrade_->isoEGPhi.at(isoEG) == l1upgrade_->egPhi.at(i) &&
				    l1upgrade_->isoEGEta.at(isoEG) == l1upgrade_->egEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoEG++;
			 } 
			 myEvt_.Isoel.push_back(isolated);

// Histogram Quantities          			 
          if(myEvt_.Bxel.at(myEvt_.Nele)==0 && !isolated) {
			    h_nIsoEG_Et-> Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_nIsoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nNIsoEG_bx0++;
			 }	 else if (myEvt_.Bxel.at(myEvt_.Nele)==0 && isolated ) {
			    h_isoEG_Et->Fill(myEvt_.Etel.at(myEvt_.Nele),puWeight);
				 h_isoEG_Eta->Fill(myEvt_.Etael.at(myEvt_.Nele),puWeight);
				 h_isoEG_Phi->Fill(myEvt_.Phiel.at(myEvt_.Nele),puWeight);
				 nIsoEG_bx0++;
			 }
			 myEvt_.Nele++;
		  }//veto EG greater than 63	 
		 }	
		 h_nIsoEG_Nele->Fill(nNIsoEG_bx0,puWeight);
		 h_isoEG_Nele->Fill(nIsoEG_bx0,puWeight);	 
		 //printf("Number of Electrons in myEvt %i \n",myEvt_.Nele);


// Note:  Taus are in the jet list.  Decide what to do with them. For now
//  leave them the there as jets (not even flagged..)

		 for(unsigned int i=0; i< l1upgrade_->nJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1upgrade_->jetBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->jetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->jetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->jetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);
			 
// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_CJet_Et ->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_CJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_CJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nCJets_bx0++;
			 }
			 myEvt_.Njet++;	 
		 }
		 h_CJet_Njet->Fill(nCJets_bx0,puWeight);
		 
		 for(unsigned int i=0; i< l1upgrade_->nFwdJets; i++) {

      	 
      	 myEvt_.Bxjet.push_back(l1upgrade_->fwdJetBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->fwdJetEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->fwdJetPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->fwdJetEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(false);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(true);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_FJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_FJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_FJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nFJets_bx0++;
			 }	 
			 myEvt_.Njet++;
		 }
		 h_FJet_Njet->Fill( nFJets_bx0,puWeight);


// NOTES:  Stage 1 has Tau Relaxed and TauIsolated.  The isolated Tau are a subset of the Relaxed.
//         so sort through the relaxed list and flag those that also appear in the isolated list.			 

		 for(unsigned int i=0; i< l1upgrade_->nTau; i++) {
         if(l1upgrade_->tauEt.at(i)<64) {
      	 myEvt_.Bxjet.push_back(l1upgrade_->tauBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->tauEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->tauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->tauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
      	 myEvt_.Fwdjet.push_back(false);

          bool isolated = false;
			 bool fnd = false;
			 int isoTau = 0;
			 while(!fnd && isoTau < l1upgrade_->nIsoTau) {
			    if(l1upgrade_->isoTauPhi.at(isoTau) == l1upgrade_->tauPhi.at(i) &&
				    l1upgrade_->isoTauEta.at(isoTau) == l1upgrade_->tauEta.at(i)) {
				    isolated = true;
					 fnd = true;
				 }	 
			    isoTau++;
			 } 
			 myEvt_.isoTaujet.push_back(isolated);



// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0 && !isolated) {
			    h_TJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_TJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_TJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nTJets_bx0++;
			 } else if (myEvt_.Bxjet.at(myEvt_.Njet)==0 && isolated) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;			 	 
			 }	 
			 myEvt_.Njet++;
         } //veto taus greater than 63
		 }	
		 h_TJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);

/*		 
		 for(unsigned int i=0; i< l1upgrade_->nIsoTau; i++) {

      	 myEvt_.Bxjet.push_back(l1upgrade_->isoTauBx.at(i));
      	 myEvt_.Etjet.push_back(l1upgrade_->isoTauEt.at(i));
      	 myEvt_.Phijet.push_back(phiINjetCoord(l1upgrade_->isoTauPhi.at(i))); //PROBLEM: real value, trigger wants bin convert with phiINjetCoord
      	 myEvt_.Etajet.push_back(etaINjetCoord(l1upgrade_->isoTauEta.at(i))); //PROBLEM: real value, trigger wants bin convert with etaINjetCoord
      	 myEvt_.Taujet.push_back(true);
			 myEvt_.isoTaujet.push_back(false);
      	 myEvt_.Fwdjet.push_back(false);

// Histogram Quantities          			 
          if(myEvt_.Bxjet.at(myEvt_.Njet)==0) {
			    h_isoTJet_Et->Fill(myEvt_.Etjet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Eta->Fill(myEvt_.Etajet.at(myEvt_.Njet),puWeight);
				 h_isoTJet_Phi->Fill(myEvt_.Phijet.at(myEvt_.Njet),puWeight);
				 nIsoTJets_bx0++;
			 }	 
			 myEvt_.Njet++;

		 }		 		 		 	 		 
		 h_isoTJet_Njet->Fill(nIsoTJets_bx0,puWeight);
		// printf("Number of Jets in myEvt %i \n",myEvt_.Njet);
*/
	   // Fill energy sums  (Are overflow flags accessible in l1extra?)
   	 for(unsigned int i=0; i< l1upgrade_->nMet; i++) {
  		     if(l1upgrade_->metBx.at(i)==0) {
			    myEvt_.ETT     = l1upgrade_->et.at(i) ; 
      	    myEvt_.ETM     = l1upgrade_->met.at(i) ; 
   	       myEvt_.PhiETM  = l1upgrade_->metPhi.at(i)  ;
				 
// Histogram Quantities
             h_Sum_ETT->Fill(myEvt_.ETT,puWeight);
				 h_Sum_ETM->Fill(myEvt_.ETM,puWeight);
				 h_Sum_PhiETM->Fill(myEvt_.PhiETM,puWeight);				 
			  }
		 }	  	  		 	
   	 myEvt_.OvETT   = 0	;//not available in l1extra
		 myEvt_.OvETM   = 0	;//not available in l1extra  
   	  
   	 for(unsigned int i=0; i< l1upgrade_->nMht; i++) {
  		     if(l1upgrade_->mhtBx.at(i)==0) {   	 
   	       myEvt_.HTT     = l1upgrade_->ht.at(i) ; 
   	       myEvt_.HTM     = l1upgrade_->mht.at(i) ; 
   	       myEvt_.PhiHTM  = l1upgrade_->mhtPhi.at(i) ;
				 
// Histogram Quantities
             h_Sum_HTT->Fill(myEvt_.HTT,puWeight);
				 h_Sum_HTM->Fill(myEvt_.HTM,puWeight);
				 h_Sum_PhiHTM->Fill(myEvt_.PhiHTM,puWeight);		
		     }
		 }	  
   	 myEvt_.OvHTM   = 0	; //not available in l1extra
	    myEvt_.OvHTT   = 0	;//not available in l1extra
	
/*	
// Get the muon information  (not currently in the L1UpgradeTree)
		 for(unsigned int i=0; i<l1upgrade_->nMuons; i++) {

      	 myEvt_.Bxmu.push_back(l1upgrade_->muonBx.at(i));
      	 myEvt_.Ptmu.push_back(l1upgrade_->muonEt.at(i));
      	 myEvt_.Phimu.push_back(l1upgrade_->muonPhi.at(i)); 
      	 myEvt_.Etamu.push_back(l1upgrade_->muonEta.at(i)); 
			 myEvt_.Qualmu.push_back(l1upgrade_->muonQuality.at(i));
          myEvt_.Isomu.push_back(l1upgrade_->muonIso.at(i));

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	
*/
// Get the muon information  from reEmul GMT
		 for(int i=0; i<gmtEmu_->N; i++) {

      	 myEvt_.Bxmu.push_back(gmtEmu_->CandBx[i]);
      	 myEvt_.Ptmu.push_back(gmtEmu_->Pt[i]);
      	 myEvt_.Phimu.push_back(gmtEmu_->Phi[i]); 
      	 myEvt_.Etamu.push_back(gmtEmu_->Eta[i]); 
			 myEvt_.Qualmu.push_back(gmtEmu_->Qual[i]);
			 myEvt_.Isomu.push_back(false);

// Histogram Quantities
			 
          if(myEvt_.Bxmu.at(myEvt_.Nmu)==0) {
             if(myEvt_.Isomu.at(myEvt_.Nmu)==1) {
			      h_isoMu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_isoMu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nIsoMu_bx0++;
             } else {
   			   h_Mu_Et-> Fill(myEvt_.Ptmu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Eta->Fill(myEvt_.Etamu.at(myEvt_.Nmu),puWeight);
				   h_Mu_Phi->Fill(myEvt_.Phimu.at(myEvt_.Nmu),puWeight);
					nMu_bx0++;				 
				 }			 				 
			 }	 
			 myEvt_.Nmu++; 
		 }	      
       h_Mu_Nmu->Fill(nMu_bx0,puWeight); //only histogram non isolated muons
		 h_isoMu_Nmu->Fill(nIsoMu_bx0,puWeight);	



       break;


		 
	 default:
       std::cout << "---Not a valid input source FULL STOP! " << std::endl;
		 
		 break;
	 }	 
	 
	 return;
}


void L1Menu2015::InsertInMenu(std::string L1name, Bool_t value) {

	Bool_t post_prescale = false;

	Int_t prescale = 1;

	std::map<std::string, int>::const_iterator it = Prescales.find(L1name);
	if (it == Prescales.end() ) {
		std::cout << " --- NO PRESCALE DEFINED FOR " << L1name << " ---  SET P = 1 " << std::endl;
	}
	else {
		prescale = Prescales[L1name];
	}

	if (prescale >0) {
		Counts[L1name] ++;
		Int_t n = Counts[L1name];
		if ( n % prescale == 0) post_prescale = value; 
	}

	insert_names[insert_ibin] = L1name;
	insert_val[insert_ibin] = post_prescale ;

	insert_ibin ++;

}

Int_t L1Menu2015::L1BitNumber(std::string l1name) {

	std::map<std::string, int>::const_iterator it = BitMapping.find(l1name);
	if (it == BitMapping.end() ) {
		std::cout << " Wrong L1 name, not in BitMapping " << l1name << std::endl;
		return -1;
	}

	return BitMapping[l1name];
}

void L1Menu2015::FilL1Bits() {
        //printf("Run %i Event %i Alg Fired:",event_->run,event_->event);
	if(!gt_) {
	  for (Int_t ibit=0; ibit < 128; ibit++) {
		  PhysicsBits[ibit] = 0;
		  if (ibit<64) {
			  PhysicsBits[ibit] = (gt_->tw1[2]>>ibit)&1;
		  }
		  else {
			  PhysicsBits[ibit] = (gt_->tw2[2]>>(ibit-64))&1;
		  }
		  //if(PhysicsBits[ibit]) printf(" %i ",ibit);
	  }
	}else {
	  PhysicsBits[0] = 1; //set zero bias on if no gt information
	}	 
	//printf("\n");
}       


void L1Menu2015::InitL1Menu(TString menuFile) {

//=================================================================================
//   Read in the menu from file.  
// 
//   Format:
//	  - Each algorithm has a single line with the following information
//	  
//	  name(string)  bitNumber(int)  prescale(int)  Primary Threshold(float)  Secondary Threshold(float) Tertiary Theshold(float) 4th Threshold(float)   etaCut(float)   minMuQual(int)  target Bandwidth(float) fixed Bandwidth(int)  
//	  
//
//   Notes:
//	  - For each algorithm in the file there must be a corresponding call in the EvalMenu routine.
//	  - If the prescale is set to -1, the algorithm will not be run.
//	  
//===================================================================================




// Open File
   printf("\n===================  L1Menu2015:  Reading L1 Menu File %s =============================\n",menuFile.Data());
   ifstream ifs( menuFile );  
	
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
	  
	  printf("Alg: %20s %2i %2i %6.2f %6.2f %6.2f %6.2f %6.2f %2i %6.2f %2i %2i\n",algName.Data(),BitMapping[algName.Data()],Prescales[algName.Data()],
	         trigParList[algName.Data()].primTh,trigParList[algName.Data()].secTh,trigParList[algName.Data()].triTh,trigParList[algName.Data()].quadTh,
				trigParList[algName.Data()].etaCut,trigParList[algName.Data()].minQual,trigParList[algName.Data()].bandwidth,trigParList[algName.Data()].scalable,trigParList[algName.Data()].locked);
	}
   printf("====================================================================================================\n\n");

for (std::map<std::string, int>::iterator it=Prescales.begin(); it != Prescales.end(); it++) {
	std::string name = it -> first;
	Counts[name] = 0;
	Biased[name] = false; 
}


}









Bool_t L1Menu2015::SingleMuEta(Float_t ptcut, Float_t etaCut, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) { 
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < qualmin) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,ptcut,etaCut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etaCut) continue;
		   if (pt >= ptcut) muon = true;
	   }
	}

	Bool_t ok = muon;
	return ok;

}

Bool_t L1Menu2015::SingleIsoMuEta(Float_t ptcut, Float_t etaCut, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) { 
		Int_t bx = myEvt_.Bxmu.at(imu);		
//		if (bx != 0 || !myEvt_.Isomu.at(imu)) continue;
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < qualmin) continue;
		Float_t eta = myEvt_.Etamu.at(imu);
		
		if(useUpgradeMuons) {
   		if( UpgradeIsolatedMuon(pt,eta,ptcut,etaCut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etaCut) continue;
		   if (pt >= ptcut) muon = true;
	   }
		
//		if(pt>16. && ptcut<23.) printf("%i Muon Pt %f  Cut %f  Pass/Fail %i \n",imu,pt,ptcut,muon);
	}


	Bool_t ok = muon;
	return ok;

}

Bool_t L1Menu2015::SingleMu(Float_t ptcut, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		//    BX = 0, +/- 1 or +/- 2
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);
		Float_t eta = myEvt_.Etamu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < qualmin) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,ptcut,5.0)      ) muon = true;
 		} else {
		   if (pt >= ptcut) muon = true;
	   }	
			
	}

	Bool_t ok = muon;
	return ok;

}

Bool_t L1Menu2015::DoubleMuHighQEtaCut(Float_t ptcut, Float_t etacut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t nmu=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		//    BX = 0, +/- 1 or +/- 2
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);		

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,ptcut,etacut)      ) nmu ++;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= ptcut) nmu ++;
	   }
		
	}

	Bool_t ok = (nmu >= 2 ) ;
	return ok;

}

Bool_t L1Menu2015::Onia(Float_t ptcut1, Float_t ptcut2, Float_t etacut, Int_t delta) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t Nmu = myEvt_.Nmu;
	Int_t n1=0;
	Int_t n2=0;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);		
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);		
		if (fabs(eta) > etacut) continue;
		if (pt >= ptcut1) n1 ++;
		if (pt >= ptcut2) n2++;
	}

	Bool_t ok = (n1 >=1 && n2 >= 2 ) ;
	if (! ok) return false;

	// -- now the CORRELATION condition
	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > etacut) continue;
		if (pt < ptcut1) continue;
		Int_t ieta1 = etaMuIdx(eta);

		for (Int_t imu2=0; imu2 < Nmu; imu2++) {
			if (imu2 == imu) continue;
			Int_t bx2 = myEvt_.Bxmu.at(imu2);		
			if (bx2 != 0) continue;
			Float_t pt2 = myEvt_.Ptmu.at(imu2);			
			Int_t qual2 = myEvt_.Qualmu.at(imu2);        
			if ( qual2 < 4) continue;
			Float_t eta2 = myEvt_.Etamu.at(imu2);        
			if (fabs(eta2) > etacut) continue;
			if (pt2 < ptcut2) continue;
			Int_t ieta2 = etaMuIdx(eta2);

			Float_t deta = ieta1 - ieta2; 
		// std::cout << "eta 1 2 delta " << ieta1 << " " << ieta2 << " " << deta << std::endl;
			if ( fabs(deta) <= delta)  CORREL = true;
		// if (fabs ( eta - eta2) <=  1.7) CORREL = true; 
		}

	}

	return CORREL;

}

Bool_t L1Menu2015::DoubleMu(Float_t cut1, Float_t cut2, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;  

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
	   if ( qual < qualmin) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,cut1,5.0)      ) n1++;
 		} else {
		   if (pt >= cut1) n1++;
	   } 
		
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,cut2,5.0)      ) n2++;
 		} else {
		   if (pt >= cut2) n2++;
	   } 		
	}

	Bool_t ok = (n1 >= 1 && n2 >= 2 );
	return ok;

}

Bool_t L1Menu2015::TripleMu(Float_t cut1, Float_t cut2, Float_t cut3, Int_t qualmin) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < qualmin) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,cut1,5.0)      ) n1++;
 		} else {
		   if (pt >= cut1) n1++;
	   } 
		
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,cut2,5.0)      ) n2++;
 		} else {
		   if (pt >= cut2) n2++;
	   } 
		
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,cut3,5.0)      ) n3++;
 		} else {
		   if (pt >= cut3) n3++;
	   } 										
	}

	Bool_t ok = ( n1 >= 1 && n2 >= 2 && n3 >= 3 );
	return ok;

}

Bool_t L1Menu2015::DoubleMuXOpen(Float_t cut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( (qual >= 5 || qual == 3 ) && pt >= cut ) n1 ++;
		if ( pt >= 0 ) n2 ++;
	}

	Bool_t ok = ( n1 >= 1 && n2 >= 2 );
	return ok;
}


Bool_t L1Menu2015::EvalMenu(double lumiWeight) {

	insert_ibin = 0;
	

	if(Prescales["L1_SingleEG"]>0)        InsertInMenu("L1_SingleEG",         SingleEG_Eta(trigParList["L1_SingleEG"].primTh,trigParList["L1_SingleEG"].etaCut) );
	if(Prescales["L1_SingleIsoEG"]>0)     InsertInMenu("L1_SingleIsoEG",      SingleIsoEG_Eta(trigParList["L1_SingleIsoEG"].primTh,trigParList["L1_SingleIsoEG"].etaCut) );
	if(Prescales["L1_SingleMu"]>0)        InsertInMenu("L1_SingleMu",         SingleMuEta(trigParList["L1_SingleMu"].primTh,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual));
	if(Prescales["L1_SingleIsoMu"]>0)     InsertInMenu("L1_SingleIsoMu",      SingleIsoMuEta(trigParList["L1_SingleIsoMu"].primTh,trigParList["L1_SingleIsoMu"].etaCut,trigParList["L1_SingleIsoMu"].minQual));
   if(Prescales["L1_SingleTau"]>0)       InsertInMenu("L1_SingleTau",        SingleTauJet(trigParList["L1_SingleTau"].primTh, trigParList["L1_SingleTau"].etaCut) );	   	   
   if(Prescales["L1_SingleIsoTau"]>0)    InsertInMenu("L1_SingleIsoTau",     SingleIsoTauJet(trigParList["L1_SingleIsoTau"].primTh, trigParList["L1_SingleIsoTau"].etaCut) );	   	   


	if(Prescales["L1_DoubleEG"]>0)        InsertInMenu("L1_DoubleEG",         DoubleEG(trigParList["L1_DoubleEG"].primTh,trigParList["L1_DoubleEG"].secTh) );
	if(Prescales["L1_isoEG_EG"]>0)        InsertInMenu("L1_isoEG_EG",         isoEG_EG(trigParList["L1_isoEG_EG"].primTh,trigParList["L1_isoEG_EG"].secTh) );
	if(Prescales["L1_EG_isoEG"]>0)        InsertInMenu("L1_EG_isoEG",         EG_isoEG(trigParList["L1_EG_isoEG"].primTh,trigParList["L1_EG_isoEG"].secTh) );
	if(Prescales["L1_DoubleIsoEG"]>0)     InsertInMenu("L1_DoubleIsoEG",      DoubleIsoEG(trigParList["L1_DoubleIsoEG"].primTh,trigParList["L1_DoubleIsoEG"].secTh) );

	if(Prescales["L1_DoubleMu"]>0)         InsertInMenu("L1_DoubleMu",         DoubleMu(trigParList["L1_DoubleMu"].primTh,trigParList["L1_DoubleMu"].secTh,trigParList["L1_DoubleMu"].minQual));
	if(Prescales["L1_isoMu_Mu"]>0)         InsertInMenu("L1_isoMu_Mu",         DoubleMu(trigParList["L1_isoMu_Mu"].primTh,trigParList["L1_isoMu_Mu"].secTh,trigParList["L1_isoMu_Mu"].minQual));
	if(Prescales["L1_DoubleIsoMu"]>0)      InsertInMenu("L1_DoubleIsoMu",      DoubleMu(trigParList["L1_DoubleIsoMu"].primTh,trigParList["L1_DoubleIsoMu"].secTh,trigParList["L1_DoubleIsoMu"].minQual));


	if(Prescales["L1_DoubleTau"]>0)         InsertInMenu("L1_DoubleTau",         DoubleTauJetEta(trigParList["L1_DoubleTau"].primTh,trigParList["L1_DoubleTau"].secTh,trigParList["L1_DoubleTau"].etaCut) );
	if(Prescales["L1_isoTau_Tau"]>0)        InsertInMenu("L1_isoTau_Tau",        isoTau_Tau(trigParList["L1_isoTau_Tau"].primTh,trigParList["L1_isoTau_Tau"].secTh, trigParList["L1_isoTau_Tau"].etaCut) );
	if(Prescales["L1_DoubleIsoTau"]>0)      InsertInMenu("L1_DoubleIsoTau",      DoubleIsoTau(trigParList["L1_DoubleIsoTau"].primTh,trigParList["L1_DoubleIsoTau"].secTh,trigParList["L1_DoubleIsoTau"].etaCut) );


	if(Prescales["L1_EG_Mu"]>0)        InsertInMenu("L1_EG_Mu",         EG_Mu(trigParList["L1_EG_Mu"].primTh,trigParList["L1_EG_Mu"].secTh, trigParList["L1_EG_Mu"].minQual ) );
	if(Prescales["L1_isoEG_Mu"]>0)     InsertInMenu("L1_isoEG_Mu",      isoEG_Mu(trigParList["L1_isoEG_Mu"].primTh,trigParList["L1_isoEG_Mu"].secTh, trigParList["L1_isoEG_Mu"].minQual ) );
	if(Prescales["L1_isoEG_isoMu"]>0)  InsertInMenu("L1_isoEG_isoMu",   isoEG_Mu(trigParList["L1_isoEG_isoMu"].primTh,trigParList["L1_isoEG_isoMu"].secTh, trigParList["L1_isoEG_isoMu"].minQual ) );


	if(Prescales["L1_Mu_EG"]>0)        InsertInMenu("L1_Mu_EG",         Mu_EG(trigParList["L1_Mu_EG"].primTh,trigParList["L1_Mu_EG"].secTh, trigParList["L1_Mu_EG"].minQual) );
	if(Prescales["L1_isoMu_EG"]>0)     InsertInMenu("L1_isoMu_EG",      Mu_EG(trigParList["L1_isoMu_EG"].primTh,trigParList["L1_isoMu_EG"].secTh, trigParList["L1_isoMu_EG"].minQual) );
	if(Prescales["L1_isoMu_isoEG"]>0)  InsertInMenu("L1_isoMu_isoEG",   Mu_isoEG(trigParList["L1_isoMu_isoEG"].primTh,trigParList["L1_isoMu_isoEG"].secTh, trigParList["L1_isoMu_isoEG"].minQual) );


	if(Prescales["L1_EG_Tau"]>0)    InsertInMenu("L1_EG_Tau",              EG_Tau(trigParList["L1_EG_Tau"].primTh,trigParList["L1_EG_Tau"].secTh, trigParList["L1_EG_Tau"].etaCut) );
	if(Prescales["L1_isoEG_Tau"]>0) InsertInMenu("L1_isoEG_Tau",           isoEG_Tau(trigParList["L1_isoEG_Tau"].primTh,trigParList["L1_isoEG_Tau"].secTh, trigParList["L1_isoEG_Tau"].etaCut) );
	if(Prescales["L1_isoEG_isoTau"]>0) InsertInMenu("L1_isoEG_isoTau",     isoEG_isoTau(trigParList["L1_isoEG_isoTau"].primTh,trigParList["L1_isoEG_isoTau"].secTh, trigParList["L1_isoEG_isoTau"].etaCut) );


	if(Prescales["L1_Mu_Tau"]>0)    InsertInMenu("L1_Mu_Tau",              Mu_Tau(trigParList["L1_Mu_Tau"].primTh,trigParList["L1_Mu_Tau"].secTh, trigParList["L1_Mu_Tau"].etaCut, trigParList["L1_Mu_Tau"].minQual) );
	if(Prescales["L1_isoMu_Tau"]>0) InsertInMenu("L1_isoMu_Tau",           Mu_Tau(trigParList["L1_isoMu_Tau"].primTh,trigParList["L1_isoMu_Tau"].secTh, trigParList["L1_isoMu_Tau"].etaCut, trigParList["L1_isoMu_Tau"].minQual) );
	if(Prescales["L1_isoMu_isoTau"]>0) InsertInMenu("L1_isoMu_isoTau",     Mu_isoTau(trigParList["L1_isoMu_isoTau"].primTh,trigParList["L1_isoMu_isoTau"].secTh, trigParList["L1_isoMu_isoTau"].etaCut, trigParList["L1_isoMu_isoTau"].minQual) );


	if(Prescales["L1_SingleJet"]>0) InsertInMenu("L1_SingleJet",      SingleJet(trigParList["L1_SingleJet"].primTh) );
	if(Prescales["L1_SingleJetC"]>0)InsertInMenu("L1_SingleJetC",     SingleJetCentral(trigParList["L1_SingleJetC"].primTh) );
   if(Prescales["L1_DoubleJet"]>0)    InsertInMenu("L1_DoubleJet",      DoubleJetCentral(trigParList["L1_DoubleJet"].primTh,trigParList["L1_DoubleJet"].secTh) );
	if(Prescales["L1_QuadJetC"]>0)     InsertInMenu("L1_QuadJetC",       QuadJetCentral(trigParList["L1_QuadJetC"].primTh,trigParList["L1_QuadJetC"].secTh,trigParList["L1_QuadJetC"].triTh,trigParList["L1_QuadJetC"].quadTh) );
	if(Prescales["L1_SixJet"]>0)       InsertInMenu("L1_SixJet",         MultiJet(6,trigParList["L1_SixJet"].primTh,trigParList["L1_SixJet"].secTh,trigParList["L1_SixJet"].triTh,trigParList["L1_SixJet"].quadTh) );


	if(Prescales["L1_SingleEG_CJet"]>0)   InsertInMenu("L1_SingleEG_CJet",    EG_JetCentral(trigParList["L1_SingleEG_CJet"].primTh,trigParList["L1_SingleEG_CJet"].secTh) );	   
	if(Prescales["L1_SingleIsoEG_CJet"]>0)InsertInMenu("L1_SingleIsoEG_CJet", IsoEG_JetCentral(trigParList["L1_SingleIsoEG_CJet"].primTh,trigParList["L1_SingleIsoEG_CJet"].secTh,trigParList["L1_SingleIsoEG_CJet"].etaCut) );	   
	if(Prescales["L1_SingleMu_CJet"]>0)    InsertInMenu("L1_SingleMu_CJet",    Muer_JetCentral(trigParList["L1_SingleMu_CJet"].primTh,trigParList["L1_SingleMu_CJet"].secTh,trigParList["L1_SingleMu_CJet"].etaCut, trigParList["L1_SingleMu_CJet"].minQual) );
	if(Prescales["L1_SingleIsoMu_CJet"]>0) InsertInMenu("L1_SingleIsoMu_CJet", IsoMuer_JetCentral(trigParList["L1_SingleIsoMu_CJet"].primTh,trigParList["L1_SingleIsoMu_CJet"].secTh,trigParList["L1_SingleIsoMu_CJet"].etaCut, trigParList["L1_SingleIsoMu_CJet"].minQual) );
// boolean methods for these need to be fixed.
//	if(Prescales["L1_SingleTau_CJet"]>0)    InsertInMenu("L1_SingleTau_CJet",    Tau_JetCentral(trigParList["L1_SingleTau_CJet"].primTh,trigParList["L1_SingleTau_CJet"].secTh,trigParList["L1_SingleTau_CJet"].etaCut) );	   
//	if(Prescales["L1_SingleIsoTau_CJet"]>0)    InsertInMenu("L1_SingleIsoTau_CJet",    IsoTau_JetCentral(trigParList["L1_SingleIsoTau_CJet"].primTh,trigParList["L1_SingleIsoTau_CJet"].secTh,trigParList["L1_SingleIsoTau_CJet"].etaCut) );	   
	if(Prescales["L1_SingleTau_TwoFJet"]>0) InsertInMenu("L1_SingleTau_TwoFJet", Tau_TwoJetForward(trigParList["L1_SingleTau_TwoFJet"].primTh,trigParList["L1_SingleTau_TwoFJet"].secTh,trigParList["L1_SingleTau_TwoFJet"].triTh) );	   
   if(Prescales["L1_DoubleFwdJet"]>0) InsertInMenu("L1_DoubleFwdJet",   DoubleJetForward(trigParList["L1_DoubleFwdJet"].primTh,trigParList["L1_DoubleFwdJet"].secTh) );



	if(Prescales["L1_SingleEG_ETM"]>0)     InsertInMenu("L1_SingleEG_ETM",     EG_ETM(trigParList["L1_SingleEG_ETM"].primTh,trigParList["L1_SingleEG_ETM"].secTh) );
	if(Prescales["L1_SingleIsoEG_ETM"]>0)  InsertInMenu("L1_SingleIsoEG_ETM",  IsoEG_ETM(trigParList["L1_SingleIsoEG_ETM"].primTh,trigParList["L1_SingleIsoEG_ETM"].secTh,trigParList["L1_SingleIsoEG_ETM"].etaCut) );
	if(Prescales["L1_SingleMu_ETM"]>0)     InsertInMenu("L1_SingleMu_ETM",     Muer_ETM(trigParList["L1_SingleMu_ETM"].primTh,trigParList["L1_SingleMu_ETM"].secTh,trigParList["L1_SingleMu_ETM"].etaCut, trigParList["L1_SingleMu_ETM"].minQual) );
	if(Prescales["L1_SingleIsoMu_ETM"]>0)  InsertInMenu("L1_SingleIsoMu_ETM",  IsoMuer_ETM(trigParList["L1_SingleIsoMu_ETM"].primTh,trigParList["L1_SingleIsoMu_ETM"].secTh,trigParList["L1_SingleIsoMu_ETM"].etaCut, trigParList["L1_SingleIsoMu_ETM"].minQual) );
	if(Prescales["L1_SingleTau_ETM"]>0)    InsertInMenu("L1_SingleTau_ETM",    Tau_ETM(trigParList["L1_SingleTau_ETM"].primTh,trigParList["L1_SingleTau_ETM"].secTh,trigParList["L1_SingleTau_ETM"].etaCut) );	   
	if(Prescales["L1_SingleIsoTau_ETM"]>0) InsertInMenu("L1_SingleIsoTau_ETM", IsoTau_ETM(trigParList["L1_SingleIsoTau_ETM"].primTh,trigParList["L1_SingleIsoTau_ETM"].secTh,trigParList["L1_SingleIsoTau_ETM"].etaCut) );	   

	if(Prescales["L1_SingleEG_HTM"]>0)     InsertInMenu("L1_SingleEG_HTM",     EG_HTM(trigParList["L1_SingleEG_HTM"].primTh,trigParList["L1_SingleEG_HTM"].secTh) );
	if(Prescales["L1_SingleIsoEG_HTM"]>0)  InsertInMenu("L1_SingleIsoEG_HTM",  IsoEG_HTM(trigParList["L1_SingleIsoEG_HTM"].primTh,trigParList["L1_SingleIsoEG_HTM"].secTh,trigParList["L1_SingleIsoEG_HTM"].etaCut) );
	if(Prescales["L1_SingleMu_HTM"]>0)     InsertInMenu("L1_SingleMu_HTM",     Muer_HTM(trigParList["L1_SingleMu_HTM"].primTh,trigParList["L1_SingleMu_HTM"].secTh,trigParList["L1_SingleMu_HTM"].etaCut, trigParList["L1_SingleMu_HTM"].minQual) );
	if(Prescales["L1_SingleIsoMu_HTM"]>0)  InsertInMenu("L1_SingleIsoMu_HTM",  IsoMuer_HTM(trigParList["L1_SingleIsoMu_HTM"].primTh,trigParList["L1_SingleIsoMu_HTM"].secTh,trigParList["L1_SingleIsoMu_HTM"].etaCut, trigParList["L1_SingleIsoMu_HTM"].minQual) );
	if(Prescales["L1_SingleTau_HTM"]>0)    InsertInMenu("L1_SingleTau_HTM",    Tau_HTM(trigParList["L1_SingleTau_HTM"].primTh,trigParList["L1_SingleTau_HTM"].secTh,trigParList["L1_SingleTau_HTM"].etaCut) );	   
	if(Prescales["L1_SingleIsoTau_HTM"]>0) InsertInMenu("L1_SingleIsoTau_HTM", IsoTau_HTM(trigParList["L1_SingleIsoTau_HTM"].primTh,trigParList["L1_SingleIsoTau_HTM"].secTh,trigParList["L1_SingleIsoTau_HTM"].etaCut) );	   



	if(Prescales["L1_HTM"]>0) InsertInMenu("L1_HTM",           HTM(trigParList["L1_HTM"].primTh) );    
	if(Prescales["L1_ETM"]>0) InsertInMenu("L1_ETM",           ETM(trigParList["L1_ETM"].primTh) );           
	if(Prescales["L1_HTT"]>0) InsertInMenu("L1_HTT",           HTT(trigParList["L1_HTT"].primTh) );


	Int_t NN = insert_ibin;

	Int_t kOFFSET_old = kOFFSET;
	for (Int_t k=0; k < NN; k++) {
		TheTriggerBits[k + kOFFSET_old] = insert_val[k];
	}
	kOFFSET += insert_ibin;

	if (first) {

		NBITS_TRIGS = NN;

		for (Int_t ibin=0; ibin < insert_ibin; ibin++) {
			TString l1name = (TString)insert_names[ibin];
			h_Trig -> GetXaxis() -> SetBinLabel(ibin+1, l1name );
		}
		h_Trig -> GetXaxis() -> SetBinLabel(NN+1, "Triggered") ;

		for (Int_t k=1; k <= kOFFSET - kOFFSET_old ; k++) {
			h_All -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
			h_Cumm -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
		   h_Corr -> GetXaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
			h_Corr -> GetYaxis() -> SetBinLabel(k +kOFFSET_old , h_Trig -> GetXaxis() -> GetBinLabel(k) );
		}
	}

	Bool_t res = false;
	for (Int_t i=0; i < NN; i++) {
		res = res || insert_val[i] ;
		if (insert_val[i]) h_Trig -> Fill(i,lumiWeight);
	}
	if (res) h_Trig -> Fill(NN,lumiWeight);

	return res;
}

void L1Menu2015::EvalThresh(Int_t calcThreshold, double lumiWeight) {
 

/* notes:

  1) Histograms are binned so that integer values are at the bin centers.
  2) For EG and Muon, the bin widths are 1 GeV (Not ideal for muons which have a non uniform distribution of values)
  3) For Jets, the bin widths are 4 GeV
  4) For EG and Jets we use the lower edge of the bin to define the threshold and avoid rounding errors
         e.g For EG:  Lower edge is 15.5 GeV so test will be whether EG > 15.5 but then this will be plotted
			             at bin center which is 16 GeV.  This avoids cases where numbers like  15.99999999999999   
							 would fail threshold or be plotted in the incorrect bin.
			For jets it is similar excepts bins step by 4 GeV.
  5) For Muons, we take the threholds from the bin centers and live with small rnding problems since muon pts are nonuniform.
  6) For non-primary thresholds, we scale by ratios. These ratio are multiplied by the bin_center values since they 
     correspond to the true threshold.
	      -> No attempt is currently made to adjust the non-primary thresholds so that they find actual integer values (or factors of 4 in the case of jets) 
  
  
*/
	//------- SingleMu -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu = h_SingleMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu+1; bin++){
//		const float bin_low_edge = h_SingleMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_byThreshold->GetBinCenter(bin);
		if(SingleMuEta(bin_center,trigParList["L1_SingleMu"].etaCut,trigParList["L1_SingleMu"].minQual)){
			h_SingleMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------- SingleIsoMu -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoMu = h_SingleIsoMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoMu+1; bin++){
//		const float bin_low_edge = h_SingleIsoMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoMu_byThreshold->GetBinCenter(bin);
		if(SingleIsoMuEta(bin_center,trigParList["L1_SingleIsoMu"].etaCut,trigParList["L1_SingleIsoMu"].minQual)){
			h_SingleIsoMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------



	//------ DoubleMu ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleMu = h_DoubleMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleMu+1; bin++){
//		const float bin_low_edge = h_DoubleMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleMu_byThreshold->GetBinCenter(bin);
		if(DoubleMu(bin_center,bin_center*(trigParList["L1_DoubleMu"].secTh/trigParList["L1_DoubleMu"].primTh),trigParList["L1_DoubleMu"].minQual )){
			h_DoubleMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleMu ---- 2-D Evaluation --Only fill half of space because symmetric (if id cuts identical) --------------
	if(calcThreshold>1) {
	  const unsigned n_bins_DoubleMu_X = h2_DoubleMu_byThreshold->GetNbinsX();
//	  const unsigned n_bins_DoubleMu_Y = h2_DoubleMu_byThreshold->GetNbinsY();
	  for(unsigned bin=1; bin <= n_bins_DoubleMu_X+1; bin++){	
//	     const float bin_low_edge_X = h2_DoubleMu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
		  const float bin_center_X = h2_DoubleMu_byThreshold->GetXaxis()->GetBinCenter(bin);
		  for(unsigned ybin=1; ybin <= bin; ybin++){	  //stop the y-axis scan at the diagonal since this is symmetric
//		     const float bin_low_edge_Y = h2_DoubleMu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		     const float bin_center_Y = h2_DoubleMu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		     if(DoubleMu(bin_center_X, bin_center_Y,trigParList["L1_DoubleMu"].minQual)){
			     h2_DoubleMu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		     }
		  }	
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ DoubleIsoMu ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleIsoMu = h_DoubleIsoMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleIsoMu+1; bin++){
//		const float bin_low_edge = h_DoubleIsoMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleIsoMu_byThreshold->GetBinCenter(bin);
		if(DoubleMu(bin_center,bin_center*(trigParList["L1_DoubleIsoMu"].secTh/trigParList["L1_DoubleIsoMu"].primTh),trigParList["L1_DoubleIsoMu"].minQual )){
			h_DoubleIsoMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleIsoMu ---- 2-D Evaluation --Only fill half of space because symmetric (if id cuts identical) --------------
	if(calcThreshold>1) {
	  const unsigned n_bins_DoubleIsoMu_X = h2_DoubleIsoMu_byThreshold->GetNbinsX();
//	  const unsigned n_bins_DoubleIsoMu_Y = h2_DoubleIsoMu_byThreshold->GetNbinsY();
	  for(unsigned bin=1; bin <= n_bins_DoubleIsoMu_X+1; bin++){	
//	     const float bin_low_edge_X = h2_DoubleIsoMu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
		  const float bin_center_X = h2_DoubleIsoMu_byThreshold->GetXaxis()->GetBinCenter(bin);
		  for(unsigned ybin=1; ybin <= bin; ybin++){	  //stop the y-axis scan at the diagonal since this is symmetric
//		     const float bin_low_edge_Y = h2_DoubleIsoMu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		     const float bin_center_Y = h2_DoubleIsoMu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		     if(DoubleMu(bin_center_X, bin_center_Y,trigParList["L1_DoubleIsoMu"].minQual)){
			     h2_DoubleIsoMu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		     }
		  }	
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------



	//------ isoMu_Mu ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoMu_Mu = h_isoMu_Mu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoMu_Mu+1; bin++){
//		const float bin_low_edge = h_isoMu_Mu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoMu_Mu_byThreshold->GetBinCenter(bin);
		if(DoubleMu(bin_center,bin_center*(trigParList["L1_isoMu_Mu"].secTh/trigParList["L1_isoMu_Mu"].primTh),trigParList["L1_isoMu_Mu"].minQual )){
			h_isoMu_Mu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ isoMu_Mu ---- 2-D Evaluation --Only fill half of space because symmetric (if id cuts identical) --------------
	if(calcThreshold>1) {
	  const unsigned n_bins_isoMu_Mu_X = h2_isoMu_Mu_byThreshold->GetNbinsX();
//	  const unsigned n_bins_isoMu_Mu_Y = h2_isoMu_Mu_byThreshold->GetNbinsY();
	  for(unsigned bin=1; bin <= n_bins_isoMu_Mu_X+1; bin++){	
//	     const float bin_low_edge_X = h2_isoMu_Mu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
		  const float bin_center_X = h2_isoMu_Mu_byThreshold->GetXaxis()->GetBinCenter(bin);
		  for(unsigned ybin=1; ybin <= bin; ybin++){	  //stop the y-axis scan at the diagonal since this is symmetric
//		     const float bin_low_edge_Y = h2_isoMu_Mu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		     const float bin_center_Y = h2_isoMu_Mu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		     if(DoubleMu(bin_center_X, bin_center_Y,trigParList["L1_isoMu_Mu"].minQual)){
			     h2_isoMu_Mu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		     }
		  }	
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------- SingleMu_ETM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu_ETM = h_SingleMu_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu_ETM+1; bin++){
//		const float bin_low_edge = h_SingleMu_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_ETM_byThreshold->GetBinCenter(bin);
		if(Muer_ETM(bin_center,bin_center*(trigParList["L1_SingleMu_ETM"].secTh/trigParList["L1_SingleMu_ETM"].primTh),trigParList["L1_SingleMu_ETM"].etaCut,trigParList["L1_SingleMu_ETM"].minQual)){
			h_SingleMu_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleMu_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
   if(calcThreshold>1) {
		const unsigned n_bins_SingleMu_ETM_X = h2_SingleMu_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleMu_ETM_Y = h2_SingleMu_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleMu_ETM_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleMu_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleMu_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleMu_ETM_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleMu_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleMu_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Muer_ETM(bin_center_X, bin_center_Y,trigParList["L1_SingleMu_ETM"].etaCut,trigParList["L1_SingleMu_ETM"].minQual)){
			   	h2_SingleMu_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------- SingleIsoMu_ETM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoMu_ETM = h_SingleIsoMu_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoMu_ETM+1; bin++){
//		const float bin_low_edge = h_SingleIsoMu_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoMu_ETM_byThreshold->GetBinCenter(bin);
		if(IsoMuer_ETM(bin_center,bin_center*(trigParList["L1_SingleIsoMu_ETM"].secTh/trigParList["L1_SingleIsoMu_ETM"].primTh),trigParList["L1_SingleIsoMu_ETM"].etaCut,trigParList["L1_SingleIsoMu_ETM"].minQual)){
			h_SingleIsoMu_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoMu_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
   if(calcThreshold>1) {
		const unsigned n_bins_SingleIsoMu_ETM_X = h2_SingleIsoMu_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoMu_ETM_Y = h2_SingleIsoMu_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoMu_ETM_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleIsoMu_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoMu_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoMu_ETM_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleIsoMu_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoMu_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoMuer_ETM(bin_center_X, bin_center_Y,trigParList["L1_SingleIsoMu_ETM"].etaCut,trigParList["L1_SingleIsoMu_ETM"].minQual)){
			   	h2_SingleIsoMu_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------
	

	//------- SingleMu_HTM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu_HTM = h_SingleMu_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu_HTM+1; bin++){
//		const float bin_low_edge = h_SingleMu_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_HTM_byThreshold->GetBinCenter(bin);
		if(Muer_HTM(bin_center,bin_center*(trigParList["L1_SingleMu_HTM"].secTh/trigParList["L1_SingleMu_HTM"].primTh),trigParList["L1_SingleMu_HTM"].etaCut,trigParList["L1_SingleMu_HTM"].minQual)){
			h_SingleMu_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleMu_HTM ---- 2-D Evaluation ------------------------------------------------------------------------------
   if(calcThreshold>1) {
		const unsigned n_bins_SingleMu_HTM_X = h2_SingleMu_HTM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleMu_HTM_Y = h2_SingleMu_HTM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleMu_HTM_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleMu_HTM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleMu_HTM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleMu_HTM_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleMu_HTM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleMu_HTM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Muer_HTM(bin_center_X, bin_center_Y,trigParList["L1_SingleMu_HTM"].etaCut,trigParList["L1_SingleMu_HTM"].minQual)){
			   	h2_SingleMu_HTM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------
	

	//------- SingleIsoMu_HTM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoMu_HTM = h_SingleIsoMu_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoMu_HTM+1; bin++){
//		const float bin_low_edge = h_SingleIsoMu_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoMu_HTM_byThreshold->GetBinCenter(bin);
		if(IsoMuer_HTM(bin_center,bin_center*(trigParList["L1_SingleIsoMu_HTM"].secTh/trigParList["L1_SingleIsoMu_HTM"].primTh),trigParList["L1_SingleIsoMu_HTM"].etaCut,trigParList["L1_SingleIsoMu_HTM"].minQual)){
			h_SingleIsoMu_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoMu_HTM ---- 2-D Evaluation ------------------------------------------------------------------------------
   if(calcThreshold>1) {
		const unsigned n_bins_SingleIsoMu_HTM_X = h2_SingleIsoMu_HTM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoMu_HTM_Y = h2_SingleIsoMu_HTM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoMu_HTM_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleIsoMu_HTM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoMu_HTM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoMu_HTM_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleIsoMu_HTM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoMu_HTM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoMuer_HTM(bin_center_X, bin_center_Y,trigParList["L1_SingleIsoMu_HTM"].etaCut,trigParList["L1_SingleIsoMu_HTM"].minQual)){
			   	h2_SingleIsoMu_HTM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	
	//------- SingleMu_CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleMu_CJet = h_SingleMu_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleMu_CJet+1; bin++){
//		const float bin_low_edge = h_SingleMu_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleMu_CJet_byThreshold->GetBinCenter(bin);
		if(Muer_JetCentral(bin_center,bin_center*(trigParList["L1_SingleMu_CJet"].secTh/trigParList["L1_SingleMu_CJet"].primTh),trigParList["L1_SingleMu_CJet"].etaCut,trigParList["L1_SingleMu_CJet"].minQual)){
			h_SingleMu_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleMu_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleMu_CJet_X = h2_SingleMu_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleMu_CJet_Y = h2_SingleMu_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleMu_CJet_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleMu_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleMu_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleMu_CJet_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleMu_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleMu_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Muer_JetCentral(bin_center_X, bin_center_Y,trigParList["L1_SingleMu_CJet"].etaCut,trigParList["L1_SingleMu_CJet"].minQual)){
			   	h2_SingleMu_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------
	
	
	//------- SingleIsoMu_CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoMu_CJet = h_SingleIsoMu_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoMu_CJet+1; bin++){
//		const float bin_low_edge = h_SingleIsoMu_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoMu_CJet_byThreshold->GetBinCenter(bin);
		if(IsoMuer_JetCentral(bin_center,bin_center*(trigParList["L1_SingleIsoMu_CJet"].secTh/trigParList["L1_SingleIsoMu_CJet"].primTh),trigParList["L1_SingleIsoMu_CJet"].etaCut,trigParList["L1_SingleIsoMu_CJet"].minQual)){
			h_SingleIsoMu_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoMu_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleIsoMu_CJet_X = h2_SingleIsoMu_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoMu_CJet_Y = h2_SingleIsoMu_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoMu_CJet_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_SingleIsoMu_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoMu_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoMu_CJet_Y+1; ybin++){	
//		   	const float bin_low_edge_Y = h2_SingleIsoMu_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoMu_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoMuer_JetCentral(bin_center_X, bin_center_Y,trigParList["L1_SingleIsoMu_CJet"].etaCut,trigParList["L1_SingleIsoMu_CJet"].minQual)){
			   	h2_SingleIsoMu_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------	


	//--------- SingleEG ---------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG = h_SingleEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG+1; bin++){
		const float bin_low_edge = h_SingleEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_byThreshold->GetBinCenter(bin);
		if(SingleEG_Eta(bin_low_edge,trigParList["L1_SingleEG"].etaCut)){
			h_SingleEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	
	//--------- SingleIsoEG ------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoEG = h_SingleIsoEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoEG+1; bin++){
		const float bin_low_edge = h_SingleIsoEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoEG_byThreshold->GetBinCenter(bin);
		if(SingleIsoEG_Eta(bin_low_edge,trigParList["L1_SingleIsoEG"].etaCut)){
			h_SingleIsoEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	

	//------ DoubleEG ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleEG = h_DoubleEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleEG+1; bin++){
		const float bin_low_edge = h_DoubleEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleEG_byThreshold->GetBinCenter(bin);
		if(DoubleEG(bin_low_edge, bin_center*(trigParList["L1_DoubleEG"].secTh/trigParList["L1_DoubleEG"].primTh))){
			h_DoubleEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleEG ---- 2-D Evaluation ----(Fill only bottom half since symmetric (if ID cuts same) -------------------
	if(calcThreshold>1) {
		const unsigned n_bins_DoubleEG_X = h2_DoubleEG_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleEG_Y = h2_DoubleEG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleEG_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleEG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleEG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	         //stop the y-axis scan at the diagonal since this is symmetric
		   	const float bin_low_edge_Y = h2_DoubleEG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleEG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleEG(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_DoubleEG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleIsoEG ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleIsoEG = h_DoubleIsoEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleIsoEG+1; bin++){
		const float bin_low_edge = h_DoubleIsoEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleIsoEG_byThreshold->GetBinCenter(bin);
		if(DoubleIsoEG(bin_low_edge, bin_center*(trigParList["L1_DoubleIsoEG"].secTh/trigParList["L1_DoubleIsoEG"].primTh))){
			h_DoubleIsoEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleIsoEG ---- 2-D Evaluation ----(Fill only bottom half since symmetric (if ID cuts same) -------------------
	if(calcThreshold>1) {
		const unsigned n_bins_DoubleIsoEG_X = h2_DoubleIsoEG_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleIsoEG_Y = h2_DoubleIsoEG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleIsoEG_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleIsoEG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleIsoEG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	         //stop the y-axis scan at the diagonal since this is symmetric
		   	const float bin_low_edge_Y = h2_DoubleIsoEG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleIsoEG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleIsoEG(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_DoubleIsoEG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//------ isoEG_EG ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoEG_EG = h_isoEG_EG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoEG_EG+1; bin++){
		const float bin_low_edge = h_isoEG_EG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoEG_EG_byThreshold->GetBinCenter(bin);
		if(isoEG_EG(bin_low_edge, bin_center*(trigParList["L1_isoEG_EG"].secTh/trigParList["L1_isoEG_EG"].primTh))){
			h_isoEG_EG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ EG_isoEG ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_EG_isoEG = h_EG_isoEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_EG_isoEG+1; bin++){
		const float bin_low_edge = h_EG_isoEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_EG_isoEG_byThreshold->GetBinCenter(bin);
		if(EG_isoEG(bin_low_edge, bin_center*(trigParList["L1_EG_isoEG"].secTh/trigParList["L1_EG_isoEG"].primTh))){
			h_EG_isoEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleEG_ETM --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG_ETM = h_SingleEG_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG_ETM+1; bin++){
		const float bin_low_edge = h_SingleEG_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_ETM_byThreshold->GetBinCenter(bin);
		if(EG_ETM(bin_low_edge, bin_center*(trigParList["L1_SingleEG_ETM"].secTh/trigParList["L1_SingleEG_ETM"].primTh))){
			h_SingleEG_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleEG_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleEG_ETM_X = h2_SingleEG_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleEG_ETM_Y = h2_SingleEG_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleEG_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleEG_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleEG_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleEG_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleEG_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleEG_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_ETM(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_SingleEG_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleIsoEG_ETM --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoEG_ETM = h_SingleIsoEG_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoEG_ETM+1; bin++){
		const float bin_low_edge = h_SingleIsoEG_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoEG_ETM_byThreshold->GetBinCenter(bin);
		if(IsoEG_ETM(bin_low_edge, bin_center*(trigParList["L1_SingleIsoEG_ETM"].secTh/trigParList["L1_SingleIsoEG_ETM"].primTh),trigParList["L1_SingleIsoEG_ETM"].etaCut)){
			h_SingleIsoEG_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoEG_ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleIsoEG_ETM_X = h2_SingleIsoEG_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoEG_ETM_Y = h2_SingleIsoEG_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoEG_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleIsoEG_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoEG_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoEG_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleIsoEG_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoEG_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoEG_ETM(bin_low_edge_X, bin_low_edge_Y,trigParList["L1_SingleIsoEG_ETM"].etaCut)){
			   	h2_SingleIsoEG_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleEG_HTM --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG_HTM = h_SingleEG_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG_HTM+1; bin++){
		const float bin_low_edge = h_SingleEG_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_HTM_byThreshold->GetBinCenter(bin);
		if(EG_HTM(bin_low_edge, bin_center*(trigParList["L1_SingleEG_HTM"].secTh/trigParList["L1_SingleEG_HTM"].primTh))){
			h_SingleEG_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleEG_HTM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleEG_HTM_X = h2_SingleEG_HTM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleEG_HTM_Y = h2_SingleEG_HTM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleEG_HTM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleEG_HTM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleEG_HTM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleEG_HTM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleEG_HTM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleEG_HTM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_HTM(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_SingleEG_HTM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleIsoEG_HTM --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoEG_HTM = h_SingleIsoEG_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoEG_HTM+1; bin++){
		const float bin_low_edge = h_SingleIsoEG_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoEG_HTM_byThreshold->GetBinCenter(bin);
		if(IsoEG_HTM(bin_low_edge, bin_center*(trigParList["L1_SingleIsoEG_HTM"].secTh/trigParList["L1_SingleIsoEG_HTM"].primTh),trigParList["L1_SingleIsoEG_HTM"].etaCut)){
			h_SingleIsoEG_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoEG_HTM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {
		const unsigned n_bins_SingleIsoEG_HTM_X = h2_SingleIsoEG_HTM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoEG_HTM_Y = h2_SingleIsoEG_HTM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoEG_HTM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleIsoEG_HTM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoEG_HTM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoEG_HTM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleIsoEG_HTM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoEG_HTM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoEG_HTM(bin_low_edge_X, bin_low_edge_Y,trigParList["L1_SingleIsoEG_HTM"].etaCut)){
			   	h2_SingleIsoEG_HTM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------



	//------ SingleEG_CJet --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleEG_CJet = h_SingleEG_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleEG_CJet+1; bin++){
		const float bin_low_edge = h_SingleEG_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleEG_CJet_byThreshold->GetBinCenter(bin);
		if(EG_JetCentral(bin_low_edge, bin_center*(trigParList["L1_SingleEG_CJet"].secTh/trigParList["L1_SingleEG_CJet"].primTh))){
			h_SingleEG_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleEG_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleEG_CJet_X = h2_SingleEG_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleEG_CJet_Y = h2_SingleEG_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleEG_CJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleEG_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleEG_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleEG_CJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleEG_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleEG_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_JetCentral(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_SingleEG_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleIsoEG_CJet --------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoEG_CJet = h_SingleIsoEG_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoEG_CJet+1; bin++){
		const float bin_low_edge = h_SingleIsoEG_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoEG_CJet_byThreshold->GetBinCenter(bin);
		if(IsoEG_JetCentral(bin_low_edge, bin_center*(trigParList["L1_SingleIsoEG_CJet"].secTh/trigParList["L1_SingleIsoEG_CJet"].primTh),trigParList["L1_SingleIsoEG_CJet"].etaCut)){
			h_SingleIsoEG_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleIsoEG_CJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleIsoEG_CJet_X = h2_SingleIsoEG_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleIsoEG_CJet_Y = h2_SingleIsoEG_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleIsoEG_CJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleIsoEG_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleIsoEG_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleIsoEG_CJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleIsoEG_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleIsoEG_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(IsoEG_JetCentral(bin_low_edge_X, bin_low_edge_Y,trigParList["L1_SingleIsoEG_CJet"].etaCut)){
			   	h2_SingleIsoEG_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------



	//-------- Mu_EG -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_Mu_EG = h_Mu_EG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_Mu_EG+1; bin++){
//		const float bin_low_edge = h_Mu_EG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_Mu_EG_byThreshold->GetBinCenter(bin);
		if(Mu_EG(bin_center, bin_center*(trigParList["L1_Mu_EG"].secTh/trigParList["L1_Mu_EG"].primTh), trigParList["L1_Mu_EG"].minQual)){
			h_Mu_EG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	
	
	
	//------- EG_Mu --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_EG_Mu = h_EG_Mu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_EG_Mu+1; bin++){
		const float bin_low_edge = h_EG_Mu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_EG_Mu_byThreshold->GetBinCenter(bin);
		if(EG_Mu(bin_low_edge,  bin_center*(trigParList["L1_EG_Mu"].secTh/trigParList["L1_EG_Mu"].primTh), trigParList["L1_EG_Mu"].minQual) ){
			h_EG_Mu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ Mu_EG ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_Mu_EG_X = h2_Mu_EG_byThreshold->GetNbinsX();
		const unsigned n_bins_Mu_EG_Y = h2_Mu_EG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_Mu_EG_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_Mu_EG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_Mu_EG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_Mu_EG_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_Mu_EG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_Mu_EG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_EG(bin_center_X, bin_low_edge_Y, trigParList["L1_Mu_EG"].minQual)){
			   	h2_Mu_EG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//-------- isoMu_EG -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoMu_EG = h_isoMu_EG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoMu_EG+1; bin++){
//		const float bin_low_edge = h_isoMu_EG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoMu_EG_byThreshold->GetBinCenter(bin);
		if(Mu_EG(bin_center, bin_center*(trigParList["L1_isoMu_EG"].secTh/trigParList["L1_isoMu_EG"].primTh), trigParList["L1_isoMu_EG"].minQual)){
			h_isoMu_EG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}


	//------ isoMu_EG ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoMu_EG_X = h2_isoMu_EG_byThreshold->GetNbinsX();
		const unsigned n_bins_isoMu_EG_Y = h2_isoMu_EG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoMu_EG_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_isoMu_EG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoMu_EG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoMu_EG_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoMu_EG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoMu_EG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_EG(bin_center_X, bin_low_edge_Y, trigParList["L1_isoMu_EG"].minQual)){
			   	h2_isoMu_EG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//-------- isoMu_isoEG -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoMu_isoEG = h_isoMu_isoEG_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoMu_isoEG+1; bin++){
//		const float bin_low_edge = h_isoMu_isoEG_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoMu_isoEG_byThreshold->GetBinCenter(bin);
		if(Mu_isoEG(bin_center, bin_center*(trigParList["L1_isoMu_isoEG"].secTh/trigParList["L1_isoMu_isoEG"].primTh), trigParList["L1_isoMu_isoEG"].minQual)){
			h_isoMu_isoEG_byThreshold->Fill(bin_center,lumiWeight);
		}
	}


	//------ isoMu_isoEG ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoMu_isoEG_X = h2_isoMu_isoEG_byThreshold->GetNbinsX();
		const unsigned n_bins_isoMu_isoEG_Y = h2_isoMu_isoEG_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoMu_isoEG_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_isoMu_isoEG_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoMu_isoEG_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoMu_isoEG_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoMu_isoEG_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoMu_isoEG_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_isoEG(bin_center_X, bin_low_edge_Y, trigParList["L1_isoMu_isoEG"].minQual)){
			   	h2_isoMu_isoEG_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//------- isoEG_Mu --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoEG_Mu = h_isoEG_Mu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoEG_Mu+1; bin++){
		const float bin_low_edge = h_isoEG_Mu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoEG_Mu_byThreshold->GetBinCenter(bin);
		if(isoEG_Mu(bin_low_edge,  bin_center*(trigParList["L1_isoEG_Mu"].secTh/trigParList["L1_isoEG_Mu"].primTh), trigParList["L1_isoEG_Mu"].minQual) ){
			h_isoEG_Mu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ isoEG_Mu ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoEG_Mu_X = h2_isoEG_Mu_byThreshold->GetNbinsX();
		const unsigned n_bins_isoEG_Mu_Y = h2_isoEG_Mu_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoEG_Mu_X+1; bin++){	
	   	const float bin_low_edge_X = h2_isoEG_Mu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoEG_Mu_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoEG_Mu_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoEG_Mu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoEG_Mu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(isoEG_Mu(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_isoEG_Mu"].minQual)){
			   	h2_isoEG_Mu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------- isoEG_isoMu --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoEG_isoMu = h_isoEG_isoMu_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoEG_isoMu+1; bin++){
		const float bin_low_edge = h_isoEG_isoMu_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoEG_isoMu_byThreshold->GetBinCenter(bin);
		if(isoEG_Mu(bin_low_edge,  bin_center*(trigParList["L1_isoEG_isoMu"].secTh/trigParList["L1_isoEG_isoMu"].primTh), trigParList["L1_isoEG_isoMu"].minQual) ){
			h_isoEG_isoMu_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ isoEG_isoMu ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoEG_isoMu_X = h2_isoEG_isoMu_byThreshold->GetNbinsX();
		const unsigned n_bins_isoEG_isoMu_Y = h2_isoEG_isoMu_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoEG_isoMu_X+1; bin++){	
	   	const float bin_low_edge_X = h2_isoEG_isoMu_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoEG_isoMu_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoEG_isoMu_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoEG_isoMu_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoEG_isoMu_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(isoEG_Mu(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_isoEG_isoMu"].minQual)){
			   	h2_isoEG_isoMu_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//-------- Mu_Tau -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_Mu_Tau = h_Mu_Tau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_Mu_Tau+1; bin++){
//		const float bin_low_edge = h_Mu_Tau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_Mu_Tau_byThreshold->GetBinCenter(bin);
		if(Mu_Tau(bin_center, bin_center*(trigParList["L1_Mu_Tau"].secTh/trigParList["L1_Mu_Tau"].primTh), trigParList["L1_Mu_Tau"].etaCut, trigParList["L1_Mu_Tau"].minQual)){
			h_Mu_Tau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	
	//------ Mu_Tau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_Mu_Tau_X = h2_Mu_Tau_byThreshold->GetNbinsX();
		const unsigned n_bins_Mu_Tau_Y = h2_Mu_Tau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_Mu_Tau_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_Mu_Tau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_Mu_Tau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_Mu_Tau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_Mu_Tau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_Mu_Tau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_Tau(bin_center_X, bin_low_edge_Y, trigParList["L1_Mu_Tau"].etaCut, trigParList["L1_Mu_Tau"].minQual)){
			   	h2_Mu_Tau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------	

	//-------- isoMu_Tau -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoMu_Tau = h_isoMu_Tau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoMu_Tau+1; bin++){
//		const float bin_low_edge = h_isoMu_Tau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoMu_Tau_byThreshold->GetBinCenter(bin);
		if(Mu_Tau(bin_center, bin_center*(trigParList["L1_isoMu_Tau"].secTh/trigParList["L1_isoMu_Tau"].primTh), trigParList["L1_isoMu_Tau"].etaCut, trigParList["L1_isoMu_Tau"].minQual)){
			h_isoMu_Tau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	
	//------ isoMu_Tau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoMu_Tau_X = h2_isoMu_Tau_byThreshold->GetNbinsX();
		const unsigned n_bins_isoMu_Tau_Y = h2_isoMu_Tau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoMu_Tau_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_isoMu_Tau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoMu_Tau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoMu_Tau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoMu_Tau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoMu_Tau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_Tau(bin_center_X, bin_low_edge_Y, trigParList["L1_isoMu_Tau"].etaCut, trigParList["L1_isoMu_Tau"].minQual)){
			   	h2_isoMu_Tau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------	
	
	//-------- isoMu_isoTau -------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoMu_isoTau = h_isoMu_isoTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoMu_isoTau+1; bin++){
//		const float bin_low_edge = h_isoMu_isoTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoMu_isoTau_byThreshold->GetBinCenter(bin);
		if(Mu_isoTau(bin_center, bin_center*(trigParList["L1_isoMu_isoTau"].secTh/trigParList["L1_isoMu_isoTau"].primTh), trigParList["L1_isoMu_isoTau"].etaCut, trigParList["L1_isoMu_isoTau"].minQual)){
			h_isoMu_isoTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	
	//------ isoMu_isoTau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoMu_isoTau_X = h2_isoMu_isoTau_byThreshold->GetNbinsX();
		const unsigned n_bins_isoMu_isoTau_Y = h2_isoMu_isoTau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoMu_isoTau_X+1; bin++){	
//	   	const float bin_low_edge_X = h2_isoMu_isoTau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoMu_isoTau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoMu_isoTau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoMu_isoTau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoMu_isoTau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Mu_isoTau(bin_center_X, bin_low_edge_Y, trigParList["L1_isoMu_isoTau"].etaCut, trigParList["L1_isoMu_isoTau"].minQual)){
			   	h2_isoMu_isoTau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------	
	

	
	
	//------- EG_Tau --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_EG_Tau = h_EG_Tau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_EG_Tau+1; bin++){
		const float bin_low_edge = h_EG_Tau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_EG_Tau_byThreshold->GetBinCenter(bin);
		if(EG_Tau(bin_low_edge, bin_center*(trigParList["L1_EG_Tau"].secTh/trigParList["L1_EG_Tau"].primTh),  trigParList["L1_EG_Tau"].etaCut) ){
			h_EG_Tau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ EG_Tau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_EG_Tau_X = h2_EG_Tau_byThreshold->GetNbinsX();
		const unsigned n_bins_EG_Tau_Y = h2_EG_Tau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_EG_Tau_X+1; bin++){	
	   	const float bin_low_edge_X = h2_EG_Tau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_EG_Tau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_EG_Tau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_EG_Tau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_EG_Tau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(EG_Tau(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_EG_Tau"].etaCut)){
			   	h2_EG_Tau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------- isoEG_Tau --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoEG_Tau = h_isoEG_Tau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoEG_Tau+1; bin++){
		const float bin_low_edge = h_isoEG_Tau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoEG_Tau_byThreshold->GetBinCenter(bin);
		if(isoEG_Tau(bin_low_edge, bin_center*(trigParList["L1_isoEG_Tau"].secTh/trigParList["L1_isoEG_Tau"].primTh),  trigParList["L1_isoEG_Tau"].etaCut) ){
			h_isoEG_Tau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ isoEG_Tau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoEG_Tau_X = h2_isoEG_Tau_byThreshold->GetNbinsX();
		const unsigned n_bins_isoEG_Tau_Y = h2_isoEG_Tau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoEG_Tau_X+1; bin++){	
	   	const float bin_low_edge_X = h2_isoEG_Tau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoEG_Tau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoEG_Tau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoEG_Tau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoEG_Tau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(isoEG_Tau(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_isoEG_Tau"].etaCut)){
			   	h2_isoEG_Tau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------- isoEG_isoTau --------------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoEG_isoTau = h_isoEG_isoTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoEG_isoTau+1; bin++){
		const float bin_low_edge = h_isoEG_isoTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoEG_isoTau_byThreshold->GetBinCenter(bin);
		if(isoEG_isoTau(bin_low_edge, bin_center*(trigParList["L1_isoEG_isoTau"].secTh/trigParList["L1_isoEG_isoTau"].primTh),  trigParList["L1_isoEG_isoTau"].etaCut) ){
			h_isoEG_isoTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	//----------------------------------------------------------------------------------------------------------------------

	//------ isoEG_isoTau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_isoEG_isoTau_X = h2_isoEG_isoTau_byThreshold->GetNbinsX();
		const unsigned n_bins_isoEG_isoTau_Y = h2_isoEG_isoTau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_isoEG_isoTau_X+1; bin++){	
	   	const float bin_low_edge_X = h2_isoEG_isoTau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_isoEG_isoTau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_isoEG_isoTau_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_isoEG_isoTau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_isoEG_isoTau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(isoEG_isoTau(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_isoEG_isoTau"].etaCut)){
			   	h2_isoEG_isoTau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//------ SingleJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleJet = h_SingleJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleJet+1; bin++){
		const float bin_low_edge = h_SingleJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleJet_byThreshold->GetBinCenter(bin);
		if(SingleJet(bin_low_edge)){
			h_SingleJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SingleJetC -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleJetC = h_SingleJetC_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleJetC+1; bin++){
		const float bin_low_edge = h_SingleJetC_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleJetC_byThreshold->GetBinCenter(bin);
		if(SingleJetCentral(bin_low_edge)){
			h_SingleJetC_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleJet = h_DoubleJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleJet+1; bin++){
		const float bin_low_edge = h_DoubleJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleJet_byThreshold->GetBinCenter(bin);
		if(DoubleJetCentral(bin_low_edge,bin_center*(trigParList["L1_DoubleJet"].secTh/trigParList["L1_DoubleJet"].primTh))){
			h_DoubleJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	

	//------ DoubleJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_DoubleJet_X = h2_DoubleJet_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleJet_Y = h2_DoubleJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin ; ybin++){  //only compute the rate below the diagonal because symmetric	
		   	const float bin_low_edge_Y = h2_DoubleJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleJetCentral(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_DoubleJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------ DoubleFwdJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleFwdJet = h_DoubleFwdJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleFwdJet+1; bin++){
		const float bin_low_edge = h_DoubleFwdJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleFwdJet_byThreshold->GetBinCenter(bin);
		if(DoubleJetForward(bin_low_edge,bin_center*(trigParList["L1_DoubleFwdJet"].secTh/trigParList["L1_DoubleFwdJet"].primTh))){
			h_DoubleFwdJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------
	

	//------ DoubleFwdJet ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_DoubleFwdJet_X = h2_DoubleFwdJet_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleFwdJet_Y = h2_DoubleFwdJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleFwdJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleFwdJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleFwdJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin ; ybin++){  //only compute the rate below the diagonal because symmetric	
		   	const float bin_low_edge_Y = h2_DoubleFwdJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleFwdJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleJetForward(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_DoubleFwdJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	
	//------- QuadJet ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_QuadJetCentral = h_QuadJetC_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_QuadJetCentral+1; bin++){
		const float bin_low_edge = h_QuadJetC_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_QuadJetC_byThreshold->GetBinCenter(bin);
		if(QuadJetCentral(bin_low_edge, bin_center*(trigParList["L1_QuadJetC"].secTh/trigParList["L1_QuadJetC"].primTh) , bin_center*(trigParList["L1_QuadJetC"].triTh/trigParList["L1_QuadJetC"].primTh), bin_center*(trigParList["L1_QuadJetC"].quadTh/trigParList["L1_QuadJetC"].primTh))){
			h_QuadJetC_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ QuadJet ---- 2-D Evaluation (A)------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_QuadJetCentral_X = h2A_QuadJetCentral_byThreshold->GetNbinsX();
//		const unsigned n_bins_QuadJetCentral_Y = h2A_QuadJetCentral_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_QuadJetCentral_X+1; bin++){	
	   	const float bin_low_edge_X = h2A_QuadJetCentral_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2A_QuadJetCentral_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){ //calculate only below the diagonal since expect Et ordering	
		   	const float bin_low_edge_Y = h2A_QuadJetCentral_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2A_QuadJetCentral_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(QuadJetCentral(bin_low_edge_X, bin_low_edge_Y,  bin_low_edge_Y,  bin_low_edge_Y)){
			   	h2A_QuadJetCentral_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------ QuadJet ---- 2-D Evaluation (B)------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_QuadJetCentral_X = h2B_QuadJetCentral_byThreshold->GetNbinsX();
//		const unsigned n_bins_QuadJetCentral_Y = h2B_QuadJetCentral_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_QuadJetCentral_X+1; bin++){	
	   	const float bin_low_edge_X = h2B_QuadJetCentral_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2B_QuadJetCentral_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	//calculate only below the diagonal since expect Et ordering
		   	const float bin_low_edge_Y = h2B_QuadJetCentral_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2B_QuadJetCentral_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(QuadJetCentral(bin_low_edge_X, bin_low_edge_X,  bin_low_edge_Y,  bin_low_edge_Y)){
			   	h2B_QuadJetCentral_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------- SixJet ------------------------------------------------------------------------------------------------------
	const unsigned n_bins_SixJet = h_SixJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SixJet+1; bin++){
		const float bin_low_edge = h_SixJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SixJet_byThreshold->GetBinCenter(bin);
		if(MultiJet(6,bin_low_edge, bin_center*(trigParList["L1_SixJet"].secTh/trigParList["L1_SixJet"].primTh) , bin_center*(trigParList["L1_SixJet"].triTh/trigParList["L1_SixJet"].primTh), bin_center*(trigParList["L1_SixJet"].quadTh/trigParList["L1_SixJet"].primTh))){
			h_SixJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------


	//------ SixJet ---- 2-D Evaluation (A: first two jets have highest threshold (x-axis) and all others are lower (y-axis) )------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SixJet_X = h2A_SixJet_byThreshold->GetNbinsX();
//		const unsigned n_bins_SixJet_Y = h2A_SixJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SixJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2A_SixJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2A_SixJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){ //calculate only below the diagonal since expect Et ordering	
		   	const float bin_low_edge_Y = h2A_SixJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2A_SixJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(MultiJet(6,bin_low_edge_X, bin_low_edge_X,  bin_low_edge_Y,  bin_low_edge_Y)){
			   	h2A_SixJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	
	//------ SingleTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau = h_SingleTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau+1; bin++){
		const float bin_low_edge = h_SingleTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_byThreshold->GetBinCenter(bin);
		if(SingleTauJet(bin_low_edge,trigParList["L1_SingleTau"].etaCut )){
			h_SingleTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	
	//------ SingleIsoTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoTau = h_SingleIsoTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoTau+1; bin++){
		const float bin_low_edge = h_SingleIsoTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoTau_byThreshold->GetBinCenter(bin);
		if(SingleIsoTauJet(bin_low_edge,trigParList["L1_SingleIsoTau"].etaCut)){
			h_SingleIsoTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}	
	
	//----------------------------------------------------------------------------------------------------------------------			
	
	//------ DoubleTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleTau = h_DoubleTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleTau+1; bin++){
		const float bin_low_edge = h_DoubleTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleTau_byThreshold->GetBinCenter(bin);
		if(DoubleTauJetEta(bin_low_edge,bin_center*(trigParList["L1_DoubleTau"].secTh/trigParList["L1_DoubleTau"].primTh),trigParList["L1_DoubleTau"].etaCut)){
			h_DoubleTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	


	//------ DoubleTau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_DoubleTau_X = h2_DoubleTau_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleTau_Y = h2_DoubleTau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleTau_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleTau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleTau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	//calculate only below the diagonal because it should be symmetric.
		   	const float bin_low_edge_Y = h2_DoubleTau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleTau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleTauJetEta(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_DoubleTau"].etaCut)){
			   	h2_DoubleTau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------ DoubleIsoTau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleIsoTau = h_DoubleIsoTau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleIsoTau+1; bin++){
		const float bin_low_edge = h_DoubleIsoTau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleIsoTau_byThreshold->GetBinCenter(bin);
		if(DoubleIsoTau(bin_low_edge,bin_center*(trigParList["L1_DoubleIsoTau"].secTh/trigParList["L1_DoubleIsoTau"].primTh),trigParList["L1_DoubleIsoTau"].etaCut)){
			h_DoubleIsoTau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	


	//------ DoubleIsoTau ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_DoubleIsoTau_X = h2_DoubleIsoTau_byThreshold->GetNbinsX();
//		const unsigned n_bins_DoubleIsoTau_Y = h2_DoubleIsoTau_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleIsoTau_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleIsoTau_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleIsoTau_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= bin; ybin++){	//calculate only below the diagonal because it should be symmetric.
		   	const float bin_low_edge_Y = h2_DoubleIsoTau_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleIsoTau_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleIsoTau(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_DoubleIsoTau"].etaCut)){
			   	h2_DoubleIsoTau_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



	//------ isoTau_Tau -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_isoTau_Tau = h_isoTau_Tau_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_isoTau_Tau+1; bin++){
		const float bin_low_edge = h_isoTau_Tau_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_isoTau_Tau_byThreshold->GetBinCenter(bin);
		if(isoTau_Tau(bin_low_edge,bin_center*(trigParList["L1_isoTau_Tau"].secTh/trigParList["L1_isoTau_Tau"].primTh),trigParList["L1_isoTau_Tau"].etaCut)){
			h_isoTau_Tau_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------	


	//------ SingleTau +ETM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau_ETM = h_SingleTau_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau_ETM+1; bin++){
		const float bin_low_edge = h_SingleTau_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_ETM_byThreshold->GetBinCenter(bin);
		if(Tau_ETM(bin_low_edge,bin_center*(trigParList["L1_SingleTau_ETM"].secTh/trigParList["L1_SingleTau_ETM"].primTh), trigParList["L1_SingleTau_ETM"].etaCut)){
			h_SingleTau_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ SingleTau + ETM---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleTau_ETM_X = h2_SingleTau_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleTau_ETM_Y = h2_SingleTau_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleTau_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleTau_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleTau_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleTau_ETM_Y+1; ybin++){	//calculate only below the diagonal because it should be symmetric.
		   	const float bin_low_edge_Y = h2_SingleTau_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleTau_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Tau_ETM(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_SingleTau_ETM"].etaCut)){
			   	h2_SingleTau_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoTau +ETM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoTau_ETM = h_SingleIsoTau_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoTau_ETM+1; bin++){
		const float bin_low_edge = h_SingleIsoTau_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoTau_ETM_byThreshold->GetBinCenter(bin);
		if(IsoTau_ETM(bin_low_edge,bin_center*(trigParList["L1_SingleIsoTau_ETM"].secTh/trigParList["L1_SingleIsoTau_ETM"].primTh), trigParList["L1_SingleIsoTau_ETM"].etaCut)){
			h_SingleIsoTau_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			


	//------ SingleTau +HTM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau_HTM = h_SingleTau_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau_HTM+1; bin++){
		const float bin_low_edge = h_SingleTau_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_HTM_byThreshold->GetBinCenter(bin);
		if(Tau_HTM(bin_low_edge,bin_center*(trigParList["L1_SingleTau_HTM"].secTh/trigParList["L1_SingleTau_HTM"].primTh), trigParList["L1_SingleTau_HTM"].etaCut)){
			h_SingleTau_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ SingleIsoTau +HTM -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoTau_HTM = h_SingleIsoTau_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoTau_HTM+1; bin++){
		const float bin_low_edge = h_SingleIsoTau_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoTau_HTM_byThreshold->GetBinCenter(bin);
		if(IsoTau_HTM(bin_low_edge,bin_center*(trigParList["L1_SingleIsoTau_HTM"].secTh/trigParList["L1_SingleIsoTau_HTM"].primTh), trigParList["L1_SingleIsoTau_HTM"].etaCut)){
			h_SingleIsoTau_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			


/* boolean for these need to be fixed
	//------ SingleTau + CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau_CJet = h_SingleTau_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau_CJet+1; bin++){
		const float bin_low_edge = h_SingleTau_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_CJet_byThreshold->GetBinCenter(bin);
		if(Tau_JetCentral(bin_low_edge,bin_center*(trigParList["L1_SingleTau_CJet"].secTh/trigParList["L1_SingleTau_CJet"].primTh), trigParList["L1_SingleTau_CJet"].etaCut)){
			h_SingleTau_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ SingleTau + CJet---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleTau_CJet_X = h2_SingleTau_CJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleTau_CJet_Y = h2_SingleTau_CJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleTau_CJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleTau_CJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleTau_CJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleTau_CJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleTau_CJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleTau_CJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Tau_JetCentral(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_SingleTau_CJet"].etaCut)){
			   	h2_SingleTau_CJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------ SingleIsoTau + CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleIsoTau_CJet = h_SingleIsoTau_CJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleIsoTau_CJet+1; bin++){
		const float bin_low_edge = h_SingleIsoTau_CJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleIsoTau_CJet_byThreshold->GetBinCenter(bin);
		if(IsoTau_JetCentral(bin_low_edge,bin_center*(trigParList["L1_SingleIsoTau_CJet"].secTh/trigParList["L1_SingleIsoTau_CJet"].primTh), trigParList["L1_SingleIsoTau_CJet"].etaCut)){
			h_SingleIsoTau_CJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			
*/


	//------ SingleTau + 2FJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleTau_TwoFJet = h_SingleTau_TwoFJet_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleTau_TwoFJet+1; bin++){
		const float bin_low_edge = h_SingleTau_TwoFJet_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleTau_TwoFJet_byThreshold->GetBinCenter(bin);
		if(Tau_TwoJetForward(bin_low_edge,bin_center*(trigParList["L1_SingleTau_TwoFJet"].secTh/trigParList["L1_SingleTau_TwoFJet"].primTh),bin_center*(trigParList["L1_SingleTau_TwoFJet"].triTh/trigParList["L1_SingleTau_TwoFJet"].primTh))){
			h_SingleTau_TwoFJet_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ SingleTau + 2FJet---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleTau_TwoFJet_X = h2_SingleTau_TwoFJet_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleTau_TwoFJet_Y = h2_SingleTau_TwoFJet_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleTau_TwoFJet_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleTau_TwoFJet_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleTau_TwoFJet_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleTau_TwoFJet_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleTau_TwoFJet_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleTau_TwoFJet_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(Tau_TwoJetForward(bin_low_edge_X, bin_low_edge_Y, bin_low_edge_Y)){
			   	h2_SingleTau_TwoFJet_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//-------- HTT ---------------------------------------------------------------------------------------------------------
	const unsigned n_bins_HTT = h_HTT_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_HTT+1; bin++){
		const float bin_low_edge = h_HTT_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_HTT_byThreshold->GetBinCenter(bin);
		if(HTT(bin_low_edge)){
			h_HTT_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//------- ETM ----------------------------------------------------------------------------------------------------------
	const unsigned n_bins_ETM = h_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_ETM+1; bin++){
		const float bin_low_edge = h_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_ETM_byThreshold->GetBinCenter(bin);
		if(ETM(bin_low_edge)){
			h_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}

	//------- HTM ----------------------------------------------------------------------------------------------------------
	const unsigned n_bins_HTM = h_HTM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_HTM+1; bin++){
		const float bin_low_edge = h_HTM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_HTM_byThreshold->GetBinCenter(bin);
		if(HTM(bin_low_edge)){
			h_HTM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}

	//-------- HTT + ETM ---------------------------------------------------------------------------------------------------------
	const unsigned n_bins_HTT_ETM = h_HTT_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_HTT_ETM+1; bin++){
		const float bin_low_edge = h_HTT_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_HTT_ETM_byThreshold->GetBinCenter(bin);
		if(HTT_ETM(bin_low_edge,bin_center*(trigParList["L1_HTT_ETM"].secTh/trigParList["L1_HTT_ETM"].primTh))){
			h_HTT_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}

	//------ HTT + ETM ---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_HTT_ETM_X = h2_HTT_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_HTT_ETM_Y = h2_HTT_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_HTT_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_HTT_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_HTT_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_HTT_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_HTT_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_HTT_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(HTT_ETM(bin_low_edge_X, bin_low_edge_Y)){
			   	h2_HTT_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------


	//------ ETM + CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_SingleCJet_ETM = h_SingleCJet_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_SingleCJet_ETM+1; bin++){
		const float bin_low_edge = h_SingleCJet_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_SingleCJet_ETM_byThreshold->GetBinCenter(bin);
		if(JetCentral_ETM(bin_low_edge, bin_center*(trigParList["L1_SingleCJet_ETM"].secTh/trigParList["L1_SingleCJet_ETM"].primTh), trigParList["L1_SingleCJet_ETM"].etaCut)){
			h_SingleCJet_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ ETM + CJet---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_SingleCJet_ETM_X = h2_SingleCJet_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_SingleCJet_ETM_Y = h2_SingleCJet_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_SingleCJet_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_SingleCJet_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_SingleCJet_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_SingleCJet_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_SingleCJet_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_SingleCJet_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(JetCentral_ETM(bin_low_edge_X, bin_low_edge_Y, trigParList["L1_SingleCJet_ETM"].etaCut)){
			   	h2_SingleCJet_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------

	//------ ETM + CJet -----------------------------------------------------------------------------------------------------
	const unsigned n_bins_DoubleCJet_ETM = h_DoubleCJet_ETM_byThreshold->GetNbinsX();
	for(unsigned bin=1; bin <= n_bins_DoubleCJet_ETM+1; bin++){
		const float bin_low_edge = h_DoubleCJet_ETM_byThreshold->GetBinLowEdge(bin);
		const float bin_center = h_DoubleCJet_ETM_byThreshold->GetBinCenter(bin);
		if(DoubleJetCentral_ETM(bin_low_edge,bin_center*(trigParList["L1_DoubleCJet_ETM"].secTh/trigParList["L1_DoubleCJet_ETM"].primTh), bin_center*(trigParList["L1_DoubleCJet_ETM"].triTh/trigParList["L1_DoubleCJet_ETM"].primTh), trigParList["L1_DoubleCJet_ETM"].etaCut)){
			h_DoubleCJet_ETM_byThreshold->Fill(bin_center,lumiWeight);
		}
	}
	//----------------------------------------------------------------------------------------------------------------------			

	//------ ETM + CJet---- 2-D Evaluation ------------------------------------------------------------------------------
	if(calcThreshold>1) {	
		const unsigned n_bins_DoubleCJet_ETM_X = h2_DoubleCJet_ETM_byThreshold->GetNbinsX();
		const unsigned n_bins_DoubleCJet_ETM_Y = h2_DoubleCJet_ETM_byThreshold->GetNbinsY();
		for(unsigned bin=1; bin <= n_bins_DoubleCJet_ETM_X+1; bin++){	
	   	const float bin_low_edge_X = h2_DoubleCJet_ETM_byThreshold->GetXaxis()->GetBinLowEdge(bin);	
			const float bin_center_X = h2_DoubleCJet_ETM_byThreshold->GetXaxis()->GetBinCenter(bin);
			for(unsigned ybin=1; ybin <= n_bins_DoubleCJet_ETM_Y+1; ybin++){	
		   	const float bin_low_edge_Y = h2_DoubleCJet_ETM_byThreshold->GetYaxis()->GetBinLowEdge(ybin);	
		   	const float bin_center_Y = h2_DoubleCJet_ETM_byThreshold->GetYaxis()->GetBinCenter(ybin);
		   	if(DoubleJetCentral_ETM(bin_low_edge_X, bin_low_edge_X, bin_low_edge_Y, trigParList["L1_DoubleCJet_ETM"].etaCut)){
			   	h2_DoubleCJet_ETM_byThreshold->Fill(bin_center_X,bin_center_Y,lumiWeight);
		   	}
			}	
		}
   }
	//----------------------------------------------------------------------------------------------------------------------



}


Bool_t L1Menu2015::Mu_EG(Float_t mucut, Float_t EGcut , Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }		
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}


Bool_t L1Menu2015::Mu_isoEG(Float_t mucut, Float_t EGcut , Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }		
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}


Bool_t L1Menu2015::EG_Mu(Float_t EGcut, Float_t mucut , Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }		
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::isoEG_Mu(Float_t EGcut, Float_t mucut , Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }		
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}


Bool_t L1Menu2015::Mu_Tau(Float_t Mucut, Float_t taucut , Float_t etaCut, Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t tau =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);
      Float_t eta = myEvt_.Etamu.at(imu);					
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,Mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= Mucut) muon = true;
	   }		
	}

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (!isTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) tau = true;
	}

	Bool_t ok = muon && tau;
	return ok;

}


Bool_t L1Menu2015::Mu_isoTau(Float_t Mucut, Float_t taucut , Float_t etaCut, Int_t minMuQual) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t tau =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {   
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);
      Float_t eta = myEvt_.Etamu.at(imu);					
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,Mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= Mucut) muon = true;
	   }		
	}

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (!isTauJet || !myEvt_.isoTaujet[ue]) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) tau = true;
	}

	Bool_t ok = muon && tau;
	return ok;

}


Bool_t L1Menu2015::EG_Tau(Float_t EGcut, Float_t taucut , Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t tau =false;
	Bool_t eg = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) {
		   eg = true;
		
// Now look for a tau that is not the same as this eg
		  Int_t Nj = myEvt_.Njet ;
		  for (Int_t uj=0; uj < Nj; uj++) {
			  bx = myEvt_.Bxjet[uj];        		
			  if (bx != 0) continue;
			  Bool_t isTauJet = myEvt_.Taujet[uj];
			  if (!isTauJet) continue;
			  Float_t eta = myEvt_.Etajet[uj];
			  if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
			  if (myEvt_.Etajet[uj] == myEvt_.Etael[ue] && myEvt_.Phijet[uj] == myEvt_.Phiel[ue] ) continue;
			  Float_t rankt = myEvt_.Etjet[uj];
			  Float_t ptt = rankt; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rank*4.,theL1JetCorrection);
			  if (ptt >= taucut) tau = true;
		  }
		
		
		}
	}  // end loop over EM objects



	Bool_t ok = eg && tau;
	return ok;

}


Bool_t L1Menu2015::isoEG_Tau(Float_t EGcut, Float_t taucut , Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t tau =false;
	Bool_t eg = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) {
		   eg = true;
		
// Now look for a tau that is not the same as this eg
		  Int_t Nj = myEvt_.Njet ;
		  for (Int_t uj=0; uj < Nj; uj++) {
			  bx = myEvt_.Bxjet[uj];        		
			  if (bx != 0) continue;
			  Bool_t isTauJet = myEvt_.Taujet[uj];
			  if (!isTauJet) continue;
			  Float_t eta = myEvt_.Etajet[uj];
			  if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
			  if (myEvt_.Etajet[uj] == myEvt_.Etael[ue] && myEvt_.Phijet[uj] == myEvt_.Phiel[ue] ) continue;
			  Float_t rankt = myEvt_.Etjet[uj];
			  Float_t ptt = rankt; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rank*4.,theL1JetCorrection);
			  if (ptt >= taucut) tau = true;
		  }
		
		
		}
	}  // end loop over EM objects



	Bool_t ok = eg && tau;
	return ok;

}


Bool_t L1Menu2015::isoEG_isoTau(Float_t EGcut, Float_t taucut , Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];    // ZeroBias
	if (! raw) return false;


	Bool_t tau =false;
	Bool_t eg = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) {
		   eg = true;
		
// Now look for a tau that is not the same as this eg
		  Int_t Nj = myEvt_.Njet ;
		  for (Int_t uj=0; uj < Nj; uj++) {
			  bx = myEvt_.Bxjet[uj];        		
			  if (bx != 0) continue;
			  Bool_t isTauJet = myEvt_.Taujet[uj];
			  if (!isTauJet || !myEvt_.isoTaujet[uj]) continue;
			  Float_t eta = myEvt_.Etajet[uj];
			  if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
			  if (myEvt_.Etajet[uj] == myEvt_.Etael[ue] && myEvt_.Phijet[uj] == myEvt_.Phiel[ue] ) continue;
			  Float_t rankt = myEvt_.Etjet[uj];
			  Float_t ptt = rankt; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rank*4.,theL1JetCorrection);
			  if (ptt >= taucut) tau = true;
		  }
		
		
		}
	}  // end loop over EM objects



	Bool_t ok = eg && tau;
	return ok;

}




Bool_t L1Menu2015::DoubleMu_EG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0]; 	// ZeroBias
	if (! raw) return false;

	Bool_t eg =false;
	Bool_t muon = false;
	Int_t  Nmuons = 0;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		// if ( qual < 4) continue;
		if (qual < 4 && qual !=3 ) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) Nmuons ++;
 		} else {
		   if (pt >= mucut) Nmuons ++;
	   }		

	}
	if (Nmuons >= 2) muon = true;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::Mu_DoubleEG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias..
	if (! raw) return false;

	Bool_t eg =false;
	Bool_t muon = false;
	Int_t  Nmuons = 0;
	Int_t Nelectrons = 0;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) Nmuons ++;
 		} else {
		   if (pt >= mucut) Nmuons ++;
	   }		
	}
	if (Nmuons >= 1) muon = true;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) Nelectrons ++;
	}  // end loop over EM objects
	if (Nelectrons >= 2) eg = true;

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::MuOpen_EG(Float_t mucut, Float_t EGcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;


	Bool_t eg =false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = muon && eg;
	return ok;

}

Bool_t L1Menu2015::Mu_JetCentral(Float_t mucut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_DoubleJetCentral(Float_t mucut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t n1 = 0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) n1 ++;
	}
	jet = ( n1 >= 2 );

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut ) {

///  CORRUPT NEEDS FIXING FOR NEW TAU DEFINITIONS.

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t central = false;
	Bool_t tau = false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;	
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) {  	// look at CentralJet
			if (pt >= jetcut) central = true;
		}
		else   {		// look at TauJets
			if (pt >= taucut) tau = true;
		}
	}
	jet = central || tau  ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Muer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

   
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::IsoMuer_JetCentral(Float_t mucut, Float_t jetcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

   
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;				
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
//		if (bx != 0 || !myEvt_.Isomu.at(imu)) continue;
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeIsolatedMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}



Bool_t L1Menu2015::Muer_JetCentral_LowerTauTh(Float_t mucut, Float_t jetcut, Float_t taucut) {

//  	CORRUPT NEEDS UPDATING FOR NEW TAU APPROACH

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t central = false;
	Bool_t tau = false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;

		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) {       // look at CentralJet
			if (pt >= jetcut) central = true;
		}
		else   {                // look at TauJets
			if (pt >= taucut) tau = true;
		}
	}
	jet = central || tau  ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,2.1)      ) muon = true;
 		} else {
		   if (fabs(eta) > 2.1) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mia(Float_t mucut, Float_t jet1, Float_t jet2) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;
	Int_t n1 = 0;
	Int_t n2 = 0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jet1) n1 ++;
		if (pt >= jet2) n2 ++;
	}       
	jet = (n1 >= 1 && n2 >= 2 );

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,2.1)      ) muon = true;
 		} else {
		   if (fabs(eta) > 2.1) continue;
		   if (pt >= mucut) muon = true;
	   }
	} 

	Bool_t ok = muon && jet;
	if (! ok) return false;

	// now the CORREL condition


	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        
		if (fabs(eta) > 2.1) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
		Float_t etamu = myEvt_.Etamu.at(imu);
		Int_t ieta_mu = etaINjetCoord(etamu);

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
		   Bool_t isTauJet = myEvt_.Taujet[ue];
		   if (isTauJet) continue;		
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jet2) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
			Float_t etajet = myEvt_.Etajet[ue];
			Int_t ieta_jet = (int)etajet;

			Bool_t corr_phi = correlateInPhi(iphi_jet, iphi_mu);
			Bool_t corr_eta = correlateInEta(ieta_jet, ieta_mu);
			Bool_t corr = corr_phi && corr_eta;
			if (corr) CORREL = true ;
		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Mu_JetCentral_delta(Float_t mucut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
		Float_t etamu = myEvt_.Etamu.at(imu);
		Int_t ieta_mu = etaINjetCoord(etamu);

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
		   Bool_t isTauJet = myEvt_.Taujet[ue];
		   if (isTauJet) continue;		
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
			Float_t etajet = myEvt_.Etajet[ue];
			Int_t ieta_jet = (int)etajet;

			Bool_t corr = correlateInPhi(iphi_jet, iphi_mu, 2) && correlateInEta(ieta_jet, ieta_mu, 2);
			if (corr) CORREL = true ;
		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Mu_JetCentral_deltaOut(Float_t mucut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;

	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		if (pt < mucut) continue;

		Float_t phimu = myEvt_.Phimu.at(imu);
		Int_t iphi_mu = phiINjetCoord(phimu);
//          Float_t etamu = myEvt_.Etamu.at(imu);
//          Int_t ieta_mu = etaINjetCoord(etamu);

// 		Int_t PhiOut[3];
// 		PhiOut[0] = iphi_mu;
// 		if (iphi_mu< 17) PhiOut[1] = iphi_mu+1;
// 		if (iphi_mu == 17) PhiOut[1] = 0;
// 		if (iphi_mu > 0) PhiOut[2] = iphi_mu - 1;
// 		if (iphi_mu == 0) PhiOut[2] = 17;

		for (Int_t ue=0; ue < Nj; ue++) {
			Int_t bxj = myEvt_.Bxjet[ue];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[ue];
			if (isFwdJet) continue;
  		   Bool_t isTauJet = myEvt_.Taujet[ue];
		   if (isTauJet) continue;		
			Float_t rank = myEvt_.Etjet[ue];
			Float_t ptj = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[ue];
			Int_t iphi_jet = (int)phijet;
//                  Float_t etajet = myEvt_.Etajet[ue];
//                  Int_t ieta_jet = (int)etajet;

			if (! correlateInPhi(iphi_jet, iphi_mu, 8)) CORREL = true;


		}
	}

	return CORREL;

}

Bool_t L1Menu2015::Muer_TripleJetCentral(Float_t mucut, Float_t jet1, Float_t jet2, Float_t jet3)  {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t muon = false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jet1) n1 ++;
		if (pt >= jet2) n2 ++;
		if (pt >= jet3) n3 ++;
	}

	jet = ( n1 >= 1 && n2 >= 2 && n3 >= 3 ) ;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;
		Float_t eta = myEvt_.Etamu.at(imu) ;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,2.1)      ) muon = true;
 		} else {
		   if (fabs(eta) > 2.1) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Bool_t ok = muon && jet;
	return ok;

}

Bool_t L1Menu2015::Mu_HTT(Float_t mucut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ht=false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
      Float_t eta = myEvt_.Etamu.at(imu);
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < 4) continue;

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,5.0)      ) muon = true;
 		} else {
		   if (pt >= mucut) muon = true;
	   }
	}

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = muon && ht;
	return ok;

}

Bool_t L1Menu2015::Muer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Bool_t ok = muon && etm;
	return ok;

}

Bool_t L1Menu2015::IsoMuer_ETM(Float_t mucut, Float_t ETMcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
//		if (bx != 0 || !myEvt_.Isomu.at(imu)) continue;
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeIsolatedMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Bool_t ok = muon && etm;
	return ok;

}


Bool_t L1Menu2015::Muer_HTM(Float_t mucut, Float_t HTMcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htm = false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Float_t adc = myEvt_.HTM ;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut);

	Bool_t ok = muon && htm;
	return ok;

}

Bool_t L1Menu2015::IsoMuer_HTM(Float_t mucut, Float_t HTMcut, Float_t etacut, Int_t minMuQual ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htm = false;
	Bool_t muon = false;

	Int_t Nmu = myEvt_.Nmu;
	for (Int_t imu=0; imu < Nmu; imu++) {
		Int_t bx = myEvt_.Bxmu.at(imu);		
//		if (bx != 0 || !myEvt_.Isomu.at(imu)) continue;
		if (bx != 0) continue;
		Float_t pt = myEvt_.Ptmu.at(imu);			
		Int_t qual = myEvt_.Qualmu.at(imu);        
		if ( qual < minMuQual) continue;
		Float_t eta = myEvt_.Etamu.at(imu);        

		if(useUpgradeMuons) {
   		if( UpgradeIsolatedMuon(pt,eta,mucut,etacut)      ) muon = true;
 		} else {
		   if (fabs(eta) > etacut) continue;
		   if (pt >= mucut) muon = true;
	   }
	}

	Float_t adc = myEvt_.HTM ;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut);

	Bool_t ok = muon && htm;
	return ok;

}


Bool_t L1Menu2015::EG_FwdJet(Float_t EGcut, Float_t FWcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {        
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (!isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= FWcut) jet = true;
	}

	Bool_t ok = ( eg && jet);
	return ok;

}

/*
Bool_t L1Menu2015::EG_JetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t eg = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = eg && jet;
	return ok;

}

Bool_t L1Menu2015::IsoEG_JetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t eg = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Bool_t ok = eg && jet;
	return ok;

}
*/

Bool_t L1Menu2015::EG_JetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t eg = false;
   Bool_t ok = false;
	
	Int_t Nj = myEvt_.Njet ;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut){
		
		    eg = true;
		
	   	 for (Int_t uj=0; uj < Nj; uj++) {
				  Int_t bxj = myEvt_.Bxjet[uj];        		
				  if (bxj != 0) continue;
				  Bool_t isFwdJet = myEvt_.Fwdjet[uj];
				  if (isFwdJet) continue;
				  Bool_t isTauJet = myEvt_.Taujet[uj];
				  if (isTauJet) continue;		
				  Float_t rankj = myEvt_.Etjet[uj];
				  Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rank*4.,theL1JetCorrection);


				  if (ptj >= jetcut &&
				     (myEvt_.Etajet[uj]!=myEvt_.Etael[ue]) &&
					  (myEvt_.Phijet[uj]!=myEvt_.Phiel[ue]) ) jet = true;		
	   	 }

		    ok = eg && jet;
		} // if good EG
	}  // end loop over EM objects


	return ok;

}

Bool_t L1Menu2015::IsoEG_JetCentral(Float_t EGcut, Float_t jetcut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;
	Bool_t eg = false;
   Bool_t ok = false;
	
	Int_t Nj = myEvt_.Njet ;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16		
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut){
		
		    eg = true;
		
	   	 for (Int_t uj=0; uj < Nj; uj++) {
				  Int_t bxj = myEvt_.Bxjet[uj];        		
				  if (bxj != 0) continue;
				  Bool_t isFwdJet = myEvt_.Fwdjet[uj];
				  if (isFwdJet) continue;
				  Bool_t isTauJet = myEvt_.Taujet[uj];
				  if (isTauJet) continue;		
				  Float_t rankj = myEvt_.Etjet[uj];
				  Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rank*4.,theL1JetCorrection);


				  if (ptj >= jetcut &&
				     (myEvt_.Etajet[uj]!=myEvt_.Etael[ue]) &&
					  (myEvt_.Phijet[uj]!=myEvt_.Phiel[ue]) ) jet = true;		
	   	 }

		    ok = eg && jet;
		} // if good EG
	}  // end loop over EM objects

	return ok;

}


Bool_t L1Menu2015::EG_DoubleJetCentral(Float_t EGcut, Float_t jetcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	Int_t njets = 0;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) njets ++;
	}
	jet = ( njets >= 2 );

	Bool_t ok = ( eg && jet);
	return ok;

}

Bool_t L1Menu2015::EG_HT(Float_t EGcut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t ht = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = ( eg && ht);
	return ok;

}


Bool_t L1Menu2015::EG_ETM(Float_t EGcut, Float_t ETMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t etm = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.ETM ;;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut) ;

	Bool_t ok = ( eg && etm);
	
	//printf("EG_ETM:  eg %i  etm %i \n",eg,etm);
	
	return ok;

}


Bool_t L1Menu2015::IsoEG_ETM(Float_t EGcut, Float_t ETMcut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t etm = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;

		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16

		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.ETM ;;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut) ;

	Bool_t ok = ( eg && etm);
	
	//printf("EG_ETM:  eg %i  etm %i \n",eg,etm);
	
	return ok;

}


Bool_t L1Menu2015::EG_HTM(Float_t EGcut, Float_t HTMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t htm = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.HTM ;;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut) ;

	Bool_t ok = ( eg && htm);
	
	//printf("EG_ETM:  eg %i  etm %i \n",eg,etm);
	
	return ok;

}


Bool_t L1Menu2015::IsoEG_HTM(Float_t EGcut, Float_t HTMcut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t htm = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0 || !myEvt_.Isoel[ue]) continue;

		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16

		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Float_t adc = myEvt_.HTM ;;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut) ;

	Bool_t ok = ( eg && htm);
	
	//printf("EG_ETM:  eg %i  etm %i \n",eg,etm);
	
	return ok;

}



Bool_t L1Menu2015::DoubleEG_HT(Float_t EGcut, Float_t HTcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Int_t n1 = 0;
	Bool_t ht = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) n1 ++;
	}  // end loop over EM objects
	eg = ( n1 >= 2 );

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;
	ht = (TheHTT >= HTcut) ;

	Bool_t ok = ( eg && ht);
	return ok;

}

Bool_t L1Menu2015::EGEta2p1_JetCentral(Float_t EGcut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;
	
	Bool_t eg = false;
	Bool_t jet = false;

	Int_t Nele = myEvt_.Nele; 
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Bool_t ok = (eg && jet);
	if (! ok) return false;


	//  -- now evaluate the delta condition :

	Bool_t CORREL = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
		   Bool_t isTauJet = myEvt_.Taujet[ue];
		   if (isTauJet) continue;		
			Float_t rankj = myEvt_.Etjet[uj];
			// Float_t ptj = rankj * 4;
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet; 

			if ( iphijet != PhiOut[0] && 
				iphijet != PhiOut[1] &&
				iphijet != PhiOut[2] ) CORREL = true;
		}  // loop over jets

	}  // end loop over EM objects

	return CORREL;
	
}

Bool_t L1Menu2015::EGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut) {

//  CORRUPT.  NEEDS UPDATING FOR NEW TAU APPROACH

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Bool_t central = false;
	Bool_t tau = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (! isTauJet) {
			if (pt >= jetcut) central = true;
		}
		else {
			if (pt >= taucut) tau = true;
		}
	}
	jet = tau || central;

	Bool_t ok = (eg && jet);
	if (! ok) return false;

	//  -- now evaluate the delta condition :

	Bool_t CORREL_CENTRAL = false;
	Bool_t CORREL_TAU = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;	
			Bool_t isTauJet = myEvt_.Taujet[uj];
			Float_t rankj = myEvt_.Etjet[uj];
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if (! isTauJet) {

				if (ptj >= jetcut) { 
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_CENTRAL = true;
				}

			}
			else {
				if (ptj >= taucut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_TAU = true;
				}

			}


		}  // loop over jets

	}  // end loop over EM objects

	Bool_t CORREL = CORREL_CENTRAL || CORREL_TAU ;
	return CORREL;

}

Bool_t L1Menu2015::IsoEGEta2p1_JetCentral_LowTauTh(Float_t EGcut, Float_t jetcut, Float_t taucut) {


// CORRUPT.  NEEDS UPDATING FOR NEW TAU APPROACH

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Bool_t central = false;
	Bool_t tau = false;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if ( ! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (! isTauJet) {
			if (pt >= jetcut) central = true;
		}
		else {
			if (pt >= taucut) tau = true;
		}
	}
	jet = tau || central;

	Bool_t ok = (eg && jet);
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL_CENTRAL = false;
	Bool_t CORREL_TAU = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if ( ! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0;   

		PhiOut[0] = iphiel; 
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue;
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
			Bool_t isTauJet = myEvt_.Taujet[uj];
			Float_t rankj = myEvt_.Etjet[uj];
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if (! isTauJet) {

				if (ptj >= jetcut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_CENTRAL = true;
				}

			}
			else {
				if (ptj >= taucut) {
					if ( iphijet != PhiOut[0] &&
						iphijet != PhiOut[1] &&
						iphijet != PhiOut[2] ) CORREL_TAU = true;
				}

			}


		}  // loop over jets

	}  // end loop over EM objects

	Bool_t CORREL = CORREL_CENTRAL || CORREL_TAU ;

	return CORREL;

}

Bool_t L1Menu2015::EGEta2p1_DoubleJetCentral(Float_t EGcut, Float_t jetcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Int_t n2=0;

	Int_t Nele = myEvt_.Nele; 
	for (Int_t ue=0; ue < Nele; ue++) { 
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true; 
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;               
	for (Int_t ue=0; ue < Nj; ue++) {      
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) n2 ++;
	}

	jet = (n2 >= 2);

	Bool_t ok = (eg && jet);
	if (! ok) return false;

		//  -- now evaluate the delta condition :

	Bool_t CORREL = false;
	Int_t PhiOut[3];

	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt < EGcut) continue;

		Float_t phiel = myEvt_.Phiel[ue];
		Int_t iphiel = (int)phiel;

		PhiOut[0]=0; PhiOut[1]=0; PhiOut[2]=0; 

		PhiOut[0] = iphiel;
		if (iphiel< 17) PhiOut[1] = iphiel+1;
		if (iphiel == 17) PhiOut[1] = 0;
		if (iphiel > 0) PhiOut[2] = iphiel - 1;
		if (iphiel == 0) PhiOut[2] = 17;

		Int_t npair = 0;

		for (Int_t uj=0; uj < Nj; uj++) {
			Int_t bxj = myEvt_.Bxjet[uj];        		
			if (bxj != 0) continue; 
			Bool_t isFwdJet = myEvt_.Fwdjet[uj];
			if (isFwdJet) continue;
		   Bool_t isTauJet = myEvt_.Taujet[ue];
		   if (isTauJet) continue;		
			Float_t rankj = myEvt_.Etjet[uj];
										// Float_t ptj = rankj * 4;
			Float_t ptj = rankj; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[uj],rankj*4.,theL1JetCorrection);
			if (ptj < jetcut) continue;
			Float_t phijet = myEvt_.Phijet[uj];
			Int_t iphijet = (int)phijet;

			if ( iphijet != PhiOut[0] &&
				iphijet != PhiOut[1] &&
				iphijet != PhiOut[2] ) npair ++;

		}  // loop over jets

		if (npair >= 2 ) CORREL = true ;

	}  // end loop over EM objects

	return CORREL;

}

Bool_t L1Menu2015::EGEta2p1_DoubleJetCentral_TripleJetCentral(Float_t EGcut, Float_t jetcut2, Float_t jetcut3) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t eg = false;
	Bool_t jet = false;
	Int_t n2=0;       
	Int_t n3=0;

	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue; 
		Float_t eta = myEvt_.Etael[ue];
		if (eta < 4.5 || eta > 16.5) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= EGcut) eg = true;  
	}  // end loop over EM objects

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut2) n2 ++;
		if (pt >= jetcut3) n3 ++;
	}

	jet = (n2 >= 2 && n3 >= 3 );

	Bool_t ok = (eg && jet);
	return ok;

}

Bool_t L1Menu2015::HTT_HTM(Float_t HTTcut, Float_t HTMcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htt = false;
	Bool_t htm = false;
	Float_t adc = myEvt_.HTT;   
	Float_t TheHTT =  adc; // / 2.   ;          
	htt = ( TheHTT >= HTTcut ) ;

	Int_t adc_HTM  = myEvt_.HTM ; 
	Float_t TheHTM = adc_HTM; // * 2.  ;           
	htm = ( TheHTM >= HTMcut );

	Bool_t ok = (htt && htm);
	return ok;

}

Bool_t L1Menu2015::HTT_ETM(Float_t HTTcut, Float_t ETMcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htt = false;
	Bool_t etm = false;
	Float_t adc = myEvt_.HTT;   
	Float_t TheHTT =  adc; // / 2.   ;          
	htt = ( TheHTT >= HTTcut ) ;

	Int_t adc_ETM  = myEvt_.ETM ; 
	Float_t TheETM = adc_ETM; // * 2.  ;           
	etm = ( TheETM >= ETMcut );

	Bool_t ok = (htt && etm);
	return ok;

}

Bool_t L1Menu2015::Tau_ETM(Float_t taucut, Float_t ETMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (!isTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) jet = true;
	}

	Bool_t ok = ( jet && etm );
	return ok;
	
}


Bool_t L1Menu2015::IsoTau_ETM(Float_t taucut, Float_t ETMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isIsoTauJet = myEvt_.isoTaujet[ue];
		if (!isIsoTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) jet = true;
	}

	Bool_t ok = ( jet && etm );
	return ok;
	
}


Bool_t L1Menu2015::Tau_HTM(Float_t taucut, Float_t HTMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.HTM ;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (!isTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) jet = true;
	}

	Bool_t ok = ( jet && htm );
	return ok;
	
}


Bool_t L1Menu2015::IsoTau_HTM(Float_t taucut, Float_t HTMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t htm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.HTM ;
	Float_t TheHTM = adc; // / 2. ;
	htm = (TheHTM >= HTMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isIsoTauJet = myEvt_.isoTaujet[ue];
		if (!isIsoTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= taucut) jet = true;
	}

	Bool_t ok = ( jet && htm );
	return ok;
	
}



Bool_t L1Menu2015::Tau_JetCentral(Float_t taucut, Float_t jetcut, Float_t etaCut) {

// CORRUPT NEEDS UPDATING FOR NEW TAU APPROACH

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t tau = false;
	Bool_t jet = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
      Float_t eta= myEvt_.Etajet[ue];
		
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) {
		  if (pt >= taucut) tau = true;
		} else {
		  if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		  if (pt >= jetcut) jet = true;
		}
		
		
	}

	Bool_t ok = ( jet && tau );
	return ok;
	
}


Bool_t L1Menu2015::IsoTau_JetCentral(Float_t taucut, Float_t jetcut, Float_t etaCut) {

// CORRUPT NEEDS UPDATING FOR NEW TAU APPROACH

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t tau = false;
	Bool_t jet = false;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
      Float_t eta= myEvt_.Etajet[ue];
		
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) {
		  if (pt >= taucut && myEvt_.isoTaujet[ue]) tau = true;
		} else {
		  if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		  if (pt >= jetcut) jet = true;
		}
		
		
	}

	Bool_t ok = ( jet && tau );
	return ok;
	
}
Bool_t L1Menu2015::Tau_TwoJetForward(Float_t taucut, Float_t jetcut1, Float_t jetcut2) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t tau = false;
	Bool_t jet = false;
	
	Int_t njcut1 = 0;
	Int_t njcut2 = 0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		if (!myEvt_.Fwdjet[ue] && !myEvt_.Taujet[ue]) continue;		
		Float_t pt = myEvt_.Etjet[ue];
//      Float_t eta= myEvt_.Etajet[ue];
		
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) {
		  if (pt >= taucut) tau = true;
		} else {
		  if (pt >= jetcut1) njcut1++;
		  if (pt >= jetcut2) njcut2++;
		}
		
		
	}
 
   if(njcut1>=1 && njcut2>=2) jet = true;  //assumes jetcut1 >= jetcut2

	Bool_t ok = ( jet && tau );
	return ok;
	
}



Bool_t L1Menu2015::JetCentral_ETM(Float_t jetcut, Float_t ETMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false;
	Bool_t jet = false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut) jet = true;
	}

	Bool_t ok = ( jet && etm );
	return ok;
	
}

Bool_t L1Menu2015::DoubleJetCentral_ETM(Float_t jetcut1, Float_t jetcut2, Float_t ETMcut, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t etm = false; 
	Bool_t jet = false;
	Int_t n1=0;
	Int_t n2=0;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;
	etm = (TheETM >= ETMcut);

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= jetcut1) n1 ++;
		if (pt >= jetcut2) n2 ++;
	}       
	jet = (n1 >= 1 && n2 >= 2);

	Bool_t ok = ( jet && etm );
	return ok;

}

Bool_t L1Menu2015::SingleJetCentral(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut) ok = true;
	} 

	return ok;

}

Bool_t L1Menu2015::SingleJet(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut) ok = true;
	}

	return ok;

}

Bool_t L1Menu2015::DoubleJetCentral(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}
	Bool_t ok = ( n1 >=1 && n2 >= 2);
	//if(ok) printf("Run %i  Event  %i \n",event_->run,event_->event);
	return ok;

}

Bool_t L1Menu2015::DoubleJetForward(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (!isFwdJet) continue;
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}
	Bool_t ok = ( n1 >=1 && n2 >= 2);
	//if(ok) printf("Run %i  Event  %i \n",event_->run,event_->event);
	return ok;

}



Bool_t L1Menu2015::DoubleJet_Eta1p7_deltaEta4(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < 5.5 || eta > 15.5) continue;  // eta = 6 - 15
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}
	Bool_t ok = ( n1 >=1 && n2 >= 2);
	if (! ok) return false;

	// -- now the correlation

	Bool_t CORREL = false;

	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta1 = myEvt_.Etajet[ue];
		if (eta1 < 5.5 || eta1 > 15.5) continue;  // eta = 6 - 15
		if (pt < cut1) continue;

		for (Int_t ve=0; ve < Nj; ve++) {
			if (ve == ue) continue;
			Int_t bx2 = myEvt_.Bxjet[ve];        		
			if (bx2 != 0) continue;
			Bool_t isFwdJet2 = myEvt_.Fwdjet[ve];
			if (isFwdJet2) continue;
			Float_t rank2 = myEvt_.Etjet[ve];
			Float_t pt2 = rank2 * 4;
			Float_t eta2 = myEvt_.Etajet[ve];
			if (eta2 < 5.5 || eta2 > 15.5) continue;  // eta = 6 - 15
			if (pt2 < cut2) continue;

			Bool_t corr = correlateInEta((int)eta1, (int)eta2, 4);
			if (corr) CORREL = true;
		}


	}

	return CORREL ;

}


Bool_t L1Menu2015::SingleTauJet(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
      if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16  // eta = 5 - 16
		if (pt >= cut) n1++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 );
	return ok;

}


Bool_t L1Menu2015::SingleIsoTauJet(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isIsoTauJet = myEvt_.isoTaujet[ue];
		if (! isIsoTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
      if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16  // eta = 5 - 16
		if (pt >= cut) n1++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 );
	return ok;

}



Bool_t L1Menu2015::DoubleTauJetEta(Float_t cut1, Float_t cut2, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 && n2 >= 2);
	
// Print out results if there was a trigger
//   if(ok) {
//	  printf("========== Tau Trigger ============\n");
//	  dumpEvent(0x41);
//	  printf("===================================\n");
//	}	
	return ok;

}


Bool_t L1Menu2015::isoTau_Tau(Float_t cut1, Float_t cut2, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		if (pt >= cut1 && myEvt_.isoTaujet[ue]) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 && n2 >= 2);
	return ok;

}


Bool_t L1Menu2015::DoubleIsoTau(Float_t cut1, Float_t cut2, Float_t etaCut) {

	Bool_t raw = PhysicsBits[0];  // ZeroBias
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue; 
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (! isTauJet) continue;
		Float_t rank = myEvt_.Etjet[ue];    // the rank of the electron
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		Float_t eta = myEvt_.Etajet[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		if (pt >= cut1 && myEvt_.isoTaujet[ue]) n1++;
		if (pt >= cut2 && myEvt_.isoTaujet[ue]) n2++;
	}  // end loop over jets

	Bool_t ok = ( n1 >=1 && n2 >= 2);
	return ok;

}




Bool_t L1Menu2015::TripleJetCentral(Float_t cut1, Float_t cut2, Float_t cut3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
	}

	Bool_t ok = ( n1 >=1 && n2 >= 2 && n3 >= 3 );
	return ok;

}

Bool_t L1Menu2015::TripleJet_VBF(Float_t jet1, Float_t jet2, Float_t jet3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t jet=false;        
	Bool_t jetf1=false;           
	Bool_t jetf2=false;   

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;

	Int_t f1=0;
	Int_t f2=0;
	Int_t f3=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);

		if (isFwdJet) {
			if (pt >= jet1) f1 ++;
			if (pt >= jet2) f2 ++;
			if (pt >= jet3) f3 ++;              
		} 
		else {
			if (pt >= jet1) n1 ++;
			if (pt >= jet2) n2 ++;
			if (pt >= jet3) n3 ++;
		}    
	}

	jet   = ( n1 >= 1 && n2 >= 2 && n3 >= 3 ) ;        
	jetf1 = ( f1 >= 1 && n2 >= 1 && n3 >= 2 ) ;  // numbers change ofcourse    
	jetf2 = ( n1 >= 1 && f2 >= 1 && n3 >= 2 ) ;  

	Bool_t ok = false;

	if( jet || jetf1 || jetf2 ) ok =true;

	return ok;
}

Bool_t L1Menu2015::QuadJetCentral(Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4 ) {

// cut1 >= cut2  >= cut3 >= cut4

	// ZeroBias
// Bool_t raw = PhysicsBits[16];  // SingleJet36
// if (! raw) return false;
	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t n4=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
		if (pt >= cut4) n4++;
	}

	Bool_t ok = ( n1 >=1 && n2 >= 2 && n3 >= 3 && n4 >= 4);
	return ok;

}

Bool_t L1Menu2015::MultiJet(Int_t nj, Float_t cut1, Float_t cut2, Float_t cut3, Float_t cut4 ) {

// cut1 >= cut2  >= cut3 >= cut4

	// ZeroBias
// Bool_t raw = PhysicsBits[16];  // SingleJet36
// if (! raw) return false;
	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;


	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t n4=0;

	Int_t Nj = myEvt_.Njet ;
	for (Int_t ue=0; ue < Nj; ue++) {
		Int_t bx = myEvt_.Bxjet[ue];        		
		if (bx != 0) continue;
//		Bool_t isFwdJet = myEvt_.Fwdjet[ue];
//		if (isFwdJet) continue;
		Bool_t isTauJet = myEvt_.Taujet[ue];
		if (isTauJet) continue;		
		Float_t rank = myEvt_.Etjet[ue];
		Float_t pt = rank; //CorrectedL1JetPtByGCTregions(myEvt_.Etajet[ue],rank*4.,theL1JetCorrection);
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
		if (pt >= cut4) n4++;
	}

	Bool_t ok = ( n1 >=1 && n2 >= 2 && n3 >= 3 && n4 >= nj);
	return ok;

}



Bool_t L1Menu2015::ETM(Float_t ETMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Float_t adc = myEvt_.ETM ;
	Float_t TheETM = adc; // / 2. ;

	if (TheETM < ETMcut) return false;
	return true;

}

Bool_t L1Menu2015::HTM(Float_t HTMcut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Float_t adc = myEvt_.HTM ;
	Float_t TheHTM = adc; // / 2. ;

	if (TheHTM < HTMcut) return false;
	return true;

}


Bool_t L1Menu2015::HTT(Float_t HTTcut) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Float_t adc = myEvt_.HTT ;
	Float_t TheHTT = adc; // / 2. ;

	if (TheHTT < HTTcut) return false;
	return true;

}

Bool_t L1Menu2015::ETT(Float_t ETTcut) {

	Float_t adc = myEvt_.ETT ;
	Float_t TheETT = adc; // / 2. ;

	if (TheETT < ETTcut) return false;

	return true;

}


Bool_t L1Menu2015::SingleEG(Float_t cut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false; 
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ; 
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok; 

}

Bool_t L1Menu2015::SingleIsoEG_Eta(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Bool_t iso = myEvt_.Isoel[ue];
		if (! iso) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok;

}

Bool_t L1Menu2015::SingleEG_Eta(Float_t cut, Float_t etaCut ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Bool_t ok=false;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t eta = myEvt_.Etael[ue];
		if (eta < etaCut || eta > 21.-etaCut) continue;  // eta = 5 - 16
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut) ok = true;
	}  // end loop over EM objects

	return ok;

}

Bool_t L1Menu2015::DoubleEG(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2) ;
	//if(ok) printf("Found doubleEG event Run %i Event %i \n",event_->run,event_->event);
	return ok;

}

Bool_t L1Menu2015::DoubleIsoEG(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1 && myEvt_.Isoel[ue]) n1++;
		if (pt >= cut2 && myEvt_.Isoel[ue]) n2++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2) ;
	//if(ok) printf("Found doubleEG event Run %i Event %i \n",event_->run,event_->event);
	return ok;

}


Bool_t L1Menu2015::isoEG_EG(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1 && myEvt_.Isoel[ue]) n1++;
		if (pt >= cut2) n2++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2) ;
	//if(ok) printf("Found doubleEG event Run %i Event %i \n",event_->run,event_->event);
	return ok;

}

Bool_t L1Menu2015::EG_isoEG(Float_t cut1, Float_t cut2 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {               
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1) n1++;
		if (pt >= cut1 && myEvt_.Isoel[ue]) n2++;
		if (pt >= cut2 && myEvt_.Isoel[ue]) n3++;
	}  // end loop over EM objects

	Bool_t ok = false;
	if(n1==n2) { 
	   ok = ( n1 >= 1 && n3 >= 2) ; //if high Pt is isolated it will also increment n3 so demand n3>=2
	} else {
	   ok = ( n1 >= 1 && n3 >= 1) ; //If n1>n2 that means there is a nonisoEG above upper threshold so only demand n3>=1.
	}	
	//if(ok) printf("Found doubleEG event Run %i Event %i \n",event_->run,event_->event);
	return ok;

}


Bool_t L1Menu2015::TripleEG(Float_t cut1, Float_t cut2, Float_t cut3 ) {

	Bool_t raw = PhysicsBits[0];   // ZeroBias  
	if (! raw) return false;

	Int_t n1=0;
	Int_t n2=0;
	Int_t n3=0;
	Int_t Nele = myEvt_.Nele;
	for (Int_t ue=0; ue < Nele; ue++) {
		Int_t bx = myEvt_.Bxel[ue];        		
		if (bx != 0) continue;
		Float_t rank = myEvt_.Etel[ue];    // the rank of the electron
		Float_t pt = rank ;
		if (pt >= cut1) n1++;
		if (pt >= cut2) n2++;
		if (pt >= cut3) n3++;
	}  // end loop over EM objects

	Bool_t ok = ( n1 >= 1 && n2 >= 2 && n3 >= 3) ;
	return ok;


}

Bool_t L1Menu2015::UpgradeMuon(Float_t pt,Float_t eta, Float_t ptcut, Float_t etaCut){

    Bool_t fndMuon = false;
	 if(pt>=ptcut && fabs(eta)<=etaCut) fndMuon = true;
	 
	 return fndMuon;
}
	
Bool_t L1Menu2015::UpgradeIsolatedMuon(Float_t pt,Float_t eta, Float_t ptcut, Float_t etaCut) {

    Bool_t isoMuon = false;

// Ad Hoc rejection factor (Note rejection factor of 1. or less means all will pass)
    double rejectionFactor = 1.0;
    if(randGen->Uniform(1.) < 1./rejectionFactor) isoMuon = true;	 
	 Bool_t fndMuon = UpgradeMuon( pt, eta, ptcut, etaCut);
	 
	 return isoMuon & fndMuon;
}


void L1Menu2015::Loop(Int_t calcThreshold, Int_t selectDataInput, TString lsRunFile, TString L1MenuFile, const int n_events_) {

	
	const Int_t nevents = (n_events_ < 0) ? GetEntries() : n_events_;

	Int_t NPASS = 0; 

// Setup Random Number generator
   randGen = new TRandom2();


//	Int_t nPAG =0;
	first = true;
	int cnt=0;
   int lastLumi = -1;

// Load luminosity information
   int ind = 0;
   double Run[450], LS[450], IntL[450], PU[450], InstL[450];
	//TString lsRunFile = "getLumi_out_pixCorrLumi_66PU_stdCorr.txt";
   ifstream ifs( lsRunFile );
	while(ifs){
		ifs >> Run[ind];
		ifs >> LS[ind];
		ifs >> IntL[ind];
		ifs >> InstL[ind];
		ifs >> PU[ind];
		//printf("Run %f  LS %f  InstL %6.3e \n",Run[ind],LS[ind],InstL[ind]);
		ind++;
	}




	for (Long64_t i=0; i<nevents; i++)
	{     
	//load the i-th event
		Long64_t ientry = LoadTree(i); if (ientry < 0) break;
		GetEntry(i);

      cnt++;
      if((cnt%100)==0) printf("Event Number %i\r",cnt); fflush(stdout);
//      if(cnt%(int)pow(10.,(double)((int)log10((double)cnt)))==0) printf("Event Number %i\r",cnt); fflush(stdout);


// Fill my event data
      fillDataStructure(selectDataInput);

// Dump events
      if(cnt<numDumpEvt) dumpEvent();


//      Fill the physics bits:
                //printf("Entry %i",i);
		FilL1Bits();

		if (first) InitL1Menu(L1MenuFile);
		
		//
		//PhysicsBits[0] = true; //force it!
		Bool_t raw = PhysicsBits[0];  // ZeroBias
		if (! raw) continue;


//  --- Reset the emulated "Trigger Bits"
		kOFFSET = 0;
		for (Int_t k=0; k < N128; k++) {
			TheTriggerBits[k] = false;
		}


// Get the instantaneous luminosity for this lumiSection (defaults at set in RunL1_HFW)
			if (event_->run == 198588) theZeroBiasPrescale = 44.;
			if (event_->run == 198603) theZeroBiasPrescale = 92.;
			if (event_->run == 198609) theZeroBiasPrescale = 92.;
		   int thisLumi = event_->lumi; //get this event's lumi
			if(thisLumi != lastLumi) {
			  // Note: If LS not found in list, value of inst. lumi will not be changed from default set in RunL1_HFW.
 			  for(int k=0; k<ind-1; k++) if( LS[k]==thisLumi && Run[k]==event_->run) theLumiForThisSetOfLumiSections = InstL[k]/1.0e32; 		   
			  lastLumi = thisLumi;
			}         

/*  Determine the Weight for this event to turn it into a contribution to a rate

  Data:
        scal = theZeroBiasPrescale * puWeight                         <-- All Prescales (L1 and any HLT if imposed) also PU weighting of MC (data puWeight = 1.)
		         (theTargetLumi/theLumiForThisSetOfLumiSections) /      <-- Ratio of luminosity
					(23.3570304 * theNumberOfUserdLumiSections) /          <-- total time
					 1000.                                                 <-- turn Hz to kHz 

  MC:
  MC:
        Kind of ad hoc at the moment.  For MC input the following:
		        ZeroBiasPrescale = 1
				  TargetLumi = LumiForThisSetOfLumiSections
				  
				  theNumberOfUserLumiSections = (num MC Events)/(23.3570304*11246.*Number of Bunches)
				  
		  		  This approach treats each  MC as a crossing
*/
// For test L1 Upgrade Triggers all are collected in one place.
      double lumiWeight = theZeroBiasPrescale * puWeight *                          
                          (theTargetLumi/theLumiForThisSetOfLumiSections) /  
                          (23.3570304 * theNumberOfUserdLumiSections) / 1000.;
		Bool_t pass = EvalMenu(lumiWeight);
		if(calcThreshold>0) EvalThresh(calcThreshold,lumiWeight);

		if (pass) NPASS ++;

// We have done the first event...used elsewhere
		first = false;

	// -- now the pure rate stuff
	// -- kOFFSET now contains the number of triggers we have calculated
      Bool_t anyTrigger = false;
		Bool_t firstTrig  = true;
		
		for (Int_t k=0; k < kOFFSET; k++) {
			if ( ! TheTriggerBits[k] ) continue;
		
		   anyTrigger = true;
			h_All -> Fill(k,lumiWeight);

		// did the event pass another trigger ?
			Bool_t nonpure = false;
			for (Int_t k2=0; k2 < kOFFSET; k2++) {
				if (k2 == k) {
				  h_Corr->Fill(k,k,lumiWeight);
				  //continue;
				} else if ( TheTriggerBits[k2] ) {
				  nonpure = true;
				  h_Corr->Fill(k,k2,lumiWeight);
				}	
			}
			Bool_t pure = !nonpure ;
			if (pure) {
			  h_Pure -> Fill(k,lumiWeight);
			  h_Pure -> Fill(kOFFSET,lumiWeight);
			} 
			
// Calculate the cummulative rates
         if(firstTrig) {
			   for(Int_t k2=k; k2<kOFFSET; k2++) h_Cumm->Fill(k2,lumiWeight);
				firstTrig = false;
			}
			
			
		}
		if(anyTrigger) h_All->Fill(kOFFSET,lumiWeight);
		h_All -> GetXaxis() -> SetBinLabel(kOFFSET+1,"Triggered Evt");
		h_Pure -> GetXaxis() -> SetBinLabel(kOFFSET+1,"Triggered Evt");

	}  // end evt loop


//	std::cout << " Prescales for: " << theTargetLumi << ", LumiForThisSetOfLumiSections = " << theLumiForThisSetOfLumiSections << ", L1NtupleFileName = " << theL1NtupleFileName << std::endl;
//	std::cout << std::endl << " --------------------------------------------------------- " << std::endl << std::endl;

/*	
// Loop over triggers and print out configuration
//=================================================
        for(std::map<std::string, trigPar>::iterator itr = trigParList.begin(); itr != trigParList.end(); itr++) {
	   std::cout  << setw(20) << itr->first << "   Parameters(" << setw(5) << (itr->second).primTh << "," 
	                                                << setw(5) << (itr->second).secTh << ","
						        << setw(5) << (itr->second).triTh << ","
						        << setw(5) << (itr->second).quadTh << ","
						        << setw(5) << (itr->second).etaCut << ","
						        << setw(5) << (itr->second).minQual << ")"  << std::endl;
	} 
	std::cout << std::endl;
*/
        	
}



Float_t L1Menu2015::RunL1Ana(
     TString L1MenuFileName,                //File name containing L1 Menu Algorithms and thresholds. (See InitL1Menu for format)
	  TString lsFileName,	                 //LumiSection Luminosity (used for data)
	  TString jobTag,                        //Job tag for storing output
	  Int_t   procNevts,                     //Number of events to run over.  If -1, run over all events in the file
	  Int_t   makeThresholdPlots,            //Flag for whether to calculate the rate vs threshold plots (0=skip all; 1=do 1-D plots; 2=do 1-d and 2-d plots)
	  Int_t  selectDataInput	                    //Flag for whether to use L1Extra values
                ) {
	
	

					 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	std::stringstream histos;
	histos << "L1RateHist_" << jobTag;
	if(makeThresholdPlots>0) histos << "Thr" << makeThresholdPlots; 
	histos << "_rates.root";
   TString outHistName = histos.str();
	TFile* outHist = new TFile(outHistName,"RECREATE");
	outHist->cd();


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Make sure we get the errors correct on histograms   
	TH1::SetDefaultSumw2();

//Cross
	h_SingleMu_ETM_byThreshold     = new TH1F("h_SingleMu_ETM_byThreshold","h_SingleMu_ETM_byThreshold",141,-0.5,140.5);
	h2_SingleMu_ETM_byThreshold    = new TH2F("h2_SingleMu_ETM_byThreshold","h2_SingleMu_ETM_byThreshold",141,-0.5,140.5,201,-0.25,100.25);
	h_SingleMu_CJet_byThreshold    = new TH1F("h_SingleMu_CJet_byThreshold","h_SingleMu_CJet_byThreshold",141,-0.5,140.5);
	h2_SingleMu_CJet_byThreshold   = new TH2F("h2_SingleMu_CJet_byThreshold","h2_SingleMu_CJet_byThreshold",141,-0.5,140.5,51,-2.,202.);
	h_SingleIsoMu_ETM_byThreshold  = new TH1F("h_SingleIsoMu_ETM_byThreshold","h_SingleIsoMu_ETM_byThreshold",141,-0.5,140.5);
	h2_SingleIsoMu_ETM_byThreshold = new TH2F("h2_SingleIsoMu_ETM_byThreshold","h2_SingleIsoMu_ETM_byThreshold",141,-0.5,140.5,201,-0.25,100.25);
	h_SingleIsoMu_CJet_byThreshold = new TH1F("h_SingleIsoMu_CJet_byThreshold","h_SingleIsoMu_CJet_byThreshold",141,-0.5,140.5);
	h2_SingleIsoMu_CJet_byThreshold= new TH2F("h2_SingleIsoMu_CJet_byThreshold","h2_SingleIsoMu_CJet_byThreshold",141,-0.5,140.5,51,-2.,202.);
	h_SingleMu_HTM_byThreshold     = new TH1F("h_SingleMu_HTM_byThreshold","h_SingleMu_HTM_byThreshold",141,-0.5,140.5);
	h_SingleIsoMu_HTM_byThreshold  = new TH1F("h_SingleIsoMu_HTM_byThreshold","h_SingleIsoMu_HTM_byThreshold",141,-0.5,140.5);
	h2_SingleMu_HTM_byThreshold    = new TH2F("h2_SingleIMu_HTM_byThreshold","h2_SingleMu_HTM_byThreshold",141,-0.5,140.5,201,-0.25,100.25);
	h2_SingleIsoMu_HTM_byThreshold = new TH2F("h2_SingleIsoMu_HTM_byThreshold","h2_SingleIsoMu_HTM_byThreshold",141,-0.5,140.5,201,-0.25,100.25);


	h_SingleEG_ETM_byThreshold     = new TH1F("h_SingleEG_ETM_byThreshold","h_SingleEG_ETM_byThreshold",64,-0.5,63.5);
	h2_SingleEG_ETM_byThreshold    = new TH2F("h2_SingleEG_ETM_byThreshold","h2_SingleEG_ETM_byThreshold",64,-0.5,63.5,201,-0.25,100.25);
	h_SingleEG_CJet_byThreshold    = new TH1F("h_SingleEG_CJet_byThreshold","h_SingleEG_CJet_byThreshold",64,-0.5,63.5);
	h2_SingleEG_CJet_byThreshold   = new TH2F("h2_SingleEG_CJet_byThreshold","h2_SingleEG_CJet_byThreshold",64,-0.5,63.5,51,-2.,202.);
	h_SingleIsoEG_ETM_byThreshold  = new TH1F("h_SingleIsoEG_ETM_byThreshold","h_SingleIsoEG_ETM_byThreshold",64,-0.5,63.5);
	h2_SingleIsoEG_ETM_byThreshold = new TH2F("h2_SingleIsoEG_ETM_byThreshold","h2_SingleIsoEG_ETM_byThreshold",64,-0.5,63.5,201,-0.25,100.25);
	h_SingleIsoEG_CJet_byThreshold = new TH1F("h_SingleIsoEG_CJet_byThreshold","h_SingleIsoEG_CJet_byThreshold",64,-0.5,63.5);
	h2_SingleIsoEG_CJet_byThreshold= new TH2F("h2_SingleIsoEG_CJet_byThreshold","h2_SingleIsoEG_CJet_byThreshold",64,-0.5,63.5,51,-2.,202.);
	h_SingleEG_HTM_byThreshold     = new TH1F("h_SingleEG_HTM_byThreshold","h_SingleEG_HTM_byThreshold",64,-0.5,63.5);
	h_SingleIsoEG_HTM_byThreshold  = new TH1F("h_SingleIsoEG_HTM_byThreshold","h_SingleIsoEG_HTM_byThreshold",64,-0.5,63.5);
	h2_SingleEG_HTM_byThreshold    = new TH2F("h2_SingleEG_HTM_byThreshold","h2_SingleEG_HTM_byThreshold",64,-0.5,63.5,201,-0.25,100.25);
	h2_SingleIsoEG_HTM_byThreshold = new TH2F("h2_SingleIsoEG_HTM_byThreshold","h2_SingleIsoEG_HTM_byThreshold",64,-0.5,63.5,201,-0.25,100.25);



	h_Mu_EG_byThreshold        = new TH1F("h_Mu_EG_byThreshold","h_Mu_EG_byThreshold",141,-0.5,140.5);
   h2_Mu_EG_byThreshold       = new TH2F("h2_Mu_EG_byThreshold","h2_Mu_EG_byThreshold",141,-0.5,140.5,64,-0.5,63.5);
	h_isoMu_EG_byThreshold     = new TH1F("h_isoMu_EG_byThreshold","h_isoMu_EG_byThreshold",141,-0.5,140.5);
   h2_isoMu_EG_byThreshold    = new TH2F("h2_isoMu_EG_byThreshold","h2_isoMu_EG_byThreshold",141,-0.5,140.5,64,-0.5,63.5);
	h_isoMu_isoEG_byThreshold  = new TH1F("h_isoMu_isoEG_byThreshold","h_isoMu_isoEG_byThreshold",141,-0.5,140.5);
   h2_isoMu_isoEG_byThreshold = new TH2F("h2_isoMu_isoEG_byThreshold","h2_isoMu_isoEG_byThreshold",141,-0.5,140.5,64,-0.5,63.5);

	h_EG_Mu_byThreshold        = new TH1F("h_EG_Mu_byThreshold","h_EG_Mu_byThreshold",64,-0.5,63.5);
	h_isoEG_Mu_byThreshold     = new TH1F("h_isoEG_Mu_byThreshold","h_isoEG_Mu_byThreshold",64,-0.5,63.5);
   h2_isoEG_Mu_byThreshold    = new TH2F("h2_isoEG_Mu_byThreshold","h2_isoEG_Mu_byThreshold",64,-0.5,63.5,64,-0.5,64.5);
	h_isoEG_isoMu_byThreshold  = new TH1F("h_isoEG_isoMu_byThreshold","h_isoEG_isoMu_byThreshold",64,-0.5,63.5);
   h2_isoEG_isoMu_byThreshold = new TH2F("h2_isoEG_isoMu_byThreshold","h2_isoEG_isoMu_byThreshold",64,-0.5,63.5,64,-0.5,64.5);

	h_EG_Tau_byThreshold       = new TH1F("h_EG_Tau_byThreshold","h_EG_Tau_byThreshold",64,-0.5,63.5);
   h2_EG_Tau_byThreshold      = new TH2F("h2_EG_Tau_byThreshold","h2_EG_Tau_byThreshold",64,-0.5,63.5,51,-2.,202.);
	h_isoEG_Tau_byThreshold    = new TH1F("h_isoEG_Tau_byThreshold","h_isoEG_Tau_byThreshold",64,-0.5,63.5);
   h2_isoEG_Tau_byThreshold   = new TH2F("h2_isoEG_Tau_byThreshold","h2_isoEG_Tau_byThreshold",64,-0.5,63.5,51,-2.,202.);
	h_isoEG_isoTau_byThreshold = new TH1F("h_isoEG_isoTau_byThreshold","h_isoEG_isoTau_byThreshold",64,-0.5,63.5);
   h2_isoEG_isoTau_byThreshold= new TH2F("h2_isoEG_isoTau_byThreshold","h2_isoEG_isoTau_byThreshold",64,-0.5,63.5,51,-2.,202.);

	h_Mu_Tau_byThreshold        = new TH1F("h_Mu_Tau_byThreshold","h_Mu_Tau_byThreshold",141,-0.5,140.5);
   h2_Mu_Tau_byThreshold       = new TH2F("h2_Mu_Tau_byThreshold","h2_Mu_Tau_byThreshold",141,-0.5,140.5,51,-2.,202.);
	h_isoMu_Tau_byThreshold     = new TH1F("h_isoMu_Tau_byThreshold","h_isoMu_Tau_byThreshold",141,-0.5,140.5);
   h2_isoMu_Tau_byThreshold    = new TH2F("h2_isoMu_Tau_byThreshold","h2_isoMu_Tau_byThreshold",141,-0.5,140.5,51,-2.,202.);
	h_isoMu_isoTau_byThreshold  = new TH1F("h_isoMu_isoTau_byThreshold","h_isoMu_isoTau_byThreshold",141,-0.5,140.5);
   h2_isoMu_isoTau_byThreshold = new TH2F("h2_isoMu_isoTau_byThreshold","h2_isoMu_isoTau_byThreshold",141,-0.5,140.5,51,-2.,202.);


	
//Jets
	h_SingleJet_byThreshold       = new TH1F("h_SingleJet_byThreshold","h_SingleJet_byThreshold",101,-2.,402.);
	h_SingleJetC_byThreshold      = new TH1F("h_SingleJetC_byThreshold","h_SingleJetC_byThreshold",101,-2.,402.);	
   h_DoubleJet_byThreshold       = new TH1F("h_DoubleJet_byThreshold","h_DoubleJet_byThreshold",101,-2.,402.);
	h2_DoubleJet_byThreshold      = new TH2F("h2_DoubleJet_byThreshold","h2_DoubleJet_byThreshold",51,-2.,202.,51,-2.,202.);
   h_DoubleFwdJet_byThreshold    = new TH1F("h_DoubleFwdJet_byThreshold","h_DoubleFwdJet_byThreshold",101,-2.,402.);
	h2_DoubleFwdJet_byThreshold   = new TH2F("h2_DoubleFwdJet_byThreshold","h2_DoubleFwdJet_byThreshold",51,-2.,202.,51,-2.,202.);
	h_QuadJetC_byThreshold        = new TH1F("h_QuadJetC_byThreshold","h_QuadJetC_byThreshold",101,-2.,402.);
	h2A_QuadJetCentral_byThreshold= new TH2F("h2A_QuadJetCentral_byThreshold","h2A_QuadJetCentral_byThreshold (1,3)",51,-2.,202.,51,-2.,202.);
	h2B_QuadJetCentral_byThreshold= new TH2F("h2B_QuadJetCentral_byThreshold","h2B_QuadJetCentral_byThreshold (2,2)",51,-2.,202.,51,-2.,202.);
	h_SixJet_byThreshold          = new TH1F("h_SixJet_byThreshold","h_SixJet_byThreshold",101,-2.,402.);
	h2A_SixJet_byThreshold        = new TH2F("h2A_SixJet_byThreshold","h2A_SixJet_byThreshold (1-2,3-6)",51,-2.,202.,51,-2.,202.);
	h_SingleCJet_ETM_byThreshold  = new TH1F("h_SingleCJet_ETM_byThreshold","h_SingleCJet_ETM_byThreshold",101,-2.,402.);
	h2_SingleCJet_ETM_byThreshold = new TH2F("h2_SingleCJet_ETM_byThreshold","h2_SingleCJet_ETM_byThreshold",101,-2.,402.,201,-0.25,100.25);
	h_DoubleCJet_ETM_byThreshold  = new TH1F("h_DoubleCJet_ETM_byThreshold","h_DoubleCJet_ETM_byThreshold",101,-2.,402.);
	h2_DoubleCJet_ETM_byThreshold = new TH2F("h2_DoubleCJet_ETM_byThreshold","h2_DoubleCJet_ETM_byThreshold",101,-2.,402.,201,-0.25,100.25);

// Taus (4 GeV Steps)
//	h_SingleTau_byThreshold          = new TH1F("h_SingleTau_byThreshold","h_SingleTau_byThreshold",101,-2.,402.);
//	h_SingleIsoTau_byThreshold       = new TH1F("h_SingleIsoTau_byThreshold","h_SingleIsoTau_byThreshold",101,-2.,402.);
//	h_DoubleTau_byThreshold          = new TH1F("h_DoubleTau_byThreshold","h_DoubleTau_byThreshold",101,-2.,402.);
//	h2_DoubleTau_byThreshold         = new TH2F("h2_DoubleTau_byThreshold","h2_DoubleTau_byThreshold",51,-2.,202.,51,-2.,202.);
//	h_SingleTau_ETM_byThreshold      = new TH1F("h_SingleTau_ETM_byThreshold","h_SingleTau_ETM_byThreshold",51,-2.,202.);
//	h2_SingleTau_ETM_byThreshold     = new TH2F("h2_SingleTau_ETM_byThreshold","h2_SingleTau_ETM_byThreshold",51,-2.,202.,201,-0.25,100.25);	
//	h_SingleTau_CJet_byThreshold     = new TH1F("h_SingleTau_CJet_byThreshold","h_SingleTau_CJet_byThreshold",101,-2.,402.);
//	h2_SingleTau_CJet_byThreshold    = new TH2F("h2_SingleTau_CJet_byThreshold","h2_SingleTau_CJet_byThreshold",51,-2.,202.,51,-2.,202.);
//	h_SingleTau_TwoFJet_byThreshold  = new TH1F("h_SingleTau_TwoFJet_byThreshold","h_SingleTau_TwoFJet_byThreshold",101,-2.,402.);
//	h2_SingleTau_TwoFJet_byThreshold = new TH2F("h2_SingleTau_TwoFJet_byThreshold","h2_SingleTau_TwoFJet_byThreshold",51,-2.,202.,51,-2.,202.);

// Taus
	h_SingleTau_byThreshold          = new TH1F("h_SingleTau_byThreshold","h_SingleTau_byThreshold",201,-0.5,200.5);
	h_SingleIsoTau_byThreshold       = new TH1F("h_SingleIsoTau_byThreshold","h_SingleIsoTau_byThreshold",201,-0.5,200.5);
	h_DoubleTau_byThreshold          = new TH1F("h_DoubleTau_byThreshold","h_DoubleTau_byThreshold",201,-0.5,200.5);
	h2_DoubleTau_byThreshold         = new TH2F("h2_DoubleTau_byThreshold","h2_DoubleTau_byThreshold",201,-0.5,200.5,201,-0.5,200.5);
	h_DoubleIsoTau_byThreshold       = new TH1F("h_DoubleIsoTau_byThreshold","h_DoubleIsoTau_byThreshold",201,-0.5,200.5);
	h2_DoubleIsoTau_byThreshold      = new TH2F("h2_DoubleIsoTau_byThreshold","h2_DoubleIsoTau_byThreshold",201,-0.5,200.5,201,-0.5,200.5);
	h_isoTau_Tau_byThreshold         = new TH1F("h_isoTau_Tau_byThreshold","h_isoTau_Tau_byThreshold",201,-0.5,200.5);
	h_SingleTau_ETM_byThreshold      = new TH1F("h_SingleTau_ETM_byThreshold","h_SingleTau_ETM_byThreshold",201,-0.5,200.5);
	h2_SingleTau_ETM_byThreshold     = new TH2F("h2_SingleTau_ETM_byThreshold","h2_SingleTau_ETM_byThreshold",201,-0.5,200.5,201,-0.25,100.25);	
	h_SingleIsoTau_ETM_byThreshold   = new TH1F("h_SingleIsoTau_ETM_byThreshold","h_SingleIsoTau_ETM_byThreshold",201,-0.5,200.5);
	h_SingleTau_HTM_byThreshold      = new TH1F("h_SingleTau_HTM_byThreshold","h_SingleTau_HTM_byThreshold",201,-0.5,200.5);
	h_SingleIsoTau_HTM_byThreshold   = new TH1F("h_SingleIsoTau_HTM_byThreshold","h_SingleIsoTau_HTM_byThreshold",201,-0.5,200.5);
	h_SingleTau_CJet_byThreshold     = new TH1F("h_SingleTau_CJet_byThreshold","h_SingleTau_CJet_byThreshold",201,-0.5,200.5);
	h2_SingleTau_CJet_byThreshold    = new TH2F("h2_SingleTau_CJet_byThreshold","h2_SingleTau_CJet_byThreshold",201,-0.5,200.5,51,-2.,202.);
	h_SingleIsoTau_CJet_byThreshold  = new TH1F("h_SingleIsoTau_CJet_byThreshold","h_SingleIsoTau_CJet_byThreshold",201,-0.5,200.5);
	h_SingleTau_TwoFJet_byThreshold  = new TH1F("h_SingleTau_TwoFJet_byThreshold","h_SingleTau_TwoFJet_byThreshold",201,-0.5,200.5);
	h2_SingleTau_TwoFJet_byThreshold = new TH2F("h2_SingleTau_TwoFJet_byThreshold","h2_SingleTau_TwoFJet_byThreshold",201,-0.5,200.5,51,-2.,202.);



//Sums 
	h_HTT_byThreshold = new TH1F("h_HTT_byThreshold","h_HTT_byThreshold",1601,-0.25,800.25);
	h_ETM_byThreshold = new TH1F("h_ETM_byThreshold","h_ETM_byThreshold",201 ,-0.25,100.25);
	h_HTM_byThreshold = new TH1F("h_HTM_byThreshold","h_HTM_byThreshold",201 ,-0.5,200.5);	
   h_HTT_ETM_byThreshold = new TH1F("h_HTT_ETM_byThreshold","h_HTT_ETM_byThreshold",1601,-0.25,800.25);
   h2_HTT_ETM_byThreshold = new TH2F("h2_HTT_ETM_byThreshold","h_HTT_ETM_byThreshold",401,-0.25,200.25,201,-0.25,100.25);

//EGamma
	h_SingleEG_byThreshold    = new TH1F("h_SingleEG_byThreshold","h_SingleEG_byThreshold",64,-0.5,63.5);
	h_SingleIsoEG_byThreshold = new TH1F("h_SingleIsoEG_byThreshold","h_SingleIsoEG_byThreshold",64,-0.5,63.5);
	h_DoubleEG_byThreshold    = new TH1F("h_DoubleEG_byThreshold","h_DoubleEG_byThreshold",64,-0.5,63.5);
	h2_DoubleEG_byThreshold   = new TH2F("h2_DoubleEG_byThreshold","h2_DoubleEG_byThreshold",64,-0.5,63.5,64,-0.5,63.5);
	h_DoubleIsoEG_byThreshold = new TH1F("h_DoubleIsoEG_byThreshold","h_DoubleIsoEG_byThreshold",64,-0.5,63.5);
	h2_DoubleIsoEG_byThreshold= new TH2F("h2_DoubleIsoEG_byThreshold","h2_DoubleIsoEG_byThreshold",64,-0.5,63.5,64,-0.5,63.5);
	h_isoEG_EG_byThreshold    = new TH1F("h_isoEG_EG_byThreshold","h_isoEG_EG_byThreshold",64,-0.5,63.5);
	h_EG_isoEG_byThreshold    = new TH1F("h_EG_isoEG_byThreshold","h_EG_isoEG_byThreshold",64,-0.5,63.5);

//Muons
	h_SingleMu_byThreshold    = new TH1F("h_SingleMu_byThreshold","h_SingleMu_byThreshold",141,-0.5,140.5);
	h_SingleIsoMu_byThreshold = new TH1F("h_SingleIsoMu_byThreshold","h_SingleIsoMu_byThreshold",141,-0.5,140.5);	
	h_DoubleMu_byThreshold    = new TH1F("h_DoubleMu_byThreshold","h_DoubleMu_byThreshold",141,-0.5,140.5);
	h2_DoubleMu_byThreshold   = new TH2F("h2_DoubleMu_byThreshold","h2_DoubleMu_byThreshold",141,-0.5,140.5,141,-0.5,140.5);
	h_DoubleIsoMu_byThreshold = new TH1F("h_DoubleIsoMu_byThreshold","h_DoubleIsoMu_byThreshold",141,-0.5,140.5);
	h2_DoubleIsoMu_byThreshold= new TH2F("h2_DoubleIsoMu_byThreshold","h2_DoubleIsoMu_byThreshold",141,-0.5,140.5,141,-0.5,140.5);
	h_isoMu_Mu_byThreshold    = new TH1F("h_isoMu_Mu_byThreshold","h_isoMu_Mu_byThreshold",141,-0.5,140.5);
	h2_isoMu_Mu_byThreshold   = new TH2F("h2_isoMu_Mu_byThreshold","h2_isoMu_Mu_byThreshold",141,-0.5,140.5,141,-0.5,140.5);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Plot of input quantities
   h_Mu_Nmu = new TH1F("h_Mu_Nmu","Num Mu",5,-0.5,4.5);
   h_Mu_Et  = new TH1F("h_Mu_Et","Mu Pt",141,-0.5,140.5);
   h_Mu_Eta = new TH1F("h_Mu_Eta","Mu Eta",100,-2.5,2.5);
   h_Mu_Phi = new TH1F("h_Mu_Phi","Mu Phi",100,0.0,TMath::TwoPi());

   h_isoMu_Nmu = new TH1F("h_isoMu_Nmu","Num isoMu",5,-0.5,4.5);
   h_isoMu_Et  = new TH1F("h_isoMu_Et","isoMu Pt",141,-0.5,140.5);
   h_isoMu_Eta = new TH1F("h_isoMu_Eta","isoMu Eta",100,-2.5,2.5);
   h_isoMu_Phi = new TH1F("h_isoMu_Phi","isoMu Phi",100,0.0,TMath::TwoPi());
	
   h_isoEG_Nele= new TH1F("h_isoEG_Nele","Num iso EG",10,-0.5,9.5);
   h_isoEG_Et  = new TH1F("h_isoEG_Et","iso EG Et",64,-0.5,63.5);
   h_isoEG_Eta = new TH1F("h_isoEG_Eta","iso EG Eta",22,-0.5,21.5);
   h_isoEG_Phi = new TH1F("h_isoEG_Phi","iso EG Phi",18,-0.5,17.5);

   h_nIsoEG_Nele= new TH1F("h_nIsoEG_Nele","Num nonIso EG",10,-0.5,9.5);
   h_nIsoEG_Et  = new TH1F("h_nIsoEG_Et","nonIso EG Et",64,-0.5,63.5);
   h_nIsoEG_Eta = new TH1F("h_nIsoEG_Eta","nonIso EG Eta",22,-0.5,21.5);
   h_nIsoEG_Phi = new TH1F("h_nIsoEG_Phi","nonIso EG Phi",18,-0.5,17.5);
	
   h_CJet_Njet= new TH1F("h_CJet_Njet","Num Central Jets",10,-0.5,9.5);
   h_CJet_Et  = new TH1F("h_CJet_Et","Central Jet Et",257,-0.5,256.5);
   h_CJet_Eta = new TH1F("h_CJet_Eta","Central Jet Eta",22,-0.5,21.5);
   h_CJet_Phi = new TH1F("h_CJet_Phi","Central Jet Phi",18,-0.5,17.5);

   h_FJet_Njet= new TH1F("h_FJet_Njet","Num Forward Jets",10,-0.5,9.5);
   h_FJet_Et  = new TH1F("h_FJet_Et","Forward Jet Et",257,-0.5,256.5);
   h_FJet_Eta = new TH1F("h_FJet_Eta","Forward Jet Eta",22,-0.5,21.5);
   h_FJet_Phi = new TH1F("h_FJet_Phi","Forward Jet Phi",18,-0.5,17.5);

   h_TJet_Njet= new TH1F("h_TJet_Njet","Num Tau Jets",10,-0.5,9.5);
	h_TJet_Et  = new TH1F("h_TJet_Et","Tau Jet Et",257,-0.5,256.5);
   h_TJet_Eta = new TH1F("h_TJet_Eta","Tau Jet Eta",22,-0.5,21.5);
   h_TJet_Phi = new TH1F("h_TJet_Phi","Tau Jet Phi",18,-0.5,17.5);

   h_isoTJet_Njet= new TH1F("h_isoTJet_Njet","Num isoTau Jets",10,-0.5,9.5);
	h_isoTJet_Et  = new TH1F("h_isoTJet_Et","isoTau Jet Et",257,-0.5,256.5);
   h_isoTJet_Eta = new TH1F("h_isoTJet_Eta","isoTau Jet Eta",22,-0.5,21.5);
   h_isoTJet_Phi = new TH1F("h_isoTJet_Phi","isoTau Jet Phi",18,-0.5,17.5);
	
   h_Sum_ETT   = new TH1F("h_Sum_ETT","ET Total",1601,-0.25,800.25);
	h_Sum_ETM   = new TH1F("h_Sum_ETM","ET Miss",201,-0.25,100.25);
	h_Sum_PhiETM= new TH1F("h_Sum_PhiETM","ET Miss Phi",100,-TMath::Pi(),TMath::Pi());	

   h_Sum_HTT   = new TH1F("h_Sum_HTT","HT Total",1601,-0.25,800.25);
	h_Sum_HTM   = new TH1F("h_Sum_HTM","HT Miss",101,-0.25,200.25);
	h_Sum_PhiHTM= new TH1F("h_Sum_PhiHTM","HT Miss Phi",100,-TMath::Pi(),TMath::Pi());


// Histograms for the Trigger Results.
   h_Trig = new TH1F("h_Trig","h_Trig",N128,-0.5,N128+0.5);
	h_All  = new TH1F("h_All","h_All",N128,-0.5,N128+0.5);
	h_Pure = new TH1F("h_Pure","h_Pure",N128,-0.5,N128+0.5);
	h_Corr = new TH2F("h_Corr","h_Corr",N128,-0.5,N128+0.5,N128,-0.5,N128+0.5);
	h_Cumm = new TH1F("h_Cumm","h_Cumm",N128,-0.5,N128+0.5);
	
	
// Histogram Event Level Quantities
   h_puWeight    = new TH1F("h_puWeight","PU Weight",100,0.,5.);
	h_Sim_meanInt = new TH1F("h_Sim_meanInt","Mean Int (sim)",100,0.,100.);
	h_lumiSec     = new TH1F("h_lumiSec","Lumi Section",120,0.,1200.);
	h_bx          = new TH1F("h_bx","Beam Crossing",100,0.,3000.);	

// Time to dump the configuration 
   std::cout << std::endl;
   std::cout << " ==============  Running L1Menu2015 with Parameters  ===================" << std::endl;
   std::cout << "  L1 menu                       = " << L1MenuFileName.Data() << std::endl; 
	std::cout << "  L1NtupleFileName              = " << theL1NtupleFileName << std::endl;	
	std::cout << "  LS File Name                  = " << lsFileName.Data() << std::endl;
	std::cout << "  Zero Bias Prescales           = " << theZeroBiasPrescale << std::endl;		
	std::cout << "  NumberOfUserdLumiSections     = " << theNumberOfUserdLumiSections << std::endl;
	std::cout << "  LumiForThisSetOfLumiSections  = " << theLumiForThisSetOfLumiSections << std::endl;
	std::cout << "  Target Luminosity             = " << theTargetLumi << std::endl;	
	std::cout << "  job Tag                       = " << jobTag.Data() << std::endl;	
	std::cout << "  Number of Events to Process   = " << procNevts << std::endl;	
   std::cout << "  Calculate Threshold Plots     = " << makeThresholdPlots  << std::endl;
	std::cout << "  Data input Selection          = " << selectDataInput << std::endl;
   std::cout << " =======================================================================" << std::endl;
   std::cout << std::endl;

//  Do the heavy lifting
//	L1Menu2015 a(usedL1MenuThr,targetlumi,NumberOfUserLumiSections,LumiForThisSetOfLumiSections,
//			       (std::string)L1NtupleFileName,AveragePU,ZeroBiasPrescale,L1JetCorrection);
//	a.Open((std::string)L1NtupleFileName);
//	a.Loop(makeThresholdPlots,selectDataInput,lsFileName,L1MenuFileName,procNevts);

//  Do the heavy lifting

//	Open(theL1NtupleFileName);
	Loop(makeThresholdPlots,selectDataInput,lsFileName,L1MenuFileName,procNevts);


// Output File for saving rates
   TString outRatesFileName = "L1Rates_"; outRatesFileName +=  jobTag; outRatesFileName += ".txt";
   FILE *outRates = fopen(outRatesFileName,"w");

// Table header	
   printf(          "================================== L1Menu2015:  Rate Measures for L1 Menu ===================================\n");
	printf(          "L1Bit      L1SeedName   pre-sc      Thresholds               rate (kHz)    cumulative (kHz)     pure (kHz)  \n");
   fprintf(outRates,"================================= L1Menu2015:  Rate Measures for L1 Menu ==================================\n");
	fprintf(outRates,"L1Bit      L1SeedName   pre-sc      Thresholds               rate (kHz)    cumulative (kHz)     pure (kHz)  \n");
		
	Float_t totalrate     = 0.;
	Float_t totalrate_err = 0.;	
   Float_t finalL1Rate = h_All -> GetBinContent(kOFFSET+1);
	Float_t finalL1Rate_err = h_All -> GetBinError(kOFFSET+1);
	Float_t pureL1Rate = h_Pure -> GetBinContent(kOFFSET+1);
	Float_t pureL1Rate_err = h_Pure -> GetBinError(kOFFSET+1);

	for (Int_t k=1; k < kOFFSET+1; k++) {  // -- kOFFSET now contains the number of triggers we have calculated
		TString name = h_All -> GetXaxis() -> GetBinLabel(k);

		Float_t rate = h_All -> GetBinContent(k);
		Float_t err_rate  = h_All -> GetBinError(k);

		Float_t cumm_rate = h_Cumm -> GetBinContent(k);
      Float_t err_cumm_rate  = h_Cumm -> GetBinError(k);

		Float_t pure_rate = h_Pure -> GetBinContent(k);
		Float_t err_pure_rate = h_Pure -> GetBinError(k);

// Extract the Threshold values for this algorithm (Just Et for now)
		std::string L1namest = (std::string)name;
		Float_t pThres = trigParList[L1namest].primTh;
		Float_t sThres = trigParList[L1namest].secTh;
		Float_t tThres = trigParList[L1namest].triTh;
		Float_t qThres = trigParList[L1namest].quadTh;

// Pull out the thresholds being used
		std::map<std::string, int>::const_iterator it = Prescales.find(L1namest);  //blw a.
		Float_t pre;
		if (it == Prescales.end() ) {  //blw a.
			std::cout << " --- SET P = 1 FOR SEED :  " << L1namest << std::endl;
			pre = 1;
		}
		else {
			pre = it -> second;
		}

      totalrate += rate;
		totalrate_err += err_rate*err_rate;
		
//print the results  //blw a.
      printf("%2i %20s   %2i   %5.1f %5.1f %5.1f %5.1f  %8.2f  %5.2f   %8.2f  %5.2f   %8.2f  %5.2f \n",L1BitNumber(L1namest),name.Data(),(int)pre,pThres,sThres,tThres,qThres,rate,err_rate,cumm_rate,err_cumm_rate,pure_rate,err_pure_rate);				
      fprintf(outRates,"%2i %20s   %2i   %5.1f %5.1f %5.1f %5.1f  %8.2f  %5.2f   %8.2f  %5.2f   %8.2f  %5.2f \n",L1BitNumber(L1namest),name.Data(),(int)pre,pThres,sThres,tThres,qThres,rate,err_rate,cumm_rate,err_cumm_rate,pure_rate,err_pure_rate);				
		
	}
   printf("---------------------------------------------------------------------------------------------------------------\n");
   printf(  " Total L1 Rate (with overlaps)    = %8.2f +/- %5.2f  kHz\n",finalL1Rate, finalL1Rate_err); 
	printf(  " Total L1 Rate (without overlaps) = %8.2f +/- %5.2f  kHz\n",totalrate,sqrt(totalrate_err));
	printf(  " Total L1 Rate (pure triggers)    = %8.2f +/- %5.2f  kHz\n",pureL1Rate, pureL1Rate_err);	
	printf("===============================================================================================================\n");

   fprintf(outRates,"---------------------------------------------------------------------------------------------------------------\n");
   fprintf(outRates,  " Total L1 Rate (with overlaps)    = %8.2f +/- %5.2f  kHz\n",finalL1Rate, finalL1Rate_err); 
	fprintf(outRates,  " Total L1 Rate (without overlaps) = %8.2f +/- %5.2f  kHz\n",totalrate,sqrt(totalrate_err));
	fprintf(outRates,  " Total L1 Rate (pure triggers)    = %8.2f +/- %5.2f  kHz\n",pureL1Rate, pureL1Rate_err);	
	fprintf(outRates,"===============================================================================================================\n");

	
	fclose(outRates);	
	outHist->Write();
	outHist->Close();

	
	return finalL1Rate;
}
					 


