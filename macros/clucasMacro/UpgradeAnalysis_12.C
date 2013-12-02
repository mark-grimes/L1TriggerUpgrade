#include "L1Ntuple.h"
#include "hist.C"
#include "Style.C"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <utility>

// --------------------------------------------------------------------
//                UpgradeAnalysis_12 macro definition
// --------------------------------------------------------------------

struct order_gt : public std::binary_function<double, double, bool> {
  bool operator()(const double& x, const double& y) {
    return ( x > y ) ;
  }
};

//struct to order on the first entry (normally pT here) in a vector of pairs
struct order_gt_pairs : public std::binary_function<pair<float, float>, pair<float, float>, bool> {
  bool operator()(const pair<float, float>& x, const pair<float, float>& y) {
    return ( x.first > y.first ) ;
  }
};


class UpgradeAnalysis_12 : public L1Ntuple
{
  public :

    //constructor    
    UpgradeAnalysis_12(std::string filename) : L1Ntuple(filename) {}
    UpgradeAnalysis_12() {}
    ~UpgradeAnalysis_12() {}
    
    int FillDistros(Long64_t nevents, int lsStart, int lsFin, int PUScen, TString runFile, int doUpgradeObj);
    void DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis);
    void setChrisStyle();
    

  private : 
	
	void BookHistos();
	
	TH1D *met_hist;
	TH1D *met_rate;
	TH1D *ett_hist;
	TH1D *ett_rate;
	TH1D *mht_hist;
	TH1D *mht_rate;
	TH1D *htt_hist;
	TH1D *htt_rate;
	TH1D *tau_hist[4];
	TH1D *tau_rate[4];
	TH1D *isoEG_hist[4];
	TH1D *isoEG_rate[4];
	TH1D *muon_hist[4];
	TH1D *muon_rate[4];
	TH1D *muon_hist_hi[4];
	TH1D *muon_rate_hi[4];
	TH1D *combEG_hist[4];
	TH1D *combEG_rate[4];
	TH1D *combJetsEt_hist[6];
	TH1D *combJetsEt_rate[6];
	TH1D *combJetsEr_hist[6];
	TH1D *combJetsEr_rate[6];
    TH1D *cenpTauJets_distro[4];
    TH1D *cenpTauJets_rate[4];
    TH1D *cenJets_distro[4];
    TH1D *cenJets_rate[4];
    TH1D *fwdJets_distro[4];
    TH1D *fwdJets_rate[4]
;    TH1D *allJets_hist;

	TH1D *diJetEr_hist;
	TH1D *diJetEr_rate;
	TH1D *sinMuEr_hist;
	TH1D *sinMuEr_rate;
	
	TH1D *doubleEG_cross_hist;
	TH1D *doubleEG_cross_rate;
	TH1D *doubleMu_cross_hist;
	TH1D *doubleMu_cross_rate;
	TH1D *EGMu_cross_hist;
	TH1D *EGMu_cross_rate;
	TH1D *MuEG_cross_hist;
	TH1D *MuEG_cross_rate;
	TH1D *MuJet_cross_hist;
	TH1D *MuJet_cross_rate;
    TH1D *TauMu_cross_hist;
    TH1D *TauMu_cross_rate;
    TH1D *TauEG_cross_hist;
    TH1D *TauEG_cross_rate;
    TH1D *IsoEGCenJet_cross_hist;
    TH1D *IsoEGCenJet_cross_rate;
    TH1D *IsoEGMET_cross_hist;
    TH1D *IsoEGMET_cross_rate;
    TH1D *TauTwoFwd_cross_hist;
    TH1D* TauTwoFwd_cross_rate;

    TH1D *h_samplePU;

	TH1I *lumi_hist;
	TH1I *lumi_distro;
	TH1I *muonQual_hist;
	
	TGraph *lumiPU;
	TGraph *lumiIntL;
	TGraph *lumiInstL;
	
	TH1D *BXnoZero;
	TH1D *BXZero;

	TH2D *trigOverlap;
	TH1D *cumul_hist;
	TH1D *cumul_rate;

	TGraphErrors *puSinJ;
	TGraphErrors *puQuadJ;
	TGraphErrors *puSinEG;
	TGraphErrors *puTripEG;
	TGraphErrors *puETM;
	TGraphErrors *puETT;
	TGraphErrors *puHTM;
	TGraphErrors *puHTT;
	TGraphErrors *puSinMu;
	TGraphErrors *puTripMu;

    // Upgrade distros
    TH1D* up_isoEmEt_distro[4];
    TH1D* up_nonIsoEmEt_distro[4];
    TH1D* up_towerJetEt_distro[6];
    TH1D* up_isoTauEt_distro[4];
    TH1D* up_nonIsoTauEt_distro[4];
    TH1D* up_combEGEt_distro[4];
    TH1D* up_combTauEt_distro[4];
    TH1D* up_doubleEG_cross_hist;
    TH1D* up_EGMu_cross_hist;
    TH1D* up_MuEG_cross_hist;


    TH1D* up_isoEmEt_rate[4];
    TH1D* up_nonIsoEmEt_rate[4];
    TH1D* up_towerJetEt_rate[6];
    TH1D* up_isoTauEt_rate[4];
    TH1D* up_nonIsoTauEt_rate[4];
    TH1D* up_combEGEt_rate[4];
    TH1D* up_combTauEt_rate[4];
    TH1D* up_doubleEG_cross_rate;
    TH1D* up_EGMu_cross_rate;
    TH1D* up_MuEG_cross_rate;

};

void UpgradeAnalysis_12::setChrisStyle()
{
	//gStyle->SetOptStats("oeRMen");
}	


// --------------------------------------------------------------------
//                         bookhistos function 
// --------------------------------------------------------------------
void UpgradeAnalysis_12::BookHistos()
{

    met_hist            = 	new TH1D("met_distro", "met_distro", 100, 0., 750.);
    met_rate            = 	new TH1D("ETM_rate", "ETM", 100, 0., 750.);
    ett_hist            = 	new TH1D("ett_distro", "ett_distro", 100, 0., 750.);
    ett_rate            = 	new TH1D("ETT_rate", "ETT", 100, 0., 750.);
    mht_hist            = 	new TH1D("mht_distro", "mht_distro", 100, 0., 750.);
    mht_rate            = 	new TH1D("HTM_rate", "HTM", 100, 0., 750.);
    htt_hist            = 	new TH1D("ht_distro", "ht_distro", 100, 0., 750.);
    htt_rate            = 	new TH1D("HTT_rate", "HTT", 100, 0., 750.);
    
    tau_hist[0]         = 	new TH1D("tau_1_distro", "tau_1_distro_er2p17", 100, 0., 255.);
    tau_hist[1]         = 	new TH1D("tau_2_distro", "tau_2_distro_er2p17", 100, 0., 255.);
    tau_hist[2]         =   new TH1D("tau_3_distro", "tau_3_distro_er2p17", 100, 0., 255.);
    tau_hist[3]         =   new TH1D("tau_4_distro", "tau_4_distro_er2p17", 100, 0., 255.);    
    tau_rate[0]         =	new TH1D("tau_1_rate", "Single Tau #eta #leq 2.17", 100, 0., 255.);
    tau_rate[1]         =	new TH1D("tau_2_rate", "Double Tau #eta #leq 2.17", 100, 0., 255.);
    tau_rate[2]         =   new TH1D("tau_3_rate", "Triple Tau #eta #leq 2.17", 100, 0., 255.);
    tau_rate[3]         =   new TH1D("tau_4_rate", "Quad Tau #eta #leq 2.17", 100, 0., 255.);
    
    isoEG_hist[0]       =	new TH1D("isoEG_1_distro", "isoEG_1_distro_er2p17",100, 0., 63.);
    isoEG_hist[1]       =	new TH1D("isoEG_2_distro", "isoEG_2_distro_er2p17",100, 0., 63.);
    isoEG_hist[2]       =	new TH1D("isoEG_3_distro", "isoEG_3_distro_er2p17",100, 0., 63.);
    isoEG_hist[3]       =	new TH1D("isoEG_4_distro", "isoEG_4_distro_er2p17",100, 0., 63.);
    isoEG_rate[0]       =	new TH1D("isoEG_1_rate", "Single IsoEG #eta #leq 2.17",100, 0., 63.);
    isoEG_rate[1]       =	new TH1D("isoEG_2_rate", "Double IsoEG #eta #leq 2.17",100, 0., 63.);
    isoEG_rate[2]       =	new TH1D("isoEG_3_rate", "Triple IsoEG #eta #leq 2.17",100, 0., 63.);
    isoEG_rate[3]       =	new TH1D("isoEG_4_rate", "Quad IsoEG #eta #leq 2.17",100, 0., 63.);
    
    muon_hist[0]        =	new TH1D("muon_1_distro", "muon_1_distro", 100, 0., 140.);
    muon_hist[1]        =	new TH1D("muon_2_distro", "muon_2_distro", 100, 0., 140.);
    muon_hist[2]        =	new TH1D("muon_3_distro", "muon_3_distro", 100, 0., 140.);
    muon_hist[3]        =	new TH1D("muon_4_distro", "muon_4_distro", 100, 0., 140.);
    muon_rate[0]        =	new TH1D("muon_1_rate", "Single Muon", 100, 0., 140.);
    muon_rate[1]        =	new TH1D("muon_2_rate", "Double Muon", 100, 0., 140.);
    muon_rate[2]        =	new TH1D("muon_3_rate", "Triple Muon", 100, 0., 140.);
    muon_rate[3]        =	new TH1D("muon_4_rate", "Quad Muon", 100, 0., 140.);
    
    muon_hist_hi[0]     =	new TH1D("muon_hi_1_distro", "muon_1_distro_hi", 100, 0., 140.);
    muon_hist_hi[1]     =	new TH1D("muon_hi_2_distro", "muon_2_distro_hi", 100, 0., 140.);
    muon_hist_hi[2]     =	new TH1D("muon_hi_3_distro", "muon_3_distro_hi", 100, 0., 140.);
    muon_hist_hi[3]     =	new TH1D("muon_hi_4_distro", "muon_4_distro_hi", 100, 0., 140.);
    muon_rate_hi[0]     =	new TH1D("muon_hi_1_rate", "Single Muon hi", 100, 0., 140.);
    muon_rate_hi[1]     =	new TH1D("muon_hi_2_rate", "Double Muon hi", 100, 0., 140.);
    muon_rate_hi[2]     =	new TH1D("muon_hi_3_rate", "Triple Muon hi", 100, 0., 140.);
    muon_rate_hi[3]     =	new TH1D("muon_hi_4_rate", "Quad Muon hi", 100, 0., 140.);
    
    
    combEG_hist[0]      =	new TH1D("combEG_1_distro", "combEG_1_distro",100, 0., 63.);
    combEG_hist[1]      =	new TH1D("combEG_2_distro", "combEG_2_distro",100, 0., 63.);
    combEG_hist[2]      =	new TH1D("combEG_3_distro", "combEG_3_distro",100, 0., 63.);
    combEG_hist[3]      =	new TH1D("combEG_4_distro", "combEG_4_distro",100, 0., 63.);
    combEG_rate[0]      =	new TH1D("combEG_1_rate", "Single EG",100, 0., 63.);
    combEG_rate[1]      =	new TH1D("combEG_2_rate", "Double EG",100, 0., 63.);
    combEG_rate[2]      =	new TH1D("combEG_3_rate", "Triple EG",100, 0., 63.);
    combEG_rate[3]      =	new TH1D("combEG_4_rate", "Quad EG",100, 0., 63.);
    
    combJetsEt_hist[0]  =	new TH1D("combJets_1_distro", "combJets_1_distro", 100, 0., 255.);
    combJetsEt_hist[1]  =	new TH1D("combJets_2_distro", "combJets_2_distro", 100, 0., 255.);
    combJetsEt_hist[2]  =	new TH1D("combJets_3_distro", "combJets_3_distro", 100, 0., 255.);
    combJetsEt_hist[3]  =	new TH1D("combJets_4_distro", "combJets_4_distro", 100, 0., 255.);
    combJetsEt_hist[4]  =   new TH1D("combJets_5_distro", "combJets_5_distro", 100, 0., 255.);
    combJetsEt_hist[5]  =   new TH1D("combJets_6_distro", "combJets_6_distro", 100, 0., 255.);    
    combJetsEt_rate[0]  =	new TH1D("combJets_1_rate", "Single Jet", 100, 0., 255.);
    combJetsEt_rate[1]  =	new TH1D("combJets_2_rate", "Double Jet", 100, 0., 255.);
    combJetsEt_rate[2]  =	new TH1D("combJets_3_rate", "Triple Jet", 100, 0., 255.);
    combJetsEt_rate[3]  =	new TH1D("combJets_4_rate", "Quad Jet", 100, 0., 255.);
    combJetsEt_rate[4]  =   new TH1D("combJets_5_rate", "Five Jet", 100, 0., 255.);
    combJetsEt_rate[5]  =   new TH1D("combJets_6_rate", "Six Jet", 100, 0., 255.);    

    combJetsEr_hist[0]  =	new TH1D("combJetsEr_1_distro", "combJetsEr_1_distro", 100, 0., 255.);
    combJetsEr_hist[1]  =	new TH1D("combJetsEr_2_distro", "combJetsEr_2_distro", 100, 0., 255.);
    combJetsEr_hist[2]  =	new TH1D("combJetsEr_3_distro", "combJetsEr_3_distro", 100, 0., 255.);
    combJetsEr_hist[3]  =	new TH1D("combJetsEr_4_distro", "combJetsEr_4_distro", 100, 0., 255.);
    combJetsEr_hist[4]  =   new TH1D("combJetsEr_5_distro", "combJetsEr_5_distro", 100, 0., 255.);
    combJetsEr_hist[5]  =   new TH1D("combJetsEr_6_distro", "combJetsEr_6_distro", 100, 0., 255.);    
    combJetsEr_rate[0]  =	new TH1D("combJetsEr_1_rate", "Single Jet #eta #leq 3", 100, 0., 255.);
    combJetsEr_rate[1]  =	new TH1D("combJetsEr_2_rate", "Double Jet #eta #leq 3", 100, 0., 255.);
    combJetsEr_rate[2]  =	new TH1D("combJetsEr_3_rate", "Triple Jet #eta #leq 3", 100, 0., 255.);
    combJetsEr_rate[3]  =	new TH1D("combJetsEr_4_rate", "Quad Jet #eta #leq 3", 100, 0., 255.);
    combJetsEr_rate[4]  =   new TH1D("combJetsEr_5_rate", "Five Jet #eta #leq 3", 100, 0., 255.);
    combJetsEr_rate[5]  =   new TH1D("combJetsEr_6_rate", "Six Jet #eta #leq 3", 100, 0., 255.);

    fwdJets_distro[0]  =   new TH1D("fwdJets_1_distro", "fwdJets_1_distro", 100, 0., 255.);
    fwdJets_distro[1]  =   new TH1D("fwdJets_2_distro", "fwdJets_2_distro", 100, 0., 255.);
    fwdJets_distro[2]  =   new TH1D("fwdJets_3_distro", "fwdJets_3_distro", 100, 0., 255.);
    fwdJets_distro[3]  =   new TH1D("fwdJets_4_distro", "fwdJets_4_distro", 100, 0., 255.);   
    fwdJets_rate[0]  =   new TH1D("fwdJets_1_rate", "fwdJets_1_rate", 100, 0., 255.);
    fwdJets_rate[1]  =   new TH1D("fwdJets_2_rate", "fwdJets_2_rate", 100, 0., 255.);
    fwdJets_rate[2]  =   new TH1D("fwdJets_3_rate", "fwdJets_3_rate", 100, 0., 255.);
    fwdJets_rate[3]  =   new TH1D("fwdJets_4_rate", "fwdJets_4_rate", 100, 0., 255.); 

    cenJets_distro[0]  =   new TH1D("cenJets_1_distro", "cenJets_1_distro", 100, 0., 255.);
    cenJets_distro[1]  =   new TH1D("cenJets_2_distro", "cenJets_2_distro", 100, 0., 255.);
    cenJets_distro[2]  =   new TH1D("cenJets_3_distro", "cenJets_3_distro", 100, 0., 255.);
    cenJets_distro[3]  =   new TH1D("cenJets_4_distro", "cenJets_4_distro", 100, 0., 255.);         
    cenJets_rate[0]  =   new TH1D("cenJets_1_rate", "cenJets_1_rate", 100, 0., 255.);
    cenJets_rate[1]  =   new TH1D("cenJets_2_rate", "cenJets_2_rate", 100, 0., 255.);
    cenJets_rate[2]  =   new TH1D("cenJets_3_rate", "cenJets_3_rate", 100, 0., 255.);
    cenJets_rate[3]  =   new TH1D("cenJets_4_rate", "cenJets_4_rate", 100, 0., 255.);   


    diJetEr_hist        =	new TH1D("DoubleJetEtar_distro", "Double Jet Distro #eta #leq 3", 100, 0., 255.);
    diJetEr_rate        =	new TH1D("DoubleJetEtar_rate", "Double Jet Rate #eta #leq 3", 100, 0., 255.);
    sinMuEr_hist        =	new TH1D("SingleMuonEtar_distro", "Single Muon Distro #eta #leq 2.1", 100, 0., 140.);
    sinMuEr_rate        =	new TH1D("SingleMuonEtar_rate", "Single Muon Rate #eta #leq 2.1", 100, 0., 140.);
    
    cenpTauJets_distro[0]   =   new TH1D("cenpTauJets_1_distro", "Single Central and Tau Jets #eta #leq 3", 100, 0, 255.);
    cenpTauJets_distro[1]   =   new TH1D("cenpTauJets_2_distro", "Double Central and Tau Jets #eta #leq 3", 100, 0, 255.);
    cenpTauJets_distro[2]   =   new TH1D("cenpTauJets_3_distro", "Triple Central and Tau Jets #eta #leq 3", 100, 0, 255.);
    cenpTauJets_distro[3]   =   new TH1D("cenpTauJets_4_distro", "Quad Central and Tau Jets #eta #leq 3", 100, 0, 255.);

    cenpTauJets_rate[0]   =   new TH1D("cenpTauJets_1_rate", "Single Central and Tau Jets Rate #eta #leq 3", 100, 0, 255.);
    cenpTauJets_rate[1]   =   new TH1D("cenpTauJets_2_rate", "Double Central and Tau Jets Rate #eta #leq 3", 100, 0, 255.);
    cenpTauJets_rate[2]   =   new TH1D("cenpTauJets_3_rate", "Triple Central and Tau Jets Rate #eta #leq 3", 100, 0, 255.);
    cenpTauJets_rate[3]   =   new TH1D("cenpTauJets_4_rate", "Quad Central and Tau Jets Rate #eta #leq 3", 100, 0, 255.);


    doubleEG_cross_hist =	new TH1D("doubleEGCross_distro", "Double EG Cross Trigger distro", 100, 0., 63.);
    doubleEG_cross_rate =	new TH1D("doubleEGCross_rate", "Double EG Cross Trigger", 100, 0., 63.);
    doubleMu_cross_hist =	new TH1D("doubleMuCross_distro", "Double Mu Cross Trigger distro", 100, 0., 140.);
    doubleMu_cross_rate =	new TH1D("doubleMuCross_rate", "Double Mu Cross Trigger", 100, 0., 140.);
    EGMu_cross_hist     =	new TH1D("EGMuCross_distro", "EG Mu Cross Trigger distro", 100, 0., 63.);
    EGMu_cross_rate     =	new TH1D("EGMuCross_rate", "EG Mu Cross Trigger", 100, 0., 63.);
    MuEG_cross_hist     =	new TH1D("MuEGCross_distro", "Mu EG Cross Trigger distro", 100, 0., 140.);
    MuEG_cross_rate     =	new TH1D("MuEGCross_rate", "Mu EG Cross Trigger", 100, 0., 140.);
    MuJet_cross_hist    =	new TH1D("MuJetCross_distro", "Mu Jet Cross Trigger distro", 100, 0., 140.);
    MuJet_cross_rate    =	new TH1D("MuJetCross_rate", "Mu Jet Cross Trigger", 100, 0., 140.);
    TauEG_cross_hist    =   new TH1D("TauEGCross_distro", "Tau EG Cross Trigger distro", 100, 0., 140.);
    TauEG_cross_rate    =   new TH1D("TauEGCross_rate", "Tau EG Cross Trigger", 100, 0., 140.);
    TauMu_cross_hist    =   new TH1D("TauMuCross_distro", "Tau Mu Cross Trigger distro", 100, 0., 140.);
    TauMu_cross_rate    =   new TH1D("TauMuCross_rate", "Tau Mu Cross Trigger", 100, 0., 140.);
    IsoEGCenJet_cross_hist = new TH1D("IsoEGCenJet_distro", "IsoEG CenJet Cross Trigger distro", 100, 0., 100.);
    IsoEGCenJet_cross_rate = new TH1D("IsoEGCenJet_rate", "IsoEG CenJet Cross Trigger", 100, 0., 100.);
    IsoEGMET_cross_hist = new TH1D("IsoEGMET_distro", "IsoEG MET Cross Trigger distro", 100, 0., 100.);
    IsoEGMET_cross_rate = new TH1D("IsoEGMET_rate", "IsoEG MET Cross Trigger", 100, 0., 100.);
    TauTwoFwd_cross_hist    =   new TH1D("TauTwoFwdCross_distro", "Tau TwoFwd Cross Trigger distro", 100, 0., 140.);
    TauTwoFwd_cross_rate    =   new TH1D("TauTwoFwdCross_rate", "Tau TwoFwd Cross Trigger", 100, 0., 140.);    

    h_samplePU          =    new TH1D("h_samplePU", "Sample PU", 5000, 0., 100.);
    trigOverlap         =	new TH2D("trigOverlap", "Trigger Overlap", 14, 0., 14., 14, 0., 14.);
    cumul_hist          =	new TH1D("cumulative_distro", "cumulative distro", 13, 0., 13.);
    cumul_rate          =	new TH1D("cumulative_rate", "Cumulative Rate", 13 , 0., 13.);
    
    allJets_hist        =	new TH1D("All Jets", "All Jets", 255, 0., 255.);
    
    lumi_hist           =	new TH1I("lumi", "lumis", 500, 0., 500.);
    lumi_distro         =	new TH1I("lumi_event_distro", "lumi_event_distro", 500, 0., 500.);
    muonQual_hist       =	new TH1I("muonQual_hist", "muonQual_hist", 20, 0., 20.);
    BXnoZero            =	new TH1D("BX before Zero Trigger", "BX before Zero Trigger", 4000, 0., 4000.);
    BXZero              =	new TH1D("BX after Zero Trigger", "BX after Zero Trigger", 4000, 0., 4000.);
   
   // upgrade histos
    up_towerJetEt_distro[0]  =    new TH1D("up_towerJet_1_distro", "Distro of Single Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_distro[1]  =    new TH1D("up_towerJet_2_distro", "Distro of Double Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_distro[2]  =    new TH1D("up_towerJet_3_distro", "Distro of Triple Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_distro[3]  =    new TH1D("up_towerJet_4_distro", "Distro of Quad Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_distro[4]  =    new TH1D("up_towerJet_5_distro", "Distro of Five Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_distro[5]  =    new TH1D("up_towerJet_6_distro", "Distro of Six Upgrade Tower Jets", 100, 0., 255.);    
    
    up_towerJetEt_rate[0]    =    new TH1D("up_towerJet_1_rate", "Single Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_rate[1]    =    new TH1D("up_towerJet_2_rate", "Double Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_rate[2]    =    new TH1D("up_towerJet_3_rate", "Triple Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_rate[3]    =    new TH1D("up_towerJet_4_rate", "Quad Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_rate[4]    =    new TH1D("up_towerJet_5_rate", "Five Upgrade Tower Jets", 100, 0., 255.);
    up_towerJetEt_rate[5]    =    new TH1D("up_towerJet_6_rate", "Six Upgrade Tower Jets", 100, 0., 255.);    

    up_nonIsoEmEt_distro[0]  =   new TH1D("up_nonIsoEm_1_distro", "Distro of Single Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_distro[1]  =   new TH1D("up_nonIsoEm_2_distro", "Distro of Double Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_distro[2]  =   new TH1D("up_nonIsoEm_3_distro", "Distro of Triple Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_distro[3]  =   new TH1D("up_nonIsoEm_4_distro", "Distro of Quad Upgrade nonIsoEG", 100, 0., 63.);
    
    up_nonIsoEmEt_rate[0]    =   new TH1D("up_nonIsoEm_1_rate", "Single Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_rate[1]    =   new TH1D("up_nonIsoEm_2_rate", "Double Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_rate[2]    =   new TH1D("up_nonIsoEm_3_rate", "Triple Upgrade nonIsoEG", 100, 0., 63.);
    up_nonIsoEmEt_rate[3]    =   new TH1D("up_nonIsoEm_4_rate", "Quad Upgrade nonIsoEG", 100, 0., 63.);
    
    up_isoEmEt_distro[0]     =   new TH1D("up_isoEm_1_distro", "Distro of Single Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_distro[1]     =   new TH1D("up_isoEm_2_distro", "Distro of Double Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_distro[2]     =   new TH1D("up_isoEm_3_distro", "Distro of Triple Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_distro[3]     =   new TH1D("up_isoEm_4_distro", "Distro of Quad Upgrade isoEG", 100, 0., 63.);
    
    up_isoEmEt_rate[0]       =   new TH1D("up_isoEm_1_rate", "Single Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_rate[1]       =   new TH1D("up_isoEm_2_rate", "Double Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_rate[2]       =   new TH1D("up_isoEm_3_rate", "Triple Upgrade isoEG", 100, 0., 63.);
    up_isoEmEt_rate[3]       =   new TH1D("up_isoEm_4_rate", "Quad Upgrade isoEG", 100, 0., 63.);
    
    up_isoTauEt_distro[0]    =   new TH1D("up_isoTau_1_distro", "Distro of Single Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_distro[1]    =   new TH1D("up_isoTau_2_distro", "Distro of Double Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_distro[2]    =   new TH1D("up_isoTau_3_distro", "Distro of Triple Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_distro[3]    =   new TH1D("up_isoTau_4_distro", "Distro of Quad Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    
    up_isoTauEt_rate[0]      =   new TH1D("up_isoTau_1_rate", "Single Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_rate[1]      =   new TH1D("up_isoTau_2_rate", "Double Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_rate[2]      =   new TH1D("up_isoTau_3_rate", "Triple Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    up_isoTauEt_rate[3]      =   new TH1D("up_isoTau_4_rate", "Quad Upgrade isoTau #eta #leq 2.17", 100, 0., 255.);
    
    up_nonIsoTauEt_distro[0] =   new TH1D("up_nonIsoTau_1_distro", "Distro of Single Upgrade nonIsoTau #eta #leq 2.17 #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_distro[1] =   new TH1D("up_nonIsoTau_2_distro", "Distro of Double Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_distro[2] =   new TH1D("up_nonIsoTau_3_distro", "Distro of Triple Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_distro[3] =   new TH1D("up_nonIsoTau_4_distro", "Distro of Quad Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    
    up_nonIsoTauEt_rate[0]   =   new TH1D("up_nonIsoTau_1_rate", "Single Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_rate[1]   =   new TH1D("up_nonIsoTau_2_rate", "Double Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_rate[2]   =   new TH1D("up_nonIsoTau_3_rate", "Triple Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);
    up_nonIsoTauEt_rate[3]   =   new TH1D("up_nonIsoTau_4_rate", "Quad Upgrade nonIsoTau #eta #leq 2.17", 100, 0., 255.);

    up_combTauEt_distro[0]     =   new TH1D("up_combTau_1_distro", "Distro of Single Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_distro[1]     =   new TH1D("up_combTau_2_distro", "Distro of Double Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_distro[2]     =   new TH1D("up_combTau_3_distro", "Distro of Triple Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_distro[3]     =   new TH1D("up_combTau_4_distro", "Distro of Quad Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    
    up_combTauEt_rate[0]     =   new TH1D("up_combTau_1_rate", "Single Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_rate[1]     =   new TH1D("up_combTau_2_rate", "Double Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_rate[2]     =   new TH1D("up_combTau_3_rate", "Triple Upgrade combTau #eta #leq 2.17", 100, 0., 255.);
    up_combTauEt_rate[3]     =   new TH1D("up_combTau_4_rate", "Quad Upgrade combTau #eta #leq 2.17", 100, 0., 255.);

    up_combEGEt_distro[0]     =   new TH1D("up_combEG_1_distro", "Distro of Single Upgrade combEG", 100, 0., 63.);
    up_combEGEt_distro[1]     =   new TH1D("up_combEG_2_distro", "Distro of Double Upgrade combEG", 100, 0., 63.);
    up_combEGEt_distro[2]     =   new TH1D("up_combEG_3_distro", "Distro of Triple Upgrade combEG", 100, 0., 63.);
    up_combEGEt_distro[3]     =   new TH1D("up_combEG_4_distro", "Distro of Quad Upgrade combEG", 100, 0., 63.);
    
    up_combEGEt_rate[0]     =   new TH1D("up_combEG_1_rate", "Single Upgrade combEG", 100, 0., 63.);
    up_combEGEt_rate[1]     =   new TH1D("up_combEG_2_rate", "Double Upgrade combEG", 100, 0., 63.);
    up_combEGEt_rate[2]     =   new TH1D("up_combEG_3_rate", "Triple Upgrade combEG", 100, 0., 63.);
    up_combEGEt_rate[3]     =   new TH1D("up_combEG_4_rate", "Quad Upgrade combEG", 100, 0., 63.);


    up_doubleEG_cross_hist =   new TH1D("up_doubleEGCross_distro", "Double EG Upgrade Cross Trigger distro", 100, 0., 63.);
    up_doubleEG_cross_rate =   new TH1D("up_doubleEGCross_rate", "Double EG Upgrade Cross Trigger", 100, 0., 63.);
    up_EGMu_cross_hist     =   new TH1D("up_EGMuCross_distro", "EG Mu Upgrade Cross Trigger distro", 100, 0., 63.);
    up_EGMu_cross_rate     =   new TH1D("up_EGMuCross_rate", "EG Mu Upgrade Cross Trigger", 100, 0., 63.);
    up_MuEG_cross_hist     =   new TH1D("up_MuEGCross_distro", "Mu EG Upgrade Cross Trigger distro", 100, 0., 140.);
    up_MuEG_cross_rate     =   new TH1D("up_MuEGCross_rate", "Mu EG Upgrade Cross Trigger", 100, 0., 140.);

}

// --------------------------------------------------------------------
//                          ratecalc function
// --------------------------------------------------------------------
void UpgradeAnalysis_12::DoRateCalc(TH1D* h1, TH1D *h2, int preScale, int nLumis) //TH1D *enScale)
{
	//short function to calculate rate plots

    Int_t nbins  = h1->GetNbinsX();

    Double_t tLumi = 23.3570304;
    Double_t targetLumi = 2.e+34; //target lumi after scaling
    Double_t nSamples = 1; //hack for more than 1 ZB sample

    Double_t normalization = (1.*targetLumi)/(nSamples*nLumis*tLumi);

    h2->GetXaxis()->SetTitle("Threshold (GeV)");
    h2->GetYaxis()->SetTitle("Rate (Hz)");

    for(Int_t i=1; i < nbins; ++i)
    {
    	double histError = 0;
        float rate = h1->IntegralAndError(i, nbins+1, histError);
        h2->SetBinContent(i, rate);
        h2->SetBinError(i, histError);
    }
    
    h2->Scale(normalization); //scale entire plot according to normalization factor


}

// --------------------------------------------------------------------
//                             run function 
// --------------------------------------------------------------------
int UpgradeAnalysis_12::FillDistros(Long64_t nevents, int lsStart, int lsFin, int PUScen, TString runFile, int doUpgradeObj)
{
	//function for making generic rate vs threshold plots

	std::cout << "\n\n**** Filling pT Distributions\n" << std::endl;


	//load TDR style
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetFillStyle(0);
	gStyle->SetFillColor(0);
	
	std::vector<int>	lumis;
	std::vector<string> instLumis;

	int counter=0, ind=0;
	int iMHT, iMET;
    double sampAvgPU=0;

	//--- Get Run details ---//

	double LS[450], IntL[450], PU[450], InstL[450];

	ifstream ifs( runFile );
	while(ifs){
		ifs >> LS[ind];
		ifs >> IntL[ind];
		ifs >> InstL[ind];
		ifs >> PU[ind];
		ind++;
	}

	//-------------------------------

	
	std::ostringstream outFileName;
	outFileName << "out_files/output_" << lsStart << "-" << lsFin << ".root";
	
	std::cout << "Creating: " << outFileName.str().c_str() << std::endl;
	TFile *outFile = TFile::Open(outFileName.str().c_str(), "RECREATE");
	outFile->cd();
	
    TH1D::SetDefaultSumw2();
	BookHistos();
	
	int ntot = 0;
    double preScale=1.;
	bool anyGoodLumi = false;
	
	
	if (nevents==-1 || nevents>GetEntries()) nevents=GetEntries();
	std::cout << nevents << " to process ..." << std::endl;
	
	//loop over the events to fill all object vectors
	for (Long64_t i=0; i<nevents; i++)
	{
		//load the i-th event 
		Long64_t ientry = LoadTree(i); if (ientry < 0) break;
		GetEntry(i);
		BXnoZero->Fill(event_->bx);

		//if (event_->run == 198588) continue; 

		//define trigger booleans
		int trigNum=14;
		bool trigRes[trigNum];

		for(int i=0; i<trigNum; i++) trigRes[i] = false; //set all to false by default

		if(i!=0 && (i%1000)==0) {std::cout << "- processing event " << i << "\r" << std::flush;}
		if ((gt_->tw1[2] & 0x0001) > 0){ // check that zero bias trigger fired
			
			BXZero->Fill(event_->bx);

			double lumiWeight=1.;
			int thisLumi = event_->lumi; //get this event's lumi
		
			// set run-specific preScales
			if (event_->run == 198588) preScale = 44.;
			if (event_->run == 198603) preScale = 92.;
			if (event_->run == 198609) preScale = 92.;
	
			for(int k=0; k<ind-1; k++){
				if( LS[k]==thisLumi ) lumiWeight = preScale/InstL[k];
            }
			//lumiWeight = 1.;

			//this is used to run over specific lumi bins
			if ((thisLumi>=lsStart) && (thisLumi<=lsFin)){
				anyGoodLumi = true;
				
				// create vectors for all jets and em's
				std::vector<double>	combEm;
                std::vector<double> combEmUp;
				std::vector< pair<double, double> > combTauUp;
                std::vector< pair<float,float> > muon_hi, muon_lo, muon_open;
				std::vector< pair<double, double> > combJets;
                std::vector< pair<double, double> > combCenTauJets;

				
				ntot++; //counter for events ran over

				//check if the current lumi is unique (used later for rate calc)
				bool newLumi = true;
				lumi_distro->Fill(thisLumi);
				

				//THIS MUST BE CHANGED IF MULTIPLE RUNS ARE USED
				for(unsigned int i=0; i<lumis.size(); ++i){
					if(lumis.at(i)==thisLumi) newLumi = false;
				}
				if(newLumi){
                    lumis.push_back(thisLumi); //enter unique lumis to vector
                    for(int k=0; k<ind-1; k++){
                        if( LS[k]==thisLumi ) sampAvgPU += PU[k];
                    }
                }
				   //------- Start Filling Objects --------//

				for(unsigned int i=0; i<l1extra_->metBx.size(); i++){
					if (l1extra_->metBx.at(i) == 0) iMET=i;
				}
				for(unsigned int i=0; i<l1extra_->mhtBx.size(); i++){
					if (l1extra_->mhtBx.at(i) == 0) iMHT=i;
				}
				
				met_hist->	Fill(l1extra_->met.at(iMET), lumiWeight);
				mht_hist->	Fill(l1extra_->mht.at(iMHT), lumiWeight);
				ett_hist->	Fill(l1extra_->et.at(iMET), lumiWeight);
				htt_hist->	Fill(l1extra_->ht.at(iMHT), lumiWeight);

				//------------------muon--------------------
				for(unsigned int i=0; i<l1extra_->nMuons; i++){
					if(l1extra_->muonEt.at(i) > 0.){
						if(i<4){
							muon_hist[i]->Fill(l1extra_->muonEt.at(i), lumiWeight);
						}
					}
				}
				
				//loop for filling hi-quality muons, with BX condition to remove cosmics
				for(int i=0; i<gmt_->N; i++){
					if(gmt_->Pt.at(i) > 0.){
                        muonQual_hist->Fill(gmt_->Qual.at(i));
						if(i<4){
							if(gmt_->CandBx.at(i)==0){ //check is the same BX as trigger
                                // quality requirements vary for different muon triggers
                                if( fabs(gmt_->Eta.at(i) <= 2.1) ){
    								if( (i==0) && (gmt_->Qual.at(i) >= 4.) ){
                                        pair<float, float> muonHiCand(gmt_->Pt.at(i), gmt_->Eta.at(i));
    									muon_hi.push_back(muonHiCand);
    								}
                                    if(gmt_->Qual.at(i) >= 3.){
                                        pair<float, float> muonLoCand(gmt_->Pt.at(i), gmt_->Eta.at(i));
                                        muon_lo.push_back(muonLoCand);
                                    }
                                    if(gmt_->Qual.at(i) >= 2.){
                                        pair<float, float> muonOpenCand(gmt_->Pt.at(i), gmt_->Eta.at(i));
                                        muon_open.push_back(muonOpenCand);                                    
                                    }
                                }
							}
						}
					}
				}

				sort(muon_hi.begin(), muon_hi.end(), order_gt_pairs());
                sort(muon_lo.begin(), muon_lo.end(), order_gt_pairs());
                sort(muon_open.begin(), muon_open.end(), order_gt_pairs());

                for(unsigned int i=0; i<muon_hi.size(); i++){
                    if (i<4) muon_hist_hi[i]->Fill(muon_hi.at(i).first, lumiWeight);
                }


				//if (muon_hi.size()>0){
                //    muon_hist_hi[0]->Fill(muon_hi.at(0).first, lumiWeight);
                //}
//
                //for(unsigned int i=1; i<muon_lo.size(); i++){
                //    if (i<4) muon_hist_hi[i]->Fill(muon_lo.at(i).first, lumiWeight);
                //}
//
                //for(unsigned int i=0; i<muon_hi.size(); i++){
				//	if(i<4){
				//		//muon_hist_hi[i]->Fill(muon_hi.at(i).first, lumiWeight);
				//		if((i==0) && (fabs(muon_hi.at(i).second) <= 2.1)) sinMuEr_hist->Fill(muon_hi.at(i).first, lumiWeight);
				//	}
				//}

				//---------------combined EG----------------
				// isoEG
				for(unsigned int i=0; i<l1extra_->nIsoEm; i++){
					if(l1extra_->isoEmEt.at(i) > 0.){
						if((l1extra_->isoEmBx.at(i))!=0) break;
						double et = l1extra_->isoEmEt.at(i);
						combEm.push_back(et);
						if( fabs(l1extra_->isoEmEta.at(i)) <= 2.17){
							if (i<4) isoEG_hist[i]->Fill(et, lumiWeight);
						}
					}
				}
				
				// nonisoEG
				for(unsigned int i=0; i<l1extra_->nNonIsoEm; i++){
					if(l1extra_->nonIsoEmEt.at(i) > 0.){
						if((l1extra_->nonIsoEmBx.at(i))!=0) break;
						double et = l1extra_->nonIsoEmEt.at(i); 
						combEm.push_back(et);
					}
				}
				
				sort(combEm.begin(), combEm.end(), order_gt()); //sort them into descending order
				
				// fill distros for combined EG objects
				for(unsigned int i=0; i<4; i++){
					if (combEm.size()>i) combEG_hist[i]->Fill(combEm.at(i), lumiWeight);
				}
				
				//---------------combined Jet----------------
				
				//tau-jets
				for(unsigned int i=0; i<l1extra_->nTauJets; i++){ 
					if(l1extra_->tauJetEt.at(i) > 0.){
						if((l1extra_->tauJetBx.at(i))!=0) break;
						double et = l1extra_->tauJetEt.at(i);
						double eta = l1extra_->tauJetEta.at(i);
						pair<double, double> candJet(et, eta);
						combJets.push_back(candJet);

                        if (fabs(eta) < 3.) combCenTauJets.push_back(candJet);

						if (fabs(eta)<=2.17 ){ //apply the tau eta restriction
							if (i<4) tau_hist[i]->Fill(et, lumiWeight);
						}

					}
				}
				
				//fwd-jets
				for(unsigned int i=0; i<l1extra_->nFwdJets; i++){
					if(l1extra_->fwdJetEt.at(i) > 0.){
						if((l1extra_->fwdJetBx.at(i))!=0) break;
						double et = l1extra_->fwdJetEt.at(i);
						double eta = l1extra_->fwdJetEta.at(i);
						pair<double, double> candJet(et, eta);

						combJets.push_back(candJet);

                        if (i<4) fwdJets_distro[i]->Fill(et, lumiWeight);

					}
				}
				
				//cen-jets
				for(unsigned int i=0; i<l1extra_->nCenJets; i++){
					if(l1extra_->cenJetEt.at(i) > 0.){
						if((l1extra_->cenJetBx.at(i))!=0) break;
						double et = l1extra_->cenJetEt.at(i);
						double eta = l1extra_->cenJetEta.at(i);

						pair<double, double> candJet(et,eta);

						combJets.push_back(candJet);
                        if (fabs(eta) < 3.) combCenTauJets.push_back(candJet);
					
                        if (i<4) cenJets_distro[i]->Fill(et, lumiWeight);
                    }
				}
				
				sort(combJets.begin(), combJets.end(), order_gt_pairs());
                sort(combCenTauJets.begin(), combCenTauJets.end(), order_gt_pairs());
				
				for(unsigned int i=0; i<combJets.size(); i++) allJets_hist->Fill(combJets.at(i).first, lumiWeight);
				
                for(unsigned int i=0; i<combCenTauJets.size(); i++) cenpTauJets_distro[i]->Fill(combCenTauJets.at(i).first, lumiWeight);

				for(unsigned int i=0; i<6; i++){
					if (combJets.size() > i){
						combJetsEt_hist[i]->Fill(combJets.at(i).first, lumiWeight);
						if ( (i==1) && (fabs(combJets.at(i).second) <= 3.) ) diJetEr_hist->Fill(combJets.at(i).first, lumiWeight);
						if ( fabs(combJets.at(i).second) <= 3. )  combJetsEr_hist[i]->Fill(combJets.at(i).first, lumiWeight );
					}
				}

				//---------------cross triggers---------------
				//doubleEG 13:7
				if (combEm.size() > 1){
				    if (combEm.at(1) > floor(0.5+7.*combEm.at(0)/13.) ) doubleEG_cross_hist->Fill(combEm.at(0), lumiWeight);
				}

				//doubleMu 12:3.5
				if (muon_hi.size() > 1){
				    if (muon_hi.at(1).first > floor(0.5+3.5*muon_hi.at(0).first/12.)) doubleMu_cross_hist->Fill(muon_hi.at(0).first, lumiWeight);
				}

				//EG+mu 12:3.5  and   mu+EG 12:7
				if ( (muon_hi.size() > 0) && (combEm.size() > 0) ){
				    if (muon_hi.at(0).first > floor(0.5 + (3.5*combEm.at(0))/12.)) EGMu_cross_hist->Fill(combEm.at(0), lumiWeight);
				    if (combEm.at(0) > floor(0.5 + (7.*muon_hi.at(0).first)/12.))  MuEG_cross_hist->Fill(muon_hi.at(0).first, lumiWeight);
				}

                //Tau+Mu 12:3.5
                if ( (muon_hi.size() > 0) && (l1extra_->tauJetEt.size() > 0) ){
                    if (l1extra_->tauJetBx.at(0) != 0) break;
                    if ( muon_hi.at(0).first > floor(0.5+3.5*l1extra_->tauJetEt.at(0)/12.) ){
                        TauMu_cross_hist->Fill(l1extra_->tauJetEt.at(0), lumiWeight);
                    }
                }

                //Tau+EG 12:7
                if ( (combEm.size()>0) && (l1extra_->tauJetEt.size() > 0) ){
                    if (l1extra_->tauJetBx.at(0) != 0) break;
                    if ( combEm.at(0) > floor(0.5+3.5*l1extra_->tauJetEt.at(0)/12.) ){
                        TauEG_cross_hist->Fill(l1extra_->tauJetEt.at(0), lumiWeight);
                    } 
                }

                //IsoEGCenJet 13:7
                for(unsigned int i=0; i<l1extra_->nIsoEm;i++){
                    if (l1extra_->nCenJets <= i) break;
                    if ( (fabs(l1extra_->isoEmEta.at(i)) <=2.17) && (l1extra_->isoEmBx.at(i)==0) && (l1extra_->cenJetBx.at(i)==0)){
                        if ( l1extra_->cenJetEt.at(i) > floor(0.5+7*l1extra_->isoEmEt.at(i)/12.) ){
                            IsoEGCenJet_cross_hist->Fill(l1extra_->isoEmEt.at(i), lumiWeight);
                            break;
                        }    
                    }
                }
//
                //IsoEGMET 13:7
                for(unsigned int i=0; i<l1extra_->nIsoEm;i++){
                    //std::cout << i << " " << iMET << " " << l1extra_->met.size() << std::endl;
                    if ( (fabs(l1extra_->isoEmEta.at(i)) <=2.17) && (l1extra_->isoEmBx.at(i)==0) && (l1extra_->metBx.at(iMET)==0) ){
                        if ( l1extra_->met.at(iMET) > floor(0.5+7*l1extra_->isoEmEt.at(i)/12.) ){
                            IsoEGMET_cross_hist->Fill(l1extra_->isoEmEt.at(i), lumiWeight);
                            break;
                        }    
                    }
                }
//
                ////TauTwoFwd
                //for(unsigned int i=0; i<l1extra_->nTauJets;i++){
                //    if ( (fabs(l1extra_->tauJetEta.at(i)) <= 2.17) && (l1extra_->tauJetBx.at(i)==0) ){
                //        //then loop over fwdJets and attach eta and pT requirements
                //    }
                //}    

                //===============UPGRADE OBJECTS===============//
                if (doUpgradeObj){
                    
                  //------Tower Jets------//
                  for(unsigned int i=0; i<l1upgrade_->nTowerJets;i++){
                    if(l1upgrade_->towerJetEt.at(i) > 0.){
                      if(l1upgrade_->towerJetBx.at(i)!=0) break;
                      if(i<4) up_towerJetEt_distro[i]->Fill(l1upgrade_->towerJetEt.at(i), lumiWeight);
                    }
                  }
                  //------nonIsoEm------//
                  for(unsigned int i=0; i<l1upgrade_->nNonIsoEm;i++){
                    double et = l1upgrade_->nonIsoEmEt.at(i);
                    if(et > 0.){
                      if(l1upgrade_->nonIsoEmBx.at(i)!=0) break;
                      if(i<4) up_nonIsoEmEt_distro[i]->Fill(et, lumiWeight);
                      
                      combEmUp.push_back(et);
                    }
                  }
                  
                  //------isoEm------//
                  for(unsigned int i=0; i<l1upgrade_->nIsoEm;i++){
                    double et = l1upgrade_->isoEmEt.at(i);
                    if(et > 0.){
                      if(l1upgrade_->isoEmBx.at(i)!=0) break;
                      if(i<4 && fabs(l1upgrade_->isoEmEta.at(i))<=2.17) up_isoEmEt_distro[i]->Fill(et, lumiWeight);
                      
                      combEmUp.push_back(et);
                    }
                  }
                  //------combEm------//
                  sort(combEmUp.begin(), combEmUp.end(), order_gt()); //sort them into descending order
                  for(unsigned int i=0; i<4; i++){
                    if(combEmUp.size() > i) up_combEGEt_distro[i]->Fill(combEmUp.at(i), lumiWeight);
                  }
        
                  //------nonIsoTau------//
                  for(unsigned int i=0; i<l1upgrade_->nNonIsoTaus;i++){
                    double et  = l1upgrade_->nonIsoTauEt.at(i);
                    double eta = l1upgrade_->nonIsoTauEta.at(i);
                    if(et > 0. && fabs(eta)<=2.17){
                      if(l1upgrade_->nonIsoTauBx.at(i)!=0) break;
                      if(i<4) up_nonIsoTauEt_distro[i]->Fill(et, lumiWeight);
                      pair<double, double> candTau(et, eta);
                      combTauUp.push_back(candTau);
                    }
                  }
                  
                  //------isoTau------//
                  for(unsigned int i=0; i<l1upgrade_->nIsoTaus;i++){
                    double et=l1upgrade_->isoTauEt.at(i);
                    double eta=l1upgrade_->isoTauEta.at(i);
                    if(et > 0. && fabs(eta) <= 2.17){
                      if(l1upgrade_->isoTauBx.at(i)!=0) break;
                      if(i<4) up_isoTauEt_distro[i]->Fill(et, lumiWeight);
                      
                      pair<double, double> candTau(et, eta);
                      combTauUp.push_back(candTau);
                    }
                  }
                  //------combTau------//
                  sort(combTauUp.begin(), combTauUp.end(), order_gt_pairs());
                  for(unsigned int i=0; i<4; i++){
                    if(combTauUp.size()>i) up_combTauEt_distro[i]->Fill(combTauUp.at(i).first, lumiWeight);
                  }
                
                  //------CrossTriggers------//
                  //doubleEG 13:7
                  if (combEmUp.size() > 1){
                     if (combEmUp.at(1) > floor(0.5+7.*combEmUp.at(0)/13.) ) up_doubleEG_cross_hist->Fill(combEmUp.at(0), lumiWeight);
                  }
  
                  //EG+mu 12:3.5  and   mu+EG 12:7
                  if ( (muon_hi.size() > 0) && (combEmUp.size() > 0) ){
                    if (muon_hi.at(0).first > floor(0.5 + (3.5*combEmUp.at(0))/12.)) up_EGMu_cross_hist->Fill(combEmUp.at(0), lumiWeight);
                    if (combEmUp.at(0) > floor(0.5 + (7.*muon_hi.at(0).first)/12.))     up_MuEG_cross_hist->Fill(muon_hi.at(0).first, lumiWeight);
                  }

                } //doUpgrade

                





				//---------------cumulative plot----------------
				//check triggers bools for filling cumulative plot
				if (PUScen == 45){
					//fill for menuv1 at 45PU
					if (combEm.size() > 0){
						if (combEm.at(0) > 34.)											trigRes[0]=true; //SingleEG34
					}
					for(unsigned int i=0; i<l1extra_->nIsoEm; i++){
						if (l1extra_->isoEmEt.at(i)>27.){
							if(fabs(l1extra_->isoEmEta.at(i)) <= 2.17) 					trigRes[1]=true; //SingleIsoEG27 er2.17
						}
					}
					if (combEm.size() > 1){
						if ((combEm.at(0) > 16.) && (combEm.at(1) > 9.))				trigRes[2]=true; //DoubleEG16,9 assym
					}
					if (muon_hi.size()>0){
						if (muon_hi.at(0).first > 120.) 									trigRes[3]=true; //SingleMu120
					}
					if (muon_hi.size() > 1){
						if ((muon_hi.at(0).first > 13.) && (muon_hi.at(1).first > 0.))	trigRes[4]=true; //DoubleMu13,open assym
					}
					if ((combEm.size()>0) && (muon_hi.size()>0)){
						if ((muon_hi.at(0).first > 5.) && (combEm.at(0) > 15.)) 		trigRes[5]=true; //EGMu 15,5. assym
						if ((muon_hi.at(0).first > 13.) && (combEm.at(0) > 8.))			trigRes[6]=true; //MuEG, 13:8 assym
					}
					if (combJets.size() > 0){
						if (combJets.at(0).first > 165.)								trigRes[7]=true; //SingleJet165
					}
					for(unsigned int i=0; i<combJets.size(); i++){
						int numFound = 0;
						if (combJets.at(i).first > 85.){
							if (fabs(combJets.at(i).second) <= 3.)	numFound++;
						}
						if (numFound >= 2) trigRes[8]=true; //DoubleJet85 er3
					}
					if (combJets.size() > 3){
						if (combJets.at(3).first > 60.)									trigRes[9]=true; //QuadJet60
					}
					for(unsigned int i=1; i<l1extra_->nTauJets; i++){ //DOUBLE
						int numFound = 0;
						if (l1extra_->tauJetEt.at(i) > 52.){
							if (fabs(l1extra_->tauJetEta.at(i)) <= 2.17) numFound++;
						}
						if (numFound >= 2) trigRes[10]=true; //DoubleTau52 er2.17
					}
					if (l1extra_->met.at(iMET) > 50.) 									trigRes[11]=true; //ETM50
					if (l1extra_->ht.at(iMHT) > 340.)									trigRes[12]=true; //HTT340

				}
				else if (PUScen == 66){
					//fill for menuv1 at 66PU
					if (combEm.size() > 0){
						if (combEm.at(0) > 37.)											trigRes[0]=true; //SingleEG37
					}
					for(unsigned int i=0; i<l1extra_->nIsoEm; i++){
						if (l1extra_->isoEmEt.at(i)>28.){
							if(fabs(l1extra_->isoEmEta.at(i)) <= 2.17) 					trigRes[1]=true; //SingleIsoEG28 er2.17
						}
					}
					if (combEm.size() > 1){
						if ((combEm.at(0) > 16.) && (combEm.at(1) > 9.))				trigRes[2]=true; //DoubleEG16,9 assym
					}
					if (muon_hi.size()>0){
						if (muon_hi.at(0).first > 120.) 									trigRes[3]=true; //SingleMu120
					}
					if (muon_hi.size() > 1){
						if ((muon_hi.at(0).first > 13.) && (muon_hi.at(1).first > 0.))	trigRes[4]=true; //DoubleMu13,open assym
					}
					if ((combEm.size()>0) && (muon_hi.size()>0)){
						if ((muon_hi.at(0).first > 5.) && (combEm.at(0) > 15.)) 		trigRes[5]=true; //EGMu, 15:5. assym
						if ((muon_hi.at(0).first > 13.) && (combEm.at(0) > 8.))			trigRes[6]=true; //MuEG, 13:8 assym
					}
					if (combJets.size() > 0){
						if (combJets.at(0).first > 192.)								trigRes[7]=true; //SingleJet192
					}
					for(unsigned int i=0; i<combJets.size(); i++){
						int numFound = 0;
						if (combJets.at(i).first > 92.){
							if (fabs(combJets.at(i).second) <= 3.)	numFound++;
						}
						if (numFound >= 2) trigRes[8]=true; //DoubleJet92 er3
					}
					if (combJets.size() > 3){
						if (combJets.at(3).first > 72.)									trigRes[9]=true; //QuadJet72
					}
					for(unsigned int i=1; i<l1extra_->nTauJets; i++){ //DOUBLE
						int numFound = 0;
						if (l1extra_->tauJetEt.at(i) > 50.){
							if (fabs(l1extra_->tauJetEta.at(i)) <= 2.17) numFound++;
						}
						if (numFound >= 2) trigRes[10]=true; //DoubleTau50 er2.17
					}
					if (l1extra_->met.at(iMET) > 60.) 									trigRes[11]=true; //ETM60
					if (l1extra_->ht.at(iMHT) > 550.)									trigRes[12]=true; //HTT550
				}

				//fill cumulative distribution
				for (int i=0; i<trigNum-1; i++){
					bool aboveFired = false;
					for(unsigned int j=0; j<i; j++){
						if(trigRes[i] && trigRes[j]) aboveFired = true;
					}
					if(aboveFired || (i==0 && trigRes[0])) cumul_rate->Fill(i, lumiWeight);
				}

				//fill the "any trigger fired" result
				if (trigRes[0] || trigRes[1] || trigRes[2] || trigRes[3] || trigRes[4] || trigRes[5] || trigRes[6] || trigRes[7] || trigRes[8] || trigRes[9] || trigRes[10] || trigRes[11] || trigRes[12]) trigRes[13] = true;

				//fill trigger results histo
				for(int i=0; i<trigNum; i++){
					for(int j=0; j<trigNum; j++){
						(trigRes[i] && trigRes[j]) ? trigOverlap->Fill(i, j, 1) : trigOverlap->Fill(i, j, 0);
					}
				}
			} //lumi binning
		} //zerobias trigger fired
	} //nevents

    
    // calculate the average PU of the sample
    sampAvgPU /= lumis.size();
    std::cout << "*** Average PU of sample: " << sampAvgPU << " ***" << std::endl;
    h_samplePU->Fill(sampAvgPU);


	std::cout << std::endl << "Total number of events: " << ntot << std::endl;

    //==============Create Rate Plots===============//

    
    for(int i=0;i<6;i++){
        if (i<4){
            DoRateCalc(isoEG_hist[i], isoEG_rate[i], preScale, lumis.size());
            DoRateCalc(combEG_hist[i], combEG_rate[i], preScale, lumis.size());
            DoRateCalc(muon_hist_hi[i], muon_rate_hi[i], preScale, lumis.size());
            DoRateCalc(tau_hist[i],tau_rate[i], preScale, lumis.size());
            DoRateCalc(cenpTauJets_distro[i], cenpTauJets_rate[i], preScale, lumis.size());
            DoRateCalc(fwdJets_distro[i], fwdJets_rate[i], preScale, lumis.size());
            DoRateCalc(cenJets_distro[i], cenJets_rate[i], preScale, lumis.size());
        }
        DoRateCalc(combJetsEt_hist[i], combJetsEt_rate[i], preScale, lumis.size());
        DoRateCalc(combJetsEr_hist[i], combJetsEr_rate[i], preScale, lumis.size());
    }
    
    DoRateCalc(met_hist, met_rate, preScale, lumis.size());
    DoRateCalc(ett_hist, ett_rate, preScale, lumis.size());
    DoRateCalc(mht_hist, mht_rate, preScale, lumis.size());
    DoRateCalc(htt_hist, htt_rate, preScale, lumis.size());
    DoRateCalc(sinMuEr_hist, sinMuEr_rate, preScale, lumis.size());
    DoRateCalc(diJetEr_hist, diJetEr_rate, preScale, lumis.size());
    DoRateCalc(doubleEG_cross_hist, doubleEG_cross_rate, preScale, lumis.size());
    DoRateCalc(doubleMu_cross_hist, doubleMu_cross_rate, preScale, lumis.size());
    DoRateCalc(EGMu_cross_hist, EGMu_cross_rate, preScale, lumis.size());
    DoRateCalc(MuEG_cross_hist, MuEG_cross_rate, preScale, lumis.size());
    DoRateCalc(TauMu_cross_hist, TauMu_cross_rate, preScale, lumis.size());
    DoRateCalc(TauEG_cross_hist, TauEG_cross_rate, preScale, lumis.size());
    DoRateCalc(IsoEGCenJet_cross_hist, IsoEGCenJet_cross_rate, preScale, lumis.size());
    DoRateCalc(IsoEGMET_cross_hist, IsoEGMET_cross_rate, preScale, lumis.size());
    DoRateCalc(TauTwoFwd_cross_hist, TauTwoFwd_cross_rate, preScale, lumis.size());

    

    if (doUpgradeObj){
        
        for(int i=0; i<4; i++){
            DoRateCalc(up_towerJetEt_distro[i], up_towerJetEt_rate[i], preScale, lumis.size());
            DoRateCalc(up_isoEmEt_distro[i], up_isoEmEt_rate[i], preScale, lumis.size());
            DoRateCalc(up_nonIsoEmEt_distro[i], up_nonIsoEmEt_rate[i], preScale, lumis.size());
            DoRateCalc(up_isoTauEt_distro[i], up_isoTauEt_rate[i], preScale, lumis.size());
            DoRateCalc(up_nonIsoTauEt_distro[i], up_nonIsoTauEt_rate[i], preScale, lumis.size());    
            DoRateCalc(up_combTauEt_distro[i], up_combTauEt_rate[i], preScale, lumis.size());
            DoRateCalc(up_combEGEt_distro[i], up_combEGEt_rate[i], preScale, lumis.size());
        }

        DoRateCalc(up_doubleEG_cross_hist, up_doubleEG_cross_rate, preScale, lumis.size());
        DoRateCalc(up_EGMu_cross_hist, up_EGMu_cross_rate, preScale, lumis.size());
        DoRateCalc(up_MuEG_cross_hist, up_MuEG_cross_rate, preScale, lumis.size());

    }


	std::cout << "Writing and closing the file." << std::endl;


	outFile->Write();
	outFile->Close();

	return lumis.size();

} 	
