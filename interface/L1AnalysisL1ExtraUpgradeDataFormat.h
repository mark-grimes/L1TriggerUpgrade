#ifndef __L1Analysis_L1AnalysisL1ExtraUpgradeDataFormat_H__
#define __L1Analysis_L1AnalysisL1ExtraUpgradeDataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>

namespace L1Analysis
{
  struct L1AnalysisL1ExtraUpgradeDataFormat
  {
    L1AnalysisL1ExtraUpgradeDataFormat(){Reset();};
    ~L1AnalysisL1ExtraUpgradeDataFormat(){};
    
    void Reset()
    {
      nEG = 0;
      egEt.clear();
      egEta.clear();
      egPhi.clear();
      egBx.clear();
      
      nIsoEG = 0;
      isoEGEt.clear();
      isoEGEta.clear();
      isoEGPhi.clear();
      isoEGBx.clear();

      nTau = 0;
      tauEt.clear();
      tauEta.clear();
      tauPhi.clear(); 
      tauBx.clear();

      nIsoTau = 0;
      isoTauEt.clear();
      isoTauEta.clear();
      isoTauPhi.clear(); 
      isoTauBx.clear();

      nJets = 0;
      jetEt.clear();
      jetEta.clear();
      jetPhi.clear();
      jetBx.clear();

      nFwdJets = 0;
      fwdJetEt.clear();
      fwdJetEta.clear();
      fwdJetPhi.clear();
      fwdJetBx.clear();

      nMuons = 0;
      muonEt.clear();
      muonEta.clear();
      muonPhi.clear();
      muonChg.clear();
      muonIso.clear();
      muonFwd.clear();
      muonMip.clear();
      muonRPC.clear();
      muonBx.clear();
      muonQuality.clear();

      nMet = 0;
      et.clear();
      met.clear();
      metPhi.clear();
      metBx.clear();

      nMht = 0;
      ht.clear();
      mht.clear();
      mhtPhi.clear();
      mhtBx.clear();

    }
   
    unsigned int nEG;
    std::vector<double> egEt;
    std::vector<double> egEta;
    std::vector<double> egPhi;
    std::vector<int>    egBx;
 
    unsigned int nIsoEG;
    std::vector<double> isoEGEt;
    std::vector<double> isoEGEta;
    std::vector<double> isoEGPhi;
    std::vector<int>    isoEGBx;
 
    unsigned int nTau;
    std::vector<double> tauEt;
    std::vector<double> tauEta;
    std::vector<double> tauPhi;
    std::vector<int>    tauBx;

    unsigned int nIsoTau;
    std::vector<double> isoTauEt;
    std::vector<double> isoTauEta;
    std::vector<double> isoTauPhi;
    std::vector<int>    isoTauBx;

    unsigned int nJets;
    std::vector<double> jetEt;
    std::vector<double> jetEta;
    std::vector<double> jetPhi;
    std::vector<int>    jetBx;
 
    unsigned int nFwdJets;
    std::vector<double> fwdJetEt;
    std::vector<double> fwdJetEta;
    std::vector<double> fwdJetPhi;
    std::vector<int>    fwdJetBx;

    unsigned int nMuons;
    std::vector<double>   muonEt;
    std::vector<double>   muonEta;
    std::vector<double>   muonPhi;
    std::vector<int>      muonChg;
    std::vector<unsigned int> muonIso;
    std::vector<unsigned int> muonFwd;
    std::vector<unsigned int> muonMip;
    std::vector<unsigned int> muonRPC;
    std::vector<int>      muonBx;
    std::vector<int>      muonQuality;
 
    unsigned int nMet;
    std::vector<double> et;
    std::vector<double> met;
    std::vector<double> metPhi;
    std::vector<double> metBx;

    unsigned int nMht;
    std::vector<double> ht;
    std::vector<double> mht;
    std::vector<double> mhtPhi;
    std::vector<double> mhtBx;

  }; 
}
#endif


