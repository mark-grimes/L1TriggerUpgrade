#ifndef __L1Analysis_L1AnalysisL1Upgrade_H__
#define __L1Analysis_L1AnalysisL1Upgrade_H__

#include <TTree.h>
#include <TMatrixD.h>

namespace L1Analysis
{
  class L1AnalysisL1Upgrade
{

  public : 
  void initTree(TTree * tree, const std::string & className);

  public:
  L1AnalysisL1Upgrade() {}
  void print(); 
  bool check();   
  
    // ---- L1AnalysisL1Upgrade information.
    unsigned nIsoEm;
    std::vector<double> isoEmEt;
    std::vector<double> isoEmEta;
    std::vector<double> isoEmPhi;
    std::vector<int>    isoEmBx;

    unsigned nNonIsoEm;
    std::vector<double> nonIsoEmEt;
    std::vector<double> nonIsoEmEta;
    std::vector<double> nonIsoEmPhi;
    std::vector<int>    nonIsoEmBx;

    unsigned nTowerJets;
    std::vector<double> towerEt;
    std::vector<double> towerEta;
    std::vector<double> towerPhi; 
    std::vector<int>    towerBx;

    unsigned nIsoTaus;
    std::vector<double> isoTauEt;
    std::vector<double> isoTauEta;
    std::vector<double> isoTauPhi;
    std::vector<int>    isoTauBx;

    unsigned nTaus;
    std::vector<double> tauEt;
    std::vector<double> tauEta;
    std::vector<double> tauPhi;
    std::vector<int>    tauBx;
    
    unsigned nMuons;
    std::vector<double>   muonEt;
    std::vector<double>   muonEta;
    std::vector<double>   muonPhi;
    std::vector<int>      muonChg;
    std::vector<unsigned int> muonIso;
    std::vector<unsigned int> muonFwd;
    std::vector<unsigned int> muonMip;
    std::vector<unsigned int> muonRPC;
    std::vector<int>          muonBx;

    std::vector<double> hfEtSum;
    std::vector<unsigned int> hfBitCnt; 
    std::vector<int>          hfBx;
 
    double met;
    double metPhi;
    int metsBx;
    double mht;
    double mhtPhi;
    int mhtsBx;
    double et;
    double ht;
};
}

#endif

#ifdef l1ntuple_cxx



void L1Analysis::L1AnalysisL1Upgrade::initTree(TTree * tree, const std::string & className)
{
  //note: names do not necessarily match due to setup of top-level config for creating ntuples
   SetBranchAddress(tree, "nIsoEm", className,  &nIsoEm);
   SetBranchAddress(tree, "isoEmEt", className,  &isoEmEt);
   SetBranchAddress(tree, "isoEmEta", className,  &isoEmEta);
   SetBranchAddress(tree, "isoEmPhi", className,  &isoEmPhi); 
   SetBranchAddress(tree, "isoEmBx", className,  &isoEmBx);
   SetBranchAddress(tree, "nNonIsoEm", className,  &nNonIsoEm);
   SetBranchAddress(tree, "nonIsoEmEt", className,  &nonIsoEmEt);
   SetBranchAddress(tree, "nonIsoEmEta", className,  &nonIsoEmEta);
   SetBranchAddress(tree, "nonIsoEmPhi", className,  &nonIsoEmPhi); 
   SetBranchAddress(tree, "nonIsoEmBx", className,  &nonIsoEmBx);
   SetBranchAddress(tree, "nCenJets", className,     &nTowerJets);
   SetBranchAddress(tree, "cenJetEt", className,     &towerJetEt);
   SetBranchAddress(tree, "cenJetEta", className,    &towerJetEta);
   SetBranchAddress(tree, "cenJetPhi", className,    &towerJetPhi);  
   SetBranchAddress(tree, "cenJetBx", className,     &towerJetBx);
   SetBranchAddress(tree, "nFwdJets", className,     &nIsoTaus);
   SetBranchAddress(tree, "fwdJetEt", className,     &isoTauEt);
   SetBranchAddress(tree, "fwdJetEta", className,    &isoTauEta);
   SetBranchAddress(tree, "fwdJetPhi", className,    &isoTauPhi); 
   SetBranchAddress(tree, "fwdJetBx",  className,    &isoTauBx);
   SetBranchAddress(tree, "nTauJets", className,     &nNonIsoTaus);
   SetBranchAddress(tree, "tauJetEt", className,     &nonIsoTauEt);
   SetBranchAddress(tree, "tauJetEta", className,    &nonIsoTauEta);
   SetBranchAddress(tree, "tauJetPhi", className,    &nonIsoTauPhi); 
   SetBranchAddress(tree, "tauJetBx",  className,    &nonIsoTauBx);
   SetBranchAddress(tree, "nMuons", className,       &nMuons);
   SetBranchAddress(tree, "muonEt", className,       &muonEt);
   SetBranchAddress(tree, "muonEta", className,      &muonEta);
   SetBranchAddress(tree, "muonPhi", className,      &muonPhi);
   SetBranchAddress(tree, "muonChg", className,      &muonChg);
   SetBranchAddress(tree, "muonIso", className,      &muonIso);
   SetBranchAddress(tree, "muonFwd", className,      &muonFwd);
   SetBranchAddress(tree, "muonMip", className,      &muonMip);
   SetBranchAddress(tree, "muonRPC", className,      &muonRPC); 
   SetBranchAddress(tree, "muonBx",  className,      &muonBx);
   SetBranchAddress(tree, "hfEtSum", className,      &hfEtSum);
   SetBranchAddress(tree, "hfBitCnt", className,     &hfBitCnt); 
   SetBranchAddress(tree, "hfBx", className,         &hfBx);
   SetBranchAddress(tree, "met", className,          &met);
   SetBranchAddress(tree, "metPhi", className,       &metPhi);
   SetBranchAddress(tree, "metsBx", className,       &metsBx);
   SetBranchAddress(tree, "mht", className,          &mht);
   SetBranchAddress(tree, "mhtPhi", className,       &mhtPhi);  
   SetBranchAddress(tree, "mhtsBx", className,       &mhtsBx);
   SetBranchAddress(tree, "et", className,           &et);
   SetBranchAddress(tree, "ht", className,           &ht);
}


void L1Analysis::L1AnalysisL1Upgrade::print()
{
}

bool L1Analysis::L1AnalysisL1Upgrade::check()
{
  bool test=true;
  return test;
}


#endif


