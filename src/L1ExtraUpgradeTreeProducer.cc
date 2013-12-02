// -*- C++ -*-
//
// Package:    UserCode/L1TriggerDPG
// Class:      L1ExtraUpgradeTreeProducer
// 
/**\class L1ExtraUpgradeTreeProducer L1ExtraUpgradeTreeProducer.cc UserCode/L1TriggerDPG/src/L1ExtraUpgradeTreeProducer.cc

Description: Produce L1 Extra tree

Implementation:
     
*/
//
// Original Author:  Alex Tapper
//         Created:  
// $Id: L1ExtraUpgradeTreeProducer.cc,v 1.5 2013/01/06 21:55:55 jbrooke Exp $
//
//


// system include files
#include <memory>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// data formats
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"


// ROOT output stuff
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisL1ExtraUpgrade.h"

//
// class declaration
//

class L1ExtraUpgradeTreeProducer : public edm::EDAnalyzer {
public:
  explicit L1ExtraUpgradeTreeProducer(const edm::ParameterSet&);
  ~L1ExtraUpgradeTreeProducer();
  
  
private:
  virtual void beginJob(void) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

public:
  
  L1Analysis::L1AnalysisL1ExtraUpgrade* l1Extra;
  L1Analysis::L1AnalysisL1ExtraUpgradeDataFormat * l1ExtraData;

private:

  unsigned maxL1Extra_;

  // output file
  edm::Service<TFileService> fs_;
  
  // tree
  TTree * tree_;
 
  // EDM input tags
  edm::InputTag egLabel_;
  edm::InputTag isoEGLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag isoTauLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag fwdJetLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag metLabel_;
  edm::InputTag mhtLabel_;

};



L1ExtraUpgradeTreeProducer::L1ExtraUpgradeTreeProducer(const edm::ParameterSet& iConfig):
  egLabel_(iConfig.getUntrackedParameter("egLabel",edm::InputTag("SLHCL1ExtraParticles:EGamma"))),
  isoEGLabel_(iConfig.getUntrackedParameter("isoEGLabel",edm::InputTag("SLHCL1ExtraParticles:IsoEGamma"))),
  tauLabel_(iConfig.getUntrackedParameter("tauLabel",edm::InputTag("SLHCL1ExtraParticles:Taus"))),
  isoTauLabel_(iConfig.getUntrackedParameter("isoTauLabel",edm::InputTag("SLHCL1ExtraParticles:IsoTaus"))),
  jetLabel_(iConfig.getUntrackedParameter("jetLabel",edm::InputTag("calibJetProducer:Tower"))),
  fwdJetLabel_(iConfig.getUntrackedParameter("fwdJetLabel",edm::InputTag(""))),
  muonLabel_(iConfig.getUntrackedParameter("muonLabel",edm::InputTag(""))),
  metLabel_(iConfig.getUntrackedParameter("metLabel",edm::InputTag(""))),
  mhtLabel_(iConfig.getUntrackedParameter("mhtLabel",edm::InputTag("")))
{
 
  maxL1Extra_ = iConfig.getParameter<unsigned int>("maxL1Extra");
 
  l1Extra     = new L1Analysis::L1AnalysisL1ExtraUpgrade();
  l1ExtraData = l1Extra->getData();
  
  // set up output
  tree_=fs_->make<TTree>("L1ExtraUpgradeTree", "L1ExtraUpgradeTree");
  tree_->Branch("L1ExtraUpgrade", "L1Analysis::L1AnalysisL1ExtraUpgradeDataFormat", &l1ExtraData, 32000, 3);

}


L1ExtraUpgradeTreeProducer::~L1ExtraUpgradeTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1ExtraUpgradeTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  l1Extra->Reset();

  edm::Handle<l1extra::L1EmParticleCollection> eg;
  edm::Handle<l1extra::L1EmParticleCollection> isoEG;
  edm::Handle<l1extra::L1JetParticleCollection> tau;
  edm::Handle<l1extra::L1JetParticleCollection> isoTau;
  edm::Handle<l1extra::L1JetParticleCollection> jet;
  edm::Handle<l1extra::L1JetParticleCollection> fwdJet;
  edm::Handle<l1extra::L1MuonParticleCollection> muon; ;
  edm::Handle<l1extra::L1EtMissParticleCollection> mets;
  edm::Handle<l1extra::L1EtMissParticleCollection> mhts;

  iEvent.getByLabel(egLabel_, eg);
  iEvent.getByLabel(isoEGLabel_, isoEG);
  iEvent.getByLabel(tauLabel_, tau);
  iEvent.getByLabel(isoTauLabel_, isoTau);
  iEvent.getByLabel(jetLabel_, jet);
  iEvent.getByLabel(fwdJetLabel_, fwdJet);
  iEvent.getByLabel(muonLabel_, muon);
  iEvent.getByLabel(metLabel_, mets);
  iEvent.getByLabel(mhtLabel_, mhts);

  if (eg.isValid()){ 
    l1Extra->SetEG(eg, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade EG not found (" << egLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (isoEG.isValid()){ 
    l1Extra->SetIsoEG(isoEG, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Iso EG not found (" << isoEGLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (tau.isValid()){ 
    l1Extra->SetTau(tau, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Tau not found (" << tauLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (isoTau.isValid()){ 
    l1Extra->SetIsoTau(isoTau, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Iso Tau not found (" << isoTauLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (jet.isValid()){ 
    l1Extra->SetJet(jet, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Jets not found (" << jetLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (fwdJet.isValid()){ 
    l1Extra->SetFwdJet(fwdJet, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Fwd Jets not found (" << fwdJetLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (muon.isValid()){ 
    l1Extra->SetMuon(muon, maxL1Extra_);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade Muons not found (" << muonLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (mets.isValid()){ 
    l1Extra->SetMet(mets);
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade MET not found (" << metLabel_ << "). Branch will not be filled" << std::endl;
  }

  if (mhts.isValid()){ 
    l1Extra->SetMht(mhts);  
  } else {
    edm::LogWarning("MissingProduct") << "L1ExtraUpgrade MHT not found (" << mhtLabel_ << "). Branch will not be filled" << std::endl;
  }

  tree_->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1ExtraUpgradeTreeProducer::beginJob(void)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1ExtraUpgradeTreeProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1ExtraUpgradeTreeProducer);
