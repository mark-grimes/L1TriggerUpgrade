#ifndef __L1Analysis_L1AnalysisL1ExtraUpgrade_H__
#define __L1Analysis_L1AnalysisL1ExtraUpgrade_H__

//-------------------------------------------------------------------------------
// Created 02/03/2010 - A.C. Le Bihan
// 
// 
// Original code : UserCode/L1TriggerDPG/L1ExtraTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"

#include "L1AnalysisL1ExtraUpgradeDataFormat.h"

namespace L1Analysis
{
  class L1AnalysisL1ExtraUpgrade 
  {
  public:
    L1AnalysisL1ExtraUpgrade();
    ~L1AnalysisL1ExtraUpgrade();
    void Reset() {l1extra_.Reset();}
    void SetEG      (const edm::Handle<l1extra::L1EmParticleCollection>   eg,     unsigned maxL1Extra);
    void SetIsoEG   (const edm::Handle<l1extra::L1EmParticleCollection>   isoEG,  unsigned maxL1Extra);
    void SetTau     (const edm::Handle<l1extra::L1JetParticleCollection>  tau,    unsigned maxL1Extra);
    void SetIsoTau  (const edm::Handle<l1extra::L1JetParticleCollection>  isoTau, unsigned maxL1Extra);

    void SetJet     (const edm::Handle<l1extra::L1JetParticleCollection>  jet,    unsigned maxL1Extra);
    void SetFwdJet  (const edm::Handle<l1extra::L1JetParticleCollection>  fwdJet, unsigned maxL1Extra);
    void SetMuon    (const edm::Handle<l1extra::L1MuonParticleCollection> muon,   unsigned maxL1Extra);
    void SetMet     (const edm::Handle<l1extra::L1EtMissParticleCollection> mets);
    void SetMht     (const edm::Handle<l1extra::L1EtMissParticleCollection> mhts);

    // Add L1TrackTriggerObjects
    void SetTkEM   (const edm::Handle<l1extra::L1TkEmParticleCollection>   tkEM,     unsigned maxL1Extra);
    void SetTkEG   (const edm::Handle<l1extra::L1TkElectronParticleCollection>   tkEG,     unsigned maxL1Extra);
    void SetTkIsoEG(const edm::Handle<l1extra::L1TkElectronParticleCollection>   tkIsoEG,     unsigned maxL1Extra);
    void SetTkMuon (const edm::Handle<l1extra::L1TkMuonParticleCollection> tkMuon,   unsigned maxL1Extra);
    void SetTkJet  (const edm::Handle<l1extra::L1TkJetParticleCollection>  tkJet,    unsigned maxL1Extra);
    void SetTkMet  (const edm::Handle<l1extra::L1TkEtMissParticleCollection> tkMets);

    L1AnalysisL1ExtraUpgradeDataFormat * getData() {return &l1extra_;}

  private :
    L1AnalysisL1ExtraUpgradeDataFormat l1extra_;
  }; 
}
#endif


