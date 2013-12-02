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

    L1AnalysisL1ExtraUpgradeDataFormat * getData() {return &l1extra_;}

  private :
    L1AnalysisL1ExtraUpgradeDataFormat l1extra_;
  }; 
}
#endif


