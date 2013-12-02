import FWCore.ParameterSet.Config as cms

# import master L1 Ntuple config
from UserCode.L1TriggerUpgrade.l1CurrentNtupleFromRAW import *

# add Stage 1 EG, tau, muon isolation
print "Adding Stage 1 EG, tau, muon iso algorithms"
process.load('UserCode.L1TriggerUpgrade.L1UCT2015_cff')
process.uct2015L1ExtraParticles.muonSource = cms.InputTag("valGmtDigis")
process.uct2015L1ExtraParticles.produceMuonParticles = cms.bool(True)

process.p += process.L1UCT2015

# add Stage 2 jets
print "Adding Stage 2 jet, sum algorithms"
process.load('SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff')
process.L1CaloTowerProducer.ECALDigis = cms.InputTag("ecalDigis:EcalTriggerPrimitives")
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("hcalDigis")

process.p += process.L1CaloTowerProducer
process.p += process.L1TowerJetProducer
process.p += process.L1TowerJetFilter1D
process.p += process.L1TowerJetFilter2D
process.p += process.L1TowerFwdJetProducer
process.p += process.L1TowerFwdJetFilter1D
process.p += process.L1TowerFwdJetFilter2D
process.p += process.L1CalibFilterTowerJetProducer

# load the upgrade tree
process.load('UserCode.L1TriggerUpgrade.l1ExtraUpgradeTreeProducer_cfi')
process.l1ExtraUpgradeTreeProducer.muonLabel = cms.untracked.InputTag("uct2015L1ExtraParticles")
process.l1ExtraUpgradeTreeProducer.egLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "Relaxed")
process.l1ExtraUpgradeTreeProducer.isoEGLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "Isolated")
process.l1ExtraUpgradeTreeProducer.tauLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "RelaxedTauLeadTrk")
process.l1ExtraUpgradeTreeProducer.isoTauLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "IsolatedTau")

process.p += process.l1ExtraUpgradeTreeProducer

print "Done"

