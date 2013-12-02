import FWCore.ParameterSet.Config as cms

# import master L1 Ntuple config
from UserCode.L1TriggerUpgrade.l1CurrentNtupleFromRAW import *

print "Adding Stage 1 algorithms"

# add Stage 1 algorithms
process.load('UserCode.L1TriggerUpgrade.L1UCT2015_cff')
process.uct2015L1ExtraParticles.muonSource = cms.InputTag("valGmtDigis")
process.uct2015L1ExtraParticles.produceMuonParticles = cms.bool(True)

# load the upgrade tree
process.load('UserCode.L1TriggerUpgrade.l1ExtraUpgradeTreeProducer_cfi')
process.l1ExtraUpgradeTreeProducer.muonLabel = cms.untracked.InputTag("uct2015L1ExtraParticles")
process.l1ExtraUpgradeTreeProducer.egLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "Relaxed")
process.l1ExtraUpgradeTreeProducer.isoEGLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "Isolated")
process.l1ExtraUpgradeTreeProducer.tauLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "RelaxedTau")
process.l1ExtraUpgradeTreeProducer.isoTauLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "IsolatedTau")
process.l1ExtraUpgradeTreeProducer.jetLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "Jets")
process.l1ExtraUpgradeTreeProducer.fwdJetLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "FwdJets")
process.l1ExtraUpgradeTreeProducer.metLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "MET")
process.l1ExtraUpgradeTreeProducer.mhtLabel = cms.untracked.InputTag("uct2015L1ExtraParticles", "MHT")

# add upgrade algorithms to the path
process.p += process.L1UCT2015
process.p += process.l1ExtraUpgradeTreeProducer    # upgrade candidates

print "Done"
