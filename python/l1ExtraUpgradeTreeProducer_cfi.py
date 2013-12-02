import FWCore.ParameterSet.Config as cms

l1ExtraUpgradeTreeProducer = cms.EDAnalyzer("L1ExtraUpgradeTreeProducer",
   egLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:EGamma"),
   isoEGLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:IsoEGamma"),
   tauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:Taus"),
   isoTauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:IsoTaus"),
   jetLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer:Cen8x8"),
   fwdJetLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer:Fwd8x8"),
   muonLabel = cms.untracked.InputTag("l1UpgradeMuonIsolator"),
   metLabel = cms.untracked.InputTag("rawSLHCL1ExtraParticles:MET"),
   mhtLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer:TowerMHT"),
   maxL1Extra = cms.uint32(20)
)
