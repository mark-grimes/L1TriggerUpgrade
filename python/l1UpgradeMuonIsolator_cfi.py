import FWCore.ParameterSet.Config as cms

l1UpgradeMuonIsolator = cms.EDProducer(
    'L1UpgradeMuonIsolator',
    muonSource = cms.InputTag("valGmtDigis"),
    regionSource = cms.InputTag("gctDigis"),
    rgnThreshold = cms.double(0.),
    regionLSB = cms.double(0.5),
    muFixedIso = cms.double(-999.),
    muRelIso = cms.double(0.35)                          
)
