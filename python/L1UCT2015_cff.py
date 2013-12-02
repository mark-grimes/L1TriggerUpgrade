import FWCore.ParameterSet.Config as cms

from L1Trigger.UCT2015.Lut import *

hackHCALMIPs = cms.EDProducer(
    "HcalTpgMipEmbedder",
    src = cms.InputTag("hcalDigis"),
    threshold = cms.double(3), # In GeV
    rawThreshold = cms.uint32(3), # In TPG rank
    cutOnRawBits = cms.bool(False), # What to cut on
)

uctDigis = cms.EDProducer(
    "L1RCTProducer",
    #hcalDigis = cms.VInputTag(cms.InputTag("hcalDigis")),
    hcalDigis = cms.VInputTag(cms.InputTag("hackHCALMIPs")),
    useEcal = cms.bool(True),
    useHcal = cms.bool(True),
    ecalDigis = cms.VInputTag(cms.InputTag("ecalDigis:EcalTriggerPrimitives")),
    BunchCrossings = cms.vint32(0),
    getFedsFromOmds = cms.bool(False),
    queryDelayInLS = cms.uint32(10),
    queryIntervalInLS = cms.uint32(100)#,
)

UCT2015Producer = cms.EDProducer(
    "UCT2015Producer",
    puCorrect = cms.bool(True),
    useUICrho = cms.bool(True),
    # All of these uint32 thresholds are in GeV.
    puETMax = cms.uint32(7),
    regionETCutForHT = cms.uint32(5),
    regionETCutForMET = cms.uint32(0),
    minGctEtaForSums = cms.uint32(4),
    maxGctEtaForSums = cms.uint32(17),
    jetSeed = cms.uint32(5),
    egtSeed = cms.uint32(5),
    relativeIsolationCut = cms.double(1.0),
    relativeJetIsolationCut = cms.double(1.0),
    egammaLSB = cms.double(1.0), # This has to correspond with the value from L1CaloEmThresholds
    regionLSB = RCTConfigProducers.jetMETLSB,
    eicIsolationThreshold = cms.double(3.)
)

from L1Trigger.UCT2015.uct2015L1ExtraParticles_cfi import *

L1UCT2015 = cms.Sequence(
    hackHCALMIPs
    + uctDigis
    + UCT2015Producer
    + uct2015L1ExtraParticles
    )
