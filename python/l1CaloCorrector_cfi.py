import FWCore.ParameterSet.Config as cms

l1CaloCorrector = cms.EDProducer(
    "L1CaloCorrector",
    hcalSource = cms.untracked.InputTag("hcalDigis"),
    hcalFile   = cms.FileInPath("UserCode/L1TriggerUpgrade/data/HCALCorrectionPU45.root")
    )
