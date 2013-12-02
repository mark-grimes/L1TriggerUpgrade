import FWCore.ParameterSet.Config as cms


# import upgrade simulation configs
from SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff import *
from SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTriggerAnalysisOnData_cfi import *

L1CaloTowerProducer.ECALDigis = cms.InputTag("ecalDigis:EcalTriggerPrimitives")
L1CaloTowerProducer.HCALDigis =  cms.InputTag("hcalDigis")

# EG/Tau algorithms
from L1TriggerConfig.RCTConfigProducers.L1RCTConfig_cff import RCTConfigProducers
RCTConfigProducers.eMaxForHoECut = cms.double(60.0)
RCTConfigProducers.hOeCut = cms.double(0.05)
RCTConfigProducers.eGammaECalScaleFactors = cms.vdouble(1.0, 1.01, 1.02, 1.02, 1.02,
                                                                1.06, 1.04, 1.04, 1.05, 1.09,
                                                                1.1, 1.1, 1.15, 1.2, 1.27,
                                                                1.29, 1.32, 1.52, 1.52, 1.48,
                                                                1.4, 1.32, 1.26, 1.21, 1.17,
                                                                1.15, 1.15, 1.15)
RCTConfigProducers.eMinForHoECut = cms.double(3.0)
RCTConfigProducers.hActivityCut = cms.double(4.0)
RCTConfigProducers.eActivityCut = cms.double(4.0)
RCTConfigProducers.jetMETHCalScaleFactors = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0)
RCTConfigProducers.eicIsolationThreshold = cms.uint32(6)
RCTConfigProducers.etMETLSB = cms.double(0.25)
RCTConfigProducers.jetMETECalScaleFactors = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0, 1.0, 1.0,
                                                                1.0, 1.0, 1.0)
RCTConfigProducers.eMinForFGCut = cms.double(100.0)
RCTConfigProducers.eGammaLSB = cms.double(0.25)



