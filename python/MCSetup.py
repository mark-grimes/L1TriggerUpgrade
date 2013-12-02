import FWCore.ParameterSet.Config as cms

def mcSetup(process, stage1=False, stage2=False):

    print "Going to run on MC"

    # adopt to MC
    process.p.remove(process.RawToDigi)
    process.load('Configuration.StandardSequences.RawToDigi_cff')
    process.p.insert(0, process.RawToDigi)

    process.l1NtupleProducer.hltSource        = cms.InputTag("none")
    process.l1NtupleProducer.verbose          = cms.untracked.bool(True)
    process.l1NtupleProducer.simulationSource = cms.InputTag("simulation")

    # bug fix for missing HCAL TPs in MC RAW
    from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
    HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
    process.p.insert(1, process.valHcalTriggerPrimitiveDigis)
    process.l1NtupleProducer.hcalSource       = cms.InputTag("valHcalTriggerPrimitiveDigis")
    process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))

    if (stage1):    
        process.hackHCALMIPs.src = cms.InputTag("valHcalTriggerPrimitiveDigis")
        from L1Trigger.UCT2015.Lut import eg_calib_v1, eg_calib_v4
        process.RCTConfigProducers.eGammaECalScaleFactors = eg_calib_v1
        process.RCTConfigProducers.jetMETECalScaleFactors = eg_calib_v1

    if (stage2):
        process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")
