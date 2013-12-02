import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1Stage1NtupleFromRAW import *

from UserCode.L1TriggerUpgrade.MCSetup import *
mcSetup(process, True, False)

# job options
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
                )

readFiles.extend( [
     '/store/mc/Summer12/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/UpgradeL1TDR-PU35_POSTLS161_V12-v2/00000/168A4B24-C53C-E211-A370-003048C68AA6.root'
] )

# Debug stuff

#process.output = cms.OutputModule(
#    "PoolOutputModule",
#    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('output.root')
#    )

#process.e = cms.EndPath(process.output)

#process.load("L1TriggerConfig.GctConfigProducers.l1GctConfigDump_cfi")
#process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
#process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
#process.MessageLogger.debugModules = cms.untracked.vstring('l1GctConfigDump')
#process.p += ( process.l1GctConfigDump )
