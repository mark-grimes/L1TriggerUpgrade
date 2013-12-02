import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1FallbackNtupleFromRAW import *

print "Going to run on data"

#process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(10000)

process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
                )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# DTTF re-emulation
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("L1MuDTTFParametersRcd"),
             tag = cms.string("L1MuDTTFParameters_dttf12_TSC_03_csc_col_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
             )
)

# RPC re-emulation
process.load("L1TriggerConfig.RPCTriggerConfig.RPCBxOrConfig_cff")
process.load("L1Trigger.RPCTrigger.RPCConeConfig_cff")
process.load("L1Trigger.RPCTrigger.l1RpcEmulDigis_cfi")
process.l1RPCBxOrConfig.lastBX = cms.int32(0)
process.l1RPCBxOrConfig.firstBX = cms.int32(0)

readFiles.extend( [
    'file:/gpfs_phys/storm/cms/data/Run2012C/ZeroBias/RAW/v1/000/198/588/1ADBB979-1ACA-E111-AA51-003048D3C932.root'
    ] )

#process.output = cms.OutputModule(
#    "PoolOutputModule",
#    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('output.root')
#    )

#process.e = cms.EndPath(process.output)
