import FWCore.ParameterSet.Config as cms

print "Adding current L1 information"

# make L1 ntuples from RAW+RECO

process = cms.Process("L1NTUPLE")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

# L1 raw to digi options
process.gctDigis.numberOfGctSamplesToUnpack = cms.uint32(5)
process.l1extraParticles.centralBxOnly = cms.bool(False)

# output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Tree.root')
)

# L1 emulator
process.load("UserCode.L1TriggerDPG.L1EmulatorTree_cff")

# L1 ntuple producers
process.load("UserCode.L1TriggerDPG.l1NtupleProducer_cfi")
process.l1NtupleProducer.hltSource = cms.InputTag("none")
process.load("UserCode.L1TriggerDPG.l1ExtraTreeProducer_cfi")
process.load("EventFilter.L1GlobalTriggerRawToDigi.l1GtTriggerMenuLite_cfi")
process.load("UserCode.L1TriggerDPG.l1MenuTreeProducer_cfi")

# remove CSCTF from ntuples due to crashes!
process.l1NtupleProducer.csctfTrkSource       = cms.InputTag("none")
process.l1NtupleProducer.csctfLCTSource       = cms.InputTag("none")
process.l1NtupleProducer.csctfStatusSource    = cms.InputTag("none")
process.l1NtupleProducer.csctfDTStubsSource   = cms.InputTag("none")

# remove RPC TT because it seg faults
process.ValL1Emulator.remove(process.valRpcTechTrigDigis)

# re-wire emulators to include upstream modules in L1Extra
#process.valGctDigis.inputLabel = 'valRctDigis'
process.valGmtDigis.DTCandidates = cms.InputTag('valDttfDigis')
process.valGmtDigis.RPCbCandidates = cms.InputTag('valRpcTriggerDigis', 'RPCb')
process.valGmtDigis.RPCfCandidates = cms.InputTag('valRpcTriggerDigis', 'RPCf')

# path
process.p = cms.Path(
    
    process.RawToDigi

    +process.l1NtupleProducer

    +process.l1GtTriggerMenuLite
    +process.l1MenuTreeProducer
    
    +process.l1extraParticles
    +process.l1ExtraTreeProducer

    +process.L1EmulatorTree

)


# global tag etc
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# GCT re-emulation
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("L1GctJetFinderParamsRcd"),
             tag = cms.string("L1GctJetFinderParams_GCTPhysics_2012_04_27_JetSeedThresh5GeV_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
             )
)

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("L1GctChannelMaskRcd"),
             tag = cms.string("L1GctChannelMask_AllEnergySumsMaskedFromHF_jetCentresToEta3Allowed_mc"),
             connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_L1T")
             )
)

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

# job options
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",
                             fileNames = readFiles,
                             secondaryFileNames = secFiles
                             )

print "Done"
