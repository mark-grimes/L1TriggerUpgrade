import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1Stage2NtupleFromRAW import *

print "Going to run on data"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
                )

readFiles.extend( [
    'file:/gpfs_phys/storm/cms/data/Run2012C/ZeroBias/RAW/v1/000/198/588/1ADBB979-1ACA-E111-AA51-003048D3C932.root'
    ] )

#process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(10000)

process.output = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                           'keep *_SLHCL1ExtraParticles_*_*',
                                           'keep *_PUsubAK5CaloJetProducer_*_*',
                                           'keep *_calibTowerJetProducer_*_*',
                                           'keep *_ak5CaloJets_*_RECO',
                                           'drop *_ak5CaloJets_rho_RECO',
                                           'keep *_ak5CaloJets_rho_RCTofflineTEST',
                                           'keep *_ak5PFJets_*_*',
                                           'keep *_l1extraParticles_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_L1TowerJetProducer_*_*', 
                                           'keep *_L1TowerJetPUcalc_*_*',
                                           'keep *_L1TowerJetFilter2D_*_*',
                                           'keep *_L1TowerFwdJetFilter2D_*_*',
                                           'keep *_ak5JetID_*_*',
                                           'keep *_L1TowerJetPU*_*_*' 
                                           ),
    fileName = cms.untracked.string('output.root')
    )

process.e = cms.EndPath(process.output)
