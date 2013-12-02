import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1CurrentNtupleFromRAW import *

print "Going to run on data"

#process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(10000)

process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1000)
            )

readFiles.extend( [
    'file:/gpfs_phys/storm/cms/data/Run2012C/ZeroBias/RAW/v1/000/198/588/1ADBB979-1ACA-E111-AA51-003048D3C932.root'
        ] )

