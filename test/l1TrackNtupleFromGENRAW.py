import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1Stage2NtupleFromRAW import *

from UserCode.L1TriggerUpgrade.MCSetup import *
mcSetup(process, False, True)

# job options
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

readFiles.extend( [
    #'/store/mc/Summer12/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/UpgradeL1TDR-PU35_POSTLS161_V12-v2/00000/168A4B24-C53C-E211-A370-003048C68AA6.root'
    '/store/mc/UpgFall13d/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94955470-F338-E311-B942-003048FF86CA.root'
    ] )

# Remove this from the sequence because it will be replaced by an SLHC version in a minute
del process.trackerNumberingGeometry
del process.trackerGeometry  # will be replaced by "trackerSLHCGeometry"

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')

process.L1TrackPhotons = cms.EDProducer("L1TrackEmParticleProducer",
	L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),
        label = cms.string("NonIsolated")
)

process.l1TrackTreeProducer = cms.EDAnalyzer("L1TrackTreeProducer",
	trackEGLabel = cms.untracked.InputTag("L1TrackPhotons","NonIsolated")
)

# Add code to do the L1 tracks and create track EM objects
process.p += process.BeamSpotFromSim
process.p += process.L1Tracks
process.p += process.L1TrackPhotons
process.p += process.l1TrackTreeProducer
