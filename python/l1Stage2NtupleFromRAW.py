import FWCore.ParameterSet.Config as cms

# import master L1 Ntuple config
from UserCode.L1TriggerUpgrade.l1CurrentNtupleFromRAW import *

print "Adding Stage 2 algorithms"

# add full upgrade objects
process.load('UserCode.L1TriggerUpgrade.L1Upgrade_EGTau_cff')
process.load('UserCode.L1TriggerUpgrade.L1Upgrade_Jet_cff')
process.load('UserCode.L1TriggerUpgrade.l1UpgradeMuonIsolator_cfi')

# load the upgrade tree
process.load('UserCode.L1TriggerUpgrade.l1ExtraUpgradeTreeProducer_cfi')

# add upgrade algorithms to the path
process.p += process.SLHCCaloTrigger
process.p += process.l1UpgradeMuonIsolator
process.p += process.l1ExtraUpgradeTreeProducer    # upgrade candidates

print "Done"
