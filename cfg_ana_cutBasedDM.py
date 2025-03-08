import os
import sys
#fileName = sys.argv[1]
fileIn = sys.argv[1:]
fileName =  " ".join(sys.argv[1:])

from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('cutbased')

import LDMX.Ecal.ecal_hardcoded_conditions
from LDMX.Ecal import EcalGeometry
from LDMX.Hcal import hcal
from LDMX.Ecal import vetos
from LDMX.Recon.simpleTrigger import TriggerProcessor
from LDMX.Recon.fiducialFlag import RecoilFiducialityProcessor

p.maxEvents = -1
p.run = 1
p.skipCorruptedInputFiles = True

# p.inputFiles = [f'events.root']
#p.inputFiles  = ["/Users/tav/Documents/1Research/LDMX/CutBasedDM/mc_v14-8gev-8.0GeV-1e-ecal_photonuclear_run24515_t1703519903_trigSkim.root"]
p.inputFiles  = fileIn
print("Input files = ", p.inputFiles)
#p.inputFiles = [f'signalIn.root']
#p.inputFiles = [f'ecalPnIn.root']
#p.histogramFile = fileName[:-5] + "_histo_v4_nonfid.root"
p.histogramFile = "/sdf/group/ldmx/users/tamasvami/CutBasedDM/"  + str(fileName.split('/')[-2]) + "/" + str(fileName.split('/')[-1][:-5]) + "_histo_v16.root"
# p.histogramFile  = str(fileName.split('/')[-1][:-5]) + "_histo_v16.root"
print("Histogram output = ", p.histogramFile)
#p.termLogLevel = 0


# Acceptance filter
signal = False
accaptance = []
if "signal" in fileName: 
    signal = True
    accaptance = [RecoilFiducialityProcessor("RecoilFiducialityProcessor")]
    print("We are running on signal\n----------------------------------")

#Ecal vetos
ecalVeto = vetos.EcalVetoProcessor()
ecalVeto.collection_name= 'EcalVetoNew'
ecalVeto.recoil_from_tracking = True

# HCAL veto
hcalVeto   = hcal.HcalVetoProcessor('hcalVeto')

# Trigger
trigger = TriggerProcessor('Trigger', 8000.)

CutBasedAna = ldmxcfg.Analyzer.from_file('CutBasedDM.cxx')
#CutBasedAna.fiducial_analysis = False
CutBasedAna.fiducial_analysis = True
CutBasedAna.trigger_name = "Trigger"
CutBasedAna.trigger_pass = "cutbased"
CutBasedAna.signal = signal
CutBasedAna.ignore_fiducial_analysis = False

p.sequence = []

from LDMX.Tracking import full_tracking_sequence



track_sqs_1 = [
    full_tracking_sequence.digi_tagger,
    full_tracking_sequence.digi_recoil,
    full_tracking_sequence.truth_tracking,
    full_tracking_sequence.seeder_tagger,
    full_tracking_sequence.seeder_recoil
]
track_sqs_2 = [
    full_tracking_sequence.tracking_tagger,
    full_tracking_sequence.tracking_recoil,
    # full_tracking_sequence.greedy_solver_tagger,
    # full_tracking_sequence.greedy_solver_recoil,
    # full_tracking_sequence.GSF_tagger,
    # full_tracking_sequence.GSF_recoil
]

p.logger.custom("TruthSeedProcessor", level = 10)
# p.logger.custom(ecalVeto, level = 10)
p.logFrequency = 1000

for seq in track_sqs_2:
    seq.input_pass_name = "cutbased"

p.sequence.extend(track_sqs_1)
p.sequence.extend(track_sqs_2)
p.sequence.extend(accaptance)

p.sequence.extend([ecalVeto, hcalVeto, trigger, CutBasedAna])


