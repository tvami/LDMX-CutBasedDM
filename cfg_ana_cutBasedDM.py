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
from LDMX.Recon.trackDeDxMassEstimator import recoilTrackMassEstimator
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
#p.histogramFile = "/sdf/group/ldmx/users/tamasvami/CutBasedDM/"  + str(fileName.split('/')[-2]) + "/" + str(fileName.split('/')[-1][:-5]) + "_histo_v16.root"
p.histogramFile = "/home/vamitamas/CutBasedDM/"  + str(fileName.split('/')[-2]) + "/" + str(fileName.split('/')[-1][:-5]) + "_histo_v18.root"
# p.histogramFile  = str(fileName.split('/')[-1][:-5]) + "_histo_v16.root"
print("Histogram output = ", p.histogramFile)


# Acceptance filter
signal = False
accaptance = []
if "signal" in fileName: 
    signal = True
    accaptance = [RecoilFiducialityProcessor("RecoilFiducialityProcessor")]
    print("We are running on signal\n----------------------------------")

if "tagetEN_8gev" in fileName: 
    sp_pass_temp = 'genie'
elif signal :
    sp_pass_temp = ''
else:
    sp_pass_temp = 'sim'

#Ecal vetos
ecalVeto = vetos.EcalVetoProcessor()
ecalVeto.collection_name= 'EcalVetoNew'
ecalVeto.sp_pass_name = sp_pass_temp
ecalVeto.recoil_from_tracking = True
#ecalVeto.track_pass_name = 'cutbased'
ecalVeto.run_lin_reg = False

# HCAL veto
# import LDMX.Hcal.HcalGeometry
# import LDMX.Hcal.hcal_hardcoded_conditions
# import LDMX.Hcal.digi as hcal_digi
# hcal_digi_reco = hcal_digi.HcalSimpleDigiAndRecProducer()
hcalVeto   = hcal.HcalVetoProcessor('hcalVeto')

# Trigger
trigger = TriggerProcessor('Trigger', 8000.)

cutBasedAna = ldmxcfg.Analyzer.from_file('CutBasedDM.cxx')
#cutBasedAna.fiducial_analysis = False
cutBasedAna.fiducial_analysis = True
cutBasedAna.trigger_name = "Trigger"
cutBasedAna.trigger_pass = "cutbased"
cutBasedAna.signal = signal
cutBasedAna.ignore_fiducial_analysis = False
cutBasedAna.sp_pass_name = sp_pass_temp
cutBasedAna.recoil_track_collection = 'RecoilTracksClean'
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
    full_tracking_sequence.greedy_solver_tagger,
    full_tracking_sequence.greedy_solver_recoil,
    # full_tracking_sequence.GSF_tagger,
    # full_tracking_sequence.GSF_recoil
]

p.termLogLevel = 10
# p.logger.custom("TruthSeedProcessor", level = 10)
# p.logger.custom(ecalVeto, level = 10)
# p.logger.custom(full_tracking_sequence.greedy_solver_recoil, level = -1)

p.logFrequency = 1000

for seq in track_sqs_1:
    seq.sp_pass_name = sp_pass_temp

for seq in track_sqs_2:
    seq.input_pass_name = "cutbased"

p.sequence.extend(track_sqs_1)
p.sequence.extend(track_sqs_2)
p.sequence.extend(accaptance)

p.sequence.extend([ecalVeto, hcalVeto, trigger])
p.sequence.extend([cutBasedAna])

