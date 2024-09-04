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

p.maxEvents = -1
p.run = 2

# p.inputFiles = [f'events.root']
#p.inputFiles  = ["/Users/tav/Documents/1Research/LDMX/CutBasedDM/mc_v14-8gev-8.0GeV-1e-ecal_photonuclear_run24515_t1703519903_trigSkim.root"]
print(fileIn)
p.inputFiles  = fileIn
#p.inputFiles = [f'signalIn.root']
#p.inputFiles = [f'ecalPnIn.root']
#p.histogramFile = fileName[:-5] + "_histo_v4_nonfid.root"
p.histogramFile = "/sdf/home/t/tamasvami/CutBasedDM/"  + str(fileName.split('/')[-2]) + "/" + str(fileName.split('/')[-1][:-5]) + "_histo_v7.root"
print("Output: " + str(fileName.split('/')[-2]) + "/" + str(fileName.split('/')[-1][:-5]) + "_histo_v7.root")
#p.termLogLevel = 0

CutBasedAna = ldmxcfg.Analyzer.from_file('/sdf/home/t/tamasvami/CutBasedDM/ldmx-sw/CutBasedDM.cxx')
#CutBasedAna.fiducial = False
CutBasedAna.fiducial = True

ecalVeto = vetos.EcalVetoProcessor()
ecalVeto.collection_name= 'EcalVetoNew'

CutBasedAna.trigger_name = "Trigger"
CutBasedAna.trigger_pass = "cutbased"
hcalVeto   =hcal.HcalVetoProcessor('hcalVeto')
p.sequence = [ecalVeto,hcalVeto,TriggerProcessor('Trigger', 8000.),CutBasedAna]


