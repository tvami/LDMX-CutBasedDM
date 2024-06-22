import os
import sys
fileName = sys.argv[1]


from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('cutbased')

import LDMX.Ecal.ecal_hardcoded_conditions
from LDMX.Ecal import EcalGeometry
from LDMX.Hcal import hcal
from LDMX.Ecal import vetos

p.maxEvents = -1
p.run = 2

# p.inputFiles = [f'events.root']
#p.inputFiles  = ["/Users/tav/Documents/1Research/LDMX/CutBasedDM/mc_v14-8gev-8.0GeV-1e-ecal_photonuclear_run24515_t1703519903_trigSkim.root"]
p.inputFiles  = [fileName]
#p.inputFiles = [f'signalIn.root']
#p.inputFiles = [f'ecalPnIn.root']
p.histogramFile = fileName[:-5] + "_histo_v3_nonfid.root"
#p.termLogLevel = 0

CutBasedAna = ldmxcfg.Analyzer.from_file('/home/vamitamas/CutBasedDM/ldmx-sw/CutBasedDM.cxx')
CutBasedAna.fiducial = False

ecalVeto = vetos.EcalVetoProcessor()
ecalVeto.collection_name= 'EcalVetoNew'

CutBasedAna.trigger_name = "TriggerSums20Layers"
CutBasedAna.trigger_pass = "signal"
hcalVeto   =hcal.HcalVetoProcessor('hcalVeto')
p.sequence = [ecalVeto,hcalVeto,CutBasedAna]



