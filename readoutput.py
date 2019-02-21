from DataFormats.FWLite import Events, Handle
from ROOT import TH1F, TCanvas, TLorentzVector, TFile
import sys
import math
#import CommonTools.CandAlgos
#import PhysicsTools.HepMCCandAlgos
#import FWCore.ParameterSet.Config as cms
from array import array
import FWCore.ParameterSet.Config as cms

inputfiles = ['/afs/cern.ch/user/r/rselvati/work/private/ECAL_sim/CMSSW_10_3_1/SingleElectronPt10_pythia8_cfi_py_GEN_SIM_DIGI.root']

events = Events(inputfiles)

simEcalDigis_h = Handle("EBDigiCollection")
simEcalDigis_l = ('simEcalDigis')
#simEcalDigis_h = Handle("vector<PCaloHit>")
#simEcalDigis_l = ('g4SimHits')

nEvents =1

for event in events:

    event.getByLabel(simEcalDigis_l,"ebDigis",simEcalDigis_h)
    #event.getByLabel(simEcalDigis_l,"EcalHitsEB",simEcalDigis_h)

    simEcalDigis = simEcalDigis_h.product()


    print "Event number ", nEvents
    print

    print "size: ", simEcalDigis.size()

    for ecalDigi in simEcalDigis:
        print ciao


    print
    print "//////////////////"
    print

    nEvents += 1 

