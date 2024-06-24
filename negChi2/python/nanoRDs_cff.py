from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *   #why do we need all these nanoaod functions?
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

# Bs chain collection
from questions.negChi2.BsToDsPhiKKPiMu_cff import *

# from PhysiscsTools.NanoAOD
nanoSequence = cms.Sequence(nanoMetadata + globalTables)

def nanoAOD_customizeBsToDsPhiKKPiMu(process):
    process.nanoBsToDsPhiKKPiMuSequence = cms.Sequence( BsToDsPhiKKPiMuSequence ) #+ CountBsToDsPhiKKPiMu)#+HighMassLowMassFlagsTables )
    return process

