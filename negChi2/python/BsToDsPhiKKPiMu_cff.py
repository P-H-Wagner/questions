import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars

BsToDsPhiKKPiMu = cms.EDProducer(
    'BsToDsPhiKKPiMuBuilder',
    muonCollection        = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
    trgResultsCollection  = cms.InputTag("TriggerResults", "", "HLT"),
    trgObjectsCollection  = cms.InputTag("slimmedPatTrigger"),
    trgPrescaleCollection = cms.InputTag("patTrigger"),
    vtxCollection         = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #cuts and selections
    trgFilterLabel = cms.string("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q"),
    muSelection    = cms.string(' &&  '.join([
    'pt > 7.0',
    'eta > -1.5',
    'eta < 1.5',
    'isPFMuon',
    'isGlobalMuon'])), # pre-selection of hadrons (k1,k2 and pion)
    hlt_7_4_p0   = cms.string("HLT_Mu7_IP4_part0_v2"), # trigger menue
    hlt_7_4_p1   = cms.string("HLT_Mu7_IP4_part1_v2"), # "
    hlt_7_4_p2   = cms.string("HLT_Mu7_IP4_part2_v2"), # "
    hlt_7_4_p3   = cms.string("HLT_Mu7_IP4_part3_v2"), # "
    hlt_7_4_p4   = cms.string("HLT_Mu7_IP4_part4_v2"), # "
    maxdR_matching = cms.double(0.05), #muon trg object matching  

    pfCand = cms.InputTag('packedPFCandidates'),
    muCand = cms.InputTag('muonTrgSelector', 'trgMuons'),
    pvCand = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tracks = cms.InputTag('packedPFCandidates'), # for isolation
    lostTracks = cms.InputTag("lostTracks"),     # consider also lost pf cands!
    hadSelection = cms.string(' &&  '.join([
'abs(pdgId) != 11',
'abs(pdgId) != 13',
'charge != 0',
'pt > 1.0', 
'eta > -2.4',
'eta < 2.4',
'hasTrackDetails'])), # pre-selection of hadrons (k1,k2 and pion)
maxdRHadMuon     = cms.double( 1.2),       # max dR between hadron and muon
mindRHadMuon     = cms.double( 0.005),     # min dR "
maxdzDiffHadMuon = cms.double( 0.5),   # difference in dz between muon/pv and had/pv
maxdxyHadPv      = cms.double( 0.6),
phiMassAllowance = cms.double( 0.015),  # allow 15 MeV when collecting candidates for phi 
dsMassAllowance  = cms.double( 0.06),   # allow 150 MeV when collecting candidates for ds
maxBsMass        = cms.double( 8.0 ),  
piMass           = cms.double( 0.13957039),      # pi mass
piMassSigma      = cms.double( 0.00000018),      # pi mass
kMass            = cms.double( 0.493677),         # kaon mass
kMassSigma       = cms.double( 0.000016),         # kaon mass
phiMass          = cms.double( 1.019461),       # phi mass
constrainPhiMass = cms.bool(False),    # constrain phi mass in the vtx fit?
minPhiVtxProb    = cms.double(0.01),
dsMass           = cms.double( 1.96834),         # ds mass
constrainDsMass  = cms.bool(False),     # constrain Ds mass in the vtx fit?
minDsVtxProb     = cms.double(0.01),
dsStarMass       = cms.double( 2.112204),    # ds star mass
muMass           = cms.double( 0.105658),        # mu ma
muMassSigma      = cms.double( 0.0000000023),        # mu mass
bsMass           = cms.double( 5.36688),         # bs mass
isoCone          = cms.double( 0.5)             # cut on dR for the mu isolation cone
)

print( " ========> Parameters used:")
print(BsToDsPhiKKPiMu.dumpPython)

BsToDsPhiKKPiMuSequence = cms.Sequence(BsToDsPhiKKPiMu)

