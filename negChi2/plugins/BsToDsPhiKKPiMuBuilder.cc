#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>    

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

// for vertex fitting (both global and sequential)
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h" 
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h" 
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for the tracks!
#include "FWCore/Framework/interface/MakerMacros.h"

#include "KinVtxFitter.h"   // 
#include "helper.h"         // helper functions
#include "TLorentzVector.h" // use this instead 
#include "TVector3.h" // for boost vector

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

// B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"

using namespace std;

class BsToDsPhiKKPiMuBuilder : public edm::global::EDProducer<> {
public:

  //define collections which dont exist by default  
  //constructor
  explicit BsToDsPhiKKPiMuBuilder(const edm::ParameterSet&);
  //destructor
  ~BsToDsPhiKKPiMuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

  reco::TransientTrack getTransientTrack(const reco::Track track) const {    
      reco::TransientTrack transientTrack(track, paramField);
      return transientTrack;
    }


private:
 
  // muon selection
  const StringCutObjectSelector<pat::Muon> muSelection_;

  //define tokens to access data later
  const edm::InputTag muonTag;
  const edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;

  const edm::InputTag triggerBitTag;
  const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  const edm::InputTag triggerObjectTag;
  const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

  const edm::InputTag triggerPrescaleTag;
  const edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  const edm::InputTag vertexSrcTag;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

  //the maximal dR you allow between pat muon and trigger muon candidate
  const string trgFilterLabel_;
  const string hlt_7_4_p0_;
  const string hlt_7_4_p1_;
  const string hlt_7_4_p2_;
  const string hlt_7_4_p3_;
  const string hlt_7_4_p4_;
  const double maxdR_;

 
  //Bfield
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //cuts 
  const StringCutObjectSelector<pat::PackedCandidate> hadSelection_; // cut on hadrons

  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double maxdxyHadPv_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double maxBsMass_;
  const double piMass_;
  const double piMassSigma_;
  const double kMass_;
  const double kMassSigma_;
  const double phiMass_;
  const bool   constrainPhiMass_;
  const double minPhiVtxProb_;
  const double dsMass_;
  const bool   constrainDsMass_;
  const double minDsVtxProb_;
  const double dsStarMass_;
  const double muMass_;
  const double muMassSigma_;
  const double bsMass_;
  const double isoCone_;
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  const edm::InputTag srcTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> src_;

  //for the muons

  const edm::InputTag trgMuonTag;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuons_;

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  // lost tracks for isolation
  const edm::InputTag tracksLostTag;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksLost_;
 
};

//define the constructor
BsToDsPhiKKPiMuBuilder::BsToDsPhiKKPiMuBuilder(const edm::ParameterSet& iConfig):

    //for the muons
    muSelection_(iConfig.getParameter<std::string>("muSelection")),
    muonTag(iConfig.getParameter<edm::InputTag>("muonCollection")),
    muonSrc_(consumes<std::vector<pat::Muon>>(muonTag)),
    //for trigger info
    triggerBitTag(iConfig.getParameter<edm::InputTag>("trgResultsCollection")),
    triggerBits_(consumes<edm::TriggerResults>(triggerBitTag)),
  
    triggerObjectTag(iConfig.getParameter<edm::InputTag>("trgObjectsCollection")),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectTag)),
  
    triggerPrescaleTag(iConfig.getParameter<edm::InputTag>("trgPrescaleCollection")),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(triggerPrescaleTag)),
    //vertex info
    vertexSrcTag(iConfig.getParameter<edm::InputTag>("vtxCollection")),
    vertexSrc_(consumes<reco::VertexCollection>(vertexSrcTag)),
  
    //parameters
  
    trgFilterLabel_(iConfig.getParameter<string>("trgFilterLabel")),
    hlt_7_4_p0_(iConfig.getParameter<string>("hlt_7_4_p0")),
    hlt_7_4_p1_(iConfig.getParameter<string>("hlt_7_4_p1")),
    hlt_7_4_p2_(iConfig.getParameter<string>("hlt_7_4_p2")),
    hlt_7_4_p3_(iConfig.getParameter<string>("hlt_7_4_p3")),
    hlt_7_4_p4_(iConfig.getParameter<string>("hlt_7_4_p4")),
    maxdR_(iConfig.getParameter<double>("maxdR_matching")),


    // f.e. hadSelection_ = cfg.getPatameter...
    hadSelection_(iConfig.getParameter<std::string>("hadSelection")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    maxdxyHadPv_(iConfig.getParameter<double>("maxdxyHadPv")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    maxBsMass_(iConfig.getParameter<double>("maxBsMass")),

    piMass_(iConfig.getParameter<double>("piMass")),
    piMassSigma_(iConfig.getParameter<double>("piMassSigma")),
    kMass_(iConfig.getParameter<double>("kMass")),
    kMassSigma_(iConfig.getParameter<double>("kMassSigma")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    constrainPhiMass_(iConfig.getParameter<bool>("constrainPhiMass")),
    minPhiVtxProb_(iConfig.getParameter<double>("minPhiVtxProb")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    constrainDsMass_(iConfig.getParameter<bool>("constrainDsMass")),
    minDsVtxProb_(iConfig.getParameter<double>("minDsVtxProb")),
    dsStarMass_(iConfig.getParameter<double>("dsStarMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    muMassSigma_(iConfig.getParameter<double>("muMassSigma")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
    isoCone_(iConfig.getParameter<double>("isoCone")),

    srcTag(iConfig.getParameter<edm::InputTag>("pfCand")),
    src_(consumes<pat::PackedCandidateCollection>(srcTag)), 

    trgMuonTag(iConfig.getParameter<edm::InputTag>("muCand")),
    trgMuons_(consumes<pat::MuonCollection>(trgMuonTag)), 

    primaryVtxTag(iConfig.getParameter<edm::InputTag>("pvCand")),
    primaryVtx_(consumes<reco::VertexCollection>(primaryVtxTag)),

    tracksLostTag(iConfig.getParameter<edm::InputTag>("lostTracks")),
    tracksLost_(consumes<pat::PackedCandidateCollection>(tracksLostTag)){

       // output collection
       produces<pat::CompositeCandidateCollection>("bs");
    }

//check const keywords 

// this starts the event loop
void BsToDsPhiKKPiMuBuilder::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //Define handles
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
 
  //is used to store information about the results of trigger decisions in a given event
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
 
  //is used to store information about the results of trigger decisions in a given event
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
 
  // for every event, get the name list of the triggers
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
 
 
  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
 
  //pat muons contain more info than reco muons, f.e. trigger info! 
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);
 
  //input
  edm::Handle<pat::PackedCandidateCollection> pcand;
  iEvent.getByToken(src_, pcand);
 
 
  edm::Handle<pat::MuonCollection> trgMuons;
  iEvent.getByToken(trgMuons_,trgMuons);
 
  edm::Handle<reco::VertexCollection> primaryVtx;
  iEvent.getByToken(primaryVtx_,primaryVtx);

  edm::Handle<pat::PackedCandidateCollection> tracksLost;
  iEvent.getByToken(tracksLost_,tracksLost);

  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> bsCandidates(new pat::CompositeCandidateCollection());

  // Getting the indices of the HLT paths
  unsigned int index_7_4_p0      = names.triggerIndex(hlt_7_4_p0_);
  unsigned int index_7_4_p1      = names.triggerIndex(hlt_7_4_p1_);
  unsigned int index_7_4_p2      = names.triggerIndex(hlt_7_4_p2_);
  unsigned int index_7_4_p3      = names.triggerIndex(hlt_7_4_p3_);
  unsigned int index_7_4_p4      = names.triggerIndex(hlt_7_4_p4_);
 
  // full lumi
  //unsigned int index_8_3      = names.triggerIndex("HLT_Mu8_IP3"); 
  //unsigned int index_8_5      = names.triggerIndex("HLT_Mu8_IP5");
  //unsigned int index_8_6      = names.triggerIndex("HLT_Mu8_IP6");
  //unsigned int index_8p5_3p5  = names.triggerIndex("HLT_Mu8p5_IP3p5");
  //unsigned int index_9_4      = names.triggerIndex("HLT_Mu9_IP4");
  //unsigned int index_9_5      = names.triggerIndex("HLT_Mu9_IP5");
  //unsigned int index_10p5_3p5 = names.triggerIndex("HLT_Mu10p5_IP3p5");
  //unsigned int index_12_6     = names.triggerIndex("HLT_Mu12_IP6");
 
  //std::cout << (index_7_4_p0 < triggerBits->size())  << "and" <<      (triggerBits->accept(index_7_4_p0)) << std::endl;
  //std::cout << (index_7_4_p1 < triggerBits->size())  << "and" <<     (triggerBits->accept(index_7_4_p1)) << std::endl;
  //std::cout << (index_7_4_p2 < triggerBits->size())  << "and" <<    (triggerBits->accept(index_7_4_p2)) << std::endl;
  //std::cout << (index_7_4_p3 < triggerBits->size())  << "and" <<     (triggerBits->accept(index_7_4_p3)) << std::endl;
  //std::cout << (index_7_4_p4 < triggerBits->size())  << "and" <<    (triggerBits->accept(index_7_4_p4)) << std::endl;
 
  //default is false  
  bool pass_7_4_p0_path      = false;
  bool pass_7_4_p1_path      = false;
  bool pass_7_4_p2_path      = false;
  bool pass_7_4_p3_path      = false;
  bool pass_7_4_p4_path      = false;
 
  // check first if the index is valid, i.e. if it is not out of range (maximum is givrn by triggerBits->size())
  // and if so, check if the trigger has fired with accept() 
 
  pass_7_4_p0_path      = ((index_7_4_p0 < triggerBits->size())   && (triggerBits->accept(index_7_4_p0)));
  pass_7_4_p1_path      = ((index_7_4_p1 < triggerBits->size())   && (triggerBits->accept(index_7_4_p1)));
  pass_7_4_p2_path      = ((index_7_4_p2 < triggerBits->size())   && (triggerBits->accept(index_7_4_p2)));
  pass_7_4_p3_path      = ((index_7_4_p3 < triggerBits->size())   && (triggerBits->accept(index_7_4_p3)));
  pass_7_4_p4_path      = ((index_7_4_p4 < triggerBits->size())   && (triggerBits->accept(index_7_4_p4)));
 
 
  //define vector out of bools
  std::vector<bool> trgMuonFrom_7_4_flag;
 
  //only continue when we the event passes the HLT_Mu7_IP4
  if (pass_7_4_p0_path || pass_7_4_p1_path || pass_7_4_p2_path || pass_7_4_p3_path || pass_7_4_p4_path){
 
  //std::cout << "i pass the pathes" << std::endl;
 
  //std::cout<<"found trigger!" << std::endl;
  // define vectors of ints of length muons->size(), all values set to 0
  std::vector<int> isTriggerMuon(muons->size(), 0);
 
  //////////////////////////////////////////////////////////////////////
  // Make sure that you can find a pat  muon matching the HLT object  //
  //////////////////////////////////////////////////////////////////////
 
  // now loop over all pat::muons
  for (unsigned int muIdx=0; muIdx<muons->size(); ++muIdx){
    
    edm::Ptr<pat::Muon> muPtr(muons, muIdx); 

    if (!muSelection_(*muPtr)) continue;
 
 
    //check if the pat muon is matched to some trigger object (by using the function triggerObjectMatchByPath()
 
    // initialize start values
    float drMuonTrgObj = -1.;
    int muonIdx        = -1;

    //////////////////////////////////////////////////////////////////////
    // Find the (pat muon, HLT)  pair with the smallest dR and save it  //
    //////////////////////////////////////////////////////////////////////
    int iTrg = 0;
    for(unsigned int objIdx=0; objIdx < triggerObjects->size(); ++objIdx){
 
      iTrg++;
 
      pat::TriggerObjectStandAlone trgObj = (*triggerObjects)[objIdx];
      //unpack trigger labels
 
      trgObj.unpackFilterLabels(iEvent, *triggerBits);
 
      std::vector<std::string> filterLabels = trgObj.filterLabels();
 
      //check if the triggermuon was actually firing the trigger
      if(!trgObj.hasFilterLabel(trgFilterLabel_)) continue;
 
      //save the dR between the triggering muon and the pat muon 
      float dr = reco::deltaR(trgObj, *muPtr);
 
 
      if( ( (dr < drMuonTrgObj) || (drMuonTrgObj == -1))  && (dr < maxdR_))
      {
              //modify dR and the reco and triggering muon index
              drMuonTrgObj = dr;
              muonIdx      = muIdx;
      }
 
    } // closing loop over trg muons
 
    // Build candidate only if we have trg matched muon 
    if(muonIdx != -1)
     {
 
        //save also the tracks
        const reco::TransientTrack ttrackTrgMuon(*(muPtr->bestTrack()), paramField);
        if (!ttrackTrgMuon.isValid()) continue;
 
        std::vector<float> dzMuPV;
      
        //Fix the primary vertex to be the one closest to the trg Muon in dz
      
        float dummy = 1.0;
        int goldenIdx = -1;
        reco::Vertex pv;
        for(size_t vtxIdx = 0; vtxIdx < primaryVtx->size(); ++vtxIdx){
          edm::Ptr<reco::Vertex> vtxPtr(primaryVtx, vtxIdx);
          float minDz = fabs(muPtr->bestTrack()->dz(vtxPtr->position())); 
          if(minDz < dummy){
            dummy = minDz;
            goldenIdx = vtxIdx;
          }
      
        }
        if (goldenIdx<0) continue;
        if (goldenIdx >= 0){
        pv = primaryVtx->at(goldenIdx);
        }
      
        //////////////////////////////////////////////////
        // Loop over k1 and select the good tracks      //
        //////////////////////////////////////////////////
      
        for(size_t k1Idx = 0; k1Idx < pcand->size() + tracksLost->size() ; ++k1Idx) {
      
          //define a pointer to the kaon at position k1Idx
          edm::Ptr<pat::PackedCandidate> k1Ptr;
          if (k1Idx < pcand->size()) k1Ptr = edm::Ptr<pat::PackedCandidate>(pcand, k1Idx); //normal tracks
          else k1Ptr = edm::Ptr<pat::PackedCandidate>(tracksLost, k1Idx - pcand->size());  //lost tracks
      
          if (!hadSelection_(*k1Ptr)) continue; 
      
          float muonK1dR = reco::deltaR(*k1Ptr,*muPtr);
      
          bool k1Sel = (( muonK1dR < maxdRHadMuon_ ) && 
          (reco::deltaR(*k1Ptr, *muPtr) > mindRHadMuon_) && 
          (abs(k1Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_)) &&
          (abs(k1Ptr->bestTrack()->dxy(pv.position()) < maxdxyHadPv_ )) ;
      
          if (!k1Sel) continue;
          //////////////////////////////////////////////////
          // Loop over k2 and select the good tracks      //
          //////////////////////////////////////////////////
      
          for(size_t k2Idx = k1Idx + 1; k2Idx < pcand->size()+ tracksLost->size() ; ++k2Idx) {
      
          //make sure k2 is not k1
          if (k2Idx == k1Idx) continue;
      
          edm::Ptr<pat::PackedCandidate> k2Ptr;
          if (k2Idx < pcand->size()) k2Ptr = edm::Ptr<pat::PackedCandidate>(pcand, k2Idx); //normal tracks
          else k2Ptr = edm::Ptr<pat::PackedCandidate>(tracksLost, k2Idx - pcand->size());  //lost tracks
        
          if(!hadSelection_(*k2Ptr)) continue;
      
          float muonK2dR = reco::deltaR(*k2Ptr,*muPtr);
          bool k2Sel = (( muonK2dR < maxdRHadMuon_ ) && 
          (reco::deltaR(*k2Ptr, *muPtr) > mindRHadMuon_) && 
          (abs(k2Ptr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_)) &&
          (abs(k2Ptr->bestTrack()->dxy(pv.position()) < maxdxyHadPv_ )) ;
      
          if (!k2Sel) continue;
      
          //////////////////////////////////////////////////
          // Loop over pi and select the good tracks      //
          //////////////////////////////////////////////////
      
          for(size_t piIdx = 0; piIdx < pcand->size()+ tracksLost->size() ; ++piIdx) {
      
            //make sure the pion is none of the kaons:
            if((piIdx == k1Idx) || (piIdx == k2Idx)) continue;
      
            edm::Ptr<pat::PackedCandidate> piPtr;
            if (piIdx < pcand->size()) piPtr = edm::Ptr<pat::PackedCandidate>(pcand, piIdx); //normal tracks
            else piPtr = edm::Ptr<pat::PackedCandidate>(tracksLost, piIdx - pcand->size());  //lost tracks
      
            // if this pion does not pass the selection, jump to the next!
            if(!hadSelection_(*piPtr)) continue;
      
            float muonPidR = reco::deltaR(*piPtr,*muPtr);
      
            bool piSel = ((muonPidR < maxdRHadMuon_) && 
            (reco::deltaR(*piPtr, *muPtr) > mindRHadMuon_) &&
            (abs(piPtr->bestTrack()->dz(pv.position()) - muPtr->bestTrack()->dz(pv.position())) < maxdzDiffHadMuon_)) &&
            (abs(piPtr->bestTrack()->dxy(pv.position()) < maxdxyHadPv_ )) ;
      
            if (!piSel) continue;
      
            //////////////////////////////////////////////////
            // Build Phi resonance                          //
            //////////////////////////////////////////////////
      
      
            //define a composite candidate pair 
            pat::CompositeCandidate kk;
      
            //PF canidates are always assigned with the pi mass, we have to force the kaon mass
            math::PtEtaPhiMLorentzVector k1P4(k1Ptr->pt(), k1Ptr->eta(), k1Ptr->phi(), kMass_);
            math::PtEtaPhiMLorentzVector k2P4(k2Ptr->pt(), k2Ptr->eta(), k2Ptr->phi(), kMass_);
            math::PtEtaPhiMLorentzVector piP4(piPtr->pt(), piPtr->eta(), piPtr->phi(), piMass_); //just to be sure lets also force the pi mass
       
            kk.setP4(k1P4 + k2P4);
      
            //only continue when they build a phi resonance, allow 15MeV:
            if (fabs(kk.mass() - phiMass_) > phiMassAllowance_) continue;     
            kk.setCharge(k1Ptr->charge() + k2Ptr->charge());
      
            //////////////////////////////////////////////////
            // Build Ds resonance                           //
            // ///////////////////////////////////////////////
      
            pat::CompositeCandidate phiPi;
            phiPi.setP4(kk.p4() + piP4); 
      
            //only continue when they build a ds resonance, allow 50MeV:
            if (fabs(phiPi.mass() - dsMass_) > dsMassAllowance_) continue;
      
            //std::cout << "we passed the ds resonance" << std::endl; 
            phiPi.setCharge(kk.charge() + piPtr->charge());
            //std::cout << "found ds resonance" << std::endl;
      
            //////////////////////////////////////////////////
            // Build Bs resonance                           //
            //////////////////////////////////////////////////
      
            pat::CompositeCandidate dsMu;
            dsMu.setP4(phiPi.p4() + muPtr->p4()); 
            dsMu.setCharge(phiPi.charge() + muPtr->charge()); //sanity check:shoould be 0
      
            if(dsMu.mass() > maxBsMass_) continue;
      
            //build bs with collinear approximation
            pat::CompositeCandidate bs;
      
            bs.setP4(dsMu.p4() * bsMass_ / dsMu.mass()); //the bs_mass will thus be fixed at 536688 (peak in the histo)
            bs.setCharge(dsMu.charge());
      
            std::cout << "\n Found mu k k pi candidate! Lets do some fits! "  << std::endl; 

            ////////////////////////////////////////////////
            // Now we do a proper fit                     //
            ////////////////////////////////////////////////
      
            //define a factory
            KinematicParticleFactoryFromTransientTrack pFactory;
      
            //define the vector for the particles to be fitted
            std::vector<RefCountedKinematicParticle> phiToFit;
            std::vector<RefCountedKinematicParticle> dsToFit;
            std::vector<RefCountedKinematicParticle> bsToFit;
      
            // add masses
            ParticleMass piMass  = piMass_;
            ParticleMass kMass   = kMass_;
            ParticleMass phiMass = phiMass_;
            ParticleMass dsMass  = dsMass_;
            ParticleMass muMass  = muMass_;
      
            
            float ndf = 0.;
            float chi = 0.;
            float kMassSigma = kMassSigma_ ;
            float piMassSigma = piMassSigma_;
            float muMassSigma = muMassSigma_;
      
            // fix the tracks to have pos def covariance matrix, taken from:
            // https://github.com/CMSKStarMuMu/miniB0KstarMuMu/blob/master/miniKstarMuMu/plugins/miniKstarMuMu.cc#L1611-L1678
      
            const reco::Track k1Track = *k1Ptr->bestTrack();
            const reco::Track k2Track = *k2Ptr->bestTrack();
            const reco::Track piTrack = *piPtr->bestTrack();
            const reco::Track muTrack = *muPtr->bestTrack();
      
            reco::Track k1TrackCorr   = correctCovMat(&k1Track, 1e-8);
            reco::Track k2TrackCorr   = correctCovMat(&k2Track, 1e-8);
            reco::Track piTrackCorr   = correctCovMat(&piTrack, 1e-8);
            reco::Track muTrackCorr   = correctCovMat(&muTrack, 1e-8);
      
            reco::TransientTrack ttK1 = ttBuilder->build(k1TrackCorr);
            reco::TransientTrack ttK2 = ttBuilder->build(k2TrackCorr);
            reco::TransientTrack ttPi = ttBuilder->build(piTrackCorr);
            reco::TransientTrack ttMu = ttBuilder->build(muTrackCorr);
      
            phiToFit.push_back( pFactory.particle(ttK1, kMass,  chi, ndf, kMassSigma ));
            phiToFit.push_back( pFactory.particle(ttK2, kMass,  chi, ndf, kMassSigma ));
            dsToFit.push_back(  pFactory.particle(ttPi, piMass, chi, ndf, piMassSigma));
            bsToFit.push_back(  pFactory.particle(ttMu, muMass, chi, ndf, muMassSigma));
      
            //////////////////////////////////
            // First we do a direct fit     //
            // of the Ds into three tracks: //
            // Ds -> K K Pi                 //
            //////////////////////////////////

            KinVtxFitter easyFitter(
            {ttK1, ttK2, getTransientTrack(piTrackCorr)},
            {kMass, kMass, piMass},
            {kMassSigma,kMassSigma,piMassSigma}
            );
            if(!easyFitter.success()) continue;
            std::cout << "\n---------- DIRECT FIT ----------\n" << std::endl;
            std::cout << "Ds fit chi2 is: " << easyFitter.chi2() << std::endl;
            std::cout << "Ds fit vtx prob is: " << easyFitter.prob() << std::endl;
            
            /////////////////////////////////
            // Now we do a multistep fit:  //
            // Start with Phi -> K K fit   // 
            /////////////////////////////////
     
            std::cout << "\n---------- MULTISTEP FIT ----------\n" << std::endl;
            //perform fit with vertexFit function (in helper.h) 
            RefCountedKinematicTree phiTree  = vertexFit(phiToFit, phiMass, constrainPhiMass_);
      
            if (!phiTree->isValid() || phiTree->isEmpty() || !phiTree->isConsistent()) continue; //check if fit result is valid
      
            //access the fitted resonance and vertex 
            phiTree->movePointerToTheTop();
            RefCountedKinematicParticle phiParticle = phiTree->currentParticle();
      
            //get vtx chi2 and ndof
      
            auto phiVtx = phiTree->currentDecayVertex();
            if (!phiVtx->vertexIsValid() || !phiParticle->currentState().isValid() ) continue; //check if fit result is valid
      
            float phiVtxChi2    = phiVtx->chiSquared();
            std::cout << "Phi fit chi2 is: " << phiVtxChi2 << std::endl;
            if (phiVtxChi2 < 0) continue;
      
            float phiVtxNDof    = phiVtx->degreesOfFreedom();
            float phiVtxProb    = ChiSquaredProbability(phiVtxChi2, phiVtxNDof); 
            std::cout << "Phi fit vtx prob is: " << phiVtxProb << std::endl;
            if (phiVtxProb < 0.01) continue;
            
      
            // access refitted children
            phiTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle phiDau1 = phiTree->currentParticle();
            phiTree->movePointerToTheNextChild();
            RefCountedKinematicParticle phiDau2 = phiTree->currentParticle();
      
            // add the phi to the list of particles (pi) to fit the ds 
            dsToFit.push_back(phiParticle);
     
            /////////////////////////////// 
            // Ds -> Phi Pi fit          //
            /////////////////////////////// 
      
            //perform fit with vertexFit function (in helper.h) 
            RefCountedKinematicTree dsTree = vertexFit(dsToFit, dsMass, constrainDsMass_);
            if (!dsTree->isValid() || dsTree->isEmpty() ) continue; //check if fit result is valid
      
            // access the fitted resonance and the refitted children
            dsTree->movePointerToTheTop();
            RefCountedKinematicParticle dsParticle = dsTree->currentParticle();
      
            // get vtx chi2 and ndof
            RefCountedKinematicVertex dsVtx = dsTree->currentDecayVertex(); //compare to the access via AlgebraicVector7
            if (!dsVtx->vertexIsValid()) continue; //check if fit result is valid
      
            float dsVtxChi2    = dsVtx->chiSquared();
            std::cout << "Ds fit chi2 is: " << dsVtxChi2 << std::endl;
            if (dsVtxChi2 < 0) continue;
            float dsVtxNDof    = dsVtx->degreesOfFreedom();
            float dsVtxProb    = ChiSquaredProbability(dsVtxChi2, dsVtxNDof); 
            std::cout << "Ds fit vtx prob is: " << dsVtxProb << std::endl;
            if (dsVtxProb < 0.01) continue;
      
            dsTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle dsDau1 = dsTree->currentParticle();
            dsTree->movePointerToTheNextChild();
            RefCountedKinematicParticle dsDau2 = dsTree->currentParticle();
      
            // add the ds to the list of particles to fit the bs
            bsToFit.push_back(dsParticle);
            
            /////////////////////////////// 
            // Bs -> Ds Mu (nu) fit      //
            /////////////////////////////// 

            //perform fit with vertexFit function (in helper.h) 
            RefCountedKinematicTree bsTree = vertexFit(bsToFit, bsMass_, false); // no constraint for bs because missing momentum
            if (!bsTree->isValid() || bsTree->isEmpty() ) continue; //check if fit result is valid 
      
            // access the fitted resonance and the refitted children
            bsTree->movePointerToTheTop();
            RefCountedKinematicParticle bsParticle = bsTree->currentParticle();
      
            // get vtx chi2 and ndof
            RefCountedKinematicVertex bsVtx = bsTree->currentDecayVertex();
            if (!bsVtx->vertexIsValid()) continue; //check if fit result is valid
      
            float bsVtxChi2    = bsVtx->chiSquared();
            std::cout << "Bs fit chi2 is: " << bsVtxChi2 << std::endl;
            if (bsVtxChi2 < 0) continue;
            float bsVtxNDof    = bsVtx->degreesOfFreedom();
            float bsVtxProb    = ChiSquaredProbability(bsVtxChi2, bsVtxNDof); 
            std::cout << "Bs fit vtx prob is: " << bsVtxProb << std::endl;
      
            bsTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle bsDau1 = bsTree->currentParticle();
            bsTree->movePointerToTheNextChild();
            RefCountedKinematicParticle bsDau2 = bsTree->currentParticle();
      
      
            //////////////////////////////////// end of global fitter /////////////////////////////////////
            
            bsCandidates->emplace_back(bs);
        } //closing pi loop
      } //closing k2 loop
    } //closing k1 loop
  } //closing trg muon condition
} // closing muon loop
} //closing HLT condition 
  iEvent.put(std::move(bsCandidates), "bs");
}//closing event loop

DEFINE_FWK_MODULE(BsToDsPhiKKPiMuBuilder);
