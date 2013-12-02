// -*- C++ -*-
//
// Package:    L1UpgradeMuonIsolator
// Class:      L1UpgradeMuonIsolator
// 
/**\class L1UpgradeMuonIsolator L1UpgradeMuonIsolator.cc UserCode/L1UpgradeMuonIsolator/src/L1UpgradeMuonIsolator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Brooke
//         Created:  Tue Jan 15 16:47:46 GMT 2013
// $Id: L1UpgradeMuonIsolator.cc,v 1.1 2013/01/19 17:26:09 jbrooke Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Event Setup
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"

// input collections
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

// output collections
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

// intermediate stuff
#include "L1Trigger/UCT2015/src/L1GObject.h"


//
// class declaration
//

class L1UpgradeMuonIsolator : public edm::EDProducer {
public:
  explicit L1UpgradeMuonIsolator(const edm::ParameterSet&);
  ~L1UpgradeMuonIsolator();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  int deltaPhi18(int phi1,
		 int phi2);
  int deltaGctPhi(const L1CaloRegion& r1,
		  const L1CaloRegion& r2);
  L1CaloRegionCollection::const_iterator matchObjectToRegion(const L1CaloGeometry* geom,
							     const edm::Handle<L1CaloRegionCollection> &newRegions,
							     const float &eta, const float &phi);
  L1GObject buildJetAtIndex(const edm::Handle<L1CaloRegionCollection> &newRegions,
			    const L1CaloRegionCollection::const_iterator &newRegion,
			    const float &regionLSB_,
			    float threshold);
  
  
  // ----------member data ---------------------------
  static double muonMassGeV_;
  
  edm::InputTag muonSource_ ;
  edm::InputTag regionSource_ ;
  
  double rgnThreshold_;
  double regionLSB_;
  double muFixedIso_;  // require 3x3 Et below this value
  double muRelIso_;    // placeholder for relative isolation

};

//
// constants, enums and typedefs
//

double L1UpgradeMuonIsolator::muonMassGeV_ = 0.105658369 ; // PDG06


//
// static data member definitions
//

//
// constructors and destructor
//
L1UpgradeMuonIsolator::L1UpgradeMuonIsolator(const edm::ParameterSet& iConfig) :
     muonSource_( iConfig.getParameter< edm::InputTag >("muonSource" ) ),
     regionSource_( iConfig.getParameter< edm::InputTag >("regionSource") ),
     rgnThreshold_( iConfig.getParameter<double>("rgnThreshold") ),
     regionLSB_( iConfig.getParameter<double>("regionLSB") ),
     muFixedIso_( iConfig.getParameter<double>("muFixedIso") ),
     muRelIso_( iConfig.getParameter<double>("muRelIso") )
{

  //register your products
  produces< l1extra::L1MuonParticleCollection >( "" ) ;
  
}


L1UpgradeMuonIsolator::~L1UpgradeMuonIsolator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1UpgradeMuonIsolator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace l1extra ;
   using namespace std ;

   // get muon/calo scales and geometry
   ESHandle< L1MuTriggerScales > muScales ;
   iSetup.get< L1MuTriggerScalesRcd >().get( muScales ) ;
   
   ESHandle< L1MuTriggerPtScale > muPtScale ;
   iSetup.get< L1MuTriggerPtScaleRcd >().get( muPtScale ) ;

   ESHandle< L1CaloGeometry > caloGeomESH ;
   iSetup.get< L1CaloGeometryRecord >().get( caloGeomESH ) ;
   const L1CaloGeometry* caloGeom = &( *caloGeomESH ) ;
   

   // get muons   
   Handle< L1MuGMTReadoutCollection > hwMuCollection ;
   iEvent.getByLabel( muonSource_, hwMuCollection ) ;
   
   vector< L1MuGMTExtendedCand > hwMuCands ;
   
   // get calo info
   edm::Handle<L1CaloRegionCollection> regions;
   iEvent.getByLabel(regionSource_, regions);

   // create output collection
   auto_ptr< L1MuonParticleCollection > muColl(
     new L1MuonParticleCollection );

   if( !hwMuCollection.isValid() ) {
     LogWarning("NoInput") << "Input muon collection not found. Not going to produce output" << std::endl;
   }
   else {

     hwMuCands = hwMuCollection->getRecord().getGMTCands() ;

     vector< L1MuGMTExtendedCand >::const_iterator muItr = hwMuCands.begin() ;
     vector< L1MuGMTExtendedCand >::const_iterator muEnd = hwMuCands.end() ;
     for( int i = 0 ; muItr != muEnd ; ++muItr, ++i ) {
	 if( !muItr->empty() ) {
	     // keep x and y components non-zero and protect against roundoff.
	     double pt =
	       muPtScale->getPtScale()->getLowEdge( muItr->ptIndex() ) + 1.e-6 ;
	     
	     // 	    cout << "L1Extra pt " << pt << endl ;
	     
	     double eta =
	       muScales->getGMTEtaScale()->getCenter( muItr->etaIndex() ) ;
	     
	     double phi =
	       muScales->getPhiScale()->getLowEdge( muItr->phiIndex() ) ;
	     
	     math::PtEtaPhiMLorentzVector p4( pt,
					      eta,
					      phi,
					      muonMassGeV_ ) ;

	     L1MuonParticle muonParticle( muItr->charge(),
					  p4,
					  *muItr,
					  muItr->bx() );
	     
	     // calculate isolation
	     bool isolated = false;
	     
	     if ( !regions.isValid() ) {
	       LogWarning("NoInput") << "Input region collection not found. Output muons will be non-isolated." << std::endl;
	     }
	     else {

	       L1CaloRegionCollection::const_iterator matchedRegion = 
		 matchObjectToRegion(caloGeom,
				     regions,
				     eta,
				     phi);
	       L1GObject seededJet = buildJetAtIndex(regions,
						     matchedRegion,
						     regionLSB_,
						     rgnThreshold_);
	       
	       isolated = (muFixedIso_ > 0.) && (seededJet.pt() < muFixedIso_);
	       
	       muonParticle.setIsolated(isolated);
	     }

	     // push to collection
	     //if (isolated) isoMuColl->push_back( muonParticle );
	     muColl->push_back( muonParticle );
	     
	   }
       }

   }

   OrphanHandle< L1MuonParticleCollection > muHandle =
     iEvent.put( muColl );

}

// ------------ method called once each job just before starting event loop  ------------
void 
L1UpgradeMuonIsolator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1UpgradeMuonIsolator::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
L1UpgradeMuonIsolator::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
L1UpgradeMuonIsolator::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
L1UpgradeMuonIsolator::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
L1UpgradeMuonIsolator::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1UpgradeMuonIsolator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int L1UpgradeMuonIsolator::deltaPhi18(int phi1, int phi2)
{   
    // Compute the difference in phi between two towers, wrapping at phi = 18
    int difference = phi1 - phi2;
    if (std::abs(phi1 - phi2) == 17) {
        difference = -difference/std::abs(difference);
    }
    return difference;
}   


int L1UpgradeMuonIsolator::deltaGctPhi(const L1CaloRegion& r1, const L1CaloRegion& r2)
{
    return deltaPhi18(r1.gctPhi(), r2.gctPhi());
}   


L1CaloRegionCollection::const_iterator L1UpgradeMuonIsolator::matchObjectToRegion(const L1CaloGeometry* geom,
										  const edm::Handle<L1CaloRegionCollection> &regions,
										  const float &eta,
										  const float &phi)
{

  // find closest region to eta value...
  L1CaloRegionCollection::const_iterator matchedRegion = regions->begin();
  float dRMin = 999.9;
  for(L1CaloRegionCollection::const_iterator region = regions->begin();
      region != regions->end(); region++)
    {
      
      double regionEta = geom->etaBinCenter( region->id() ) ;
      double regionPhi = geom->emJetPhiBinCenter( region->id() ) ;
      float dEta = (eta - regionEta);
      float dPhi = (phi - regionPhi);
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      if (dR > dRMin) continue;
      matchedRegion = region;
      dRMin = dR;
    }
  
  return matchedRegion;
  
}


L1GObject L1UpgradeMuonIsolator::buildJetAtIndex(const edm::Handle<L1CaloRegionCollection> &newRegions,
						 const L1CaloRegionCollection::const_iterator &newRegion, 
						 const float &regionLSB_, 
						 float threshold) {

  double regionET = newRegion->et() * regionLSB_;
  if (regionET < threshold) regionET = 0.0;
  
  double neighborN_et = 0;
  double neighborS_et = 0;
  double neighborE_et = 0;
  double neighborW_et = 0;
  double neighborNE_et = 0;
  double neighborSW_et = 0;
  double neighborNW_et = 0;
  double neighborSE_et = 0;
  unsigned int nNeighbors = 0;
  
  //std::cout << "Looking for seed @ " << newRegion->gctPhi() << " " << newRegion->gctEta() << std::endl;
  
  for(L1CaloRegionCollection::const_iterator neighbor = newRegions->begin();
      neighbor != newRegions->end(); neighbor++)
    {
      
      double neighborET = neighbor->et() * regionLSB_;
      if (neighborET < threshold) neighborET = 0.0;
      
      if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	 (newRegion->gctEta()    ) == neighbor->gctEta()) {
	neighborN_et = neighborET;
	nNeighbors++;
	//debug("N", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta()    ) == neighbor->gctEta()) {
	neighborS_et = neighborET;
	nNeighbors++;
	//debug("S", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborE_et = neighborET;
	nNeighbors++;
	//debug("E", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 0 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborW_et = neighborET;
	nNeighbors++;
	//debug("W", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborNE_et = neighborET;
	nNeighbors++;
	//debug("NE", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborSW_et = neighborET;
	nNeighbors++;
	//debug("SW", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == 1 &&
	      (newRegion->gctEta() - 1) == neighbor->gctEta()) {
	neighborNW_et = neighborET;
	nNeighbors++;
	//debug("NW", *newRegion, *neighbor);
	continue;
      }
      else if(deltaGctPhi(*newRegion, *neighbor) == -1 &&
	      (newRegion->gctEta() + 1) == neighbor->gctEta()) {
	neighborSE_et = neighborET;
	nNeighbors++;
	//debug("SE", *newRegion, *neighbor);
	continue;
      }
    }
  
  //   ---- this should probably be removed for constructing "isolation" jets
  
  //    if(regionET > neighborN_et &&
  //            regionET > neighborNW_et &&
  //            regionET > neighborW_et &&
  //            regionET > neighborSW_et &&
  //            regionET >= neighborNE_et &&
  //            regionET >= neighborE_et &&
  //            regionET >= neighborSE_et &&
  //            regionET >= neighborS_et) 
  //    {

  unsigned int jetET = regionET +
    neighborN_et + neighborS_et + neighborE_et + neighborW_et +
    neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;
  
  // Temporarily use the region granularity -- we will try to improve as above when code is debugged
  int jetPhi = newRegion->gctPhi();
  int jetEta = newRegion->gctEta();
  
  bool neighborCheck = (nNeighbors == 8);
  // On the eta edge we only expect 5 neighbors
  if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
    neighborCheck = true;
  if (!neighborCheck) {
    std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
    assert(false);
  }
  
  L1GObject jet(jetET, jetEta, jetPhi, "MuonSeededJet");
  jet.associatedRegionEt_ = regionET;
  return jet;
  //    }
  
  std::cout << "it's all gone a little bit wrong" << std::endl;
  assert(false);
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(L1UpgradeMuonIsolator);
