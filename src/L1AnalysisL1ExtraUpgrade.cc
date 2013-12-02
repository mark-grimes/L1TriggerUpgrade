#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisL1ExtraUpgrade.h"

L1Analysis::L1AnalysisL1ExtraUpgrade::L1AnalysisL1ExtraUpgrade()
{
}

L1Analysis::L1AnalysisL1ExtraUpgrade::~L1AnalysisL1ExtraUpgrade()
{

}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetEG(const edm::Handle<l1extra::L1EmParticleCollection> eg, unsigned maxL1Extra)
{
  for(l1extra::L1EmParticleCollection::const_iterator it=eg->begin(); it!=eg->end() && l1extra_.nEG<maxL1Extra; it++){  
      l1extra_.egEt .push_back(it->et());
      l1extra_.egEta.push_back(it->eta());
      l1extra_.egPhi.push_back(it->phi());
      l1extra_.egBx .push_back(it->bx());
      l1extra_.nEG++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetIsoEG(const edm::Handle<l1extra::L1EmParticleCollection> isoEG, unsigned maxL1Extra)
{
  for(l1extra::L1EmParticleCollection::const_iterator it=isoEG->begin(); it!=isoEG->end() && l1extra_.nIsoEG<maxL1Extra; it++){
      l1extra_.isoEGEt .push_back(it->et());
      l1extra_.isoEGEta.push_back(it->eta());
      l1extra_.isoEGPhi.push_back(it->phi());
      l1extra_.isoEGBx .push_back(it->bx());
      l1extra_.nIsoEG++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetTau(const edm::Handle<l1extra::L1JetParticleCollection> tau, unsigned maxL1Extra)
{
      //std::cout << "Filling L1 Extra tauJets" << std::endl;      
   for(l1extra::L1JetParticleCollection::const_iterator it=tau->begin(); it!=tau->end() && l1extra_.nTau<maxL1Extra; it++){
      
     // printf("L1tauJet (et,eta,phi,bx,) (%f,%f,%f,%d)\n",it->et(),it->eta(),it->phi(),it->bx() );
      l1extra_.tauEt .push_back(it->et());
      l1extra_.tauEta.push_back(it->eta());
      l1extra_.tauPhi.push_back(it->phi());
      l1extra_.tauBx .push_back(it->bx());
      l1extra_.nTau++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetIsoTau(const edm::Handle<l1extra::L1JetParticleCollection> isoTau, unsigned maxL1Extra)
{
      //std::cout << "Filling L1 Extra tauJets" << std::endl;      
   for(l1extra::L1JetParticleCollection::const_iterator it=isoTau->begin(); it!=isoTau->end() && l1extra_.nIsoTau<maxL1Extra; it++){
      
     // printf("L1tauJet (et,eta,phi,bx,) (%f,%f,%f,%d)\n",it->et(),it->eta(),it->phi(),it->bx() );
      l1extra_.isoTauEt .push_back(it->et());
      l1extra_.isoTauEta.push_back(it->eta());
      l1extra_.isoTauPhi.push_back(it->phi());
      l1extra_.isoTauBx .push_back(it->bx());
      l1extra_.nIsoTau++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetJet(const edm::Handle<l1extra::L1JetParticleCollection> jet, unsigned maxL1Extra)
{
//      std::cout << "Filling L1 Extra cenJets" << maxL1Extra << " " << cenJet->size() << std::endl;      
 
  for(l1extra::L1JetParticleCollection::const_iterator it=jet->begin(); it!=jet->end() && l1extra_.nJets<maxL1Extra; it++){
      //printf("L1CenJet (et,eta,phi,bx,) (%f,%f,%f,%d) \n",it->et(),it->eta(),it->phi(),it->bx() );
//      std::cout << "L1 CenJets et,eta,phi,bx = " << it->et() << ", " << it->eta() <<", " <<it->phi() <<", " << it->bx() << std::endl;
      l1extra_.jetEt .push_back(it->et());
      l1extra_.jetEta.push_back(it->eta());
      l1extra_.jetPhi.push_back(it->phi());
      l1extra_.jetBx .push_back(it->bx());
      l1extra_.nJets++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetFwdJet(const edm::Handle<l1extra::L1JetParticleCollection> fwdJet, unsigned maxL1Extra)
{
      //std::cout << "Filling L1 Extra fwdJets" << std::endl;      
   for(l1extra::L1JetParticleCollection::const_iterator it=fwdJet->begin(); it!=fwdJet->end() && l1extra_.nFwdJets<maxL1Extra; it++){ 
      //printf("L1fwdJet (et,eta,phi,bx,) (%f,%f,%f,%d)\n",it->et(),it->eta(),it->phi(),it->bx() );
      l1extra_.fwdJetEt .push_back(it->et());
      l1extra_.fwdJetEta.push_back(it->eta());
      l1extra_.fwdJetPhi.push_back(it->phi());
      l1extra_.fwdJetBx .push_back(it->bx());
      l1extra_.nFwdJets++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetMuon(const edm::Handle<l1extra::L1MuonParticleCollection> muon, unsigned maxL1Extra)
{
  for(l1extra::L1MuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && l1extra_.nMuons<maxL1Extra; it++){
      
      l1extra_.muonEt .push_back( it->et());
      l1extra_.muonEta.push_back(it->eta());
      l1extra_.muonPhi.push_back(it->phi());
      l1extra_.muonChg.push_back(it->charge());
      l1extra_.muonIso.push_back(it->isIsolated());
      l1extra_.muonMip.push_back(it->isMip());
      l1extra_.muonFwd.push_back(it->isForward());
      l1extra_.muonRPC.push_back(it->isRPC());
      l1extra_.muonBx .push_back(it->bx());
      l1extra_.muonQuality .push_back(it->gmtMuonCand().quality());
		
//		std::cout << "gmtmuon cand: pt " << it->gmtMuonCand().ptValue() 
//					<< "; ptExtra " << it->et() 
//					<< "; qual " << it->gmtMuonCand().quality() 
//					<< std::endl;
      l1extra_.nMuons++;
    }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetMet(const edm::Handle<l1extra::L1EtMissParticleCollection> mets)
{
  for(l1extra::L1EtMissParticleCollection::const_iterator it=mets->begin(); it!=mets->end(); it++) {
    l1extra_.et.    push_back( it->etTotal() ); 
    l1extra_.met.   push_back( it->et() );
    l1extra_.metPhi.push_back( it->phi() );
    l1extra_.metBx. push_back( it->bx() );
    l1extra_.nMet++;
  }
}

void L1Analysis::L1AnalysisL1ExtraUpgrade::SetMht(const edm::Handle<l1extra::L1EtMissParticleCollection> mhts)
{
  for(l1extra::L1EtMissParticleCollection::const_iterator it=mhts->begin(); it!=mhts->end(); it++) {
    l1extra_.ht.    push_back( it->etTotal() );
    l1extra_.mht.   push_back( it->et() );
    l1extra_.mhtPhi.push_back( it->phi() );
    l1extra_.mhtBx. push_back( it->bx() );
    l1extra_.nMht++;
  }
}
