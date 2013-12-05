#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisL1Track.h"


class L1TrackTreeProducer : public edm::EDAnalyzer
{
public:
	explicit L1TrackTreeProducer( const edm::ParameterSet& config );
	~L1TrackTreeProducer();

private:
	virtual void beginJob( void );
	virtual void analyze( const edm::Event& event, const edm::EventSetup& eventSetup );
	virtual void endJob();

	std::auto_ptr<L1Analysis::L1AnalysisL1Track> pL1Track_;
	L1Analysis::L1AnalysisL1TrackDataFormat* pL1TrackDataFormat_;

	// tree
	TTree* pTree_;

	// EDM input tags
	edm::InputTag trackEGLabel_;
};

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackTreeProducer);


L1TrackTreeProducer::L1TrackTreeProducer( const edm::ParameterSet& config )
	: trackEGLabel_( config.getUntrackedParameter("trackEGLabel",edm::InputTag("")) )
{
	pL1Track_.reset( new L1Analysis::L1AnalysisL1Track );
	pL1TrackDataFormat_=pL1Track_->getData();

	edm::Service<TFileService> fileService;
	pTree_=fileService->make<TTree>("L1TrackTree", "L1TrackTree");
	pTree_->Branch("L1Track", "L1Analysis::L1AnalysisL1TrackDataFormat", &pL1TrackDataFormat_, 32000, 3);
}

L1TrackTreeProducer::~L1TrackTreeProducer()
{

}

void L1TrackTreeProducer::beginJob( void )
{

}

void L1TrackTreeProducer::analyze( const edm::Event& event, const edm::EventSetup& eventSetup )
{
	pL1Track_->Reset();

	edm::Handle<l1extra::L1TrackEmParticleCollection> hTrackEGs;

	event.getByLabel( trackEGLabel_, hTrackEGs );

	if( hTrackEGs.isValid() )
	{
		pL1Track_->SetTrackEG( hTrackEGs );
	}
	else
	{
		edm::LogWarning( "MissingProduct" ) << "L1ExtraUpgrade Track EG objects not found (" << trackEGLabel_ << "). Branch will not be filled" << std::endl;
	}

	pTree_->Fill();
}

void L1TrackTreeProducer::endJob()
{

}
