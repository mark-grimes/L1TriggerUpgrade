#include "UserCode/L1TriggerUpgrade/interface/L1AnalysisL1Track.h"

// Includes for things that were only forward declared
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticle.h"

L1Analysis::L1AnalysisL1Track::L1AnalysisL1Track()
{
}

L1Analysis::L1AnalysisL1Track::~L1AnalysisL1Track()
{

}

void L1Analysis::L1AnalysisL1Track::SetTrackEG( const edm::Handle<l1extra::L1TrackEmParticleCollection> trackEGs )
{
	for( const l1extra::L1TrackEmParticle& trackEMParticle : *trackEGs )
	{
		l1TrackData_.trackEG_eT.push_back( trackEMParticle.et() );
		l1TrackData_.trackEG_eta.push_back( trackEMParticle.eta() );
		l1TrackData_.trackEG_phi.push_back( trackEMParticle.phi() );
		l1TrackData_.trackEG_bunchCrossing.push_back( trackEMParticle.getEGRef()->bx() );
		l1TrackData_.trackEG_trackIsolation.push_back( trackEMParticle.getTrkIsol() );
		l1TrackData_.trackEG_emObjectIndex.push_back( trackEMParticle.getEGRef().index() );
		l1TrackData_.nTrackEG++;

		// This is extra debug information I'm storing about the referenced
		// L1EmParticle, so that I can check I'm looking it up correctly
		// using trackEG_emObjectIndex.
		l1TrackData_.trackEG_debug_EmParticleEt.push_back( trackEMParticle.getEGRef()->et() );
		l1TrackData_.trackEG_debug_EmParticleEta.push_back( trackEMParticle.getEGRef()->eta() );
		l1TrackData_.trackEG_debug_EmParticlePhi.push_back( trackEMParticle.getEGRef()->phi() );
	}
}
