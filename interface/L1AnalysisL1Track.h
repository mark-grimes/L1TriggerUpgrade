#ifndef L1Analysis_L1AnalysisL1Track_h
#define L1Analysis_L1AnalysisL1Track_h

#include "L1AnalysisL1TrackDataFormat.h"
#include "DataFormats/L1TrackTrigger/interface/L1TrackEmParticleFwd.h"

namespace L1Analysis
{
	/** @brief Class that fills a L1AnalysisL1TrackDataFormat structure given the edm information
	 *
	 * Form copied from L1AnalysisL1ExtraUpgrade. Just trying to match whatever is the norm for trigger
	 * analysis work practices.
	 *
	 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
	 * @date 03/Dec/2013
	 */
	class L1AnalysisL1Track
	{
	public:
		L1AnalysisL1Track();
		~L1AnalysisL1Track();
		void Reset() { l1TrackData_.Reset(); }

	    /** @brief Adds the electron collection with L1 track information to the tree. */
	    void SetTrackEG( const edm::Handle<l1extra::L1TrackEmParticleCollection> trackEGs );

	    L1AnalysisL1TrackDataFormat* getData() { return &l1TrackData_; }
	  private :
	    L1AnalysisL1TrackDataFormat l1TrackData_;
	};

} // end of namespace L1Analysis

#endif
