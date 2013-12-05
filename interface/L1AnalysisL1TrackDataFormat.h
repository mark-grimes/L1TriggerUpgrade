#ifndef L1Analysis_L1AnalysisL1TrackDataFormat_h
#define L1Analysis_L1AnalysisL1TrackDataFormat_h

#include <vector>

namespace L1Analysis
{
	/** @brief Data format to hold L1 track trigger information in ntuples.
	 *
	 * Can't say I'm keen on this way of working - if a change needs to be made I've got to
	 * change it in about 5 places. Just following what seems to be the norm for trigger studies.
	 * I basically just copied what was done for L1AnalysisL1ExtraUpgradeDataFormat.
	 *
	 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
	 * @date 03/Dec/2013
	 */
	struct L1AnalysisL1TrackDataFormat
	{
		L1AnalysisL1TrackDataFormat() { Reset(); }
		~L1AnalysisL1TrackDataFormat() {}

		void Reset()
		{
			nTrackEG=0;
			trackEG_eT.clear();
			trackEG_eta.clear();
			trackEG_phi.clear();
			trackEG_bunchCrossing.clear();
			trackEG_trackIsolation.clear();
			trackEG_emObjectIndex.clear();
			trackEG_emType.clear();
			trackEG_debug_EmParticleEt.clear();
			trackEG_debug_EmParticleEta.clear();
			trackEG_debug_EmParticlePhi.clear();
		}

		unsigned int nTrackEG;
		std::vector<float> trackEG_eT;
		std::vector<float> trackEG_eta;
		std::vector<float> trackEG_phi;
		std::vector<int> trackEG_bunchCrossing;
		std::vector<float> trackEG_trackIsolation;
		std::vector<int> trackEG_emObjectIndex;
		std::vector<int> trackEG_emType;

		// I'll temporarily add some extra information while debugging
		// to make sure that I'm indexing the L1EmParticleCollection
		// properly (record this information here and then check it
		// against what I get when I look up the EmParticle later).
		std::vector<float> trackEG_debug_EmParticleEt;
		std::vector<float> trackEG_debug_EmParticleEta;
		std::vector<float> trackEG_debug_EmParticlePhi;
	};

} // end of namespace L1Analysis
#endif
