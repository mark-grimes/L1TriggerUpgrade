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

// data formats
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

class L1CaloCorrector : public edm::EDProducer {
public:
  L1CaloCorrector(const edm::ParameterSet& pset);
  virtual ~L1CaloCorrector();
  void produce(edm::Event& iEvent, const edm::EventSetup& iEventSetup);
private:
  
  edm::InputTag hcalSource_;
  //    edm::InputTag ecalSource_;
  
  // HCAL corrections
  TH1D* reassignment_likelihoods[32][500];//each iphi bin, each possible value of adc count

};

L1CaloCorrector::L1CaloCorrector(const edm::ParameterSet& pset) :
  hcalSource_(pset.getParameter<edm::InputTag>("hcalSource"))
				//  ecalTPSrc_(pset.getParameter<edm::InputTag>("ecalSource")),
{

  edm::FileInPath fip = pset.getParameter<edm::FileInPath>("hcalFile");
  std::string iFile = fip.fullPath();

  edm::LogWarning("InputFile") << "Taking input from file : " << iFile << std::endl;

  TFile down_shift_file(iFile.c_str());

  assert(down_shift_file.IsOpen());

  for(int ieta=1; ieta <=32; ieta++ ){
    TGraph* tmp_graph=(TGraph*) down_shift_file.Get(Form("graph_%d",ieta));
    assert(tmp_graph != NULL);//This is used for interpolation for empty slices below
    const int n_graph_points = tmp_graph->GetN();
    double last_x,last_y;
    tmp_graph->GetPoint(n_graph_points-1,last_x,last_y);//so I know how far to trust the interp
    
    TH2D* transfer_function = (TH2D*)down_shift_file.Get(Form("transfer_function_%d",ieta));
    assert(transfer_function != NULL);
    assert(transfer_function->GetNbinsX()==500+1);//I started it at -1 on purpose.  Should have an extra bin.
    for(int adc_value=0; adc_value < 500; adc_value++){
      const int adc_bin_number=adc_value+2;//the first bin is centered on -1.
      reassignment_likelihoods[ieta-1][adc_value] = 
	(TH1D*)transfer_function->ProjectionY(Form("ieta%d_py%d",ieta,adc_value),adc_bin_number,adc_bin_number);//vertical slice
      assert(reassignment_likelihoods[ieta-1][adc_value] != NULL);
      reassignment_likelihoods[ieta-1][adc_value]->SetDirectory(0);//protect from closing files
      if(reassignment_likelihoods[ieta-1][adc_value]->Integral()<0.5){//empty
	if( adc_value < last_x ){//I need an interpolation and can trust the graph.
	  const float interpolation = tmp_graph->Eval(adc_value);//Just linear, no extrapolation. Chill.
	  reassignment_likelihoods[ieta-1][adc_value]->Fill(interpolation);
	}else{//no idea, just keep the current value
	  reassignment_likelihoods[ieta-1][adc_value]->Fill(adc_value);
	}
      }//make sure there's no empty histograms
    }//loop over adc_counts
  }//loop over iphi
  down_shift_file.Close();

  produces<HcalTrigPrimDigiCollection>("");

}


L1CaloCorrector::~L1CaloCorrector() {

  for (int i=0; i<32; ++i) {
    for (int j=0; j<500; ++j) {
      if (reassignment_likelihoods[i][j] != 0) {
	delete reassignment_likelihoods[i][j];
      }
    }
  }

}


void L1CaloCorrector::produce(edm::Event& iEvent, const edm::EventSetup& iEventSetup) {

  // get input
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPin;
  iEvent.getByLabel(hcalSource_, hcalTPin);

  // create output collections
  std::auto_ptr<HcalTrigPrimDigiCollection> hcalTPout(new HcalTrigPrimDigiCollection);
  
  if (hcalTPin.isValid()) {
	
    hcalTPout->reserve(hcalTPin->size());
    
    // apply correction
    edm::SortedCollection<HcalTriggerPrimitiveDigi>::const_iterator old_digi;
    for (old_digi = hcalTPin->begin(); old_digi != hcalTPin->end(); old_digi++) {
      
      HcalTriggerPrimitiveDigi new_digi(*old_digi);//start by total copy
      
      //-------------------------------------------
      //scaling part
      const int presample = old_digi->presamples();//which presample do we rescale
      
      int old_et = old_digi->sample(presample).compressedEt();
      const bool old_finegrain = old_digi->sample(presample).fineGrain();
      const int old_slb = old_digi->sample(presample).slb();
      const int old_slbchan = old_digi->sample(presample).slbChan();
      
      if(old_et >= 500) old_et = 500-1;//for binning of reassignment_likelihoods
      
      const int ieta = new_digi.id().ieta();
      const float tmp_et = reassignment_likelihoods[abs(ieta)-1][old_et]->GetRandom();//Boosh!
      const int new_et = (int)floor(tmp_et+0.5);//make sure only the bin centers are generated
      
      HcalTriggerPrimitiveSample new_sample  //use the new et
	= HcalTriggerPrimitiveSample(new_et, old_finegrain, old_slb, old_slbchan);
      new_digi.setSample(presample,new_sample);
      
      if(new_et >= 0) {//can't do this with HcalTriggerPrimitiveDigi::SOI_compressedEt()?  Never returns -1....
	hcalTPout->push_back(new_digi);//add to the list
      }//otherwise, pretend it never existed
    }//for each old digi
    
  }
  else {
    edm::LogWarning("MissingProduct") << "Missing input " << hcalSource_ << " Not going to generate output" << std::endl;
  }
  
  // store the output
  iEvent.put(hcalTPout);
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1CaloCorrector);

