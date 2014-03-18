#ifndef __L1Analysis_L1AnalysisDataFormat_H__
#define __L1Analysis_L1AnalysisDataFormat_H__

#include <vector>
// #include <inttypes.h>
#include <TROOT.h>

namespace L1Analysis
{
  struct L1AnalysisDataFormat
  {
    L1AnalysisDataFormat(){Reset();};
    ~L1AnalysisDataFormat(){};
    
    void Reset()
    { 

      //Event Info
		 Run = -1;
		 Event = -1;
		 LS = -1;
		 

    	//EG info
       Nele = 0;
       Bxel.clear();
       Etel.clear();
       Phiel.clear();
       Etael.clear();
       Isoel.clear();

    	//tkEG info
       NTkele = 0;
       BxTkel.clear();
       EtTkel.clear();
       PhiTkel.clear();
       EtaTkel.clear();
       IsoTkel.clear();

    	//track Taus info
       NTktau = 0;
       BxTktau.clear();
       EtTktau.clear();
       PhiTktau.clear();
       EtaTktau.clear();
       IsoTktau.clear();

       
       // calo Jets
       Njet = 0;
       Bxjet.clear();
       Etjet.clear();
       Phijet.clear();
       Etajet.clear();
       Taujet.clear();
       isoTaujet.clear();
       Fwdjet.clear();  

       // Trk Jets
       NTkjet = 0;
       BxTkjet.clear();
       EtTkjet.clear();
       PhiTkjet.clear();
       EtaTkjet.clear();
       
	
	// muons	 
       Nmu = 0;
       Bxmu.clear();
       Ptmu.clear();
       Phimu.clear();
       Etamu.clear();
       Qualmu.clear();
       Isomu.clear();

	// track muons	 
       NTkmu = 0;
       BxTkmu.clear();
       PtTkmu.clear();
       PhiTkmu.clear();
       EtaTkmu.clear();
       QualTkmu.clear();
       IsoTkmu.clear();
		 

// ------ ETT, ETM, HTT and HTM from PSB14:

	ETT = -1.;
	OvETT = false;
	HTT = -1.;
	OvHTT = false;
	ETM = -1.;
	PhiETM = -1;
	OvETM = false;
	HTM = -1.;
	PhiHTM = -1;
	OvHTM = false;
	
	
// ------- Track base "Energy Sum" Quantities
        TkETT   = -1.;
	TkETM = -1.;
	TkETMPhi = -1.;
	
	
	TkHTT = -1.;
        TkHTM = -1.;
        TkHTMPhi = -1.;
		
    }
   
    // ---- L1AnalysisDataFormat information.
   
    //
	 int Run;
	 int Event;
	 int LS;
	 
    // info
    int            Nele;
    std::vector<int>    Bxel;
    std::vector<float>  Etel;
    std::vector<float>  Phiel;
    std::vector<float>  Etael;
    std::vector<bool>   Isoel;

    int            NTkele;
    std::vector<int>    BxTkel;
    std::vector<float>  EtTkel;
    std::vector<float>  PhiTkel;
    std::vector<float>  EtaTkel;
    std::vector<float>  zVtxTkel;
    std::vector<bool>   IsoTkel;
    std::vector<float>  tIsoTkel; 

    int            NTktau;
    std::vector<int>    BxTktau;
    std::vector<float>  EtTktau;
    std::vector<float>  PhiTktau;
    std::vector<float>  EtaTktau;
    std::vector<float>  zVtxTktau;
    std::vector<bool>   IsoTktau;
    std::vector<float>  tIsoTktau; 

    
    int            NTkjet;
    std::vector<int>    BxTkjet;
    std::vector<float>  EtTkjet;
    std::vector<float>  PhiTkjet;
    std::vector<float>  EtaTkjet;
    std::vector<float>  zVtxTkjet;

    int            Njet;
    std::vector<int>    Bxjet;
    std::vector<float>  Etjet;
    std::vector<float>  Phijet;
    std::vector<float>  Etajet;
    std::vector<bool>   Taujet;
    std::vector<bool>   isoTaujet;
    std::vector<bool>   Fwdjet;



    int            Nmu;
    std::vector<int>    Bxmu;
    std::vector<float>  Ptmu;
    std::vector<float>  Phimu;
    std::vector<float>  Etamu;
    std::vector<int>    Qualmu;
    std::vector<float>  Isomu;
    
    int            NTkmu;
    std::vector<int>    BxTkmu;
    std::vector<float>  PtTkmu;
    std::vector<float>  PhiTkmu;
    std::vector<float>  EtaTkmu;
    std::vector<int>    QualTkmu;
    std::vector<float>  IsoTkmu;  
    std::vector<float>  zVtxTkmu; 
    std::vector<float>  tIsoTkmu;  

// ------ ETT, ETM, HTT and HTM :

    float ETT;
    bool OvETT;

    float HTT;
    bool OvHTT;

    float ETM;
    int PhiETM;
    bool OvETM;

    float HTM;
    int PhiHTM;
    bool OvHTM;

// ------- Track base "Energy Sum" Quantities
   float TkETT;
   float TkETM;
   float TkETMPhi;

   float TkHTT;
   float TkHTM;
   float TkHTMPhi;
    

  }; 
} 
#endif


