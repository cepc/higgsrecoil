#ifndef _HiggsRecoil_hh_
#define _HiggsRecoil_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class HiggsRecoil  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new HiggsRecoil ; }

		HiggsRecoil();

		~HiggsRecoil() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;

		int _overwrite, _leptonID;
		float _cmsE; 
		TTree *_outputTree;

		int _NMuP, _NMuM, _NChP, _NChM, _NCh;
		float _P_MuP[4], _P_MuM[4], _P_DL[4];
		int _EventType; 
		float _InvMass;
                float _RecoilMass[11]; 

		int _PID1, _PID2;
		float _PL1[4], _PL2[4], _RPL1[4], _RPL2[4], _SM[4], _P_allCharged[4], _P_allNeutral[4], _P_Higgs[4], _P_allReco[4], _P_Photon[4];
		float _Hmass; 
		int _Num;
		int _NHDaug; 
		int _HdaughterPID, _HDaughterPID[2]; 
		int _ZdaughterPID; 
		float _Pz[4], _Ph[4], _PzD1[4], _PzD2[4], _PhD1[4], _PhD2[4], _RPzD1[4], _RPzD2[4], _RPhD1[4], _RPhD2[4];
		float _P[4], _SumP[4], _VisP[4], _MissP[4];
		int _PID, _NFMCP, _MotherFlag, _NNeutrino, _NZNeutrino; 
		float _ENeutrino, _DiPhMass, _DiPhMassCorr; 
//		float _CosTheta, _Phi, _Charge;
		float _Mz, _Mrecoil, _MzReco, _MhReco, _MrecoilReco; 
                float _Pt_Z, _Pt_photon, _D_phi, _cosZ, _cosmis, _acol, _acop, _DeltaPt;
		float KthEn[7][9];
		unsigned int _eventNr;

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


