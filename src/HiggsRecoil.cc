#include <HiggsRecoil.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"

using namespace std;

const float sqrts = 250.0; 	//GeV

HiggsRecoil a_HiggsRecoil_instance;

HiggsRecoil::HiggsRecoil()
	: Processor("HiggsRecoil"),
	_output(0)
{
	_description = "Print MC Truth" ;

	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	_colName="MCParticle";
	registerProcessorParameter( "MCObjects" ,
			"The name of the PFOs" ,
			_colName ,
			_colName);

	_treeName="MCPart";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	_leptonID = 13;
	registerProcessorParameter( "LeptonIDTag" ,
                        "Lepton ID that will be used in this analysis." ,
                        _leptonID ,
                        _leptonID);
	
	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

}

void HiggsRecoil::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("Num", &_Num, "Num/I");

	_outputTree->Branch("ZdaughterPID",&_ZdaughterPID,"ZDauPID/I");
	_outputTree->Branch("HdaughterPID",&_HdaughterPID,"HDauPID/I");
	_outputTree->Branch("HDaughterPID",&_HDaughterPID,"HDAUPID[2]/I");

	_outputTree->Branch("PzD1", _PzD1, "PzD1[4]/F");
	_outputTree->Branch("PzD2", _PzD2, "PzD2[4]/F");

	_outputTree->Branch("PhD1", _PhD1, "PhD1[4]/F");
	_outputTree->Branch("PhD2", _PhD2, "PhD2[4]/F");

	_outputTree->Branch("Pz",_Pz,"Pz[4]/F");
	_outputTree->Branch("Mz", &_Mz, "Mz/F");
	_outputTree->Branch("Ph",_Ph,"Ph[4]/F");
	_outputTree->Branch("NHDaug", &_NHDaug, "NHDaug/I");

	_outputTree->Branch("NNeutrino", &_NNeutrino, "NNeutrino/I");
	_outputTree->Branch("NZNeutrino", &_NZNeutrino, "NZNeutrino/I");
	_outputTree->Branch("ENeutrino", &_ENeutrino, "ENeutrino/F");


	_outputTree->Branch("NMuP", &_NMuP, "NMuP/I");
	_outputTree->Branch("NMuM", &_NMuM, "NMuM/I");

	_outputTree->Branch("NChP", &_NChP, "NChP/I");
        _outputTree->Branch("NChM", &_NChM, "NChM/I");
        _outputTree->Branch("NCh",  &_NCh,  "NCh/I");

	_outputTree->Branch("EventType", &_EventType, "EventType/I");

	_outputTree->Branch("P_MuP", _P_MuP, "P_MuP[4]/F");
        _outputTree->Branch("P_MuM", _P_MuM, "P_MuM[4]/F");
        _outputTree->Branch("P_Photon", _P_Photon, "P_Photon[4]/F");
        _outputTree->Branch("Pt_Z", &_Pt_Z, "Pt_Z/F");
        _outputTree->Branch("Pt_photon", &_Pt_photon, "Pt_photon/F");
        _outputTree->Branch("DeltaPt", &_DeltaPt, "DeltaPt/F");
        _outputTree->Branch("cosZ", &_cosZ, "cosZ/F");
        _outputTree->Branch("cosmis", &_cosmis, "cosmis/F");
        _outputTree->Branch("acop", &_acop, "acop/F");
        _outputTree->Branch("acol", &_acol, "acol/F");
        _outputTree->Branch("D_phi", &_D_phi, "D_phi/F");
	_outputTree->Branch("P_DL", _P_DL, "P_DL[4]/F");
	_outputTree->Branch("P_allCharged", _P_allCharged, "P_allCharged[4]/F");
	_outputTree->Branch("P_allNeutral", _P_allNeutral, "P_allNeutral[4]/F");
	_outputTree->Branch("P_allReco", _P_allReco, "P_allReco[4]/F");
	_outputTree->Branch("P_Higgs", _P_Higgs, "P_Higgs[4]/F");
	_outputTree->Branch("Hmass", &_Hmass, "Hmass/F");

	_outputTree->Branch("DiPhMass", &_DiPhMass, "DiPhMass/F");
	_outputTree->Branch("DiPhMassCorr", &_DiPhMassCorr, "DiPhMassCorr/F");

	_outputTree->Branch("InvMass", &_InvMass, "InvMass/F");
	_outputTree->Branch("RecoilMass", _RecoilMass, "RecoilMass[11]/F");

	_Num = 0;
}

void HiggsRecoil::processEvent( LCEvent * evtP ) 
{	

	if (evtP) 								
	{		
		try 	
		{    
			LCCollection* col_MCP = evtP->getCollection( _colName ) ;
			LCCollection* col_RecoP = evtP->getCollection( "ArborPFOs" );

			int _nMCP=col_MCP->getNumberOfElements();
			int _nRecoP = col_RecoP->getNumberOfElements();

			_eventNr=evtP->getEventNumber();

			int tmpPID = 0; 
			int NParent = 0; 
			int NDaughter = 0;
			float MCPEn = 0; 

			_Hmass = 0; 
			_Mz = 0; 
			_NMuP = 0;
			_NMuM = 0; 
			_NChP = 0;
			_NChM = 0;
			_NCh = 0;
			_EventType = -1; 
			_InvMass = 0; 
			_NNeutrino = 0; 
			_NZNeutrino = 0; 
			_ENeutrino = 0; 
			_DiPhMass = 0; 
			_DiPhMassCorr = 0; 
                        _Pt_Z = -100;
                        _D_phi = -100;
                        _cosZ = -100;

                        TLorentzVector ecms(0,0,0,250);

			for(int i0 = 0; i0 < 4; i0++)
			{
				_P_allCharged[i0] = 0;
				_P_allNeutral[i0] = 0;
				_P_allReco[i0] = 0;
				_P_Higgs[i0] = 0;
				_P_MuP[i0] = 0;
				_P_MuM[i0] = 0;
				_P_DL[i0] = 0;
				_PzD1[i0] = -100;
				_PzD2[i0] = -100;
				_PhD1[i0] = -100;
				_PhD2[i0] = -100;
				_Pz[i0] = 0;
				_Ph[i0] = 0;
				_P[i0] = 0;	//tmp
			}
			for(int i0 = 0; i0 < 11; i0++)
			{
                           _RecoilMass[i0]=0;
                        }

			TVector3 CurrPMomentum; 

			int RecoPID = 0; 
			int DIndex = 0;
			float RecoE = 0; 
			float RecoP[3] = {0, 0, 0};

			std::vector<TLorentzVector> FourMom_MuonP; 
			std::vector<TLorentzVector> FourMom_MuonM;
			std::vector<TLorentzVector> P_ChP; 
			std::vector<TLorentzVector> P_ChM; 
			std::vector<TLorentzVector> P_Ph; 
			
			TLorentzVector P_P, P_M, P_T[11];
			TLorentzVector P_totalPhoton(0, 0, 0, 0);
                        TLorentzVector beamp(0,0,125.0,125.0);
                        TLorentzVector beamm(0,0,-125.0,125.0);
                        float spreadfactor[11]={0.0,0.0004,0.0008,0.0012,0.0016,0.0020,0.0024,0.0028,0.0032,0.0036,0.0040};
                        for(int i=0; i<11; i++)
                        {
                           float scale1 = (gRandom->Gaus(1, spreadfactor[i]));
                           float scale2 = (gRandom->Gaus(1, spreadfactor[i]));
            		   P_T[i]=scale1*beamp+scale2*beamm;
                        }
			float currInvMass(0), currRecoilMass(0);
			float MinZThrDis = 1.0E10;

	
                        float photone=-1.;
			for(int i0 = 0; i0 < _nRecoP; i0++)
			{
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoP->getElementAt(i0));	//add some Photon ID as you like...
                                if(a_RecoP->getCharge()!=0) continue;
				_P_allNeutral[3] += a_RecoP->getEnergy();
				_P_allNeutral[0] += a_RecoP->getMomentum()[0];
				_P_allNeutral[1] += a_RecoP->getMomentum()[1];
				_P_allNeutral[2] += a_RecoP->getMomentum()[2];
				TLorentzVector currP( a_RecoP->getMomentum()[0], a_RecoP->getMomentum()[1], a_RecoP->getMomentum()[2], a_RecoP->getEnergy());
				if(_P_allNeutral[3] > 20)	//
				{
					P_Ph.push_back(currP);
					P_totalPhoton += currP;
				}
                               
                                if(a_RecoP->getType()!=22) continue;
                                if(fabs(a_RecoP->getMomentum()[2]/a_RecoP->getEnergy())>0.995) continue; // acceptance of the detector 
                                if(a_RecoP->getEnergy()>photone)
                                {
                                   photone = a_RecoP->getEnergy();
                                   _P_Photon[3] = a_RecoP->getEnergy();
                                   _P_Photon[0] = a_RecoP->getMomentum()[0];
                                   _P_Photon[1] = a_RecoP->getMomentum()[1];
                                   _P_Photon[2] = a_RecoP->getMomentum()[2];
                                   _Pt_photon = sqrt(_P_Photon[0]*_P_Photon[0]+_P_Photon[1]*_P_Photon[1]);
                                }
			}


			for(int j = 0; j < _nRecoP; j++)
			{
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoP->getElementAt(j));
                                if(a_RecoP->getCharge()==0) continue;
				RecoPID = a_RecoP->getType();
				RecoE = a_RecoP->getEnergy();
				RecoP[0] = a_RecoP->getMomentum()[0];
				RecoP[1] = a_RecoP->getMomentum()[1];
				RecoP[2] = a_RecoP->getMomentum()[2];

				TLorentzVector currP(RecoP[0], RecoP[1], RecoP[2], RecoE);
                               
                                if(RecoE>2.0) _NCh++;


				for(int s = 0; s < 4; s++)
				{
					_P_allCharged[s] += currP[s];
				}

				if( RecoE > 10 && RecoE < 100 )	//0.4*sqrt(s)
				{	
					if(abs(RecoPID) == _leptonID )	//Put by hand... guess enough
					{
						if(RecoPID == _leptonID )	//Got swapped...gosh!
						{
							FourMom_MuonM.push_back(currP);
						}
						else
						{
							FourMom_MuonP.push_back(currP);
						}
					}
					else if( a_RecoP->getCharge() > 0.5 )
					{
						P_ChP.push_back(currP);
					}
					else if( a_RecoP->getCharge() < -0.5 )
					{
						P_ChM.push_back(currP);

					}
				}
			}			

			_NMuP = FourMom_MuonP.size();
			_NMuM = FourMom_MuonM.size();
			_NChP = P_ChP.size();
			_NChM = P_ChM.size();

			std::vector<TLorentzVector> CandiP; 
			std::vector<TLorentzVector> CandiM;

			if(_NMuP && _NMuM)
			{
				_EventType = 0;
				CandiP = FourMom_MuonP;
				CandiM = FourMom_MuonM;
			}
			else if( _NMuP == 0 && _NChP && _NMuM)
			{
				_EventType = 1;
				CandiP = P_ChP;
				CandiM = FourMom_MuonM;
			}
			else if( _NMuM == 0 && _NChM && _NMuP )
			{
				_EventType = 2;
				CandiM = P_ChM;
				CandiP = FourMom_MuonP;
			}
			else if( _NMuM == 0 && _NChM && _NChP && _NMuP == 0 )
			{
				_EventType = 3;
				CandiM = P_ChM;
				CandiP = P_ChP;
			}

			int NCandiP = CandiP.size();
			int NCandiM = CandiM.size();
//                        float scale1 = (gRandom->Gaus(1, 0.0016))-1;
//                        float scale2 = (gRandom->Gaus(1, 0.0016))-1;

			if( NCandiP > 0 && NCandiM > 0 )
			{
				for(int p = 0; p < NCandiP; p++)
				{
					P_P = CandiP[p];

					for(int m = 0; m < NCandiM; m++)
					{
						P_M = CandiM[m];

						currInvMass = (P_P + P_M).M();
                                                
						if(fabs(currInvMass - 91.2) < MinZThrDis)
						{
                                                        MinZThrDis = fabs(currInvMass - 91.2);
							_InvMass = currInvMass; 
                                                        for(int i=0; i<11; i++)
                                                        { 
						            currRecoilMass = (P_T[i] - P_P - P_M).M();
                                                            _RecoilMass[i] = currRecoilMass;
                                                        } 
							for(int s = 0; s < 4; s++)
							{
								_P_MuP[s] = P_P[s];
								_P_MuM[s] = P_M[s];
								_P_DL[s] = _P_MuP[s] + _P_MuM[s];
							}
                                                        _acop = fabs(P_P.Phi()-P_M.Phi());
                                                        TLorentzVector miss = ecms - _P_DL;
                                                        _cosmis = miss.CosTheta();
                                                        _acol = P_P.Angle(P_M.Vect())*180./3.1415926;
                                                        _Pt_Z = sqrt(_P_DL[0]*_P_DL[0]+_P_DL[1]*_P_DL[1]);
                                                        _DeltaPt = _Pt_Z - _Pt_photon;
                                                        _cosZ = _P_DL[2]/sqrt(_P_DL[0]*_P_DL[0]+_P_DL[1]*_P_DL[1]+_P_DL[2]*_P_DL[2]);
                                                        float phi_p_tmp = atan2(_P_MuP[1],_P_MuP[0])*180./3.14159265;
                                                        float phi_m_tmp = atan2(_P_MuM[1],_P_MuM[0])*180./3.14159265;
                                                        if(_P_MuP[1] < 0) phi_p_tmp = phi_p_tmp + 360.;
                                                        if(_P_MuM[1] < 0) phi_m_tmp = phi_m_tmp + 360.;
                                                        _D_phi = fabs(phi_p_tmp - phi_m_tmp);
                                                        if (_D_phi > 180) _D_phi = 360. - _D_phi;
						}
					}
				}
			}

			//cout<<"NMu+ : "<<FourMom_MuonP.size()<<" NMu- : "<<FourMom_MuonM.size()<<"  Inv/Recoil Mass "<<_InvMass<<" : "<<_RecoilMass<<endl;
			//cout<<_NChP<<" : "<<_NChM<<endl; 

			for(int i2  = 0; i2 < 4; i2++)
			{
				_P_allReco[i2] = _P_allCharged[i2] + _P_allNeutral[i2];
				_P_Higgs[i2] = _P_allReco[i2] - _P_DL[i2]; 
			}

			_Hmass = sqrt( _P_Higgs[3]* _P_Higgs[3] - _P_Higgs[0]* _P_Higgs[0] - _P_Higgs[1]* _P_Higgs[1] - _P_Higgs[2]* _P_Higgs[2] );


			for(int i = 0; i < _nMCP; i++)
			{
				MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(col_MCP->getElementAt(i));
				tmpPID = a1_MCP->getPDG();
				NParent = a1_MCP->getParents().size();
				NDaughter = a1_MCP->getDaughters().size();
				MCPEn = a1_MCP->getEnergy();

				if(NParent == 0 && NDaughter == 1 && abs(tmpPID) < 20 )	//Including all Z decay
				{
					_ZdaughterPID = abs(tmpPID);
					_Pz[3] += MCPEn;
					_Pz[0] += a1_MCP->getMomentum()[0];
					_Pz[1] += a1_MCP->getMomentum()[1];
					_Pz[2] += a1_MCP->getMomentum()[2];
					CurrPMomentum = a1_MCP->getMomentum();

					if(_ZdaughterPID < 6 || _ZdaughterPID == 13 || _ZdaughterPID == 11 || _ZdaughterPID == 15)
					{
						if(tmpPID > 0)
						{
							_PzD1[0] = a1_MCP->getMomentum()[0];
							_PzD1[1] = a1_MCP->getMomentum()[1];
							_PzD1[2] = a1_MCP->getMomentum()[2];
							_PzD1[3] = MCPEn;
						}
						else if(tmpPID < 0)
						{
							_PzD2[3] = MCPEn;
							_PzD2[0] = a1_MCP->getMomentum()[0];
							_PzD2[1] = a1_MCP->getMomentum()[1];
							_PzD2[2] = a1_MCP->getMomentum()[2];
						}
					}
				}

				if(tmpPID == 25 && NDaughter > 1 && NParent !=0 ) //Higgs
				{
					_NHDaug = NDaughter;
					_Ph[3] = MCPEn;
					_Ph[0] = a1_MCP->getMomentum()[0];
					_Ph[1] = a1_MCP->getMomentum()[1];
					_Ph[2] = a1_MCP->getMomentum()[2];

					for(int i1 = 0; i1 < NDaughter; i1++)   //NDaughter == 2
					{
						MCParticle *b_MCP = a1_MCP->getDaughters()[i1];
						_HdaughterPID = b_MCP->getPDG();
						CurrPMomentum = b_MCP->getMomentum();
                                                _HDaughterPID[DIndex] = _HdaughterPID;
						
						if(_HdaughterPID < 6 || _HdaughterPID == 13 || _HdaughterPID == 11 || _HdaughterPID == 15 )
						{
							if(_HdaughterPID > 0)
							{
								_PhD1[3] = b_MCP->getEnergy();
								_PhD1[0] = b_MCP->getMomentum()[0];
								_PhD1[1] = b_MCP->getMomentum()[1];
								_PhD1[2] = b_MCP->getMomentum()[2];
							}
							else if(_HdaughterPID < 0)
							{
								_PhD2[3] = b_MCP->getEnergy();
								_PhD2[0] = b_MCP->getMomentum()[0];
								_PhD2[1] = b_MCP->getMomentum()[1];
								_PhD2[2] = b_MCP->getMomentum()[2];
							}
						}
						else if(DIndex == 0)
						{	
							_PhD1[3] = b_MCP->getEnergy();
							_PhD1[0] = b_MCP->getMomentum()[0];
							_PhD1[1] = b_MCP->getMomentum()[1];
							_PhD1[2] = b_MCP->getMomentum()[2];
							DIndex = 1; 
						}
						else if(DIndex == 1)
						{
							_PhD2[3] = b_MCP->getEnergy();
							_PhD2[0] = b_MCP->getMomentum()[0];
							_PhD2[1] = b_MCP->getMomentum()[1];
							_PhD2[2] = b_MCP->getMomentum()[2];
						}
                                                if(_HdaughterPID==23)
                                                {
                                                        NDaughter = b_MCP->getDaughters().size();
                                                        for(int i=0; i<NDaughter; i++)
                                                        {
                                                            MCParticle *c_MCP = b_MCP->getDaughters()[i];
                                                            if(abs(c_MCP->getPDG())==12 || abs(c_MCP->getPDG())==14 || abs(c_MCP->getPDG())==16)
                                                            _NZNeutrino++;
                                                        }
                                                }
						_HdaughterPID = abs(_HdaughterPID);
					}
				}

				if(NDaughter == 0)
				{
					if(abs(tmpPID) == 12 || abs(tmpPID) == 14 || abs(tmpPID) == 16)
					{
						_NNeutrino++;
						_ENeutrino += MCPEn;  
					}
				}	
			}
			_Mz = sqrt(_Pz[3]*_Pz[3] - _Pz[0]*_Pz[0]- _Pz[1]*_Pz[1]- _Pz[2]*_Pz[2]);

			_outputTree->Fill();
			_Num++;

		}		
		catch (lcio::DataNotAvailableException err) { }

	}  	  

}	

void HiggsRecoil::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}



