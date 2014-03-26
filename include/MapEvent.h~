#ifndef _RAWEVENT_CLASS
#define _RAWEVENT_CLASS


//Root
#include "Riostream.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TH1.h" 
#include "TH3D.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
using namespace ROOT::Math;

//c++
#include "string"
#include "vector"
#include "iostream"
#include "fstream"
#include "sstream"
#include "math.h"
#include "map"
#include "utility"
#include "algorithm"
#include "stdexcept"
#include "exception"
#include "stdio.h"
#include "stdlib.h"
using namespace std;
              
              

class MapEvent
{

 private:

//From geant, HC stands for hit collection
vector<Int_t> fHCPdg; // type of detected particles
vector<Double_t> fHCEnergy; // depositied energy
vector<TVector3> fHCPosition; // position of the hit, we are intersted with the first position

vector<Int_t> fHCPrimaryID; // the ID number of the mother particle generated in the source
vector<Int_t> fHCPrimaryPdg; // the type of the particle generated in the source
vector<Double_t> fHCPrimaryEnergy; // usefull for sorting and for analysis
vector<TVector3> fHCPrimaryMoment; // Momement of the particle, we are intersted with the first position

// Sorted vectors
vector<Int_t>    fDetectedPrimaryParticleID ; // If we have several primary particle ending in the same pad
vector<Double_t> fDetectedPrimaryEnergy; // the DETECTED energy of the primary particle, this will be the sum of all the energy from the same ID

map<int,int> fMapPrimaryPdg;
map<double,int> fMapPrimaryEnergy;
	
						/************************************************
							  			M E T H O D S 
						************************************************/

/*************************
C O N  (D E) -  S T R U C T E R S 
************************/
 public:
MapEvent(void );
~MapEvent(void);


/*************************
      S O R T E R S  
************************/
 public:
void Treat(void);



/*************************
      S E T T E R S 
************************/
 public:
 
void FillVectors(int pdg, // particle pdg
 				 double Energy, // depositid energy
				 double Px, double Py, double Pz, // position vector
				 
				 int ID, // original(primary) ID
				 int PrimPdg,// primary particle definition        PDG encoding				 
				 double PrimEnergy,//original(primary) energy
				 double Mx, double My, double Mz); // primary particle momentum vector

/*************************
      G E T T E R S 
************************/
 public:
 
// from Hit collections
Int_t       GetHCPrimaryPdg(int i);// say you have in one pad two primary particles  1 1 1 1 1 1 0 0 0  2 2 2 2 2 => this function will return  1 for i=0, 1 for i=1 and  1 for i=2 ...
Double_t    GetHCPrimaryEnergy(int i);

// after treatment
Double_t 	GetFullEnergy(void) ; // full energy in pad
TVector3 	GetFirstHitPosition(void) ; // first hit position
                                  // the reason for this treatment is because our basic element is a pad not one hit from the collection of Hit
Int_t    	GetPrimaryPdgMult(void); // say you have in one pad two primary particles  1 1 1 1 1 1 0 0 0  2 2 2 2 2 => this function will return '3'
Int_t    	GetPrimaryPdg(int i); // say you have in one pad two primary particles 1 1 1 1 1 1 0 0 0  2 2 2 2 2 => 1 - 0 - 2 => this function will return  1 for i=0, 0 for i=1 and 2 for i=2

Int_t    	GetPrimaryEnergyMult(void); 
Double_t    GetPrimaryEnergy(int i); 




// getters sorted per primary ID
Int_t    GetPrimaryPdgForID(int ID) ; // for a specific type of particle PDG
Double_t GetPrimaryEnergyForID(int ID) ; // for a specific type of particle PDG
Double_t GetDetectedEnergyForID(int ID) ; // full energy in pad for a specific type of particle PDG
Double_t GetAngleOfEmissionForID(int ID) ; // for a specific type of particle PDG // from the moment

/*

// Per particle type (PDG)
Double_t GetInputEnergyForPdg(int pdg) ; // full energy in pad for a specific type of particle PDG
//Int_t    GetPrimaryIDForPdg(int pdg) ; // for a specific type of particle PDG          <---- could produce conflict
Int_t    GetPrimaryParticleDefForPdg(int pdg) ; // for a specific type of particle PDG
Double_t GetPrimaryEnergyForPdg(int pdg) ; // for a specific type of particle PDG
Double_t GetDetectedEnergyForPdg(int pdg) ; // full energy in pad for a specific type of particle PDG
Double_t GetAngleOfEmissionForPdg(int pdg) ; // for a specific type of particle PDG // from the moment

// per index
Double_t GetPDG(int index) ; // full energy in pad for a specific type of particle PDG
Int_t    GetPrimaryID(int index) ; // for a specific type of particle PDG
Int_t    GetPrimaryParticleDef(int index) ; // for a specific type of particle PDG
Double_t GetPrimaryEnergy(int index) ; // for a specific type of particle PDG
Double_t GetDetectedEnergy(int index) ; // full energy in pad for a specific type of particle PDG
Double_t GetAngleOfEmission(int index) ; // for a specific type of particle PDG // from the moment

*/


/*************************
      H E L P E R S 
************************/
void ShowVectorsContent(void) ;

 };
 

#endif



