#ifndef _RAWEVENT_CLASS
#define _RAWEVENT_CLASS

//Root
#include "TVector3.h"

//c++
#include "vector"
#include "iostream"
#include "map"
#include "utility"
using namespace std;

//g4              
#include "globals.hh"             

class RawG4Event
{

 private:

//From geant, HC stands for hit collection
vector<Int_t> 		fHCPdg; // type of detected particles
vector<Double_t> 	fHCEnergy; // depositied energy
vector<Int_t> 		fHCDetector; // Detector
vector<Int_t> 		fHCCrystal; // Crystal
vector<TVector3> 	fHCPosition; // position of the hit, we are intersted with the first position

vector<Int_t> 		fHCPrimaryID; // the ID number of the mother particle generated in the source
vector<Int_t> 		fHCPrimaryPdg; // the type of the particle generated in the source
vector<Double_t> 	fHCPrimaryEnergy; // usefull for sorting and for analysis
vector<TVector3> 	fHCPrimaryMomentum; // Momentum of the particle,
vector<Double_t> 	fHCPrimaryTheta; // Angle theta with respect to the beam of the particle,
vector<Double_t> 	fHCPrimaryPhi; // Angle Phi the particle,

// Sorted vectors
vector<Int_t>    	fDetectedPrimaryParticleID ; // If we have several primary particle ending in the same pad
vector<Double_t> 	fDetectedPrimaryEnergy; // the DETECTED energy of the primary particle, this will be the sum of all the energy from the same ID

map<int,int> 		fMapPrimaryPdg;
map<double,int> 	fMapPrimaryEnergy;
map<double,int> 	fMapPrimaryTheta;
map<double,int> 	fMapPrimaryPhi;


/********************************************
C O N S T R U C T O R -  D E S T R U C T E R 
*********************************************/
 public:
RawG4Event(void );
~RawG4Event(void);


/*************************
      S O R T E R S  
************************/
 public:
void SortPrimary(void);



/*************************
      S E T T E R S 
************************/
 public:
 
void FillVectors(int pdg, // particle pdg
 				 double Energy, // depositid energy
 				 int detector, //detector
 				 int crystal, //crystal
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
TVector3    GetHCPrimaryMomentum(int i);
Double_t    GetHCPrimaryTheta(int i);
Double_t    GetHCPrimaryPhi(int i);

// after treatment
Double_t 	GetFullEnergy(void) ; // full energy in pad
TVector3 	GetFirstHitPosition(void) ; // first hit position
TVector3 	GetSecondHitPosition(void) ; // first hit position
Int_t 	 	GetDetector(void) ; // detector
Int_t 	 	GetCrystal(void) ; // crystal
                                  // the reason for this treatment is because our basic element is a pad not one hit from the collection of Hit
Int_t    	GetPrimaryPdgMult(void); // say you have in one pad two primary particles  1 1 1 1 1 1 0 0 0  2 2 2 2 2 => this function will return '3'
Int_t    	GetPrimaryPdg(int i); // say you have in one pad two primary particles 1 1 1 1 1 1 0 0 0  2 2 2 2 2 => 1 - 0 - 2 => this function will return  1 for i=0, 0 for i=1 and 2 for i=2

Int_t    	GetPrimaryEnergyMult(void); 
Double_t    GetPrimaryEnergy(int i); 

Int_t    	GetPrimaryThetaMult(void); 
Double_t    GetPrimaryTheta(int i);
Double_t    GetPrimaryPhi(int i); 

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



