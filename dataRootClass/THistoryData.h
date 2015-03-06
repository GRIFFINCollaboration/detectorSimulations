#ifndef __HistoryDATA__
#define __HistoryDATA__
/*****************************************************************************
 * Original Author: Mhd Moukaddam  contact address: moukaddam@triumf.ca      *
 *---------------------------------------------------------------------------*
 * Decription: This class stores selected information of the G4 simulation   * 
 *             to help in analysing the final simulated physics spectra.	 *	
 *             Information such as, the number of generation in a cascade,   *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/


//c++
#include <vector>

//ROOT
#include "TObject.h"
#include "TVector3.h"

class THistoryData : public TObject {

 private:
vector<Double_t>	fHistoryPrimaryEnergy; // the type of the particle generated in the source e.g. energy of gamma (primary) in  0 --> gamma (primary) --> e-   --> gamma  --> e-
vector<Int_t> 		fHistoryPrimaryPdg; // the type of the particle generated in the source e.g. pdg = 22 (gamma) in  0 --> gamma (primary) --> e-   --> gamma  --> e-
vector<Int_t> 		fHistoryCurrentPdg; //the type of the last particle produced  e.g. pdg = 11 (e-) in  0 --> gamma (primary) --> e-   --> gamma  --> e-
vector<Int_t> 		fHistoryGNumber; // Number of generations leading to the last particle e.g 4 in : 0 --> gamma (primary) --> e-   --> gamma  --> e- 
vector<TString> 	fHistoryG2BirthVolume; // Birth volume of the 2nd generation e.g magnets in : 0 --> gamma (primary from source) --> e- (magnets) --> gamma (magnets) --> ...
vector<TString> 	fHistoryGLastBirthVolume; // Birth volume of the last generat. e.g SiLi in : 0 --> gamma (primary from source) --> e- (magnets) --> gamma (magnets) --> e- (SiLi)
vector<TString> 	fHistoryG2Process; // The process at the second generation e.g  Compton in this example :  0---> gamma (primary) ---> electron by Compton ----> gamma by Bremsstrahlung ....

   
 public:
   THistoryData();
   virtual ~THistoryData();

   void   Clear();
   void   Dump() const;



   /////////////////////           GETTERS           ////////////////////////
inline Int_t    	GetHistoryMult(void)		    { return fHistoryPrimaryEnergy.size(); } 
inline Double_t    	GetHistoryPrimaryEnergyAt(int i){ return fHistoryPrimaryEnergy.at(i); }
inline Int_t 		GetHistoryPrimaryPdgAt(int i){ return fHistoryPrimaryPdg.at(i); } 
inline Int_t		GetHistoryCurrentPdgAt(int i){ return fHistoryCurrentPdg.at(i); } 
inline Int_t		GetHistoryGNumberAt(int i){ return fHistoryGNumber.at(i); } 
inline TString		GetHistoryG2BirthVolumeAt(int i){ return fHistoryG2BirthVolume.at(i); } 
inline TString		GetHistoryGLastBirthVolumeAt(int i){ return fHistoryGLastBirthVolume.at(i); } 
inline TString		GetHistoryG2ProcessAt(int i){ return fHistoryG2Process.at(i); } 

inline vector<Double_t>    	GetHistoryPrimaryEnergy(void){ return fHistoryPrimaryEnergy; }
inline vector<Int_t>		GetHistoryPrimaryPdg(void){ return fHistoryPrimaryPdg; } 
inline vector<Int_t>		GetHistoryCurrentPdg(void){ return fHistoryCurrentPdg; } 
inline vector<Int_t>		GetHistoryGNumber(void){ return fHistoryGNumber; } 
inline vector<TString>	GetHistoryG2BirthVolume(void){ return fHistoryG2BirthVolume; } 
inline vector<TString>	GetHistoryGLastBirthVolume(void){ return fHistoryGLastBirthVolume; } 
inline vector<TString>	GetHistoryG2Process(void){ return fHistoryG2Process; } 

   
   /////////////////////           SETTERS           ////////////////////////
inline void		SetHistoryPrimaryEnergy(Double_t i){  fHistoryPrimaryEnergy.push_back(i); }
inline void 	SetHistoryPrimaryPdg(Int_t i){  fHistoryPrimaryPdg.push_back(i); } 
inline void		SetHistoryCurrentPdg(Int_t i){  fHistoryCurrentPdg.push_back(i); } 
inline void		SetHistoryGNumber(Int_t i){  fHistoryGNumber.push_back(i); } 
inline void		SetHistoryG2BirthVolume(TString i){  fHistoryG2BirthVolume.push_back(i); } 
inline void		SetHistoryGLastBirthVolume(TString i){  fHistoryGLastBirthVolume.push_back(i); } 
inline void 	SetHistoryG2Process(TString i){ fHistoryG2Process.push_back(i); } 
                 
   ClassDef(THistoryData,1)  // HistoryData structure
};

#endif
