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
#include <iostream>

//ROOT
#include "TObject.h"
#include "TVector3.h"

class THistoryData : public TObject {

 private:
vector<Int_t>		fHistoryPrimaryID; 
vector<Int_t> 		fHistoryPrimaryPdg; // the type of the particle generated in the source e.g. pdg = 22 (gamma) in  0 --> gamma (primary) --> e-   --> gamma  --> e-
vector<Double_t>	fHistoryPrimaryEnergy; 
vector<TVector3> 	fHistoryPrimaryTrajectory; // NOTA BENE : Only for the first emitted primary charged particle! useful for e-/e+ in a magnetic field.

vector<TVector3> 	fHistoryPrimaryPositionVertex;
vector<TVector3> 	fHistoryPrimaryMomentumVertex;
vector<TVector3> 	fHistoryPrimaryPosition1stImpact;
vector<TVector3> 	fHistoryPrimaryMomentum1stImpact; 
vector<TString> 	fHistoryPrimaryVolume1stImpact; 

vector<Int_t>		fHistoryParentID; 

vector<Int_t> 		fHistoryCurrentID; 
vector<Int_t> 		fHistoryCurrentPdg; //the type of the last particle produced  e.g. pdg = 11 (e-) in  0 --> gamma (primary) --> e-   --> gamma  --> e-
vector<Double_t>	fHistoryCurrentEnergy;
vector<Double_t>	fHistoryCurrentTime; // time when the detector is hit
vector<TString> 	fHistoryCurrentCreatorProcess; 
vector<TString> 	fHistoryCurrentVolume1stImpact; 

vector<TVector3> 	fHistoryCurrentPositionVertex;
vector<TVector3> 	fHistoryCurrentPosition1stImpact; 
vector<TVector3> 	fHistoryCurrentPositionDetector; 
vector<TVector3> 	fHistoryCurrentPositionDeath; 

vector<TVector3> 	fHistoryCurrentMomentumVertex;
vector<TVector3> 	fHistoryCurrentMomentum1stImpact; 
vector<TVector3> 	fHistoryCurrentMomentumDetector; 
vector<TVector3> 	fHistoryCurrentMomentumDeath; 

   
 public:
   THistoryData();
   virtual ~THistoryData();

   void   ClearVariables();
   void   Dump() const;



   /////////////////////           GETTERS           ////////////////////////
inline Int_t    	GetHistoryMult(void)		    { return fHistoryPrimaryEnergy.size(); } 
inline Double_t		GetHistoryPrimaryID(int i)		{  return fHistoryPrimaryID.at(i); }
inline Int_t 		GetHistoryPrimaryPdg(int i)		{  return fHistoryPrimaryPdg.at(i); } 
inline Double_t		GetHistoryPrimaryEnergy(int i)	{  return fHistoryPrimaryEnergy.at(i); } 
inline TVector3		GetHistoryPrimaryTrajectory(int i)			{  return fHistoryPrimaryTrajectory.at(i); } 

inline TVector3	GetHistoryPrimaryPositionVertex(int i)			{  return fHistoryPrimaryPositionVertex.at(i); } 
inline TVector3	GetHistoryPrimaryMomentumVertex(int i)			{  return fHistoryPrimaryMomentumVertex.at(i); }
inline TVector3	GetHistoryPrimaryMomentum1stImpact(int i)		{  return fHistoryPrimaryMomentum1stImpact.at(i); }
inline TVector3	GetHistoryPrimaryPosition1stImpact(int i)		{  return fHistoryPrimaryPosition1stImpact.at(i); }
inline TString	GetHistoryPrimaryVolume1stImpact(int i)			{  return fHistoryPrimaryVolume1stImpact.at(i); } 

inline Int_t	GetHistoryParentID(int i)	{  return fHistoryParentID.at(i); }

inline Double_t	GetHistoryCurrentID(int i)				{  return fHistoryCurrentID.at(i); }
inline Int_t 	GetHistoryCurrentPdg(int i)				{  return fHistoryCurrentPdg.at(i); } 
inline Double_t	GetHistoryCurrentEnergy(int i)			{  return fHistoryCurrentEnergy.at(i); } 
inline Double_t	GetHistoryCurrentTime(int i)			{  return fHistoryCurrentTime.at(i); } 
inline TString	GetHistoryCurrentCreatorProcess(int i)	{  return fHistoryCurrentCreatorProcess.at(i); } 

inline TVector3	GetHistoryCurrentPositionVertex(int i)		{  return fHistoryCurrentPositionVertex.at(i); } 
inline TVector3	GetHistoryCurrentPosition1stImpact(int i)	{  return fHistoryCurrentPosition1stImpact.at(i); }
inline TVector3	GetHistoryCurrentPositionDetector(int i)	{  return fHistoryCurrentPositionDetector.at(i); } 
inline TVector3	GetHistoryCurrentPositionDeath(int i)		{  return fHistoryCurrentPositionDeath.at(i); } 

inline TVector3	GetHistoryCurrentMomentumVertex(int i)		{  return fHistoryCurrentMomentumVertex.at(i); } 
inline TVector3	GetHistoryCurrentMomentum1stImpact(int i)	{  return fHistoryCurrentMomentum1stImpact.at(i); }
inline TVector3	GetHistoryCurrentMomentumDetector(int i)	{  return fHistoryCurrentMomentumDetector.at(i); } 
inline TVector3	GetHistoryCurrentMomentumDeath(int i)		{  return fHistoryCurrentMomentumDeath.at(i); } 


   /////////////////////           SETTERS           ////////////////////////
inline void		SetHistoryPrimaryID(Double_t i)		{  fHistoryPrimaryID.push_back(i); }
inline void 	SetHistoryPrimaryPdg(Int_t i)		{  fHistoryPrimaryPdg.push_back(i); } 
inline void		SetHistoryPrimaryEnergy(Double_t i)	{  fHistoryPrimaryEnergy.push_back(i); } 
inline void		SetHistoryPrimaryTrajectory(Double_t x, Double_t y, Double_t z)	{ fHistoryPrimaryTrajectory.push_back(TVector3(x,y,z)); } 

inline void		SetHistoryPrimaryPositionVertex(Double_t x, Double_t y, Double_t z)			{   fHistoryPrimaryPositionVertex.push_back(TVector3(x,y,z)); } 
inline void		SetHistoryPrimaryMomentumVertex(Double_t x, Double_t y, Double_t z)			{   fHistoryPrimaryMomentumVertex.push_back(TVector3(x,y,z)); }
inline void		SetHistoryPrimaryMomentum1stImpact(Double_t x, Double_t y, Double_t z)		{   fHistoryPrimaryMomentum1stImpact.push_back(TVector3(x,y,z)); }
inline void		SetHistoryPrimaryPosition1stImpact(Double_t x, Double_t y, Double_t z)		{  	fHistoryPrimaryPosition1stImpact.push_back(TVector3(x,y,z)); }
inline void		SetHistoryPrimaryVolume1stImpact(TString i)									{   fHistoryPrimaryVolume1stImpact.push_back(i); } 

inline void		SetHistoryParentID(Int_t i)	{  fHistoryParentID.push_back(i); }

inline void		SetHistoryCurrentCreatorProcess(TString i)			{  fHistoryCurrentCreatorProcess.push_back(i); } 
inline void		SetHistoryCurrentID(Double_t i)						{  fHistoryCurrentID.push_back(i); }
inline void 	SetHistoryCurrentPdg(Int_t i)						{  fHistoryCurrentPdg.push_back(i); } 
inline void		SetHistoryCurrentEnergy(Double_t i)					{  fHistoryCurrentEnergy.push_back(i); } 
inline void		SetHistoryCurrentTime(Double_t i)					{  fHistoryCurrentTime.push_back(i); } // Time whe the detector is hit 
inline void		SetHistoryCurrentVolume1stImpact(TString i)			{  fHistoryCurrentVolume1stImpact.push_back(i); } 


inline void		SetHistoryCurrentPositionVertex(Double_t x, Double_t y, Double_t z)		{  fHistoryCurrentPositionVertex.push_back(TVector3(x,y,z)); } 
inline void		SetHistoryCurrentPosition1stImpact(Double_t x, Double_t y, Double_t z)	{  fHistoryCurrentPosition1stImpact.push_back(TVector3(x,y,z)); }
inline void 	SetHistoryCurrentPositionDetector(Double_t x, Double_t y, Double_t z)	{  fHistoryCurrentPositionDetector.push_back(TVector3(x,y,z)); } 
inline void		SetHistoryCurrentPositionDeath(Double_t x, Double_t y, Double_t z)		{  fHistoryCurrentPositionDeath.push_back(TVector3(x,y,z)); } 

inline void		SetHistoryCurrentMomentumVertex(Double_t x, Double_t y, Double_t z)		{  fHistoryCurrentMomentumVertex.push_back(TVector3(x,y,z)); } 
inline void		SetHistoryCurrentMomentum1stImpact(Double_t x, Double_t y, Double_t z)	{  fHistoryCurrentMomentum1stImpact.push_back(TVector3(x,y,z)); }
inline void 	SetHistoryCurrentMomentumDetector(Double_t x, Double_t y, Double_t z)	{  fHistoryCurrentMomentumDetector.push_back(TVector3(x,y,z)); } 
inline void		SetHistoryCurrentMomentumDeath(Double_t x, Double_t y, Double_t z)		{  fHistoryCurrentMomentumDeath.push_back(TVector3(x,y,z)); } 
                 
   ClassDef(THistoryData,1)  // HistoryData structure
};

#endif
