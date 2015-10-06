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

#include <iostream>
using namespace std;

#include "THistoryData.h"


ClassImp(THistoryData)

THistoryData::THistoryData() {
   ClearVariables();
}

THistoryData::~THistoryData() {}


void THistoryData::ClearVariables(){
	
	fHistoryPrimaryID.clear(); 
	fHistoryPrimaryPdg.clear(); 
	fHistoryPrimaryEnergy.clear(); 

	fHistoryPrimaryPositionVertex.clear();
	fHistoryPrimaryMomentumVertex.clear();
	fHistoryPrimaryPosition1stImpact.clear();
	fHistoryPrimaryMomentum1stImpact.clear(); 
	fHistoryPrimaryVolume1stImpact.clear(); 

	fHistoryParentID.clear(); 

	fHistoryCurrentCreatorProcess.clear(); 
	fHistoryCurrentID.clear(); 
	fHistoryCurrentPdg.clear(); 
	fHistoryCurrentEnergy.clear(); 
	fHistoryCurrentTime.clear(); 
	fHistoryCurrentVolume1stImpact.clear(); 
	
	fHistoryCurrentPositionVertex.clear();
	fHistoryCurrentPosition1stImpact.clear(); 
	fHistoryCurrentPositionDetector.clear(); 
	fHistoryCurrentPositionDeath.clear(); 

	fHistoryCurrentMomentumVertex.clear();
	fHistoryCurrentMomentum1stImpact.clear(); 
	fHistoryCurrentMomentumDetector.clear(); 
	fHistoryCurrentMomentumDeath.clear(); 

}


void THistoryData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX HISTORY XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

   cout << "History_Mult = " <<fHistoryPrimaryID.size() << endl;
   
   for (UShort_t i = 0; i < fHistoryPrimaryID.size(); i++)
      cout << "Primary ID: " << fHistoryPrimaryID[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryPrimaryPdg.size(); i++)
      cout << "Primary Pdg: " << fHistoryPrimaryPdg[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryPrimaryEnergy.size(); i++)
      cout << "Primary Energy: " << fHistoryPrimaryEnergy[i] << endl;
 
 	for (UShort_t i = 0; i < fHistoryPrimaryPositionVertex.size(); i++)
		cout << " History Primary Position Vertex    : " << fHistoryPrimaryPositionVertex[i].X() << " " << fHistoryPrimaryPositionVertex[i].Y() << " " << fHistoryPrimaryPositionVertex[i].Z() << endl;	
	for (UShort_t i = 0; i < fHistoryPrimaryMomentumVertex.size(); i++)
		cout << " History Primary Momentum Vertex    : " << fHistoryPrimaryMomentumVertex[i].X() << " " << fHistoryPrimaryMomentumVertex[i].Y() << " " << fHistoryPrimaryMomentumVertex[i].Z() << endl;

	for (UShort_t i = 0; i < fHistoryPrimaryVolume1stImpact.size(); i++)
		cout << " History Primary volume 1st Impact : " << fHistoryPrimaryVolume1stImpact[i] << endl;
	for (UShort_t i = 0; i < fHistoryPrimaryPosition1stImpact.size(); i++)
		cout << " History Primary Position 1st Impact : " << fHistoryPrimaryPosition1stImpact[i].X() << " " << fHistoryPrimaryPosition1stImpact[i].Y() << " " << fHistoryPrimaryPosition1stImpact[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryPrimaryMomentum1stImpact.size(); i++)
		cout << " History Primary Momentum 1st Impact : " << fHistoryPrimaryMomentum1stImpact[i].X() << " " << fHistoryPrimaryMomentum1stImpact[i].Y() << " " << fHistoryPrimaryMomentum1stImpact[i].Z() << endl;

  for (UShort_t i = 0; i < fHistoryCurrentCreatorProcess.size(); i++)
  	 cout << "Current Creator Proc : " << fHistoryCurrentCreatorProcess[i] << endl;
  	
   for (UShort_t i = 0; i < fHistoryParentID.size(); i++)
      cout << "Parent ID: " << fHistoryParentID[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryCurrentID.size(); i++)
      cout << "Current ID: " << fHistoryCurrentID[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryCurrentPdg.size(); i++)
      cout << "Current Pdg: " << fHistoryCurrentPdg[i] << endl;
     
   for (UShort_t i = 0; i < fHistoryCurrentEnergy.size(); i++)
  	 cout << "Current Energy: " << fHistoryCurrentEnergy[i] << endl;
   
   for (UShort_t i = 0; i < fHistoryCurrentTime.size(); i++)
  	 cout << "Current Energy: " << fHistoryCurrentTime[i] << endl;

	for (UShort_t i = 0; i < fHistoryCurrentVolume1stImpact.size(); i++)
		cout << " History Current volume 1st Impact : " << fHistoryCurrentVolume1stImpact[i] << endl;
		  	    

	for (UShort_t i = 0; i < fHistoryCurrentPositionVertex.size(); i++)
		cout << " History Current Position Vertex    : " << fHistoryCurrentPositionVertex[i].X() << " " << fHistoryCurrentPositionVertex[i].Y() << " " << fHistoryCurrentPositionVertex[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentPosition1stImpact.size(); i++)
		cout << " History Current Position 1st Impact : " << fHistoryCurrentPosition1stImpact[i].X() << " " << fHistoryCurrentPosition1stImpact[i].Y() << " " << fHistoryCurrentPosition1stImpact[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentPositionDetector.size(); i++)
		cout << " History Current Position Detector  : " << fHistoryCurrentPositionDetector[i].X() << " " << fHistoryCurrentPositionDetector[i].Y() << " " << fHistoryCurrentPositionDetector[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentPositionDeath.size(); i++)
		cout << " History Current Position Death     : " << fHistoryCurrentPositionDeath[i].X() << " " << fHistoryCurrentPositionDeath[i].Y() << " " << fHistoryCurrentPositionDeath[i].Z() << endl;
 
	for (UShort_t i = 0; i < fHistoryCurrentMomentumVertex.size(); i++)
		cout << " History Current Momentum Vertex    : " << fHistoryCurrentMomentumVertex[i].X() << " " << fHistoryCurrentMomentumVertex[i].Y() << " " << fHistoryCurrentMomentumVertex[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentMomentum1stImpact.size(); i++)
		cout << " History Current Momentum 1st Impact : " << fHistoryCurrentMomentum1stImpact[i].X() << " " << fHistoryCurrentMomentum1stImpact[i].Y() << " " << fHistoryCurrentMomentum1stImpact[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentMomentumDetector.size(); i++)
		cout << " History Current Momentum Detector  : " << fHistoryCurrentMomentumDetector[i].X() << " " << fHistoryCurrentMomentumDetector[i].Y() << " " << fHistoryCurrentMomentumDetector[i].Z() << endl;
	for (UShort_t i = 0; i < fHistoryCurrentMomentumDeath.size(); i++)
		cout << " History Current Momentum Death     : " << fHistoryCurrentMomentumDeath[i].X() << " " << fHistoryCurrentMomentumDeath[i].Y() << " " << fHistoryCurrentMomentumDeath[i].Z() << endl;
  
  	 
  	 
   
}
