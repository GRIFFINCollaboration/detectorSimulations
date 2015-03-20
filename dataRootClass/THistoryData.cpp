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
	fHistoryPrimaryEnergy.clear(); 
	fHistoryPrimaryPdg.clear(); 
	fHistoryCurrentPdg.clear(); 
	fHistoryGNumber.clear(); 
	fHistoryG2BirthVolume.clear(); 
	fHistoryGLastBirthVolume.clear(); 
	fHistoryG2Process.clear(); 
}


void THistoryData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX HISTORY XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

   cout << "History_Mult = " <<fHistoryPrimaryEnergy.size() << endl;
   
   for (UShort_t i = 0; i < fHistoryPrimaryEnergy.size(); i++)
      cout << "Primary Energy: " << fHistoryPrimaryEnergy[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryPrimaryPdg.size(); i++)
      cout << "Primary Pdg: " << fHistoryPrimaryPdg[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryCurrentPdg.size(); i++)
      cout << "Current Pdg: " << fHistoryCurrentPdg[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryGNumber.size(); i++)
      cout << "Number of Generations : " << fHistoryGNumber[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryG2BirthVolume.size(); i++)
      cout << "Birth volume at G-2: " << fHistoryG2BirthVolume[i] << endl;
      
   for (UShort_t i = 0; i < fHistoryGLastBirthVolume.size(); i++)
      cout << "Genrating Process at G-2 : " << fHistoryGLastBirthVolume[i] << endl;
     
   for (UShort_t i = 0; i < fHistoryG2Process.size(); i++)
  	 cout << "Birth volume at G-last : " << fHistoryG2Process[i] << endl;
   
   
}
