/*****************************************************************************
 * Original Author    : Mhd Moukaddam  contact address: moukaddam@triumf.ca  *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             PACESdetector. And was adapted from S1 detector Class.        *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "TPacesData.h"


ClassImp(TPacesData)

TPacesData::TPacesData() {
   Clear();
}

TPacesData::~TPacesData() {
}


void TPacesData::Clear() {
 
   fPaces_DetNbr.clear();
   fPaces_Energy.clear();

   fPrimaryTheta.clear(); 
   fPrimaryPhi.clear();  
   fPrimaryEnergy.clear();   
   fPrimaryPdg.clear();
   fPdg.clear();
   fPositionFirstHit.clear();
   fEventNumber = -1 ;
}


void TPacesData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event : " << fEventNumber <<  " XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

   // DSSD
   // (Th,E)
   cout << "Paces_MultE = " << fPaces_DetNbr.size() << endl;
   for (UShort_t i = 0; i < fPaces_DetNbr.size(); i++)
      cout << "DetNbr: " << fPaces_DetNbr[i] << endl;
   for (UShort_t i = 0; i < fPaces_Energy.size(); i++)
      cout << "Energy: " << fPaces_Energy[i] << endl;
      
   // First hit
   cout << "Position of First Hit Mult = " << fPositionFirstHit.size() << endl;
   for (UShort_t i = 0; i < fPositionFirstHit.size(); i++)   
   cout << "    X - " << fPositionFirstHit.at(i).X() << "  Y - " <<fPositionFirstHit.at(i).Y() << " Z - " << fPositionFirstHit.at(i).Z() << endl;
   
   // Primary angles
   cout << "Primary Particle angle Mult = " << fPrimaryTheta.size() << endl;
   for (UShort_t i = 0; i < fPrimaryTheta.size(); i++)   
   cout << "    Theta - " << fPrimaryTheta.at(i) << "  Phi - " <<fPrimaryPhi.at(i) << endl;
   // Primary Pdg
   cout << "Primary Particle Pdg Mult = " << fPrimaryPdg.size() << endl;
   for (UShort_t i = 0; i < fPrimaryPdg.size(); i++)   
   cout << "   Pdg - " << fPrimaryPdg.at(i) << endl;
   // Primary energy
   cout << "Primary Particle Energy Mult = " << fPrimaryEnergy.size() << endl;
   for (UShort_t i = 0; i < fPrimaryEnergy.size(); i++)   
   cout << "   Energy - " << fPrimaryEnergy.at(i) << endl;
   
   cout << "Particle Pdg Mult = " << fPdg.size() << endl;
   for (UShort_t i = 0; i < fPdg.size(); i++)   
   cout << "   Pdg - " << fPdg.at(i) << endl;
   
}
