/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Mhd Moukaddam  contact address: moukaddam@triumf.ca      *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             Spice detector. And was adapted from S1 detector Class.        *
 *             in th NPTOOL project                                          *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "TSpiceData.h"


ClassImp(TSpiceData)

TSpiceData::TSpiceData()
{
   // Default constructor
   Clear();
}

TSpiceData::~TSpiceData() {}


void TSpiceData::Clear()
{
   // DSSD
   // (Th,E)
   fSpice_Theta_DetNbr.clear();
   fSpice_Theta_StripNbr.clear();
   fSpice_Theta_Energy.clear();
   fSpice_Theta_ResEnergy.clear();

   // (Ph,E)
   fSpice_Phi_DetNbr.clear();
   fSpice_Phi_StripNbr.clear();
   fSpice_Phi_Energy.clear();
   fSpice_Phi_ResEnergy.clear();

   fPrimaryTheta.clear(); 
   fPrimaryPhi.clear();  
   fPrimaryEnergy.clear();   
   fPrimaryPdg.clear();
   fPdg.clear();
   fPositionFirstHit.clear();
   fEventNumber = -1 ;
}


void TSpiceData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event : " << fEventNumber <<  "XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

   // DSSD
   // (Th,E)
   cout << "Spice_MultThE = " << fSpice_Theta_DetNbr.size() << endl;
   for (UShort_t i = 0; i < fSpice_Theta_DetNbr.size(); i++)
      cout << "DetThE: " << fSpice_Theta_DetNbr[i] << " StripThE: " << fSpice_Theta_StripNbr[i] 
      		 << " EnergyTh: " << fSpice_Theta_Energy[i] << " ResEnergyTh: " << fSpice_Theta_ResEnergy[i] << endl;

   // (Ph,E)
   cout << "Spice_MultPhE = " << fSpice_Phi_DetNbr.size() << endl;
   for (UShort_t i = 0; i < fSpice_Phi_DetNbr.size(); i++)
      cout << "DetPhE: " << fSpice_Phi_DetNbr[i] << " StripPhE: " << fSpice_Phi_StripNbr[i] 
      		 << " EnergyPh: " << fSpice_Phi_Energy[i] << " ResEnergyPh: " << fSpice_Phi_ResEnergy[i] << endl;

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
