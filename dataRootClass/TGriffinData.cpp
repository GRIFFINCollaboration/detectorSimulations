/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Mhd Moukaddam  contact address: moukaddam@triumf.ca      *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             Griffin detector. And was adapted from S1 detector Class.        *
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

#include "TGriffinData.h"


ClassImp(TGriffinData)

TGriffinData::TGriffinData()
{
   // Default constructor
   Clear();
}

TGriffinData::~TGriffinData() {}


void TGriffinData::Clear()
{
   // DSSD
   // (Th,E)
   fGriffin_Theta_DetNbr.clear();
   fGriffin_Theta_StripNbr.clear();
   fGriffin_Theta_Energy.clear();

   // (Ph,E)
   fGriffin_Phi_DetNbr.clear();
   fGriffin_Phi_StripNbr.clear();
   fGriffin_Phi_Energy.clear();

   fPrimaryTheta.clear(); 
   fPrimaryPhi.clear();  
   fPrimaryEnergy.clear();   
   fPrimaryPdg.clear();
   fPdg.clear();
   fPositionFirstHit.clear();
   fEventNumber = -1 ;
}


void TGriffinData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event : " << fEventNumber <<  "XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

   // DSSD
   // (Th,E)
   cout << "Griffin_MultThE = " << fGriffin_Theta_DetNbr.size() << endl;
   for (UShort_t i = 0; i < fGriffin_Theta_DetNbr.size(); i++)
      cout << "DetThE: " << fGriffin_Theta_DetNbr[i] << " StripThE: " << fGriffin_Theta_StripNbr[i] << " EnergyTh: " << fGriffin_Theta_Energy[i] << endl;

   // (Ph,E)
   cout << "Griffin_MultPhE = " << fGriffin_Phi_DetNbr.size() << endl;
   for (UShort_t i = 0; i < fGriffin_Phi_DetNbr.size(); i++)
      cout << "DetPhE: " << fGriffin_Phi_DetNbr[i] << " StripPhE: " << fGriffin_Phi_StripNbr[i] << " EnergyPh: " << fGriffin_Phi_Energy[i] << endl;

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
