/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Ryan Dunlop  contact address: rdunlop@uoguelph.ca        *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             Griffin detector. And was adapted from S1 detector Class.     *
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
        fDetNbr.clear();
        fPrimaryEnergy.clear();
        fPrimaryPdg.clear();
        fPdg.clear();
        fPositionFirstHit.clear();
        fEventNumber = -1;
}


void TGriffinData::Dump() const
{
        cout << "Stuff will go here" << endl;
   
}
