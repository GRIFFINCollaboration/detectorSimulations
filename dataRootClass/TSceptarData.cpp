

#include <iostream>
using namespace std;

#include "TSceptarData.h"


ClassImp(TSceptarData)

TSceptarData::TSceptarData() 
{
    ClearVariables();
}

TSceptarData::~TSceptarData()
{
}

void TSceptarData::ClearVariables()
{
    fSceptar_DetNbr.clear();
    fSceptar_Energy.clear();

    fPrimaryTheta.clear();
    fPrimaryPhi.clear();
    fPrimaryEnergy.clear();
    fPrimaryPdg.clear();
    fPdg.clear();
    fPositionFirstHit.clear();
    fEventNumber = -1 ;
}

void TSceptarData::Dump() const
{
    
    cout << "XXXXXXXXXXXXXXXXXXXXXXXX Sceptar Event : " << fEventNumber << "XXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    // DSSD
    // (Th,E)
    cout << "Sceptar_MultE = " << fSceptar_DetNbr.size() << endl;
    for (UShort_t i = 0; i < fSceptar_DetNbr.size(); i++)
    cout << "DetNbr: " << fSceptar_DetNbr[i] << "Energy: " << fSceptar_Energy[i] << endl;

    // First hit
    cout << "Position of First Hit Mult = " << fPositionFirstHit.size() << endl;
    for (UShort_t i = 0; i < fPositionFirstHit.size(); i++)
    cout << " X - " << fPositionFirstHit.at(i).X() << " Y - " <<fPositionFirstHit.at(i).Y() << " Z - " << fPositionFirstHit.at(i).Z() << endl;

    // Primary angles
    cout << "Primary Particle angle Mult = " << fPrimaryTheta.size() << endl;
    for (UShort_t i = 0; i < fPrimaryTheta.size(); i++)
    cout << " Theta - " << fPrimaryTheta.at(i) << " Phi - " <<fPrimaryPhi.at(i) << endl;
    // Primary Pdg
    cout << "Primary Particle Pdg Mult = " << fPrimaryPdg.size() << endl;
    for (UShort_t i = 0; i < fPrimaryPdg.size(); i++)
    cout << " Pdg - " << fPrimaryPdg.at(i) << endl;
    // Primary energy
    cout << "Primary Particle Energy Mult = " << fPrimaryEnergy.size() << endl;
    for (UShort_t i = 0; i < fPrimaryEnergy.size(); i++)
    cout << " Energy - " << fPrimaryEnergy.at(i) << endl;

    cout << "Particle Pdg Mult = " << fPdg.size() << endl;
    for (UShort_t i = 0; i < fPdg.size(); i++)
    cout << " Pdg - " << fPdg.at(i) << endl;

}
