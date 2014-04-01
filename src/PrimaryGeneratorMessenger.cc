//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: PrimaryGeneratorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
:Action(Gun)
{
  gunDir = new G4UIdirectory("/DetSys/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  energyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/energy",this);
  energyCmd->SetGuidance("Sets energy of the photons.");
  energyCmd->SetUnitCategory("Energy");
  energyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  particleCmd = new G4UIcmdWithAString("/DetSys/gun/particle",this);
  particleCmd->SetGuidance("Set particle type.");
  particleCmd->SetParameterName("choice",false);
  particleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ionCmd = new G4UIcommand("/DetSys/gun/ion",this);
  ionCmd->SetGuidance("Set ion type - atomic number, mass number, excitation energy (units in keV).");
  G4UIparameter *parameter1,*parameter2, *parameter3;
  G4bool omitable;
  parameter1 = new G4UIparameter ("Z", 'i', omitable = false);
  ionCmd->SetParameter(parameter1);
  parameter2 = new G4UIparameter ("A", 'i', omitable = false);
  ionCmd->SetParameter(parameter2);
  parameter3 = new G4UIparameter ("E", 'd', omitable = false);
  ionCmd->SetParameter(parameter3);
  ionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  betaPlusEmissionCmd = new G4UIcmdWithAString("/DetSys/gun/betaPlusEmission",this);
  betaPlusEmissionCmd->SetGuidance("Simulate a beta plus decay with an energy distribution input file");
  betaPlusEmissionCmd->SetParameterName("choice",false);
  betaPlusEmissionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  betaMinusEmissionCmd = new G4UIcmdWithAString("/DetSys/gun/betaMinusEmission",this);
  betaMinusEmissionCmd->SetGuidance("Simulate a beta negative decay with an energy distribution input file");
  betaMinusEmissionCmd->SetParameterName("choice",false);
  betaMinusEmissionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  radioactiveBetaDecayCmd = new G4UIcmdWithAString("/DetSys/gun/radioactiveBetaDecay",this);
  radioactiveBetaDecayCmd->SetGuidance("Simulate a complete radioactive beta negative decay with a simulation directory");
  radioactiveBetaDecayCmd->SetParameterName("choice",false);
  radioactiveBetaDecayCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  radioactiveSourceDecayCmd = new G4UIcmdWithAString("/DetSys/gun/radioactiveSourceDecay",this);
  radioactiveSourceDecayCmd->SetGuidance("Simulate a radioactive source decay with a file containing the decay data");
  radioactiveSourceDecayCmd->SetParameterName("choice",false);
  radioactiveSourceDecayCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  polarizationCmd = new G4UIcmdWithADouble("/DetSys/gun/polarization",this);
  polarizationCmd->SetGuidance("Set Polarization of Nuclei");
  polarizationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  emitBetaParticleCmd = new G4UIcmdWithABool("/DetSys/gun/emitBetaParticle",this);
  emitBetaParticleCmd->SetGuidance("Emit Beta Particle? True/False");
  emitBetaParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  includeXRayInputFileKShellCmd = new G4UIcmdWithABool("/DetSys/gun/includeXRayInputFileKShell",this);
  includeXRayInputFileKShellCmd->SetGuidance("Emit X-rays from K-shell vacancies using input file? True/False");
  includeXRayInputFileKShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  includeXRayInputFileLShellCmd = new G4UIcmdWithABool("/DetSys/gun/includeXRayInputFileLShell",this);
  includeXRayInputFileLShellCmd->SetGuidance("Emit X-rays from L-shell vacancies using input file? True/False");
  includeXRayInputFileLShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  includeXRayInputFileMShellCmd = new G4UIcmdWithABool("/DetSys/gun/includeXRayInputFileMShell",this);
  includeXRayInputFileMShellCmd->SetGuidance("Emit X-rays from M-shell vacancies using input file? True/False");
  includeXRayInputFileMShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  radioactiveDecayHalflifeCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/radioactiveDecayHalflife",this);
  radioactiveDecayHalflifeCmd->SetGuidance("Half-life of radioactive isotope simulation");
  radioactiveDecayHalflifeCmd->SetUnitCategory("Time");
  radioactiveDecayHalflifeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  numberOfRadioactiveNucleiCmd = new G4UIcmdWithAnInteger("/DetSys/gun/numberOfRadioactiveNuclei",this);
  numberOfRadioactiveNucleiCmd->SetGuidance("Set the number of radioactive nuclei");
  numberOfRadioactiveNucleiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  directionCmd = new G4UIcmdWith3Vector("/DetSys/gun/direction",this);
  directionCmd->SetGuidance("Set momentum direction.");
  directionCmd->SetGuidance("Direction needs not to be a unit vector.");
  directionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  positionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/gun/position",this);
  positionCmd->SetGuidance("Set particle position.");
  positionCmd->SetUnitCategory("Length");
  positionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  radiusCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/radius",this);
  radiusCmd->SetGuidance("Set source active radius.");
  radiusCmd->SetUnitCategory("Length");
  radiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  rangeCmd = new G4UIcommand("/DetSys/gun/energyrange",this);
  rangeCmd->SetGuidance("Set energy range - minimum, maximum, step size (all in keV)");
  G4UIparameter *parameter4, *parameter5, *parameter6;
	parameter4 = new G4UIparameter ("minimum", 'd', omitable = false);
  rangeCmd->SetParameter(parameter4);
  parameter5 = new G4UIparameter ("maximum", 'd', omitable = false);
  rangeCmd->SetParameter(parameter5);
  parameter6 = new G4UIparameter ("step", 'd', omitable = false);
  rangeCmd->SetParameter(parameter6);
  rangeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  simKinematicsCmd = new G4UIcmdWithABool("/DetSys/gun/simulateKinematics",this);
  simKinematicsCmd->SetGuidance("Choose to simulate kinematic shift/broadening of particles");
  simKinematicsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  simKinematicsBetaValueCmd = new G4UIcmdWithADouble("/DetSys/gun/kinematicsBetaValue",this);
  simKinematicsBetaValueCmd->SetGuidance("Set beta value of heavy ion");
  simKinematicsBetaValueCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  simKinematicsIonEnergyCmd = new G4UIcmdWithADoubleAndUnit("/DetSys/gun/kinematicsIonEnergy",this);
	simKinematicsIonEnergyCmd->SetGuidance("Sets energy of the ions in kinematic simulation.");
  simKinematicsIonEnergyCmd->SetUnitCategory("Energy");
  simKinematicsIonEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete energyCmd;
  delete directionCmd;
  delete positionCmd;
  delete radiusCmd;
  delete rangeCmd;
  delete particleCmd;
  delete ionCmd;
  delete gunDir;
  delete betaPlusEmissionCmd;
  delete betaMinusEmissionCmd;
  delete radioactiveBetaDecayCmd;
  delete radioactiveSourceDecayCmd;
  delete polarizationCmd;
  delete emitBetaParticleCmd;
  delete includeXRayInputFileKShellCmd;
  delete includeXRayInputFileLShellCmd;
  delete includeXRayInputFileMShellCmd;
  delete radioactiveDecayHalflifeCmd;
  delete numberOfRadioactiveNucleiCmd;
  delete simKinematicsCmd;
  delete simKinematicsBetaValueCmd;
  delete simKinematicsIonEnergyCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue( G4UIcommand* command, G4String newValue )
{ 
  if( command == energyCmd ) { 
    Action->SetEnergy(energyCmd->GetNewDoubleValue(newValue));
  }
  if( command == particleCmd ) { 
    Action->SetParticleType(newValue);
  }
  if( command == ionCmd ) { 
    G4int Z,A;
    G4double E;
    const char* s = newValue;
    std::istringstream is ((char*)s);
    is>>Z>>A>>E;
    Action->SetIonType(Z,A,E);
  }
  if( command == betaPlusEmissionCmd ) { 
    Action->SetBetaPlusEmission(newValue);
  }
  if( command == betaMinusEmissionCmd ) { 
    Action->SetBetaMinusEmission(newValue);
  }
  if( command == radioactiveSourceDecayCmd ) {
    Action->SetRadioactiveSourceDecay(newValue);
  }
  if( command == radioactiveBetaDecayCmd ) {
    Action->SetRadioactiveBetaDecay(newValue);
  }
  if( command == polarizationCmd ) {
    Action->SetPolarization(polarizationCmd->GetNewDoubleValue(newValue));
  }
  if( command == emitBetaParticleCmd ) {
    Action->SetEmitBetaParticle(emitBetaParticleCmd->GetNewBoolValue(newValue));
  }
  if( command == includeXRayInputFileKShellCmd ) {
    Action->IncludeXRayInputFileKShell(includeXRayInputFileKShellCmd->GetNewBoolValue(newValue));
  }
  if( command == includeXRayInputFileLShellCmd ) {
    Action->IncludeXRayInputFileLShell(includeXRayInputFileLShellCmd->GetNewBoolValue(newValue));
  }
  if( command == includeXRayInputFileMShellCmd ) {
    Action->IncludeXRayInputFileMShell(includeXRayInputFileMShellCmd->GetNewBoolValue(newValue));
  }
  if( command == radioactiveDecayHalflifeCmd ) {
    Action->SetRadioactiveDecayHalflife(radioactiveDecayHalflifeCmd->GetNewDoubleValue(newValue));
  }
  if( command == numberOfRadioactiveNucleiCmd ) {
    Action->SetNumberOfRadioactiveNuclei(numberOfRadioactiveNucleiCmd->GetNewIntValue(newValue));
  }
  if( command == directionCmd ) { 
    Action->SetDirection(directionCmd->GetNew3VectorValue(newValue));
  }
  if( command == positionCmd ) { 
    Action->SetPosition(positionCmd->GetNew3VectorValue(newValue));
  }
	if( command == radiusCmd ) { 
    Action->SetSourceRadius(radiusCmd->GetNewDoubleValue(newValue));
  }
  if( command == rangeCmd ) { 
    G4double minimum, maximum, step;
    const char* s = newValue;
    std::istringstream is ((char*)s);
    is>>minimum>>maximum>>step;
    Action->SetEnergyRange(minimum,maximum,step);
  }
  if( command == simKinematicsCmd) {
  	Action->SetKinematicsActive(simKinematicsCmd->GetNewBoolValue(newValue));
  }
  if( command == simKinematicsBetaValueCmd ) {
  	Action->SetKinematicsBetaValue(simKinematicsBetaValueCmd->GetNewDoubleValue(newValue));
  }
  if( command == simKinematicsIonEnergyCmd ) { 
    Action->SetKinematicsIonEnergy(simKinematicsIonEnergyCmd->GetNewDoubleValue(newValue));
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

