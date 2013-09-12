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
// $Id: DetectorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:Detector(Det)
{ 
  DetSysDir = new G4UIdirectory("/DetSys/");
  DetSysDir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/DetSys/det/");
  detDir->SetGuidance("detector control");

  appDir = new G4UIdirectory("/DetSys/app/");
  appDir->SetGuidance("apparatus control");

  worldDir = new G4UIdirectory("/DetSys/world/");
  worldDir->SetGuidance("world control");

  WorldMaterialCmd = new G4UIcmdWithAString("/DetSys/world/material",this);
  WorldMaterialCmd->SetGuidance("Select material for the world.");
  WorldMaterialCmd->SetParameterName("choice",false);
  WorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/world/dimensions",this);
  WorldDimensionsCmd->SetGuidance("Set world dimensions - x y z unit.");
  WorldDimensionsCmd->SetUnitCategory("Length");
  WorldDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldVisCmd = new G4UIcmdWithABool("/DetSys/world/vis",this);
  WorldVisCmd->SetGuidance("Set the visulization of the world");
  WorldVisCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldMagneticFieldCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/world/magneticField",this);
  WorldMagneticFieldCmd->SetGuidance("Set world magnetic field - x y z unit.");
  WorldMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
  WorldMagneticFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  GenericTargetCmd = new G4UIcmdWithAString("/DetSys/app/genericTarget",this);
  GenericTargetCmd->SetGuidance("Select material of the target.");
  GenericTargetCmd->SetParameterName("choice",false);
  GenericTargetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  GenericTargetDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetDimensions",this);
  GenericTargetDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
  GenericTargetDimensionsCmd->SetUnitCategory("Length");
  GenericTargetDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  GenericTargetPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/genericTargetPosition",this);
  GenericTargetPositionCmd->SetGuidance("Set target position - x y z unit.");
  GenericTargetPositionCmd->SetUnitCategory("Length");
  GenericTargetPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FieldBoxMaterialCmd = new G4UIcmdWithAString("/DetSys/app/fieldBoxMaterial",this);
  FieldBoxMaterialCmd->SetGuidance("Select material of the target.");
  FieldBoxMaterialCmd->SetParameterName("choice",false);
  FieldBoxMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FieldBoxDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxDimensions",this);
  FieldBoxDimensionsCmd->SetGuidance("Set target dimensions - x y z unit.");
  FieldBoxDimensionsCmd->SetUnitCategory("Length");
  FieldBoxDimensionsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FieldBoxPositionCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxPosition",this);
  FieldBoxPositionCmd->SetGuidance("Set target position - x y z unit.");
  FieldBoxPositionCmd->SetUnitCategory("Length");
  FieldBoxPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FieldBoxMagneticFieldCmd = new G4UIcmdWith3VectorAndUnit("/DetSys/app/fieldBoxMagneticField",this);
  FieldBoxMagneticFieldCmd->SetGuidance("Set magnetic field - x y z unit.");
  FieldBoxMagneticFieldCmd->SetUnitCategory("Magnetic flux density");
  FieldBoxMagneticFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddApparatusSpiceTargetChamberCmd = new G4UIcmdWithoutParameter("/DetSys/app/addSpiceTargetChamber",this);
  AddApparatusSpiceTargetChamberCmd->SetGuidance("Add SPICE target chamber.");
  AddApparatusSpiceTargetChamberCmd->AvailableForStates(G4State_Idle);

  AddApparatus8piVacuumChamberCmd = new G4UIcmdWithoutParameter("/DetSys/app/add8piVacuumChamber",this);
  AddApparatus8piVacuumChamberCmd->SetGuidance("Add 8pi vacuum chamber.");
  AddApparatus8piVacuumChamberCmd->AvailableForStates(G4State_Idle);

  AddApparatus8piVacuumChamberAuxMatShellCmd = new G4UIcmdWithAnInteger("/DetSys/app/add8piVacuumChamberAuxMatShell",this);
  AddApparatus8piVacuumChamberAuxMatShellCmd->SetGuidance("Add AuxMat shell around 8pi vacuum chamber");
  AddApparatus8piVacuumChamberAuxMatShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/DetSys/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  AddDetectionSystemGammaTrackingCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGammaTracking",this);
  AddDetectionSystemGammaTrackingCmd->SetGuidance("Add Detection System GammaTracking");
  AddDetectionSystemGammaTrackingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemBrillance380V1Cmd = new G4UIcmdWithAnInteger("/DetSys/det/addBrillance380V1",this);
  AddDetectionSystemBrillance380V1Cmd->SetGuidance("Add Detection System Brillance380V1");
  AddDetectionSystemBrillance380V1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemSodiumIodideCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSodiumIodide",this);
  AddDetectionSystemSodiumIodideCmd->SetGuidance("Add Detection System SodiumIodide");
  AddDetectionSystemSodiumIodideCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystem8piCmd = new G4UIcmdWithAnInteger("/DetSys/det/add8pi",this);
  AddDetectionSystem8piCmd->SetGuidance("Add Detection System 8pi");
  AddDetectionSystem8piCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystem8piDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/add8piDetector",this);
  AddDetectionSystem8piDetectorCmd->SetGuidance("Add 8pi Detector");
  AddDetectionSystem8piDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemGriffinForwardCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForward",this);
  AddDetectionSystemGriffinForwardCmd->SetGuidance("Add Detection System GriffinForward");
  AddDetectionSystemGriffinForwardCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemGriffinForwardDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinForwardDetector",this);
  AddDetectionSystemGriffinForwardDetectorCmd->SetGuidance("Add GriffinForward Detector");
  AddDetectionSystemGriffinForwardDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemGriffinBackCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBack",this);
  AddDetectionSystemGriffinBackCmd->SetGuidance("Add Detection System GriffinBack");
  AddDetectionSystemGriffinBackCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemGriffinBackDetectorCmd = new G4UIcmdWithAnInteger("/DetSys/det/addGriffinBackDetector",this);
  AddDetectionSystemGriffinBackDetectorCmd->SetGuidance("Add GriffinBack Detector");
  AddDetectionSystemGriffinBackDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemSceptarCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSceptar",this);
  AddDetectionSystemSceptarCmd->SetGuidance("Add Detection System Sceptar");
  AddDetectionSystemSceptarCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemSpiceCmd = new G4UIcmdWithAnInteger("/DetSys/det/addSpice",this);
  AddDetectionSystemSpiceCmd->SetGuidance("Add Detection System Spice");
  AddDetectionSystemSpiceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemSpiceV02Cmd = new G4UIcmdWithAnInteger("/DetSys/det/addSpiceV02",this);
  AddDetectionSystemSpiceV02Cmd->SetGuidance("Add Detection System SpiceV02");
  AddDetectionSystemSpiceV02Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddDetectionSystemPacesCmd = new G4UIcmdWithAnInteger("/DetSys/det/addPaces",this);
  AddDetectionSystemPacesCmd->SetGuidance("Add Detection System Paces");
  AddDetectionSystemPacesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete detDir;
  delete appDir;
  delete worldDir;
  delete WorldMaterialCmd;
  delete WorldDimensionsCmd;
  delete WorldVisCmd;
  delete WorldMagneticFieldCmd;
  delete DetSysDir; 
  delete UpdateCmd;
  delete GenericTargetCmd;
  delete GenericTargetDimensionsCmd;
  delete GenericTargetPositionCmd;
  delete FieldBoxMaterialCmd;
  delete FieldBoxDimensionsCmd;
  delete FieldBoxPositionCmd;
  delete FieldBoxMagneticFieldCmd;
  delete AddApparatusSpiceTargetChamberCmd;
  delete AddDetectionSystemGammaTrackingCmd;
  delete AddApparatus8piVacuumChamberCmd;
  delete AddApparatus8piVacuumChamberAuxMatShellCmd;
  delete AddDetectionSystemBrillance380V1Cmd;
  delete AddDetectionSystemSodiumIodideCmd;
  delete AddDetectionSystem8piCmd;
  delete AddDetectionSystem8piDetectorCmd;
  delete AddDetectionSystemSceptarCmd;
  delete AddDetectionSystemGriffinForwardCmd;
  delete AddDetectionSystemGriffinForwardDetectorCmd;
  delete AddDetectionSystemGriffinBackCmd;
  delete AddDetectionSystemGriffinBackDetectorCmd;
  delete AddDetectionSystemSpiceCmd;
  delete AddDetectionSystemSpiceV02Cmd;
  delete AddDetectionSystemPacesCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == WorldMaterialCmd ) { 
    Detector->SetWorldMaterial(newValue);
  } 
  if( command == WorldDimensionsCmd ) { 
    Detector->SetWorldDimensions(WorldDimensionsCmd->GetNew3VectorValue(newValue));
  } 
  if( command == WorldVisCmd ) { 
    Detector->SetWorldVis(WorldVisCmd->GetNewBoolValue(newValue)); 
  }
  if( command == WorldMagneticFieldCmd ) {
    Detector->SetWorldMagneticField(WorldMagneticFieldCmd->GetNew3VectorValue(newValue));
  }
  if( command == UpdateCmd ) { 
    Detector->UpdateGeometry(); 
  }
  if( command == GenericTargetCmd ) { 
    Detector->SetGenericTargetMaterial(newValue);
  }
  if( command == GenericTargetDimensionsCmd ) { 
    Detector->SetGenericTargetDimensions(GenericTargetDimensionsCmd->GetNew3VectorValue(newValue));
  }
  if( command == GenericTargetPositionCmd ) { 
    Detector->SetGenericTargetPosition(GenericTargetPositionCmd->GetNew3VectorValue(newValue));
  }
  if( command == FieldBoxMaterialCmd ) {
    Detector->SetFieldBoxMaterial(newValue);
  }
  if( command == FieldBoxDimensionsCmd ) {
    Detector->SetFieldBoxDimensions(FieldBoxDimensionsCmd->GetNew3VectorValue(newValue));
  }
  if( command == FieldBoxPositionCmd ) {
    Detector->SetFieldBoxPosition(FieldBoxPositionCmd->GetNew3VectorValue(newValue));
  }
  if( command == FieldBoxMagneticFieldCmd ) {
    Detector->SetFieldBoxMagneticField(FieldBoxMagneticFieldCmd->GetNew3VectorValue(newValue));
  }
  if( command == AddApparatusSpiceTargetChamberCmd ) { 
    Detector->AddApparatusSpiceTargetChamber(); 
  }
  if( command == AddApparatus8piVacuumChamberCmd ) {
    Detector->AddApparatus8piVacuumChamber();
  }
  if( command == AddApparatus8piVacuumChamberAuxMatShellCmd ) {
    Detector->AddApparatus8piVacuumChamberAuxMatShell(AddApparatus8piVacuumChamberAuxMatShellCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemGammaTrackingCmd ) {
    Detector->AddDetectionSystemGammaTracking(AddDetectionSystemGammaTrackingCmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystemBrillance380V1Cmd ) { 
    Detector->AddDetectionSystemBrillance380V1(AddDetectionSystemBrillance380V1Cmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystemSodiumIodideCmd ) {
    Detector->AddDetectionSystemSodiumIodide(AddDetectionSystemSodiumIodideCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystem8piCmd ) { 
    Detector->AddDetectionSystem8pi(AddDetectionSystem8piCmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystem8piDetectorCmd ) { 
    Detector->AddDetectionSystem8piDetector(AddDetectionSystem8piDetectorCmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystemSceptarCmd ) { 
    Detector->AddDetectionSystemSceptar(AddDetectionSystemSceptarCmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystemGriffinForwardCmd ) {
    Detector->AddDetectionSystemGriffinForward(AddDetectionSystemGriffinForwardCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemGriffinForwardDetectorCmd ) {
    Detector->AddDetectionSystemGriffinForwardDetector(AddDetectionSystemGriffinForwardDetectorCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemGriffinBackCmd ) {
    Detector->AddDetectionSystemGriffinBack(AddDetectionSystemGriffinBackCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemGriffinBackDetectorCmd ) {
    Detector->AddDetectionSystemGriffinBackDetector(AddDetectionSystemGriffinBackDetectorCmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemSpiceCmd ) { 
    Detector->AddDetectionSystemSpice(AddDetectionSystemSpiceCmd->GetNewIntValue(newValue)); 
  }
  if( command == AddDetectionSystemSpiceV02Cmd ) {
    Detector->AddDetectionSystemSpiceV02(AddDetectionSystemSpiceV02Cmd->GetNewIntValue(newValue));
  }
  if( command == AddDetectionSystemPacesCmd ) {
    Detector->AddDetectionSystemPaces(AddDetectionSystemPacesCmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
