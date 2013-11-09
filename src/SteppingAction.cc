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
// $Id: SteppingAction.cc,v 1.1 2010-11-08 10:38:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                                         EventAction* evt)
:detector(det), eventaction(evt)					 
{
    griffinDetectorMapSet = false;
    numberOfAssemblyVols = 13;
    
  // List of assembly volumes just to keep track:
  // 
  // assembly 
  // leftSuppressorCasingAssembly
  // rightSuppressorCasingAssembly
  // leftSuppressorExtensionAssembly
  // rightSuppressorExtensionAssembly
  // suppressorBackAssembly 
  // extensionSuppressorShellAssembly   
  // backAndSideSuppressorShellAssembly
  // hevimetAssembly 
  // germaniumAssemblyCry[0]
  // germaniumAssemblyCry[1]
  // germaniumAssemblyCry[2]
  // germaniumAssemblyCry[3]

    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4int particleType = 0;
  G4int volNameOver9;

  // get volume of the current step
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4String volname = volume->GetName();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();

  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();

  // Track particle type in EVERY step
  //G4cout << "Particle name = " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << G4endl;
  if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma")         particleType = 1;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-")       particleType = 2;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e+")       particleType = 3;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "proton")   particleType = 4;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "neutron")  particleType = 5;
  else                                                                                  particleType = 0;

  eventaction->AddParticleType(particleType);

  //G4cout << "Found Edep = " << edep/keV << " keV in " << volname << G4endl;
  // example volname
  //volname = av_1_impr_6_sodium_iodide_crystal_block_log_pv_0

  size_t found;
  G4String search;
  G4int searchLength;

  // Grid Cell
  found = volname.find("gridcell");
  if (ekin != 0 && found!=G4String::npos && particleType == 1) {
      SetDetNumberForGenericDetector(volname);
      eventaction->SetGridEKinGammaDet(ekin,stepl,det-1);
  }

  found = volname.find("gridcell");
  if (ekin != 0 && found!=G4String::npos && particleType == 2) {
      SetDetNumberForGenericDetector(volname);
      eventaction->SetGridEKinElectronDet(ekin,stepl,det-1);
  }

  // Griffin energy deposits
  found = volname.find("germanium_block1");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("back_quarter_suppressor");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorBackDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("left_suppressor_extension");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      //G4cout << "left_suppressor_extension Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinSuppressorLeftExtensionDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("right_suppressor_extension");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorRightExtensionDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("left_suppressor_casing");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      //G4cout << "left_suppressor_casing_log Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinSuppressorLeftSideDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("right_suppressor_casing");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorRightSideDet(edep,stepl,det-1,cry-1);
  }

  // Dead layer specific code
  found = volname.find("germanium_dls_block1");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(volname);
      //G4cout << "germanium_dls_block1 Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
  }

  // LaBr
  found = volname.find("brillance380_crystal_block");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddLaBrCrystDet(edep,stepl,det-1);
  }

  // Sodium Iodide
  found = volname.find("sodium_iodide_crystal_block");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSodiumIodideCrystDet(edep,stepl,det-1);
  }

  // Sceptar
  found = volname.find("sceptar_square_scintillator_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSceptarSquareCrystDet(edep,stepl,det-1);
  }

  found = volname.find("sceptar_angled_scintillator_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSceptarAngledCrystDet(edep,stepl,det-1);
  }



}

void SteppingAction::SetDetAndCryNumberForGriffinComponent(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if( avOver9 == 47 ) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if( avOver99 == 47 ) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')+10+(cstr[5]-'0');
        impr = cstr[12]-'0';
    }

    det = (G4int)(ceil(av/numberOfAssemblyVols));
    cry = impr;

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if(avOver99 == 47) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')+10+(cstr[5]-'0');
        impr = cstr[12]-'0';
    }

    det = (G4int)(ceil(av/numberOfAssemblyVols))+1;
    cry = av - numberOfAssemblyVols*(det-1);

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetNumberForGenericDetector(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int volNameOver9 = cstr[11]-'0';
    if(volNameOver9 == 47) {
        // if volNameOver9 = 47 then this character is "_" and the imprint is between 0 - 9, else between 10 - 99
        det = cstr[10]-'0';
    }
    else {
        det = ((cstr[10]-'0')*10)+volNameOver9 ;
    }
    //G4cout << "Found electron ekin in " << volname << " det = " << det << G4endl;
}

G4int SteppingAction::FindTrueGriffinDetector(G4int det)
{
    G4int trueDet;
    trueDet = detector->griffinDetectorsMap[det-1];

    return trueDet;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
