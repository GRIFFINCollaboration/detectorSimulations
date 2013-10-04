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
{ }

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
  found = volname.find("germanium_block1_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("back_quarter_suppressor_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinSuppressorBackDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("left_suppressor_extension_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinSuppressorLeftExtensionDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("right_suppressor_extension_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinSuppressorRightExtensionDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("left_suppressor_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinSuppressorLeftSideDet(edep,stepl,det-1,cry-1);
  }

  found = volname.find("right_suppressor_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinDetector(volname);
      eventaction->AddGriffinSuppressorRightSideDet(edep,stepl,det-1,cry-1);
  }

  // LaBr
  found = volname.find("brillance380_crystal_block_log");
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

}

void SteppingAction::SetDetAndCryNumberForGriffinDetector(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int volNameOver9 = cstr[11]-'0';
    if(volNameOver9 == 47) {
        // if volNameOver9 = 47 then this character is "_" and the imprint is between 0 - 9, else between 10 - 99
        cry = cstr[10]-'0';
        det = (G4int)(ceil(cry/4.0));
    }
    else {
        cry = ((cstr[10]-'0')*10)+volNameOver9 ;
        det = (G4int)(ceil(cry/4.0));
    }
    cry = cry - 4*(det-1);
    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << G4endl;
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
