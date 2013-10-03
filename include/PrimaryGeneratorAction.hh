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
// $Id: PrimaryGeneratorAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class EventAction;
class G4ParticleTable;
class G4ParticleDefinition;
class BeamRequestBetaParticle;
class BeamRequestGammaAndIC;
class BeamRequestXRay;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction( DetectorConstruction* );    
    virtual ~PrimaryGeneratorAction();

    void GeneratePrimaries( G4Event* ) ;

    inline G4double GetParticleEnergy() {return energy;};
  
    void SetEnergy( G4double );
    void SetParticleType( G4String );
    void SetIonType( G4int Z, G4int A, G4double E );
    void SetDirection( G4ThreeVector );
    void SetPosition( G4ThreeVector );
    void DefineIsotropicRadOnBox( G4ThreeVector );
    
    void SetBetaPlusEmission( G4String );
    void SetBetaMinusEmission( G4String );
    void SetRadioactiveBetaDecay( G4String );
    void SetPolarization( G4double );
    void SetEmitBetaParticle( G4bool );
    void IncludeXRayInputFileKShell( G4bool );
    void IncludeXRayInputFileLShell( G4bool );
    void IncludeXRayInputFileMShell( G4bool );
    void SetRadioactiveDecayHalflife( G4double );
    void SetNumberOfRadioactiveNuclei( G4int );

    void ReadEnergyDistribution( G4String );

     G4String PrepareLine(); 
  
  private:
    G4double UniformRand48();
  
    G4ParticleDefinition* 		particle;
    G4ParticleTable* 					particleTable;
    DetectorConstruction* 		Detector ;

    BeamRequestBetaParticle* 	myBetaParticle;
    BeamRequestGammaAndIC* 		myGammaAndICParticle;
    BeamRequestXRay* 					myXRay;

    G4double      						energy;
    G4double      						halflife;
    G4double      						polarization;
    G4int         						numberOfNuclei;
    G4bool        						emit_beta_flag;
    G4bool        						emit_gamma_ic_flag;
    G4bool        						emit_xray_flag;
    G4bool        						xray_input_kShell;
    G4bool        						xray_input_lShell;
    G4bool        						xray_input_mShell;
    G4int         						Z;
    
    // isotropicRadOnBox
    G4bool 										isoRadOnBox;
    G4double 									totalLengthOfBox_x;
    G4double 									totalLengthOfBox_y;
    G4double 									totalLengthOfBox_z;
    G4double 									randBox_x;
    G4double 									randBox_y;
    G4double 									randBox_z;    
    G4double 									sourceShellRadius;

    G4double 									outputTimeInSeconds;
    G4double 									previousTimeInSeconds;

    G4String      						particleType;
    G4String      						simulationDir;

    G4ThreeVector 						position;
    G4ThreeVector 						direction;
    G4bool        						directionSpecified;
    G4bool        						newParticleType;
    G4bool        						emissionSimulation;
    G4bool        						radioactiveDecaySimulation;
    G4int         						eventID;  
    G4int         						numberOfEvents;  
    G4int         						eventSum;

    std::vector <double> 			energyDist;
    std::vector <double> 			weightDist;
    std::vector <double> 			monteCarlo;

    G4ParticleGun*           particleGun;	 //pointer a to G4  class
    
    PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
