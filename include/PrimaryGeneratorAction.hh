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

//c++  // MHD 1 May 2013 
#include <map>
#include <utility>  // pair 
#include <string>
#include <vector>
#include <iostream>

using namespace std;
using namespace CLHEP;

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
    void GetRandomDirection();

    inline G4double GetParticleEnergy() {return energy;};
    
    G4int GetPrimaryParticleType(); // MHD : 12 April 2013
  
  	
    void SetEnergy( G4double );
    void SetParticleType( G4String );
    void SetIonType( G4int Z, G4int A, G4double E );
    void SetDirection( G4ThreeVector );
    void SetPosition( G4ThreeVector );
    void SetSourceRadius( G4double );
    void SetEnergyRange( G4double minimum, G4double maximum, G4double step);
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
    
    void SetSourceRecord( G4double/*energy*/,string/*process name*/);
    void WriteSourceRecord(void);
    
 		// methods used for radioactive beta decay as done for spice
  	void SetRadioactiveSourceDecay( G4String myInputFile );
  	void EmissionForRadioactiveSourceDecay( G4Event* myEvent );
  	void LevelSchemeReader( const char* filename );
  	void EmitBetaForSourceDecay(G4double myEndPointEnergy, G4Event* myEvent);
  	void EmitParticleForSourceDecay(G4double energy, G4Event* myEvent);
  	void EmissionForVacantShell(int shell, G4Event* myEvent);
  	
  	void SetKinematicsActive( G4bool tf );
  	void SetKinematicsBetaValue( G4double beta );
  	void SetKinematicsIonEnergy( G4double value );
  	G4double KinematicEnergyBroadening( G4double energy, G4Event* anEvent );
  	void EmitIon(G4int ionZ, G4int ionA, G4double ionE, G4Event* anEvent);
  	
  	G4ThreeVector GetEmissionPositionAtSource( G4ThreeVector value );

    void ReadEnergyDistribution( G4String );

    G4String PrepareLine(); 
  
  private:
    G4double UniformRand48();
  
    G4ParticleDefinition* 		particle;
    G4ParticleTable* 					particleTable;
    DetectorConstruction* 		Detector ;
    
    G4bool radioactiveSourceDecaySimulation; 

  	// initialise variables used in the source decay generator
  	//***
  	// in method EmissionForRadioactiveSourceDecay
  	//***
  	G4double betaDecayRandomiser;
  	G4double betaEmissionRandomiser;
  	G4double levelRandomiser;
  	G4double levelProbSum;
  	G4double transitionConversionK;
  	G4double transitionConversionL1;
  	G4double transitionConversionL2;
  	G4double transitionConversionL3;
  	G4double particleEnergy;
  	G4double conversionTest;
  	G4double xrayTest;
  	//***
  	// in method levelSchemeReader
  	//***
  	G4bool levelSchemeRead;
  	vector<vector<double> > levelScheme;
  	// header information of level scheme table
  	// (Number of levels, number of parameters, binding energies)
  	G4int nLevels, nTransPerLevel, nParam; 
  	G4double bindingEnergyK, bindingEnergyL1, bindingEnergyL2, bindingEnergyL3;

  	// arrays to hold x-ray and auger information
		vector<double> fKXRayEnergy;
		vector<double> fKXRayIntensity;
		vector<double> fKXRayOrigin;

		vector<double> fL1XRayEnergy;
		vector<double> fL1XRayIntensity;
		vector<double> fL1XRayOrigin;

		vector<double> fL2XRayEnergy;
		vector<double> fL2XRayIntensity;
		vector<double> fL2XRayOrigin;

		vector<double> fL3XRayEnergy;
		vector<double> fL3XRayIntensity;
		vector<double> fL3XRayOrigin;

		vector<double> fAugerIntensity;
		vector<double> fAugerEnergy;
		vector<double> fAugerRecFrom;
		vector<double> fAugerEjecFrom;
                         
  	G4double KXRayEnergy[10],L1XRayEnergy[4],L2XRayEnergy[4],L3XRayEnergy[6];
  	G4int KXRayOrigin[10],L1XRayOrigin[4],L2XRayOrigin[4],L3XRayOrigin[6];
  	G4double KXRayIntensity[10],L1XRayIntensity[4],L2XRayIntensity[4],L3XRayIntensity[6];
  	G4double augerIntensity[48], augerEnergy[48];
  
  	// level information from level scheme table
  	G4double levelID, levelEnergy, levelSpin, levelParity, betaDecayProb;
    G4double betaParticleType, betaEmissionProb,endPointEnergy;
  	vector<double> daughterID;
  	vector<double> daughterBranching;
  	vector<double> ICProbabilityK;
  	vector<double> ICProbabilityL1;
  	vector<double> ICProbabilityL2;
  	vector<double> ICProbabilityL3;
  	//***
  	// in method EmitBetaForSourceDecay and EmitParticleForSourceDecay
  	//***
  	G4double direction_x, direction_y, direction_z;
  	G4double direction_norm, direction_norm2;
  
  	// source geometry
  	G4double position_x, position_y, position_z, radius;
  
  	//***
  	// in method EmitBetaForSourceDecay
  	//***
  	G4double betaRandomiser;
  	G4int betaEnergy;
  	G4double betaEnergyDouble;
  	vector<double> betaEnergyProbabilitySum;
  	//***
  	// in method EmissionForVacantShell
  	//***
  	G4double totalXRayIntensity, totalAugerIntensity, xOrARandom, addedIntensity, xRayRandom;

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
    
    G4double			minimumEnergy;
    G4double			maximumEnergy;
    G4double			stepSizeEnergy;
    G4double			previousEnergy;

    G4ThreeVector 						position;
    G4ThreeVector 						direction;
    G4bool        						directionSpecified;
    G4bool										energyRange;
    G4bool        						newParticleType;
    G4bool        						emissionSimulation;
    G4bool        						radioactiveDecaySimulation;
    G4int         						eventID;  
    G4int         						numberOfEvents;  
    G4int         						eventSum;
    
    G4double                  fSpeedOfLight;
	  G4double                  fRestMassOfElectron;
	  G4double                  fAtomicMassUnit;
	  
	  G4double                  fBetaHeavyIon;
    G4bool                    fBetaDefinedByUser;
    G4bool										fSimulateKinematics;
    G4bool										ionDefined;
    G4bool                    ionEmittedThisEvent;
    G4int											ionDefinitionZ;
    G4int											ionDefinitionA;
    G4double									ionDefinitionE;
    G4ThreeVector							ionDirection;
    G4double									fIonKineticEnergy;
    G4bool                    fIonEnergyDefinedByUser;
    
    
    
    std::vector <double> 			energyDist;
    std::vector <double> 			weightDist;
    std::vector <double> 			monteCarlo;
    
    // map to keep a record on the energies the primary 
		map<double, pair < string,int > > fMapSourceRecord ;

    G4ParticleGun*           particleGun;	 //pointer a to G4  class
    
    PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
