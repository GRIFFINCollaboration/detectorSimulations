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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define DEBUG  0

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "BeamRequestBetaParticle.hh"
#include "BeamRequestGammaAndIC.hh"
#include "BeamRequestXRay.hh"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:Detector(DC)
{
  srand48(time(NULL));
  
  

  directionSpecified         = false;
  newParticleType            = false;
  emissionSimulation         = false;
  radioactiveDecaySimulation = false;
  radioactiveDecaySimulation = false;
  energyRange                = false;
  previousEnergy             = 0*keV;
  ionDefined                 = false;
  
  fSpeedOfLight              = 299.8*km/s;
	fRestMassOfElectron        = 511*keV;
	fAtomicMassUnit            = 1.660539E-27*kg;

  eventSum = 0;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // defaults
  energy = 1000.*keV;
  
  fBetaHeavyIon               = 0.0;
  fBetaDefinedByUser          = false;
  fSimulateKinematics         = false;
  fIonEnergyDefinedByUser     = false;
  

  this->halflife              = 0.0; //seconds
  this->polarization          = 0.0;
  this->numberOfNuclei        = 0;
  this->outputTimeInSeconds   = 0.0;
  this->previousTimeInSeconds = 0.0;

  this->emit_beta_flag     = true;
  this->emit_gamma_ic_flag = true;
  this->emit_xray_flag     = true;

  this->xray_input_kShell  = false;
  this->xray_input_lShell  = false;
  this->xray_input_mShell  = false;
  
  this->radioactiveSourceDecaySimulation = false;

  // Messenger options
  isoRadOnBox               = false;

  // default particle position
  position = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);

  G4ParticleTable* myParticleTable = G4ParticleTable::GetParticleTable();
  this->particleTable = myParticleTable;
  G4String particleName = "gamma";
  G4ParticleDefinition* myParticle = this->particleTable->FindParticle(particleName);
  this->particle = myParticle;
  particleGun->SetParticleDefinition(this->particle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	WriteSourceRecord();
	delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double randomClock = 0.0, lamda = 0.0;
  G4double deltaTimeInSeconds = 0.0;
  G4int theBetaBranch, numberOfParticles, numberOfXRays;  
  ionEmittedThisEvent = false;  

	// if /DetSys/gun/radioactiveBetaDecay is called
  if(emissionSimulation) {
    G4double selectedEnergy = 0.00*keV;
    G4double monteCarloRand = UniformRand48();

		// energy distribution
    for(G4int i = 0 ; i < energyDist.size() ; i++ )
    {
      if(monteCarloRand <= monteCarlo[i])
      {
        selectedEnergy = energyDist[i]*1.0*keV;
        break;
      }
    } // end for(G4int i
    
    energy = selectedEnergy;
    
    if( isoRadOnBox || !directionSpecified) GetRandomDirection();
    
    // If the particle is an electron/positron and the user has decided to simulate kinematic
    // broadening, do so.
    if (this->particle->GetParticleName() == "e-" || this->particle->GetParticleName() == "e-") 
      if (fSimulateKinematics) energy = KinematicEnergyBroadening(energy, anEvent);
    
    particleGun->SetParticleDefinition(this->particle);
    particleGun->SetParticlePosition(position);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleTime(randomClock);

    eventID = anEvent->GetEventID();
    eventSum++;
    particleGun->GeneratePrimaryVertex(anEvent);
  } // end if /DetSys/gun/radioactiveBetaDecay
  else if(this->radioactiveDecaySimulation) {
    if(this->halflife > 0 && this->numberOfNuclei > 0)
    {
      randomClock = UniformRand48();
      lamda = log(2)/( this->halflife );    // lamda in seconds
      deltaTimeInSeconds = (-1/(( (double)(this->numberOfNuclei - eventID) )*lamda))*(log(randomClock)); // calculate delta time
      this->outputTimeInSeconds = deltaTimeInSeconds + this->previousTimeInSeconds;
      this->numberOfNuclei = this->numberOfNuclei - 1;
    }
    else
    {
      this->outputTimeInSeconds = 0.0;
    }
    this->previousTimeInSeconds = this->outputTimeInSeconds;

    //==================== Simulate Decay ============================//
    theBetaBranch = myGammaAndICParticle->GetBetaBranch();  // get beta branch from gamma data
    if(emit_beta_flag)
    { 
        if(!(myGammaAndICParticle->IsECDecay())) {
            myBetaParticle->AddBeamParticle(particleGun, theBetaBranch, this->outputTimeInSeconds);  // load particleGun with beta
            particleGun->SetParticleTime( this->outputTimeInSeconds*second );  // tag beta-particle with decay time
            particleGun->GeneratePrimaryVertex(anEvent);  // fire beta event
        }
    }
    myGammaAndICParticle->GenerateParticles(this->outputTimeInSeconds);  // generate gamma and IC partile decay
    //myGammaAndICParticle->PrintAllGeneratedParticles();
    numberOfParticles = myGammaAndICParticle->GetNumberOfGeneratedParticles();
    for( G4int i = 0 ; i < numberOfParticles ; i++ )  // get the number of particles in decay and loop over them
    {
      if(emit_gamma_ic_flag)
      {
        myGammaAndICParticle->AddBeamParticle(particleGun, i, this->outputTimeInSeconds);  // load particleGun with first particle in decay
        particleGun->SetParticleTime( this->outputTimeInSeconds*second );  // tag particle with decay time
        particleGun->GeneratePrimaryVertex(anEvent);  // fire particle event
      }
      if( emit_xray_flag && myGammaAndICParticle->IsParticleAnElectron(i) )  // if particle was an electron we need to emit an x-ray
      {
        myXRay->GenerateParticles(myGammaAndICParticle->GetParticleShell(i));
        numberOfXRays = myXRay->GetNumberOfGeneratedParticles();
        for( G4int j = 0 ; j < numberOfXRays ; j++ )  // get the number of x-rays and loop over them
        {
          myXRay->AddBeamParticle(particleGun, j);  // get x-ray from empty shell
          particleGun->SetParticleTime( this->outputTimeInSeconds*second ); // tag x-ray with decay time
          particleGun->GeneratePrimaryVertex(anEvent);  // fire x-ray event
        }
        //myXRay->PrintAllGeneratedParticles();
      }
    }
    //================================================================//

  }
  else if(radioactiveSourceDecaySimulation)
  {
      EmissionForRadioactiveSourceDecay(anEvent);
  } //end if(radioactiveSourceDecaySimulation)
  else if(energyRange)
  {
		if( isoRadOnBox || !directionSpecified) GetRandomDirection();
	
		// Stepping
		if (previousEnergy/keV == 0 || (abs(previousEnergy-maximumEnergy+stepSizeEnergy))/keV < 0.001)
		{
			energy = minimumEnergy;
			previousEnergy = minimumEnergy;
		}
		else
		{
			previousEnergy = energy;
			energy = energy + stepSizeEnergy;
		}
		
		// In this instance a new variable has to be created with the energy, as the above if statement
		// relies on the energy prior to it being corrected by kinematics.
		G4double correctedEnergy = energy;
    if (this->particle->GetParticleName() == "e-" || this->particle->GetParticleName() == "e-") 
      if (fSimulateKinematics) correctedEnergy = KinematicEnergyBroadening(correctedEnergy, anEvent);
    
		particleGun->SetParticleDefinition(this->particle);
    particleGun->SetParticlePosition(position);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->SetParticleEnergy(correctedEnergy);
    particleGun->SetParticleTime(randomClock);
    
  	eventID = anEvent->GetEventID();
  	eventSum++;
  	particleGun->GeneratePrimaryVertex(anEvent);	
  	  
  } // end if(energyRange)
  else 
  {
    
    if( isoRadOnBox || !directionSpecified) GetRandomDirection();
    
    // If the particle is an electron/positron and the user has decided to simulate kinematic
    // broadening, do so.
    if (this->particle->GetParticleName() == "e-" || this->particle->GetParticleName() == "e-") 
      if (fSimulateKinematics) energy = KinematicEnergyBroadening(energy, anEvent);

    particleGun->SetParticleDefinition(this->particle);
    particleGun->SetParticlePosition(position);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleTime(randomClock);

    eventID = anEvent->GetEventID();
    eventSum++;
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GetRandomDirection()
{
	G4double posx, posy, posz ;
	
	// if DefineIsotropicRadOnBox() has been called
	if( isoRadOnBox )
	{
     posx = position.x();
     posy = position.y();
     posz = position.z();

     randBox_x = ( UniformRand48() * totalLengthOfBox_x ) - ( totalLengthOfBox_x / 2.0 ) ;
     randBox_y = ( UniformRand48() * totalLengthOfBox_y ) - ( totalLengthOfBox_y / 2.0 ) ;
     randBox_z = ( UniformRand48() * totalLengthOfBox_z ) - ( totalLengthOfBox_z / 2.0 ) ;
      
     direction = G4ThreeVector( randBox_x-posx, randBox_y-posy, randBox_z-posz) ;
	 }	// end if( isoRadOnBox )

	// if direction has not been specified, emit in random direction
   else if(!directionSpecified)
   {
     G4double costheta = 2.*UniformRand48()-1.0;
     G4double sintheta = sqrt( 1. - costheta*costheta );
     G4double phi      = (360.*deg)*UniformRand48();
     direction = G4ThreeVector(sintheta*cos(phi), sintheta*sin(phi), costheta);
     
   } // end if(!directionSpecified)
   
   	
}

void PrimaryGeneratorAction::SetEnergy( G4double value )
{
  if( value > 0. ) {
    energy = value;
    G4cout << " --> Particle energy has been set to " << energy/keV << " keV" << G4endl;
  }
  else {
    G4cout << " --> Invalid value, keeping previous energy (" << energy/keV << " keV)" << G4endl;
  }
}

void PrimaryGeneratorAction::SetEnergyRange( G4double minimum, G4double maximum, G4double step)
{
	if (minimum + maximum + step > 0.)
	{
		minimumEnergy = minimum*keV;
    maximumEnergy = maximum*keV;
    stepSizeEnergy = step*keV;
    G4cout << "Minimum is " << minimum/keV << G4endl;
    G4cout << "Maximum is " << maximum/keV << G4endl;
    G4cout << "Step size is " << stepSizeEnergy/keV << G4endl;
    energyRange = true;
	}
	else
	{
		minimumEnergy = 100*keV;
    maximumEnergy = 1000*keV;
    stepSizeEnergy = 100*keV;
		G4cout << " --> Invalid value, set to defaults" << G4endl;
	}
}

void PrimaryGeneratorAction::SetSourceRadius( G4double value )
{
  if( value > 0. )
  		{
		radius = value;
		G4cout << " --> Source radius has been set to " << radius/mm << " mm" << G4endl;
	 	}
  else
  	{
    G4cout << " --> Invalid value, keeping default radius (" << radius/mm << " mm)" << G4endl;
  	}
}

void PrimaryGeneratorAction::SetParticleType( G4String name )
{
  if( name == "gamma" || name == "Gamma" ) 
  {
    particleType = "gamma";
  }
  else if( name == "electron" || name == "Electron" || name == "e-" ) 
  {
    particleType = "e-";
  }
  else if( name == "positron" || name == "Positron" || name == "e+" ) 
  {
    particleType = "e+";
  }
  else if( name == "proton" || name == "Proton") 
  {
    particleType = "proton";
  }
  else if( name == "alpha" || name == "Alpha") 
  {
    particleType = "alpha";
  }
  else
  {
    G4cout << " --> Can't find particle type! Defaulting to gamma!" << G4endl;
    particleType = "gamma";

  }

  this->particle = this->particleTable->FindParticle(particleType);
  particleGun->SetParticleDefinition(this->particle);
  //G4cout << " --> Particle type has been set to " << particleType << G4endl;
}

void PrimaryGeneratorAction::SetIonType( G4int Z, G4int A, G4double E )
{
	// If ion hasn't been defined previously, assign Z, A & E to variables outside 
	// the scope of this function... this allows the ion to be emitted repeatedly 
	// at each new event along with gammas, electrons, etc 
	// (e.g. in the KinematicEnergyBroadening function)
	if (!ionDefined) {
		ionDefinitionZ = Z;
		ionDefinitionA = A;
		ionDefinitionE = E;
		ionDefined = true;
	}
  this->particle = this->particleTable->GetIon(Z,A,E*keV);
  particleGun->SetParticleDefinition(this->particle);
  G4cout << " --> Particle type has been set to an ion with Z=" << Z << " and A=" << A 
  			 << " and an excitation energy of " << E << " keV." << G4endl;
}

void PrimaryGeneratorAction::DefineIsotropicRadOnBox( G4ThreeVector boxSize )
{
    G4double maxLength;

    isoRadOnBox = true;

    totalLengthOfBox_x = boxSize.x();
    totalLengthOfBox_y = boxSize.y();
    totalLengthOfBox_z = boxSize.z();

    if(totalLengthOfBox_x > totalLengthOfBox_y) {
        maxLength = totalLengthOfBox_x;
    }
    else {
        maxLength = totalLengthOfBox_y;
    }
    if(maxLength < totalLengthOfBox_z) {
        maxLength = totalLengthOfBox_z;
    }

    sourceShellRadius =  maxLength;
}

void PrimaryGeneratorAction::SetDirection( G4ThreeVector value )
{
  G4ThreeVector null = G4ThreeVector(0.0*mm,0.0*mm,0.0*mm);
  if( value != null) {
    direction = value;
    directionSpecified = true;
    G4cout << " --> Particle direction has been set to (" << value.x() << ", " << value.y() << ", " << value.z() << ")" << G4endl;
  }
  else {
    directionSpecified = false;
    G4cout << " --> Invalid value, using isotropic distribution" << G4endl;
  }
}

void PrimaryGeneratorAction::SetPosition( G4ThreeVector value )
{

  if(isoRadOnBox) 
		{
	    // random position on "shell" source
	    G4double costheta = 2. * UniformRand48() - 1.0;
	    G4double sintheta = sqrt( 1. - costheta * costheta );
	    G4double phi      = (360.*deg) * UniformRand48();
	    position = G4ThreeVector( sourceShellRadius*sintheta*cos(phi), sourceShellRadius*sintheta*sin(phi), sourceShellRadius*costheta);
		}
	else
  	position = value;
  
  G4cout << " --> Particle position has been set to (" << value.x()/mm << ", " << value.y()/mm << ", " << value.z()/mm << ") in mm" << G4endl;

}

G4ThreeVector PrimaryGeneratorAction::GetEmissionPositionAtSource( G4ThreeVector value )
{
    //initiate the new position
    G4ThreeVector new_position = G4ThreeVector(0,0,0);
    G4double x_offset=0;
    G4double y_offset=0;
   
    //randomize in the source radius
    do
    {
     x_offset = RandFlat::shoot(-radius,+radius) ; //cout << x_offset <<endl;
     y_offset = RandFlat::shoot(-radius,+radius) ; //cout << y_offset <<endl;
    } while ( (x_offset*x_offset + y_offset*y_offset) > radius*radius) ;
    
    //set the new position
    new_position.set(value.x() + x_offset, value.y() + y_offset, value.z());

   // G4cout << " --> Particle position Emitted at source (" << new_position.x()/mm << ", " << new_position.y()/mm << ", " << new_position.z()/mm << ") in mm" << G4endl;
    
    
return new_position;

}

G4String PrimaryGeneratorAction::PrepareLine()
{
  char line[256];
  
  sprintf(line, " -1 %7.3f %7.3f %7.3f %7.3f %10d\n", energy/keV, direction.x(), direction.y(), direction.z(), eventID );
  
  return G4String(line);  
}

void PrimaryGeneratorAction::SetBetaPlusEmission( G4String file )
{
  emissionSimulation = true;
  G4double sumWeight = 0;
  G4double percentage = 0;

  particleType = "e+";
  this->particle = this->particleTable->FindParticle(particleType);
  particleGun->SetParticleDefinition(this->particle);
  G4cout << " --> Particle type has been set to " << particleType << G4endl;

  ReadEnergyDistribution(file);
  G4cout << " --> Number of events listed in file is " << numberOfEvents << G4endl;

  if(energyDist.size() != weightDist.size())
  {
    G4cout << "Error: 9829384" << G4endl;
    exit (1);
  } 

  for(G4int i = 0 ; i < weightDist.size() ; i++ )
  {
    sumWeight = sumWeight + weightDist[i];
  }

  for(G4int i = 0 ; i < weightDist.size() ; i++ )
  {
    percentage = weightDist[i]/sumWeight;
    if(i == 0)
    {
      monteCarlo.push_back(percentage);
    }
    else
    {
      monteCarlo.push_back(percentage + monteCarlo[i-1]);
    }
  }
}

void PrimaryGeneratorAction::SetBetaMinusEmission( G4String file )
{
  emissionSimulation = true;
  G4double sumWeight = 0;
  G4double percentage = 0;

  particleType = "e-";
  this->particle = this->particleTable->FindParticle(particleType);
  particleGun->SetParticleDefinition(this->particle);
  G4cout << " --> Particle type has been set to " << particleType << G4endl;

  ReadEnergyDistribution(file);
  G4cout << " --> Number of events listed in file is " << numberOfEvents << G4endl;

  if(energyDist.size() != weightDist.size())
  {
    G4cout << "Error: 9829384" << G4endl;
    exit (1);
  } 

  for(G4int i = 0 ; i < weightDist.size() ; i++ )
  {
    sumWeight = sumWeight + weightDist[i];
  }

  for(G4int i = 0 ; i < weightDist.size() ; i++ )
  {
    percentage = weightDist[i]/sumWeight;
    if(i == 0)
    {
      monteCarlo.push_back(percentage);
    }
    else
    {
      monteCarlo.push_back(percentage + monteCarlo[i-1]);
    }
  }
}

// -----------------------------------------------------
// The next several functions are used in simulating a decay scheme (text file),
// primarily used for SPICE
// -----------------------------------------------------

// Messenger set-up
void PrimaryGeneratorAction::SetRadioactiveSourceDecay( G4String myInputFile )
{
  radioactiveSourceDecaySimulation = true;
  PrimaryGeneratorAction::LevelSchemeReader(myInputFile);
}

// Superordinate function for decay
void PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay( G4Event* myEvent)
{
  // Variables used in control-loops
  G4int i = 0;
  G4int j = 0;
  G4int k = 0;
  
  // Find starting point of the decay by comparing beta decay probability to random number
    	
  // RandFlat is a random-number generating class from CLHEP
  // RandFlat::shoot(m,n) excludes m and n itself!
  betaDecayRandomiser = RandFlat::shoot(0.,100.);
  levelProbSum=0;
  while( (levelProbSum += levelScheme[i][4]) < betaDecayRandomiser )
    {
      i++;
      // make sure that i is running only through table elements
      if( i==nLevels )
		{
		  G4cout << "Total beta decay probability found = " <<levelProbSum<< endl;
		  i=0;
		  levelProbSum=0;
		  betaDecayRandomiser = RandFlat::shoot(0.,100.);
		  G4cout << "Total beta decay probability below 100%." << endl;
		  G4cout << "Repeat finding beta." << endl;
		}
    }
  
  // check for the emission probability of the generated beta particle
  // always 100 for e-
  // <100 for beta+/EC competition
  // there is no difference between e- and e+
  // is the probability normalised to 100?
  // !!CHECK!!
  // betaEmissionRandomiser = RandFlat::shoot(0.,100.);
  // if(betaEmissionRandomiser <= levelScheme[i][6])
  // EmitBetaForSourceDecay(levelScheme[i][7],myEvent);
  
  while(i!=0)
    {     
     /*
	Find the right decay branch
	rename initial level to j
	compare random number to branching ratio
	to get final ID i, and conversion coefficients
      */
      j = i;
      k = 0;
      levelProbSum = 0.;
      levelRandomiser = RandFlat::shoot(0.,100.);
      while( (levelProbSum += levelScheme[j][9+k*6]) < levelRandomiser )
	{
	  k++;
	  if( k==nTransPerLevel )
	    {
	      k=0;
	      levelProbSum=0;
	      levelRandomiser = RandFlat::shoot(0.,100.);
	      G4cout << "Total branching ratio of level " << j << " below 100%." << endl;
	      G4cout << "Repeat finding decay branch." << endl;
	    }
	}
      i=(int)levelScheme[j][8+k*6];
      transitionConversionK=levelScheme[j][10+k*6];
      transitionConversionL1=levelScheme[j][11+k*6];
      transitionConversionL2=levelScheme[j][12+k*6];
      transitionConversionL3=levelScheme[j][13+k*6];
      
      /*
	compute transition energy
	and emit gamma or electron after conversion test
      */      
      particleEnergy = levelScheme[j][1]-levelScheme[i][1];
      conversionTest = RandFlat::shoot(0.,1.);
      if(transitionConversionK > conversionTest)
		{
		  SetParticleType( "e-" ) ; // MHD 12 April 2013
		  if (DEBUG) G4cout << "PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay K-IC " << (particleEnergy-bindingEnergyK)*keV << endl;
		  EmitParticleForSourceDecay((particleEnergy-bindingEnergyK)*keV,myEvent); // MHD 12 April 2013
		  SetSourceRecord( (particleEnergy-bindingEnergyK)*keV , "IC-k"); // MHD 01 May 2013
		  EmissionForVacantShell(0, myEvent);
		}
      else if((transitionConversionK + transitionConversionL1) > conversionTest)
		{
		  SetParticleType( "e-" ) ; // MHD 12 April 2013
		  if (DEBUG) G4cout << "PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay L1-IC " <<  (particleEnergy-bindingEnergyL1)*keV << endl;
		  EmitParticleForSourceDecay((particleEnergy-bindingEnergyL1)*keV,myEvent); // MHD 12 April 2013
		  SetSourceRecord( (particleEnergy-bindingEnergyL1)*keV , "IC-L1"); // MHD 01 May 2013
		  EmissionForVacantShell(1, myEvent);
		}
      else if((transitionConversionK + transitionConversionL1 + transitionConversionL2) > conversionTest)
		{
		  SetParticleType( "e-" ) ; // MHD 12 April 2013
		  if (DEBUG) G4cout << "PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay L2-IC " <<(particleEnergy-bindingEnergyL2)*keV<< endl;
		  EmitParticleForSourceDecay((particleEnergy-bindingEnergyL2)*keV,myEvent); // MHD 12 April 2013
		  SetSourceRecord( (particleEnergy-bindingEnergyL2)*keV , "IC-L2"); // MHD 01 May 2013
		  EmissionForVacantShell(2, myEvent);
		}
      else if((transitionConversionK + transitionConversionL1 + transitionConversionL2 + transitionConversionL3) > conversionTest)
		{
		  SetParticleType( "e-" ) ; // MHD 12 April 2013
		  if (DEBUG) G4cout << "PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay L3-IC " <<(particleEnergy-bindingEnergyL3)*keV<< endl;
		  EmitParticleForSourceDecay((particleEnergy-bindingEnergyL3)*keV,myEvent); // MHD 12 April 2013
		  SetSourceRecord( (particleEnergy-bindingEnergyL3)*keV , "IC-L3"); // MHD 01 May 2013
		  EmissionForVacantShell(3, myEvent);
		}
      else
		{
		 SetParticleType( "gamma" ) ; // MHD 12 April 2013
		 if (DEBUG) G4cout << "PrimaryGeneratorAction::EmissionForRadioactiveSourceDecay gamma"<< (particleEnergy)*keV << endl;
		 EmitParticleForSourceDecay((particleEnergy)*keV,myEvent); // MHD 12 April 2013
		 SetSourceRecord( (particleEnergy)*keV , "gamma"); // MHD 01 May 2013
		} 
    } // end while(i!=0) 
  
}

// -------------------------------------------------------
// leve scheme reader for radioactive source decay
// -------------------------------------------------------

void PrimaryGeneratorAction::LevelSchemeReader( const char* filename )
{
  G4cout << "\n-----------------------------------------------------"
         << "\n    Level Scheme"
         << "\n-----------------------------------------------------";
  
  G4cout << "\n ----> Reading the level scheme from "
	 << filename << "..." << endl;
	 
	//G4cout << "\n Enter to continue! "<< endl;	G4cin.get();
	 	   
  ifstream file( filename ); // open file for reading
  
  // ignore first blank line
  string buffer;
  getline(file,buffer); G4cout <<"skipping this line :" <<buffer << endl; 
  
  // read in number of levels
  file >> nLevels >> nTransPerLevel >> nParam >> bindingEnergyK >> bindingEnergyL1 >> bindingEnergyL2 >> bindingEnergyL3;
  G4cout << " [Number of levels "
	 << nLevels << "]"
	 << endl;
  G4cout << " [Number of transition per level "
	 << nTransPerLevel << "]"
	 << endl;
  G4cout << " [binding energies "<< bindingEnergyK<< " " << bindingEnergyL1 << " " <<bindingEnergyL2<< " " << bindingEnergyL3 << "]" << endl;
  	
  // jump two lines to line 5
  getline(file,buffer);  	G4cout <<"skipping this line :" <<buffer << endl; 
  getline(file,buffer);     G4cout <<"skipping this line :" <<buffer << endl; 
  // Read in the energies & intensities of the K, L lines
  G4int shellFrom [24] = {3,2,1,6,10,5,12,7,17,26,8,7,7,12,6,5,19,9,12,10,11,19,4,4};
  G4int shellTo [24] = {0,0,0,0,0,0,0,0,0,0,3,3,2,3,1,1,3,3,2,1,1,2,2,3};
  G4int k = 0;
  G4int l1 = 0;
  G4int l2 = 0;
  G4int l3 = 0;
  
  G4double dumDouble = 0;

  for(G4int i=0;i<24;i++)
    {
      if(shellTo[i]==0)
	{
	  file >> KXRayEnergy[k];
	  KXRayOrigin[k]=shellFrom[i];
	  
	  fKXRayEnergy.push_back(KXRayEnergy[k]);
	  fKXRayOrigin.push_back(shellFrom[i]);
	  
	  k++;
	}
      else if(shellTo[i]==1)
	{
	  file >> L1XRayEnergy[l1];
	  L1XRayOrigin[l1]=shellFrom[i];
	  
	  fL1XRayEnergy.push_back(L1XRayEnergy[l1]);
	  fL1XRayOrigin.push_back(shellFrom[i]);
	  
	  l1++;
	}
      else if(shellTo[i]==2)
	{
	  file >> L2XRayEnergy[l2];
	  L2XRayOrigin[l2]=shellFrom[i];
	  
	  fL2XRayEnergy.push_back(L2XRayEnergy[l2]);
	  fL2XRayOrigin.push_back(shellFrom[i]);
	  
	  l2++;
	}
      else if(shellTo[i]==3)
	{
	  file >> L3XRayEnergy[l3];
	  L3XRayOrigin[l3]=shellFrom[i];
	  
	  fL3XRayEnergy.push_back(L3XRayEnergy[l3]);
	  fL3XRayOrigin.push_back(shellFrom[i]);
	  
	  l3++;
	}
	
    }
  getline(file,buffer); G4cout <<"skipping this line :" <<buffer << endl; 
  

  for(G4int i=0;i<24;i++)
    {
    	file >> dumDouble; 
      if(shellTo[i]==0)
		{ 
		fKXRayIntensity.push_back(dumDouble);
		}   
      else if(shellTo[i]==1)
		{ 
		fL1XRayIntensity.push_back(dumDouble);
		}
      else if(shellTo[i]==2)
		{ 
		fL2XRayIntensity.push_back(dumDouble);
		 }
      else if(shellTo[i]==3)
		{ 
		fL3XRayIntensity.push_back(dumDouble);
		}	
    } 
   
    
  // Skip two more lines and read in the Auger intensities	
  getline(file,buffer);  G4cout <<"skipping this line :" <<buffer << endl; 
  getline(file,buffer);  G4cout <<"skipping this line :" <<buffer << endl; 
  // table of Auger electrons resulting from recombination to K shell:
  // augerRecFrom gives the level, from which the recombination happens
  // augerEjecFrom gives the level, from which the Auger electron is ejected
  // The numbering goes L1 - 1, L2 - 3, L3 - 3, M1 - 4, etc.
  for(G4int i = 0; i<48; i++)
    	{ file >> dumDouble;  fAugerIntensity.push_back(dumDouble) ; }
 
  //skip line  
  getline(file,buffer); G4cout <<"skipping this line :" <<buffer << endl;
  for(G4int i = 0; i<48; i++)
   	    { file >> dumDouble;  fAugerEnergy.push_back(dumDouble) ; }
 
   	  
  // set up storage space for table for magnetic field
  levelScheme.resize( nLevels );
  
  G4int ix = 0;
  G4int iy = 0;
  for(ix=0; ix<nLevels; ix++)
    {
      levelScheme[ix].resize(nParam);
    }
  G4cout << "Level scheme storage set." << endl;

  // skip another 3 lines
  getline(file,buffer);  	G4cout <<"skipping this line :" <<buffer << endl;
  getline(file,buffer);     G4cout <<"skipping this line :" <<buffer << endl;
  getline(file,buffer);     G4cout <<"skipping this line :" <<buffer << endl;
  getline(file,buffer);     G4cout <<"skipping this line :" <<buffer << endl;
  G4cout << "reading in data" << endl;

  // resize vectors
  daughterID.resize(nTransPerLevel);
  daughterBranching.resize(nTransPerLevel);
  ICProbabilityK.resize(nTransPerLevel);
  ICProbabilityL1.resize(nTransPerLevel);
  ICProbabilityL2.resize(nTransPerLevel);
  ICProbabilityL3.resize(nTransPerLevel);

  /*
     for each level, insert the level information to the levelScheme array
     then loop through the daughter levels and insert the ID, branching and
     internal conversion probability for each level
  */
  for(ix=0;ix<nLevels;ix++)
    {
      file >> levelID >> levelEnergy >> levelSpin >> levelParity
	   >> betaDecayProb >> betaParticleType >> betaEmissionProb
	   >> endPointEnergy;

      levelScheme[ix][0]	= levelID;
      levelScheme[ix][1]	= levelEnergy;
      levelScheme[ix][2]	= levelSpin;
      levelScheme[ix][3]	= levelParity;
      levelScheme[ix][4]	= betaDecayProb;
      levelScheme[ix][5]	= betaParticleType;
      levelScheme[ix][6]	= betaEmissionProb;
      levelScheme[ix][7]	= endPointEnergy;
/*
	G4cout << "levelID            " << levelScheme[ix][0] << endl;
	G4cout << "levelEnergy        " << levelScheme[ix][1] << endl;
	G4cout << "levelSpin          " << levelScheme[ix][2] << endl;
	G4cout << "levelParity        " << levelScheme[ix][3] << endl;
	G4cout << "betaDecayProb      " << levelScheme[ix][4] << endl;
	G4cout << "betaParticleType   " << levelScheme[ix][5] << endl;
	G4cout << "betaEmissionProb   " << levelScheme[ix][6] << endl;
	G4cout << "endPointEnergy     " << levelScheme[ix][7] << endl;
*/
	
      for(iy=0;iy<nTransPerLevel;iy++)
	{
		file >> daughterID[iy]
			>> daughterBranching[iy]
			>> ICProbabilityK[iy]
			>> ICProbabilityL1[iy]
			>> ICProbabilityL2[iy]
			>> ICProbabilityL3[iy];

		levelScheme[ix][6*iy+8]   = daughterID[iy];
		levelScheme[ix][6*iy+8+1] = daughterBranching[iy];
		levelScheme[ix][6*iy+8+2] = ICProbabilityK[iy];
		levelScheme[ix][6*iy+8+3] = ICProbabilityL1[iy];
		levelScheme[ix][6*iy+8+4] = ICProbabilityL2[iy];
		levelScheme[ix][6*iy+8+5] = ICProbabilityL3[iy];
/*		
		G4cout << "daughterID[" << iy << "]       " << levelScheme[ix][6*iy+8+0] << endl;
		G4cout << "daughterBranching[" << iy << "]" << levelScheme[ix][6*iy+8+1] << endl;
		G4cout << "ICProbabilityK[" << iy << "]   " << levelScheme[ix][6*iy+8+2] << endl;
		G4cout << "ICProbabilityL1[" << iy << "]  " << levelScheme[ix][6*iy+8+3] << endl;
		G4cout << "ICProbabilityL2[" << iy << "]  " << levelScheme[ix][6*iy+8+4] << endl;
		G4cout << "ICProbabilityL3[" << iy << "]  " << levelScheme[ix][6*iy+8+5] << endl;
	*/
	}

      //G4cout << "Level ID: " <<  "\t" <<  levelScheme[ix][0] << G4endl;
      //G4cout << "First daughter ID: " <<  "\t" << levelScheme[ix][8] << G4endl;
    }
  file.close();
  
  G4cout << "---> ...finished reading " << endl;
} // end levelSchemeReader

// ------------------------------------------
// Emit Beta for Decay of Radioactive Source
// ------------------------------------------
void PrimaryGeneratorAction::EmitBetaForSourceDecay(G4double myEndPointEnergy, G4Event* myEvent)
{  
  // Get random direction for emission
  if( isoRadOnBox || !directionSpecified) GetRandomDirection();
  
  // compare a random number to the energy distribution of
  // the beta decay to get correct energy
  // maximum energy is endpoint energy given from table
  G4int numberOfChannels = (G4int)(myEndPointEnergy);
  betaEnergyProbabilitySum.resize(numberOfChannels);
  betaEnergyProbabilitySum[0] = 0.;
  for(int i=1;i<numberOfChannels;i++)
    {
      betaEnergyProbabilitySum[i] = 0;
      betaEnergyProbabilitySum[i] = sqrt((i+1022)*i)
	*(numberOfChannels-i)*(numberOfChannels-i)*(i+511)
	+ betaEnergyProbabilitySum[i-1];
    }
  
  betaRandomiser = 0.;
  betaRandomiser = RandFlat::shoot(0.,1.);
  betaEnergy = 0;
  while(betaEnergyProbabilitySum[betaEnergy]<betaRandomiser*betaEnergyProbabilitySum[numberOfChannels-1])
    {
      betaEnergy++;
    }
  betaEnergyDouble = (G4double)betaEnergy*keV;
  
  SetParticleType( "e-" ) ; // MHD 12 April 2013
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));  // MHD 12 April 2013
  particleGun->SetParticleMomentumDirection(direction);  // MHD 12 April 2013
  particleGun->SetParticleEnergy(betaEnergyDouble);  // MHD 12 April 2013
  particleGun->GeneratePrimaryVertex(myEvent);  // MHD 12 April 2013
  
} // end EmitBetaForSourceDecay(...)

// Function called each time a particle (e.g. K-electron, x-ray or Auger) needs to be fired
void PrimaryGeneratorAction::EmitParticleForSourceDecay(G4double energy, G4Event* myEvent)
{
  // Get random direction for emission
  if( isoRadOnBox || !directionSpecified) GetRandomDirection();
  
  // If the particle is an electron/positron and the user has decided to simulate kinematic
  // broadening, do so.
  if (this->particle->GetParticleName() == "e-" || this->particle->GetParticleName() == "e-") 
    if (fSimulateKinematics) energy = KinematicEnergyBroadening(energy, myEvent);

  particleGun->SetParticlePosition(GetEmissionPositionAtSource(position));
  particleGun->SetParticleMomentumDirection(direction);
  particleGun->SetParticleEnergy(energy);
  particleGun->GeneratePrimaryVertex(myEvent);
  
} // end EmitParticleForSourceDecay(...)

// --------------------------------------------------------
// emission of x rays and auger electrons for vacant shell
// --------------------------------------------------------
void PrimaryGeneratorAction::EmissionForVacantShell(int shell, G4Event* myEvent)
{
  /*
    for the shell vacated, get the X-ray and Auger arrays
    calculate the total X-ray and Auger intensity
    and find the process by comparing to a random number
    (no process is possible, as we do not have complete data)
    update the vacated shell and repeat the process if it is K, L1-3
  */
  // Initialise the auger arrays
  G4int augerRecFrom[48]  = {1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,6,6,6,6,6,6};
  G4int augerEjecFrom[48] = {1,2,3,4,5,6,7,8,9,10,11,2,3,4,5,6,7,8,9,10,11,13,3,4,5,6,7,8,9,10,11,12,13,4,5,6,9,10,11,6,9,11,6,7,8,9,10,11};

  // K-Shell vacancy
  if(shell == 0)
    {
      // sum up the X-ray intensities and Auger intensities
      totalXRayIntensity = 0.;     
      for(G4int i=0;i<fKXRayIntensity.size();i++)
	totalXRayIntensity += fKXRayIntensity.at(i);
	
      totalAugerIntensity = 0.;
      for(G4int i = 0; i<fAugerIntensity.size(); i++)
	totalAugerIntensity += fAugerIntensity.at(i);
      
      // use random number to determine process
      xOrARandom = RandFlat::shoot(0.,100);
      if(totalXRayIntensity >= xOrARandom)
	{ 	  	  	  
	  G4int i = 0;
	  addedIntensity = fKXRayIntensity.at(0);
	  while( addedIntensity < xOrARandom )
	    {
	      i++;
	      addedIntensity += fKXRayIntensity.at(i);
	      if(i>fKXRayIntensity.size())
			{
			  G4cout << "Reached end of array for K X ray intensities." << endl;
			  exit (1);
			}
	    }

	 SetParticleType( "gamma" ) ; // MHD 12 April 2013
	 //G4cout << "PrimaryGeneratorAction::EmissionForVacantShell, K-XRay" << endl;
	 EmitParticleForSourceDecay((fKXRayEnergy.at(i))*keV,myEvent); // MHD 12 April 2013
	 SetSourceRecord((fKXRayEnergy.at(i))*keV,"K-XRay");
	 EmissionForVacantShell(fKXRayOrigin.at(i), myEvent);	  
	}
      else if( (totalXRayIntensity+totalAugerIntensity) >= xOrARandom)
	{
	  G4int i = 0;
	  addedIntensity = totalXRayIntensity+fAugerIntensity.at(0);
	  while( addedIntensity < xOrARandom )
	    {
	      i++;
	      addedIntensity += fAugerIntensity.at(i);
	      if(i>fAugerEnergy.size())
			{
			  G4cout << "Reached end of array for auger intensities." << endl;
			  exit (1);
			}
	    }
	    
	  SetParticleType( "e-" ) ; // MHD 12 April 2013
	  //G4cout << "PrimaryGeneratorAction::EmissionForVacantShell Auger" << endl;
	  SetSourceRecord((fAugerEnergy.at(i))*keV,"K-Auger");
      EmitParticleForSourceDecay((fAugerEnergy.at(i))*keV,myEvent); // MHD 12 April 2013
      
	  EmissionForVacantShell(augerRecFrom[i], myEvent);	  
	  EmissionForVacantShell(augerEjecFrom[i], myEvent);	  
	}
      //      else
      //	G4cout << "K shell vacancy does not lead to X ray or Auger electron." << endl;
    }
  // For vacancies in shells L1, L2, L3 we only have data for x-ray emission
  // L1-Shell vacancy
  else if (shell == 1)
    {
      // sum up the X-ray intensities and Auger intensities
      totalXRayIntensity = 0.;      
      for(G4int i=0;i<fL1XRayIntensity.size();i++)
	  totalXRayIntensity += fL1XRayIntensity.at(i);

      // use random number to determine process
      xRayRandom = RandFlat::shoot(0.,100);
      if(totalXRayIntensity > xRayRandom)
	{   	    
	  G4int i = 0;
	  addedIntensity = fL1XRayIntensity.at(0);
	  while( addedIntensity < xRayRandom )
	    {
			i++;
			addedIntensity += fL1XRayIntensity.at(i);
	      if(i>fL1XRayIntensity.size())
			{
			  G4cout << "Reached end of array for L1 X ray intensities." << endl;
			  exit (1);
			}
	    }

	SetParticleType( "gamma" ) ; // MHD 12 April 2013
    //G4cout << "PrimaryGeneratorAction::EmissionForVacantShell L1-XRay" << endl;
    SetSourceRecord((fL1XRayEnergy.at(i))*keV,"L1-Xray");
	EmitParticleForSourceDecay((fL1XRayEnergy.at(i))*keV,myEvent); // MHD 12 April 2013
	
	EmissionForVacantShell(fL1XRayOrigin.at(i), myEvent);	  
	}
      //else
      //G4cout << "L1 shell vacancy does not lead to X ray or Auger electron." << endl;
    }
  // L2-Shell vacancy
  else if (shell == 2)
    {
      // sum up the X-ray intensities and Auger intensities
      totalXRayIntensity = 0.;      
      for(G4int i=0;i<fL2XRayIntensity.size();i++)
	totalXRayIntensity += fL2XRayIntensity.at(i);

      // use random number to determine process
      xRayRandom = RandFlat::shoot(0.,100);
      if(totalXRayIntensity > xRayRandom)
	{
	  G4int i = 0;
	  addedIntensity = fL2XRayIntensity.at(0);
	  while( addedIntensity < xRayRandom )
	    {
	      i++;
	      addedIntensity += fL2XRayIntensity.at(i);
	      if(i>fL2XRayIntensity.size())
			{
			  G4cout << "Reached end of array for L2 X ray intensities." << endl;
			  exit (1);
			}
	    }
	  SetParticleType( "gamma" ) ; // MHD 12 April 2013
	  //G4cout << "PrimaryGeneratorAction::EmissionForVacantShell L2-XRay" << endl;	  
	  EmitParticleForSourceDecay((fL2XRayEnergy.at(i))*keV,myEvent); // MHD 12 April 2013
	  SetSourceRecord((fL2XRayEnergy.at(i))*keV,"L2-Xray");
	  EmissionForVacantShell(fL2XRayOrigin.at(i), myEvent);	  
	}
      //      else
      //	G4cout << "L2 shell vacancy does not lead to X ray or Auger electron." << endl;
    }
  // L3-Shell vacancy
  else if (shell == 3)
    {
      // sum up the X-ray intensities and Auger intensities
      totalXRayIntensity = 0.;      
      for(G4int i=0;i<fL3XRayIntensity.size();i++)
	totalXRayIntensity += fL3XRayIntensity.at(i);

      // use random number to determine process
      xRayRandom = RandFlat::shoot(0.,100);
      if(totalXRayIntensity > xRayRandom)
	{
	  G4int i = 0;
	  addedIntensity = fL3XRayIntensity.at(0);
	  while( addedIntensity < xRayRandom )
	    {
	      i++;
	      addedIntensity += fL3XRayIntensity.at(i);
	      if(i>fL3XRayIntensity.size())
			{
			  G4cout << "Reached end of array for L3 X ray intensities." << endl;
			  exit (1);
			}
	    }

	  SetParticleType( "gamma" ) ; // MHD 12 April 2013
	  //G4cout << "PrimaryGeneratorAction::EmissionForVacantShell L3-XRay" << endl;
	  SetSourceRecord((fL3XRayEnergy.at(i))*keV,"L3-Xray");
	  EmitParticleForSourceDecay((fL3XRayEnergy.at(i))*keV,myEvent); // MHD 12 April 2013
	  
	  EmissionForVacantShell(fL3XRayOrigin.at(i), myEvent);	  
	}
      //      else
      //	G4cout << "L3 shell vacancy does not lead to X ray or Auger electron." << endl;
    }
  // No data for all other shells
  //  else
    // G4cout << "Vacancy left in shell " << shell << ", no data to continue." << G4endl;
   

}
// -----------------------------------------------------
// End of functions used for simulating a decay scheme (text file),
// primarily used for SPICE
// -----------------------------------------------------

void PrimaryGeneratorAction::SetRadioactiveBetaDecay( G4String path )
{
  particleType = "e-";
  this->particle = this->particleTable->FindParticle(particleType);
  particleGun->SetParticleDefinition(this->particle);
  //G4cout << " --> Particle type has been set to " << particleType << G4endl;

  this->radioactiveDecaySimulation = true;
  this->simulationDir = path;
  G4cout << " --> Radioactive simulation path: " << path << G4endl;

  this->myBetaParticle       = new BeamRequestBetaParticle(path);
  this->Z                    = myBetaParticle->GetZ();
  this->myGammaAndICParticle = new BeamRequestGammaAndIC(path, this->polarization);
  this->myXRay               = new BeamRequestXRay(path, this->Z, xray_input_kShell, xray_input_lShell, xray_input_mShell);
}

void PrimaryGeneratorAction::SetPolarization( G4double polar )
{
  // Polarization needs to be set BEFORE SetRadioactiveBetaDecay
  this->polarization = polar;
}

void PrimaryGeneratorAction::SetEmitBetaParticle( G4bool tf )
{
  this->emit_beta_flag = tf;
}

void PrimaryGeneratorAction::IncludeXRayInputFileKShell( G4bool tf )
{
  this->xray_input_kShell = tf;
}

void PrimaryGeneratorAction::IncludeXRayInputFileLShell( G4bool tf )
{
  this->xray_input_lShell = tf;
}

void PrimaryGeneratorAction::IncludeXRayInputFileMShell( G4bool tf )
{
  this->xray_input_mShell = tf;
}

void PrimaryGeneratorAction::SetRadioactiveDecayHalflife( G4double hl )
{
  this->halflife = hl/second;
  G4cout << " --> halflife = " << this->halflife << " seconds" << G4endl;
}

void PrimaryGeneratorAction::SetNumberOfRadioactiveNuclei( G4int num )
{
  this->numberOfNuclei = num;
}

void PrimaryGeneratorAction::ReadEnergyDistribution( G4String file )
{
	G4int ncolumns = 0;
	G4int nlines = 0;
	G4int lastline = -1;
	G4double doubleBuf;
	G4String buf, currentline;

	numberOfEvents = 0;
	std::ifstream infile;
	G4String inputfile;
	infile.open(file.c_str());
	if (!infile)
	{
		G4cout << "unable to find/open file" << G4endl;
		exit (1);
	}

	while(getline(infile, currentline))
	{
		if (!currentline.empty())
		{
			std::stringstream ss(currentline);
			ncolumns=0;
			while (ss >> buf)
			{
				doubleBuf = atof( buf.c_str() );
				if(ncolumns == 0)
				{
					if(buf == "#")
					{
						lastline = nlines;
					}
					else
					{
						energyDist.push_back( doubleBuf );					
					}
				}  
				else if(ncolumns == 1)
				{
					if(nlines == lastline)
					{
						numberOfEvents = (int)(doubleBuf);
					}
					else
					{
						weightDist.push_back( doubleBuf );
					}
				}
				else
				{
					G4cout << "Error: 1324431" << G4endl;
					exit (1);
				}
				ncolumns++;
			}
			nlines++;
		}
	}

	infile.close(); 
}

G4double PrimaryGeneratorAction::UniformRand48()
{
  G4double rand = drand48();
  return rand;
}

G4int PrimaryGeneratorAction::GetPrimaryParticleType()  // MHD : 12 April 2013
 {
 
 G4int type=-100 ;
 

   if( particleType == "gamma") 
  {
    type = 0; // A*10+Z
  }
  else if( particleType == "e-" ) 
	  {
		type = -1;  // 0*10 -1
	  }
	  else if( particleType == "e+" ) 
		  {
			type = +1; // 0*10 +1
		  }
		  else if( particleType == "proton") 
			  {
				type = +11; // 1*10 + 1
			  }
			  else if( particleType == "alpha") 
				  {
					type = +42; // 4*10 + 2
				  }
				  else
					  {
						type = -100; // default
					  } 
					  
	return type ; 				  
 }

// Messenger command to activate the kinematic effects
void PrimaryGeneratorAction::SetKinematicsActive( G4bool tf )
{
  fSimulateKinematics = tf;
}

// Messenger command to set beta value of heavy ion
void PrimaryGeneratorAction::SetKinematicsBetaValue( G4double beta )
{
  fBetaHeavyIon = beta;
  fBetaDefinedByUser = true;
}

// Messenger command to set kinetic energy of ion used in kinematic simulations
void PrimaryGeneratorAction::SetKinematicsIonEnergy( G4double value )
{
  fIonKineticEnergy = value;
  fIonEnergyDefinedByUser = true;
}


G4double PrimaryGeneratorAction::KinematicEnergyBroadening( G4double energy, G4Event* anEvent )
{
	// If the ion has been declared with '/DetSys/gun/ion' ...
	if ( ionDefined)
	{
	  // ... and no ion has been emitted this event, emit one
	  // At the very least, the direction of the ion is required in order to simulate kinematics
		if ( !ionEmittedThisEvent ) EmitIon( ionDefinitionZ, ionDefinitionA, ionDefinitionE, anEvent ); 
	}
	else 
	{
	  // Return an error if the ion has not been declared and return original energy
		G4cout << "ERROR: You need to define an ion with '/DetSys/gun/ion Z A E*' prior to using"
		       << " the KinematicEnergyBroadening function" << G4endl;
		return energy;
	}
	
	// Retrieve the theta and phi of the emitted particle & ion
  G4double thetaParticleIonFrame = direction.getTheta();
  G4double phiParticleIonFrame = direction.getPhi();
  if (phiParticleIonFrame < 0) phiParticleIonFrame += 2*3.141592654; // So that range is 0 to 2pi
  G4double thetaIonLabFrame = ionDirection.getTheta();
  G4double phiIonLabFrame = ionDirection.getPhi();
  if (phiIonLabFrame < 0) phiIonLabFrame += 2*3.141592654;
	
	// Using kinematics, theta in the lab frame is calculated starting from the
	// energy, theta & beta in the ion (source) frame.
	G4double gammaHeavyIon = 1 / sqrt(1 - pow(fBetaHeavyIon,2)); 
	
	G4double totalEnergyIonFrame = energy + fRestMassOfElectron;
  G4double gammaIonFrame = totalEnergyIonFrame / fRestMassOfElectron;
  G4double betaIonFrame = sqrt(1 - pow(gammaIonFrame, -2));
  G4double momentumIonFrame = sqrt(pow(energy,2) + ( 2 * energy * fRestMassOfElectron));	
  
  G4double totalEnergyLabFrame = 
           gammaHeavyIon * (totalEnergyIonFrame + fBetaHeavyIon * cos(thetaParticleIonFrame) * momentumIonFrame);
  G4double kineticEnergyLabFrame = totalEnergyLabFrame - fRestMassOfElectron;
  G4double gammaLabFrame = totalEnergyLabFrame / fRestMassOfElectron;
  G4double betaLabFrame = sqrt(1 - pow(gammaLabFrame, -2));
  G4double momentumLabFrame = 
           sqrt(pow(kineticEnergyLabFrame,2) + ( 2 * kineticEnergyLabFrame * fRestMassOfElectron));
           
	G4double betaLabThetaLab = (betaIonFrame * cos(thetaParticleIonFrame) + fBetaHeavyIon)
	                           /(1 + fBetaHeavyIon * betaIonFrame * cos(thetaParticleIonFrame));
  G4double thetaParticleLabFrame = acos(betaLabThetaLab/betaLabFrame);
   
  // Beautiful rotation matrices to find the new (x,y,z) components of the particle emission
  direction = G4ThreeVector(
                cos(phiIonLabFrame) * cos(thetaIonLabFrame) * cos(phiParticleIonFrame) * sin(thetaParticleLabFrame) - 
								sin(phiIonLabFrame) * sin(phiParticleIonFrame) * sin(thetaParticleLabFrame) +
								cos(phiIonLabFrame) * sin(thetaIonLabFrame) * cos(thetaParticleLabFrame),
                
                sin(phiIonLabFrame) * cos(thetaIonLabFrame) * cos(phiParticleIonFrame) * sin(thetaParticleLabFrame) +
								cos(phiIonLabFrame) * sin(phiParticleIonFrame) * sin(thetaParticleLabFrame) +
								sin(phiIonLabFrame) * sin(thetaIonLabFrame) * cos(thetaParticleLabFrame),
								
                cos(thetaIonLabFrame) * cos(thetaParticleLabFrame) - 
	              sin(thetaIonLabFrame) * cos(phiParticleIonFrame) * sin(thetaParticleLabFrame)
	            );
  
  return kineticEnergyLabFrame;      
}

// Emit ion for Kinematic Broadening simulation
void PrimaryGeneratorAction::EmitIon(G4int ionZ, G4int ionA, G4double ionE, G4Event* anEvent)
{
  // Random generator decides direction of ion
  // TODO The direction of the ion will actually form a distribution but for now it's directed at the S3
  G4double theta = 0;
  while (theta < 0.523598775 || theta > 1.047197551)
  {
    G4double costheta = 2.*UniformRand48()-1.0;
    G4double sintheta = sqrt( 1. - costheta*costheta );
    G4double phi      = (360.*deg)*UniformRand48();
    ionDirection = G4ThreeVector(sintheta*cos(phi), sintheta*sin(phi), costheta);
    theta = ionDirection.getTheta();
  }
  
  // Create a gun to fire the ion
  G4ParticleGun* ionGun = new G4ParticleGun(1);
  
  // Set ion attributes (direction, definition, position)
  ionGun->SetParticleMomentumDirection(ionDirection);  
  ionGun->SetParticleDefinition(this->particleTable->GetIon(ionZ,ionA,ionE));
	ionGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  
  // If beta has been defined by user, calculate appropriate kinetic energy of ion
  // or if kinetic energy has been defined by user, calculate beta
  if (fBetaDefinedByUser) fIonKineticEnergy = 100000 * ionA * fAtomicMassUnit * fBetaHeavyIon * pow(fSpeedOfLight, 2)/2;
  else if (fIonEnergyDefinedByUser) fBetaHeavyIon = sqrt( 2 * fIonKineticEnergy / (ionA * fAtomicMassUnit * pow(fSpeedOfLight, 2)) ) / 1000;
  else fIonKineticEnergy = 100*MeV; // Default
  
  // Fire ion  
  ionGun->SetParticleEnergy(fIonKineticEnergy);
  ionGun->GeneratePrimaryVertex(anEvent);
  ionEmittedThisEvent = true;
  delete ionGun;
}
 
       
void PrimaryGeneratorAction::SetSourceRecord( G4double energy , string process_name)
{
//G4cout<< energy << "   " << process_name << " size before "<< fMapSourceRecord.size()<<G4endl;  
fMapSourceRecord[energy].first = process_name;
fMapSourceRecord[energy].second++;
//G4cout<< "                                 size after "<< fMapSourceRecord.size()<<G4endl;  
}

void PrimaryGeneratorAction::WriteSourceRecord( void )
{

ofstream source_record ;
source_record.open("SourceRecord.txt");
source_record << "Energy (keV)    process    counts "<<endl;

for ( map<double,pair<string,int> >::iterator it = fMapSourceRecord.begin() ; it!=fMapSourceRecord.end() ; ++it )
    {
    //G4cout         <<setw(7)<<std::setprecision(5)<<it->first/keV<< "  " <<setw(7)<< it->second.first << "  "<<setw(7)<< it->second.second<<G4endl;
    source_record  <<setw(7)<<std::setprecision(5)<<it->first/keV<< "  " <<setw(7)<< it->second.first << "  "<<setw(7)<< it->second.second<<endl;
    }

	source_record.close();

}        
