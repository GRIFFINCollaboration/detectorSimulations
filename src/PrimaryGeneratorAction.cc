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

  eventSum = 0;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // defaults
  energy = 1000.*keV;

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
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double randomClock = 0.0, lamda = 0.0;
  G4double deltaTimeInSeconds = 0.0;
  G4double posx, posy, posz ;
  G4int theBetaBranch, numberOfParticles, numberOfXRays;    

  //G4cout << "GeneratePrimaries" << G4endl;

  if(emissionSimulation) {
    G4double selectedEnergy = 0.00*keV;
    G4double monteCarloRand = UniformRand48();

    for(G4int i = 0 ; i < energyDist.size() ; i++ )
    {
      if(monteCarloRand <= monteCarlo[i])
      {
        selectedEnergy = energyDist[i]*1.0*keV;
        break;
      }
    }
    
    energy = selectedEnergy;

		if( isoRadOnBox )
		{
      posx = position.x();
      posy = position.y();
      posz = position.z();

      randBox_x = ( UniformRand48() * totalLengthOfBox_x ) - ( totalLengthOfBox_x / 2.0 ) ;
      randBox_y = ( UniformRand48() * totalLengthOfBox_y ) - ( totalLengthOfBox_y / 2.0 ) ;
      randBox_z = ( UniformRand48() * totalLengthOfBox_z ) - ( totalLengthOfBox_z / 2.0 ) ;
      
      direction = G4ThreeVector( randBox_x-posx, randBox_y-posy, randBox_z-posz) ;
	  }		

    else if(!directionSpecified)
    {
      // random direction
      G4double costheta = 2.*UniformRand48()-1.0;
      G4double sintheta = sqrt( 1. - costheta*costheta );
      G4double phi      = (360.*deg)*UniformRand48();
      direction = G4ThreeVector(sintheta*cos(phi), sintheta*sin(phi), costheta);
    }

    particleGun->SetParticleDefinition(this->particle);
    particleGun->SetParticlePosition(position);
    particleGun->SetParticleMomentumDirection(direction);
    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleTime(randomClock);

    eventID = anEvent->GetEventID();
    eventSum++;
    particleGun->GeneratePrimaryVertex(anEvent);
  }
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
  else 
  {
    
    if( isoRadOnBox )
		{
      posx = position.x();
      posy = position.y();
      posz = position.z();

      randBox_x = ( UniformRand48() * totalLengthOfBox_x ) - ( totalLengthOfBox_x / 2.0 ) ;
      randBox_y = ( UniformRand48() * totalLengthOfBox_y ) - ( totalLengthOfBox_y / 2.0 ) ;
      randBox_z = ( UniformRand48() * totalLengthOfBox_z ) - ( totalLengthOfBox_z / 2.0 ) ;
      
      direction = G4ThreeVector( randBox_x-posx, randBox_y-posy, randBox_z-posz) ;
	  }	
    else if(!directionSpecified)
    {
      // random direction
      //G4double costheta = 2.*UniformRand48()-1.0;
      G4double costheta = 2. * ( CLHEP::RandFlat::shoot() ) - 1.0 ;
      G4double sintheta = sqrt( 1. - costheta*costheta );
      //G4double phi      = (360.*deg)*UniformRand48();
      G4double phi      = (360.*deg) * ( CLHEP::RandFlat::shoot() ) ;
      direction = G4ThreeVector( sintheta * cos(phi) , sintheta * sin(phi) , costheta);
    }

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
  G4cout << " --> Particle type has been set to " << particleType << G4endl;
}

void PrimaryGeneratorAction::SetIonType( G4int Z, G4int A, G4double E )
{
  this->particle = this->particleTable->GetIon(Z,A,E*keV);
  particleGun->SetParticleDefinition(this->particle);
  G4cout << " --> Particle type has been set to an ion with Z=" << Z << " and A=" << A << " and an excitation energy of " << E << " keV." << G4endl;
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
