#include "BeamRequestXRay.hh"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits.h>

using namespace std;

BeamRequestXRay::BeamRequestXRay(G4String path, G4int Z, G4bool xray_input_kShell, G4bool xray_input_lShell, G4bool xray_input_mShell)
{
  srand48(time(NULL));

  this->path = path;
  this->k_shell_filename = "kShell.dat";
  this->l_shell_filename = "lShell.dat";
  this->m_shell_filename = "mShell.dat";

  this->Z = Z; // Z of daughter   
  this->angularDistribution = "Isotropic";

  inputXray[0].bools = xray_input_kShell;
  inputXray[1].bools = xray_input_lShell;
  inputXray[2].bools = xray_input_mShell;

  inputXray[0].vac = const_cast<char *>("Kvac");
  inputXray[1].vac = const_cast<char *>("Lvac");
  inputXray[2].vac = const_cast<char *>("Mvac");
  inputXray[3].vac = const_cast<char *>("Nvac");    

  ///////////////////////////////////////////////////////////////////
  // X-ray input tables
  ///////////////////////////////////////////////////////////////////
  if(xray_input_kShell) {
    ReadXRayData(this->k_shell_filename); 
    bool_xrays = true;
  }
  if(xray_input_lShell) {
    ReadXRayData(this->l_shell_filename); 
    bool_xrays = true;      
  }
  if(xray_input_mShell) {
    ReadXRayData(this->m_shell_filename); 
    bool_xrays = true;      
  }
       
  if (bool_xrays) {
    for(G4int i = 0 ; i < ARRAYDIM ; i++) {
      for(G4int j = 0 ; j < COLSGDATA ; j++) {
        for(G4int k = 0 ; k < MAXSHELLS ; k++) {
          this->inputXrayFileData[i][j][k] = ReturnXRayData(i,j,k);
        }
      }
    }
  }     
  ClearMyParticleStruct();
  this->particle.clear();  
}

BeamRequestXRay::~BeamRequestXRay()
{
}

void BeamRequestXRay::ReadXRayData(G4String filename)
{
  G4String bufString = this->path + filename;
  G4int ncolumns = 0;
  G4int numberLines = 0;
  G4int slice = -1;
  G4double doubleBuf;
  G4String buf, currentline;

  ifstream infile;
  string inputfile;
  infile.open(bufString.c_str());

  for(G4int i = 0 ; i < ARRAYDIM ; i++) {     // set inputXrayFileData to zero
    for(G4int j = 0 ; j < COLSGDATA ; j++) {
      for(G4int k = 0 ; k < MAXSHELLS ; k++) {
        this->inputXrayFileData[i][j][k] = 0.0;
      }
    }
  }

  numberLines = 0;

  if (!infile) {
    cerr << "unable to find/open inputFile" << endl;
    exit(1);
  }  

  while(getline(infile, currentline)) {
    if (!currentline.empty()) {
      stringstream ss(currentline);
      ncolumns=0;
      while (ss >> buf) {
        if (numberLines == 0) { // looks at the header
          slice = GetXRayShellNumberFromFile(buf);    // Get Slice 
        }
        else { // gets data
          if (numberLines > 0 ) {
            doubleBuf = atof( buf.c_str() );
            this->inputXrayFileData[numberLines-1][ncolumns][slice] = doubleBuf;
          }       
          else {
            cout << "Error in X-ray energy and intensity input file!" << endl; 
          }  
        }     
        ncolumns++;
      }
      numberLines++;
    }
  }

  infile.close();   

  if (numberLines >= ARRAYDIM) {  // Check
    G4cout << "Error in length of Data, exceeded the expected arrays" << G4endl;
    exit(1);
  }
}

G4int BeamRequestXRay::GetXRayShellNumberFromFile(G4String buf)
{
  G4int output;
  char *compare;
  char *pBuf = const_cast<char*> ( buf.c_str() );   
  char *BUF1 = strtok(pBuf, "_");
  char *BUF2 = strtok(NULL, "_");    
  char *BUF3 = strtok(NULL, "\0");

  compare = const_cast<char *>(BUF1);

  for (int i = 0; i < MAXSHELLS; i++) {
    if(!strcmp(inputXray[i].vac,compare)) {
      output = i;
      break;
    }
    else {
      if(i == MAXSHELLS - 1) {
        G4cout << "Error finding Shell Vacancies" << G4endl;
        output = -1; 
        exit(1);       
      }
    }
  }
  return output;
}

G4int BeamRequestXRay::GetXRayShellNumber(G4String shell)
{
  G4int output;

  if(shell == "K") {
    output = 1;
  }
  else if(shell == "L") {
    output = 2;  
  }
  else if(shell == "M") {
    output = 3;  
  }
  else if(shell == "N") {
    output = 4;  
  }   
  else {
    G4cout << "Error: 5646815 finding Shell Vacancies" << G4endl;
    exit(1);
  }   
  return output;
}

void BeamRequestXRay::GenerateParticles(G4String shell)
{
  G4double value, random, outputEnergy;
  
  this->particle.clear();  
  ClearMyParticleStruct();  

  G4int shellNumber = GetXRayShellNumber(shell);

  if(shellNumber < 1 || shellNumber > 4) {
    G4cout << "Error: 3214893" << G4endl;
    G4cout << "shell value for xray energy is out_uLong of range... " << G4endl;
    exit(1);  
  }
  if(inputXray[shellNumber-1].bools && inputXrayFileData[0][0][shellNumber-1] != 0) {
    for(G4int j = 0 ; j < COLSGDATA ; j=j+2) {
      random = drand48();
      value = 0;
      outputEnergy = 0;
      for(G4int i = 0 ; i < ARRAYDIM ; i++) {
        if (inputXrayFileData[i][j][shellNumber-1] == 0) {
          outputEnergy = 0;
          break;
        }
        else {
          value = value + double(inputXrayFileData[i][j+1][shellNumber-1]/100.0);
          if(random <= value) {
            outputEnergy = inputXrayFileData[i][j][shellNumber-1];
            break;
          }
        }
      }
      if(outputEnergy != 0) {
        this->my_particle.energy = outputEnergy;
        this->my_particle.shellNumber = shellNumber;  
        particle.push_back(my_particle);  
        ClearMyParticleStruct();  
      }
    }    
  }
  else { // Use Average X-Ray Energies
    if(shellNumber == 4) {
      outputEnergy = 0.0;
    }
    else {
      for (G4int z = 0; z < XRAYROWS ; z++) {
        if ( averageXRayEnergies[z][0] == this->Z ) {
          outputEnergy = averageXRayEnergies[z][shellNumber];
          //G4cout << "Z = " << this->Z << " shellNumber = " << shellNumber << " outputEnergy = " << outputEnergy << G4endl; 
          break;
        }
        if ( z == XRAYROWS-1 ) {
          G4cout << "Error: 32423933" << G4endl;
          G4cout << "can't find xray energy... " << G4endl;
          exit(1);  
        }   
      }  
    }  
    this->my_particle.energy = outputEnergy;
    this->my_particle.shellNumber = shellNumber;  
    particle.push_back(my_particle);  
    ClearMyParticleStruct();  
  }
}

void BeamRequestXRay::ClearMyParticleStruct()
{
  this->my_particle.energy = 0;
  this->my_particle.shellNumber = 0;
}

void BeamRequestXRay::AddBeamParticle(G4ParticleGun* beamParticle, G4int pNum)
{
  beamParticle -> SetParticleDefinition(G4ParticleTable::GetParticleTable() -> FindParticle( GetParticleName() ));
  beamParticle -> SetParticleEnergy( (this->particle[pNum].energy)/(1000.0) );  // in MeV!
  beamParticle -> SetParticlePosition( GetParticlePosition() );
  beamParticle -> SetParticleMomentumDirection( GetParticleDirection() ); 
}

G4String BeamRequestXRay::GetParticleName()
{
  G4String particleName = "gamma";
  return particleName;
}

G4ThreeVector BeamRequestXRay::GetParticlePosition()
{
  G4ThreeVector particlePosition = G4ThreeVector(0, 0, 0);
  return particlePosition;
}

G4ThreeVector BeamRequestXRay::GetParticleDirection()
{
  G4ThreeVector particleDirection;

  if(this->angularDistribution == "Isotropic") {
    // Randomize the direction over 4pi
    G4double rand1, rand2, phi, theta;
    rand1 = drand48(); 
    rand2 = drand48();

    phi = rand1*2.0*M_PI;
    G4double angle = M_PI;   //for isotropic scattering
    theta = acos(1.0 - (rand2*(1.0 - cos(angle))));

    particleDirection = G4ThreeVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  }
  else {
    G4cout << "Error 12312200 : Xray Angular Distribution" << G4endl;
    exit(1);
  }

  return particleDirection;
}

void BeamRequestXRay::PrintAllGeneratedParticles()
{
  G4cout << "==============================================New Decay======================================================" << G4endl;

  for( int i = 0 ; i < particle.size() ; i++ ) {
    G4cout << "Energy       = " << particle[i].energy << "\t\t";
    G4cout << "Shell Number = " << particle[i].shellNumber << G4endl;
    G4cout << "-------------------------------------------------------------------------------------------------------------" << G4endl;  
  }
  G4cout << "=============================================================================================================" << G4endl;
}
