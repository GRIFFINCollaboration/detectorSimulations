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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BeamRequestBetaParticle_H
#define BeamRequestBetaParticle_H 1

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

const int MAXLEVELS       = 200;
const int MAXQVALUE       = 100000;
const int MAXBETABRANCHES = 200;
const int VECSIZE         = 10000;

using namespace std;

class BeamRequestBetaParticle 
{
  public:
    BeamRequestBetaParticle(G4String path);
    ~BeamRequestBetaParticle();

    void AddBeamParticle(G4ParticleGun* beamParticle, G4int betaBranch, G4double timeSeconds);
    const G4int GetZ() {return this->myZ;};

  private:
    G4String GetParticleName();
    G4double GetParticleEnergy(G4int betaBranch);
    G4ThreeVector GetParticlePosition();
    G4ThreeVector GetParticleDirection(G4double timeSeconds);

    void BetaParticleDistrib();
    void GenerateWThetaAngularDistribution();
    G4double UniformRand48();

    G4String inputFile;
    G4String path;
    G4String angularDistribution;
    G4String betaFile;
    G4String stringT1;
    G4String stringT2;

    G4int numberOfBetaBranches;
    G4int betaBranch;
    G4int numberOfEvents;
    G4int betaPlusMinusType;
    G4int maxBranch;
    G4int myZ;

    G4double degreeOfPolarization;
    G4double betaAsymmetry;
    G4double spinTempBeta;
    G4double timeT1;
    G4double timeT2;
    G4double betaDistrib[MAXQVALUE][2][MAXLEVELS];
    G4double wTheta[VECSIZE];
};

#endif
