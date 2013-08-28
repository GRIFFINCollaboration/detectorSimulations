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

#ifndef BeamRequestXRay_H
#define BeamRequestXRay_H 1

//#include "globals.hh"
//#include "GlobalDefines.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
//#include "SearchInputFile.hh"
#include "TableOfIsotopesData.hh"

#include <vector>

#define ARRAYDIM 512
#define MAXDECAYS 10
#define MAXSHELLS 4
#define MAXLEVELS 200
#define MAXQVALUE 100000
#define GCVTMAG 10
#define COLSGDATA 16
#define CHARLENGTH 1000
#define MAXBETABRANCHES 200
#define VECSIZE 10000

#define GAMMADIM 10
#define SPINDIM 10
#define THETASTEP 100

using namespace std;

class BeamRequestXRay 
{    
  public:
    BeamRequestXRay(G4String path, G4int Z, G4bool xray_input_kShell, G4bool xray_input_lShell, G4bool xray_input_mShell);
    ~BeamRequestXRay();

    void GenerateParticles(G4String shell);
    void AddBeamParticle(G4ParticleGun* beamParticle, G4int pNum);
    void PrintAllGeneratedParticles();
    G4int GetNumberOfGeneratedParticles() {return particle.size();}

  private:
    G4String GetParticleName();
    G4ThreeVector GetParticlePosition();
    G4ThreeVector GetParticleDirection();
    G4double ReturnXRayData(int i, int j, int k) {return inputXrayFileData[i][j][k];}
    void ReadXRayData(G4String filename);
    G4int GetXRayShellNumber(G4String shell);
    G4int GetXRayShellNumberFromFile(G4String buf);
    void ClearMyParticleStruct();
    G4double inputXrayFileData[ARRAYDIM][COLSGDATA][MAXSHELLS];    

    struct shellLabelsStruct 
    {
      char *shell;   
      char *vac;     
      G4bool bools;    
    } inputXray[MAXSHELLS]; 
    
    struct particleStruct 
    {
      G4double energy;
      G4int shellNumber;
    } my_particle;     
   
    vector<particleStruct> particle;         

    G4String angularDistribution;
    G4String path;
    G4String k_shell_filename;
    G4String l_shell_filename;
    G4String m_shell_filename;
    G4bool bool_xrays;    
    G4int Z;
};

#endif
