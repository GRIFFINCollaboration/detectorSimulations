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

#ifndef RunAction_h
#define RunAction_h 1

//#include <fstream>
//#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4UserRunAction.hh"

//using namespace std;

class G4Run;
class HistoManager ;
//class RunActionMessenger;

class RunAction : public G4UserRunAction
{
  public:
    RunAction( HistoManager* );
    virtual ~RunAction();

    void BeginOfRunAction (const G4Run*);
    void EndOfRunAction   (const G4Run*);

//    void EnableWrite      ( G4bool );
//    void EnableWriteHex   ( G4bool );

//    void SetRunNumber     ( G4int );

//    inline void  LMwrite(const G4String &theString) { outFileLMD << theString;   };
//    inline void  LMwrite(const char     *theChars)  { outFileLMD << theChars;    };
//    void  HEXwrite(G4int ndet, G4double ener, G4double time, G4ThreeVector pos1, G4int nseg, G4String process, G4String collection);

//    void  HEXwriteStartFlag();
//    void  HEXwriteEndFlag();

//    void  TimeFlagStart();
//    void  TimeFlagEnd();
//    void  TimeEntry(G4double time_value);

//    void  ParticleInfoFlagStart();
//    void  ParticleInfoFlagEnd();
//    void  ParticleInfoEntry(G4double time_value);

//    void  SetTimeUnit(G4String newTimeUnit);

//    inline G4bool DoWriteLMD()                         { return writeLMD;           };
//    inline G4bool DoWriteHEX()                         { return writeHEX;           };
//    inline G4int  GetRunNumber()                       { return runNumber;          };
    
  private:
  
  	HistoManager* histoManager ;
//    void  WriteHeader ();
//    G4int OpenLMFile  ();  
    
//    std::ofstream outFileLMD;
    
//    G4int    runNumber;
//    G4bool   writeLMD;
//    G4bool   writeHEX;
//    G4String outFileName;
//    G4String timeUnit;

//    RunActionMessenger* myMessenger;
};

#endif





