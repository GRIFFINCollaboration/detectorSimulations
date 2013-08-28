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

#ifndef DetectorHit_h
#define DetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include <string>

#include "G4Types.hh"

using namespace std;

class DetectorHit : public G4VHit
{
  public:
      DetectorHit();
     ~DetectorHit();
      DetectorHit(const DetectorHit&);

  public:
      inline void* operator new(size_t);
      inline void  operator delete(void*);

  public:
      void Draw();
      void Print();

  public:
      void SetDetNb     (const G4int det)          { detNb = det;       };  
      void SetSegNb     (const G4int seg)          { segNb = seg;       };
      void SetEdep      (const G4double de)        { edep = de;         };
      void SetTime      (const G4double ti)        { time = ti;         };
      void SetProcess   (const G4String pr)        { process = pr;      };
      void SetCollect   (const G4String col)       { collection = col;  };
      void SetPos       (const G4ThreeVector &xyz) { pos = xyz;         }; // position in the frame of the world
      
  public:
      inline G4int         GetDetNb()     const    { return detNb;      };
      inline G4int         GetSegNb()     const    { return segNb;      };
      inline G4double      GetEdep()      const    { return edep;       };
      inline G4double      GetTime()      const    { return time;       };  
      inline G4String      GetProcess()   const    { return process;    };
      inline G4String      GetCollect()   const    { return collection; };  
      inline G4ThreeVector GetPos()       const    { return pos;        };
      
  private:
      G4int         detNb;
      G4int         segNb;
      G4double      edep;
      G4double      time;
      G4String      process;
      G4String      collection;
      G4ThreeVector pos;
};

typedef G4THitsCollection<DetectorHit> DetectorHitCollection;

extern G4Allocator<DetectorHit> DetectorHitAllocator;

inline void* DetectorHit::operator new(size_t)
{
  void *aHit = (void *) DetectorHitAllocator.MallocSingle();
  return aHit;
}

inline void DetectorHit::operator delete(void *aHit)
{
  DetectorHitAllocator.FreeSingle((DetectorHit*) aHit);
}

#endif


