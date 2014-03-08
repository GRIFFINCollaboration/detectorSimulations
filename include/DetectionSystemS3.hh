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

#ifndef DetectionSystemS3_h
#define DetectionSystemS3_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define AL_COL 0.5,0.5,0.5

class DetectionSystemS3
{
public:
  DetectionSystemS3();
  ~DetectionSystemS3();
  
  //------------------------------------------------//
  // logical and physical volumes
  //------------------------------------------------//
private:
  G4AssemblyVolume* assembly;
  G4AssemblyVolume* assemblyS3Ring01;
  G4AssemblyVolume* assemblyS3Ring02;
  G4AssemblyVolume* assemblyS3Ring03;
  G4AssemblyVolume* assemblyS3Ring04;
  G4AssemblyVolume* assemblyS3Ring05;
  G4AssemblyVolume* assemblyS3Ring06;
  G4AssemblyVolume* assemblyS3Ring07;
  G4AssemblyVolume* assemblyS3Ring08;
  G4AssemblyVolume* assemblyS3Ring09;
  G4AssemblyVolume* assemblyS3Ring10;
  G4AssemblyVolume* assemblyS3Ring11;
  G4AssemblyVolume* assemblyS3Ring12;
  G4AssemblyVolume* assemblyS3Ring13;
  G4AssemblyVolume* assemblyS3Ring14;
  G4AssemblyVolume* assemblyS3Ring15;
  G4AssemblyVolume* assemblyS3Ring16;
  G4AssemblyVolume* assemblyS3Ring17;
  G4AssemblyVolume* assemblyS3Ring18;
  G4AssemblyVolume* assemblyS3Ring19;
  G4AssemblyVolume* assemblyS3Ring20;
  G4AssemblyVolume* assemblyS3Ring21;
  G4AssemblyVolume* assemblyS3Ring22;
  G4AssemblyVolume* assemblyS3Ring23;
  G4AssemblyVolume* assemblyS3Ring24;
  SensitiveDetector* siDetS3Ring01_SD;
  SensitiveDetector* siDetS3Ring02_SD;
  SensitiveDetector* siDetS3Ring03_SD;
  SensitiveDetector* siDetS3Ring04_SD;
  SensitiveDetector* siDetS3Ring05_SD;
  SensitiveDetector* siDetS3Ring06_SD;
  SensitiveDetector* siDetS3Ring07_SD;
  SensitiveDetector* siDetS3Ring08_SD;
  SensitiveDetector* siDetS3Ring09_SD;
  SensitiveDetector* siDetS3Ring10_SD;
  SensitiveDetector* siDetS3Ring11_SD;
  SensitiveDetector* siDetS3Ring12_SD;
  SensitiveDetector* siDetS3Ring13_SD;
  SensitiveDetector* siDetS3Ring14_SD;
  SensitiveDetector* siDetS3Ring15_SD;
  SensitiveDetector* siDetS3Ring16_SD;
  SensitiveDetector* siDetS3Ring17_SD;
  SensitiveDetector* siDetS3Ring18_SD;
  SensitiveDetector* siDetS3Ring19_SD;
  SensitiveDetector* siDetS3Ring20_SD;
  SensitiveDetector* siDetS3Ring21_SD;
  SensitiveDetector* siDetS3Ring22_SD;
  SensitiveDetector* siDetS3Ring23_SD;
  SensitiveDetector* siDetS3Ring24_SD;
  
public:
  G4int Build(G4SDManager* mySDman);
  G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move,
		      G4int ringNumber, G4int nRadSeg, G4int detectorNumber);
  G4int PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move);
  
private:
  G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);
  
  G4LogicalVolume* S3InnerGuardRing_log;
  G4LogicalVolume* S3OuterGuardRing_log;
  G4LogicalVolume* siDetS3Ring01_log;
  G4LogicalVolume* siDetS3Ring02_log;
  G4LogicalVolume* siDetS3Ring03_log;
  G4LogicalVolume* siDetS3Ring04_log;
  G4LogicalVolume* siDetS3Ring05_log;
  G4LogicalVolume* siDetS3Ring06_log;
  G4LogicalVolume* siDetS3Ring07_log;
  G4LogicalVolume* siDetS3Ring08_log;
  G4LogicalVolume* siDetS3Ring09_log;
  G4LogicalVolume* siDetS3Ring10_log;
  G4LogicalVolume* siDetS3Ring11_log;
  G4LogicalVolume* siDetS3Ring12_log;
  G4LogicalVolume* siDetS3Ring13_log;
  G4LogicalVolume* siDetS3Ring14_log;
  G4LogicalVolume* siDetS3Ring15_log;
  G4LogicalVolume* siDetS3Ring16_log;
  G4LogicalVolume* siDetS3Ring17_log;
  G4LogicalVolume* siDetS3Ring18_log;
  G4LogicalVolume* siDetS3Ring19_log;
  G4LogicalVolume* siDetS3Ring20_log;
  G4LogicalVolume* siDetS3Ring21_log;
  G4LogicalVolume* siDetS3Ring22_log;
  G4LogicalVolume* siDetS3Ring23_log;
  G4LogicalVolume* siDetS3Ring24_log;
  
  //--------------------------------------------------------//
  // SPICE physical properties
  // OBS: crystal properties are public, others are private
  //--------------------------------------------------------//
private:
  G4String wafer_material;
  
  //-----------------------------//
  // parameters for the annular  //
  // planar detector crystal     //
  //-----------------------------//
public:
  G4double S3DetCrystalThickness;
  G4double S3DetCrystalOuterDiameter;
  G4double S3DetCrystalInnerDiameter;
  G4double S3DetRadialSegments;
  G4double S3DetPhiSegments;

  //-------------------------------//
  // parameters for the guard ring //
  //-------------------------------//
private:
  G4double S3DetGuardRingInnerDiameter;
  G4double S3DetGuardRingOuterDiameter;
   
    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//
private:
  G4int BuildSiliconWafer(G4int ringID);
  G4int BuildInnerGuardRing();
  G4int BuildOuterGuardRing();
  
  G4Tubs*             BuildCrystal(G4int myRingID);
};

#endif
