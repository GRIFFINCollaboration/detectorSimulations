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
#define PEEK_COL 0.5, 0.5, 0.0


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
	G4AssemblyVolume* assemblyS3Ring[24];
	G4VPhysicalVolume* s3_mount_phys;
  

private:
  G4ThreeVector GetDirectionXYZ(G4double theta, G4double phi);
  
  G4LogicalVolume* S3InnerGuardRing_log;
  G4LogicalVolume* S3OuterGuardRing_log;
  G4LogicalVolume* siDetS3Ring_log[24];
  G4LogicalVolume* s3_mount_log; 
   
  //--------------------------------------------------------//
  // SPICE physical properties
  // OBS: crystal properties are public, others are private
  //--------------------------------------------------------//
private:
  G4String wafer_material;
  G4String s3_mount_material; 
  
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

  // -------------------------
  // Dimensions of Si-CD Mount
  // -------------------------
  G4double s3_mount_length;
  G4double s3_mount_thickness;
  G4double s3_active_radius;
  G4double s3_mount_chamfer;
  G4double s3_mount_centre_to_chamfer;
  
  // -------------------------
  // Copy number
  // -------------------------
  G4int s3MountCopyNumber;
       
    //------------------------------------------------//
    // internal methods in Build()
    //------------------------------------------------//
private:
  G4int BuildSiliconWafer(G4int ringID);
  G4int BuildInnerGuardRing();
  G4int BuildOuterGuardRing();
  void BuildS3Mount();  

	   
  G4Tubs*  BuildCrystal(G4int myRingID);
 
 public:
  G4int Build(); 
  G4int PlaceDetector( G4LogicalVolume* exp_hall_log, G4ThreeVector move,  G4double angle_offset, G4int ringNumber, G4int nRadSeg, G4int detectorNumber);
  G4int PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move);
  void  PlaceS3Mount(G4LogicalVolume* exp_hall_log, G4ThreeVector move,  G4double angle_offset); 

};

#endif
