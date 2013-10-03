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

#ifndef DetectionSystemSpice_h
#define DetectionSystemSpice_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectionSystemSpice
{
  public:
    DetectionSystemSpice();
    ~DetectionSystemSpice();

    G4double GetDetector2Origin(DetectionSystemSpice*);
    G4double GetCrystalWidth();

    G4int Build() ; //G4SDManager* mySDman); 
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number); 

  private:
    G4ThreeVector getDirection(G4int);
    G4ThreeVector transTriangle(G4int);
    G4double getTheta(G4int);
    G4double getPhi(G4int);
    G4ThreeVector getstripDirection(G4int);

    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblySi;

//    SensitiveDetector* crystal_block_SD;

    // LogicalVolumes used in SPICEDetection()
    G4LogicalVolume* casing_log;
    G4LogicalVolume* crystal_block_log;

    G4AssemblyVolume* spiceDetector;

    ///////////////////////////////////////////////////////////////////
    // SPICE Properties
    ///////////////////////////////////////////////////////////////////

    G4int number_of_detectors;

    G4double octagon2[18][2];

    G4double crystal_length_x;   
    G4double crystal_length_y;   
    G4double crystal_length_z; 	 
    G4double crystal_dist_from_origin;

    G4double crystal_tube_length;
    G4double crystal_inner_radius;
    G4double crystal_outer_radius;
    G4double tube_to_target;

    G4double crystal_strip_thickness;
    G4double crystal_strip_width;
    G4double crystal_strip_length;
    G4double strip_to_target;
    G4double strip_rad_position;
    
    G4double detect_dist;

    G4String casing_material; 
    G4String wafer_material;

    G4double casing_width;    
    G4double casing_thickness;
    G4double casing_threshold;

    G4double detector2target;

    G4double theta;
    G4double phi;

    G4ThreeVector finalPlace;
  
    G4int BuildSiliconWafer();

};

#endif
