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

#ifndef ApparatusSpiceTargetChamber_h
#define ApparatusSpiceTargetChamber_h 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class MagneticField;

#define DELRIN_COL 0.4,0.2,1.0
#define AL_COL 0.5,0.5,0.5
#define PB_COL 0.3,0.3,0.3
#define NDFEB_COL 0.7,0.3,0.3

///////////////////////////////////////////////////////////////////////
// ApparatusSpiceTargetChamber
///////////////////////////////////////////////////////////////////////
class ApparatusSpiceTargetChamber
{
  public:
    ApparatusSpiceTargetChamber();
    ~ApparatusSpiceTargetChamber();

  public:
    void Build(G4LogicalVolume*);
  
  private:
    G4LogicalVolume* expHallLog;

  private:
//    MagneticField* solMagField;   	// pointer to the magnetic field

  private:
    // LogicalVolumes used in ApparatusSpiceTargetChamber()

    G4LogicalVolume* target_chamber_cylinder_log;
    G4LogicalVolume* target_chamber_sphere_log;
    G4LogicalVolume* transport_chamber_al_plate_log;
    G4LogicalVolume* transport_chamber_wide_cylinder_log;
    G4LogicalVolume* transport_chamber_narrow_cylinder_log;
    G4LogicalVolume* photon_shield_log;
    G4LogicalVolume* transport_magnet_log;
    G4LogicalVolume* inner_transport_magnet_log;
    G4LogicalVolume* inner_transport_magnet_cover_log;
    G4LogicalVolume* collector_magnet_log;

    G4LogicalVolume* det_cham_end_plate_log;
    G4LogicalVolume* det_cham_plate_log;
    G4LogicalVolume* det_cham_tube_log;

    G4LogicalVolume* target_chamber_vacuum_log;
    G4LogicalVolume* target_chamber_sph_vacuum_log;
    G4LogicalVolume* transport_wide_vacuum_log;
    G4LogicalVolume* transport_narrow_vacuum_log;
    G4LogicalVolume* al_plate_vacuum_log;
    G4LogicalVolume* det_cham_plate_vacuum_log;
    G4LogicalVolume* det_cham_tube_vacuum_log;


    
  private:
    // Physical Volumes used in ApparatusSpiceTargetChamber()

    G4VPhysicalVolume* target_chamber_cylinder_phys;
    G4VPhysicalVolume* target_chamber_sphere_phys;
    G4VPhysicalVolume* transport_chamber_al_plate_phys;
    G4VPhysicalVolume* transport_chamber_wide_cylinder_phys;
    G4VPhysicalVolume* transport_chamber_narrow_cylinder_phys;
    G4VPhysicalVolume* photon_shield_phys;
    G4VPhysicalVolume* transport_magnet_phys;
    G4VPhysicalVolume* inner_transport_magnet_phys;
    G4VPhysicalVolume* inner_transport_magnet_cover_phys;
    G4VPhysicalVolume* collector_magnet_phys;

    G4VPhysicalVolume* det_cham_end_plate_phys;
    G4VPhysicalVolume* det_cham_plate_phys;
    G4VPhysicalVolume* det_cham_tube_phys;

    G4VPhysicalVolume* target_chamber_vacuum_phys;
    G4VPhysicalVolume* target_chamber_sph_vacuum_phys;
    G4VPhysicalVolume* transport_wide_vacuum_phys;
    G4VPhysicalVolume* transport_narrow_vacuum_phys;
    G4VPhysicalVolume* al_plate_vacuum_phys;
    G4VPhysicalVolume* det_cham_plate_vacuum_phys;
    G4VPhysicalVolume* det_cham_tube_vacuum_phys;

  private: 
    ///////////////////////////////////////////////////////////////////
    // ApparatusSpiceTargetChamberCylinder Properties 
    ///////////////////////////////////////////////////////////////////

    // Materials
    G4String target_chamber_cylinder_material;
    G4String target_chamber_sphere_material;
    G4String transport_chamber_al_plate_material;
    G4String transport_chamber_wide_cylinder_material;
    G4String transport_chamber_narrow_cylinder_material;
    G4String photon_shield_material;
    G4String transport_magnet_material;
    G4String inner_transport_magnet_material;
    G4String inner_transport_magnet_cover_material;
    G4String collector_magnet_material;
    G4String vacuum_material;
    G4String det_cham_end_plate_material;
    G4String det_cham_plate_material;
    G4String det_cham_tube_material;

    // Dimensions
    G4double target_chamber_cylinder_length;
    G4double target_chamber_inner_radius;
    G4double target_chamber_outer_radius;
    G4double Al_plate_thickness;
    G4double Al_plate_inner_radius;
    G4double Al_plate_outer_radius;
    G4double trans_chamber_wide_cylinder_length;
    G4double trans_chamber_wide_cylinder_inner_radius;   
    G4double trans_chamber_wide_cylinder_outer_radius;
    G4double trans_chamber_narrow_cylinder_length;           
    G4double trans_chamber_narrow_cylinder_inner_radius;    
    G4double trans_chamber_narrow_cylinder_outer_radius;
    G4double photon_shield_inner_radius;   
    G4double photon_shield_tip_outer_radius;
    G4double photon_shield_base_outer_radius;
    G4double photon_shield_length;
    G4double transport_magnet_length;
    G4double transport_magnet_inner_radius;
    G4double transport_magnet_outer_radius;
    G4double inner_transport_magnet_length;
    G4double inner_transport_magnet_inner_radius;
    G4double inner_transport_magnet_outer_radius;
    G4double inner_transport_magnet_cover_length;
    G4double inner_transport_magnet_cover_outer_radius;
    G4double inner_transport_magnet_cover_inner_radius;
    G4double collector_magnet_length_x;
    G4double collector_magnet_length_y_inner;
    G4double collector_magnet_length_y_outer;
    G4double collector_magnet_length_z;
    G4double det_cham_end_plate_thickness;
    G4double det_cham_end_plate_inner_radius;
    G4double det_cham_end_plate_outer_radius;
    G4double det_cham_plate_thickness;
    G4double det_cham_plate_inner_radius;
    G4double det_cham_plate_outer_radius;
    G4double det_cham_tube_thickness;
    G4double det_cham_tube_inner_radius;
    G4double det_cham_tube_outer_radius;

    // Positions
    G4double target_chamber_start_position;
    G4double al_plate_position;
    G4double trans_chamber_wide_cylinder_position;
    G4double trans_chamber_narrow_cylinder_position;
    G4double photon_shield_position;
    G4double transport_magnet_position;
    G4double inner_transport_magnet_position;
    G4double inner_transport_magnet_cover_position;
    G4double collector_magnet_position_z;
    G4double collector_magnet_position_y;
    G4double collector_magnet_axial_separation;

    G4double det_cham_end_plate_position;
    G4double det_cham_plate_position;
    G4double det_cham_tube_position;


    //Cross Section
    G4double target_chamber_cross_section;
    
    // Field Strengths
    G4double solMagFieldStrength;


  private: 
    // internal methods for ConstructApparatusSpiceTargetChamber() for Building
    void BuildApparatusSpiceTargetChamberCylinder();
    void BuildApparatusSpiceTargetChamberSphere();
    void BuildTransportChamberAlPlate();
    void BuildTransportChamberWideCylinder();
    void BuildTransportChamberNarrowCylinder();
    void BuildPhotonShield();
    void BuildTransportMagnets();
    void BuildInnerTransportMagnets();
    void BuildInnerTransportMagnetsCover();
    void BuildCollectorMagnet();

    void BuildDetectorChamberAlEndPlate();
    void BuildDetectorChamberAlPlate();
    void BuildDetectorChamberAlTube();   

  private:
    // Internal Functions 
    G4RotationMatrix* RotateMagnets(G4int);
    G4ThreeVector TranslateMagnets(G4int, G4double, G4double);
    
  private: 
    // internal methods for ConstructApparatusSpiceTargetChamber() for Placing 
    void PlaceApparatusSpiceTargetChamberCylinder();   
    void PlaceApparatusSpiceTargetChamberSphere();
    void PlaceTransportChamberAlPlate();
    void PlaceTransportChamberWideCylinder();
    void PlaceTransportChamberNarrowCylinder();
    void PlacePhotonShield();
    void PlaceTransportMagnets();
    void PlaceInnerTransportMagnets();
    void PlaceInnerTransportMagnetsCover();
    void PlaceCollectorMagnet(G4int);

    void PlaceDetectorChamberAlEndPlate();
    void PlaceDetectorChamberAlPlate();
    void PlaceDetectorChamberAlTube();

    
};

#endif
