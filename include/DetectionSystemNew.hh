
#ifndef DetectionSystemNew_h
#define DetectionSystemNew_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
//class MagneticField;



// define custom colours for visualisation
#define CHAMBER_COL 0.5,0.5,0.0
#define LAYER1_COL 0.7,0.0,0.8
#define LAYER2_COL 0.7,0.2,0.8
#define LAYER3_COL 0.7,0.4,0.8
#define NDFEB_COL 0.7,0.3,0.3
#define MAGCOV_COL 0.3,0.3,0.3



class DetectionSystemNew
{


    public:
    DetectionSystemNew();
    ~DetectionSystemNew();


    public:
    void Build(G4LogicalVolume*);


    private:
    G4LogicalVolume* expHallLog;


    private:

    ////////////////////////////////////////////////
    // Logical Volumes used in NewDetectionSystem //
    ////////////////////////////////////////////////

    G4LogicalVolume* final_chamber_log;
    G4LogicalVolume* magnet_log;
    G4LogicalVolume* magnet_cover_log;
    G4LogicalVolume* photon_shield_layer_one_log;
    G4LogicalVolume* photon_shield_layer_two_log;
    G4LogicalVolume* photon_shield_layer_three_log;

    private:

    /////////////////////////////////////////////////
    // Physical Volumes used in NewDetectionSystem //
    /////////////////////////////////////////////////

    G4VPhysicalVolume* final_chamber_phys;
    G4VPhysicalVolume* magnet_phys;
    G4VPhysicalVolume* magnet_cover_phys;
    G4VPhysicalVolume* photon_shield_layer_one_phys;
    G4VPhysicalVolume* photon_shield_layer_two_phys;
    G4VPhysicalVolume* photon_shield_layer_three_phys;

    private:

    ///////////////////////////////////////////
    // Properties used in NewDetectionSystem //
    ///////////////////////////////////////////

    ////////////////////////
    // Magnet properties: //
    ////////////////////////


    G4int NUMBER_OF_MAGNETS;
//    G4double zOffset;

    ////////////////
    // Materials: //
    ////////////////

    G4String magnet_material;
    G4String magnet_cover_material;
    G4String chamber_material;
    G4String photon_shield_layer_one_material;
    G4String photon_shield_layer_two_material;
    G4String photon_shield_layer_three_material;

    /////////////////
    // Dimensions: //
    /////////////////

    //////////////
    // Chamber: //
    //////////////
    
    G4double chamber_hemisphere_inner_radius;
    G4double chamber_hemisphere_outer_radius;
    G4double chamber_cylinder_inner_radius;
    G4double chamber_cylinder_outer_radius;
    G4double chamber_cylinder_height;
    G4double chamber_cut_out_inner_radius;
    G4double chamber_cut_out_outer_radius;
    G4double chamber_cut_out_height;

    ////////////////////
    // Photon Shield: //
    ////////////////////

    G4double photon_shield_layer_one_thickness;
    G4double photon_shield_layer_two_thickness;
    G4double photon_shield_layer_three_thickness;

    G4double photon_shield_target_distance ;
	G4double detector_full_width ;
	G4double detector_target_distance;

    G4double photon_shield_beam_cut_outer_radius;
    G4double photon_shield_beam_cut_inner_radius;
    G4double photon_shield_beam_cut_height;


    //////////////
    // Magnets: //
    //////////////

    G4double plate_one_one_thickness;
    G4double plate_one_one_length;
    G4double plate_one_one_height;
    G4double plate_one_one_lower_height;
    G4double plate_one_edge_x;
    G4double plate_one_two_thickness;
    G4double plate_one_two_length;
    G4double plate_one_two_height;

    G4double plate_two_one_thickness;
    G4double plate_two_one_length;
    G4double plate_two_one_height;
    G4double plate_two_one_edge_x;
    G4double plate_two_two_thickness;
    G4double plate_two_two_length;
    G4double plate_two_two_height;

    G4double plate_three_one_thickness;
    G4double plate_three_one_length;
    G4double plate_three_one_height;
    G4double plate_three_one_edge_x;
    G4double plate_three_two_thickness;
    G4double plate_three_two_length;
    G4double plate_three_two_height;

    G4double distance_from_target;
    G4double cutting_box_angle;
	G4int no_magnet_layer;
	G4double magnet_cover_thickness;
  

    private:

    ///////////////////////////////////////////////////////////////////
    // internal methods and functions in NewDetectionSystem::Build() //
    ///////////////////////////////////////////////////////////////////

    // methods
    void BuildChamber();
    void BuildPhotonShield();
    void BuildMagnet();
    void BuildMagnetCover();

    void PlaceChamber();
    void PlacePhotonShield();
    void PlaceMagnet(G4int);
    void PlaceMagnetCover(G4int);

    // functions
    G4RotationMatrix* RotateMagnets(G4int);
    G4ThreeVector TranslateMagnets(G4int,G4double,G4double);


};

#endif
