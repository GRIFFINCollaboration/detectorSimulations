#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
//#include "SensitiveDetector.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

//#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemSpice.hh"

///////////////////////////////////////////////////////////////////////
// The ::DetectionSystemSpice constructor instatiates all the Logical 
// and Physical Volumes used in the SPICE geometery, and the 
//  ::DetectionSystemSpice  destructor deletes them from the stack 
// when they go out of scope
///////////////////////////////////////////////////////////////////////
DetectionSystemSpice::DetectionSystemSpice() :
   // LogicalVolumes
  crystal_block_log(0), 
  casing_log(0)
{ 
  /////////////////////////////////////////////////////////////////////
  // SPICE Physical Properties 
  /////////////////////////////////////////////////////////////////////
  this->number_of_detectors                        = 4;

  this->casing_material                            = "Aluminum";
  this->wafer_material                             = "Silicon"; 

  //Careful, these are cartesian dimensions of a triangle NOT length of sides (used for hexagonal plane
  this->crystal_length_x 	                   = 5.00*mm;   
  this->crystal_length_y 	                   = 10.*mm; 
  this->crystal_length_z 	                   = 10.*mm; 

  this->crystal_dist_from_origin                   = 11*mm;  //From calculation based on spacing between detector triplette groups
  this->casing_width                               = 1.0*mm;
  this->casing_thickness                           = 1.0*mm;
  this->casing_threshold                           = 0.0*mm;       //Typically not used, but implemented in code just in case
  this->detector2target                            = -50*mm;

  //These Dimensions are used for cylindrical detector.
  this->crystal_tube_length                        =30.*mm;
  this->crystal_inner_radius                       =13.*mm;
  this->crystal_outer_radius                       =18.*mm;
  this->tube_to_target                             =-50.*mm;

  this->octagon2[0][0] = 0.0;
  this->octagon2[0][1] = 0.0;
  this->octagon2[1][0] = 180.0;
  this->octagon2[1][1] = 0.0;
  this->octagon2[2][0] = 0.0;
  this->octagon2[2][1] = 0.0;
  this->octagon2[3][0] = 0.0;
  this->octagon2[3][1] = 60.0;
  this->octagon2[4][0] = 180.0;
  this->octagon2[4][1] = 60.0;
  this->octagon2[5][0] = 0.0;
  this->octagon2[5][1] = 60.0;
  this->octagon2[6][0] = 0.0;
  this->octagon2[6][1] = 120.0;
  this->octagon2[7][0] = 180.0;
  this->octagon2[7][1] = 120.0;
  this->octagon2[8][0] = 0.0;
  this->octagon2[8][1] = 120.0;
  this->octagon2[9][0] = 0.0;
  this->octagon2[9][1] = 180.0;
  this->octagon2[10][0] = 180.0;
  this->octagon2[10][1] = 180.0;
  this->octagon2[11][0] = 0.0;
  this->octagon2[11][1] = 180.0;
  this->octagon2[12][0] = 0.0;
  this->octagon2[12][1] = 240.0;
  this->octagon2[13][0] = 180.0;
  this->octagon2[13][1] = 240.0;
  this->octagon2[14][0] = 0.0;
  this->octagon2[14][1] = 240.0;
  this->octagon2[15][0] = 0.0;
  this->octagon2[15][1] = 300.0;
  this->octagon2[16][0] = 180.0;
  this->octagon2[16][1] = 300.0;
  this->octagon2[17][0] = 0.0;
  this->octagon2[17][1] = 300.0;

//  this->octagon2 = {
//    { 0.000000,   0.000000},
//    { 180.000000, 0.000000},
//    { 0.000000,   0.00000},
//    { 0.000000,   60.0000},
//    { 180.000000, 60.0000},
//    { 0.000000,   60.0000},
//    { 0.000000,   120.000000},
//    { 180.000000, 120.00000},
//    { 0.000000,   120.00000},
//    { 0.000000,   180.0000},
//    { 180.000000, 180.0000},
//    { 0.000000,   180.0000},
//    { 0.000000,   240.000000},
//    { 180.000000, 240.00000},
//    { 0.000000,   240.00000},
//    { 0.000000,   300.0000},
//    { 180.000000, 300.0000},
//    { 0.000000,   300.0000}
//  };

}// end ::DetectionSystemSpice

DetectionSystemSpice::~DetectionSystemSpice()
{
  // LogicalVolumes in ConstructDetectionSystemSpice
  delete casing_log;
  delete crystal_block_log;

//  delete crystal_block_SD;

}// end ::~DetectionSystemSpice



G4int DetectionSystemSpice::Build()//G4SDManager* mySDman)
{ 
//  if( !crystal_block_SD ) {
//    crystal_block_SD = new SensitiveDetector("/sd/allSpice", "CollectionSpice");
//    mySDman->AddNewDetector( crystal_block_SD );
//  }

  // Build assembly volume
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;
  G4AssemblyVolume* myAssemblySi = new G4AssemblyVolume();
  this->assemblySi = myAssemblySi;

  BuildSiliconWafer();

//  G4RotationMatrix* finalRotate = new G4RotationMatrix();
//  finalRotate->rotateX(90.*deg);finalRotate->rotateY(0.*deg);finalRotate->rotateZ(0);

//   G4ThreeVector finalPlace;
//   finalPlace.setX(0); finalPlace.setY(0); finalPlace.setZ(this->detector2target);
  
  // Sensitive Detector
//  crystal_block_log->SetSensitiveDetector( crystal_block_SD );  

  return 1;
}

G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number)
{
  // Create Ring of Detectors
  for(G4int detector_number=0; detector_number < number_of_detectors; detector_number++) {
    G4ThreeVector move =  getDirection(detector_number); //direction vector from origin

    //Crystal Rotation
    G4RotationMatrix* rotate = new G4RotationMatrix; //rotation matrix corresponding to direction vector
    rotate->rotateX(-theta);
    if(detector_number == 1 || detector_number ==3) {
      rotate->rotateZ(-45.*deg);
      rotate->rotateY(phi);
      rotate->rotateX(90.*deg);
    }
    if(detector_number == 2) {
      rotate->rotateZ(90*deg);
      rotate->rotateY(phi);
      rotate->rotateX(-45.*deg);
    }
    if(detector_number == 0) {
      rotate->rotateZ(90*deg);
      rotate->rotateY(phi);
      rotate->rotateX(45.*deg);
    }
    assemblySi->MakeImprint(exp_hall_log, move, rotate, detector_number);
  }  

  return 1;
}

///////////////////////////////////////////////////////////////////////
//methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4int DetectionSystemSpice::BuildSiliconWafer()
{
  //vis attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  // silicon detector measurement  
  G4double half_length_x = (this->crystal_length_x)/2.0;
  G4double half_length_y2 = (this->crystal_length_y)/2.0;
  G4double half_length_z = (this->crystal_length_z)/2.0;
  G4double half_length_y1 = (this->crystal_length_y)/2.0;	
 //primtive volume
  G4VSolid* crystal_block = new G4Trd("crystal_block", half_length_x, half_length_x, half_length_y1, half_length_y2, half_length_z );

 //Establish Logical Volumes

  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }

  // Define rotation and movement objects
  G4ThreeVector move_null 	= G4ThreeVector(0,0,0);
  G4RotationMatrix* rotate_null = new G4RotationMatrix;
  
  //logical volume
  if( crystal_block_log == NULL )
  {
    crystal_block_log = new G4LogicalVolume(crystal_block, material, "crystal_block_log", 0, 0, 0);
    crystal_block_log->SetVisAttributes(vis_att);
  }

  this->assemblySi->AddPlacedVolume(crystal_block_log, move_null, rotate_null);

}//end ::BuildSiliconWafer


//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemSpice::getDirection(G4int detector_number)
{
  theta = getTheta(detector_number);
  phi = getPhi(detector_number);
  G4double transx, transy, transz;
  G4double distance = (this->crystal_dist_from_origin + (this->crystal_length_z/2));

//  if (detector_number%3 == 0 )  //Creates Side wafer
//  {
//    transx = -this->crystal_length_y/2. * cos(theta);
//    transz =  this->crystal_length_y/2. * sin(theta);
//  }
//  else if (detector_number%3 == 1)  //Create Center Wafer
//  {
    transx = 0.;
    transz = 0.;
//  }
//  else if (detector_number%3 == 2)
//  {
//    transx =  this->crystal_length_y/2. * cos(theta);
//    transz = -this->crystal_length_y/2. * sin(theta);
//  }

  G4double x,y,z;
  x = sin(theta) * distance + transx; // in xz plane
  y = -casing_threshold;
  z = cos(theta) * distance + transz;

  x = sin(theta) * distance + transx; // in xy plane
  z = -casing_threshold + this->detector2target;
  y = cos(theta) * distance + transz;


  G4ThreeVector direction(x,y,z);
	
  return direction;
}//end ::end getDirection


G4double DetectionSystemSpice::getTheta(G4int detector_number)
{
   theta = detector_number*90.*deg;
   
   return theta;
}

G4double DetectionSystemSpice::getPhi(G4int detector_number)
{
   phi = octagon2[detector_number][0]*deg;
   
   return phi;
}

 G4double DetectionSystemSpice::GetDetector2Origin(DetectionSystemSpice* spiceObj)
{  
  return 25.8212;
}

G4double DetectionSystemSpice::GetCrystalWidth()
{
  return 23.358*mm;
}

//void DetectionSystemSpice::BuildSiliconTube()
//{
//  //vis attributes
//  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
//  vis_att->SetVisibility(true);

//  // silicon detector measurement  
//  G4double half_length = (this->crystal_tube_length)/2.0;
//  G4double inner_rad = (this->crystal_inner_radius);
//  G4double outer_rad = (this->crystal_outer_radius);
//  G4Tubs* crystal_block = new G4Tubs("crystal_block", inner_rad, outer_rad, half_length, 0.*deg, 360*deg);

// //Establish Logical Volumes
//  G4Material* wafer_material = G4Material::GetMaterial(this->wafer_material);
//    crystal_block_log = new G4LogicalVolume(crystal_block, wafer_material, "crystal_block_log", 0, 0, 0);
//    crystal_block_log->SetVisAttributes(vis_att);
//  

//}//end::BuildSiliconTube


//void DetectionSystemSpice::BuildAluminumCasing()
//{

//   // Set visualization attributes
//  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
//  vis_att->SetVisibility(true);

//  // silicon detector measurements
//  G4double half_length_x = this->crystal_length_x/2.0;
//  G4double half_length_y = this->crystal_length_y/2.0;
//  G4double half_length_z = this->crystal_length_z/2.0;	

//  //primtive volume
//  G4VSolid* casing = new G4Trd("casing",half_length_x + casing_width, half_length_x + casing_width, half_length_y + casing_width*sin(60.*deg), 2*half_length_y + casing_width/cos(60.*deg), half_length_z + casing_width);
//  G4VSolid* detector_block = new G4Trd("detector_block", half_length_x, half_length_x, half_length_y, 2*half_length_y, half_length_z );
//  G4VSolid* top_opening = new G4Trd("top_opening", 0.6*mm, 0.6*mm, half_length_y - casing_width*sin(60.*deg), 2*half_length_y - casing_width/cos(60.*deg), half_length_z - 1.0*mm);
// 
//  G4VSolid* vacuum_gap = new G4UnionSolid("vacuum_gap", detector_block, top_opening, 0, G4ThreeVector(-half_length_x - 0.6*mm + 0.01*mm,0.,0.) );
//  
//  G4VSolid* casingshell = new G4SubtractionSolid("Casing_Shell", casing, vacuum_gap,0, G4ThreeVector(0.,0.,0.));  

//  // Establish logical volumes
//  G4Material* casing_material = G4Material::GetMaterial(this->casing_material);
//  casing_log = new G4LogicalVolume(casingshell, casing_material, "casing_log", 0, 0, 0);
//  casing_log->SetVisAttributes(vis_att);

//}//end ::BuildAluminumCasing

