#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemBrillance380V1.hh"

DetectionSystemBrillance380V1::DetectionSystemBrillance380V1() :
    // LogicalVolumes 
    detector_volume_log(0),    
    crystal_block_log(0),
    can_cylinder_log(0),
    can_front_lid_log(0),
    can_back_lid_log(0),
    can_vacuum_cylinder_log(0),
    can_vacuum_front_lid_log(0),
    can_vacuum_back_lid_log(0)
{ 
  
  /////////////////////////////////////////////////////////////////////
  // Brillance380 Crystal and Can Physical Properties 
  /////////////////////////////////////////////////////////////////////
  this->detail_view_end_angle	 	= 360.0*deg;     
  this->number_of_segments          = 20;  // huh?
  this->crystal_material            = "Cerium_Doped_Lanthanum_Bromide";
  this->can_material                = "G4_Al";
  this->vacuum_material             = "G4_Galactic";
  this->crystal_length_x            = 0.00*mm;
  this->crystal_length_y            = 0.00*mm;
  this->crystal_length_z            = 5.08*cm;
  this->crystal_inner_radius 		= 0.00*cm;
  this->crystal_outer_radius 		= 2.54*cm;
  this->can_thickness               = 0.05*cm;
  this->can_inner_radius            = 2.64*cm;
  this->can_lid_inner_radius 		= 0*cm;
  this->can_lid_outer_radius 		= 2.64*cm;	// sits flush inside can cylinder
  this->can_front_lid_thickness		= 0.05*cm;
  this->can_back_lid_thickness 		= 0.1*cm;
//  this->can_face_dist_from_origin	= this->inner_radius*cm;
  this->crystal_dist_from_can_face 	= 0.05*cm;  	// filled with vacuum
  this->crystal_dist_from_can_back 	= 0.1*cm;   	// filled with vacuum

  this->can_length_z			= this->crystal_length_z 
					+ this->can_front_lid_thickness
					+ this->can_back_lid_thickness 
					+ this->crystal_dist_from_can_face
					+ this->crystal_dist_from_can_back;	    

  //G4double triangleThetaAngle = (180/M_PI)*(atan((1/sqrt(3))/sqrt((11/12) + (1/sqrt(2))) )+atan((sqrt(2))/(1+sqrt(2))))*deg;
  G4double triangleThetaAngle = 54.735610317245360*deg;

  // theta
  this->detectorAngles[0][0] 	= triangleThetaAngle;
  this->detectorAngles[1][0] 	= triangleThetaAngle;
  this->detectorAngles[2][0] 	= triangleThetaAngle;
  this->detectorAngles[3][0] 	= triangleThetaAngle;
  this->detectorAngles[4][0] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[5][0] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[6][0] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[7][0] 	= 180.0*deg - triangleThetaAngle;
  // phi
  this->detectorAngles[0][1] 	= 22.5*deg;
  this->detectorAngles[1][1] 	= 112.5*deg;
  this->detectorAngles[2][1] 	= 202.5*deg;
  this->detectorAngles[3][1] 	= 292.5*deg;
  this->detectorAngles[4][1] 	= 22.5*deg;
  this->detectorAngles[5][1] 	= 112.5*deg;
  this->detectorAngles[6][1] 	= 202.5*deg;
  this->detectorAngles[7][1] 	= 292.5*deg;
  // yaw (alpha)
  this->detectorAngles[0][2] 	= 0.0*deg;
  this->detectorAngles[1][2] 	= 0.0*deg;
  this->detectorAngles[2][2] 	= 0.0*deg;
  this->detectorAngles[3][2] 	= 0.0*deg;
  this->detectorAngles[4][2] 	= 0.0*deg;
  this->detectorAngles[5][2] 	= 0.0*deg;
  this->detectorAngles[6][2] 	= 0.0*deg;
  this->detectorAngles[7][2] 	= 0.0*deg;
  // pitch (beta)
  this->detectorAngles[0][3] 	= triangleThetaAngle;
  this->detectorAngles[1][3] 	= triangleThetaAngle;
  this->detectorAngles[2][3] 	= triangleThetaAngle;
  this->detectorAngles[3][3] 	= triangleThetaAngle;
  this->detectorAngles[4][3] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[5][3] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[6][3] 	= 180.0*deg - triangleThetaAngle;
  this->detectorAngles[7][3] 	= 180.0*deg - triangleThetaAngle;
  // roll (gamma)
  this->detectorAngles[0][4] 	= 22.5*deg;
  this->detectorAngles[1][4] 	= 112.5*deg;
  this->detectorAngles[2][4] 	= 202.5*deg;
  this->detectorAngles[3][4] 	= 292.5*deg;
  this->detectorAngles[4][4] 	= 22.5*deg;
  this->detectorAngles[5][4] 	= 112.5*deg;
  this->detectorAngles[6][4] 	= 202.5*deg;
  this->detectorAngles[7][4] 	= 292.5*deg;
}

DetectionSystemBrillance380V1::~DetectionSystemBrillance380V1()
{
    // LogicalVolumes 
    delete detector_volume_log;    
    delete crystal_block_log;
    delete can_cylinder_log;
    delete can_front_lid_log;
    delete can_back_lid_log;
    delete can_vacuum_cylinder_log;
    delete can_vacuum_front_lid_log;
    delete can_vacuum_back_lid_log;

//    delete crystal_block_SD;
}

//G4int DetectionSystemBrillance380V1::Build(G4SDManager* mySDman)
G4int DetectionSystemBrillance380V1::Build()
{ 
//  if( !crystal_block_SD ) {
//    crystal_block_SD = new SensitiveDetector("/sd/allBrillance380V1", "CollectionBrillance380V1");
//    mySDman->AddNewDetector( crystal_block_SD );
//  }

  // Build assembly volume
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;

  G4cout << "BuildCrystalVolume" << G4endl;
  BuildCrystalVolume();      
  G4cout << "BuildAluminumCanVolume" << G4endl;
  BuildAluminumCanVolume(); 
  G4cout << "BuildCanVacuumVolume" << G4endl;
  BuildCanVacuumVolume();   

  // Sensitive Detector
//  crystal_block_log->SetSensitiveDetector( crystal_block_SD );  

  return 1;
}

G4int DetectionSystemBrillance380V1::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector_number)
{
  G4int detector_copy_ID = 0;

  G4cout << "Brillance380V1 Detector Number = " << detector_number << G4endl;

  G4int copy_number = detector_copy_ID + detector_number;

  G4double position = 12.5*cm + ( this->can_length_z / 2.0 ) ;

  G4double theta  = this->detectorAngles[detector_number][0];
  G4double phi    = this->detectorAngles[detector_number][1];
  G4double alpha  = this->detectorAngles[detector_number][2]; // yaw
  G4double beta   = this->detectorAngles[detector_number][3]; // pitch
  G4double gamma  = this->detectorAngles[detector_number][4]; // roll

  G4double x = 0;
  G4double y = 0;
  G4double z = position;

  G4RotationMatrix* rotate = new G4RotationMatrix;    // rotation matrix corresponding to direction vector
  rotate->rotateY(M_PI);
  rotate->rotateY(beta);
  rotate->rotateZ(gamma);

  G4ThreeVector move(transX(x,y,z,theta,phi), transY(x,y,z,theta,phi), transZ(x,y,z,theta,phi));

  assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);

  return 1;
}

G4int DetectionSystemBrillance380V1::BuildCrystalVolume()
{
  G4Material* material = G4Material::GetMaterial(this->crystal_material);
  if( !material ) {
    G4cout << " ----> Material " << this->crystal_material << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.2,1.0,0.3));
  vis_att->SetVisibility(true);

  G4Tubs* crystal_block = BuildCrystal();

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= ( (can_back_lid_thickness + crystal_dist_from_can_back) - (can_front_lid_thickness + crystal_dist_from_can_face) )/2.0;
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  
  //logical volume
  if( crystal_block_log == NULL )
  {
    crystal_block_log = new G4LogicalVolume(crystal_block, material, "brillance380_crystal_block", 0, 0, 0);
    crystal_block_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(crystal_block_log, move, rotate);

  return 1;
}

G4int DetectionSystemBrillance380V1::BuildAluminumCanVolume()
{
  G4Material* material = G4Material::GetMaterial(this->can_material);
  if( !material ) {
    G4cout << " ----> Material " << this->can_material << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  vis_att->SetVisibility(true);  

  G4ThreeVector direction =  G4ThreeVector(0,0,1);
  G4double z_position;
  G4ThreeVector move; 
  G4RotationMatrix* rotate = new G4RotationMatrix;

  /////////////////////////////////////////////////////////////////////
  // Build and Place Aluminum Can
  /////////////////////////////////////////////////////////////////////
  G4Tubs* can_cylinder = BuildAluminumCan(); 

  //logical volume
  if( can_cylinder_log == NULL )
  {
    can_cylinder_log = new G4LogicalVolume(can_cylinder, material, "can_cylinder_log", 0, 0, 0);
    can_cylinder_log->SetVisAttributes(vis_att);
  }

  // place front can_lid
  z_position 	= 0;
  move 		= z_position * direction;
  
  //add physical cylinder
  this->assembly->AddPlacedVolume(can_cylinder_log, move, rotate);

  /////////////////////////////////////////////////////////////////////
  // Build and Place Aluminum Front Lid
  /////////////////////////////////////////////////////////////////////  
  G4Tubs* can_front_lid = BuildAluminumCanFrontLid();
  
  // logical volume
  if( can_front_lid_log == NULL )
  {
    can_front_lid_log = new G4LogicalVolume(can_front_lid, material, "can_front_lid_log", 0, 0, 0);
    can_front_lid_log->SetVisAttributes(vis_att);
  }

  // place front can_lid
  z_position 	= (can_length_z/2.0) - (can_front_lid_thickness/2.0);
  move 		= z_position * direction;
  
  //add physical front can_lid
  this->assembly->AddPlacedVolume(can_front_lid_log, move, rotate);
    
  /////////////////////////////////////////////////////////////////////
  // Build and Place Aluminum Back Lid
  /////////////////////////////////////////////////////////////////////
  G4Tubs* can_back_lid = BuildAluminumCanBackLid();
    
  // logical volume
  if( can_back_lid_log == NULL )
  {
    can_back_lid_log = new G4LogicalVolume(can_back_lid, material, "can_back_lid_log", 0, 0, 0);
    can_back_lid_log->SetVisAttributes(vis_att);
  }

  // place back can_lid
  z_position 	= -(can_length_z/2.0) + (can_back_lid_thickness/2.0);
  move 		= z_position * direction;

  // add physical back can_lid
  this->assembly->AddPlacedVolume(can_back_lid_log, move, rotate);
  
  return 1;
}


G4int DetectionSystemBrillance380V1::BuildCanVacuumVolume()
{
  G4Material* material = G4Material::GetMaterial(this->vacuum_material);
  if( !material ) {
    G4cout << " ----> Material " << this->vacuum_material << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
 
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.1,0.0,0.9));
  vis_att->SetVisibility(true);

  G4ThreeVector direction =  G4ThreeVector(0,0,1);
  G4double z_position;
  G4ThreeVector move; 
  G4RotationMatrix* rotate = new G4RotationMatrix;

  /////////////////////////////////////////////////////////////////////
  // Build and Place Vacuum Can
  /////////////////////////////////////////////////////////////////////
  G4Tubs* can_vacuum_cylinder = BuildCanVacuum(); 

  // Define rotation and movement objects
  z_position	= ( (can_back_lid_thickness + crystal_dist_from_can_back) - (can_front_lid_thickness + crystal_dist_from_can_face) )/2.0;
  move 		= z_position * direction;

  // logical volume
  if( can_vacuum_cylinder_log == NULL )
  {
    can_vacuum_cylinder_log = new G4LogicalVolume(can_vacuum_cylinder, material, "can_vacuum_cylinder_log", 0, 0, 0);
    can_vacuum_cylinder_log->SetVisAttributes(vis_att);
  }
	  
  // add physical can cylinder
  this->assembly->AddPlacedVolume(can_vacuum_cylinder_log, move, rotate); 
        
  /////////////////////////////////////////////////////////////////////
  // Build and Place Vacuum Front Lid
  /////////////////////////////////////////////////////////////////////  
  G4Tubs* can_vacuum_front_lid = BuildCanVacuumFrontLid();

  // logical volume
  if( can_vacuum_front_lid_log == NULL )
  {
    can_vacuum_front_lid_log = new G4LogicalVolume(can_vacuum_front_lid, material, "can_vacuum_front_lid_log", 0, 0, 0);
    can_vacuum_front_lid_log->SetVisAttributes(vis_att);
  }

  // place can vacuum front lid
  z_position 	= (can_length_z/2.0) - (can_front_lid_thickness) - (crystal_dist_from_can_face/2.0);
  move 		= z_position * direction;
  
  // add physical BACK outer_can_lid
  this->assembly->AddPlacedVolume(can_vacuum_front_lid_log, move, rotate); 
  
  /////////////////////////////////////////////////////////////////////
  // Build and Place Vacuum Back Lid
  /////////////////////////////////////////////////////////////////////
  G4Tubs* can_vacuum_back_lid = BuildCanVacuumBackLid();

  // logical volume
  if( can_vacuum_back_lid_log == NULL )
  {
    can_vacuum_back_lid_log = new G4LogicalVolume(can_vacuum_back_lid, material, "can_vacuum_back_lid_log", 0, 0, 0);
    can_vacuum_back_lid_log->SetVisAttributes(vis_att);
  }

  // place can vacuum front lid
  z_position 	= - (can_length_z/2.0) + (can_back_lid_thickness) + (crystal_dist_from_can_back/2.0);
  move 		= z_position * direction;
  
  // add physical BACK outer_can_lid
  this->assembly->AddPlacedVolume(can_vacuum_back_lid_log, move, rotate); 
  
  return 1;
}//end ::BuildCanVacuumVolume

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemBrillance380V1::BuildCrystal()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = crystal_inner_radius;
  G4double outer_radius = crystal_outer_radius;
  G4double half_length_z = (crystal_length_z)/2.0;

  G4Tubs* crystal_block = new G4Tubs("crystal_block", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return crystal_block;
}//end ::BuildCrystal

G4Tubs* DetectionSystemBrillance380V1::BuildAluminumCan()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius 	= can_inner_radius;  			
  G4double outer_radius 	= can_inner_radius + can_thickness;
  G4double half_length_z 	= can_length_z/2.0;

  G4Tubs* can_cylinder = new G4Tubs("can_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_cylinder;
}//end ::BuildAluminumCan

G4Tubs* DetectionSystemBrillance380V1::BuildAluminumCanFrontLid()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = can_lid_inner_radius;
  G4double outer_radius = can_lid_outer_radius;
  G4double half_length_z = can_front_lid_thickness/2.0;
  
  G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_lid;
}//end ::BuildAluminumFrontLid

G4Tubs* DetectionSystemBrillance380V1::BuildAluminumCanBackLid()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = can_lid_inner_radius;
  G4double outer_radius = can_lid_outer_radius;
  G4double half_length_z = can_back_lid_thickness/2.0;
  
  G4Tubs* can_lid = new G4Tubs("can_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_lid;
}//end ::BuildAluminumBackLid

G4Tubs* DetectionSystemBrillance380V1::BuildCanVacuum()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = crystal_outer_radius;	
  G4double outer_radius = can_inner_radius;
  G4double half_length_z = crystal_length_z/2.0;
  
  G4Tubs* can_vacuum_cylinder = new G4Tubs("can_vacuum_cylinder", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_vacuum_cylinder;
}//end ::BuildCanVacuum

G4Tubs* DetectionSystemBrillance380V1::BuildCanVacuumFrontLid()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = can_lid_inner_radius;
  G4double outer_radius = can_lid_outer_radius;
  G4double half_length_z = crystal_dist_from_can_face/2.0; 
  
  G4Tubs* can_vacuum_cylinder_lid = new G4Tubs("can_vacuum_cylinder_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_vacuum_cylinder_lid;
}//end ::BuildCanVacuumFrontLid

G4Tubs* DetectionSystemBrillance380V1::BuildCanVacuumBackLid()
{
  G4double start_phi = 0.0;
  G4double end_phi = this->detail_view_end_angle;

  G4double inner_radius = can_lid_inner_radius;
  G4double outer_radius = can_lid_outer_radius;
  G4double half_length_z = crystal_dist_from_can_back/2.0; 
  
  G4Tubs* can_vacuum_cylinder_lid = new G4Tubs("can_vacuum_cylinder_lid", inner_radius, outer_radius, half_length_z, start_phi, end_phi);

  return can_vacuum_cylinder_lid;
}//end ::BuildCanVacuumBackLid

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemBrillance380V1::GetDirectionXYZ(G4double theta, G4double phi)
{
  G4double x,y,z;
  x = sin(theta) * cos(phi);
  y = sin(theta) * sin(phi);
  z = cos(theta);
	
  G4ThreeVector direction = G4ThreeVector(x,y,z);
	
  return direction;
}//end ::GetDirection

G4double DetectionSystemBrillance380V1::transX(G4double x, G4double y, G4double z, G4double theta, G4double phi){
  return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*cos(phi) );
}

G4double DetectionSystemBrillance380V1::transY(G4double x, G4double y, G4double z, G4double theta, G4double phi){
  return ( pow(x*x+y*y+z*z,0.5)*sin(theta)*sin(phi) );
}

G4double DetectionSystemBrillance380V1::transZ(G4double x, G4double y, G4double z, G4double theta, G4double phi){
  return ( pow(x*x+y*y+z*z,0.5)*cos(theta) );
}
