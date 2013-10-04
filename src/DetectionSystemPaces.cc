#include "DetectorMessenger.hh"
//#include "SensitiveDetector.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

//#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemPaces.hh"

DetectionSystemPaces::DetectionSystemPaces() :
  // Logical Volumes
  aluminum_hemisphere_log(0),
  aluminum_annulus_top_log(0),
  aluminum_annulus_bot_log(0),
  canister_log(0),
  silicon_block_log(0),
  silicon_dead_layer_log(0),
  teflon_annulus_top_log(0),
  teflon_annulus_bot_log(0),
  delrin_hemisphere_log(0)
{
  /*Measurement status*/
  //~   estimate
  //OK  confirmed measurement
  
  // Cut clearance
  this->cut_clearance = 0.001*mm;
  
  // Aluminum hemisphere
  this->aluminum_hemisphere_inner_radius =         (82.50 / 2)*mm; //OK
  this->aluminum_hemisphere_outer_radius =         (88.86 / 2)*mm; //OK
  this->aluminum_hemisphere_beam_hole_radius =     (30.00 / 2)*mm; //~
  this->aluminum_hemisphere_beam_hole_rim_height = (2.0)*mm; //~
  
  // Annulus that sits in front of silicon crystal
  this->aluminum_annulus_top_inner_radius =        (15.86 / 2)*mm; //OK
  this->aluminum_annulus_top_outer_radius =        (31.71 / 2)*mm; //OK
  this->aluminum_annulus_top_thickness =           (2.37)*mm; //OK
  
  // Annulus that sits in back of silicon crystal
  this->aluminum_annulus_bot_inner_radius =        (20.00 / 2)*mm; //~
  this->aluminum_annulus_bot_outer_radius =        (25.00 / 2)*mm; //~
  this->aluminum_annulus_bot_thickness =           (5.00)*mm; //~
  
  // Canister holding crystal
  this->canister_inner_radius =                    (20.00 / 2)*mm; //~
  this->canister_outer_radius =                    (21.00 / 2)*mm; //~
  this->canister_thickness =                       (5.00)*mm; //~
  
  // Silicon crystal detector (block parameters includes dead layer thickness)
  this->silicon_block_radius =                     (20.00 / 2)*mm; //~
  this->silicon_block_thickness =                  (5.00)*mm; //~
  this->silicon_dead_layer_thickness =             (0.10)*mm; //~
  
  // Screws
  this->screw_radius =                             (1.00 / 2)*mm; //~
  this->screw_placement_radius =                   (22.50 / 2)*mm; //~
  
  // Teflon bases
  this->teflon_annulus_top_inner_radius =          (21.00 / 2)*mm; //~
  this->teflon_annulus_top_outer_radius =          (25.00 / 2)*mm; //~
  this->teflon_annulus_top_thickness =             (2.22)*mm; //OK
  this->teflon_annulus_bot_inner_radius =          (25.00 / 2)*mm; //~
  this->teflon_annulus_bot_outer_radius =          (32.00 / 2)*mm; //~
  this->teflon_annulus_bot_thickness =             (5.00)*mm; //~
  
  // delrin hemisphere
  this->delrin_hemisphere_inner_radius =           (84.4)*mm; //~
  this->delrin_hemisphere_outer_radius =           (89.4)*mm; //~
  this->delrin_hemisphere_beam_hole_radius =       (20.00 / 2)*mm; //~
  
  // placement distances
  this->aluminum_hemisphere_dist =     (0.)*mm;
  this->delrin_hemisphere_dist =       (0.)*mm;
  this->aluminum_annulus_top_dist =    (aluminum_annulus_top_thickness / 2);
  this->canister_dist =                (aluminum_annulus_top_thickness) + (canister_thickness / 2) ;
  this->silicon_block_dist =           (aluminum_annulus_top_thickness) + (silicon_block_thickness / 2) ;
  this->teflon_annulus_top_dist =      (aluminum_annulus_top_thickness) + (silicon_block_thickness) - (teflon_annulus_top_thickness / 2);
  this->teflon_annulus_bot_dist =      (teflon_annulus_top_dist) + (teflon_annulus_top_thickness / 2) + (teflon_annulus_bot_thickness / 2);
  this->aluminum_annulus_bot_dist =    (teflon_annulus_top_dist) + (teflon_annulus_top_thickness / 2) + (aluminum_annulus_bot_thickness / 2);
  
  // Placement and Orientation scalars
  /*paces_placement_distance[5] = {  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm ,  30.00*mm };
  paces_placement_theta[5] =    {   0.0*deg ,  72.0*deg , 144.0*deg , 216.0*deg , 288.0*deg };
  paces_placement_phi[5] =      {  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg ,  60.0*deg };
  paces_orientation_theta[5] =  {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };
  paces_orientation_phi[5] =    {   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg ,   0.0*deg };*/
  paces_placement_distance[0] = 30.00*mm;
  paces_placement_distance[1] = 30.00*mm;
  paces_placement_distance[2] = 30.00*mm;
  paces_placement_distance[3] = 30.00*mm;
  paces_placement_distance[4] = 30.00*mm;
  paces_placement_theta[0] =   0.0000*deg;
  paces_placement_theta[1] =  68.7575*deg;
  paces_placement_theta[2] = 138.4082*deg;
  paces_placement_theta[3] = 212.1280*deg;
  paces_placement_theta[4] = 287.4275*deg;
  paces_placement_phi[0] = 66.0*deg;
  paces_placement_phi[1] = 66.0*deg;
  paces_placement_phi[2] = 66.0*deg;
  paces_placement_phi[3] = 66.0*deg;
  paces_placement_phi[4] = 66.0*deg;
  paces_orientation_theta[0] = 0.0*deg;
  paces_orientation_theta[1] = 0.0*deg;
  paces_orientation_theta[2] = 0.0*deg;
  paces_orientation_theta[3] = 0.0*deg;
  paces_orientation_theta[4] = 0.0*deg;
  paces_orientation_phi[0] = 0.1780*deg;
  paces_orientation_phi[1] = 0.8275*deg;
  paces_orientation_phi[2] =-0.2579*deg;
  paces_orientation_phi[3] = 0.2992*deg;
  paces_orientation_phi[4] = 0.1934*deg;
  
  // assmebly volume sizes
  this->detector_assembly_radius =    (32.00 / 2)*mm;
  this->detector_assembly_thickness = (15.00)*mm;
}

DetectionSystemPaces::~DetectionSystemPaces()
{
  // logical volumes
  delete exp_hall_log;
  delete aluminum_hemisphere_log;
  delete aluminum_annulus_top_log;
  delete aluminum_annulus_bot_log;
  delete canister_log;
  delete silicon_block_log;
  delete silicon_dead_layer_log;
  delete teflon_annulus_top_log;
  delete teflon_annulus_bot_log;
  delete delrin_hemisphere_log;
  // sensitive detectors
//  delete silicon_block_SD;
}

G4int DetectionSystemPaces::Build()// G4SDManager* mySDman)
{
  // Add sensitive detectors
//  if( !silicon_block_SD ) {
//    silicon_block_SD = new SensitiveDetector("/sd/allPaces", "CollectionPacesSilicon");
//    mySDman->AddNewDetector( silicon_block_SD );
//  }
  
  // Build assembly volumes
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;
  G4AssemblyVolume* myAssemblyDetector = new G4AssemblyVolume();
  this->assemblyDetector = myAssemblyDetector;
  G4AssemblyVolume* myAssemblySilicon = new G4AssemblyVolume();
  this->assemblySilicon = myAssemblySilicon;
  
  // Add silicon assembly
  AddSiliconBlock();
  AddSiliconDeadLayer();
  
  // Add detector assembly
  AddCanister();
  AddAluminumAnnulusTop();
  AddTeflonAnnulusTop();
  AddTeflonAnnulusBot();
  AddAluminumAnnulusBot();
  AddScrews();
  
  // Add total assembly
  AddAluminumHemisphere();
  
  // Assemble the ensemble
  CombineAssemblySilicon();
  CombineAssemblyDetector();
  CombineAssembly();
  
//  silicon_block_log->SetSensitiveDetector( silicon_block_SD );
  
  return 1;
}//end ::Build

G4int DetectionSystemPaces::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int ndet)
{
  //place detectors
  G4int d_i;
  G4double* ptr_pd = this->paces_placement_distance;
  G4double* ptr_pt = this->paces_placement_theta;
  G4double* ptr_pp = this->paces_placement_phi;
  G4double* ptr_ot = this->paces_orientation_theta;
  G4double* ptr_op = this->paces_orientation_phi;
  if (ndet > 5 || ndet < 0) ndet = 5;
  for (G4int i=0; i<ndet; i++)
  {
    d_i = i;
    G4double d_dist = ptr_pd[d_i], d_theta = ptr_pt[d_i], d_phi = ptr_pp[d_i];
    G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
    Ta.setX( d_dist * cos(d_theta) * sin(d_phi) );
    Ta.setY( d_dist * sin(d_theta) * sin(d_phi) );
    Ta.setZ( d_dist *      1.0     * cos(d_phi) );
    G4double ori_theta = d_theta + ptr_ot[d_i] + pi/2; //plus 90 deg
    G4double ori_phi = d_phi + ptr_op[d_i];
    G4ThreeVector yprimeaxis = G4ThreeVector(cos(ori_theta), sin(ori_theta), 0);
    Ra->set(yprimeaxis, ori_phi);
    //G4cout << "----------- d_i = " << d_i << G4endl;
    this->assemblySilicon->MakeImprint(exp_hall_log, Ta, Ra, d_i+1);
    this->assemblyDetector->MakeImprint(exp_hall_log, Ta, Ra, d_i*20);
  }

  //place hemisphere
  G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
  this->assembly->MakeImprint(exp_hall_log, T0, R0, 5*40);
  
  return 1;
}//end ::PlaceDetector

G4int DetectionSystemPaces::CombineAssemblySilicon()
{
  G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
  Ta.setZ( this->silicon_block_dist );
  //this->assemblySilicon->AddPlacedVolume(silicon_block_log, Ta, Ra);
  this->assemblySilicon->AddPlacedVolume(silicon_block_log, Ta, Ra);
  
  return 1;
}//end ::CombineAssemblySilicon


G4int DetectionSystemPaces::CombineAssemblyDetector()
{
  G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
  Ta.setZ( this->silicon_block_dist );
  //this->assemblyDetector->AddPlacedAssembly(assemblySilicon, Ta, Ra);
  //this->assemblyDetector->AddPlacedVolume(silicon_block_log, Ta, Ra);
  this->assemblyDetector->AddPlacedVolume(silicon_dead_layer_log, Ta, Ra);
  Ta.setZ( this->canister_dist );
  this->assemblyDetector->AddPlacedVolume(canister_log, Ta, Ra);
  Ta.setZ( this->aluminum_annulus_top_dist );
  this->assemblyDetector->AddPlacedVolume(aluminum_annulus_top_log, Ta, Ra);
  Ta.setZ( this->aluminum_annulus_bot_dist );
  this->assemblyDetector->AddPlacedVolume(aluminum_annulus_bot_log, Ta, Ra);
  Ta.setZ( this->teflon_annulus_top_dist );
  this->assemblyDetector->AddPlacedVolume(teflon_annulus_top_log, Ta, Ra);
  Ta.setZ( this->teflon_annulus_bot_dist );
  this->assemblyDetector->AddPlacedVolume(teflon_annulus_bot_log, Ta, Ra);
  
  return 1;
}//end ::CombineAssemblyDetector

G4int DetectionSystemPaces::CombineAssembly()
{
  //place aluminum hemisphere
  G4RotationMatrix* Ra = new G4RotationMatrix; G4ThreeVector Ta;
  Ta.setZ( this->aluminum_hemisphere_dist );
  this->assembly->AddPlacedVolume(aluminum_hemisphere_log, Ta, Ra);
  
  return 1;
}//end ::CombineAssembly

//Add silicon detector
G4int DetectionSystemPaces::AddSiliconBlock()
{
  //material
  G4Material* material = G4Material::GetMaterial("Silicon");
  if( !material ) {
    G4cout << " ----> Material " << "Silicon" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* silicon_block_vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  silicon_block_vis_att->SetVisibility(true);

  //silicon block (without dead layer, smaller than actual)
  G4double dead_layer_cut = this->silicon_dead_layer_thickness + this->cut_clearance;
  G4double inner_radius = 0.*mm;
  G4double outer_radius = this->silicon_block_radius - dead_layer_cut;
  G4double half_length_z = this->silicon_block_thickness/2.0 - dead_layer_cut;
  
  //primitive volume
  G4Tubs* silicon_block = new G4Tubs("silicon_block", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  
  //logical volume
  if( silicon_block_log == NULL ) {
    silicon_block_log = new G4LogicalVolume(silicon_block, material, "silicon_block_log", 0, 0, 0);
    silicon_block_log->SetVisAttributes(silicon_block_vis_att);
  }
  
  return 1;
}//end ::end AddSiliconBlock

//Add detector dead layer
G4int DetectionSystemPaces::AddSiliconDeadLayer()
{
  //material
  G4Material* material = G4Material::GetMaterial("Silicon");
  if( !material ) {
    G4cout << " ----> Material " << "Silicon" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* silicon_dead_layer_vis_att = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  silicon_dead_layer_vis_att->SetVisibility(true);
  
  //dead layer (outside of silicon block)
  G4double dl_thickness = this->silicon_dead_layer_thickness;
  G4double inner_radius = this->silicon_block_radius - dl_thickness;
  G4double outer_radius = this->silicon_block_radius;
  G4double half_length_z = this->silicon_block_thickness/2.0;
  //cut of original block
  G4double cut_inner_radius = 0.*mm;
  G4double cut_outer_radius = inner_radius;
  G4double cut_half_length_z = half_length_z - dl_thickness;
  
  //primitive volume
  G4Tubs* silicon_dead_layer_with_block = new G4Tubs("silicon_dead_layer_with_block", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_silicon_block = new G4Tubs("cut_silicon_block", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0*pi, 2*pi);
  
  //cut out original block
  G4ThreeVector move_cut = G4ThreeVector(0, 0, 0);
  G4SubtractionSolid* silicon_dead_layer = new G4SubtractionSolid("silicon_dead_layer", silicon_dead_layer_with_block, cut_silicon_block, 0, move_cut);
  
  //logical volume
  if( silicon_dead_layer_log == NULL ) {
    silicon_dead_layer_log = new G4LogicalVolume(silicon_dead_layer, material, "silicon_dead_layer_log", 0, 0, 0);
    silicon_dead_layer_log->SetVisAttributes(silicon_dead_layer_vis_att);
  }
  
  return 1;
}//end ::end AddSiliconDeadLayer

//Add detector canister
G4int DetectionSystemPaces::AddCanister()
{
  //material
  G4Material* material = G4Material::GetMaterial("Mylar");
  if( !material ) {
    G4cout << " ----> Material " << "Mylar" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* canister_vis_att = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  canister_vis_att->SetVisibility(true);
  
  //main annulus
  G4double inner_radius = this->canister_inner_radius + this->cut_clearance;
  G4double outer_radius = this->canister_outer_radius;
  G4double half_length_z = this->canister_thickness/2.0;
  
  //primitive volume
  G4Tubs* canister = new G4Tubs("canister", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  
  //logical volume
  if( canister_log == NULL ) {
    canister_log = new G4LogicalVolume(canister, material, "canister_log", 0, 0, 0);
    canister_log->SetVisAttributes(canister_vis_att);
  }
  
  return 1;
}//end ::end AddCanister

//Add annulus that is on top of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusTop()
{
  //material
  G4Material* material = G4Material::GetMaterial("Aluminum");
  if( !material ) {
    G4cout << " ----> Material " << "Aluminum" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* aluminum_annulus_top_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  aluminum_annulus_top_vis_att->SetVisibility(true);
  
  //main annulus
  G4double inner_radius = this->aluminum_annulus_top_inner_radius;
  G4double outer_radius = this->aluminum_annulus_top_outer_radius;
  G4double half_length_z = this->aluminum_annulus_top_thickness/2.0;
  //screw cuts
  G4double cut_inner_radius = 0.*mm;
  G4double cut_outer_radius = this->screw_radius;
  G4double cut_half_length_z = 20.*mm;
  
  //primitive volume
  G4Tubs* al_an_to_without_cuts = new G4Tubs("al_an_to_without_cuts", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_screw_hole = new G4Tubs("cut_screw_hole", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0*pi, 2*pi);
  
  //cut out screw holes
  G4double c_i;
  G4double c_r = this->screw_placement_radius;
  G4double c_dt = 45.*deg;
  G4ThreeVector move_cut;
  c_i = 0; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_0 = new G4SubtractionSolid("al_an_to_0", al_an_to_without_cuts, cut_screw_hole, 0, move_cut);
  c_i = 1; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_1 = new G4SubtractionSolid("al_an_to_1", al_an_to_0, cut_screw_hole, 0, move_cut);
  c_i = 2; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_2 = new G4SubtractionSolid("al_an_to_2", al_an_to_1, cut_screw_hole, 0, move_cut);
  c_i = 3; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_3 = new G4SubtractionSolid("al_an_to_3", al_an_to_2, cut_screw_hole, 0, move_cut);
  c_i = 4; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_4 = new G4SubtractionSolid("al_an_to_4", al_an_to_3, cut_screw_hole, 0, move_cut);
  c_i = 5; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_5 = new G4SubtractionSolid("al_an_to_5", al_an_to_4, cut_screw_hole, 0, move_cut);
  c_i = 6; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_to_6 = new G4SubtractionSolid("al_an_to_6", al_an_to_5, cut_screw_hole, 0, move_cut);
  //final cut
  c_i = 7; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* aluminum_annulus_top = new G4SubtractionSolid("aluminum_annulus_top", al_an_to_6, cut_screw_hole, 0, move_cut);
  
  //logical volume
  if( aluminum_annulus_top_log == NULL ) {
    aluminum_annulus_top_log = new G4LogicalVolume(aluminum_annulus_top, material, "aluminum_annulus_top_log", 0, 0, 0);
    aluminum_annulus_top_log->SetVisAttributes(aluminum_annulus_top_vis_att);
  }
  
  return 1;
}//end ::end AddAluminumAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddAluminumAnnulusBot()
{
  //material
  G4Material* material = G4Material::GetMaterial("Aluminum");
  if( !material ) {
    G4cout << " ----> Material " << "Aluminum" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* aluminum_annulus_bot_vis_att = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  aluminum_annulus_bot_vis_att->SetVisibility(true);
  
  //main annulus
  G4double inner_radius = this->aluminum_annulus_bot_inner_radius;
  G4double outer_radius = this->aluminum_annulus_bot_outer_radius;
  G4double half_length_z = this->aluminum_annulus_bot_thickness/2.0;
  //screw cuts
  G4double cut_inner_radius = 0.*mm;
  G4double cut_outer_radius = this->screw_radius;
  G4double cut_half_length_z = 20.*mm;
  
  //primitive volume
  G4Tubs* al_an_bo_without_cuts = new G4Tubs("al_an_bo_without_cuts", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_screw_hole = new G4Tubs("cut_screw_hole", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0*pi, 2*pi);
  
  //cut out screw holes
  G4double c_i;
  G4double c_r = this->screw_placement_radius;
  G4double c_dt = 45.*deg;
  G4ThreeVector move_cut;
  c_i = 0; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_0 = new G4SubtractionSolid("al_an_bo_0", al_an_bo_without_cuts, cut_screw_hole, 0, move_cut);
  c_i = 1; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_1 = new G4SubtractionSolid("al_an_bo_1", al_an_bo_0, cut_screw_hole, 0, move_cut);
  c_i = 2; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_2 = new G4SubtractionSolid("al_an_bo_2", al_an_bo_1, cut_screw_hole, 0, move_cut);
  c_i = 3; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_3 = new G4SubtractionSolid("al_an_bo_3", al_an_bo_2, cut_screw_hole, 0, move_cut);
  c_i = 4; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_4 = new G4SubtractionSolid("al_an_bo_4", al_an_bo_3, cut_screw_hole, 0, move_cut);
  c_i = 5; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_5 = new G4SubtractionSolid("al_an_bo_5", al_an_bo_4, cut_screw_hole, 0, move_cut);
  c_i = 6; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* al_an_bo_6 = new G4SubtractionSolid("al_an_bo_6", al_an_bo_5, cut_screw_hole, 0, move_cut);
  //final cut
  c_i = 7; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* aluminum_annulus_bot = new G4SubtractionSolid("aluminum_annulus_bot", al_an_bo_6, cut_screw_hole, 0, move_cut);
  
  //logical volume
  if( aluminum_annulus_bot_log == NULL ) {
    aluminum_annulus_bot_log = new G4LogicalVolume(aluminum_annulus_bot, material, "aluminum_annulus_bot_log", 0, 0, 0);
    aluminum_annulus_bot_log->SetVisAttributes(aluminum_annulus_bot_vis_att);
  }
  
  return 1;
}//end ::end AddAluminumAnnulusBot

//Add annulus that is on side of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusTop()
{
  //material
  G4Material* material = G4Material::GetMaterial("Teflon");
  if( !material ) {
    G4cout << " ----> Material " << "Teflon" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* teflon_annulus_top_vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  teflon_annulus_top_vis_att->SetVisibility(true);
  
  //main annulus
  G4double inner_radius = this->teflon_annulus_top_inner_radius;
  G4double outer_radius = this->teflon_annulus_top_outer_radius;
  G4double half_length_z = this->teflon_annulus_top_thickness/2.0;
  //screw cuts
  G4double cut_inner_radius = 0.*mm;
  G4double cut_outer_radius = this->screw_radius;
  G4double cut_half_length_z = 20.*mm;
  
  //primitive volume
  G4Tubs* te_an_to_without_cuts = new G4Tubs("te_an_to_without_cuts", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_screw_hole = new G4Tubs("cut_screw_hole", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0*pi, 2*pi);
  
  //cut out screw holes
  G4double c_i;
  G4double c_r = this->screw_placement_radius;
  G4double c_dt = 45.*deg;
  G4ThreeVector move_cut;
  c_i = 0; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_0 = new G4SubtractionSolid("te_an_to_0", te_an_to_without_cuts, cut_screw_hole, 0, move_cut);
  c_i = 1; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_1 = new G4SubtractionSolid("te_an_to_1", te_an_to_0, cut_screw_hole, 0, move_cut);
  c_i = 2; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_2 = new G4SubtractionSolid("te_an_to_2", te_an_to_1, cut_screw_hole, 0, move_cut);
  c_i = 3; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_3 = new G4SubtractionSolid("te_an_to_3", te_an_to_2, cut_screw_hole, 0, move_cut);
  c_i = 4; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_4 = new G4SubtractionSolid("te_an_to_4", te_an_to_3, cut_screw_hole, 0, move_cut);
  c_i = 5; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_5 = new G4SubtractionSolid("te_an_to_5", te_an_to_4, cut_screw_hole, 0, move_cut);
  c_i = 6; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_to_6 = new G4SubtractionSolid("te_an_to_6", te_an_to_5, cut_screw_hole, 0, move_cut);
  //final cut
  c_i = 7; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* teflon_annulus_top = new G4SubtractionSolid("teflon_annulus_top", te_an_to_6, cut_screw_hole, 0, move_cut);
  
  //logical volume
  if( teflon_annulus_top_log == NULL ) {
    teflon_annulus_top_log = new G4LogicalVolume(teflon_annulus_top, material, "teflon_annulus_top_log", 0, 0, 0);
    teflon_annulus_top_log->SetVisAttributes(teflon_annulus_top_vis_att);
  }
  
  return 1;
}//end ::end AddTeflonAnnulusTop

//Add annulus that is on bottom of silicon detector
G4int DetectionSystemPaces::AddTeflonAnnulusBot()
{
  //material
  G4Material* material = G4Material::GetMaterial("Teflon");
  if( !material ) {
    G4cout << " ----> Material " << "Teflon" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* teflon_annulus_bot_vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  teflon_annulus_bot_vis_att->SetVisibility(true);
  
  //main annulus
  G4double inner_radius = this->teflon_annulus_bot_inner_radius;
  G4double outer_radius = this->teflon_annulus_bot_outer_radius;
  G4double half_length_z = this->teflon_annulus_bot_thickness/2.0;
  //screw cuts
  G4double cut_inner_radius = 0.*mm;
  G4double cut_outer_radius = this->screw_radius;
  G4double cut_half_length_z = 20.*mm;
  
  //primitive volume
  G4Tubs* te_an_bo_without_cuts = new G4Tubs("te_an_bo_without_cuts", inner_radius, outer_radius, half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_screw_hole = new G4Tubs("cut_screw_hole", cut_inner_radius, cut_outer_radius, cut_half_length_z, 0*pi, 2*pi);
  
  //cut out screw holes
  G4double c_i;
  G4double c_r = this->screw_placement_radius;
  G4double c_dt = 45.*deg;
  G4ThreeVector move_cut;
  c_i = 0; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_0 = new G4SubtractionSolid("te_an_bo_0", te_an_bo_without_cuts, cut_screw_hole, 0, move_cut);
  c_i = 1; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_1 = new G4SubtractionSolid("te_an_bo_1", te_an_bo_0, cut_screw_hole, 0, move_cut);
  c_i = 2; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_2 = new G4SubtractionSolid("te_an_bo_2", te_an_bo_1, cut_screw_hole, 0, move_cut);
  c_i = 3; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_3 = new G4SubtractionSolid("te_an_bo_3", te_an_bo_2, cut_screw_hole, 0, move_cut);
  c_i = 4; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_4 = new G4SubtractionSolid("te_an_bo_4", te_an_bo_3, cut_screw_hole, 0, move_cut);
  c_i = 5; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_5 = new G4SubtractionSolid("te_an_bo_5", te_an_bo_4, cut_screw_hole, 0, move_cut);
  c_i = 6; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* te_an_bo_6 = new G4SubtractionSolid("te_an_bo_6", te_an_bo_5, cut_screw_hole, 0, move_cut);
  //final cut
  c_i = 7; move_cut.setX( c_r*cos(c_i*c_dt) ); move_cut.setY( c_r*sin(c_i*c_dt) );
  G4SubtractionSolid* teflon_annulus_bot = new G4SubtractionSolid("teflon_annulus_bot", te_an_bo_6, cut_screw_hole, 0, move_cut);
  
  //logical volume
  if( teflon_annulus_bot_log == NULL ) {
    teflon_annulus_bot_log = new G4LogicalVolume(teflon_annulus_bot, material, "teflon_annulus_bot_log", 0, 0, 0);
    teflon_annulus_bot_log->SetVisAttributes(teflon_annulus_bot_vis_att);
  }
  
  return 1;
}//end ::end AddTeflonAnnulusBot

//Add the aluminum hemisphere, housing for PACES
G4int DetectionSystemPaces::AddAluminumHemisphere()
{
  //material
  G4Material* material = G4Material::GetMaterial("Aluminum");
  if( !material ) {
    G4cout << " ----> Material " << "Aluminum" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* aluminum_hemisphere_vis_att = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  aluminum_hemisphere_vis_att->SetVisibility(true);
  
  //main hemisphere shell
  G4double inner_radius = this->aluminum_hemisphere_inner_radius;
  G4double outer_radius = this->aluminum_hemisphere_outer_radius;
  G4double thickness = abs(outer_radius - inner_radius);
  G4double start_theta = 0.*pi;
  G4double end_theta = pi/2 - start_theta;
  //beam hole rim
  G4double rim_inner_radius = abs(inner_radius - this->aluminum_hemisphere_beam_hole_rim_height);
  G4double rim_outer_radius = inner_radius;
  G4double rim_start_theta = 0.*pi;
  G4double rim_end_theta = asin( (this->aluminum_hemisphere_beam_hole_radius + thickness) / inner_radius );
  //beam hole cut
  G4double beam_hole_inner_radius = 0.*mm;
  G4double beam_hole_outer_radius = this->aluminum_hemisphere_beam_hole_radius;
  G4double beam_hole_half_length_z = outer_radius;
  //cuts where detectors go
  G4double detector_inner_radius = 0.*mm;
  G4double detector_outer_radius = this->detector_assembly_radius + this->cut_clearance;
  G4double detector_half_length_z = this->detector_assembly_thickness/2.0;
  
  //primitive volumes
  G4Sphere* al_he_shell = new G4Sphere("al_he_shell", inner_radius, outer_radius, 0*pi, 2*pi, start_theta, end_theta);
  G4Sphere* al_he_rim = new G4Sphere("al_he_rim", rim_inner_radius, rim_outer_radius, 0*pi, 2*pi, rim_start_theta, rim_end_theta);
  G4Tubs* cut_beam_hole = new G4Tubs("cut_beam_hole", beam_hole_inner_radius, beam_hole_outer_radius, beam_hole_half_length_z, 0*pi, 2*pi);
  G4Tubs* cut_detector = new G4Tubs("cut_detector", detector_inner_radius, detector_outer_radius, detector_half_length_z, 0*pi, 2*pi);
  
  //add rim and beam hole
  G4RotationMatrix* R0 = new G4RotationMatrix; G4ThreeVector T0;
  G4UnionSolid* al_he_shell_rim = new G4UnionSolid("al_he_shell_rim", al_he_shell, al_he_rim, R0, T0);
  G4SubtractionSolid* al_he_shell_rim_hole = new G4SubtractionSolid("al_he_shell_rim_hole", al_he_shell_rim, cut_beam_hole, R0, T0);
  //cut out detectors
  G4int d_i;
  G4double* ptr_pd = this->paces_placement_distance;
  G4double* ptr_pt = this->paces_placement_theta;
  G4double* ptr_pp = this->paces_placement_phi;
  G4double* ptr_ot = this->paces_orientation_theta;
  G4double* ptr_op = this->paces_orientation_phi;
  /*
  G4RotationMatrix** rotate_cut[5];
  G4ThreeVector* move_cut[5];
  for (int i=0; i<5; i++)
  {
    d_i = i;
    G4double d_dist = ptr_pd[d_i] + detector_half_length_z, d_theta = ptr_pt[d_i], d_phi = ptr_pp[d_i];
    (*rotate_cut[i]) = new G4RotationMatrix;
    (*move_cut[i]).setX( d_dist * cos(d_theta) * sin(d_phi) );
    (*move_cut[i]).setY( d_dist * sin(d_theta) * sin(d_phi) );
    (*move_cut[i]).setZ( d_dist *      1.0     * cos(d_phi) );
    G4double ori_theta = d_theta + ptr_ot[d_i] - pi/2; //minus 90 deg
    G4double ori_phi = d_phi + ptr_op[d_i];
    G4ThreeVector yprimeaxis = G4ThreeVector(cos(ori_theta), sin(ori_theta), 0);
    (*rotate_cut[i])->set(yprimeaxis, ori_phi);
  }
  G4SubtractionSolid* al_he_0 = new G4SubtractionSolid("al_he_0", al_he_with_rim, cut_detector, *rotate_cut[0], *move_cut[0]);
  G4SubtractionSolid* al_he_1 = new G4SubtractionSolid("al_he_1", al_he_0, cut_detector, *rotate_cut[1], *move_cut[1]);
  G4SubtractionSolid* al_he_2 = new G4SubtractionSolid("al_he_2", al_he_1, cut_detector, *rotate_cut[2], *move_cut[2]);
  G4SubtractionSolid* al_he_3 = new G4SubtractionSolid("al_he_3", al_he_2, cut_detector, *rotate_cut[3], *move_cut[3]);
  G4SubtractionSolid* aluminum_hemisphere = new G4SubtractionSolid("aluminum_hemisphere", al_he_3, cut_detector, *rotate_cut[4], *move_cut[4]);
  for (int i=0; i<5; i++) delete rotate_cut[i]; //safety
  */
  G4ThreeVector move_cut, yprimeaxis;
  G4RotationMatrix* rotate_cut = new G4RotationMatrix;
  G4double d_dist, d_theta, d_phi;
  G4double ori_theta, ori_phi;
  //one
  d_i = 0;
  d_dist = ptr_pd[d_i] + detector_half_length_z; d_theta = ptr_pt[d_i]; d_phi = ptr_pp[d_i];
  move_cut.setX( d_dist * cos(d_theta) * sin(d_phi) );
  move_cut.setY( d_dist * sin(d_theta) * sin(d_phi) );
  move_cut.setZ( d_dist *      1.0     * cos(d_phi) );
  ori_theta = d_theta + ptr_ot[d_i] - pi/2;
  ori_phi = d_phi + ptr_op[d_i];
  yprimeaxis.set(cos(ori_theta), sin(ori_theta), 0);
  rotate_cut->set(yprimeaxis, ori_phi);
  G4SubtractionSolid* al_he_0 = new G4SubtractionSolid("al_he_0", al_he_shell_rim_hole, cut_detector, rotate_cut, move_cut);
  //two
  d_i = 1;
  d_dist = ptr_pd[d_i] + detector_half_length_z; d_theta = ptr_pt[d_i]; d_phi = ptr_pp[d_i];
  move_cut.setX( d_dist * cos(d_theta) * sin(d_phi) );
  move_cut.setY( d_dist * sin(d_theta) * sin(d_phi) );
  move_cut.setZ( d_dist *      1.0     * cos(d_phi) );
  ori_theta = d_theta + ptr_ot[d_i] - pi/2;
  ori_phi = d_phi + ptr_op[d_i];
  yprimeaxis.set(cos(ori_theta), sin(ori_theta), 0);
  rotate_cut->set(yprimeaxis, ori_phi);
  G4SubtractionSolid* al_he_1 = new G4SubtractionSolid("al_he_1", al_he_0, cut_detector, rotate_cut, move_cut);
  //three
  d_i = 2;
  d_dist = ptr_pd[d_i] + detector_half_length_z; d_theta = ptr_pt[d_i]; d_phi = ptr_pp[d_i];
  move_cut.setX( d_dist * cos(d_theta) * sin(d_phi) );
  move_cut.setY( d_dist * sin(d_theta) * sin(d_phi) );
  move_cut.setZ( d_dist *      1.0     * cos(d_phi) );
  ori_theta = d_theta + ptr_ot[d_i] - pi/2;
  ori_phi = d_phi + ptr_op[d_i];
  yprimeaxis.set(cos(ori_theta), sin(ori_theta), 0);
  rotate_cut->set(yprimeaxis, ori_phi);
  G4SubtractionSolid* al_he_2 = new G4SubtractionSolid("al_he_2", al_he_1, cut_detector, rotate_cut, move_cut);
  //four
  d_i = 3;
  d_dist = ptr_pd[d_i] + detector_half_length_z; d_theta = ptr_pt[d_i]; d_phi = ptr_pp[d_i];
  move_cut.setX( d_dist * cos(d_theta) * sin(d_phi) );
  move_cut.setY( d_dist * sin(d_theta) * sin(d_phi) );
  move_cut.setZ( d_dist *      1.0     * cos(d_phi) );
  ori_theta = d_theta + ptr_ot[d_i] - pi/2;
  ori_phi = d_phi + ptr_op[d_i];
  yprimeaxis.set(cos(ori_theta), sin(ori_theta), 0);
  rotate_cut->set(yprimeaxis, ori_phi);
  G4SubtractionSolid* al_he_3 = new G4SubtractionSolid("al_he_3", al_he_2, cut_detector, rotate_cut, move_cut);
  //five
  d_i = 4;
  d_dist = ptr_pd[d_i] + detector_half_length_z; d_theta = ptr_pt[d_i]; d_phi = ptr_pp[d_i];
  move_cut.setX( d_dist * cos(d_theta) * sin(d_phi) );
  move_cut.setY( d_dist * sin(d_theta) * sin(d_phi) );
  move_cut.setZ( d_dist *      1.0     * cos(d_phi) );
  ori_theta = d_theta + ptr_ot[d_i] - pi/2;
  ori_phi = d_phi + ptr_op[d_i];
  yprimeaxis.set(cos(ori_theta), sin(ori_theta), 0);
  rotate_cut->set(yprimeaxis, ori_phi);
  G4SubtractionSolid* aluminum_hemisphere = new G4SubtractionSolid("aluminum_hemisphere", al_he_3, cut_detector, rotate_cut, move_cut);
  //end
  
  //logical volume
  if( aluminum_hemisphere_log == NULL ) {
    aluminum_hemisphere_log = new G4LogicalVolume(aluminum_hemisphere, material, "aluminum_hemisphere_log", 0, 0, 0);
    aluminum_hemisphere_log->SetVisAttributes(aluminum_hemisphere_vis_att);
  }
  
  return 1;
}//end ::end AddAluminumHemisphere

//Add the delrin hemisphere, housing for the hemisphere
G4int DetectionSystemPaces::AddDelrinHemisphere()
{
  //material
  G4Material* material = G4Material::GetMaterial("Delrin");
  if( !material ) {
    G4cout << " ----> Material " << "Delrin" << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }
  
  //vis attributes
  G4VisAttributes* delrin_hemisphere_vis_att = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  delrin_hemisphere_vis_att->SetVisibility(true);
  
  //main hemisphere
  G4double inner_radius = this->delrin_hemisphere_inner_radius;
  G4double outer_radius = this->delrin_hemisphere_outer_radius;
  G4double start_theta = asin( this->delrin_hemisphere_beam_hole_radius / inner_radius );
  G4double end_theta = pi/2 - start_theta;
  
  //primitive volumes
  G4Sphere* delrin_hemisphere = new G4Sphere("delrin_hemisphere", inner_radius, outer_radius, 0*pi, 2*pi, start_theta, end_theta);
  
  //logical volume
  if( delrin_hemisphere_log == NULL ) {
    delrin_hemisphere_log = new G4LogicalVolume(delrin_hemisphere, material, "delrin_hemisphere_log", 0, 0, 0);
    delrin_hemisphere_log->SetVisAttributes(delrin_hemisphere_vis_att);
  }
  
  return 1;
}//end ::end AddDelrinHemisphere

G4int DetectionSystemPaces::AddScrews()
{
  return 1;
}//end ::end AddScrews
