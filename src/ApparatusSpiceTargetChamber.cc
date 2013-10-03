#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
//#include "SensitiveDetector.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
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

//#include "DetectorConstruction.hh"

//#include "G4SDManager.hh"
//#include "DetectableSD.hh"

//#include "G4FieldManager.hh"
//#include "G4UniformMagField.hh"
//#include "MagneticField.hh"
//#include "G4TransportationManager.hh"

#include "ApparatusSpiceTargetChamber.hh"

///////////////////////////////////////////////////////////////////////
// The ::ApparatusSpiceTargetChamber constructor initiates all the Logical 
// and Physical Volumes used in the Target Chamber geometery, and the 
//  ::ApparatusSpiceTargetChamber  destructor deletes them from the stack 
// when they go out of scope
///////////////////////////////////////////////////////////////////////
ApparatusSpiceTargetChamber::ApparatusSpiceTargetChamber() :
   // LogicalVolumes
   target_chamber_cylinder_log(0), 
   target_chamber_sphere_log(0), 
   transport_chamber_al_plate_log(0), 
   transport_chamber_wide_cylinder_log(0),
   transport_chamber_narrow_cylinder_log(0),
   photon_shield_log(0), 
   transport_magnet_log(0), 
   inner_transport_magnet_log(0), 
   collector_magnet_log(0),
   det_cham_end_plate_log(0), 
   det_cham_plate_log(0),
   det_cham_tube_log(0),
   target_chamber_vacuum_log(0), 
   target_chamber_sph_vacuum_log(0), 
   transport_wide_vacuum_log(0), 
   transport_narrow_vacuum_log(0), 
   al_plate_vacuum_log(0),
   det_cham_plate_vacuum_log(0), 
   det_cham_tube_vacuum_log(0), 
   inner_transport_magnet_cover_log(0),

   // PhysicalVolumes
   target_chamber_cylinder_phys(0),
   target_chamber_sphere_phys(0), 
   transport_chamber_al_plate_phys(0), 
   transport_chamber_wide_cylinder_phys(0),
   transport_chamber_narrow_cylinder_phys(0), 
   photon_shield_phys(0), 
   transport_magnet_phys(0),
   inner_transport_magnet_phys(0),
   collector_magnet_phys(0), 
   det_cham_end_plate_phys(0),
   det_cham_plate_phys(0), 
   det_cham_tube_phys(0), 
   inner_transport_magnet_cover_phys(0),
   target_chamber_vacuum_phys(0),
   target_chamber_sph_vacuum_phys(0), 
   transport_wide_vacuum_phys(0), 
   transport_narrow_vacuum_phys(0), 
   al_plate_vacuum_phys(0),
   det_cham_plate_vacuum_phys(0), 
   det_cham_tube_vacuum_phys(0)//,

   // Fields
//   solMagField(0)

{ 
//  solMagField = new MagneticField();		// Global field is set to zero
//  cout << "after MagneticField();" << endl;

  /////////////////////////////////////////////////////////////////////
  // ApparatusSpiceTargetChamber Physical Properties 
  /////////////////////////////////////////////////////////////////////

  this->vacuum_material                                 = "Vacuum";

  this->target_chamber_cylinder_material 	       				= "Delrin";
  this->target_chamber_sphere_material 	                = "Delrin";
  this->transport_chamber_al_plate_material             = "Aluminum";
  this->transport_chamber_wide_cylinder_material        = "Delrin";
  this->transport_chamber_narrow_cylinder_material      = "Delrin";
  this->photon_shield_material                          = "Tungsten";
  this->transport_magnet_material                       = "NdFeB";
  this->inner_transport_magnet_material                 = "NdFeB";
  this->inner_transport_magnet_cover_material           = "Delrin";
  this->collector_magnet_material                       = "NdFeB";
  
  this->det_cham_end_plate_material                     = "Aluminum";
  this->det_cham_plate_material                         = "Aluminum";
  this->det_cham_tube_material                          = "Aluminum";

  this->target_chamber_cylinder_length                  = 100.*mm;
  this->target_chamber_inner_radius                     = 92.*mm;
  this->target_chamber_outer_radius                     = 100.*mm;
  this->Al_plate_thickness                              = 10.*mm;
  this->Al_plate_inner_radius                           = 50.*mm;
  this->Al_plate_outer_radius                           = 120.*mm;
  this->trans_chamber_wide_cylinder_length              = 150.*mm;
  this->trans_chamber_wide_cylinder_inner_radius        = 60.*mm;
  this->trans_chamber_wide_cylinder_outer_radius        = 70.*mm;
  this->trans_chamber_narrow_cylinder_length            = 290.*mm;
  this->trans_chamber_narrow_cylinder_inner_radius      = 5.*mm;
  this->trans_chamber_narrow_cylinder_outer_radius      = 8.*mm;
  this->photon_shield_inner_radius                      = 5.*mm;
  this->photon_shield_tip_outer_radius                  = 9.8*mm;
  this->photon_shield_base_outer_radius                 = 21.*mm;
  this->photon_shield_length                            = 24.*mm;
  this->transport_magnet_length                         = 150.*mm;
  this->transport_magnet_inner_radius                   = 63.*mm;
  this->transport_magnet_outer_radius                   = 89.*mm;
  this->inner_transport_magnet_length                   = 150.*mm;
  this->inner_transport_magnet_inner_radius             = 11.*mm;
  this->inner_transport_magnet_outer_radius             = 18.*mm;
  this->inner_transport_magnet_cover_length             = 150.*mm;
  this->inner_transport_magnet_cover_inner_radius       = 18.*mm;
  this->inner_transport_magnet_cover_outer_radius       = 19.*mm;
  this->collector_magnet_length_x                       = 40.*mm;
  this->collector_magnet_length_y_inner                 = 6.*mm;
  this->collector_magnet_length_y_outer                 = 6.*mm;
  this->collector_magnet_length_z                       = 90.*mm;
  this->det_cham_end_plate_thickness                    = 10.*mm;
  this->det_cham_end_plate_inner_radius                 = this->trans_chamber_narrow_cylinder_outer_radius;
  this->det_cham_end_plate_outer_radius                 = 120.*mm;
  this->det_cham_plate_thickness                        = 10.*mm;
  this->det_cham_plate_inner_radius                     = 50.*mm;
  this->det_cham_plate_outer_radius                     = 120.*mm;
  this->det_cham_tube_thickness                         = 30.*mm;
  this->det_cham_tube_inner_radius                      = 110.*mm;
  this->det_cham_tube_outer_radius                      = 120.*mm;


  this->target_chamber_start_position                   = -100.*mm;
  this->al_plate_position                               = -100.*mm;
  this->trans_chamber_wide_cylinder_position            = -110.*mm;
  this->trans_chamber_narrow_cylinder_position          = -45.*mm;
  this->photon_shield_position                          = -21.*mm;
  this->transport_magnet_position                       = -110.*mm;
  this->inner_transport_magnet_position                 = -110.*mm;
  this->inner_transport_magnet_cover_position           = -110.*mm;
  this->collector_magnet_position_z                     = 0.*mm;
  this->collector_magnet_position_y                     = 26.0*mm;//45.0*mm;
  this->det_cham_end_plate_position                     =-300.*mm;
  this->det_cham_plate_position                         =-260.*mm;
  this->det_cham_tube_position                          =-270.*mm;
 
  this->target_chamber_cross_section                    = 360.*deg;

  this->solMagFieldStrength                             = 0.0*tesla;


}// end ::ApparatusSpiceTargetChamber

ApparatusSpiceTargetChamber::~ApparatusSpiceTargetChamber()
{
   // LogicalVolumes in ConstructApparatusSpiceTargetChamber
   delete target_chamber_cylinder_log;
   delete target_chamber_sphere_log;
   delete transport_chamber_al_plate_log;
   delete transport_chamber_wide_cylinder_log;
   delete transport_chamber_narrow_cylinder_log;
   delete photon_shield_log;
   delete transport_magnet_log;
   delete inner_transport_magnet_log;
   delete inner_transport_magnet_cover_log;
   delete collector_magnet_log;
   delete det_cham_end_plate_log;
   delete det_cham_plate_log;
   delete det_cham_tube_log;

   delete target_chamber_vacuum_log;
   delete target_chamber_sph_vacuum_log;
   delete det_cham_plate_vacuum_log;
   delete det_cham_tube_vacuum_log;
  // delete expHallLog;
   
   // PhysicalVolumes in ConstructApparatusSpiceTargetChamber
   delete target_chamber_cylinder_phys;
   delete target_chamber_sphere_phys;
   delete transport_chamber_al_plate_phys;
   delete transport_chamber_wide_cylinder_phys;
   delete transport_chamber_narrow_cylinder_phys;
   delete photon_shield_phys;
   delete transport_magnet_phys;
   delete inner_transport_magnet_phys;
   delete inner_transport_magnet_cover_phys;
   delete collector_magnet_phys;
   delete det_cham_end_plate_phys;
   delete det_cham_plate_phys;
   delete det_cham_tube_phys;

   delete target_chamber_vacuum_phys;
   delete target_chamber_sph_vacuum_phys;
   delete det_cham_plate_vacuum_phys;
   delete det_cham_tube_vacuum_phys;

}// end ::~ApparatusSpiceTargetChamber

///////////////////////////////////////////////////////////////////////
//ConstructApparatusSpiceTargetChamber builds the Target Chamber at the origin
///////////////////////////////////////////////////////////////////////
void ApparatusSpiceTargetChamber::Build(G4LogicalVolume* exp_hall_log)
{ 
   this->expHallLog = exp_hall_log;

   BuildApparatusSpiceTargetChamberCylinder();             //Includes a vacuum 
   BuildApparatusSpiceTargetChamberSphere();               //Includes a vacuum 
   BuildTransportChamberAlPlate();           //Includes a vacuum   
   BuildTransportChamberWideCylinder();      //Includes a vacuum 
   BuildTransportChamberNarrowCylinder();    //Includes a vacuum 

 //  BuildPhotonShield();
 //  BuildTransportMagnets();                //As per new lens design
 //  BuildInnerTransportMagnets();           //As per new lens design 
//   BuildInnerTransportMagnetsCover();      //As per new lens design

   BuildCollectorMagnet();

   BuildDetectorChamberAlEndPlate();         //Includes a vacuum 
   BuildDetectorChamberAlPlate();            //Includes a vacuum 
   BuildDetectorChamberAlTube();             //Includes a vacuum 

   PlaceApparatusSpiceTargetChamberCylinder();             
   PlaceApparatusSpiceTargetChamberSphere();               
   PlaceTransportChamberAlPlate();           //Places a vacuum   Must include these first as vacuum material overwrites existing materials
   PlaceTransportChamberWideCylinder();
   PlaceTransportChamberNarrowCylinder();   
   PlaceDetectorChamberAlEndPlate();
   PlaceDetectorChamberAlPlate();
   PlaceDetectorChamberAlTube();

 //  PlacePhotonShield();
  // PlaceTransportMagnets();                //As per new lens design
  // PlaceInnerTransportMagnets();           //As per new lens design
  // PlaceInnerTransportMagnetsCover();      //As per new lens design

		for (G4int copyID = 0; copyID < 4; copyID++)
   		PlaceCollectorMagnet( copyID ) ;

}//end ::Build


///////////////////////////////////////////////////////////////////////
//methods used to build shapes
///////////////////////////////////////////////////////////////////////
void ApparatusSpiceTargetChamber::BuildApparatusSpiceTargetChamberCylinder()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(DELRIN_COL));
   vis_att->SetVisibility(true);

   G4double length = this->target_chamber_cylinder_length; 

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;
                 
   G4double inner_radius = this->target_chamber_inner_radius;
   G4double outer_radius = this->target_chamber_outer_radius;

   G4Tubs* target_chamber_cylinder = new G4Tubs("target_chamber_cylinder", inner_radius, outer_radius, length/2., startPhi, endPhi);

   // Establish logical volumes
   G4Material* target_chamber_cylinder_material = G4Material::GetMaterial(this->target_chamber_cylinder_material);
   target_chamber_cylinder_log = new G4LogicalVolume(target_chamber_cylinder, target_chamber_cylinder_material, "target_chamber_cylinder_log", 0, 0, 0);
   target_chamber_cylinder_log->SetVisAttributes(vis_att);

//   //Create Vacuum to fill Cylinder
//   G4Tubs* target_chamber_vacuum = new G4Tubs("target_chamber_vacuum", 0.*mm, inner_radius, length/2., startPhi, endPhi);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   target_chamber_vacuum_log = new G4LogicalVolume(target_chamber_vacuum, vacuum_material, "target_chamber_vacuum_log",0,0,0);
//   target_chamber_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);	

}//end ::ApparatusSpiceTargetChamberCylinder



void ApparatusSpiceTargetChamber::BuildApparatusSpiceTargetChamberSphere()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(DELRIN_COL));
   vis_att->SetVisibility(true);

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double startTheta = 0;
   G4double endTheta = 1./2.*pi;

   G4double inner_radius = this->target_chamber_inner_radius;
   G4double outer_radius = this->target_chamber_outer_radius;

   G4Sphere* target_chamber_sphere = new G4Sphere("target_chamber_sphere",  inner_radius, outer_radius, startPhi, endPhi, startTheta, endTheta);

   // Establish logical volumes
   G4Material* target_chamber_sphere_material = G4Material::GetMaterial(this->target_chamber_sphere_material);
   target_chamber_sphere_log = new G4LogicalVolume(target_chamber_sphere, target_chamber_sphere_material, "target_chamber_sphere_log", 0, 0, 0);
   target_chamber_sphere_log->SetVisAttributes(vis_att);

//   //Create Vacuum to fill Cylinder
//   G4Sphere* target_chamber_sph_vacuum = new G4Sphere("target_chamber_sph_vacuum", 0., inner_radius, startPhi, endPhi, startTheta, endTheta);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   target_chamber_sph_vacuum_log = new G4LogicalVolume(target_chamber_sph_vacuum, vacuum_material, "target_chamber_vacuum_log",0,0,0);
//   target_chamber_sph_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);	

}//end ::ApparatusSpiceTargetChamberCylinder



void ApparatusSpiceTargetChamber::BuildTransportChamberAlPlate()
{

// Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
   vis_att->SetVisibility(true);

   G4double length = this->Al_plate_thickness;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;


   G4double inner_radius = this->Al_plate_inner_radius;
   G4double outer_radius = this->Al_plate_outer_radius;

   G4Tubs* transport_chamber_al_plate = new G4Tubs("transport_chamber_al_plate", inner_radius, outer_radius, length/2., startPhi, endPhi);
   
   // Establish logical volumes
   G4Material* transport_chamber_al_plate_material = G4Material::GetMaterial(this->transport_chamber_al_plate_material);
   transport_chamber_al_plate_log = new G4LogicalVolume(transport_chamber_al_plate, transport_chamber_al_plate_material, "transport_chamber_al_plate_log", 0, 0, 0);
   transport_chamber_al_plate_log->SetVisAttributes(vis_att);

//   //Create Vacuum to fill Cylinder
//   G4Tubs* al_plate_vacuum = new G4Tubs("al_plate_vacuum", 0., inner_radius, length/2., startPhi, endPhi);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   al_plate_vacuum_log = new G4LogicalVolume(al_plate_vacuum, vacuum_material, "al_plate_vacuum_log",0,0,0);
//   al_plate_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);	

}//end ::TransportChamberAlPlate


void ApparatusSpiceTargetChamber::BuildTransportChamberWideCylinder()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(DELRIN_COL));
   vis_att->SetVisibility(true);


   G4double length = this->trans_chamber_wide_cylinder_length;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double inner_radius = this->trans_chamber_wide_cylinder_inner_radius;   
   G4double outer_radius = this->trans_chamber_wide_cylinder_outer_radius;

   G4Tubs* transport_chamber_wide_cylinder = new G4Tubs("transport_chamber_wide_cylinder", inner_radius, outer_radius, length/2., startPhi, endPhi);    
   // Establish logical volumes
   G4Material* transport_chamber_wide_cylinder_material = G4Material::GetMaterial(this->transport_chamber_wide_cylinder_material);
   transport_chamber_wide_cylinder_log = new G4LogicalVolume(transport_chamber_wide_cylinder, transport_chamber_wide_cylinder_material, "transport_chamber_wide_cylinder_log", 0, 0, 0);
   transport_chamber_wide_cylinder_log->SetVisAttributes(vis_att);

//   //Create Vacuum to fill Cylinder
//   G4Tubs* transport_wide_vacuum = new G4Tubs("transport_wide_vacuum", 0.*mm, inner_radius, length/2., startPhi, endPhi);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   transport_wide_vacuum_log = new G4LogicalVolume(transport_wide_vacuum, vacuum_material, "transport_wide_vacuum_log",0,0,0);
//   transport_wide_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);

//   //Fill Vacuum with B-Field
//   G4double zstrength = this->solMagFieldStrength;

//   G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0., zstrength)); //Create Field
//   G4FieldManager*fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();//Set as default field
//   fieldMgr->SetDetectorField(magField);
//   fieldMgr->CreateChordFinder(magField);

//   transport_wide_vacuum_log->SetFieldManager(fieldMgr,true);


}//end ::TransportChamberWideCylinder


void ApparatusSpiceTargetChamber::BuildTransportChamberNarrowCylinder()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(DELRIN_COL));
   vis_att->SetVisibility(true);

   G4double length = this->trans_chamber_narrow_cylinder_length;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double inner_radius = this->trans_chamber_narrow_cylinder_inner_radius;   
   G4double outer_radius = this->trans_chamber_narrow_cylinder_outer_radius;

   G4Tubs* transport_chamber_narrow_cylinder = new G4Tubs("transport_chamber_narrow_cylinder", inner_radius, outer_radius, length/2., startPhi, endPhi);

   // Establish logical volumes
   G4Material* transport_chamber_narrow_cylinder_material = G4Material::GetMaterial(this->transport_chamber_narrow_cylinder_material);
   transport_chamber_narrow_cylinder_log = new G4LogicalVolume(transport_chamber_narrow_cylinder, transport_chamber_narrow_cylinder_material, "transport_chamber_narrow_cylinder_log", 0, 0, 0);
   transport_chamber_narrow_cylinder_log->SetVisAttributes(vis_att);
 
//  //Create Vacuum to fill Cylinder
//   G4Tubs* transport_narrow_vacuum = new G4Tubs("transport_narrow_vacuum", 0., inner_radius, length/2., startPhi, endPhi);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   transport_narrow_vacuum_log = new G4LogicalVolume(transport_narrow_vacuum, vacuum_material, "transport_narrow_vacuum_log",0,0,0);
//   transport_narrow_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);

}//end ::TransportChamberNarrowCylinder


void ApparatusSpiceTargetChamber::BuildPhotonShield()
{
    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(PB_COL));
    vis_att->SetVisibility(true);

    G4double length = this->photon_shield_length;

    G4double startPhi = 0;
    G4double endPhi = this->target_chamber_cross_section;

    G4double tip_inner_radius = this->photon_shield_inner_radius;
    G4double base_inner_radius = this->photon_shield_inner_radius;
    G4double tip_outer_radius = this->photon_shield_tip_outer_radius;
    G4double base_outer_radius = this->photon_shield_base_outer_radius;

    G4Cons* photon_shield = new G4Cons("photon_shield", base_inner_radius, base_outer_radius, tip_inner_radius, tip_outer_radius, length/2., startPhi, endPhi);

    // Establish logical volumes
    G4Material* photon_shield_material = G4Material::GetMaterial(this->photon_shield_material);
    photon_shield_log = new G4LogicalVolume(photon_shield, photon_shield_material, "photon_shield_log", 0, 0, 0);
    photon_shield_log->SetVisAttributes(vis_att);

}//end ::TransportChamberNarrowCylinder


void ApparatusSpiceTargetChamber::BuildTransportMagnets()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(NDFEB_COL));
   vis_att->SetVisibility(true);
    
   G4double length = this->transport_magnet_length;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double inner_radius = this->transport_magnet_inner_radius;   
   G4double outer_radius = this->transport_magnet_outer_radius;

   G4Tubs* transport_magnet = new G4Tubs("transport_magnet", inner_radius, outer_radius, length/2., startPhi, endPhi);

   // Establish logical volumes
   G4Material* transport_magnet_material = G4Material::GetMaterial(this->transport_magnet_material);
   transport_magnet_log = new G4LogicalVolume(transport_magnet, transport_magnet_material, "transport_magnet_log", 0, 0, 0);
   transport_magnet_log->SetVisAttributes(vis_att);


}//end ::TransportChamberNarrowCylinder

void ApparatusSpiceTargetChamber::BuildInnerTransportMagnets()
{
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Color(NDFEB_COL));
  vis_att->SetVisibility(true);

  G4double length = this->inner_transport_magnet_length;

  G4double startPhi = 0;
  G4double endPhi = 360.*deg;

  G4double inner_radius = this->inner_transport_magnet_inner_radius;
  G4double outer_radius = this->inner_transport_magnet_outer_radius;

  G4Tubs* inner_transport_magnet = new G4Tubs("inner_transport_magnet", inner_radius, outer_radius, length/2., startPhi, endPhi);

  //Establish Logical Volumes
  G4Material* inner_transport_magnet_material = G4Material::GetMaterial(this->inner_transport_magnet_material);
  inner_transport_magnet_log = new G4LogicalVolume(inner_transport_magnet, inner_transport_magnet_material, "inner_transport_magnet_log", 0, 0, 0);
  inner_transport_magnet_log->SetVisAttributes(vis_att);


}//end ::InnerTransportMagnets


void ApparatusSpiceTargetChamber::BuildInnerTransportMagnetsCover()
{
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Color(DELRIN_COL));
  vis_att->SetVisibility(true);

  G4double length = this->inner_transport_magnet_cover_length;

  G4double startPhi = 0;
  G4double endPhi = 360.*deg;

  G4double inner_radius = this->inner_transport_magnet_cover_inner_radius;
  G4double outer_radius = this->inner_transport_magnet_cover_outer_radius;

  G4Tubs* inner_transport_magnet_cover = new G4Tubs("inner_transport_magnet_cover", inner_radius, outer_radius, length/2., startPhi, endPhi);

  //Establish Logical Volumes
  G4Material* inner_transport_magnet_cover_material = G4Material::GetMaterial(this->inner_transport_magnet_cover_material);
  inner_transport_magnet_cover_log = new G4LogicalVolume(inner_transport_magnet_cover, inner_transport_magnet_cover_material, "inner_transport_magnet_cover_log", 0, 0, 0);
  inner_transport_magnet_cover_log->SetVisAttributes(vis_att);


}//end ::InnerTransportMagnetsCover

void ApparatusSpiceTargetChamber::BuildCollectorMagnet()
{
   // Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(NDFEB_COL));
   vis_att->SetVisibility(true);

   G4double length_x       = this->collector_magnet_length_x/2.;
   G4double length_y_inner = this->collector_magnet_length_y_inner/2.;
   G4double length_y_outer = this->collector_magnet_length_y_outer/2.;
   G4double length_z       = this->collector_magnet_length_z/2.;

   G4Trd* collector_magnet = new G4Trd("collector_magnet", length_z, length_z, length_y_inner, length_y_outer, length_x);

   // Establish logical volumes
   G4Material* collector_magnet_material = G4Material::GetMaterial(this->collector_magnet_material);
   collector_magnet_log = new G4LogicalVolume(collector_magnet, collector_magnet_material, "collector_magnet_log", 0, 0, 0);
   collector_magnet_log->SetVisAttributes(vis_att);


}//end ::CollectorMagnet


void ApparatusSpiceTargetChamber::BuildDetectorChamberAlEndPlate()
{

// Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
   vis_att->SetVisibility(true);

   G4double length = this->det_cham_end_plate_thickness;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;


   G4double inner_radius = this->det_cham_end_plate_inner_radius;
   G4double outer_radius = this->det_cham_end_plate_outer_radius;

   G4Tubs* det_cham_end_plate = new G4Tubs("det_cham_end_plate", inner_radius, outer_radius, length/2., startPhi, endPhi);
   
   // Establish logical volumes
   G4Material* det_cham_end_plate_material = G4Material::GetMaterial(this->det_cham_end_plate_material);
   det_cham_end_plate_log = new G4LogicalVolume(det_cham_end_plate, det_cham_end_plate_material, "det_cham_end_plate_log", 0, 0, 0);
   det_cham_end_plate_log->SetVisAttributes(vis_att);

}//end ::DetectorChamberAlEndPlate()


void ApparatusSpiceTargetChamber::BuildDetectorChamberAlPlate()
{

// Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
   vis_att->SetVisibility(true);

   G4double length = this->det_cham_plate_thickness;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double inner_radius = this->det_cham_plate_inner_radius;
   G4double outer_radius = this->det_cham_plate_outer_radius;

   G4Tubs* det_cham_plate = new G4Tubs("det_cham_plate", inner_radius, outer_radius, length/2., startPhi, endPhi);
   
   // Establish logical volumes
   G4Material* det_cham_plate_material = G4Material::GetMaterial(this->det_cham_plate_material);
   det_cham_plate_log = new G4LogicalVolume(det_cham_plate, det_cham_plate_material, "det_cham_plate_log", 0, 0, 0);
   det_cham_plate_log->SetVisAttributes(vis_att);

//  //Create Vacuum to fill Cylinder
//   G4Tubs* det_cham_plate_vacuum = new G4Tubs("det_cham_plate_vacuum", 0., inner_radius, length/2., startPhi, endPhi);
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   det_cham_plate_vacuum_log = new G4LogicalVolume(det_cham_plate_vacuum, vacuum_material, "det_cham_plate_vacuum_log",0,0,0);
//   det_cham_plate_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);

}//end ::DetectorChamberAlPlate()


void ApparatusSpiceTargetChamber::BuildDetectorChamberAlTube()
{
// Set visualization attributes
   G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(AL_COL));
   vis_att->SetVisibility(true);

   G4double length = this->det_cham_tube_thickness;

   G4double startPhi = 0;
   G4double endPhi = this->target_chamber_cross_section;

   G4double inner_radius = this->det_cham_tube_inner_radius;
   G4double outer_radius = this->det_cham_tube_outer_radius;

   G4Tubs* det_cham_tube = new G4Tubs("det_cham_tube", inner_radius, outer_radius, length/2., startPhi, endPhi);
   
   // Establish logical volumes
   G4Material* det_cham_tube_material = G4Material::GetMaterial(this->det_cham_tube_material);
   det_cham_tube_log = new G4LogicalVolume(det_cham_tube, det_cham_tube_material, "det_cham_tube_log", 0, 0, 0);
   det_cham_tube_log->SetVisAttributes(vis_att);

//  //Create Vacuum to fill Cylinder
//   G4Tubs* det_cham_tube_vacuum = new G4Tubs("det_cham_tube_vacuum",  this->trans_chamber_narrow_cylinder_outer_radius, inner_radius, length/2., startPhi, endPhi); //inner radius avoids transport cylinder
//   G4Material* vacuum_material = G4Material::GetMaterial(this->vacuum_material);
//   det_cham_tube_vacuum_log = new G4LogicalVolume(det_cham_tube_vacuum, vacuum_material, "det_cham_tube_vacuum_log",0,0,0);
//   det_cham_tube_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);

}//end ::DetectorChamberAlTube()



///////////////////////////////////////////////////////////////////////
//methods used in Build()
///////////////////////////////////////////////////////////////////////

void ApparatusSpiceTargetChamber::PlaceApparatusSpiceTargetChamberCylinder()
{

   G4double z_position = this->target_chamber_start_position + this->target_chamber_cylinder_length/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   target_chamber_cylinder_phys = new G4PVPlacement(0, move, target_chamber_cylinder_log, "target_chamber_cylinder_phys", expHallLog, false, 0);
   //Place Vacuum
   //target_chamber_vacuum_phys = new G4PVPlacement(0, move, target_chamber_vacuum_log, "target_chamber_vacuum_phys", expHallLog, false, 0);

}//end ::PlaceApparatusSpiceTargetChamberCylinder()



void ApparatusSpiceTargetChamber::PlaceApparatusSpiceTargetChamberSphere()
{

   G4double z_position = 0.0 ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   target_chamber_sphere_phys = new G4PVPlacement(0, move, target_chamber_sphere_log, "target_chamber_sphere_phys", expHallLog, false, 0);
   //Place Vacuum
   //target_chamber_sph_vacuum_phys = new G4PVPlacement(0, move, target_chamber_sph_vacuum_log, "target_chamber_sph_vacuum_phys", expHallLog, false, 0);


}//end ::PlaceApparatusSpiceTargetChamberCylinder()




void ApparatusSpiceTargetChamber::PlaceTransportChamberAlPlate()
{

   G4double z_position = this->al_plate_position - this->Al_plate_thickness/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   transport_chamber_al_plate_phys = new G4PVPlacement(0, move, transport_chamber_al_plate_log, "transport_chamber_al_plate_phys", expHallLog, false, 0);

   //Place Vacuum
   //al_plate_vacuum_phys = new G4PVPlacement(0, move, al_plate_vacuum_log, "al_plate_vacuum_phys", expHallLog, false, 0);


}//end ::PlaceTransportChamberAlPlate()



void ApparatusSpiceTargetChamber::PlaceTransportChamberWideCylinder()
{

   G4double z_position = this->trans_chamber_wide_cylinder_position - this->trans_chamber_wide_cylinder_length/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   transport_chamber_wide_cylinder_phys = new G4PVPlacement(0, move, transport_chamber_wide_cylinder_log, "transport_chamber_wide_cylinder_phys", expHallLog, false, 0);
   //Place Vacuum
   //transport_wide_vacuum_phys = new G4PVPlacement(0, move, transport_wide_vacuum_log, "transport_wide_vacuum_phys", expHallLog, false, 0);

}//end ::PlaceTransportChamberWideCylinder()


void ApparatusSpiceTargetChamber::PlaceTransportChamberNarrowCylinder()
{

   G4double z_position = this->trans_chamber_narrow_cylinder_position - this->trans_chamber_narrow_cylinder_length/2.;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   transport_chamber_narrow_cylinder_phys = new G4PVPlacement(0, move, transport_chamber_narrow_cylinder_log, "transport_chamber_narrow_cylinder_phys", expHallLog, false, 0);
   //Place Vacuum
   //transport_narrow_vacuum_phys = new G4PVPlacement(0, move, transport_narrow_vacuum_log, "transport_narrow_vacuum_phys", expHallLog, false, 0);

}//end ::PlaceTransportChamberWideCylinder()


void ApparatusSpiceTargetChamber::PlacePhotonShield()
{
   G4double z_position = this->photon_shield_position - this->photon_shield_length/2. ;

   G4ThreeVector move(0, 0, z_position);
 
   // Establish physical volumes
   photon_shield_phys = new G4PVPlacement(0, move, photon_shield_log, "photon_shield_phys",expHallLog, false, 0);

}//end ::PlaceTransportChamberWideCylinder()


void ApparatusSpiceTargetChamber::PlaceTransportMagnets()
{
   G4double x_position = 0.0 ;
   G4double y_position = 0.0 ;
   G4double z_position = this->transport_magnet_position - this->transport_magnet_length/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   transport_magnet_phys = new G4PVPlacement(0, move, transport_magnet_log, "transport_magnet_phys", expHallLog, false, 0);


}//end ::PlaceTransportMagnets()


void ApparatusSpiceTargetChamber::PlaceInnerTransportMagnets()
{
    G4double x_position = 0.0;
    G4double y_position = 0.0;
    G4double z_position = this->inner_transport_magnet_position - this->inner_transport_magnet_length/2.;

    G4ThreeVector move(x_position, y_position, z_position);

    //Establish Physical Volume
    inner_transport_magnet_phys = new G4PVPlacement(0, move, inner_transport_magnet_log, "inner_transport_magnet_phys", expHallLog, false, 0);

}//end ::PlaceInnerTransportMagnets()


void ApparatusSpiceTargetChamber::PlaceInnerTransportMagnetsCover()
{
    G4double x_position = 0.0;
    G4double y_position = 0.0;
    G4double z_position = this->inner_transport_magnet_cover_position - this->inner_transport_magnet_cover_length/2.;

    G4ThreeVector move(x_position, y_position, z_position);

    //Establish Physical Volume
    inner_transport_magnet_cover_phys = new G4PVPlacement(0, move, inner_transport_magnet_cover_log, "inner_transport_magnet_cover_phys", expHallLog, false, 0);

}//end ::PlaceInnerTransportMagnets()


void ApparatusSpiceTargetChamber::PlaceCollectorMagnet(G4int copyID)
{
   G4double z_position = this->collector_magnet_position_z - this->collector_magnet_length_z/2. ;
 
   G4RotationMatrix* rotation;

   rotation = RotateMagnets(copyID);

   G4double radial_position = this->collector_magnet_position_y + this->collector_magnet_length_x/2.; //Convert to Cylindricals      

   G4ThreeVector move = TranslateMagnets(copyID, radial_position, z_position);
      
   // Establish physical volumes
   collector_magnet_phys = new G4PVPlacement(rotation, move, collector_magnet_log, "collector_magnet_phys", expHallLog, false, copyID);

}//end ::PlaceCollectorMagnet()


void ApparatusSpiceTargetChamber::PlaceDetectorChamberAlEndPlate()
{

   G4double z_position = this->det_cham_end_plate_position - this->det_cham_end_plate_thickness/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   det_cham_end_plate_phys = new G4PVPlacement(0, move, det_cham_end_plate_log, "det_cham_end_plate_phys", expHallLog, false, 0);


}//end ::PlaceDetectorChamberAlEndPlate()


void ApparatusSpiceTargetChamber::PlaceDetectorChamberAlPlate()
{

   G4double z_position = this->det_cham_plate_position - this->det_cham_plate_thickness/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   det_cham_plate_phys = new G4PVPlacement(0, move, det_cham_plate_log, "det_cham_plate_phys", expHallLog, false, 0);

 //Place Vacuum
 //  det_cham_plate_vacuum_phys = new G4PVPlacement(0, move, det_cham_plate_vacuum_log, "det_cham_plate_vacuum_phys", expHallLog, false, 0);

}//end ::PlaceDetectorChamberAlPlate()


void ApparatusSpiceTargetChamber::PlaceDetectorChamberAlTube()
{

   G4double z_position = this->det_cham_tube_position - this->det_cham_tube_thickness/2. ;

   G4ThreeVector move(0, 0, z_position);
      
   // Establish physical volumes
   det_cham_tube_phys = new G4PVPlacement(0, move, det_cham_tube_log, "det_cham_tube_phys", expHallLog, false, 0);

 //Place Vacuum
 //  det_cham_tube_vacuum_phys = new G4PVPlacement(0, move, det_cham_tube_vacuum_log, "det_cham_tube_vacuum_phys", expHallLog, false, 0);

}//end ::PlaceDetectorChamberAlTube()



///////////////////////////////////////////////////////
//Functions Within Target Chamber
//////////////////////////////////////////////////////

G4RotationMatrix* ApparatusSpiceTargetChamber::RotateMagnets(G4int copyID)
 {

   G4RotationMatrix* rotate = new G4RotationMatrix;
   rotate->rotateZ(-copyID*90.*deg);
   rotate->rotateY(-90*deg);

   return rotate;

 }

G4ThreeVector ApparatusSpiceTargetChamber::TranslateMagnets(G4int copyID, G4double radial_position, G4double z_position)
{

  G4double x_position = radial_position*cos(copyID * 90.*deg);
  G4double y_position = radial_position*sin(copyID * 90.*deg);

  return G4ThreeVector(x_position, y_position, z_position);

}

