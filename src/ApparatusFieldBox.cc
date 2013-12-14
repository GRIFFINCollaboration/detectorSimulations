#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
//#include "LocalMagneticField.hh"
//#include "MagneticField.hh"
//#include "G4TransportationManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ApparatusFieldBox.hh"

ApparatusFieldBox::ApparatusFieldBox() :
    // LogicalVolumes 
    box_log(0)//, magneticField(0)
{ 
  this->box_material = "Vacuum";
  this->box_length_x = 0.0*cm;
  this->box_length_y = 0.0*cm;
  this->box_length_z = 0.0*cm;
}

ApparatusFieldBox::~ApparatusFieldBox()
{
    // LogicalVolumes 
    delete box_log;
//    delete magneticField;
}

G4int ApparatusFieldBox::Build(G4String box_material, G4double box_length_x, G4double box_length_y, G4double box_length_z, G4ThreeVector fieldVector)
{

  //magneticField = new LocalMagneticField(fieldVector);

  this->fieldVector = fieldVector;

  this->box_material = box_material;
  this->box_length_x = box_length_x*mm;
  this->box_length_y = box_length_y*mm;
  this->box_length_z = box_length_z*mm;

  // Build assembly volume
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;

  G4cout << "BuildboxVolume" << G4endl;
  BuildboxVolume();

  return 1;
}

G4int ApparatusFieldBox::PlaceApparatus(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate)
{
  G4int copy_ID = 0;

  this->assembly->MakeImprint(exp_hall_log, move, rotate, copy_ID);

  return 1;
}

G4int ApparatusFieldBox::BuildboxVolume()
{
  G4Material* material = G4Material::GetMaterial(this->box_material);
  if( !material ) {
    G4cout << " ----> Material " << this->box_material << " not found, cannot build!" << G4endl;
    return 0;
  }
  
  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  vis_att->SetVisibility(true);  

  G4ThreeVector move; 
  G4RotationMatrix* rotate = new G4RotationMatrix;

  G4Box* box = Buildbox();

  //logical volume
  if( box_log == NULL )
  {
    box_log = new G4LogicalVolume(box, material, "box_log", 0, 0, 0);

    //box_log->SetFieldManager( magneticField->GetLocalFieldManager(), true ) ;

//    G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0., 1.0*tesla, 0.)); //Create Field
//    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();//Set as default field
//    fieldMgr->SetDetectorField(magField);
//    fieldMgr->CreateChordFinder(magField);

//    box_log->SetFieldManager(fieldMgr,true);


    box_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(box_log, move, rotate);

  return 1;
}

G4Box* ApparatusFieldBox::Buildbox()
{
  G4double half_length_x = this->box_length_x/2.0;
  G4double half_length_y = this->box_length_y/2.0;
  G4double half_length_z = this->box_length_z/2.0;

  G4Box* box = new G4Box("box", half_length_x, half_length_y, half_length_z);
  return box;
}

//Calculate a direction vector from spherical theta & phi components
G4ThreeVector ApparatusFieldBox::GetDirectionXYZ(G4double theta, G4double phi)
{
  G4double x,y,z;
  x = sin(theta) * cos(phi);
  y = sin(theta) * sin(phi);
  z = cos(theta);
	
  G4ThreeVector direction = G4ThreeVector(x,y,z);
	
  return direction;
}
