#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

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

#include "DetectionSystemGammaTracking.hh"

DetectionSystemGammaTracking::DetectionSystemGammaTracking() :
    // Logical Volumes 
    logicShell(0)  
{}

DetectionSystemGammaTracking::~DetectionSystemGammaTracking()
{
    // Logical Volumes 
    delete logicShell;    
}

///////////////////////////////////////////////////////////////////////
// ConstructDetectionSystemGammaTracking builds the DetectionSystemGammaTracking 
// at the origin
///////////////////////////////////////////////////////////////////////
G4int DetectionSystemGammaTracking::Build()//G4SDManager* mySDman)
{ 
  // Build assembly volume
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;

  // detector shell
  
  // search the shell material by its name
  G4String matShellName = "Germanium"; // Should use G4_Ge if it is defined.
  G4Material* matShell = G4Material::GetMaterial(matShellName);
  
  if( !matShell ) {
    G4cout << " ----> Material " << matShellName << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }

  // default parameter values of the detector shell
  shellRmin  = 7.*cm;
  shellThick = 2.*cm;
  
  phi_in     = 0.;
  d_phi      = 360.*deg;
  th_in      = 0.;
  d_th       = 180.*deg;
  
  if( phi_in < 0. ) phi_in = 0.;
  if( phi_in > 360.*deg ) phi_in = 0.;
  
  G4double phi_max = phi_in + d_phi;
  if( phi_max > 360.*deg ) phi_max = 360.*deg;
  
  if( th_in < 0. ) th_in = 0.;
  if( th_in > 180.*deg ) th_in = 0.;
  
  G4double th_max = th_in + d_th;
  if( th_max > 180.*deg ) th_max = 180.*deg;
    
  G4Sphere        *solidShell = new G4Sphere( "shell", shellRmin, shellRmin+shellThick, phi_in, phi_max, th_in, th_max);
  G4LogicalVolume *logicShell = new G4LogicalVolume( solidShell, matShell, "Shell", 0, 0, 0 );

  // Vis Attributes
  G4VisAttributes *pVA = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
  logicShell->SetVisAttributes( pVA );
  
  // Describe Orientation
  G4ThreeVector direction = G4ThreeVector(0,0,1);
  G4ThreeVector move = 0.0 * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateX(0.0);
  rotate->rotateY(0.0);
  rotate->rotateZ(0.0); 

  this->assembly->AddPlacedVolume(logicShell, move, rotate);

  return 1;
}

G4int DetectionSystemGammaTracking::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4RotationMatrix* rotate, G4int detector_number)
{
  G4int detector_copy_ID = 0;

  G4cout << "Gamma Tracking Detector Number = " << detector_number << G4endl;

  G4int copy_number = detector_copy_ID + detector_number;

  assembly->MakeImprint(exp_hall_log, move, rotate, copy_number);

  return 1;
}

// Calculate a direction vector from spherical theta & phi components
G4ThreeVector DetectionSystemGammaTracking::GetDirectionXYZ(G4double theta, G4double phi)
{
  G4double x,y,z;
  x = sin(theta) * cos(phi);
  y = sin(theta) * sin(phi);
  z = cos(theta);
	
  G4ThreeVector direction = G4ThreeVector(x,y,z);
	
  return direction;
}
