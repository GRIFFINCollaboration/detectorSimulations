#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SensitiveDetector.hh"

#include "G4Material.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"

#include "DetectionSystemSpice.hh"

DetectionSystemSpice::DetectionSystemSpice() :
  // LogicalVolumes
  siInnerGuardRing_log(0),
  siOuterGuardRing_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    this->wafer_material             = "Silicon";

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    this->siDetCrystalOuterDiameter = 94.*mm;
    this->siDetCrystalInnerDiameter = 16.*mm;
    this->siDetCrystalThickness = 6.15*mm;
    this->siDetRadialSegments = 10.;
    this->siDetPhiSegments = 12.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    this->siDetGuardRingOuterDiameter = 102*mm;
    this->siDetGuardRingInnerDiameter = 10*mm;
}

DetectionSystemSpice::~DetectionSystemSpice()
{    // LogicalVolumes in ConstructSPICEDetectionSystem
    delete [] siDetSpiceRing_log;

    delete siInnerGuardRing_log;
    delete siOuterGuardRing_log;

		delete [] siDetSpiceRing_SD;

}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build(G4SDManager* mySDman)
{
	
	G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;
   
	// Loop through each ring ...
	for(int ringID=0; ringID<10; ringID++) {
		
		// ... if the sensitive detector does not exist, create it and add to sd manager
		if( !siDetSpiceRing_SD[ringID] ) {
			G4String ringName = "/sd/allSpiceRing";
			ringName += ringID+1;
			G4String HCname = "CollectionSpiceRing";
			HCname += ringID+1;
    	siDetSpiceRing_SD[ringID] = new SensitiveDetector(ringName, HCname);
    	mySDman->AddNewDetector( siDetSpiceRing_SD[ringID] );
    } // end if( !siDetSpiceRing_SD[ringID] )
    
    // Build assembly volumes
    G4AssemblyVolume* myAssemblySiRing = new G4AssemblyVolume();    
  	this->assemblySiRing[ringID] = myAssemblySiRing;

		// Build Silicon Ring
  	BuildSiliconWafer(ringID+1);
  	
  	// Set Sensitive Detector
  	siDetSpiceRing_log[ringID]->SetSensitiveDetector( siDetSpiceRing_SD[ringID] );
  	
  	
  } // end for(int ringID)
  
  BuildInnerGuardRing();
  BuildOuterGuardRing();

  return 1;
} // end Build

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemSpice::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4int ringNumber, G4int nRadSeg, G4int detectorNumber)
{
  G4RotationMatrix* rotate = new G4RotationMatrix;
  G4int nRadSegTot = (G4int)this->siDetPhiSegments;
  G4double angle = (360./nRadSegTot)*(nRadSeg-0.5)*deg;
  rotate->rotateZ(angle);
  
  assemblySiRing[ringNumber-1]->MakeImprint(exp_hall_log, move, rotate, detectorNumber);

  return 1;
}

G4int DetectionSystemSpice::PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move)
{
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);
  assembly->MakeImprint(exp_hall_log, move, rotate, 0);

  return 1;
}

//---------------------------------------------------------//
// build functions for different parts                     //
// called in main build function                           //
//---------------------------------------------------------//
G4int DetectionSystemSpice::BuildSiliconWafer(G4int myRingID)
{
	// Define the material, return error if not found
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
  	G4cout << " ----> Material " << this->wafer_material 
  				 << " not found, cannot build the detector shell! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  // construct solid
  G4Tubs* siDetSpiceRingSec = BuildCrystal(myRingID);
  // construct logical volume
  if( !siDetSpiceRing_log[myRingID-1] )
	{
		G4String ringName = "siDetSpiceRing";
		ringName += myRingID;
		
		siDetSpiceRing_log[myRingID-1] = new G4LogicalVolume(siDetSpiceRingSec, material, ringName, 0, 0, 0);
		siDetSpiceRing_log[myRingID-1]->SetVisAttributes(vis_att);
	}
	this->assemblySiRing[myRingID-1]->AddPlacedVolume(siDetSpiceRing_log[myRingID-1], move, rotate);


  return 1;
}

G4int DetectionSystemSpice::BuildInnerGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the inner guard ring of Spice! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
					 this->siDetGuardRingInnerDiameter/2.,
					 this->siDetCrystalInnerDiameter/2.,
					 this->siDetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( siInnerGuardRing_log == NULL )
  {
    siInnerGuardRing_log = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
    siInnerGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(siInnerGuardRing_log, move, rotate);

  return 1;
}

G4int DetectionSystemSpice::BuildOuterGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the outer guard ring of Spice! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
					 this->siDetCrystalOuterDiameter/2.,
					 this->siDetGuardRingOuterDiameter/2.,
					 this->siDetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->siDetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( siOuterGuardRing_log == NULL )
  {
    siOuterGuardRing_log = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
    siOuterGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(siOuterGuardRing_log, move, rotate);

  return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemSpice::BuildCrystal(G4int myRingID)
{
  // define angle, length, thickness, and inner and outer diameter
  // of silicon detector segment
  G4double tube_element_length = (this->siDetCrystalOuterDiameter - this->siDetCrystalInnerDiameter)/(2*(this->siDetRadialSegments));
  G4double tube_element_angular_width = (360./this->siDetPhiSegments)*deg;
  G4double tube_element_outer_radius = ((G4double)this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(myRingID);
  G4double tube_element_inner_radius = (this->siDetCrystalInnerDiameter)/2.0 + tube_element_length*(myRingID-1);
  G4double tube_element_half_thickness = (this->siDetCrystalThickness)/2.0;

  // establish solid
  G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

  return crystal_block;
}//end ::BuildCrystal

