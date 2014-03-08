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
  siOuterGuardRing_log(0),
  siDetSpiceRing01_log(0),
  siDetSpiceRing02_log(0),
  siDetSpiceRing03_log(0),
  siDetSpiceRing04_log(0),
  siDetSpiceRing05_log(0),
  siDetSpiceRing06_log(0),
  siDetSpiceRing07_log(0),
  siDetSpiceRing08_log(0),
  siDetSpiceRing09_log(0),
  siDetSpiceRing10_log(0)
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
    delete siDetSpiceRing01_log;
    delete siDetSpiceRing02_log;
    delete siDetSpiceRing03_log;
    delete siDetSpiceRing04_log;
    delete siDetSpiceRing05_log;
    delete siDetSpiceRing06_log;
    delete siDetSpiceRing07_log;
    delete siDetSpiceRing08_log;
    delete siDetSpiceRing09_log;
    delete siDetSpiceRing10_log;
    delete siInnerGuardRing_log;
    delete siOuterGuardRing_log;

		delete siDetSpiceRing01_SD;
    delete siDetSpiceRing02_SD;
    delete siDetSpiceRing03_SD;
    delete siDetSpiceRing04_SD;
    delete siDetSpiceRing05_SD;
    delete siDetSpiceRing06_SD;
    delete siDetSpiceRing07_SD;
    delete siDetSpiceRing08_SD;
    delete siDetSpiceRing09_SD;
    delete siDetSpiceRing10_SD;
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemSpice::Build(G4SDManager* mySDman)
{
 	if( !siDetSpiceRing01_SD ) {
    siDetSpiceRing01_SD = new SensitiveDetector("/sd/allSpiceRing01", "CollectionSpiceRing01");
    mySDman->AddNewDetector( siDetSpiceRing01_SD );
  }

  if( !siDetSpiceRing02_SD ) {
    siDetSpiceRing02_SD = new SensitiveDetector("/sd/allSpiceRing02", "CollectionSpiceRing02");
    mySDman->AddNewDetector( siDetSpiceRing02_SD );
  }

  if( !siDetSpiceRing03_SD ) {
    siDetSpiceRing03_SD = new SensitiveDetector("/sd/allSpiceRing03", "CollectionSpiceRing03");
    mySDman->AddNewDetector( siDetSpiceRing03_SD );
  }

  if( !siDetSpiceRing04_SD ) {
    siDetSpiceRing04_SD = new SensitiveDetector("/sd/allSpiceRing04", "CollectionSpiceRing04");
    mySDman->AddNewDetector( siDetSpiceRing04_SD );
  }

  if( !siDetSpiceRing05_SD ) {
    siDetSpiceRing05_SD = new SensitiveDetector("/sd/allSpiceRing05", "CollectionSpiceRing05");
    mySDman->AddNewDetector( siDetSpiceRing05_SD );
  }

  if( !siDetSpiceRing06_SD ) {
    siDetSpiceRing06_SD = new SensitiveDetector("/sd/allSpiceRing06", "CollectionSpiceRing06");
    mySDman->AddNewDetector( siDetSpiceRing06_SD );
  }

  if( !siDetSpiceRing07_SD ) {
    siDetSpiceRing07_SD = new SensitiveDetector("/sd/allSpiceRing07", "CollectionSpiceRing07");
    mySDman->AddNewDetector( siDetSpiceRing07_SD );
  }

  if( !siDetSpiceRing08_SD ) {
    siDetSpiceRing08_SD = new SensitiveDetector("/sd/allSpiceRing08", "CollectionSpiceRing08");
    mySDman->AddNewDetector( siDetSpiceRing08_SD );
  }

  if( !siDetSpiceRing09_SD ) {
    siDetSpiceRing09_SD = new SensitiveDetector("/sd/allSpiceRing09", "CollectionSpiceRing09");
    mySDman->AddNewDetector( siDetSpiceRing09_SD );
  }

  if( !siDetSpiceRing10_SD ) {
    siDetSpiceRing10_SD = new SensitiveDetector("/sd/allSpiceRing10", "CollectionSpiceRing10");
    mySDman->AddNewDetector( siDetSpiceRing10_SD );
  }
  
  // Build assembly volumes
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;
  G4AssemblyVolume* myAssemblySiRing01 = new G4AssemblyVolume();    
  this->assemblySiRing01 = myAssemblySiRing01;
  G4AssemblyVolume* myAssemblySiRing02 = new G4AssemblyVolume();
  this->assemblySiRing02 = myAssemblySiRing02;
  G4AssemblyVolume* myAssemblySiRing03 = new G4AssemblyVolume();
  this->assemblySiRing03 = myAssemblySiRing03;
  G4AssemblyVolume* myAssemblySiRing04 = new G4AssemblyVolume();
  this->assemblySiRing04 = myAssemblySiRing04;
  G4AssemblyVolume* myAssemblySiRing05 = new G4AssemblyVolume();
  this->assemblySiRing05 = myAssemblySiRing05;
  G4AssemblyVolume* myAssemblySiRing06 = new G4AssemblyVolume();
  this->assemblySiRing06 = myAssemblySiRing06;
  G4AssemblyVolume* myAssemblySiRing07 = new G4AssemblyVolume();
  this->assemblySiRing07 = myAssemblySiRing07;
  G4AssemblyVolume* myAssemblySiRing08 = new G4AssemblyVolume();
  this->assemblySiRing08 = myAssemblySiRing08;
  G4AssemblyVolume* myAssemblySiRing09 = new G4AssemblyVolume();
  this->assemblySiRing09 = myAssemblySiRing09;
  G4AssemblyVolume* myAssemblySiRing10 = new G4AssemblyVolume();
  this->assemblySiRing10 = myAssemblySiRing10;

  G4cout << "BuildSiliconWafer" << G4endl;
  for(int ringID=0; ringID<10; ringID++)
    BuildSiliconWafer(ringID+1);
  BuildInnerGuardRing();
  BuildOuterGuardRing();
  
  // Sensitive Detector
  siDetSpiceRing01_log->SetSensitiveDetector( siDetSpiceRing01_SD );
  siDetSpiceRing02_log->SetSensitiveDetector( siDetSpiceRing02_SD );
  siDetSpiceRing03_log->SetSensitiveDetector( siDetSpiceRing03_SD );
  siDetSpiceRing04_log->SetSensitiveDetector( siDetSpiceRing04_SD );
  siDetSpiceRing05_log->SetSensitiveDetector( siDetSpiceRing05_SD );
  siDetSpiceRing06_log->SetSensitiveDetector( siDetSpiceRing06_SD );
  siDetSpiceRing07_log->SetSensitiveDetector( siDetSpiceRing07_SD );
  siDetSpiceRing08_log->SetSensitiveDetector( siDetSpiceRing08_SD );
  siDetSpiceRing09_log->SetSensitiveDetector( siDetSpiceRing09_SD );
  siDetSpiceRing10_log->SetSensitiveDetector( siDetSpiceRing10_SD );

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
  if(ringNumber == 1)
    assemblySiRing01->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 2)
    assemblySiRing02->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 3)
    assemblySiRing03->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 4)
    assemblySiRing04->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 5)
    assemblySiRing05->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 6)
    assemblySiRing06->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 7)
    assemblySiRing07->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 8)
    assemblySiRing08->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 9)
    assemblySiRing09->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 10)
    assemblySiRing10->MakeImprint(exp_hall_log, move, rotate, detectorNumber);

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
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the detector shell! " << G4endl;
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

  if(myRingID == 1)
    {
      // construct solid
      G4Tubs* siDetSpiceRing01Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing01_log == NULL )
		{
		  siDetSpiceRing01_log = new G4LogicalVolume(siDetSpiceRing01Sec, material, "siDetSpiceRing01", 0, 0, 0);
		  siDetSpiceRing01_log->SetVisAttributes(vis_att);
		}
      this->assemblySiRing01->AddPlacedVolume(siDetSpiceRing01_log, move, rotate);
    } // end if(myRingID == 1)
    
  if(myRingID == 2)
    {
      // construct solid
      G4Tubs* siDetSpiceRing02Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing02_log == NULL )
		{
		  siDetSpiceRing02_log = new G4LogicalVolume(siDetSpiceRing02Sec, material, "siDetSpiceRing02", 0, 0, 0);
		  siDetSpiceRing02_log->SetVisAttributes(vis_att);
		}
      this->assemblySiRing02->AddPlacedVolume(siDetSpiceRing02_log, move, rotate);
    } // end if(myRingID == 2)
  if(myRingID == 3)
    {
      // construct solid
      G4Tubs* siDetSpiceRing03Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing03_log == NULL )
		{
		  siDetSpiceRing03_log = new G4LogicalVolume(siDetSpiceRing03Sec, material, "siDetSpiceRing03", 0, 0, 0);
		  siDetSpiceRing03_log->SetVisAttributes(vis_att);
		}
      this->assemblySiRing03->AddPlacedVolume(siDetSpiceRing03_log, move, rotate);
    } // end if(myRingID == 3)
  if(myRingID == 4)
    {
      // construct solid
      G4Tubs* siDetSpiceRing04Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing04_log == NULL )
	{
	  siDetSpiceRing04_log = new G4LogicalVolume(siDetSpiceRing04Sec, material, "siDetSpiceRing04", 0, 0, 0);
	  siDetSpiceRing04_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing04->AddPlacedVolume(siDetSpiceRing04_log, move, rotate);
    } // end if(myRingID == 4)
  if(myRingID == 5)
    {
      // construct solid
      G4Tubs* siDetSpiceRing05Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing05_log == NULL )
	{
	  siDetSpiceRing05_log = new G4LogicalVolume(siDetSpiceRing05Sec, material, "siDetSpiceRing05", 0, 0, 0);
	  siDetSpiceRing05_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing05->AddPlacedVolume(siDetSpiceRing05_log, move, rotate);
    } // end if(myRingID == 5)
  if(myRingID == 6)
    {
      // construct solid
      G4Tubs* siDetSpiceRing06Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing06_log == NULL )
	{
	  siDetSpiceRing06_log = new G4LogicalVolume(siDetSpiceRing06Sec, material, "siDetSpiceRing06", 0, 0, 0);
	  siDetSpiceRing06_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing06->AddPlacedVolume(siDetSpiceRing06_log, move, rotate);
    } // end if(myRingID == 6)
  if(myRingID == 7)
    {
      // construct solid
      G4Tubs* siDetSpiceRing07Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing07_log == NULL )
	{
	  siDetSpiceRing07_log = new G4LogicalVolume(siDetSpiceRing07Sec, material, "siDetSpiceRing07", 0, 0, 0);
	  siDetSpiceRing07_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing07->AddPlacedVolume(siDetSpiceRing07_log, move, rotate);
    } // end if(myRingID == 7)
  if(myRingID == 8)
    {
      // construct solid
      G4Tubs* siDetSpiceRing08Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing08_log == NULL )
	{
	  siDetSpiceRing08_log = new G4LogicalVolume(siDetSpiceRing08Sec, material, "siDetSpiceRing08", 0, 0, 0);
	  siDetSpiceRing08_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing08->AddPlacedVolume(siDetSpiceRing08_log, move, rotate);
    } // end if(myRingID == 8)
  if(myRingID == 9)
    {
      // construct solid
      G4Tubs* siDetSpiceRing09Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing09_log == NULL )
	{
	  siDetSpiceRing09_log = new G4LogicalVolume(siDetSpiceRing09Sec, material, "siDetSpiceRing09", 0, 0, 0);
	  siDetSpiceRing09_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing09->AddPlacedVolume(siDetSpiceRing09_log, move, rotate);
    } // end if(myRingID == 9)
  if(myRingID == 10)
    {
      // construct solid
      G4Tubs* siDetSpiceRing10Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetSpiceRing10_log == NULL )
	{
	  siDetSpiceRing10_log = new G4LogicalVolume(siDetSpiceRing10Sec, material, "siDetSpiceRing10", 0, 0, 0);
	  siDetSpiceRing10_log->SetVisAttributes(vis_att);
	}
      this->assemblySiRing10->AddPlacedVolume(siDetSpiceRing10_log, move, rotate);
    } // end if(myRingID == 10)

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

