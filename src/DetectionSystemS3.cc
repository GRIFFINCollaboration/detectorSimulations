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

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemS3.hh"

DetectionSystemS3::DetectionSystemS3() :
  // LogicalVolumes
  S3InnerGuardRing_log(0),
  S3OuterGuardRing_log(0),
  siDetS3Ring01_log(0),
  siDetS3Ring02_log(0),
  siDetS3Ring03_log(0),
  siDetS3Ring04_log(0),
  siDetS3Ring05_log(0),
  siDetS3Ring06_log(0),
  siDetS3Ring07_log(0),
  siDetS3Ring08_log(0),
  siDetS3Ring09_log(0),
  siDetS3Ring10_log(0),
  siDetS3Ring11_log(0),
  siDetS3Ring12_log(0),
  siDetS3Ring13_log(0),
  siDetS3Ring14_log(0),
  siDetS3Ring15_log(0),
  siDetS3Ring16_log(0),
  siDetS3Ring17_log(0),
  siDetS3Ring18_log(0),
  siDetS3Ring19_log(0),
  siDetS3Ring20_log(0),
  siDetS3Ring21_log(0),
  siDetS3Ring22_log(0),
  siDetS3Ring23_log(0),
  siDetS3Ring24_log(0)
{
    /////////////////////////////////////////////////////////////////////
    // SPICE Physical Properties
    /////////////////////////////////////////////////////////////////////

    this->wafer_material             = "Silicon";

    //-----------------------------//
    // parameters for the annular  //
    // planar detector crystal     //
    //-----------------------------//
    this->S3DetCrystalOuterDiameter = 70.*mm;
    this->S3DetCrystalInnerDiameter = 22.*mm;
    this->S3DetCrystalThickness = .15*mm;
    this->S3DetRadialSegments = 24.;
    this->S3DetPhiSegments = 32.;

    //-----------------------------//
    // parameters for guard ring   //
    //-----------------------------//
    this->S3DetGuardRingOuterDiameter = 76*mm;
    this->S3DetGuardRingInnerDiameter = 20*mm;
}

DetectionSystemS3::~DetectionSystemS3()
{    // LogicalVolumes in ConstructSPICEDetectionSystem
    delete siDetS3Ring01_log;
    delete siDetS3Ring02_log;
    delete siDetS3Ring03_log;
    delete siDetS3Ring04_log;
    delete siDetS3Ring05_log;
    delete siDetS3Ring06_log;
    delete siDetS3Ring07_log;
    delete siDetS3Ring08_log;
    delete siDetS3Ring09_log;
    delete siDetS3Ring10_log;
    delete siDetS3Ring11_log;
    delete siDetS3Ring12_log;
    delete siDetS3Ring13_log;
    delete siDetS3Ring14_log;
    delete siDetS3Ring15_log;
    delete siDetS3Ring16_log;
    delete siDetS3Ring17_log;
    delete siDetS3Ring18_log;
    delete siDetS3Ring19_log;
    delete siDetS3Ring20_log;
    delete siDetS3Ring21_log;
    delete siDetS3Ring22_log;
    delete siDetS3Ring23_log;
    delete siDetS3Ring24_log;
    delete S3InnerGuardRing_log;
    delete S3OuterGuardRing_log;

    delete siDetS3Ring01_SD;
    delete siDetS3Ring02_SD;
    delete siDetS3Ring03_SD;
    delete siDetS3Ring04_SD;
    delete siDetS3Ring05_SD;
    delete siDetS3Ring06_SD;
    delete siDetS3Ring07_SD;
    delete siDetS3Ring08_SD;
    delete siDetS3Ring09_SD;
    delete siDetS3Ring10_SD;
    delete siDetS3Ring11_SD;
    delete siDetS3Ring12_SD;
    delete siDetS3Ring13_SD;
    delete siDetS3Ring14_SD;
    delete siDetS3Ring15_SD;
    delete siDetS3Ring16_SD;
    delete siDetS3Ring17_SD;
    delete siDetS3Ring18_SD;
    delete siDetS3Ring19_SD;
    delete siDetS3Ring20_SD;
    delete siDetS3Ring21_SD;
    delete siDetS3Ring22_SD;
    delete siDetS3Ring23_SD;
    delete siDetS3Ring24_SD;
}

//---------------------------------------------------------//
// main build function called in DetectorConstruction      //
// when detector is constructed                            //
//---------------------------------------------------------//
G4int DetectionSystemS3::Build(G4SDManager* mySDman)
{
  if( !siDetS3Ring01_SD ) {
    siDetS3Ring01_SD = new SensitiveDetector("/sd/allS3Ring01", "CollectionS3Ring01");
    mySDman->AddNewDetector( siDetS3Ring01_SD );
  }

  if( !siDetS3Ring02_SD ) {
    siDetS3Ring02_SD = new SensitiveDetector("/sd/allS3Ring02", "CollectionS3Ring02");
    mySDman->AddNewDetector( siDetS3Ring02_SD );
  }

  if( !siDetS3Ring03_SD ) {
    siDetS3Ring03_SD = new SensitiveDetector("/sd/allS3Ring03", "CollectionS3Ring03");
    mySDman->AddNewDetector( siDetS3Ring03_SD );
  }

  if( !siDetS3Ring04_SD ) {
    siDetS3Ring04_SD = new SensitiveDetector("/sd/allS3Ring04", "CollectionS3Ring04");
    mySDman->AddNewDetector( siDetS3Ring04_SD );
  }

  if( !siDetS3Ring05_SD ) {
    siDetS3Ring05_SD = new SensitiveDetector("/sd/allS3Ring05", "CollectionS3Ring05");
    mySDman->AddNewDetector( siDetS3Ring05_SD );
  }

  if( !siDetS3Ring06_SD ) {
    siDetS3Ring06_SD = new SensitiveDetector("/sd/allS3Ring06", "CollectionS3Ring06");
    mySDman->AddNewDetector( siDetS3Ring06_SD );
  }

  if( !siDetS3Ring07_SD ) {
    siDetS3Ring07_SD = new SensitiveDetector("/sd/allS3Ring07", "CollectionS3Ring07");
    mySDman->AddNewDetector( siDetS3Ring07_SD );
  }

  if( !siDetS3Ring08_SD ) {
    siDetS3Ring08_SD = new SensitiveDetector("/sd/allS3Ring08", "CollectionS3Ring08");
    mySDman->AddNewDetector( siDetS3Ring08_SD );
  }

  if( !siDetS3Ring09_SD ) {
    siDetS3Ring09_SD = new SensitiveDetector("/sd/allS3Ring09", "CollectionS3Ring09");
    mySDman->AddNewDetector( siDetS3Ring09_SD );
  }

  if( !siDetS3Ring10_SD ) {
    siDetS3Ring10_SD = new SensitiveDetector("/sd/allS3Ring10", "CollectionS3Ring10");
    mySDman->AddNewDetector( siDetS3Ring10_SD );
  }
  if( !siDetS3Ring11_SD ) {
    siDetS3Ring11_SD = new SensitiveDetector("/sd/allS3Ring11", "CollectionS3Ring11");
    mySDman->AddNewDetector( siDetS3Ring11_SD );
  }

  if( !siDetS3Ring12_SD ) {
    siDetS3Ring12_SD = new SensitiveDetector("/sd/allS3Ring12", "CollectionS3Ring12");
    mySDman->AddNewDetector( siDetS3Ring12_SD );
  }

  if( !siDetS3Ring13_SD ) {
    siDetS3Ring13_SD = new SensitiveDetector("/sd/allS3Ring13", "CollectionS3Ring13");
    mySDman->AddNewDetector( siDetS3Ring13_SD );
  }

  if( !siDetS3Ring14_SD ) {
    siDetS3Ring14_SD = new SensitiveDetector("/sd/allS3Ring14", "CollectionS3Ring14");
    mySDman->AddNewDetector( siDetS3Ring14_SD );
  }

  if( !siDetS3Ring15_SD ) {
    siDetS3Ring15_SD = new SensitiveDetector("/sd/allS3Ring15", "CollectionS3Ring15");
    mySDman->AddNewDetector( siDetS3Ring15_SD );
  }

  if( !siDetS3Ring16_SD ) {
    siDetS3Ring16_SD = new SensitiveDetector("/sd/allS3Ring16", "CollectionS3Ring16");
    mySDman->AddNewDetector( siDetS3Ring16_SD );
  }

  if( !siDetS3Ring17_SD ) {
    siDetS3Ring17_SD = new SensitiveDetector("/sd/allS3Ring17", "CollectionS3Ring17");
    mySDman->AddNewDetector( siDetS3Ring17_SD );
  }

  if( !siDetS3Ring18_SD ) {
    siDetS3Ring18_SD = new SensitiveDetector("/sd/allS3Ring18", "CollectionS3Ring18");
    mySDman->AddNewDetector( siDetS3Ring18_SD );
  }

  if( !siDetS3Ring19_SD ) {
    siDetS3Ring19_SD = new SensitiveDetector("/sd/allS3Ring19", "CollectionS3Ring19");
    mySDman->AddNewDetector( siDetS3Ring19_SD );
  }

  if( !siDetS3Ring20_SD ) {
    siDetS3Ring20_SD = new SensitiveDetector("/sd/allS3Ring20", "CollectionS3Ring20");
    mySDman->AddNewDetector( siDetS3Ring20_SD );
  }

  if( !siDetS3Ring21_SD ) {
    siDetS3Ring21_SD = new SensitiveDetector("/sd/allS3Ring21", "CollectionS3Ring21");
    mySDman->AddNewDetector( siDetS3Ring21_SD );
  }

  if( !siDetS3Ring22_SD ) {
    siDetS3Ring22_SD = new SensitiveDetector("/sd/allS3Ring22", "CollectionS3Ring22");
    mySDman->AddNewDetector( siDetS3Ring22_SD );
  }

  if( !siDetS3Ring23_SD ) {
    siDetS3Ring23_SD = new SensitiveDetector("/sd/allS3Ring23", "CollectionS3Ring23");
    mySDman->AddNewDetector( siDetS3Ring23_SD );
  }

  if( !siDetS3Ring24_SD ) {
    siDetS3Ring24_SD = new SensitiveDetector("/sd/allS3Ring24", "CollectionS3Ring24");
    mySDman->AddNewDetector( siDetS3Ring24_SD );
  }

  // Build assembly volumes
  G4AssemblyVolume* myAssembly = new G4AssemblyVolume();
  this->assembly = myAssembly;
  G4AssemblyVolume* myAssemblySiRing01 = new G4AssemblyVolume();
  this->assemblyS3Ring01 = myAssemblySiRing01;
  G4AssemblyVolume* myAssemblySiRing02 = new G4AssemblyVolume();
  this->assemblyS3Ring02 = myAssemblySiRing02;
  G4AssemblyVolume* myAssemblySiRing03 = new G4AssemblyVolume();
  this->assemblyS3Ring03 = myAssemblySiRing03;
  G4AssemblyVolume* myAssemblySiRing04 = new G4AssemblyVolume();
  this->assemblyS3Ring04 = myAssemblySiRing04;
  G4AssemblyVolume* myAssemblySiRing05 = new G4AssemblyVolume();
  this->assemblyS3Ring05 = myAssemblySiRing05;
  G4AssemblyVolume* myAssemblySiRing06 = new G4AssemblyVolume();
  this->assemblyS3Ring06 = myAssemblySiRing06;
  G4AssemblyVolume* myAssemblySiRing07 = new G4AssemblyVolume();
  this->assemblyS3Ring07 = myAssemblySiRing07;
  G4AssemblyVolume* myAssemblySiRing08 = new G4AssemblyVolume();
  this->assemblyS3Ring08 = myAssemblySiRing08;
  G4AssemblyVolume* myAssemblySiRing09 = new G4AssemblyVolume();
  this->assemblyS3Ring09 = myAssemblySiRing09;
  G4AssemblyVolume* myAssemblySiRing10 = new G4AssemblyVolume();
  this->assemblyS3Ring10 = myAssemblySiRing10;
  G4AssemblyVolume* myAssemblySiRing11 = new G4AssemblyVolume();
  this->assemblyS3Ring11 = myAssemblySiRing11;
  G4AssemblyVolume* myAssemblySiRing12 = new G4AssemblyVolume();
  this->assemblyS3Ring12 = myAssemblySiRing12;
  G4AssemblyVolume* myAssemblySiRing13 = new G4AssemblyVolume();
  this->assemblyS3Ring13 = myAssemblySiRing13;
  G4AssemblyVolume* myAssemblySiRing14 = new G4AssemblyVolume();
  this->assemblyS3Ring14 = myAssemblySiRing14;
  G4AssemblyVolume* myAssemblySiRing15 = new G4AssemblyVolume();
  this->assemblyS3Ring15 = myAssemblySiRing15;
  G4AssemblyVolume* myAssemblySiRing16 = new G4AssemblyVolume();
  this->assemblyS3Ring16 = myAssemblySiRing16;
  G4AssemblyVolume* myAssemblySiRing17 = new G4AssemblyVolume();
  this->assemblyS3Ring17 = myAssemblySiRing17;
  G4AssemblyVolume* myAssemblySiRing18 = new G4AssemblyVolume();
  this->assemblyS3Ring18 = myAssemblySiRing18;
  G4AssemblyVolume* myAssemblySiRing19 = new G4AssemblyVolume();
  this->assemblyS3Ring19 = myAssemblySiRing19;
  G4AssemblyVolume* myAssemblySiRing20 = new G4AssemblyVolume();
  this->assemblyS3Ring20 = myAssemblySiRing20;
  G4AssemblyVolume* myAssemblySiRing21 = new G4AssemblyVolume();
  this->assemblyS3Ring21 = myAssemblySiRing21;
  G4AssemblyVolume* myAssemblySiRing22 = new G4AssemblyVolume();
  this->assemblyS3Ring22 = myAssemblySiRing22;
  G4AssemblyVolume* myAssemblySiRing23 = new G4AssemblyVolume();
  this->assemblyS3Ring23 = myAssemblySiRing23;
  G4AssemblyVolume* myAssemblySiRing24 = new G4AssemblyVolume();
  this->assemblyS3Ring24 = myAssemblySiRing24;

  G4cout << "BuildSiliconWafer" << G4endl;
  for(int ringID=0; ringID<24; ringID++)
    BuildSiliconWafer(ringID+1);
  BuildInnerGuardRing();
  BuildOuterGuardRing();

  // Sensitive Detector
  siDetS3Ring01_log->SetSensitiveDetector( siDetS3Ring01_SD );
  siDetS3Ring02_log->SetSensitiveDetector( siDetS3Ring02_SD );
  siDetS3Ring03_log->SetSensitiveDetector( siDetS3Ring03_SD );
  siDetS3Ring04_log->SetSensitiveDetector( siDetS3Ring04_SD );
  siDetS3Ring05_log->SetSensitiveDetector( siDetS3Ring05_SD );
  siDetS3Ring06_log->SetSensitiveDetector( siDetS3Ring06_SD );
  siDetS3Ring07_log->SetSensitiveDetector( siDetS3Ring07_SD );
  siDetS3Ring08_log->SetSensitiveDetector( siDetS3Ring08_SD );
  siDetS3Ring09_log->SetSensitiveDetector( siDetS3Ring09_SD );
  siDetS3Ring10_log->SetSensitiveDetector( siDetS3Ring10_SD );
  siDetS3Ring11_log->SetSensitiveDetector( siDetS3Ring11_SD );
  siDetS3Ring12_log->SetSensitiveDetector( siDetS3Ring12_SD );
  siDetS3Ring13_log->SetSensitiveDetector( siDetS3Ring13_SD );
  siDetS3Ring14_log->SetSensitiveDetector( siDetS3Ring14_SD );
  siDetS3Ring15_log->SetSensitiveDetector( siDetS3Ring15_SD );
  siDetS3Ring16_log->SetSensitiveDetector( siDetS3Ring16_SD );
  siDetS3Ring17_log->SetSensitiveDetector( siDetS3Ring17_SD );
  siDetS3Ring18_log->SetSensitiveDetector( siDetS3Ring18_SD );
  siDetS3Ring19_log->SetSensitiveDetector( siDetS3Ring19_SD );
  siDetS3Ring20_log->SetSensitiveDetector( siDetS3Ring20_SD );
  siDetS3Ring21_log->SetSensitiveDetector( siDetS3Ring21_SD );
  siDetS3Ring22_log->SetSensitiveDetector( siDetS3Ring22_SD );
  siDetS3Ring23_log->SetSensitiveDetector( siDetS3Ring23_SD );
  siDetS3Ring24_log->SetSensitiveDetector( siDetS3Ring24_SD );

  return 1;
}

//---------------------------------------------------------//
// "place" function called in DetectorMessenger            //
// if detector is added                                    //
//---------------------------------------------------------//
G4int DetectionSystemS3::PlaceDetector(G4LogicalVolume* exp_hall_log, G4ThreeVector move, G4int ringNumber, G4int nRadSeg, G4int detectorNumber)
{
  G4RotationMatrix* rotate = new G4RotationMatrix;
  G4int nRadSegTot = (G4int)this->S3DetPhiSegments;
  G4double angle = (360./nRadSegTot)*(nRadSeg-0.5)*deg;
  rotate->rotateZ(angle);
  if(ringNumber == 1)
    assemblyS3Ring01->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 2)
    assemblyS3Ring02->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 3)
    assemblyS3Ring03->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 4)
    assemblyS3Ring04->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 5)
    assemblyS3Ring05->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 6)
    assemblyS3Ring06->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 7)
    assemblyS3Ring07->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 8)
    assemblyS3Ring08->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 9)
    assemblyS3Ring09->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 10)
    assemblyS3Ring10->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 11)
    assemblyS3Ring11->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 12)
    assemblyS3Ring12->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 13)
    assemblyS3Ring13->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 14)
    assemblyS3Ring14->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 15)
    assemblyS3Ring15->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 16)
    assemblyS3Ring16->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 17)
    assemblyS3Ring17->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 18)
    assemblyS3Ring18->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 19)
    assemblyS3Ring19->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 20)
    assemblyS3Ring20->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 21)
    assemblyS3Ring21->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 22)
    assemblyS3Ring22->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 23)
    assemblyS3Ring23->MakeImprint(exp_hall_log, move, rotate, detectorNumber);
  if(ringNumber == 24)
    assemblyS3Ring24->MakeImprint(exp_hall_log, move, rotate, detectorNumber);

  return 1;
}

G4int DetectionSystemS3::PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4ThreeVector move)
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
G4int DetectionSystemS3::BuildSiliconWafer(G4int myRingID)
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
  G4double z_position		= -(this->S3DetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  if(myRingID == 1)
    {
      // construct solid
      G4Tubs* siDetS3Ring01Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring01_log == NULL )
	{
	  siDetS3Ring01_log = new G4LogicalVolume(siDetS3Ring01Sec, material, "siDetS3Ring01", 0, 0, 0);
	  siDetS3Ring01_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring01->AddPlacedVolume(siDetS3Ring01_log, move, rotate);
    } // end if(myRingID == 1)
  if(myRingID == 2)
    {
      // construct solid
      G4Tubs* siDetS3Ring02Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring02_log == NULL )
	{
	  siDetS3Ring02_log = new G4LogicalVolume(siDetS3Ring02Sec, material, "siDetS3Ring02", 0, 0, 0);
	  siDetS3Ring02_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring02->AddPlacedVolume(siDetS3Ring02_log, move, rotate);
    } // end if(myRingID == 2)
  if(myRingID == 3)
    {
      // construct solid
      G4Tubs* siDetS3Ring03Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring03_log == NULL )
	{
	  siDetS3Ring03_log = new G4LogicalVolume(siDetS3Ring03Sec, material, "siDetS3Ring03", 0, 0, 0);
	  siDetS3Ring03_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring03->AddPlacedVolume(siDetS3Ring03_log, move, rotate);
    } // end if(myRingID == 3)
  if(myRingID == 4)
    {
      // construct solid
      G4Tubs* siDetS3Ring04Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring04_log == NULL )
	{
	  siDetS3Ring04_log = new G4LogicalVolume(siDetS3Ring04Sec, material, "siDetS3Ring04", 0, 0, 0);
	  siDetS3Ring04_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring04->AddPlacedVolume(siDetS3Ring04_log, move, rotate);
    } // end if(myRingID == 4)
  if(myRingID == 5)
    {
      // construct solid
      G4Tubs* siDetS3Ring05Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring05_log == NULL )
	{
	  siDetS3Ring05_log = new G4LogicalVolume(siDetS3Ring05Sec, material, "siDetS3Ring05", 0, 0, 0);
	  siDetS3Ring05_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring05->AddPlacedVolume(siDetS3Ring05_log, move, rotate);
    } // end if(myRingID == 5)
  if(myRingID == 6)
    {
      // construct solid
      G4Tubs* siDetS3Ring06Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring06_log == NULL )
	{
	  siDetS3Ring06_log = new G4LogicalVolume(siDetS3Ring06Sec, material, "siDetS3Ring06", 0, 0, 0);
	  siDetS3Ring06_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring06->AddPlacedVolume(siDetS3Ring06_log, move, rotate);
    } // end if(myRingID == 6)
  if(myRingID == 7)
    {
      // construct solid
      G4Tubs* siDetS3Ring07Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring07_log == NULL )
	{
	  siDetS3Ring07_log = new G4LogicalVolume(siDetS3Ring07Sec, material, "siDetS3Ring07", 0, 0, 0);
	  siDetS3Ring07_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring07->AddPlacedVolume(siDetS3Ring07_log, move, rotate);
    } // end if(myRingID == 7)
  if(myRingID == 8)
    {
      // construct solid
      G4Tubs* siDetS3Ring08Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring08_log == NULL )
	{
	  siDetS3Ring08_log = new G4LogicalVolume(siDetS3Ring08Sec, material, "siDetS3Ring08", 0, 0, 0);
	  siDetS3Ring08_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring08->AddPlacedVolume(siDetS3Ring08_log, move, rotate);
    } // end if(myRingID == 8)
  if(myRingID == 9)
    {
      // construct solid
      G4Tubs* siDetS3Ring09Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring09_log == NULL )
	{
	  siDetS3Ring09_log = new G4LogicalVolume(siDetS3Ring09Sec, material, "siDetS3Ring09", 0, 0, 0);
	  siDetS3Ring09_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring09->AddPlacedVolume(siDetS3Ring09_log, move, rotate);
    } // end if(myRingID == 9)
  if(myRingID == 10)
    {
      // construct solid
      G4Tubs* siDetS3Ring10Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring10_log == NULL )
	{
	  siDetS3Ring10_log = new G4LogicalVolume(siDetS3Ring10Sec, material, "siDetS3Ring10", 0, 0, 0);
	  siDetS3Ring10_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring10->AddPlacedVolume(siDetS3Ring10_log, move, rotate);
    } // end if(myRingID == 10)
  if(myRingID == 11)
    {
      // construct solid
      G4Tubs* siDetS3Ring11Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring11_log == NULL )
	{
	  siDetS3Ring11_log = new G4LogicalVolume(siDetS3Ring11Sec, material, "siDetS3Ring11", 0, 0, 0);
	  siDetS3Ring11_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring11->AddPlacedVolume(siDetS3Ring11_log, move, rotate);
    } // end if(myRingID == 11)
  if(myRingID == 12)
    {
      // construct solid
      G4Tubs* siDetS3Ring12Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring12_log == NULL )
	{
	  siDetS3Ring12_log = new G4LogicalVolume(siDetS3Ring12Sec, material, "siDetS3Ring12", 0, 0, 0);
	  siDetS3Ring12_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring12->AddPlacedVolume(siDetS3Ring12_log, move, rotate);
    } // end if(myRingID == 12)
  if(myRingID == 13)
    {
      // construct solid
      G4Tubs* siDetS3Ring13Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring13_log == NULL )
	{
	  siDetS3Ring13_log = new G4LogicalVolume(siDetS3Ring13Sec, material, "siDetS3Ring13", 0, 0, 0);
	  siDetS3Ring13_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring13->AddPlacedVolume(siDetS3Ring13_log, move, rotate);
    } // end if(myRingID == 13)
  if(myRingID == 14)
    {
      // construct solid
      G4Tubs* siDetS3Ring14Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring14_log == NULL )
	{
	  siDetS3Ring14_log = new G4LogicalVolume(siDetS3Ring14Sec, material, "siDetS3Ring14", 0, 0, 0);
	  siDetS3Ring14_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring14->AddPlacedVolume(siDetS3Ring14_log, move, rotate);
    } // end if(myRingID == 14)
  if(myRingID == 15)
    {
      // construct solid
      G4Tubs* siDetS3Ring15Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring15_log == NULL )
	{
	  siDetS3Ring15_log = new G4LogicalVolume(siDetS3Ring15Sec, material, "siDetS3Ring15", 0, 0, 0);
	  siDetS3Ring15_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring15->AddPlacedVolume(siDetS3Ring15_log, move, rotate);
    } // end if(myRingID == 15)
  if(myRingID == 16)
    {
      // construct solid
      G4Tubs* siDetS3Ring16Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring16_log == NULL )
	{
	  siDetS3Ring16_log = new G4LogicalVolume(siDetS3Ring16Sec, material, "siDetS3Ring16", 0, 0, 0);
	  siDetS3Ring16_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring16->AddPlacedVolume(siDetS3Ring16_log, move, rotate);
    } // end if(myRingID == 16)
  if(myRingID == 17)
    {
      // construct solid
      G4Tubs* siDetS3Ring17Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring17_log == NULL )
	{
	  siDetS3Ring17_log = new G4LogicalVolume(siDetS3Ring17Sec, material, "siDetS3Ring17", 0, 0, 0);
	  siDetS3Ring17_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring17->AddPlacedVolume(siDetS3Ring17_log, move, rotate);
    } // end if(myRingID == 17)
  if(myRingID == 18)
    {
      // construct solid
      G4Tubs* siDetS3Ring18Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring18_log == NULL )
	{
	  siDetS3Ring18_log = new G4LogicalVolume(siDetS3Ring18Sec, material, "siDetS3Ring18", 0, 0, 0);
	  siDetS3Ring18_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring18->AddPlacedVolume(siDetS3Ring18_log, move, rotate);
    } // end if(myRingID == 18)
  if(myRingID == 19)
    {
      // construct solid
      G4Tubs* siDetS3Ring19Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring19_log == NULL )
	{
	  siDetS3Ring19_log = new G4LogicalVolume(siDetS3Ring19Sec, material, "siDetS3Ring19", 0, 0, 0);
	  siDetS3Ring19_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring19->AddPlacedVolume(siDetS3Ring19_log, move, rotate);
    } // end if(myRingID == 19)
  if(myRingID == 20)
    {
      // construct solid
      G4Tubs* siDetS3Ring20Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring20_log == NULL )
	{
	  siDetS3Ring20_log = new G4LogicalVolume(siDetS3Ring20Sec, material, "siDetS3Ring20", 0, 0, 0);
	  siDetS3Ring20_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring20->AddPlacedVolume(siDetS3Ring20_log, move, rotate);
    } // end if(myRingID == 20)
  if(myRingID == 21)
    {
      // construct solid
      G4Tubs* siDetS3Ring21Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring21_log == NULL )
	{
	  siDetS3Ring21_log = new G4LogicalVolume(siDetS3Ring21Sec, material, "siDetS3Ring21", 0, 0, 0);
	  siDetS3Ring21_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring21->AddPlacedVolume(siDetS3Ring21_log, move, rotate);
    } // end if(myRingID == 21)
  if(myRingID == 22)
    {
      // construct solid
      G4Tubs* siDetS3Ring22Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring22_log == NULL )
	{
	  siDetS3Ring22_log = new G4LogicalVolume(siDetS3Ring22Sec, material, "siDetS3Ring22", 0, 0, 0);
	  siDetS3Ring22_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring22->AddPlacedVolume(siDetS3Ring22_log, move, rotate);
    } // end if(myRingID == 22)
  if(myRingID == 23)
    {
      // construct solid
      G4Tubs* siDetS3Ring23Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring23_log == NULL )
	{
	  siDetS3Ring23_log = new G4LogicalVolume(siDetS3Ring23Sec, material, "siDetS3Ring23", 0, 0, 0);
	  siDetS3Ring23_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring23->AddPlacedVolume(siDetS3Ring23_log, move, rotate);
    } // end if(myRingID == 23)
  if(myRingID == 24)
    {
      // construct solid
      G4Tubs* siDetS3Ring24Sec = BuildCrystal(myRingID);
      // construct logical volume
      if( siDetS3Ring24_log == NULL )
	{
	  siDetS3Ring24_log = new G4LogicalVolume(siDetS3Ring24Sec, material, "siDetS3Ring24", 0, 0, 0);
	  siDetS3Ring24_log->SetVisAttributes(vis_att);
	}
      this->assemblyS3Ring24->AddPlacedVolume(siDetS3Ring24_log, move, rotate);
    } // end if(myRingID == 24)

  return 1;
}

G4int DetectionSystemS3::BuildInnerGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the inner guard ring of the S3 detector! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* innerGuardRing = new G4Tubs("innerGuardRing",
					 this->S3DetGuardRingInnerDiameter/2.,
					 this->S3DetCrystalInnerDiameter/2.,
					 this->S3DetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->S3DetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( S3InnerGuardRing_log == NULL )
  {
    S3InnerGuardRing_log = new G4LogicalVolume(innerGuardRing, material, "innerGuardRing", 0,0,0);
    S3InnerGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(S3InnerGuardRing_log, move, rotate);

  return 1;
}

G4int DetectionSystemS3::BuildOuterGuardRing()
{
  G4Material* material = G4Material::GetMaterial(this->wafer_material);
  if( !material ) {
    G4cout << " ----> Material " << this->wafer_material << " not found, cannot build the outer guard ring of the S3 detector! " << G4endl;
    return 0;
  }

  // Set visualization attributes
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  vis_att->SetVisibility(true);

  G4Tubs* outerGuardRing = new G4Tubs("outerGuardRing",
					 this->S3DetCrystalOuterDiameter/2.,
					 this->S3DetGuardRingOuterDiameter/2.,
					 this->S3DetCrystalThickness/2.,0,360);

  // Define rotation and movement objects
  G4ThreeVector direction 	= G4ThreeVector(0,0,1);
  G4double z_position		= -(this->S3DetCrystalThickness/2.);
  G4ThreeVector move 		= z_position * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateZ(0*deg);

  //logical volume
  if( S3OuterGuardRing_log == NULL )
  {
    S3OuterGuardRing_log = new G4LogicalVolume(outerGuardRing, material, "outerGuardRing", 0,0,0);
    S3OuterGuardRing_log->SetVisAttributes(vis_att);
  }

  this->assembly->AddPlacedVolume(S3OuterGuardRing_log, move, rotate);

  return 1;
}

///////////////////////////////////////////////////////////////////////
// Methods used to build shapes
///////////////////////////////////////////////////////////////////////
G4Tubs* DetectionSystemS3::BuildCrystal(G4int myRingID)
{
  // define angle, length, thickness, and inner and outer diameter
  // of silicon detector segment
  G4double tube_element_length = (this->S3DetCrystalOuterDiameter - this->S3DetCrystalInnerDiameter)/(2*(this->S3DetRadialSegments));
  G4double tube_element_angular_width = (360./this->S3DetPhiSegments)*deg;
  G4double tube_element_outer_radius = ((G4double)this->S3DetCrystalInnerDiameter)/2.0 + tube_element_length*(myRingID);
  G4double tube_element_inner_radius = (this->S3DetCrystalInnerDiameter)/2.0 + tube_element_length*(myRingID-1);
  G4double tube_element_half_thickness = (this->S3DetCrystalThickness)/2.0;

  // establish solid
  G4Tubs* crystal_block = new G4Tubs("crystal_block",tube_element_inner_radius,tube_element_outer_radius,tube_element_half_thickness,0,tube_element_angular_width);

  return crystal_block;
}//end ::BuildCrystal

