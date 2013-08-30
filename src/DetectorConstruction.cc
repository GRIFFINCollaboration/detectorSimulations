//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SensitiveDetector.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

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

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//#include "G4FieldManager.hh"
//#include "G4UniformMagField.hh"
//#include "MagneticField.hh"
//#include "G4TransportationManager.hh"
//#include "Field.hh"
//#include "GlobalField.hh"

#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

#include "DetectionSystemGammaTracking.hh"
#include "DetectionSystemBrillance380V1.hh"
#include "DetectionSystem8pi.hh"
#include "DetectionSystemGriffin.hh"
#include "DetectionSystemSceptar.hh"
#include "DetectionSystemSpice.hh"
#include "DetectionSystemSpiceV02.hh"
#include "DetectionSystemPaces.hh"
#include "DetectionSystemSodiumIodide.hh"

#include "ApparatusGenericTarget.hh"
#include "ApparatusSpiceTargetChamber.hh"
#include "Apparatus8piVacuumChamber.hh"
#include "Apparatus8piVacuumChamberDelrinShell.hh"
#include "ApparatusFieldBox.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
    // Fields
    expHallMagField(0)
{
  // materials
  DefineMaterials();
  this->builtDetectors = false;
  this->matWorldName = "G4_AIR";

  this->hall_x = 10.0*m;
  this->hall_y = 10.0*m;
  this->hall_z = 10.0*m;

  this->hall_vis = false;

  // Generic Target Apparatus
  this->setGenericTargetMaterial   = false;
  this->setGenericTargetDimensions = false;
  this->setGenericTargetPosition   = false;

  // Field Box
  this->setFieldBoxMaterial= false;
  this->setFieldBoxDimensions= false;
  this->setFieldBoxPosition= false;
  this->setFieldBoxMagneticField= false;

  // parameters to suppress:
  DefineSuppressedParameters();

  // create commands for interactive definition
  detectorMessenger = new DetectorMessenger(this);

  // ensure the global field is initialized
  //(void)GlobalField::getObject();

  //expHallMagField = new MagneticField(); // Global field is set to zero
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Experimental hall (world volume)
  
  // search the world material by its name
  G4Material* matWorld = G4Material::GetMaterial(matWorldName);
  
  if( !matWorld ) {
    G4cout << " ----> Material " << matWorldName << " not found, cannot build the experimental hall! " << G4endl;
    return 0;
  }
  
  G4Box*             hallBox  = new G4Box("hallBox", this->hall_x, this->hall_y, this->hall_z);
  G4LogicalVolume*   hallLog  = new G4LogicalVolume(hallBox, matWorld, "hallLog", 0, 0, 0);
  G4VPhysicalVolume* hallPhys = new G4PVPlacement(0, G4ThreeVector(), "hallPhys", hallLog, 0, false, 0);
  
  G4VisAttributes* hallVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  hallVisAtt->SetForceWireframe(true);
  hallVisAtt->SetVisibility(this->hall_vis);
  hallLog->SetVisAttributes(hallVisAtt);
  this->hallLog = hallLog;  

  if(this->builtDetectors) {
      G4cout << "Already Built Detectors!" << G4endl;
  }
  else {
      // sensitive detector manager
      G4SDManager* mySDman = G4SDManager::GetSDMpointer();

      //  DetectionSystemGammaTracking* pGammaTracking = new DetectionSystemGammaTracking();
      //  this->myGammaTracking = pGammaTracking;
      //  this->myGammaTracking->Build(mySDman);

      DetectionSystemBrillance380V1* pBrillance380V1 = new DetectionSystemBrillance380V1();
      this->myBrillance380V1 = pBrillance380V1;
      this->myBrillance380V1->Build(mySDman);

      DetectionSystemSodiumIodide* pSodiumIodide = new DetectionSystemSodiumIodide();
      this->mySodiumIodide = pSodiumIodide;
      this->mySodiumIodide->Build(mySDman);

      DetectionSystemGriffin* pGriffinForward = new DetectionSystemGriffin(0); // Select Forward (0) or Back (1)
      this->myGriffinForward = pGriffinForward;
      this->myGriffinForward->Build(mySDman);

      DetectionSystemGriffin* pGriffinBack = new DetectionSystemGriffin(1); // Select Forward (0) or Back (1)
      this->myGriffinBack = pGriffinBack;
      this->myGriffinBack->Build(mySDman);

      DetectionSystem8pi* p8pi = new DetectionSystem8pi();
      this->my8pi = p8pi;
      this->my8pi->Build(mySDman);

      DetectionSystemSceptar* pSceptar = new DetectionSystemSceptar();
      this->mySceptar = pSceptar;
      this->mySceptar->Build(mySDman);

      DetectionSystemSpice* pSpice = new DetectionSystemSpice();
      this->mySpice = pSpice;
      this->mySpice->Build(mySDman);

      DetectionSystemSpiceV02* pSpiceV02 = new DetectionSystemSpiceV02();
      this->mySpiceV02 = pSpiceV02;
      this->mySpiceV02->Build(mySDman);

      DetectionSystemPaces* pPaces = new DetectionSystemPaces();
      this->myPaces = pPaces;
      this->myPaces->Build(mySDman);

      G4cout << "Built Detectors!" << G4endl;
  }

  this->builtDetectors = true;
  return hallPhys;
}

void DetectorConstruction::SetWorldMaterial( G4String name )
{
  this->matWorldName = name;
  UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldDimensions( G4ThreeVector vec )
{
  this->hall_x = vec.x();
  this->hall_y = vec.y();
  this->hall_z = vec.z();
  UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldVis( G4bool vis )
{
  this->hall_vis = vis;
  UpdateGeometry(); // auto update
}

void DetectorConstruction::SetWorldMagneticField( G4ThreeVector vec )
{
    //expHallMagField->SetFieldValue(G4ThreeVector(vec.x(),vec.y(),vec.z()));
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void DetectorConstruction::SetGenericTargetMaterial( G4String name )
{
  this->setGenericTargetMaterial = true;
  this->genericTargetMaterial = name;
  SetGenericTarget();
}

void DetectorConstruction::SetGenericTargetDimensions( G4ThreeVector vec )
{
  this->setGenericTargetDimensions = true;
  this->genericTargetDimensions = vec;
  SetGenericTarget();
}

void DetectorConstruction::SetGenericTargetPosition( G4ThreeVector vec )
{
  this->setGenericTargetPosition = true;
  this->genericTargetPosition = vec;
  SetGenericTarget();
}

void DetectorConstruction::SetGenericTarget( )
{
  if(this->setGenericTargetMaterial) {
    if(this->setGenericTargetDimensions) {
      if(this->setGenericTargetPosition) {
        G4String name = this->genericTargetMaterial;
        G4double vec_x = this->genericTargetDimensions.x()/mm;
        G4double vec_y = this->genericTargetDimensions.y()/mm;
        G4double vec_z = this->genericTargetDimensions.z()/mm;
        ApparatusGenericTarget* pApparatusGenericTarget = new ApparatusGenericTarget();
        pApparatusGenericTarget->Build(name, vec_x, vec_y, vec_z);
        G4RotationMatrix* rotate = new G4RotationMatrix;
        pApparatusGenericTarget->PlaceApparatus(this->hallLog, this->genericTargetPosition, rotate);
      }
      else {
      }
    }
    else {
    }
  }
  else {
  }
}

void DetectorConstruction::SetFieldBoxMaterial( G4String name )
{
  this->setFieldBoxMaterial = true;
  this->fieldBoxMaterial = name;
  SetFieldBox();
}

void DetectorConstruction::SetFieldBoxDimensions( G4ThreeVector vec )
{
  this->setFieldBoxDimensions = true;
  this->fieldBoxDimensions = vec;
  SetFieldBox();
}

void DetectorConstruction::SetFieldBoxPosition( G4ThreeVector vec )
{
  this->setFieldBoxPosition = true;
  this->fieldBoxPosition = vec;
  SetFieldBox();
}

void DetectorConstruction::SetFieldBoxMagneticField( G4ThreeVector vec )
{
  this->setFieldBoxMagneticField = true;
  this->fieldBoxMagneticField = vec;
  SetFieldBox();
}

void DetectorConstruction::SetFieldBox( )
{
  if(this->setFieldBoxMagneticField) {
    if(this->setFieldBoxMaterial) {
      if(this->setFieldBoxDimensions) {
        if(this->setFieldBoxPosition) {
          G4String name = this->fieldBoxMaterial;
          G4double vec_x = this->fieldBoxDimensions.x()/mm;
          G4double vec_y = this->fieldBoxDimensions.y()/mm;
          G4double vec_z = this->fieldBoxDimensions.z()/mm;
          ApparatusFieldBox* pApparatusFieldBox = new ApparatusFieldBox();
          pApparatusFieldBox->Build(name, vec_x, vec_y, vec_z, this->fieldBoxMagneticField);
          G4RotationMatrix* rotate = new G4RotationMatrix;
          pApparatusFieldBox->PlaceApparatus(this->hallLog, this->fieldBoxPosition, rotate);
        }
        else {
        }
      }
      else {
      }
    }
    else {
    }
  }
  else {
  }
}


void DetectorConstruction::AddApparatusSpiceTargetChamber()
{
   //Create Target Chamber
   ApparatusSpiceTargetChamber* myApparatusSpiceTargetChamber = new ApparatusSpiceTargetChamber();
   myApparatusSpiceTargetChamber->Build(this->hallLog); 
}

void DetectorConstruction::AddApparatus8piVacuumChamber()
{
   //Create Vacuum Chamber
   Apparatus8piVacuumChamber* myApparatus8piVacuumChamber = new Apparatus8piVacuumChamber();
   myApparatus8piVacuumChamber->Build(this->hallLog);
}

void DetectorConstruction::AddApparatus8piVacuumChamberDelrinShell(G4int thickness)
{
   //Create Shell Around Vacuum Chamber
   Apparatus8piVacuumChamberDelrinShell* myApparatus8piVacuumChamberDelrinShell = new Apparatus8piVacuumChamberDelrinShell();
   myApparatus8piVacuumChamberDelrinShell->Build(this->hallLog, thickness);
}

void DetectorConstruction::AddDetectionSystemGammaTracking(G4int ndet)
{
  // Describe Placement
  G4ThreeVector direction = G4ThreeVector(0,0,1);
  G4ThreeVector move = 0.0 * direction;
  G4RotationMatrix* rotate = new G4RotationMatrix;
  rotate->rotateX(0.0);
  rotate->rotateY(0.0);
  rotate->rotateZ(0.0);

  G4int detector_number = 0;

  this->myGammaTracking->PlaceDetector(this->hallLog, move, rotate, detector_number);
}

void DetectorConstruction::AddDetectionSystemBrillance380V1(G4int ndet)
{
  // Describe Placement
  G4double detectorAngles[8][2] = {0};
  G4double theta,phi,position;
  G4ThreeVector move,direction;
	
  detectorAngles[0][0] 	= 0.0;
  detectorAngles[1][0] 	= 45.0;  
  detectorAngles[2][0] 	= 90.0;  
  detectorAngles[3][0] 	= 135.0;  
  detectorAngles[4][0] 	= 180.0;  
  detectorAngles[5][0] 	= 225.0;  
  detectorAngles[6][0] 	= 270.0;  
  detectorAngles[7][0] 	= 315.0;  
  detectorAngles[0][1] 	= 90.0;
  detectorAngles[1][1] 	= 90.0;  
  detectorAngles[2][1] 	= 90.0;  
  detectorAngles[3][1] 	= 90.0;  
  detectorAngles[4][1] 	= 90.0;  
  detectorAngles[5][1] 	= 90.0;  
  detectorAngles[6][1] 	= 90.0;  
  detectorAngles[7][1] 	= 90.0;

  for(G4int detector_number = 0; detector_number < ndet; detector_number++)
  {
    phi = detectorAngles[detector_number][0]*deg; // Creates a ring in phi plane
    theta = detectorAngles[detector_number][1]*deg;     

    direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    position = 11.0*cm + (this->myBrillance380V1->GetDetectorLengthOfUnitsCM()/2.0);
    move = position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector
    rotate->rotateX(theta);
    rotate->rotateY(0);
    rotate->rotateZ(phi+0.5*M_PI); 
      
    this->myBrillance380V1->PlaceDetector(this->hallLog, move, rotate, detector_number);
  }
}

void DetectorConstruction::AddDetectionSystemSodiumIodide(G4int ndet)
{
  // Describe Placement
  G4double detectorAngles[8][2] = {0};
  G4double theta,phi,position;
  G4ThreeVector move,direction;

  detectorAngles[0][0] 	= 0.0;
  detectorAngles[1][0] 	= 45.0;
  detectorAngles[2][0] 	= 90.0;
  detectorAngles[3][0] 	= 135.0;
  detectorAngles[4][0] 	= 180.0;
  detectorAngles[5][0] 	= 225.0;
  detectorAngles[6][0] 	= 270.0;
  detectorAngles[7][0] 	= 315.0;
  detectorAngles[0][1] 	= 90.0;
  detectorAngles[1][1] 	= 90.0;
  detectorAngles[2][1] 	= 90.0;
  detectorAngles[3][1] 	= 90.0;
  detectorAngles[4][1] 	= 90.0;
  detectorAngles[5][1] 	= 90.0;
  detectorAngles[6][1] 	= 90.0;
  detectorAngles[7][1] 	= 90.0;

  for(G4int detector_number = 0; detector_number < ndet; detector_number++)
  {
    phi = detectorAngles[detector_number][0]*deg; // Creates a ring in phi plane
    theta = detectorAngles[detector_number][1]*deg;

    direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    position = 25.0*cm + (this->mySodiumIodide->GetDetectorLengthOfUnitsCM()/2.0);
    move = position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector
    rotate->rotateX(theta);
    rotate->rotateY(0);
    rotate->rotateZ(phi+0.5*M_PI);

    this->mySodiumIodide->PlaceDetector(this->hallLog, move, rotate, detector_number);
  }
}

void DetectorConstruction::AddDetectionSystemGriffinForward(G4int ndet)
{
  G4double theta,phi,position;
  G4ThreeVector move,direction;

  for(G4int detector_number = 0; detector_number < ndet; detector_number++)
  {
    direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    position = this->griffinFwdBackPosition;
    move = position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector

    this->myGriffinForward->PlaceDetector(this->hallLog, move, rotate, detector_number);
  }
}

void DetectorConstruction::AddDetectionSystemGriffinForwardDetector(G4int ndet)
{
  G4double theta,phi,position;
  G4ThreeVector move,direction;

  direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  position = this->griffinFwdBackPosition;
  move = position * direction;

  G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector

  this->myGriffinForward->PlaceDetector(this->hallLog, move, rotate, ndet);
}

void DetectorConstruction::AddDetectionSystemGriffinBack(G4int ndet)
{
  G4double theta,phi,position;
  G4ThreeVector move,direction;

  for(G4int detector_number = 0; detector_number < ndet; detector_number++)
  {

    direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    position = this->griffinFwdBackPosition;
    move = position * direction;

    G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector

    this->myGriffinBack->PlaceDetector(this->hallLog, move, rotate, detector_number);
  }
}

void DetectorConstruction::AddDetectionSystemGriffinBackDetector(G4int ndet)
{
  G4double theta,phi,position;
  G4ThreeVector move,direction;

  direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  position = this->griffinFwdBackPosition;
  move = position * direction;

  G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector

  this->myGriffinBack->PlaceDetector(this->hallLog, move, rotate, ndet);
}

void DetectorConstruction::AddDetectionSystemSceptar(G4int ndet)
{
  this->mySceptar->PlaceDetector(this->hallLog, ndet);
}

void DetectorConstruction::AddDetectionSystemSpice(G4int ndet)
{
  this->mySpice->PlaceDetector(this->hallLog, ndet);
}

void DetectorConstruction::AddDetectionSystemSpiceV02(G4int ndet)
{
  // Place in world !
  G4double phi = 0.0*deg;
  G4double theta = 0.0*deg;

  G4ThreeVector direction = G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  G4double position = 0.0*mm;
  G4ThreeVector move = position * direction;

  G4RotationMatrix* rotate = new G4RotationMatrix; 		//rotation matrix corresponding to direction vector
  rotate->rotateX(0);
  rotate->rotateY(0);
  rotate->rotateZ(0);

  this->mySpiceV02->PlaceDetector(this->hallLog, move, rotate, ndet);

}

void DetectorConstruction::AddDetectionSystemPaces(G4int ndet)
{
  this->myPaces->PlaceDetector(this->hallLog, ndet);
}
