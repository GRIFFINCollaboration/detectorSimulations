#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

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

#include "NewSquareDetector.hh"


NewSquareDetector::NewSquareDetector()


{

    ///////////////////////////////
    // Square Detector Materials //
    ///////////////////////////////

    this->detector_material = "Silicon";

    ////////////////////////////
    // Dimensions for Crystal //
    ////////////////////////////

    this->square_segment_element_length = 50*mm ; //55*mm for optimum detector size, 30mm for in pipe detector
    this->square_segment_element_width = 36*mm ;    // 45*mm for optimum detector size , 22mm for in pipe detector
    this->square_segment_thickness = 5*mm ; //thickness in X-Y plane

    //////////////////////////////////////
    // No. of Segments in x-y dimension //
    //////////////////////////////////////

    this->no_x_segments = 4. ;
    this->no_y_segments = 3. ;
    this->no_detectors = 4. ;

    //////////////////////////////
    // Dimensions of Guard Ring //
    //////////////////////////////

    this->square_guard_ring_depth = 2.*mm ;
    this->square_guard_ring_thickness = 5.*mm ;

    //////////////////////////////////////////////////
    // How far back the detector is from the target //
    //////////////////////////////////////////////////

    this->SquareDetectorDistance = - 75*mm; // 65*mm for optimum detector size, 120mm for in pipe detector

    /////////////////////////////////////////////////////
    // Relative position in X-Y Plane around Beam Axis //
    /////////////////////////////////////////////////////
    
    //Guard Ring of Detector
    this->BeamPipeXDistanceGR = -7*mm;
    this->BeamPipeYDistanceGR = -7*mm;

    // Active Area of Detector
    this->BeamPipeXDistanceD = (BeamPipeXDistanceGR + square_guard_ring_depth)*mm;
    this->BeamPipeYDistanceD = (BeamPipeYDistanceGR - square_guard_ring_depth)*mm;

}


////////////////////////////
// Delete Logical Volumes //
////////////////////////////

NewSquareDetector::~NewSquareDetector()
{

    delete [] SquareDetect_log;
    delete [] square_guard_ring_log;

}


/////////////////////////////////////////
// Main Build Function in Construction //
/////////////////////////////////////////

G4int NewSquareDetector::Build()
{

    //Build assembly volumes
    this->assembly = new G4AssemblyVolume();

    for(int detectorID=1; detectorID<5; detectorID++)
    {

        this->assemblySquareDet[detectorID] = new G4AssemblyVolume();
        this->assemblyGuardRing[detectorID] = new G4AssemblyVolume();
        //Build logical segment of square detector
        BuildDetectorFace(detectorID);
        BuildGuardRing(detectorID);
    }

    return 1;
} // end Build


/////////////////////////////////////////////////////////////
// Building a Singular Segment of a Singular Detector Face //
/////////////////////////////////////////////////////////////

G4Box* NewSquareDetector::BuildSegment()
{

    //Define the geometry of a single segment
    G4double segment_element_length = ((this->square_segment_element_length - (this->square_guard_ring_depth*2))/no_x_segments)/2;
    G4double segment_element_width = ((this->square_segment_element_width - (this->square_guard_ring_depth*2))/no_y_segments)/2;
    G4double segment_thickness = this->square_segment_thickness/2;

    //cout << "seglength" << segment_element_width <<endl;

    //Creating the logical volume of a single segment
    G4Box* crystal_segment = new G4Box("cystal_segment", segment_element_length, segment_element_width, segment_thickness);

    return crystal_segment;

}//end Build Segment


////////////////////////////////////////////////////////////////
// Building a Singular Guard Ring of a Singular Detector Face //
////////////////////////////////////////////////////////////////

G4int NewSquareDetector::BuildGuardRing(G4int detectorID)
{

    // Set visualization attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(GUARD_COL));
    vis_att->SetVisibility(true);

    //Define the material
    G4Material* material = G4Material::GetMaterial(this->detector_material);

    //Define Outer Edge Geometry of the Guard Ring
    G4double guard_ring_outer_length = this->square_segment_element_length/2 ;
    G4double guard_ring_outer_width = this->square_segment_element_width /2;
    G4double guard_ring_outer_thickness = this->square_segment_thickness /2;

    //Define Inner Edge Geometry of the Guard Ring
    G4double guard_ring_inner_length = (this->square_segment_element_length - (this->square_guard_ring_depth *2))/2 ;
    G4double guard_ring_inner_width = (this->square_segment_element_width - (this->square_guard_ring_depth *2))/2 ;
    G4double guard_ring_inner_thickness = (this->square_guard_ring_thickness)/2 ;

    //Create the two Sections of a Single Guard Ring
    G4Box* outer_guard_ring = new G4Box("outer_guard_ring", guard_ring_outer_length, guard_ring_outer_width, guard_ring_outer_thickness);
    G4Box* inner_guard_ring = new G4Box("inner_guard_ring", guard_ring_inner_length, guard_ring_inner_width, guard_ring_inner_thickness);

    //Defining the Combination of the two Sections forming a singular Guard Ring
    G4ThreeVector Trans = G4ThreeVector(0,0,0);

    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateZ(0*deg);
    G4SubtractionSolid* square_guard_ring = new G4SubtractionSolid("squard_guard_ring", outer_guard_ring, inner_guard_ring,rot, Trans);

    square_guard_ring_log[detectorID] = new G4LogicalVolume(square_guard_ring, material, "square_guard_ring",0,0,0);
    square_guard_ring_log[detectorID]->SetVisAttributes(vis_att);

    // Define rotation and movement of  a singular Guard Ring before creating logical volume
    G4double BeamPipeXDistance = this->BeamPipeXDistanceGR;
    G4double BeamPipeYDistance = this->BeamPipeYDistanceGR;

    G4ThreeVector xdirection 	= G4ThreeVector(1,0,0);
    G4double x_position		    = (this->square_segment_element_length/2) +BeamPipeXDistanceGR ;
    G4ThreeVector xmove 		= x_position * xdirection;

    G4ThreeVector ydirection 	= G4ThreeVector(0,1,0);
    G4double y_position		    = -(this->square_segment_element_width/2.)+BeamPipeYDistanceGR;
    G4ThreeVector ymove 		= y_position * ydirection;

    G4ThreeVector zdirection 	= G4ThreeVector(0,0,1);
    G4double z_position		    = -(this->square_segment_thickness/2.);
    G4ThreeVector zmove 		= z_position * zdirection;

    G4ThreeVector move 		    = (xmove + ymove + zmove);
    G4RotationMatrix* rotate    = new G4RotationMatrix;
    rotate->rotateZ(0*deg);

    //cout << "moveGR" << move << endl;

    //Creating the logical volume of singular guard ring

    this->assemblyGuardRing[detectorID]->AddPlacedVolume(square_guard_ring_log[detectorID], move, rotate);

    return 1;

}//end Build Guard Ring 


///////////////////////////////////////
// Building a Singular Detector Face //
///////////////////////////////////////

G4int NewSquareDetector::BuildDetectorFace(G4int detectorID)
{

    //Define the material to be used
    G4Material* material = G4Material::GetMaterial(this->detector_material);

    //Set Visualisation Attributes
    G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(SIL_COL));
    vis_att->SetVisibility(true);

    // Define rotation and movement of the first segment away from the origin
    G4double BeamPipeXDistanceD = this->BeamPipeXDistanceD;
    G4double BeamPipeYDistanceD = this->BeamPipeYDistanceD;

    G4ThreeVector xdirection 	= G4ThreeVector(1,0,0);
    G4double x_position		    = (((this->square_segment_element_length - (this->square_guard_ring_depth*2))/this->no_x_segments)/2) + BeamPipeXDistanceD;
    G4ThreeVector xrealign         = x_position * xdirection;

    G4ThreeVector ydirection 	= G4ThreeVector(0,1,0);
    G4double y_position		    =(-(((this->square_segment_element_width- (this->square_guard_ring_depth*2))/no_y_segments)/2 )) +  BeamPipeYDistanceD;
    G4ThreeVector yrealign         = y_position * ydirection;

    G4ThreeVector zdirection 	= G4ThreeVector(0,0,1);
    G4double z_position		    = -(this->square_segment_thickness/2.);
    G4ThreeVector zrealign         = z_position * zdirection;

    G4ThreeVector realign 		    = (xrealign + yrealign + zrealign);
    G4RotationMatrix* rotate    = new G4RotationMatrix();
    rotate->rotateZ(0*deg);

    //Build the logical volume of a single segment

    G4Box* SquareDetect = BuildSegment();

    // construct logical volume
                
	G4String detectName = "SquareDetect_";
    G4String detectorNo = G4UIcommand::ConvertToString(detectorID);
	detectName += detectorNo;
    detectName += "_Log";
		
    SquareDetect_log[detectorID] = new G4LogicalVolume(SquareDetect, material, detectName, 0, 0, 0);
	SquareDetect_log[detectorID]->SetVisAttributes(vis_att);
        
    //Using for loops to define the array of segments for one detector crystal    

    for(int xseg = 0; xseg<this->no_x_segments; xseg++)
	{
		for(int yseg=0; yseg<this->no_y_segments; yseg++)
		{

            G4ThreeVector xdirection = G4ThreeVector(1,0,0);
            G4double xposition = ((this->square_segment_element_length - (this->square_guard_ring_depth *2))/this->no_x_segments)*mm;
            G4ThreeVector xmove =  xposition * xdirection ;

            G4ThreeVector ydirection = G4ThreeVector(0,-1,0);
            G4double yposition = ((this->square_segment_element_width - (this->square_guard_ring_depth *2))/this->no_y_segments)*mm;
            G4ThreeVector ymove = yposition * ydirection;

            G4ThreeVector allmove = ((xmove*xseg) + (ymove*yseg));

            G4ThreeVector move = realign + allmove;

            //Create the logical volume of the array
            this->assemblySquareDet[detectorID]->AddPlacedVolume(SquareDetect_log[detectorID], move, rotate);

        }
    }
    
    return 1;

} //end build
                

//////////////////////////////////////////////////
// "place" function called in DetectorMessenger //
//////////////////////////////////////////////////

G4int NewSquareDetector::PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector)
{

    //Define movement of Logical volume in relation to target/source (z-axis)
    G4double SquareDetectorTargetDistance = this->SquareDetectorDistance;
    G4ThreeVector move(0,0,SquareDetectorTargetDistance); 

    //Define the rotation to create the detector plates
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ((360.*deg/this->no_detectors) * detector);

    // Build physical volume of the detectors
    assemblySquareDet[detector]->MakeImprint(exp_hall_log, move, rotate, detector);

  return 1;
} //end place

G4int NewSquareDetector::PlaceGuardRing(G4LogicalVolume* exp_hall_log , G4int detector)
{

    //Define movement of Logical volume in relation to target/source (z-axis)
    G4double SquareDetectorTargetDistance = this->SquareDetectorDistance;
    G4ThreeVector move(0,0,SquareDetectorTargetDistance); 

    //Define the rotation to create the detector guard rings
    G4RotationMatrix* rotate = new G4RotationMatrix;
    rotate->rotateZ((360.*deg/this->no_detectors)* detector);

    // Build physical volume of the guard rings
    assemblyGuardRing[detector]->MakeImprint(exp_hall_log, move, rotate, detector);

  return 1;
} //end place












