#ifndef NewSquareDetector_h
#define NewSquareDetector_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#define SIL_COL 0.1,0.9,0.9
#define GUARD_COL 0.1,0.1,0.9

class NewSquareDetector
{

    public:
    NewSquareDetector();
    ~NewSquareDetector();
  
    //////////////////////////////////
    // logical and physical volumes //
    //////////////////////////////////

    private:
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* assemblySquareDet[5];
    G4AssemblyVolume* assemblyGuardRing[5];

  
    public:
    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* exp_hall_log, G4int detector);
    G4int PlaceGuardRing(G4LogicalVolume* exp_hall_log, G4int detector);
  

    private:
    G4LogicalVolume* expHallLog;
  
    private:  
    G4LogicalVolume* square_guard_ring_log[5];
    G4LogicalVolume* SquareDetect_log[5];

  
    ////////////////////////////////////////////////////////////
    // OBS: crystal properties are public, others are private //
    ////////////////////////////////////////////////////////////

    private:
    G4String 	detector_material;
  
    ////////////////////////////////////////////////
    // parameters for the planar detector crystal //
    ////////////////////////////////////////////////

    public:
    G4double 	square_segment_element_length;
    G4double 	square_segment_element_width;
    G4double 	square_segment_thickness;
    G4double  no_x_segments;
    G4double  no_y_segments;
    G4double no_detectors;
    G4double SquareDetectorDistance;
    G4double BeamPipeXDistanceGR;
    G4double BeamPipeYDistanceGR;
    G4double BeamPipeXDistanceD;
    G4double BeamPipeYDistanceD;
	G4int detector_alignment;
	G4double SquareDetectorRotationXY;
	G4double SquareDetectorRotationZ;
	G4double SquareDetectorCorrection;

	G4ThreeVector TranslateDetectors(G4int,G4double,G4double);


    ///////////////////////////////////
    // parameters for the guard ring //
    ///////////////////////////////////

    private:
    G4double 	square_guard_ring_depth;
    G4double 	square_guard_ring_thickness;
   
    /////////////////////////////////
    // internal methods in Build() //
    /////////////////////////////////

    private:
    G4int 	BuildDetectorFace1(G4int detectorID);
    G4int BuildGuardRing1(G4int detectorID);
    G4Box*	BuildSegment1();
    G4int 	BuildDetectorFace2(G4int detectorID);
    G4int BuildGuardRing2(G4int detectorID);
    G4Box*	BuildSegment2();

};

#endif
