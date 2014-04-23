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

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "HistoManager.hh"

class RunAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

static const int MAXSTEPS = 1000;
static const int NUMSTEPVARS = 16;

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*, HistoManager*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);

  G4int GetEventNumber(){return evtNb;};

  void AddStepTracker(G4double eventNumber, G4double stepNumber, G4String volume, 
  						G4double cryNumber, G4double detNumber, 
  						G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time, 
  						G4double originDirectionX, G4double originDirectionY, G4double originDirectionZ, 
  						G4double originEnergy, G4int originPdg, G4int originID, G4int trackID) {
  						
					  if(histoManager->GetStepTrackerBool())   {
						  	stepTracker[0][stepIndex] = eventNumber; 
						  	stepTracker[1][stepIndex] = stepNumber; 
						  	stepTracker[2][stepIndex] = cryNumber; 
						  	stepTracker[3][stepIndex] = detNumber; 
						  	stepTracker[4][stepIndex] = depEnergy; 
						  	stepTracker[5][stepIndex] = posx; 
						  	stepTracker[6][stepIndex] = posy; 
						  	stepTracker[7][stepIndex] = posz; 
						  	stepTracker[8][stepIndex] = time; 
						  	stepTracker[9][stepIndex] = originDirectionX;
							stepTracker[10][stepIndex] = originDirectionY;
							stepTracker[11][stepIndex] = originDirectionZ;
						  	stepTracker[12][stepIndex] = originEnergy;
						  	stepTracker[13][stepIndex] = originPdg;
						  	stepTracker[14][stepIndex] = originID;
						  	stepTracker[15][stepIndex] = trackID;
      
						  	stepVolume[stepIndex] = volume ; 	
						  	stepIndex++; 
						  	if(stepIndex == MAXSTEPS) 	{
						  		G4cout << "\n ----> error 13423549 \n" << G4endl; 
						  		exit(1);
						  		}
						  }; 
 	};

  // particle types
  void AddParticleType(G4int index) {particleTypes[index] += 1;};

  // Grid kinetic energy of gammas and electrons
  void SetGridEKinElectronDet(G4double de, G4double dl, G4int det) { if(gridEKinElectronDet[det] < de) gridEKinElectronDet[det] = de; gridTrackElectronDet[det] += dl;};
  void SetGridEKinGammaDet(G4double de, G4double dl, G4int det) { if(gridEKinGammaDet[det] < de) gridEKinGammaDet[det] = de; gridTrackGammaDet[det] += dl;};
  // Energy deposit in detection systems
  void AddGriffinCrystDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinCrystEnergyDet[det][cry] += de; GriffinCrystTrackDet[det][cry] += dl;};
  void AddGriffinSuppressorBackDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorBackEnergyDet[det][cry] += de; GriffinSuppressorBackTrackDet[det][cry] += dl;};
  void AddGriffinSuppressorLeftExtensionDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorLeftExtensionEnergyDet[det][cry] += de; GriffinSuppressorLeftExtensionTrackDet[det][cry] += dl;};
  void AddGriffinSuppressorLeftSideDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorLeftSideEnergyDet[det][cry] += de; GriffinSuppressorLeftSideTrackDet[det][cry] += dl;};
  void AddGriffinSuppressorRightExtensionDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorRightExtensionEnergyDet[det][cry] += de; GriffinSuppressorRightExtensionTrackDet[det][cry] += dl;};
  void AddGriffinSuppressorRightSideDet(G4double de, G4double dl, G4int det, G4int cry) {GriffinSuppressorRightSideEnergyDet[det][cry] += de; GriffinSuppressorRightSideTrackDet[det][cry] += dl;};

  void AddSodiumIodideCrystDet(G4double de, G4double dl, G4int det) {SodiumIodideCrystEnergyDet[det] += de; SodiumIodideCrystTrackDet[det] += dl;};

  void AddLaBrCrystDet(G4double de, G4double dl, G4int det) {LaBrCrystEnergyDet[det] += de; LaBrCrystTrackDet[det] += dl;};

  void AddSceptarSquareCrystDet(G4double de, G4double dl, G4int det) {SceptarSquareCrystEnergyDet[det] += de; SceptarSquareCrystTrackDet[det] += dl;} ;
  void AddSceptarAngledCrystDet(G4double de, G4double dl, G4int det) {SceptarAngledCrystEnergyDet[det] += de; SceptarAngledCrystTrackDet[det] += dl;} ;

	// NOTE: I am initializing and strucuring these based on AddSodiumIodideCrystDet and AddLaBrCrystDet
	// please correct them if they are wrong. 
  void Add8piCrystDet(G4double de, G4double dl, G4int det) {EightPiCrystEnergyDet[det] += de; EightPiCrystTrackDet[det] += dl;} ;  
  void AddSpiceCrystDet(G4double de, G4double dl, G4int det) {SpiceCrystEnergyDet[det] += de; SpiceCrystTrackDet[det] += dl;} ;
  void AddS3CrystDet(G4double de, G4double dl, G4int det) {SpiceCrystEnergyDet[det] += de; SpiceCrystTrackDet[det] += dl;} ;
  void AddPacesCrystDet(G4double de, G4double dl, G4int det) {PacesCrystEnergyDet[det] += de; PacesCrystTrackDet[det] += dl;} ;



private:

	void ClearVariables();

	void FillParticleType();
	void FillGridEkin();
	void FillGriffinCryst();
	void FillLaBrCryst();
	void FillSodiumIodideCryst();
	void FillSceptarCryst() ;
	void FillSpiceCryst() ;
	void FillS3Cryst();
	void FillPacesCryst() ; 
	void Fill8piCryst() ;

	RunAction*    runAct;
	HistoManager* histoManager;

	G4int     printModulo;
    G4int     evtNb;
    G4bool    stepTrackerBool;

    // tracking info
    G4double stepTracker[NUMSTEPVARS][MAXSTEPS];
    G4String stepVolume[MAXSTEPS]; // volume at each step 
    G4int    stepIndex;

	// Particle types in simulation
	G4int particleTypes[NUMPARTICLETYPES];


	// Grid kinetic energy / track length of gamma and electon
	G4double gridEKinElectronDet[MAXNUMDET];
	G4double gridTrackElectronDet[MAXNUMDET];
	G4double gridEKinGammaDet[MAXNUMDET];
	G4double gridTrackGammaDet[MAXNUMDET];

	// Energy deposit in detection systems
	G4double GriffinCrystEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinCrystTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorBackEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorBackTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

	G4double GriffinSuppressorLeftExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorLeftExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorLeftSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorLeftSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

	G4double GriffinSuppressorRightExtensionEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorRightExtensionTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorRightSideEnergyDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];
	G4double GriffinSuppressorRightSideTrackDet[MAXNUMDETGRIFFIN][MAXNUMCRYGRIFFIN];

	G4double LaBrCrystEnergyDet[MAXNUMDET];
	G4double LaBrCrystTrackDet[MAXNUMDET];

	G4double SodiumIodideCrystEnergyDet[MAXNUMDET];
	G4double SodiumIodideCrystTrackDet[MAXNUMDET];


	// NOTE: I am initializing these arrays based on the LaBr and SodiumIodide definitions. If
	// they are not correct please correct them. 
    G4double SceptarSquareCrystEnergyDet[MAXNUMDET] ;
    G4double SceptarSquareCrystTrackDet[MAXNUMDET] ;

    G4double SceptarAngledCrystEnergyDet[MAXNUMDET] ;
    G4double SceptarAngledCrystTrackDet[MAXNUMDET] ;

	G4double EightPiCrystEnergyDet[MAXNUMDET] ;
	G4double EightPiCrystTrackDet[MAXNUMDET] ; 	

	G4double SpiceCrystEnergyDet[MAXNUMDET] ;
	G4double SpiceCrystTrackDet[MAXNUMDET] ;	
	
	G4double S3CrystEnergyDet[MAXNUMDET] ;
	G4double S3CrystTrackDet[MAXNUMDET] ;	
	
	G4double PacesCrystEnergyDet[MAXNUMDET] ;
	G4double PacesCrystTrackDet[MAXNUMDET] ;
	
};

#endif




// *****************************************************************************
// *  Old code.
// *	
// *	This section was used in the original (ie. early September 2013) 
// *  codebase but has since been replaced with a newer, more efficient version. 
// *  I have left it here in case I screw anything up while trying to merge the 
// *  two bases together. 
// * 
// *	If this section is to be used, everything above this comment must be removed. 
// *	
// *  Tyler Ballast
// ***************************************************************************** 

//#ifndef EventAction_h
//#define EventAction_h 1

//#include "G4UserEventAction.hh"

//using namespace std;

//class G4Event;

//class EventAction : public G4UserEventAction
//{
//  public:
//    EventAction();
//   ~EventAction();

//  public:
//    void BeginOfEventAction(const G4Event*);
//    void EndOfEventAction(const G4Event*);

//  private:
//    G4String PrepareLine( G4int, G4double, G4double, G4ThreeVector, G4int, G4String, G4String);

//  private:
//    G4int  eventID;

//};
//#endif


    
