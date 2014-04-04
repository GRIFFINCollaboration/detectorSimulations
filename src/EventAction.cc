// *****************************************************************************
// *  New Code. 
// *	
// *	This is a more efficient version of EventAction, implemented 
// *	on September 26, 2013. Event action is used to control each event and  	
// *	stores the data in root files as opposed to the old version which wrote 
// *	the data to a HEX file. SteppingAction is used to process each step in any
// *	given event, and RunAction is used after all events have taken place to 
// * 	write the root file to disk. 
// *  Tyler Ballast
// ***************************************************************************** 

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
// $Id: EventAction.cc,v 1.1 2010-11-08 10:38:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4
#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"

//Root 
#include "RootManager.hh"
#include "HistoManager.hh"

//c++
#include <sstream>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run, HistoManager* histo)
:runAct(run),histoManager(histo)
	{
        evtNb = 0;
		printModulo = 1000; 
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) 
//    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    printf( " ---> Ev.# %5d\r", evtNb);
    G4cout.flush();

	RootManager::instance()->SetEventNumber(evtNb); 
    ClearVariables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
        G4double eventNumber = 0 ;
        G4double stepNumber = 0;
        G4int	 cryNumber  = 0;
        G4int	 detNumber  = 0;
        G4double depEnergy  = 0;
        G4double posx       = 0;
        G4double posy       = 0;
        G4double posz       = 0;
        G4double time       = 0;
        G4double initialDirectionX       = 0;
        G4double initialDirectionY       = 0;
        G4double initialDirectionZ       = 0;
        G4double initialEnergy       	 = 0;
        G4int trackID       = 0;
        G4String volume     = "";
                        
    for(G4int i = 0; i < MAXSTEPS; i++) {
         
        if(stepTracker[1][i] != 0 && histoManager->GetStepTrackerBool()) {

        eventNumber = stepTracker[0][i] ;
        stepNumber = stepTracker[1][i];
        cryNumber  = stepTracker[2][i];
        detNumber  = stepTracker[3][i];
        depEnergy  = stepTracker[4][i];
        posx       = stepTracker[5][i]/mm;
        posy       = stepTracker[6][i]/mm;
        posz       = stepTracker[7][i]/mm;
        time       = stepTracker[8][i]/second;
        initialDirectionX = stepTracker[9][i];
        initialDirectionY = stepTracker[10][i];
        initialDirectionZ = stepTracker[11][i];
        initialEnergy = stepTracker[12][i];
        trackID    = stepTracker[13][i];
        volume    = stepVolume[i];
        
        histoManager->FillNtuple(eventNumber, stepNumber, cryNumber, detNumber, depEnergy, posx, posy, posz, time );     
        RootManager::instance()->FillG4Hit(volume, detNumber, cryNumber, 11, depEnergy, posx, posy, posz, trackID, 11, initialEnergy, initialDirectionX, initialDirectionY, initialDirectionZ);	   // this is one hit of a Hit Collection
        //RootManager::instance()->FillHist(1000/keV);		//optional
        }
		
    }
    
    if (depEnergy>0.0) { // if condition satisfied Sort the HitCollection and make a physical event 
		RootManager::instance()->SortEvent();
		}
    
  FillParticleType() ; 
  FillGridEkin() ;
  FillGriffinCryst() ;
  FillSodiumIodideCryst() ;
  FillLaBrCryst() ;	
  
  // I included 'Cryst' on the following to match your naming convention above. If that 
  // doesn't fit with reality then please do change them. 
  
  FillSceptarCryst() ;
  FillSpiceCryst() ;
  FillPacesCryst() ; 
  Fill8piCryst() ;


  //G4int i =0;
  //histoManager->FillNtuple(stepTracker[0][i], stepTracker[1][i], stepTracker[2][i], stepTracker[3][i], stepTracker[4][i]/keV, stepTracker[5][i]/mm, stepTracker[6][i]/mm, stepTracker[7][i]/mm, stepTracker[8][i]/second );

  //histoManager->FillNtuple(stepTracker[0][i], stepTracker[1][i], stepTracker[2][i], stepTracker[3][i], stepTracker[4][i]/keV, stepTracker[5][i]/mm, stepTracker[6][i]/mm, stepTracker[7][i]/mm, stepTracker[8][i]/second );

  ClearVariables() ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//void EventAction::SetVec3D(G4double x, G4double y, G4double z)
//{
//  //histoManager->FillHisto2D(1, x, y);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void EventAction::ClearVariables()
{
  if(histoManager->GetStepTrackerBool()) {
      stepIndex = 0;
      for (G4int i = 0 ; i < MAXSTEPS; i++) {
        for (G4int j = 0 ; j < NUMSTEPVARS; j++) {
            stepTracker[j][i] = 0;
        }
      }
  }

  for (G4int i = 0 ; i < NUMPARTICLETYPES; i++) {
      particleTypes[i]                  = 0;
  }
  for (G4int i = 0 ; i < MAXNUMDET; i++) {
      gridEKinGammaDet[i]               = 0 ;
      gridTrackGammaDet[i]              = 0 ;
      gridEKinElectronDet[i]            = 0 ;
      gridTrackElectronDet[i]           = 0 ;
      LaBrCrystEnergyDet[i]             = 0 ;
      LaBrCrystTrackDet[i]              = 0 ;
      SodiumIodideCrystEnergyDet[i]     = 0 ;
      SodiumIodideCrystTrackDet[i]      = 0 ;
      
      
      SceptarSquareCrystEnergyDet[i]					= 0 ;
      SceptarSquareCrystTrackDet[i] 					= 0 ;
      SceptarAngledCrystEnergyDet[i]					= 0 ;
      SceptarAngledCrystTrackDet[i] 					= 0 ;

      EightPiCrystEnergyDet[i]					= 0 ; 
      EightPiCrystTrackDet[i]						= 0 ;
      SpiceCrystEnergyDet[i]						= 0 ;
      SpiceCrystTrackDet[i]							= 0 ;
      PacesCrystEnergyDet[i]						= 0 ;
      PacesCrystTrackDet[i]							= 0 ;
      
      
  }
  for (G4int i = 0 ; i < MAXNUMDETGRIFFIN; i++) {
    for (G4int j = 0 ; j < MAXNUMCRYGRIFFIN; j++) {

        GriffinCrystEnergyDet[i][j]                 = 0;
        GriffinCrystTrackDet[i][j]                  = 0;
        GriffinSuppressorBackEnergyDet[i][j]        = 0;
        GriffinSuppressorBackTrackDet[i][j]         = 0;
        GriffinSuppressorLeftExtensionEnergyDet[i][j]   = 0;
        GriffinSuppressorLeftExtensionTrackDet[i][j]    = 0;
        GriffinSuppressorLeftSideEnergyDet[i][j]        = 0;
        GriffinSuppressorLeftSideTrackDet[i][j]         = 0;
        GriffinSuppressorRightExtensionEnergyDet[i][j]   = 0;
        GriffinSuppressorRightExtensionTrackDet[i][j]    = 0;
        GriffinSuppressorRightSideEnergyDet[i][j]        = 0;
        GriffinSuppressorRightSideTrackDet[i][j]         = 0;
    }
  }
  
  // NOTE: Clear the variables from the new Fill___Cryst functions.   
  
}


void EventAction::FillParticleType()
{
    G4int numParticleTypes = 0;
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++) 
		  {
		    if (particleTypes[i] != 0) 
				  { // if particle type 'i' has non-zero counts
				    for (G4int j = 0 ; j< particleTypes[i]; j++) 
						  { // loop over the number of time we saw it
						    //G4cout << "particleTypes[" << i << "] = " << particleTypes[i] << G4endl;
						    histoManager->FillHisto(astats_particle_type_in_each_step, i);
						  }
				  }
		  }	

    // Fill the number of particle types in the event
    for (G4int i = 0 ; i < NUMPARTICLETYPES; i++)
	  {
	    if (particleTypes[i] != 0) 
	    	numParticleTypes++;
	  }
    histoManager->FillHisto(astats_particle_type_in_each_event, numParticleTypes);
}

void EventAction::FillGridEkin()
{
    for (G4int i = 0 ; i < MAXNUMDET; i++)
    {
      if(gridEKinGammaDet[i] > MINENERGYTHRES) 
      {
        // fill energies in each detector
        if(WRITEEKINHISTOS)   histoManager->FillHisto(gridcell_gamma_ekin_det0+i, gridEKinGammaDet[i]);
        if(WRITETRACKLHISTOS) histoManager->FillHisto(gridcell_gamma_trackl_det0+i, gridTrackGammaDet[i]);
      }
    }
    for (G4int i = 0 ; i < MAXNUMDET; i++) 
    {
      if(gridEKinElectronDet[i] > MINENERGYTHRES) 
      {
        // fill energies in each detector
        if(WRITEEKINHISTOS)   histoManager->FillHisto(gridcell_electron_ekin_det0+i, gridEKinElectronDet[i]);
        if(WRITETRACKLHISTOS) histoManager->FillHisto(gridcell_electron_trackl_det0+i, gridTrackElectronDet[i]);
      }
    }
}

void EventAction::FillGriffinCryst()
{
    G4double  energySum = 0;
    G4double  energySumDet = 0;
    G4bool SuppressorBackFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorExtensionFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorSideFired[MAXNUMDETGRIFFIN] = {0};
    G4bool SuppressorFired = false;

    // Fill Griffin Histos
    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        energySumDet = 0;
        // Find if any suppressors were fired
        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if (GriffinSuppressorBackEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorBackFired[i] = true;
            }
            if (GriffinSuppressorLeftExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorRightExtensionEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorExtensionFired[i] = true;
            }
            if (GriffinSuppressorLeftSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if (GriffinSuppressorRightSideEnergyDet[i][j] > MINENERGYTHRES) {
                SuppressorSideFired[i] = true;
            }
            if ( !SuppressorFired && ( SuppressorBackFired[i] || SuppressorExtensionFired[i] || SuppressorSideFired[i] ) ) {
                SuppressorFired = true;
            }
        }

        for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
            if(GriffinCrystEnergyDet[i][j] > MINENERGYTHRES) {
                // fill energies in each crystal
                if(WRITEEDEPHISTOS)     histoManager->FillHisto((griffin_crystal_unsup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                if(WRITEEDEPHISTOS)     histoManager->FillHisto(griffin_crystal_unsup_edep_cry, GriffinCrystEnergyDet[i][j]);
                if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) { // Suppressor fired?
                    if(WRITEEDEPHISTOS) histoManager->FillHisto((griffin_crystal_sup_edep_det0_cry0+(MAXNUMDETGRIFFIN*j))+i, GriffinCrystEnergyDet[i][j]);
                    if(WRITEEDEPHISTOS) histoManager->FillHisto(griffin_crystal_sup_edep_cry, GriffinCrystEnergyDet[i][j]);
                }
                energySumDet += GriffinCrystEnergyDet[i][j];
            }
        }
        if(energySumDet > MINENERGYTHRES) {
            // fill energies in each detector
            if(WRITEEDEPHISTOS)     histoManager->FillHisto(griffin_crystal_unsup_edep_det0+i, energySumDet);
            // fill standard energy and track spectra
            if(WRITEEDEPHISTOS)     histoManager->FillHisto(griffin_crystal_unsup_edep, energySumDet);
            if(!SuppressorBackFired[i] && !SuppressorExtensionFired[i] && !SuppressorSideFired[i]) {
                // fill energies in each detector
                if(WRITEEDEPHISTOS) histoManager->FillHisto(griffin_crystal_sup_edep_det0+i, energySumDet);
                // fill standard energy and track spectra
                if(WRITEEDEPHISTOS) histoManager->FillHisto(griffin_crystal_sup_edep, energySumDet);
            }
        }
        energySum += energySumDet;
    }

    if(energySum > MINENERGYTHRES) {
        if(WRITEEDEPHISTOS)     histoManager->FillHisto(griffin_crystal_unsup_edep_sum, energySum);
        if(!SuppressorFired) {
            if(WRITEEDEPHISTOS) histoManager->FillHisto(griffin_crystal_sup_edep_sum, energySum);
        }
    }
}

void EventAction::FillLaBrCryst()
{
    G4double  energySum = 0, trackSum = 0;
    for (G4int i=0; i < MAXNUMDET; i++) {
      if(LaBrCrystEnergyDet[i] > MINENERGYTHRES) {
        // fill energies in each detector
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(labr_crystal_edep_det0+i, LaBrCrystEnergyDet[i]);
        //if(WRITETRACKLHISTOS) histoManager->FillHisto(labr_crystal_trackl_det0+i, LaBrCrystTrackDet[i]);
        // fill standard energy and track spectra
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(labr_crystal_edep, LaBrCrystEnergyDet[i]);
        //if(WRITETRACKLHISTOS) histoManager->FillHisto(labr_crystal_trackl, LaBrCrystTrackDet[i]);
        // add sum energies
        energySum    += LaBrCrystEnergyDet[i];
        trackSum     += LaBrCrystTrackDet[i];
      }
    }
    if(energySum > MINENERGYTHRES) {
      if(WRITEEDEPHISTOS)     histoManager->FillHisto(labr_crystal_edep_sum, energySum);
      //if(WRITETRACKLHISTOS)   histoManager->FillHisto(labr_crystal_trackl_sum, trackSum);
    }
}

void EventAction::FillSodiumIodideCryst()
{
    G4double  energySum = 0, trackSum = 0;
    for (G4int i=0; i < MAXNUMDET; i++) {
      if(SodiumIodideCrystEnergyDet[i] > MINENERGYTHRES) {
        // fill energies in each detector
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sodiumIodide_crystal_edep_det0+i, SodiumIodideCrystEnergyDet[i]);
        //if(WRITETRACKLHISTOS) histoManager->FillHisto(sodiumIodide_crystal_trackl_det0+i, SodiumIodideCrystTrackDet[i]);
        // fill standard energy and track spectra
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sodiumIodide_crystal_edep, SodiumIodideCrystEnergyDet[i]);
        //if(WRITETRACKLHISTOS) histoManager->FillHisto(sodiumIodide_crystal_trackl, SodiumIodideCrystTrackDet[i]);
        // add sum energies
        energySum    += SodiumIodideCrystEnergyDet[i];
        trackSum     += SodiumIodideCrystTrackDet[i];
      }
    }
    if(energySum > MINENERGYTHRES) {
      if(WRITEEDEPHISTOS)     histoManager->FillHisto(sodiumIodide_crystal_edep_sum, energySum);
      //if(WRITETRACKLHISTOS)   histoManager->FillHisto(sodiumIodide_crystal_trackl_sum, trackSum);
    }
}

// These still need to be filled to carry out the appropriate tasks. I am unsure what exactly needs
// to be done in each so I will initialize them and leave it to someone who knows. 
// NOTE: The arrays defined here also need to be cleared in ClearVariables() above. 

void EventAction::FillSceptarCryst() 
{
    G4double  energySum = 0, trackSum = 0;
    for (G4int i=0; i < MAXNUMDET; i++) {
      if(SceptarSquareCrystEnergyDet[i] > MINENERGYTHRES) {
        // fill energies in each detector
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sceptar_square_edep_det0+i, SceptarSquareCrystEnergyDet[i]);
        // fill standard energy and track spectra
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sceptar_square_edep, SceptarSquareCrystEnergyDet[i]);
        // add sum energies
        energySum    += SceptarSquareCrystEnergyDet[i];
        trackSum     += SceptarSquareCrystTrackDet[i];
      }
    }
    for (G4int i=0; i < MAXNUMDET; i++) {
      if(SceptarAngledCrystEnergyDet[i] > MINENERGYTHRES) {
        // fill energies in each detector
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sceptar_angled_edep_det0+i, SceptarAngledCrystEnergyDet[i]);
        // fill standard energy and track spectra
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(sceptar_angled_edep, SceptarAngledCrystEnergyDet[i]);
        // add sum energies
        energySum    += SceptarAngledCrystEnergyDet[i];
        trackSum     += SceptarAngledCrystTrackDet[i];
      }
    }
    if(energySum > MINENERGYTHRES) {
      if(WRITEEDEPHISTOS)     histoManager->FillHisto(sceptar_edep_sum, energySum);
    }
}

void EventAction::Fill8piCryst() 
{

}

void EventAction::FillSpiceCryst() 
{

}

void EventAction::FillS3Cryst() 
{

}

void EventAction::FillPacesCryst() 
{
		G4double  energySum = 0, trackSum = 0;
    for (G4int i=0; i < MAXNUMDET; i++) {
      if(PacesCrystEnergyDet[i] > MINENERGYTHRES) {
        // fill energies in each detector
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(paces_crystal_edep_det0+i, PacesCrystEnergyDet[i]);
        // fill standard energy and track spectra
        if(WRITEEDEPHISTOS)   histoManager->FillHisto(paces_crystal_edep, PacesCrystEnergyDet[i]);
        // add sum energies
        energySum    += PacesCrystEnergyDet[i];
        trackSum     += PacesCrystTrackDet[i];
      }
    }
    if(energySum > MINENERGYTHRES) {
      if(WRITEEDEPHISTOS)     histoManager->FillHisto(paces_crystal_edep_sum, energySum);
    }     
}

//void AddStepTracker(G4int eventNumber, G4int stepNumber, G4int cryNumber, G4int detNumber, G4double depEnergy, G4double posx, G4double posy, G4double posz, G4double time)
//{
//    stepTracker[0][stepIndex] = (G4double)(eventNumber);
//    stepTracker[1][stepIndex] = (G4double)(stepNumber);
//    stepTracker[2][stepIndex] = (G4double)(cryNumber);
//    stepTracker[3][stepIndex] = (G4double)(detNumber);
//    stepTracker[4][stepIndex] = depEnergy;
//    stepTracker[5][stepIndex] = posx;
//    stepTracker[6][stepIndex] = posy;
//    stepTracker[7][stepIndex] = posz;
//    stepTracker[8][stepIndex] = time;
//    stepIndex++;



//}

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


//#include "G4Event.hh"
//#include "G4RunManager.hh"
//#include "G4EventManager.hh"
//#include "G4Point3D.hh"
//#include "G4TrajectoryContainer.hh"
//#include "G4Trajectory.hh"
//#include "G4VVisManager.hh"
//#include "G4ios.hh"
//#include <cstdio>

//#include "SensitiveDetector.hh"
//#include "EventAction.hh"
//#include "RunAction.hh"
//#include "PrimaryGeneratorAction.hh"
//#include "DetectorHit.hh"

//using namespace std;

//EventAction::EventAction()
//{
//  eventID = 0;
//}

//EventAction::~EventAction()
//{}

//void EventAction::BeginOfEventAction(const G4Event* evt)
//{
//  //G4cout << "BeginOfEventAction" << G4endl;
//  eventID = evt->GetEventID();
//  // periodic printing
//  if (eventID < 100 || eventID%500 == 0)  {
//    printf( " >>> Ev.# %5d\r", eventID);
//    G4cout.flush(); 
//  }
//}

//void EventAction::EndOfEventAction(const G4Event* evt)
//{

//	// Writes data to HEX files
//	
//  //G4cout << "EndOfEventAction" << G4endl;
//  G4HCofThisEvent* HCEV = evt->GetHCofThisEvent();
//  if( HCEV==NULL ) {
//    G4cout << " --> Return NULL G4HCofThisEvent" << G4endl;
//    return;
//  }

//  G4int nHCEV = HCEV->GetNumberOfCollections();
//  
//  G4RunManager*           runManager = G4RunManager::GetRunManager();
//  RunAction*              theRun     = (RunAction *)runManager->GetUserRunAction();
//  PrimaryGeneratorAction* theAction  = (PrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction();

//  //////////////////////////
//  // Write to List Mode File
//  //////////////////////////
//  if( theRun->DoWriteLMD() || theRun->DoWriteHEX() ) {

//    if( theRun->DoWriteHEX() ) {theRun->HEXwriteStartFlag();}
////    if( theRun->DoWriteLMD() ) {theRun->LMwrite( theAction->PrepareLine() );} // this writes out initial condition of particle
//    
//    //////////////////////////
//    // write the event hits
//    
//    G4int nn, nhits;
//    DetectorHit  *aHit;
//    G4double      edep;
//    G4double      time;
//    G4double      parte;
//    G4Point3D     posA;
//    G4int         detnum;
//    G4int         segnum;
//    G4String      process;
//    G4String      collection;

//    G4bool        writeTime = true;
//    
//    char          evh[128];
//    G4String      eventheader;

//    for(nn = 0; nn < nHCEV; nn++) 
//    {
//      DetectorHitCollection* theHits = (DetectorHitCollection*)HCEV->GetHC(nn);
//      nhits = theHits->entries();

//      parte      = theAction->GetParticleEnergy(); // Initial Particle Info
//      if(nhits > 0) { // this will only write ONCE per enent, unlike code below.
//         sprintf(evh,"Hits: %d\n",nhits);
//         eventheader = G4String(evh);
//         if( theRun->DoWriteLMD() ) {theRun->LMwrite(eventheader);}
////      if( theRun->DoWriteHEX() ) {
////          theRun->ParticleInfoFlagStart();
////          theRun->ParticleInfoEntry(parte);
////          theRun->ParticleInfoFlagEnd();
////        }
//      }

//      for(G4int i = 0; i < nhits; i++)
//      {

//        aHit       = (*theHits)[i];
//        edep       = aHit->GetEdep();
//        time       = aHit->GetTime();
//        process    = aHit->GetProcess();
//        collection = aHit->GetCollect();
//        posA       = aHit->GetPos();
//        detnum     = aHit->GetDetNb();
//        segnum     = aHit->GetSegNb();

//        if( theRun->DoWriteHEX() && writeTime ) { // Hit Time Info
//          // The first hit stamps the event with a time.
//          theRun->TimeFlagStart();
//          theRun->TimeEntry(time);
//          theRun->TimeFlagEnd();
//          writeTime = false;
//        }

//        if( theRun->DoWriteLMD() ) {theRun->LMwrite( PrepareLine(detnum, edep, time, posA, segnum, process, collection) );}
//        if( theRun->DoWriteHEX() ) {theRun->HEXwrite(detnum, edep, time, posA, segnum, process, collection);}
//      }
//    }
//    if( theRun->DoWriteHEX() ) {theRun->HEXwriteEndFlag();}
//  }  
//}

//G4String EventAction::PrepareLine( G4int ndet, G4double ener, G4double time, G4ThreeVector pos1, G4int nseg, G4String process, G4String collection)
//{
//  G4String dummy1 = "";
//  char     dummy2[128];
//  
//  sprintf( dummy2, "%5d ", ndet );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%2.2d ", nseg );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%9.3f ", ener/keV );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%9.3f ", time/second );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%8.3f %8.3f %8.3f ", pos1.x()/mm, pos1.y()/mm, pos1.z()/mm );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%s ",process.c_str() );
//  dummy1 += G4String(dummy2);

//  sprintf( dummy2, "%s ",collection.c_str() );
//  dummy1 += G4String(dummy2);

//  dummy1.replace( dummy1.size()-1, 1, "\n" );
//  return dummy1;
//}

