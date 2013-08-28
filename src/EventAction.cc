#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4Point3D.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include <cstdio>

#include "SensitiveDetector.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorHit.hh"

using namespace std;

EventAction::EventAction()
{
  eventID = 0;
}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  //G4cout << "BeginOfEventAction" << G4endl;
  eventID = evt->GetEventID();
  // periodic printing
  if (eventID < 100 || eventID%500 == 0)  {
    printf( " >>> Ev.# %5d\r", eventID);
    G4cout.flush(); 
  }
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //G4cout << "EndOfEventAction" << G4endl;
  G4HCofThisEvent* HCEV = evt->GetHCofThisEvent();
  if( HCEV==NULL ) {
    G4cout << " --> Return NULL G4HCofThisEvent" << G4endl;
    return;
  }

  G4int nHCEV = HCEV->GetNumberOfCollections();
  
  G4RunManager*           runManager = G4RunManager::GetRunManager();
  RunAction*              theRun     = (RunAction *)runManager->GetUserRunAction();
  PrimaryGeneratorAction* theAction  = (PrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction();

  //////////////////////////
  // Write to List Mode File
  //////////////////////////
  if( theRun->DoWriteLMD() || theRun->DoWriteHEX() ) {

    if( theRun->DoWriteHEX() ) {theRun->HEXwriteStartFlag();}
//    if( theRun->DoWriteLMD() ) {theRun->LMwrite( theAction->PrepareLine() );} // this writes out initial condition of particle
    
    //////////////////////////
    // write the event hits
    
    G4int nn, nhits;
    DetectorHit  *aHit;
    G4double      edep;
    G4double      time;
    G4double      parte;
    G4Point3D     posA;
    G4int         detnum;
    G4int         segnum;
    G4String      process;
    G4String      collection;

    G4bool        writeTime = true;

    for(nn = 0; nn < nHCEV; nn++) 
    {
      DetectorHitCollection* theHits = (DetectorHitCollection*)HCEV->GetHC(nn);
      nhits = theHits->entries();

      parte      = theAction->GetParticleEnergy(); // Initial Particle Info
      if(nhits > 0) { // this will only write ONCE per enent, unlike code below.
//      if( theRun->DoWriteHEX() ) {
//          theRun->ParticleInfoFlagStart();
//          theRun->ParticleInfoEntry(parte);
//          theRun->ParticleInfoFlagEnd();
//        }
      }

      for(G4int i = 0; i < nhits; i++)
      {

        aHit       = (*theHits)[i];
        edep       = aHit->GetEdep();
        time       = aHit->GetTime();
        process    = aHit->GetProcess();
        collection = aHit->GetCollect();
        posA       = aHit->GetPos();
        detnum     = aHit->GetDetNb();
        segnum     = aHit->GetSegNb();

        if( theRun->DoWriteHEX() && writeTime ) { // Hit Time Info
          // The first hit stamps the event with a time.
          theRun->TimeFlagStart();
          theRun->TimeEntry(time);
          theRun->TimeFlagEnd();
          writeTime = false;
        }

        if( theRun->DoWriteLMD() ) {theRun->LMwrite( PrepareLine(detnum, edep, time, posA, segnum, process, collection) );}
        if( theRun->DoWriteHEX() ) {theRun->HEXwrite(detnum, edep, time, posA, segnum, process, collection);}
      }
    }
    if( theRun->DoWriteHEX() ) {theRun->HEXwriteEndFlag();}
  }  
}

G4String EventAction::PrepareLine( G4int ndet, G4double ener, G4double time, G4ThreeVector pos1, G4int nseg, G4String process, G4String collection)
{
  G4String dummy1 = "";
  char     dummy2[128];
  
  sprintf( dummy2, "%5d ", ndet );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%2.2d ", nseg );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%9.3f ", ener/keV );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%9.3f ", time/second );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%8.3f %8.3f %8.3f ", pos1.x()/mm, pos1.y()/mm, pos1.z()/mm );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%s ",process.c_str() );
  dummy1 += G4String(dummy2);

  sprintf( dummy2, "%s ",collection.c_str() );
  dummy1 += G4String(dummy2);

  dummy1.replace( dummy1.size()-1, 1, "\n" );
  return dummy1;
}
