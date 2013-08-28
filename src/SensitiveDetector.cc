#include "SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <string>
#include <cmath>
#include "globals.hh"

using namespace std;

SensitiveDetector::SensitiveDetector(G4String name, G4String HCname)
:G4VSensitiveDetector(name)
{
  HCID = -1;
  collectionName.insert(HCname);
}

SensitiveDetector::~SensitiveDetector()
{}

void SensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  if( !HCE ) return;
  theHits     = new DetectorHitCollection(SensitiveDetectorName,collectionName[0]);
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  G4int detNum, segNum, bufferCopyID;

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.)
    return false;

  G4double time = aStep->GetTrack()->GetGlobalTime();

  G4String creator_process = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

//  G4ThreeVector preStepPosition = aStep->GetPreStepPoint()->GetPosition();
//  G4String name1 = aStep->GetTrack()->GetVolume()->GetName();
//  G4String name2 = aStep->GetTrack()->GetNextVolume()->GetName();
//  G4cout << "--> " << preStepPosition.x()/mm << " " << preStepPosition.y()/mm << " " << preStepPosition.z()/mm << G4endl;
//  G4cout << "--> " << name1 << " to " << name2 << G4endl;

  // position of interaction point
  G4ThreeVector position  = aStep->GetPostStepPoint()->GetPosition();
//  G4String volname        = aStep->GetTrack()->GetVolume()->GetName();

  // The following code could be done a better way.
  // Instead of looking for copyIDs (and keeping track of them), 
  // we could just search the collection string for supressors etc.

  G4String findstr("Griffin");
  size_t found;
  found=collectionName[0].find(findstr);
  if (found!=string::npos)
  {  // GRIFFIN // G4cout << "found at: " << int(found) << endl;
    
    detNum = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(depth);

    if(detNum > 4000 && detNum < 4100 ) 	// germanium
    {
      bufferCopyID = detNum - 4000;
    }
    else if(detNum > 4100 && detNum < 4200 ) 	// left_suppressor_side
    {
      bufferCopyID = detNum - 4100;
    }
    else if(detNum > 4200 && detNum < 4300 ) 	// right_suppressor_side
    {
      bufferCopyID = detNum - 4200;
    }
    else if(detNum > 4300 && detNum < 4400 ) 	// left_suppressor_extension
    {
      bufferCopyID = detNum - 4300;
    }
    else if(detNum > 4400 && detNum < 4500 ) 	// right_suppressor_extension
    {
      bufferCopyID = detNum - 4400;
    }      
    else if(detNum > 4500 && detNum < 4600 ) 	// back_suppressor
    {
      bufferCopyID = detNum - 4500;
    }  
    else if(detNum > 5000 && detNum < 5100 ) 	// germanium
    {
      bufferCopyID = detNum - 5000;
    }
    else if(detNum > 5100 && detNum < 5200 ) 	// left_suppressor_side
    {
      bufferCopyID = detNum - 5100;
    }
    else if(detNum > 5200 && detNum < 5300 ) 	// right_suppressor_side
    {
      bufferCopyID = detNum - 5200;
    }
    else if(detNum > 5300 && detNum < 5400 ) 	// left_suppressor_extension
    {
      bufferCopyID = detNum - 5300;
    }
    else if(detNum > 5400 && detNum < 5500 ) 	// right_suppressor_extension
    {
      bufferCopyID = detNum - 5400;
    }
    else if(detNum > 5500 && detNum < 5600 ) 	// back_suppressor
    {
      bufferCopyID = detNum - 5500;
    }
    else
    {
      G4cout << "Error, unknown copy ID for GRIFFIN" << G4endl;
      G4cout << "detNum = " << detNum << G4endl;
      exit(1);
    }

    G4double mod = bufferCopyID % 4;
    G4double div = static_cast<double>(bufferCopyID/4.0);
    G4int ceiling = static_cast<int>(ceil(div));

    detNum = ceiling;
    segNum = mod + 1;
  }
  else
  { // Not GRIFFFIN, so just count
    detNum = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(depth);
    segNum = 1;
  }

//  G4cout << "--> " << volname << "\t" << "det num = " << detNum << "\t" << "bufferCopyID = " << bufferCopyID << "\t" << G4endl;

  DetectorHit* newHit = new DetectorHit();
  
  newHit->SetDetNb( detNum ); 
  newHit->SetSegNb( segNum );                             // segment number 
  newHit->SetEdep( edep );                                // energy release
  newHit->SetTime( time );                                // global time
  newHit->SetProcess( creator_process );                  // the process
  newHit->SetCollect( collectionName[0] );                // the collection name
  newHit->SetPos( position );                             // position of the interaction

  theHits->insert( newHit );

  return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
  if(HCID<0) { 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
  }
  HCE->AddHitsCollection( HCID, theHits );
}  

void SensitiveDetector::DeleteHits()
{
  delete theHits;
}

