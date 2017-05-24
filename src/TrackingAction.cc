#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "TrackInformation.hh"
#include "G4Material.hh"

TrackingAction::TrackingAction()
{;}

TrackingAction::~TrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  
  if(aTrack->GetParentID()==0 && aTrack->GetUserInformation()==0)
  {
    TrackInformation* anInfo = new TrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation(anInfo);
  }
  
    // Mhd 02 May 2013 
    if(aTrack->GetParentID()!=0) {
  
		TrackInformation* oldInfo = (TrackInformation*)aTrack->GetUserInformation();
		TrackInformation* newInfo = oldInfo;
		G4Track* theTrack = (G4Track*)aTrack;

		//Set (append) the Pdg of the ancestor particle
		//newInfo->SetCurrentParentID(newInfo->GetCurrentTrackID());
		newInfo->SetCurrentParentID(theTrack->GetParentID()); 
		
		//Set (append) the Pdg of the ancestor particle
		newInfo->SetCurrentTrackID(theTrack->GetTrackID());
			
		//Set (append) the Pdg of the ancestor particle
		newInfo->SetSecondariesPdgElement(theTrack->GetDefinition()->GetPDGEncoding());

		//Set (append) the process name at birth of the ancestor particle 
	  	if( theTrack->GetNextVolume() != 0 ) 	{
		newInfo->SetSecondariesProcessElement(theTrack->GetCreatorProcess()->GetProcessName());
		} 
		else {
		newInfo->SetSecondariesProcessElement("OOW/Dead"); // OOW = OutOfWorld
		}
				
		//Set (append) the Birth volume of the ancestor particle 
	  	if( theTrack->GetNextVolume() != 0 ) 	{
		newInfo->SetSecondariesBirthVolumeElement(theTrack->GetNextVolume()->GetName());
		} 
		else {
		newInfo->SetSecondariesBirthVolumeElement("OOW/Dead");
		}
	
		// set (append) the new user information 
    	theTrack->SetUserInformation(newInfo);
  }
  
  
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

   // G4cout << " ----------------> in PostUserTrackingAction  "<<G4endl;
	G4Track* theTrack = (G4Track*)aTrack;    
    TrackInformation* newInfo = (TrackInformation*)aTrack->GetUserInformation(); // here it represents the old 

    
// Set (append) the death volume for this track 
  if( theTrack->GetNextVolume() != 0 ) 
	{
	//G4cout << " Adding the next DEATH volume "<< theTrack->GetNextVolume()->GetName() <<G4endl;
	newInfo->SetSecondariesDeathVolumeElement(theTrack->GetNextVolume()->GetName());
	} 
else 
	{
	//G4cout << std::setw(11) << "OutOfWorld" << " "<<G4endl;
	newInfo->SetSecondariesDeathVolumeElement("OOW/Dead");
	}
	
	// set (append) the new user information 
    theTrack->SetUserInformation(newInfo);

     // G4cout << " ----------------> New info After DEATH in Secondary loop"<<G4endl;
  // copying the info for the descendancies 
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    TrackInformation* info = (TrackInformation*)(aTrack->GetUserInformation());
	//info->Print();
    
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      { 
        TrackInformation* infoNew = new TrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
  
}

