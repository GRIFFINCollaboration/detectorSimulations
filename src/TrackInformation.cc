#include "TrackInformation.hh"
#include "G4ios.hh"

G4Allocator<TrackInformation> aTrackInformationAllocator;

TrackInformation::TrackInformation()
  : G4VUserTrackInformation()
{
    originalTrackID = 0;
    originalPdg = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
    
    AncestorsBirthVolume.clear();
    AncestorsDeathVolume.clear(); 
    AncestorsPdg.clear(); 
}

TrackInformation::TrackInformation(const G4Track* aTrack)
  : G4VUserTrackInformation()
{

    //G4cout<< " this is a track with a parent ID = " <<aTrack->GetParentID()<<G4endl;
    originalTrackID = aTrack->GetTrackID();
    originalPdg = aTrack->GetDefinition()->GetPDGEncoding();
    originalPosition = aTrack->GetPosition();
    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetKineticEnergy();
    originalTime = aTrack->GetGlobalTime();
   
	// MHD : 2 May 2013 ;
	// Not sure if I need to add something here.. 
	// this is only used once when the Track informations are retrived from the primary track 
	AncestorsBirthVolume.push_back(aTrack->GetNextVolume()->GetName()) ;
	// Mhd : Redundant information for the first element but i'll keep it for the sake of coherence'
	AncestorsPdg.push_back(aTrack->GetDefinition()->GetPDGEncoding()); 
	// We cant do this step since the particle isn't dead yet  
	// AncestorsBirthVolume.push_back(aTrack->GetNextVolume()->GetName()) ;   

	/*
	// MHD 25 april 2013
	// detect change of material and make sure the information object is already
	// const G4Step* GetStep() const;   aTrack->GetStep() ; // Not working 
	if(  ( aTrack->GetStep()->GetPreStepPoint()->GetMaterial() !=  aTrack->GetStep()->GetPostStepPoint()->GetMaterial() )  && aTrack->GetUserInformation()!=0 )
	{
	G4cout<<aTrack->GetStep()->GetPreStepPoint()->GetMaterial() << " -> " <<aTrack->GetStep()->GetPreStepPoint()->GetMaterial()<<G4endl;
	flagDensity.push_back(aTrack->GetStep()->GetPostStepPoint()->GetMaterial()->GetDensity());
	}   
	*/
    
}
                       
                       
TrackInformation::TrackInformation(const TrackInformation* aTrackInfo)
  : G4VUserTrackInformation()
{
    originalTrackID = aTrackInfo->originalTrackID;
    originalPdg = aTrackInfo->originalPdg;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    
    AncestorsBirthVolume = aTrackInfo->AncestorsBirthVolume;
    AncestorsDeathVolume = aTrackInfo->AncestorsDeathVolume;
    AncestorsPdg         = aTrackInfo->AncestorsPdg;   
   
}


TrackInformation::~TrackInformation(){;}


void TrackInformation::Print() const
{
    G4cout
			<< " at " << originalPosition
			<< "Original track ID " << originalTrackID  << G4endl
			<< "Original track Pdg " << originalPdg << G4endl
			<< "Original track Energy " << originalEnergy << G4endl
			<< "Original track Time " << originalTime << G4endl ;

    
     for(unsigned i = 0 ; i < AncestorsBirthVolume.size() ; i++ )
    G4cout      
     << " AncestorsBirthVolume at " <<i  <<" : "<< AncestorsBirthVolume.at(i) << G4endl;
   
     for(unsigned i = 0 ; i < AncestorsDeathVolume.size() ; i++ )
    G4cout      
     << " AncestorsDeathVolume at " <<i  <<" : "<< AncestorsDeathVolume.at(i) << G4endl;
     
     for(unsigned i = 0 ; i < AncestorsPdg.size() ; i++ )
    G4cout      
     << " AncestorsPdg at " <<i  <<" : "<< AncestorsPdg.at(i) << G4endl;
          
}

