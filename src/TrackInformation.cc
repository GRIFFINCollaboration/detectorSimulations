#include "TrackInformation.hh"
#include "G4ios.hh"

G4Allocator<TrackInformation> aTrackInformationAllocator;

TrackInformation::TrackInformation()
  : G4VUserTrackInformation()
{

Clear() ; 

}


void TrackInformation::Clear()
{
    Tagged = false;
    
    originalTrackID = -1;
    originalPdg = -1;
    originalPosition = G4ThreeVector(-1.,-1.,-1.);
    originalMomentum = G4ThreeVector(-1.,-1.,-1.);
    originalEnergy = -1;
    originalTime = -1;
 
 	originalImpactVolume = "None";
 	originalImpactPosition = G4ThreeVector(-1.,-1.,-1.);
 	originalImpactMomentum = G4ThreeVector(-1.,-1.,-1.);

	currentProcess = "None";	   
	currentTrackID = -1;
    currentparentTrackID = -1 ;
    currentPdg = -1 ;
    
    currentPositionAtVertex= G4ThreeVector(-1.,-1.,-1.);
    currentPositionAtDetector= G4ThreeVector(-1.,-1.,-1.);
    currentPositionAtDeath= G4ThreeVector(-1.,-1.,-1.); // death or out of world 
    
    currentMomentumAtVertex= G4ThreeVector(-1.,-1.,-1.); // contains the angle at the emission
    currentMomentumAtDetector= G4ThreeVector(-1.,-1.,-1.); // contains angle at the impact on the detector 
    currentMomentumAtDeath= G4ThreeVector(-1.,-1.,-1.); // if != 0, then the particle has escaped 
        
    currentEnergyAtVertex = -1;
    currentTimeAtVertex = -1; 
  
   	currentImpactVolume = "None";
 	currentImpactPosition = G4ThreeVector(-1.,-1.,-1.);
 	currentImpactMomentum = G4ThreeVector(-1.,-1.,-1.);
 	  
    SecondariesProcess.clear();
    SecondariesBirthVolume.clear();
    SecondariesDeathVolume.clear(); 
    SecondariesPdg.clear();
 
}

void TrackInformation::PartialClear()
{
    Tagged = false;
 
    currentProcess = -1;
	currentTrackID = -1;
    currentparentTrackID = -1 ;
    currentPdg = -1 ;
    
    currentPositionAtVertex= G4ThreeVector(-1.,-1.,-1.);
    currentPositionAtDetector= G4ThreeVector(-1.,-1.,-1.);
    currentPositionAtDeath= G4ThreeVector(-1.,-1.,-1.); // death or out of world 
    
    currentMomentumAtVertex= G4ThreeVector(-1.,-1.,-1.); // contains the angle at the emission
    currentMomentumAtDetector= G4ThreeVector(-1.,-1.,-1.); // contains angle at the impact on the detector 
    currentMomentumAtDeath= G4ThreeVector(-1.,-1.,-1.); // death or out of world   
    
    currentEnergyAtVertex = -1;
    currentTimeAtVertex = -1;
        
    SecondariesProcess.clear();
    SecondariesBirthVolume.clear();
    SecondariesDeathVolume.clear(); 
    SecondariesPdg.clear();
 
}


TrackInformation::TrackInformation(const G4Track* aTrack)
  : G4VUserTrackInformation()
{

    //G4cout<< " this is a track with a parent ID = " <<aTrack->GetParentID()<<G4endl;
    Tagged = false; 
    originalParentID = aTrack->GetParentID();
    originalTrackID = aTrack->GetTrackID();
    originalPdg = aTrack->GetDefinition()->GetPDGEncoding();
    originalPosition = aTrack->GetVertexPosition();
    originalMomentum = aTrack->GetVertexMomentumDirection();
    originalEnergy = aTrack->GetVertexKineticEnergy();
    originalTime = aTrack->GetGlobalTime();
   
    currentTrackID = aTrack->GetTrackID(); 
    currentparentTrackID = aTrack->GetParentID();
    
    // Process name at birth
    if (aTrack->GetCreatorProcess() == 0){
    SecondariesProcess.push_back("Source") ; 
    }
    else {
    SecondariesProcess.push_back(aTrack->GetCreatorProcess()->GetProcessName()) ; 
    }
	// this is only used once when the Track informations are retrived from the primary track 
	SecondariesBirthVolume.push_back(aTrack->GetNextVolume()->GetName()) ;
	// Mhd : Redundant information for the first element but i'll keep it for the sake of coherence
	SecondariesPdg.push_back(aTrack->GetDefinition()->GetPDGEncoding()); 
	// We cant do this step since the particle isn't dead yet  
	// SecondariesDeathVolume.push_back(aTrack->GetNextVolume()->GetName()) ;   
    
}
                       
                       
TrackInformation::TrackInformation(const TrackInformation* aTrackInfo)
  : G4VUserTrackInformation()
{
    Tagged = aTrackInfo->Tagged;
    
    //original at source     
    originalParentID = aTrackInfo->originalParentID;
    originalTrackID = aTrackInfo->originalTrackID;
    originalPdg = aTrackInfo->originalPdg;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    // after impact 
    originalImpactVolume = aTrackInfo->originalImpactVolume;
	originalImpactPosition = aTrackInfo->originalImpactPosition;
	originalImpactMomentum = aTrackInfo->originalImpactMomentum;
 
    // secondaries
    currentProcess = aTrackInfo->currentProcess;
  	currentTrackID = aTrackInfo->currentTrackID;
    currentparentTrackID = aTrackInfo->currentparentTrackID ;
    currentPdg = aTrackInfo->currentPdg ;
    
    currentPositionAtVertex= aTrackInfo->currentPositionAtVertex;
    currentPositionAtDetector= aTrackInfo->currentPositionAtDetector;
    currentPositionAtDeath= aTrackInfo->currentPositionAtDeath; // death or out of world 
    
    currentMomentumAtVertex= aTrackInfo->currentMomentumAtVertex; // contains the angle at the emission
    currentMomentumAtDetector= aTrackInfo->currentMomentumAtDetector; // contains angle at the impact on the detector 
    currentMomentumAtDeath= aTrackInfo->currentMomentumAtDeath; // death or out of world 
        
    currentEnergyAtVertex = aTrackInfo->currentEnergyAtVertex;
    currentTimeAtVertex = aTrackInfo->currentTimeAtVertex; 
      
    SecondariesProcess 	 = aTrackInfo->SecondariesProcess;
    SecondariesBirthVolume = aTrackInfo->SecondariesBirthVolume;
    SecondariesDeathVolume = aTrackInfo->SecondariesDeathVolume;
    SecondariesPdg         = aTrackInfo->SecondariesPdg;   
   
}


TrackInformation::~TrackInformation(){;}


void TrackInformation::Print() const
{
  
    G4cout	<< " ============================================================= " << G4endl ;			
   
       G4cout	<< " +++++ Original at source +++++ " << G4endl ;			
   if (Tagged)
   	G4cout	<< " This Track is tagged !! " << Tagged  << G4endl ; 
   	else 
   	   	G4cout	<< " No tagging is set " << Tagged  << G4endl ; 
   	   	
    G4cout	<< " at Source " << G4endl
    		<< "    Position          " << originalPosition << G4endl  
    		<< "    Momentum          " << originalMomentum << G4endl  
    		<< "Original parent ID    " << originalParentID  << G4endl
			<< "Original track ID     " << originalTrackID  << G4endl
			<< "Original track Pdg    " << originalPdg << G4endl
			<< "Original track Energy " << originalEnergy << G4endl
			<< "Original track Time   " << originalTime << G4endl ;

    G4cout	<< " ----- Original after first impact ----- " << G4endl ;			
    G4cout	<< " at Volume            " << originalImpactVolume  << G4endl		
    		<< "    Position          " << originalImpactPosition << G4endl ; 
    G4cout	<< "    Momentum          " << originalImpactMomentum  << G4endl ;

    G4cout	<< " +++++ Secondaries +++++ " << G4endl ;			
    G4cout	<< "Current track Process " << currentProcess  << G4endl  
    		<< "Current parent ID     " << currentparentTrackID  << G4endl
			<< "Current track ID      " << currentTrackID  << G4endl  
			<< "Current track Pdg     " << currentPdg << G4endl 
			<< "Current track Energy  " << currentEnergyAtVertex << G4endl 
			<< "Current track time    " << currentTimeAtVertex << G4endl ;
    G4cout	<< " ----- Secondary key impact ----- " << G4endl ;			
    G4cout	<< " at Vertex            " << G4endl
    		<< "    Position          " << currentPositionAtVertex << G4endl  
    		<< "    Momentum          " << currentMomentumAtVertex << G4endl ;  

    G4cout	<< " at 1st impact        " << currentImpactVolume  << G4endl		
    		<< "    Position          " << currentImpactPosition << G4endl ; 
    G4cout	<< "    Momentum          " << currentImpactMomentum  << G4endl ;
    
    G4cout	<< " at Detector " << G4endl
    		<< "    Position          " << currentPositionAtDetector << G4endl  
    		<< "    Momentum          " << currentMomentumAtDetector << G4endl ;  

    G4cout	<< " at Death " << G4endl
    		<< "    Position          " << currentPositionAtDeath << G4endl 
    		<< "    Momentum          " << currentMomentumAtDeath << G4endl ;
    		  
    for(unsigned i = 0 ; i < SecondariesProcess.size() ; i++ )
    G4cout      
     << " SecondariesProcess at       " <<i  <<" : "<< SecondariesProcess.at(i) << G4endl;
        
     for(unsigned i = 0 ; i < SecondariesBirthVolume.size() ; i++ )
    G4cout      
     << " SecondariesBirthVolume at   " <<i  <<" : "<< SecondariesBirthVolume.at(i) << G4endl;
   
     for(unsigned i = 0 ; i < SecondariesDeathVolume.size() ; i++ )
    G4cout      
     << " SecondariesDeathVolume at   " <<i  <<" : "<< SecondariesDeathVolume.at(i) << G4endl;
     
     for(unsigned i = 0 ; i < SecondariesPdg.size() ; i++ )
    G4cout      
     << " SecondariesPdg at           " <<i  <<" : "<< SecondariesPdg.at(i) << G4endl;
         
}

