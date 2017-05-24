#ifndef TrackInformation_h
#define TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include "G4VProcess.hh"

//c++
#include <vector>
using namespace std; 

class TrackInformation : public G4VUserTrackInformation 
{
  public:
    TrackInformation();
    TrackInformation(const G4Track* aTrack);
    TrackInformation(const TrackInformation* aTrackInfo);
           
    virtual ~TrackInformation();
   
    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);
    inline int operator ==(const TrackInformation& right) const
    {return (this==&right);}

    void Print() const;
   	void Clear() ;
   	void PartialClear() ;
   	
  private:
  // Original particle at source 
    G4int                 originalParentID;
    G4int                 originalTrackID;
    G4int				  originalPdg;
    G4ThreeVector         originalPosition;
    vector<G4ThreeVector> originalTrajectory; // till first impact 
    G4ThreeVector         originalMomentum;
    G4double              originalEnergy;
    G4double              originalTime;
    // At first impact 
    G4String			  originalImpactVolume ; 
    G4ThreeVector         originalImpactPosition;
    G4ThreeVector		  originalImpactMomentum;

// secondaries 
    G4String              currentProcess;
    G4int                 currentTrackID;
    G4int                 currentparentTrackID;
    G4int				  currentPdg;
    G4ThreeVector         currentPositionAtVertex;
    G4ThreeVector         currentPositionAtDetector;
    G4ThreeVector         currentPositionAtDeath; // death or out of world
    G4ThreeVector         currentMomentumAtDeath; // death or out of world 
   
    //At first impact 
    G4String			  currentImpactVolume ; 
    G4ThreeVector         currentImpactPosition;
    G4ThreeVector		  currentImpactMomentum;   
    G4ThreeVector         currentMomentumAtVertex; // contains the angle at the emission
    G4ThreeVector         currentMomentumAtDetector; // contains angle at the impact on the detector 
    
    G4double              currentEnergyAtVertex;
    G4double              currentTimeAtVertex;  
    
    vector<G4String>      SecondariesProcess; 
    vector<G4String>      SecondariesBirthVolume; 
    vector<G4String>      SecondariesDeathVolume; 
    vector<G4int>         SecondariesPdg;         

	G4bool                Tagged;  
	G4bool                TagOriginalImpact; // used to tag the storage for first impact 
	G4bool                TagCurrentImpact;  // used to tag the storage for first impact
      
  public:
     
  //Getters
  //Original particle at source 
    inline G4bool GetTagged(void) const {return Tagged;}
    inline G4bool GetTagOriginalImpact(void) const {return TagOriginalImpact;}
    inline G4bool GetTagCurrentImpact(void) const {return TagCurrentImpact;}
    
    inline G4int GetOriginalParentID(void) const {return originalParentID;}
    inline G4int GetOriginalTrackID(void) const {return originalTrackID;}
    inline G4int GetOriginalPdg(void) const {return originalPdg;}
    inline G4ThreeVector GetOriginalPosition(void) const {return originalPosition;}
    inline G4ThreeVector GetOriginalMomentum(void) const {return originalMomentum;}
    inline G4double GetOriginalEnergy(void) const {return originalEnergy;}
    inline G4double GetOriginalTime(void) const {return originalTime;}

    inline G4String      GetOriginalImpactVolume(void) { return originalImpactVolume; }
    inline G4ThreeVector GetOriginalImpactMomentum(void) { return originalImpactMomentum; }
    inline G4ThreeVector GetOriginalImpactPosition(void) { return originalImpactPosition; }
    inline vector<G4ThreeVector> GetOriginalTrajectory(void) { return originalTrajectory;}
// Secondaries
    inline G4String GetCurrentProcess(void) const {return currentProcess;}
    inline G4int GetCurrentParentID(void) const {return currentparentTrackID;}
    inline G4int GetCurrentTrackID(void) const {return currentTrackID;}
    inline G4int GetCurrentPdg(void) const {return currentPdg;}
    inline G4ThreeVector GetCurrentPositionAtVertex(void) const {return currentPositionAtVertex;}
    inline G4ThreeVector GetCurrentMomentumAtVertex(void) const {return currentMomentumAtVertex;}
	inline G4ThreeVector GetCurrentPositionAtDetector(void) const {return currentPositionAtDetector;}
	inline G4ThreeVector GetCurrentMomentumAtDetector(void) const {return currentMomentumAtDetector;}
	inline G4ThreeVector GetCurrentPositionAtDeath(void) const {return currentPositionAtDeath;}
	inline G4ThreeVector GetCurrentMomentumAtDeath(void) const {return currentMomentumAtDeath;}
    inline G4double GetCurrentEnergyAtVertex(void) const {return currentEnergyAtVertex;}
    inline G4double GetCurrentTimeAtVertex(void) const {return currentTimeAtVertex;}
		// all the vector 
    inline vector<G4String> GetSecondariesProcess(void) const {return SecondariesProcess;}    
    inline vector<G4String> GetSecondariesBirthVolume(void) const {return SecondariesBirthVolume;}
    inline vector<G4String> GetSecondariesDeathVolume(void) const {return SecondariesDeathVolume;}
    inline vector<G4int>    GetSecondariesPdg(void)         const {return SecondariesPdg;}
    	// 1 element
    inline G4String GetSecondariesProcessAt(unsigned i) const {return SecondariesProcess.at(i);}    
    inline G4String GetSecondariesBirthVolumeAt(unsigned i) const {return SecondariesBirthVolume.at(i);}
    inline G4String GetSecondariesDeathVolumeAt(unsigned i) const {return SecondariesDeathVolume.at(i);}
    inline G4int    GetSecondariesPdgAt(unsigned i)         const {return SecondariesPdg.at(i);}
		// size of vectors 
    inline G4int GetSecondariesProcessSize(void) const {return SecondariesProcess.size();}
    inline G4int GetSecondariesBirthVolumeSize(void) const {return SecondariesBirthVolume.size();}
    inline G4int GetSecondariesDeathVolumeSize(void) const {return SecondariesDeathVolume.size();}
    inline G4int GetSecondariesPdgSize(void)         const {return SecondariesPdg.size();}            
  
    inline G4String      GetCurrentImpactVolume(void) { return currentImpactVolume; }
    inline G4ThreeVector GetCurrentImpactMomentum(void) { return currentImpactMomentum; }
    inline G4ThreeVector GetCurrentImpactPosition(void) { return currentImpactPosition; }
     
    //Setters
    inline void   	SetTagged( G4bool tag) { Tagged = tag ;} 
    inline void		SetTagOriginalImpact(G4bool tag) { TagOriginalImpact = tag ;}
    inline void 	SetTagCurrentImpact(G4bool tag)  { TagCurrentImpact = tag ;}    
        
        // original 
    inline void SetOriginalImpactVolume(G4String vol ) { originalImpactVolume = vol ; }
    inline void SetOriginalImpactMomentum(G4ThreeVector momentum) { originalImpactMomentum = momentum;}
    inline void SetOriginalImpactPosition(G4ThreeVector position) { originalImpactPosition = position;}
    inline void SetOriginalTrajectoryElement(G4ThreeVector position) { originalTrajectory.push_back(position);}
	   // secondaries
	inline void SetCurrentProcess(G4String process) { currentProcess = process ;}
	inline void SetCurrentParentID(int id) { currentparentTrackID = id ;}
    inline void SetCurrentTrackID(int id) { currentTrackID = id ;}
    inline void SetCurrentPdg(int pdg) { currentPdg = pdg ;}
    inline void SetCurrentPositionAtVertex(G4ThreeVector position) { currentPositionAtVertex = position;}
    inline void SetCurrentMomentumAtVertex(G4ThreeVector momentum) { currentMomentumAtVertex = momentum;}
	inline void SetCurrentPositionAtDetector(G4ThreeVector position) { currentPositionAtDetector = position ;}
	inline void SetCurrentMomentumAtDetector(G4ThreeVector momentum) { currentMomentumAtDetector = momentum;}
	inline void SetCurrentPositionAtDeath(G4ThreeVector position) { currentPositionAtDeath = position ;}
	inline void SetCurrentMomentumAtDeath(G4ThreeVector position) { currentMomentumAtDeath = position ;}
    inline void SetCurrentEnergyAtVertex(G4double energy ) {currentEnergyAtVertex = energy ;}
    inline void SetCurrentTimeAtVertex(G4double time ) {currentTimeAtVertex = time ;}

    inline void SetCurrentImpactVolume(G4String vol ) { currentImpactVolume = vol ; }
    inline void SetCurrentImpactMomentum(G4ThreeVector momentum) { currentImpactMomentum = momentum;}
    inline void SetCurrentImpactPosition(G4ThreeVector position) { currentImpactPosition = position;}
     
    inline void SetSecondariesProcessElement(G4String ans_birth_vol) { SecondariesProcess.push_back(ans_birth_vol);}
    inline void SetSecondariesBirthVolumeElement(G4String ans_birth_vol) { SecondariesBirthVolume.push_back(ans_birth_vol);}
	inline void SetSecondariesDeathVolumeElement(G4String ans_death_vol) { SecondariesDeathVolume.push_back(ans_death_vol);} 
	inline void SetSecondariesPdgElement(G4int ans_pdg) { SecondariesPdg.push_back(ans_pdg);} 

};

extern G4Allocator<TrackInformation> aTrackInformationAllocator;

inline void* TrackInformation::operator new(size_t){ 
	void* aTrackInfo;
	aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
	return aTrackInfo;
}

inline void TrackInformation::operator delete(void *aTrackInfo){ 
	aTrackInformationAllocator.FreeSingle((TrackInformation*)aTrackInfo);}

#endif

  


    
