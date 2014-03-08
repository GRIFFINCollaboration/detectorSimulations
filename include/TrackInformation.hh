#ifndef TrackInformation_h
#define TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

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

  private:
    G4int                 originalTrackID;
    G4int				  originalPdg;
    G4ThreeVector         originalPosition;
    G4ThreeVector         originalMomentum;
    G4double              originalEnergy;
    G4double              originalTime;
    
    vector<G4String>      AncestorsBirthVolume; // MHD 02 May 2013
    vector<G4String>      AncestorsDeathVolume; // MHD 02 May 2013
    vector<G4int>         AncestorsPdg;         // MHD 02 May 2013
        
  public:
  //Getters
    inline G4int GetOriginalTrackID() const {return originalTrackID;}
    inline G4int GetOriginalPdg() const {return originalPdg;}
    inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
    inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
    inline G4double GetOriginalEnergy() const {return originalEnergy;}
    inline G4double GetOriginalTime() const {return originalTime;}
    inline vector<G4String> GetAncestorsBirthVolume() const {return AncestorsBirthVolume;}
    inline vector<G4String> GetAncestorsDeathVolume() const {return AncestorsDeathVolume;}
    inline vector<G4int>    GetAncestorsPdg()         const {return AncestorsPdg;}
            
    //Setters
    void   SetAncestorsBirthVolumeElement(G4String ans_birth_vol) { AncestorsBirthVolume.push_back(ans_birth_vol);}
	void   SetAncestorsDeathVolumeElement(G4String ans_death_vol) { AncestorsDeathVolume.push_back(ans_death_vol);} 
	void   SetAncestorsPdgElement(G4int ans_pdg) { AncestorsPdg.push_back(ans_pdg);} 
};

extern G4Allocator<TrackInformation> aTrackInformationAllocator;

inline void* TrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void TrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((TrackInformation*)aTrackInfo);}

#endif

  


    
