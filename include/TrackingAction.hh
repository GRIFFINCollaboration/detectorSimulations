#ifndef TrackingAction_h
#define TrackingAction_h 1

class G4Track;
#include "G4UserTrackingAction.hh"
#include "globals.hh"

class TrackingAction : public G4UserTrackingAction
{
  public:
  
// Constructor & Destructor
   TrackingAction();
   virtual ~TrackingAction();

// Member functions
   //void SetTrackingManagerPointer(G4TrackingManager* pValue);
   virtual void PreUserTrackingAction(const G4Track*) ;
   virtual void PostUserTrackingAction(const G4Track*) ;


};


#endif


/*
I got this version from another Tip by Makoto Asai
http://geant4.slac.stanford.edu/Tips/event/5.html

#ifndef T01TrackingAction_h
#define T01TrackingAction_h 1

class G4Track;
#include "G4UserTrackingAction.hh"
#include "globals.hh"

class T01TrackingAction : public G4UserTrackingAction 
{
  public:
    T01TrackingAction();
    virtual ~T01TrackingAction();
   
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

  private:
    G4String processName;
};

#endif
*/
