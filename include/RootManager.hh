#ifndef ROOTMANAGER_H 
#define ROOTMANAGER_H 1


//c++
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <sstream>
#include <ctime>
#include <vector>

using namespace std ;

// CLHEP
#include "CLHEP/Random/RandGauss.h"

//G4 
#include "TrackInformation.hh"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

//User
#include "RawG4Event.hh"
#include "../dataRootClass/TTigFragment.h"
#include "../dataRootClass/TSpiceData.h"
#include "../dataRootClass/TS3Data.h"
#include "../dataRootClass/TGriffinData.h"
#include "../dataRootClass/TPacesData.h"
#include "../dataRootClass/THistoryData.h"


class RootManager   {
    
    public:
        static RootManager *instance();
        
    private:
        static RootManager *fRootManager;
        RootManager();
        ~RootManager();
                
        //File to stor te tree
        TFile *fRootfile;
        
        // Output Root Tree
        TTree *fOutputTree;
        
        //Some histograms 
        TH1F *fHist;
        
        //G4 Event
        map<string,RawG4Event> fGeantEvent;
        
        //Writing Class for detectors goes here
        TSpiceData* fSpiceData;
        TS3Data*    fS3Data;   
        TPacesData* fPacesData;
        TGriffinData* fGriffinData;  
        
        THistoryData* fHistoryData;
        TTigFragment* 	fFragment;    

    public:
    // fill the histograms 
        void FillHist(double);
        
    //fill the map hit by hit    
      void FillG4Hit(string , // Word representing the key used to identify the detector, and to build mnemonics
					int , // detector 
					int ,  // crystal 	
					int, // particle pdg 
					double, // particle depositid energy
					double, double, double, // particle position vector
					int, // original particle Track ID
					int,// primary particle pdg encoding
					double,//original (primary) particle energy
					double, double, double);// primary particle momentum vector
       
       // clear all vectors  
       void Clear(void);
              
       //Build the mnemonic used in TRIUMF  	
	   string BuildMnemonic(string volume, int detector, int crystal);
       
       // Set the branches on the tree
       void SetTree();
       
       //SortEvent
       void SortEvent(int eventNb);
       
       //Set the data in Spice writing Class
       void SetHistory( vector <TrackInformation*> info );  // History of the LAST particles in a cascade of events (Whether a part or all of this cascade ended in the detector or )
       void SetSpiceEvent(int eventNb, string mnemonic, int Ring, int Seg);
       void SetS3Event(int eventNb, string mnemonic, int Ring, int Seg);
       void SetPacesEvent(int eventNb, string mnemonic, int Ring, int Seg);

       
       void SetGriffinEvent(int key);
       void SetFragmentEvent(string key);
       
       // Close the root Manager        							
       void Close();  
       
  public:
    static double SpiceResolution[2];

};


#endif
