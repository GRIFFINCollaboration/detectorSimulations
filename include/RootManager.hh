#ifndef ROOTMANAGER_H 
#define ROOTMANAGER_H 1


//c++
#include <iostream>
#include <map>
#include <ctime>

using namespace std ;

//ROOT
#include<TFile.h>
#include<TH1.h>
#include "TTree.h"

//Geant4
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4Tubs.hh"

//User
#include "RawG4Event.hh"
#include "DetectionSystemSpice.hh"
#include "../dataRootClass/TSpiceData.h"
#include "../dataRootClass/TS3Data.h"
#include "../dataRootClass/TGriffinData.h"


class DetectionSystemSpice;

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
        map<Int_t,RawG4Event> fGeantEvent;
        
        //Writing Class for detectors goes here
        TSpiceData* fSpiceData;
        TS3Data*    fS3Data;      
     
        DetectionSystemSpice *fDetectorSpice; 
                     
       
    public:
    // fill the histograms 
        void FillHist(double);
        
    //fill the map hit by hit    
      void FillG4Hit(int, // key of the detector
					 int, // particle pdg 
					 double, // particle depositid energy
					 double, double, double, // particle position vector
					 int, // original particle Track ID
					 int,// primary particle pdg encoding
					 double,//original (primary) particle energy
					 double, double, double);// primary particle momentum vector

       //Set event number  						 
       void SetEventNumber(int) ;
       
       // Set the branches on the tree
       void SetTree();
       
       //SortEvent
       void SortEvent();
       
       //Set the data in Spice writing Class
       void SetSpiceEvent(int key);
       void SetS3Event(int key);
       
       // Close the root Manager        							
       void Close();  

};


#endif
