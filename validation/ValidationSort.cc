// ValidationSort.cc - 5 Sept 2013 - C Unsworth
// To Compile:
// g++ ValidationSort.cc --std=c++0x -o ValidationSort -O2 `root-config --cflags --libs`
//
// Quick sort to process TIGRESS/GRIFFIN GEANT4 simulation output and produce spectra
// Which can be compared with experiment.

// C/C++ libraries:
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

// ROOT libraries:
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

// Definitions
#define CRYSTALS 64   // Total number of HPGe crystals
#define CLOVERS  16   // Total number of clovers in array
#define COLOURS   4   // Number of crystals per clover (Blue,Green,Red,White)

// File pointers:
static TFile *outfile = 0;

// Spectra pointers here:
static TH1F *CrystalEnergy[CRYSTALS];


int main(int argc, char **argv){
   
   // General stuff...
   int i, crystal;
   char name[256], title[256];  // For spectra 
   string str;
   // For holding the data items from each interaction....
   int DetN, SegN; 
   float EkeV, Tsec, Xmm, Ymm, Zmm; 
   char Process[24], Collection[64];
   // For holding data items for the whole event....   
   float CrysEn[CRYSTALS];
   
   //  Open the input file
   ifstream infile;
   if(argc>1) {
      cout << "Opening input file: " << argv[1] << "...." << endl;
      infile.open(argv[1]);
      if(infile) {
         cout << "Success!" << endl;
      }
      else {
         cout << "Could not open file!" << endl;
         return -1;
      }   
   }
       
   // Open output file                                       
   outfile = new TFile("ValidationOut.root","RECREATE"); 
   
   // Initialise spectra...
   for(i=0;i<=CRYSTALS;i++) {
      sprintf(name,"Crystal%d",i); 
      sprintf(title,"Crystal %d Energy (keV)",i); 
      CrystalEnergy[i] = new TH1F(name,title,8192,0,2048);
   }
   
   // Start main loop....
   int line = 0;  // count lines to handle header
   while(!infile.eof()) {  
      getline(infile,str);
      line += 1;
      if(line<10) {continue;} // skip output header, should also skip first occurance of event header "Hits: x:"
      
      if(str.substr(0,4)=="Hits") {  // If this is a new event then fill spectra and reset event data
         for(i=0;i<CRYSTALS;i++) {
            if(CrysEn[i] > 0.0) {CrystalEnergy[i]->Fill(CrysEn[i]);}  
         }   
         memset(CrysEn,0.0,(CRYSTALS*sizeof(float)));
         continue;
      }
      else {  // This is an interaction to be added to event record...
         // Scan input string for expected data items:
         sscanf(str.c_str(),"%i %i %f %f %f %f %f %s %s", &DetN, &SegN, &EkeV, &Tsec, &Xmm, &Ymm, &Zmm, Process, Collection); 
         
         crystal = ((DetN -1) * 4) + SegN -1; // Get crystal number from detector and segment
                     // note: "segment" in GEANT output means crystal in real world 
         
         CrysEn[crystal] += EkeV; // sum energy in this crystal
                     // note: At the moment this is summing HPGe and BGO energy for each colour, need to parse
                     // "Collection" string to fix this        
      } // close secttion for processing interactions   
   }  // Close main while loop
   
   // Now write the spectra....  (I think I can do these all at once but not sure how just yet.
   for(i=0;i<CRYSTALS;i++) {
      CrystalEnergy[i]->Write();
   }
   
   outfile->Close();

   return 0;
}

