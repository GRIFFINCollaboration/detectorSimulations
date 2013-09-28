//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: RunAction.cc,v 1.1 2010-11-08 10:38:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo)
:histoManager(histo)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
    
  //histograms
  //
  histoManager->book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //save histograms
  //
  histoManager->PrintStatistic();
  histoManager->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// *****************************************************************************
// *  Old code.
// *	
// *	This section was used in the original (ie. early September 2013) 
// *  codebase but has since been replaced with a newer, more efficient version. 
// *  I have left it here in case I screw anything up while trying to merge the 
// *  two bases together. 
// * 
// *	If this section is to be used, everything above this comment must be removed. 
// *	
// *  Tyler Ballast
// ***************************************************************************** 

//#include "G4Run.hh"
//#include "G4RunManager.hh"
//#include "G4UImanager.hh"
//#include "G4VVisManager.hh"
//#include "G4ios.hh"
//#include <cmath>
//#include <time.h>
//#include <sys/time.h>
//#include <sys/stat.h>
//#include <sys/types.h>
//#include <unistd.h>
//#include "Randomize.hh"

//#include "RunAction.hh"
//#include "RunMessenger.hh"

//#define  FORMAT_VERSION "7.7.7"

//using namespace std;

//RunAction::RunAction()
//{
//  runNumber       = 0;
//  writeLMD        = false;   // don't write list mode file
//  writeHEX        = false;

//  myMessenger     = new RunActionMessenger(this);

//  this->timeUnit = "NoUnit";
//}

//RunAction::~RunAction()
//{
//  delete myMessenger;  
//}

//void RunAction::BeginOfRunAction(const G4Run* aRun)
//{
//  runNumber  = aRun->GetRunID();

//  G4cout << G4endl  << "### Run " << runNumber << " start." << G4endl << G4endl;
//  
//  // Output file in list-mode
//  if( writeLMD ) {
//    if( OpenLMFile() ) {
//      G4cout << " Could not open list-mode file, disabling output! " << G4endl;
//      writeLMD = false;
//    }
//    else
//      G4cout << outFileName << " opened" << G4endl;
//  }
//  
//  if (G4VVisManager::GetConcreteInstance()) {
//    G4UImanager* UI = G4UImanager::GetUIpointer();
//    UI->ApplyCommand("/vis/scene/notifyHandlers");
//  }
//}

//void RunAction::EndOfRunAction(const G4Run* aRun)
//{
//  if( !aRun ) return;
//  G4cout << "\n" << G4endl;
//  if (G4VVisManager::GetConcreteInstance()) {
//    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
//  }
//    
//  if( writeLMD ) {
//    outFileLMD.close();
//    G4cout << outFileName << " file closed" << G4endl;
//  }

//  G4cout << "\n### Run " << runNumber << " stop." << G4endl;
//}

//void RunAction::WriteHeader()
//{
//  time_t vTime;
//  time( &vTime );

//  outFileLMD << "AGATA " << FORMAT_VERSION << G4endl;
//  
//  outFileLMD << "DATE " << ctime(&vTime);  // date has already end-of-line character!

//  outFileLMD << "GEOMETERY" << G4endl << "SHELL" << G4endl << "150 240" << G4endl << "ENDGEOMETRY" << G4endl;  

//  outFileLMD << "GAMMA 1" << G4endl << "1000.500" << G4endl;

//  outFileLMD << "$" << G4endl;    
//}

//G4int RunAction::OpenLMFile()
//{
//  char evname[32];
//  
//  sprintf(evname, "events.%4.4d", runNumber);
//  
//  outFileName = G4String(evname);
//  
//  outFileLMD.open(outFileName);
//  if( !outFileLMD.is_open() )
//    return 1;
//  
//  WriteHeader();
//  return 0;
//}

//void RunAction::HEXwrite(G4int ndet, G4double ener, G4double time, G4ThreeVector pos1, G4int nseg, G4String process, G4String collection)
//{
//  G4double remain;
//  G4double roll;
//  unsigned short sd_type = 0;

//  remain = fmod(((ener/eV)/10),60000); // 0.1 eV's
//  roll = floor(((ener/eV)/10)/60000);

//  if (collection == "CollectionBrillance380V1")
//  {
//    sd_type = 1001;    
//  }
//  else if (collection == "CollectionGammaTracking")
//  {
//    sd_type = 2001;    
//  }
//  else if (collection == "Collection8piGermanium")
//  {
//    sd_type = 3001;    
//    ndet = ndet - (20*0);
//  }
//  else if (collection == "Collection8piInnerBGO")
//  {
//    sd_type = 3002;    
//    ndet = ndet - (20*1);
//  }
//  else if (collection == "Collection8piOuterLowerBGO")
//  {
//    sd_type = 3003;    
//    ndet = ndet - (20*2);
//  }
//  else if (collection == "Collection8piOuterUpperBGO")
//  {
//    sd_type = 3004;    
//    ndet = ndet - (20*3);
//  }
//  else if (collection == "CollectionGriffinForwardGe" || collection == "CollectionGriffinBackGe")
//  {
//    sd_type = 4001;    
//  }
//  else if (collection == "CollectionGriffinForwardLeftCasing" || collection == "CollectionGriffinBackLeftCasing")
//  {
//    sd_type = 4002;    
//  }
//  else if (collection == "CollectionGriffinForwardRightCasing" || collection == "CollectionGriffinBackRightCasing")
//  {
//    sd_type = 4003;    
//  }
//  else if (collection == "CollectionGriffinForwardLeftExtension" || collection == "CollectionGriffinBackLeftExtension")
//  {
//    sd_type = 4004;    
//  }
//  else if (collection == "CollectionGriffinForwardRightExtension" || collection == "CollectionGriffinBackRightExtension")
//  {
//    sd_type = 4005;    
//  }
//  else if (collection == "CollectionGriffinForwardBackPlug" || collection == "CollectionGriffinBackBackPlug")
//  {
//    sd_type = 4006;    
//  }
//  else if (collection == "CollectionSquareScint")
//  {
//    sd_type = 5001;
//    ndet = ndet - (20*5);
//  }
//  else if (collection == "CollectionAngledScint")
//  {
//    sd_type = 5001;
//    ndet = ndet - (20*6);
//  }
//  else if (collection == "CollectionSpice" || collection == "CollectionSpiceV02")
//  {
//    sd_type = 6001;    
//  }
//  else if (collection == "CollectionPacesSilicon")
//  {
//    sd_type = 7001;
//    ndet = ndet -1;
//  }
//  else if (collection == "CollectionSodiumIodide")
//  {
//    sd_type = 3001;
//  }
//  else
//  {
//    sd_type = 0;
//    G4cout << "No Collection Name Found!" << G4endl;
//    exit(1);
//  }

//  if (process == "eBrem")
//  {
//	// do something
//  }
//  else if (process == "annihil")
//  {
//	// do something
//  }
//  else
//  {

//  }  

////  if(ndet == 1)
////  {
////    G4cout << "ndet = " << ndet << " collection = " << collection << G4endl;
////  }



//  unsigned short detector = static_cast<unsigned short> (ndet);
//  unsigned short segment = static_cast<unsigned short> (nseg);
//  unsigned short rollover = static_cast<unsigned short> (roll);
//  unsigned short energy = static_cast<unsigned short> (remain+0.5); // outputs to 0.1 eV's (to avoid rounding later on) (goes up to 600keV)
//  //unsigned short itime = static_cast<unsigned short> (((time+0.5*second)/second));
//  //unsigned short posx = static_cast<unsigned short> (pos1.x()/mm);
//  //unsigned short posy = static_cast<unsigned short> (pos1.y()/mm);
//  //unsigned short posz = static_cast<unsigned short> (pos1.z()/mm);




//  if(energy != 0)
//  {
//    ofstream output_file_pointer;
//    output_file_pointer.open("output.dat", ios::app);
//    if (!output_file_pointer.is_open())
//    {
//      exit(1); // error
//    }

//    output_file_pointer.write(reinterpret_cast<char *> (&sd_type), sizeof(unsigned short));
//    output_file_pointer.write(reinterpret_cast<char *> (&detector), sizeof(unsigned short));
//    output_file_pointer.write(reinterpret_cast<char *> (&segment), sizeof(unsigned short));
//    output_file_pointer.write(reinterpret_cast<char *> (&rollover), sizeof(unsigned short));
//    output_file_pointer.write(reinterpret_cast<char *> (&energy), sizeof(unsigned short));

//    output_file_pointer.close();
//  }
//}

//void RunAction::HEXwriteStartFlag()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFA0;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::HEXwriteEndFlag()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFFA;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::TimeFlagStart()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFB0;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::TimeFlagEnd()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFFB;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::TimeEntry(G4double time_value)
//{
//  unsigned long time_output_value;
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open()) {
//    exit(1); // error
//  }

//  if(this->timeUnit == "NoUnit") {
//    //G4cout << "Error 2342342 : No Time Unit Specified. Time = 0" << G4endl;
//    time_output_value = static_cast<unsigned long> (0.0);
//  }
//  else {
//    if(this->timeUnit == "picosecond") {
//      time_output_value = static_cast<unsigned long> (time_value/picosecond + 0.5);
//     }
//    else if(this->timeUnit == "nanosecond") {
//      time_output_value = static_cast<unsigned long> (time_value/nanosecond + 0.5);
//    }
//    else if(this->timeUnit == "microsecond") {
//      time_output_value = static_cast<unsigned long> (time_value/microsecond + 0.5);
//    }
//    else if(this->timeUnit == "millisecond") {
//      time_output_value = static_cast<unsigned long> (time_value/millisecond + 0.5);
//    }
//    else if(this->timeUnit == "second") {
//        time_output_value = static_cast<unsigned long> (time_value/second + 0.5);
//    }
//    else {
//      G4cout << "Error 18181281 : Not a Valid Time Unit" << G4endl;
//      exit(1);
//    }
//  }

//  output_file_pointer.write(reinterpret_cast<char *> (&time_output_value), sizeof(unsigned long));                  

//  output_file_pointer.close();
//}

//void RunAction::ParticleInfoFlagStart()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFC0;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::ParticleInfoFlagEnd()
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open())
//  {
//    exit(1); // error
//  }

//  unsigned short start_flag = 0xFFFC;
//  
//  output_file_pointer.write(reinterpret_cast<char *> (&start_flag), sizeof(unsigned short));
//  output_file_pointer.close();
//}

//void RunAction::ParticleInfoEntry(G4double value)
//{
//  ofstream output_file_pointer;
//  output_file_pointer.open("output.dat", ios::app);
//  if (!output_file_pointer.is_open()) {
//    exit(1); // error
//  }
//  unsigned long output_value = static_cast<unsigned long> ( (value/eV)*10.0 + 0.5); // tens of eV

//  output_file_pointer.write(reinterpret_cast<char *> (&output_value), sizeof(unsigned long));               

//  output_file_pointer.close();
//}



///////////////////////////////////////////////////////////////////////////
//// Methods for the Messenger

//void RunAction::EnableWrite( G4bool enable )
//{
//  writeLMD = enable;

//  if(writeLMD)
//    G4cout << " ---> Event file will be written!" << G4endl; 
//  else
//    G4cout << " ---> Event file will not be written!" << G4endl; 
//}

//void RunAction::EnableWriteHex( G4bool enable )
//{
//  writeHEX = enable;

//  if(writeHEX)
//    G4cout << " ---> Hex Event file will be written!" << G4endl; 
//  else
//    G4cout << " ---> Hex Event file will not be written!" << G4endl; 
//}

//void RunAction::SetRunNumber( G4int number )
//{
//  G4bool doWrite    = false;

//  if( number < 0 ) {
//    G4cout << " Invalid value, could not change run number! " << G4endl;
//  }
//  else {
//    // save status
//    doWrite = writeLMD;
//    if( doWrite )
//      writeLMD = false;
//  
//    G4RunManager * runManager = G4RunManager::GetRunManager();
//    runManager->BeamOn(0);  // this is needed to initialize the run number
//    runManager->SetRunIDCounter( number );
//    G4cout << " ----> Run number has been set to " << number << G4endl;
//    // restore status
//    if( doWrite )
//      writeLMD = true;
//  }
//}

//void RunAction::SetTimeUnit( G4String newTimeUnit )
//{
//  this->timeUnit = newTimeUnit;
//}
