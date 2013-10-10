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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fileName[0] = "g4out";
  factoryOn = false;

  makeHistoIndex = 0;

  // histograms
  for (G4int k=0; k<MAXHISTO; k++) {
    fHistId[k] = 0;
    fHistPt[k] = 0;    
  }
  // ntuple
  for (G4int k=0; k<MAXNTCOL; k++) {
    fNtColId[k] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
  G4String filename;
  G4String title;
  G4String detString;
  G4String cryString;
  G4double xmin;
  G4double xmax;
  G4int    nbins;
//  G4int	 	 detStringNum ; 
	
	// determine the maximum number of file names we will need
	// this corresponds to the maximum number of detectors.
//	if( ( MAXNUMDET >= MAXNUMDETGRIFFIN ) && MAXNUMDET != 0 )
//		detStringNum = MAXNUMDET ; 
//	else
//		detStringNum = MAXNUMDETGRIFFIN ; 
//	
//	G4String detString[detStringNum] ; 
//	
//	// This should reduce calls to G4intToG4String in the future. 
//	for( G4int i = 0 ; i < detStringNum ; i++ ) 
//		detString[i] = G4intToG4String(i);
	
	
	
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  G4String extension = analysisManager->GetFileType();
  fileName[1] = fileName[0] + "." + extension;
      
  // Create directories
  //analysisManager->SetHistoDirectoryName("histo");
  //analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  G4bool fileOpen = analysisManager->OpenFile(fileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName[1] 
           << G4endl;
    return;
  }
  
  // create selected histograms
  analysisManager->SetFirstHistoId(1);

  filename  = "astats_particle_type_in_each_step";
  title     = "Particle Creation";
  nbins     = 20;
  xmin      = 0.;
  xmax      = 20.;
  MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

  filename  = "astats_particle_type_in_each_event";
  title     = "Number of Particle Types in Event";
  nbins     = 100;
  xmin      = 0.;
  xmax      = 100.;
  MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

  if(WRITEEKINHISTOS) {
    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);
        filename  = "gridcell_electron_ekin_det" + detString;
        title     = "EKin in cell (keV)";
        nbins     = EKINNBINS;
        xmin      = EKINXMIN;
        xmax      = EKINXMAX;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
  }

  if(WRITETRACKLHISTOS) {
    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);
        filename  = "gridcell_electron_trackl_det" + detString;
        title     = "Trackl in cell (keV)";
        nbins     = TRACKLNBINS;
        xmin      = TRACKLXMIN;
        xmax      = TRACKLXMAX;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
  }

  if(WRITEEKINHISTOS) {
    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);
        filename  = "gridcell_gamma_ekin_det" + detString;
        title     = "EKin in cell (keV)";
        nbins     = EKINNBINS;
        xmin      = EKINXMIN;
        xmax      = EKINXMAX;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
  }

  if(WRITETRACKLHISTOS) {
    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);
        filename  = "gridcell_gamma_trackl_det" + detString;
        title     = "Trackl in cell (keV)";
        nbins     = TRACKLNBINS;
        xmin      = TRACKLXMIN;
        xmax      = TRACKLXMAX;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
  }


  if(WRITEEDEPHISTOS)
  {
  	// Variables and title used for all detectors
    nbins     = EDEPNBINS;
    xmin      = EDEPXMIN;
    xmax      = EDEPXMAX;
    title     = "Edep in crystal (keV)";
    
  	
		// Griffin Suppressors  
    filename  = "griffin_crystal_sup_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "griffin_crystal_sup_edep_cry";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "griffin_crystal_sup_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    
    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
      detString = G4intToG4String(i);
      
      filename  = "griffin_crystal_sup_edep_det" + detString ;
      MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }

    for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
      for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        detString = G4intToG4String(i);
        cryString = G4intToG4String(j);

        filename  = "griffin_crystal_sup_edep_det" + detString + "_cry" + cryString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
      }
    }
  
    // Griffin Crystal Unsuppressed
    filename  = "griffin_crystal_unsup_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "griffin_crystal_unsup_edep_cry";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "griffin_crystal_unsup_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
      detString = G4intToG4String(i);

      filename  = "griffin_crystal_unsup_edep_det" + detString;
      MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }

    for (G4int j=0; j < MAXNUMCRYGRIFFIN; j++) {
      for (G4int i=0; i < MAXNUMDETGRIFFIN; i++) {
        detString = G4intToG4String(i);
        cryString = G4intToG4String(j);

        filename  = "griffin_crystal_unsup_edep_det" + detString + "_cry" + cryString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
      }
    }

		// Brilliance Detector
    filename  = "labr_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "labr_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "labr_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }

		// Sodium Iodide detector
    filename  = "sodiumIodide_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "sodiumIodide_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "sodiumIodide_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
        
    // Sceptar detector
    filename  = "sceptar_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "sceptar_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "sceptar_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
    
    // 8pi detector
    filename  = "Eightpi_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "Eightpi_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "Eightpi_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }
    
    // spice detector
    filename  = "spice_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "spice_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "spice_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    } 
       
    // paces detector
    filename  = "paces_crystal_edep";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    filename  = "paces_crystal_edep_sum";
    MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);

    for (G4int i=0; i < MAXNUMDET; i++) {
        detString = G4intToG4String(i);

        filename  = "paces_crystal_edep_det" + detString;
        MakeHisto(analysisManager, filename,  title, xmin, xmax, nbins);
    }    
	}

  /////////////////////////////////////////////////////////////////////
  // Create 1 ntuple
  //    
//  analysisManager->CreateNtuple("101", "Edep and TrackL");
//  fNtColId[0] = analysisManager->CreateNtupleDColumn("Eabs");
//  fNtColId[1] = analysisManager->CreateNtupleDColumn("Egap");
//  fNtColId[2] = analysisManager->CreateNtupleDColumn("Labs");
//  fNtColId[3] = analysisManager->CreateNtupleDColumn("Lgap");
//  analysisManager->FinishNtuple();
  
  factoryOn = true;       
  G4cout << "\n----> Histogram Tree is opened in " << fileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  if (factoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
    G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
  }                    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::MakeHisto(G4AnalysisManager* analysisManager, G4String filename,  G4String title, G4double xmin, G4double xmax, G4int nbins)
{
  makeHistoIndex++;
  if (makeHistoIndex >= MAXHISTO) {
    G4cout << "---> Exceeded maximum number of histograms. Increase MAXHISTO in HistoManager.hh" << G4endl;
    exit(1);
  }
  fHistId[makeHistoIndex] = analysisManager->CreateH1(filename, title, nbins, xmin, xmax);
  fHistPt[makeHistoIndex] = analysisManager->GetH1(fHistId[makeHistoIndex]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih > MAXHISTO) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does note xist; xbin= " << xbin << " w= " << weight << G4endl;
    return;
  }
  if (fHistPt[ih]) fHistPt[ih]->fill(xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MAXHISTO) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << "  fac= " << fac << G4endl;
    return;
  }
  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs, G4double trackLGap)
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(fNtColId[0], energyAbs);
  analysisManager->FillNtupleDColumn(fNtColId[1], energyGap);
  analysisManager->FillNtupleDColumn(fNtColId[2], trackLAbs);
  analysisManager->FillNtupleDColumn(fNtColId[2], trackLGap);
  analysisManager->AddNtupleRow();  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if(factoryOn) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " EAbs : mean = " << G4BestUnit(fHistPt[1]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[1]->rms(),  "Energy") 
               << G4endl;
    G4cout                
       << " EGap : mean = " << G4BestUnit(fHistPt[2]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[2]->rms(),  "Energy") 
               << G4endl;
    G4cout 
       << " LAbs : mean = " << G4BestUnit(fHistPt[3]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[3]->rms(),  "Length") 
               << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(fHistPt[4]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[4]->rms(),  "Length") 
               << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String HistoManager::G4intToG4String(G4int value)
{
  G4String theString;
  std::stringstream out;
  out << value;
  theString = out.str();
  return theString;
}
