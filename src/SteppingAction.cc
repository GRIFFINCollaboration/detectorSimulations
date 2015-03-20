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
// $Id: SteppingAction.cc,v 1.1 2010-11-08 10:38:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "TrackInformation.hh" 
#include "TrackingAction.hh" 

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                                         EventAction* evt)
:detector(det), eventaction(evt)					 
{
    griffinDetectorMapSet = false;
    numberOfAssemblyVols = 13;
    
  // List of assembly volumes just to keep track:
  // 
  // assembly 
  // leftSuppressorCasingAssembly
  // rightSuppressorCasingAssembly
  // leftSuppressorExtensionAssembly
  // rightSuppressorExtensionAssembly
  // suppressorBackAssembly 
  // extensionSuppressorShellAssembly   
  // backAndSideSuppressorShellAssembly
  // hevimetAssembly 
  // germaniumAssemblyCry[0]
  // germaniumAssemblyCry[1]
  // germaniumAssemblyCry[2]
  // germaniumAssemblyCry[3]

  stepNumber = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4int particleType = 0;
  G4int evntNb;

  det = 0;
  cry = 0;

  // Get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String volname = volume->GetName();
   
  // collect energy and track length step by step
  // As it's called more than once, get the Track and assign to variable
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();

	G4Track* theTrack = aStep->GetTrack();
	G4double stepl = 0.;
	if (theTrack->GetDefinition()->GetPDGCharge() != 0.)
		stepl = aStep->GetStepLength();

  // Track particle type in EVERY step
  //G4cout << "Particle name = " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << G4endl;
  if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma")         particleType = 1;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-")       particleType = 2;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e+")       particleType = 3;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "proton")   particleType = 4;
  else if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "neutron")  particleType = 5;
  else                                                                                  particleType = 0;

  eventaction->AddParticleType(particleType);
  evntNb =  eventaction->GetEventNumber();
   
  // Get the track (extra) information , mainly about the primary particle, but could be used for other stuff
  TrackInformation* info = (TrackInformation*)(aStep->GetTrack()->GetUserInformation()); 
  //info->SetTagged(false);


// ------------- add information to the track info --------------------

// After the very first step
	if (aStep->GetTrack()->GetCurrentStepNumber()== 1 ){ 
        
		if ( aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()==0) { // if out of world 
			
			if (aStep->GetTrack()->GetCreatorProcess() == 0 ) {  // if original particle  (this information will be past by to all the descendants)
				info->SetTagged(true);
				info->SetOriginalImpactVolume("OOW");// Out Of World
				info->SetOriginalImpactPosition(aStep->GetTrack()->GetPosition()) ;  
				info->SetOriginalImpactMomentum(aStep->GetTrack()->GetMomentum()) ;
				}
			else {	// if other particle  (this information will change with every new track)
				info->SetCurrentImpactVolume("OOW");
				info->SetCurrentImpactPosition(aStep->GetTrack()->GetPosition()) ;  
				info->SetCurrentImpactMomentum(aStep->GetTrack()->GetMomentum()) ;
				//cin.get() ;
				}
			}
		else {  // if NOT out of world 
		    // if original particle (this information will be past by to all the descendants)
			if (aStep->GetTrack()->GetCreatorProcess() == 0 ) {
				info->SetOriginalImpactVolume(aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName());
				info->SetOriginalImpactPosition(aStep->GetTrack()->GetPosition()) ;  
				info->SetOriginalImpactMomentum(aStep->GetTrack()->GetMomentum()) ;
				}
			else {// if other particle (this information will change with every new track)
				info->SetCurrentImpactVolume(aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName());
				info->SetCurrentImpactPosition(aStep->GetTrack()->GetPosition()) ;  
				info->SetCurrentImpactMomentum(aStep->GetTrack()->GetMomentum()) ;  
				}
			}
			
		} 


// Tag this track if the particle will hit the detector
	int anchor = 0 ; 
	
// Get volume before and after for each step
	G4String after = "VolAfter" ;
	G4String before = "VolBefore" ;
	if ( aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()!=0) {
		before = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ; 
		after  = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ; 
		}
	else {
		before = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() ; 
		after = "OOW" ;	
		}
	
// If the detector is touched get position and momentum (and time?)
	if( (before.contains("siDetS3") || after.contains("siDetS3")) && !info->GetTagged() ){
				
		// Get position and momentum
		info->SetCurrentPositionAtDetector(aStep->GetTrack()->GetPosition()) ;  
		info->SetCurrentMomentumAtDetector(aStep->GetTrack()->GetMomentum()) ;
		info->SetTagged(true);
		}

// this is the last step, in case the particle is tagged => add some new info and keep the track 
	if (aStep->GetTrack()->GetTrackStatus() == 2 && info->GetTagged() ) {
	
		// Get position 
		info->SetCurrentPositionAtDeath(aStep->GetTrack()->GetPosition()) ;  
		info->SetCurrentMomentumAtDeath(aStep->GetTrack()->GetMomentum()) ;  
		// Get Pdg 
		info->SetCurrentPdg(aStep->GetTrack()->GetDefinition()->GetPDGEncoding());
		// Get Energy and time of the track at vertex 
		info->SetCurrentEnergyAtVertex(aStep->GetTrack()->GetVertexKineticEnergy());
		info->SetCurrentTimeAtVertex(aStep->GetTrack()->GetGlobalTime()); // in ns
		// Get Vertex of the track 
		info->SetCurrentPositionAtVertex(aStep->GetTrack()->GetVertexPosition());
		// Get momentum at the vertex 
		info->SetCurrentMomentumAtVertex(aStep->GetTrack()->GetVertexMomentumDirection());                
		// Get the creator process
		if(aStep->GetTrack()->GetCreatorProcess()!=0)
			info->SetCurrentProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
		else
			info->SetCurrentProcess("source")  ;                     
		
		// Keep the track in the record
		//cout << " <<<<<< Setting info >>>>>> " << endl ; 
		eventaction->SetPrimaryInfo( info ) ;
		//cout << " <<<<<< <<<<<< >>>>>> >>>>>> " << endl ; 
		
		//clear the parts unrelated to the original particle information, and untag for the new track 
		info->PartialClear();  // info->SetTagged(false); included 

		}


  G4int 	OriginID = info->GetOriginalTrackID() ;  
  G4int 	OriginPdg = info->GetOriginalPdg() ;     
  G4double 	OriginEnergy = info->GetOriginalEnergy() ;  // Kinetic Energy                                     
  G4ThreeVector OriginMoment = info->GetOriginalMomentum() ;    
 
  G4StepPoint* point1 = aStep->GetPreStepPoint();
  G4StepPoint* point2 = aStep->GetPostStepPoint();

  G4ThreeVector pos1 = point1->GetPosition();
  G4ThreeVector pos2 = point2->GetPosition();

  //G4double time1 = point1->GetGlobalTime();
  G4double time2 = point2->GetGlobalTime();

  G4int trackID = theTrack->GetTrackID();
  
  size_t found;
  G4String search;

  // Grid Cell
  found = volname.find("gridcell");
  if (ekin != 0 && found!=G4String::npos && particleType == 1) {
      SetDetNumberForGenericDetector(volname);
      eventaction->SetGridEKinGammaDet(ekin,stepl,det-1);
  }

  found = volname.find("gridcell");
  if (ekin != 0 && found!=G4String::npos && particleType == 2) {
      SetDetNumberForGenericDetector(volname);
      eventaction->SetGridEKinElectronDet(ekin,stepl,det-1);
  }

  // Griffin energy deposits
  found = volname.find("germanium_block1");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  found = volname.find("back_quarter_suppressor");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorBackDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  found = volname.find("left_suppressor_extension");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      //G4cout << "left_suppressor_extension Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinSuppressorLeftExtensionDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  found = volname.find("right_suppressor_extension");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorRightExtensionDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  found = volname.find("left_suppressor_casing");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      //G4cout << "left_suppressor_casing_log Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinSuppressorLeftSideDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  found = volname.find("right_suppressor_casing");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForGriffinComponent(volname);
      eventaction->AddGriffinSuppressorRightSideDet(edep,stepl,det-1,cry-1);
      //eventaction->AddStepTracker(evntNb, stepNumber, volname, cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  // Dead layer specific code
  found = volname.find("germanium_dls_block1");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(volname);
      //G4cout << "germanium_dls_block1 Found Edep = " << edep/keV << " keV in det = " << det << " cry = " << cry << " found = " << found << " volname = " << volname << G4endl;
      eventaction->AddGriffinCrystDet(edep,stepl,det-1,cry-1);
      eventaction->AddStepTracker(evntNb, stepNumber, "GRG", cry-1, det-1, edep, pos2.x(), pos2.y(), pos2.z(), time2, OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), OriginEnergy, OriginPdg, OriginID, trackID);
  }

  // LaBr
  found = volname.find("lanthanum_bromide_crystal_block");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddLaBrCrystDet(edep,stepl,det-1);
  }

  // Sodium Iodide
  found = volname.find("sodium_iodide_crystal_block");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSodiumIodideCrystDet(edep,stepl,det-1);
  }

  // Sceptar
  found = volname.find("sceptar_square_scintillator_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSceptarSquareCrystDet(edep,stepl,det-1);
  }

  found = volname.find("sceptar_angled_scintillator_log");
  if (edep != 0 && found!=G4String::npos) {
      SetDetNumberForGenericDetector(volname);
      eventaction->AddSceptarAngledCrystDet(edep,stepl,det-1);
  }

  // Paces
  found = volname.find("paces_silicon_block_log");
  if (edep != 0 && found!=G4String::npos) {   
      SetDetAndCryNumberForPacesDetector(volname);
      eventaction->AddPacesCrystDet(edep,stepl,det);
      eventaction->AddStepTracker(evntNb, stepNumber, "PAC", cry, det, edep, 
      pos2.x(), pos2.y(), pos2.z(), time2, 
      OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), 
      OriginEnergy, OriginPdg, OriginID, trackID);  
  }
  
  //SPICE  
    found = volname.find("siDetSpiceRing");
  if (edep != 0 && found!=G4String::npos) {
      SetDetAndCryNumberForSpiceDetector(volname);
      eventaction->AddSpiceCrystDet(edep,stepl,det);
      eventaction->AddStepTracker(evntNb, stepNumber, "SPI", cry, det, edep, 
      pos2.x(), pos2.y(), pos2.z(), time2, 
      OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), 
      OriginEnergy, OriginPdg, OriginID, trackID);
  }
  
  //S3 of SPICE
    found = volname.find("siDetS3Ring");
  if (edep != 0 && found!=G4String::npos) {
	SetDetAndCryNumberForS3Detector(volname);
	eventaction->AddSpiceCrystDet(edep,stepl,det);
	eventaction->AddStepTracker(evntNb, stepNumber, "SPE", cry, det, edep, 
	pos2.x(), pos2.y(), pos2.z(), time2, 
	OriginMoment.getX(), OriginMoment.getY(), OriginMoment.getZ(), 
	OriginEnergy, OriginPdg, OriginID, trackID);
  }
 
 
   stepNumber++;
 
 
}

void SteppingAction::SetDetAndCryNumberForGriffinComponent(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if( avOver9 == 47 ) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if( avOver99 == 47 ) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')*10+(cstr[5]-'0');  // This was fixed
        impr = cstr[12]-'0';
    }

    det = (G4int)(ceil(((G4double)(av)-5.0)/(G4double)(numberOfAssemblyVols))); // This was fixed
    cry = impr;

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetAndCryNumberForDeadLayerSpecificGriffinCrystal(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int av;
    G4int impr;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        av = cstr[3]-'0';
        impr = cstr[10]-'0';
    }
    else if(avOver99 == 47) { // under 100
        av = (cstr[3]-'0')*10+(cstr[4]-'0');
        impr = cstr[11]-'0';
    }
    else { // OVER 100
        av = (cstr[3]-'0')*100+(cstr[4]-'0')*10+(cstr[5]-'0');
        impr = cstr[12]-'0';
    }

    det = (G4int)(ceil((G4double)(av)/(G4double)(numberOfAssemblyVols)));
    cry = av - numberOfAssemblyVols*(det-1);

    det = FindTrueGriffinDetector(det);

    //G4cout << "Found Edep in " << volname <<  " cry = " << cry << " det = " << det << " av = " << av << G4endl;
}

void SteppingAction::SetDetAndCryNumberForSpiceDetector(G4String volname)
{
    // the volume name contains five underscrores : av_xxx_impr_SegmentID_siDetSpiceRing_RingID_etc...
    G4String dummy="";                          
    size_t UnderScoreIndex[6];
    size_t old = -1 ;  
    for (int i = 0 ; i < 6 ; i++ ){
    UnderScoreIndex[i] = volname.find_first_of("_",old+1);
    old = UnderScoreIndex[i] ;
    }

   dummy = volname.substr (UnderScoreIndex[2]+1,UnderScoreIndex[3]-UnderScoreIndex[2]-1); // select the substring between the underscores 
   cry = atoi(dummy.c_str()) - 1; // subtract one : In spice we start counting ring or sectors from zero 
   
   dummy = volname.substr (UnderScoreIndex[4]+1,UnderScoreIndex[5]-UnderScoreIndex[4]-1);
   det = atoi(dummy.c_str()); // ring 

    //G4cout << " (Stepping action) in " << volname <<  " segment = " << cry << " ring = " << det << G4endl;
    //G4cout << " in " << volname <<  " segment = " << cry << " ring = " << det << G4endl;
    //G4cin.get();
}

void SteppingAction::SetDetAndCryNumberForS3Detector(G4String volname)
{
    // the volume name contains five underscrores : av_xxx_impr_SegmentID_siDetSpiceRing_RingID_etc...
    G4String dummy="";                          
    size_t UnderScoreIndex[6];
    size_t old = -1 ;  
    for (int i = 0 ; i < 6 ; i++ ){
    UnderScoreIndex[i] = volname.find_first_of("_",old+1);
    old = UnderScoreIndex[i] ;
    }

   dummy = volname.substr (UnderScoreIndex[2]+1,UnderScoreIndex[3]-UnderScoreIndex[2]-1); // select the substring between the underscores 
   cry = atoi(dummy.c_str()) - 1; // subtract one : In spice we start counting ring or sectors from zero 
   
   dummy = volname.substr (UnderScoreIndex[4]+1,UnderScoreIndex[5]-UnderScoreIndex[4]-1);
   det = atoi(dummy.c_str()); // ring 

    //G4cout << " in " << volname <<  " segment = " << cry << " ring = " << det << G4endl;
    //G4cin.get();
}


void SteppingAction::SetDetAndCryNumberForPacesDetector(G4String volname)
{
    // the volume name contains five underscrores : av_xxx_impr_SegmentID_siDetSpiceRing_RingID_etc...
    G4String dummy="";                          
    size_t UnderScoreIndex[6];
    size_t old = -1 ;  
    for (int i = 0 ; i < 6 ; i++ ){
		UnderScoreIndex[i] = volname.find_first_of("_",old+1);
		old = UnderScoreIndex[i] ;
		}

   dummy = volname.substr (UnderScoreIndex[2]+1,UnderScoreIndex[3]-UnderScoreIndex[2]-1); // select the substring between the underscores 
   cry = atoi(dummy.c_str()) ; // in paces, start counting from 1 
   
   dummy = volname.substr (UnderScoreIndex[4]+1,UnderScoreIndex[5]-UnderScoreIndex[4]-1);
   det = atoi(dummy.c_str()); // ring 
   
    //G4cout << " (Stepping action) in " << volname <<  " segment = " << cry << " ring = " << det << G4endl;
    //G4cin.get();
}

void SteppingAction::SetDetNumberForGenericDetector(G4String volname)
{
    const char *cstr = volname.c_str();
    G4int volNameOver9;
    G4int avOver9 = cstr[4]-'0';
    G4int avOver99 = cstr[5]-'0';
    if(avOver9 == 47) { // under 10
        volNameOver9 = cstr[11]-'0';
        if(volNameOver9 == 47) {
            det = cstr[10]-'0';
        }
        else {
            det = ((cstr[10]-'0')*10)+volNameOver9 ;
        }
    }
    else if(avOver99 == 47) { // under 100
        volNameOver9 = cstr[12]-'0';
        if(volNameOver9 == 47) {
            det = cstr[11]-'0';
        }
        else {
            det = ((cstr[11]-'0')*10)+volNameOver9 ;
        }
    }
    else { // OVER 100
        volNameOver9 = cstr[13]-'0';
        if(volNameOver9 == 47) {
            det = cstr[12]-'0';
        }
        else {
            det = ((cstr[12]-'0')*10)+volNameOver9 ;
        }
    }

    G4cout << "Stepping Action :: Found electron ekin in " << volname << " det = " << det << G4endl;
}

G4int SteppingAction::FindTrueGriffinDetector(G4int det)
{
    G4int trueDet;
    trueDet = detector->griffinDetectorsMap[det-1];

    return trueDet;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
