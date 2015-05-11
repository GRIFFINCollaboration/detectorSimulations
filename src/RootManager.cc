

#include "RootManager.hh"

double RootManager::SpiceResolution[2];

RootManager *RootManager::fRootManager = NULL;

RootManager *RootManager::instance()    {
    if(!fRootManager)
        fRootManager = new RootManager();   
    return fRootManager;
}



RootManager::RootManager() 
{
	printf("RootManager has been created.\n");    
	fflush(stdout);

	// open root file
	int t = (int) time(0);   // get time now
	fRootfile = new TFile(Form("output_at_%d_seconds.root",t ),"RECREATE");
	
	if(!fRootfile) {
		cout << "\nCould not open file " <<endl;
		exit(-1);
		}


	// Create objects to hold the data
	//spice event
	fSpiceData = new TSpiceData();
	//fS3Data = new TS3Data();
	//
	//fDetectorSpice = new DetectionSystemSpice();

    //Paces event 
	//fPacesData = new TPacesData();

    //New event
    //fNewData = new TNewData();

	//SCEPTAR
	//fSceptarData = new TSceptarData();
	
	//GRIFFIN Event
	//fGriffinData = new TGriffinData();

	//Fragment Event
	//fFragment = new TTigFragment();

	//History
	fHistoryData = new THistoryData();
		
	//histograms 
	//fHist = new TH1F("h","h",500,0,1800);

	//Attach detector branches to the tree
	SetTree();

}

//Destructor
RootManager::~RootManager()  {}


void RootManager::SetTree(){


	//Creating the tree ;
	printf("RootManager : Setting the Tree.\n");    
	fOutputTree = new TTree("Simulated_Data","Simulated Data Tree");
	if(!fOutputTree) {
		cout << "\nCould not create Simulated Data Tree in root file" << endl;
		exit(-1);
		}
	fOutputTree->SetAutoSave(100000);
	
 	/*
	At this stage you can define what branches are written in the tree
	*/
	//fOutputTree->Branch("S3Branch","TS3Data",&fS3Data);
	fOutputTree->Branch("SpiceBranch","TSpiceData",&fSpiceData); 
	//----------------
	//fOutputTree->Branch("TTigFragment","TTigFragment",&fFragment, 1000, 99);
	//----------------
	//fOutputTree->Branch("PacesBranch","TPacesData",&fPacesData);
	//----------------
	//fOutputTree->Branch("NewBranch", "TNewData", &fNewData);
	//----------------
	//fOutputTree->Branch("SceptarBranch","TSceptarData",&fSceptarData);
	//----------------
	// fOutputTree->Branch("GriffinBranch","TGriffinData",&fGriffinData);
	//----------------
	fOutputTree->Branch("HistoryBranch","THistoryData",&fHistoryData);
	//----------------
		
	/*
	Other detector branches goes here
	*/
 	
 } 

void RootManager::FillHist(double temp) {
		float tempf = (float)temp; 
		fHist->Fill(tempf);
	}
					
void RootManager::FillG4Hit(string volume, // Word representing the key used to identify the detector, and to build mnemonics
							int detector,
							int crystal,  	
							int pdg,			// integer representing the type of particle
							double Energy, 			// depositid energy in a step
							double Px, double Py, double Pz,// position vector
							int Id, 			// original(primary) ID
							int PrimPdg,			// primary particle definition        PDG encoding
							double PrimEnergy,		// original(primary) energy
							double Mx, double My, double Mz) { // primary particle momentum vector 
	
	/*			
	Some other functions 
	*/
	//Get the menmonic
	string mnemonic = BuildMnemonic(volume,detector,crystal);
	//Fill the event, calculate the full deposited energy etc...
	fGeantEvent[mnemonic].FillVectors( pdg, Energy, detector,crystal, Px, Py, Pz,  Id,  PrimPdg, PrimEnergy,  Mx, My, Mz);
}




void RootManager::SortEvent(int eventNb) {

	// clear the Fragment 
	//fFragment->ClearVariables();

	// clear the SpiceData object
	//fSpiceData->ClearVariables();
	//fS3Data->ClearVariables();

	// clear the PacesData object
	//fPacesData->ClearVariables();

	//clear the New Object
	//fNewData->ClearVariables();

	//clear the Sceptar Object
	//fSceptarData->ClearVariables();

	// clear the GriffinData object
	//	fGriffinData->ClearVariables();

	// Sort Data from the map second element by getters and set them
	std::map<string,RawG4Event>::iterator it;
	for (std::map<string,RawG4Event>::iterator it=fGeantEvent.begin(); it!=fGeantEvent.end(); ++it) {
		
		string system (it->first, 0, 3); // first three letters of the mnemonic (it->first) defines the system
		//		cout << " RootManager::SortEvent -- system " << system << endl ; 
		//Fragment
		//if (1) SetFragmentEvent(it->first); // take all the event in the fragment tree

		//Spice
		//if (system=="SPI") SetSpiceEvent(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());
		//if (system=="SPE") SetS3Event(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());
		
		//Paces
		//if (system=="PAC") {
		//	SetPacesEvent(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());
		//	fPacesData->Dump();
		//	}

        //NEW
        if (system=="NEW") SetNewEvent(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());

        //Sceptar
        //if (system=="SEP") SetSceptarEvent(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());
		
		//Griffin
		//Need to put Griffin key here
		
		//New Detector
		//if (system==XYZ) SetOtherDetectorEvent(key);
		}

 
// fill the tree 
	fOutputTree->Fill();

// ClearVariables the map
	fGeantEvent.clear();
	
}


string RootManager::BuildMnemonic(string volume, int detector, int crystal) {

	std::string system (volume, 0, 2);
	std::string sub_system (volume, 2, 1);
	std::string number = "00" ; 

	//Build the mnemonic for spice 
	//---------------------  
	if (system == "SP"){ // SPICE SiLi
		if (sub_system == "I") {
			detector = 9 - detector ; // (detector here = ring) spice is upstream the target, according to the mnemonics outer ring is '0' inner ring is 9.
			ostringstream convert;   // stream used for the conversion	
			convert << std::setw(3) << std::setfill('0') << (detector * 12 + crystal);      //  set the width to 3 (xyz), fill the blanks by zeros
			number = convert.str(); // set 'Result' to the contents of the stream
			return ( system + sub_system  + "00XN" + number) ; 
			}
		else if (sub_system == "E") { // SPICE S3
				return system + sub_system  +"00"+"XN"+ number ; 
				}
	 }		
    //--------------------- PACES
	else if (system == "PA" && sub_system == "C") {
		// (detector here = ring), Paces has only one "ring",  =>  detector = 0 
		ostringstream convert;   // stream used for the conversion	
		convert << std::setw(2) << std::setfill('0') << (crystal);      //  set the width to 2, i.e. (xy), fill the blanks by zeros = { 01, ... , 05 }
		number = convert.str(); // set 'Result' to the contents of the stream
		//cout << " RootManager::BuildMnemonic : " << system + sub_system << number <<"XN00" << "  |  " << system + sub_system + number+"XN"+"00" << endl ; 
		return system + sub_system  + number + "XN" + "00"; 
		}
    //--------------------- New PACES  
	else if (system == "NE" && sub_system == "W") {
			detector = 4 - detector ; // New is upstream the target
			ostringstream convert;   // stream used for the conversion	
			convert << std::setw(3) << std::setfill('0') << (detector * 12 + crystal);      //  set the width to 3 (xyz), fill the blanks by zeros
			number = convert.str(); // set 'Result' to the contents of the stream
			return (volume + "00XN" + number) ; 
			}

    //--------------------- Sceptar
	else if (system == "SE"){
		if (sub_system == "P") {
			ostringstream convert;   // stream used for the conversion	
			convert << std::setw(3) << std::setfill('0') << (crystal );      //  set the width to 3 (xyz), fill the blanks by zeros
			number = convert.str(); // set 'Result' to the contents of the stream
			return (volume + "00XN" + number) ; 
			}
		}

    //--------------------- 
	//Griffin 	

    //--------------------- 
	//Other detectors
	
	//------------------ 
	//Default
	 return volume;  

}

	

void RootManager::SetFragmentEvent(string mnemonic) {	
	
		// treat
		fGeantEvent.at(mnemonic).SortPrimary();			
		// get the Energy      
		double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
		// Set the data into the fragment , NB : TTigFragment Class has the all the memebers "public", no setters are used, instead we could assigne the values immediatly. 
		fFragment->ChannelName	= mnemonic ;
		fFragment->ChargeCal = energy ;
	}

void RootManager::SetSpiceEvent(int eventNb, string mnemonic, int Ring, int Seg) {	
	// treat
	fGeantEvent.at(mnemonic).SortPrimary();			

	fSpiceData->SetEventNumber(eventNb) ;
	// get primary
	// Pdg	
	int mult = fGeantEvent.at(mnemonic).GetPrimaryPdgMult(); // inside this particular pad 
    for (int i = 0 ; i<mult ;  i++ )
	fSpiceData->SetPrimaryPdg( fGeantEvent.at(mnemonic).GetPrimaryPdg(i) ) ;

	// Energy      
    mult = fGeantEvent.at(mnemonic).GetPrimaryEnergyMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ )	{
	fSpiceData->SetPrimaryEnergy( fGeantEvent.at(mnemonic).GetPrimaryEnergy(i) ) ;
	}

	// Momentum
	    mult = fGeantEvent.at(mnemonic).GetPrimaryThetaMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ )	{
	fSpiceData->SetPrimaryTheta(fGeantEvent.at(mnemonic).GetPrimaryTheta(i) ) ;
	fSpiceData->SetPrimaryPhi( fGeantEvent.at(mnemonic).GetPrimaryPhi(i) ) ;
	}

	// get the energy 			
	double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
	double stDev = SpiceResolution[1] * energy + SpiceResolution[0];
	double applied_resolution = CLHEP::RandGauss::shoot(energy, stDev); 
	TVector3 pos = fGeantEvent.at(mnemonic).GetSecondHitPosition() ;
	
	// fill the SpiceData object
	// (Th,E)
	fSpiceData->SetSpiceThetaEDetectorNbr(1) ; 
	fSpiceData->SetSpiceThetaEStripNbr(Ring) ;    
	fSpiceData->SetSpiceThetaEEnergy(energy/keV) ;     
	fSpiceData->SetSpiceThetaEResEnergy(applied_resolution) ;


	// (Ph,E)
	fSpiceData->SetSpicePhiEDetectorNbr(1) ;
	fSpiceData->SetSpicePhiEStripNbr(Seg)  ;
	fSpiceData->SetSpicePhiEEnergy(energy/keV) ;
	fSpiceData->SetSpicePhiEResEnergy(applied_resolution) ;

	fSpiceData->SetPositionFirstHit(pos) ;				
}


void RootManager::SetS3Event(int eventNb, string mnemonic, int Ring, int Seg) {	
	// treat
	fGeantEvent.at(mnemonic).SortPrimary();			

	fS3Data->SetEventNumber(eventNb) ;
	// get primary
	// Pdg	
	int mult = fGeantEvent.at(mnemonic).GetPrimaryPdgMult(); // inside this particular pad 
    for (int i = 0 ; i<mult ;  i++ )
	fS3Data->SetPrimaryPdg( fGeantEvent.at(mnemonic).GetPrimaryPdg(i) ) ;

	// Energy      
    mult = fGeantEvent.at(mnemonic).GetPrimaryEnergyMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ )	{
	fS3Data->SetPrimaryEnergy( fGeantEvent.at(mnemonic).GetPrimaryEnergy(i) ) ;
	}

	// Momentum
	    mult = fGeantEvent.at(mnemonic).GetPrimaryThetaMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ )	{
	fS3Data->SetPrimaryTheta(fGeantEvent.at(mnemonic).GetPrimaryTheta(i) ) ;
	fS3Data->SetPrimaryPhi( fGeantEvent.at(mnemonic).GetPrimaryPhi(i) ) ;
	}

	// get the energy 			
	double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
	double stDev = (SpiceResolution[1]*keV) * energy  + (SpiceResolution[0]*keV); // the result is in MeV
	double applied_resolution = CLHEP::RandGauss::shoot(energy, stDev);  
	TVector3 pos = fGeantEvent.at(mnemonic).GetSecondHitPosition() ;

	// fill the S3Data object
	// (Th,E)
	fS3Data->SetS3ThetaEDetectorNbr(1) ; 
	fS3Data->SetS3ThetaEStripNbr(Ring) ;	
	fS3Data->SetS3ThetaEEnergy(energy/keV) ;     
 
	// (Ph,E)
	fS3Data->SetS3PhiEDetectorNbr(1) ;
	fS3Data->SetS3PhiEStripNbr(Seg)  ;
	fS3Data->SetS3PhiEEnergy(applied_resolution/keV) ;

	fS3Data->SetPositionFirstHit(pos) ;		

}


void RootManager::SetPacesEvent(int eventNb, string mnemonic, int Ring, int Seg)
{	

	// treat
	fGeantEvent.at(mnemonic).SortPrimary();	
    
    //cout << " mnemonic " << mnemonic << endl;
    //cout << " Seg " << Seg << endl;
    //cout << " Ring " << Ring << endl;
    
	//fGeantEvent.at(mnemonic).ShowVectorsContent();   	
	//G4cin.get(); 
	  
	fPacesData->SetEventNumber(eventNb) ;
	// get primary
	// Pdg	
	int mult = fGeantEvent.at(mnemonic).GetPrimaryPdgMult(); // inside this particular pad 
    for (int i = 0 ; i<mult ;  i++ )
		fPacesData->SetPrimaryPdg( fGeantEvent.at(mnemonic).GetPrimaryPdg(i) ) ;

	// Energy      
    mult = fGeantEvent.at(mnemonic).GetPrimaryEnergyMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ ) {
		fPacesData->SetPrimaryEnergy( fGeantEvent.at(mnemonic).GetPrimaryEnergy(i) ) ;
		}

	// Momentum
	mult = fGeantEvent.at(mnemonic).GetPrimaryThetaMult(); // this should be the same as above
	for (int i = 0 ; i<mult ;  i++ ) {
		fPacesData->SetPrimaryTheta( fGeantEvent.at(mnemonic).GetPrimaryTheta(i) ) ;
		fPacesData->SetPrimaryPhi( fGeantEvent.at(mnemonic).GetPrimaryPhi(i) ) ;
		}

	// get the energy 			
	double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
	TVector3 pos = fGeantEvent.at(mnemonic).GetFirstHitPosition() ;
    
	// fill the PacesData object
	Ring = Ring ; // irrelevant, in PACES we only have one "ring", should take off this parameter at some point.
	fPacesData->SetPacesDetectorNbr(Seg) ; 
	fPacesData->SetPacesEnergy(energy/keV) ;    
    
    //pos.Dump();
    //getchar();
	fPacesData->SetPositionFirstHit(pos) ;	

}


void RootManager::SetNewEvent(int eventNb, string mnemonic, int detector, int Seg)
{

    //treat
    fGeantEvent.at(mnemonic).SortPrimary();

	fNewData->SetEventNumber(eventNb) ;
    //get primary
    //Pdg
    int mult = fGeantEvent.at(mnemonic).GetPrimaryPdgMult();
    for( int i = 0 ; i<mult ; i++) {
        fNewData->SetPrimaryPdg(fGeantEvent.at(mnemonic).GetPrimaryPdg(i) ) ;
    }

    //Energy
    mult = fGeantEvent.at(mnemonic).GetPrimaryEnergyMult();
    for(int i = 0 ; i<mult ; i++){
        fNewData->SetPrimaryEnergy(fGeantEvent.at(mnemonic).GetPrimaryEnergy(i) );
    }

    //Momentum
    mult = fGeantEvent.at(mnemonic).GetPrimaryThetaMult();
    for(int i = 0 ; i<mult ; i++){
        fNewData->SetPrimaryTheta(fGeantEvent.at(mnemonic).GetPrimaryTheta(i) ) ;
        fNewData->SetPrimaryPhi(fGeantEvent.at(mnemonic).GetPrimaryPhi(i) );
    }

    //get the energy
    double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
    TVector3 pos = fGeantEvent.at(mnemonic).GetFirstHitPosition();

    //fill the NewData Object
    fNewData->SetNewDetectorNbr(detector);
    fNewData->SetNewEnergy(energy);
    fNewData->SetPositionFirstHit(pos);

	//fNewData->Dump();

}


void RootManager::SetSceptarEvent(int eventNb, string mnemonic, int detector, int Seg)
{

    //treat
    fGeantEvent.at(mnemonic).SortPrimary();

	fSceptarData->SetEventNumber(eventNb) ;
    //get primary
    //Pdg
    int mult = fGeantEvent.at(mnemonic).GetPrimaryPdgMult();
    for( int i = 0 ; i<mult ; i++) {
        fSceptarData->SetPrimaryPdg(fGeantEvent.at(mnemonic).GetPrimaryPdg(i) ) ;
    }

    //Energy
    mult = fGeantEvent.at(mnemonic).GetPrimaryEnergyMult();
    for(int i = 0 ; i<mult ; i++){
        fSceptarData->SetPrimaryEnergy(fGeantEvent.at(mnemonic).GetPrimaryEnergy(i) );
    }

    //Momentum
    mult = fGeantEvent.at(mnemonic).GetPrimaryThetaMult();
    for(int i = 0 ; i<mult ; i++){
        fSceptarData->SetPrimaryTheta(fGeantEvent.at(mnemonic).GetPrimaryTheta(i) ) ;
        fSceptarData->SetPrimaryPhi(fGeantEvent.at(mnemonic).GetPrimaryPhi(i) );
    }

    //get the energy
    double energy = fGeantEvent.at(mnemonic).GetFullEnergy();
    TVector3 pos = fGeantEvent.at(mnemonic).GetFirstHitPosition();

    //fill the SceptarData Object
    fSceptarData->SetSceptarDetectorNbr(detector);
    fSceptarData->SetSceptarEnergy(energy);
    fSceptarData->SetPositionFirstHit(pos);

	//fSceptarData->Dump();

}

     void RootManager::ClearVariables(void){
     
		//printf("RootManager : ClearVariables .\n");    

		// clear the Fragment 
		//fFragment->ClearVariables();

		// clear the HistoryData object
		fHistoryData->ClearVariables();

		// ClearVariables the SpiceData object
		fSpiceData->ClearVariables();
		//fS3Data->ClearVariables();

		// ClearVariables the PacesData object
		//fPacesData->ClearVariables();

		// clear the NewPaces 
		//fNewData->ClearVariables();

		// clear Sceptar 
		//fSceptarData->ClearVariables();

		// ClearVariables the GriffinData object
		//fGriffinData->ClearVariables();
     }

void RootManager::SetHistory( vector <TrackInformation*> info ){

bool filled = false ; 

	for (unsigned iInfo = 0 ; iInfo < info.size() ;  iInfo ++ ) {
		
		//info.at(iInfo)->Print();
		double x, y, z ;// dummy variabales  
		vector<G4ThreeVector> trajectory ; // for the first primary charged particle, typically the emitted electron from the source  

		//primary information 		
		fHistoryData->SetHistoryPrimaryID(info.at(iInfo)->GetOriginalTrackID())		;
		fHistoryData->SetHistoryPrimaryPdg(info.at(iInfo)->GetOriginalPdg())		; 
		fHistoryData->SetHistoryPrimaryEnergy(info.at(iInfo)->GetOriginalEnergy())  ; 
		
		if (info.at(iInfo)->GetOriginalPdg()==11 && info.at(iInfo)->GetCurrentParentID() == 0 && (!filled) ){ // make sure it's a charged electron and a primary particle
			trajectory = info.at(iInfo)->GetOriginalTrajectory() ; 
			for (unsigned iTer = 0 ; iTer < trajectory.size() ; iTer++ ) {
				x = trajectory.at(iTer).getX(); 
				y = trajectory.at(iTer).getY(); 
				z = trajectory.at(iTer).getZ();
				//G4cout << iTer << " (setting) " << trajectory.at(iTer) << G4endl ; 	
				fHistoryData->SetHistoryPrimaryTrajectory(x,y,z)  ;
				//G4cin.get(); 
				if (trajectory.at(iTer).mag()> 300*mm) { /* G4cout << " maximum mag reached  " << trajectory.at(iTer).mag() << G4endl ; G4cin.get() ;*/ break; } // only store the trajectory with a sphere radius of 150 mm 
				if(iTer>999) { /*G4cout << " maximum iteration reached  " << iTer << " " << trajectory.at(iTer).mag() << G4endl ; G4cin.get() ; */ break ; } // only store 100 steps
				}
			//G4cin.get(); 
			filled = true ; 
			}
		
		x = info.at(iInfo)->GetOriginalPosition().getX(); 
		y = info.at(iInfo)->GetOriginalPosition().getY(); 
		z = info.at(iInfo)->GetOriginalPosition().getZ(); 
		fHistoryData->SetHistoryPrimaryPositionVertex(x,y,z)  ;
		 
		x = info.at(iInfo)->GetOriginalMomentum().getX(); 
		y = info.at(iInfo)->GetOriginalMomentum().getY(); 
		z = info.at(iInfo)->GetOriginalMomentum().getZ(); 
		fHistoryData->SetHistoryPrimaryMomentumVertex(x,y,z)  ; 	

        // 1s timpact of the primary 
		x = info.at(iInfo)->GetOriginalImpactMomentum().getX(); 
		y = info.at(iInfo)->GetOriginalImpactMomentum().getY(); 
		z = info.at(iInfo)->GetOriginalImpactMomentum().getZ(); 
		fHistoryData->SetHistoryPrimaryMomentum1stImpact(x,y,z)	; 
		
		x = info.at(iInfo)->GetOriginalImpactPosition().getX(); 
		y = info.at(iInfo)->GetOriginalImpactPosition().getY(); 
		z = info.at(iInfo)->GetOriginalImpactPosition().getZ(); 
		fHistoryData->SetHistoryPrimaryPosition1stImpact(x,y,z)	; 
		fHistoryData->SetHistoryPrimaryVolume1stImpact(info.at(iInfo)->GetOriginalImpactVolume())		; 
		
		// Parent of the current track hitting the target
		fHistoryData->SetHistoryParentID(info.at(iInfo)->GetCurrentParentID())	;

		// the current track hitting the target information 
		fHistoryData->SetHistoryCurrentID(info.at(iInfo)->GetCurrentTrackID())	;
		fHistoryData->SetHistoryCurrentPdg(info.at(iInfo)->GetCurrentPdg())		; 
		fHistoryData->SetHistoryCurrentEnergy(info.at(iInfo)->GetCurrentEnergyAtVertex())	; 
		fHistoryData->SetHistoryCurrentCreatorProcess(info.at(iInfo)->GetCurrentProcess())	; 
		fHistoryData->SetHistoryCurrentTime(info.at(iInfo)->GetCurrentTimeAtVertex())	; 
		fHistoryData->SetHistoryCurrentVolume1stImpact(info.at(iInfo)->GetCurrentImpactVolume())		; 
  

		x = info.at(iInfo)->GetCurrentPositionAtVertex().getX(); 
		y = info.at(iInfo)->GetCurrentPositionAtVertex().getY(); 
		z = info.at(iInfo)->GetCurrentPositionAtVertex().getZ(); 
		fHistoryData->SetHistoryCurrentPositionVertex(x,y,z)	; 
		
		x = info.at(iInfo)->GetCurrentImpactPosition().getX(); 
		y = info.at(iInfo)->GetCurrentImpactPosition().getY(); 
		z = info.at(iInfo)->GetCurrentImpactPosition().getZ(); 
		fHistoryData->SetHistoryCurrentPosition1stImpact(x,y,z) ;
		
		x = info.at(iInfo)->GetCurrentPositionAtDetector().getX(); 
		y = info.at(iInfo)->GetCurrentPositionAtDetector().getY(); 
		z = info.at(iInfo)->GetCurrentPositionAtDetector().getZ(); 
		fHistoryData->SetHistoryCurrentPositionDetector(x,y,z)	; 

		x = info.at(iInfo)->GetCurrentPositionAtDeath().getX(); 
		y = info.at(iInfo)->GetCurrentPositionAtDeath().getY(); 
		z = info.at(iInfo)->GetCurrentPositionAtDeath().getZ();  
		fHistoryData->SetHistoryCurrentPositionDeath(x,y,z)	;  



		x = info.at(iInfo)->GetCurrentMomentumAtVertex().getX(); 
		y = info.at(iInfo)->GetCurrentMomentumAtVertex().getY(); 
		z = info.at(iInfo)->GetCurrentMomentumAtVertex().getZ(); 
		fHistoryData->SetHistoryCurrentMomentumVertex(x,y,z)	; 

		x = info.at(iInfo)->GetCurrentImpactMomentum().getX(); 
		y = info.at(iInfo)->GetCurrentImpactMomentum().getY(); 
		z = info.at(iInfo)->GetCurrentImpactMomentum().getZ(); 		
		fHistoryData->SetHistoryCurrentMomentum1stImpact(x,y,z)	;
		
		x = info.at(iInfo)->GetCurrentMomentumAtDetector().getX(); 
		y = info.at(iInfo)->GetCurrentMomentumAtDetector().getY(); 
		z = info.at(iInfo)->GetCurrentMomentumAtDetector().getZ(); 
		fHistoryData->SetHistoryCurrentMomentumDetector(x,y,z)	 ; 
		
		x = info.at(iInfo)->GetCurrentMomentumAtDeath().getX(); 
		y = info.at(iInfo)->GetCurrentMomentumAtDeath().getY(); 
		z = info.at(iInfo)->GetCurrentMomentumAtDeath().getZ(); 
		fHistoryData->SetHistoryCurrentMomentumDeath(x,y,z)	; 

		}
		
	//fHistoryData->Dump();    
	//cin.get(); 
}

void RootManager::SetGriffinEvent(int Crystal)
{	
	Crystal = Crystal; 
	//Stuff goes here
	//cout << " " << endl;		

}



void RootManager::Close()
{
    fRootfile->Write();
    fRootfile->Close();    
}






