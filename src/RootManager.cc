

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
	fS3Data = new TS3Data();

    //Paces event 
	fPacesData = new TPacesData();
	
	//GRIFFIN Event
	//fGriffinData = new TGriffinData();

	//Fragment Event
	//fFragment = new TTigFragment();
	
	//histograms 
	//fHist = new TH1F("h","h",500,0,1800);

	//Attach detector branches to the tree
	SetTree();

}

//Destructor
RootManager::~RootManager()  {}


void RootManager::SetTree(){

	//Creating the tree ;
	fOutputTree = new TTree("Simulated_Data","Simulated Data Tree");
	if(!fOutputTree) {
		cout << "\nCould not create Simulated Data Tree in root file" << endl;
		exit(-1);
		}
	fOutputTree->SetAutoSave(100000);
	
 	/*
	At this stage you can define what branches are written in the tree
	*/
		fOutputTree->Branch("S3Branch","TS3Data",&fS3Data);
	//fOutputTree->Branch("SpiceBranch","TSpiceData",&fSpiceData); 
	//----------------
	//fOutputTree->Branch("TTigFragment","TTigFragment",&fFragment, 1000, 99);
	//----------------
	fOutputTree->Branch("PacesBranch","TPacesData",&fPacesData);
	//----------------


	//----------------
	// fOutputTree->Branch("GriffinBranch","TGriffinData",&fGriffinData);
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
		//fFragment->Clear();
	
	// clear the SpiceData object
		//fSpiceData->Clear();
		//fS3Data->Clear();
	
	// clear the PacesData object
		fPacesData->Clear();

	// clear the GriffinData object
	//	fGriffinData->Clear();

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
		if (system=="PAC") {
			SetPacesEvent(eventNb, it->first, it->second.GetDetector(), it->second.GetCrystal());
			fPacesData->Dump();
			}
		//Griffin
		//Need to put Griffin key here
		//New Detector
		//if (system==XYZ) SetOtherDetectorEvent(key);
		}

 
// fill the tree 
	fOutputTree->Fill();

// clear the map
	fGeantEvent.clear();
	
}


string RootManager::BuildMnemonic(string volume, int detector, int crystal) {

	std::string system (volume, 0, 2);
	std::string sub_system (volume, 2, 1);
	std::string number = "00" ; 

	//Build the mnemonic for spice 
	//---------------------  
	if (system == "SP")
		if (sub_system == "I") {
			detector = 9 - detector ; // (detector here = ring) spice is upstream the target, according to the mnemonics outer ring is '0' inner ring is 9.
			ostringstream convert;   // stream used for the conversion	
			convert << std::setw(3) << std::setfill('0') << (detector * 12 + crystal);      //  set the width to 3 (xyz), fill the blanks by zeros
			number = convert.str(); // set 'Result' to the contents of the stream
			return ( system + sub_system  + "00XN" + number) ; 
			}
		else if (sub_system == "E") {
				return system + sub_system  +"00"+"XN"+ number ; 
				}
				
    //--------------------- 
	if (system == "PA" && sub_system == "C") {
		// (detector here = ring), Paces has only one "ring",  =>  detector = 0 
		ostringstream convert;   // stream used for the conversion	
		convert << std::setw(2) << std::setfill('0') << (crystal);      //  set the width to 2, i.e. (xy), fill the blanks by zeros = { 01, ... , 05 }
		number = convert.str(); // set 'Result' to the contents of the stream
		//cout << " RootManager::BuildMnemonic : " << system + sub_system << number <<"XN00" << "  |  " << system + sub_system + number+"XN"+"00" << endl ; 
		return system + sub_system  + number + "XN" + "00"; 
		}

    //--------------------- 
	//Griffin 	

    //--------------------- 
	//Other detectors 	

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
	fSpiceData->SetSpiceThetaEEnergy(energy) ;     
	fSpiceData->SetSpiceThetaEResEnergy(applied_resolution) ;


	// (Ph,E)
	fSpiceData->SetSpicePhiEDetectorNbr(1) ;
	fSpiceData->SetSpicePhiEStripNbr(Seg)  ;
	fSpiceData->SetSpicePhiEEnergy(energy) ;
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
	TVector3 pos = fGeantEvent.at(mnemonic).GetSecondHitPosition() ;

	// fill the S3Data object
	// (Th,E)
	fS3Data->SetS3ThetaEDetectorNbr(1) ; 
	fS3Data->SetS3ThetaEStripNbr(Ring) ;    
	fS3Data->SetS3ThetaEEnergy(energy) ;     

	// (Ph,E)
	fS3Data->SetS3PhiEDetectorNbr(1) ;
	fS3Data->SetS3PhiEStripNbr(Seg)  ;
	fS3Data->SetS3PhiEEnergy(energy) ;

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
	fPacesData->SetPacesDetectorNbr(Seg) ; 
	fPacesData->SetPacesEnergy(energy) ;    
    
    //pos.Dump();
    //getchar();
	fPacesData->SetPositionFirstHit(pos) ;	

}



void RootManager::SetGriffinEvent(int Crystal)
{	
	//Stuff goes here
	cout << " " << endl;		

}



void RootManager::Close()
{
    fRootfile->Write();
    fRootfile->Close();    
}




