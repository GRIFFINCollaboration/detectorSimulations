

#include "RootManager.hh"


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
	fGriffinData = new TGriffinData();
<<<<<<< HEAD
	fS3Data = new TS3Data();

=======
	fDetectorSpice = new DetectionSystemSpice();
	
>>>>>>> griffin/master
	//histograms 
	fHist = new TH1F("h","h",500,0,1800);


	//Creating the tree ;
	fOutputTree = new TTree("Simulated_Data","Simulated Data Tree");
	if(!fOutputTree) {
		cout << "\nCould not create Simulated Data Tree in root file" << endl;
		exit(-1);
		}
	fOutputTree->SetAutoSave(100000);

	//Attach detector branches to the tree
	SetTree();

}

//Destructor
RootManager::~RootManager()  {}


void RootManager::SetTree()
{
 fOutputTree->Branch("SpiceBranch","TSpiceData",&fSpiceData);
 fOutputTree->Branch("S3Branch","TS3Data",&fS3Data);
 
 /*
Other detector branches goes here
*/

} 

void RootManager::FillHist(double temp)
{
    float tempf = (float)temp; 
    fHist->Fill(tempf);
}


void RootManager::FillG4Hit(int key,  						// integer representing the key, could be used to identify the detector
							int pdg,			 				// integer representing the type of particle
							double Energy, 						// depositid energy in a step
							double Px, double Py, double Pz, 	// position vector
							int Id, 							// original(primary) ID
							int PrimPdg,						// primary particle definition        PDG encoding
							double PrimEnergy,					// original(primary) energy
							double Mx, double My, double Mz)	// primary particle momentum vector
							{
							/*			
							Some other functions 
							*/
							//Fill the event, calculate the full deposited energy etc...
							fGeantEvent[key].FillVectors( pdg, Energy,  Px, Py, Pz,  Id,  PrimPdg, PrimEnergy,  Mx, My, Mz);
							}



void RootManager::SortEvent(void)
{

// Sort Data from the map second element by getters and set them
       std::map<Int_t,RawG4Event>::iterator it;
  for (std::map<Int_t,RawG4Event>::iterator it=fGeantEvent.begin(); it!=fGeantEvent.end(); ++it) {
  
			int key = it->first ;

			//Spice
			if (key>1 && key <50) SetSpiceEvent(key);
			if (key>50 && key <1000) SetS3Event(key);

			//other
			//if (key>1 && key <1000) SetOtherDetectorEvent(key);
			}

// fill the tree by SpiceData
	fOutputTree->Fill();
 
// clear the map
	fGeantEvent.clear();

// clear the SpiceData object
	fSpiceData->Clear();
	fS3Data->Clear();
}


void RootManager::SetSpiceEvent(int RingSeg)
{	
			// treat
			fGeantEvent.at(RingSeg).SortPrimary();			

			// get the segment and ring
			int Seg = (RingSeg%100) ;   //?? 100 should be 12... CHECK, MHD 20Dec2013 
			int Ring = (RingSeg-Seg)/100;

			// get primary
			// Pdg	
			int mult = fGeantEvent.at(RingSeg).GetPrimaryPdgMult(); // inside this particular pad 
		    for (int i = 0 ; i<mult ;  i++ )
			fSpiceData->SetPrimaryPdg( fGeantEvent.at(RingSeg).GetPrimaryPdg(i) ) ;

		  	// Energy      
		    mult = fGeantEvent.at(RingSeg).GetPrimaryEnergyMult(); // this should be the same as above
			for (int i = 0 ; i<mult ;  i++ )
			{
			fSpiceData->SetPrimaryEnergy( fGeantEvent.at(RingSeg).GetPrimaryEnergy(i) ) ;
			}

			// Momentum
			    mult = fGeantEvent.at(RingSeg).GetPrimaryThetaMult(); // this should be the same as above
			for (int i = 0 ; i<mult ;  i++ )
			{
			fSpiceData->SetPrimaryTheta(fGeantEvent.at(RingSeg).GetPrimaryTheta(i) ) ;
			fSpiceData->SetPrimaryPhi( fGeantEvent.at(RingSeg).GetPrimaryPhi(i) ) ;
			}

			// get the energy 			
			double energy = fGeantEvent.at(RingSeg).GetFullEnergy();
			double applied_resolution = fDetectorSpice->ApplySpiceResolution(energy);
			TVector3 pos = fGeantEvent.at(RingSeg).GetFirstHitPosition() ;

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
<<<<<<< HEAD

			fSpiceData->SetPositionFirstHit(pos) ;				
}


void RootManager::SetS3Event(int RingSeg)
{	
			// treat
			fGeantEvent.at(RingSeg).SortPrimary();			

			// get the segment and ring
			int Seg = (RingSeg%100) ;   //?? 100 should be 12... CHECK, MHD 20Dec2013 
			int Ring = (RingSeg-Seg)/100;

			// get primary
			// Pdg	
			int mult = fGeantEvent.at(RingSeg).GetPrimaryPdgMult(); // inside this particular pad 
		    for (int i = 0 ; i<mult ;  i++ )
			fS3Data->SetPrimaryPdg( fGeantEvent.at(RingSeg).GetPrimaryPdg(i) ) ;

		  	// Energy      
		    mult = fGeantEvent.at(RingSeg).GetPrimaryEnergyMult(); // this should be the same as above
			for (int i = 0 ; i<mult ;  i++ )
			{
			fS3Data->SetPrimaryEnergy( fGeantEvent.at(RingSeg).GetPrimaryEnergy(i) ) ;
			}

			// Momentum
			    mult = fGeantEvent.at(RingSeg).GetPrimaryThetaMult(); // this should be the same as above
			for (int i = 0 ; i<mult ;  i++ )
			{
			fS3Data->SetPrimaryTheta(fGeantEvent.at(RingSeg).GetPrimaryTheta(i) ) ;
			fS3Data->SetPrimaryPhi( fGeantEvent.at(RingSeg).GetPrimaryPhi(i) ) ;
			}

			// get the energy 			
			double energy = fGeantEvent.at(RingSeg).GetFullEnergy();
			TVector3 pos = fGeantEvent.at(RingSeg).GetFirstHitPosition() ;

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

=======
			fSpiceData->SetSpicePhiEResEnergy(applied_resolution) ;
			
			fSpiceData->SetPositionFirstHit(pos) ;		
		
>>>>>>> griffin/master
}



// set event number
void RootManager::SetEventNumber(int i )
{
fSpiceData->SetEventNumber(i) ;
fS3Data->SetEventNumber(i) ;
}		


void RootManager::Close()
{
    fRootfile->Write();
    fRootfile->Close();    
}


