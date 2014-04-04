#include "RawG4Event.hh"

RawG4Event::RawG4Event(void){
}

RawG4Event::~RawG4Event(void){
}
			 

void RawG4Event::FillVectors( 
							int pdg,
							double Energy, // depositid energy
							int detector, // detector 
							int crystal, // crystal 
							double Px, double Py, double Pz, // position vector

							int ID, // original(primary) TrackID
							int PrimPdg,// primary particle pdg        PDG encoding
							double PrimEnergy,//original(primary) energy
							double Mx, double My, double Mz) { // primary particle momentum vector
							
							fHCPdg.push_back(pdg); 
							fHCEnergy.push_back(Energy) ;
							fHCPosition.push_back(TVector3(Px,Py,Pz)) ;
							fHCDetector.push_back(detector) ;
							fHCCrystal.push_back(crystal) ;
							
							
							fHCPrimaryID.push_back(ID) ;
							fHCPrimaryPdg.push_back(PrimPdg) ;
							fHCPrimaryEnergy.push_back(PrimEnergy) ;
							fHCPrimaryMomentum.push_back(TVector3(Mx,My,Mz)) ;
							fHCPrimaryTheta.push_back( TVector3(Mx,My,Mz).Angle(TVector3(0,0,1)) / deg ); // TVector3(0,0,1) is beam direction  in degrees
							fHCPrimaryPhi.push_back(TVector3(Mx,My,Mz).Phi()  / deg );
							}



void RawG4Event::SortPrimary(void) {

	for(unsigned i=0 ; i < fHCPrimaryPdg.size() ; i++) 	{
		fMapPrimaryPdg[fHCPrimaryPdg.at(i)]++;
		fMapPrimaryEnergy[fHCPrimaryEnergy.at(i)]++;
       	fMapPrimaryTheta[fHCPrimaryTheta.at(i)]++; 
		fMapPrimaryPhi[fHCPrimaryPhi.at(i)]++;
	}

}
	

// ------------------- GETTERS --------------------


TVector3 RawG4Event::GetFirstHitPosition(void) 	{  // generates artificially a spurious peak in the interstrip region 
	return fHCPosition.at(0) ; // returns the first position registered in this vector
	}

TVector3 RawG4Event::GetSecondHitPosition(void)	{
	if (fHCPosition.size()>1) return fHCPosition.at(1) ; // returns the second position registered in this vector if the vector size > 1 
	else return fHCPosition.at(0) ; // else returns the first position 
	}
	
	
Int_t RawG4Event::GetDetector(void) {
	return fHCDetector.at(0) ; // returns the detector number, since the key of the map is the mnemonic, all the values should be the same 
	}
	
Int_t RawG4Event::GetCrystal(void) 	{
	return fHCCrystal.at(0) ; // returns the crystal number, since the key of the map is the mnemonic, all the values should be the same 
	}
	
Double_t RawG4Event::GetFullEnergy(void) {

	Double_t energy=0 ;
	for(unsigned i=0 ; i< fHCEnergy.size() ; i++)	 energy += fHCEnergy.at(i);

	return energy ;
	}


// Get values for an element of the hit collection
Int_t RawG4Event::GetHCPrimaryPdg(int i) 
	{
	return fHCPrimaryPdg.at(i);
	}
		
Double_t RawG4Event::GetHCPrimaryEnergy(int i) 
	{
	return  fHCPrimaryEnergy.at(i);
	}

TVector3 RawG4Event::GetHCPrimaryMomentum(int i) 
	{
	return  fHCPrimaryMomentum.at(i);
	}

Double_t RawG4Event::GetHCPrimaryTheta(int i) 
	{
	return  fHCPrimaryTheta.at(i);
	}
	
Double_t RawG4Event::GetHCPrimaryPhi(int i) 
	{
	return  fHCPrimaryPhi.at(i);
	}

				
Int_t RawG4Event::GetPrimaryPdgForID(int ID) 
	{
	int pdg = -1 ;
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { pdg = fHCPrimaryPdg.at(i); break; }
	return pdg ;
	}
	

// Get the multiplicities 	
Int_t RawG4Event::GetPrimaryPdgMult(void) 
	{
	return (int) fMapPrimaryPdg.size() ;	
	}
	
Int_t RawG4Event::GetPrimaryEnergyMult(void) 
	{
	return (int) fMapPrimaryEnergy.size() ;	// Pdg Map 
	}
	
Int_t RawG4Event::GetPrimaryThetaMult(void)  // this will dictate the multiplicity theta and Phi of the primary particle
	{
	return (int) fMapPrimaryTheta.size() ;	
	}	

//Get the specific values of an element in a collection hit 
Int_t RawG4Event::GetPrimaryPdg(int i)  // This is always used  with the multiplicity GetPrimaryPdgMult 
	{
	std::map<int,int>::iterator it ; // create iterator
	it=fMapPrimaryPdg.begin();  // set it on the begining of the map
	advance(it,i); // advance "it" to the argument position
    //std::cout << it->first << " => " << it->second << '\n';
	return it->first ;
	}	
Double_t RawG4Event::GetPrimaryEnergy(int i)  
	{
	map<double,int>::iterator it ; 
	it=fMapPrimaryEnergy.begin(); 
	advance(it,i); 
	return it->first ;
	}	
Double_t RawG4Event::GetPrimaryTheta(int i)  
	{
	map<Double_t,int>::iterator it ; 
	it=fMapPrimaryTheta.begin(); 
	advance(it,i); 
	return it->first ;
	}	
Double_t RawG4Event::GetPrimaryPhi(int i)  
	{
	map<Double_t,int>::iterator it ; 
	it=fMapPrimaryPhi.begin(); 
	advance(it,i);
	return it->first ;
	}
	

Double_t RawG4Event::GetPrimaryEnergyForID(int ID) 
	{

	int energy = -1 ;
	
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { energy = fHCPrimaryEnergy.at(i); break; }

	return energy ;

	}


double RawG4Event::GetDetectedEnergyForID(int ID) 
	{

	double energy = 0 ;
	
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { energy = energy + fHCPrimaryEnergy.at(i);}

	if (energy == 0) energy = -1 ;
	return energy ;

	}

double RawG4Event::GetAngleOfEmissionForID(int ID) 
{

double angle = -1 ;

for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
if (fHCPrimaryID.at(i) ==  ID ) { angle = fHCPrimaryTheta.at(i); }

return angle ;

}


void RawG4Event::ShowVectorsContent(void)
{

cout<<" * * * * * * * * * * *  "<<endl;
cout<<" Particle Informations "<<endl;
cout<<" * * * * * * * * * * *  "<<endl;

cout<<" Particle Definition (PDG)"<<endl;
for(unsigned i=0 ; i< fHCPdg.size() ; i++)
cout<<" "<<fHCPdg.at(i)<<endl;

cout<<" Particle Deposited Energy "<<endl;
for(unsigned i=0 ; i< fHCEnergy.size() ; i++)
cout<<" "<<fHCEnergy.at(i)<<endl;

cout<<" Particle Position  "<<endl;
for(unsigned i=0 ; i< fHCPosition.size() ; i++)
	{
	cout<<" - X "<<fHCPosition.at(i).X();
	cout<<" - Y "<<fHCPosition.at(i).Y();
	cout<<" - Z "<<fHCPosition.at(i).Z()<<endl;
	}
	
cout<<" * * * * * * * * * * *  "<<endl;
cout<<" Primary Informations "<<endl;
cout<<" * * * * * * * * * * *  "<<endl;

cout<<" Particle ID (// time of emission) "<<endl;
for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
cout<<" "<<fHCPrimaryID.at(i)<<endl;

cout<<" Particle Definition (PDG)"<<endl;
for(unsigned i=0 ; i< fHCPrimaryPdg.size() ; i++)
cout<<" "<<fHCPrimaryPdg.at(i)<<endl;

cout<<" Particle Full Energy  "<<endl;
for(unsigned i=0 ; i< fHCPrimaryEnergy.size() ; i++)
cout<<" "<<fHCPrimaryEnergy.at(i)<<endl;

cout<<" Particle Moment  "<<endl;
for(unsigned i=0 ; i< fHCPrimaryMomentum.size() ; i++)
	{
	cout<<" - Theta "<<fHCPrimaryTheta.at(i); // angle with beam 
	cout<<" - Phi "<<fHCPrimaryPhi.at(i);		
	cout<<" - X "<<fHCPrimaryMomentum.at(i).X();
	cout<<" - Y "<<fHCPrimaryMomentum.at(i).Y();
	cout<<" - Z "<<fHCPrimaryMomentum.at(i).Z()<<endl;
	}

}



