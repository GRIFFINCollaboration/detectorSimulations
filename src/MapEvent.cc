#include "MapEvent.h"

MapEvent::MapEvent(void)
{
//cout << "Inside MapEvent::MapEvent()" << endl;
}

MapEvent::~MapEvent(void)
{
//cout << "Inside MapEvent::~MapEvent()" << endl;
}
			 

void MapEvent::FillVectors( int pdg,
							double Energy, // depositid energy
							 double Px, double Py, double Pz, // position vector
							 
							 int ID, // original(primary) TrackID
							 int PrimPdg,// primary particle pdg        PDG encoding
							 double PrimEnergy,//original(primary) energy
							 double Mx, double My, double Mz) // primary particle momentum vector
							
							{
							    fHCPdg.push_back(pdg); 
								fHCEnergy.push_back(Energy) ;
								fHCPosition.push_back(TVector3(Px,Py,Pz)) ;
								
								fHCPrimaryID.push_back(ID) ;
								fHCPrimaryPdg.push_back(PrimPdg) ;
								fHCPrimaryEnergy.push_back(PrimEnergy) ;
								fHCPrimaryMoment.push_back(TVector3(Mx,My,Mz)) ;
							}




void MapEvent::Treat(void)
{

// treat map of Pdg  < Key = Primary Pdg,  Value = multiplicity hits in the collection >
	for(unsigned i=0 ; i < fHCPrimaryPdg.size() ; i++)
	{
		fMapPrimaryPdg[fHCPrimaryPdg.at(i)]++;
		fMapPrimaryEnergy[fHCPrimaryEnergy.at(i)]++;
		//std::cout <<i<< " : " << fHCPrimaryPdg.at(i) << '\n';
	}


}
	


TVector3 MapEvent::GetFirstHitPosition(void) 
	{
	return fHCPosition.at(0) ; // returns the first position registered in this vector
	}


Double_t MapEvent::GetFullEnergy(void) 
	{

	Double_t energy=0 ;
	
	for(unsigned i=0 ; i< fHCEnergy.size() ; i++)
	 energy += fHCEnergy.at(i);

	return energy ;
	}



Int_t MapEvent::GetHCPrimaryPdg(int i) 
	{
	return fHCPrimaryPdg.at(i);
	}
	
	
Double_t MapEvent::GetHCPrimaryEnergy(int i) 
	{
	return  fHCPrimaryEnergy.at(i);
	}
	
	
	
Int_t MapEvent::GetPrimaryPdgForID(int ID) 
	{

	int pdg = -1 ;
	
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { pdg = fHCPrimaryPdg.at(i); break; }

	return pdg ;

	}
	
	
Int_t MapEvent::GetPrimaryPdgMult(void) 
	{
	return (int) fMapPrimaryPdg.size() ;	
	}
	
Int_t MapEvent::GetPrimaryEnergyMult(void) 
	{
	return (int) fMapPrimaryEnergy.size() ;	// Pdg Map 
	}
	
	
Int_t MapEvent::GetPrimaryPdg(int i)  // This is always used  with the multiplicity GetPrimaryPdgMult 
	{
	
	std::map<int,int>::iterator it ; // create iterator
	it=fMapPrimaryPdg.begin();  // set it on the begining of the map
	advance(it,i); // advance it to the argument position
    //std::cout << it->first << " => " << it->second << '\n';
	return it->first ;
	
	}
	
	
Double_t MapEvent::GetPrimaryEnergy(int i)  // This is always used  with the multiplicity GetPrimaryPdgMult 
	{

	map<double,int>::iterator it ; // create iterator
	it=fMapPrimaryEnergy.begin();  // set it on the begining of the map
	advance(it,i); // advance it to the argument position
    //cout << " i : "<< i << " - "<< it->first << " => " << it->second << '\n';
	return it->first ;
	
	}

	

Double_t MapEvent::GetPrimaryEnergyForID(int ID) 
	{

	int energy = -1 ;
	
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { energy = fHCPrimaryEnergy.at(i); break; }

	return energy ;

	}


double MapEvent::GetDetectedEnergyForID(int ID) 
	{

	double energy = 0 ;
	
	for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
	if (fHCPrimaryID.at(i) ==  ID ) { energy = energy + fHCPrimaryEnergy.at(i);}

	if (energy == 0) energy = -1 ;
	return energy ;

	}


double MapEvent::GetAngleOfEmissionForID(int ID) 
{

double angle = -1 ;
TVector3 beamDirection(0,0,1) ;

for(unsigned i=0 ; i< fHCPrimaryID.size() ; i++)
if (fHCPrimaryID.at(i) ==  ID ) { angle = fHCPrimaryMoment.at(i).Angle(beamDirection); }

return angle ;

}


void MapEvent::ShowVectorsContent(void)
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
for(unsigned i=0 ; i< fHCPrimaryMoment.size() ; i++)
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
for(unsigned i=0 ; i< fHCPrimaryMoment.size() ; i++)
	{
	cout<<" - X "<<fHCPrimaryMoment.at(i).X();
	cout<<" - Y "<<fHCPrimaryMoment.at(i).Y();
	cout<<" - Z "<<fHCPrimaryMoment.at(i).Z()<<endl;
	}

}



