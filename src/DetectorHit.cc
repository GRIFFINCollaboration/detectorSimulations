#include "DetectorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include <cstdio>
#include <iomanip>

using namespace std;

G4Allocator<DetectorHit> DetectorHitAllocator;

DetectorHit::DetectorHit() {}

DetectorHit::~DetectorHit() {}

DetectorHit::DetectorHit(const DetectorHit& theHit) : G4VHit()
{

	OriginID       		= theHit.OriginID;
  OriginEnergy          = theHit.OriginEnergy;
  OriginPdg          	= theHit.OriginPdg;
  OriginMoment          = theHit.OriginMoment; 		
  
  AncestorsBirthVolume  = theHit.AncestorsBirthVolume;  // MHD 01 May 2013
  AncestorsDeathVolume  = theHit.AncestorsDeathVolume;  // MHD 01 May 2013
  AncestorsPdg          = theHit.AncestorsPdg;          // MHD 01 May 2013
  
  detNb          = theHit.detNb;
  segNb          = theHit.segNb;
  pdg            = theHit.pdg;
  edep           = theHit.edep;
  time           = theHit.time;
  process        = theHit.process;
  collection     = theHit.collection;
  pos            = theHit.pos;
}

void DetectorHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(1.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,1.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void DetectorHit::Print()
{

  G4int prec = G4cout.precision(2);
  G4cout.setf(ios::fixed);

  G4cout  << G4endl << "================ Print Detector Hit ================== " << G4endl ;
  G4cout  << G4endl << "---------------- Primary Ancestor particle  ---------------- " << G4endl ;
  G4cout  << " OriginID "     		<< std::setw(6) << OriginID
  		  << " OriginEnergy "     	<< std::setw(6) << OriginEnergy
  		  << " OriginPdg "     	    << std::setw(6) << OriginPdg 
  		  << " OriginMoment ( "     << std::setw(7) << OriginMoment.x()
                            		<< std::setw(7) << OriginMoment.y()
                            		<< std::setw(7) << OriginMoment.z() << ")" << G4endl;
                                   
  G4cout  << G4endl << "---------------- Secondary Ancestors particles  ---------------- " << G4endl ;
  
  G4cout  <<G4endl<< " >>>> Birth volumes "<<G4endl;                     
  for (unsigned i = 0 ; i < AncestorsBirthVolume.size(); i++ )                
  G4cout  << i<<"-"<< AncestorsBirthVolume.at(i)<< "     " ;
  G4cout  << G4endl;
  
  G4cout  <<G4endl<< " >>>> Death volumes "<<G4endl;           
  for (unsigned i = 0 ; i < AncestorsDeathVolume.size(); i++ )                    
  G4cout  << i<<"-"<< AncestorsDeathVolume.at(i)<< "     " ;
  G4cout  << G4endl;
  
  G4cout  <<G4endl<< " >>>> Pdg "<<G4endl;   		  
  for (unsigned i = 0 ; i < AncestorsPdg.size(); i++ )   
  G4cout  << i<<"-"<< AncestorsPdg.at(i)<< "     " ;
  G4cout  << G4endl;
                
                   
  G4cout  << G4endl << "---------------- grand daughter particle ---------------- " << G4endl ;                            		                            		   		   		     
  G4cout  << "   det# "     		<< std::setw(6) << detNb
          << "   seg# "     		<< std::setw(6) << segNb
          << "   Pdg  " 		    << std::setw(6) << pdg
          << "   Edep "     		<< std::setw(7) << edep/keV << " keV"
          << "   Time "     		<< std::setw(7) << time/second << " sec"
          << "   Process "  		<< process << " "
          << "   collection "  		<< collection << " "
          << "   pos ("     		<< std::setw(7) << pos.x()/cm
                            		<< std::setw(7) << pos.y()/cm
                            		<< std::setw(7) << pos.z()/cm << ") cm" << G4endl;
  
  G4cout.unsetf(ios::fixed);
  G4cout.precision(prec);
}

