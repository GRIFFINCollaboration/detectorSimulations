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
  detNb          = theHit.detNb;
  segNb          = theHit.segNb;
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

  G4cout  << "   det# "     << std::setw(6) << detNb
          << "   Edep "     << std::setw(7) << edep/keV << " keV"
          << "   Time "     << std::setw(7) << time/second << " sec"
          << "   Process "  << process << " "
          << "   pos ("     << std::setw(7) << pos.x()/cm
                            << std::setw(7) << pos.y()/cm
                            << std::setw(7) << pos.z()/cm << ") cm" << G4endl;
  
  G4cout.unsetf(ios::fixed);
  G4cout.precision(prec);
}

