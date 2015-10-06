#ifndef _PACESDATA_
#define _PACESDATA_
/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Mhd Moukaddam  contact address: moukaddam@triumf.ca      *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             PACES detector. And was adapted from S1 detector Class.       *
 *             This class derives from TObject (ROOT) and its aim is to be   *
 *             stored in the output TTree of the G4 simulation               *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/


//c++
#include <vector>

//ROOT
#include "TObject.h"
#include "TVector3.h"

class TPacesData : public TObject {
 private:
   // SiLi
   vector<Int_t>   fPaces_DetNbr;
   vector<Double_t>   fPaces_Energy;

   vector<TVector3> fPositionFirstHit;

   vector<Double_t> fPrimaryTheta;
   vector<Double_t> fPrimaryPhi;
   vector<Double_t> fPrimaryEnergy;
   vector<Int_t>   fPrimaryPdg;
   vector<Int_t>   fPdg;
   
   Int_t fEventNumber;
   
 public:
   TPacesData();
   virtual ~TPacesData();

   void   ClearVariables();
   void   Dump() const;

   /////////////////////           GETTERS           ////////////////////////
   Int_t   GetPacesMult()                 {return fPaces_DetNbr.size();}
   Int_t   GetPacesDetectorNbr(Int_t i)   {return fPaces_DetNbr.at(i);}
   Double_t   GetPacesEnergy(Int_t i)     {return fPaces_Energy.at(i);}


   TVector3 GetPositionFirstHit(Int_t i)  {return fPositionFirstHit.at(i);}    
   Int_t GetEventNumber(void)             {return fEventNumber;}
   
   Double_t GetPrimaryTheta(Int_t i)    {return fPrimaryTheta.at(i);}
   Double_t GetPrimaryPhi(Int_t i)    {return fPrimaryPhi.at(i);}
   Double_t GetPrimaryEnergy(Int_t i)    {return fPrimaryEnergy.at(i);}
   Int_t   GetPrimaryPdg(Int_t i)    {return fPrimaryPdg.at(i);}
   Int_t   GetPdg(Int_t i)           {return fPdg.at(i);}
   
   /////////////////////           SETTERS           ////////////////////////
   void   SetPacesDetectorNbr(Int_t det)  {fPaces_DetNbr.push_back(det);}
   void   SetPacesEnergy(Double_t E)         {fPaces_Energy.push_back(E);}
   
   void   SetPositionFirstHit(TVector3 position)  {fPositionFirstHit.push_back(position);}
   void   SetEventNumber(Int_t i)          {fEventNumber = i;}

   void   SetPrimaryTheta(Double_t theta)  {fPrimaryTheta.push_back(theta);}
   void   SetPrimaryPhi(Double_t phi)  {fPrimaryPhi.push_back(phi);}
   void   SetPrimaryEnergy(Double_t energy)       {fPrimaryEnergy.push_back(energy);}
   void   SetPrimaryPdg(Int_t pdg)       {fPrimaryPdg.push_back(pdg);}
   void   SetPdg(Int_t pdg)              {fPdg.push_back(pdg);}
     
   ClassDef(TPacesData,1)  // PacesData structure
};

#endif
