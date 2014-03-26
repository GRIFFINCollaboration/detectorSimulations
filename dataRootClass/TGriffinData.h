#ifndef __GriffinDATA__
#define __GriffinDATA__
/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Mhd Moukaddam  contact address: moukaddam@triumf.ca      *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             SiLi detector. And was adapted from S1 detector Class.        *
 *             The format is the same as the one which is used for the GANIL *
 *             experiments after conversion of the raw data with GRU.        *
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

class TGriffinData : public TObject {
 private:
   // DSSD
   // Theta strips
   vector<Int_t>   fGriffin_Theta_DetNbr;
   vector<Int_t>   fGriffin_Theta_StripNbr;
   vector<Double_t>   fGriffin_Theta_Energy;

   // Phi strips
   vector<Int_t>   fGriffin_Phi_DetNbr;
   vector<Int_t>   fGriffin_Phi_StripNbr;
   vector<Double_t>   fGriffin_Phi_Energy;

   vector<TVector3> fPositionFirstHit;

   vector<Double_t> fPrimaryTheta;
   vector<Double_t> fPrimaryPhi;
   vector<Double_t> fPrimaryEnergy;
   vector<Int_t>   fPrimaryPdg;
   vector<Int_t>   fPdg;
   
   Int_t fEventNumber;
   
 public:
   TGriffinData();
   virtual ~TGriffinData();

   void   Clear();
   void   Dump() const;



   /////////////////////           GETTERS           ////////////////////////
   // (Th,E)
   Int_t   GeTGriffinThetaEMult()                 {return fGriffin_Theta_StripNbr.size();}
   Int_t   GeTGriffinThetaEDetectorNbr(Int_t i)   {return fGriffin_Theta_DetNbr.at(i);}
   Int_t   GeTGriffinThetaEStripNbr(Int_t i)      {return fGriffin_Theta_StripNbr.at(i);}
   Double_t   GeTGriffinThetaEEnergy(Int_t i)     {return fGriffin_Theta_Energy.at(i);}

   // (Ph,E)
   Int_t   GeTGriffinPhiEMult()                {return fGriffin_Phi_StripNbr.size();}
   Int_t   GeTGriffinPhiEDetectorNbr(Int_t i)  {return fGriffin_Phi_DetNbr.at(i);}
   Int_t   GeTGriffinPhiEStripNbr(Int_t i)     {return fGriffin_Phi_StripNbr.at(i);}
   Double_t   GeTGriffinPhiEEnergy(Int_t i)       {return fGriffin_Phi_Energy.at(i);}

   TVector3 GetPositionFirstHit(Int_t i)  {return fPositionFirstHit.at(i);}    
   Int_t GetEventNumber(void)             {return fEventNumber;}
   
   Double_t GetPrimaryTheta(Int_t i)    {return fPrimaryTheta.at(i);}
   Double_t GetPrimaryPhi(Int_t i)    {return fPrimaryPhi.at(i);}
   Double_t GetPrimaryEnergy(Int_t i)    {return fPrimaryEnergy.at(i);}
   Int_t   GetPrimaryPdg(Int_t i)    {return fPrimaryPdg.at(i);}
   Int_t   GetPdg(Int_t i)           {return fPdg.at(i);}
   
   /////////////////////           SETTERS           ////////////////////////
   // (Th,E)
   void   SeTGriffinThetaEDetectorNbr(Int_t det)  {fGriffin_Theta_DetNbr.push_back(det);}
   void   SeTGriffinThetaEStripNbr(Int_t Nr)      {fGriffin_Theta_StripNbr.push_back(Nr);}
   void   SeTGriffinThetaEEnergy(Double_t E)         {fGriffin_Theta_Energy.push_back(E);}

   // (Ph,E)
   void   SeTGriffinPhiEDetectorNbr(Int_t det) {fGriffin_Phi_DetNbr.push_back(det);}
   void   SeTGriffinPhiEStripNbr(Int_t Nr)     {fGriffin_Phi_StripNbr.push_back(Nr);}
   void   SeTGriffinPhiEEnergy(Double_t E)        {fGriffin_Phi_Energy.push_back(E);}
   
   void   SetPositionFirstHit(TVector3 position)  {fPositionFirstHit.push_back(position);}
   void   SetEventNumber(Int_t i)          {fEventNumber = i;}

   void   SetPrimaryTheta(Double_t theta)  {fPrimaryTheta.push_back(theta);}
   void   SetPrimaryPhi(Double_t phi)  {fPrimaryPhi.push_back(phi);}
   void   SetPrimaryEnergy(Double_t energy)       {fPrimaryEnergy.push_back(energy);}
   void   SetPrimaryPdg(Int_t pdg)       {fPrimaryPdg.push_back(pdg);}
   void   SetPdg(Int_t pdg)              {fPdg.push_back(pdg);}
     
   ClassDef(TGriffinData,1)  // GriffinData structure
};

#endif
