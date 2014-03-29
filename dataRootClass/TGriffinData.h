#ifndef __GriffinDATA__
#define __GriffinDATA__
/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 * Modified by    : Ryan Dunlop  contact address: rdunlop@uoguelph.ca        *
 *---------------------------------------------------------------------------*
 * Decription: This class stores the results of the G4 simulation for the    *
 *             GRIFFIN detector. And was adapted from S1 detector Class.     *
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

        vector<Int_t>   fDetNbr;
        vector<TVector3> fPositionFirstHit;

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
        Int_t   GetGriffinEDetectorNbr(Int_t i)   {return fDetNbr.at(i);}
        Int_t   GetGriffinEMult(Int_t i)        {return fDetNbr.size();}

        TVector3 GetPositionFirstHit(Int_t i)  {return fPositionFirstHit.at(i);}    
        Int_t GetEventNumber(void)             {return fEventNumber;}

        Double_t GetPrimaryEnergy(Int_t i)    {return fPrimaryEnergy.at(i);}
        Int_t   GetPrimaryPdg(Int_t i)    {return fPrimaryPdg.at(i);}
        Int_t   GetPdg(Int_t i)           {return fPdg.at(i);}

        /////////////////////           SETTERS           ////////////////////////
        void   SetGriffinEDetectorNbr(Int_t det)  {fDetNbr.push_back(det);}

        void   SetPositionFirstHit(TVector3 position)  {fPositionFirstHit.push_back(position);}
        void   SetEventNumber(Int_t i)          {fEventNumber = i;}

        void   SetPrimaryEnergy(Double_t energy)       {fPrimaryEnergy.push_back(energy);}
        void   SetPrimaryPdg(Int_t pdg)         {fPrimaryPdg.push_back(pdg);}
        void   SetPdg(Int_t pdg)              {fPdg.push_back(pdg);}
     
   ClassDef(TGriffinData,1)  // GriffinData structure
};

#endif
        
