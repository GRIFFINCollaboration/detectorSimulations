
#ifndef _NEWDATA_
#define _NEWDATA_


//c++
#include <vector>

//ROOT
#include "TObject.h"
#include "TVector3.h"

class TNewData : public TObject 
{

private:
    
    vector<Int_t> fNew_DetNbr;
    vector<Double_t> fNew_Energy;

    vector<TVector3> fPositionFirstHit;

    vector<Double_t> fPrimaryTheta; 
    vector<Double_t> fPrimaryPhi;
    vector<Double_t> fPrimaryEnergy;
    vector<Double_t> fPrimaryPdg;
    vector<Double_t> fPdg;

    Int_t fEventNumber;

public:

    TNewData();
    virtual ~TNewData();

    void ClearVariables();
    void Dump() const;

    /////////////////////             GETTERS             ////////////////////////
    Int_t       GetNewMult()                            {return fNew_DetNbr.size();}
    Int_t       GetNewDetectorNbr(Int_t i)              {return fNew_DetNbr.at(i);}
    Double_t    GetNewEnergy(Int_t i)                   {return fNew_Energy.at(i);}


    TVector3    GetPositionFirstHit(Int_t i)            {return fPositionFirstHit.at(i);}
    Int_t       GetEventNumber(void)                    {return fEventNumber;}

    Double_t    GetPrimaryTheta(Int_t i)                {return fPrimaryTheta.at(i);}
    Double_t    GetPrimaryPhi(Int_t i)                  {return fPrimaryPhi.at(i);}
    Double_t    GetPrimaryEnergy(Int_t i)               {return fPrimaryEnergy.at(i);}
    Int_t       GetPrimaryPdg(Int_t i)                  {return fPrimaryPdg.at(i);}
    Int_t       GetPdg(Int_t i)                         {return fPdg.at(i);}

    /////////////////////           SETTERS                  ////////////////////////
    void        SetNewDetectorNbr(Int_t det)            {fNew_DetNbr.push_back(det);}
    void        SetNewEnergy(Double_t E)                {fNew_Energy.push_back(E);}

    void        SetPrimaryTheta(Double_t theta)         {fPrimaryTheta.push_back(theta);}
    void        SetPrimaryPhi(Double_t phi)             {fPrimaryPhi.push_back(phi);}
    void        SetPrimaryEnergy(Double_t energy)       {fPrimaryEnergy.push_back(energy);}
    void        SetPrimaryPdg(Int_t pdg)                {fPrimaryPdg.push_back(pdg);}
    void        SetPdg(Int_t pdg)                       {fPdg.push_back(pdg);}

    void        SetPositionFirstHit(TVector3 position)  {fPositionFirstHit.push_back(position);}
    void        SetEventNumber(Int_t i)                 {fEventNumber = i;}

    ClassDef(TNewData,1)

};

#endif






