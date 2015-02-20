#ifndef TTIGFRAGMENT_H
#define TTIGFRAGMENT_H

#include<string>

#include "TFragment.h"

class TTigFragment:public TFragment {
   // nothing yet
 public:
   TTigFragment(void);
   ~TTigFragment();
   
   std::string DigitizerType;  //->      // currently tig10 or tig64
   int ODBType;                 //->

   std::string ChannelName;    //->
   int ChannelNumber;           //->

   int Led;                     //->
   int TimeToTrig;              //->

   float ChargeCal;             //->

   int TimeStampLow;            //->
   int TimeStampHigh;           //->
   int TimeStampLive;           //->
   int TimeStampTR;             //->       // triggers requested
   int TimeStampTA;             //->         // triggers accepted

   bool SlowRiseTime;           //->
   bool PileUp;                 //->


   void Clear(const Option_t * /* option */="");                //!
   void Print(const Option_t * /* option */="") const;                //!

   void SetTimeStamp();         //!

   //in order to use root sort function when the fragments are in arrays 
   //puts newset id at the beginning of the array, older ids at the end.
   bool IsSortable() const {
      return kTRUE;
   };                           //!
   int Compare(const TObject * obj) const;      //! compare by trigger and fragment id 

    ClassDef(TTigFragment, 13); // TTigFragment structure
};

#endif                          // TTIGFRAGMENT_H
