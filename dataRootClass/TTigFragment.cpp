#include "TTigFragment.h"

ClassImp(TTigFragment);

TTigFragment::TTigFragment()
{
   Clear();
}

TTigFragment::~TTigFragment()
{
   //Clear();
}

void TTigFragment::Clear(Option_t *)
{
   TFragment::Clear();

   ODBType = 0;

   ChannelName = "XXXXXXXXXX";
   ChannelNumber = 0;

   Led = -1;
   TimeToTrig = -9999;

   ChargeCal = 0.0;

   TimeStampLow = -1;
   TimeStampHigh = -1;
   TimeStampLive = 0;
   TimeStampTR = 0;             // triggers requested
   TimeStampTA = 0;             // triggers accepted

   //waveform.SetBit(TH1::kCanRebin);

   SlowRiseTime = false;
   PileUp = false;

}

void TTigFragment::Print(Option_t *) const
{
   printf("%s Event at	%i:\n", DigitizerType.c_str(), MidasId);
   printf("MidasId    	%i\n", MidasId);
   printf("TriggerId: 	%i\n", TriggerId);
   printf("FragmentId:   %i\n", FragmentId);
   printf("TriggerBit:	0x%08x\n", TriggerBitPattern);
   printf("Channel: %i\tName: %s\n", ChannelNumber, ChannelName.c_str());
   printf("\tChannel Address: 0x%08x", ChannelAddress);
   printf("\t\tChannel Num:      %i\n", ChannelNumber);
   printf("\tCharge		   0x%08x", Charge);
   printf("\t\tEnergy:		      %f\n", ChargeCal);
   printf("\tLED:         0x%08x", Led);
   printf("\t\tTimeStamp High: 0x%08x\n", TimeStampHigh);
   printf("\tCFD:         0x%08x", Cfd);
   printf("\t\tTimeStamp Low:    0x%08x\n", TimeStampLow);
   printf("\tTimeToTrig:  %i\n", TimeToTrig);
   unsigned short temptime = (TimeStampLow & 0x0000ffff) - ((Cfd >> 4) & 0x0000ffff);   //TimeStampLow&0x0000ffff; 
   printf("\ttime from timestamp(to the nearest 10ns):    0x%04x\t%ins\n", temptime, temptime * 10);
   if (!wavebuffer.empty())
      printf("Has a wave form stored.\n");
   else
      printf("Does Not have a wave form stored.\n");
}


void TTigFragment::SetTimeStamp()
{
   TimeStamp = (((Long64_t) TimeStampHigh) << 24) | (TimeStampLow & 0x00ffffff); // [MHD : 1 April 2014 ] replaced int64_t by Long64_t
}


int TTigFragment::Compare(const TObject * obj) const
{
   TTigFragment *other = (TTigFragment *) obj;
   if (TriggerId > other->TriggerId) {
      return -1;
   } else if (TriggerId == other->TriggerId) {
      if (FragmentId > other->FragmentId) {
         return -1;
      } else if (FragmentId == other->FragmentId) {
         return 0;              // this shouldn't happen
      } else {
         return 1;
      }
   } else {
      return 1;
   }
}
