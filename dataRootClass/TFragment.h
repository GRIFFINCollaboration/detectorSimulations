#ifndef TFRAGMENT_H
#define TFRAGMENT_H

#include<vector>

#include<time.h>

#include<TObject.h>
#include<TNamed.h>

//#ifndef __CINT__
//#include "Globals.h"  // commented [MHD : 1 April 2014]
//#endif

//using namespace std;

class TFragment : public TNamed	{

  public:
    TFragment();
    // Create the class and reserve the vector size for MaxChannel
    // Limit the number of dynamic allocation, so better perf but more memory used
    ~TFragment(); 

    time_t MidasTimeStamp;	//->
    unsigned int MidasId;	//->
    int TriggerId;	//->
    int FragmentId;	//->
    int TriggerBitPattern;		//->

    int ChannelAddress;	//->
    int Cfd;	//->
    int Charge;	//->

    unsigned long TimeStamp;
    
    std::vector<int>  wavebuffer;	//->	
      
    virtual void	Clear(const Option_t * /* option */ =""); //!
    virtual void 	Print(const Option_t * /* option */ ="") const; //!
        
    ClassDef(TFragment,2);  // TFragment structure
};

#endif // TFRAGMENT_H
