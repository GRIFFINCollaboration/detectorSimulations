use :
 ./detectionSystems.exe run.mac to run the simulation

Notes :
The code is modular, we should keep it this way, A detector <=> separate writing class 
For each new detector a "Setter" root-class must be written.
For now we have one for Spice. However Spice is not yet integrated completely, we are working on it.

To add another  "Setter" class for another detector :
- Declare the root class in the CMakeLists on top
- The files corresponding to the user-writing-class goes in the folder dataRootClass/
- Change the makefile in dataRootClass/ to compile the new files, 'make clean' and 'make' to create the libraries
- Finally, in ./src/RootManager.hh,  add the header of the new class and use it 

#include "RawG4Event.hh"
#include "path/to/dataRootClass/TSpiceData.h" 
>>> #include "path/to/dataRootClass/TNewData.h" 
( or Update the LD_LIBRARY_PATH to avoid writing the full path)
