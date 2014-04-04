#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)


#----------------------------------------------------------------------------
# Adding a variable to use root through out the G4 project. 
#
set(mytarget    			detectionSystems)  			# equivalent to "name := myexample" in G4 GNUmakefiles
set(useROOT     			true            ) 		 	# use true or false (or comment to set to false)
set(myROOTclassDir 			dataRootClass	) 			# directory of the root classes
set(myROOTclass 			true			)   		# comment if none (please see http://root.cern.ch/phpBB3//viewtopic.php?t=6172)
set(myROOTS3class 			TS3Data.cpp		)   		# comment if none 
set(myROOTFragmentclass 	TFragment.cpp	) 			# comment if none 
set(myROOTTigFragmentclass 	TTigFragment.cpp) 			# comment if none  
set(myROOTSpiceclass 		TSpiceData.cpp	)   		# comment if none  
set(myROOTGriffclass 		TGriffinData.cpp)  			# comment if none
#set(myROOTNewclass 		TNewData.cpp) 				# comment if none  


# http://www.cmake.org/cmake/help/cmake_tutorial.html
# http://www.cmake.org/cmake/help/cmake2.6docs.html
# http://www.cmake.org/Wiki/CMake_FAQ

set(mymain ${mytarget}.cc)
set(myexe  ${mytarget})
set(myproj ${mytarget}.proj)  # any name as long as it is different from myexe

project($myproj)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../shared/include 
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Use root in the project 
#                   
if(useROOT)
	if(myROOTclass)
	
		if(myROOTFragmentclass)
		message(" --- compiling Fragment --- ")
		EXECUTE_PROCESS(COMMAND rootcint -f dictFragment.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTFragmentclass} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/Fragment_linkdef.h )
		endif(myROOTFragmentclass)
			
		if(myROOTTigFragmentclass)
		message(" --- compiling TigFragment --- ")
		EXECUTE_PROCESS(COMMAND rootcint -f dictTigFragment.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTTigFragmentclass} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/TigFragment_linkdef.h )
		endif(myROOTTigFragmentclass)
		
		if(myROOTS3class)
		message(" --- compiling S3 --- ")
		EXECUTE_PROCESS(COMMAND rootcint -f dictS3.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTS3class} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/S3_linkdef.h )
		endif(myROOTS3class)

		if(myROOTSpiceclass)
		message(" --- compiling Spice --- ")
		EXECUTE_PROCESS(COMMAND rootcint -f dictSpice.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTSpiceclass} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/Spice_linkdef.h )
		endif(myROOTSpiceclass)

		if(myROOTGriffclass)
		message(" --- compiling Griffin --- ")
		EXECUTE_PROCESS(COMMAND rootcint -f dictGriffin.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTGriffclass} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/Griffin_linkdef.h )
		endif(myROOTGriffclass)

		#if(myROOTNewclass)
		#message(" --- compiling New --- ")
		#EXECUTE_PROCESS(COMMAND rootcint -f dictNew.cxx -c ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/${myROOTNewclass} ${PROJECT_SOURCE_DIR}/${myROOTclassDir}/New_linkdef.h)
		#endif(myROOTNewclass)
		
	endif(myROOTclass)
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(useROOT) 
                    
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Get examples sources from common and shared
#
#list(APPEND sources 
#  ${PROJECT_SOURCE_DIR}/../shared/src/DetectorConstruction.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/DetectorMessenger.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/EventAction.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/PhysicsList.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/PrimaryGeneratorAction.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/RunAction.cc
#  ${PROJECT_SOURCE_DIR}/../shared/src/SteppingAction.cc)

#----------------------------------------------------------------------------
# Find HBOOK (optional package)
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../common/analysis/cmake ${CMAKE_MODULE_PATH})
find_package(HBOOK QUIET)
if(NOT HBOOK_FOUND)
  message(STATUS "G4 Examples: HBOOK package not found. --> g4tools/hbook analysis disabled")  
else()
  message(STATUS "G4 Examples: HBOOK package found. --> g4tools/hbook analysis enabled")  
  add_definitions(-DG4_USE_HBOOK)  
endif()

#----------------------------------------------------------------------------
# Add g4tools/hbook sources
#

if (HBOOK_FOUND)
  include_directories(
       ${PROJECT_SOURCE_DIR}/../../common/analysis/include)
  list(APPEND sources 
       ${PROJECT_SOURCE_DIR}/../../common/analysis/src/ExG4HbookAnalysisManager.cc)

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/close.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/close.f )

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/setpawc.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/setpawc.f )

  add_custom_command(
       OUTPUT  ${PROJECT_BINARY_DIR}/setntuc.o
       COMMAND gfortran
       ARGS  -c ${Geant4_INCLUDE_DIR}/tools/hbook/setntuc.f )
  set(TOOLS_FORTRAN_OBJECTS close.o setpawc.o setntuc.o)       
endif()       

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${myexe} ${mymain} ${sources} ${headers} ${TOOLS_FORTRAN_OBJECTS})
if(useROOT)
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(myROOTclass)

		if(myROOTFragmentclass)
		message(" --- Linking :  dictFragment.cxx --- ")
		set(myROOTDict ${myROOTDict} dictFragment.cxx)
		endif(myROOTFragmentclass)
		
		if(myROOTTigFragmentclass)
		message(" --- Linking :  dictTigFragment.cxx --- ")
		set(myROOTDict ${myROOTDict} dictTigFragment.cxx)
		endif(myROOTTigFragmentclass)
		
		if(myROOTS3class)
		message(" --- Linking :  dictS3.cxx --- ")
		set(myROOTDict ${myROOTDict} dictS3.cxx)
		endif(myROOTS3class)

		if(myROOTSpiceclass)
			message(" --- Linking :  dictSpice.cxx --- ")
		set(myROOTDict ${myROOTDict} dictSpice.cxx)
		endif(myROOTSpiceclass)

		if(myROOTGriffclass)
			message(" --- Linking :  dictGriffin.cxx --- ")
		set(myROOTDict ${myROOTDict} dictGriffin.cxx)
		endif(myROOTGriffclass)
				
		#if(myROOTNewclass)
		#	message(" --- Linking :  dictNew.cxx --- ")
		#set(myROOTDict ${myROOTDict} dictNew.cxx)
		#endif(myROOTNewclass)
		
		# This will replace the semicolons in a cmake list ";" to spaces " "	 
		foreach(arg ${myROOTDict})
			set(myROOTDictSPACED "${myROOTDictSPACED} ${arg}")
		endforeach(arg ${myROOTDict})

		message(" --- Dictionaries linked are : ${myROOTDictSPACED} --- ")

	  	set(CMAKE_EXE_LINKER_FLAGS "${ROOT_LD_FLAGS} ${myROOTDictSPACED}")

	else(myROOTclass)
		set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
	endif(myROOTclass)
endif(useROOT)
target_link_libraries(${myexe} ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build detectionSystems. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(detectionSystems_SCRIPTS
    detectionSystems.in vis.mac 
  )

foreach(_script ${detectionSystems_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${myexe} DESTINATION bin)
