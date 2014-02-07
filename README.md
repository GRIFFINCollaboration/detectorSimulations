detectorSimulations
===================

The detectorSimulations package contains the Geant4 simulations for GRIFFIN, TIGRESS, and all of their auxilary detectors.  Please note that in order to respect our non-disclosure agreements, all source containing third party IP have been omitted from this repo, and can be obtained from your colleagues directly at the lab.

8pi support is available in the 8pi branch.

#Setup

To setup the simulation package on a computer with GEANT4 already present, just copy the code to your machine:

    git clone https://github.com/GRIFFINCollaboration/detectorSimulations.git
    
Then you'll need to get the files containing our NDA-protected parameters from one of your colleagues at the lab; place these in the /src directory, and everything should compile and run as expected. 

###Setup FAQ

- Yes, you need both the secret suppressed files AND their unsuppressed equivalents, not just either / or.
