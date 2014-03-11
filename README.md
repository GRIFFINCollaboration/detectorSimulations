detectorSimulations
===================

The detectorSimulations package contains the Geant4 simulations for GRIFFIN, TIGRESS, and all of their auxilary detectors.  Please note that in order to respect our non-disclosure agreements, all source containing third party IP have been omitted from this repo, and can be obtained from your colleagues directly at the lab.


#Setup

To setup the simulation package on a computer with GEANT4 already present, just copy the code to your machine:

    git clone https://github.com/GRIFFINCollaboration/detectorSimulations.git
    
Then you'll need to get the files containing our NDA-protected parameters from one of your colleagues at the lab; place these in the /src directory, and everything should compile and run as expected. 

###Setup FAQ

- Yes, you need both the secret suppressed files AND their unsuppressed equivalents, not just either / or.


#Usage

###SPICE

In order to use SPICE, you must first place/build it with:

```
/DetSys/app/addSpiceTargetChamber
/DetSys/det/addSpice 10
/DetSys/det/addS3 24
```

But you also need to define the tabulated magnetic field.  These files are 50Mb each and are different for each lens so until a better solution comes along they will be kept on the network at TRIUMF (email Mohamad or Lee for precise location):

    /DetSys/world/tabMagneticField SPICEField3D.TABLE

Then you can run the simulation as normal.

There are now additional functions in the PrimaryGeneratorAction (some of which may be merged in the near future):

    /DetSys/gun/radioactiveSourceDecay filename

will read in a decay scheme that is written previously by the user
(not to be confused with the similar /DetSys/gun/radioactiveBetaDecay)

    /DetSys/gun/radius  2.75 mm

 ... sets the radius of the source
 
    /DetSys/gun/energyrange 100 2000 50

... will emit particles every 50keV between 100 and 2000keV 
i.e. 100keV, 150keV, 200keV ... 2000keV, 100keV...
