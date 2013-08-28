//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BeamRequestGammaAndIC_H
#define BeamRequestGammaAndIC_H 1

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

#define NUM_PARTICLES 1

#define ARRAYDIM 512
#define MAXDECAYS 10
#define MAXSHELLS 4
#define MAXLEVELS 200
#define MAXQVALUE 100000
#define GCVTMAG 10
#define COLSGDATA 16
#define CHARLENGTH 1000
#define MAXBETABRANCHES 200
#define VECSIZE 10000

#define GAMMADIM 12
#define SPINDIM 12
#define THETASTEP 100

using namespace std;

class BeamRequestGammaAndIC 
{    
  public:
    BeamRequestGammaAndIC(G4String path, G4double polar);
    ~BeamRequestGammaAndIC();
    void AddBeamParticle(G4ParticleGun* beamParticle, G4int pNum, G4double timeSeconds);

    G4int    GetNumberOfGeneratedParticles() {return particle.size();}
    G4int    GetBetaBranch();
    G4bool   get_bool() {return betaBool;}
    G4bool   IsParticleAnElectron(G4int pNum);
    G4bool   IsECDecay() {return ecDecay;}
    G4double get_energy() {return energy;}
    G4double get_betaType() {return betaPlusMinusType;}
    G4double GetParticleEnergy(G4int i) {return particle[i].energy;}
    G4String GetParticleName(G4int i) {return particle[i].name;}
    G4String GetParticleShell(G4int i) {return particle[i].shell;}

    void GenerateParticles(G4double timeSeconds);
    void GenerateEfficiencyParticles();
    void PrintSurface();
    void PrintAllGeneratedParticles();
    //void PrintAllGeneratedParticlesCheck();
    void PrintLevels(G4bool printGammasBool);
    void PrintBetaBranch();

  private:
    G4ThreeVector GetParticlePosition();
    G4ThreeVector GetParticleDirection(G4int pNum, G4double timeSeconds);
    G4double      FindElectronBindingEnergy(G4int Z, G4int shell);
    G4double      findXrayEnergy(G4int Z, G4int shell);
    G4double      FindGammaIntensity(G4int energyID, G4int levelIndex);
    G4double      FindGammaIntensityWithIntConv(G4int energyID, G4int levelIndex);
    G4double      FindGammaIntConv(G4int energyID, G4int levelIndex);
    G4int         FindGammaIndex(G4int energyID, G4int levelIndex);
    G4int         FindLFromMultipolarity(const char *input);
    G4int         FindEMFromMultipolarity(const char *input); // 1 = E , 2 = M
    G4String      IntToG4String(G4int i);
    G4String      GetDate();
    G4String      GetParticleName();

    G4double W_CompleteAndPartialAlignment(G4int pNum, G4double theta, G4double sigma);
    G4double LegendreP(G4int n, G4double x);
    G4double Amax(G4int k, G4double ji, G4int L1, G4int L2, G4double jf, G4double delta);
    G4double f(G4int k, G4double jf, G4int L1, G4int L2, G4double ji);
    G4double B(G4int k, G4double j);
    G4double F(G4int k, G4double jf, G4int L1, G4int L2, G4double ji);
    G4double ClebschGordan(G4double j1, G4double m1, G4double j2, G4double m2, G4double j, G4double m);
    G4double RacahW(G4double a, G4double b, G4double c, G4double d, G4double e, G4double f);
    G4double Wigner6j(G4double J1, G4double J2, G4double J3, G4double J4, G4double J5, G4double J6);
    G4double Wigner3j(G4double j1, G4double j2, G4double j3, G4double m1, G4double m2, G4double m3);
    G4double factorial(G4double value);
    G4double Alpha(G4int k, G4double ji, G4double sigma);
    //G4double MyAlpha(G4int k, G4int pNum);
    G4double MyAlphaM(G4int k, G4double ji, G4double mi);
    G4double Rho(G4int k, G4double j, G4double sigma);
    G4double P(G4double m, G4double j, G4double sigma);

    void ReadGammaData();
    void ReadBetaData();
    void ReadPolarData();

    void CalculateBetaBranchProbabilities();
    void ClearMyParticleStruct();
    void ClearMyLevelStruct();
    void ClearMyGammaStruct();
    void ClearMyAngularStruct();
    void CalculateCoefficientsForAngularDistributions();
    void GetAllGammaTransitions();
    void GenerateWThetaAngularDistribution();
    void WriteOutWThetaData(G4int gDex, G4int jDex, G4double m);

    G4double UniformRand48();


    G4bool writeOutWThetaFiles;
    G4bool firstWriteOfWThetaFiles;
    int betaPlusMinusType;
    double inputXrayTable[ARRAYDIM][COLSGDATA][MAXSHELLS];

    G4bool betaBool;
    G4bool includeInitialBetaDepolarization;
    G4bool includeGammaDepolarization;
    G4double energy;



    G4bool eff_calib_flag;
    G4bool outputBeta;
    G4bool emit_beta_flag;
    G4bool emit_gamma_flag;
    G4bool emit_ic_flag;
    G4bool emit_xray_flag;
    G4bool angularDistribFlag;
    G4bool precessionFlag;
    G4bool ecDecay;

    G4double precessionFreq;
    G4double polarization;
    G4double spinTempBeta;
    G4double timeT1;
    G4double timeT2;

    G4int numberOfBetaBranches;
    G4int betaBranch;
    G4int myBetaBranch;
    G4int numberOfEvents;
    G4int numberOfGammaTransitions;

    G4String path;
    G4String gamma_filename;
    G4String beta_filename;
    G4String polar_filename;
    G4String matrix_filename;
    G4String k_shell_filename;
    G4String l_shell_filename;
    G4String m_shell_filename;
    G4String stringT1;
    G4String stringT2;
    G4String gammaAngularDistribution;
    G4String icAngularDistribution;
    G4String inputFile;

    G4bool k_shell_bool;
    G4bool l_shell_bool;
    G4bool m_shell_bool;

    // Level Info
    G4int    daughterZ;
    G4int    daughterA;
    G4int    parentParity;
    G4double parentSpin;
    G4double maxSpin;
    G4double qValue;
    G4double betaBranchToGS;
    G4double betaBranchRelativeSum;

    G4double ClebschGordanIntJK2[SPINDIM][SPINDIM];
    G4double ClebschGordanHalfJK2[SPINDIM][SPINDIM];
    G4double ClebschGordanIntJK4[SPINDIM][SPINDIM];
    G4double ClebschGordanHalfJK4[SPINDIM][SPINDIM];
    G4double parentMStateProbs[SPINDIM-1][SPINDIM][SPINDIM];
    G4double wTheta[ARRAYDIM][VECSIZE][SPINDIM];
    G4double wThetaSum[ARRAYDIM][VECSIZE][SPINDIM];

    struct shellLabelsStruct
    {
      char *shell;
      char *vac;
      bool bools;
    } inputXray[MAXSHELLS];

    struct effDataStruct
    {
      G4bool   isotropic;
      G4double energy;
      G4double x_momentum;
      G4double y_momentum;
      G4double z_momentum;
      G4String energyUnits;
      G4String particle;
    } eff;

    struct particleStruct
    {
      G4int    gammaIndex;
      G4int    levelEnergy;
      G4int    levelIndex;
      G4int    levelInitial;
      G4int    levelFinal;
      G4double delta;
      G4double energy;
      G4double m;
      G4double mSum[SPINDIM];
      G4String name;
      G4String shell;
    } my_particle;

    vector<particleStruct> particle;

    // Gamma Info
    struct gammaStruct
    {
      G4int    energyID;
      G4int    level;
      G4int    L1;
      G4int    L2;
      G4int    EM1;
      G4int    EM2;
      G4int    mAve[SPINDIM];
      G4double energy;
      G4double intensity;
      G4double value;
      G4double delta;
      G4double alpha;
      G4double ji;
      G4double jf;
      G4double k;
      G4double l;
      G4double m;
      G4double n;
      G4double o;
      G4double p;
      G4double q;
      G4double A2;
      G4double A4;
      G4double MyAlpha2;
      G4double MyAlpha4;
      G4double mStateProbs[SPINDIM][SPINDIM];
      G4String multipolarity;
    } my_gamma;

    vector<gammaStruct> gamma;

    struct levelStruct
    {
      G4int    energyID;
      G4int    parity;
      G4int    gammaID[GAMMADIM];
      G4double spin;
      G4double energy;
      G4double ec2betaplus; // if 0 100% beta, if 1 100% EC
      G4double betaBranchPercent;
      G4double betaBranchRelative;
    } my_level;

    vector<levelStruct> level;
};

#endif
