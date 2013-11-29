#include "BeamRequestGammaAndIC.hh"
#include "TableOfIsotopesData.hh"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits.h>
#include <cmath>
#include <sys/stat.h>   // Needed for defines used in the mkdir() call.

#include<iomanip>

using namespace std;

BeamRequestGammaAndIC::BeamRequestGammaAndIC(G4String path, G4double polar)
{
  this->path = path;
  this->gamma_filename  = "gammaData.dat";
  this->beta_filename   = "betaData.dat";
  this->polar_filename  = "polarData.dat";

  this->emit_beta_flag  = true;
  this->emit_gamma_flag = true;
  this->emit_ic_flag    = true;
  this->emit_xray_flag  = true;

//  this->eff_calib_flag  = false;
//  this->eff.particle    = "gamma";
//  this->eff.energy      = 1000.0;
//  this->eff.isotropic   = true;
//  this->eff.x_momentum  = 0.0;
//  this->eff.y_momentum  = 0.0;
//  this->eff.z_momentum  = 1.0;
//  this->eff.energyUnits = "keV";

//  this->gammaAngularDistribution = "Isotropic"; // "Isotropic" or "W(theta)"
//  this->icAngularDistribution    = "Isotropic"; // "Isotropic" or "W(theta)"
//  this->precessionFreq           = 0.0;

//  this->spinTempBeta             = 0.0;
//  this->timeT1                   = 1.0;
//  this->timeT2                   = 1.0;
//  this->stringT1                 = "INF";
//  this->stringT2                 = "INF";

  this->polarization             = polar;

  this->includeInitialBetaDepolarization = true;
  this->includeGammaDepolarization       = true;
  this->writeOutWThetaFiles              = false;
  this->firstWriteOfWThetaFiles          = true;

  ////////////////////////////////////////////////////////////////////
  // X-Ray File Properties
  ////////////////////////////////////////////////////////////////////
  this->k_shell_filename   = "kshell.dat";
  this->l_shell_filename   = "lshell.dat";
  this->m_shell_filename   = "mshell.dat";
  this->k_shell_bool       = false;
  this->l_shell_bool       = false;
  this->m_shell_bool       = false;
  this->inputXray[0].shell = const_cast<char *>("K");
  this->inputXray[1].shell = const_cast<char *>("L");
  this->inputXray[2].shell = const_cast<char *>("M");
  this->inputXray[3].shell = const_cast<char *>("N");
  this->inputXray[0].vac   = const_cast<char *>("Kvac");
  this->inputXray[1].vac   = const_cast<char *>("Lvac");
  this->inputXray[2].vac   = const_cast<char *>("Mvac");
  this->inputXray[3].vac   = const_cast<char *>("Nvac");

  srand48(time(NULL));

  ClearMyGammaStruct();
  ClearMyParticleStruct();
  ClearMyLevelStruct();

  if(!this->eff_calib_flag) {
    ReadGammaData();   // read in gamma data
    ReadBetaData();    // read in beta data
    if(polar != 0.0) {
        ReadPolarData();
    }
    else {
        this->gammaAngularDistribution = "Isotropic";
        this->icAngularDistribution = "Isotropic";
    }
    CalculateBetaBranchProbabilities();    
  }

  if(!(this->eff_calib_flag) && this->emit_gamma_flag) {
    if(gammaAngularDistribution == "W(theta)") {
      CalculateCoefficientsForAngularDistributions();  
      GenerateWThetaAngularDistribution();
    }
  }

  this->particle.clear();  
  
  PrintLevels(true);  // print to screen the levels with beta branch probabilites.
}

BeamRequestGammaAndIC::~BeamRequestGammaAndIC()
{
}

void BeamRequestGammaAndIC::ReadGammaData()
{
  // This function is straight forward, reads in data from gamma_filename (gammaData.dat)
  ClearMyGammaStruct();
  
  G4String bufString = this->path + this->gamma_filename;
  G4int ncolumns = 0;
  G4int nlines   = 0;
  G4int intBuf;
  G4double doubleBuf;
  string buf, currentline;
  G4int L1, L2, EM1, EM2;
  char *pBuf;
  char *BUF1;
  char *BUF2;

  ifstream infile;
  infile.open(bufString.c_str());

  if(!infile) {
    cerr << "unable to find/open iFile_gData" << endl;
    exit (1);
  }

  while(getline(infile, currentline)) {
    if(!currentline.empty()) {
      stringstream ss(currentline);
      ncolumns = 0;
      while (ss >> buf) {
        if(nlines == 0) { // looks at the header:gammaID,gammaEnergy,gammaIntensity,initialLevel,nEscape,Trans,MixRatio,IntConv
          if(ncolumns == 0 && buf == "gammaID"){}
          else if(ncolumns == 1 && buf == "gammaEnergy"){}
          else if(ncolumns == 2 && buf == "gammaIntensity"){}
          else if(ncolumns == 3 && buf == "initialLevel"){}
          else if(ncolumns == 4 && buf == "nEscape"){}
          else if(ncolumns == 5 && buf == "Trans"){}
          else if(ncolumns == 6 && buf == "MixRatio"){}
          else if(ncolumns == 7 && buf == "IntConv"){}
          // ...
          //else { cout << "Error in Header!" << endl; }             
        }
        else {
          doubleBuf = atof( buf.c_str() );
          intBuf = atoi( buf.c_str() );
          if(ncolumns == 0)
            this->my_gamma.energyID = intBuf;
          else if(ncolumns == 1)
            this->my_gamma.energy = doubleBuf;
          else if(ncolumns == 2)
            this->my_gamma.intensity = doubleBuf;
          else if(ncolumns == 3)
            this->my_gamma.level = intBuf;
          else if(ncolumns == 4)
            this->my_gamma.value = doubleBuf;
          else if(ncolumns == 5) {
            if(this->my_gamma.energy == 0.0) {
              this->my_gamma.multipolarity = "";
            }
            else {
              this->my_gamma.multipolarity = buf.c_str();
              if(buf.find ('+') != string::npos) {  // If String contains a '+' --> Mixing! ie. E1+M2
                pBuf = const_cast<char*> ( buf.c_str() );   
                BUF1  = strtok(pBuf, "+");
                BUF2 = strtok(NULL, "\0");
                L1 = FindLFromMultipolarity(BUF1);
                L2 = FindLFromMultipolarity(BUF2);    
                EM1 = FindEMFromMultipolarity(BUF1);
                EM2 = FindEMFromMultipolarity(BUF2);                              
              }
              else { // No Mixing ie. E1
                L1 = FindLFromMultipolarity(buf.c_str());
                L2 = 0;  
                EM1 = FindEMFromMultipolarity(buf.c_str());
                EM2 = 0;                                    
              }             
              this->my_gamma.L1 = L1;   
              this->my_gamma.L2 = L2;    
              this->my_gamma.EM1 = EM1;   
              this->my_gamma.EM2 = EM2;                  
            }
          }
          else if(ncolumns == 6)
            this->my_gamma.delta = doubleBuf;
          else if(ncolumns == 7)
            this->my_gamma.alpha = doubleBuf;
          else if(ncolumns == 8)
            this->my_gamma.k = doubleBuf;
          else if(ncolumns == 9)
            this->my_gamma.l = doubleBuf;
          else if(ncolumns == 10)
            this->my_gamma.m = doubleBuf;
          else if(ncolumns == 11)
            this->my_gamma.n = doubleBuf;
          else if(ncolumns == 12)
            this->my_gamma.o = doubleBuf;
          else if(ncolumns == 13)
            this->my_gamma.p = doubleBuf;
          else if(ncolumns == 14)
            this->my_gamma.q = doubleBuf;
          else {
            G4cout << "Error: 889234" << G4endl; 
            exit(1);
          }
        }     
        ncolumns++;
      }
      if(nlines > 0) {
        gamma.push_back(my_gamma);
        ClearMyGammaStruct();       
      }          
      nlines++;
    }
  }

  infile.close();   

  if(nlines >= ARRAYDIM) { // Check
    cout << "Error in length of Data, exceeded the expected arrays" << endl;
    exit (1);
  }
}

void BeamRequestGammaAndIC::ReadBetaData()
{
  // This function is straight forward, reads in data from beta_filename (betaData.dat)
  ClearMyLevelStruct();
  
  G4String bufString = this->path + this->beta_filename;
  G4int ncolumns = 0;
  G4int nlines = 0;
  G4int intBuf;
  G4double doubleBuf;
  string buf, currentline;

  ifstream infile;
  infile.open(bufString.c_str());
  
  this->maxSpin = 0;

  if(!infile) {
    cerr << "unable to find/open beta file" << endl;
    exit (1);
  }

  while(getline(infile, currentline)) {
    if(!currentline.empty()) {
      stringstream ss(currentline);
      ncolumns=0;
      while (ss >> buf) {
        intBuf = atoi( buf.c_str() );
        doubleBuf = atof( buf.c_str() );
        if(nlines == 0) {
          if(ncolumns == 0) {
            if(intBuf < 0) {
                this->betaPlusMinusType = 1;
            }
            else {
                this->betaPlusMinusType = 0;
            }
            this->daughterZ = abs(intBuf);
            }
            else if(ncolumns == 1) {
              this->daughterA = intBuf;
            }
        }  
        else if(nlines == 1) {
          if(ncolumns == 0)
            this->parentSpin = doubleBuf;
          else if(ncolumns == 1)
            this->parentParity = intBuf;          
        }
        else if(nlines == 2) {
          if(ncolumns == 0)
            this->qValue = doubleBuf;
          if(ncolumns == 1)
            this->betaBranchToGS = doubleBuf;
        }        
        else if(nlines > 2) {
          if(ncolumns == 0)
            this->my_level.energyID = intBuf;
          else if(ncolumns == 1)
            this->my_level.spin = doubleBuf;
            if(this->maxSpin < this->my_level.spin)
              this->maxSpin = this->my_level.spin;
          else if(ncolumns == 2)
            this->my_level.parity = intBuf; 
          else if(ncolumns == 3)
            this->my_level.energy = doubleBuf; 
          else if(ncolumns == 4)
            this->my_level.ec2betaplus = doubleBuf;
          else if(ncolumns > 4)
            this->my_level.gammaID[ncolumns-5] = intBuf; 
          else if(ncolumns > 4+GAMMADIM) {
            cout << "Error: 7798064" << endl;
            exit (1);
          }                       
        }
        ncolumns++;
      }
      if(nlines > 2) {
        level.push_back(my_level);
        ClearMyLevelStruct();       
      }          
      nlines++;
    }
  }
  infile.close(); 
}

void BeamRequestGammaAndIC::ReadPolarData()
{
  // This function is straight forward, reads in data from polar_filename (polarData.dat)
  G4String bufString = this->path + this->polar_filename;
  G4int ncolumns = 0;
  G4int nlines = 0;
  G4int intBuf;
  G4double doubleBuf;
  string buf, currentline;

  ifstream infile;
  infile.open(bufString.c_str());

  if(!infile) {
    cerr << "unable to find/open beta file" << endl;
    exit (1);
  }

  while(getline(infile, currentline)) {
    if(!currentline.empty()) {
      stringstream ss(currentline);
      ncolumns=0;
      while (ss >> buf) {
        intBuf = atoi( buf.c_str() );
        doubleBuf = atof( buf.c_str() );
        if(nlines == 0) {
          if(ncolumns == 0) {
            if(intBuf == 0) {
              gammaAngularDistribution = "Isotropic";
              cout << "Error: 9073522" << endl;
              exit(1);
            }
            else if(intBuf == 1) {
              gammaAngularDistribution = "W(theta)";
            }
            else {
              cout << "Error: 9073522" << endl;
              exit (1);
            }
          }
        }
        else if(nlines == 1) {
          if(ncolumns == 0) {
            if(intBuf == 0) {
              icAngularDistribution = "Isotropic";
            }
            else if(intBuf == 1) {
              icAngularDistribution = "W(theta)";
            }
            else {
              cout << "Error: 0349237" << endl;
              exit (1);
            }
          }
        }
        else if(nlines == 2) {
          if(ncolumns == 0)
            this->spinTempBeta = doubleBuf;
        }
        else if(nlines == 3) {
          if(ncolumns == 0)
            this->precessionFreq = doubleBuf;
        }
        else if(nlines == 4) {
          if(ncolumns == 0)
            this->timeT1 = doubleBuf;
          if(ncolumns == 1)
            this->stringT1 = buf.c_str();
        }
        else if(nlines == 5) {
          if(ncolumns == 0)
            this->timeT2 = doubleBuf;
          if(ncolumns == 1)
            this->stringT2 = buf.c_str();
        }
        else {
          cout << "Error: 3546665" << endl;
          cout << "Too many lines in " << bufString << endl;
          exit (1);
        }
        ncolumns++;
      }
      nlines++;
    }
  }
  infile.close();
}

void BeamRequestGammaAndIC::CalculateBetaBranchProbabilities()
{
  G4double betaAdd, betaSub, groundStateZeroIntensity, groundStateIntensity, sumRelative;
  G4int aboveLevelEnergy,  diffLevelEnergies, newGamma;
  
  // loop over all levels and find groundStateIntensity
  sumRelative = 0;
  for(G4int i = 0; i < this->level.size() ; i++) {
    betaAdd = 0;
    betaSub = 0;
    for(G4int j = 0; j < GAMMADIM; j++) {
      if(level[i].gammaID[j] == 0)
        break;
      betaAdd += FindGammaIntensityWithIntConv(level[i].gammaID[j], i);    
    }
    for(G4int j = i+1; j < this->level.size(); j++) {
      aboveLevelEnergy = this->level[j].energyID;
      diffLevelEnergies = aboveLevelEnergy - this->level[i].energyID;
      for(G4int k = 0; k < GAMMADIM; k++) {
        newGamma = this->level[j].gammaID[k];
        if( newGamma == 0 )
          break;
        if( newGamma == diffLevelEnergies)
          betaSub += FindGammaIntensityWithIntConv(level[j].gammaID[k], j);  
      }     
    }
    if(i==0) {
      betaAdd = betaSub;  // this normalizes g.s. to zero beta branch probability
      groundStateZeroIntensity = betaSub;
    }
    
    this->level[i].betaBranchRelative = betaAdd - betaSub;
    if(this->level[i].betaBranchRelative < 0) {
      this->level[i].betaBranchRelative = 0;
    }
    sumRelative += this->level[i].betaBranchRelative;
  }
  
  // now using correct beta branch probability to g.s. (from data file) we can 
  // use sumRelative to find the "fake" ground state intensity.
  groundStateIntensity = (this->betaBranchToGS*sumRelative)/(1-this->betaBranchToGS) + groundStateZeroIntensity;

  sumRelative = 0;
  // now do the same loop over all levels and use the "fake" ground state 
  // intensity. This will give correct beta branches.
  for(G4int i = 0; i < this->level.size() ; i++) {
    betaAdd = 0;
    betaSub = 0;
    for(G4int j = 0; j < GAMMADIM; j++) {
      if(level[i].gammaID[j] == 0)
        break;
      betaAdd += FindGammaIntensityWithIntConv(level[i].gammaID[j], i);    
    }
    for(G4int j = i+1; j < this->level.size(); j++) {
      aboveLevelEnergy = this->level[j].energyID;
      diffLevelEnergies = aboveLevelEnergy - this->level[i].energyID;
      for(G4int k = 0; k < GAMMADIM; k++) {
        newGamma = this->level[j].gammaID[k];
        if( newGamma == 0 )
          break;
        if( newGamma == diffLevelEnergies)
          betaSub += FindGammaIntensityWithIntConv(level[j].gammaID[k], j);  
      }     
    }
    if(i==0) {
      betaAdd = groundStateIntensity;
    }
    this->level[i].betaBranchRelative = betaAdd - betaSub;
    if(this->level[i].betaBranchRelative < 0) {
      this->level[i].betaBranchRelative = 0;
    }    
    sumRelative += this->level[i].betaBranchRelative;
  }
  
  // normalize to give percentages i.e. 0.8324. (less than 1)
  for(G4int i = 0; i < this->level.size() ; i++) {
    this->level[i].betaBranchPercent = this->level[i].betaBranchRelative/sumRelative;
  }
}

G4double BeamRequestGammaAndIC::FindElectronBindingEnergy(G4int Z, G4int shell)
{
  G4double energy;
  if(shell < 1 || shell > 4) {
    G4cout << "Error: 99202" << G4endl;
    exit(1);  
  }
  for(G4int z = 0; z < EBEROWS ; z++) {
    if( averageAtomicElectronBindingEnergies[z][0] == Z ) {
      energy = averageAtomicElectronBindingEnergies[z][shell];
      break;
    }
    if( z == EBEROWS-1 ) {
      G4cout << "Error: 115279" << G4endl;      
      exit(1);  
    }
  }  
  return energy;
}

G4double BeamRequestGammaAndIC::FindGammaIntConv(G4int energyID, G4int levelIndex)
{
  G4double intConv;
  G4int gammaIndex;
  
  gammaIndex = FindGammaIndex(energyID, levelIndex);
  intConv    = this->gamma[gammaIndex].alpha;

  return intConv;
}

G4double BeamRequestGammaAndIC::FindGammaIntensity(G4int energyID, G4int levelIndex)
{
  G4double intensity;
  G4int gammaIndex;

  gammaIndex = FindGammaIndex(energyID, levelIndex);
  intensity  = this->gamma[gammaIndex].intensity;

  return intensity;
}

G4double BeamRequestGammaAndIC::FindGammaIntensityWithIntConv(G4int energyID, G4int levelIndex)
{
  G4double intConv;
  G4double intensity;
  G4double out;
  G4int gammaIndex;

  gammaIndex = FindGammaIndex(energyID, levelIndex);
  intensity  = this->gamma[gammaIndex].intensity;
  intConv    = this->gamma[gammaIndex].alpha;
  out        = (1 + intConv)*intensity;

  return out;
}

G4int BeamRequestGammaAndIC::FindGammaIndex(G4int energyID, G4int levelIndex)
{
  G4int gammaIndex = -1;
  for(G4int i = 0; i < this->gamma.size(); i++) {
    if(this->gamma[i].energyID == energyID) {
      if(this->gamma[i].level == 0) { // if the gammaID is not repeated for other transitions
        gammaIndex = i;
        break;
      }
      else { // else if there are more than one transition with this gammaID, check that it's from the right level.
        if(this->gamma[i].level == this->level[levelIndex].energyID) {
          gammaIndex = i;
          break;
        }
      }
    }
  }
  if(gammaIndex == -1) {
    cout << "Error: 2344444" << endl;
    cout << "energyID = " << energyID << endl;
    cout << "levelIndex = " << levelIndex << endl;
    exit (1);
  }
  
  return gammaIndex;
}

G4int BeamRequestGammaAndIC::GetBetaBranch()
{
  G4double doubRandom, doubRandomEC, sum = 0;
  G4int betaBranch = -1;
  doubRandom = UniformRand48();

  
  for(G4int i = 0; i < this->level.size(); i++) {
    sum += this->level[i].betaBranchPercent;
    if(doubRandom <= sum) {
      betaBranch = i;
      if( get_betaType() == 1) {
          doubRandomEC = UniformRand48();
          if(doubRandomEC <= this->level[i].ec2betaplus) {
              ecDecay = true;
          }
          else {
              ecDecay = false;
          }
      }
      break;
    }
  }
  if(betaBranch == -1) {
    cout << "Error: 8893039" << endl;
    exit (1);
  }

  this->myBetaBranch = betaBranch;

  return this->myBetaBranch;
}

void BeamRequestGammaAndIC::GenerateParticles(G4double timeSeconds)
{
  G4int    nextLevelEnergyID, numGammasAtLevel, index, miindex, mfindex, betavar, shellNumber, mdex, gammaID;
  G4double ji, jf, sum, probGamma, probIntConv, alpha, probShellK, probShellL, probShellM, electronBindingEnergy;
  G4double pSpin, polarization, rhoo, expSum, beta, rand1, doubRandom, doubRandomShell, sumIntensity;
  G4double rhom[SPINDIM]         = {0.0};
  G4double temp_mSum[SPINDIM]    = {0.0};
  G4double probOfGamma[GAMMADIM] = {0.0};
  G4String shell;

  this->particle.clear(); // reset particle vector

  G4int dex = this->myBetaBranch; // get the beta branch

  for(G4int x = 0; x < ARRAYDIM; x++) {
    miindex = -1;
    mfindex = -1;
    if( dex < 0 )
      break;
    sumIntensity = 0;

    for(G4int y = 0; y < GAMMADIM; y++) {
      gammaID = this->level[dex].gammaID[y];
      if(gammaID == 0) {
        numGammasAtLevel = y;
        break;
      }
      sumIntensity += FindGammaIntensityWithIntConv(level[dex].gammaID[y], dex); // calculate total transition intensity
    }
    
    if(numGammasAtLevel == 0) // we must be at ground state.
      break;
    for(G4int y = 0; y < numGammasAtLevel; y++)
      probOfGamma[y] = (FindGammaIntensityWithIntConv(level[dex].gammaID[y], dex))/(sumIntensity);

    gammaID = -1;
    sum = 0;
    doubRandom = UniformRand48();

    for(G4int j = 0; j < GAMMADIM ; j++) {
      sum += probOfGamma[j];
      if(doubRandom <= sum) {
        gammaID = level[dex].gammaID[j];
        break;
      }
    }
    if(gammaID == -1) {
      G4cout << "Error: 886212" << G4endl;
      exit(1);
    }
    // found transition

    index = FindGammaIndex(gammaID, dex); 

    if(gammaAngularDistribution == "W(theta)") {
      beta = this->spinTempBeta;
      pSpin = this->level[dex].spin;
      polarization = (exp(beta)-1)/(exp(beta)+1);
      if(stringT1 == "INF") {
      }
      else if(stringT1 == "second" || stringT1 == "seconds" || stringT1 == "sec" || stringT1 == "secs" || stringT1 == "s") {
        polarization = polarization * exp(-1*(timeSeconds/this->timeT1));
      }
      else {
        G4cout << "Error: 32426566" << G4endl;
        G4cout << "Please input T1 time in seconds." << G4endl;
        exit(1);
      }
      if(stringT2 == "INF") {
      }
      else if(stringT2 == "second" || stringT2 == "seconds" || stringT2 == "sec" || stringT2 == "secs" || stringT2 == "s") {
        polarization = polarization * exp(-1*(timeSeconds/this->timeT2));
      }
      else {
        G4cout << "Error: 32426567" << G4endl;
        G4cout << "Please input T2 time in seconds." << G4endl;
        exit(1);
      }

      if(polarization == 1.0) {
        beta = this->spinTempBeta;
      }
      else {
        beta = log((1+polarization)/(1-polarization));
      }
      
      if(this->particle.size() == 0) { // first decay in series. In other words, we just beta decayed to this new level.
        expSum = 1;
        for(G4int i = 1; i <= 2*this->parentSpin ; i++ )
          expSum = expSum + exp(-1*i*beta);

        rhoo = 1/expSum;
        rhom[0] = rhoo;
        for(G4int i = 1; i <= 2*this->parentSpin ; i++ )
          rhom[i] = rhom[i-1]*exp(-1*beta);
  
        sum = 0;
        for(G4int k = 0; k <= 2*this->parentSpin ; k++ ) {  // mSum goes from -J...+J
          temp_mSum[k] = rhom[int(2*this->parentSpin) - k]; // this is just from pure polarization
          sum += temp_mSum[k];
        }
        for(G4int k = 0; k <= 2*this->parentSpin ; k++ ) {  // mSum goes from -J...+J
          temp_mSum[k] = temp_mSum[k]/sum;
        }

        if(this->includeInitialBetaDepolarization) { // Now we include the de-polarization via the Beta decay
          ji = this->parentSpin;
          jf = this->gamma[index].ji;    
          betavar = (jf-ji+4);
          if(  betavar >= SPINDIM ) {
            G4cout << "Error: 234211434" << G4endl;
            G4cout << "Beta Decay Transition Too Large! Steps over 4 integers!" << G4endl;
            exit(1);
          }

          for(G4int l = 0; l <= 2*jf ; l++ ) {
            sum = 0;
            for(G4int k = 0; k <= 2*ji ; k++ ) {
              sum += (this->parentMStateProbs[betavar][k][l] * temp_mSum[k]);
            }
            this->my_particle.mSum[l] = sum;
          }
        } // End if(this->includeInitialBetaDepolarization) {
        else { // if we don't include beta decay depolarization then put 100% m state probability into m = +jf, ie. perfect polariztion.
          ji = this->parentSpin;
          jf = this->gamma[index].ji;          
          if(polarization == 1.0) {
            for(G4int l = 0; l <= 2*jf ; l++ ) {
              this->my_particle.mSum[l] = 0.0;
              if(l == 2*jf)
                this->my_particle.mSum[l] = 1.0;
            }
          }
          else {
            G4cout << "Error: 123214323" << G4endl;
            G4cout << "Need More Info" << G4endl;
            exit(1);
          }
        }
      } // End if(this->particle.size() == 0) { //first decay
      else { // If this is NOT the first decay in the series, then we can depolarize the system from gamma decay.
             // If we don't include gamma depolarization, then put the mState Prob at 100% for m = +jf,
             // but this only makes sense if out polarization is 1.0, ie perfect polarization which implies
             // that if we don't include gamma depolarization, we  shouldn't include beta depolarization either.
        if(this->includeGammaDepolarization) { // include gamma decay depolarization physics
          ji = this->gamma[this->particle[this->particle.size()-1].gammaIndex].ji;
          jf = this->gamma[this->particle[this->particle.size()-1].gammaIndex].jf;
          for(G4int l = 0; l <= 2*jf ; l++ ) {
            sum = 0;
            for(G4int k = 0; k <= 2*ji ; k++ ) {
              sum += (this->gamma[this->particle[this->particle.size()-1].gammaIndex].mStateProbs[k][l] * this->particle[this->particle.size()-1].mSum[k]);
            }
            this->my_particle.mSum[l] = sum;
          }
        }
        else { // do NOT include gamma decay depolarization physics
          jf = this->gamma[this->particle[this->particle.size()-1].gammaIndex].jf;        
          if(polarization == 1.0) {
            for(G4int l = 0; l <= 2*jf ; l++ ) {
              this->my_particle.mSum[l] = 0.0;
              if(l == 2*jf)
                this->my_particle.mSum[l] = 1.0; //put 100% m state probability into m = +jf;
            }
          }
          else {
            G4cout << "Error: 2131222" << G4endl;
            G4cout << "Need More Info" << G4endl;
            exit(1);
          }
        }        
      }
      
      sum = 0;
      for(G4int l = 0; l <= 2*jf ; l++ ) {
        sum += this->my_particle.mSum[l];
      }
      if(sum < 0.9999999999 || sum > 1.0000000001) {
        // Make sure the sum of mState Probs equal 100%.
        // Since the sum depends on calculations of Clebsch-Gordan coefficients, it probably will not EXACTLY equal 100%,
        // due to rounding in the calculations. It should however, be good to at least 10 decimal places.
        G4cout << "the sum = " << setprecision(12) << sum << G4endl;
        exit(1);
      }
            
      // Pick an m state!!!
      rand1 = UniformRand48();
      mdex = -1;
      sum = 0;
      for(G4int l = 0; l <= 2*jf ; l++ ) {
        sum += this->my_particle.mSum[l];
        if(rand1 < sum) {
          mdex = l;
          this->my_particle.m = (-1*jf)+double(mdex);
          this->gamma[index].mAve[mdex]++;
          break;
        }
      }
      if(mdex == -1) {
        G4cout << "Error: 21312322" << G4endl;
        G4cout << "Couldn't find an m state, this error should be impossible!" << G4endl;
        exit(1);
      }
      // Now that we have picked an m state, let's "collapse" the m-state probabilities to reflect our chosen m state.
      for(G4int l = 0; l <= 2*jf ; l++ ) {
        this->my_particle.mSum[l] = 0.0;
        if(l == mdex)
          this->my_particle.mSum[l] = 1.0;
      }
    } // End if(gammaAngularDistribution == "W(theta)") {
    
    alpha       = FindGammaIntConv(gammaID, dex); // total internal conversion coefficient
    probGamma   = (1/(1+alpha));
    probIntConv = (alpha/(1+alpha));

    doubRandom = UniformRand48();

    if( doubRandom >= 0 && doubRandom <= probGamma) {
      // Emitt a Gamma!
      this->my_particle.name         = "gamma";
      this->my_particle.energy       = this->gamma[index].energy;
      this->my_particle.levelEnergy  = this->level[dex].energy;
      this->my_particle.gammaIndex   = index;
      this->my_particle.levelIndex   = dex;
      this->my_particle.levelInitial = this->level[dex].energyID;
      this->my_particle.levelFinal   = this->level[dex].energyID - gammaID;
      this->my_particle.delta        = this->gamma[index].delta;
    }
    else if(doubRandom > probGamma && doubRandom <= (probGamma+probIntConv)) {
      // Internal Conversion!
      // Figure out what shell to emit from.
      probShellK = this->gamma[index].k/alpha;
      probShellL = this->gamma[index].l/alpha;      
      probShellM = this->gamma[index].m/alpha; 

      doubRandomShell = UniformRand48();

      if(doubRandomShell >= 0 && doubRandomShell <= probShellK) {
        shellNumber = 1; // K Shell
        shell = "K";
      }
      else if(doubRandomShell > probShellK && doubRandomShell <= (probShellL+probShellK) ) {
        shellNumber = 2; // L Shell
        shell = "L";        
      }
      else if(doubRandomShell > (probShellL+probShellK) && doubRandomShell <= (probShellL+probShellK+probShellM) ) {
        shellNumber = 3; // M Shell
        shell = "M";        
      }
      else { // if not K, L, or M shell, emit from N shell.
             // we can not go any higher since we don't have electron binding energies for higher shells.
        shellNumber = 4; // N Shell
        shell = "N";        
      }

      electronBindingEnergy = FindElectronBindingEnergy(daughterZ, shellNumber); // Find electron binding energy

      this->my_particle.name         = "e-";
      this->my_particle.energy       = this->gamma[index].energy - electronBindingEnergy;
      this->my_particle.levelEnergy  = this->level[dex].energy;
      this->my_particle.gammaIndex   = index;
      this->my_particle.levelIndex   = dex;
      this->my_particle.levelInitial = this->level[dex].energyID;
      this->my_particle.levelFinal   = this->level[dex].energyID - gammaID;
      this->my_particle.delta        = this->gamma[index].delta;
      this->my_particle.shell        = shell;
    }
    else {
      G4cout << "Error! Could not select gamma decay or internal conversion" << G4endl;
      G4cout << "Error: 992302" << G4endl;
      exit(1);
    }        

    nextLevelEnergyID = this->level[dex].energyID - gammaID; // Determind next energy level ID

    if(nextLevelEnergyID == 0) { // if nextLevelEnergyID = 0, we must have hit the ground state. No more particles, break!
      particle.push_back(my_particle);       
      ClearMyParticleStruct();         
      break;
    }

    dex = -1; // Find the next level!
    for(G4int y = 0; y < this->level.size() ; y++) {
      if(this->level[y].energyID == nextLevelEnergyID ) {
        dex = y;
        particle.push_back(my_particle);         
        ClearMyParticleStruct();           
        break;    
      }
    }
    if(dex == -1) {
      G4cout << "Error, can't find next energy level!" << G4endl;   
      G4cout << "nextLevelEnergyID = " << nextLevelEnergyID << G4endl;   
      G4cout << "Error: 9123211" << G4endl;   
      exit(1);
    }    
  }
}

void BeamRequestGammaAndIC::GenerateEfficiencyParticles()
{
  this->particle.clear();

  if( this->eff.energyUnits == "keV" ) {
    if(this->eff.particle == "gamma")
      this->my_particle.name = "gamma";
    else if(this->eff.particle == "neutron")
      this->my_particle.name = "neutron";
    else if(this->eff.particle == "electron")
      this->my_particle.name = "e-";
    else if(this->eff.particle == "positron")
      this->my_particle.name = "e+";
    else {
      G4cout << "Error 4534534" << G4endl;
      exit(1);
    }
    this->my_particle.gammaIndex = 0;
    this->my_particle.levelIndex = 0;    
    this->my_particle.levelInitial = 0;    
    this->my_particle.levelFinal = 0;        
    this->my_particle.delta = 0;     
    this->my_particle.energy = this->eff.energy;
    this->my_particle.m = 0;     
    particle.push_back(my_particle);
    ClearMyParticleStruct();
  }
  else {
    G4cout << "Error 444511" << G4endl;
    exit(1);
  }
}

G4int BeamRequestGammaAndIC::FindLFromMultipolarity(const char *input)
{  
  G4int buf = -1;
  G4int L = -1;  
  const int NUMBER_MULTIPOLES = 10;

  struct emMultipoleMap {
    char *sym;
  } emMultipole[NUMBER_MULTIPOLES];

  emMultipole[0].sym=const_cast<char *>("E1");
  emMultipole[1].sym=const_cast<char *>("M1");
  emMultipole[2].sym=const_cast<char *>("E2");
  emMultipole[3].sym=const_cast<char *>("M2");
  emMultipole[4].sym=const_cast<char *>("E3");
  emMultipole[5].sym=const_cast<char *>("M3");
  emMultipole[6].sym=const_cast<char *>("E4");
  emMultipole[7].sym=const_cast<char *>("M4");     
  emMultipole[8].sym=const_cast<char *>("E5");
  emMultipole[9].sym=const_cast<char *>("M5"); 

  for(G4int i = 0; i < NUMBER_MULTIPOLES; i++) {
    if(!strcmp(emMultipole[i].sym,input)) {
      buf = i;
      break;
    }
  }
  if(buf == -1) {
    G4cout << "input = " << input << G4endl;
    G4cout << "Error: 6734534" << G4endl;          
    exit(1);
  }
  if(buf == 0 || buf == 1)
    L = 1;
  else if(buf == 2 || buf == 3)
    L = 2;
  else if(buf == 4 || buf == 5)
    L = 3;
  else if(buf == 6 || buf == 7)
    L = 4;  
  else if(buf == 8 || buf == 9)
    L = 5;  
  else {
    G4cout << "Error: 878965" << G4endl;          
    exit(1);
  }            
  return L; 
}

G4int BeamRequestGammaAndIC::FindEMFromMultipolarity(const char *input)
{  // E=1, M=2
  G4int buf = -1;
  G4int EM = -1;  
  const int NUMBER_MULTIPOLES = 10;

  struct emMultipoleMap {
    char *sym;
  } emMultipole[NUMBER_MULTIPOLES];

  emMultipole[0].sym=const_cast<char *>("E1");
  emMultipole[1].sym=const_cast<char *>("M1");
  emMultipole[2].sym=const_cast<char *>("E2");
  emMultipole[3].sym=const_cast<char *>("M2");
  emMultipole[4].sym=const_cast<char *>("E3");
  emMultipole[5].sym=const_cast<char *>("M3");
  emMultipole[6].sym=const_cast<char *>("E4");
  emMultipole[7].sym=const_cast<char *>("M4");     
  emMultipole[8].sym=const_cast<char *>("E5");
  emMultipole[9].sym=const_cast<char *>("M5"); 

  for(G4int i = 0; i < NUMBER_MULTIPOLES; i++) {
    if(!strcmp(emMultipole[i].sym,input)) {
      buf = i;
      break;
    }
  }
  if(buf == -1) {
    G4cout << "Error: 2342333" << G4endl;          
    exit(1);
  }
  if(buf == 0 || buf == 2 || buf == 4 || buf == 6 || buf == 8)
    EM = 1;
  else if(buf == 1 || buf == 3 || buf == 5 || buf == 7 || buf == 9)
    EM = 2;
  else {
    G4cout << "Error: 778744" << G4endl;          
    exit(1);
  }            
  return EM; 
}

void BeamRequestGammaAndIC::AddBeamParticle(G4ParticleGun* beamParticle, G4int pNum, G4double timeSeconds)
{
  beamParticle -> SetParticleDefinition(G4ParticleTable::GetParticleTable() -> FindParticle( particle[pNum].name ));
  beamParticle -> SetParticleEnergy( (particle[pNum].energy)/1000.0 );  // in MeV!
  beamParticle -> SetParticlePosition( GetParticlePosition() );
  beamParticle -> SetParticleMomentumDirection( GetParticleDirection(pNum, timeSeconds) ); 
}

G4bool BeamRequestGammaAndIC::IsParticleAnElectron(G4int pNum)
{
  G4bool booli;

  if(particle[pNum].name == "e-")
    booli = true;
  else if(particle[pNum].name == "gamma")
    booli = false;
  else {
    G4cout << "Error, not expected particle type!" << G4endl;
    G4cout << "Error: 9927251" << G4endl;  
    exit(1);
  }

  return booli;
}


G4ThreeVector BeamRequestGammaAndIC::GetParticlePosition()
{
  G4ThreeVector particlePosition = G4ThreeVector(0.0, 0.0, 0.0);
  return particlePosition;
}

G4ThreeVector BeamRequestGammaAndIC::GetParticleDirection(G4int pNum, G4double timeSeconds)
{
  G4ThreeVector particleDirection;
  G4double x,y,z;
//  G4double rand1, rand2, phi, theta, A2, A4, Alpha2, Alpha4, B, a, b, func, funcSum, norm, r, val, mytheta, myphi;
  G4double rand1, rand2, phi, theta, funcSum, r ;
  G4double x_prime,y_prime,z_prime;
  G4double rotate;
  G4bool breakOut = false;
    
  // G4double P2;
  // G4double P4;    
  // G4int mdex,th,ph;
  G4int mdex ;  

  funcSum = 0;

  if(this->eff_calib_flag) {
    if(this->eff.isotropic) {
      // Randomize the direction over 4pi
      rand1 = UniformRand48();
      rand2 = UniformRand48();
      phi = rand1*2.0*M_PI;
      G4double angle = M_PI;   //for isotropic scattering
      theta = acos(1.0 - (rand2*(1.0 - cos(angle))));
    
      x = sin(theta) * cos(phi);
      y = sin(theta) * sin(phi);
      z = cos(theta);
    }
    else {
      x = this->eff.x_momentum;
      y = this->eff.y_momentum;
      z = this->eff.z_momentum;            
    }
  }
  else if(this->particle[pNum].name == "gamma") {
    if(gammaAngularDistribution == "W(theta)" && particle[pNum].energy >= 0) {
    
      mdex = int(this->particle[pNum].m + this->gamma[this->particle[pNum].gammaIndex].ji);
      r = UniformRand48();
      
      theta = -1.0;
      funcSum = 0;
      for(G4int y = 0; y < VECSIZE; y++) {
        if(r < this->wThetaSum[this->particle[pNum].gammaIndex][y][mdex]) {
          theta = (double(y)/(double(VECSIZE-1)))*M_PI;
          break;
        }
      }
      if(theta == -1.0) {
        G4cout << "this->particle[pNum].gammaIndex = " << this->particle[pNum].gammaIndex << G4endl;
        G4cout << "mdex = " << mdex << G4endl;
        for(G4int y = 0; y < VECSIZE; y++) {
          G4cout << this->wThetaSum[this->particle[pNum].gammaIndex][y][mdex] << "\t";// << G4endl;
        }
        G4cout << G4endl;
        G4cout << "Error 425234232 : No Theta Found" << G4endl;
        exit(1);    
      }      

      rand2 = UniformRand48();
      phi = rand2*2.0*M_PI;

      x = sin(theta) * cos(phi);
      y = sin(theta) * sin(phi);
      z = cos(theta);

      if( this->precessionFreq != 0.0 ) {

        // Right now axis of quantization is along z-axis.
        // Let's rotate our frame such that it is along the x-axis
        rotate = 90.0*(M_PI/180.0);

        x_prime = x*cos(rotate) + z*sin(rotate);
        y_prime = y;
        z_prime = -x*sin(rotate) + z*cos(rotate);

        x = x_prime;
        y = y_prime;
        z = z_prime;

        // Now rotate
        rotate = this->precessionFreq * (2.0*M_PI) * timeSeconds;

        x_prime = x*cos(rotate) - y*sin(rotate);
        y_prime = x*sin(rotate) + y*cos(rotate);
        z_prime = z;

        x = x_prime;
        y = y_prime;
        z = z_prime;
      }
    }
    else if(gammaAngularDistribution == "Isotropic") {
      // Randomize the direction over 4pi
      rand1 = UniformRand48();
      rand2 = UniformRand48();
      phi = rand1*2.0*M_PI;
      G4double angle = M_PI;   //for isotropic scattering
      theta = acos(1.0 - (rand2*(1.0 - cos(angle))));

      x = sin(theta) * cos(phi);
      y = sin(theta) * sin(phi);
      z = cos(theta);
    }
    else {
      G4cout << "Error 324324333 : gamma angular distributions" << G4endl;
      exit(1);
    }
  }
  else if(this->particle[pNum].name == "e-") {
    if(icAngularDistribution == "Isotropic") {
      // Randomize the direction over 4pi
      rand1 = UniformRand48();
      rand2 = UniformRand48();
      phi = rand1*2.0*M_PI;
      G4double angle = M_PI;   //for isotropic scattering
      theta = acos(1.0 - (rand2*(1.0 - cos(angle))));
      x = sin(theta) * cos(phi);
      y = sin(theta) * sin(phi);
      z = cos(theta);
    }
    else {
      G4cout << "Error 555232 : ic angular distributions" << G4endl;
      exit(1);
    }   
  }
  
  particleDirection = G4ThreeVector(x,y,z);

  return particleDirection;
}

void BeamRequestGammaAndIC::ClearMyParticleStruct()
{
  this->my_particle.name = "";
  this->my_particle.energy = 0;
  this->my_particle.gammaIndex = 0;
  this->my_particle.levelIndex = 0;
  this->my_particle.levelInitial = 0;    
  this->my_particle.levelFinal = 0;    
  this->my_particle.delta = 0; 
  this->my_particle.m = 0;          
  this->my_particle.shell = "";  
  for(G4int i = 0 ; i < SPINDIM ; i++ ) {
    this->my_particle.mSum[i] = 0.0;
  }    
}

void BeamRequestGammaAndIC::ClearMyLevelStruct()
{
  this->my_level.energyID = 0;
  this->my_level.spin = 0;
  this->my_level.parity = 0;
  this->my_level.energy = 0;
  this->my_level.ec2betaplus = 0;
  this->my_level.betaBranchPercent = 0;
  this->my_level.betaBranchRelative = 0;    
  
  for(G4int i = 0 ; i < GAMMADIM ; i++ ) {
    this->my_level.gammaID[i] = 0;
  }
}

void BeamRequestGammaAndIC::ClearMyGammaStruct()
{
  this->my_gamma.energyID = 0;
  this->my_gamma.energy = 0;
  this->my_gamma.intensity = 0;
  this->my_gamma.level = 0;
  this->my_gamma.value = 0;  
  this->my_gamma.multipolarity = "";
  this->my_gamma.L1 = 0;
  this->my_gamma.L2 = 0;  
  this->my_gamma.EM1 = 0;  
  this->my_gamma.EM2 = 0;      
  this->my_gamma.delta = 0;
  this->my_gamma.alpha = 0;
  this->my_gamma.ji = 0;
  this->my_gamma.jf = 0;
  this->my_gamma.k = 0;    
  this->my_gamma.l = 0;    
  this->my_gamma.m = 0;    
  this->my_gamma.n = 0;    
  this->my_gamma.o = 0;    
  this->my_gamma.p = 0;    
  this->my_gamma.q = 0;    
  this->my_gamma.A2 = 0;    
  this->my_gamma.A4 = 0;
  this->my_gamma.MyAlpha2 = 0;    
  this->my_gamma.MyAlpha4 = 0;  
  for(G4int i = 0 ; i < SPINDIM ; i++ ) {
    for(G4int j = 0 ; j < SPINDIM ; j++ ) {
      this->my_gamma.mStateProbs[i][j] = 0;
      for(G4int k = 0 ; k < SPINDIM-1 ; k++ ) {
        this->parentMStateProbs[k][i][j] = 0;
      }      
    }
  }
  for(G4int i = 0 ; i < SPINDIM ; i++ ) {
    this->my_gamma.mAve[i] = 0;
  }        
}

void BeamRequestGammaAndIC::PrintSurface()
{  
//  G4cout << "========================================== surface counts = " << surfaceCounts << " ===========================================================" << G4endl;
//  for(G4int i = 0; i < 181; i++)
//  {  
//    for(G4int j = 0; j < 361; j++)
//    {
//      G4cout << this->surface[i][j] << " ";
//    }
//    G4cout << G4endl;
//  }
//  G4cout << "=====================================================================================================" << G4endl;    
}


void BeamRequestGammaAndIC::PrintAllGeneratedParticles()
{
  G4cout << "==============================================New Decay======================================================" << G4endl;
  //G4cout << "Number of Particles \t" << particle.size() << "\t" << G4endl;

  for( G4int i = 0 ; i < particle.size() ; i++ ) {
    G4cout << "Name   = " << particle[i].name << "\t\t";
    G4cout << "Energy = " << setprecision(2) << fixed << particle[i].energy << "\t\t";
    G4cout << "Leveli = " <<particle[i].levelInitial << "\t\t";    
    G4cout << "Levelf = " <<particle[i].levelFinal << "\t\t";   
    G4cout << "Delta  = " << setprecision(2) << fixed <<particle[i].delta << "\t\t" << G4endl;    
    G4cout << "Ji     = " << this->gamma[particle[i].gammaIndex].ji << "\t\t";
    G4cout << "Jf     = " << this->gamma[particle[i].gammaIndex].jf << "\t\t"; 
    G4cout << "L1     = " << this->gamma[particle[i].gammaIndex].L1 << "\t\t"; 
    G4cout << "L2     = " << this->gamma[particle[i].gammaIndex].L2 << "\t\t";     
    G4cout << "m      = " << particle[i].m << "\t\t" << G4endl;         
    G4cout << "A2     = " << setprecision(3) << fixed << this->gamma[particle[i].gammaIndex].A2 << "\t\t";    
    G4cout << "A4     = " << setprecision(3) << fixed << this->gamma[particle[i].gammaIndex].A4 << "\t\t"; 
    G4cout << "Alpha2 = " << setprecision(3) << fixed << this->gamma[particle[i].gammaIndex].MyAlpha2 << "\t\t";    
    G4cout << "Alpha4 = " << setprecision(3) << fixed << this->gamma[particle[i].gammaIndex].MyAlpha4 << "\t\t" << G4endl; 
    G4cout << " pops = [ "; 
    for( G4int j = 0 ; j <= (2*(this->gamma[particle[i].gammaIndex].ji)) ; j++ ) {
      G4cout << setprecision(4) << fixed << particle[i].mSum[j];  
      if( j!= (2*(this->gamma[particle[i].gammaIndex].ji))) {
        G4cout << "\t\t";
      }
    }       
    G4cout << " ]" << G4endl; 
    
    G4double sum = 0;
    G4cout << "Apops = [ "; 
    for( G4int j = 0 ; j <= (2*(this->gamma[particle[i].gammaIndex].ji)) ; j++ ) {
      sum = sum + this->gamma[particle[i].gammaIndex].mAve[j];
    } 
    for( G4int j = 0 ; j <= (2*(this->gamma[particle[i].gammaIndex].ji)) ; j++ ) {
      G4cout << setprecision(4) << fixed << (this->gamma[particle[i].gammaIndex].mAve[j]/sum) ;  
      if( j!= (2*(this->gamma[particle[i].gammaIndex].ji))) {
        G4cout << "\t\t"; 
      }
    }            
    G4cout << " ]" << G4endl;     
    G4cout << "-------------------------------------------------------------------------------------------------------------" << G4endl;  
  }
  G4cout << "=============================================================================================================" << G4endl;      
}

void BeamRequestGammaAndIC::PrintLevels(G4bool printGammasBool)
{
  G4cout << "********************************************************************************" << G4endl;
  G4cout << "********\t\t\tLevel Information\t\t\t********" << G4endl;
  G4cout << "********************************************************************************" << G4endl;
  G4cout << "Number of Levels \t" << level.size() << "\t" << G4endl;
  G4cout << "********************************************************************************" << G4endl;
  G4cout << "Level \tEnergyID \tbetaBranch \tEnergy \tSpin \tParity \tGammaID" << G4endl;
  G4cout << "--------------------------------------------------------------------------------" << G4endl;
  for( int i = level.size()-1 ; i >= 0 ; i-- )
  {
    G4cout << i << "\t" << level[i].energyID << "\t\t" << setprecision(4) << fixed << level[i].betaBranchPercent << "\t\t" << setprecision(2) << fixed << level[i].energy << "\t" << level[i].spin << "\t" << level[i].parity << G4endl;
    if(printGammasBool)
    {
      for( int j = 0 ; j < GAMMADIM ; j++ )
      {
        if(level[i].gammaID[j] == 0)
        {
          break;
        }
        G4cout << "\t\t\t\t\t\t\t\t"<< level[i].gammaID[j] << G4endl;
      }
    }
    G4cout << "--------------------------------------------------------------------------------" << G4endl;
  }
  G4cout << "********************************************************************************" << G4endl;
}

void BeamRequestGammaAndIC::PrintBetaBranch()
{
  G4cout << "betaBranch = \t" << this->myBetaBranch << "\t" << G4endl;
}

G4String BeamRequestGammaAndIC::IntToG4String(G4int i) 
{
  std::string s;
  std::stringstream out;
  out << i;
  s = out.str();
  return s;
}

void BeamRequestGammaAndIC::CalculateCoefficientsForAngularDistributions()
{
  G4double ji,jf,mi,mf,delta,L1,L2,m;
  G4int gammaIndex, gammaID, energyID, energyDiff;
  G4int levelf = -1;
  G4int leveli = -1;  
  G4int maxSpinInt;  
  G4double sum = 0;

  for(G4int i = 0; i < SPINDIM ; i++ )
  {
    for(G4int j = 0; j < SPINDIM ; j++ )
    {
      ClebschGordanIntJK2[i][j] = 0.0;
      ClebschGordanHalfJK2[i][j] = 0.0;
      ClebschGordanIntJK4[i][j] = 0.0;
      ClebschGordanHalfJK4[i][j]= 0.0;
    }
  }

  // Depolarization due to beta
  ji = this->parentSpin;
  for(G4int i = 0; i < SPINDIM-1 ; i++ )
  {  
    jf = ji-4+i;
    for(G4int k = 0; k <= 2*ji ; k++ )
    {
      mi = (-1*ji)+k;
      sum = 0;
      for(G4int l = 0; l <= 2*jf ; l++ )
      {
        mf = (-1*jf)+l;
        this->parentMStateProbs[i][k][l] = pow(ClebschGordan(ji, mi, jf, mf, ji+jf, mi+mf),2);
        sum += this->parentMStateProbs[i][k][l];
      }
      for(G4int l = 0; l <= 2*jf ; l++ ) // Normalize
      {
        this->parentMStateProbs[i][k][l] = this->parentMStateProbs[i][k][l] / sum;
      }  
    }
  }    
  
  for(G4int i = 1; i < this->level.size(); i++ )
  {
    leveli = i;
    for(G4int j = 0; j < GAMMADIM; j++ )
    {
      gammaID = this->level[leveli].gammaID[j];
      if(gammaID == 0)
        break;
      energyID = this->level[leveli].energyID;
      energyDiff = energyID - gammaID;
      for(G4int k = 0; k < i; k++ )
      {
        if(this->level[k].energyID == energyDiff)
        {
          levelf = k;
          break;
        }
      }
      if(levelf == -1)
      {
        G4cout << "Error: 8956321" << G4endl;
        exit(1);
      }
      ji = this->level[leveli].spin;
      jf = this->level[levelf].spin;
      gammaIndex = FindGammaIndex(gammaID, leveli);
      this->gamma[gammaIndex].ji = ji;
      this->gamma[gammaIndex].jf = jf;
        
      L1 = this->gamma[gammaIndex].L1;
      L2 = this->gamma[gammaIndex].L2;
      delta = this->gamma[gammaIndex].delta;        
          
      this->gamma[gammaIndex].A2 = Amax(2,ji,L1,L2,jf,delta);
      this->gamma[gammaIndex].A4 = Amax(4,ji,L1,L2,jf,delta);  
      
      for(G4int k = 0; k <= 2*ji ; k++ )
      {
        mi = (-1*ji)+k;
        sum = 0;
        for(G4int l = 0; l <= 2*jf ; l++ )
        {
          mf = (-1*jf)+l;
          this->gamma[gammaIndex].mStateProbs[k][l] = pow(ClebschGordan(ji, mi, jf, mf, ji+jf, mi+mf),2);
          sum += this->gamma[gammaIndex].mStateProbs[k][l];
        }
        for(G4int l = 0; l <= 2*jf ; l++ ) // Normalize
        {
          this->gamma[gammaIndex].mStateProbs[k][l] = this->gamma[gammaIndex].mStateProbs[k][l] / sum;
        }        
        
      }        
    }
  }

  maxSpinInt = ceil(abs(this->maxSpin));
  for(G4int j = 0; j <= maxSpinInt ; j++ )
  {
    for(G4int k = 0; k <= 2*j ; k++ )
    {
      m = -j+k;
      ClebschGordanIntJK2[j][k] = ClebschGordan(j, m, j, -1*m, 2, 0);
    }
  }
  for(G4int j = 0; j <= maxSpinInt ; j++ )
  {
    for(G4int k = 0; k <= 2*j ; k++ )
    {
      m = -j+k;
      ClebschGordanIntJK4[j][k] = ClebschGordan(j, m, j, -1*m, 4, 0);
    }
  }  
  for(G4int j = 0; j <= maxSpinInt ; j++ )
  {
    for(G4int k = 0; k <= 2*j+1 ; k++ )
    {
      m = -1*(j+0.5)+k;
      ClebschGordanHalfJK2[j][k] = ClebschGordan(j+0.5, m, j+0.5, -1*m, 2, 0);
    }
  }
  for(G4int j = 0; j <= maxSpinInt ; j++ )
  {
    for(G4int k = 0; k <= 2*j+1 ; k++ )
    {
      m = -1*(j+0.5)+k;
      ClebschGordanHalfJK4[j][k] = ClebschGordan(j+0.5, m, j+0.5, -1*m, 4, 0);
    }
  }    
}

G4String BeamRequestGammaAndIC::GetDate()
{
  string date;
  char *strctime;
  time_t curr;
  curr=time(0);
  strctime = ctime(&curr);

  // Example of strctime: "Mon Apr 27 15:53:41 2009"
  char *pBuf = const_cast<char*> ( strctime );   
  char *BUF1  = strtok(pBuf, " ");
  char *BUF2 = strtok(NULL, " "); 
  char *BUF3 = strtok(NULL, " ");  
  char *BUF4 = strtok(NULL, " ");         
  char *BUF5 = strtok(NULL, "\0");    

  date = "";
  date.append(BUF1,3);
  date.append(1,'_');
  date.append(BUF2,3);
  date.append(1,'_'); 
  date.append(BUF3,2);
  date.append(1,'_');  
  date.append(BUF4,8);   
  date.append(1,'_');   
  date.append(BUF5,4);   

  return date ;      
}

////////////////////////////////////////////////////////////////////////
// General Expression for gamma-distribution from T. Yamazaki
// Nuclear Data, Section A, Vol 3, Number 1, August 1967
// For Functions from H.A. Tolhoek and J.A.M. Cox Papers
// Physica XIX (1953) 101-119 and
// Physica XX (1954) 1310-1313 see minotaur1.50
////////////////////////////////////////////////////////////////////////
G4double BeamRequestGammaAndIC::W_CompleteAndPartialAlignment(G4int pNum, G4double theta, G4double sigma)  
{
  G4int L1 = 1;   
  G4int L2 = 1;
  G4double ji = 1;
  G4double jf = 1;
  G4double delta = 0;
  
  G4double A2;
  G4double A4;
  G4double P2 = LegendreP(2,cos(theta));
  G4double P4 = LegendreP(4,cos(theta));
    
  if(sigma == 0)
  {
    A2 = Amax(2,ji,L1,L2,jf,delta);
    A4 = Amax(4,ji,L1,L2,jf,delta);
  }
  else
  {
    A2 = Alpha(2,ji,sigma)*Amax(2,ji,L1,L2,jf,delta);
    A4 = Alpha(4,ji,sigma)*Amax(4,ji,L1,L2,jf,delta);
  }
  
  G4double W = 1 + A2*P2 + A4*P4;

  return W;
}

G4double BeamRequestGammaAndIC::LegendreP(G4int n, G4double x)
{
  G4double P;
  if(n == 0)
  {
    P = 1;
  }
  else if(n == 2)
  {
    P = (1.0/2.0)*(3*(pow(x,2))-1);
  }
  else if(n == 4)
  {
    P = (1.0/8.0)*(35*(pow(x,4))-30*(pow(x,2))+3);
  }
  else
  {
    G4cout << "Error 344546 : Legendre Polynomial Not Found" << G4endl;
    exit(1);
  }

  return P;
}

G4double BeamRequestGammaAndIC::Amax(G4int k, G4double ji, G4int L1, G4int L2, G4double jf, G4double delta)
{
  G4double out;
  if(delta == 0)
  {
    out = (1.0/(1+pow(delta,2)))*(f(k,jf,L1,L1,ji));
  }
  else
  {
    out = (1.0/(1+pow(delta,2)))*(f(k,jf,L1,L1,ji)+2*delta*f(k,jf,L1,L2,ji)+(pow(delta,2))*f(k,jf,L2,L2,ji));
  }
  return out;
}

G4double BeamRequestGammaAndIC::f(G4int k, G4double jf, G4int L1, G4int L2, G4double ji)
{
  G4double out;
  G4double b = B(k,ji);
  if(b == 0)
  {
    return 0; 
  }
  out = b*F(k,jf,L1,L2,ji);
  return out;
}

G4double BeamRequestGammaAndIC::F(G4int k, G4double jf, G4int L1, G4int L2, G4double ji)
{
  G4double out;
  G4double CG = ClebschGordan(L1,1,L2,-1,k,0);
  if(CG == 0)
  {
    return 0;
  }
  G4double W = RacahW(ji,ji,L1,L2,k,jf);
  if(W == 0)
  {
    return 0;
  }
  out = pow((-1),(jf-ji-1))*(pow((2*L1+1)*(2*L2+1)*(2*ji+1),(1.0/2.0)))*CG*W;
  return out;
}

G4double BeamRequestGammaAndIC::B(G4int k, G4double j)
{
  G4double out;
  if(fmod(2*j,2) == 0 || j == 1)        // integral spin
  {    
    out = pow((2*j+1),(1.0/2.0))*(pow(-1,j))*ClebschGordan(j,0,j,0,k,0);
  }
  else                   // half-integral spin
  {
    out = pow((2*j+1),(1.0/2.0))*(pow(-1,(j-1.0/2.0)))*ClebschGordan(j,1.0/2.0,j,-1.0/2.0,k,0);
  } 
  return out;
}


G4double BeamRequestGammaAndIC::ClebschGordan(G4double j1, G4double m1, G4double j2, G4double m2, G4double j, G4double m)
{
  // Check Conditions ////////////////////////////////////////////////////////
  if( 2*j1 != floor(2*j1) || 2*j2 != floor(2*j2) || 2*j != floor(2*j) || 2*m1 != floor(2*m1) || 2*m2 != floor(2*m2) || 2*m != floor(2*m) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 != m)
  {
    //G4cout << "m1 + m2 must equal m." << G4endl;
    return 0;
  }
  if( j1 - m1 != floor ( j1 - m1 ) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if( j2 - m2 != floor ( j2 - m2 ) )
  {
    //G4cout << "2*j2 and 2*m2 must have the same parity" << G4endl;
    return 0;
  }
  if( j - m != floor ( j - m ) )
  {
    //G4cout << "2*j and 2*m must have the same parity" << G4endl;
    return 0;
  }
  if(j > j1 + j2 || j < abs(j1 - j2))
  {
    //G4cout << "j is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m) > j)
  {
    //warning('m is out of bounds." << G4endl;
    return 0 ;
  }
  ////////////////////////////////////////////////////////////////////////////
  G4double term, cg;
  G4double term1 = pow((((2*j+1)/factorial(j1+j2+j+1))*factorial(j2+j-j1)*factorial(j+j1-j2)*factorial(j1+j2-j)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j+m)*factorial(j-m)),(0.5));
  G4double sum = 0;

  for(G4int k = 0 ; k <= 99 ; k++ )
  {
    if( (j1+j2-j-k < 0) || (j-j1-m2+k < 0) || (j-j2+m1+k < 0) || (j1-m1-k < 0) || (j2+m2-k < 0) )
    {
    }
    else
    {
      term = factorial(j1+j2-j-k)*factorial(j-j1-m2+k)*factorial(j-j2+m1+k)*factorial(j1-m1-k)*factorial(j2+m2-k)*factorial(k);
      if((k%2) == 1)
      {
        term = -1*term;
      }
      sum = sum + 1.0/term;
    }
  }

  cg = term1*sum;
  return cg;
  // Reference: An Effective Algorithm for Calculation of the C.G.
  // Coefficients Liang Zuo, et. al.
  // J. Appl. Cryst. (1993). 26, 302-304
}

G4double BeamRequestGammaAndIC::RacahW(G4double a, G4double b, G4double c, G4double d, G4double e, G4double f)
{
  G4double out = pow((-1),(a+b+d+c))*Wigner6j(a,b,e,d,c,f);
  return out;
}

G4double BeamRequestGammaAndIC::Wigner6j(G4double J1, G4double J2, G4double J3, G4double J4, G4double J5, G4double J6)
{
  // error checking //////////////////////////////////////////////////////////
  if(J3 > J1 + J2 || J3 < abs(J1 - J2))
  {
    //G4cout << "first J3 triange condition not satisfied. J3 > J1 + J2 || J3 < abs(J1 - J2)" << G4endl;
    return 0;
  }
  if(J3 > J4 + J5 || J3 < abs(J4 - J5))
  {
    //G4cout << "second J3 triange condition not satisfied. J3 > J4 + J5 || J3 < abs(J4 - J5)" << G4endl;    
    return 0;
  }
  if(J6 > J2 + J4 || J6 < abs(J2 - J4))
  {
    //G4cout << "first J6 triange condition not satisfied. J6 > J2 + J4 || J6 < abs(J2 - J4)" << G4endl;    
    return 0;
  }
  if(J6 > J1 + J5 || J6 < abs(J1 - J5))
  {
    //G4cout << "second J6 triange condition not satisfied. J6 > J1 + J5 || J6 < abs(J1 - J5)" << G4endl;    
    return 0;
  }

  G4double j1 = J1;
  G4double j2 = J2;
  G4double j12 = J3;
  G4double j3 = J4;
  G4double j = J5;
  G4double j23 = J6;
  G4double sum = 0;

  for(G4double m1 = -j1 ; m1 <= j1 ; m1++ )
  {
    for(G4double m2 = -j2 ; m2 <= j2 ; m2++ )
    {
      for(G4double m3 = -j3 ; m3 <= j3 ; m3++ )
      {
        for(G4double m12 = -j12 ; m12 <= j12 ; m12++ )
        {
          for(G4double m23 = -j23 ; m23 <= j23 ; m23++ )
          {
            for(G4double m = -j ; m <= j ; m++ )
            {
              sum = sum + pow((-1),(j3+j+j23-m3-m-m23))*Wigner3j(j1,j2,j12,m1,m2,m12)*Wigner3j(j1,j,j23,m1,-m,m23)*Wigner3j(j3,j2,j23,m3,m2,-m23)*Wigner3j(j3,j,j12,-m3,m,m12);
            }
          }
        }
      }
    }
  }
  return sum;
}

G4double BeamRequestGammaAndIC::Wigner3j(G4double j1, G4double j2, G4double j3, G4double m1, G4double m2, G4double m3)
{
  // error checking //////////////////////////////////////////////////////////
  if( 2*j1 != floor(2*j1) || 2*j2 != floor(2*j2) || 2*j3 != floor(2*j3) || 2*m1 != floor(2*m1) || 2*m2 != floor(2*m2) || 2*m3 != floor(2*m3) )
  {
    G4cout << "All arguments must be integers or half-integers." << G4endl;
    return 0;
  }
  if(m1 + m2 + m3 != 0)
  {
    //G4cout << "m1 + m2 + m3 must equal zero." << G4endl;
    return 0;
  }
  if( j1 + j2 + j3 != floor(j1 + j2 + j3) )
  {
    //G4cout << "2*j1 and 2*m1 must have the same parity" << G4endl;
    return 0;
  }
  if(j3 > j1 + j2 || j3 < abs(j1 - j2))
  {
    //G4cout << "j3 is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m1) > j1)
  {
    //G4cout << "m1 is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m2) > j2)
  {
    //G4cout << "m2 is out of bounds." << G4endl;
    return 0;
  }
  if(abs(m3) > j3)
  {
    //G4cout << "m3 is out of bounds." << G4endl;
    return 0;
  }
  ///////////////////////////////////////////////////////////////////////////

  G4double out = (pow((-1),(j1-j2-m3)))/(pow((2*j3+1),(1.0/2.0)))*ClebschGordan(j1,m1,j2,m2,j3,-1*m3);
  return out;
}

G4double BeamRequestGammaAndIC::factorial(G4double value)
{
  G4double fac;
  if(value > 1)
  {
    fac = value*factorial(value-1);
  }
  else
  {
    fac = 1;
  }
  return fac;
}

G4double BeamRequestGammaAndIC::Alpha(G4int k, G4double ji, G4double sigma)
{
  G4double out;
  G4double rho = Rho(k,ji,sigma);
  if(rho == 0)
  {
    out = 0;
  }
  else
  {
    out = rho/B(k,ji);
  }
  return out;
}

G4double BeamRequestGammaAndIC::MyAlphaM(G4int k, G4double ji, G4double mi)
{
  G4double sum = 0;
  G4double out,b,j;
  G4int mdex, jdex, pop;
  j = abs(ji);
  jdex = floor(j);
  pop = 0;
  if(fmod(2*j,2) == 0 || j == 1)        // integral spin
  {
    if( k == 2)
    {
      sum = 0;
      for(G4double m = -j ; m <= j ; m++)
      {  
        mdex = m+j;
        if(m == mi)
          pop = 1.0;
        else
          pop = 0.0;  
        sum = sum + pow((-1),(j-m))*ClebschGordanIntJK2[jdex][mdex]*pop;
      }
      mdex = 2*j - j;
      b = (pow(-1,j))*ClebschGordanIntJK2[jdex][mdex];      
    }
    else if( k == 4 )
    {
      sum = 0;
      for(G4double m = -j ; m <= j ; m++)
      {  
        mdex = m+j;
        if(m == mi)
          pop = 1.0;
        else
          pop = 0.0;            
        sum = sum + pow((-1),(j-m))*ClebschGordanIntJK4[jdex][mdex]*pop;
      }  
      mdex = 2*j - j;
      b = (pow(-1,j))*ClebschGordanIntJK4[jdex][mdex];        
    }
    else
    {
      G4cout << "Error: 12382388" << G4endl;
      exit(1);
    }
  }
  else                   // half-integral spin
  {
      
    jdex = floor(j);
    if( k == 2)
    {   
      sum = 0; 
      for(G4double m = -j ; m <= j ; m++)
      {  
        mdex = m+j;
        if(m == mi)
          pop = 1.0;
        else
          pop = 0.0;            
        sum = sum + pow((-1),(j-m))*ClebschGordanHalfJK2[jdex][mdex]*pop;
      }
      mdex = ceil(j);
      b = (pow(-1,j-0.5))*ClebschGordanHalfJK2[jdex][mdex];
    }
    else if( k == 4 )
    {
      sum = 0;
      for(G4double m = -j ; m <= j ; m++)
      {  
        mdex = m+j;  
        if(m == mi)
          pop = 1.0;
        else
          pop = 0.0;          
        sum = sum + pow((-1),(j-m))*ClebschGordanHalfJK4[jdex][mdex]*pop;
      }
      mdex = ceil(j);
      b = (pow(-1,j-0.5))*ClebschGordanHalfJK4[jdex][mdex];  
    }
    else
    {
      G4cout << "Error: 342222" << G4endl;
      exit(1);
    }
  } 
  
  
  if(b == 0.0 && sum != 0.0)
  {
    G4cout << "Error: 234324322 b == 0.0 and sum != 0.0" << G4endl;
    exit(1);
  }
  if(b == 0.0)
  {
    out = 0.0;
  }
  else
  {
    out = double(sum/b);
  }  
  
  return out;
}

G4double BeamRequestGammaAndIC::Rho(G4int k, G4double j, G4double sigma)
{
  G4double sum = 0;
  for(G4double m = -j ; m <= j ; m++) {  
    sum = sum + pow((-1),(j-m))*ClebschGordan(j,m,j,-m,k,0)*P(m,j,sigma);
  }
  G4double out = (pow((2*j+1),(1/2)))*sum;
  return out;
}

G4double BeamRequestGammaAndIC::P(G4double m, G4double j, G4double sigma)
{
  G4double sum = 0;
  for(G4double mp = -j ; mp <= j ; mp++) {
    sum = sum + exp((-1)*((pow(mp,2))/(2*pow(sigma,2))));
  }
  G4double out = exp(-(pow(m,2))/(2*pow(sigma,2)))/sum;
  return out;
}

void BeamRequestGammaAndIC::GenerateWThetaAngularDistribution()
{
  G4double ji, jf, m, delta, L1, L2, A2, A4, P2, P4, Alpha2, Alpha4, value, func, funcSum;
  G4double sum = 0;
  for(G4int i = 0; i < ARRAYDIM ; i++) {
    for(G4int j = 0; j < VECSIZE ; j++) {
      for(G4int k = 0; k < SPINDIM ; k++) {
        this->wThetaSum[i][j][k] = 0.0;
        this->wTheta[i][j][k] = 0.0;
      }
    }
  }

  for(G4int i = 0; i < this->gamma.size() ; i++) {
    ji = this->gamma[i].ji;
    jf = this->gamma[i].jf;
    L1 = this->gamma[i].L1;
    L2 = this->gamma[i].L2;
    delta = this->gamma[i].delta;
    A2 = Amax(2,ji,L1,L2,jf,delta);
    A4 = Amax(4,ji,L1,L2,jf,delta);

    if(2*ji > SPINDIM) {
      G4cout << "Error 344112123: 2*ji too big!" << G4endl;
      exit(1);
    }


    for(G4int k = 0; k <= 2*ji ; k++ ) {
      m = (-1*ji)+k;
      Alpha2 = MyAlphaM(2, ji, m);
      Alpha4 = MyAlphaM(4, ji, m);
      funcSum = 0;
      for(G4int l = 0; l < VECSIZE ; l++ ) {
        value = (double(l)/(double(VECSIZE-1)))*M_PI;
        P2 = LegendreP(2,cos(value));
        P4 = LegendreP(4,cos(value));
        func = (M_PI/2.0)*(1 +  Alpha2*A2*P2 +  Alpha4*A4*P4)*sin(value)*(1.0/double(VECSIZE-1));
        funcSum = funcSum + func;

        this->wThetaSum[i][l][k] = funcSum;
        this->wTheta[i][l][k] = func;
        if(l == VECSIZE-1) {
          if(funcSum < 0.9999999999 || funcSum > 1.0000000001) {
            //G4cout << "gamma funcSum =" << setprecision(10) << funcSum << G4endl;
            //G4cout << "Need to renormalize gamma W(theta)" << G4endl;
            sum = funcSum;
            for(G4int n = 0; n < VECSIZE ; n++ )
            {
              this->wThetaSum[i][n][k] = this->wThetaSum[i][n][k]/sum;
              this->wTheta[i][n][k] = this->wTheta[i][n][k]/sum;
            }
          }
        }        
      }
    
      if(writeOutWThetaFiles)
        WriteOutWThetaData(i,k,m);
    }           
  }
}
void BeamRequestGammaAndIC::WriteOutWThetaData(G4int gDex, G4int jDex, G4double m)
{
  G4double theta, ji, jf, L1, L2, delta, A2, A4, Alpha2, Alpha4, energyID;
  string s;
  stringstream out;
  out << gDex;
  s = "wTheta_gamma";
  s.append(out.str());
  out.str("");
  s.append("_j");
  out << jDex;
  s.append(out.str());
  out.str("");
  s.append(".dat");      

  ofstream oo;

  oo.open(s.c_str());
  if( !oo ) { // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    exit(1);
  }
  for(G4int n = 0; n < VECSIZE ; n++ ) {
    theta = (double(n)/(double(VECSIZE-1)))*M_PI;
    oo << theta << "\t" << this->wTheta[gDex][n][jDex] << endl;
  }        
  oo.close(); 

  ofstream app;
  if(this->firstWriteOfWThetaFiles) {
    app.open("gammaKey.dat");
    firstWriteOfWThetaFiles = false;
  }
  else {
    app.open("gammaKey.dat", ios::app);
  }

  if( !app ) { // file couldn't be opened
    cerr << "Error: file could not be opened" << endl;
    exit(1);
  }
  energyID = this->gamma[gDex].energyID;
  ji = this->gamma[gDex].ji;
  jf = this->gamma[gDex].jf;    
  L1 = this->gamma[gDex].L1;
  L2 = this->gamma[gDex].L2;
  delta = this->gamma[gDex].delta;  
  A2 = Amax(2,ji,L1,L2,jf,delta);
  A4 = Amax(4,ji,L1,L2,jf,delta);
  Alpha2 = MyAlphaM(2, ji, m);
  Alpha4 = MyAlphaM(4, ji, m);
  app << "index = " << gDex << "\t" << "energyID = " << energyID << "\t" << "ji = " << ji << "\t" << "jf = " << jf << "\t" << "m = " << m << "\t" << "L1 = " << L1 << "\t" << "L2 = " << L2 << "\t" << "delta = " << setprecision(4) << delta << "\t" << "A2max = " << A2 << "\t" << "A4max = " << A4 << "\t" << "Alpha2 = " << Alpha2 << "\t" << "Alpha4 = " << Alpha4 << "\t" << endl;
  oo.close(); 
}

G4double BeamRequestGammaAndIC::UniformRand48()
{
  G4double rand = drand48();
  return rand;
}
