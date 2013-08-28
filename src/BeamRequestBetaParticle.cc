#include "BeamRequestBetaParticle.hh"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits.h>

using namespace std;

// Defines used in Carl's Beta Particle Energy Distrib. Code
// Energies are in keV!
const int MaxN     = 16384;
const int MaxB     = 200;
const double hBarc = 197329.0;
const double Tstep = 1.0;
const double me    = 511.003;
const double ro    = 1.25;
const double alpha = 0.007297;
const double Norm  = 1.0;

BeamRequestBetaParticle::BeamRequestBetaParticle(G4String path)
{
  this->path                 = path;
  this->betaFile             = "betaData.dat";

  // These variables are added if someone wants to hard code
  // in a angular distribution to the beta particles. Like the
  // gamma code, this could be integrated into the messenger.
  this->angularDistribution  = "Isotropic"; // "Isotropic" or "W(theta)"
  this->degreeOfPolarization = 0.0;
  this->betaAsymmetry        = 0.0;
  this->timeT1               = 1.0; // Hz
  this->timeT2               = 1.0; // Hz
  this->stringT1             = "INF";
  this->stringT2             = "INF";
  this->spinTempBeta         = 0.0;

  for(int i = 0 ; i < MAXQVALUE ; i++) {
    for(int j = 0 ; j < 2 ; j++) {
      for(int k = 0 ; k < MAXLEVELS ; k++) {
        this->betaDistrib[i][j][k] = 0.0;
      }
    }
  }
  for(int i = 0 ; i < VECSIZE ; i++) {
    wTheta[i] = 0;
  }  
  GenerateWThetaAngularDistribution();
  BetaParticleDistrib();
}

BeamRequestBetaParticle::~BeamRequestBetaParticle()
{
}

void BeamRequestBetaParticle::BetaParticleDistrib()
{
  // BetaParticleDistrib(); contains Carl's orginal beta energy 
  // distribution code with minimal changes.
  G4int N, i, j, k, Nbranch, Nbig, ansInt;
  G4double Ii, Pi, Idif, Pc, Qvalue, Z, betaBranchToGS;
  G4double Id[MaxB], If[MaxB], Pf[MaxB], Ex[MaxB], Br[MaxB];
  G4double spec[MaxN], specT[MaxN];
  G4double Q, T, E, p, q, gam, y, phi, R, F, pi, Area, Shape;
  G4double variBranches;    // for probability distrib.
  G4double variTotal;       // for probability distrib.
  G4double A;
  
  G4int ncolumns = 0;
  G4int nlines = 0;
  G4int intBuf;
  G4double floatBuf;
  string buf, currentline;

  G4String myFile = this->path + this->betaFile;

  Nbranch = 0;

  ifstream infile;
  infile.open(myFile.c_str());

  if (!infile) {
    cerr << "unable to find/open file" << endl;
    exit (1);
  }

  while(getline(infile, currentline)) {
    if (!currentline.empty()) {
      stringstream ss(currentline);
      ncolumns=0;
      while (ss >> buf) {
        intBuf = atoi( buf.c_str() );
        floatBuf = atof( buf.c_str() );
        if(nlines == 0) {
          if(ncolumns == 0) {
            Z = intBuf;
          }
          else if(ncolumns == 1) {
            A = intBuf;
          }                              
        }  
        else if(nlines == 1) {
          if(ncolumns == 0) {
            Ii = floatBuf;
          }
          else if(ncolumns == 1) {
            Pi = intBuf;          
          }                              
        }
        else if(nlines == 2) {
          if(ncolumns == 0) {
            Qvalue = floatBuf;
          }
          if(ncolumns == 1) {
            betaBranchToGS = floatBuf;
          }          
        }        
        else if(nlines > 2) {
          if(ncolumns == 0) {
            Id[Nbranch] = floatBuf;
          }
          else if(ncolumns == 1) {
            If[Nbranch] = floatBuf;
          }
          else if(ncolumns == 2) {
            Pf[Nbranch] = floatBuf; 
          }     
          else if(ncolumns == 3) {
            Ex[Nbranch] = floatBuf; 
          } 
          else if(ncolumns == 5) {
            Br[Nbranch] = floatBuf; 
          }                     
        }          
        else {
          G4cout << "Error: 1324431" << G4endl;
          exit (1);
        }
        ncolumns++;
      }
      if(nlines > 2) {
        Nbranch++;
      }          
      nlines++;
    }
  }

  this->myZ = abs(Z);

  infile.close(); 

  if (Z < 0) {
    this->betaPlusMinusType = 1;  // positrons
  }
  else {
    this->betaPlusMinusType = 0;  // electrons
  }

  if(Nbranch >= MaxB) {
    G4cout << "Sorry, maximum number of branches is " << MaxB-1 << G4endl;
    exit(1);
  }

  pi = acos(-1.0);
  Nbig = 0;

  for(i=0; i < MaxN; i++) {
    specT[i] = 0.0;
  }

  for(j=0; j < Nbranch; j++) {
    Q = Qvalue - Ex[j];
    Idif = fabs(Ii - If[j]);
    Pc = Pi*Pf[j];

    for(i=0; i < MaxN; i++) {
      spec[i] = 0.0;
    }   

    i = 0;
    T = 0.5*Tstep;
    Area = 0.0;   

    while(T < Q) {
      E = T + me;
      p = sqrt(E*E - me*me);
      q = Q - T;

      // Basic phase space contribution to spectrum shape 

      spec[i] = p*E*q*q;

      // Fermi Function approximation given on page 682 of Blatt and Weisskopf    
      // With a basic size correction L(Z,E) correction: 2*(1+gam) rather than 4. 
      // For exact complex gamma-function evaluation see Numerical Recipes.       
      // No electronic screening correction is applied.  This makes little        
      // difference for electrons, where there is only a 5% error for T = 10 keV  
      // and Z = 95.  However, it is a huge effect for low-energy positrons,      
      // changing F by a factor of 735 for T = 10 keV and Z = -95.  See Table 2   
      // in Buhring, NPA 61, 110 (1965) and Fig 1 of Wilkinson NPA 262, 58 (1974) 
      // for the effect on the integrated f-value as a function of beta endpoint. 
      // One can either get the numerical tables of Behrens and Janecke,         
      // Landolt-Bornstein Group I, Vol 4, 1969 and fit the tabulated data,       
      // QC 61.L332 GRP.1.BD.4 in the Reference Section), or                     
      // once I have the exact magnitude of the complex gamma-function evaluated, 
      // implement all of the analytic corrections described in Wilkinson Part II 
      // NIM A 335, 12 (1993).  For now, do not trust the positron spectrum shape 
      // below about 50 keV for intermediate mass nuclei.                         

      gam = sqrt(1.0 - alpha*alpha*Z*Z);
      y = alpha*Z*E/p;
      phi = atan(gam/y);
      if(phi < 0) {phi += pi;}
      R   = ro*pow(A,(1.0/3.0));
      F   = 2.0*(gam + 1.0)*pow((2.0*gam + 1.0),-1.0*(4.0*gam + 1.0))
          * pow((gam*gam + y*y),(gam - 0.5))
          * pow((2.0*p*R/hBarc),(2.0*gam - 2.0))
          * exp(2.0*phi*y + 2.0*(gam + 1.0) + (1.0/6.0)*( (gam/(gam*gam + y*y))
          - (1.0 / (2.0*gam + 1.0)) ) );
          
      spec[i] *= F;

      // Determine the Shape Factor 
      // These assume that the centrifugal barrier is small compared to the          
      // Coulomb barrier, valid for (Ze^2/hbarv)^2 << 1, i.e. all but very low       
      // energy betas.  For the allowed and unique transitons the shape factors      
      // are trustworthy.  For the non-unique transitions, the analytic expression   
      // also requires, Ze^2/R << 2p, a dubious assumption for high Q-value decays   
      // and does NOT account for any interference between the opposite-parity       
      // l-1 relativitic matrix elements and the normal-parity order-l matrix        
      // elements.  The relative signs and magnitudes of these matrix elements can   
      // only be determined by detailed calculation, and I have arbitrarily          
      // assumed that the order l normal contributions dominate.  These cases should 
      // thus not be trusted if these contribute significant decay branches.         

      if(Idif < 1.1) {
        if(Pc > 0.0) { // Allowed
          Shape = 1.0;
          //G4cout << "Shape - Allowed" << G4endl;
        }
        else {
          Shape = 1.0;
          if(i == 0) {
            G4cout << "Shape - Warning, Branch " << j+1 << " is First-Forbidden Non-Unique!" << G4endl;
          }
        }
      }
      else if(Idif < 2.1) {
        if(Pc < 0.0) { // First-Forbidden Unique
          Shape = p*p + q*q;
          //G4cout << "Shape - First-Forbidden Unique" << G4endl;
        }
        else {
          Shape = p*p + 4.0*q*q;
          if(i == 0) {
            G4cout << "Shape - Warning, Branch " << j+1 << " is Second-Forbidden Non-Unique!" << G4endl;
          }
        }
      }
      else if(Idif < 3.1) {
        if(Pc > 0.0) { // Second-Forbidden Unique
          Shape = 3.0*p*p*p*p + 10.0*p*p*q*q  + 3.0*q*q*q*q;
          //G4cout << "Shape - Second-Forbidden Unique" << G4endl;
        } 
        else {
          Shape = p*p*p*p + 15.0*p*p*q*q + 18.0*q*q*q*q;
          if(i == 0) {
            G4cout << "Shape - Warning, Branch " << j+1 << " is Third-Forbidden Non-Unique!" << G4endl;
          }
        }
      }
      else if(Idif < 4.1) {
        if(Pc < 0.0) { // Third-Forbidden Unique
          Shape = 4.0*p*p*p*p*p*p +  28.0*p*p*p*p*q*q +  28.0*p*p*q*q*q*q +   4.0*q*q*q*q*q*q;
          //G4cout << "Shape - Third-Forbidden Unique" << G4endl;
        } 
        else {
          Shape = 9.0*p*p*p*p*p*p + 112.0*p*p*p*p*q*q + 252.0*p*p*q*q*q*q + 144.0*q*q*q*q*q*q;
          if(i == 0)
          {
            G4cout << "Shape - Warning, Branch " << j+1 << " is Fourth-Forbidden Non-Unique!" << G4endl;
          }
        }
      }
      else if(Idif < 5.1) {      
        if(Pc > 0.0) { // Fourth-Forbidden Unique
          Shape = 5.0*p*p*p*p*p*p*p*p + 60.0*p*p*p*p*p*p*q*q + 126.0*p*p*p*p*q*q*q*q + 60.0*p*p*q*q*q*q*q*q + 5.0*q*q*q*q*q*q*q*q;
          //G4cout << "Shape - Fourth-Forbidden Unique" << G4endl;
        } 
        else {
            Shape = 4.0*p*p*p*p*p*p*p*p + 75.0*p*p*p*p*p*p*q*q + 280.0*p*p*p*p*q*q*q*q + 300.0*p*p*q*q*q*q*q*q + 100.0*q*q*q*q*q*q*q*q;
            if(i == 0)
            {
              G4cout << "Shape - Warning, Branch " << j+1 << " is Sixth-Forbidden Non-Unique!" << G4endl;
            }
        }
      }
      else if(Idif < 6.1) {
        if(Pc < 0.0) { // Fifth-Forbidden Unique
          Shape = 6.0*p*p*p*p*p*p*p*p*p*p + 110.0*p*p*p*p*p*p*p*p*q*q + 396.0*p*p*p*p*p*p*q*q*q*q + 396.0*p*p*p*p*q*q*q*q*q*q + 110.0*p*p*q*q*q*q*q*q*q*q + 6.0*q*q*q*q*q*q*q*q*q*q;
          //G4cout << "Shape - Fifth-Forbidden Unique" << G4endl;
        }
        else {
            Shape = 10.0*p*p*p*p*p*p*p*p*p*p + 264.0*p*p*p*p*p*p*p*p*q*q + 1485.0*p*p*p*p*p*p*q*q*q*q + 2640.0*p*p*p*p*q*q*q*q*q*q + 1650.0*p*p*q*q*q*q*q*q*q*q + 360.0*q*q*q*q*q*q*q*q*q*q;
            if(i == 0)
            {
              G4cout << "Shape - Warning, Branch " << j+1 << " is Seventh-Forbidden Non-Unique!" << G4endl;
            }
        }
      }
      else if(Idif < 7.1) {
          if(Pc > 0.0) { // Sixth-Forbidden Unique
          Shape = 7.0*p*p*p*p*p*p*p*p*p*p*p*p + 182.0*p*p*p*p*p*p*p*p*p*p*q*q + 1001.0*p*p*p*p*p*p*p*p*q*q*q*q + 1716.0*p*p*p*p*p*p*q*q*q*q*q*q + 1001.0*p*p*p*p*q*q*q*q*q*q*q*q + 182.0*p*p*q*q*q*q*q*q*q*q*q*q + 7.0*q*q*q*q*q*q*q*q*q*q*q*q;
          //G4cout << "Shape - Sixth-Forbidden Unique" << G4endl;
        }
        else {
            Shape = 900.0*p*p*p*p*p*p*p*p*p*p*p*p + 31850.0*p*p*p*p*p*p*p*p*p*p*q*q + 252252.0*p*p*p*p*p*p*p*p*q*q*q*q + 675675.0*p*p*p*p*p*p*q*q*q*q*q*q + 700700.0*p*p*p*p*q*q*q*q*q*q*q*q + 286650.0*p*p*q*q*q*q*q*q*q*q*q*q + 44100.0*q*q*q*q*q*q*q*q*q*q*q*q;
            if(i == 0)
            {
              G4cout << "Shape - Warning, Branch " << j+1 << " is Eighth-Forbidden Non-Unique!" << G4endl;
            }
        }
      }
      else {
            G4cout << "Shape - Sorry, Branch " << j+1 << " has Idif = " << Idif << ", Pc = " << Pc << G4endl;
        exit(1);
      }

      spec[i] *= Shape;
      // Not multiplied by Tstep here to give Area = 10^6 x Tstep
      Area += spec[i];

      T += Tstep;
      i += 1;
    }

    N = i;
    if(N > Nbig) {Nbig = N;}
    //Area /= Br[j]*Norm;

    for(i = 0; i < N; i++) {
      T = ((float)i + 0.5)*Tstep;
      variBranches = 0;
      for(k = 0; k < i; k++) {
        variBranches += spec[k]/Area;
      }

      // if(maxBranch < j){maxBranch = j;}
      this->betaDistrib[i][0][j] = T;
      this->betaDistrib[i][1][j] = variBranches;
      specT[i] += spec[i]/Area;
    }
    //fclose(ofp);
  }
  this->maxBranch            = Nbranch;
  this->numberOfBetaBranches = Nbranch;

  if(Nbranch > MAXBETABRANCHES) {
    G4cout << "Error 2407 : too many beta branches" << G4endl;
    exit(1);
  }

  for(i = 0; i < Nbig; i++) {
    T = ((float)i + 0.5)*Tstep;
    variTotal = 0;
    for(k = 0; k < i; k++) {
      variTotal += specT[j];
    }
  }
}

void BeamRequestBetaParticle::AddBeamParticle(G4ParticleGun* beamParticle, G4int betaBranch, G4double timeSeconds)
{
  beamParticle -> SetParticleDefinition(G4ParticleTable::GetParticleTable() -> FindParticle( GetParticleName() ));
  beamParticle -> SetParticleEnergy( (GetParticleEnergy(betaBranch)/1000.0) );  // in MeV!
  beamParticle -> SetParticlePosition( GetParticlePosition() );
  beamParticle -> SetParticleMomentumDirection( GetParticleDirection(timeSeconds) ); 
}

G4String BeamRequestBetaParticle::GetParticleName()
{
  G4String particleName;

  if(this->betaPlusMinusType == 0) {
    particleName = "e-";    // electrons
  }
  else if(this->betaPlusMinusType == 1) {
    particleName = "e+";  // positrons
  }
  else {
    G4cout << "Error 17930 : can't determine beta type" << G4endl;
    exit(1);
  }

  return particleName;
}


G4double BeamRequestBetaParticle::GetParticleEnergy(G4int betaBranch)
{
  G4int select;
  G4double particleEnergy;
  G4double rand = UniformRand48();

  if (betaBranch == -1) {
    G4cout << "Error 09374 : beta branch never passed" << G4endl;
    exit(1);
  }

  if (betaBranch <= this->numberOfBetaBranches) {
    for (int y = 0; y < MAXQVALUE; y++) {
      if (y == 0) {
        if (rand > 0 && rand < this->betaDistrib[y][1][betaBranch]) {
          select = y;
          particleEnergy = this->betaDistrib[y][0][betaBranch];
          break;
        }
      }
      else {
        if (rand > this->betaDistrib[y-1][1][betaBranch]  && rand < this->betaDistrib[y][1][betaBranch]) {
          select = y;
          particleEnergy = this->betaDistrib[y][0][betaBranch];
          break;
        }
      }
    }        
  }
  else {
    G4cout << "Error 11034 : beta branch is too large" << G4endl;
    exit(1);
  }

  return particleEnergy;
}

G4ThreeVector BeamRequestBetaParticle::GetParticlePosition()
{
  G4ThreeVector particlePosition = G4ThreeVector(0, 0, 0);
  return particlePosition;
}

G4ThreeVector BeamRequestBetaParticle::GetParticleDirection(G4double timeSeconds)
{
  G4ThreeVector particleDirection;
  G4double rand, phi, theta, r, polarization, beta;
    
  if(this->angularDistribution == "Isotropic") {
    // Randomize the direction over 4pi
    G4double costheta = 2.*UniformRand48()-1.0;
    G4double sintheta = sqrt( 1. - costheta*costheta );
    G4double phi      = (360.*deg)*UniformRand48();
    particleDirection = G4ThreeVector(sintheta*cos(phi), sintheta*sin(phi), costheta);
  }
  else if(this->angularDistribution == "W(theta)") {
    beta = this->spinTempBeta;
    //pSpin = this->level[dex].spin;
    polarization = (exp(beta)-1)/(exp(beta)+1);
    if(stringT1 == "INF") {
      polarization = 1.0;
    }
    else {
      polarization = polarization * exp(-1*(timeSeconds/this->timeT1));
    }
    if(stringT2 == "INF") {
      polarization = 1.0;
    }
    else {
      polarization = polarization * exp(-1*(timeSeconds/this->timeT2));
    }      
    
    if(polarization != 1.0) {
      beta = log((1+polarization)/(1-polarization));
    }

    rand = UniformRand48();
    r    = UniformRand48();
    if(rand < polarization*this->betaAsymmetry) { // wTheta Component
      theta = -1.0;
      for (G4int y = 0; y < VECSIZE; y++) {
        if(r < wTheta[y]) {
          theta = (double(y)/(double(VECSIZE-1)))*M_PI;
          break;
        }
      }
      if(theta == -1.0) {
        G4cout << "Error 32432999 : No Theta Found" << G4endl;
        exit(1);    
      }      
    }
    else { // Isotropic Component
      theta = acos(1.0 - (r*(1.0 - cos(M_PI))));    
    }

    phi = UniformRand48()*2.0*M_PI;
    
    particleDirection = G4ThreeVector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  }
  else {
    G4cout << "Error 23123213 : Beta Angular Distribution" << G4endl;
    exit(1);
  }

  return particleDirection;
}

void BeamRequestBetaParticle::GenerateWThetaAngularDistribution()
{
  G4double value, func, funcSum, sum;
  
  funcSum = 0;
  for (G4int i = 0; i < VECSIZE; i++) {
    value = (double(i)/(double(VECSIZE-1)))*M_PI;
    func = (M_PI/2.0)*(1 + cos(value) )*sin(value)*(1.0/double(VECSIZE-1));
    funcSum = funcSum + func;
    wTheta[i] = funcSum;
    if(i == VECSIZE-1) {
      if (funcSum < 0.9999999 || funcSum > 1.0000001) {
        G4cout << "beta funcSum =" << funcSum << G4endl;
        G4cout << "NEED TO RENORMALIZE BETA W(THETA)! CHECK THIS!" << G4endl;
        sum = funcSum;
        for (G4int k = 0; k < VECSIZE ; k++ ) {
          wTheta[k] = wTheta[k]/sum;
        }
      }  
    }
  }
}

G4double BeamRequestBetaParticle::UniformRand48()
{
  G4double rand = drand48();
  return rand;
}

