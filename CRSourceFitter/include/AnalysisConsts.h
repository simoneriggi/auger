#ifndef _AnalysisConsts_H
#define _AnalysisConsts_H

#include "TMath.h"


namespace CRSourceFitter_ns {

static const int fNmassAtSource = 4;
static const int fNmassAtEarth = 4;
//static const int fNmassAtEarth = 26;

static const int fNmassPropData= 7;


//#####################
//### MASS GROUP   ####
//#####################
// SOURCE 
// 1) Z0=1 proton
// 2) Z0=2 He
// 3) Z0=6,8,12  CNO
// 4) Z0=14,26   heavy
  
// EARTH
// 1) Z=0,1
// 2) Z=2
// 3) 4<=Z<=12
// 4) 12<Z<=26

static const float fMinCharge[fNmassAtSource]= {0.0,1.5,3.5,12.5};
static const float fMaxCharge[fNmassAtSource]= {1.5,3.5,12.5,26.0};
static const float fMinMass[fNmassAtSource]= {0.0,3.5,4.5,24.5};
static const float fMaxMass[fNmassAtSource]= {1.5,4.5,24.5,56.0};

static const float fMinChargeAtEarth[fNmassAtEarth]= {0.0,1.5,3.5,12.5};
static const float fMaxChargeAtEarth[fNmassAtEarth]= {1.5,3.5,12.5,26.0};
static const float fMinMassAtEarth[fNmassAtEarth]= {0.0,3.5,4.5,24.5};
static const float fMaxMassAtEarth[fNmassAtEarth]= {1.5,4.5,24.5,56.0};



/*
//## NO GROUPING - FULL MASS/CHARGE RANGE
static const float fMinChargeAtEarth[fNmassAtEarth]= { 0.0,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,
																										13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5};
static const float fMaxChargeAtEarth[fNmassAtEarth]= { 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,
																										14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5};

static const float fMinMassAtEarth[fNmassAtEarth]= { 0.5,3.5,6.5,8.5,10.5,11.5,13.5,15.5,18.5,19.5,22.5,23.5,26.5,27.5,
																										30.5,31.5,34.5,38.5,39.5,39.5,44.5,47.5,50.5,51.5,54.5,55.5};

static const float fMaxMassAtEarth[fNmassAtEarth]= { 1.5,4.5,7.5,9.5,11.5,12.5,14.5,16.5,19.5,20.5,23.5,24.5,27.5,28.5,
																										31.5,32.5,35.5,39.5,40.5,40.5,45.5,48.5,51.5,52.5,55.5,56.5};
*/

//## Z-dependent maximum energy
static const float fMaxGenEnergyAtSourceInSimulation[fNmassAtSource]={log10(pow(10,21.)),log10(2.*pow(10,21.)),log10(12.*pow(10,21.)),log10(26.*pow(10,21.))};

//## A-dependent maximum energy
//static const float fMaxGenEnergyAtSourceInSimulation[fNmassAtSource]={log10(pow(10,21.)),log10(4.*pow(10,21.)),log10(24.*pow(10,21.)),log10(56.*pow(10,21.))};

//############################
//###   Xmax BINNING 
//############################

//## EG binning ## PRL PAPER
static const int fNbins_Xmax= 7;
static const double fEmin_Xmax[fNbins_Xmax]={18.6,18.7,18.8,18.9,19.0,19.2,19.4};
static const double fEmax_Xmax[fNbins_Xmax]={18.7,18.8,18.9,19.0,19.2,19.4,19.8};


//## ENERGY SYST SHIFT +20% (x1.2)
//static const int fNbins_Xmax= 7;
//static const double fEmin_Xmax[fNbins_Xmax]={18.6792,18.7792,18.8792,18.9792,19.0792,19.2792,19.4792};
//static const double fEmax_Xmax[fNbins_Xmax]={18.7792,18.8792,18.9792,19.0792,19.2792,19.4792,19.8792};


//## ENERGY SYST SHIFT -20% (x0.8)
//static const int fNbins_Xmax= 7;
//static const double fEmin_Xmax[fNbins_Xmax]={18.5031,18.6031,18.7031,18.8031,18.9031,19.1031,19.3031};
//static const double fEmax_Xmax[fNbins_Xmax]={18.6031,18.7031,18.8031,18.9031,19.1031,19.3031,19.7031};



//############################
//###   SPECTRUM  BINNING 
//############################
//## Spectrum binning ABOVE ANKLE
static const int fNbins_Spectrum= 19;
static const double fEmin_Spectrum[fNbins_Spectrum]={18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4};
static const double fEmax_Spectrum[fNbins_Spectrum]={18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4,20.5};


//### ENERGY SYST SHIFT +20% (x1.2)
//static const int fNbins_Spectrum= 19;
//static const double fEmin_Spectrum[fNbins_Spectrum]={18.6792,18.7792,18.8792,18.9792,19.0792,19.1792,19.2792,19.3792,19.4792,19.5792,19.6792,19.7792,19.8792,19.9792,20.0792,20.1792,20.2792,20.3792,20.4792};
//static const double fEmax_Spectrum[fNbins_Spectrum]={18.7792,18.8792,18.9792,19.0792,19.1792,19.2792,19.3792,19.4792,19.5792,19.6792,19.7792,19.8792,19.9792,20.0792,20.1792,20.2792,20.3792,20.4792,20.5792};
//##########


//### ENERGY SYST SHIFT -20% (x0.8)
//static const int fNbins_Spectrum= 19;
//static const double fEmin_Spectrum[fNbins_Spectrum]={18.5031,18.6031,18.7031,18.8031,18.9031,19.0031,19.1031,19.2031,19.3031,19.4031,19.5031,19.6031,19.7031,19.8031,19.9031,20.0031,20.1031,20.2031,20.3031};
//static const double fEmax_Spectrum[fNbins_Spectrum]={18.6031,18.7031,18.8031,18.9031,19.0031,19.1031,19.2031,19.3031,19.4031,19.5031,19.6031,19.7031,19.8031,19.9031,20.0031,20.1031,20.2031,20.3031,20.4031};
//##########


//##################################
//###   PROPAGATION MC  BINNING 
//##################################
//## Propagation data binning 
static const int fNbins_PropData= 45;
static const double fEmin_PropData[fNbins_PropData]={18.0,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4,20.5,20.6,20.7,20.8,20.9,21.0,21.1,21.2,21.3,21.4,21.5,21.6,21.7,21.8,21.9,22.0,22.1,22.2,22.3,22.4};
static const double fEmax_PropData[fNbins_PropData]={18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4,20.5,20.6,20.7,20.8,20.9,21.0,21.1,21.2,21.3,21.4,21.5,21.6,21.7,21.8,21.9,22.0,22.1,22.2,22.3,22.4,22.5};


//### ENERGY SYST SHIFT +20% (x1.2)
//static const int fNbins_PropData= 45;
//static const double fEmin_PropData[fNbins_PropData]={18.0792,18.1792,18.2792,18.3792,18.4792,18.5792,18.6792,18.7792,18.8792,18.9792,19.0792,19.1792,19.2792,19.3792,19.4792,19.5792,19.6792,19.7792,19.8792,19.9792,20.0792,20.1792,20.2792,20.3792,20.4792,20.5792,20.6792,20.7792,20.8792,20.9792,21.0792,21.1792,21.2792,21.3792,21.4792,21.5792,21.6792,21.7792,21.8792,21.9792,22.0792,22.1792,22.2792,22.3792,22.4792};
//static const double fEmax_PropData[fNbins_PropData]={18.1792,18.2792,18.3792,18.4792,18.5792,18.6792,18.7792,18.8792,18.9792,19.0792,19.1792,19.2792,19.3792,19.4792,19.5792,19.6792,19.7792,19.8792,19.9792,20.0792,20.1792,20.2792,20.3792,20.4792,20.5792,20.6792,20.7792,20.8792,20.9792,21.0792,21.1792,21.2792,21.3792,21.4792,21.5792,21.6792,21.7792,21.8792,21.9792,22.0792,22.1792,22.2792,22.3792,22.4792,22.5792};
//###

//### ENERGY SYST SHIFT -20% (x0.8)
//static const int fNbins_PropData= 45;
//static const double fEmin_PropData[fNbins_PropData]={17.9031,18.0031,18.1031,18.2031,18.3031,18.4031,18.5031,18.6031,18.7031,18.8031,18.9031,19.0031,19.1031,19.2031,19.3031,19.4031,19.5031,19.6031,19.7031,19.8031,19.9031,20.0031,20.1031,20.2031,20.3031,20.4031,20.5031,20.6031,20.7031,20.8031,20.9031,21.0031,21.1031,21.2031,21.3031,21.4031,21.5031,21.6031,21.7031,21.8031,21.9031,22.0031,22.1031,22.2031,22.3031};
//static const double fEmax_PropData[fNbins_PropData]={18.0031,18.1031,18.2031,18.3031,18.4031,18.5031,18.6031,18.7031,18.8031,18.9031,19.0031,19.1031,19.2031,19.3031,19.4031,19.5031,19.6031,19.7031,19.8031,19.9031,20.0031,20.1031,20.2031,20.3031,20.4031,20.5031,20.6031,20.7031,20.8031,20.9031,21.0031,21.1031,21.2031,21.3031,21.4031,21.5031,21.6031,21.7031,21.8031,21.9031,22.0031,22.1031,22.2031,22.3031,22.4031};
//###


}//close namespace

#endif

