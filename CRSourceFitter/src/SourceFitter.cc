/**
* @file SourceFitter.cc
* @class SourceFitter
* @brief Perform a likelihood fit of spectrum&composition data, determining the source parameters (gamma,Emax,abundances,...) 
* 
* @author S. Riggi
* @date 22/04/2010
*/

#include "SourceFitter.h"
#include "DataReader.h"
#include "PropagationMCReader.h"
#include "AnalysisConsts.h"
#include "ConfigParser.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TRandom.h>

#include<TFractionFitter.h>
#include<TPaveText.h>
#include<TVirtualFitter.h>
#include<TObjArray.h>
#include<TMatrixD.h>

#include <TMinuit.h>


#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <cmath>


#include<vector>


using namespace std;


ClassImp(CRSourceFitter_ns::SourceFitter)

namespace CRSourceFitter_ns {

int SourceFitter::fVerbosity;
int SourceFitter::fNmass;

std::vector<TH1D*> SourceFitter::fXmaxData;
std::vector<TH1D*> SourceFitter::fEnergyData;

std::vector< std::vector<TH1D*> > SourceFitter::fXmaxMC;
std::vector< std::vector<TH1D*> > SourceFitter::fXmaxMC_scaled;
std::vector< std::vector<TH1D*> > SourceFitter::fXmaxGenMC;
	
std::vector< std::vector<TH1D*> > SourceFitter::fEnergyMC;
std::vector< std::vector<TH1D*> > SourceFitter::fGenEnergyMC;

TH1D* SourceFitter::fAugerEnergySpectrum;
TH1D* SourceFitter::fAugerEnergySpectrumE3;
TH1D* SourceFitter::fAugerEnergySpectrumEvents;
std::vector<double> SourceFitter::fExposure;

std::vector<TH1D*> SourceFitter::fPropEnergySpectrum_spectrumfit;
std::vector<TH1D*> SourceFitter::fPropEnergySpectrumE3_spectrumfit;
std::vector<TH1D*> SourceFitter::fPropDirectEnergySpectrum_spectrumfit;
std::vector<TH1D*> SourceFitter::fPropDirectEnergySpectrumE3_spectrumfit;
std::vector<TH1D*> SourceFitter::fPropEnergySpectrum_Xmaxfit;
std::vector<TH1D*> SourceFitter::fPropEnergySpectrumE3_Xmaxfit;
TH1D* SourceFitter::fTotPropEnergySpectrum_spectrumfit;
TH1D* SourceFitter::fTotPropEnergySpectrumE3_spectrumfit;
TH1D* SourceFitter::fGalacticSpectrum_spectrumfit;

std::vector<double> SourceFitter::fStartFractParamsAtSource;
std::vector<double> SourceFitter::fFractionsAtSource;
std::vector<double> SourceFitter::fFractionsAtSourceLowErr;
std::vector<double> SourceFitter::fFractionsAtSourceHighErr;
std::vector< std::vector<double> > SourceFitter::fFractionsAtEarth;
std::vector< std::vector<double> > SourceFitter::fFractionsAtEarth_spectrumfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtEarth_Xmaxfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtEarthHighErr_Xmaxfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtEarthLowErr_Xmaxfit;
std::vector< std::vector<double> > SourceFitter::fDirectFractionsAtEarth_spectrumfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtDetector_Xmaxfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtDetectorHighErr_Xmaxfit;
std::vector< std::vector<double> > SourceFitter::fFractionsAtDetectorLowErr_Xmaxfit;


double SourceFitter::fInjIndex;
double SourceFitter::fEmax;
double SourceFitter::fRedShiftEvolution;
double SourceFitter::fSpectrumLikelihood;
double SourceFitter::fXmaxLikelihood;
double SourceFitter::fTotLikelihood;

PropagationMCReader* SourceFitter::fPropMCReader;
DataReader* SourceFitter::fDataReader;

double SourceFitter::fInjIndex_Xmax;
double SourceFitter::fInjIndexErr_Xmax;  
double SourceFitter::fInjIndex_spectrum;
double SourceFitter::fInjIndexErr_spectrum;

double SourceFitter::fEmax_Xmax;
double SourceFitter::fEmaxErr_Xmax; 
double SourceFitter::fEmax_spectrum;
double SourceFitter::fEmaxErr_spectrum;

double SourceFitter::fEmax_start;
double SourceFitter::fInjIndex_start; 
double SourceFitter::fFluxNormFactor_start; 
double SourceFitter::fGalacticEmax_start;
double SourceFitter::fGalacticInjIndex_start; 
double SourceFitter::fGalacticFluxNormFactor_start; 

double SourceFitter::fFluxNormFactor_Xmax;
double SourceFitter::fFluxNormFactorErr_Xmax;
double SourceFitter::fFluxNormFactor_spectrum;
double SourceFitter::fFluxNormFactorErr_spectrum;
double SourceFitter::fFluxNormFactor;
double SourceFitter::fFluxNormFactorErr;

double SourceFitter::fGalacticInjIndex_spectrum;
double SourceFitter::fGalacticInjIndexErr_spectrum;
double SourceFitter::fGalacticEmax_spectrum;
double SourceFitter::fGalacticEmaxErr_spectrum;
double SourceFitter::fGalacticFluxNormFactor_spectrum;
double SourceFitter::fGalacticFluxNormFactorErr_spectrum;

int SourceFitter::fSourceCutoffMode;
double SourceFitter::fSourceCutoffShape;
double SourceFitter::fSourceCutoffShapeErr;
bool SourceFitter::fEmaxADependent;
int SourceFitter::fNdof_spectrum;
int SourceFitter::fNdof_Xmax;
double SourceFitter::fNevents_spectrum;
double SourceFitter::fNevents_Xmax;


TTree* SourceFitter::fXmaxMinimizationInfo;
TTree* SourceFitter::fSpectrumMinimizationInfo;

int SourceFitter::fFitStatus;

//### FIT OPTIONS ###
bool SourceFitter::fIsCombinedFit;
bool SourceFitter::fFitExtragalFractions; 
bool SourceFitter::fFitGalSpectrum;
bool SourceFitter::fFitGalFractions;
bool SourceFitter::fFitSourceCutoffShape;

bool SourceFitter::fFixShift;
bool SourceFitter::fFixMassFractions;
bool SourceFitter::fFixEmax;
bool SourceFitter::fFixInjIndex;
bool SourceFitter::fFixSourceCutoffShape;

bool SourceFitter::fWeightLikelihood;
double SourceFitter::fSpectrumLikelihoodWeight;
double SourceFitter::fXmaxLikelihoodWeight;


SourceFitter::SourceFitter(){
  
  
	fVerbosity= 2;
 
	
 	//initialize fit options
	fEmaxADependent= false;

	fIsCombinedFit= false;	
	fFitExtragalFractions= false;
	fFitGalFractions= false;
	fFitGalSpectrum= false;
	fFitSourceCutoffShape= false;

	
	fFixShift= false;
	fFixMassFractions= false;
	fFixEmax= false;
	fFixInjIndex= false;
	fFixSourceCutoffShape= false;
	
	fWeightLikelihood= false;
	fSpectrumLikelihoodWeight= 1.;
	fXmaxLikelihoodWeight= 1.;  

 
	fNdof_spectrum=0;
	fNdof_Xmax=0;
  fNevents_spectrum=0;
	fNevents_Xmax=0;

  
	//initialize pointers
  fDataReader= NULL;
  
	theMassFractStruct.fSourceFraction = NULL;
	theStartMassFractStruct.fSourceFraction = NULL;

}//close constructor



SourceFitter::~SourceFitter(){
	
	cout<<"~SourceFitter()"<<endl;
	
	
	if(fPropMCReader) {
		cout<<"deleting fPropMCReader"<<endl;
		delete fPropMCReader;
		cout<<"done"<<endl;
	}
	if(fDataReader) {
		cout<<"deleting fDataReader"<<endl;
		delete fDataReader;
		cout<<"done"<<endl;
	}

	
	

	cout<<"deleting theMassFractStruct.fSourceFraction"<<endl;
	delete [] theMassFractStruct.fSourceFraction;
	theMassFractStruct.fSourceFraction = NULL;
	cout<<"done"<<endl;

	cout<<"deleting theStartMassFractStruct.fSourceFraction"<<endl;
  delete [] theStartMassFractStruct.fSourceFraction;
	theStartMassFractStruct.fSourceFraction = NULL;
	cout<<"done"<<endl;
	
}//close destructor


void SourceFitter::GetConfig(){

	//## Get steering info from ConfigParser

	//## FileNames
	fDataFileName= ConfigParser::fDataFileName;
	fSpectrumTableFileName= ConfigParser::fSpectrumTableFileName;
	fMCFileName= ConfigParser::fMCFileName;
	fMatrixFileName= ConfigParser::fMatrixFileName;
	fOutputFileName= ConfigParser::fOutputFileName;

	fNmass= fNmassAtEarth;
 
	//## Source vars
  fInjIndex= ConfigParser::fInjIndex;
	fEmax= ConfigParser::fEmax;
	fFluxNormFactor= ConfigParser::fFluxNorm; 
  fRedShiftEvolution= ConfigParser::fRedshiftEvolution;
	fFractionsAtSource= ConfigParser::fFract;

	fStartFractParamsAtSource.clear();
	fStartFractParamsAtSource.resize(fFractionsAtSource.size()-1);	
	double denom= 1;
	fStartFractParamsAtSource[0]= fFractionsAtSource[0];
	for(unsigned k=1;k<fStartFractParamsAtSource.size();k++) {
		denom-= fFractionsAtSource[k-1];
		fStartFractParamsAtSource[k]= fFractionsAtSource[k]/denom;
	}

	fEmaxADependent= ConfigParser::fEmaxADependent;
	
  //## SET FIT OPTIONS
	fIsCombinedFit= ConfigParser::fCombinedFit;
  fFitExtragalFractions= ConfigParser::fFitExtragalFractions; 
  fFitGalFractions= ConfigParser::fFitGalFractions;  
  fFitGalSpectrum= ConfigParser::fFitGalSpectrum;
	fFitSourceCutoffShape= ConfigParser::fFitSourceCutoffShape;
	
	fFixMassFractions= ConfigParser::fFixMassFractions;
	fFixEmax= ConfigParser::fFixEmax;
	fFixInjIndex= ConfigParser::fFixInjIndex;
	fFixSourceCutoffShape= ConfigParser::fFixSourceCutoffShape;

	fWeightLikelihood= ConfigParser::fWeightLikelihood;
	fSpectrumLikelihoodWeight= ConfigParser::fSpectrumLikelihoodWeight;
	fXmaxLikelihoodWeight= ConfigParser::fXmaxLikelihoodWeight;
	
	fSourceCutoffMode= ConfigParser::fSourceCutoffMode;
	fSourceCutoffShape= ConfigParser::fSourceCutoffShape;


}//close SourceFitter::SetConfig()


void SourceFitter::SetData(){

	cout << "SourceFitter::SetData() - Reading data from " << fDataFileName.c_str() << endl;
    
	//## Create Data Reader class
	fDataReader= new DataReader();		
	fDataReader->SetDataFileName(fDataFileName);
	fDataReader->SetMCFileName(fMCFileName);
	fDataReader->SetSpectrumTableFileName(fSpectrumTableFileName);
	fDataReader->SetOutputFileName(std::string("OutputData.root"));

	fDataReader->Init();//allocate (just once!) all histograms
	fDataReader->ReadSpectrumData();//Read spectrum table
	fDataReader->ReadData();//Read data histos
	fDataReader->ReadMC();//Read data histos
	fDataReader->Save();//Save to file
		
	//## Create PropagationMCReader class
	fPropMCReader= new PropagationMCReader();
	fPropMCReader->SetOutputFileName(std::string("ProcessedPropData.root"));
	
	fPropMCReader->SetInjectionIndex(fInjIndex);
	fPropMCReader->SetMaxInjEnergyAtSource(fEmax);
	fPropMCReader->SetExtraGalSpectrumNormalization(pow(10,fFluxNormFactor));
	fPropMCReader->SetAbundanceAtSource(fFractionsAtSource);//how many masses? more than the fitted at Earth (p,He,CNO,Fe)?
	fPropMCReader->SetRedShiftEvolution(fRedShiftEvolution);
	fPropMCReader->SetSourceCutoffMode(fSourceCutoffMode);
	fPropMCReader->SetSourceCutoffShape(fSourceCutoffShape);
 	fPropMCReader->SetEmaxADependence(fEmaxADependent); 	 	
	
	cout<<"SourceFitter::SetData(): Setting propagation matrix from file "<<fMatrixFileName.c_str()<<endl;
	TFile* MatrixFile= new TFile(fMatrixFileName.c_str(),"READ");
	TMatrixD* TransferMatrix= (TMatrixD*)MatrixFile->Get("TMatrixT<double>");
	fPropMCReader->SetTransferMatrix(TransferMatrix);
	
	fPropMCReader->Init();
	fPropMCReader->CalculateExpectedSpectrumAtEarth();
	fPropMCReader->GenerateGalacticSpectrum();
	fPropMCReader->CalculateMassFractions();

	
	//## Get data histo
	fXmaxData= fDataReader->GetXmaxData();
	fEnergyData= fDataReader->GetEnergyData();
	
	fAugerEnergySpectrum= fDataReader->GetAugerEnergySpectrum();
  fAugerEnergySpectrumE3= fDataReader->GetAugerEnergySpectrumE3();
	fAugerEnergySpectrumEvents= fDataReader->GetAugerEnergySpectrumEvents();
	fExposure= fDataReader->GetExposure();

	//## Get MC histo
	for(int i=0;i<fNmassAtEarth;i++){
		fXmaxMC[i]= fDataReader->GetXmaxMC(i);
		fXmaxGenMC[i]= fDataReader->GetGenXmaxMC(i);
		fEnergyMC[i]= fDataReader->GetEnergyMC(i);	
		fGenEnergyMC[i]= fDataReader->GetGenEnergyMC(i);
	}//end loop masses at Earth

	//## Get PropagationMC histo
	fPropEnergySpectrum_spectrumfit= fPropMCReader->GetExpEnergySpectrumAtEarth_Spectrum();
	fPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpEnergySpectrumE3AtEarth_Spectrum();
  fPropDirectEnergySpectrum_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumAtEarth_Spectrum();
	fPropDirectEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumE3AtEarth_Spectrum(); 
	fGalacticSpectrum_spectrumfit= fPropMCReader->GetGalacticSpectrum_Spectrum();

	fTotPropEnergySpectrum_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumAtEarth_Spectrum();
	fTotPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumE3AtEarth_Spectrum();

	fFractionsAtEarth_spectrumfit= fPropMCReader->GetExpFractionsAtEarth_Spectrum();
	fFractionsAtEarth_Xmaxfit= fPropMCReader->GetExpFractionsAtEarth_Xmax();	
	fDirectFractionsAtEarth_spectrumfit= fPropMCReader->GetExpDirectFractionsAtEarth_Spectrum();


	cout<<"SourceFitter::SetData ==> "<<fXmaxData.size()<<" fit histograms stored "<<endl;

}//close SourceFitter::SetData()



void SourceFitter::Init(){

	//## Get steering config from ConfigParser
	GetConfig();

	//## Inizialize vector histo
	for(int i=0;i<fNmassAtEarth;i++) {
		fXmaxMC.push_back ( std::vector<TH1D*>() );
		fXmaxMC_scaled.push_back ( std::vector<TH1D*>() );
		fXmaxGenMC.push_back ( std::vector<TH1D*>() );
    
		fEnergyMC.push_back ( std::vector<TH1D*>() );
		fGenEnergyMC.push_back ( std::vector<TH1D*>() );
		    
    fFractionsAtEarth.push_back ( std::vector<double>() );
		fFractionsAtEarth_spectrumfit.push_back ( std::vector<double>() );
		fDirectFractionsAtEarth_spectrumfit.push_back ( std::vector<double>() );
		fFractionsAtEarth_Xmaxfit.push_back ( std::vector<double>() );
		fFractionsAtEarthHighErr_Xmaxfit.push_back ( std::vector<double>() );
		fFractionsAtEarthLowErr_Xmaxfit.push_back ( std::vector<double>() );		

		fFractionsAtDetector_Xmaxfit.push_back ( std::vector<double>() );
		fFractionsAtDetectorLowErr_Xmaxfit.push_back ( std::vector<double>() );
		fFractionsAtDetectorHighErr_Xmaxfit.push_back ( std::vector<double>() );

		fFractionsAtSourceLowErr.push_back(0.);
		fFractionsAtSourceHighErr.push_back(0.);

		fMeanXmaxMC.push_back ( vector<double>() );
  	fMeanXmaxErrMC.push_back ( vector<double>() );
  	fRMSMC.push_back ( vector<double>() );
  	fRMSErrMC.push_back ( vector<double>() );
  	fMeanEnergyMC.push_back ( vector<double>() );
  	fMeanEnergyErrMC.push_back ( vector<double>() );
		
		fMeanXmaxConexMC.push_back ( vector<double>() );
		fMeanXmaxErrConexMC.push_back ( vector<double>() );
		fXmaxRMSConexMC.push_back ( vector<double>() );
		fXmaxRMSErrConexMC.push_back ( vector<double>() );
  }//end loop masses at Earth
    
	//init with zeros
	for(int s=0;s<fNbins_Xmax;s++) {
		fMeanXmaxData.push_back(0.);
		fMeanXmaxErrData.push_back(0.);
		fRMSData.push_back(0.);
		fRMSErrData.push_back(0.);
		fMeanEnergyData.push_back(0.);
		fMeanEnergyErrData.push_back(0.);

		fMeanXmaxFit.push_back(0.);
		fMeanXmaxErrFit.push_back(0.);
		fRMSFit.push_back(0.);
		fRMSErrFit.push_back(0.);
		fRMSTrueFit.push_back(0.);
		fRMSTrueErrFit.push_back(0.);

		fERRaw.push_back(0.);
		fERRawErr.push_back(0.);
		fERTrue.push_back(0.);
		fERTrueErr.push_back(0.);
  }//end loop Xmax energy bins

	for(int i=0;i<fNmassAtEarth;i++) {
		for(int s=0;s<fNbins_Xmax;s++) {
			fFractionsAtEarth_Xmaxfit[i].push_back(0.);
			fFractionsAtEarthLowErr_Xmaxfit[i].push_back(0.);
			fFractionsAtEarthHighErr_Xmaxfit[i].push_back(0.);

			fFractionsAtDetector_Xmaxfit[i].push_back(0.);
			fFractionsAtDetectorLowErr_Xmaxfit[i].push_back(0.);
			fFractionsAtDetectorHighErr_Xmaxfit[i].push_back(0.);

			fMeanXmaxMC[i].push_back(0.);
			fMeanXmaxErrMC[i].push_back(0.);
			fRMSMC[i].push_back(0.);
			fRMSErrMC[i].push_back(0.);
			fMeanEnergyMC[i].push_back(0.);
  		fMeanEnergyErrMC[i].push_back(0.);	
			fMeanXmaxConexMC[i].push_back(0.);
			fMeanXmaxErrConexMC[i].push_back(0.);		
	    fXmaxRMSConexMC[i].push_back(0.);
			fXmaxRMSErrConexMC[i].push_back(0.);		
		}//end loop Xmax energy bins
	}//end loop masses at Earth


	//#################################
	//## SET INITIAL PARAMETER VALUES
	//#################################
	fInjIndex_Xmax= fInjIndex;
	fEmax_Xmax= fEmax;
	fFluxNormFactor_Xmax= fFluxNormFactor;

	fInjIndex_spectrum= fInjIndex;
	fEmax_spectrum= fEmax;
	fFluxNormFactor_spectrum= fFluxNormFactor;

  fInjIndex_start= fInjIndex;
	fEmax_start= fEmax;
	fFluxNormFactor_start= fFluxNormFactor;

	fGalacticInjIndex_start= 2.0;
	fGalacticEmax_start= 18.00;
	fGalacticFluxNormFactor_start= 20.0;

	fSourceCutoffShape_start= fSourceCutoffShape;
	
	cout<<"*********************************"<<endl;
	cout<<"** INITIAL PASSED PARAMETERS ****"<<endl;
	cout<<"*********************************"<<endl;
	cout<<"fInjIndex_start="<<fInjIndex_start<<endl;
	cout<<"fEmax_start="<<fEmax_start<<endl;
	cout<<"fFluxNormFactor_start="<<fFluxNormFactor_start<<endl;
	

	//Allocate source fraction
	theMassFractStruct.fSourceFraction= new double[fNmassAtSource];
  theStartMassFractStruct.fSourceFraction= new double[fNmassAtSource];

	cout<<endl;
	
}//close SourceFitter::Init()


void SourceFitter::ResetData(){

	//Inizialize vector2D MC histo
	for(int i=0;i<fNmassAtEarth;i++) {
		fXmaxMC[i].clear();
		fXmaxMC[i].resize(0);
		
		fXmaxMC_scaled[i].clear();
		fXmaxMC_scaled[i].resize(0);

		fEnergyMC[i].clear();
		fEnergyMC[i].resize(0);

		fGenEnergyMC[i].clear();
		fGenEnergyMC[i].resize(0);
		
		fXmaxGenMC[i].clear();
		fXmaxGenMC[i].resize(0);
	}//close for i

}//close SourceFitter::ResetData()



void SourceFitter::MaxLikelihoodFcn_Xmax(int& nPar, double* const grad,double& value, double* const par,const int iFlag){

 
	value = 0.;

	//*******************************
  //****  SET FIT PARAMETER
  //*******************************
	int par_counter= 0;

	//## MASS FRACTIONS
	if(fFitExtragalFractions){
  	//***************************************************
  	//****  METHOD TO COSTRAIN THE FRACTIONS IN THE FIT
  	//***************************************************
  	// Unconstrained fractions are only constrained in [0,1] in Minuit set parameters
  	// Real fractions must sum to unity. Use method in "Statistical Methods in Data Analysis - W. J. Metzger"
  	// fract0= fract0_uncostr
  	// fract1= fract1_uncostr * (1-fract0_uncostr)
  	// fract2= fract2_uncostr * (1-fract0_uncostr) * (1-fract1_uncostr)
  	// ...
  	// fractN= (1-fract0_uncostr) * (1-fract1_uncostr) * ... * (1-fractN-1_uncostr)
  	double UnconstrainedFraction[fNmassAtSource-1];

  	for(int i=0; i< fNmassAtSource-1; i++) UnconstrainedFraction[i]= 0.;

  	UnconstrainedFraction[0]= *(par+0);
  	fFractionsAtSource[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
  	for(int i=1; i< fNmassAtSource-1; i++){
      UnconstrainedFraction[i]= *(par+i);
      fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
      ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}
  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;		

		par_counter+= fNmassAtSource-1;

	}//close if fFitExtragalFractions

	
	//## Injection index
	fInjIndex_Xmax= *(par+par_counter);
	par_counter++;
  
  //## Emax
	fEmax_Xmax= *(par+par_counter);
	par_counter++;
  
	//## Source cutoff shape
	if(fFitSourceCutoffShape){	
		fSourceCutoffShape = *(par+par_counter);
		par_counter++;
	}//close if add DATA shift in fit pars

	

	//###########################################
  //###  RESET PROP MC AND SET CURRENT VALUES
  //###########################################
  fPropMCReader->Reset();
 	
  fPropMCReader->SetInjectionIndex(fInjIndex_Xmax);
  fPropMCReader->SetMaxInjEnergyAtSource(fEmax_Xmax);
  fPropMCReader->SetRedShiftEvolution(fRedShiftEvolution);
  fPropMCReader->SetAbundanceAtSource(fFractionsAtSource);
  fPropMCReader->SetExtraGalSpectrumNormalization(fFluxNormFactor_Xmax);
	fPropMCReader->SetSourceCutoffMode(fSourceCutoffMode);
	fPropMCReader->SetSourceCutoffShape(fSourceCutoffShape);
	
  fPropMCReader->CalculateExpectedSpectrumAtEarth();
  fPropMCReader->CalculateMassFractions();
	
  fFractionsAtEarth_Xmaxfit= fPropMCReader->GetExpFractionsAtEarth_Xmax();
	

  //## Calculate current Xmax likelihood
  double theXmaxLikelihood= GetXmaxLikelihood();
  value= theXmaxLikelihood;

  if(fVerbosity>1){
  	cout<<"*****************************"<<endl;
  	cout<<"*** Xmax LL fit results ****"<<endl;
  	cout<<"*****************************"<<endl;
 	 	cout<<"current fractions @ Source "<<endl;
  	for(unsigned int i=0;i<fFractionsAtSource.size();i++){
			cout<<"FRACT @ Source"<<i+1<<": "<<fFractionsAtSource[i]<<endl;
  	}	  
  	cout<<"current InjIndex: "<<fInjIndex_Xmax<<endl;
  	cout<<"current Emax: "<<fEmax_Xmax<<endl;
		if(fFitSourceCutoffShape) cout<<"current SourceCutoffShape: "<<fSourceCutoffShape<<endl; 
		cout<<"current Xmax LL ==> "<<setprecision(5)<<value<<endl;
		cout<<"*****************************"<<endl;
	}
	fXmaxLikelihood= theXmaxLikelihood;
	
}//close SourceFitter::MaxLikelihoodFcn_Xmax()


double SourceFitter::GetXmaxLikelihood(){

	double TotalLogLikelihood= 0.;
	double PartialLogLikelihood= 0.;	
	
	fNdof_Xmax=0;
  fNevents_Xmax=0; 
	
	for(int s=0;s<fNbins_Xmax;s++){

		double Ndata = fXmaxData[s]->Integral();
		int    nbins = fXmaxData[s]->GetNbinsX();
		fNevents_Xmax+= Ndata;		

		PartialLogLikelihood = 0.;
  	
  	//Calculate Xmax LL 
  	for(int i=0;i<nbins;i++){
  		double f= 0.;
  		double n= fXmaxData[s]->GetBinContent(i+1);
  	
  		for(int j=0;j<fNmassAtEarth;j++){ 		
  			double Nmc= fXmaxMC[j][s]->Integral();
				double aji=0.;
				aji= fXmaxMC[j][s]->GetBinContent(i+1);
				f+= Ndata * fFractionsAtEarth_Xmaxfit[j][s]*aji/Nmc;		
			}//end loop masses at Earth
  	 
  		double LogLikelihood = 0.;
  	
  		if(f>0){
  			LogLikelihood= n*log(f)-f;
				fNdof_Xmax++;
			}
  	  	
    	PartialLogLikelihood+= -LogLikelihood; 
    	
		}//end loop Xmax bins

		TotalLogLikelihood+= PartialLogLikelihood;
		
	}//end loop Xmax energy bins

	return TotalLogLikelihood;

}//close SourceFitter::GetXmaxLikelihood()


//---------------------------------------------------------------------------
bool SourceFitter::FitXmax(int verbosity) {

  if (verbosity>0) cout << "\n Xmax fit " << endl;

	
  int ierflag;
  double arglist[2];  
  int nPar= 0;
	int nFittedFract= fNmassAtSource-1;

	nPar+= 2;//add InjIndex & Emax
	if(fFitExtragalFractions) nPar+= nFittedFract;//add fractions
	if(fFitSourceCutoffShape) nPar++;//add source cutoff shape

	cout << " XmaxFit: Number of parameters ==> " <<nPar<< endl;

  TMinuit theMinuit(nPar); 
  theMinuit.SetPrintLevel(0);//0 normal; -1 quiet; 1 verbose
  theMinuit.SetMaxIterations(10000);
  theMinuit.SetFCN(SourceFitter::MaxLikelihoodFcn_Xmax);
  theMinuit.mnexcm("SET NOW",arglist,0,ierflag);
  
  const double stepSize    = 0.1;
  const double minFraction = 0.;
  const double maxFraction = 1.;
  const double startValue  = 1./(nFittedFract+1);//initialize to equal values for all fractions  
	

  const double startEmax   = fEmax_start;
  const double stepEmax    = 0.1;
  const double minEmax     = 18.0;
  const double maxEmax     = 21.0;
    
  const double startInjIndex  = fInjIndex_start;
  const double stepInjIndex   = 0.1;
  const double minInjIndex	  = 0.0;
  const double maxInjIndex	  = 4.0;

	const double startSourceCutoffShape = fSourceCutoffShape_start;
  const double stepSourceCutoffShape  = 0.1;
  const double minSourceCutoffShape	  = 0.;
  const double maxSourceCutoffShape  = 0.;


	//## Save initial starting values
	theStartMassFractStruct.fSourceFraction[0]= fStartFractParamsAtSource[0];
  double ConstrFactor= (1.-fStartFractParamsAtSource[0]);
  for(int i=1; i< fNmassAtSource-1; i++){
      theStartMassFractStruct.fSourceFraction[i]= fStartFractParamsAtSource[i]*ConstrFactor;
      ConstrFactor*= (1.-fStartFractParamsAtSource[i]);
  }
  theStartMassFractStruct.fSourceFraction[fNmassAtSource-1]= ConstrFactor;	

	
	//## Convert start fraction values to start fraction params
	double startFractParamValue[fNmassAtSource-1];
	startFractParamValue[0]= startValue;
	ConstrFactor= 1.-startFractParamValue[0];
	for(int k=1;k<fNmassAtSource-1;k++){		
		startFractParamValue[k]= startValue/ConstrFactor;
		ConstrFactor*= (1.-startFractParamValue[k]);
	}

	
	//########################################
	//# SET FIT PARAMS
	//########################################
	int par_counter=0;

	//## ADD EXTRAGALACTIC FRACTIONS
	if(fFitExtragalFractions){
  	int thisPar=0;
  	for(int iPar=0;iPar<nFittedFract;iPar++){

  		ostringstream parName;
    	parName << "f"<<"_" << iPar;
     
			theMinuit.mnparm(thisPar,parName.str().c_str(),startFractParamValue[iPar],stepSize,minFraction,maxFraction,ierflag);
			fStartFractFitParam[iPar]= startFractParamValue[iPar]; 	
		
			//# fix the fractions?
			if(fFixMassFractions){
				theMinuit.FixParameter(thisPar);// fix mass fractions
			}

    	thisPar++;
			par_counter++;
  	}//end loop extragal fractions		 
	}//close if fFitExtragalFractions

	//## ADD INJ INDEX
	theMinuit.mnparm(par_counter,"InjIndex",startInjIndex,stepInjIndex,minInjIndex,maxInjIndex,ierflag);
	if(fFixInjIndex) theMinuit.FixParameter(par_counter);//fix gamma
	par_counter++;

	//## ADD EMAX
  theMinuit.mnparm(par_counter,"Emax",startEmax,stepEmax,minEmax,maxEmax,ierflag); 
	if(fFixEmax) theMinuit.FixParameter(par_counter);//fix Emax 
	par_counter++;

	//## ADD SOURCE CUTOFF SHAPE
	if(fFitSourceCutoffShape) {
		theMinuit.mnparm(par_counter,"SourceCutoffShape",startSourceCutoffShape,stepSourceCutoffShape,minSourceCutoffShape,maxSourceCutoffShape,ierflag);
		if(fFixSourceCutoffShape) theMinuit.FixParameter(par_counter);//fix source cutoff shape
		par_counter++;  
	}

	
  
	//## SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
 	arglist[0]= 0.5;//0.5 likelihood, 1 ChiSquare
  theMinuit.mnexcm("SET ERR",arglist,1,ierflag);

	//## SET VERBOSITY
	//arglist[0]= 3;
  //theMinuit.mnexcm("SET PRI",arglist,1,ierflag);

	
	//## SET MINUIT COMMAND
	arglist[0]= 10000;
  arglist[1]= 0.1;
  theMinuit.mnexcm("MINIMIZE", arglist ,2,ierflag);//use MINIMIZE for scan
	//theMinuit.mnexcm("MIGRAD", arglist ,2,ierflag);//use MIGRAD for scan
	

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  theMinuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

	//## Get covariance matrix
	double covMatrix[nPar][nPar];		
  theMinuit.mnemat(&covMatrix[0][0],nPar);
  	
  cout<<"** COV MATRIX ***"<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
		
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i

	//Store cov matrix
	fCovMatrix.clear();
	fCovMatrix.resize(0);
	fCovMatrix.push_back(CovarianceMatrix);  


  //#################################
	//####  GET FIT PARAMS    #########
	//#################################
	double p,ep;
	par_counter=0;
	fNmassFree= fNmassAtSource-1;

	fitParameterList.clear();
	fitParameterList.resize(0);
	fitParameterErrorList.clear();
	fitParameterErrorList.resize(0);

	for(int i=0;i<nPar;i++){
		theMinuit.GetParameter(i,p,ep);	
		fitParameterList.push_back(p);
		fitParameterErrorList.push_back(ep);
	}


	if(fFitExtragalFractions){
  	theMinuit.GetParameter(0,p,ep);
		fFractFitParam[0]= p;
		fFractFitParamErr[0]= ep;

		double UnconstrainedFraction[fNmassAtSource-1];
  	UnconstrainedFraction[0]= p;
		
  	fFractionsAtSource[0]= UnconstrainedFraction[0];
		fSourceFract[0]= UnconstrainedFraction[0];

  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
	
  	for(int i=1;i<nFittedFract;i++){
			theMinuit.GetParameter(i,p,ep);
			fFractFitParam[i]= p;
			fFractFitParamErr[i]= ep;

			UnconstrainedFraction[i]= p;
    	fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
			fSourceFract[i]= UnconstrainedFraction[i]*ConstrFactor;

    	ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}//end loop fitted fractions

  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;
		fSourceFract[fNmassAtSource-1]= ConstrFactor;

		for(int i=0;i<fNmassAtSource;i++) theMassFractStruct.fSourceFraction[i]= fFractionsAtSource[i];

		//Calculate mass fraction uncertainties
		fitParameterErrorList= GetParamUncertainty(fitParameterList,CovarianceMatrix,nPar+1);
		for(unsigned int k=0;k<fitParameterErrorList.size();k++) cout<<"ParErr "<<k<<" ==> "<<fitParameterErrorList[k]<<endl;

		
		//Calculate low and upper bound
		double ErrLowBound=0.;
  	double ErrHighBound=0.;
		for(int i=0;i<fNmassAtSource;i++){
	  	ErrHighBound= fmin(1.,fFractionsAtSource[i]+fitParameterErrorList[i])-fFractionsAtSource[i];
 	  	ErrLowBound= fFractionsAtSource[i]-fmax(0.,fFractionsAtSource[i]-fitParameterErrorList[i]);
 	  	fFractionsAtSourceHighErr[i]= ErrHighBound; 
 	  	fFractionsAtSourceLowErr[i]= ErrLowBound;
			fSourceFractHighErr[i]= ErrHighBound;
			fSourceFractLowErr[i]= ErrLowBound;
			cout<<"Mass Fractions ==> "<<i+1<<"  "<<fFractionsAtSource[i]<<" +- "<< fFractionsAtSourceHighErr[i]<<"  "<<fFractionsAtSourceLowErr[i]<<endl;
  	}//end loop fractions at source
			
		par_counter+= nFittedFract;
	}//close if fFitExtragalFractions

	//get InjIndex
  theMinuit.GetParameter(par_counter,p,ep);
  fInjIndex_Xmax= p;
  fInjIndexErr_Xmax= ep;
	par_counter++;
  
	//get Emax
  theMinuit.GetParameter(par_counter,p,ep);
  fEmax_Xmax= p;
  fEmaxErr_Xmax= ep; 
	par_counter++; 

	if(fFitSourceCutoffShape) {
		theMinuit.GetParameter(par_counter,p,ep);
  	fSourceCutoffShape= p; 
		fSourceCutoffShapeErr= ep;
		par_counter++;
	}

	fFitStatus= theMinuit.GetStatus();
  

  cout<<"***************************************"<<endl;
  cout<<"*** Xmax LL fit final results     ****"<<endl;
  cout<<"***************************************"<<endl; 
  cout<<"Fractions @ Source "<<endl;
  for(unsigned int i=0;i<fFractionsAtSource.size();i++) cout<<"FRACT @ Source"<<i+1<<" ==> "<<fFractionsAtSource[i]<<endl;
  cout<<"InjIndex Xmax= "<<fInjIndex_Xmax<<endl;
  cout<<"Emax Xmax= "<<fEmax_Xmax<<endl;
  cout<<"Xmax LL= "<<fXmaxLikelihood<<endl;
  if(fFitSourceCutoffShape) cout<<"SourceCutoffShape="<<fSourceCutoffShape<<endl;
	
  fXmaxMinimizationInfo->Fill();


	
	//###########################
	//#####   STORE HISTO    ####
	//###########################
	StoreXmaxFitResults();
	StoreElongationRate();
	StoreRMS();
	StoreMassFraction();
	
  return true;

}//close SourceFitter::FitXmax


void SourceFitter::MaxLikelihoodFcn_EnergySpectrum(int& nPar, double* const grad,double& value, double* const par,const int iFlag){

	value= 0.;

	//############################
	//##  SET FIT PARAMS
	//###########################
	int par_counter=0;

	if(fFitExtragalFractions){
  	//***************************************************
  	//****  METHOD TO COSTRAIN THE FRACTIONS IN THE FIT
  	//***************************************************
  	// Unconstrained fractions are only constrained in [0,1] in Minuit set parameters
  	// Real fractions must sum to unity. Use method in "Statistical Methods in Data Analysis - W. J. Metzger"
  	// fract0= fract0_uncostr
  	// fract1= fract1_uncostr * (1-fract0_uncostr)
  	// fract2= fract2_uncostr * (1-fract0_uncostr) * (1-fract1_uncostr)
  	// ...
  	// fractN= (1-fract0_uncostr) * (1-fract1_uncostr) * ... * (1-fractN-1_uncostr)
  	double UnconstrainedFraction[fNmassAtSource-1];
  	
		for(int i=0; i< fNmassAtSource-1; i++) UnconstrainedFraction[i]=0.;

  	UnconstrainedFraction[0]= *(par+0);
  	fFractionsAtSource[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
  	for(int i=1; i< fNmassAtSource-1; i++){
      UnconstrainedFraction[i]= *(par+i);
      fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
      ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}
  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;	

		par_counter+= fNmassAtSource-1;
	}//close if fFitExtragalFractions


	//set gamma injection index
  fInjIndex_spectrum= *(par+par_counter);
	par_counter++;
  
  //set Emax
  fEmax_spectrum= *(par+par_counter);
	par_counter++;

	//set normalization
  fFluxNormFactor_spectrum= pow(10,*(par+par_counter));
	par_counter++;

	//set source cutoff shape
	if(fFitSourceCutoffShape){
		fSourceCutoffShape= *(par+par_counter);
		par_counter++;
	}

	if(fFitGalSpectrum){
		//set gamma galactic injection index
		fGalacticInjIndex_spectrum= *(par+par_counter);
		par_counter++;
  
		//set Emax
		fGalacticEmax_spectrum= *(par+par_counter);
		par_counter++;

		//set normalization
		fGalacticFluxNormFactor_spectrum= pow(10,*(par+par_counter));
		par_counter++;
	}//close if fFitGalSpectrum
	
	
	//###########################################
 	//###  RESET PROP MC AND READ WITH CURRENT VALUES
 	//###########################################
	fPropMCReader->Reset();
 	
	//## EXTRAGALACTIC
	fPropMCReader->SetInjectionIndex(fInjIndex_spectrum);
	fPropMCReader->SetMaxInjEnergyAtSource(fEmax_spectrum);
	fPropMCReader->SetRedShiftEvolution(fRedShiftEvolution);
	fPropMCReader->SetAbundanceAtSource(fFractionsAtSource);
	fPropMCReader->SetExtraGalSpectrumNormalization(fFluxNormFactor_spectrum);
	fPropMCReader->SetSourceCutoffMode(fSourceCutoffMode);
	fPropMCReader->SetSourceCutoffShape(fSourceCutoffShape);
		
	//## GALACTIC
	fPropMCReader->SetGalacticInjectionIndex(fGalacticInjIndex_spectrum);
	fPropMCReader->SetGalacticMaxInjEnergyAtSource(fGalacticEmax_spectrum);
	fPropMCReader->SetGalSpectrumNormalization(fGalacticFluxNormFactor_spectrum);
	if(fFitGalFractions) fPropMCReader->IncludeGalacticComposition(true);//add also galactic mass to the fractions

	fPropMCReader->CalculateExpectedSpectrumAtEarth();
	fPropMCReader->GenerateGalacticSpectrum();
	fPropMCReader->CalculateMassFractions();
	
	fPropEnergySpectrum_spectrumfit= fPropMCReader->GetExpEnergySpectrumAtEarth_Spectrum();
	fPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpEnergySpectrumE3AtEarth_Spectrum();
  fPropDirectEnergySpectrum_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumAtEarth_Spectrum();
	fPropDirectEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumE3AtEarth_Spectrum();
  
	fTotPropEnergySpectrum_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumAtEarth_Spectrum();
	fTotPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumE3AtEarth_Spectrum();
	fGalacticSpectrum_spectrumfit= fPropMCReader->GetGalacticSpectrum_Spectrum();

	fFractionsAtEarth_spectrumfit= fPropMCReader->GetExpFractionsAtEarth_Spectrum();
	fFractionsAtEarth_Xmaxfit= fPropMCReader->GetExpFractionsAtEarth_Xmax();	
	fDirectFractionsAtEarth_spectrumfit= fPropMCReader->GetExpDirectFractionsAtEarth_Spectrum();


  //Get LL spectrum
  double theSpectrumLikelihood= GetEnergySpectrumLikelihood();
  value= theSpectrumLikelihood;
  
	if(fVerbosity>1){
		cout<<"*********************************"<<endl;
  	cout<<"*** Spectrum LL fit results ****"<<endl;
  	cout<<"*********************************"<<endl;
  	cout<<"current fractions @ Source "<<endl;
  	for(unsigned int i=0;i<fFractionsAtSource.size();i++){
			cout<<"FRACT @ Source "<<i+1<<" ==> "<<fFractionsAtSource[i]<<endl;
  	}	  
  	cout<<"GammaEG= "<<fInjIndex_spectrum<<endl;
  	cout<<"EmaxEG= "<<fEmax_spectrum<<endl;
  	cout<<"NormEG= "<<fFluxNormFactor_spectrum<<endl;	
		if(fFitGalSpectrum){
  		cout<<"GammaGAL= "<<fGalacticInjIndex_spectrum<<endl;
  		cout<<"EmaxGAL= "<<fGalacticEmax_spectrum<<endl;
  		cout<<"NormGAL= "<<fGalacticFluxNormFactor_spectrum<<endl;
		}
		if(fFitSourceCutoffShape) cout<<"SourceCutoffShape= "<<fSourceCutoffShape<<endl;	
		cout<<"current LL ==> "<<value<<endl;
		cout<<"*********************************"<<endl;
	}
	fSpectrumLikelihood= theSpectrumLikelihood;
	
}//close SourceFitter::MaxLikelihoodFcn_EnergySpectrum()


double SourceFitter::GetEnergySpectrumLikelihood(){

	double TotalLogLikelihood=0.;
	double PartialLogLikelihood=0.;	
	
	fNdof_spectrum=0;
	fNevents_spectrum=0;	

	const double Ndata  = fAugerEnergySpectrumEvents->Integral();
	const int    nbins = fAugerEnergySpectrumEvents->GetNbinsX();
  	
	fNevents_spectrum= Ndata;

  //Calculate LogLike
  for(int i=0;i<nbins;i++){
  	double f=0.;
  		
		const double n = fAugerEnergySpectrumEvents->GetBinContent(i+1);		
		double aji_EG= fTotPropEnergySpectrum_spectrumfit->GetBinContent(i+1);
		double aji_GAL= fGalacticSpectrum_spectrumfit->GetBinContent(i+1);
		
  	if(fFitGalSpectrum) f= fExposure[i]*(aji_EG + aji_GAL);
		else f= fExposure[i]*aji_EG;
    
  	double LogLikelihood = 0.;
  	
  	if(f>0){
  		LogLikelihood= n*log(f)-f;
			fNdof_spectrum++;
		}
  	 		
    PartialLogLikelihood+= -LogLikelihood; 
   	
	}//end loop spectrum bins

	TotalLogLikelihood= PartialLogLikelihood;
	
	return TotalLogLikelihood;

}//close SourceFitter::GetEnergySpectrumLikelihood()


bool SourceFitter::FitEnergySpectrum(int verbosity) {

  if (verbosity>0) cout << "\n EnergySpectrum fit " << endl;

	
  int ierflag;
  double arglist[2]; 
	int nPar;
	int nFittedFract= fNmassAtSource-1;

	nPar= 3; //add EG InjIndex + EG Emax + EG Norm
	if(fFitExtragalFractions) nPar+= nFittedFract; //add fractions
	if(fFitSourceCutoffShape) nPar++;//add source cutoff shape
	if(fFitGalSpectrum) nPar+= 3; //add Gal InjIndex + Gal Emax + Gal Norm
	
  cout << " SpectrumFit: Number of parameters ==> " <<nPar<< endl;
 
  TMinuit theMinuit(nPar); 
  theMinuit.SetPrintLevel(0);//0 normal; -1 quiet; 1 verbose
  theMinuit.SetFCN(SourceFitter::MaxLikelihoodFcn_EnergySpectrum);
  theMinuit.mnexcm("SET NOW",arglist,0,ierflag);
  theMinuit.SetMaxIterations(10000);

	const double stepSize    = 0.1;
  const double minFraction = 0.;
  const double maxFraction = 1.;
  const double startValue  = 1./(nFittedFract+1);

  const double startEmax	 = fEmax_start;
  const double stepEmax    = 0.1;
  const double minEmax     = 18.0;
  const double maxEmax     = 21.0;
    
	const double startInjIndex  = fInjIndex_start;
	const double stepInjIndex   = 0.1;
	const double minInjIndex	  = 0.0;
	const double maxInjIndex	  = 4.0;

	double startSourceCutoffShape = fSourceCutoffShape_start;
  const double stepSourceCutoffShape  = 0.1;
  const double minSourceCutoffShape	  = 0.;
  const double maxSourceCutoffShape  = 0.;

  double NormStartValue  = 10;
  const double NormStepSize    = 1; 
  const double minNorm         = 0.001;
  const double maxNorm         = 100.;
 
  //change norm according to index
  if(fInjIndex_spectrum<=1.5) {
		NormStartValue = -8;	
  }
  if(fInjIndex_spectrum>1.5&&fInjIndex_spectrum<2.) {
	  NormStartValue = 1;	
  }
  if(fInjIndex_spectrum>=2.&&fInjIndex_spectrum<2.5) {
	  NormStartValue = 10;	
  }
  if(fInjIndex_spectrum>=2.5&&fInjIndex_spectrum<3.) {
		NormStartValue = 20;	
  }
  if(fInjIndex_spectrum>=3&&fInjIndex_spectrum<3.5) {
		NormStartValue = 29;	
  }
  if(fInjIndex_spectrum>=3.5&&fInjIndex_spectrum<4.) {
  	NormStartValue = 37;	
  }

	fFluxNormFactor_start= NormStartValue;	

  
  //### GALACTIC PARAMETERS ###
  fGalacticEmax_spectrum= 18.2;
  const double startGalEmax	 = fGalacticEmax_start;
  const double stepGalEmax    = 0.1;
  const double minGalEmax     = 18.0;
  const double maxGalEmax     = 18.5;
    
  fGalacticInjIndex_spectrum= 3.0;
  const double startGalInjIndex  = fGalacticInjIndex_start;
  const double stepGalInjIndex   = 0.1;
  const double minGalInjIndex	 = 1.;
  const double maxGalInjIndex	 = 5.0;

  fGalacticFluxNormFactor_spectrum= 20;
  const double GalNormStartValue  = fGalacticFluxNormFactor_start;
  const double GalNormStepSize    = 1; 
  const double minGalNorm         = 0.001;
  const double maxGalNorm         = 100;


	//########################################
	//# Store initial values of the fractions
	//########################################
	theStartMassFractStruct.fSourceFraction[0]= fStartFractParamsAtSource[0];
  double ConstrFactor= (1.-fStartFractParamsAtSource[0]);
  for(int i=1; i< fNmassAtSource-1; i++){
  	theStartMassFractStruct.fSourceFraction[i]= fStartFractParamsAtSource[i]*ConstrFactor;
    ConstrFactor*= (1.-fStartFractParamsAtSource[i]);
  }
  theStartMassFractStruct.fSourceFraction[fNmassAtSource-1]= ConstrFactor;	

	//convert into start values for fraction params
	double startFractParamValue[fNmassAtSource-1];

	startFractParamValue[0]= startValue;
	ConstrFactor= 1.-startFractParamValue[0];
	for(int k=1;k<fNmassAtSource-1;k++){				
		startFractParamValue[k]= startValue/ConstrFactor;
		ConstrFactor*= (1.-startFractParamValue[k]);
	}
	
	//########################################
	//# SET FIT PARAMS
	//########################################
	int par_counter=0;

	//## ADD EXTRAGALACTIC PARAMETERS	
	if(fFitExtragalFractions){
  	int thisPar=0;
  	for(int iPar=0;iPar<nFittedFract;iPar++){
  		ostringstream parName;
    	parName << "f"<<"_" << iPar;
     
    	theMinuit.mnparm(thisPar,parName.str().c_str(),startFractParamValue[iPar],stepSize,minFraction,maxFraction,ierflag);
			fStartFractFitParam[iPar]= startFractParamValue[iPar];
			
			//# fix the fractions?
			if(fFixMassFractions){
				theMinuit.FixParameter(thisPar);// fix mass fractions
			}

    	thisPar++;
			par_counter++;
  	}//end loop extragal fractions
		
	}//close if fFitExtragalFractions
	
	//add inj index
  theMinuit.mnparm(par_counter,"InjIndex",startInjIndex,stepInjIndex,minInjIndex,maxInjIndex,ierflag);
	if(fFixInjIndex) theMinuit.FixParameter(par_counter);// fix InjIndex
	par_counter++;

  //add Emax
  theMinuit.mnparm(par_counter,"Emax",startEmax,stepEmax,minEmax,maxEmax,ierflag); 
	if(fFixEmax) theMinuit.FixParameter(par_counter);// fix Emax
	par_counter++;
 
  //add normalization to fit parameters
  theMinuit.mnparm(par_counter,"Flux Norm",NormStartValue,NormStepSize,-60,60,ierflag);
	par_counter++;  

	//add cutoff shape
	if(fFitSourceCutoffShape) {
		theMinuit.mnparm(par_counter,"SourceCutoffShape",startSourceCutoffShape,stepSourceCutoffShape,minSourceCutoffShape,maxSourceCutoffShape,ierflag);
		if(fFixSourceCutoffShape) theMinuit.FixParameter(par_counter);//fix source cutoff shape
		par_counter++;  
	}



	//## ADD GALACTIC PARAMETERS	
	if(fFitGalSpectrum){
		//add inj index
  	theMinuit.mnparm(par_counter,"InjIndexGalactic",startGalInjIndex,stepGalInjIndex,minGalInjIndex,maxGalInjIndex,ierflag);
		par_counter++;

  	//add Emax @ source
  	theMinuit.mnparm(par_counter,"EmaxGalactic",startGalEmax,stepGalEmax,minGalEmax,maxGalEmax,ierflag);
		par_counter++;
  
  	//add normalization to fit parameters
  	theMinuit.mnparm(par_counter,"Gal Flux Norm",GalNormStartValue,GalNormStepSize,-30,60,ierflag);
		par_counter++;
	}//close if



	//## SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
  arglist[0]= 1;
  theMinuit.mnexcm("SET STR",arglist,1,ierflag);

	//## SET MINUIT ERROR STRATEGY
  arglist[0]= 0.5;//0.5 likelihood, 1 ChiSquare
  theMinuit.mnexcm("SET ERR",arglist,1,ierflag);


	//## SET MINUIT VERBOSITY
  //arglist[0]= 3;
  //theMinuit.mnexcm("SET PRI",arglist,1,ierflag);


	//## SET MINUIT COMMAND
	arglist[0]= 10000;
  arglist[1]= 0.1;
  theMinuit.mnexcm("MINIMIZE", arglist ,2,ierflag);//use MIGRAD for scan
	//theMinuit.mnexcm("MIGRAD", arglist ,2,ierflag);
  //theMinuit.mnexcm("HESSE", arglist ,0,ierflag);
  //improve local minimum ==> try to exit from bad local minima (hopefully!)
  //theMinuit.mnexcm("IMPROVE", arglist ,0,ierflag);

  double LikelihoodMin, edm, errdef;
  int nvpar, nparx, icstat;
  theMinuit.mnstat(LikelihoodMin, edm, errdef, nvpar, nparx, icstat);
  
	//Get covariance matrix
	double covMatrix[nPar][nPar];		
  theMinuit.mnemat(&covMatrix[0][0],nPar);
  	
  cout<<"*** COV MATRIX ***"<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
		
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i

	//Store cov matrix
	fCovMatrix.clear();
	fCovMatrix.resize(0);
	fCovMatrix.push_back(CovarianceMatrix);  


	//##########################
	//### GET FIT PARAMS     ###
	//##########################
  double p,ep;
	par_counter= 0;
	fNmassFree= fNmassAtSource-1;

	fitParameterList.clear();
	fitParameterList.resize(0);
	fitParameterErrorList.clear();
	fitParameterErrorList.resize(0);

	for(int i=0;i<nPar;i++){
		theMinuit.GetParameter(i,p,ep);	
		fitParameterList.push_back(p);
		fitParameterErrorList.push_back(ep);
	}

  
	if(fFitExtragalFractions){
  	double UnconstrainedFraction[fNmassAtSource-1];

  	theMinuit.GetParameter(0,p,ep);
  	UnconstrainedFraction[0]= p;
		fFractFitParam[0]= p;
		fFractFitParamErr[0]= ep;

  	fFractionsAtSource[0]= UnconstrainedFraction[0];
		fSourceFract[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
	
  	for(int i=1;i<nFittedFract;i++){
			theMinuit.GetParameter(i,p,ep);
			fFractFitParam[i]= p;
			fFractFitParamErr[i]= ep;

			UnconstrainedFraction[i]= p;
    	fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
			fSourceFract[i]= UnconstrainedFraction[i]*ConstrFactor;
    	ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}//close for i
  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;
		fSourceFract[fNmassAtSource-1]= ConstrFactor;

		//Calculate mass fraction uncertainties
		fitParameterErrorList= GetParamUncertainty(fitParameterList,CovarianceMatrix,nPar+1);
		for(unsigned int k=0;k<fitParameterErrorList.size();k++) cout<<"ParErr "<<k<<" ==> "<<fitParameterErrorList[k]<<endl;

		
		//Calculate low and upper bound
		double ErrLowBound=0.;
  	double ErrHighBound=0.;
		for(int i=0;i<fNmassAtSource;i++){
	  	ErrHighBound= fmin(1.,fFractionsAtSource[i]+fitParameterErrorList[i])-fFractionsAtSource[i];
 	  	ErrLowBound= fFractionsAtSource[i]-fmax(0.,fFractionsAtSource[i]-fitParameterErrorList[i]);
 	  	fFractionsAtSourceHighErr[i]= ErrHighBound; 
 	  	fFractionsAtSourceLowErr[i]= ErrLowBound;
			fSourceFractHighErr[i]= ErrHighBound;
			fSourceFractLowErr[i]= ErrLowBound;
			cout<<"Mass Fractions ==> "<<i+1<<"  "<<fFractionsAtSource[i]<<" +- "<< fFractionsAtSourceHighErr[i]<<"  "<<fFractionsAtSourceLowErr[i]<<endl;
  	}//end loop masses at source 

		par_counter+= nFittedFract;
	}//close if fFitExtragalFractions 



	//get InjIndex
  theMinuit.GetParameter(par_counter,p,ep);
  fInjIndex_spectrum= p;
  fInjIndexErr_spectrum= ep;
	par_counter++;
 
	//get Emax
  theMinuit.GetParameter(par_counter,p,ep);
  fEmax_spectrum= p;
  fEmaxErr_spectrum= ep;  
	par_counter++;

	//get norm factor
	theMinuit.GetParameter(par_counter,p,ep);
  fFluxNormFactor_spectrum= p;
  fFluxNormFactorErr_spectrum= ep;
	par_counter++;

	//get source cutoff shape
	if(fFitSourceCutoffShape) {
		theMinuit.GetParameter(par_counter,p,ep);
  	fSourceCutoffShape= p; 
		fSourceCutoffShapeErr= ep;
		par_counter++;
	} 


	//get gal parameter
	if(fFitGalSpectrum){
		theMinuit.GetParameter(par_counter,p,ep);
  	fGalacticInjIndex_spectrum= p;
		fGalacticInjIndexErr_spectrum= ep;
		par_counter++;
  
  	theMinuit.GetParameter(par_counter,p,ep);
  	fGalacticEmax_spectrum= p;
		fGalacticEmaxErr_spectrum= ep;
		par_counter++;
    
  	theMinuit.GetParameter(par_counter,p,ep);
  	fGalacticFluxNormFactor_spectrum= p;
		fGalacticFluxNormFactorErr_spectrum= ep;
		par_counter++;
	}//close if

	
  fFitStatus= theMinuit.GetStatus();
 
  
	cout<<"***************************************"<<endl;
  cout<<"*** Spectrum LL fit final results ****"<<endl;
  cout<<"***************************************"<<endl; 
	cout<<"Fractions @ Source "<<endl;
  for(unsigned int i=0;i<fFractionsAtSource.size();i++) cout<<"FRACT @ Source"<<i+1<<" ==> "<<fFractionsAtSource[i]<<endl; 
  cout<<"InjIndex= "<<fInjIndex_spectrum<<endl;
  cout<<"Emax= "<<fEmax_spectrum<<endl;
  cout<<"FluxNorm= "<<fFluxNormFactor_spectrum<<endl;
	if(fFitSourceCutoffShape) cout<<"SourceCutoffShape= "<<fSourceCutoffShape<<endl;	
	if(fFitGalSpectrum){
		cout<<"Gal InjIndex= "<<fGalacticInjIndex_spectrum<<endl;
  	cout<<"Gal Emax= "<<fGalacticEmax_spectrum<<endl;
  	cout<<"Gal FluxNorm= "<<fGalacticFluxNormFactor_spectrum<<endl;
	}
	for(int i=0;i<fNmassAtSource;i++) theMassFractStruct.fSourceFraction[i]= fFractionsAtSource[i];
	
  fSpectrumMinimizationInfo->Fill();
	

	//###########################
	//#####   STORE HISTO    ####
	//###########################
	StoreSpectrumFitResults();

  return true;

}//close SourceFitter::FitEnergySpectrum()



void SourceFitter::MaxLikelihoodFcnCombined(int& nPar,double* const grad, double& value,double* const par,const int iFlag) {
  
  value  = 0;
  
  //############################
	//##  SET FIT PARAMS
	//###########################
	int par_counter=0;
  
	if(fFitExtragalFractions){
	
  	//***************************************************
  	//****  METHOD TO COSTRAIN THE FRACTIONS IN THE FIT
  	//***************************************************
  	// Unconstrained fractions are only constrained in [0,1] in Minuit set parameters
  	// Real fractions must sum to unity. Use method in "Statistical Methods in Data Analysis - W. J. Metzger"
 	 	// fract0= fract0_uncostr
  	// fract1= fract1_uncostr * (1-fract0_uncostr)
  	// fract2= fract2_uncostr * (1-fract0_uncostr) * (1-fract1_uncostr)
  	// ...
  	// fractN= (1-fract0_uncostr) * (1-fract1_uncostr) * ... * (1-fractN-1_uncostr)
  	double UnconstrainedFraction[fNmassAtSource-1];

  	for(int i=0; i< fNmassAtSource-1; i++) UnconstrainedFraction[i]=0.;

  	UnconstrainedFraction[0]= *(par+0);
  	fFractionsAtSource[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
  	for(int i=1; i< fNmassAtSource-1; i++){
      UnconstrainedFraction[i]= *(par+i);
      fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
      ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}
  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;		

		par_counter+= fNmassAtSource-1;
	}//close if fFitExtragalFractions


  //set gamma injection index
  fInjIndex= *(par+par_counter);
	fInjIndex_spectrum= fInjIndex;
	par_counter++;
  
  //set Emax
  fEmax= *(par+par_counter);
	fEmax_spectrum= fEmax;
  par_counter++;

  //Set flux normalization
  fFluxNormFactor= pow(10,*(par+par_counter));
	fFluxNormFactor_spectrum= fFluxNormFactor;
	par_counter++;

	//Set source cutoff shape
	if(fFitSourceCutoffShape){
		fSourceCutoffShape = *(par+par_counter);
		par_counter++;
	}

  	
	if(fFitGalSpectrum){
 	 	//set galactic injection index
  	fGalacticInjIndex_spectrum= *(par+par_counter);
		par_counter++;
  
  	//set galactic Emax
  	fGalacticEmax_spectrum= *(par+par_counter);
		par_counter++;
  
  	//Set flux normalization
  	fGalacticFluxNormFactor_spectrum= pow(10,*(par+par_counter));
		par_counter++;

  }//close if fFitGalSpectrum


  //###########################################
 	//###  RESET PROP MC AND READ WITH CURRENT VALUES
 	//###########################################
	fPropMCReader->Reset();
 	
	//## EXTRAGALACTIC
	fPropMCReader->SetInjectionIndex(fInjIndex);
	fPropMCReader->SetMaxInjEnergyAtSource(fEmax);
	fPropMCReader->SetRedShiftEvolution(fRedShiftEvolution);
	fPropMCReader->SetAbundanceAtSource(fFractionsAtSource);//how many masses? more than the fitted at Earth (p,He,CNO,Fe)?
	fPropMCReader->SetExtraGalSpectrumNormalization(fFluxNormFactor);
	fPropMCReader->SetSourceCutoffMode(fSourceCutoffMode);
	fPropMCReader->SetSourceCutoffShape(fSourceCutoffShape);
	
	//## GALACTIC
	fPropMCReader->SetGalacticInjectionIndex(fGalacticInjIndex_spectrum);
	fPropMCReader->SetGalacticMaxInjEnergyAtSource(fGalacticEmax_spectrum);
	fPropMCReader->SetGalSpectrumNormalization(fGalacticFluxNormFactor_spectrum);
	if(fFitGalFractions) fPropMCReader->IncludeGalacticComposition(true);//add also galactic mass to the fractions

	fPropMCReader->CalculateExpectedSpectrumAtEarth();
	fPropMCReader->GenerateGalacticSpectrum();
	fPropMCReader->CalculateMassFractions();
	
	fPropEnergySpectrum_spectrumfit= fPropMCReader->GetExpEnergySpectrumAtEarth_Spectrum();
	fPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpEnergySpectrumE3AtEarth_Spectrum();
  fPropDirectEnergySpectrum_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumAtEarth_Spectrum();
	fPropDirectEnergySpectrumE3_spectrumfit= fPropMCReader->GetExpDirectEnergySpectrumE3AtEarth_Spectrum();
  
	fTotPropEnergySpectrum_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumAtEarth_Spectrum();
	fTotPropEnergySpectrumE3_spectrumfit= fPropMCReader->GetTotExpEnergySpectrumE3AtEarth_Spectrum();
	fGalacticSpectrum_spectrumfit= fPropMCReader->GetGalacticSpectrum_Spectrum();

	fFractionsAtEarth_spectrumfit= fPropMCReader->GetExpFractionsAtEarth_Spectrum();
	fFractionsAtEarth_Xmaxfit= fPropMCReader->GetExpFractionsAtEarth_Xmax();	
	fDirectFractionsAtEarth_spectrumfit= fPropMCReader->GetExpDirectFractionsAtEarth_Spectrum();


  //## Get spectrum LL
  double theSpectrumLikelihood= GetEnergySpectrumLikelihood();
  
  //## Get Xmax LL
  double theXmaxLikelihood= GetXmaxLikelihood();
  
	cout<<"NdofSpectrum="<<fNdof_spectrum<<"  NdofXmax="<<fNdof_Xmax<<endl;
	cout<<"NeventsSpectrum="<<fNevents_spectrum<<"  NeventsXmax="<<fNevents_Xmax<<endl;
  
	//## use weights?
	if(fWeightLikelihood) {
		theSpectrumLikelihood/= fSpectrumLikelihoodWeight; 
		theXmaxLikelihood/= fXmaxLikelihoodWeight; 
		//theSpectrumLikelihood/= fNevents_spectrum; //weights given by number of events
		//theXmaxLikelihood/= fNevents_Xmax; //weights given by number of events
	}
	
	value= theSpectrumLikelihood + theXmaxLikelihood;

	fSpectrumLikelihood= theSpectrumLikelihood;
	fXmaxLikelihood= theXmaxLikelihood;
	fTotLikelihood= theSpectrumLikelihood + theXmaxLikelihood;

	if(fVerbosity>1){
  	cout<<"****************************************"<<endl;
  	cout<<"****  COMBINED LL FIT RESULTS    *******"<<endl;
  	cout<<"****************************************"<<endl;
  	cout<<"current fractions @ Source "<<endl;
  	for(int i=0;i<fNmassAtSource;i++){
			cout<<"FRACT "<<i+1<<" ==> "<<fFractionsAtSource[i]<<endl;
  	}	  
  	cout<<"current InjIndex ==> "<<fInjIndex<<endl;
  	cout<<"current Emax ==> "<<fEmax<<endl;
  	cout<<"current Flux Normalization ==> "<<fFluxNormFactor<<endl;
		if(fFitSourceCutoffShape) cout<<"current SourceCutoffShape ==> "<<fSourceCutoffShape<<endl;
		if(fFitGalSpectrum){
  		cout<<"current Gal InjIndex ==> "<<fGalacticInjIndex_spectrum<<endl;
  		cout<<"current Gal Emax ==> "<<fGalacticEmax_spectrum<<endl;
  		cout<<"current Gal Flux Normalization ==> "<<fGalacticFluxNormFactor_spectrum<<endl; 
		}	
		cout<<"current SpectrumLikelihood ==> "<<theSpectrumLikelihood<<endl;
  	cout<<"current XmaxLikelihood ==> "<<theXmaxLikelihood<<endl;
  	cout<<"current Likelihood ==> "<<value<<endl;
	}

}//close SourceFitter::MaxLikelihoodFcnCombined()


bool SourceFitter::CombinedFit(int verbosity) {

  if(verbosity>0 ) cout << " Spectrum-Xmax combined fit " << endl;
  int ierflag;
  double arglist[2];  
  int nPar;
	int nFittedFract= fNmassAtSource-1;

	nPar= 3;// EG InjIndex + EG Emax + EG Norm 
	if(fFitExtragalFractions) nPar+= nFittedFract;
	if(fFitSourceCutoffShape) nPar++;
	if(fFitGalSpectrum) nPar+= 3;// Gal InjIndex + Gal Emax + Gal Norm
	
  cout << " CombinedFit: Number of parameters ==> " <<nPar<< endl;
 
  TMinuit theMinuit(nPar); 
  theMinuit.SetPrintLevel(0);
  theMinuit.SetMaxIterations(10000);
  theMinuit.SetFCN(SourceFitter::MaxLikelihoodFcnCombined);
  
  
  theMinuit.mnexcm("SET NOW", arglist ,0,ierflag);
  
  const double stepSize    = 0.1;
  const double minFraction = 0.;
  const double maxFraction = 1.;
  const double startValue  = 1./(nFittedFract+1);

  const double startEmax	 = fEmax_start;
  const double stepEmax    = 0.1;
	const double minEmax     = 18.0;
  const double maxEmax     = 21.0;
    
  const double startInjIndex  = fInjIndex_start;
  const double stepInjIndex   = 0.1;
	const double minInjIndex	  = 0.0;
  const double maxInjIndex	  = 4.0;

	const double startSourceCutoffShape = fSourceCutoffShape_start;
	const double stepSourceCutoffShape  = 1;
  const double minSourceCutoffShape	  = 0.0001;
  const double maxSourceCutoffShape   = 100.;


	double NormStartValue  = fFluxNormFactor_start;
  const double NormStepSize    = 1; 
  const double minNorm         = 0.001;
  const double maxNorm         = 100.;

	
  //change norm according to index
  if(fInjIndex>=0.0&&fInjIndex<0.3) {
      NormStartValue = -35;
  }
  if(fInjIndex>=0.3&&fInjIndex<0.5) {   
      NormStartValue = -30;
  }
  if(fInjIndex>=0.5&&fInjIndex<0.7) {
      NormStartValue = -25;
  }
  if(fInjIndex>=0.7&&fInjIndex<0.9) {
      NormStartValue = -20;
  }
  if(fInjIndex>=0.9&&fInjIndex<1.1) {
      NormStartValue = -15;     
  }
  if(fInjIndex>=1.1&&fInjIndex<1.3) {
      NormStartValue = -10;
  }
  if(fInjIndex>=1.3&&fInjIndex<1.5) {
      NormStartValue = -5; 
  }     
  if(fInjIndex>=1.5&&fInjIndex<1.7) {
      NormStartValue = -2; 
  }
	if(fInjIndex>=1.7&&fInjIndex<1.9) {
      NormStartValue = 2;  
  }
  if(fInjIndex>=1.9&&fInjIndex<2.1) {
      NormStartValue = 5;       
  }
  if(fInjIndex>=2.1&&fInjIndex<2.3) {
      NormStartValue = 8;  
  }
  if(fInjIndex>=2.3&&fInjIndex<2.5) {
      NormStartValue = 11; 
  }     
  if(fInjIndex>=2.5&&fInjIndex<3.) { 
      NormStartValue = 20; 
  }
	if(fInjIndex>=3&&fInjIndex<3.5) {  
      NormStartValue = 29; 
  }
  if(fInjIndex>=3.5&&fInjIndex<4.) { 
      NormStartValue = 37;      
  }
	
	
  fFluxNormFactor_start= NormStartValue;
	

  //### GALACTIC PARAMETERS ###
  fGalacticEmax_spectrum= 18.0;
  const double startGalEmax	 = fGalacticEmax_start;
  const double stepGalEmax    = 0.1;
	const double Zgal           = 26.;
  const double minGalEmax     = 16.00;
  const double maxGalEmax     = 17.10;
    
  fGalacticInjIndex_spectrum= 2.5;
  const double startGalInjIndex  = fGalacticInjIndex_start;
  const double stepGalInjIndex   = 0.1;
  const double minGalInjIndex	 = 1.;
  const double maxGalInjIndex	 = 5.0;

  fGalacticFluxNormFactor_spectrum= 17.;
	const double GalNormStartValue  = fGalacticFluxNormFactor_start;
	const double GalNormStepSize    = 1; 
  const double minGalNorm         = 0.001;
  const double maxGalNorm         = 100;


	
	//########################################
	//# Store initial values of the fractions
	//########################################
	theStartMassFractStruct.fSourceFraction[0]= fStartFractParamsAtSource[0];
  double ConstrFactor= (1.-fStartFractParamsAtSource[0]);
  for(int i=1; i< fNmassAtSource-1; i++){
      theStartMassFractStruct.fSourceFraction[i]= fStartFractParamsAtSource[i]*ConstrFactor;
      ConstrFactor*= (1.-fStartFractParamsAtSource[i]);
  }
  theStartMassFractStruct.fSourceFraction[fNmassAtSource-1]= ConstrFactor;	

	
  //convert into start values for fraction params
	double startFractParamValue[fNmassAtSource-1];
	
	
	startFractParamValue[0]= startValue;
	ConstrFactor= 1.-startFractParamValue[0];
	for(int k=1;k<fNmassAtSource-1;k++){					
		startFractParamValue[k]= startValue/ConstrFactor;//equal values
		ConstrFactor*= (1.-startFractParamValue[k]);
	}
	

	//########################################
	//# SET FIT PARAMS
	//########################################
	int par_counter=0;

  if(fFitExtragalFractions){
  	int thisPar=0;
  	for(int iPar=0;iPar<nFittedFract;iPar++){
  		ostringstream parName;
    	parName << "f"<<"_" << iPar;
     
			theMinuit.mnparm(thisPar,parName.str().c_str(),startFractParamValue[iPar],stepSize,minFraction,maxFraction,ierflag);
			fStartFractFitParam[iPar]= fStartFractParamsAtSource[iPar]; 

			//# fix the fractions?
			if(fFixMassFractions){
				theMinuit.FixParameter(thisPar);// fix mass fractions
			}//close if FixMassFractions

    	thisPar++;
			par_counter++;
  	}//end loop extragal fractions
	}//close if fFitExtragalFractions		
	

	//## ADD EXTRAGALACTIC PARAMETERS	
 	//add inj index
  theMinuit.mnparm(par_counter,"InjIndex",startInjIndex,stepInjIndex,minInjIndex,maxInjIndex,ierflag);
	if(fFixInjIndex) theMinuit.FixParameter(par_counter);// fix InjIndex
	par_counter++;


  //add Emax
  theMinuit.mnparm(par_counter,"Emax",startEmax,stepEmax,minEmax,maxEmax,ierflag); 
	if(fFixEmax) theMinuit.FixParameter(par_counter);// fix Emax
	par_counter++;

 
 	//add normalization to fit parameters
  theMinuit.mnparm(par_counter,"Flux Norm",NormStartValue,NormStepSize,-60,60,ierflag);
	par_counter++;

	//add cutoff shape
	if(fFitSourceCutoffShape) {
		theMinuit.mnparm(par_counter,"SourceCutoffShape",startSourceCutoffShape,stepSourceCutoffShape,minSourceCutoffShape,maxSourceCutoffShape,ierflag);
		if(fFixSourceCutoffShape) theMinuit.FixParameter(par_counter);//fix source cutoff width  		
		par_counter++;  
	}


  //## ADD GALACTIC PARAMETERS
	if(fFitGalSpectrum){			
  	//add inj index
  	theMinuit.mnparm(par_counter,"InjIndexGalactic",startGalInjIndex,stepGalInjIndex,minGalInjIndex,maxGalInjIndex,ierflag);
		par_counter++;  		

		//add Emax @ source
  	theMinuit.mnparm(par_counter,"EmaxGalactic",startGalEmax,stepGalEmax,minGalEmax,maxGalEmax,ierflag);  
		par_counter++;  		

		//add normalization to fit parameters
  	theMinuit.mnparm(par_counter,"Gal Flux Norm",GalNormStartValue,GalNormStepSize,5,100,ierflag);//range 5-100
		par_counter++;
	}//close if 

		
	
	
	//## SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
  arglist[0]= 1;
  theMinuit.mnexcm("SET STR",arglist,1,ierflag);

	//## SET MINUIT ERROR
  arglist[0]= 0.5;//0.5 likelihood, 1 ChiSquare
  theMinuit.mnexcm("SET ERR",arglist,1,ierflag);


	//## SET MINUIT VERBOSITY
  //arglist[0]= 3;
  //theMinuit.mnexcm("SET PRI",arglist,1,ierflag);

	//SET MINUIT COMMAND
	arglist[0]= 10000;
	//arglist[1]= 0.5;
  arglist[1]= 0.1;

	theMinuit.mnexcm("MINIMIZE", arglist, 2, ierflag);//use MINIMIZE
	//theMinuit.mnexcm("MIGRAD", arglist ,2,ierflag);//use MIGRAD 
	//theMinuit.mnexcm("HESS", arglist ,2,ierflag);//use HESS
		
	//improve local minimum ==> try to exit from bad local minima (hopefully!)
  //theMinuit.mnexcm("IMPROVE", arglist ,0,ierflag);
 
	fFitStatus= theMinuit.GetStatus();
	fFitStatus= ierflag;

  double LikelihoodMin, edm, errdef;
  int nvpar, nparx, icstat;
  theMinuit.mnstat(LikelihoodMin, edm, errdef, nvpar, nparx, icstat);
  
	//Get covariance matrix
	double covMatrix[nPar][nPar];		
  theMinuit.mnemat(&covMatrix[0][0],nPar);
  	
  cout<<"*** COV MATRIX ***"<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
		
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i

	//Store cov matrix
	fCovMatrix.clear();
	fCovMatrix.resize(0);
	fCovMatrix.push_back(CovarianceMatrix);


	//######################
	//### GET PARAMS     ###
	//######################  
  double p,ep;
  par_counter= 0;
	fNmassFree= fNmassAtSource-1;

	fitParameterList.clear();
	fitParameterList.resize(0);
	fitParameterErrorList.clear();
	fitParameterErrorList.resize(0);

	for(int i=0;i<nPar;i++){
		theMinuit.GetParameter(i,p,ep);	
		fitParameterList.push_back(p);
		fitParameterErrorList.push_back(ep);
	}


	if(fFitExtragalFractions){
  	theMinuit.GetParameter(0,p,ep);
		fFractFitParam[0]= p;
		fFractFitParamErr[0]= ep;

		double UnconstrainedFraction[fNmassAtSource-1];
  	UnconstrainedFraction[0]= p;

  	fFractionsAtSource[0]= UnconstrainedFraction[0];
		fSourceFract[0]= UnconstrainedFraction[0];
  	double ConstrFactor= (1.-UnconstrainedFraction[0]);
	
  	for(int i=1;i<nFittedFract;i++){
			theMinuit.GetParameter(i,p,ep);
			fFractFitParam[i]= p;
			fFractFitParamErr[i]= ep;
			UnconstrainedFraction[i]= p;
    	fFractionsAtSource[i]= UnconstrainedFraction[i]*ConstrFactor;
			fSourceFract[i]= UnconstrainedFraction[i]*ConstrFactor;
    	ConstrFactor*= (1.-UnconstrainedFraction[i]);
  	}//close for i
  	fFractionsAtSource[fNmassAtSource-1]= ConstrFactor;
		fSourceFract[fNmassAtSource-1]= ConstrFactor;

		for(int i=0;i<fNmassAtSource;i++) theMassFractStruct.fSourceFraction[i]= fFractionsAtSource[i];
	
		//Calculate mass fraction uncertainties
		fitParameterErrorList= GetParamUncertainty(fitParameterList,CovarianceMatrix,nPar+1);
		for(unsigned int k=0;k<fitParameterErrorList.size();k++) cout<<"ParErr "<<k<<" ==> "<<fitParameterErrorList[k]<<endl;

		
		//Calculate low and upper bound
		double ErrLowBound=0.;
  	double ErrHighBound=0.;
		for(int i=0;i<fNmassAtSource;i++){
	  	ErrHighBound= fmin(1.,fFractionsAtSource[i]+fitParameterErrorList[i])-fFractionsAtSource[i];
 	  	ErrLowBound= fFractionsAtSource[i]-fmax(0.,fFractionsAtSource[i]-fitParameterErrorList[i]);
 	  	fFractionsAtSourceHighErr[i]= ErrHighBound; 
 	  	fFractionsAtSourceLowErr[i]= ErrLowBound;
			fSourceFractHighErr[i]= ErrHighBound;
			fSourceFractLowErr[i]= ErrLowBound;
			cout<<"Mass Fractions ==> "<<i+1<<"  "<<fFractionsAtSource[i]<<" +- "<< fFractionsAtSourceHighErr[i]<<"  "<<fFractionsAtSourceLowErr[i]<<endl;
  	}//close for i 

		par_counter+= nFittedFract;
	}//close if fFitExtragalFractions


	//get final InjIndex
  theMinuit.GetParameter(par_counter,p,ep);
  fInjIndex= p;
	fInjIndex_spectrum= p;
	fInjIndexErr_spectrum= ep;
	fInjIndex_Xmax= p;
	fInjIndexErr_Xmax= ep;
	par_counter++;

  //get final Emax
	theMinuit.GetParameter(par_counter,p,ep);
	fEmax= p;
	fEmax_spectrum= p;
	fEmaxErr_spectrum= ep;
	fEmax_Xmax= p;
	fEmaxErr_Xmax= ep;
	par_counter++;  

	//get final Norm
  theMinuit.GetParameter(par_counter,p,ep);
	fFluxNormFactor= p;
	fFluxNormFactor_spectrum= p;
	fFluxNormFactorErr_spectrum= ep;
	par_counter++;

	if(fFitGalSpectrum){
		//get final gal InjIndex
    theMinuit.GetParameter(par_counter,p,ep);
    fGalacticInjIndex_spectrum= p;
		fGalacticInjIndexErr_spectrum= ep;
		par_counter++;  

		//get final gal Emax
    theMinuit.GetParameter(par_counter,p,ep);
    fGalacticEmax_spectrum= p;
		fGalacticEmaxErr_spectrum= ep;
		par_counter++;
    
		//get final gal Norm
    theMinuit.GetParameter(par_counter,p,ep);
    fGalacticFluxNormFactor_spectrum= p;
		fGalacticFluxNormFactorErr_spectrum= ep;
		par_counter++;
	}//close if

	
  //get source cutoff width
	if(fFitSourceCutoffShape) {
		theMinuit.GetParameter(par_counter,p,ep);
  	fSourceCutoffShape= p; 
		fSourceCutoffShapeErr= ep;
		par_counter++;
	} 


  cout<<"FIT STATUS ==> "<<fFitStatus<<endl;

  
  cout<<"****************************************"<<endl;
  cout<<"****  COMBINED FIT RESULTS    **********"<<endl;
  cout<<"****************************************"<<endl;
  cout<<"final fractions @ Source "<<endl;
  for(unsigned int i=0;i<fFractionsAtSource.size();i++){
		cout<<"FRACT "<<i+1<<" ==> "<<fFractionsAtSource[i]<<endl;
  }	  
  cout<<"final InjIndex ==> "<<fInjIndex<<endl;
  cout<<"final Emax ==> "<<fEmax<<endl;
  cout<<"final Flux Normalization ==> "<<fFluxNormFactor<<endl;
	if(fFitSourceCutoffShape) cout<<"final SourceCutoffShape ==> "<<fSourceCutoffShape<<endl; 	
	if(fFitGalSpectrum){
  	cout<<"final Gal InjIndex ==> "<<fGalacticInjIndex_spectrum<<endl;
  	cout<<"final Gal Emax ==> "<<fGalacticEmax_spectrum<<endl;
  	cout<<"final Gal Flux Normalization ==> "<<fGalacticFluxNormFactor_spectrum<<endl;	
	}
	cout<<"final SpectrumLikelihood ==> "<<fSpectrumLikelihood<<endl;
  cout<<"final XmaxLikelihood ==> "<<fXmaxLikelihood<<endl;
  cout<<"final Likelihood ==> "<<fTotLikelihood<<endl;

 
  fSpectrumMinimizationInfo->Fill();
  fXmaxMinimizationInfo->Fill();
  

  //###########################
  //#####   STORE HISTO    ####
  //###########################
	StoreXmaxFitResults();
	StoreElongationRate();
	StoreRMS();
	StoreMassFraction();
	StoreSpectrumFitResults();
	
  return true;

}//close function




void SourceFitter::DoSourceFit(){
	 
	Init();
	SetData();

	fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
	fOutputFile->cd();
  
	fSpectrumMinimizationInfo= new TTree("SpectrumMinimizationInfo","SpectrumMinimizationInfo");
	fSpectrumMinimizationInfo->Branch("Nmass",&fNmass,"Nmass/I");
	fSpectrumMinimizationInfo->Branch("NmassFree",&fNmassFree,"NmassFree/I");
	fSpectrumMinimizationInfo->Branch("SourceFract",fSourceFract,"SourceFract[Nmass]/D");
	fSpectrumMinimizationInfo->Branch("SourceFractLowErr",fSourceFractLowErr,"SourceFractLowErr[Nmass]/D");  
	fSpectrumMinimizationInfo->Branch("SourceFractHighErr",fSourceFractHighErr,"SourceFractHighErr[Nmass]/D");  
	fSpectrumMinimizationInfo->Branch("FractFitParam",fFractFitParam,"FractFitParam[NmassFree]/D");
	fSpectrumMinimizationInfo->Branch("FractFitParamErr",fFractFitParamErr,"FractFitParamErr[NmassFree]/D");   	
	
	fSpectrumMinimizationInfo->Branch("InjIndex",&fInjIndex_spectrum,"InjIndex/D");
  fSpectrumMinimizationInfo->Branch("InjIndexErr",&fInjIndexErr_spectrum,"InjIndexErr/D");
	fSpectrumMinimizationInfo->Branch("Emax",&fEmax_spectrum,"Emax/D");
	fSpectrumMinimizationInfo->Branch("EmaxErr",&fEmaxErr_spectrum,"EmaxErr/D");	
	fSpectrumMinimizationInfo->Branch("FluxNormalization",&fFluxNormFactor_spectrum,"FluxNormalization/D");
	fSpectrumMinimizationInfo->Branch("FluxNormalizationErr",&fFluxNormFactorErr_spectrum,"FluxNormalizationErr/D");
	fSpectrumMinimizationInfo->Branch("SourceCutoffShape",&fSourceCutoffShape,"SourceCutoffShape/D");
	fSpectrumMinimizationInfo->Branch("SourceCutoffShapeErr",&fSourceCutoffShapeErr,"SourceCutoffShapeErr/D");
	
	fSpectrumMinimizationInfo->Branch("RedShiftEv",&fRedShiftEvolution,"RedShiftEv/D");

	fSpectrumMinimizationInfo->Branch("GalInjIndex",&fGalacticInjIndex_spectrum,"GalInjIndex/D");	
	fSpectrumMinimizationInfo->Branch("GalInjIndexErr",&fGalacticInjIndexErr_spectrum,"GalInjIndexErr/D");
	fSpectrumMinimizationInfo->Branch("GalEmax",&fGalacticEmax_spectrum,"GalEmax/D");
	fSpectrumMinimizationInfo->Branch("GalEmaxErr",&fGalacticEmaxErr_spectrum,"GalEmaxErr/D");
	fSpectrumMinimizationInfo->Branch("GalFluxNormalization",&fGalacticFluxNormFactor_spectrum,"GalFluxNormalization/D");	
	
	fSpectrumMinimizationInfo->Branch("SourceCutoffShape",&fSourceCutoffShape,"fSourceCutoffShape/D");
	fSpectrumMinimizationInfo->Branch("SourceCutoffShapeErr",&fSourceCutoffShapeErr,"fSourceCutoffShapeErr/D");

			
	fSpectrumMinimizationInfo->Branch("GalFluxNormalizationErr",&fGalacticFluxNormFactorErr_spectrum,"GalFluxNormalizationErr/D");
	fSpectrumMinimizationInfo->Branch("SpectrumLikelihood",&fSpectrumLikelihood,"SpectrumLikelihood/D");
	fSpectrumMinimizationInfo->Branch("SpectrumLikelihoodWeight",&fSpectrumLikelihoodWeight,"SpectrumLikelihoodWeight/D"); 	
	fSpectrumMinimizationInfo->Branch("FitStatus",&fFitStatus,"FitStatus/I");

	fSpectrumMinimizationInfo->Branch("fSourceFraction",theMassFractStruct.fSourceFraction,"theMassFractStruct.fSourceFraction[Nmass]/D"); 
  fSpectrumMinimizationInfo->Branch("fStartSourceFraction",theStartMassFractStruct.fSourceFraction,"theStartMassFractStruct.fSourceFraction[Nmass]/D"); 

	fSpectrumMinimizationInfo->Branch("StartFractFitParam",fStartFractFitParam,"StartFractFitParam[NmassFree]/D");
	fSpectrumMinimizationInfo->Branch("StartInjIndex",&fInjIndex_start,"StartInjIndex/D");
  fSpectrumMinimizationInfo->Branch("StartEmax",&fEmax_start,"StartEmax/D");
  fSpectrumMinimizationInfo->Branch("StartFluxNorm",&fFluxNormFactor_start,"StartFluxNorm/D");
  fSpectrumMinimizationInfo->Branch("StartGalInjIndex",&fGalacticInjIndex_start,"StartGalInjIndex/D");
	fSpectrumMinimizationInfo->Branch("StartGalEmax",&fGalacticEmax_start,"StartGalEmax/D");
  fSpectrumMinimizationInfo->Branch("StartGalFluxNorm",&fGalacticFluxNormFactor_start,"StartGalFluxNorm/D");


	fXmaxMinimizationInfo= new TTree("XmaxMinimizationInfo","XmaxMinimizationInfo");
	fXmaxMinimizationInfo->Branch("Nmass",&fNmass,"Nmass/I");
	fXmaxMinimizationInfo->Branch("NmassFree",&fNmassFree,"NmassFree/I");
	fXmaxMinimizationInfo->Branch("SourceFract",fSourceFract,"SourceFract[Nmass]/D");
	fXmaxMinimizationInfo->Branch("SourceFractLowErr",fSourceFractLowErr,"SourceFractLowErr[Nmass]/D");  
	fXmaxMinimizationInfo->Branch("SourceFractHighErr",fSourceFractHighErr,"SourceFractHighErr[Nmass]/D");  
	fXmaxMinimizationInfo->Branch("FractFitParam",fFractFitParam,"FractFitParam[NmassFree]/D");
	fXmaxMinimizationInfo->Branch("FractFitParamErr",fFractFitParamErr,"FractFitParamErr[NmassFree]/D");   	
	
	fXmaxMinimizationInfo->Branch("InjIndex",&fInjIndex_Xmax,"InjIndex/D");
	fXmaxMinimizationInfo->Branch("InjIndexErr",&fInjIndexErr_Xmax,"InjIndexErr/D");
	fXmaxMinimizationInfo->Branch("Emax",&fEmax_Xmax,"Emax/D");
	fXmaxMinimizationInfo->Branch("EmaxErr",&fEmaxErr_Xmax,"EmaxErr/D");
	fXmaxMinimizationInfo->Branch("RedShiftEv",&fRedShiftEvolution,"RedShiftEv/D");
	fXmaxMinimizationInfo->Branch("SourceCutoffShape",&fSourceCutoffShape,"fSourceCutoffShape/D");
	fXmaxMinimizationInfo->Branch("SourceCutoffShapeErr",&fSourceCutoffShapeErr,"fSourceCutoffShapeErr/D");
		
	fXmaxMinimizationInfo->Branch("XmaxLikelihood",&fXmaxLikelihood,"XmaxLikelihood/D");
	fXmaxMinimizationInfo->Branch("XmaxLikelihoodWeight",&fXmaxLikelihoodWeight,"XmaxLikelihoodWeight/D");  
  fXmaxMinimizationInfo->Branch("FitStatus",&fFitStatus,"FitStatus/I");
    
  fXmaxMinimizationInfo->Branch("fStartSourceFraction",theStartMassFractStruct.fSourceFraction,"theStartMassFractStruct.fSourceFraction[Nmass]/D"); 

	fXmaxMinimizationInfo->Branch("StartInjIndex",&fInjIndex_start,"StartInjIndex/D");
	fXmaxMinimizationInfo->Branch("StartEmax",&fEmax_start,"StartEmax/D");
   
	if(fIsCombinedFit){
		CombinedFit(3);
	}
	else{
		FitEnergySpectrum(3);
		//FitXmax(3);
	}
		
	cout<<"Writing trees in file"<<endl;
	fOutputFile->cd();
	fSpectrumMinimizationInfo->Write();
	if(fIsCombinedFit) fXmaxMinimizationInfo->Write();
 			
	fOutputFile->Close();

	cout<<"*** END SOURCE FITTER ***"<<endl;

}//close SourceFitter::DoSourceFit()


double SourceFitter::GetDerivativeMatrixElement(int i,int j,std::vector<double> parameters){

	//index i is relative to real mass fractions
	//index j is relative to fit mass fraction params
	//int si= (int)(i/fNmassAtSource); 
	//int ii= i-si*fNmassAtSource;
	//int sj= (int)(j/(fNmassAtSource-1));
	//int jj= j-sj*(fNmassAtSource-1);

	//calculate matrix elements according to the general indexes
	double matrixValue= 0.;

	//mass fit params
	std::vector<double> MassParams;
	MassParams.clear();
	MassParams.resize(fNmassAtSource-1);
	std::vector<double> MassFract;
	MassFract.clear();
	MassFract.resize(fNmassAtSource);
	double ConstrFactor;

	//copy from full vector
	for(unsigned int k=0;k<MassParams.size();k++){		
		MassParams[k]= parameters[k];
	}

	MassFract[0]= MassParams[0];
  ConstrFactor= (1.-MassParams[0]);

	for(int k=1;k<fNmassAtSource-1;k++){   
    MassFract[k]= MassParams[k]*ConstrFactor;
    ConstrFactor*= (1.-MassParams[k]);   
  }//end loop masses at source
	MassFract[fNmassAtSource-1]= ConstrFactor;
	

	if(i<fNmassAtSource && j<fNmassAtSource-1){
		//fract vs fract params
		if(j==i && j==0) matrixValue= 1.;
		else if(j==i && j!=0) matrixValue= MassFract[i]/MassParams[j];
		else if(j>i) matrixValue= 0.;	
		else if(j<i) matrixValue= -MassFract[i]/(1.-MassParams[j]); 
	}//close if
	else if(i<fNmassAtSource && j>=fNmassAtSource-1){
		//fract vs other params
		matrixValue= 0.;	
	}
	else if(i>=fNmassAtSource && j<fNmassAtSource-1){
		//other params vs fract params
		matrixValue= 0.;	
	}
	else if(i>=fNmassAtSource && j>=fNmassAtSource-1){
		//other params vs other params
		if(j==i-1) matrixValue= 1.;
		else matrixValue= 0.;	
	}
	else{
		cerr<<"Invalid index for derivative calculation!..please check!"<<endl;
		exit(1);
	}

	return matrixValue;

}//close SourceFitter::GetDerivativeMatrixElement()


std::vector<double> SourceFitter::GetParamUncertainty(std::vector<double> parameters,TMatrixD CovarianceMatrix,int Npar){

	unsigned int Npar_fit= parameters.size();	

	int Ndim= CovarianceMatrix.GetNrows();
	//check dimensions
	if(Npar_fit!=(unsigned int)(Ndim)){
		cerr<<"SourceFitter::GetParamUncertainty(): Error Npar_fit!=Ndim ... exit!"<<endl;
		exit(1);
	}//close if

	
	std::vector<double> parErrors;
	parErrors.resize(Npar);
	for(unsigned int i=0;i<parErrors.size();i++) parErrors[i]=0.;

	//calculate derivative matrix
	TMatrixD DerMatrix(Npar,Npar_fit);
	for(int i=0;i<Npar;i++){
		for(unsigned int j=0;j<Npar_fit;j++){
			DerMatrix(i,j)= GetDerivativeMatrixElement(i,j,parameters);
		}//close for j
	}//close for i

	//print der matrix
	//cout<<"Printing Derivative Matrix"<<endl;
	//DerMatrix.Print();
	//cout<<endl;

	//cout<<"Printing Covariance Matrix"<<endl;
	//CovarianceMatrix.Print();
	//cout<<endl;

	TMatrixD DerMatrixTransp(Npar_fit,Npar);
	DerMatrixTransp.Transpose(DerMatrix);//transpose of DerivativeMatrix


	TMatrixD tmp(Npar,Npar_fit);
	tmp.Mult(DerMatrix,CovarianceMatrix);

	TMatrixD FinalErrMatrix(Npar,Npar);
	FinalErrMatrix.Mult(tmp,DerMatrixTransp); //store Final Correct Covariance Matrix

  //cout<<"Printing Error Matrix"<<endl;
	//FinalErrMatrix.Print();
	//cout<<endl;

	//get diagonal
	for(int i=0;i<Npar;i++) parErrors[i]= sqrt(FinalErrMatrix(i,i));

	return parErrors;

}//close GetParamUncertainty()


void SourceFitter::StoreXmaxFitResults(){

	fOutputFile->cd();

	for(int i=0;i<fNmassAtEarth;i++) {   
		fXmaxMC[i]= fDataReader->GetXmaxMC(i);
	  fXmaxGenMC[i]= fDataReader->GetGenXmaxMC(i);		
	  fEnergyMC[i]= fDataReader->GetEnergyMC(i);
		fGenEnergyMC[i]= fDataReader->GetGenEnergyMC(i);
	}//end loop masses at Earth

	fFractionsAtEarth_spectrumfit= fPropMCReader->GetExpFractionsAtEarth_Spectrum();
	fFractionsAtEarth_Xmaxfit= fPropMCReader->GetExpFractionsAtEarth_Xmax();
	
  double XmaxDataKurtosis[fNbins_Xmax];
  double XmaxKurtosis[fNbins_Xmax];
  double XmaxTrueKurtosis[fNbins_Xmax];
  double NevXmax[fNbins_Xmax];

	TH1D* XmaxMC_scaled;
 	TH1D* XmaxFit;
	fXmaxFit.clear();
	fXmaxFit.resize(0);
	fXmaxTrueFit.clear();
	fXmaxTrueFit.resize(0);
	for(int i=0;i<fNmassAtEarth;i++){
		fXmaxMC_scaled[i].clear();
		fXmaxMC_scaled[i].resize(0);		
	}
	

	// *********************
  // ** STORE Xmax histo
  // *********************
  for(int s=0;s<fNbins_Xmax;s++){
		//### DATA INFO
		fXmaxData[s]->SetName(Form("XmaxData_bin%d",s+1));
		fXmaxData[s]->SetLineColor(kBlack);
		fXmaxData[s]->SetMarkerColor(kBlack);
		fXmaxData[s]->SetMarkerStyle(8);
		fXmaxData[s]->SetDrawOption("E1p");
		fXmaxData[s]->Write();
		
		NevXmax[s]= fXmaxData[s]->GetEntries();
	  
		XmaxDataKurtosis[s]= fXmaxData[s]->GetKurtosis();
		double datarms= fXmaxData[s]->GetRMS();
		double datam4= (XmaxDataKurtosis[s]+3)*pow(datarms,4);
		double datarmserr= sqrt((1/NevXmax[s])*(datam4-pow(datarms,4)*(NevXmax[s]-3)/(NevXmax[s]-1)))/(2.*datarms);
				
		
		//Fill data vectors
		fMeanXmaxData[s]= fXmaxData[s]->GetMean();
		fMeanXmaxErrData[s]= fXmaxData[s]->GetMeanError();
		fRMSData[s]= datarms;
		fRMSErrData[s]= datarmserr;
		fMeanEnergyData[s]= fEnergyData[s]->GetMean();
		fMeanEnergyErrData[s]= fEnergyData[s]->GetMeanError();

	    
	  double NormFactor=0.;  
		double ERRaw=0.;
		double ERRawErr=0.;
		double ERTrue=0.;
		double ERTrueErr=0.;
		
		//MassFractions are @ top of atmosphere
		//Calculate fractions @ detector level
		for(int i=0;i<fNmassAtEarth;i++){
			fFractionsAtDetector_Xmaxfit[i][s]= fFractionsAtEarth_Xmaxfit[i][s];
			fFractionsAtDetectorLowErr_Xmaxfit[i][s]= 0.;
			fFractionsAtDetectorHighErr_Xmaxfit[i][s]= 0.;
		}
		
		
  
		//### MC INFO
		for(int i=0;i<fNmassAtEarth;i++){
			fXmaxMC[i][s]->SetName(Form("XmaxMC%d_bin%d",i+1,s+1));
			TH1D* XmaxMC_scaled= (TH1D*)fXmaxMC[i][s]->Clone(Form("XmaxMC%d_bin%d_scaled",i+1,s+1));

			double scale= (fXmaxData[s]->Integral()/fXmaxMC[i][s]->Integral())* fFractionsAtEarth_Xmaxfit[i][s];
			XmaxMC_scaled->Scale(scale);			
			fXmaxMC_scaled[i].push_back(XmaxMC_scaled);//add to vector
			

			if(i==0){
				fXmaxMC[i][s]->SetLineColor(kRed);
				fXmaxMC[i][s]->SetMarkerColor(kRed);
				fXmaxMC_scaled[i][s]->SetLineColor(kRed);
				fXmaxMC_scaled[i][s]->SetMarkerColor(kRed);
			}
			else if(i==1){
				fXmaxMC[i][s]->SetLineColor(kCyan);
				fXmaxMC[i][s]->SetMarkerColor(kCyan);
				fXmaxMC_scaled[i][s]->SetLineColor(kCyan);
				fXmaxMC_scaled[i][s]->SetMarkerColor(kCyan);
			}
			else if(i==2){
				fXmaxMC[i][s]->SetLineColor(kGreen);
				fXmaxMC[i][s]->SetMarkerColor(kGreen);
				fXmaxMC_scaled[i][s]->SetLineColor(kGreen);
				fXmaxMC_scaled[i][s]->SetMarkerColor(kGreen);
			}
			else if(i==3){
				fXmaxMC[i][s]->SetLineColor(kBlue);
				fXmaxMC[i][s]->SetMarkerColor(kBlue);
				fXmaxMC_scaled[i][s]->SetLineColor(kBlue);
				fXmaxMC_scaled[i][s]->SetMarkerColor(kBlue);
			}			
			else{
				fXmaxMC[i][s]->SetLineColor(kBlack);
				fXmaxMC[i][s]->SetMarkerColor(kBlack);
				fXmaxMC_scaled[i][s]->SetLineColor(kBlack);
				fXmaxMC_scaled[i][s]->SetMarkerColor(kBlack);
			}
			
			//Write to file
			fXmaxMC[i][s]->SetDrawOption("hist");
			fXmaxMC[i][s]->Write();		
			fXmaxMC_scaled[i][s]->SetDrawOption("hist");
			fXmaxMC_scaled[i][s]->Write();

				
			//Store MC vectors
			fMeanXmaxMC[i][s]= fXmaxMC[i][s]->GetMean();
			fMeanXmaxErrMC[i][s]= fXmaxMC[i][s]->GetMeanError();
			fRMSMC[i][s]= fXmaxMC[i][s]->GetRMS();
			fRMSErrMC[i][s]= fXmaxMC[i][s]->GetRMSError();
			fMeanEnergyMC[i][s]= fEnergyMC[i][s]->GetMean();
			fMeanEnergyErrMC[i][s]= fEnergyMC[i][s]->GetMeanError();
				
			//Store pure conex info
			fMeanXmaxConexMC[i][s]= fXmaxGenMC[i][s]->GetMean();
			fMeanXmaxErrConexMC[i][s]= fXmaxGenMC[i][s]->GetMeanError();
	    fXmaxRMSConexMC[i][s]= fXmaxGenMC[i][s]->GetRMS();
			fXmaxRMSErrConexMC[i][s]= fXmaxGenMC[i][s]->GetRMSError();
						      
			//Calculate ER
			ERRaw+= fMeanXmaxMC[i][s]*fFractionsAtDetector_Xmaxfit[i][s];//RAW ER      
			ERTrue+= fMeanXmaxConexMC[i][s]*fFractionsAtEarth_Xmaxfit[i][s];//TRUE ER
     			
		}//end loop masses at Earth
		
		//Fill ER vectors
		fERRaw[s]= ERRaw;
		fERTrue[s]= ERTrue;

		//## Calculate fit histo
		TH1D* fit= (TH1D*)fXmaxMC_scaled[0][s]->Clone("fit");
		fit->SetName(Form("XmaxFit_bin%d",s+1));
		for(int i=1;i<fNmassAtEarth;i++) fit->Add(fXmaxMC_scaled[i][s]);	
		fit->SetLineColor(kBlack);
		fit->SetMarkerColor(kBlack);
		fit->SetDrawOption("hist");
		fit->Write();		
		fXmaxFit.push_back(fit);
		
		XmaxKurtosis[s]= fit->GetKurtosis();//get Xmax distr 4th moment
		
		
		//## Calculate true fit histo
		TH1D* truefit= (TH1D*)fXmaxGenMC[0][s]->Clone("truefit");;
		double scale= fFractionsAtEarth_Xmaxfit[0][s]*fXmaxData[s]->Integral()/fXmaxGenMC[0][s]->Integral();	
		truefit->Scale(scale);
		for(int i=1;i<fNmassAtEarth;i++) 
			truefit->Add(fXmaxGenMC[i][s],fFractionsAtEarth_Xmaxfit[i][s]*fXmaxData[s]->Integral()/fXmaxGenMC[i][s]->Integral());

		truefit->SetName(Form("XmaxTrueFit_bin%d",s+1));	
		truefit->SetLineColor(kBlack);
		truefit->SetMarkerColor(kBlack);
		truefit->SetDrawOption("hist");
		truefit->Write();		
		fXmaxTrueFit.push_back(truefit);			
		XmaxTrueKurtosis[s]= truefit->GetKurtosis();//get Xmax distr 4th moment
				
		//######################################
		//###  METHOD TO CALCULATE RMS ERROR ###
		//######################################
		//## to be done after Earth fraction err calculation		
		//## Calculate RMS uncertainty ==> RMSErr= sqr{(1/N)*[m4-RMS^4*(N-3)/(N-1)]}/2RMS		
		double rms= fXmaxFit[s]->GetRMS();
		double m4= (XmaxKurtosis[s]+3)*pow(rms,4);
		double truerms= fXmaxTrueFit[s]->GetRMS();
		double truem4= (XmaxTrueKurtosis[s]+3)*pow(truerms,4);
		double rmserr= sqrt((1/NevXmax[s])*(m4-pow(rms,4)*(NevXmax[s]-3)/(NevXmax[s]-1)))/(2.*rms);
		double rmstrueerr= sqrt((1/NevXmax[s])*(truem4-pow(truerms,4)*(NevXmax[s]-3)/(NevXmax[s]-1)))/(2.*truerms);		
									
		//Store fit vectors
		fMeanXmaxFit[s]= fXmaxFit[s]->GetMean();
		fMeanXmaxErrFit[s]= fXmaxFit[s]->GetMeanError();
		fRMSFit[s]= fXmaxFit[s]->GetRMS();
		//fRMSErrFit[s]= fXmaxFit[s]->GetRMSError();//Gaussian approximation
		fRMSErrFit[s]= rmserr;
				
		fRMSTrueFit[s]= fXmaxTrueFit[s]->GetRMS();
		//fRMSTrueErrFit[s]= fXmaxTrueFit[s]->GetRMSError(); //Gaussian approximation
		fRMSTrueErrFit[s]= rmstrueerr;
			
	}//end loop Xmax energy bins
	
	//## Calculate ER errors
	//## to be done after Earth fraction err calculation		
	for(int s=0;s<fNbins_Xmax;s++){		
		fERRawErr[s]= 0.;
		fERTrueErr[s]= 0.;	
	}//end loop energy bins

	
}//close SourceFitter::StoreFitResults()


void SourceFitter::StoreSpectrumFitResults(){

	fOutputFile->cd();

	//## Store Spectrum Data
	fAugerEnergySpectrumEvents->SetNameTitle("AugerEnergySpectrumEvents","AugerEnergySpectrumEvents");
	fAugerEnergySpectrumEvents->SetMarkerStyle(8);
	fAugerEnergySpectrumEvents->SetMarkerColor(kBlack);
	fAugerEnergySpectrumEvents->Write();

	fAugerEnergySpectrumE3->SetNameTitle("AugerEnergySpectrumE3","AugerEnergySpectrumE3");
	fAugerEnergySpectrumE3->SetMarkerStyle(8);
	fAugerEnergySpectrumE3->SetMarkerColor(kBlack);    
	fAugerEnergySpectrumE3->Write(); 

	fAugerEnergySpectrum->SetNameTitle("AugerEnergySpectrum","AugerEnergySpectrum");  
	fAugerEnergySpectrum->SetMarkerStyle(8);
	fAugerEnergySpectrum->SetMarkerColor(kBlack);
	fAugerEnergySpectrum->Write();
  
	//## Store fit histo
	TH1D* fExtragalSpectrumE3Fit= (TH1D*)fTotPropEnergySpectrum_spectrumfit->Clone("fExtragalSpectrumE3Fit"); 
  TH1D* fGalSpectrumE3Fit= (TH1D*)fGalacticSpectrum_spectrumfit->Clone("fGalSpectrumE3Fit");

  std::vector<TH1D*> fExtragalSpectrumE3Component(fPropMCReader->GetExpEnergySpectrumAtEarth_Spectrum());
	std::vector<TH1D*> fExtragalSpectrumE3DirectComponent(fPropMCReader->GetExpDirectEnergySpectrumAtEarth_Spectrum());

	//## Multiply dN/dlgE x E^2 to get E^3 dN/dE
	//## E^3 dN/dE= E^3 x 1/E x dN/dlgE= E^2 x dN/dlgE 
	for(int s=0;s<fExtragalSpectrumE3Fit->GetNbinsX();s++){   
		double bincontent= fExtragalSpectrumE3Fit->GetBinContent(s+1);	
		double bincenter= fExtragalSpectrumE3Fit->GetBinCenter(s+1);		
		double weight= pow(pow(10,bincenter),2);		
		fExtragalSpectrumE3Fit->SetBinContent(s+1,weight*bincontent);	
		
		bincontent= fGalSpectrumE3Fit->GetBinContent(s+1);	
		bincenter= fGalSpectrumE3Fit->GetBinCenter(s+1);	
		weight= pow(pow(10,bincenter),2);		
		fGalSpectrumE3Fit->SetBinContent(s+1,weight*bincontent);

		for(unsigned int i=0;i<fExtragalSpectrumE3Component.size();i++){
			bincontent= fExtragalSpectrumE3Component[i]->GetBinContent(s+1);			
			bincenter= fExtragalSpectrumE3Component[i]->GetBinCenter(s+1);		
			weight= pow(pow(10,bincenter),2);	
			fExtragalSpectrumE3Component[i]->SetBinContent(s+1,weight*bincontent);
			
			bincontent= fExtragalSpectrumE3DirectComponent[i]->GetBinContent(s+1);			
			bincenter= fExtragalSpectrumE3DirectComponent[i]->GetBinCenter(s+1);		
			weight= pow(pow(10,bincenter),2);	
			fExtragalSpectrumE3DirectComponent[i]->SetBinContent(s+1,weight*bincontent);
			
		}//end loop extragalactic masses
  }//end loop spectrum bins
 
  fExtragalSpectrumE3Fit->SetNameTitle("ExtragalSpectrumE3Fit","ExtragalSpectrumE3Fit");
  fExtragalSpectrumE3Fit->SetMarkerStyle(8);
  fExtragalSpectrumE3Fit->SetMarkerColor(kBlack);
  fExtragalSpectrumE3Fit->Write();

  fGalSpectrumE3Fit->SetNameTitle("GalSpectrumE3Fit","GalSpectrumE3Fit");
  fGalSpectrumE3Fit->SetMarkerStyle(8);
  fGalSpectrumE3Fit->SetMarkerColor(kBlack);
  fGalSpectrumE3Fit->Write();

	TH1D* fTotSpectrumE3Fit= (TH1D*)fExtragalSpectrumE3Fit->Clone("fTotSpectrumE3Fit");
  if(fFitGalSpectrum) {
	  fTotSpectrumE3Fit->Add(fGalSpectrumE3Fit);
  }
  
  fTotSpectrumE3Fit->SetNameTitle("TotSpectrumE3Fit","TotSpectrumE3Fit");  
  fTotSpectrumE3Fit->Write();
  

	for(unsigned int i=0;i<fExtragalSpectrumE3Component.size();i++) {
		fExtragalSpectrumE3Component[i]->SetName(Form("ExtragalSpectrumE3Component%d",i+1));
    fExtragalSpectrumE3Component[i]->SetTitle(Form("ExtragalSpectrumE3Component%d",i+1));
		fExtragalSpectrumE3Component[i]->Write();

		fExtragalSpectrumE3DirectComponent[i]->SetName(Form("ExtragalSpectrumE3DirectComponent%d",i+1));
    fExtragalSpectrumE3DirectComponent[i]->SetTitle(Form("ExtragalSpectrumE3DirectComponent%d",i+1));
		fExtragalSpectrumE3DirectComponent[i]->Write();
	}//end loop masses at Earth


}//close SourceFitter::StoreSpectrumFitResults()



void SourceFitter::StoreElongationRate(){

	double MeanEnergy[fNbins_Xmax];
	double MeanEnergyErr[fNbins_Xmax];
	double MeanXmaxData[fNbins_Xmax];
	double MeanXmaxDataErr[fNbins_Xmax];
	double ERRaw[fNbins_Xmax];
	double ERRawErr[fNbins_Xmax];
	double ERTrue[fNbins_Xmax];
	double ERTrueErr[fNbins_Xmax];

	for(int s=0;s<fNbins_Xmax;s++){
		MeanEnergy[s]= fMeanEnergyData[s];
		MeanEnergyErr[s]= fMeanEnergyErrData[s];
		MeanXmaxData[s]= fMeanXmaxData[s];
		MeanXmaxDataErr[s]= fMeanXmaxErrData[s];

		ERRaw[s]= fERRaw[s];
		ERRawErr[s]= fERRawErr[s];
		ERTrue[s]= fERTrue[s];
		ERTrueErr[s]= fERTrueErr[s];
	}

	TGraphErrors* ERData= new TGraphErrors(fNbins_Xmax,MeanEnergy,MeanXmaxData,MeanEnergyErr,MeanXmaxDataErr);
	ERData->SetMarkerStyle(8);
  ERData->SetMarkerSize(1.5);
  ERData->SetMarkerColor(kBlack);
  ERData->SetLineColor(kBlack);
	ERData->SetName("ERData");

	TGraphErrors* ERRawGraph= new TGraphErrors(fNbins_Xmax,MeanEnergy,ERRaw,MeanEnergyErr,ERRawErr);
	ERRawGraph->SetMarkerStyle(29);
  ERRawGraph->SetMarkerSize(2.);
  ERRawGraph->SetMarkerColor(kRed);
  ERRawGraph->SetLineColor(kRed);
	ERRawGraph->SetName("ERRaw");
	

	TGraphErrors* ERTrueGraph= new TGraphErrors(fNbins_Xmax,MeanEnergy,ERTrue,MeanEnergyErr,ERTrueErr);
	ERTrueGraph->SetMarkerStyle(8);
  ERTrueGraph->SetMarkerSize(2.);
  ERTrueGraph->SetMarkerColor(kBlack);
  ERTrueGraph->SetLineColor(kBlack);
	ERTrueGraph->SetName("ERTrue");

  
	fOutputFile->cd();
	ERData->Write();
	ERRawGraph->Write();
	ERTrueGraph->Write();

}//close SourceFitter::StoreElongationRate()


void SourceFitter::StoreRMS(){

	double MeanEnergy[fNbins_Xmax];
	double MeanEnergyErr[fNbins_Xmax];
	double XmaxRMSData[fNbins_Xmax];
	double XmaxRMSDataErr[fNbins_Xmax];

	double XmaxRMSFit[fNbins_Xmax];
	double XmaxRMSFitErr[fNbins_Xmax];
	
	double XmaxTrueRMSFit[fNbins_Xmax];
	double XmaxTrueRMSFitErr[fNbins_Xmax];
	
	for(int s=0;s<fNbins_Xmax;s++){
		MeanEnergy[s]= fMeanEnergyData[s];
		MeanEnergyErr[s]= fMeanEnergyErrData[s];
		XmaxRMSData[s]= fRMSData[s];
		XmaxRMSDataErr[s]= fRMSErrData[s];
		XmaxRMSFit[s]= fRMSFit[s];
		XmaxRMSFitErr[s]= fRMSErrFit[s];
		
		XmaxTrueRMSFit[s]= fRMSTrueFit[s];
		XmaxTrueRMSFitErr[s]= fRMSTrueErrFit[s];
	}

	TGraphErrors* RMSData= new TGraphErrors(fNbins_Xmax,MeanEnergy,XmaxRMSData,MeanEnergyErr,XmaxRMSDataErr);
	RMSData->SetName("RMSData");
	RMSData->SetMarkerStyle(8);
  RMSData->SetMarkerSize(1.5);
  RMSData->SetMarkerColor(kBlack);
  RMSData->SetLineColor(kBlack);

	TGraphErrors* RMSFit= new TGraphErrors(fNbins_Xmax,MeanEnergy,XmaxRMSFit,MeanEnergyErr,XmaxRMSFitErr);
	RMSFit->SetName("RMSFit");
	RMSFit->SetMarkerStyle(29);
  RMSFit->SetMarkerSize(2.);
  RMSFit->SetMarkerColor(kRed);
  RMSFit->SetLineColor(kRed);

	TGraphErrors* RMSTrueFit= new TGraphErrors(fNbins_Xmax,MeanEnergy,XmaxTrueRMSFit,MeanEnergyErr,XmaxTrueRMSFitErr);
	RMSTrueFit->SetName("RMSTrueFit");
	RMSTrueFit->SetMarkerStyle(8);
  RMSTrueFit->SetMarkerSize(2.);
  RMSTrueFit->SetMarkerColor(kBlack);
  RMSTrueFit->SetLineColor(kBlack);
 
  fOutputFile->cd();
	RMSData->Write();
	RMSFit->Write();
	RMSTrueFit->Write();
	

}//close SourceFitter::StoreRMS()


void SourceFitter::StoreMassFraction(){

	double Fract[fNbins_Xmax];
	double FractErrHigh[fNbins_Xmax];
	double FractErrLow[fNbins_Xmax];
	
	double Fract_DetLevel[fNbins_Xmax];
	double FractErrLow_DetLevel[fNbins_Xmax];
	double FractErrHigh_DetLevel[fNbins_Xmax];
	
	double MeanEnergy[fNbins_Xmax];
	double MeanEnergyErr[fNbins_Xmax];
	
  for(int i=0;i<fNmassAtEarth;i++){

		for(int s=0;s<fNbins_Xmax;s++){
			MeanEnergy[s]= fMeanEnergyData[s];
		  MeanEnergyErr[s]= fMeanEnergyErrData[s];
	
			Fract[s]= fFractionsAtEarth_Xmaxfit[i][s];
			FractErrHigh[s]= fFractionsAtEarthHighErr_Xmaxfit[i][s];
			FractErrLow[s]= fFractionsAtEarthLowErr_Xmaxfit[i][s];
			
			Fract_DetLevel[s]= fFractionsAtDetector_Xmaxfit[i][s];
			FractErrHigh_DetLevel[s]= fFractionsAtDetectorHighErr_Xmaxfit[i][s];
			FractErrLow_DetLevel[s]= fFractionsAtDetectorLowErr_Xmaxfit[i][s];
	 	}//end loop Xmax energy bins

	 	TGraphAsymmErrors* MassFractGraph= new TGraphAsymmErrors(fNbins_Xmax,MeanEnergy,Fract,0,0,FractErrLow,FractErrHigh);
	 	TString currentGraphName= Form("MassFraction%d",i+1);
	 	MassFractGraph->SetName(currentGraphName);
		
		TGraphAsymmErrors* MassFractDetGraph= new TGraphAsymmErrors(fNbins_Xmax,MeanEnergy,Fract_DetLevel,0,0,FractErrLow_DetLevel,FractErrHigh_DetLevel);
    currentGraphName=Form("MassDetFraction%d",i+1);
		MassFractDetGraph->SetName(currentGraphName);
		
		if(i==0){
			MassFractGraph->SetMarkerColor(kRed);
  		MassFractGraph->SetLineColor(kRed);
  		MassFractGraph->SetMarkerSize(1.5);
  		MassFractGraph->SetMarkerStyle(8);

			MassFractDetGraph->SetMarkerColor(kRed+2);
  		MassFractDetGraph->SetLineColor(kRed+2);
  		MassFractDetGraph->SetMarkerSize(1.3);
  		MassFractDetGraph->SetMarkerStyle(24);
		}
		else if(i==1){
			MassFractGraph->SetMarkerColor(kBlue);
  		MassFractGraph->SetLineColor(kBlue);
  		MassFractGraph->SetMarkerSize(2);
  		MassFractGraph->SetMarkerStyle(22);

			MassFractDetGraph->SetMarkerColor(kBlue+2);
  		MassFractDetGraph->SetLineColor(kBlue+2);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(26);
		}
		else if(i==2){
			MassFractGraph->SetMarkerColor(kGreen);
  		MassFractGraph->SetLineColor(kGreen);
  		MassFractGraph->SetMarkerSize(2);
  		MassFractGraph->SetMarkerStyle(29);

			MassFractDetGraph->SetMarkerColor(kGreen+2);
  		MassFractDetGraph->SetLineColor(kGreen+2);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(26);
		}
		else if(i==3){
			MassFractGraph->SetMarkerColor(kCyan);
  		MassFractGraph->SetLineColor(kCyan);
  		MassFractGraph->SetMarkerSize(2);
  		MassFractGraph->SetMarkerStyle(29);

			MassFractDetGraph->SetMarkerColor(kCyan+2);
  		MassFractDetGraph->SetLineColor(kCyan+2);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(26);
		}
		else {			
			MassFractGraph->SetMarkerColor(kBlack);
  		MassFractGraph->SetLineColor(kBlack);
  		MassFractGraph->SetMarkerSize(2);
  		MassFractGraph->SetMarkerStyle(29);

			MassFractDetGraph->SetMarkerColor(kBlack);
  		MassFractDetGraph->SetLineColor(kBlack);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(26);
		}
		
    fOutputFile->cd();
	  MassFractGraph->Write();
		MassFractDetGraph->Write();
		
	}//end loop masses at Earth

	
}//close SourceFitter::StoreMassFraction()

}//close namespace

