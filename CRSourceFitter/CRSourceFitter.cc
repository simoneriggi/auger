/**
* @mainpage SourceCompositionFitter tool

* <img src="SourceCompositionFitter_screenshot.jpg" alt="Screenshot">
* <b>Description of the project </b> <br>
* This project is developed to perform a combined fit of spectrum and Xmax data. <br>
The project settings can be specified in a configuration file which is parsed by the ConfigParser class.<br>
The MCReader and DataReader classes are provided to read and get data from ROOT files containing the data (Xmax, spectrum) and MC information (both reconstructed and generated). A reader is provided also to get data from a ROOT file containing the simulated propagation data. Its primary function is to build the transfer matrix used to propagate the spectrum from the source to the Earth, and to calculate spectra observed @ Earth. <br>
A global model with a set of parameters is assumed and fitted to spectrum+Xmax data. The fitted parameters are: injection index gamma, maximum acceleration energy Emax @ source, flux normalization, mass fractions. Eventually the shift of Xmax distribution can be included as fit parameter. The likelihood fit of Xmax, spectrum, Xmax+spectrum data is done by the SourceFitter class. Eventually the EarthFitter class perform the Xmax fit @ Earth level <br>
As utilities the Drawer class plus some other ROOT macros are provided to draw the obtained fit results and other operations.
* 
* @author S. Riggi <br> Dept. of Physics and Astronomy, University of Catania <br> INFN, Section of Catania  <br> Italy
*/

/**
* @file main.cc
* @class main
* @brief main class used to handle the entire project
*
* ...
* @author S. Riggi
* @date 22/04/2010
*/

#include "ConfigParser.h"
#include "SourceFitter.h"
#include "DrawFitResults.h"
#include "PropagationMCReader.h"
#include "AnalysisConsts.h"
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
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>

#include<TPaveText.h>
#include<TVirtualFitter.h>
#include<TObjArray.h>
#include<TMatrixD.h>
#include<TColor.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>


using namespace std;


//FileNames
std::string DataFileName;
std::string MCFileName;
std::string MCPropFileName;
std::string OutputFileName;
std::string configFileName;
std::string inputFileName;


//Source params
int Nmass;
double InjIndexAtSource;
double EmaxAtSource;
double RedShiftEvolution;
double FluxNorm;
int SourceEnergyCutoffMode;
double SourceEnergyCutoffShape;
bool EmaxADependent;
std::vector<double> FractAtSource;
std::vector<double> FractParamsAtSource;

//Fit options
bool FitGalSpectrum;
bool FitExtragalFractions;
bool FitGalFractions;
bool FitSourceCutoffShape;
bool CombinedFit;
bool StoreFitResults;
bool FixMassFractions;
bool FixGamma;
bool FixEmax;
bool FixCutoffShape;
bool WeightLikelihood;
double SpectrumLikelihoodWeight;
double XmaxLikelihoodWeight;

void usage() {
  cout<<endl;
	cout<<"*** PROGRAM USAGE ***"<<endl;
	cout<<"CompositionFitter --[runmode] --config=[path-to-configfile] --input=[path-to-inputfile]"<<endl;
	cout<<"  [runmode]"<<endl;
	cout<<"      sourcefit: run source composition fit"<<endl;
	cout<<"      draw: draw fitted results"<<endl;
	cout<<"**********************"<<endl;
} //close usage() 

void DoSourceFit();
void Draw();

int main(int argc, char **argv) {

	if(argc<3){
		cout<<endl;
		cerr<< "ERROR: Incorrect number of arguments...see program usage"<<endl;
		usage();		
		exit(1);
	}
	
	//## handle arguments
	std::string config_arg_format="--config=";
	configFileName="";

	//get path of config macro
	if(!strstr(argv[2],config_arg_format.c_str())){
		//mistake in giving the config argument 
		std::string errMsg = "ERROR: Give the config argument in form --config=[path-to-configmacro]";
    throw std::runtime_error(errMsg);	
		exit(1);
	}
	else{
		//config argument specified correctly => get the file
		std::string config_arg= string(argv[2]);
		int readpos = config_arg.find("="); 
		configFileName = config_arg.substr(readpos+1);// get from "--config=" to the end
		//cout<<"INFO: Read settings from config file "<<configFileName<<endl;
	}

	
	//## parsing info
	CRSourceFitter_ns::ConfigParser* theConfigParser= new CRSourceFitter_ns::ConfigParser(configFileName.c_str());
	theConfigParser->SetVerbosity(1);
	theConfigParser->ReadConfig();

	
	//## handle --input arguments
	if(argv[3]!=NULL){
		std::string input_arg_format="--input=";
		inputFileName="";

		//get path of config macro
		if(!strstr(argv[3],input_arg_format.c_str())){
			//mistake when giving the config argument 
			std::string errMsg = "ERROR: Give the config argument in form --input=[path-to-inputfile] ###";
    	throw std::runtime_error(errMsg);	
			exit(1);
		}
		else{
			//input argument specified correctly => get the file
			std::string input_arg= string(argv[3]);
			int readpos = input_arg.find("="); 
			inputFileName = input_arg.substr(readpos+1);// get from "--input=" to the end
			cout<<"INFO: Drawing fit info from file "<<inputFileName<<endl;
		}
	}//close if
	
	cout<<"get run mode"<<endl;
	//## get run mode
	if(strcmp(argv[1],"--fit")==0) {	
		cout<<"************************************"<<endl;
		cout<<"** SOURCE FITTING                ***"<<endl;
		cout<<"************************************"<<endl;
		DoSourceFit();
	}
	
	else if(strcmp(argv[1],"--draw")==0){
		cout<<"************************************"<<endl;
		cout<<"** DRAWING                       ***"<<endl;
		cout<<"************************************"<<endl;
		Draw();
	}	
	else{
		//mistake in giving the run mode argument
		std::string errMsg = "ERROR: Execution mode " + string(argv[1])
      + " not defined!";
    throw std::runtime_error(errMsg);	 
		exit(1);
	}  
  
  return 0;
  
}//close main



void DoSourceFit(){

  //Create a CompositionFitter class
  CRSourceFitter_ns::SourceFitter theSourceFitter;
	theSourceFitter.SetVerbosity(2);
  theSourceFitter.DoSourceFit();	
  
}//close DoSourceFit()


void Draw(){

	//Create a CompositionFitter class
  CRSourceFitter_ns::DrawFitResults theDrawer;
  
  //Set some class variables
  theDrawer.SetNumberOfComponents(Nmass);
  theDrawer.SetDataFileName(inputFileName.c_str());
	theDrawer.SetOutputFileName("Plot.root");	

	theDrawer.Draw();

}//close Draw()


