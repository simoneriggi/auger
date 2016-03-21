/**
* @file ConfigParser.cc
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#include "ConfigParser.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

ClassImp(CRSourceFitter_ns::ConfigParser)

namespace CRSourceFitter_ns {

//## FileNames
std::string ConfigParser::fDataFileName;
std::string ConfigParser::fMCFileName;
std::string ConfigParser::fSpectrumTableFileName;
std::string ConfigParser::fMatrixFileName;
std::string ConfigParser::fOutputFileName;

//## Source vars
double ConfigParser::fInjIndex;
double ConfigParser::fEmax;
double ConfigParser::fFluxNorm;
double ConfigParser::fRedshiftEvolution;
std::vector<double> ConfigParser::fFract;	
	
bool ConfigParser::fEmaxADependent;
int ConfigParser::fSourceCutoffMode;	
double ConfigParser::fSourceCutoffShape;

	
//## Fit options
bool ConfigParser::fFitGalSpectrum;
bool ConfigParser::fFitExtragalFractions;
bool ConfigParser::fFitGalFractions;
bool ConfigParser::fFitXmaxSyst;
bool ConfigParser::fCombinedFit;
bool ConfigParser::fFitSourceCutoffShape;
bool ConfigParser::fFixMassFractions;	
bool ConfigParser::fFixInjIndex;
bool ConfigParser::fFixEmax;
bool ConfigParser::fFixShift;
bool ConfigParser::fFixSourceCutoffShape;
bool ConfigParser::fWeightLikelihood;
double ConfigParser::fSpectrumLikelihoodWeight;
double ConfigParser::fXmaxLikelihoodWeight;
bool ConfigParser::fStoreFitResults;	

using namespace std;


ConfigParser::ConfigParser(const char* filename){

	//## init params
	verbosity= 1;
	configFile= filename;

	//## FileNames
	fDataFileName= "";
	fSpectrumTableFileName= "";
	fMCFileName= "";
	fMatrixFileName= "";
	fOutputFileName= "Output.root";


	//## Source vars
	fInjIndex= 2.0;
	fEmax= 20.0;
	fFluxNorm= 10.0;
		
	fEmaxADependent= false;
	fSourceCutoffMode= eExpo;	
	fSourceCutoffShape= 1;
	
	//## Fit options
	fFitGalSpectrum= false;
	fFitExtragalFractions= true;
	fFitGalFractions= false;
	fFitXmaxSyst= false;
	
	fCombinedFit= true;
	fFitSourceCutoffShape= true;
	fFixMassFractions= false;	
	fFixInjIndex= false;
	fFixEmax= false;
	fFixShift= false;
	fFixSourceCutoffShape= false;
	fWeightLikelihood= false;
	fSpectrumLikelihoodWeight= 1;
	fXmaxLikelihoodWeight= 1;



}//close costructor

ConfigParser::~ConfigParser(){

}//close destructor


void ConfigParser::ReadConfig() {

	// read configuration from file
	cout<<"ConfigParser::ReadConfig(): Parsing config file "<<configFile<<endl;

  ifstream in;  
  in.open(configFile);
  if(!in.good()) {
    string errMsg = "ConfigParser::ReadConfig(): Error reading config file " + string(configFile)
      + " \n  ********* no configuration settings!!! ********* ";
    throw std::runtime_error(errMsg);
		exit(1);
  }

	//Start parsing the config file
	char buffer[1000];//container for a full line
  std::string descriptor;//container for config descriptor 
  in.getline(buffer,1000);//get the full line
	std::string parsedline;

	
	while(std::getline(in,parsedline)) {
		char first_char= *(parsedline.c_str());
		
		if(first_char!='#' && first_char!='\n' && first_char!=' '){
			stringstream line(parsedline);
			stringstream line_copy(parsedline);
      
			line_copy >> descriptor;
      
			if(verbosity>1) cout<<"descriptor="<<descriptor<<endl;

      if(descriptor!="\n" && descriptor!=""){
				//get all config parameters

				//########################
				//####  DATA FILES
				//########################
				if(descriptor.compare("DataFile")==0){
					std::string thisFileName;
		    	line >> descriptor >> thisFileName;			
					fDataFileName= thisFileName;
		  	}//close else if
				else if(descriptor.compare("MCFile")==0){
					std::string thisFileName;
		    	line >> descriptor >> thisFileName;
					fMCFileName= thisFileName;			
				}//close else if
				else if(descriptor.compare("SpectrumTableFile")==0){
					std::string thisFileName;
		    	line >> descriptor >> thisFileName;			
					fSpectrumTableFileName= thisFileName;
		  	}//close else if						
				else if(descriptor.compare("MatrixFile")==0){
					std::string thisFileName;
		    	line >> descriptor >> thisFileName;			
					fMatrixFileName= thisFileName;
		  	}//close else if
				else if(descriptor.compare("OutputFile")==0){
					std::string thisFileName;
		    	line >> descriptor >> thisFileName;			
					fOutputFileName= thisFileName;
		  	}//close else if
				
				//########################
				//####  FIT PARAMETERS
				//########################
				else if(descriptor.compare("InjIndex")==0){
		    	line >> descriptor >> fInjIndex;			
		  	}//close else if
				else if(descriptor.compare("Emax")==0){
		    	line >> descriptor >> fEmax;			
		  	}//close else if
				else if(descriptor.compare("FluxNormalization")==0){
		    	line >> descriptor >> fFluxNorm;			
		  	}//close else if
				else if(descriptor.compare("Fraction")==0){
					double thisEntry;
		    	line >> descriptor >> thisEntry;
					fFract.push_back(thisEntry);	
				}//close else if
				else if(descriptor.compare("SourceCutoffMode")==0){
		    	line >> descriptor >> fSourceCutoffMode;			
		  	}//close else if
				else if(descriptor.compare("SourceCutoffShape")==0){
		    	line >> descriptor >> fSourceCutoffShape;			
		  	}//close else if			
				else if(descriptor.compare("EmaxADependent")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fEmaxADependent= true;
					else if(thisFlagValue.compare("F")==0) fEmaxADependent= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				
				//########################
				//####  FIT OPTIONS
				//########################
				else if(descriptor.compare("FitGalSpectrum")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFitGalSpectrum= true;
					else if(thisFlagValue.compare("F")==0) fFitGalSpectrum= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FitExtragalFractions")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFitExtragalFractions= true;
					else if(thisFlagValue.compare("F")==0) fFitExtragalFractions= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FitGalFractions")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fFitGalFractions= true;
					else if(thisFlagValue.compare("F")==0) fFitGalFractions= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("CombinedFit")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;			
					if(thisFlagValue.compare("T")==0) fCombinedFit= true;
					else if(thisFlagValue.compare("F")==0) fCombinedFit= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FitSourceCutoffShape")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;	
					if(thisFlagValue.compare("T")==0) fFitSourceCutoffShape= true;
					else if(thisFlagValue.compare("F")==0) fFitSourceCutoffShape= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FixMassFractions")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;	
					if(thisFlagValue.compare("T")==0) fFixMassFractions= true;
					else if(thisFlagValue.compare("F")==0) fFixMassFractions= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FixGamma")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;	
					if(thisFlagValue.compare("T")==0) fFixInjIndex= true;
					else if(thisFlagValue.compare("F")==0) fFixInjIndex= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FixEmax")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;	
					if(thisFlagValue.compare("T")==0) fFixEmax= true;
					else if(thisFlagValue.compare("F")==0) fFixEmax= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("FixSourceCutoffShape")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue;		
					if(thisFlagValue.compare("T")==0) fFixSourceCutoffShape= true;
					else if(thisFlagValue.compare("F")==0) fFixSourceCutoffShape= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if
				else if(descriptor.compare("WeightLikelihood")==0){
					std::string thisFlagValue;
		    	line >> descriptor >> thisFlagValue >> fSpectrumLikelihoodWeight >> fXmaxLikelihoodWeight;	
					if(thisFlagValue.compare("T")==0) fWeightLikelihood= true;
					else if(thisFlagValue.compare("F")==0) fWeightLikelihood= false;	
					else{
						string errMsg = "ConfigParser::ReadConfig(): invalid setting, use T or F";
    				throw std::runtime_error(errMsg);
						exit(1);
					}
		  	}//close else if

				else{
					//config setting not defined
					line >> descriptor;
					string errMsg = "ConfigParser::ReadConfig(): Descriptor " + descriptor
      	                + " not defined \n  ********* bad settings!!! ********* ";
    			throw std::runtime_error(errMsg);
					exit(1);
				}//close else

			}//close if descriptor
		}//close if buffer

		if (!in.good()) break;
	}//close while

	in.close();

	//now check the integrity of parsed info
	bool status= Check();
	
	//eventually print parsed info
	if(verbosity>0){
		Dump();
	}

}//close ConfigParser::ReadConfig()


bool ConfigParser::Check(){

	//perform a check of integrity of parsed information from the config file
	//e.g. check vector size, throw errors or warnings
	if(fDataFileName.c_str()==""){
		string errMsg = "ConfigParser::Check(): Empty data filename...exit!";
    throw std::runtime_error(errMsg);
		exit(1);
	}
	if(fMCFileName.c_str()==""){
		string errMsg = "ConfigParser::Check(): Empty MC filename...exit!";
    throw std::runtime_error(errMsg);
		exit(1);
	}
	if(fMatrixFileName.c_str()==""){
		string errMsg = "ConfigParser::Check(): Empty matrix filename...exit!";
    throw std::runtime_error(errMsg);
		exit(1);
	}
	if(fSpectrumTableFileName.c_str()==""){
		string errMsg = "ConfigParser::Check(): Empty spectrum table filename...exit!";
    throw std::runtime_error(errMsg);
		exit(1);
	}

	
	return true;
	
}//close ConfigParser::Check()


void ConfigParser::Dump(){

	//print parsed information
	cout<<"#############################"<<endl;
	cout<<"###     RUN SETTINGS      ###"<<endl;
	cout<<"#############################"<<endl;
	cout<<"****** DATA FILE ******"<<endl;
	cout<<"Data file: "<<fDataFileName<<endl;
	cout<<"MC file: "<<fMCFileName<<endl;
	cout<<"Spectrum table file: "<<fSpectrumTableFileName<<endl;
	cout<<"Matrix file: "<<fMatrixFileName<<endl;
	cout<<"Output file: "<<fOutputFileName<<endl;

	cout<<"****** FIT OPTIONS ******"<<endl;
	cout<<"InjIndex: "<<fInjIndex<<endl;
	cout<<"Emax: "<<fEmax<<endl;
	cout<<"FluxNorm: "<<fFluxNorm<<endl;
	for(unsigned int k=0;k<fFract.size();k++) cout<<"Fract["<<k<<"]: "<<fFract[k]<<endl;
	
	cout<<"SourceCutoffMode: "<<fSourceCutoffMode<<endl;
	cout<<"SourceCutoffShape: "<<fSourceCutoffShape<<endl;
	cout<<"EmaxADependent? "<<fEmaxADependent<<endl;
	cout<<"CombinedFit? "<<fCombinedFit<<endl;
	
	cout<<"FitGalSpectrum? "<<fFitGalSpectrum<<endl;
	cout<<"FitExtragalFractions? "<<fFitExtragalFractions<<endl;
	cout<<"FitGalFractions? "<<fFitGalFractions<<endl;
	cout<<"FitXmaxSyst? "<<fFitXmaxSyst<<endl;
	cout<<"FitCutoffShape? "<<fFitSourceCutoffShape<<endl;
	
	cout<<"FixMassFractions? "<<fFixMassFractions<<endl;
	cout<<"FixGamma? "<<fFixInjIndex<<endl;
	cout<<"FixEmax? "<<fFixEmax<<endl;
	cout<<"FixCutoffShape? "<<fFixSourceCutoffShape<<endl;
	cout<<"WeightLikelihood? "<<fWeightLikelihood<<"  WSpectrum="<<fXmaxLikelihoodWeight<<"  WXmax="<<fXmaxLikelihoodWeight<<endl;


}//close ConfigParser::Dump()

}//close namespace
