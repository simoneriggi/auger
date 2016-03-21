/**
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/
#ifndef _CONFIG_PARSER_H_
#define _CONFIG_PARSER_H_

#include <TObject.h>

#include <vector>
#include <string>

namespace CRSourceFitter_ns {

class ConfigParser : public TObject {
  
	public:
  
		/** 
		\brief Class constructor
 		*/
  	ConfigParser(const char* filename);

		/** 
		\brief Class destructor
 		*/
  	virtual ~ConfigParser();

	public:
		
		enum SourceCutoffModel { eExpo=1, eFermi= 2, eKink= 3, eKink2= 4};
		
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		void ReadConfig();
		/** 
		\brief Set the printing verbosity
 		*/
  	void SetVerbosity(int value){verbosity= value;}
	
	private:

		/** 
		\brief Check integrity of parsed information
 		*/
		bool Check();
		/** 
		\brief Print parsed information
 		*/
		void Dump();

	public:	

		//## FileNames
		static std::string fDataFileName;
		static std::string fSpectrumTableFileName;
		static std::string fMCFileName;
		static std::string fMatrixFileName;
		static std::string fOutputFileName;


		//## Source vars
		static double fInjIndex;
		static double fEmax;
		static double fRedshiftEvolution;
		static double fFluxNorm;
		static std::vector<double> fFract;	
	
		static bool fEmaxADependent;
		static double fSourceDistanceCutoff;
		static int fSourceCutoffMode;	
		static double fSourceCutoffShape;
	
	
		//## Fit options
		static bool fFitGalSpectrum;
		static bool fFitExtragalFractions;
		static bool fFitGalFractions;
		static bool fFitXmaxSyst;
	
		static bool fStoreFitResults;

		static bool fCombinedFit;
		static bool fFitSourceCutoffShape;
		static bool fFixMassFractions;	
		static bool fFixInjIndex;
		static bool fFixEmax;
		static bool fFixShift;
		static bool fFixSourceCutoffShape;
		static bool fWeightLikelihood;
		static double fSpectrumLikelihoodWeight;
		static double fXmaxLikelihoodWeight;
	
	private:

		const char* configFile;
		int verbosity;

	ClassDef(ConfigParser,1)
	
};

}//close namespace

#endif
