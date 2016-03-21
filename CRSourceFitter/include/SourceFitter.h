/**
* @file SourceFitter.h
* @class SourceFitter
* @brief Perform a likelihood fit of spectrum&composition data, determining the source parameters (gamma,Emax,abundances,...) 
* 
* @author S. Riggi
* @date 22/04/2010
*/

#ifndef _SOURCEFITTER_H_
#define _SOURCEFITTER_H_

#include <PropagationMCReader.h>
#include <DataReader.h>

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>

#include <vector>
#include <string>

namespace CRSourceFitter_ns {

class SourceFitter : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/	
		SourceFitter();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceFitter();
		/**
		* \brief Set the verbosity
		*/
		void SetVerbosity(int value) { fVerbosity = value;}
	

		/**
		* \brief Set the number of mass components to be fitted
		*/
		void SetNumberOfComponents(int comp) { fNmass = comp;}
	
		/**
		* \brief Set the input filename containing data info 
		*/
		void SetDataFileName(std::string & fileName) { fDataFileName = fileName;}
		/**
		* \brief Set the input filename containing MC info 
		*/
		void SetMCFileName(std::string & fileName) { fMCFileName = fileName;}
	
		/**
		* \brief Set the input filename containing MC propagation info 
		*/
		void SetPropMCFileName(std::string & fileName) { fPropMCFileName = fileName;}
	
		/**
		* \brief Set the output filename containing the fitted results
		*/
		void SetOutputFileName(std::string & fileName) { fOutputFileName = fileName;}
     
    
  	//####################################
  	//### SET PROP MC OPTIONS     ########
		//####################################
		/**
		* \brief Set the starting value of injection index parameter
		*/
  	void SetInjIndex(double val) {fInjIndex= val;}
		/**
		* \brief Set the starting value of Emax parameter
		*/
  	void SetEmax(double val) {fEmax= val;}
		/**
		* \brief Set the starting value of flux normalization parameter
		*/
		void SetFluxNormalization(double val) {fFluxNormFactor= val;}
		/**
		* \brief Set the starting value of redshift evolution (this is fixed, not included as parameter in the fit)
		*/
  	void SetRedShiftEvolution(double val) {fRedShiftEvolution= val;} 
		/**
		* \brief Set the starting value of mass fractions parameters
		*/         
  	void SetFractionsAtSource(std::vector<double> vect) {fFractionsAtSource= vect;} 
	
		/**
		* \brief Set the starting value of mass fractions parameters. These represent a linear transformation of the mass fractions and are effectively used as fit parameters instead of the mass fractions. This is introduced to guarantee the constraint sum_of_fractions=1 in the fit.
		*/  
		void SetStartFractParamsAtSource(std::vector<double> vect) {fStartFractParamsAtSource= vect;} 
	      
		/**
		* \brief Set the source cutoff mode
		*/       
		void SetSourceCutoffMode(int val) {fSourceCutoffMode= val;}  
		/**
		* \brief Set the source cutoff shape width
		*/       
		void SetSourceCutoffShape(double val) {fSourceCutoffShape= val;}  
		/**
		* \brief Set the A-dependent maximum energy
		*/       
		void SetADependentEmax(bool val) {fEmaxADependent= val;}

    
  	/**
		* \brief Include the source cutoff shape as fit parameter
		*/
  	void FitSourceCutoffShape(bool choice) {fFitSourceCutoffShape= choice;}
	
		/**
		* \brief Include the mass fractions as fit parameters
		*/
		void FitExtragalFractions(bool choice) {fFitExtragalFractions= choice;} 
		/**
		* \brief Include the galactic mass fractions in the fit
		*/ 
		void FitGalFractions(bool choice) {fFitGalFractions= choice;} 
	 
		/**
		* \brief Include the galactic spectrum component in the fit
		*/  
  	void FitGalSpectrum(bool choice) {fFitGalSpectrum= choice;} 
		/**
		* \brief Perform a combined fit of spectrum+Xmax, instead of spectrum alone
		*/  
  	void IsCombinedFit(bool choice) {fIsCombinedFit= choice;} 
	

		/**
		* \brief Fix the shift in the fit?
		*/
		void IsShiftFixed(bool choice) {fFixShift= choice;} 
		/**
		* \brief Fix the fractions in the fit?
		*/
		void IsMassFractionFixed(bool choice) {fFixMassFractions= choice;} 
		/**
		* \brief Fix Emax in the fit?
		*/
		void IsEmaxFixed(bool choice) {fFixEmax= choice;} 
		/**
		* \brief Fix Gamma in the fit?
		*/
		void IsInjIndexFixed(bool choice) {fFixInjIndex= choice;} 
		/**
		* \brief Fix the source cutoff width in the fit?
		*/
		void IsCutoffShapeFixed(bool choice) {fFixSourceCutoffShape= choice;} 
	
		/**
		* \brief Is the likelihood weighted in the fit?
		*/
		void IsLikelihoodWeighted(bool choice) {fWeightLikelihood= choice;} 
		/**
		* \brief Set spectrum likelihood weight
		*/
		void SetSpectrumLikelihoodWeight(double value) {fSpectrumLikelihoodWeight= value;} 
		/**
		* \brief Get Xmax likelihood weight
		*/
		void SetXmaxLikelihoodWeight(double value) {fXmaxLikelihoodWeight= value;}
	
  	/**
		* \brief Execute the source fit
		*/
  	void DoSourceFit();

	private:
  
  	/**
		* \brief Get steering from ConfigParser
		*/
  	void GetConfig();
		/**
		* \brief Init data structures
		*/
  	void Init();
		/**
		* \brief Set the data structures (e.g. histograms,...)
		*/
  	void SetData();

		/**
		* \brief Reset the data structures (e.g. empty histograms,...)
		*/
  	void ResetData();			              
		/**
		* \brief Get the Xmax likelihood
		*/
  	static double GetXmaxLikelihood();
		/**
		* \brief Get the energy spectrum likelihood
		*/
  	static double GetEnergySpectrumLikelihood();
 
		/**
		* \brief Xmax Likelihood definition for MINUIT use
		*/
		static void MaxLikelihoodFcn_Xmax(int& nPar, double* const grad, 
				    		    	double& value, double* const par, 
				    					const int iFlag);
		/**
		* \brief Spectrum Likelihood definition for MINUIT use
		*/
  	static void MaxLikelihoodFcn_EnergySpectrum(int& nPar, double* const grad, 
				    		    	double& value, double* const par, 
				    					const int iFlag);
		/**
		* \brief Combined Likelihood definition for MINUIT use
		*/			    					
  	static void MaxLikelihoodFcnCombined(int& nPar, double* const grad, 
				    		    	double& value, double* const par, 
				    					const int iFlag);				    					
		/**
		* \brief Do the energy spectrum fit
		*/
  	bool FitEnergySpectrum(int verbosity);
		/**
		* \brief Do the Xmax fit
		*/
  	bool FitXmax(int verbosity);
  	/**
		* \brief Do the combined fit
		*/
  	bool CombinedFit(int verbosity);

		/**
		* \brief Get uncertainty of mass fractions 
		*/
		std::vector<double> GetParamUncertainty(std::vector<double> fitparams,TMatrixD CovarianceMatrix,int Npar);
		/**
		* \brief Get element of derivative matrix of mass fractions vs fit params
		*/
		double GetDerivativeMatrixElement(int i,int j,std::vector<double> fitparams);
		/**
		* \brief Store Xmax fit results in output file
		*/
		void StoreXmaxFitResults();
		/**
		* \brief Store spectrum fit results in output file
		*/
		void StoreSpectrumFitResults();
		/**
		* \brief Store elongation rate in output file
		*/
		void StoreElongationRate();
		/**
		* \brief Store RMS in output file
		*/
		void StoreRMS();
		/**
		* \brief Store mass fraction in output file
		*/
		void StoreMassFraction();

	private:

  	static int fNmass;
  
		//## FILENAMES
		std::string fDataFileName;
		std::string fSpectrumTableFileName;
		std::string fMCFileName;
		std::string fPropMCFileName;
		std::string fMatrixFileName;
		std::string fOutputFileName;

		//## FIT VARS
		static DataReader* fDataReader; 
		static PropagationMCReader* fPropMCReader; 
  
  	static double fInjIndex;
  	static double fEmax;
  	static double fRedShiftEvolution;
  	static double fSpectrumLikelihood;
  	static double fXmaxLikelihood;
  	static double fTotLikelihood;
  	static int fFitStatus;
  	static double fEmax_start;
  	static double fInjIndex_start; 
		static double fFluxNormFactor_start;
  	static double fGalacticEmax_start;
  	static double fGalacticInjIndex_start; 
		static double fGalacticFluxNormFactor_start;

		static int fSourceCutoffMode;
		static double fSourceCutoffShape;
		static double fSourceCutoffShapeErr;
		double fSourceCutoffShape_start;
		static bool fEmaxADependent;
		static int fNdof_spectrum;
		static int fNdof_Xmax;
		static double fNevents_spectrum;
		static double fNevents_Xmax;


		//## FIT OPTS
		static bool fFitExtragalFractions;
		static bool fFitGalFractions;	
		static bool fFitGalSpectrum;
		static bool fFitSourceCutoffShape;
	
  	static bool fFixShift;
		static bool fFixMassFractions;
		static bool fFixEmax;
		static bool fFixInjIndex;
		static bool fFixSourceCutoffShape;
		static bool fWeightLikelihood;
		static double fSpectrumLikelihoodWeight;
		static double fXmaxLikelihoodWeight;

		static bool fAddSystematicsInFit;
  
		static bool fIsCombinedFit;	
	
		//## DATA VECTORS
		static std::vector<TH1D*> fXmaxData;
		static std::vector<TH1D*> fEnergyData;		

		std::vector<double> fMeanXmaxData;
  	std::vector<double> fMeanXmaxErrData;
  	std::vector<double> fRMSData;
  	std::vector<double> fRMSErrData;
  	std::vector<double> fMeanEnergyData;
  	std::vector<double> fMeanEnergyErrData;

		static TH1D* fAugerEnergySpectrum;
		static TH1D* fAugerEnergySpectrumE3;
		static TH1D* fAugerEnergySpectrumEvents;
  	static std::vector<double> fExposure;
  
  	static std::vector<double> fFractionsAtSource;
		static std::vector<double> fFractionsAtSourceLowErr;
		static std::vector<double> fFractionsAtSourceHighErr;
  	static std::vector<double> fStartFractParamsAtSource;
  	static std::vector< std::vector<double> > fFractionsAtEarth;
  	static std::vector< std::vector<double> > fFractionsAtEarth_spectrumfit;
  	static std::vector< std::vector<double> > fFractionsAtEarth_Xmaxfit;
		static std::vector< std::vector<double> > fFractionsAtEarthHighErr_Xmaxfit;
		static std::vector< std::vector<double> > fFractionsAtEarthLowErr_Xmaxfit;
		static std::vector< std::vector<double> > fDirectFractionsAtEarth_spectrumfit;
		static std::vector< std::vector<double> > fFractionsAtDetector_Xmaxfit;
  	static std::vector< std::vector<double> > fFractionsAtDetectorHighErr_Xmaxfit;
  	static std::vector< std::vector<double> > fFractionsAtDetectorLowErr_Xmaxfit;

	  struct MassFractionStruct{
	    double* fSourceFraction;
	  };

	  MassFractionStruct theMassFractStruct;
	  MassFractionStruct theStartMassFractStruct;
  
  	static std::vector< std::vector<TH1D*> > fXmaxMC;
		static std::vector< std::vector<TH1D*> > fXmaxMC_scaled; 
		static std::vector< std::vector<TH1D*> > fXmaxGenMC;	
		static std::vector< std::vector<TH1D*> > fEnergyMC;
		static std::vector< std::vector<TH1D*> > fGenEnergyMC;
	  
  	static TH1D* fTotPropEnergySpectrum_spectrumfit;
  	static TH1D* fTotPropEnergySpectrumE3_spectrumfit;
  	static TH1D* fGalacticSpectrum_spectrumfit;
		static std::vector<TH1D*> fPropEnergySpectrum_spectrumfit;
		static std::vector<TH1D*> fPropDirectEnergySpectrum_spectrumfit;
  	static std::vector<TH1D*> fPropEnergySpectrumE3_spectrumfit;
  	static std::vector<TH1D*> fPropDirectEnergySpectrumE3_spectrumfit;
  	static std::vector<TH1D*> fPropEnergySpectrum_Xmaxfit;
  	static std::vector<TH1D*> fPropEnergySpectrumE3_Xmaxfit;
  
		std::vector< std::vector<double> > fMeanXmaxMC;
  	std::vector< std::vector<double> > fMeanXmaxErrMC;
  	std::vector< std::vector<double> > fRMSMC;
  	std::vector< std::vector<double> > fRMSErrMC;
  	std::vector< std::vector<double> > fMeanEnergyMC;
  	std::vector< std::vector<double> > fMeanEnergyErrMC;
  
  	std::vector< std::vector<double> > fMeanXmaxConexMC;
  	std::vector< std::vector<double> > fMeanXmaxErrConexMC;
  	std::vector< std::vector<double> > fXmaxRMSConexMC;
  	std::vector< std::vector<double> > fXmaxRMSErrConexMC;

		std::vector<TH1D*> fXmaxFit;
  	std::vector<TH1D*> fXmaxTrueFit;
  	std::vector<double> fMeanXmaxFit;
  	std::vector<double> fMeanXmaxErrFit;
  	std::vector<double> fRMSFit;
  	std::vector<double> fRMSErrFit;
  	std::vector<double> fRMSTrueFit;
  	std::vector<double> fRMSTrueErrFit;
		//ER vectors
		std::vector<double> fERRaw;
		std::vector<double> fERRawErr;
		std::vector<double> fERTrue;
		std::vector<double> fERTrueErr;

		static double fInjIndex_Xmax;
	  static double fInjIndexErr_Xmax;
	  
	  static double fInjIndex_spectrum;
	  static double fInjIndexErr_spectrum;
	
	  static double fEmax_Xmax;
	  static double fEmaxErr_Xmax;
  
	  static double fEmax_spectrum;
	  static double fEmaxErr_spectrum;

	  static double fFluxNormFactor_Xmax;
	  static double fFluxNormFactorErr_Xmax;
  
	  static double fFluxNormFactor_spectrum;
	  static double fFluxNormFactorErr_spectrum;

	  static double fFluxNormFactor;
	  static double fFluxNormFactorErr;

		std::vector<TMatrixD> fCovMatrix;
		std::vector<double> fitParameterList;
		std::vector<double> fitParameterErrorList;

	  //### GALACTIC PARAMETERS
	  static double fGalacticInjIndex_spectrum;
	  static double fGalacticInjIndexErr_spectrum;
	  static double fGalacticEmax_spectrum;
	  static double fGalacticEmaxErr_spectrum;
	  static double fGalacticFluxNormFactor_spectrum;
	  static double fGalacticFluxNormFactorErr_spectrum;
	
		//## OUTPUT ROOT TREE VARS
		TFile* fOutputFile;
		static TTree* fXmaxMinimizationInfo;
		static TTree* fSpectrumMinimizationInfo;
	
		static const int NmaxComp= 10;
		int fNmassFree;
		double fSourceFract[NmaxComp];
		double fSourceFractLowErr[NmaxComp];
		double fSourceFractHighErr[NmaxComp];
		double fStartFractFitParam[NmaxComp-1];	
		double fFractFitParam[NmaxComp-1];
		double fFractFitParamErr[NmaxComp-1];
		double fXmaxShift;
		double fXmaxShiftErr;

		static int fVerbosity;

	ClassDef(SourceFitter,1)
	
};

}//close namespace

#endif
