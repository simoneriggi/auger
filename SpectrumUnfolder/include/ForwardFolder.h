/**
* @file ForwardFolder.h
* @class ForwardFolder
* @brief Forward folding of the energy spectrum
* 
* @author S. Riggi
* @date 20/09/2011
*/

#ifndef _FORWARD_FOLDER_H_
#define _FORWARD_FOLDER_H_

#include <MathUtils.h>

#include <TObject.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TMatrixD.h>
#include <TTree.h>

namespace Unfolder_ns {

class ForwardFolder : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    ForwardFolder();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~ForwardFolder();

		
	public:
		int RunUnfold(TH1D* spectrumHisto,TH2D* responseMatrix,SpectrumPars& initFitPars,bool computeUncertainties=true,int nRandomSamples=100,bool useFitRange=false,double fitMin=-1,double fitMax=-1);

		TH1D* GetUnfoldedSpectrum(){return fUnfoldedSpectrum;}
		TH1D* GetTrueSpectrum(){return fTrueSpectrum;}
		TH1D* GetRecSpectrum(){return fRecSpectrum;}
		TH1D* GetForwardFoldedSpectrum(){return fForwardFoldedSpectrum;}

	private:
		int Init();
		int Clear();
		
		TH1D* UnfoldSpectrum(SpectrumPars& initFitPars,std::string runMode="FIT");
		static void LogLikelihoodFcn(int& nPar, double* gin, double &f, double* par, int iflag);		
		int ComputeUncertainties(SpectrumPars& initFitPars,int nRandomSamples=100);

	private:
		
		//Data
		static TF1* fTrueSpectrumModelFcn;
		TH1D* fRecSpectrum;
		static TH1D* fCurrentRecSpectrum;		
		static TH1D* fTrueSpectrum;
		static TH1D* fForwardFoldedSpectrum;
		static TH1D* fCurrentForwardFoldedSpectrum;
		TH1D* fUnfoldedSpectrum;
		static TH2D* fResponseMatrix;
		int fNTrueBins;
		double fLgEMin_true;
		double fLgEMax_true;
		std::vector<double> fLgEBins_true;

		int fNRecBins;
		double fLgEMin_rec;
		double fLgEMax_rec;
		std::vector<double> fLgEBins_rec;
		

		//Options
		static bool fUseFitRange;
		static double fLgEMin_fit;
		static double fLgEMax_fit;
		
		//Fit results
		int fFitStatus;
		static const int NMAXFITTEDPARS= 10;
		TMatrixD* fCovarianceMatrix;	
		double fFitPar[NMAXFITTEDPARS];
		double fFitParErr[NMAXFITTEDPARS];

		static SpectrumPars* fInitFitPars;

		ClassDef(ForwardFolder,1)

};//close class


#ifdef __MAKECINT__
#pragma link C++ class ForwardFolder+; 
#endif

}//close namespace


#endif

