/**
* @file DrawFitResults.h
* @class DrawFitResults
* @brief Draw the results of the fit (fitted data and MC, elongation rate, rms, chi square, ...)
* 
* @author S. Riggi
* @date 22/04/2010
*/

#ifndef _DRAWFITRESULTS_H_
#define _DRAWFITRESULTS_H_

#include <vector>
#include <string>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TStyle.h>
#include <TApplication.h>

namespace CRSourceFitter_ns {

class DrawFitResults : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		DrawFitResults();
		/** 
		\brief Class destructor
 		*/
		virtual ~DrawFitResults();

	public:
		/** 
		\brief Draw all
 		*/
		void Draw();

		/**
		* \brief Set the number of mass components to be drawn
		*/
		void SetNumberOfComponents(int comp) { fNmass = comp;}
		
		/**
		* \brief Set the input filename containing data info to be drawn
		*/
		void SetDataFileName(const char * fileName) { fDataFileName = fileName;}
		/**
		* \brief Set the output filename containing plots of the fit results
		*/
		void SetOutputFileName(const char * fileName) { fOutputFileName = fileName;}

		/**
		* \brief Set the Xmax MC shift
		*/
  	void SetXmaxShift(double val) {fXmaxShift= val;}
	

	private:
  
  
		/**
		* \brief Inititialize data structures
		*/
		void Init();
		/**
		* \brief Draw the Xmax fit results (distributions, chisquare,...)
		*/
		void DrawXmaxFit();
		/**
		* \brief Draw the fitted mass fractions
		*/
		void DrawMassFractions();
		/**
		* \brief Draw the spectrum fit results
		*/
		void DrawSpectrumFit();
		/**
		* \brief Draw the fitted elongation rate
		*/
		void DrawElongationRate();
		/**
		* \brief Draw the fitted Xmax RMS
		*/
		void DrawRMS();

		/**
		* \brief Function to draw the ER reference given by the models 
		*/
		TF1* defineElong(double* p, int color, int style);
		/**
		* \brief Function to draw the RMS reference given by the models
		*/
		TF1* defineRMSref(double* p, int color, int style);

	
	private:

  	int fNmass; 
		const char * fDataFileName;
		const char * fOutputFileName;

		TFile* fInputFile;
		TFile* fPlotFile;
		TApplication* fApplication;
		TStyle* myStyle;
		double fXmaxShift;

	ClassDef(DrawFitResults,1)
	
};

}//close namespace

#endif

