/**
* @file SpectrumUtils.h
* @class SpectrumUtils
* @brief Spectrum utilities
* 
* @author S. Riggi
* @date 20/09/2011
*/

#ifndef _SPECTRUM_UTILS_H_
#define _SPECTRUM_UTILS_H_

#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>


namespace Unfolder_ns {

//===================================================
//==        SPECTRUM UTILS
//===================================================
class SpectrumUtils : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		SpectrumUtils();
		/** 
		\brief Class destructor
 		*/
		virtual ~SpectrumUtils();
		

	public:
		static TH1D* GetFluctuatedSpectrum(TH1D*);
		static TH1D* GetFoldedSpectrum(TH1D* trueSpectrum,TH2D* responseMatrix);
		static int CopySpectrumContent(TH1D* h1,TH1D* h2);
		static int GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel);		
		static bool CheckSpectrumBinnings(const TAxis* a1, const TAxis* a2);
		
	ClassDef(SpectrumUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class SpectrumUtils+; 
#endif


}//close namespace 

#endif
