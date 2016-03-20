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

#include <MathUtils.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TF2.h>


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
		static int GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel,bool integrateBins=true);		
		static bool CheckSpectrumBinnings(const TAxis* a1, const TAxis* a2);
		static TF1* ComputeSpectrumModel(SpectrumPars& pars,double xmin,double xmax,int npts=1000);
		static TF2* ComputeResponseModel(SpectrumPars& spectrumPars,ResoPars& biasPars,ResoPars& sigmaPars,TriggerPars& triggerPars,double xmin,double xmax,double ymin,double ymax,int npts=1000);
		static TH2D* ComputeParametricResponse(SpectrumPars& spectrumPars,ResoPars& biasPars,ResoPars& sigmaPars,TriggerPars& triggerPars,std::vector<double>& TrueBins, std::vector<double>& RecBins);

	ClassDef(SpectrumUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class SpectrumUtils+; 
#endif


}//close namespace 

#endif
