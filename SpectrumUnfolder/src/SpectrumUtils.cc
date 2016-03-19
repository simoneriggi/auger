/**
* @file SpectrumUtils.cc
* @class SpectrumUtils
* @brief Spectrum utilities
* 
* @author S. Riggi
* @date 20/09/2011
*/


#include <SpectrumUtils.h>
#include <MathUtils.h>

#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>


#include <Math/GSLRndmEngines.h>
#include <TDecompChol.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <cmath>


#include <vector>

using namespace std;
using namespace ROOT::Math;


ClassImp(Unfolder_ns::SpectrumUtils)

namespace Unfolder_ns {


SpectrumUtils::SpectrumUtils(){
	
}//close constructor


SpectrumUtils::~SpectrumUtils(){


}//close destructor


TH1D* SpectrumUtils::GetFoldedSpectrum(TH1D* trueSpectrum,TH2D* responseMatrix){

	//## Forward-folding of the input histogram with the response matrix

	//Check input args
	if(!trueSpectrum || !responseMatrix){
		cerr<<"SpectrumUtil::GetFoldedSpectrum(): ERROR: Null ptr to given input spectrum and/or response matrix!"<<endl;	
		return 0;
	}

	int nTrueBins= responseMatrix->GetNbinsX();
	int nRecBins= responseMatrix->GetNbinsY();

	if(trueSpectrum->GetNbinsX() != nTrueBins){
		cerr<<"SpectrumUtil::GetFoldedSpectrum(): ERROR: Number of bins of passed histogram differs from that of the response matrix!"<<endl;	
		return 0;
	}

	//Check consistency between the true input spectrum and response matrix (xaxis)
	if(!CheckSpectrumBinnings(trueSpectrum->GetXaxis(),responseMatrix->GetXaxis())){
		cerr<<"SpectrumUtil::GetFoldedSpectrum(): ERROR: Bin mismatch between response matrix and true spectrum to be folded!"<<endl;		
		return 0;
	}

	//Create the folded spectrum
	bool isVariableBin= responseMatrix->GetYaxis()->IsVariableBinSize();
	TH1D* foldSpectrum= 0;
	if(isVariableBin){
		const double* ybins= responseMatrix->GetYaxis()->GetXbins()->GetArray();
		foldSpectrum= new TH1D("foldSpectrum","foldSpectrum",nRecBins,ybins);
	}
	else{
		foldSpectrum= new TH1D("foldSpectrum","foldSpectrum",nRecBins,responseMatrix->GetYaxis()->GetXmin(),responseMatrix->GetYaxis()->GetXmax());
	}
	foldSpectrum->Reset();

	//Fill the folded spectrum
	for(int i=0;i<nRecBins;i++){
		double nFolded= 0.;		

		for(int j=0;j<nTrueBins;j++){
			double binCenterX= responseMatrix->GetXaxis()->GetBinCenter(j+1);
			double w= responseMatrix->GetBinContent(j+1,i+1);

			//Find corresponding bin in true spectrum 
			int binId= trueSpectrum->FindBin(binCenterX);
			double nTrue= trueSpectrum->GetBinContent(binId);
			
			nFolded+= w*nTrue;
		}//end loop rec bins

		foldSpectrum->SetBinContent(i+1,nFolded);
		foldSpectrum->SetBinError(i+1,0.);
	}//end loop true bins

	return foldSpectrum;

}//close GetFoldedSpectrum()


TH1D* SpectrumUtils::GetFluctuatedSpectrum(TH1D* inputSpectrum){

	//## Takes an histogram as input and returns an histogram with entries
	//## fluctuated according to a multinomial distribution 
	if(!inputSpectrum){
		cerr<<"SpectrumUtils::GetFluctuatedSpectrum(): ERROR: Null pointer to input histogram...exit!"<<endl;	
		return 0;
	}
	
	if(inputSpectrum->GetEntries()<=0){
		cerr<<"SpectrumUtils::GetFluctuatedSpectrum(): WARNING: Input histogram is empty...returning same histo!"<<endl;	
		return ((TH1D*)inputSpectrum->Clone("fluctuatedSpectrum"));
	}

	//### INIT RANDOM GENERATOR
	//## Generate a suitable seed
	unsigned long int seed1= time(0);
	unsigned long int seed2= clock();
	timespec timeStruct;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timeStruct);
	
	unsigned long int seed3= timeStruct.tv_nsec;
	
	ROOT::Math::GSLRandomEngine randomEngine;
	randomEngine.Initialize(); 
	randomEngine.SetSeed(seed1+seed2+seed3);

	//## Compute vector with probabilities
	//## e.g. given the input histogram, regarded as the expected counts, p= n_i/n_tot 
	std::vector<double> multinomialProbability;
	multinomialProbability.clear();
	multinomialProbability.assign(inputSpectrum->GetNbinsX(),0.);//init

	std::vector<double> poissonNEvents;
	poissonNEvents.clear();
	poissonNEvents.assign(inputSpectrum->GetNbinsX(),0.);//init
	double ntot= inputSpectrum->Integral();
	for(int i=0;i<inputSpectrum->GetNbinsX();i++){
		double n= inputSpectrum->GetBinContent(i+1);
		multinomialProbability[i]= n/ntot;

		poissonNEvents[i]= randomEngine.Poisson(n);
	}

	std::vector<unsigned int> fluctuatedEntries = randomEngine.Multinomial(inputSpectrum->GetEntries(), multinomialProbability);
	
	TH1D* fluctuatedSpectrum= (TH1D*)inputSpectrum->Clone("fluctuatedSpectrum");
	fluctuatedSpectrum->Reset();
	
	for(int i=0;i<fluctuatedSpectrum->GetNbinsX();i++){
		fluctuatedSpectrum->SetBinContent(i+1,fluctuatedEntries[i]);
		fluctuatedSpectrum->SetBinError(i+1,sqrt(fluctuatedEntries[i]));
	}
	
	return fluctuatedSpectrum;
	
}//close GetFluctuatedSpectrum()


int SpectrumUtils::CopySpectrumContent(TH1D* h1,TH1D* h2){

	if(!h1){
		cerr<<"SpectrumUtils::CopySpectrumContent(): ERROR: Passing null histogram...exit!"<<endl;
	  return -1;
	}

	h2->Reset();
	for(int i=0;i<h1->GetNbinsX();i++){
		h2->SetBinContent(i+1,h1->GetBinContent(i+1));
		h2->SetBinError(i+1,h1->GetBinError(i+1));
	}
	
	return 0;

}//close SpectrumUtils::CopyHistoContent()


bool SpectrumUtils::CheckSpectrumBinnings(const TAxis* a1, const TAxis* a2){

	//Axis 1 shall be contained inside the Axis 2 AND binning shall match
	int Nbins_1= a1->GetNbins();
	int Nbins_2= a2->GetNbins();
	double xmin_1= a1->GetXmin();
	double xmax_1= a1->GetXmax();
	double xmin_2= a2->GetXmin();
	double xmax_2= a2->GetXmax();

	if(xmin_1<xmin_2 || xmax_1>xmax_2){
		cout<<"SpectrumUtils::CheckSpectrumBinnings(): WARN: First histo axis is not contained inside the second histo!"<<endl;
		return false;
	}

	//Check binnings
	bool isFailed= false;
	for(int i=0;i<Nbins_1;i++){
		double binCenter_1= a1->GetBinCenter(i+1);
		double binLowEdge_1= a1->GetBinLowEdge(i+1);
		double binUpEdge_1= binLowEdge_1 + a1->GetBinWidth(i+1);
		
		//Find bin id in axis 2
		int binId_2= a2->FindBin(binCenter_1);
		double binCenter_2= a2->GetBinCenter(binId_2);
		double binLowEdge_2= a2->GetBinLowEdge(binId_2);
		double binUpEdge_2= binLowEdge_2+a2->GetBinWidth(binId_2);
		
		if( !TMath::AreEqualRel(binLowEdge_1,binLowEdge_2,1E-10) || !TMath::AreEqualRel(binUpEdge_1,binUpEdge_2,1E-10)){
			isFailed= true;
			break;
		}
	}//end loop internal histo bins

	if(isFailed){
		cout<<"SpectrumUtils::CheckSpectrumBinnings(): WARN: Bin mismatch between the two histograms!"<<endl;
		return false;
	}

	return true;

}//close CheckSpectrumBinnings()


int SpectrumUtils::GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel){

	if(!spectrumModel) return -1;
	spectrum.Reset();
	for(int i=0;i<spectrum.GetNbinsX();i++){
		double x= spectrum.GetBinCenter(i+1);
		double xmin= spectrum.GetBinLowEdge(i+1);
		double xmax= xmin+spectrum.GetBinWidth(i+1);
		double w= spectrumModel->Integral(xmin,xmax)/(xmax-xmin);
		//double w= spectrumModel->Eval(x);
		spectrum.SetBinContent(i+1,w);
		spectrum.SetBinError(i+1,0);
	}
	
	return 0;

}//close GetModelSpectrum()





}//close namespace 

