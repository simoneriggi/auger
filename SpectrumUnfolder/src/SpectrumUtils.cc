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


int SpectrumUtils::GetModelSpectrum(TH1D& spectrum,TF1* spectrumModel,bool integrateBins){

	if(!spectrumModel) return -1;
	spectrum.Reset();
	if(integrateBins){
		for(int i=0;i<spectrum.GetNbinsX();i++){
			double xmin= spectrum.GetBinLowEdge(i+1);
			double xmax= xmin+spectrum.GetBinWidth(i+1);
			double w= spectrumModel->Integral(xmin,xmax)/(xmax-xmin);
			if(!TMath::IsNaN(w) && fabs(w)!=TMath::Infinity() ){
				spectrum.SetBinContent(i+1,w);
				spectrum.SetBinError(i+1,0);
			}
		}
	}//close if	
	else{
		for(int i=0;i<spectrum.GetNbinsX();i++){
			double x= spectrum.GetBinCenter(i+1);
			double w= spectrumModel->Eval(x);
			if(!TMath::IsNaN(w) && fabs(w)!=TMath::Infinity() ){
				spectrum.SetBinContent(i+1,w);
				spectrum.SetBinError(i+1,0);
			}
		}
	}//close else

	return 0;

}//close GetModelSpectrum()


TF1* SpectrumUtils::ComputeSpectrumModel(SpectrumPars& pars,double xmin,double xmax,int npts){

	//Create model function
	TF1* SpectrumModel= 0;
	int nSpectrumPars= pars.GetNPars();
	int spectrumModel= pars.GetModel();

	if(spectrumModel==ePowerLaw){
		SpectrumModel= new TF1("SpectrumModel",MathUtils::PowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==eFlat){
		SpectrumModel= new TF1("SpectrumModel","[0]",xmin,xmax);
	}
	else if(spectrumModel==eBrokenPowerLaws){
		SpectrumModel= new TF1("SpectrumModel",MathUtils::BrokenPowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		SpectrumModel= new TF1("SpectrumModel",MathUtils::SmoothCutoffPowerLawSpectrum,xmin,xmax,nSpectrumPars);
	}	
	else{
		cerr<<"SpectrumUtils::ComputeSpectrumModel(): INvalid model selected!"<<endl;
		return 0;
	}

	//Set model parameters
	SpectrumModel->SetNpx(npts);
	for(int i=0;i<nSpectrumPars;i++){
		double parValue= pars.GetPar(i)->GetValue();
		SpectrumModel->SetParameter(i,parValue);
	}	
	
	//Compute integral and normalize to 1
	double integral= pars.GetIntegral(xmin,xmax);
	if(spectrumModel==eFlat){
		integral= (xmax-xmin);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		integral= SpectrumModel->Integral(xmin,xmax);
	}
	SpectrumModel->SetParameter(0,1./integral);

	return SpectrumModel;

}//close ComputeSpectrumModel()



TF2* SpectrumUtils::ComputeResponseModel(SpectrumPars& spectrumPars,ResoPars& biasPars,ResoPars& sigmaPars,TriggerPars& triggerPars,double xmin,double xmax,double ymin,double ymax,int npts){

	//Get spectrum pars
	int nSpectrumPars= spectrumPars.GetNPars();
	int spectrumModel= spectrumPars.GetModel();
	FitPars spectrumParList= spectrumPars.GetPars();

	//Get bias pars
	int nBiasPars= biasPars.GetNPars();
	int biasModel= biasPars.GetModel();
	std::vector<double> biasParList= biasPars.GetPars();

	//Get sigma pars
	int nSigmaPars= sigmaPars.GetNPars();
	int sigmaModel= sigmaPars.GetModel();
	std::vector<double> sigmaParList= sigmaPars.GetPars();

	//Get trigger pars
	int nTriggerPars= triggerPars.GetNPars();
	int triggerModel= triggerPars.GetModel();
	std::vector<double> triggerParList= triggerPars.GetPars();

	int nTotPars= (nBiasPars+1) + (nSigmaPars+1) + (nTriggerPars+1) + (nSpectrumPars+1);
	cout<<"SpectrumUtils::ComputeResponseModel(): INFO: nBiasPars="<<nBiasPars<<", nSigmaPars="<<nSigmaPars<<", nTriggerPars="<<nTriggerPars<<", nSpectrumPars="<<nSpectrumPars<<", nTotPars="<<nTotPars<<endl;

	//Init response fcn
	TF2* ResponseModelFcn= new TF2("ResponseModel",MathUtils::ResponseModel,xmin,xmax,ymin,ymax,nTotPars);
	ResponseModelFcn->SetNpx(npts);
	ResponseModelFcn->SetNpy(npts);

	//Set spectrum pars
	cout<<"SpectrumUtils::ComputeResponseModel(): INFO: spectrumModel: "<<spectrumModel<<endl;

	TF1* spectrumModelFcn= ComputeSpectrumModel(spectrumPars,xmin,xmax);
	double integral= spectrumPars.GetIntegral(xmin,xmax);
	if(spectrumModel==eFlat){
		integral= (xmax-xmin);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		integral= spectrumModelFcn->Integral(xmin,xmax);
	}	

	int par_counter= 0;
	ResponseModelFcn->SetParameter(par_counter,spectrumModel);
	par_counter++;
	for(int i=0;i<spectrumParList.GetNPars();i++){
		double parValue= spectrumParList.GetPar(i)->GetValue();
		if(i==0) parValue= 1./integral;
		cout<<"SpectrumUtils::ComputeResponseModel(): INFO: Spectrum Par no. "<<i<<"="<<parValue<<endl;
		ResponseModelFcn->SetParameter(par_counter,parValue);
		par_counter++;
	}

	

	//Set bias pars
	cout<<"SpectrumUtils::ComputeResponseModel(): INFO: BiasModel: "<<biasModel<<endl;

	ResponseModelFcn->SetParameter(par_counter,biasModel);
	par_counter++;
	for(unsigned int i=0;i<biasParList.size();i++){
		cout<<"SpectrumUtils::ComputeResponseModel(): INFO: Bias Par no. "<<i<<"="<<biasParList[i]<<endl;
		ResponseModelFcn->SetParameter(par_counter,biasParList[i]);
		par_counter++;
	}

	//Set sigma pars
	cout<<"SpectrumUtils::ComputeResponseModel(): INFO: SigmaModel="<<sigmaModel<<endl;
	ResponseModelFcn->SetParameter(par_counter,sigmaModel);
	par_counter++;
	for(unsigned int i=0;i<sigmaParList.size();i++){
		cout<<"SpectrumUtils::ComputeResponseModel(): INFO: Sigma Par no. "<<i<<"="<<sigmaParList[i]<<endl;
		
		ResponseModelFcn->SetParameter(par_counter,sigmaParList[i]);
		par_counter++;
	}

	
	//Set trigger pars
	cout<<"SpectrumUtils::ComputeResponseModel(): INFO: TriggerModel="<<triggerModel<<endl;
	ResponseModelFcn->SetParameter(par_counter,triggerModel);
	par_counter++;
	for(unsigned int i=0;i<triggerParList.size();i++){
		cout<<"SpectrumUtils::ComputeResponseModel(): INFO: Trigger Par no. "<<i<<"="<<triggerParList[i]<<endl;
		ResponseModelFcn->SetParameter(par_counter,triggerParList[i]);
		par_counter++;
	}


	//Check given pars
	for(int i=0;i<ResponseModelFcn->GetNpar();i++){
		double parValue= ResponseModelFcn->GetParameter(i);
		cout<<"SpectrumUtils::ComputeResponseModel(): INFO: Response par["<<i<<"]="<<parValue<<endl;
	}

	if(spectrumModelFcn) spectrumModelFcn->Delete();

	return ResponseModelFcn;

}//close ComputeResponseModel()


TH2D* SpectrumUtils::ComputeParametricResponse(SpectrumPars& spectrumPars,ResoPars& biasPars,ResoPars& sigmaPars,TriggerPars& triggerPars,std::vector<double>& TrueBins, std::vector<double>& RecBins){

	//## Check bins
	int NTrueBins= (int)TrueBins.size()-1;
	int NRecBins= (int)RecBins.size()-1;
	if(NTrueBins<=0 || NRecBins<=0){
		cerr<<"SpectrumUtils::ComputeParametricResponse(): ERROR: Invalid number of bins specified!"<<endl;
		return 0;
	}
	double BinEdge_Rec[NRecBins+1];
	double BinEdge_True[NTrueBins+1];

	for(int s=0;s<NRecBins;s++) {
		BinEdge_Rec[s]= RecBins[s];
		cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: BinEdge_Rec["<<s<<"]="<<BinEdge_Rec[s]<<endl;
	}
	BinEdge_Rec[NRecBins]= RecBins[NRecBins];		
	cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: BinEdge_Rec["<<NRecBins<<"]="<<BinEdge_Rec[NRecBins]<<endl;

	for(int s=0;s<NTrueBins;s++) {	
		BinEdge_True[s]= TrueBins[s];
		cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: BinEdge_True["<<s<<"]="<<BinEdge_True[s]<<endl;
	}
	BinEdge_True[NTrueBins]= TrueBins[NTrueBins];
	cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: BinEdge_True["<<NTrueBins<<"]="<<BinEdge_True[NTrueBins]<<endl;

	double LgEmin_true= BinEdge_True[0];
	double LgEmax_true= BinEdge_True[NTrueBins];
	double LgEmin_rec= BinEdge_Rec[0];
	double LgEmax_rec= BinEdge_Rec[NRecBins];


	//## Init spectrum model fcn
	TF1* SpectrumModelFcn= ComputeSpectrumModel(spectrumPars,LgEmin_true,LgEmax_true);
	if(!SpectrumModelFcn){
		cerr<<"SpectrumUtils::ComputeParametricResponse(): ERROR: Failed to compute the spectrum model!"<<endl;
		return 0;
	}

	//## Init response model
	TF2* ResponseModelFcn= ComputeResponseModel(spectrumPars,biasPars,sigmaPars,triggerPars,LgEmin_true,LgEmax_true,LgEmin_rec,LgEmax_rec);
	if(!ResponseModelFcn){
		cerr<<"SpectrumUtils::ComputeParametricResponse(): ERROR: Failed to compute the response model!"<<endl;
		if(SpectrumModelFcn) SpectrumModelFcn->Delete();
		return 0;
	}

	//## Init response matrix
	TH2D* ResponseMat= new TH2D("ResponseMat","ResponseMat",NTrueBins,BinEdge_True,NRecBins,BinEdge_Rec);
	ResponseMat->Sumw2();


	//## Fill matrix
	for(int i=0;i<ResponseMat->GetNbinsX();i++){
		double binWidth_true= ResponseMat->GetXaxis()->GetBinWidth(i+1);	
		double lgEMin_true= ResponseMat->GetXaxis()->GetBinLowEdge(i+1);
		double lgEMax_true= lgEMin_true + binWidth_true;
		double ProbNorm= SpectrumModelFcn->Integral(lgEMin_true,lgEMax_true);
		
		for(int j=0;j<ResponseMat->GetNbinsY();j++){
			double binWidth_rec= ResponseMat->GetYaxis()->GetBinWidth(j+1);	
			double lgEMin_rec= ResponseMat->GetYaxis()->GetBinLowEdge(j+1);
			double lgEMax_rec= lgEMin_rec + binWidth_rec;

			double Rji_noNorm= ResponseModelFcn->Integral(lgEMin_true,lgEMax_true,lgEMin_rec,lgEMax_rec);
			double Rji= Rji_noNorm/ProbNorm;

			cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: Etrue("<<lgEMin_true<<","<<lgEMax_true<<"), Erec("<<lgEMin_rec<<","<<lgEMax_rec<<") ProbNorm="<<ProbNorm<<", Rji_noNorm="<<Rji_noNorm<<" Rji="<<Rji<<endl;
			
			ResponseMat->SetBinContent(i+1,j+1,Rji);
			ResponseMat->SetBinError(i+1,j+1,0.);
			
		}//end loop rec bins
	}//end loop true bins

	
	
	
	/*
	//## Check matrix normalization
	double sumOfRecBins[NTrueBins];//this should be in [0,1] 
  double sumOfTrueBins[NRecBins];//this should sum to 1
  for(int i=0;i<NTrueBins;i++){
		sumOfRecBins[i]= 0;
		for(int j=0;j<NRecBins;j++){
			sumOfRecBins[i]+= ResponseMat->GetBinContent(i+1,j+1);
		}
	}

  for(int i=0;i<NRecBins;i++){
		sumOfTrueBins[i]= 0;
		for(int j=0;j<NTrueBins;j++){
			sumOfTrueBins[i]+= ResponseMat->GetBinContent(j+1,i+1);
		}
	}
	
	for(int i=0;i<NTrueBins;i++){
		cout<<"TRUE BIN "<<i+1<<"  sumRecBins="<<	sumOfRecBins[i]<<endl;
	}
	for(int i=0;i<NRecBins;i++){
		cout<<"REC BIN "<<i+1<<"  sumTrueBins="<<	sumOfTrueBins[i]<<endl;
	}
	*/

	cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: Deleting spectrum model..."<<endl;
	if(SpectrumModelFcn) SpectrumModelFcn->Delete();
	cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: Deleting response model..."<<endl;
	if(ResponseModelFcn) ResponseModelFcn->Delete();	
	cout<<"SpectrumUtils::ComputeParametricResponse(): INFO: done!"<<endl;
	
	return ResponseMat;

}//close SpectrumUtils::ComputeParametricResponse()


}//close namespace 

