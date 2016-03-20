/**
* @file ForwardFolder.cc
* @class ForwardFolder
* @brief Forward folding of the energy spectrum
* 
* @author S. Riggi
* @date 20/09/2011
*/

#include <ForwardFolder.h>
#include <MathUtils.h>
#include <SpectrumUtils.h>

#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TDecompChol.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Unfolder_ns::ForwardFolder)

namespace Unfolder_ns {


bool ForwardFolder::fUseFitRange;
double ForwardFolder::fLgEMin_fit;
double ForwardFolder::fLgEMax_fit;
SpectrumPars* ForwardFolder::fInitFitPars;

TF1* ForwardFolder::fTrueSpectrumModelFcn;
TH1D* ForwardFolder::fCurrentRecSpectrum;		
TH1D* ForwardFolder::fTrueSpectrum;
TH1D* ForwardFolder::fForwardFoldedSpectrum;
TH1D* ForwardFolder::fCurrentForwardFoldedSpectrum;
TH2D* ForwardFolder::fResponseMatrix;

ForwardFolder::ForwardFolder(){

	fTrueSpectrumModelFcn= 0;
	fUnfoldedSpectrum= 0;
	fTrueSpectrum= 0;
	fForwardFoldedSpectrum= 0;
	fCurrentForwardFoldedSpectrum= 0;
	fCurrentRecSpectrum= 0;
	fCovarianceMatrix= 0;

}//close costructor


ForwardFolder::~ForwardFolder(){

	Clear();

}//close destructor

int ForwardFolder::Clear(){
		
	if(fTrueSpectrumModelFcn) fTrueSpectrumModelFcn->Delete();
	if(fUnfoldedSpectrum) fUnfoldedSpectrum->Delete();
	if(fTrueSpectrum) fTrueSpectrum->Delete();
	if(fForwardFoldedSpectrum) fForwardFoldedSpectrum->Delete();
	if(fCurrentForwardFoldedSpectrum) fCurrentForwardFoldedSpectrum->Delete();
	if(fCurrentRecSpectrum) fCurrentRecSpectrum->Delete();
	if(fCovarianceMatrix) fCovarianceMatrix->Delete();

	return 0;

}//close Clear()

int ForwardFolder::Init(){

	//## Clear existing data
	if(Clear()<0) {
		cerr<<"ForwardFolder::Init(): ERROR: Failed to clear data!"<<endl;
		return -1;
	}

	//## Set spectrum & matrix bins
	if(!fRecSpectrum || !fResponseMatrix){
		return -1;
	}

	fNTrueBins= fResponseMatrix->GetNbinsX();
	fNRecBins= fResponseMatrix->GetNbinsY();
	
	fLgEBins_true.assign(fNTrueBins+1,0);
	fLgEBins_rec.assign(fNRecBins+1,0);
	
	double BinEdge_RecSpectrum[fNRecBins+1];
	double BinEdge_TrueSpectrum[fNTrueBins+1];

	for(int s=0;s<fNTrueBins;s++) {
		fLgEBins_true[s]= fResponseMatrix->GetXaxis()->GetBinLowEdge(s+1);
		BinEdge_TrueSpectrum[s]= fLgEBins_true[s];
	}
	fLgEBins_true[fNTrueBins]= fLgEBins_true[fNTrueBins-1]+fResponseMatrix->GetXaxis()->GetBinWidth(fNTrueBins);
	BinEdge_TrueSpectrum[fNTrueBins]= fLgEBins_true[fNTrueBins];

	for(int s=0;s<fNRecBins;s++) {
		fLgEBins_rec[s]= fResponseMatrix->GetYaxis()->GetBinLowEdge(s+1);
		BinEdge_RecSpectrum[s]= fLgEBins_rec[s];
	}
	fLgEBins_rec[fNRecBins]= fLgEBins_rec[fNRecBins-1]+fResponseMatrix->GetXaxis()->GetBinWidth(fNRecBins);
	BinEdge_RecSpectrum[fNRecBins]= fLgEBins_rec[fNRecBins];


	//## Ensure rec spectrum and response matrix have no bin offsets
	bool hasSameBinLimits= SpectrumUtils::CheckSpectrumBinnings(fRecSpectrum->GetXaxis(),fResponseMatrix->GetXaxis());
	if(!hasSameBinLimits){
		cerr<<"ForwardFolder::Init(): ERROR: Response matrix and rec spectrum have different bin limits!"<<endl;
		return -1;
	}
	

	//Spectrum model
	int npts= 1000;
	if(fTrueSpectrumModelFcn) fTrueSpectrumModelFcn->Delete();
	int spectrumModel= fInitFitPars->GetModel();
	int nSpectrumPars= fInitFitPars->GetNPars();
	if(spectrumModel==ePowerLaw){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",MathUtils::PowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eBrokenPowerLaws){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",MathUtils::BrokenPowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eSmoothBrokenPowerLaws){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn",MathUtils::SmoothCutoffPowerLawSpectrum,fLgEMin_true,fLgEMax_true,nSpectrumPars);
	}
	else if(spectrumModel==eFlat){
		fTrueSpectrumModelFcn= new TF1("TrueSpectrumModelFcn","[0]",fLgEMin_true,fLgEMax_true);
	}
	else{
		cerr<<"ForwardFolder::Init(): ERROR: INvalid spectrum model selected!"<<endl;
		return -1;
	}

	for(int i=0;i<nSpectrumPars;i++){
		double parValue= fInitFitPars->GetPar(i)->GetValue();
		fTrueSpectrumModelFcn->SetParameter(i,parValue);	
	}

	fTrueSpectrumModelFcn->SetNpx(npts);
	
	if(fUnfoldedSpectrum) fUnfoldedSpectrum->Delete();
	fUnfoldedSpectrum= new TH1D("UnfoldedSpectrum","UnfoldedSpectrum",fNTrueBins,BinEdge_TrueSpectrum);
	fUnfoldedSpectrum->Sumw2();
	fUnfoldedSpectrum->SetMarkerStyle(24);
	fUnfoldedSpectrum->SetMarkerSize(1.1);
	fUnfoldedSpectrum->SetMarkerColor(kBlack);
	fUnfoldedSpectrum->SetLineColor(kBlack);
	
	if(fTrueSpectrum) fTrueSpectrum->Delete();
	fTrueSpectrum= new TH1D("TrueSpectrum","TrueSpectrum",fNTrueBins,BinEdge_TrueSpectrum);
	fTrueSpectrum->Sumw2();
	fTrueSpectrum->SetMarkerStyle(8);
	fTrueSpectrum->SetMarkerSize(1.1);
	fTrueSpectrum->SetMarkerColor(kBlack);
	fTrueSpectrum->SetLineColor(kBlack);

	if(fForwardFoldedSpectrum) fForwardFoldedSpectrum->Delete();
	fForwardFoldedSpectrum= new TH1D("ForwardFoldedSpectrum","ForwardFoldedSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fForwardFoldedSpectrum->Sumw2();
	fForwardFoldedSpectrum->SetMarkerStyle(8);
	fForwardFoldedSpectrum->SetMarkerSize(1.1);
	fForwardFoldedSpectrum->SetMarkerColor(kBlack);
	fForwardFoldedSpectrum->SetLineColor(kBlack);
	
	if(fCurrentForwardFoldedSpectrum) fCurrentForwardFoldedSpectrum->Delete();
	fCurrentForwardFoldedSpectrum= new TH1D("CurrentForwardFoldedSpectrum","CurrentForwardFoldedSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fCurrentForwardFoldedSpectrum->Sumw2();
	fCurrentForwardFoldedSpectrum->SetMarkerStyle(8);
	fCurrentForwardFoldedSpectrum->SetMarkerSize(1.1);
	fCurrentForwardFoldedSpectrum->SetMarkerColor(kBlack);
	fCurrentForwardFoldedSpectrum->SetLineColor(kBlack);
	
	if(fCurrentRecSpectrum) fCurrentRecSpectrum->Delete();
	fCurrentRecSpectrum= new TH1D("CurrentRecSpectrum","CurrentRecSpectrum",fNRecBins,BinEdge_RecSpectrum);
	fCurrentRecSpectrum->Sumw2();
	fCurrentRecSpectrum->SetMarkerStyle(8);
	fCurrentRecSpectrum->SetMarkerSize(1.1);
	fCurrentRecSpectrum->SetMarkerColor(kBlack);
	fCurrentRecSpectrum->SetLineColor(kBlack);

	fCurrentRecSpectrum->Reset();
	for(int i=0;i<fRecSpectrum->GetNbinsX();i++) {
		fCurrentRecSpectrum->SetBinContent(i+1,fRecSpectrum->GetBinContent(i+1));
		fCurrentRecSpectrum->SetBinError(i+1,fRecSpectrum->GetBinError(i+1));
	}

	//Fit data
	if(fCovarianceMatrix) fCovarianceMatrix->Delete();
	fCovarianceMatrix= new TMatrixD(nSpectrumPars,nSpectrumPars);

	//Init random numbers
	delete gRandom;
	gRandom= new TRandom3;

	return 0;

}//close Init()

int ForwardFolder::RunUnfold(TH1D* recSpectrum,TH2D* responseMatrix,SpectrumPars& initFitPars,bool computeUncertainties,int nRandomSamples,bool useFitRange,double fitMin,double fitMax){

	//## Check input spectrum/matrix
	if(!recSpectrum || !responseMatrix){
		cerr<<"ForwardFolder::RunUnfold(): ERROR: Null ptr to input spectrum given!"<<endl;
		return -1;
	}
	fRecSpectrum= recSpectrum;
	fResponseMatrix= responseMatrix;
	fInitFitPars= initFitPars.Clone();

	fUseFitRange= useFitRange;
	fLgEMin_fit= fRecSpectrum->GetXaxis()->GetXmin();
	fLgEMax_fit= fRecSpectrum->GetXaxis()->GetXmax();
	if(useFitRange){
		fLgEMin_fit= fitMin;
		fLgEMax_fit= fitMax;
	}

	//## Init
	if(Init()<0){
		cerr<<"ForwardFolder::RunUnfold(): ERROR: Failed to init!"<<endl;
		return -1;
	}

	//## Unfold spectrum
	fUnfoldedSpectrum= UnfoldSpectrum(initFitPars,"FIT");

	//## Copy 
	SpectrumUtils::CopySpectrumContent(fCurrentForwardFoldedSpectrum,fForwardFoldedSpectrum);

	//## Calculate syst uncertainties
	if(computeUncertainties) {
		cout<<"ForwardFolder::Run(): INFO: Calculating systematics on the unfolded flux..."<<endl;
		ComputeUncertainties(initFitPars,nRandomSamples);
	}
	return 0;

}//close Unfold()


TH1D* ForwardFolder::UnfoldSpectrum(SpectrumPars& initFitPars,std::string runMode){

	//#######################
	//##    INIT FIT
	//#######################	
	double arglist[10];
  int ierflag = 0;
	double amin,edm,errdef;
  int nvpar,nparx,icstat;
	int nPar= initFitPars.GetNPars();	
	
	TMinuit* gMinuit = new TMinuit(nPar);
  gMinuit->SetPrintLevel(0);
	gMinuit->SetMaxIterations(100000);
	gMinuit->SetFCN(ForwardFolder::LogLikelihoodFcn);
	gMinuit->mnexcm("SET NOW",arglist,0,ierflag);

	arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflag);

	//## Set starting values and step sizes for parameters
	for(int i=0;i<nPar;i++){
		FitPar* thisFitPar= initFitPars.GetPar(i);
		std::string parName= thisFitPar->GetName();
		double parValue= thisFitPar->GetValue();
		double stepValueRel= thisFitPar->GetStepSize();
		double stepValue= fabs(stepValueRel*parValue);
		bool isLimited= thisFitPar->IsLimited();
		double parMinValue= thisFitPar->GetMinValue();
		double parMaxValue= thisFitPar->GetMaxValue();
		double isFixed= thisFitPar->IsFixed();
		if(isLimited) gMinuit->mnparm(i,parName.c_str(), parValue, stepValue, parMinValue, parMaxValue, ierflag);
		else gMinuit->mnparm(i,parName.c_str(), parValue, stepValue, 0, 0, ierflag);

		if(isFixed) gMinuit->FixParameter(i);
	}

	
	//## SET MINUIT STRATEGY
  // 0 ==> low level minimization but small number of FCN calls
  // 1 ==> intermediate
  // 2 ==> max level but many FCN calls
  arglist[0]= 1;
  gMinuit->mnexcm("SET STR",arglist,1,ierflag);

	arglist[0]= 0.5;//0.5 likelihood, 1 ChiSquare
  gMinuit->mnexcm("SET ERR",arglist,1,ierflag);

	if(runMode=="CALL"){ //## Call Chi2 with initial values
		gMinuit->mnexcm("CALL FCN", arglist ,0,ierflag);
	}
	else{ //## Fit
  	arglist[0] = 100000;
  	arglist[1] = 0.1;
  	gMinuit->mnexcm("MINIMIZE", arglist ,2,ierflag);
	}

	//## Print results
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  
	fFitStatus= gMinuit->GetStatus();

	//## Get fitted parameters
	double p,ep;
	double FinalPars[nPar];
	double FinalParsErr[nPar];

	cout<<"*** Fit ***"<<endl;
	for(int i=0;i<nPar;i++){
		gMinuit->GetParameter(i,p,ep);	
		FinalPars[i]= p;
		FinalParsErr[i]= ep;
		fFitPar[i]= p;
		fFitParErr[i]= ep;
		cout<<"FinalPar"<<i+1<<"="<<FinalPars[i]<<"  +- "<<FinalParsErr[i]<<endl;
	}
	

	//## Get covariance matrix
	double covMatrix[nPar][nPar];		
  gMinuit->mnemat(&covMatrix[0][0],nPar);
  	
  cout<<"*** COV MATRIX ***"<<endl;  	
	TMatrixD CovarianceMatrix(nPar,nPar);
		
	for(int i=0;i<nPar;i++){
		for(int j=0;j<nPar;j++){
		 	CovarianceMatrix(i,j)= covMatrix[i][j];		
			(*fCovarianceMatrix)(i,j)= covMatrix[i][j];	
		  if(j==nPar-1) cout<<covMatrix[i][j]<<endl;
		  else cout<<covMatrix[i][j]<<"  ";
		}//close for j
	}//close for i
	cout<<endl;

	//## Correction factor for rec flux
	//## Create unfolded histo and fill with results
	TH1D* aUnfoldedSpectrumHisto= (TH1D*)fUnfoldedSpectrum->Clone(); 
	aUnfoldedSpectrumHisto->Reset();
	fTrueSpectrumModelFcn->SetParameters(FinalPars);

	double nTot= fCurrentRecSpectrum->Integral();
	double nTrueTot= fTrueSpectrum->Integral();

	for(int j=0;j<fNRecBins;j++){
		double lgERec= fCurrentRecSpectrum->GetBinCenter(j+1);
		double nRec= fCurrentRecSpectrum->GetBinContent(j+1);
		double nRecErr= fCurrentRecSpectrum->GetBinError(j+1);
		int recBinId= fTrueSpectrum->FindBin(lgERec);
				
		double nTrue= fTrueSpectrum->GetBinContent(recBinId);
		double nTrueErr= fTrueSpectrum->GetBinError(recBinId);
		fTrueSpectrum->SetBinContent(recBinId,nTrue);
		fTrueSpectrum->SetBinError(recBinId,nTrueErr);

		//double nFF= fForwardFoldedSpectrum->GetBinContent(j+1);
		double nFF= fCurrentForwardFoldedSpectrum->GetBinContent(j+1);
		
		double correctionFactor= 0;
		if(nFF!=0) correctionFactor= (double)(nTrue)/(double)(nFF);
		double nUnfold= nRec*correctionFactor;
		double nUnfoldErr= nRecErr*correctionFactor;
		aUnfoldedSpectrumHisto->SetBinContent(recBinId,nUnfold);
		aUnfoldedSpectrumHisto->SetBinError(recBinId,nUnfoldErr);

		cout<<"ForwardFolder::UnfoldSpectrum(): lgERec="<<lgERec<<"  nRec="<<nRec<<"  nTrue="<<nTrue<<"  nFF="<<nFF<<"  corr="<<correctionFactor<<"  nUnfold="<<nUnfold<<"  Unfold/nRec-1="<<nUnfold/nRec-1<<endl;
	}//end loop rec bins

	return aUnfoldedSpectrumHisto;

}//close UnfoldSpectrum()


void ForwardFolder::LogLikelihoodFcn(int& nPar, double* gin, double &f, double* par, int iflag){

	//## The parameters are the correction factors to be applied to the measured spectrum 
	//## Set correction in current true spectrum
	f = 0.;
	double LogLikelihood= 0.;
	double Chi2= 0.;
	double NDF= 0.;
	int nTotPar= fInitFitPars->GetNPars();//total number of parameters (nPar is the free number of pars)

	//## Update "true" spectrum given current fit parameters
	fTrueSpectrumModelFcn->SetParameters(par);
	fTrueSpectrumModelFcn->SetParameter(0,pow(10,par[0]));
	
	SpectrumUtils::GetModelSpectrum(*fTrueSpectrum,fTrueSpectrumModelFcn,false);

	//## Forward fold the "true" spectrum with resolution model
	TH1D* FoldedSpectrum= SpectrumUtils::GetFoldedSpectrum(fTrueSpectrum,fResponseMatrix);


	//fForwardFoldedSpectrum->Reset();
	fCurrentForwardFoldedSpectrum->Reset();

	for(int i=0;i<fCurrentRecSpectrum->GetNbinsX();i++){
		double lgERec= fCurrentRecSpectrum->GetBinCenter(i+1);
		double lgERecMin= fCurrentRecSpectrum->GetBinLowEdge(i+1);
		double lgERecMax= lgERecMin + fCurrentRecSpectrum->GetBinWidth(i+1);
		
		double nData= fCurrentRecSpectrum->GetBinContent(i+1);

		//Find corresponding bin in folded spectrum
		int binId= FoldedSpectrum->FindBin(lgERec);
		double nExp= FoldedSpectrum->GetBinContent(binId);
		
		//fForwardFoldedSpectrum->SetBinContent(i+1,nExp);
		//fForwardFoldedSpectrum->SetBinError(i+1,sqrt(nExp));
		fCurrentForwardFoldedSpectrum->SetBinContent(i+1,nExp);
		fCurrentForwardFoldedSpectrum->SetBinError(i+1,sqrt(nExp));

		//Skip bins if requested
		if(fUseFitRange && (lgERecMin<fLgEMin_fit || lgERecMax>fLgEMax_fit) ) continue;

		double thisChi2= 0;
		if(nExp>0){
  		LogLikelihood+= -(nData*log(nExp)-nExp);
			if(nData>0){
				Chi2+= nExp-nData + nData*log(nData/nExp);
				thisChi2= nExp-nData + nData*log(nData/nExp);
			}
			else{
				Chi2+= nExp;
				thisChi2= nExp;
			}
			//cout<<"INFO: lgERec="<<lgERec<<"  nData="<<nData<<"  nExp="<<nExp<<"  LL="<<LogLikelihood<<"  Chi2="<<Chi2<<"  DeltaChi2="<<thisChi2<<endl;
		}
	
		NDF++;

	}//end loop rec bins

	
	f = LogLikelihood;

	Chi2*= 2;
	NDF-= nPar-2;

	/*
	cout<<"*** CURRENT FIT ***"<<endl;
	cout<<"nPar="<<nPar<<" nTotPar="<<nTotPar<<endl;
	cout<<"LL="<<LogLikelihood<<endl;
	cout<<"Chi2="<<Chi2<<"  Ndf="<<NDF<<endl;
	cout<<"Chi2/Ndf="<<Chi2/NDF<<endl;
	cout<<"== CURRENT PARS =="<<endl;
	for(int i=0;i<nTotPar;i++){
		cout<<"Par["<<i<<"]="<<par[i]<<endl;
	}
	cout<<"*******************"<<endl;
	cout<<endl;
	*/
	FoldedSpectrum->Delete();

}//close ForwardFolder::LogLikelihoodFcn()


int ForwardFolder::ComputeUncertainties(SpectrumPars& initFitPars,int nRandomSamples){

	//## Check covariance matrix
	if(!fCovarianceMatrix){
		cerr<<"ForwardFolder::ComputeUncertainties(): ERROR: Null pointer to covariance matrix, cannot compute uncertainties!"<<endl;
		return -1;
	}
	if(nRandomSamples<=0){
		cerr<<"ForwardFolder::ComputeUncertainties(): ERROR: Invalid number of random sampled given!"<<endl;
		return -1;
	}

	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Covariance matrix:"<<endl;
	fCovarianceMatrix->Print();

	int nPar= fCovarianceMatrix->GetNrows();
	int nFreePar= 0;
	std::vector<bool> isFreeFlags;
	std::vector<int> freeParIndex;
	for(int i=0;i<initFitPars.GetNPars();i++){
		FitPar* thisFitPar= initFitPars.GetPar(i);
		bool isFixed= thisFitPar->IsFixed();
		if(!isFixed) {
			nFreePar++;
			isFreeFlags.push_back(true);
		}	
		else{
			isFreeFlags.push_back(false);
			freeParIndex.push_back(i);
		}
	}

	TMatrixD C(nFreePar,nFreePar);
	double FittedPar[nPar];
	int ix= 0;
	int iy= 0;

	for(int i=0;i<nPar;i++){
		FittedPar[i]= fFitPar[i];

		for(int j=0;j<nPar;j++){
			if(isFreeFlags[i] && isFreeFlags[j]) {
				C(ix,iy)= (*fCovarianceMatrix)(i,j);
				if(j==nPar-1) iy=0;
				else iy++;
			}
		}//end loop y
		if(isFreeFlags[i]) ix++;
	}//end loop x

	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Covariance matrix for free pars (N="<<nFreePar<<"):"<<endl;
	C.Print();


	/*
	TMatrixD C(nPar,nPar);
	double FittedPar[nPar];
	for(int i=0;i<nPar;i++){
		FittedPar[i]= fFitPar[i];
		for(int j=0;j<nPar;j++){
			C(i,j)= (*fCovarianceMatrix)(i,j);
		}
	}
	*/

	//## Find a Cholesky decomposition of the covariance matrix for error propagation
	cout<<"|C|="<<C.Determinant()<<endl;
	
	TDecompChol CCholDecomp(C);
  CCholDecomp.Decompose();
	cout<<"--> Choleski Decomposition..."<<endl;	
	CCholDecomp.Print();
	
	//TMatrixD* CCholDecompTriang= new TMatrixD(nPar,nPar);
	TMatrixD* CCholDecompTriang= new TMatrixD(nFreePar,nFreePar);
	CCholDecompTriang= (TMatrixD*)(&CCholDecomp.GetU());
	cout<<"--> Choleski Decomposition Triangular..."<<endl;	
	CCholDecompTriang->Print();

	//TMatrixD* CCholDecompTriangTransp= new TMatrixD(nPar,nPar);
	TMatrixD* CCholDecompTriangTransp= new TMatrixD(nFreePar,nFreePar);
	CCholDecompTriangTransp->Transpose(*CCholDecompTriang);


	//################################
	//##  STATISTICAL UNCERTAINTIES 
	//################################
	double StatSigma[fNTrueBins];
	std::vector< std::vector<double> > nUnfold_Stat;
	//double nUnfold_Stat[nRandomSamples][fNTrueBins];
	double nUnfoldMean_Stat[fNTrueBins];
	int nToys_Stat= 0;
	for(int i=0;i<fNTrueBins;i++) {
		StatSigma[i]= 0.;
		nUnfoldMean_Stat[i]= 0.;
		nUnfold_Stat.push_back( std::vector<double>() );
	}

	
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Propagate rec spectrum statistical uncertainties into the unfolded flux ..."<<endl;

	for(int k=0;k<nRandomSamples;k++){
		cout<<"== RAND SPECTRUM NO. "<<k+1<<" =="<<endl;	

		//Get a fluctuated version of the rec spectrum
		TH1D* thisFluctuatedHisto= SpectrumUtils::GetFluctuatedSpectrum(fRecSpectrum);
		if(!thisFluctuatedHisto) continue;

		//Copy content to tmp rec spectrum (used in the unfolding fit)
		fCurrentRecSpectrum->Reset();
		if(SpectrumUtils::CopySpectrumContent(thisFluctuatedHisto,fCurrentRecSpectrum)<0) continue;
			
		//Run the unfolding over the fluctuated spectrum
		TH1D* FluctUnfoldedSpectrum_Stat= UnfoldSpectrum(initFitPars,"FIT");
		if(fFitStatus!=0) {
			cout<<"ForwardFolder::ComputeUncertainties(): WARN: Unfolding run not converged for rndom sample no. "<<k+1<<", skip it!"<<endl;
			continue;
		}

		//Update unfolded stats counts fr uncertainty calculation
		cout<<"ForwardFolder::ComputeUncertainties(): INFO: Updating mean counts for forward-folded random spectra ..."<<endl;	
		nToys_Stat++;

		for(int i=0;i<fNTrueBins;i++) {
			double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
			double nUnfoldToyMC_Stat= FluctUnfoldedSpectrum_Stat->GetBinContent(i+1);
			if(TMath::IsNaN(nUnfoldToyMC_Stat)){
				cerr<<"ForwardFolder::ComputeUncertainties(): WARN: Unfolded bin "<<i+1<<" for random spectrum sample no. "<<k<<" is invalid (NAN)!"<<endl;
			}
			nUnfoldMean_Stat[i]+= nUnfoldToyMC_Stat;
			//nUnfold_Stat[k][i]= nUnfoldToyMC_Stat;
			nUnfold_Stat[i].push_back(nUnfoldToyMC_Stat);
		}//end loop true bins
				
		if(FluctUnfoldedSpectrum_Stat) FluctUnfoldedSpectrum_Stat->Delete();
	}//end loop toys
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: #"<<nToys_Stat<<" unfolding over the rec spectrum performed..."<<endl;
	
	if(nToys_Stat<=1){
		cerr<<"ForwardFolder::ComputeUncertainties(): INFO: No random unfolding succeeded, cannot compute uncertainty!"<<endl;
		if(CCholDecompTriang) CCholDecompTriang->Delete();
		if(CCholDecompTriangTransp) CCholDecompTriangTransp->Delete();
		return -1;
	}
	

	//## Calculate sigma syst true model
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Computing sigma stat ..."<<endl;
	for(int i=0;i<fNTrueBins;i++) {
		double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
		nUnfoldMean_Stat[i]/= (double)(nToys_Stat);
			
		//for(int k=0;k<nRandomSamples;k++){
		for(int k=0;k<nUnfold_Stat[i].size();k++){
			//StatSigma[i]+= pow(nUnfold_Stat[k][i]-nUnfoldMean_Stat[i],2);	
			StatSigma[i]+= pow(nUnfold_Stat[i][k]-nUnfoldMean_Stat[i],2);
		}
		StatSigma[i]/= (double)(nToys_Stat-1);
		StatSigma[i]= sqrt(StatSigma[i]);
		
		double relBias_Stat= 0.;
		double relSyst_Stat= 0.;
		if(nUnfold>0) {
			relBias_Stat= nUnfoldMean_Stat[i]/nUnfold-1;	
		}
		if(nUnfoldMean_Stat[i]>0){
			relSyst_Stat= StatSigma[i]/nUnfoldMean_Stat[i];
		}
		cout<<"ForwardFolder::ComputeUncertainties(): INFO: bin "<<i<<" StatErr="<<StatSigma[i]<<"  Bias="<<relBias_Stat<<"  RelStatErr="<<relSyst_Stat<<endl;
	}//end loop true bins
		



	//##############################
	//##   FIT PARS SYSTEMATICS 
	//##############################
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Propagate true model systematics into the unfolded flux ..."<<endl;
	
	double SystSigmaTrueModel[fNTrueBins];
	double nUnfold_SystTrueModel[nRandomSamples][fNTrueBins];
	double nUnfoldMean_SystTrueModel[fNTrueBins];
	for(int i=0;i<fNTrueBins;i++) {
		SystSigmaTrueModel[i]= 0.;
		nUnfoldMean_SystTrueModel[i]= 0.;
	}

	for(int k=0;k<nRandomSamples;k++){
		cout<<"== RANDOM SPECTRUM NO. "<<k+1<<" =="<<endl;		

		//## Generate new true model around uncertainties		
		/*
		TMatrixD RandModelPar(nPar,1);
		for(int i=0;i<nPar;i++){
			RandModelPar(i,0)= gRandom->Gaus(0,1);
		}
		TMatrixD FluctModelPar(nPar,1);
		FluctModelPar= (*CCholDecompTriangTransp)*RandModelPar;

		double newModelPar[nPar];
		for(int i=0;i<nPar;i++){
			newModelPar[i]= FittedPar[i] + FluctModelPar(i,0);
		}
		*/

		TMatrixD RandModelPar(nFreePar,1);
		for(int i=0;i<nFreePar;i++){
			RandModelPar(i,0)= gRandom->Gaus(0,1);
		}
		TMatrixD FluctModelPar(nFreePar,1);
		FluctModelPar= (*CCholDecompTriangTransp)*RandModelPar;
		

		double newModelPar[nPar];
		for(int i=0;i<nPar;i++) newModelPar[i]= FittedPar[i];
		for(unsigned int i=0;i<freeParIndex.size();i++){
			int parIndex= freeParIndex[i];	
			newModelPar[parIndex]+= FluctModelPar(i,0);
		}
		

		//Update pars
		SpectrumPars* modelFitPars= initFitPars.Clone();
		for(int i=0;i<nPar;i++){
			modelFitPars->SetParValue(i,newModelPar[i]);
		}

		TH1D* FluctUnfoldedSpectrum= UnfoldSpectrum(*modelFitPars,"CALL");

		cout<<"ForwardFolder::ComputeUncertainties(): INFO: Calculate mean counts for forward-folded toy spectra..."<<endl;	
		for(int i=0;i<fNTrueBins;i++) {
			double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
			double nUnfoldToyMC= FluctUnfoldedSpectrum->GetBinContent(i+1);
			//cout<<"bin "<<i<<" nUnfold="<<nUnfold<<"  nUnfoldToyMC="<<nUnfoldToyMC<<endl;
			nUnfoldMean_SystTrueModel[i]+= nUnfoldToyMC;
			nUnfold_SystTrueModel[k][i]= nUnfoldToyMC;
			
		}//end loop true bins

		if(modelFitPars){
			delete modelFitPars;
			modelFitPars= 0;
		}

	}//end loop toys


	//## Calculate sigma syst true model
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Calculate sigma syst for true model..."<<endl;
	for(int i=0;i<fNTrueBins;i++) {
		double nUnfold= fUnfoldedSpectrum->GetBinContent(i+1);
		nUnfoldMean_SystTrueModel[i]/= (double)(nRandomSamples);
		for(int k=0;k<nRandomSamples;k++){
			SystSigmaTrueModel[i]+= pow(nUnfold_SystTrueModel[k][i]-nUnfoldMean_SystTrueModel[i],2); 	
		}
		SystSigmaTrueModel[i]/= (double)(nRandomSamples-1);
		SystSigmaTrueModel[i]= sqrt(SystSigmaTrueModel[i]);
		double relBias= 0.;
		double relSyst= 0.;
		if(nUnfold>0) {
			relBias= nUnfoldMean_SystTrueModel[i]/nUnfold-1;
		}
		if(nUnfoldMean_SystTrueModel[i]>0){
			relSyst= SystSigmaTrueModel[i]/nUnfoldMean_SystTrueModel[i];
		}
		cout<<"ForwardFolder::ComputeUncertainties(): INFO: bin "<<i<<" Syst(TrueModel)="<<SystSigmaTrueModel[i]<<"  Bias(TrueModel)="<<relBias<<"  RelSyst(TrueModel)="<<relSyst<<endl;
	}
	
	
	//## Forward-folding results
	cout<<"ForwardFolder::ComputeUncertainties(): INFO: Final unfolding results..."<<endl;
	
	for(int i=0;i<fUnfoldedSpectrum->GetNbinsX();i++){
		double binCenter= fUnfoldedSpectrum->GetBinCenter(i+1);
		double binContent= fUnfoldedSpectrum->GetBinContent(i+1);
		int recBinId= fRecSpectrum->FindBin(binCenter);
		double nRec= fRecSpectrum->GetBinContent(recBinId);			
		double nRecErr= fRecSpectrum->GetBinError(recBinId);			
		
		double totErrorSqr= 0;

		//Add stat uncertainties
		double statError_rec= fRecSpectrum->GetBinError(recBinId);
		double statError= StatSigma[i];
		totErrorSqr+= statError*statError;
	
		//Add syst uncertainty due to unfolding fit
		double systError_UnfoldingAlgo= SystSigmaTrueModel[i];
		totErrorSqr+= systError_UnfoldingAlgo*systError_UnfoldingAlgo;

		/*
		double systError_Exposure= binContent*fExposureSystUncertainty;
		double systError_ECal= SystSigmaECal[i];
		double systError_EnergyFD= SystSigmaEnergyFD[i];
		double systError_EnergySDReso= SystSigmaEnergySDReso[i];
		double systError_UnfoldingMatrix= SystSigmaMatrix[i];
		double systError_UnfoldingProtonMatrix= SystSigmaProtonMatrix[i];
		double systError_UnfoldingIronMatrix= SystSigmaIronMatrix[i];

		double totError= sqrt(statError*statError + systError_ECal*systError_ECal + systError_UnfoldingMatrix*systError_UnfoldingMatrix + systError_UnfoldingAlgo*systError_UnfoldingAlgo);
		*/

		double totError= sqrt(totErrorSqr);
		double binError= statError;	
			
		fUnfoldedSpectrum->SetBinContent(i+1,binContent);
		fUnfoldedSpectrum->SetBinError(i+1,binError);	
		cout<<"bin "<<i+1<<"  nRec="<<nRec<<"  +- "<<nRecErr<<"("<<nRecErr/nRec<<")  nUnfold="<<binContent<<"  +- "<<binError<<"("<<binError/binContent<<")  (recStatErr="<<statError_rec<<"), syst(unfoldAlgo)="<<systError_UnfoldingAlgo<<", totError="<<totError<<endl;
	
		/*
		fLogE= fUnfoldedSpectrum->GetBinCenter(i+1);
		double E= pow(10,fLogE);
		double binWidth= fUnfoldedSpectrum->GetBinWidth(i+1);
		fNUEvent= fUnfoldedSpectrum->GetBinContent(i+1);
		fUFlux= fNUEvent/(fExposure*binWidth)* 1/(E*log(10)); 
		fUFluxE3= fUFlux* pow(E,3);
		fUFluxE26= fUFlux* pow(E,2.6);

		fUFlux_StatErr= statError/(fExposure*binWidth)* 1/(E*log(10));
		fUFlux_ECalSystErr= systError_ECal/(fExposure*binWidth)* 1/(E*log(10));
		fUFlux_EnergyFDSystErr= systError_EnergyFD/(fExposure*binWidth)* 1/(E*log(10));
		fUFlux_EnergySDResoSystErr= systError_EnergySDReso/(fExposure*binWidth)* 1/(E*log(10));

		fUFlux_ExpSystErr= fUFlux*fExposureSystUncertainty;
		fUFlux_UnfoldingAlgoSystErr= fUFlux*fUnfoldingAlgoSystUncertainty;
		fUFlux_UnfoldingMatrixSystErr= systError_UnfoldingMatrix/(fExposure*binWidth)* 1/(E*log(10));
		fUFlux_UnfoldingProtonMatrixSystErr= systError_UnfoldingProtonMatrix/(fExposure*binWidth)* 1/(E*log(10));
		fUFlux_UnfoldingIronMatrixSystErr= systError_UnfoldingIronMatrix/(fExposure*binWidth)* 1/(E*log(10));
	
	
		double nRec= fRecSpectrum->GetBinContent(recBinId);			
		double nRecErr= fRecSpectrum->GetBinError(recBinId);			
		binWidth= fRecSpectrum->GetBinWidth(recBinId);
		fNEvent= fRecSpectrum->GetBinContent(recBinId); 			
	 	fFlux= fNEvent/(fExposure*binWidth)* 1/(E*log(10));
		fFlux_StatErr= fRecSpectrum->GetBinError(recBinId)/(fExposure*binWidth)* 1/(E*log(10)); 
		fFluxE3= fFlux* pow(E,3);
		fFluxE26= fFlux* pow(E,2.6);
			

		cout<<"bin "<<i+1<<"  nRec="<<nRec<<"  +- "<<nRecErr<<"("<<nRecErr/nRec<<")  nUnfold="<<binContent<<"  +- "<<binError<<"("<<binError/binContent<<")  (recStatErr="<<statError_rec<<")  syst(ECal)="<<systError_ECal<<"  syst(EFD)="<<systError_EnergyFD<<"  syst(ESD)="<<systError_EnergySDReso<< "  syst(matrix)="<<systError_UnfoldingMatrix<<"  syst(unfoldAlgo)="<<systError_UnfoldingAlgo<<"  systError_UnfoldingProtonMatrix="<<systError_UnfoldingProtonMatrix<<"  systError_UnfoldingIronMatrix="<<systError_UnfoldingIronMatrix<<"  totError="<<totError<<endl;
	
		fUnfoldedSpectrumTable->Fill();	
		*/
	}//end loop true bins
	
	if(CCholDecompTriang) CCholDecompTriang->Delete();
	if(CCholDecompTriangTransp) CCholDecompTriangTransp->Delete();
	
	return 0;

}//close ComputeUncertainties()

}//close namespace
