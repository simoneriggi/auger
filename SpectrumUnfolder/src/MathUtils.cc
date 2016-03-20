/**
* @file MathUtils.cc
* @class MathUtils
* @brief Math utilities
* 
* @author S. Riggi
* @date 20/09/2011
*/


#include <MathUtils.h>

#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TH2.h>
#include <TMath.h>

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

//Reso Model classes
ClassImp(Unfolder_ns::ResoPars)
ClassImp(Unfolder_ns::ConstResoPars)
ClassImp(Unfolder_ns::LinearResoPars)
ClassImp(Unfolder_ns::Pol2ResoPars)

//Trigger Model classes
ClassImp(Unfolder_ns::TriggerPars)
ClassImp(Unfolder_ns::ConstTriggerPars)
ClassImp(Unfolder_ns::SigmoidTriggerPars)

//Spectrum Model classes
ClassImp(Unfolder_ns::FitPar)
ClassImp(Unfolder_ns::FitPars)
ClassImp(Unfolder_ns::SpectrumPars)
ClassImp(Unfolder_ns::PowerLawPars)
ClassImp(Unfolder_ns::BrokenPowerLawsPars)
ClassImp(Unfolder_ns::SmoothCutoffPowerLaws)

//Mth util class
ClassImp(Unfolder_ns::MathUtils)

namespace Unfolder_ns {

int ResoPars::nPars;
//int ConstResoPars::nPars;
//int LinearResoPars::nPars;
//int Pol2ResoPars::nPars;
int TriggerPars::nPars;
int SpectrumPars::nPars;


MathUtils::MathUtils(){
	
}//close constructor


MathUtils::~MathUtils(){

}//close destructor


double MathUtils::ResponseModel(double* x, double* par){

	double lgE_true= x[0];
	double lgE_rec= x[1];
	int par_counter= 0;

	//## Set spectrum pars
	int spectrumModel= par[par_counter++];
	int nSpectrumPars= 0;
	
	if(spectrumModel==Unfolder_ns::ePowerLaw){
		nSpectrumPars= Unfolder_ns::PowerLawPars::GetParNumber();
	}
	else if(spectrumModel==Unfolder_ns::eBrokenPowerLaws){
		nSpectrumPars= Unfolder_ns::BrokenPowerLawsPars::GetParNumber();
	}
	else if(spectrumModel==Unfolder_ns::eFlat){
		nSpectrumPars= Unfolder_ns::FlatSpectrumPars::GetParNumber();
	}
	else if(spectrumModel==Unfolder_ns::eSmoothBrokenPowerLaws){
		nSpectrumPars= Unfolder_ns::SmoothCutoffPowerLaws::GetParNumber();
	}
	else{
		cerr<<"MathUtils::ResponseModel(): ERROR: Invalid spectrum model given!"<<endl;	
		return 0;
	}

	
	double spectrumPars[nSpectrumPars];
	for(int i=0;i<nSpectrumPars;i++){
		spectrumPars[i]= par[par_counter];
		//cout<<"spectrumPars["<<i<<"]="<<spectrumPars[i]<<endl;
		par_counter++;
	}

	double spectrum= 0;
	
	if(spectrumModel==Unfolder_ns::ePowerLaw){
		spectrum= PowerLawSpectrum(&lgE_true,spectrumPars);
	}
	else if(spectrumModel==Unfolder_ns::eBrokenPowerLaws){
		spectrum= BrokenPowerLawSpectrum(&lgE_true,spectrumPars);
	}	
	else if(spectrumModel==Unfolder_ns::eSmoothBrokenPowerLaws){
		spectrum= SmoothCutoffPowerLawSpectrum(&lgE_true,spectrumPars);
	}	
	else if(spectrumModel==Unfolder_ns::eFlat){
		spectrum= 1;
	}

	//## Set reco bias
	int biasModel= par[par_counter++];
	int nBiasPars= 0;
	if(biasModel==Unfolder_ns::eConst){
		nBiasPars= Unfolder_ns::ConstResoPars::GetParNumber();
	}	
	else if(biasModel==Unfolder_ns::eLinear){
		nBiasPars= Unfolder_ns::LinearResoPars::GetParNumber();
	}
	else if(biasModel==Unfolder_ns::ePol2){
		nBiasPars= Unfolder_ns::Pol2ResoPars::GetParNumber();
	}

	double biasPars[nBiasPars];
	
	for(int i=0;i<nBiasPars;i++){
		biasPars[i]= par[par_counter];	
		//cout<<"biasPars["<<i<<"]="<<biasPars[i]<<endl;
		par_counter++;
	}

	double bias= BiasModel(&lgE_true,biasPars);
	
	//## Set reco sigma
	int sigmaModel= par[par_counter++];
	int nSigmaPars= 0;
	if(sigmaModel==Unfolder_ns::eConst){
		nSigmaPars= Unfolder_ns::ConstResoPars::GetParNumber();
	}	
	else if(sigmaModel==Unfolder_ns::eLinear){
		nSigmaPars= Unfolder_ns::LinearResoPars::GetParNumber();
	}
	else if(sigmaModel==Unfolder_ns::ePol2){
		nSigmaPars= Unfolder_ns::Pol2ResoPars::GetParNumber();
	}

	//cout<<"nSpectrumPars="<<nSpectrumPars<<" nSigmaPars="<<nSigmaPars<<" nBiasPars"<<nBiasPars<<endl;



	double sigmaPars[nSigmaPars];
	
	for(int i=0;i<nSigmaPars;i++){
		sigmaPars[i]= par[par_counter];
		//cout<<"sigmaPars["<<i<<"]="<<sigmaPars[i]<<endl;
		par_counter++;
	}

	double sigma= ResolutionModel(&lgE_true,sigmaPars);//Relative sigma 
	
	//## Set trigger eff
	int triggerModel= par[par_counter++];
	int nTriggerPars= 0;
	if(triggerModel==Unfolder_ns::eConstTrigger){
		nTriggerPars= Unfolder_ns::ConstTriggerPars::GetParNumber();
	}	
	else if(triggerModel==Unfolder_ns::eSigmoidTrigger){
		nTriggerPars= Unfolder_ns::SigmoidTriggerPars::GetParNumber();
	}

	double triggerPars[nTriggerPars];
	for(int i=0;i<nTriggerPars;i++){
		triggerPars[i]= par[par_counter];
		par_counter++;
	}

	double triggerEfficiency= 1;
	if(triggerModel==Unfolder_ns::eConstTrigger) triggerEfficiency= triggerPars[0];
	else if(triggerModel==Unfolder_ns::eSigmoidTrigger) triggerEfficiency= TriggerEfficiencyModel(&lgE_rec,triggerPars);
	
	
	//## Set energy response 
  double arg = (lgE_rec - (lgE_true+bias) )/sigma;
	double gaussNorm= 1./(sigma*sqrt(2.*TMath::Pi()));
	double gaussResponse= gaussNorm*TMath::Exp(-0.5*arg*arg);

  double response = spectrum* gaussResponse * triggerEfficiency;
	 
	return response;

}//close ResponseModel()


double MathUtils::ResolutionModel(double* x, double* par){

	double lgE= x[0];

	int nPars= sizeof(par)/sizeof(double);
	double lgEFact[nPars];
	for(int i=0;i<nPars;i++){
		if(i==0) lgEFact[i]= 1;
		else lgEFact[i]= lgE*lgEFact[i-1];
	}

	double sigma= 0.;
	for(int i=0;i<nPars;i++){
		sigma+= par[i]*lgEFact[i];
	}
	
	return sigma;

}//close ResolutionModel()


double MathUtils::BiasModel(double* x, double* par){

	double lgE= x[0];

	int nPars= sizeof(par)/sizeof(double);
	double lgEFact[nPars];
	for(int i=0;i<nPars;i++){
		if(i==0) lgEFact[i]= 1;
		else lgEFact[i]= lgE*lgEFact[i-1];
	}

	double bias= 0.;
	for(int i=0;i<nPars;i++){
		bias+= par[i]*lgEFact[i];
	}
	
	return bias;

}//close BiasModel()


double MathUtils::TriggerEfficiencyModel(double* x, double* par){

	double lgE= x[0];
	double norm= par[0];
	double lgE0= par[1];
	double k= par[2];

	double fval= norm*0.5*(1. + TMath::Erf((lgE-lgE0)/k));
	return fval;

}//close TriggerEfficiencyModel()


double MathUtils::PowerLawSpectrum(double* x, double* par){

	double lgE= x[0];
	double E= pow(10,lgE);
	double norm= par[0];
	double gamma= par[1]; 
	
	double spectrum= norm*log(10)*pow(E,-gamma+1);
	
	return spectrum;

}//close PowerLawSpectrum()


double MathUtils::BrokenPowerLawSpectrum(double* x, double* par){

	double lgE= x[0];	
	
	double norm= par[0];
	double gamma1= par[1];
	double gamma2= par[2];
	double gamma3= par[3];
	double LgEb1= par[4];
	double LgEb2= par[5];
	
	double norm1= 1;
	double norm2= norm1* pow(pow(10,LgEb1),gamma2-gamma1);
	double norm3= norm2* pow(pow(10,LgEb2),gamma3-gamma2);
	

	double spectrum1Par[2]= {norm1,gamma1};
	double spectrum2Par[2]= {norm2,gamma2};
	double spectrum3Par[2]= {norm3,gamma3};
		
	double spectrum1= PowerLawSpectrum(&lgE,spectrum1Par);
	double spectrum2= PowerLawSpectrum(&lgE,spectrum2Par);
	double spectrum3= PowerLawSpectrum(&lgE,spectrum3Par);
	

	double spectrum= 0.;
	if(lgE<LgEb1){
		spectrum= spectrum1;
	}
	if(lgE>=LgEb1 && lgE<LgEb2){
		spectrum= spectrum2;
	}
	if(lgE>=LgEb2){
		spectrum= spectrum3;
	}

	return norm*spectrum;

}//close BrokenPowerLawSpectrum()


double MathUtils::SmoothCutoffPowerLawSpectrum(double* x, double* par){

	double lgE= x[0];
	double E= pow(10,lgE);

	double norm= pow(10,par[0]);
	double gamma1= par[1];
	double gamma2= par[2];
	double LgEb1= par[3];
	double LgEb2= par[4];
	double Wc= par[5];
	
	double cutoffArg= (lgE-LgEb2)/Wc;
	double Cutoff= 1./(1.+TMath::Exp(cutoffArg));

	double Eb1= pow(10,LgEb1);
	//double Eb2= pow(10,LgEb2);
	
	double NormBeforeBreak= norm* pow(Eb1,gamma1-gamma2)*Cutoff;
	double SpectrumAfterBreak= norm* log(10)*pow(E,-gamma2+1)*Cutoff;
	double SpectrumBeforeBreak= NormBeforeBreak* log(10)*pow(E,-gamma1+1);
	
	double Spectrum= 0.;	

	if(lgE<LgEb1){
		Spectrum= SpectrumBeforeBreak;
	}
	else{
		Spectrum= SpectrumAfterBreak;
	}
	
	double fval= Spectrum;

	return fval;

}//close SmoothCutoffPowerLawSpectrum()


double MathUtils::GetPowerLawIntegral(double gamma,double lgEMin, double lgEMax){

	double EMin= pow(10,lgEMin);
	double EMax= pow(10,lgEMax);
	
	double SpectrumInt= (pow(EMax,-gamma+1)-pow(EMin,-gamma+1))/(-gamma+1);

	return SpectrumInt;

}//close GetPowerLawIntegral()



double MathUtils::GetBrokenPowerLawIntegral(double Gamma,double Gamma2,double Gamma3,double Break,double Cutoff,double lgEMin, double lgEMax) {
	double norm1= 1;
	double norm2= norm1* pow(pow(10,Break),Gamma2-Gamma);
	double norm3= norm2* pow(pow(10,Cutoff),Gamma3-Gamma2);
	double SpectrumInt1= norm1*MathUtils::GetPowerLawIntegral(Gamma,lgEMin,Break);
	double SpectrumInt2= norm2*MathUtils::GetPowerLawIntegral(Gamma2,Break,Cutoff);
	double SpectrumInt3= norm3*MathUtils::GetPowerLawIntegral(Gamma3,Cutoff,lgEMax);
	double SpectrumInt= SpectrumInt1+SpectrumInt2+SpectrumInt3;
	return SpectrumInt;
}

}//close namespace

