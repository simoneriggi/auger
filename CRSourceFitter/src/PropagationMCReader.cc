/**
* @file PropagationMCReader.cc
* @class PropagationMCReader
* @brief Read the propagation MC data entries, create the transfer matrix, calculate the propagated spectra and composition.
* 
* @author S. Riggi
* @date 22/04/2010
*/

#include "PropagationMCReader.h"
#include "ConfigParser.h"
#include "AnalysisConsts.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TRandom3.h>

#include <TFractionFitter.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>
#include <TObjArray.h>
#include <TMatrixD.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;


ClassImp(CRSourceFitter_ns::PropagationMCReader)

namespace CRSourceFitter_ns {

int PropagationMCReader::fSourceCutoffMode;
double PropagationMCReader::fSourceCutoffShape;


PropagationMCReader::PropagationMCReader(){
    
  for(int j=0;j<fNmassAtSource;j++) 
		fAbundanceAtSource.push_back(1./(double)(fNmassAtSource));//initialize with equal fractions at source
  
  fInjIndexAtSource= 1.1;//initialize with same values used in Denis simulations ==> no reweighting with this value
  fRedShiftIndex= 0.;
  // redshift prop to (1+z)^alfa  
  // uniform scenario ==> no evolution alfa=0? original simulations
  // SFR model ==> alfa= 3
  // StrongEvolutionModel ==> alfa= 4
    
  fEmaxAtSource= 20.5;	
  fNormExtraGalSpectrum= 1.e+23;
  						
  //## GALACTIC
  fInjIndexAtSource_Galactic= 4.0;
  fEmaxAtSource_Galactic= 18.5;
  fNormGalSpectrum= 1.e+24;
  fIncludeGalacticComposition= false;

	fGalacticAbundance.resize(fNmassAtEarth);
  fGalacticAbundance[0]= 0.;
	fGalacticAbundance[1]= 0.;
	fGalacticAbundance[2]= 0.;
  fGalacticAbundance[3]= 1.;

	fCallDiagonalMatrix= false;
	fIsEmaxADependent= false;
	
	fSourceCutoffMode= 1;// exponential cutoff	
	fSourceCutoffShape= 1.0;//cutoff shape
	
	fIsGroupingInCharge= true;

}//close constructor


PropagationMCReader::~PropagationMCReader(){
	
  for(int k=0;k<fNmassAtEarth;k++){
		delete fExpEnergySpectrumAtEarth[k];
		delete fExpEnergySpectrumAtEarth_diag[k];
		delete fExpEnergySpectrumAtEarth_Xmax[k];
		delete fExpEnergySpectrumAtEarth_Xmax_diag[k];
		delete fExpEnergySpectrumAtEarth_Spectrum[k];
		delete fExpEnergySpectrumAtEarth_Spectrum_diag[k];
		delete fExpEnergySpectrumE3AtEarth_Spectrum[k];
		delete fExpEnergySpectrumE3AtEarth_Spectrum_diag[k];	  
  }
	
  delete fTotExpEnergySpectrumAtEarth;
  delete fTotExpEnergySpectrumE3AtEarth;
  delete fTotExpEnergySpectrumAtEarth_Xmax;
  delete fTotExpEnergySpectrumAtEarth_Spectrum;
  delete fTotExpEnergySpectrumE3AtEarth_Spectrum;
  delete fGalacticSpectrumHisto_Spectrum;

  for(int k=0;k<fNmassAtSource;k++){
	  delete fInjectedEnergySpectrumAtSource[k]; 
  }
 
}//close destructor


void PropagationMCReader::Init(){

	//##############################
	//##       INIT HISTO
	//##############################
	double BinEdge_Xmax[fNbins_Xmax+1];
	for(int s=0;s<fNbins_Xmax;s++) BinEdge_Xmax[s]= fEmin_Xmax[s];
	BinEdge_Xmax[fNbins_Xmax]= fEmax_Xmax[fNbins_Xmax-1];

	double BinEdge_PropData[fNbins_PropData+1];
	for(int s=0;s<fNbins_PropData;s++) BinEdge_PropData[s]= fEmin_PropData[s];
	BinEdge_PropData[fNbins_PropData]= fEmax_PropData[fNbins_PropData-1];

	double BinEdge_Spectrum[fNbins_Spectrum+1];
	for(int s=0;s<fNbins_Spectrum;s++) BinEdge_Spectrum[s]= fEmin_Spectrum[s];
	BinEdge_Spectrum[fNbins_Spectrum]= fEmax_Spectrum[fNbins_Spectrum-1];

	//Set binning energy for each nucleus according to their maximum energy
	fGenEminNucleus.resize(fNmassAtSource);
	fGenEmaxNucleus.resize(fNmassAtSource);
	for(int k=0;k<fNmassAtSource;k++) {
		for(int s=0;s<fNbins_PropData;s++){
			if(fEmin_PropData[s]<fMaxGenEnergyAtSourceInSimulation[k]){
				fGenEminNucleus[k].push_back(fEmin_PropData[s]);
				fGenEmaxNucleus[k].push_back(fEmax_PropData[s]);				
			}	
		}//end loop energy bins	
	}//end loop masses at source

	
  TString currentHistoName;
	fFractionsAtEarth_Spectrum.resize(fNmassAtEarth);
	fFractionsAtEarth_Xmax.resize(fNmassAtEarth);
	fFractionsAtEarth_Spectrum_diag.resize(fNmassAtEarth);
	fFractionsAtEarth_Xmax_diag.resize(fNmassAtEarth);

	
  //Initialize vectors
	for(int i=0;i<fNmassAtEarth;i++) {        
  	for(int s=0;s<fNbins_Xmax;s++){
			fFractionsAtEarth_Xmax[i].push_back(0.);//initialize with all 0
			fFractionsAtEarth_Xmax_diag[i].push_back(0.);//initialize with all 0
		}

		for(int s=0;s<fNbins_Spectrum;s++) {
			fFractionsAtEarth_Spectrum[i].push_back(0.);//initialize with all 0
			fFractionsAtEarth_Spectrum_diag[i].push_back(0.);//initialize with all 0
		}
  	  
	  currentHistoName= Form("ExpEnergySpectrumAtEarthHisto_%d",i+1);
	  fExpEnergySpectrumAtEarthHisto= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
	  fExpEnergySpectrumAtEarthHisto->Sumw2();
	  fExpEnergySpectrumAtEarth.push_back(fExpEnergySpectrumAtEarthHisto);

		currentHistoName= Form("ExpEnergySpectrumE3AtEarthHisto_%d",i+1);
	  fExpEnergySpectrumE3AtEarthHisto= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
	  fExpEnergySpectrumE3AtEarthHisto->Sumw2();
	  fExpEnergySpectrumE3AtEarth.push_back(fExpEnergySpectrumE3AtEarthHisto);

		currentHistoName= Form("ExpEnergySpectrumAtEarthDiagHisto_%d",i+1);
	  fExpEnergySpectrumAtEarthDiagHisto= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
	  fExpEnergySpectrumAtEarthDiagHisto->Sumw2();
	  fExpEnergySpectrumAtEarth_diag.push_back(fExpEnergySpectrumAtEarthDiagHisto);
 
	  currentHistoName= Form("ExpEnergySpectrumAtEarthHisto_Xmax_%d",i+1);
	  fExpEnergySpectrumAtEarthHisto_Xmax= new TH1D(currentHistoName,currentHistoName,fNbins_Xmax,BinEdge_Xmax);
	  fExpEnergySpectrumAtEarthHisto_Xmax->Sumw2();
	  fExpEnergySpectrumAtEarth_Xmax.push_back(fExpEnergySpectrumAtEarthHisto_Xmax);

		currentHistoName= Form("ExpEnergySpectrumAtEarthDiagHisto_Xmax_%d",i+1);
	  fExpEnergySpectrumAtEarthDiagHisto_Xmax= new TH1D(currentHistoName,currentHistoName,fNbins_Xmax,BinEdge_Xmax);
	  fExpEnergySpectrumAtEarthDiagHisto_Xmax->Sumw2();
	  fExpEnergySpectrumAtEarth_Xmax_diag.push_back(fExpEnergySpectrumAtEarthDiagHisto_Xmax);

	  currentHistoName= Form("ExpEnergySpectrumAtEarthHisto_Spectrum_%d",i+1);
	  fExpEnergySpectrumAtEarthHisto_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
	  fExpEnergySpectrumAtEarthHisto_Spectrum->Sumw2();
	  fExpEnergySpectrumAtEarth_Spectrum.push_back(fExpEnergySpectrumAtEarthHisto_Spectrum);

		currentHistoName= Form("ExpEnergySpectrumAtEarthDiagHisto_Spectrum_%d",i+1);
	  fExpEnergySpectrumAtEarthDiagHisto_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
	  fExpEnergySpectrumAtEarthDiagHisto_Spectrum->Sumw2();
	  fExpEnergySpectrumAtEarth_Spectrum_diag.push_back(fExpEnergySpectrumAtEarthDiagHisto_Spectrum);
	  
	  currentHistoName= Form("ExpEnergySpectrumE3AtEarthHisto_Spectrum_%d",i+1);
	  fExpEnergySpectrumE3AtEarthHisto_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
	  fExpEnergySpectrumE3AtEarthHisto_Spectrum->Sumw2();
	  fExpEnergySpectrumE3AtEarth_Spectrum.push_back(fExpEnergySpectrumE3AtEarthHisto_Spectrum);

		currentHistoName= Form("ExpEnergySpectrumE3AtEarthDiagHisto_Spectrum_%d",i+1);
	  fExpEnergySpectrumE3AtEarthDiagHisto_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
	  fExpEnergySpectrumE3AtEarthDiagHisto_Spectrum->Sumw2();
	  fExpEnergySpectrumE3AtEarth_Spectrum_diag.push_back(fExpEnergySpectrumE3AtEarthDiagHisto_Spectrum);

  }//end loop masses at Earth

  
  currentHistoName= "TotExpEnergySpectrumAtEarthHisto";
  fTotExpEnergySpectrumAtEarth= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
  fTotExpEnergySpectrumAtEarth->Sumw2();

  currentHistoName= "TotExpEnergySpectrumE3AtEarthHisto";
  fTotExpEnergySpectrumE3AtEarth= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
  fTotExpEnergySpectrumE3AtEarth->Sumw2();
  
 
  currentHistoName=" TotExpEnergySpectrumAtEarthHisto_Xmax";
  fTotExpEnergySpectrumAtEarth_Xmax= new TH1D(currentHistoName,currentHistoName,fNbins_Xmax,BinEdge_Xmax);
  fTotExpEnergySpectrumAtEarth_Xmax->Sumw2();

 
  currentHistoName= "TotExpEnergySpectrumAtEarthHisto_Spectrum";
  fTotExpEnergySpectrumAtEarth_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
  fTotExpEnergySpectrumAtEarth_Spectrum->Sumw2();

  currentHistoName= "TotExpEnergySpectrumE3AtEarthHisto_Spectrum";
  fTotExpEnergySpectrumE3AtEarth_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
  fTotExpEnergySpectrumE3AtEarth_Spectrum->Sumw2();

  currentHistoName= "GalacticSpectrum_Spectrum";
  fGalacticSpectrumHisto_Spectrum= new TH1D(currentHistoName,currentHistoName,fNbins_Spectrum,BinEdge_Spectrum);
  fGalacticSpectrumHisto_Spectrum->Sumw2();

  for(int j=0;j<fNmassAtSource;j++) {
  	currentHistoName=Form("InjectedEnergySpectrumAtSourceHisto_%d",j+1);
	  fInjectedEnergySpectrumAtSourceHisto= new TH1D(currentHistoName,currentHistoName,fNbins_PropData,BinEdge_PropData);
	  fInjectedEnergySpectrumAtSourceHisto->Sumw2();
	  fInjectedEnergySpectrumAtSource.push_back(fInjectedEnergySpectrumAtSourceHisto); 
  }//end loop mass at source

  //Initialize migration matrix
  //Set dimensions (Nmass@Earth x NEnergyBins) x (Nmass@Source x NEnergyBins)
  fMatrixSize= 0;
  for(int k=0;k<fNmassAtSource;k++) {
		fMatrixSize+= fGenEminNucleus[k].size();
  }
	

  //Initialize expected spectrum @ Earth
  fExpSpectrumAtEarth.resize(fNmassAtEarth);
	fExpSpectrumAtEarth_diag.resize(fNmassAtEarth);
  for(int k=0;k<fNmassAtEarth;k++) {
		fExpSpectrumAtEarth[k].resize(fNbins_PropData);
		fExpSpectrumAtEarth_diag[k].resize(fNbins_PropData);  
	}

  //Initialize injected spectrum @ source
  fInjectedSpectrumAtSource.resize(fNmassAtSource);  	  
  for(int k=0;k<fNmassAtSource;k++) fInjectedSpectrumAtSource[k].resize(fGenEminNucleus[k].size());		
	
  //init with zeros
  for(int k=0;k<fNmassAtEarth;k++){
		for(int s=0;s<fNbins_PropData;s++) {
			fExpSpectrumAtEarth[k][s]= 0.;   	
			fExpSpectrumAtEarth_diag[k][s]= 0.;   	
		}	
  }//end loop masses at Earth

  for(int k=0;k<fNmassAtSource;k++){
		for(unsigned int s=0;s<fGenEminNucleus[k].size();s++) fInjectedSpectrumAtSource[k][s]=0.; 
  }//end loop masses at source


	fSpectrumAtEarthMatrix= new TMatrixD(1,fNmassAtEarth*fNbins_PropData);	

}//close PropagationMCReader::Init()




void PropagationMCReader::Reset(){

	//##############################
	//##       RESET HISTO
	//##############################
	for(int i=0;i<fNmassAtEarth;i++){	
		fExpEnergySpectrumAtEarth[i]->Reset();	
		fExpEnergySpectrumE3AtEarth[i]->Reset();	
		fExpEnergySpectrumAtEarth_diag[i]->Reset();			
		fExpEnergySpectrumAtEarth_Xmax[i]->Reset();		
		fExpEnergySpectrumAtEarth_Xmax_diag[i]->Reset();				
		fExpEnergySpectrumAtEarth_Spectrum[i]->Reset();	
		fExpEnergySpectrumAtEarth_Spectrum_diag[i]->Reset();			
		fExpEnergySpectrumE3AtEarth_Spectrum[i]->Reset();	
		fExpEnergySpectrumE3AtEarth_Spectrum_diag[i]->Reset();	
	}

	for(int i=0;i<fNmassAtSource;i++){  
	  fInjectedEnergySpectrumAtSource[i]->Reset();
	}

  fTotExpEnergySpectrumAtEarth->Reset();
	fTotExpEnergySpectrumE3AtEarth->Reset();
	fTotExpEnergySpectrumAtEarth_Xmax->Reset();

	fTotExpEnergySpectrumAtEarth_Spectrum->Reset();
	fTotExpEnergySpectrumE3AtEarth_Spectrum->Reset();

	fGalacticSpectrumHisto_Spectrum->Reset();

}//close PropagationMCReader::ResetMCHisto()




double PropagationMCReader::InjSpectrumFcn(double* x, double* par){
   
	//## Injection spectrum at source dN/dlgE= log(10) x E x dN/dE 
	//## dN/dE= E^(-gamma) x cutoff  ==> dN/dlgE= log(10) x E^(-gamma+1) x cutoff
	//## SOURCE ENERGY CUTOFF
	//# 1) exponential cutoff flux = uncutFlux*exp(-E/ZEmax)
	//# 2) Fermi cutoff       flux = uncutFlux/(1+exp((lgE-log10(Emax*Z))/Wc)); 
	//# 3) kink cutoff        flux= uncutFlux(Z*EMax)*pow(E/(Z*Emax),-alpha) for E>=Z*Emax
	//# 4) kink2 cutoff       flux= uncutFlux * pow(E/(Z*Emax),-alpha+delta) for E>=Z*Emax
	double fval=0.;
	double lgE= x[0];
	double E= pow(10,lgE);

	double NormFactor= par[0];
  double InjIndex= par[1];
  double Z= par[2];
  double Emax= par[3];
  
  double powlaw = 0;
  powlaw= pow(E,-InjIndex+1);

  fval= NormFactor*log(10.)*powlaw;

	switch (fSourceCutoffMode) {
  	case ConfigParser::eExpo://exponential cutoff
			{
    	if(E/Z>pow(10,Emax))
				fval*= TMath::Exp(-E/(Z*pow(10,Emax)))/TMath::Exp(-1);
			}
    	break;
		case ConfigParser::eFermi://Fermi cutoff
			{
			double fermiExpArg= (lgE-log10(Z*pow(10,Emax)))/fSourceCutoffShape;
			double fermiCutoff= 1./(1.+TMath::Exp(fermiExpArg));
			fermiCutoff/= 0.5;
			fval*= fermiCutoff;	
			}
			break;
		case ConfigParser::eKink://kink cutoff
			{
			if(E/Z>=pow(10,Emax)) {
				double Kink= NormFactor*log(10.)*pow(Z*pow(10,Emax),-InjIndex+1);  
				fval= Kink * pow( E/(Z*pow(10,Emax)), -fSourceCutoffShape);//kink cutoff
			}
			}
			break;
		case ConfigParser::eKink2://kink cutoff
			{
			if(E/Z>=pow(10,Emax)) {
				double Kink= NormFactor*log(10.)*pow(Z*pow(10,Emax),-InjIndex+1);  
				fval= Kink * pow( E/(Z*pow(10,Emax)), -InjIndex+1-fSourceCutoffShape );//kink cutoff 2
			}
			}
			break;
		default:
    	cerr<<"PropagationMCReader::InjSpectrumFcn: Invalid source cutoff mode...exit!"<<endl;
			exit(1);
    break;
 	}//close switch


  return fval;

}//close PropagationMCReader::InjSpectrumFcn



void PropagationMCReader::GenerateGalacticSpectrum(){

	//## The galactic spectrum is modeled as a simple power law (index InjIndex_Gal)
	//## with a maximum energy cutoff Emax_Gal 
	double Zgal= 1;

	TF1* InjEnergySpectrumGalactic = new TF1("InjEnergySpectrumGalactic",PropagationMCReader::InjSpectrumFcn,fEmin_Spectrum[0],fEmax_Spectrum[fNbins_Spectrum-1],4);
	InjEnergySpectrumGalactic->SetParameters(fNormGalSpectrum,fInjIndexAtSource_Galactic,Zgal,fEmaxAtSource_Galactic);

	for(int s=0;s<fNbins_Spectrum;s++){
		double binWidth= fGalacticSpectrumHisto_Spectrum->GetBinWidth(s+1);
		double binContent= InjEnergySpectrumGalactic->Integral(fEmin_Spectrum[s],fEmax_Spectrum[s])/binWidth;		
		fGalacticSpectrumHisto_Spectrum->SetBinContent(s+1,binContent);
	}//end loop energy bins

	InjEnergySpectrumGalactic->Delete();

}//close PropagationMCReader::GenerateGalacticSpectrum()



void PropagationMCReader::CalculateMassFractions(){

	//## Calculate expected mass fraction at the Earth
	//## using spectrum energy binning 
	double fGalacticMassContent[fNmassAtEarth][fNbins_Spectrum];
	double fExtraGalacticMassContent[fNmassAtEarth][fNbins_Spectrum];
	double fExtraGalacticMassContent_diag[fNmassAtEarth][fNbins_Spectrum];
	double fTotContent[fNbins_Spectrum];

	//## Initialize with zeros
	for(int s=0;s<fNbins_Spectrum;s++){
		fTotContent[s]=0.; 
		for(int i=0;i<fNmassAtEarth;i++){
  		fGalacticMassContent[i][s]=0.; 
			fExtraGalacticMassContent[i][s]=0.; 
			fExtraGalacticMassContent_diag[i][s]=0.; 
  	}//end loop masses at Earth
	}//end loop energy bins

	//## Calculate spectrum contribution from extragalactic and galactic part
	for(int i=0;i<fNmassAtEarth;i++){		 
		for(int s=0;s<fNbins_Spectrum;s++){
			double aji_GAL= fGalacticSpectrumHisto_Spectrum->GetBinContent(s+1);
  		fGalacticMassContent[i][s]= aji_GAL*fGalacticAbundance[i];
			
			double aji_EG= fExpEnergySpectrumAtEarth_Spectrum[i]->GetBinContent(s+1);
			double aii_EG= fExpEnergySpectrumAtEarth_Spectrum_diag[i]->GetBinContent(s+1);//diagonal elements

			fExtraGalacticMassContent[i][s]= aji_EG;
			fExtraGalacticMassContent_diag[i][s]= aii_EG;
  	}//end loop energy bins		
	}//end loop mass at Earth
	
	//## Calculate total spectrum contribution in each bin
	for(int s=0;s<fNbins_Spectrum;s++){
		for(int i=0;i<fNmassAtEarth;i++){			
			if(fIncludeGalacticComposition) fTotContent[s]+= fGalacticMassContent[i][s] + fExtraGalacticMassContent[i][s];
			else fTotContent[s]+= fExtraGalacticMassContent[i][s];
		}//end loop masses at Earth
	}//end loop energy bins
		
	
	//## Calculate expected fractions
	for(int i=0;i<fNmassAtEarth;i++){
		for(int s=0;s<fNbins_Spectrum;s++){		
			if(fTotContent[s]!=0) {
				if(fIncludeGalacticComposition) 
					fFractionsAtEarth_Spectrum[i][s]= (fGalacticMassContent[i][s]+fExtraGalacticMassContent[i][s])/fTotContent[s];
				else 
					fFractionsAtEarth_Spectrum[i][s]= fExtraGalacticMassContent[i][s]/fTotContent[s];

				fFractionsAtEarth_Spectrum_diag[i][s]= fExtraGalacticMassContent_diag[i][s]/fTotContent[s];
			}//close if
			else {
				cerr<<"PropagationMCReader::CalculateMassFractions(): Total flux= 0 ...exit!"<<endl;
				fFractionsAtEarth_Spectrum[i][s]= 0.;
				fFractionsAtEarth_Spectrum_diag[i][s]= 0.;
				exit(1);
			}
		}//end loop energy bins 		
	}//end loop masses at Earth
	

	
	//## Calculate expected mass fraction at the Earth
	//## using Xmax energy binning 
	double GalContent_Xmax[fNmassAtEarth][fNbins_Xmax];
	double ExtraGalContent_Xmax[fNmassAtEarth][fNbins_Xmax];	
	double ExtraGalContent_Xmax_diag[fNmassAtEarth][fNbins_Xmax];
	double TotContent_Xmax[fNbins_Xmax];

	//## Initialize with zeros
	for(int j=0;j<fNbins_Xmax;j++) {
		TotContent_Xmax[j]=0.;
		for(int k=0;k<fNmassAtEarth;k++){
			GalContent_Xmax[k][j]=0;
			ExtraGalContent_Xmax[k][j]=0;		
			ExtraGalContent_Xmax_diag[k][j]=0;
		}//close for k
	}//close for j

	//## Calculate spectrum contribution from extragalactic and galactic part
	for(int k=0;k<fNmassAtEarth;k++){
		for(int j=0;j<fNbins_Xmax;j++) {
			GalContent_Xmax[k][j]=0;
			ExtraGalContent_Xmax[k][j]=0;
			ExtraGalContent_Xmax_diag[k][j]=0;
		}//end loop Xmax energy bins
		
		for(int s=0;s<fNbins_Spectrum;s++){
			double BinCenter= fExpEnergySpectrumAtEarth_Spectrum[k]->GetBinCenter(s+1);	
			double BinCenter_diag= fExpEnergySpectrumAtEarth_Spectrum_diag[k]->GetBinCenter(s+1);	
			
			for(int j=0;j<fNbins_Xmax;j++){
				if(BinCenter>= fEmin_Xmax[j] && BinCenter<=fEmax_Xmax[j]) {
					GalContent_Xmax[k][j]+= fGalacticMassContent[k][s];
					ExtraGalContent_Xmax[k][j]+= fExtraGalacticMassContent[k][s];
				}

				if(BinCenter_diag>= fEmin_Xmax[j] && BinCenter_diag<=fEmax_Xmax[j]) {
					ExtraGalContent_Xmax_diag[k][j]+= fExtraGalacticMassContent_diag[k][s];
				}
			}//end loop Xmax energy bins
		}//end loop spectrum energy bins
	}//end loop masses at Earth

	//## Calculate total spectrum contribution in each bin
	for(int s=0;s<fNbins_Xmax;s++){
		for(int i=0;i<fNmassAtEarth;i++){
			if(fIncludeGalacticComposition) TotContent_Xmax[s]+= GalContent_Xmax[i][s]+ ExtraGalContent_Xmax[i][s];
			else TotContent_Xmax[s]+= ExtraGalContent_Xmax[i][s];
		}//end loop masses at Earth
	}//end loop Xmax energy bins


	//## Calculate expected fractions
	for(int i=0;i<fNmassAtEarth;i++){
		for(int s=0;s<fNbins_Xmax;s++){
			if(fIncludeGalacticComposition) 
				fFractionsAtEarth_Xmax[i][s]= (GalContent_Xmax[i][s]+ExtraGalContent_Xmax[i][s])/TotContent_Xmax[s];
			else 
				fFractionsAtEarth_Xmax[i][s]= ExtraGalContent_Xmax[i][s]/TotContent_Xmax[s];
			fFractionsAtEarth_Xmax_diag[i][s]= ExtraGalContent_Xmax_diag[i][s]/TotContent_Xmax[s];
		}//end loop Xmax energy bins
	}//end loop masses at Earth

	
}//close PropagationMCReader::CalculateMassFractions()


void PropagationMCReader::CalculateDiagonalMatrix(){

	cout<<"PropagationMCReader::CalculateDiagonalMatrix()"<<endl;
	//########################################
  //##  Calculate diagonal transfer matrix
	//##  p->p, He->He, CNO->CNO, Fe->Fe   
  //########################################
	fTransferDiagMatrix= new TMatrixD(fNmassAtEarth*fNbins_PropData,fMatrixSize);	
	
  for(int i=0;i<fNmassAtEarth*fNbins_PropData;i++) {
		int i_block= (int)(i/fNbins_PropData);
		int j_block= -1;	
		for(int j=0;j<fMatrixSize;j++){
			int low_edge=0;
			int up_edge=0;
			for(int k=0;k<fNmassAtSource;k++) {
				up_edge= low_edge + fGenEminNucleus[k].size();
				if(j>=low_edge && j<up_edge){
					j_block= k;	
					break;
				}
				low_edge+= fGenEminNucleus[k].size();
			}//end loop masses at source
	
			if(i_block==j_block) (*fTransferDiagMatrix)(i,j)= (*fTransferMatrix)(i,j);		
			else (*fTransferDiagMatrix)(i,j)= 0.;
		}//close for j
	}//close for i

	fCallDiagonalMatrix= true;

}//close PropagationMCReader::CalculateDiagonalMatrix()


void PropagationMCReader::CalculateExpectedSpectrumAtEarth(){

	//## Generated number of events are calculated for each component
	//## up to the maximum energy Zx10^21 (different for each mass).
	
	//## Init with zeros
  for(int k=0;k<fNmassAtSource;k++){
		for(unsigned int s=0;s<fGenEminNucleus[k].size();s++){
  		fInjectedSpectrumAtSource[k][s]=0.; 
  	}//end loop energy bins
	}//end loop masses at source
	
  double EmaxInSimulation[7];
  double EminInSimulation= pow(10,18.);
  
	fGeneratedEnergySpectrum.resize(7);


	double Z;
	double A;
	int SourceGroupIndex=0;

	//## Loop over the 7 masses generated by Denis
	//## join some of them according to chosen grouping at source	
	for(int k=0;k<7;k++){
		SourceGroupIndex= -999;

		//## Assign charge and mass
    if(k==0) {//proton
			Z= eProtonZ;
			A= eProtonA;
		}
    else if(k==1) {//He
			Z= eHeliumZ;
			A= eHeliumA;
		}
    else if(k==2) {//C
			Z= eCarbonZ;
			A= eCarbonA;
		}
    else if(k==3) {//O
			Z= eOxygenZ;
			A= eOxygenA;
		}
    else if(k==4) {//Mg
			Z= eMagnesiumZ;
			A= eMagnesiumA;	
		}
    else if(k==5) {//Si
			Z= eSiliconZ;
			A= eSiliconA;
		}
    else if(k==6) {//Fe
			Z= eIronZ;
			A= eIronA;
		}
    else{
    	cerr<<"PropagationMCReader::CalculateExpectedSpectrumAtEarth(): Error assigning charge of parent nucleus...exit!"<<endl;
    	exit(1);
    }

	
		if(fIsEmaxADependent) EmaxInSimulation[k]= A*pow(10,21.);
		else EmaxInSimulation[k]= Z*pow(10,21.);

		//## Find mass group at source
		//## Group in charge or in mass?
		if(fIsGroupingInCharge){
			for(int j=0;j<fNmassAtSource;j++){
				if(Z>=fMinCharge[j] && Z<=fMaxCharge[j]){				
					SourceGroupIndex= j;
					break;
				}	
			}//end loop source group
		}
		else{
			for(int j=0;j<fNmassAtSource;j++){
				if(A>=fMinMass[j] && A<=fMaxMass[j]){				
					SourceGroupIndex= j;
					break;
				}
			}//end loop source group
		}
		

		if(SourceGroupIndex<0){
			cerr<<"PropagationMCReader::CalculateExpectedSpectrumAtEarth(): Error assigning charge of son nucleus...exit!"<<endl;
			exit(1);
		}	
	
		InjEnergySpectrum = new TF1(Form("InjEnergySpectrum_%d",k+1),PropagationMCReader::InjSpectrumFcn,log10(EminInSimulation),log10(EmaxInSimulation[k]),4);
		InjEnergySpectrum->SetParameters(fNormExtraGalSpectrum,fInjIndexAtSource,Z,fEmaxAtSource);

		fGeneratedEnergySpectrum[k]= InjEnergySpectrum;

		for(unsigned int s=0;s<fGenEminNucleus[SourceGroupIndex].size();s++){
			if(fGenEmaxNucleus[SourceGroupIndex][s]> log10(EmaxInSimulation[k])){
				if(log10(EmaxInSimulation[k])>fGenEminNucleus[SourceGroupIndex][s]) fInjectedSpectrumAtSource[SourceGroupIndex][s]+= InjEnergySpectrum->Integral(fGenEminNucleus[SourceGroupIndex][s],log10(EmaxInSimulation[k]));
			}//close if
			else fInjectedSpectrumAtSource[SourceGroupIndex][s]+= InjEnergySpectrum->Integral(fGenEminNucleus[SourceGroupIndex][s],fGenEmaxNucleus[SourceGroupIndex][s]);

		}//end loop energy bins
		
		//delete InjEnergySpectrum;

	}//end loop simulated nuclei
	

	//## Multiply spectrum group at source by the relative abundance
	for(int k=0;k<fNmassAtSource;k++){
		for(unsigned int s=0;s<fGenEminNucleus[k].size();s++){
			fInjectedSpectrumAtSource[k][s]*= fAbundanceAtSource[k];
			fInjectedEnergySpectrumAtSource[k]->SetBinContent(s+1,fInjectedSpectrumAtSource[k][s]);
		}//end loop energy bins
	}//end loop masses at source


	//Create also TMatrixD with generated spectrum
	TMatrixD* theGenSpectrumMatrix= new TMatrixD(fMatrixSize,1);	
	
  //fill the matrix 
	int index=0;

  for(int k=0;k<fNmassAtSource;k++){
		for(unsigned int s=0;s<fGenEminNucleus[k].size();s++){
			(*theGenSpectrumMatrix)(index,0)= fInjectedSpectrumAtSource[k][s];
			index++;
		}//end loop energy bins
	}//end loop masses at source

	if(!fCallDiagonalMatrix) CalculateDiagonalMatrix();
	

  //## Calculate spectrum @ Earth
	TMatrixD* thePropSpectrumMatrix= new TMatrixD(1,fNmassAtEarth*fNbins_PropData);	
	TMatrixD* theDiagPropSpectrumMatrix= new TMatrixD(1,fNmassAtEarth*fNbins_PropData);	
	
	thePropSpectrumMatrix->Mult(*fTransferMatrix,*theGenSpectrumMatrix); //calculate prop spectrum @ Earth
	theDiagPropSpectrumMatrix->Mult(*fTransferDiagMatrix,*theGenSpectrumMatrix); //calculate prop spectrum @ Earth for not fragmented primaries (diagonal elements)

	fSpectrumAtEarthMatrix->Mult(*fTransferMatrix,*theGenSpectrumMatrix); //calculate prop spectrum @ Earth

  TMatrixD thePropSpectrumMatrix_Transp(fNmassAtEarth*fNbins_PropData,1);
	thePropSpectrumMatrix_Transp.Transpose(*thePropSpectrumMatrix);//transpose of thePropSpectrumMatrix

	TMatrixD theDiagPropSpectrumMatrix_Transp(fNmassAtEarth*fNbins_PropData,1);
	theDiagPropSpectrumMatrix_Transp.Transpose(*theDiagPropSpectrumMatrix);//transpose of theDiagPropSpectrumMatrix

	
	index= 0; 
  int blockIndex= 0;	
  for(int k=0;k<fNmassAtEarth;k++) {
  	for(int s=0;s<fNbins_PropData;s++) {
  		index= s+blockIndex;
  		fExpSpectrumAtEarth[k][s]= (*thePropSpectrumMatrix)(0,index);
			fExpSpectrumAtEarth_diag[k][s]= (*theDiagPropSpectrumMatrix)(0,index);
  	}//end loop energy bins
  	blockIndex+= fNbins_PropData;
	}//end loop masses at Earth
  	

	//## Fill spectrum histos
	for(int k=0;k<fNmassAtEarth;k++){
		for(int s=0;s<fNbins_PropData;++s){
			fExpEnergySpectrumAtEarth[k]->SetBinContent(s+1,fExpSpectrumAtEarth[k][s]);
			fExpEnergySpectrumAtEarth_diag[k]->SetBinContent(s+1,fExpSpectrumAtEarth_diag[k][s]);
		}//end loop energy bins
	}//end loop masses at Earth
	

  for(int k=0;k<fNmassAtEarth;k++){
		fTotExpEnergySpectrumAtEarth->Add(fExpEnergySpectrumAtEarth[k]);
  }//end loop masses at Earth
	
	//## Calculate exp total spectrum x E^3
  for(int s=0;s<fNbins_PropData;s++){
		double BinCenter= fTotExpEnergySpectrumAtEarth->GetBinCenter(s+1);
		double BinContent= fTotExpEnergySpectrumAtEarth->GetBinContent(s+1);
		double weightE3= pow(pow(10,BinCenter),2);

		double BinContentE3= BinContent*weightE3;
		fTotExpEnergySpectrumE3AtEarth->Fill(BinCenter,BinContentE3);
	}//end loop energy bins

	
	//## Fill mass & total spectrum @ Earth using Xmax binning
	double SumBinContentXmax[fNbins_Xmax];
	double SumBinContentXmax_diag[fNbins_Xmax];

	for(int k=0;k<fNmassAtEarth;k++){
		for(int j=0;j<fNbins_Xmax;j++) {
			SumBinContentXmax[j]=0;
			SumBinContentXmax_diag[j]=0;
		}

		for(int s=0;s<fNbins_PropData;s++){
			double BinCenter= fExpEnergySpectrumAtEarth[k]->GetBinCenter(s+1);	
			double BinContent= fExpEnergySpectrumAtEarth[k]->GetBinContent(s+1);	

			double BinCenter_diag= fExpEnergySpectrumAtEarth_diag[k]->GetBinCenter(s+1);	
			double BinContent_diag= fExpEnergySpectrumAtEarth_diag[k]->GetBinContent(s+1);	
			
			for(int j=0;j<fNbins_Xmax;j++){
				if(BinCenter>=fEmin_Xmax[j] && BinCenter<=fEmax_Xmax[j]) SumBinContentXmax[j]+= BinContent;
				if(BinCenter_diag>=fEmin_Xmax[j] && BinCenter_diag<=fEmax_Xmax[j]) SumBinContentXmax_diag[j]+= BinContent_diag;
			}//end loop Xmax energy bins
		}//end loop energy bins

		for(int j=0;j<fNbins_Xmax;j++) {
			fExpEnergySpectrumAtEarth_Xmax[k]->SetBinContent(j+1,SumBinContentXmax[j]);
			fExpEnergySpectrumAtEarth_Xmax_diag[k]->SetBinContent(j+1,SumBinContentXmax_diag[j]);
		}//end loop Xmax energy bins

	}//end loop masses at Earth

	for(int k=0;k<fNmassAtEarth;k++) fTotExpEnergySpectrumAtEarth_Xmax->Add(fExpEnergySpectrumAtEarth_Xmax[k]);


	//## Fill mass & total spectrum @ Earth using Spectrum binning	
	double SumBinContent[fNbins_Spectrum];
	double SumBinContent_diag[fNbins_Spectrum];

	for(int k=0;k<fNmassAtEarth;k++){
		for(int j=0;j<fNbins_Spectrum;j++) {
			SumBinContent[j]=0;
			SumBinContent_diag[j]=0;
		}

		for(int s=0;s<fNbins_PropData;s++){
			double BinCenter= fExpEnergySpectrumAtEarth[k]->GetBinCenter(s+1);	
			double BinContent= fExpEnergySpectrumAtEarth[k]->GetBinContent(s+1);

			double BinCenter_diag= fExpEnergySpectrumAtEarth_diag[k]->GetBinCenter(s+1);	
			double BinContent_diag= fExpEnergySpectrumAtEarth_diag[k]->GetBinContent(s+1);
			
			for(int j=0;j<fNbins_Spectrum;j++){
				if(BinCenter>=fEmin_Spectrum[j] && BinCenter<=fEmax_Spectrum[j]) SumBinContent[j]+= BinContent;
				if(BinCenter_diag>=fEmin_Spectrum[j] && BinCenter_diag<=fEmax_Spectrum[j]) SumBinContent_diag[j]+= BinContent_diag;
			}//end loop spectrum energy bins
		}//end loop energy bins

		for(int j=0;j<fNbins_Spectrum;j++) {
			fExpEnergySpectrumAtEarth_Spectrum[k]->SetBinContent(j+1,SumBinContent[j]);
			fExpEnergySpectrumAtEarth_Spectrum_diag[k]->SetBinContent(j+1,SumBinContent_diag[j]);
		}//end loop spectrum energy bins
	
	}//end loop masses at Earth

  for(int k=0;k<fNmassAtEarth;k++) fTotExpEnergySpectrumAtEarth_Spectrum->Add(fExpEnergySpectrumAtEarth_Spectrum[k]);

		
	for(int k=0;k<fNmassAtEarth;k++){
		for(int s=0;s<fNbins_Spectrum;s++) {
			double BinCenter= fExpEnergySpectrumAtEarth_Spectrum[k]->GetBinCenter(s+1);
			double BinContent= fExpEnergySpectrumAtEarth_Spectrum[k]->GetBinContent(s+1);

			double BinCenter_diag= fExpEnergySpectrumAtEarth_Spectrum_diag[k]->GetBinCenter(s+1);
			double BinContent_diag= fExpEnergySpectrumAtEarth_Spectrum_diag[k]->GetBinContent(s+1);

			double weightE3= pow(pow(10,BinCenter),2);
			double BinContentE3= BinContent*weightE3;

			double weightE3_diag= pow(pow(10,BinCenter_diag),2);
			double BinContentE3_diag= BinContent_diag*weightE3_diag;

			fExpEnergySpectrumE3AtEarth_Spectrum[k]->Fill(BinCenter,BinContentE3);
			fExpEnergySpectrumE3AtEarth_Spectrum_diag[k]->Fill(BinCenter_diag,BinContentE3_diag);
		}//end loop spectrum energy bins


		for(int s=0;s<fNbins_PropData;s++) {
			double BinCenter= fExpEnergySpectrumAtEarth[k]->GetBinCenter(s+1);
			double BinContent= fExpEnergySpectrumAtEarth[k]->GetBinContent(s+1);
			double weightE3= pow(pow(10,BinCenter),2);
			double BinContentE3= BinContent*weightE3;

			fExpEnergySpectrumE3AtEarth[k]->Fill(BinCenter,BinContentE3);
		}//end loop energy bins
	}//end loop masses at Earth

	
	//Calculate exp total spectrum x E^3
  for(int s=0;s<fNbins_Spectrum;s++){
		double BinCenter= fTotExpEnergySpectrumAtEarth_Spectrum->GetBinCenter(s+1);
		double BinContent= fTotExpEnergySpectrumAtEarth_Spectrum->GetBinContent(s+1);
		double weightE3= pow(pow(10,BinCenter),2);
		double BinContentE3= BinContent*weightE3;
		fTotExpEnergySpectrumE3AtEarth_Spectrum->Fill(BinCenter,BinContentE3);
	}//end loop spectrum energy bins
	

	theGenSpectrumMatrix->Delete();	
	thePropSpectrumMatrix->Delete();
	
}//close function


void PropagationMCReader::StoreHisto(){

	
 	fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
 	fOutputFile->cd();

	cout<<"Writing spectrum histo @Earth"<<endl;
	for(int i=0;i<fNmassAtEarth;i++){
		fExpEnergySpectrumAtEarth[i]->Write();
		fExpEnergySpectrumAtEarth_diag[i]->Write();	
		fExpEnergySpectrumE3AtEarth[i]->Write();	

		fExpEnergySpectrumAtEarth_Xmax[i]->Write();
		fExpEnergySpectrumAtEarth_Xmax_diag[i]->Write();

		fExpEnergySpectrumAtEarth_Spectrum[i]->Write();
		fExpEnergySpectrumAtEarth_Spectrum_diag[i]->Write();
		fExpEnergySpectrumE3AtEarth_Spectrum[i]->Write();
		fExpEnergySpectrumE3AtEarth_Spectrum_diag[i]->Write();
	}
	
	fTotExpEnergySpectrumAtEarth->Write();
	fTotExpEnergySpectrumE3AtEarth->Write();

	
	fTotExpEnergySpectrumAtEarth_Xmax->Write();
	fTotExpEnergySpectrumAtEarth_Spectrum->Write();
	fTotExpEnergySpectrumE3AtEarth_Spectrum->Write();

	fGalacticSpectrumHisto_Spectrum->Write();
		
	cout<<"Writing gen spectrum function"<<endl;
	for(int j=0;j<7;j++){	
	  fGeneratedEnergySpectrum[j]->Write();	
	}

	cout<<"Writing gen spectrum hist @ source"<<endl;
	for(int j=0;j<fNmassAtSource;j++){	 
    fInjectedEnergySpectrumAtSource[j]->Write();
	}
	
	cout<<"Writing transfer matrix"<<endl;
	fTransferMatrix->Write();
	

 	fOutputFile->Close();
 
}//close function

}//close namespace

