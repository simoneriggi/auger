/**
* @file DataReader.cc
* @class DataReader
* @brief Read the hybrid data entries and create histograms for the composition fit
* 
* @author Dr. Simone Riggi
* @date 22/04/2010
*/

#include "DataReader.h"
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
#include <vector>

using namespace std;

ClassImp(CRSourceFitter_ns::DataReader)

namespace CRSourceFitter_ns {

DataReader::DataReader(){
  
	fDataFileName= "Data.root";
	fSpectrumTableFileName= "data/TableSpectrumAuger_ICRC09_fit.txt";
  fOutputFileName= "DataOutput.root";
   
}//close constructor

DataReader::~DataReader(){
  
	cout<<"Free DataReader() allocated memory"<<endl;
  
  for(int s=0;s<fNbins_Xmax;s++){
		//if(fXmaxData_fit[s]) delete fXmaxData_fit[s];
	 	//if(fEnergyData_fit[s]) delete fEnergyData_fit[s];
		if(fXmaxData_fit[s]) fXmaxData_fit[s]->Delete();
	 	if(fEnergyData_fit[s]) fEnergyData_fit[s]->Delete();

		for(int j=0;j<fNmassAtEarth;j++){
			//if(fXmaxMC_fit[j][s]) delete fXmaxMC_fit[j][s];
			//if(fEnergyMC_fit[j][s]) delete fEnergyMC_fit[j][s];
			if(fXmaxMC_fit[j][s]) fXmaxMC_fit[j][s]->Delete();
			if(fEnergyMC_fit[j][s]) fEnergyMC_fit[j][s]->Delete();
			if(fGenXmaxMC_fit[j][s]) fGenXmaxMC_fit[j][s]->Delete();
			if(fGenEnergyMC_fit[j][s]) fGenEnergyMC_fit[j][s]->Delete();		
		}//close for masses @ Earth
  }//end loop energy bin

	if(fEnergyHistoData) fEnergyHistoData->Delete();
	if(fXmaxHistoData) fXmaxHistoData->Delete();
	if(fEnergyHistoMC) fEnergyHistoMC->Delete();
	if(fXmaxHistoMC) fXmaxHistoMC->Delete();
	if(fGenEnergyHistoMC) fGenEnergyHistoMC->Delete();
	if(fGenXmaxHistoMC) fGenXmaxHistoMC->Delete();

  //if(fAugerEnergySpectrum) delete fAugerEnergySpectrum;
  //if(fAugerEnergySpectrumE3) delete fAugerEnergySpectrumE3;
  //if(fAugerEnergySpectrumEvents) delete fAugerEnergySpectrumEvents;
	if(fAugerEnergySpectrum) fAugerEnergySpectrum->Delete();
  if(fAugerEnergySpectrumE3) fAugerEnergySpectrumE3->Delete();
  if(fAugerEnergySpectrumEvents) fAugerEnergySpectrumEvents->Delete();

}//close destructor


void DataReader::Init(){

	//Initialize histograms
	int Nbins_Energy= 100;
  int Nbins_Xmax= 100;
	
  double XmaxMin= 0.;
  double XmaxMax= 2000.;
  
  for(int j=0;j<fNmassAtEarth;j++){
		//MC histo vector
		fXmaxMC_fit.push_back ( std::vector<TH1D*>() );
    fEnergyMC_fit.push_back ( std::vector<TH1D*>() );
		fGenXmaxMC_fit.push_back ( std::vector<TH1D*>() );
    fGenEnergyMC_fit.push_back ( std::vector<TH1D*>() );
	}

  TString currentHistoName;
	
  for(int s=0;s<fNbins_Xmax;s++){
		//XmaxHisto
		currentHistoName= Form("fXmaxHistoData_%d",s);
  	fXmaxHistoData= new TH1D(currentHistoName,currentHistoName,Nbins_Xmax,XmaxMin,XmaxMax);
  	fXmaxHistoData->Sumw2();
  	fXmaxData_fit.push_back(fXmaxHistoData);
  	 	
		//EnergyHisto
		currentHistoName= Form("fEnergyHistoData_%d",s);
  	fEnergyHistoData= new TH1D(currentHistoName,currentHistoName,Nbins_Energy,pow(10,fEmin_Xmax[s]),pow(10,fEmax_Xmax[s]));
  	fEnergyHistoData->Sumw2();
  	fEnergyData_fit.push_back(fEnergyHistoData);


		for(int j=0;j<fNmassAtEarth;j++){	
  		//XmaxHisto
			currentHistoName= Form("fXmaxHistoMC_%d_%d",j+1,s+1);
  		fXmaxHistoMC= new TH1D(currentHistoName,currentHistoName,Nbins_Xmax,XmaxMin,XmaxMax);
  		fXmaxHistoMC->Sumw2();
  		fXmaxMC_fit[j].push_back(fXmaxHistoMC);
  		
			//GenXmaxHisto
			currentHistoName= Form("fGenXmaxHistoMC_%d_%d",j+1,s+1);
  		fGenXmaxHistoMC= new TH1D(currentHistoName,currentHistoName,Nbins_Xmax,XmaxMin,XmaxMax);
  		fGenXmaxHistoMC->Sumw2();
  		fGenXmaxMC_fit[j].push_back(fGenXmaxHistoMC);

  		//EnergyHisto
			currentHistoName=Form("fEnergyHistoMC_%d_%d",j+1,s+1);
  		fEnergyHistoMC= new TH1D(currentHistoName,currentHistoName,Nbins_Energy,pow(10,fEmin_Xmax[s]),pow(10,fEmax_Xmax[s]));
  		fEnergyHistoMC->Sumw2();
  		fEnergyMC_fit[j].push_back(fEnergyHistoMC);
  	 	  	
			//GenEnergyHisto
			currentHistoName= Form("fGenEnergyHistoMC_%d_%d",j+1,s+1);
  		fGenEnergyHistoMC= new TH1D(currentHistoName,currentHistoName,Nbins_Energy,pow(10,fEmin_Xmax[s]),pow(10,fEmax_Xmax[s]));
  		fGenEnergyHistoMC->Sumw2();
  		fGenEnergyMC_fit[j].push_back(fGenEnergyHistoMC);
  	}//close for masses @ Earth
  	   	  	
  }//close for energy bins

	double BinEdge_Spectrum[fNbins_Spectrum+1];
	for(int s=0;s<fNbins_Spectrum;s++) BinEdge_Spectrum[s]= fEmin_Spectrum[s];
	BinEdge_Spectrum[fNbins_Spectrum]= fEmax_Spectrum[fNbins_Spectrum-1];

	fAugerEnergySpectrum= new TH1D("fAugerEnergySpectrum","fAugerEnergySpectrum",fNbins_Spectrum,BinEdge_Spectrum);
  fAugerEnergySpectrum->Sumw2();

  fAugerEnergySpectrumE3= new TH1D("fAugerEnergySpectrumE3","fAugerEnergySpectrumE3",fNbins_Spectrum,BinEdge_Spectrum);
  fAugerEnergySpectrumE3->Sumw2();

  fAugerEnergySpectrumEvents= new TH1D("fAugerEnergySpectrumEvents","fAugerEnergySpectrumEvents",fNbins_Spectrum,BinEdge_Spectrum);
  fAugerEnergySpectrumEvents->Sumw2();


}//close Init()

void DataReader::ReadSpectrumData(){

	cout<<"DataReader::ReadSpectrumData()"<<endl;

  ifstream in_spectrumdata;
  
  in_spectrumdata.open(fSpectrumTableFileName.c_str());
  if ( ! in_spectrumdata.good() ) {
  	string errMsg = " DataReader::ReadSpectrumData() - Error reading spectrum data table " + fSpectrumTableFileName
                  + " \n  ********* check file path or spectrum data table!!! ********* ";
    throw std::runtime_error(errMsg);
	  exit(1);
  }//close if
	
	double EnergyFlux;
	double EnergyFluxUpErr;
	double EnergyFluxLowErr;
	double EnergyFluxE;
	double EnergyFluxE3;
	double Nev_FD;
	double Nev_SD;
	double Nev_comb;
	double Exposure;
	double BinEnergy;

	fExposure.clear();
	fExposure.resize(0);

	cout<<"*** SPECTRUM TABLE DUMP ***"<<endl;
	for(int s=0;s<fNbins_Spectrum;s++){
		in_spectrumdata >> BinEnergy >> EnergyFlux >> EnergyFluxUpErr >> EnergyFluxLowErr >> EnergyFluxE >> EnergyFluxE3 >> Nev_FD >> Nev_SD >> Nev_comb >> Exposure;
		EnergyFlux*= 1.e-30;
		EnergyFluxUpErr*= 1.e-30;

		fAugerEnergySpectrum->SetBinContent(s+1,EnergyFlux);
		fAugerEnergySpectrum->SetBinError(s+1,EnergyFluxUpErr);	
		fAugerEnergySpectrumE3->SetBinContent(s+1,EnergyFluxE3);
		fAugerEnergySpectrumE3->SetBinError(s+1,EnergyFluxUpErr * pow(pow(10,BinEnergy),3));
		fAugerEnergySpectrumEvents->SetBinContent(s+1,Nev_comb);
		fExposure.push_back(Exposure);

		cout<<"E="<< BinEnergy << " Phi=" << EnergyFlux << " Nev="<< Nev_comb <<"  Exposure="<< Exposure<<endl;
	}//close for energy bins spectrum


	in_spectrumdata.close();

}//close ReadSpectrumData()


void DataReader::ResetData(){

	//Reset all histograms
	for(int s=0;s<fNbins_Xmax;s++){
		fXmaxData_fit[s]->Reset();
		fEnergyData_fit[s]->Reset();
		
    for(int j=0;j<fNmassAtEarth;j++){
			fEnergyMC_fit[j][s]->Reset();
			fXmaxMC_fit[j][s]->Reset();
			fGenEnergyMC_fit[j][s]->Reset();
			fGenXmaxMC_fit[j][s]->Reset();
		}//end loop masses
	}//end loop energy bins

}//close ResetDataHisto()



void DataReader::ReadData(){

	cout<<"DataReader::ReadData(): Reading data file "<<fDataFileName.c_str()<<endl;

	TFile* inputFile= new TFile(fDataFileName.c_str());
  if ( inputFile->IsZombie() ) {
    cerr << " DataReader::ReadData() - error opening " << fDataFileName.c_str() << endl;
    throw std::runtime_error("error opening data file");
  }

	//read all data histograms from file
	for(int s=0;s<fNbins_Xmax;s++){
		fEnergyData_fit[s]= (TH1D*)inputFile->Get(Form("Energy_bin%d",s+1));
   	fXmaxData_fit[s]= (TH1D*)inputFile->Get(Form("Xmax_bin%d",s+1));
 	}//end loop energy bins
 
	
}//close DataReader::ReadData()


void DataReader::ReadMC(){

  cout<<"DataReader::ReadMC(): Reading MC data file "<<fMCFileName.c_str()<<endl;

	TFile* inputFile= new TFile(fMCFileName.c_str());

  if ( inputFile->IsZombie() ) {
    cerr << "DataReader::ReadMC() - error opening " << fMCFileName.c_str() << endl;
    throw std::runtime_error("error opening data file");
  }

	//read all histograms from file
	for(int j=0;j<fNmassAtEarth;j++){	
		for(int s=0;s<fNbins_Xmax;s++){
			fEnergyMC_fit[j][s]= (TH1D*)inputFile->Get(Form("EnergyMC%d_bin%d",j+1,s+1));
  		fXmaxMC_fit[j][s]= (TH1D*)inputFile->Get(Form("XmaxMC%d_bin%d",j+1,s+1));
			fGenEnergyMC_fit[j][s]= (TH1D*)inputFile->Get(Form("GenEnergyMC%d_bin%d",j+1,s+1));
  		fGenXmaxMC_fit[j][s]= (TH1D*)inputFile->Get(Form("GenXmaxMC%d_bin%d",j+1,s+1));      		
		}//end loop energy bins
	}//end loop masses
	

}//close DataReader::ReadMC()

void DataReader::Save(){
   
	cout<<"DataReader::Save(): Saving data to file "<<fOutputFileName.c_str()<<endl;

	TFile* fOutputFile= new TFile(fOutputFileName.c_str(),"RECREATE");
 	fOutputFile->cd();
	
 	for(int s=0;s<fNbins_Xmax;s++){
		fEnergyData_fit[s]->SetNameTitle(Form("Energy_bin%d",s+1),Form("Energy bin%d",s+1));
   	fEnergyData_fit[s]->SetDrawOption("ep");
   	fEnergyData_fit[s]->Write();

   	fXmaxData_fit[s]->SetNameTitle(Form("Xmax_bin%d",s+1),Form("Xmax bin%d",s+1));
   	fXmaxData_fit[s]->SetDrawOption("ep");
   	fXmaxData_fit[s]->Write();

    for(int j=0;j<fNmassAtEarth;j++){	
			fEnergyMC_fit[j][s]->SetNameTitle(Form("EnergyMC%d_bin%d",j+1,s+1),Form("EnergyMC%d_bin%d",j+1,s+1));
   		fEnergyMC_fit[j][s]->SetDrawOption("hist");
   		fEnergyMC_fit[j][s]->Write();

			fXmaxMC_fit[j][s]->SetNameTitle(Form("XmaxMC%d_bin%d",j+1,s+1),Form("XmaxMC%d_bin%d",j+1,s+1));
   		fXmaxMC_fit[j][s]->SetDrawOption("hist");
   		fXmaxMC_fit[j][s]->Write();

			fGenEnergyMC_fit[j][s]->SetNameTitle(Form("GenEnergyMC%d_bin%d",j+1,s+1),Form("GenEnergyMC%d_bin%d",j+1,s+1));
   		fGenEnergyMC_fit[j][s]->SetDrawOption("hist");
   		fGenEnergyMC_fit[j][s]->Write();

			fGenXmaxMC_fit[j][s]->SetNameTitle(Form("GenXmaxMC%d_bin%d",j+1,s+1),Form("GenXmaxMC%d_bin%d",j+1,s+1));
   		fGenXmaxMC_fit[j][s]->SetDrawOption("hist");
   		fGenXmaxMC_fit[j][s]->Write();
		}//end loop masses
 }//end loop energy bins
 
 fOutputFile->Close();
 
}//close Save()

}//close namespace

