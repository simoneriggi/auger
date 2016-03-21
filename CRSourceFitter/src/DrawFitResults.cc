/**
* @file DrawFitResults.cc
* @class DrawFitResults
* @brief Draw the results of the fit (fitted data and MC, elongation rate, rms, chi square, ...)
* 
* @author S. Riggi
* @date 22/04/2010
*/

#include "DrawFitResults.h"
#include "PropagationMCReader.h"
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
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TRandom.h>

#include<TFractionFitter.h>
#include<TPaveText.h>
#include<TVirtualFitter.h>
#include<TObjArray.h>
#include<TMatrixD.h>

#include <TMinuit.h>
#include <TApplication.h>


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


ClassImp(CRSourceFitter_ns::DrawFitResults)

namespace CRSourceFitter_ns {

DrawFitResults::DrawFitResults(){

  fNmass= 3;
	fXmaxShift= 0.;
	fApplication= NULL;
		
}//close constructor


DrawFitResults::~DrawFitResults(){
	
}//close destructor


void DrawFitResults::Init(){
	
	//open input file
	fInputFile= new TFile(fDataFileName);
  if ( fInputFile->IsZombie() ) {
    cerr << " DrawFitResults::Init() - error opening " << fDataFileName << endl;
    throw std::runtime_error("error opening data file");
  }

	//open output file
	fPlotFile= new TFile(fOutputFileName,"RECREATE");

	//init graphics options
	gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1110); 
	gStyle->SetOptTitle(0); 
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);

	myStyle= new TStyle("myStyle","left centered");
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 

	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(1);
  myStyle->SetStatColor(0);
  myStyle->SetStatBorderSize(1);
  myStyle->SetOptTitle(0);
  myStyle->SetOptStat(0);
  myStyle->SetOptFit(0);
	myStyle->SetOptLogx(1);
  myStyle->SetPalette(1,0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
  myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);
  myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.2);                
  myStyle->SetStatH(0.15);                
  myStyle->SetTitleXOffset(1);
  myStyle->SetTitleYOffset(1);
  //myStyle->SetMarkerStyle(20);
  //myStyle->SetMarkerSize(0.5);
  myStyle->SetFuncWidth(1.);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetPadRightMargin(0.05);
  myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
	myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");
	myStyle->SetErrorX(0.);
	myStyle->cd();

	//create the interactive application
	fApplication= new TApplication("Application", 0, 0);

}//close destructor


void DrawFitResults::Draw(){

	cout<<"DrawFitResults::Draw()"<<endl;
	Init();//open input file
	
	fPlotFile->cd();

	DrawSpectrumFit();
	DrawXmaxFit();
	DrawMassFractions();
	DrawElongationRate();
	DrawRMS();
		
	//interactive draw of all canvas
	if(fApplication) fApplication->Run();


	fPlotFile->Close();	

}//close Draw()


void DrawFitResults::DrawSpectrumFit(){
		
	cout<<"DrawFitResults::DrawSpectrumFit(): Drawing spectrum fit results"<<endl;
	
	myStyle->SetOptStat(0); 
	myStyle->SetOptTitle(0); 
	myStyle->SetOptLogy(1);
	myStyle->SetOptLogx(0);
  //myStyle->SetCanvasBorderMode(0);
  //myStyle->SetPadBorderMode(0);
	
	
 	TCanvas* SpectrumFitPlot=new TCanvas("SpectrumFitPlot","SpectrumFitPlot");
 	SpectrumFitPlot->cd();

 	TH2D* SpectrumFitPlotBkg= new TH2D("SpectrumFitPlotBkg","",100,18.0,20.5,100,1.e+22,1.e+25);
	SpectrumFitPlotBkg->GetXaxis()->SetTitle("log_{10}[E/eV]");
 	SpectrumFitPlotBkg->GetXaxis()->SetTitleSize(0.06);
	SpectrumFitPlotBkg->GetXaxis()->SetLabelSize(0.05);
 	SpectrumFitPlotBkg->GetXaxis()->SetTitleOffset(0.85);
 	SpectrumFitPlotBkg->GetYaxis()->SetTitle("E^{3}J(E) [m^{-2}s^{-1}sr^{-1}eV^{2}]");
 	SpectrumFitPlotBkg->GetYaxis()->SetTitleSize(0.06);
	SpectrumFitPlotBkg->GetYaxis()->SetLabelSize(0.05);
 	SpectrumFitPlotBkg->GetYaxis()->SetTitleOffset(1.2);
 	SpectrumFitPlotBkg->Draw();

	TH1D* AugerEnergySpectrumE3= (TH1D*)fInputFile->Get("AugerEnergySpectrumE3");
	AugerEnergySpectrumE3->SetMarkerStyle(8);
	AugerEnergySpectrumE3->SetMarkerColor(kBlack);
	AugerEnergySpectrumE3->SetLineColor(kBlack);
	AugerEnergySpectrumE3->SetMarkerSize(1.3);
 	AugerEnergySpectrumE3->Draw("ep same");

 
	TH1D* ExtragalSpectrumE3Fit= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Fit");
	ExtragalSpectrumE3Fit->SetLineColor(kBlack);
	ExtragalSpectrumE3Fit->Draw("l same");

 	TH1D* GalSpectrumE3Fit= (TH1D*)fInputFile->Get("GalSpectrumE3Fit");
 	GalSpectrumE3Fit->SetLineColor(kBlack);
 	GalSpectrumE3Fit->SetLineStyle(1);
	GalSpectrumE3Fit->SetLineStyle(kDashed);
 	GalSpectrumE3Fit->Draw("l same");
 
 	TH1D* TotSpectrumE3Fit= (TH1D*)fInputFile->Get("TotSpectrumE3Fit");
 	TotSpectrumE3Fit->SetLineWidth(2);
	TotSpectrumE3Fit->Draw("l same");

	
	/*	
	TH1D* ExtragalSpectrumE3ComponentVector[fNmassAtEarth];
	TH1D* ExtragalSpectrumE3SumComponentVector[fNmassAtSource];
	for(int i=0;i<fNmassAtSource;i++) {
		ExtragalSpectrumE3SumComponentVector[i]= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Component1");
		ExtragalSpectrumE3SumComponentVector[i]->Reset();
	}

	for(int i=0;i<fNmassAtEarth;i++){
		cout<<"Drawing component "<<i+1<<endl;

		ExtragalSpectrumE3ComponentVector[i]= (TH1D*)fInputFile->Get(Form("ExtragalSpectrumE3Component%d",i+1));
		int color;
		double A= fMinMassAtEarth[i] + 0.5*(fMaxMassAtEarth[i]-fMinMassAtEarth[i]);
		cout<<"A="<<A<<endl;
		if(A == PropagationMCReader::eProtonA){
			color= kRed;

			ExtragalSpectrumE3SumComponentVector[0]->Add(ExtragalSpectrumE3ComponentVector[i]]);
			ExtragalSpectrumE3SumComponentVector[0] ->SetLineColor(color);
		}
		else if(A == PropagationMCReader::eHeliumA){
			color= kCyan;

			ExtragalSpectrumE3SumComponentVector[1]->Add(ExtragalSpectrumE3ComponentVector[i]]);
			ExtragalSpectrumE3SumComponentVector[1] ->SetLineColor(color);
		}
		else if(A>PropagationMCReader::eHeliumA && A<=PropagationMCReader::eMagnesiumA){
			color= kGreen;

			ExtragalSpectrumE3SumComponentVector[2]->Add(ExtragalSpectrumE3ComponentVector[i]]);
			ExtragalSpectrumE3SumComponentVector[2] ->SetLineColor(color);
		}	
		else if(A>PropagationMCReader::eMagnesiumA){
			color= kBlue;
			
			ExtragalSpectrumE3SumComponentVector[3]->Add(ExtragalSpectrumE3ComponentVector[i]]);
			ExtragalSpectrumE3SumComponentVector[3] ->SetLineColor(color);
		}
						
		fExtragalSpectrumE2Fit[i]->SetLineColor(color);
		fExtragalSpectrumE2Fit[i]->Draw("l same"); 
	}

	
	fExtragalSpectrumE2Fit_SumComponents[0]->SetLineColor(kRed);
	fExtragalSpectrumE2Fit_SumComponents[0]->SetLineWidth(2);
	fExtragalSpectrumE2Fit_SumComponents[0]->Draw("l same");
	fExtragalSpectrumE2Fit_SumComponents[1]->SetLineColor(kCyan);
	fExtragalSpectrumE2Fit_SumComponents[1]->SetLineWidth(2);
	fExtragalSpectrumE2Fit_SumComponents[1]->Draw("l same");
	fExtragalSpectrumE2Fit_SumComponents[2]->SetLineColor(kGreen);	
	fExtragalSpectrumE2Fit_SumComponents[2]->SetLineWidth(2);	
	fExtragalSpectrumE2Fit_SumComponents[2]->Draw("l same");
	fExtragalSpectrumE2Fit_SumComponents[3]->SetLineColor(kBlue);
	fExtragalSpectrumE2Fit_SumComponents[3]->SetLineWidth(2);
	fExtragalSpectrumE2Fit_SumComponents[3]->Draw("l same");
	*/

	
 	TH1D* ExtragalSpectrumE3Component1= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Component1"); 
 	ExtragalSpectrumE3Component1->SetLineColor(kRed);
	ExtragalSpectrumE3Component1->Draw("l same");

	TH1D* ExtragalSpectrumE3Component2= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Component2"); 
 	ExtragalSpectrumE3Component2->SetLineColor(kCyan);
	ExtragalSpectrumE3Component2->Draw("l same");

  TH1D* ExtragalSpectrumE3Component3= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Component3"); 
 	ExtragalSpectrumE3Component3->SetLineColor(kGreen);
	ExtragalSpectrumE3Component3->Draw("l same");

  TH1D* ExtragalSpectrumE3Component4= (TH1D*)fInputFile->Get("ExtragalSpectrumE3Component4"); 
 	ExtragalSpectrumE3Component4->SetLineColor(kBlue);
	ExtragalSpectrumE3Component4->Draw("l same");
	
	//## Get direct proton component
	TH1D* ExtragalSpectrumE3DirectComponent1= (TH1D*)fInputFile->Get("ExtragalSpectrumE3DirectComponent1"); 
	ExtragalSpectrumE3DirectComponent1->SetLineColor(kRed);
	ExtragalSpectrumE3DirectComponent1->SetLineStyle(kDashed);
	//ExtragalSpectrumE3DirectComponent1->Draw("l same");
	
	cout<<"*** DIRECT PROTON FRACTION ***"<<endl;
	for(int s=0;s<ExtragalSpectrumE3DirectComponent1->GetNbinsX();s++){
		double thisDirectProtonFraction= ExtragalSpectrumE3DirectComponent1->GetBinContent(s+1)/TotSpectrumE3Fit->GetBinContent(s+1);
		cout<<"bin "<<s+1<<"  lgE="<<ExtragalSpectrumE3DirectComponent1->GetBinCenter(s+1)<<"  pFract="<<thisDirectProtonFraction<<endl;  
	}//end loop energy bins

	
	gPad->RedrawAxis();	
	

	TLegend* legend = new TLegend(0.6,0.5,0.7,0.6,"","brNDC");
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetTextSize(0.035);
	legend->AddEntry(AugerEnergySpectrumE3,"ICRC09","p");
	legend->AddEntry(TotSpectrumE3Fit,"LL fit","l");
	//legend->AddEntry(ExtragalSpectrumE3Fit,"EG","l");
	//legend->AddEntry(GalSpectrumE3Fit,"GAL","l");	
	legend->AddEntry(ExtragalSpectrumE3Component1,"EG - p","l");
	legend->AddEntry(ExtragalSpectrumE3Component2,"EG - He","l");
	legend->AddEntry(ExtragalSpectrumE3Component3,"EG - CNO","l");
	legend->AddEntry(ExtragalSpectrumE3Component4,"EG - Fe","l");
	legend->Draw("same");
 
	fPlotFile->cd();	
	SpectrumFitPlot->Write();

}//close DrawFitResults::DrawSpectrumFit()


void DrawFitResults::DrawXmaxFit(){

	cout<<"DrawFitResults::DrawXmaxFit(): Drawing Xmax fit results"<<endl;

	TCanvas* XmaxFitPlot;
	TH1D* data;
	TH1D* mc;
	TH1D* fit;
	TLegend* legend;
	
	
	for(int s=0;s<fNbins_Xmax;s++){
		//myStyle->SetOptStat(1110);
		myStyle->SetOptStat(0);  
		myStyle->SetOptTitle(0); 
		myStyle->SetOptLogy(0);
		myStyle->SetOptLogx(0);
  	myStyle->SetCanvasBorderMode(0);
  	myStyle->SetPadBorderMode(0);

		//single fit canvas
		XmaxFitPlot= new TCanvas(Form("XmaxFitPlot_%d",s+1),Form("XmaxFitPlot_%d",s+1));
		XmaxFitPlot->SetBorderMode(0);
  	XmaxFitPlot->SetFrameFillColor(0);
  	XmaxFitPlot->SetFillColor(0);
  	XmaxFitPlot->SetFrameBorderMode(0);
		

		//Create Legenda
		legend= new TLegend(0.6,0.5,0.7,0.75,"","brNDC");
  	legend->SetFillColor(0);
		legend->SetBorderSize(0);
  	legend->SetTextSize(0.04);
		

		//Get and draw data
		data= (TH1D*)fInputFile->Get(Form("XmaxData_bin%d",s+1));
		data->SetMarkerSize(1.1);
		data->SetMarkerStyle(8);
		data->SetMarkerColor(kBlack);
		data->SetLineColor(kBlack);
		data->GetXaxis()->SetTitleFont(52);
		data->GetXaxis()->SetLabelFont(42);
		data->GetXaxis()->SetTitleSize(0.06);
  	data->GetXaxis()->SetTitleOffset(0.7);
  	data->GetXaxis()->SetLabelSize(0.05);
		data->GetXaxis()->SetTitle("X_{max} [g/cm^{2}]");
		data->GetYaxis()->SetTitle("entries");
		data->GetYaxis()->SetTitleFont(52);
		data->GetYaxis()->SetLabelFont(42);
		data->GetYaxis()->SetTitleSize(0.06);
		data->GetYaxis()->SetTitleOffset(1.2);
		data->GetYaxis()->SetLabelSize(0.05);
		data->SetTitle(0);
		data->Draw("ep");
		legend->AddEntry(data,"Data","p");	
		
		//get and draw fit
		fit= (TH1D*)fInputFile->Get(Form("XmaxFit_bin%d",s+1));
		fit->SetLineColor(kBlack);
		legend->AddEntry(fit,"LL fit","l");
		
		XmaxFitPlot->cd();
		fit->Draw("hist same");
	
		//get and draw mass components
		for(int i=0;i<fNmassAtEarth;i++){
			mc= (TH1D*)fInputFile->Get(Form("XmaxMC%d_bin%d_scaled",i+1,s+1));
			mc->SetName(Form("XmaxMC%d_bin%d_scaled",i+1,s+1));

			if(i==0) {
				mc->SetLineColor(kRed);
				legend->AddEntry(mc,"MC Proton","l");	
			}
			else if(i==1) {
				mc->SetLineColor(kCyan);
				legend->AddEntry(mc,"MC Helium","l");
			}
			else if(i==2) {
				mc->SetLineColor(kGreen);
				legend->AddEntry(mc,"MC Oxygen","l");	
			}
			else if(i==3) {
				mc->SetLineColor(kBlue);
				legend->AddEntry(mc,"MC Iron","l");
			}
			else {
				mc->SetLineColor(kBlack);			
				legend->AddEntry(mc,"MC XXX","l");
			}
			
			XmaxFitPlot->cd();
			mc->Draw("hist same");	
		}//close loop mass components
		
		
		XmaxFitPlot->cd();
		gPad->RedrawAxis();		
		
		
  	legend->Draw("same");
		
		fPlotFile->cd();
		XmaxFitPlot->Write();		
	}//close loop energy bins

	
}//close DrawFitResults::DrawXmaxFit()


void DrawFitResults::DrawMassFractions(){

	cout<<"DrawFitResults::DrawMassFractions(): Drawing fitted fractions"<<endl;

	myStyle->SetOptStat(0); 
	myStyle->SetOptTitle(0); 
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
	myStyle->SetOptLogx(1);
	myStyle->SetOptLogy(0);

	TCanvas* FitFractionPlot= new TCanvas("FitFractionPlot","FitFractionPlot");
	FitFractionPlot->cd(); 	
	FitFractionPlot->SetBorderMode(0);
  FitFractionPlot->SetFrameFillColor(0);
  FitFractionPlot->SetFillColor(0);
  FitFractionPlot->SetFrameBorderMode(0);
	FitFractionPlot->cd();

	TH2D * backGround = new TH2D("backMassPlot","",100,pow(10,17.6),pow(10,19.8),100,0.,1.);
  backGround->GetXaxis()->SetTitle("E/eV");
	backGround->GetXaxis()->SetTitleSize(0.05);
	backGround->GetXaxis()->SetTitleOffset(0.7);
  backGround->GetYaxis()->SetTitle("Mass Fractions");
	backGround->GetYaxis()->SetTitleSize(0.05);
  backGround->Draw();

	TLegend* legend = new TLegend(0.6,0.5,0.7,0.8,"","brNDC");
  legend->SetFillColor(0);
	legend->SetBorderSize(0);
  legend->SetTextSize(0.03);
  
	TGraphAsymmErrors* MassFractDetGraph;
	TGraphAsymmErrors* MassFractEarthGraph;


	for(int i=0;i<fNmassAtEarth;i++){

	 	MassFractDetGraph= (TGraphAsymmErrors*)fInputFile->Get(Form("MassDetFraction%d",i+1));
		MassFractEarthGraph=(TGraphAsymmErrors*)fInputFile->Get(Form("MassFraction%d",i+1));		

		if(i==0){
			cout<<"MassFractDetGraph i="<<i<<endl;
			MassFractDetGraph->SetMarkerColor(kRed+2);
  		MassFractDetGraph->SetLineColor(kRed+2);
  		MassFractDetGraph->SetMarkerSize(1.5);
  		MassFractDetGraph->SetMarkerStyle(8);
			//legend->AddEntry(MassFractDetGraph,"p @ det","p");

			MassFractEarthGraph->SetMarkerColor(kRed);
  		MassFractEarthGraph->SetLineColor(kRed);
  		MassFractEarthGraph->SetMarkerSize(1.5);
  		MassFractEarthGraph->SetMarkerStyle(8);
  		MassFractEarthGraph->Draw("P");
			legend->AddEntry(MassFractEarthGraph,"p","p");
		}
		else if(i==1){
			MassFractDetGraph->SetMarkerColor(kCyan+2);
  		MassFractDetGraph->SetLineColor(kCyan+2);
			//legend->AddEntry(MassFractDetGraph,"He @ det","p");

			MassFractEarthGraph->SetMarkerColor(kCyan);
  		MassFractEarthGraph->SetLineColor(kCyan);
			legend->AddEntry(MassFractEarthGraph,"He","p");
		}	
		else if(i==2){
			MassFractDetGraph->SetMarkerColor(kGreen+2);
  		MassFractDetGraph->SetLineColor(kGreen+2);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(29);
			//legend->AddEntry(MassFractDetGraph,"O @ det","p");

			MassFractEarthGraph->SetMarkerColor(kGreen);
  		MassFractEarthGraph->SetLineColor(kGreen);
  		MassFractEarthGraph->SetMarkerSize(2);
  		MassFractEarthGraph->SetMarkerStyle(29);
  		MassFractEarthGraph->Draw("P");
			legend->AddEntry(MassFractEarthGraph,"O","p");
		}
		else if(i==3){
			MassFractDetGraph->SetMarkerColor(kBlue+2);
  		MassFractDetGraph->SetLineColor(kBlue+2);
			//legend->AddEntry(MassFractDetGraph,"Fe @ det","p");

			MassFractEarthGraph->SetMarkerColor(kBlue);
  		MassFractEarthGraph->SetLineColor(kBlue);
			legend->AddEntry(MassFractEarthGraph,"Fe","p");
		}
			
		else {
			MassFractDetGraph->SetMarkerColor(kOrange+2);
  		MassFractDetGraph->SetLineColor(kOrange+2);
  		MassFractDetGraph->SetMarkerSize(2);
  		MassFractDetGraph->SetMarkerStyle(29);
  		
			legend->AddEntry(MassFractDetGraph,"XXX @ det","p");

			MassFractEarthGraph->SetMarkerColor(kOrange);
  		MassFractEarthGraph->SetLineColor(kOrange);
  		MassFractEarthGraph->SetMarkerSize(2);
  		MassFractEarthGraph->SetMarkerStyle(29);
  		
			legend->AddEntry(MassFractEarthGraph,"XXX","p");
		}

		MassFractDetGraph->SetMarkerSize(2);
  	MassFractDetGraph->SetMarkerStyle(21);
  	//MassFractDetGraph->Draw("P");
			
  	MassFractEarthGraph->SetMarkerSize(2);
  	MassFractEarthGraph->SetMarkerStyle(29);
  	MassFractEarthGraph->Draw("P");
		
		legend->Draw("same");
  
	}//end loop masses at Earth

	fPlotFile->cd();
	FitFractionPlot->Write();


}//close DrawFitResults::DrawMassFractions()


TF1* DrawFitResults::defineElong(double* p, int color, int style) {

  TF1* elong = new TF1("elong","[3]+[0]+[1]*(18.-log10(x))+[2]*(18.-log10(x))*(18.-log10(x))",1.e17,1.e20);
  elong->SetParameters(p); 
  elong->SetLineColor(color);
  elong->SetLineStyle(style); 
  elong->SetLineWidth(2);
  elong->SetNpx(500); 
  
  return elong;

}//close DrawFitResults::defineElong()


void DrawFitResults::DrawElongationRate(){
	
	cout<<"DrawFitResults::DrawElongationRate(): Drawing fitted elongation rate"<<endl;

	myStyle->SetOptStat(0); 
	myStyle->SetOptTitle(0); 
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
	myStyle->SetOptLogx(1);
	myStyle->SetOptLogy(0);

	//Draw ER
	TCanvas* ERPlot= new TCanvas("ERPlot","ERPlot");
	ERPlot->cd();	
	ERPlot->SetBorderMode(0);
  ERPlot->SetFrameFillColor(0);
  ERPlot->SetFillColor(0);
  ERPlot->SetFrameBorderMode(0);
 
  double epos100[4] = {735.816,-62.6457,1.51601, 0};   // eponEnergy_100.root 
  double epos5600[4]= {634.933,-59.3276,-0.476072, 0}; // eponEnergy_5600.root 
  double qgII100[4] = {737.823,-47.9116,-0.989034, 0}; // qgsIInEnergy_100.root 
  double qgII5600[4]= {651.467,-56.0661,-2.67465, 0};  // qgsIInEnergy_5600.root 
  double qg01100[4] = {728.141,-51.6506,-1.24245, 0};  // qgsnEnergy_100.root 
  double qg015600[4]= {634.926,-59.2864,-1.44251, 0};  // qgsnEnergy_5600.root 
  double sibl100[4] = {740.219,-56.3145,1.2876, 0};    // sibnEnergy_100.root
  double sibl5600[4]= {644.855,-58.6166,-0.935736, 0}; // sibnEnergy_5600.root 

	//assign shift as last parameter
	epos100[3]= fXmaxShift;
	epos5600[3]= fXmaxShift;
	qgII100[3]= fXmaxShift;
	qgII5600[3]= fXmaxShift;
	qg01100[3]= fXmaxShift;
	qg015600[3]= fXmaxShift;
  sibl100[3]= fXmaxShift;
	sibl5600[3]= fXmaxShift;

  const int protonColor=kRed;
  const int ironColor=kBlue;

	//define background
  TH2F * backGround = new TH2F("backGround","",100,pow(10,17.6),pow(10,19.8),100,520.,880.);
  backGround->GetXaxis()->SetTitle("E [eV]");
	backGround->GetXaxis()->SetTitleSize(0.06);
	backGround->GetXaxis()->SetLabelSize(0.05);
	backGround->GetXaxis()->SetTitleOffset(0.7);
  backGround->GetYaxis()->SetTitle("<X_{max}> [g/cm^{2}]");	
	backGround->GetYaxis()->SetTitleSize(0.06);
	backGround->GetYaxis()->SetLabelSize(0.05);
	backGround->GetYaxis()->SetTitleOffset(1.2);
  backGround->Draw();

	//define MC reference lines
	/*
  TF1* ERpredict_pSibyll= defineElong(sibl100,protonColor,2);
	ERpredict_pSibyll->Draw("same");
  TF1* ERpredict_FeSibyll= defineElong(sibl5600,ironColor,2);
	ERpredict_FeSibyll->Draw("same");
  TF1* ERpredict_pEPOS= defineElong(epos100,protonColor,3);
	ERpredict_pEPOS->Draw("same");
  TF1* ERpredict_FeEPOS= defineElong(epos5600,ironColor,3);
	ERpredict_FeEPOS->Draw("same");
  TF1* ERpredict_pQGSJETII= defineElong(qgII100,protonColor,4);
	ERpredict_pQGSJETII->Draw("same");
  TF1* ERpredict_FeQGSJETII= defineElong(qgII5600,ironColor,4);
	ERpredict_FeQGSJETII->Draw("same");
	TF1* ERpredict_pQGSJET01= defineElong(qg01100,protonColor,5);
	ERpredict_pQGSJET01->Draw("same");
  TF1* ERpredict_FeQGSJET01= defineElong(qg015600,ironColor,5);
	ERpredict_FeQGSJET01->Draw("same");
	*/

	//add QGSJET01 prediction
	const int N_pQGSJET01= 47;
	double Energy_pQGSJET01[]={5.341e+17,5.9e+17,6.321e+17,6.618e+17,7.144e+17,7.952e+17,8.851e+17,1.031e+18,1.139e+18,1.23e+18,1.39e+18,1.583e+18,1.803e+18,2.022e+18,2.217e+18,2.393e+18,2.603e+18,2.875e+18,3.128e+18,3.481e+18,4.183e+18,4.516e+18,5.262e+18,6.42e+18,6.722e+18,7.312e+18,7.833e+18,8.392e+18,9.778e+18,1.105e+19,1.268e+19,1.359e+19,1.524e+19,1.683e+19,1.803e+19,2.007e+19,2.411e+19,2.789e+19,3.176e+19,3.376e+19,4.056e+19,4.346e+19,4.8e+19,5.383e+19,6.37e+19,7.654e+19,8.074e+19};
	double MeanXmax_pQGSJET01[]={709.2,711.9,713.3,714.2,716.5,718.8,721.1,724.3,726.6,728.4,731.2,734.4,737.2,739.9,741.7,743.6,745.4,747.2,749.1,751.8,755.5,757.8,761,765.1,766.5,768.3,769.7,771.6,774.3,777.1,780.3,781.7,783.9,786.2,787.6,789.9,794,797.2,799.1,800.9,804.6,806,808.3,810.6,814.2,817.9,818.8};

	const int N_FeQGSJET01= 48;
	double Energy_FeQGSJET01[]={7.259e+17,7.777e+17,8.142e+17,8.591e+17,9.134e+17,1.017e+18,1.185e+18,1.36e+18,1.537e+18,1.876e+18,2.202e+18,2.452e+18,2.647e+18,2.901e+18,3.156e+18,3.54e+18,4.063e+18,4.386e+18,4.808e+18,5.352e+18,5.957e+18,6.431e+18,6.89e+18,7.553e+18,8.091e+18,9.075e+18,9.649e+18,1.066e+19,1.168e+19,1.261e+19,1.47e+19,1.611e+19,1.739e+19,1.921e+19,2.239e+19,2.649e+19,3.017e+19,3.257e+19,3.626e+19,3.945e+19,4.194e+19,4.492e+19,5.235e+19,5.695e+19,6.291e+19,6.948e+19,7.559e+19,7.854e+19};
	double MeanXmax_FeQGSJET01[]={621.1,622.5,623.9,625.7,627.5,630.3,633.5,637.2,640.4,645.9,649.5,652.8,654.1,656.4,658.7,661.9,664.7,667.4,669.3,672.5,674.8,676.6,678.9,680.7,682.6,685.3,687.2,689,691.3,693.6,697.2,699.5,701.4,703.7,706.9,711,714.2,716.1,719.3,721.1,722.5,723.9,727.5,729.4,732.1,733.9,736.2,737.2};

	TGraph* ERpredict_pQGSJET01= new TGraph(N_pQGSJET01,Energy_pQGSJET01,MeanXmax_pQGSJET01);
	ERpredict_pQGSJET01->SetLineColor(kRed);
	ERpredict_pQGSJET01->SetLineStyle(1);
	ERpredict_pQGSJET01->SetLineWidth(2);
	ERpredict_pQGSJET01->Draw("l same");

	TGraph* ERpredict_FeQGSJET01= new TGraph(N_FeQGSJET01,Energy_FeQGSJET01,MeanXmax_FeQGSJET01);
	ERpredict_FeQGSJET01->SetLineColor(kBlue);
	ERpredict_FeQGSJET01->SetLineStyle(1);
	ERpredict_FeQGSJET01->SetLineWidth(2);
	ERpredict_FeQGSJET01->Draw("l same");

	//add QGSJETII prediction
	const int N_pQGSJETII= 39;
	double Energy_pQGSJETII[]={5.27e+17,5.603e+17,6.002e+17,6.381e+17,6.836e+17,7.156e+17,7.965e+17,9.141e+17,1.033e+18,1.177e+18,1.33e+18,1.526e+18,2.237e+18,2.567e+18,2.924e+18,3.356e+18,3.763e+18,4.352e+18,4.919e+18,5.559e+18,6.428e+18,8.596e+18,9.865e+18,1.106e+19,1.27e+19,1.468e+19,1.647e+19,1.861e+19,2.136e+19,2.451e+19,2.749e+19,3.154e+19,3.592e+19,4.029e+19,4.695e+19,5.346e+19,5.951e+19,6.881e+19,7.896e+19};
	double MeanXmax_pQGSJETII[]={721.1,722.5,724.3,725.2,726.6,727.5,729.4,732.6,735.3,738.5,740.8,743.6,751.8,755,757.3,760.6,762.8,765.6,768.3,770.6,773.9,779.8,782.1,784.9,787.2,789.9,792.7,795,797.2,800.5,802.8,805,807.3,809.6,812.8,815.1,817.4,820.2,822.5};

	const int N_FeQGSJETII= 40;
  double Energy_FeQGSJETII[]={5.237e+17,6.01e+17,6.898e+17,7.737e+17,8.812e+17,1.004e+18,1.143e+18,1.302e+18,1.494e+18,1.689e+18,1.924e+18,2.158e+18,2.495e+18,2.863e+18,3.212e+18,3.686e+18,4.072e+18,4.745e+18,5.322e+18,6.155e+18,7.339e+18,8.487e+18,9.447e+18,1.084e+19,1.235e+19,1.417e+19,1.614e+19,1.824e+19,2.078e+19,2.384e+19,2.674e+19,3.092e+19,3.364e+19,3.95e+19,4.499e+19,5.163e+19,5.879e+19,6.645e+19,7.567e+19,7.923e+19};
	double MeanXmax_FeQGSJETII[]={629.4,633,636.2,639,642.2,645.9,648.6,651.4,655.5,658.3,661.5,664.2,667.4,670.6,673.4,676.6,679.4,682.6,684.9,688.5,692.7,695.9,698.6,701.8,704.1,707.3,710.6,712.8,716.5,718.8,722,724.8,726.6,730.3,733,735.8,738.1,741.3,743.6,745};

	TGraph* ERpredict_pQGSJETII= new TGraph(N_pQGSJETII,Energy_pQGSJETII,MeanXmax_pQGSJETII);
	ERpredict_pQGSJETII->SetLineColor(kRed);
	ERpredict_pQGSJETII->SetLineStyle(7);
	ERpredict_pQGSJETII->SetLineWidth(2);
	ERpredict_pQGSJETII->Draw("l same");

	TGraph* ERpredict_FeQGSJETII= new TGraph(N_FeQGSJETII,Energy_FeQGSJETII,MeanXmax_FeQGSJETII);
	ERpredict_FeQGSJETII->SetLineColor(kBlue);
	ERpredict_FeQGSJETII->SetLineStyle(7);
	ERpredict_FeQGSJETII->SetLineWidth(2);
	ERpredict_FeQGSJETII->Draw("l same");

	//add Sibyll2.1 prediction
	const int N_pSibyll= 41;
	double Energy_pSibyll[]={5.392e+17,6.095e+17,6.43e+17,7.268e+17,8.214e+17,8.868e+17,9.428e+17,1.066e+18,1.205e+18,1.372e+18,1.551e+18,1.766e+18,2.011e+18,2.273e+18,2.589e+18,2.949e+18,3.333e+18,3.826e+18,4.324e+18,4.963e+18,5.567e+18,6.897e+18,7.736e+18,8.879e+18,1.004e+19,1.143e+19,1.302e+19,1.472e+19,1.676e+19,1.909e+19,2.141e+19,2.439e+19,2.778e+19,3.164e+19,3.576e+19,4.073e+19,4.603e+19,5.243e+19,6.349e+19,7.287e+19,7.927e+19};
	double MeanXmax_pSibyll[]={721.1,724.3,725.7,728.9,731.7,733.5,735.3,738.1,741.7,745,747.7,751.4,754.6,757.3,761,763.8,767,770.6,773.4,777.1,780.7,785.8,788.5,791.7,795,798.2,801.4,805,807.8,811.5,814.7,817.9,821.1,824.3,828,830.7,833.9,837.2,842.7,845.9,847.7};

	const int N_FeSibyll= 39;
	double Energy_FeSibyll[]={5.518e+17,6.334e+17,7.104e+17,8.153e+17,9.216e+17,1.05e+18,1.186e+18,1.341e+18,1.539e+18,1.74e+18,1.981e+18,2.257e+18,2.57e+18,2.883e+18,3.309e+18,3.74e+18,4.227e+18,4.814e+18,5.484e+18,6.246e+18,7.169e+18,8.289e+18,9.228e+18,1.059e+19,1.197e+19,1.353e+19,1.541e+19,1.755e+19,1.984e+19,2.26e+19,2.573e+19,3.092e+19,3.549e+19,4.012e+19,4.604e+19,5.203e+19,5.882e+19,6.699e+19,7.747e+19};
  double MeanXmax_FeSibyll[]={622.9,627.1,629.4,632.6,636.2,639.9,642.7,646.3,649.5,652.8,656,658.7,662.4,665.6,668.3,671.6,675.2,678,681.7,684.9,688.5,691.3,695,698.2,701.4,704.6,707.8,711,714.2,717.4,720.2,725.2,728.4,731.7,734.9,737.6,741.3,744.5,747.2};

	TGraph* ERpredict_pSibyll= new TGraph(N_pSibyll,Energy_pSibyll,MeanXmax_pSibyll);
	ERpredict_pSibyll->SetLineColor(kRed);
	ERpredict_pSibyll->SetLineStyle(2);
	ERpredict_pSibyll->SetLineWidth(2);
	ERpredict_pSibyll->Draw("l same");

	TGraph* ERpredict_FeSibyll= new TGraph(N_FeSibyll,Energy_FeSibyll,MeanXmax_FeSibyll);
	ERpredict_FeSibyll->SetLineColor(kBlue);
	ERpredict_FeSibyll->SetLineStyle(2);
	ERpredict_FeSibyll->SetLineWidth(2);
	ERpredict_FeSibyll->Draw("l same");


	//add EPOS 1.99 prediction
	const int N_pEPOS199= 58;
	double Energy_pEPOS199[]={5.392e+17,6.046e+17,6.676e+17,6.988e+17,7.776e+17,8.718e+17,8.988e+17,1.008e+18,1.156e+18,1.286e+18,1.442e+18,1.487e+18,1.654e+18,1.854e+18,1.926e+18,2.144e+18,2.385e+18,2.478e+18,2.757e+18,3.067e+18,3.186e+18,3.545e+18,3.945e+18,4.098e+18,4.56e+18,5.073e+18,5.311e+18,5.909e+18,6.281e+18,6.426e+18,6.625e+18,7.315e+18,8.078e+18,8.392e+18,9.337e+18,1.039e+19,1.087e+19,1.21e+19,1.357e+19,1.399e+19,1.544e+19,1.745e+19,1.785e+19,2.001e+19,2.244e+19,2.331e+19,2.554e+19,2.842e+19,2.953e+19,3.31e+19,3.683e+19,3.856e+19,4.701e+19,4.847e+19,5.475e+19,6.046e+19,6.233e+19,7.042e+19};
	double MeanXmax_pEPOS199[]={732.6,735.3,738.5,739.9,742.2,745,745.9,749.1,752.8,755.5,758.3,759.6,761.9,765.1,766.1,768.8,771.6,772.9,775.2,778.4,778.9,782.1,784.9,786.2,788.5,791.7,792.7,795.4,797.2,798.2,798.6,801.4,804.1,805,808.3,811,811.9,814.7,817.9,818.8,821.6,824.8,825.7,828.4,831.2,832.6,835.3,838.1,839,841.7,845,845.4,850.9,852.3,855.5,858.3,859.2,861.9};

	const int N_FeEPOS199= 60;
	double Energy_FeEPOS199[]={5.23e+17,6e+17,6.935e+17,7.716e+17,8.651e+17,8.988e+17,1e+18,1.113e+18,1.156e+18,1.276e+18,1.42e+18,1.475e+18,1.629e+18,1.84e+18,1.897e+18,2.111e+18,2.349e+18,2.422e+18,2.715e+18,2.998e+18,3.138e+18,3.519e+18,3.655e+18,4.005e+18,4.491e+18,4.997e+18,5.191e+18,5.775e+18,6.329e+18,6.575e+18,7.26e+18,8.139e+18,8.328e+18,9.337e+18,1.039e+19,1.087e+19,1.201e+19,1.336e+19,1.399e+19,1.544e+19,1.718e+19,1.799e+19,1.986e+19,2.21e+19,2.313e+19,2.842e+19,2.975e+19,3.31e+19,3.655e+19,3.826e+19,4.257e+19,4.701e+19,4.81e+19,4.959e+19,5.475e+19,6.046e+19,6.281e+19,7.042e+19,7.6e+19,7.835e+19};
	double MeanXmax_FeEPOS199[]={621.6,625.7,629.8,633,636.2,637.2,640.4,643.6,644.5,647.2,650.5,651.4,654.1,657.8,658.7,661.5,664.2,665.1,668.8,670.6,672,675.7,676.6,679.4,682.1,684.9,686.7,689,691.7,692.7,695.4,698.6,699.5,702.3,705.5,706.4,709.2,711.9,713.3,716.5,719.3,720.2,722.9,725.7,727.1,732.6,733.9,736.7,739.4,740.8,743.1,745.9,746.8,747.2,750,753.2,754.1,756.9,759.6,760.1};

	TGraph* ERpredict_pEPOS199= new TGraph(N_pEPOS199,Energy_pEPOS199,MeanXmax_pEPOS199);
	ERpredict_pEPOS199->SetLineColor(kRed);
	ERpredict_pEPOS199->SetLineStyle(4);	
	ERpredict_pEPOS199->SetLineWidth(2);
	ERpredict_pEPOS199->Draw("l same");

	TGraph* ERpredict_FeEPOS199= new TGraph(N_FeEPOS199,Energy_FeEPOS199,MeanXmax_FeEPOS199);
	ERpredict_FeEPOS199->SetLineColor(kBlue);
	ERpredict_FeEPOS199->SetLineStyle(4);
	ERpredict_FeEPOS199->SetLineWidth(2);
	ERpredict_FeEPOS199->Draw("l same");

 	//define Auger ER 
	const int N= 13;
	const double ERSystUncertainty= 13.;
  double ER_ICRC07[N]={688.5,699.2,703.6,710.1,719.2,729.5,729.8,734.4,738.7,751.1,767.3,771.8,762.6};
  double ERErr_ICRC07[N]={2.8,3.0,2.8,3.2,3.3,3.7,3.6,3.3,3.7,4.1,4.9,7.8,11.8};
  double Energy_ICRC07[N]={17.85,17.95,18.05,18.15,18.25,18.35,18.44,18.59,18.79,19.00,19.19,19.39,19.62};

	double ER_ICRC09[N]={714.5,725.0,734.7,735.7,745.8,737.3,745.7,741.3,750.7,756.0,755.5,759.0,765.8};
  double ERErr_ICRC09[N]={2.3,2.6,2.6,2.7,3.0,3.2,3.7,4.0,4.6,4.8,3.4,4.0,5.5};
  double Energy_ICRC09[N]={18.051,18.150,18.248,18.349,18.449,18.549,18.650,18.748,18.850,18.951,19.095,19.278,19.541};
	

  TGraphErrors* ElongationRateICRC07=new TGraphErrors();
	TGraphErrors* ElongationRateICRC09=new TGraphErrors();
  
  for(int i=0;i<N;i++){
    ElongationRateICRC07->SetPoint(i,pow(10,Energy_ICRC07[i]),ER_ICRC07[i]);
    ElongationRateICRC07->SetPointError(i,0.,ERErr_ICRC07[i]);

		ElongationRateICRC09->SetPoint(i,pow(10,Energy_ICRC09[i]),ER_ICRC09[i]);
    ElongationRateICRC09->SetPointError(i,0.,ERErr_ICRC09[i]);
  }
  
  ElongationRateICRC07->SetMarkerStyle(26);
  ElongationRateICRC07->SetMarkerSize(1.5);
  ElongationRateICRC07->SetMarkerColor(1);
  ElongationRateICRC07->SetLineColor(1);
  //ElongationRateICRC07->Draw("P");

	//ElongationRateICRC09->SetMarkerStyle(26);
	ElongationRateICRC09->SetMarkerStyle(8);
  //ElongationRateICRC09->SetMarkerSize(1.3);
	ElongationRateICRC09->SetMarkerSize(1.5);
  ElongationRateICRC09->SetMarkerColor(1);
  ElongationRateICRC09->SetLineColor(1);
  ElongationRateICRC09->Draw("P");

	//copy graph and set systematics bar
	TGraphErrors* ElongationRate_syst=(TGraphErrors*)ElongationRateICRC09->Clone("ElongationRate_syst");
	for(int i=0;i<N;i++){
		ElongationRate_syst->SetPointError(i,0.,ERSystUncertainty);
	}
	ElongationRate_syst->SetMarkerSize(1.5);	
	ElongationRate_syst->Draw(">");
	
	
	//get and draw fitted ER
	TGraphErrors* ElongationRateFit= (TGraphErrors*)fInputFile->Get("ERTrue");
	ElongationRateFit->SetMarkerStyle(8);
  ElongationRateFit->SetMarkerSize(1.3);
  ElongationRateFit->SetMarkerColor(kGreen);
  ElongationRateFit->SetLineColor(kGreen);
	ElongationRateFit->SetLineStyle(1);
	ElongationRateFit->SetLineWidth(2);
	ElongationRateFit->SetFillColor(kGreen);
	ElongationRateFit->SetFillStyle(3001);
  //ElongationRateFit->Draw("3l");
	//ElongationRateFit->Draw("l");

	TGraphErrors* ElongationRateFit_noErr= (TGraphErrors*)fInputFile->Get("ERTrue");
	for(int i=0;i<ElongationRateFit_noErr->GetN();i++){
		ElongationRateFit_noErr->SetPointError(i,0.,0.);
	}
	ElongationRateFit_noErr->SetMarkerStyle(8);
  ElongationRateFit_noErr->SetMarkerSize(1.3);
  ElongationRateFit_noErr->SetMarkerColor(kGreen);
  ElongationRateFit_noErr->SetLineColor(kGreen);
	ElongationRateFit_noErr->SetLineStyle(1);
	ElongationRateFit_noErr->SetLineWidth(2);
	ElongationRateFit_noErr->SetFillColor(kGreen);
	ElongationRateFit_noErr->SetFillStyle(3001);
  //ElongationRateFit_noErr->Draw("3l");
	ElongationRateFit_noErr->Draw("l");


	TLine* lineStyle_QGSJET01=new TLine(0,1,0,1);
	lineStyle_QGSJET01->SetLineStyle(1);
	lineStyle_QGSJET01->SetLineWidth(2);

	TLine* lineStyle_QGSJETII=new TLine(0,1,0,1);
	lineStyle_QGSJETII->SetLineStyle(7);
	lineStyle_QGSJETII->SetLineWidth(2);

	TLine* lineStyle_Sibyll=new TLine(0,1,0,1);
	lineStyle_Sibyll->SetLineStyle(2);
	lineStyle_Sibyll->SetLineWidth(2);

	TLine* lineStyle_EPOS199=new TLine(0,1,0,1);
	lineStyle_EPOS199->SetLineStyle(4);
	lineStyle_EPOS199->SetLineWidth(2);
  
	TLegend* legendRef = new TLegend(0.6,0.5,0.7,0.75,"","brNDC");
  legendRef->SetFillColor(0);
  legendRef->SetTextSize(0.04);
	legendRef->SetBorderSize(0);
	legendRef->AddEntry(lineStyle_QGSJET01,"QGSJET01","l");
	legendRef->AddEntry(lineStyle_QGSJETII,"QGSJETII","l");
	legendRef->AddEntry(lineStyle_Sibyll,"Sibyll 2.1","l");
	legendRef->AddEntry(lineStyle_EPOS199,"EPOS 1.99","l");
  legendRef->Draw("same");

  TLegend* legend = new TLegend(0.6,0.5,0.7,0.6,"","brNDC");
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(ElongationRateICRC09,"PRL 2010","p");
  //legend->AddEntry(ElongationRateICRC07,"ICRC07","p");
  //legend->AddEntry(ElongationRateFit,"LL fit","3l");
	legend->AddEntry(ElongationRateFit_noErr,"LL fit","l");
	//legend->AddEntry(ERGal_noErr,"GAL inferred","l");
  legend->Draw("same");
	
	
	fPlotFile->cd();
	ERPlot->Write();

}//close DrawFitResults::DrawElongationRate()


TF1* DrawFitResults::defineRMSref(double* p, int color, int style) {

  TF1* RMSref = new TF1("RMSref","[0]+[1]*(18.-log10(x))+[2]*(18.-log10(x))*(18.-log10(x))",1.e17,1.e20);
  RMSref->SetParameters(p); 
  RMSref->SetLineColor(color);
  RMSref->SetLineStyle(style); 
  RMSref->SetLineWidth(2);
  RMSref->SetNpx(500); 

  return RMSref;
}

void DrawFitResults::DrawRMS(){
	
	cout<<"DrawFitResults::DrawRMSRate(): Drawing fitted Xmax RMS"<<endl;

	myStyle->SetOptStat(0); 
	myStyle->SetOptTitle(0); 
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
	myStyle->SetOptLogx(1);
	myStyle->SetOptLogy(0);

	TCanvas* RMSPlot= new TCanvas("RMSPlot","RMSPlot");
	RMSPlot->cd();
	RMSPlot->SetBorderMode(0);
  RMSPlot->SetFrameFillColor(0);
  RMSPlot->SetFillColor(0);
  RMSPlot->SetFrameBorderMode(0);
	

	//define background
	TH2F * backGround = new TH2F("backGround2","",100,pow(10,17.8),pow(10,19.8),100,10.,80.);
  backGround->GetXaxis()->SetTitle("E [eV]");
	backGround->GetXaxis()->SetTitleSize(0.06);
	backGround->GetXaxis()->SetLabelSize(0.05);
	backGround->GetXaxis()->SetTitleOffset(0.7);
  backGround->GetYaxis()->SetTitle("RMS(X_{max}) [g/cm^{2}]");	
	backGround->GetYaxis()->SetTitleSize(0.06);
	backGround->GetYaxis()->SetLabelSize(0.05);
	backGround->GetYaxis()->SetTitleOffset(1.2);
  backGround->Draw();

	//add QGSJET01 prediction
	const int N_pQGSJET01= 17;
	double Energy_pQGSJET01[]={5.155e+17,6.051e+17,7.102e+17,8.529e+17,9.784e+17,1.338e+18,1.72e+18,2.098e+18,3.191e+18,4.637e+18,6.388e+18,9.643e+18,2.052e+19,3.368e+19,4.433e+19,6.743e+19,7.734e+19};
	double RMS_pQGSJET01[]={67.04,66.67,66.54,66.17,66.05,65.56,65.19,64.94,64.44,63.83,63.58,62.96,62.1,61.73,61.48,60.99,60.99};

	const int N_FeQGSJET01= 23;
	double Energy_FeQGSJET01[]={5.293e+17,6.119e+17,8.558e+17,1.068e+18,1.384e+18,1.726e+18,2.306e+18,2.77e+18,3.508e+18,4.41e+18,5.586e+18,6.759e+18,8.627e+18,1.052e+19,1.637e+19,2.027e+19,2.452e+19,2.967e+19,3.404e+19,4.797e+19,5.939e+19,6.971e+19,7.815e+19};
	double RMS_FeQGSJET01[]={23.09,22.84,22.72,22.47,22.22,22.1,21.98,21.85,21.6,21.6,21.36,21.23,21.11,20.99,20.74,20.62,20.62,20.49,20.37,20.25,20.12,20,20};


	TGraph* RMSpredict_pQGSJET01= new TGraph(N_pQGSJET01,Energy_pQGSJET01,RMS_pQGSJET01);
	RMSpredict_pQGSJET01->SetLineColor(kRed);
	RMSpredict_pQGSJET01->SetLineStyle(1);
	RMSpredict_pQGSJET01->SetLineWidth(2);
	RMSpredict_pQGSJET01->Draw("l same");

	TGraph* RMSpredict_FeQGSJET01= new TGraph(N_FeQGSJET01,Energy_FeQGSJET01,RMS_FeQGSJET01);
	RMSpredict_FeQGSJET01->SetLineColor(kBlue);
	RMSpredict_FeQGSJET01->SetLineStyle(1);
	RMSpredict_FeQGSJET01->SetLineWidth(2);
	RMSpredict_FeQGSJET01->Draw("l same");

	//add QGSJETII prediction
	const int N_pQGSJETII= 20;
	double Energy_pQGSJETII[]={5.242e+17,6.013e+17,6.951e+17,9.148e+17,1.195e+18,1.56e+18,1.763e+18,2.641e+18,3.03e+18,3.957e+18,5.207e+18,8.356e+18,1.091e+19,1.414e+19,1.861e+19,2.788e+19,3.669e+19,4.828e+19,6.258e+19,7.748e+19};
	double RMS_pQGSJETII[]={61.48,61.23,60.99,60.49,60.12,59.63,59.26,58.77,58.52,58.15,57.78,56.91,56.67,56.42,55.93,55.43,55.06,54.81,54.57,54.2};

	const int N_FeQGSJETII= 21;
	double Energy_FeQGSJETII[]={5.173e+17,6.025e+17,7.989e+17,1.035e+18,1.362e+18,1.779e+18,2.341e+18,3.057e+18,3.962e+18,4.58e+18,5.294e+18,6.309e+18,8.366e+18,1.101e+19,1.438e+19,1.892e+19,2.489e+19,3.301e+19,4.245e+19,5.629e+19,7.407e+19};
	double RMS_FeQGSJETII[]={23.83,23.7,23.7,23.46,23.33,23.09,22.96,22.84,22.72,22.59,22.47,22.35,22.22,22.1,21.98,21.85,21.6,21.6,21.48,21.23,21.11};


	TGraph* RMSpredict_pQGSJETII= new TGraph(N_pQGSJETII,Energy_pQGSJETII,RMS_pQGSJETII);
	RMSpredict_pQGSJETII->SetLineColor(kRed);
	RMSpredict_pQGSJETII->SetLineStyle(7);
	RMSpredict_pQGSJETII->SetLineWidth(2);
	RMSpredict_pQGSJETII->Draw("l same");

	TGraph* RMSpredict_FeQGSJETII= new TGraph(N_FeQGSJETII,Energy_FeQGSJETII,RMS_FeQGSJETII);
	RMSpredict_FeQGSJETII->SetLineColor(kBlue);
	RMSpredict_FeQGSJETII->SetLineStyle(7);
	RMSpredict_FeQGSJETII->SetLineWidth(2);
	RMSpredict_FeQGSJETII->Draw("l same");

	
	//add Sibyll prediction
	const int N_pSibyll= 21;
	double Energy_pSibyll[]={5.088e+17,6.252e+17,7.623e+17,9.956e+17,1.134e+18,1.31e+18,1.598e+18,2.087e+18,2.544e+18,3.126e+18,3.871e+18,4.648e+18,7.071e+18,8.622e+18,1.312e+19,1.807e+19,2.547e+19,3.326e+19,4.343e+19,5.759e+19,7.464e+19};
	double RMS_pSibyll[]={58.77,58.27,57.9,57.41,57.04,56.79,56.42,55.8,55.43,55.19,54.69,54.44,53.7,53.46,52.84,52.47,52.1,51.73,51.6,51.23,50.99};

	const int N_FeSibyll= 20;
	double Energy_FeSibyll[]={5.171e+17,5.887e+17,6.305e+17,7.178e+17,9.446e+17,1.331e+18,2.12e+18,2.585e+18,3.401e+18,4.211e+18,5.845e+18,8.756e+18,1.322e+19,1.713e+19,2.396e+19,2.586e+19,4.41e+19,5.059e+19,5.759e+19,7.695e+19};
	double RMS_FeSibyll[]={24.94,24.94,24.81,24.69,24.44,24.2,23.83,23.7,23.33,23.21,23.09,22.72,22.35,22.1,21.85,21.85,21.36,21.23,21.11,20.99};


	TGraph* RMSpredict_pSibyll= new TGraph(N_pSibyll,Energy_pSibyll,RMS_pSibyll);
	RMSpredict_pSibyll->SetLineColor(kRed);
	RMSpredict_pSibyll->SetLineStyle(2);
	RMSpredict_pSibyll->SetLineWidth(2);
	RMSpredict_pSibyll->Draw("l same");

	TGraph* RMSpredict_FeSibyll= new TGraph(N_FeSibyll,Energy_FeSibyll,RMS_FeSibyll);
	RMSpredict_FeSibyll->SetLineColor(kBlue);
	RMSpredict_FeSibyll->SetLineStyle(2);
	RMSpredict_FeSibyll->SetLineWidth(2);
	RMSpredict_FeSibyll->Draw("l same");


	//add EPOS1.99 prediction
	const int N_pEPOS199= 19;
	double Energy_pEPOS199[]={5.163e+17,6.742e+17,1.177e+18,1.548e+18,1.749e+18,2.302e+18,3.006e+18,3.956e+18,5.166e+18,7.281e+18,9.508e+18,1.251e+19,1.621e+19,2.133e+19,2.807e+19,3.638e+19,4.786e+19,6.25e+19,7.561e+19};
	double RMS_pEPOS199[]={61.48,61.11,60.12,59.75,59.51,59.26,59.01,58.89,58.77,58.64,58.77,58.64,58.77,58.77,58.89,59.14,59.51,59.75,60};

	const int N_FeEPOS199= 18;
	double Energy_FeEPOS199[]={5.299e+17,7.878e+17,1.037e+18,1.354e+18,1.781e+18,2.309e+18,3.968e+18,5.221e+18,7.358e+18,9.536e+18,1.264e+19,1.639e+19,2.156e+19,2.816e+19,3.705e+19,4.875e+19,6.318e+19,7.882e+19};
	double RMS_FeEPOS199[]={19.01,18.77,18.4,18.27,18.02,17.9,17.53,17.28,17.16,17.16,17.04,16.91,17.04,16.79,16.67,16.79,16.67,16.67};


	TGraph* RMSpredict_pEPOS199= new TGraph(N_pEPOS199,Energy_pEPOS199,RMS_pEPOS199);
	RMSpredict_pEPOS199->SetLineColor(kRed);
	RMSpredict_pEPOS199->SetLineStyle(2);
	RMSpredict_pEPOS199->SetLineWidth(2);
	RMSpredict_pEPOS199->Draw("l same");

	TGraph* RMSpredict_FeEPOS199= new TGraph(N_FeEPOS199,Energy_FeEPOS199,RMS_FeEPOS199);
	RMSpredict_FeEPOS199->SetLineColor(kBlue);
	RMSpredict_FeEPOS199->SetLineStyle(2);
	RMSpredict_FeEPOS199->SetLineWidth(2);
	RMSpredict_FeEPOS199->Draw("l same");

	//Define Auger RMS
	const int N= 13;
	const double RMSSystUncertainty= 6;
	double RMS_ICRC09[N]= {55.2,57.7,56.1,53.7,53.7,49.2,48.0,43.9,48.7,42.9,35.5,28.9,26.2};
	double RMSErr_ICRC09[N]= {3.0,3.4,3.2,3.6,4.0,3.7,4.8,3.8,5.1,6.1,3.4,4.5,10.0};
 	double Energy_ICRC09[N]={18.051,18.150,18.248,18.349,18.449,18.549,18.650,18.748,18.850,18.951,19.095,19.278,19.541};

	TGraphErrors* RMSICRC09= new TGraphErrors();

	for(int i=0;i<N;i++){
		RMSICRC09->SetPoint(i,pow(10,Energy_ICRC09[i]),RMS_ICRC09[i]);
		RMSICRC09->SetPointError(i,0.,RMSErr_ICRC09[i]);
	}//close for s
	
	//RMSICRC09->SetMarkerStyle(26);
	RMSICRC09->SetMarkerStyle(8);
  RMSICRC09->SetMarkerSize(1.5);
	//RMSICRC09->SetMarkerSize(1.1);
  RMSICRC09->SetMarkerColor(kBlack);
  RMSICRC09->SetLineColor(kBlack);
  RMSICRC09->Draw("P");

	TGraphErrors* RMSICRC09_syst= (TGraphErrors*)RMSICRC09->Clone("RMSICRC09_syst");
	for(int i=0;i<N;i++){
		RMSICRC09_syst->SetPointError(i,0.,RMSSystUncertainty);	
	}
	RMSICRC09_syst->SetMarkerSize(1.5);	
	RMSICRC09_syst->Draw(">");

	TGraphErrors* RMSFit= (TGraphErrors*)fInputFile->Get("RMSTrueFit");
	RMSFit->SetMarkerStyle(8);
  RMSFit->SetMarkerSize(1.3);
  RMSFit->SetMarkerColor(kGreen);
  RMSFit->SetLineColor(kGreen);
	RMSFit->SetLineStyle(1);
	RMSFit->SetLineWidth(2);
	RMSFit->SetFillColor(kGreen);
	RMSFit->SetFillStyle(3001);
  //RMSFit->Draw("3l");
	//RMSFit->Draw("l");
  
	TGraphErrors* RMSFit_noErr= (TGraphErrors*)RMSFit->Clone("RMSFit_noErr");
	for(int i=0;i<RMSFit_noErr->GetN();i++){
		RMSFit_noErr->SetPointError(i,0.,0.);	
	}
	RMSFit_noErr->SetMarkerStyle(8);
  RMSFit_noErr->SetMarkerSize(1.3);
  RMSFit_noErr->SetMarkerColor(kGreen);
  RMSFit_noErr->SetLineColor(kGreen);
	RMSFit_noErr->SetLineStyle(1);
	RMSFit_noErr->SetLineWidth(2);
	RMSFit_noErr->SetFillColor(kGreen);
	RMSFit_noErr->SetFillStyle(3001);
	RMSFit_noErr->Draw("l");

	TLine* lineStyle_QGSJET01=new TLine(0,1,0,1);
	lineStyle_QGSJET01->SetLineStyle(1);
	lineStyle_QGSJET01->SetLineWidth(2);

	TLine* lineStyle_QGSJETII=new TLine(0,1,0,1);
	lineStyle_QGSJETII->SetLineStyle(7);
	lineStyle_QGSJETII->SetLineWidth(2);

	TLine* lineStyle_Sibyll=new TLine(0,1,0,1);
	lineStyle_Sibyll->SetLineStyle(2);
	lineStyle_Sibyll->SetLineWidth(2);

	TLine* lineStyle_EPOS199=new TLine(0,1,0,1);
	lineStyle_EPOS199->SetLineStyle(4);
	lineStyle_EPOS199->SetLineWidth(2);

	TLegend* legendRef = new TLegend(0.6,0.5,0.7,0.75,"","brNDC");
  legendRef->SetFillColor(0);
  legendRef->SetTextSize(0.04);
	legendRef->SetBorderSize(0);

	legendRef->AddEntry(lineStyle_QGSJET01,"QGSJET01","l");
	legendRef->AddEntry(lineStyle_QGSJETII,"QGSJETII","l");
	legendRef->AddEntry(lineStyle_Sibyll,"Sibyll 2.1","l");
	legendRef->AddEntry(lineStyle_EPOS199,"EPOS 1.99","l");
	//legendRef->Draw("same");

  TLegend* rmsLegend = new TLegend(0.6,0.5,0.7,0.6,"","brNDC");
  rmsLegend->SetFillColor(0);
  rmsLegend->SetTextSize(0.04);
	rmsLegend->SetBorderSize(0);
  rmsLegend->AddEntry(RMSICRC09,"PRL 2010","p");
  //rmsLegend->AddEntry(RMSFit,"LL fit","l"); 
	rmsLegend->AddEntry(RMSFit_noErr,"LL fit","l");    
	//rmsLegend->AddEntry(RMSGal_noErr,"GAL inferred","l"); 
  rmsLegend->Draw("same");

	fPlotFile->cd();
  RMSPlot->Write();
  
  
}//close DrawFitResults::DrawRMS()

}//close namespace


