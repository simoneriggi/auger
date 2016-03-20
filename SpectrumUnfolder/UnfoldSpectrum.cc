#include <ForwardFolder.h>
#include <SpectrumUtils.h>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TApplication.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;

void Usage(){
	cout<<"*** USAGE ***"<<endl;
	cout<<"[EXE] --config=[CONFIG-FILE]"<<endl;
	cout<<"*************"<<endl;
}

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
  {(char*)0, (int)0, (int*)0, (int)0}
};

void SetStyle();

int main(int argc, char *argv[]){

	TApplication* app= new TApplication("App",&argc,argv);
	
	SetStyle();

	//====================================================
	//==         REC SPECTRUM READ
	//=====================================================
	cout<<"INFO: Reading rec spectrum from file..."<<endl;
	TFile* inputFile= new TFile("data/AugerInclinedSpectrum.root","READ");
	if(!inputFile || inputFile->IsZombie()){
		cerr<<"ERROR: Failed to read input file with rec spectrum!"<<endl;
		return -1;
	}
	TH1D* recSpectrum= (TH1D*)inputFile->Get("RecSpectrum");
	if(!recSpectrum){
		cerr<<"ERROR: Failed to read rec spectrum from input file!"<<endl;
		return -1;
	}

	
	//Set rec binning
	std::vector<double> RecBins;
	double binLowEdge= 0;
	for(int i=0;i<recSpectrum->GetNbinsX();i++){
		binLowEdge= recSpectrum->GetBinLowEdge(i+1);
		RecBins.push_back(binLowEdge);
	}//end loop bins
	RecBins.push_back( binLowEdge+recSpectrum->GetBinWidth( recSpectrum->GetNbinsX() ) );

	//Set true binning (it should be equal or larger than RecBins)
	//Assume here same rec binning for simplicity
	std::vector<double> TrueBins;
	TrueBins.assign(RecBins.begin(),RecBins.end());

	cout<<"INFO: RecBins(";
	for(unsigned int i=0;i<RecBins.size();i++){
		cout<<RecBins[i]<<",";
	}
	cout<<endl;
	

	
	double MinTrueBin= TrueBins[0];
	double MaxTrueBin= TrueBins[TrueBins.size()-1];
	double MinRecBin= RecBins[0];
	double MaxRecBin= RecBins[RecBins.size()-1];
	cout<<"INFO: MinTrueBin="<<MinTrueBin<<" MaxTrueBin="<<MaxTrueBin<<endl;
	


	/*
	double TrueBinsArray[]= {17.5,17.6,17.7,17.8,17.9,18,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4,20.5};
	//double RecBinsArray[]= {18,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0};	
	double RecBinsArray[]= {17.5,17.6,17.7,17.8,17.9,18,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.0,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.0,20.1,20.2,20.3,20.4,20.5};
	int NTrueBins= sizeof(TrueBinsArray)/sizeof(double);
	int NRecBins= sizeof(RecBinsArray)/sizeof(double);
	std::vector<double> TrueBins;
	std::vector<double> RecBins;
	TrueBins.assign(TrueBinsArray, TrueBinsArray + NTrueBins);
	RecBins.assign(RecBinsArray, RecBinsArray + NRecBins);

	double MinTrueBin= TrueBinsArray[0];
	double MaxTrueBin= TrueBinsArray[NTrueBins-1];
	double MinRecBin= RecBinsArray[0];
	double MaxRecBin= RecBinsArray[NRecBins-1];
	cout<<"NTrueBins="<<NTrueBins<<" NRecBins="<<NRecBins<<endl;
	*/

	
	//====================================================
	//==         SET RESPONSE MODEL PARS
	//=====================================================	
	double gamma1= 3.27;
	double gamma2= 2.68;
	double gamma3= 4.2;
	double LgEBreak= 18.61;
	double LgECutoff= 19.41; 
	Unfolder_ns::FitPar normPar= Unfolder_ns::FitPar("Norm",1);
	Unfolder_ns::FitPar gammaPar1= Unfolder_ns::FitPar("Gamma1",gamma1,2,5);
	Unfolder_ns::FitPar gammaPar2= Unfolder_ns::FitPar("Gamma2",gamma2,2,5);
	Unfolder_ns::FitPar gammaPar3= Unfolder_ns::FitPar("Gamma3",gamma3,2,20);
	Unfolder_ns::FitPar breakPar= Unfolder_ns::FitPar("Break",LgEBreak,18,19);
	Unfolder_ns::FitPar cutoffPar= Unfolder_ns::FitPar("Cutoff",LgECutoff,19,20);
	Unfolder_ns::BrokenPowerLawsPars brokenPowerLawPars= Unfolder_ns::BrokenPowerLawsPars(normPar,gammaPar1,gammaPar2,gammaPar3,breakPar,cutoffPar);
	Unfolder_ns::PowerLawPars powerLawPars= Unfolder_ns::PowerLawPars(normPar,gammaPar2);
	Unfolder_ns::ConstResoPars biasPars= Unfolder_ns::ConstResoPars(0);
	Unfolder_ns::ConstResoPars sigmaPars= Unfolder_ns::ConstResoPars(0.1);
	Unfolder_ns::ConstTriggerPars triggerPars= Unfolder_ns::ConstTriggerPars(1);
	
	
	//Fix spectrum pars?
	brokenPowerLawPars.FixPar("Gamma1",gamma1);
	brokenPowerLawPars.FixPar("Break",LgEBreak);
	

	//====================================================
	//==         BUILD RESPONSE MATRIX
	//=====================================================
	TH2D* responseMatrix= Unfolder_ns::SpectrumUtils::ComputeParametricResponse(brokenPowerLawPars,biasPars,sigmaPars,triggerPars,TrueBins,RecBins);
	

	//====================================================
	//==         RUN UNFOLDING
	//=====================================================
	bool useFitRange= true;
	double fitMin= 18.5;
	double fitMax= 20.1;
	bool computeUncertainties= true;//false;
	int nRandSamples= 100;
	Unfolder_ns::ForwardFolder* ff= new Unfolder_ns::ForwardFolder;
	int status= ff->RunUnfold(recSpectrum,responseMatrix,brokenPowerLawPars,computeUncertainties,nRandSamples,useFitRange,fitMin,fitMax);	
	if(status<0) {
		cerr<<"ERROR: Unfolding run failed!"<<endl;
		return -1;
	}
	TH1D* UnfoldedSpectrum= ff->GetUnfoldedSpectrum();
	TH1D* FFSpectrum= ff->GetForwardFoldedSpectrum();
		
	
	
	//====================================================
	//==         DRAW PLOTS
	//=====================================================
	//--> Response Matrix
	TCanvas* ResponseMatrixPlot= new TCanvas("ResponseMatrixPlot","ResponseMatrixPlot");
	ResponseMatrixPlot->cd();
	responseMatrix->SetStats(0);
	responseMatrix->Draw("COLZ");

	TF1* diagFcn= new TF1("diagFcn","x",MinTrueBin,MaxTrueBin);
	diagFcn->Draw("lsame");
	
	//--> Results Plot
	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();

	recSpectrum->SetTitle(0);
	recSpectrum->SetStats(0);
	recSpectrum->GetXaxis()->SetTitle("lg(E/eV)");
	recSpectrum->GetXaxis()->SetTitleSize(0.06);
	recSpectrum->GetXaxis()->SetTitleOffset(0.8);
	recSpectrum->GetYaxis()->SetTitle("nentries");
	recSpectrum->GetYaxis()->SetTitleSize(0.06);
	recSpectrum->GetYaxis()->SetTitleOffset(1.3);
	recSpectrum->SetMarkerColor(kBlack);
	recSpectrum->SetLineColor(kBlack);
	recSpectrum->SetMarkerStyle(8);
	recSpectrum->Draw("ep");
	
	UnfoldedSpectrum->SetMarkerColor(kRed);
	UnfoldedSpectrum->SetLineColor(kRed);
	UnfoldedSpectrum->SetMarkerStyle(21);
	UnfoldedSpectrum->Draw("ep same");
	
	FFSpectrum->SetMarkerColor(kGreen+1);
	FFSpectrum->SetLineColor(kGreen+1);
	FFSpectrum->SetMarkerStyle(23);
	FFSpectrum->Draw("ep same");

	TLegend* PlotLegend= new TLegend(0.6,0.6,0.8,0.8);	
	PlotLegend->SetTextSize(0.045);
	PlotLegend->SetTextFont(52);
	PlotLegend->AddEntry(recSpectrum,"rec","PL");
	PlotLegend->AddEntry(UnfoldedSpectrum,"unfolded","PL");
	PlotLegend->AddEntry(FFSpectrum,"fitted","PL");
	PlotLegend->Draw("same");


	//inputFile->Close();
	
	cout<<"INFO: End application run"<<endl;
	if(app) app->Run();
	
	return 0;

}//close main


void SetStyle(){

	TStyle* myStyle= new TStyle("myStyle","myStyle");

	//## CANVAS & PAD
	myStyle->SetCanvasDefH(700); 
  myStyle->SetCanvasDefW(700); 
	myStyle->SetFrameBorderMode(0);
	myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
	myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadBottomMargin(0.12);
  myStyle->SetPadLeftMargin(0.16);
  myStyle->SetPadRightMargin(0.1);

	//## TITLE
	myStyle->SetOptTitle(0);
	myStyle->SetTitleX(0.1f);
	myStyle->SetTitleW(0.8f);                 
  myStyle->SetTitleXOffset(0.8);
  myStyle->SetTitleYOffset(1.1);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(0);//border size of Title PavelLabel
	myStyle->SetTitleSize(0.06,"X");
  myStyle->SetTitleSize(0.06,"Y");
	myStyle->SetTitleSize(0.06,"Z");
	
	//## STAT
	myStyle->SetOptStat("eMR");
	//myStyle->SetOptStat(1);
  myStyle->SetStatColor(0);
	myStyle->SetStatY(0.975);                
  myStyle->SetStatX(0.95);                
  myStyle->SetStatW(0.35);//0.2                
  myStyle->SetStatH(0.10);//0.15
  myStyle->SetStatBorderSize(1);

	myStyle->SetTitleFont(52,"X");
  myStyle->SetTitleFont(52,"Y");
  myStyle->SetTitleFont(52,"Z");
  myStyle->SetLabelFont(42,"X");
  myStyle->SetLabelFont(42,"Y");
  myStyle->SetLabelFont(42,"Z");   
	
	//## OTHER
  myStyle->SetOptFit(1);
	myStyle->SetOptLogx(0);
	myStyle->SetOptLogy(0);
  //myStyle->SetPalette(1,0);
  myStyle->SetMarkerStyle(8);
  myStyle->SetMarkerSize(0.6);
  myStyle->SetFuncWidth(1.); 
  myStyle->SetErrorX(0.);

	myStyle->SetNumberContours(999);
	myStyle->SetPalette(55);

	gROOT->SetStyle("myStyle");
	gStyle= myStyle;
	myStyle->cd();

}//close SetStyle()



