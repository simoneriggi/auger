/**
* @file DataReader.h
* @class DataReader
* @brief Read the hybrid data entries and create histograms for the composition fit
* 
* @author S. Riggi
* @date 22/04/2010
*/

#ifndef _DATAREADER_H_
#define _DATAREADER_H_

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TF1.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <vector>

namespace CRSourceFitter_ns {

class DataReader : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		DataReader();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~DataReader();
 
	public:
		
		void Set();
		/**
 		* \brief Set the data filename
 		*/
		void SetDataFileName(std::string fileName) { fDataFileName = fileName;}
		/**
 		* \brief Set the MC filename
 		*/
		void SetMCFileName(std::string fileName) { fMCFileName = fileName;}	
		/**
 		* \brief Set the output filename where to save data
 		*/
		void SetOutputFileName(std::string fileName) { fOutputFileName = fileName;}
		/**
 		* \brief Set the ascii data table containing spectrum information
 		*/
		void SetSpectrumTableFileName(std::string fileName) { fSpectrumTableFileName = fileName;}
	
		/**
 		* \brief Get list of data energy histograms
 		*/
		std::vector<TH1D*> GetEnergyData() { return fEnergyData_fit;}
		/**
 		* \brief Get list of data Xmax histograms
 		*/
		std::vector<TH1D*> GetXmaxData() { return fXmaxData_fit;}

		/**
 		* \brief Get list of MC energy histograms
 		*/
		std::vector<TH1D*> GetEnergyMC(int mcComponent) { 
			if(mcComponent>=0 && mcComponent<fEnergyMC_fit.size()) 
				return fEnergyMC_fit[mcComponent];
			else{
				std::string errMsg = "DataReader::GetEnergyMC(): Error, requested index exceeds vector size...exit!";
  	  	throw std::runtime_error(errMsg);
			}	
		}
		std::vector<TH1D*> GetGenEnergyMC(int mcComponent) { 
			if(mcComponent>=0 && mcComponent<fGenEnergyMC_fit.size()) 
				return fGenEnergyMC_fit[mcComponent];
			else{
				std::string errMsg = "DataReader::GetGenEnergyMC(): Error, requested index exceeds vector size...exit!";
    		throw std::runtime_error(errMsg);
			}	
		}

		/**
 		* \brief Get list of MC Xmax histograms
 		*/
		std::vector<TH1D*> GetXmaxMC(int mcComponent) { 
			if(mcComponent>=0 && mcComponent<fXmaxMC_fit.size()) 
				return fXmaxMC_fit[mcComponent];
			else{
				std::string errMsg = "DataReader::GetXmaxMC(): Error, requested index exceeds vector size...exit!";
  	  	throw std::runtime_error(errMsg);
			}	
		}
		std::vector<TH1D*> GetGenXmaxMC(int mcComponent) { 
			if(mcComponent>=0 && mcComponent<fGenXmaxMC_fit.size()) 
				return fGenXmaxMC_fit[mcComponent];
			else{
				std::string errMsg = "DataReader::GetGenXmaxMC(): Error, requested index exceeds vector size...exit!";
  	  	throw std::runtime_error(errMsg);
			}	
		}
	
		/**
 		* \brief Get the Auger energy spectrum histogram
 		*/
		TH1D* GetAugerEnergySpectrum() { return fAugerEnergySpectrum;}
		/**
 		* \brief Get the Auger energy spectrum histogram multiplied by E^3
 		*/
		TH1D* GetAugerEnergySpectrumE3() { return fAugerEnergySpectrumE3;}
		/**
 		* \brief Get the Auger energy spectrum number of events
 		*/
		TH1D* GetAugerEnergySpectrumEvents() { return fAugerEnergySpectrumEvents;}
		/**
 		* \brief Get the Auger energy spectrum exposure
 		*/
		std::vector<double> GetExposure() { return fExposure;}
	
 	
		/**
 		* \brief Initialize the data structures 
 		*/
 		void Init();//allocate all histograms
		/**
 		* \brief Reset the data structures 
 		*/
 		void ResetData();//reset histograms
	
		/**
 		* \brief Read spectrum data table
 		*/
		void ReadSpectrumData();
		/**
 		* \brief Read data histograms from file
 		*/
		void ReadData();
		/**
 		* \brief Read MC data histograms from file
 		*/
		void ReadMC();
		/**
 		* \brief Store data structures in output ROOT file
 		*/
 		void Save();
 	
	private:

		/**
 		* \brief Filename of data histogram
 		*/
		std::string fDataFileName;
		/**
 		* \brief Filename of MC histogram
 		*/
		std::string fMCFileName;
		/**
 		* \brief Filename of output ROOT file with processed info
 		*/
		std::string fOutputFileName;
		/**
 		* \brief Filename of input ascii table with spectrum information
 		*/
		std::string fSpectrumTableFileName;


  	//### SPECTRUM VARS
		/**
 		* \brief Vector of energy spectrum exposure
 		*/
  	std::vector<double> fExposure;
		/**
 		* \brief Energy spectrum histogram
 		*/
  	TH1D* fAugerEnergySpectrum;
		/**
 		* \brief Energy spectrum histogram multiplied by E^3
 		*/
  	TH1D* fAugerEnergySpectrumE3;
		/**
 		* \brief Energy spectrum number of events
 		*/
  	TH1D* fAugerEnergySpectrumEvents;

	
  	//## DATA HISTOGRAMS
		/**
 		* \brief Histogram of energy data
 		*/
  	TH1D* fEnergyHistoData;
		/**
 		* \brief Histogram of Xmax data
 		*/
 		TH1D* fXmaxHistoData;
	
		/**
 		* \brief List of histograms of energy data
 		*/
		std::vector<TH1D*> fEnergyData_fit;
		/**
 		* \brief List of histograms of Xmax data
 		*/
 		std::vector<TH1D*> fXmaxData_fit;

		//## MC HISTOGRAMS
		/**
 		* \brief Histogram of Energy MC data
 		*/
  	TH1D* fEnergyHistoMC;
		TH1D* fGenEnergyHistoMC;
		/**
 		* \brief Histogram of Xmax MC data
 		*/
 		TH1D* fXmaxHistoMC;
		TH1D* fGenXmaxHistoMC;
		/**
 		* \brief List of histograms of Xmax MC data
 		*/
		std::vector< std::vector<TH1D*> > fXmaxMC_fit;
		std::vector< std::vector<TH1D*> > fGenXmaxMC_fit;
		/**
 		* \brief List of histograms of Energy MC data
 		*/
		std::vector< std::vector<TH1D*> > fEnergyMC_fit;
		std::vector< std::vector<TH1D*> > fGenEnergyMC_fit;


	ClassDef(DataReader,1)

};

}//close namespace 

#endif

