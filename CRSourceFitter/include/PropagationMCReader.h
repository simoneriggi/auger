/**
* @file PropagationMCReader.h
* @class PropagationMCReader
* @brief Read the propagation MC data entries, create the transfer matrix, calculate the propagated spectra and composition.
* 
* @author S. Riggi
* @date 22/04/2010
*/
#ifndef _PROPAGATIONMCREADER_H_
#define _PROPAGATIONMCREADER_H_


#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TGraph.h>

#include <vector>
#include <string>

namespace CRSourceFitter_ns {

class PropagationMCReader : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		PropagationMCReader();
		/** 
		\brief Class destructor
 		*/
		virtual ~PropagationMCReader();
	
	public:
		/**
 		* \brief Set the output filename containing processed data
 		*/
		void SetOutputFileName(std::string fileName) { fOutputFileName = fileName;}
	
		//Set Propagation resampling info
		//### EXTRAGALACTIC
		/**
 		* \brief Set the cosmological redshift evolution of sources
 		*/
		void SetRedShiftEvolution(double value){ fRedShiftIndex = value;}
		/**
 		* \brief Set the injection index of sources
 		*/
		void SetInjectionIndex(double value){ fInjIndexAtSource = value;}
	
		/**
 		* \brief Set the maximum energy of sources
 		*/
		void SetMaxInjEnergyAtSource(double value){ fEmaxAtSource = value;}
		/**
 		* \brief Set the flux normalization of the extragalactic spectrum component
 		*/
		void SetExtraGalSpectrumNormalization(double value){ fNormExtraGalSpectrum = value;} 
		/**
 		* \brief Set the mass abundances at source
 		*/
		void SetAbundanceAtSource(std::vector<double> vect) { fAbundanceAtSource = vect;}
		/**
 		* \brief Set the transfer matrix from source to Earth
 		*/
		void SetTransferMatrix(TMatrixD* matrix) { fTransferMatrix = matrix;}
		/**
 		* \brief Set how the maximum injected energy is calculated (false=Z x Emax, true=A x Emax) 
 		*/
		void SetEmaxADependence(bool choice) { fIsEmaxADependent = choice;}
		/**
 		* \brief Set the source cutoff mode
 		*/
		void SetSourceCutoffMode(int value) { fSourceCutoffMode= value;}
		/**
 		* \brief Set the source cutoff shape
 		*/
		void SetSourceCutoffShape(double value) { fSourceCutoffShape= value;}
	

		//## GALACTIC
		/**
 		* \brief Set the injection index of galactic sources
 		*/
		void SetGalacticInjectionIndex(double value){ fInjIndexAtSource_Galactic = value;}
		/**
 		* \brief Set the maximum energy of galactic sources
 		*/
		void SetGalacticMaxInjEnergyAtSource(double value){ fEmaxAtSource_Galactic = value;}
		/**
 		* \brief Set the flux normalization of the galactic spectrum component
 		*/
		void SetGalSpectrumNormalization(double value){ fNormGalSpectrum = value;} 	
		/**
 		* \brief Set the galactic mass abundances at Earth
 		*/
		void SetGalacticComposition(std::vector<double> vect) { fGalacticAbundance = vect;}
		/**
 		* \brief Enable/disable the galactic composition in the fit (e.g. the total mass composition in a given energy energy bin is calculated by summing the extragalactic and galactic contribution)
 		*/
		void IncludeGalacticComposition(bool choice) { fIncludeGalacticComposition = choice;}

	
	
	
		//Getters
		/**
 		* \brief Get the source cutoff mode
 		*/
		int GetSourceCutoffMode() { return fSourceCutoffMode;}
		/**
 		* \brief Get the source cutoff shape
 		*/
		double GetSourceCutoffShape() { return fSourceCutoffShape;}
	

		//get prop spectra with same auger spectrum binning
		/**
 		* \brief Get list of histograms with propagated spectrum @ Earth for the primary mass groups (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpEnergySpectrumAtEarth() { return fExpEnergySpectrumAtEarth;}
		/**
 		* \brief Get total (summed up over all mass groups) propagated spectrum @ Earth (calculation with matrix method)
 		*/
		TH1D* GetTotExpEnergySpectrumAtEarth() { return fTotExpEnergySpectrumAtEarth;}
		/**
 		* \brief Get list of histograms with propagated spectrum @ Earth for the primary mass groups with energy binning of Xmax analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpEnergySpectrumAtEarth_Xmax() { return fExpEnergySpectrumAtEarth_Xmax;}
		/**
 		* \brief Get list of histograms with direct (i.e. p->p,He->He,...) propagated spectrum @ Earth for the primary mass groups with energy binning of Xmax analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpDirectEnergySpectrumAtEarth_Xmax() { return fExpEnergySpectrumAtEarth_Xmax_diag;}
		/**
 		* \brief Get total (summed up over all mass groups) propagated spectrum @ Earth with energy binning of Xmax analysis (calculation with matrix method)
 		*/
		TH1D* GetTotExpEnergySpectrumAtEarth_Xmax() { return fTotExpEnergySpectrumAtEarth_Xmax;}

		/**
 		* \brief Get list of histograms with propagated spectrum @ Earth for the primary mass groups with energy binning of spectrum analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpEnergySpectrumAtEarth_Spectrum() { return fExpEnergySpectrumAtEarth_Spectrum;}
		/**
 		* \brief Get list of histograms with direct (i.e. p->p,He->He,...) propagated spectrum @ Earth for the primary mass groups with energy binning of spectrum analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpDirectEnergySpectrumAtEarth_Spectrum() { return fExpEnergySpectrumAtEarth_Spectrum_diag;}
		/**
 		* \brief Get total (summed up over all mass groups) propagated spectrum @ Earth with energy binning of spectrum analysis (exact calculation with weight method)
 		*/
		TH1D* GetTotExpEnergySpectrumAtEarth_Spectrum() { return fTotExpEnergySpectrumAtEarth_Spectrum;}
		/**
 		* \brief Get list of histograms with propagated spectrum @ Earth multiplied by E^3 for the primary mass groups with energy binning of spectrum analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpEnergySpectrumE3AtEarth_Spectrum() { return fExpEnergySpectrumE3AtEarth_Spectrum;}
		/**
 		* \brief Get list of histograms with direct (i.e. p->p,He->He,...) propagated spectrum @ Earth multiplied by E^3 for the primary mass groups with energy binning of spectrum analysis (calculation with matrix method)
 		*/
		std::vector<TH1D*> GetExpDirectEnergySpectrumE3AtEarth_Spectrum() { return fExpEnergySpectrumE3AtEarth_Spectrum_diag;}
		/**
 		* \brief Get total (summed up over all mass groups) propagated spectrum @ Earth multiplied by E^3 with energy binning of spectrum analysis (exact calculation with weight method)
 		*/
		TH1D* GetTotExpEnergySpectrumE3AtEarth_Spectrum() { return fTotExpEnergySpectrumE3AtEarth_Spectrum;}
		/**
 		* \brief Get total (summed up over all mass groups) propagated galactic spectrum @ Earth with energy binning of spectrum analysis
 		*/
		TH1D* GetGalacticSpectrum_Spectrum() { return fGalacticSpectrumHisto_Spectrum;} 

		/**
 		* \brief Get mass fractions observed @ Earth in each energy bin with energy binning of spectrum analysis
 		*/
		std::vector< std::vector<double> > GetExpFractionsAtEarth_Spectrum() { return fFractionsAtEarth_Spectrum;}
		/**
 		* \brief Get mass fractions observed @ Earth in each energy bin with energy binning of Xmax analysis
 		*/
		std::vector< std::vector<double> > GetExpFractionsAtEarth_Xmax() { return fFractionsAtEarth_Xmax;}
		/**
 		* \brief Get direct (i.e. p->p,He->He,...) mass fractions observed @ Earth in each energy bin with energy binning of spectrum analysis
 		*/
		std::vector< std::vector<double> > GetExpDirectFractionsAtEarth_Spectrum() { return fFractionsAtEarth_Spectrum_diag;}
		/**
 		* \brief Get direct (i.e. p->p,He->He,...) mass fractions observed @ Earth in each energy bin with energy binning of Xmax analysis
 		*/
		std::vector< std::vector<double> > GetExpDirectFractionsAtEarth_Xmax() { return fFractionsAtEarth_Xmax_diag;}

		/**
 		* \brief Get the transfer matrix
 		*/
		TMatrixD* GetTransferMatrix() { return fTransferMatrix;}

		/**
 		* \brief Initialize data structures
 		*/
 		void Init();
		/**
 		* \brief Reset data structures
 		*/
 		void Reset();
	
		/**
 		* \brief Store data structures in output file
 		*/
  	void StoreHisto();
  	/**
 		* \brief Calculate the propagated spectrum @ Earth from the injected one with the matrix approach
 		*/
  	void CalculateExpectedSpectrumAtEarth();
		/**
 		* \brief Generate the galactic spectrum according to a power-law function
 		*/
  	void GenerateGalacticSpectrum();
		/**
 		* \brief Calculate the observed mass fractions @ Earth in each energy bin
 		*/
  	void CalculateMassFractions();
		/**
 		* \brief Store the diagonal part of the transfer matrix for direct (i.e. p->p,He->He,...) spectra calculation
 		*/
		void CalculateDiagonalMatrix();
		/**
 		* \brief Injection spectrum function definition: power-law * exponential cutoff
 		*/
  	static double InjSpectrumFcn(double* x, double* par);

	
		enum NuclearMassPrimary { eNeutronA = 1, eProtonA = 1, eHeliumA = 4, eLithiumA= 7, eBerylliumA= 9, eBoronA= 11, eCarbonA= 12, 
														eNitrogenA= 14, eOxygenA= 16, eFluorineA= 19, eNeonA= 20, eSodiumA= 23, eMagnesiumA= 24, 
														eAluminiumA= 27, eSiliconA= 28, ePhosphorusA= 31, eSulphurA= 32, eChlorineA= 35,
														eArgonA= 40, ePotassiumA= 39, eCalciumA= 40, eScandiumA= 45, eTitaniumA= 48, 
														eVanadiumA= 51, eChromiumA= 52, eManganeseA= 55, eIronA= 56 };

		enum NuclearChargePrimary { eNeutronZ= 0, eProtonZ = 1, eHeliumZ = 2, eLithiumZ= 3, eBerylliumZ= 4, eBoronZ= 5, eCarbonZ= 6, 
															eNitrogenZ= 7, eOxygenZ= 8, eFluorineZ= 9, eNeonZ= 10, eSodiumZ= 11, eMagnesiumZ= 12, 
															eAluminiumZ= 13, eSiliconZ= 14, ePhosphorusZ= 15, eSulphurZ= 16, eChlorineZ= 17,
															eArgonZ= 18, ePotassiumZ= 19, eCalciumZ= 20, eScandiumZ= 21, eTitaniumZ= 22, 
															eVanadiumZ= 23, eChromiumZ= 24, eManganeseZ= 25, eIronZ= 26 };

	private:

		TFile* fOutputFile;
		std::string fOutputFileName;
  	int fMatrixSize;
  
  	//Histograms
  	TH1D* fExpEnergySpectrumAtEarthHisto;
		TH1D* fExpEnergySpectrumE3AtEarthHisto;
		TH1D* fExpEnergySpectrumAtEarthDiagHisto;
  	TH1D* fTotExpEnergySpectrumAtEarth;
  	TH1D* fTotExpEnergySpectrumE3AtEarth;

  	TH1D* fExpEnergySpectrumAtEarthHisto_Xmax;
		TH1D* fExpEnergySpectrumAtEarthDiagHisto_Xmax;
  	TH1D* fTotExpEnergySpectrumAtEarth_Xmax;

	  TH1D* fTotEnergySpectrumAtEarth_Spectrum;
  	TH1D* fExpEnergySpectrumAtEarthHisto_Spectrum;
		TH1D* fExpEnergySpectrumAtEarthDiagHisto_Spectrum;
  	TH1D* fExpEnergySpectrumE3AtEarthHisto_Spectrum;
		TH1D* fExpEnergySpectrumE3AtEarthDiagHisto_Spectrum;
  	TH1D* fTotExpEnergySpectrumAtEarth_Spectrum;
  	TH1D* fTotExpEnergySpectrumE3AtEarth_Spectrum;
   
  	TH1D* fGalacticSpectrumHisto_Spectrum;
  	TH1D* fEnergySpectrumAtSourceHisto;
  	TH1D* fInjectedEnergySpectrumAtSourceHisto;

 	
  	TMatrixD* fTransferMatrix;
		TMatrixD* fTransferDiagMatrix;
		TMatrixD* fSpectrumAtEarthMatrix;

  	//vectors of histograms	
  	std::vector<TH1D*> fExpEnergySpectrumAtEarth;
		std::vector<TH1D*> fExpEnergySpectrumE3AtEarth;
		std::vector<TH1D*> fExpEnergySpectrumAtEarth_diag;
  	std::vector<TH1D*> fInjectedEnergySpectrumAtSource;
  	std::vector<TH1D*> fExpEnergySpectrumAtEarth_Xmax;
		std::vector<TH1D*> fExpEnergySpectrumAtEarth_Xmax_diag;
  	std::vector<TH1D*> fExpEnergySpectrumAtEarth_Spectrum;
		std::vector<TH1D*> fExpEnergySpectrumAtEarth_Spectrum_diag;
  	std::vector<TH1D*> fExpEnergySpectrumE3AtEarth_Spectrum;
		std::vector<TH1D*> fExpEnergySpectrumE3AtEarth_Spectrum_diag;

  	std::vector<TF1*> fGeneratedEnergySpectrum;
 		TF1* InjEnergySpectrum;


		//## EXTRAGALACTIC COMPONENT SPECTRUM
		double fRedShiftIndex;
		double fInjIndexAtSource;
		double fEmaxAtSource;
		double fNormExtraGalSpectrum;
	
		bool fIsEmaxADependent;
		bool fDoMatrixCalculation;
		static int fSourceCutoffMode;	
		static double fSourceCutoffShape;
	
		//## GALACTIC COMPONENT SPECTRUM
  	double fInjIndexAtSource_Galactic;
		double fEmaxAtSource_Galactic;
		double fNormGalSpectrum;
		bool fIncludeGalacticComposition;
		bool fCallDiagonalMatrix;
		bool fIsGroupingInCharge;

		std::vector<double> fAbundanceAtSource;
		std::vector<double> fGalacticAbundance;
	
  	std::vector< std::vector<double> > fFractionsAtEarth_Spectrum;
		std::vector< std::vector<double> > fFractionsAtEarth_Xmax;
		std::vector< std::vector<double> > fFractionsAtEarth_Spectrum_diag;
		std::vector< std::vector<double> > fFractionsAtEarth_Xmax_diag;
		std::vector< std::vector<double> > fGenEminNucleus;
		std::vector< std::vector<double> > fGenEmaxNucleus;

  	std::vector< std::vector<double> > fExpSpectrumAtEarth;
		std::vector< std::vector<double> > fExpSpectrumAtEarth_diag;
	  std::vector< std::vector<double> > fInjectedSpectrumAtSource;
	 
	ClassDef(PropagationMCReader,1)

};//close class


}//close namespace

#endif
