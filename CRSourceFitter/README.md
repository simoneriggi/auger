# CRSourceFitter
A tool to fit a cosmic ray source model scenario to joint spectrum and composition data

## Instructions
- Compile and run the program passing the run steerings in a config file   

  Example: make clean  
           make   
           SourceFitter --fit --config=myConfig.dat  

- Draw fitted results stored in a file:  

  Example: SourceFitter --draw --config=myConfig.dat --input=myFitResults.root  

- Xmax Data, spectrum tables and the propagation matrix are in data/ and tab/ directory.
  Xmax data and MC distributions are pre-calculated and stored in histograms.  
	Set the desired data paths in the config file.   
   
  When applying an energy shift to data, the energy binning changes so **YOU MUST CHANGE**   
  the energy binning in the include/AnalysisConsts.h.  
  (...would be better to set binning directly in config file...)  
  You must also set the corresponding data and MC, spectrum table and propagation matrix file, in the config file.  

  Example: I want to fit the data applying an energy shift of +20%  
					 --> put HybridData_ICRC09_EnergyShift20.root as DataFile in the config file  
					 --> put MCData_qgsjetII_EnergyShift20.root as MCFile in the config file  
					 --> put MatrixFile_SourceEv0_EnergySyst20.root as MatrixFile steering in the config file  
           --> put AugerSpectrumTable_ICRC09_EnergyShift20.dat as SpectrumTableFile in the config file  
           --> replace default AnalysisConsts.h to AnalysisConsts.h.ENERGYSHIFT20  
           make clean  
					 make  

- If you want to run many fits, i.e. make scans with (Gamma-Emax) fixed, you can use the  
  script script/GammaEmaxScanner.sh. This generates the needed config files and run scripts:  

  Example: GammaEmaxScanner.sh [runNo] [GammaStart] [GammaEnd] [GammaStep] [EmaxStart] [EmaxEnd] [EmaxStep]  
           GammaEmaxScanner.sh 1 1.00 4.00 0.05 18.00 21.00 0.05  
           This produces lots of run scripts and corresponding config files making a grid in Gamma-Emax  
           with step size 0.05, from Gamma=1 to Gamma=4 and from Emax=18 to Emax=21  
           
