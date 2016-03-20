/**
* @file MathUtils.h
* @class MathUtils
* @brief Math utilities
* 
* @author S. Riggi
* @date 20/09/2011
*/

#ifndef _MATHUTILS_H_
#define _MATHUTILS_H_

#include <TH1D.h>
#include <TH2D.h>

#include <map>

namespace Unfolder_ns {



//===================================================
//==        MATH UTILS
//===================================================
class MathUtils : public TObject {

	public:

		/** 
		\brief Class constructor
 		*/
		MathUtils();
		/** 
		\brief Class destructor
 		*/
		virtual ~MathUtils();
		

	public:
		//Spectrum models
		static double PowerLawSpectrum(double* x, double* par);	
		static double BrokenPowerLawSpectrum(double* x, double* par);
		static double SmoothCutoffPowerLawSpectrum(double* x, double* par);
		static double GetPowerLawIntegral(double gamma,double lgEMin, double lgEMax);
		static double GetBrokenPowerLawIntegral(double Gamma,double Gamma2,double Gamma3,double Break,double Cutoff,double lgEMin, double lgEMax);
		
		//Response models
		static double TriggerEfficiencyModel(double* x, double* par);
		static double ResponseModel(double* x, double* par);
		static double ResolutionModel(double* x, double* par);
		static double BiasModel(double* x, double* par);

		
	ClassDef(MathUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class MathUtils+; 
#endif




//=================================================
//===             BIAS/RESO MODEL 
//==================================================
enum ResoModel {
	eConst= 0,
	eLinear= 1,
	ePol2= 2
};


class ResoPars : public TObject {
	public:
		virtual ~ResoPars() {};
	public:
		//Pure virtual
		virtual double GetA() const = 0;	
		virtual double GetB() const = 0;	
		virtual double GetC() const = 0;	
		virtual std::vector<double> GetPars() const = 0;
		
		//Standard
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return nPars;}
		static int GetParNumber() {return 0;}
		
	protected:
		static int nPars;
		int model;
		double A;
		double B;
		double C;
	ClassDef(ResoPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class ResoPars+; 
#endif

class ConstResoPars : public ResoPars {
	public:		
		ConstResoPars(double m_A)	{ 
			A= m_A;
			model= eConst;
			nPars= 1;	
		};
		virtual ~ConstResoPars() {};
	public:	
		static int GetParNumber() {return 1;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return 0;}
		virtual double GetC() const {return 0;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			return pars;
		}
	protected:
		//static int nPars;

	ClassDef(ConstResoPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class ConstResoPars+; 
#endif

class LinearResoPars : public ResoPars {
	public:		
		LinearResoPars(double m_A,double m_B)	{
			A= m_A;
			B= m_B;
			model= eLinear;
			nPars= 2;	
		};
		virtual ~LinearResoPars() {};
	public:
		static int GetParNumber() {return 2;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return 0;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			return pars;
		}
	protected:
		//static int nPars;
		
	ClassDef(LinearResoPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class LinearResoPars+; 
#endif

class Pol2ResoPars : public ResoPars {
	public:		
		Pol2ResoPars(double m_A,double m_B,double m_C){ 
			A= m_A; 
			B= m_B; 
			C= m_C;
			model= ePol2;		
			nPars= 3;
		};
		virtual ~Pol2ResoPars() {};
	public:
		static int GetParNumber() {return 3;}
		virtual double GetA() const {return A;}
		virtual double GetB() const {return B;}
		virtual double GetC() const {return C;}
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(A);
			pars.push_back(B);
			pars.push_back(C);
			return pars;
		}
	protected:
		//static int nPars;
		
	ClassDef(Pol2ResoPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class Pol2ResoPars+; 
#endif
//==================================================


//=================================================
//===             TRIGGER MODEL 
//==================================================
enum TriggerModel {
	eConstTrigger= 0,
	eSigmoidTrigger= 1
};
class TriggerPars : public TObject {
	public:
		virtual ~TriggerPars() {};
	public:
		//Pure virtual
		virtual double GetNorm() const = 0;
		virtual double GetBreak() const = 0;
		virtual double GetSmooth() const = 0;
		virtual std::vector<double> GetPars() const = 0;

		//Default
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return nPars;}
		static int GetParNumber() {return nPars;}

	protected:
		int model;
		static int nPars;
		double Norm;
		double Break;
		double Smooth;

	ClassDef(TriggerPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class TriggerPars+; 
#endif

class ConstTriggerPars : public TriggerPars {
	public:
		ConstTriggerPars(double m_Norm){ 
			Norm= m_Norm;
			model= eConstTrigger;		
			nPars= 1;
		};
		virtual ~ConstTriggerPars() {};
	public:
		static int GetParNumber() {return 1;}
		virtual double GetNorm() const {return Norm;}
		virtual double GetBreak() const {return 0;}	
		virtual double GetSmooth() const {return 0;}	
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(Norm);
			return pars;
		}
	ClassDef(ConstTriggerPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class ConstTriggerPars+; 
#endif

class SigmoidTriggerPars : public TriggerPars {
	public:
		SigmoidTriggerPars(double m_Norm,double m_Break,double m_Smooth){ 
			Norm= m_Norm; 
			Break= m_Break; 	
			Smooth= m_Smooth;
			model= eSigmoidTrigger;		
			nPars= 3;
		};
		virtual ~SigmoidTriggerPars() {};
	public:
		static int GetParNumber() {return 3;}
		virtual double GetNorm() const {return Norm;}
		virtual double GetBreak() const {return Break;}	
		virtual double GetSmooth() const {return Smooth;}	
		virtual std::vector<double> GetPars() const {
			std::vector<double> pars;
			pars.push_back(Norm);
			pars.push_back(Break);
			pars.push_back(Smooth);
			return pars;
		}
	ClassDef(SigmoidTriggerPars,1)
};
#ifdef __MAKECINT__
#pragma link C++ class SigmoidTriggerPars+; 
#endif


//======================================================
//==        SPECTRUM FITTING PARS
//======================================================
class FitPar : public TObject {
	public:
		FitPar(){};
		FitPar(std::string name,double val)
			: m_Name(name), m_Value(val)
		{
			m_StepSize= 0.1;//10% of value
			m_IsFixed= false;
			m_IsLimited= false;
			m_MinValue= -1;
			m_MaxValue= -1;
			m_ValueError= 0;
		}
		FitPar(std::string name,double val,double minval,double maxval)
			: m_Name(name), m_Value(val), m_MinValue(minval), m_MaxValue(maxval)
		{
			m_StepSize= 0.1;//10% of value
			m_IsFixed= false;
			m_IsLimited= true;
			m_ValueError= 0;
		}
		virtual ~FitPar() {};
	public:
	
		//Getters
		double GetValue() const {return m_Value;}
		double GetValueError() const {return m_ValueError;}
		std::string GetParName() const {return m_Name;}
		bool IsFixed() const {return m_IsFixed;}
		bool IsLimited() const {return m_IsLimited;}
		double GetMinValue() const {return m_MinValue;}
		double GetMaxValue() const {return m_MaxValue;}
		double GetStepSize() const {return m_StepSize;}

		//Setters
		void SetParName(std::string name){m_Name=name;}
		void SetValue(double val){m_Value=val;}
		void SetValueError(double val){m_ValueError=val;}
		void Fix(){m_IsFixed=true;}
		void Free(){m_IsFixed=false;}
		void SetLimits(double minval,double maxval){
			m_MinValue= minval;
			m_MaxValue= maxval;
			m_IsLimited= true;
		}
		void SetStep(double val){m_StepSize=val;}

	protected:
		std::string m_Name;
		double m_Value;
		double m_ValueError;
		bool m_IsFixed;
		bool m_IsLimited;
		double m_MinValue;
		double m_MaxValue;
		double m_StepSize;

	ClassDef(FitPar,1)
};
#ifdef __MAKECINT__
#pragma link C++ class FitPar+; 
#endif


class FitPars : public TObject {

	public:
		FitPars(){
			
		};		
		virtual ~FitPars() {
			Clear();	
		};

	public: 
		struct MatchName {
			MatchName(const std::string& name) : m_name(name) {}
 			bool operator()(const FitPar& obj) const {
   			return obj.GetParName() == m_name;
 			}
 			private:
   			const std::string& m_name;
		};
	
		int GetNPars() const {return (int)(m_FitPars.size());}
		
		void AddPar(FitPar& par){
			m_FitPars.push_back(par);
		}

		FitPar* GetPar(int index) {
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return 0;
			return &m_FitPars[index];
		}

		FitPar* FindPar(int& index,std::string name){
			if(m_FitPars.empty()) return 0;
			std::vector<FitPar>::iterator it = find_if(m_FitPars.begin(), m_FitPars.end(), MatchName(name));
			if(it==m_FitPars.end()) return 0;
			size_t pos= it-m_FitPars.begin();
			index= pos;
			return GetPar(index);
		}

		int SetParValue(int index,double value){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetValue(value);
			return 0;
		}
		int SetParValue(std::string name,double value){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].SetValue(value);
			return 0;
		}
		int FixPar(int index,double value){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetValue(value);
			m_FitPars[index].Fix();
			return 0;
		}
		int FixPar(std::string name,double value){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].SetValue(value);
			m_FitPars[index].Fix();
			return 0;
		}
		int FreePar(int index){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].Free();
			return 0;
		}
		int FreePar(std::string name){
			int index= -1;
			FitPar* par= FindPar(index,name);
			if(!par) return -1;
			m_FitPars[index].Free();
			return 0;
		}
		int LimitPar(int index,double xmin,double xmax){
			if(m_FitPars.empty() || index<0 || index>=GetNPars() ) return -1;
			m_FitPars[index].SetLimits(xmin,xmax);
			return 0;
		}
		int LimitPar(std::string name,double xmin,double xmax){
			int index= -1;
			FitPar* par= FindPar(index,name);	
			if(!par) return -1;
			m_FitPars[index].SetLimits(xmin,xmax);
			return 0;
		}
	private:
		void Clear(){
			m_FitPars.clear();	
		}
	protected:
		std::vector<FitPar> m_FitPars;
	
	ClassDef(FitPars,1)

};//close FitPars()

#ifdef __MAKECINT__
#pragma link C++ class FitPars+; 
#endif



class SpectrumPars : public TObject {
	public:
		virtual ~SpectrumPars() {};
	public:
		//Pure virtual
		virtual double GetIntegral(double xmin,double xmax) = 0;
		virtual SpectrumPars* Clone() const = 0;

		//Standard fcn
		virtual int GetModel() const {return model;}
		virtual int GetNPars() const {return pars.GetNPars();}
		static int GetParNumber() {return 0;}
		virtual FitPar* GetPar(int index){return pars.GetPar(index);}
		virtual FitPars GetPars() const {return pars;}
		virtual int SetParValue(int index,double value){return pars.SetParValue(index,value);}
		virtual int FixPar(int index,double val) {
			if(SetParValue(index,val)<0) return -1;
			return pars.FixPar(index,val);
		}
		virtual int FreePar(int index) {return pars.FreePar(index);}
		virtual int LimitPar(int index,double xmin,double xmax) {return pars.LimitPar(index,xmin,xmax);}


		virtual int SetParValue(std::string name,double value){return pars.SetParValue(name,value);}
		virtual int FixPar(std::string name,double val) {
			if(SetParValue(name,val)<0) return -1;
			return pars.FixPar(name,val);
		}
		virtual int FreePar(std::string name) {return pars.FreePar(name);}
		virtual int LimitPar(std::string name,double xmin,double xmax) {return pars.LimitPar(name,xmin,xmax);}

	protected:
		FitPars pars;			
		int model;
		static int nPars;

	ClassDef(SpectrumPars,1)		
		
};
#ifdef __MAKECINT__
#pragma link C++ class SpectrumPars+; 
#endif


enum SpectrumModel {
	ePowerLaw= 0,
	eBrokenPowerLaws= 1,
	eFlat= 2,
	eSmoothBrokenPowerLaws= 3
};

class FlatSpectrumPars : public SpectrumPars {

	public:		
		FlatSpectrumPars(FitPar normPar){ 
			model= eFlat;	
			nPars= 1;
			pars.AddPar(normPar);
		};
		virtual ~FlatSpectrumPars() {};
	public:
		//Overridden method
		static int GetParNumber() {return 1;}
		virtual double GetIntegral(double xmin,double xmax) {
			double Gamma= (pars.GetPar(1))->GetValue();
			return MathUtils::GetPowerLawIntegral(Gamma,xmin,xmax);
		}
		SpectrumPars* Clone() const { 
			return new FlatSpectrumPars(*this); 
		}

		
	protected:
		
	
	ClassDef(FlatSpectrumPars,1)	
};


class PowerLawPars : public SpectrumPars {
	
	public:
		PowerLawPars(FitPar normPar,FitPar gammaPar) 
		{
			model= ePowerLaw;	
			nPars= 2;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
		};
		
		virtual ~PowerLawPars() {};
	public:
		//Overridden method
		static int GetParNumber() {return 2;}
		virtual double GetIntegral(double xmin,double xmax) {
			double Gamma= (pars.GetPar(1))->GetValue();
			return MathUtils::GetPowerLawIntegral(Gamma,xmin,xmax);
		}
		SpectrumPars* Clone() const { 
			return new PowerLawPars(*this); 
		}

		
	protected:
		
	ClassDef(PowerLawPars,1)		
		
};
#ifdef __MAKECINT__
#pragma link C++ class PowerLawPars+; 
#endif


class BrokenPowerLawsPars : public SpectrumPars {
	
	public:
		BrokenPowerLawsPars(FitPar normPar,FitPar gammaPar,FitPar gamma2Par,FitPar gamma3Par,FitPar breakPar,FitPar cutoffPar) 
		{
			model= eBrokenPowerLaws;	
			nPars= 6;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(gamma2Par);
			pars.AddPar(gamma3Par);
			pars.AddPar(breakPar);
			pars.AddPar(cutoffPar);
		};
		
		virtual ~BrokenPowerLawsPars() {};
	public:
		//Overridden method
		static int GetParNumber() {return 6;}
		virtual double GetIntegral(double xmin,double xmax) {
			double Gamma= pars.GetPar(1)->GetValue();
			double Gamma2= pars.GetPar(2)->GetValue();
			double Gamma3= pars.GetPar(3)->GetValue();
			double Break= pars.GetPar(4)->GetValue();
			double Cutoff= pars.GetPar(5)->GetValue();
			return MathUtils::GetBrokenPowerLawIntegral(Gamma,Gamma2,Gamma3,Break,Cutoff,xmin,xmax);
		}
		SpectrumPars* Clone() const { 
			return new BrokenPowerLawsPars(*this); 
		}

		
	protected:
		
	ClassDef(BrokenPowerLawsPars,1)		
		
};
#ifdef __MAKECINT__
#pragma link C++ class BrokenPowerLawsPars+; 
#endif

class SmoothCutoffPowerLaws : public SpectrumPars {
	
	public:
		SmoothCutoffPowerLaws(FitPar normPar,FitPar gammaPar,FitPar gamma2Par,FitPar breakPar,FitPar cutoffPar,FitPar smoothCutoffPar) 
		{
			model= eSmoothBrokenPowerLaws;	
			nPars= 6;
			pars.AddPar(normPar);
			pars.AddPar(gammaPar);
			pars.AddPar(gamma2Par);
			pars.AddPar(breakPar);
			pars.AddPar(cutoffPar);
			pars.AddPar(smoothCutoffPar);
		};
		
		virtual ~SmoothCutoffPowerLaws() {};
	public:
		//Overridden method
		static int GetParNumber() {return 6;}
		virtual double GetIntegral(double xmin,double xmax) {
			return 1;
		}
		SpectrumPars* Clone() const { 
			return new SmoothCutoffPowerLaws(*this); 
		}
		
		
	protected:
		
	ClassDef(SmoothCutoffPowerLaws,1)		
		
};
#ifdef __MAKECINT__
#pragma link C++ class SmoothCutoffPowerLaws+; 
#endif


}//close namespace


#endif
