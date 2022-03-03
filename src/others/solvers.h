#pragma once
#include <vector>
// #include <valarray>
using namespace std;

#include "./mesh.h"  
#include "./controls.h" 
#include "./variables.h" 
#include "./math.h" 

class MASCH_Control;
class MASCH_Mesh;
class MASCH_Variables;
class MASCH_Math;

class MASCH_NVD {
private:
public:
	double Minmod(double phiUU, double phiU, double phiD);
};

class MASCH_Gradient {
private:
public:
	double n_weight = 1.0;

	void init(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void gaussGreen(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int inp, string inp_bc, int oup0, int oup1, int oup2);
	void leastSquare(
		MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, string inp_s);

	void leastSquare(
		MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, vector<string>& inp_s);
};

class MASCH_Solver {
private:
public:
	MASCH_NVD NVD;
	MASCH_Math math;
	MASCH_Gradient calcGradient;
	// double Minmod(double phiUU, double phiU, double phiD);

	// using functType = function<double(double* inp)>;
	// using functType2 = function<double(double* inpL, double* inpR, double* inpF)>;
	// using functTypeRecon = function<double(double* inpL, double* inpR)>;
	// valarray<d
	// vector<int> gradientInpStr;
	// vector<int> gradientCellInp;
	// vector<string> gradientCellInpBC;
	// vector<vector<int>> gradientCellOut;
	// vector<functType> gradient;
	
	// boundary face 관련
	using Bound_Funct_type = function<int(double* cells, double* faces)>;
	vector<Bound_Funct_type> calcBoundFacePrimVal;
	
	// additional values 관련
	using Bound_Face_Funct_type = function<int(double* faces)>;
	using Bound_Cell_Funct_type = function<int(double* cells)>;
	vector<Bound_Face_Funct_type> calcFaceAddiVal;
	vector<Bound_Cell_Funct_type> calcCellAddiVal;
	
	// 올드값 관련
	vector<unsigned short> saveCurrIds, saveOld1Ids, saveOld2Ids;
	
	// gradient 관련
	vector<string> gradLSIds_name;
	// vector<unsigned short> gradGGIds_inp, gradGGIds_x_out, gradGGIds_y_out, gradGGIds_z_out;
	// vector<unsigned short> gradLSIds_inp, gradLSIds_x_out, gradLSIds_y_out, gradLSIds_z_out;
	
	// high-order recon. 관련
	using HO_Funct_type = function<int(double* cellsL, double* cellsR, double* faces)>;
	vector<HO_Funct_type> calcHO_FaceVal;
	
	// convective flux 관련
	using conv_Funct_type = 
	function<int(double* cellsL, double* cellsR, 
	double* faces, double* fluxA, double* fluxB)>;
	vector<conv_Funct_type> calcConvFlux;
	vector<bool> checkImplicitConvFlux;
	
	// diffusive flux 관련
	using diff_Funct_type = 
	function<int(double* cellsL, double* cellsR, 
	double* faces, double* fluxA, double* fluxB)>;
	vector<diff_Funct_type> calcLaplFlux, calcNLaplFlux;
	vector<bool> checkImplicitLaplFlux;
	
	// source 관련
	using source_Funct_type = 
	function<int(double* cells, double* fluxA, double* fluxB)>;
	vector<source_Funct_type> calcSource;
	vector<bool> checkSource;
	
	// temporal 관련
	using temporal_Funct_type = 
	function<int(double* cells, double* fields, double* fluxA, double* fluxB)>;
	vector<temporal_Funct_type> calcTemporal;
	
	// time-step 관련
	using dtstep_Cell_Funct_type = 
	function<int(double* faces, double* fields)>;
	vector<dtstep_Cell_Funct_type> calcTempStepCell;
	using dtstep_Face_Funct_type = 
	function<int(double* cellsL, double* cellsR, double* faces, double* fields)>;
	vector<dtstep_Face_Funct_type> calcTempStepFace;
	
	
	// vector<function<int(double* inpL, double* inpR, double* inpF)>> reconstruction;
	// vector<function<double(double* inpL, double* inpR, double* inpF)>> faceValue;
	// vector<function<int(double* inpL, double* inpR, double* inpF, double* out)>> convective;
	// vector<function<int(double* inpL, double* inpR, double* inpF, double* out)>> laplacian;
	// vector<function<int(double* inpL, double* inpR, double* inpF, double* out)>> nonLaplacian;
	// vector<function<double(double* inp)>> source;
	
	// vector<int> temporalInpStr;
	// vector<int> temporalInp;
	// vector<function<double(double* inp)>> temporal;

	void setBoundaryFunctions(
		MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void setFunctions(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setTempFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setRecoFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setFValFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setConvFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setDiffFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setSourFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setTempStepFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	
	void fvm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void fvm_inline(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void updateCellAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	// void updateFaceAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateProcRightCellPrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateProcRightCellAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateBoundaryFacePrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateBoundaryFaceAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateBoundaryFaceValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void initOldValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateOldValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void gradientTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void highOrderTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void cellLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void faceLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void linearSystem(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateCellPrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	
	void calcTempSteps(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	void eosIdeal(
		double cv, double gam,
		double& P, double& U, double& V, double& W, double& T,
		double& rho, double& C, double& Ht,
		double& drhodP, double& drhodT, double& dhdP, double& dhdT);
	void eosNASG(
		double pinf, double cv, double gam, double bNASG, double q,
		double& P, double& U, double& V, double& W, double& T,
		double& rho, double& C, double& Ht,
		double& drhodP, double& drhodT, double& dhdP, double& dhdT);
	
};

