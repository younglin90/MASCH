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
	double getHO_MSTACS(double phiUU, double phiU, double phiD, double coDD, double gamF);
};

class MASCH_Gradient {
private:
public:
	double n_weight = 1.0;

	void init(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void gaussGreen(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int inp, string inp_bc, int oup0, int oup1, int oup2);
	void leastSquare(
		MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, string inp_s);

	void leastSquare(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
		vector<string>& inp_cell, vector<string>& inp_bcFace);
	void leastSquare_zeroGradient(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
		vector<string>& inp_cell);
};

// class MASCH_DPM {
// private:
// public:
	// // void calcMechanicalForceModel(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	// // void calcBreakupModel(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	// // void calcHeatAndMassTransferModel(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
// };

class MASCH_Solver {
private:
public:
	MASCH_NVD NVD;
	MASCH_Math math;
	MASCH_Gradient calcGradient;
	// MASCH_DPM dpm;
	
	void dpm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	int calcTimeStepParcel(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void addParcelModel(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void parcelLoop(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void searchLocationParcelsToOutside(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateProcRightParcels(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void refreshParcels(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	// int get_icell(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
		// double parcel_x, double parcel_y, double parcel_z);
	
	
	void calcCurvature(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
		vector<string>& inp_cell);
		
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
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
	vector<vector<Bound_Funct_type>> calcBoundFacePrimVal;
	
	// additional values 관련
	using Bound_Face_Funct_type = function<int(double* faces)>;
	using Bound_Cell_Funct_type = function<int(double* cells)>;
	vector<Bound_Face_Funct_type> calcFaceAddiVal;
	vector<Bound_Cell_Funct_type> calcCellAddiVal;
	
	// 올드값 관련
	vector<string> saveCurrValueNames, saveOld1ValueNames, saveOld2ValueNames;
	
	// gradient 관련
	vector<vector<string>> gradLSIds_cell_name, gradLSIds_bcFace_name;
	
	// 곡률 관련
	vector<vector<string>> curvatureIds_cell_name, curvatureIds_bcFace_name;
	// vector<unsigned short> gradGGIds_inp, gradGGIds_x_out, gradGGIds_y_out, gradGGIds_z_out;
	// vector<unsigned short> gradLSIds_inp, gradLSIds_x_out, gradLSIds_y_out, gradLSIds_z_out;
	
	// high-order recon. 관련
	using HO_Funct_type = function<int(double* fields, double* cellsL, double* cellsR, double* faces)>;
	vector<HO_Funct_type> calcHO_FaceVal;
	
	// implicit 인지
	vector<bool> checkImplicit;
	// convective flux 관련
	using conv_diff_Funct_type = 
		function<int(
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR)>;
	vector<conv_diff_Funct_type> calcConvFlux;
	using conv_diff_Funct_BC_type = 
		function<int(
		double* fields, double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB)>;
	vector<vector<vector<conv_diff_Funct_BC_type>>> calcConvFlux_BC;
	// diffusive flux 관련
	vector<conv_diff_Funct_type> calcLaplFlux, calcNLaplFlux;
	
	// source 관련
	using source_Funct_type = 
		function<int(double* cells, double* fluxA, double* fluxB)>;
	vector<source_Funct_type> calcSource;
	vector<bool> checkImplicitSource;
	
	// temporal 관련
	using temporal_Funct_type = 
		function<int(double* cells, double* fields, double* fluxA, double* fluxB)>;
	vector<temporal_Funct_type> calcTemporal;
	
	// time-step 관련
	using dtstep_Cell_Funct_type = 
		function<int(double* cells, double* fields)>;
	vector<dtstep_Cell_Funct_type> calcTempStepCell;
	using dtstep_Face_Funct_type = 
		function<int(double* cellsL, double* cellsR, double* faces, double* fields)>;
	vector<dtstep_Face_Funct_type> calcTempStepFace;
	
	using update_prim_Funct_type = 
		function<int(double* fields, double* cells, double* Xvalues)>;
	vector<update_prim_Funct_type> calcUpdatePrim;
	
	
	// DPM 관련
	vector<function<int(double* fields, vector<vector<double>>&)>> calcDPM_injection;
	vector<function<int(double* cells, double* fields, double* parcels, double* fluxB)>> calcDPM_parcelLoop;
	
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
	void setCurvatureFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	// void setFValFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setTermsCellLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setTermsFaceLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	// void setDiffFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	// void setSourFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setTimeStepFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setUpdatePrimFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	void setDPMFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls);
	
	void fvm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void fvm_inline(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	// 초기 셋팅
	void updateCellAddiValues_All(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateBoundaryFacePrimValues_All(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void gradientTerms_All(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void curvatureTerms_All(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void updateProcRightCellValues_All(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	void initOldValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	// 아웃터 펑션
	void updateOldValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var);
	
	// 계산
	void updateCellAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateProcRightCellPrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateProcRightCellAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateBoundaryFacePrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateBoundaryFaceAddiValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateBoundaryFaceValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void gradientTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void curvatureTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateProcRightCellGradValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void highOrderTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void cellLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void faceLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void linearSystem(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void updateCellPrimValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	void calcTempSteps(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq);
	
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
	
	// 선형솔버 관련
	void solveAMGCL(
		vector<int>& istr_CSR, vector<int>& j_CSR, vector<double>& Aval, 
		vector<double>& Bval, vector<double>& Xval);
	
	
};

