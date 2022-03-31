#pragma once
#include <vector>
#include <dlfcn.h>
using namespace std;

#include "mesh.h"  
#include "controls.h" 
#include "mpi.h" 
#include "math.h" 


class MASCH_Control;
class MASCH_MPI;
class MASCH_Mesh;
class MASCH_Math;

// ==============================
// 변수 클래스
class MASCH_Variables {
private:
public:
	void addParcel(int nParcelVar, vector<double> vars){
		(*this).parcels.push_back(vector<double>());
		(*this).parcels.back().resize(nParcelVar);
		for(int j=0; j<nParcelVar; ++j){
			(*this).parcels.back()[j] = vars[j];
		}
	}
				
	vector<vector<double>> cells;
	vector<vector<double>> procRightCells;
	vector<vector<double>> faces;
	vector<vector<double>> points;
	vector<vector<double>> edges;
	vector<vector<double>> boundaries;
	vector<double> fields;
	// double** procRightCells;
	// double** faces;
	// double** points;
	// double** edges;
	// double** boundaries;
	// double* fields;
	vector<vector<double>> parcels;
	
	
	
	vector<vector<int>> j_displ_CSR;
	vector<vector<int>> i_str_CSR;
	vector<vector<double>> Avalues;
	vector<vector<double>> Xvalues;
	vector<vector<double>> Bvalues;
	// double* Avar;
	// double* Xvar;
	// double* Bvar;
	vector<int> cellcell_displ;
	vector<int> face_LR_displ;
	vector<int> face_RL_displ;
	vector<int> cellnFaces;
	// int* cellcell_displ;
	// int* face_LR_displ;
	// int* face_RL_displ;
	// int* cellnFaces;
	int nCells;
	
	void setSparCSR(MASCH_Mesh& mesh, MASCH_Control& controls);
	
	void accumSparD(int iSegEq, int i, int iEq, int jEq, double inp);
	double getSparD(int iSegEq, int i, int iEq, int jEq);
	void accumSparLR(int iSegEq, int i, int iL, int iEq, int jEq, double inp);
	void accumSparRL(int iSegEq, int i, int iR, int iEq, int jEq, double inp);
	void clearLinearSystems(int iSegEq);
	
	void accumB(int iSegEq, int i, int iEq, double inp);
	double getB(int iSegEq, int i, int iEq);
	
	void setX(int iSegEq, int i, int iEq, double inp);
	double getX(int iSegEq, int i, int iEq);
	
	
	// 임시변수
	// vector<double> tmp_sparA;
	// vector<double> tmp_B;
	// vector<double> tmp_X;
	
};

