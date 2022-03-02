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
	vector<vector<double>> parcel;
	
	
	
	vector<int> j_displ_CSR;
	vector<int> i_str_CSR;
	vector<double> Avar;
	vector<double> Xvar;
	vector<double> Bvar;
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
	
	void accumSparD(int i, int iEq, int jEq, double inp);
	double getSparD(int i, int iEq, int jEq);
	void accumSparLR(int i, int iL, int iEq, int jEq, double inp);
	void accumSparRL(int i, int iR, int iEq, int jEq, double inp);
	
	void accumB(int i, int iEq, double inp);
	double getB(int i, int iEq);
	
	void setX(int i, int iEq, double inp);
	double getX(int i, int iEq);
	
};

