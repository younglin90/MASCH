#pragma once
#include <iostream>
#include <numeric>
#include <vector>
// #include <math.h>
#include <cmath>
#include <mpi.h>
using namespace std;

// ==============================
class MASCH_Math_Sparse_Matrix {
private:
public:
	vector<int> CSR_i;
	vector<int> CSR_displ;
	vector<double> val;
};


class MASCH_Math_Matrix {
private:
public:
	vector<double> val;
};


class MASCH_Math {
private:
public:

	// void GaussSeidel(double* a, double* b, double* x, int n);
	// void GaussSeidel(vector<vector<double>>& AA, int N);
	void GaussSeidel(double* AA, int N);
	void GaussSeidelSOR(double* A, int N);
	// void GaussSeidelSOR(vector<double>& A, int N);
	void GaussSeidelSOR(vector<vector<double>>& A);
	
	void calcUnitNormals_Area3dPolygon(
		int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
		vector<double>& unitNormals, double& area,
		double& x, double& y, double& z,
		double& VSn, vector<double>& cellCentroid );
};

