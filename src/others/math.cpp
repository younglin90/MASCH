
#include "./math.h"

/*
new, 2021-10-15
http://paulbourke.net/geometry/polygonmesh/
*/
void MASCH_Math::calcUnitNormals_Area3dPolygon(
	int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	vector<double>& unitNormals, double& area,
	double& x, double& y, double& z,
	double& VSn, vector<double>& cellCentroid ){
		
	double x_tmp = 0.0;  double y_tmp = 0.0; double z_tmp = 0.0;
	double Vx0_old = 0.0;  double Vy0_old = 0.0; double Vz0_old = 0.0;
	vector<double> unitNormals_tmp(3,0.0);
	vector<double> cellCentroid_tmp(3,0.0);
	double area_tmp = 0.0;
	double VSn_tmp = 0.0;
	
	double x_avg, y_avg, z_avg;
	x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)Vx.size();
	y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)Vy.size();
	z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)Vz.size();
	
	for(int two=0; two<5; ++two){
		
		unitNormals_tmp[0] = 0.0;
		unitNormals_tmp[1] = 0.0;
		unitNormals_tmp[2] = 0.0;
		
		cellCentroid_tmp[0] = 0.0;
		cellCentroid_tmp[1] = 0.0;
		cellCentroid_tmp[2] = 0.0;
		
		area_tmp = 0.0;
		VSn_tmp = 0.0;
		double VSn_x = 0.0; double VSn_y = 0.0; double VSn_z = 0.0;
		double Vx0 = 0.0; double Vy0 = 0.0; double Vz0 = 0.0;
		
		if(two==0){
			// Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			
			Vx0 = x_avg; Vy0 = y_avg; Vz0 = z_avg;
		}
		else{
			Vx0 = Vx0_old + 1.0*(x_tmp-Vx0_old); 
			Vy0 = Vy0_old + 1.0*(y_tmp-Vy0_old); 
			Vz0 = Vz0_old + 1.0*(z_tmp-Vz0_old); 
		}
		Vx0_old = Vx0; Vy0_old = Vy0; Vz0_old = Vz0;
			// Vx0 = 0; Vy0 = 0; Vz0 = 0;
			// Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			// Vx0 = Vx.back(); Vy0 = Vy.back(); Vz0 = Vz.back();
			
		x_tmp = 0.0;  y_tmp = 0.0; z_tmp = 0.0;
		
		double totArea = 0.0;
		
		for(int i = 0; i < n; ++i ) {
			int b=i;
			int c=i+1;
			if(i==n-1) c=0;
			
			double vect_A[3];
			double vect_B[3];
			
			vect_A[0] = Vx[b]-Vx0; vect_A[1] = Vy[b]-Vy0; vect_A[2] = Vz[b]-Vz0;
			vect_B[0] = Vx[c]-Vx0; vect_B[1] = Vy[c]-Vy0; vect_B[2] = Vz[c]-Vz0;
			// vect_A[0] = Vx[b]; vect_A[1] = Vy[b]; vect_A[2] = Vz[b];
			// vect_B[0] = Vx[c]; vect_B[1] = Vy[c]; vect_B[2] = Vz[c];
		  
			double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
			double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
			double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

			unitNormals_tmp[0] += Nx_tmp;
			unitNormals_tmp[1] += Ny_tmp;
			unitNormals_tmp[2] += Nz_tmp;
			
			double tmp_area = sqrtl(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp);
			// double tmp_area = pow(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp,0.5);
			double rx = (Vx0+Vx[b]+Vx[c]) / 3.0;
			double ry = (Vy0+Vy[b]+Vy[c]) / 3.0;
			double rz = (Vz0+Vz[b]+Vz[c]) / 3.0;
			x_tmp += rx*tmp_area; y_tmp += ry*tmp_area; z_tmp += rz*tmp_area;
			
			totArea += tmp_area;
			
			VSn_tmp += Vx0*Nx_tmp; VSn_tmp += Vy0*Ny_tmp; VSn_tmp += Vz0*Nz_tmp;
			
			cellCentroid_tmp[0] += 2.0*Nx_tmp*(
				(Vx0+Vx[b])*(Vx0+Vx[b]) + (Vx[b]+Vx[c])*(Vx[b]+Vx[c]) + (Vx[c]+Vx0)*(Vx[c]+Vx0));
			cellCentroid_tmp[1] += 2.0*Ny_tmp*(
				(Vy0+Vy[b])*(Vy0+Vy[b]) + (Vy[b]+Vy[c])*(Vy[b]+Vy[c]) + (Vy[c]+Vy0)*(Vy[c]+Vy0));
			cellCentroid_tmp[2] += 2.0*Nz_tmp*(
				(Vz0+Vz[b])*(Vz0+Vz[b]) + (Vz[b]+Vz[c])*(Vz[b]+Vz[c]) + (Vz[c]+Vz0)*(Vz[c]+Vz0));
			
			
			
		}
		
		
		double mag_unitNormals = sqrtl(
			unitNormals_tmp[0]*unitNormals_tmp[0]+
			unitNormals_tmp[1]*unitNormals_tmp[1]+
			unitNormals_tmp[2]*unitNormals_tmp[2]);
		unitNormals_tmp[0] /= mag_unitNormals;
		unitNormals_tmp[1] /= mag_unitNormals;
		unitNormals_tmp[2] /= mag_unitNormals;
		
		area_tmp = mag_unitNormals;
		
		x_tmp /= totArea;
		y_tmp /= totArea;
		z_tmp /= totArea;
	
	}
	
	unitNormals.clear();
	unitNormals.resize(3,0.0);
	cellCentroid.clear();
	cellCentroid.resize(3,0.0);
	
	for(int i=0; i<3; ++i){
		unitNormals[i] = unitNormals_tmp[i];
		cellCentroid[i] = cellCentroid_tmp[i];
	}
	
	area = area_tmp;
	
	VSn = VSn_tmp;
		
	x = x_tmp; y = y_tmp; z = z_tmp;
	
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
}





void MASCH_Math::GaussSeidelSOR(double* A, int N) { 


	// double N = A.size();
	
	// double** L = new double*[N];
	// double** U = new double*[N];
	// for(int i=0; i<N; ++i){
		// L[i] = new double[N];
		// U[i] = new double[N];
	// }
	// double* B = new double[N];
	// double* D = new double[N];
	// double* X = new double[N];
	
	double L[N][N];
	double U[N][N];
	double B[N];
	double D[N];
	double X[N];
	
	
	for(int k=1; k<=N-1; ++k){
		for(int i=k+1; i<=N; ++i){
			// double coeff = A[i-1][k-1]/A[k-1][k-1];
			double coeff = A[(i-1)*N+k-1]/A[(k-1)*N+k-1];
			L[i-1][k-1] = coeff;
			for(int j=k+1; j<=N; ++j){
				// A[i-1][j-1] = A[i-1][j-1] - coeff*A[k-1][j-1];
				A[(i-1)*N+j-1] = A[(i-1)*N+j-1] - coeff*A[(k-1)*N+j-1];
			}
		}
	}
	
	for(int i=1; i<=N; ++i){
		L[i-1][i-1] = 1.0;
	}
	
	for(int j=1; j<=N; ++j){
		for(int i=1; i<=j; ++i){
			// U[i-1][j-1] = A[i-1][j-1];
			U[i-1][j-1] = A[(i-1)*N+j-1];
		}
	}
	
	for(int k=1; k<=N; ++k){
		B[k-1] = 1.0;
		D[0] = B[0];
		for(int i=2; i<=N; ++i){
			D[i-1] = B[i-1];
			for(int j=1; j<=i-1; ++j){
				D[i-1] = D[i-1] - L[i-1][j-1]*D[j-1];
			}
		}
		
		X[N-1] = D[N-1] / U[N-1][N-1];
		
		// cout << U[N-1][N-1] << endl;
		
		for(int i=N-1; i>=1; --i){
			X[i-1] = D[i-1];
			for(int j=N; j>=i+1; --j){
				X[i-1] = X[i-1] - U[i-1][j-1]*X[j-1];
			}
			X[i-1] = X[i-1]/U[i-1][i-1];
		}
		
		for(int i=1; i<=N; ++i){
			// A[i-1][k-1] = X[i-1];
			A[(i-1)*N+k-1] = X[i-1];
		}
		
		B[k-1] = 0.0;
	}


	// for(int i=0; i<N; ++i){
		// delete L[i];
		// delete U[i];
	// }
	// delete L;
	// delete U;
	// delete B;
	// delete D;
	// delete X;

}
