#include <iostream>
#include <vector>
using namespace std;

// extern "C" void setFunctionPressure(double time, double x, double y, double z, double& fP);
// extern "C" void setFunctionUvelocity(double time, double x, double y, double z, double& fU);
// extern "C" void setFunctionVvelocity(double time, double x, double y, double z, double& fV);
// extern "C" void setFunctionWvelocity(double time, double x, double y, double z, double& fW);
// extern "C" void setFunctionTemperature(double time, double x, double y, double z, double& fT);
// extern "C" void setFunctionVolumeFractions(double time, double x, double y, double z, double& fVF);
// extern "C" void setFunctionMassFractions(double time, double x, double y, double z, double& fMF);


//extern "C" int p_boundary(double time, double x, double y, double z, int* inp_id, double* cells);
extern "C" int U_boundary(
double time, double x, double y, double z,
int id_u, int id_v, int id_w, int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces);

int U_boundary(
double time, double x, double y, double z, 
int id_u, int id_v, int id_w, int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces){

	double r = 0.00017;
	if( (y-0.6)*(y-0.6)+(z-0.6)*(z-0.6) <= r*r ){
        faces[id_uF] = 4.6;
	faces[id_vF] = 0.0;
	faces[id_wF] = 0.0;
	}
	else{
	faces[id_uF] = 0.0;
	faces[id_vF] = 0.0;
	faces[id_wF] = 0.0;
	}

	return 0;
}







extern "C" int Y_boundary(
double time, double x, double y, double z,
int id_Y, int id_YF,double* cells, double* faces);

int Y_boundary(
double time, double x, double y, double z, 
int id_Y, int id_YF, double* cells, double* faces){

	double r = 0.00017;
	if( (y-0.6)*(y-0.6)+(z-0.6)*(z-0.6) <= r*r ){
        faces[id_YF] = 1.0;
	}
	else{
	faces[id_YF] = 0.0;
	}

	return 0;
}
