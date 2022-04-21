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


extern "C" int p_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int p_initial(double time, double x, double y, double z, int* inp_id, double* cells){
	cells[inp_id[0]] = 1.e5;
	if(x<0.017){
		cells[inp_id[0]] = 129680.0;
	}
	
	return 0;
}

extern "C" int U_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int U_initial(double time, double x, double y, double z, int* inp_id, double* cells){

	cells[inp_id[0]] = 0.0;
	cells[inp_id[1]] = 0.0;
	cells[inp_id[2]] = 0.0;
	if(x<0.017){
        cells[inp_id[0]] = 65.7;
        cells[inp_id[1]] = 0.0;
        cells[inp_id[2]] = 0.0;
	}
	
	return 0;
}

extern "C" int T_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int T_initial(double time, double x, double y, double z, int* inp_id, double* cells){
	cells[inp_id[0]] = 300.0;
	if(x<0.017){
		cells[inp_id[0]] = 323.29;
	}
	
	return 0;
}

extern "C" int Y_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int Y_initial(double time, double x, double y, double z, int* inp_id, double* cells){

	double r = 2.5*0.001*0.5;

	double dx = x-19.0*0.001;
	double dy = y-15.0*0.001;
	double dz = z-15.0*0.001;


	cells[inp_id[0]] = 0.0;

	if(dx*dx+dy*dy+dz*dz<=r*r){

		cells[inp_id[0]] = 1.0;
	}


	
	return 0;
}

// void setFunctionUvelocity(double time, double x, double y, double z, double& fU){
	// // cout << "BBB" << endl;
	// fU=0.0;
// }
// void setFunctionVvelocity(double time, double x, double y, double z, double& fV){
	// // cout << "BBB" << endl;
	// double radius = 0.01/2.0;
	// if( (x-0.01)*(x-0.01)+z*z <= radius*radius ){
		// fV = 22.6;
	// }
	// else{
		// fV = 0.0;
	// }
// }
// void setFunctionWvelocity(double time, double x, double y, double z, double& fW){
	// // cout << "BBB" << endl;
	// fW=0.0;
// }
// void setFunctionTemperature(double time, double x, double y, double z, double& fT)
// {
	// // cout << "CCC" << endl;
// }
// void setFunctionMassFractions(double time, double x, double y, double z, double& fMF)
// {
	// double radius = 0.01/2.0;
	// if( (x-0.01)*(x-0.01)+z*z <= radius*radius ){
		// fMF = 1.0;
	// }
	// else{
		// fMF = 0.0;
	// }
	// // cout << "EEE" << endl;
// }

// void setFunctionVolumeFractions(double time, double x, double y, double z, double& fVF)
// {
	// // cout << "DDD" << endl;
// }

