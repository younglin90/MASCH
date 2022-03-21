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


extern "C" int setBoundaryFunctionPressure(double time, double x, double y, double z, int* inp_id, double* cells);
int setBoundaryFunctionPressure(double time, double x, double y, double z, int* inp_id, double* cells){
cout << "AA" << endl;
	cells[inp_id[0]] = 55.0;
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

