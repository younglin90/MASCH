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
	cells[inp_id[0]] = 29000.0;
	if(y<0){
		cells[inp_id[0]] = 29000.0;
	}
	
	return 0;
}

extern "C" int U_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int U_initial(double time, double x, double y, double z, int* inp_id, double* cells){



        double Uuniform = 674.0;

                vector<double> ynBd(13,0.0);
                vector<double> Ubd(13,0.0);

                ynBd[0] = -100.0;
                ynBd[1] = 1.e-12;
                ynBd[2] = 0.3359173126615005;
                ynBd[3] = 0.6976744186046524;
                ynBd[4] = 1.1369509043927657;
                ynBd[5] = 1.757105943152455;
                ynBd[6] = 2.4806201550387605;
                ynBd[7] = 3.2816537467700266;
                ynBd[8] = 4.211886304909561;
                ynBd[9] = 5.116279069767443;
                ynBd[10] = 6.124031007751938;
                ynBd[11] = 6.511627906976744;
                ynBd[12] = 7.881136950904393;

                Ubd[0] = 0.0;
                Ubd[1] = 0.0;
                Ubd[2] = 0.2320108557119938;
                Ubd[3] = 0.40308578111645904;
                Ubd[4] = 0.5654341535093637;
                Ubd[5] = 0.7161246307708455;
                Ubd[6] = 0.8174621024695246;
                Ubd[7] = 0.8839608354699895;
                Ubd[8] = 0.9330139594860032;
                Ubd[9] = 0.9704691646799513;
                Ubd[10] = 0.9904862579281182;
                Ubd[11] = 0.9932751563132562;
                Ubd[12] = 0.9957791672289447;

                double Unormal = 0.0;
                double yn = y/0.001;

                for(int i=0; i<12; ++i){
                        if(ynBd[i]<yn && yn<=ynBd[i+1]){
                           Unormal = Ubd[i] + (Ubd[i+1]-Ubd[i])/(ynBd[i+1]-ynBd[i])*(yn-ynBd[i]);
                        }
                }
                if(ynBd[12]<yn) Unormal = 1.0;

                 cells[inp_id[0]] = Unormal * Uuniform;
                cells[inp_id[1]] = 0.0;
                cells[inp_id[2]] = 0.0;

                //if( yn <= 0.3359 ) cout << fU << endl;
                //                //        }
                //                                //
                //

		cells[inp_id[0]] = 105.2;

	if(y<0){

		cells[inp_id[0]] = 0.0;
	}
	

	return 0;
}

extern "C" int T_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int T_initial(double time, double x, double y, double z, int* inp_id, double* cells){
	cells[inp_id[0]] = 300.0;
	if(y<0){
		cells[inp_id[0]] = 300.0;
	}
	
	return 0;
}

extern "C" int Y_initial(double time, double x, double y, double z, int* inp_id, double* cells);
int Y_initial(double time, double x, double y, double z, int* inp_id, double* cells){
	cells[inp_id[0]] = 0.0;
	if(y<0){
		cells[inp_id[0]] = 0.0;
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

