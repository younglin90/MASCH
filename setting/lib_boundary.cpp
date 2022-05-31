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
int id_u, int id_v, int id_w, 
int id_p, int id_T,
int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces);

int U_boundary(
double time, double x, double y, double z, 
int id_u, int id_v, int id_w, 
int id_p, int id_T,
int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces){




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

                faces[id_uF] = Unormal * Uuniform;
		faces[id_vF] = 0.0;
		faces[id_wF] = 0.0;



		//faces[id_uF] = 0.5*(faces[id_uF]+cells[id_u]);
		//faces[id_vF] = 0.5*(faces[id_vF]+cells[id_v]);
		//faces[id_wF] = 0.5*(faces[id_wF]+cells[id_w]);

                //if( yn <= 0.3359 ) cout << fU << endl;
                //        }
                //

	return 0;
}





extern "C" int U_liq_boundary(
double time, double x, double y, double z,
int id_u, int id_v, int id_w, 
int id_p, int id_T,
int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces);

int U_liq_boundary(
double time, double x, double y, double z, 
int id_u, int id_v, int id_w, 
int id_p, int id_T,
int id_uF, int id_vF, int id_wF, 
int id_nx, int id_ny, int id_nz, double* cells, double* faces){

	if((x-0.01)*(x-0.01)+z*z<0.0005*0.0005){
		faces[id_uF] = 0.0;
		faces[id_vF] = 22.6;
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

	if((x-0.01)*(x-0.01)+z*z<0.0005*0.0005){
		faces[id_YF] = 1.0;
	}
	else{
		faces[id_YF] = 0.0;
	}

	return 0;
}
