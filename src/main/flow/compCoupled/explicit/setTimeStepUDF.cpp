
#include "../../../../others/solvers.h"


void MASCH_Solver::setTimeStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// for(int i=0; i<3; ++i)
	{
		int id_c = controls.getId_cellVar("speed-of-sound");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_maxA = controls.getId_cellVar("maximum area");
		int id_minA = controls.getId_cellVar("minimum area");
		int id_vol = controls.getId_cellVar("volume");
		int id_mu = controls.getId_cellVar("viscosity");
		int id_rho = controls.getId_cellVar("density");
		int id_dt = controls.getId_fieldVar("time-step");
		double maxCFL = stod(controls.controlDictMap["maxCFL"]);
		double maxdt = stod(controls.controlDictMap["maxTimeStep"]);
		double dt = controls.timeStep;
		
		solver.calcTempStepCell.push_back(
		[id_dt,dt,id_vol,id_u,id_v,id_w,maxCFL,maxdt,id_c,id_mu,id_maxA,id_rho,
		id_minA](
		double* cells, double* fields) ->int {
			// double vol = cells[id_vol];
			// double u = cells[id_u]; double v = cells[id_v]; double w = cells[id_w];
			// // double char_speed = (sqrt(u*u+v*v+w*w)+cells[id_c]);
			// double char_speed = sqrt(u*u+v*v+w*w);
			// double tmp_dt = maxCFL*pow(vol,0.3)/(char_speed+1.e-100);
			// tmp_dt = min(tmp_dt,maxdt);
			// tmp_dt = min(tmp_dt,dt);
			
			// fields[id_dt] = min(fields[id_dt],tmp_dt);
			
			double maxA = cells[id_maxA];
			double minA = cells[id_minA];
			double vol = cells[id_vol];
			double c = cells[id_c];
			double u = cells[id_u]; double v = cells[id_v]; double w = cells[id_w];
			double char_speed = sqrt(u*u+v*v+w*w) + c;
			double convTime = maxCFL*vol/maxA/char_speed;
			double mu = cells[id_mu];
			double viscTime = maxCFL* minA / ( mu + 1.e-200 ) / 6.0 * id_rho; // in 3D
			fields[id_dt] = min(fields[id_dt],min(convTime,viscTime));
			
			return 0;
		});
	}
	
}


