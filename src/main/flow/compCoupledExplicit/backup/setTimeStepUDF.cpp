
#include "../../../others/solvers.h"


void MASCH_Solver::setTimeStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// for(int i=0; i<3; ++i)
	{
		int id_c = controls.getId_cellVar("speed-of-sound");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_vol = controls.getId_cellVar("volume");
		int id_dt = controls.getId_fieldVar("time-step");
		double maxCFL = stod(controls.controlDictMap["maxCFL"]);
		double maxdt = stod(controls.controlDictMap["maxTimeStep"]);
		double dt = controls.timeStep;
		
		solver.calcTempStepCell.push_back(
		[id_dt,dt,id_vol,id_u,id_v,id_w,maxCFL,maxdt,id_c](
		double* cells, double* fields) ->int {
			double vol = cells[id_vol];
			double c = cells[id_c];
			double u = cells[id_u]; double v = cells[id_v]; double w = cells[id_w];
			// double char_speed = (sqrt(u*u+v*v+w*w)+cells[id_c]);
			double char_speed = sqrt(u*u+v*v+w*w);
			double tmp_dt = maxCFL*pow(vol,0.3)/(char_speed+c);
			fields[id_dt] = min(fields[id_dt],tmp_dt);
			
			return 0;
		});
	}
	
}
