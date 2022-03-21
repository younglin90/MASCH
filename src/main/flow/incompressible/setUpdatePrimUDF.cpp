
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["segregated.relaxationFactors"]);
	{
		
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		
		// double relaxation_factor = stod(controls.fvSolutionMap["segregated"]);
		double relaxation_factor = stod(s_iter[0]);
		
		calcUpdatePrim.push_back(
		[id_u,id_v,id_w,relaxation_factor](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_u] += relaxation_factor*Xvalues[0];
			cells[id_v] += relaxation_factor*Xvalues[1];
			cells[id_w] += relaxation_factor*Xvalues[2];
		}); 
		
	}
	
	{
		int id_p = controls.getId_cellVar("pressure");
		int id_dp = controls.getId_cellVar("delta-pressure");
		
		// double relaxation_factor = stod(controls.fvSolutionMap["segregated"]);
		double relaxation_factor = stod(s_iter[1]);
		
		calcUpdatePrim.push_back(
		[id_p,id_dp,relaxation_factor](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_dp] = relaxation_factor*Xvalues[0];
			cells[id_p] += relaxation_factor*Xvalues[0];
		}); 
		
	}
	
	{
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_xgrad_dp = controls.getId_cellVar("x-gradient delta-pressure");
		int id_ygrad_dp = controls.getId_cellVar("y-gradient delta-pressure");
		int id_zgrad_dp = controls.getId_cellVar("z-gradient delta-pressure");
		int id_vol = controls.getId_cellVar("volume");
		int id_dt = controls.getId_fieldVar("time-step");
		
		// double relaxation_factor = stod(controls.fvSolutionMap["segregated"]);
		double relaxation_factor = stod(s_iter[2]);
		
		calcUpdatePrim.push_back(
		[id_u,id_v,id_w,relaxation_factor,id_xgrad_dp,id_ygrad_dp,id_zgrad_dp,
		id_vol,id_dt](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_u] += Xvalues[0];
			cells[id_v] += Xvalues[1];
			cells[id_w] += Xvalues[2];
			// cells[id_u] += 1.e-5/cells[id_vol]*Xvalues[0];
			// cells[id_v] += 1.e-5/cells[id_vol]*Xvalues[1];
			// cells[id_w] += 1.e-5/cells[id_vol]*Xvalues[2];
			// cells[id_u] -= fields[id_dt]*cells[id_xgrad_dp];
			// cells[id_v] -= fields[id_dt]*cells[id_ygrad_dp];
			// cells[id_w] -= fields[id_dt]*cells[id_zgrad_dp];
		}); 
		
	}
}
