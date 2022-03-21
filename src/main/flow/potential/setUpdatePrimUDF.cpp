
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		int id_p = controls.getId_cellVar("pressure");
		int id_dp = controls.getId_cellVar("delta-pressure");
		
		double relaxation_factor = stod(controls.mapArgv["-r"]);
		
		calcUpdatePrim.push_back(
		[id_p,id_dp,relaxation_factor](
		double* cells, double* Xvalues) ->int {
			cells[id_dp] = relaxation_factor*Xvalues[0];
			cells[id_p] += relaxation_factor*Xvalues[0];
		}); 
		
	}
	
	{
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_rho = controls.getId_cellVar("density");
		
		calcUpdatePrim.push_back(
		[id_u,id_v,id_w,id_rho](
		double* cells, double* Xvalues) ->int {
			cells[id_u] -= 1.0/cells[id_rho]*Xvalues[0];
			cells[id_v] -= 1.0/cells[id_rho]*Xvalues[1];
			cells[id_w] -= 1.0/cells[id_rho]*Xvalues[2];
		}); 
		
	}
	
}
