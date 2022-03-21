
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		int id_p = controls.getId_cellVar("unknown");
		
		double relaxation_factor = stod(controls.mapArgv["-r"]);
		
		calcUpdatePrim.push_back(
		[id_p,relaxation_factor](
		double* cells, double* Xvalues) ->int {
			cells[id_p] += relaxation_factor*Xvalues[0];
		}); 
		
	}
}
