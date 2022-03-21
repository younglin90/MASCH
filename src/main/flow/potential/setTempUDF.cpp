
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTempFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	{
		solver.calcTemporal.push_back(
		[](double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			fluxA[0] = 1.0;
			
			fluxB[0] = 0.0;
			
			return 0;
		}); 
	}
	{
		int id_vol = controls.getId_cellVar("volume");
		solver.calcTemporal.push_back(
		[id_vol](double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			double volume = cells[id_vol];
			fluxA[0] = volume;
			fluxA[1] = 0.0;
			fluxA[2] = 0.0;
			
			fluxA[3] = 0.0;
			fluxA[4] = volume;
			fluxA[5] = 0.0;
			
			fluxA[6] = 0.0;
			fluxA[7] = 0.0;
			fluxA[8] = volume;
			
			fluxB[0] = 0.0;
			fluxB[1] = 0.0;
			fluxB[2] = 0.0;
			
			return 0;
		}); 
	}
	
	
}