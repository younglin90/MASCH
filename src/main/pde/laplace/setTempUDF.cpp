
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTempFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	{
		solver.calcTemporal.push_back(
		[](double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			fluxA[0] = 0.0;
			fluxB[0] = 0.0;
			return 0;
		}); 
	}
	
}