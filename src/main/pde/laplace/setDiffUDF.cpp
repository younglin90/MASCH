
#include "../../../others/solvers.h"

void MASCH_Solver::setDiffFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		solver.calcLaplFlux.push_back(
		[](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			return 0;
		});
	}
	
	{
		solver.calcNLaplFlux.push_back(
		[](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			return 0;
		});
	}
	
	
}