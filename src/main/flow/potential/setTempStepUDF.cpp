
#include "../../../others/solvers.h"


void MASCH_Solver::setTempStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		solver.calcTempStepCell.push_back(
		[](
		double* cells, double* fields) ->int {
			return 0;
		});
	}
	{
		solver.calcTempStepCell.push_back(
		[](
		double* cells, double* fields) ->int {
			return 0;
		});
	}
	
	
}
