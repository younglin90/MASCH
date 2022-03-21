
#include "../../../others/solvers.h"


void MASCH_Solver::setTimeStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	for(int i=0; i<3; ++i)
	{
		int id_dt = controls.getId_fieldVar("time-step");
		double dt = controls.timeStep;
		solver.calcTempStepCell.push_back(
		[id_dt,dt](
		double* cells, double* fields) ->int {
			fields[id_dt] = dt;
			return 0;
		});
	}
	
}
