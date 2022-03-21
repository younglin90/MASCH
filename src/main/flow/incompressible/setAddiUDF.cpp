
#include "../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	for(int i=0; i<3; ++i)
	{
		calcCellAddiVal.push_back(
			[&solver] (
			double* cells) ->int {
				return 0;
			}
		);
		calcFaceAddiVal.push_back( 
			[&solver] (
			double* faces) ->int {
				return 0;
			}
		);
	}
	
	
}

