
#include "../../../../others/solvers.h"

void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	
	gradLSIds_name.push_back("pressure");
	gradLSIds_name.push_back("x-velocity");
	gradLSIds_name.push_back("y-velocity");
	gradLSIds_name.push_back("z-velocity");
	
}