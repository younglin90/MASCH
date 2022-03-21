
#include "../../../others/solvers.h"

void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	saveCurrValueNames.push_back("x-velocity");
	saveOld1ValueNames.push_back("old x-velocity");
	// saveOld2ValueNames.push_back("old2 x-velocity");
	
	saveCurrValueNames.push_back("y-velocity");
	saveOld1ValueNames.push_back("old y-velocity");
	// saveOld2ValueNames.push_back("old2 y-velocity");
	
	saveCurrValueNames.push_back("z-velocity");
	saveOld1ValueNames.push_back("old z-velocity");
	// saveOld2ValueNames.push_back("old2 z-velocity");
	
}
