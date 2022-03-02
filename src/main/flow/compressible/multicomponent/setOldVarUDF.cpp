
#include "../../../../others/solvers.h"

void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	// 올드값 관련
	saveCurrIds.push_back(controls.cellVar["pressure"].id);
	saveCurrIds.push_back(controls.cellVar["x-velocity"].id);
	saveCurrIds.push_back(controls.cellVar["y-velocity"].id);
	saveCurrIds.push_back(controls.cellVar["z-velocity"].id);
	saveCurrIds.push_back(controls.cellVar["density"].id);
	saveCurrIds.push_back(controls.cellVar["total enthalpy"].id);
	
	saveOld1Ids.push_back(controls.cellVar["old pressure"].id);
	saveOld1Ids.push_back(controls.cellVar["old x-velocity"].id);
	saveOld1Ids.push_back(controls.cellVar["old y-velocity"].id);
	saveOld1Ids.push_back(controls.cellVar["old z-velocity"].id);
	saveOld1Ids.push_back(controls.cellVar["old density"].id);
	saveOld1Ids.push_back(controls.cellVar["old total enthalpy"].id);
	
	saveOld2Ids.push_back(controls.cellVar["old2 pressure"].id);
	saveOld2Ids.push_back(controls.cellVar["old2 x-velocity"].id);
	saveOld2Ids.push_back(controls.cellVar["old2 y-velocity"].id);
	saveOld2Ids.push_back(controls.cellVar["old2 z-velocity"].id);
	saveOld2Ids.push_back(controls.cellVar["old2 density"].id);
	saveOld2Ids.push_back(controls.cellVar["old2 total enthalpy"].id);
}
