
#include "../../../others/solvers.h"

void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	saveCurrValueNames.push_back("x-velocity");
	saveOld1ValueNames.push_back("old x-velocity");
	
	saveCurrValueNames.push_back("y-velocity");
	saveOld1ValueNames.push_back("old y-velocity");
	
	saveCurrValueNames.push_back("z-velocity");
	saveOld1ValueNames.push_back("old z-velocity");
	
	saveCurrValueNames.push_back("charge-density");
	saveOld1ValueNames.push_back("old charge-density");
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		saveCurrValueNames.push_back("volume-fraction-"+controls.spName[i]);
		saveOld1ValueNames.push_back("old volume-fraction-"+controls.spName[i]);
	}
	
	
}
