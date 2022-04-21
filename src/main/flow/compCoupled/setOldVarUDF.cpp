
#include "../../../others/solvers.h"

void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	saveCurrValueNames.push_back("pressure");
	saveOld1ValueNames.push_back("old pressure");
	
	saveCurrValueNames.push_back("x-velocity");
	saveOld1ValueNames.push_back("old x-velocity");
	
	saveCurrValueNames.push_back("y-velocity");
	saveOld1ValueNames.push_back("old y-velocity");
	
	saveCurrValueNames.push_back("z-velocity");
	saveOld1ValueNames.push_back("old z-velocity");
    
	saveCurrValueNames.push_back("temperature");
	saveOld1ValueNames.push_back("old temperature");
	
	saveCurrValueNames.push_back("total-enthalpy");
	saveOld1ValueNames.push_back("old total-enthalpy");
	
	saveCurrValueNames.push_back("density");
	saveOld1ValueNames.push_back("old density");
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		saveCurrValueNames.push_back(tmp_name);
		saveOld1ValueNames.push_back("old "+tmp_name);
	}
	
	
}
