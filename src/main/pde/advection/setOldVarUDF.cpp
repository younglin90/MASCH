
#include "../../../others/solvers.h"

void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	saveCurrValueNames.push_back("pressure");
	saveOld1ValueNames.push_back("old pressure");
	saveOld2ValueNames.push_back("old2 pressure");
	
	saveCurrValueNames.push_back("x-velocity");
	saveOld1ValueNames.push_back("old x-velocity");
	saveOld2ValueNames.push_back("old2 x-velocity");
	
	saveCurrValueNames.push_back("y-velocity");
	saveOld1ValueNames.push_back("old y-velocity");
	saveOld2ValueNames.push_back("old2 y-velocity");
	
	saveCurrValueNames.push_back("z-velocity");
	saveOld1ValueNames.push_back("old z-velocity");
	saveOld2ValueNames.push_back("old2 z-velocity");
    
	saveCurrValueNames.push_back("temperature");
	saveOld1ValueNames.push_back("old temperature");
	saveOld2ValueNames.push_back("old2 temperature");
	
	saveCurrValueNames.push_back("total-enthalpy");
	saveOld1ValueNames.push_back("old total-enthalpy");
	saveOld2ValueNames.push_back("old2 total-enthalpy");
	
	saveCurrValueNames.push_back("density");
	saveOld1ValueNames.push_back("old density");
	saveOld2ValueNames.push_back("old2 density");
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		saveCurrValueNames.push_back(tmp_name);
		saveOld1ValueNames.push_back("old "+tmp_name);
		saveOld2ValueNames.push_back("old2 "+tmp_name);
	}
	
	
}
