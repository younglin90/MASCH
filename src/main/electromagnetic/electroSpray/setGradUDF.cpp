
#include "../../../others/solvers.h"

void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	gradLSIds_cell_name.resize(8);
	gradLSIds_bcFace_name.resize(8);
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		gradLSIds_cell_name[0].push_back("volume-fraction-"+controls.spName[i]);
		gradLSIds_bcFace_name[0].push_back("volume-fraction-"+controls.spName[i]);
	}
	
	gradLSIds_cell_name[0].push_back("density");
	gradLSIds_bcFace_name[0].push_back("density"); 
	
	gradLSIds_cell_name[6].push_back("pressure");
	gradLSIds_bcFace_name[6].push_back("pressure");
	
	
	
}