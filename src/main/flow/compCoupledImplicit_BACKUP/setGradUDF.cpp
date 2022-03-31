
#include "../../../others/solvers.h"

void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	// gradLSIds_name.clear();
	gradLSIds_cell_name.resize(1);
	gradLSIds_bcFace_name.resize(1);
	
	gradLSIds_cell_name[0].push_back("pressure");
	gradLSIds_bcFace_name[0].push_back("pressure");
	
	gradLSIds_cell_name[0].push_back("x-velocity");
	gradLSIds_bcFace_name[0].push_back("x-velocity");
	
	gradLSIds_cell_name[0].push_back("y-velocity");
	gradLSIds_bcFace_name[0].push_back("y-velocity");
	
	gradLSIds_cell_name[0].push_back("z-velocity");
	gradLSIds_bcFace_name[0].push_back("z-velocity");
	
	gradLSIds_cell_name[0].push_back("temperature");
	gradLSIds_bcFace_name[0].push_back("temperature");
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		gradLSIds_cell_name[0].push_back(tmp_name);
		gradLSIds_bcFace_name[0].push_back(tmp_name);
	}
	
	gradLSIds_cell_name[0].push_back("density");
	gradLSIds_bcFace_name[0].push_back("density");
	
	
	
}