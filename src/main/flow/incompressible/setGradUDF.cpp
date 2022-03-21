
#include "../../../others/solvers.h"

void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	// gradLSIds_name.clear();
	gradLSIds_cell_name.resize(3);
	gradLSIds_bcFace_name.resize(3);
	
	gradLSIds_cell_name[0].push_back("pressure");
	gradLSIds_bcFace_name[0].push_back("pressure");
	
	gradLSIds_cell_name[1].push_back("pressure");
	gradLSIds_bcFace_name[1].push_back("pressure");
	
	gradLSIds_cell_name[2].push_back("delta-pressure");
	gradLSIds_bcFace_name[2].push_back("delta-pressure");
	
}