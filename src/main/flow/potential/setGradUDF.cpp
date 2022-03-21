
#include "../../../others/solvers.h"

void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	gradLSIds_name.clear();
	gradLSIds_name.resize(2);
	
	// gradLSIds_name.back().push_back("pressure");
	// gradLSIds_name.back().push_back("x-velocity");
	// gradLSIds_name.back().push_back("y-velocity");
	// gradLSIds_name.back().push_back("z-velocity");
	// gradLSIds_name.back().push_back("temperature");

	// int nSp = controls.faceVar["left mass-fraction"].sub_name.size();
	// for(int i=0; i<nSp-1; ++i){
		// string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
		// gradLSIds_name.back().push_back("mass-fraction-"+tmp_name);
	// }
	
	
}