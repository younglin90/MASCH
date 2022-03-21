
#include "../../../../others/solvers.h"

void MASCH_Solver::setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// controls.nEq.resize(5+controls.nSp);
	
	controls.nEq.clear();
	controls.nEq.push_back(5+controls.nSp-1);
	
	
	controls.primIdsSegEq.clear();
	controls.primIdsSegEq.push_back(vector<int>());
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("pressure"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("x-velocity"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("y-velocity"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("z-velocity"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("temperature"));
	
	int nSp = controls.nSp;
	for(int i=0; i<nSp-1; ++i){
		string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
		controls.primIdsSegEq.back().push_back(controls.getId_cellVar("mass-fraction-"+tmp_name));
	}
	
}