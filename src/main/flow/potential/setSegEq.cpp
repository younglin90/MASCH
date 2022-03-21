
#include "../../../others/solvers.h"

void MASCH_Solver::setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// controls.nEq.resize(5+controls.nSp);
	
	controls.nEq.clear();
	controls.nEq.push_back(1);
	controls.nEq.push_back(3);
	
	checkImplicit.push_back(true);
	checkImplicit.push_back(false);
	
	
	controls.primIdsSegEq.clear();
	controls.primIdsSegEq.push_back(vector<int>());
	// controls.primIdsSegEq.back().push_back(controls.getId_cellVar("pressure"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("delta-pressure"));
	
	controls.primIdsSegEq.push_back(vector<int>());
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("x-velocity"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("y-velocity"));
	controls.primIdsSegEq.back().push_back(controls.getId_cellVar("z-velocity"));
	
}