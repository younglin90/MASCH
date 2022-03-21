
#include "../../../others/solvers.h"

void MASCH_Solver::setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// controls.nEq.resize(5+controls.nSp);
	
	controls.nEq.clear();
	controls.sendProcValueNames.clear();
	checkImplicit.clear();
	
	int segSize = 3;
	
	controls.nEq.resize(segSize);
	controls.sendProcValueNames.resize(segSize);
	checkImplicit.resize(segSize);
	
	controls.nEq[0] = 3;
	controls.nEq[1] = 1;
	controls.nEq[2] = 3;
	
	checkImplicit[0] = true;
	checkImplicit[1] = true;
	checkImplicit[2] = false;
	
	controls.sendProcValueNames[0].push_back("x-velocity");
	controls.sendProcValueNames[0].push_back("y-velocity");
	controls.sendProcValueNames[0].push_back("z-velocity");
	
	controls.sendProcValueNames[1].push_back("delta-pressure");
	controls.sendProcValueNames[1].push_back("pressure");
	
	controls.sendProcValueNames[2].push_back("x-velocity");
	controls.sendProcValueNames[2].push_back("y-velocity");
	controls.sendProcValueNames[2].push_back("z-velocity");
}