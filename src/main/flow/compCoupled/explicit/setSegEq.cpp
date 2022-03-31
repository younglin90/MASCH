
#include "../../../../others/solvers.h"

void MASCH_Solver::setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// controls.nEq.resize(5+controls.nSp);
	
	controls.nEq.clear();
	controls.sendProcValueNames.clear();
	checkImplicit.clear();
	
	int segSize = 1;
	
	controls.nEq.resize(segSize);
	controls.sendProcValueNames.resize(segSize);
	checkImplicit.resize(segSize);
	
	controls.nEq[0] = 5+controls.nSp-1;
	
	// checkImplicit[0] = true;
	checkImplicit[0] = false;
	
	controls.sendProcValueNames[0].push_back("pressure");
	controls.sendProcValueNames[0].push_back("x-velocity");
	controls.sendProcValueNames[0].push_back("y-velocity");
	controls.sendProcValueNames[0].push_back("z-velocity");
	controls.sendProcValueNames[0].push_back("temperature");
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		controls.sendProcValueNames[0].push_back(tmp_name);
	}
	
}