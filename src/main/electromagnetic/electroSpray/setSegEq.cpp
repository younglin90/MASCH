
#include "../../../../others/solvers.h"

void MASCH_Solver::setSegEqUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// controls.nEq.resize(5+controls.nSp);
	
	controls.nEq.clear();
	controls.sendProcValueNames.clear();
	checkImplicit.clear();
	
	int segSize = 8;
	
	controls.nEq.resize(segSize);
	controls.sendProcValueNames.resize(segSize);
	checkImplicit.resize(segSize);
	
	controls.nEq[0] = controls.nSp-1; // volume-fraction eq.
	controls.nEq[1] = 1; // charge-conservation eq.
	controls.nEq[2] = 1; // electric-potential eq.
	controls.nEq[3] = 3; // electric-field eq.
	controls.nEq[4] = 3; // electric-force eq.
	controls.nEq[5] = 3; // navier-stokes momentum eq.
	controls.nEq[6] = 1; // pressure eq.
	controls.nEq[7] = 3; // velocity-correction eq.
	
	checkImplicit[0] = true;
	checkImplicit[1] = true;
	checkImplicit[2] = true;
	checkImplicit[3] = false;
	checkImplicit[4] = false;
	checkImplicit[5] = true;
	checkImplicit[6] = true;
	checkImplicit[7] = false;
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		controls.sendProcValueNames[0].push_back("volume-fraction-"+controls.spName[i]);
	}
	
	controls.sendProcValueNames[1].push_back("charge-density");
	
	controls.sendProcValueNames[2].push_back("electric-potential");
	
	controls.sendProcValueNames[3].push_back("x-electric-field");
	controls.sendProcValueNames[3].push_back("y-electric-field");
	controls.sendProcValueNames[3].push_back("z-electric-field");
	
	controls.sendProcValueNames[4].push_back("x-electric-force");
	controls.sendProcValueNames[4].push_back("y-electric-force");
	controls.sendProcValueNames[4].push_back("z-electric-force");
	
	controls.sendProcValueNames[5].push_back("x-velocity");
	controls.sendProcValueNames[5].push_back("y-velocity");
	controls.sendProcValueNames[5].push_back("z-velocity");
	
	controls.sendProcValueNames[6].push_back("pressure");
	controls.sendProcValueNames[6].push_back("delta-pressure");
	
	controls.sendProcValueNames[7].push_back("x-velocity");
	controls.sendProcValueNames[7].push_back("y-velocity");
	controls.sendProcValueNames[7].push_back("z-velocity");
	
}