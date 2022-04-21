
#include "../../../others/solvers.h"

void MASCH_Solver::setCurvatureFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	curvatureIds_cell_name.resize(8);
	
	// 곡률 관련
	for(int i=0; i<controls.spName.size(); ++i){
		string name = controls.spName[i];
		string type = controls.thermophysicalProperties[name+".transport.sigma.type"];
		if(type=="constant"){
			double value = stod(controls.thermophysicalProperties[name+".transport.sigma.value"]);
			if(value<1.e-200) continue;
			curvatureIds_cell_name[0].push_back(name);
		}
	}
	
	
}