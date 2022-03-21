
#include "../../../others/solvers.h"

void MASCH_Solver::setRecoFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		vector<int> phiC, phiFL, phiFR;
		
		phiC.push_back(controls.getId_cellVar("unknown")); 
		phiFL.push_back(controls.getId_faceVar("left unknown")); 
		phiFR.push_back(controls.getId_faceVar("right unknown"));
		
		// phiC.push_back(controls.getId_cellVar("delta-unknown")); 
		// phiFL.push_back(controls.getId_faceVar("left delta-unknown")); 
		// phiFR.push_back(controls.getId_faceVar("right delta-unknown"));
		
		calcHO_FaceVal.push_back(
			[&solver,phiC,phiFL,phiFR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				faces[phiFL[0]] = cellsL[phiC[0]]; faces[phiFR[0]] = cellsR[phiC[0]];
				// faces[phiFL[1]] = cellsL[phiC[1]]; faces[phiFR[1]] = cellsR[phiC[1]];
				
				return 0;
			}
		);
	}
	
}

