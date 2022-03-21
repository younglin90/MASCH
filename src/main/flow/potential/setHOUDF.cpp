
#include "../../../others/solvers.h"

void MASCH_Solver::setRecoFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		vector<int> phiC, phiFL, phiFR;
		
		phiC.push_back(controls.getId_cellVar("x-velocity")); 
		phiFL.push_back(controls.getId_faceVar("left x-velocity")); 
		phiFR.push_back(controls.getId_faceVar("right x-velocity"));
		
		phiC.push_back(controls.getId_cellVar("y-velocity")); 
		phiFL.push_back(controls.getId_faceVar("left y-velocity")); 
		phiFR.push_back(controls.getId_faceVar("right y-velocity"));
		
		phiC.push_back(controls.getId_cellVar("z-velocity")); 
		phiFL.push_back(controls.getId_faceVar("left z-velocity")); 
		phiFR.push_back(controls.getId_faceVar("right z-velocity"));
		
		calcHO_FaceVal.push_back(
			[&solver,phiC,phiFL,phiFR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
			
				double uL, uR, vL, vR, wL, wR;
				uL = cellsL[phiC[0]]; uR = cellsR[phiC[0]];
				vL = cellsL[phiC[1]]; vR = cellsR[phiC[1]];
				wL = cellsL[phiC[2]]; wR = cellsR[phiC[2]];
				
				faces[phiFL[0]] = uL; faces[phiFR[0]] = uR;
				faces[phiFL[1]] = vL; faces[phiFR[1]] = vR;
				faces[phiFL[2]] = wL; faces[phiFR[2]] = wR;
				
				return 0;
			}
		);
	}
	
	
	
	
	
	{
		vector<int> phiC, phiFL, phiFR;
		
		phiC.push_back(controls.getId_cellVar("pressure")); 
		phiFL.push_back(controls.getId_faceVar("left pressure")); 
		phiFR.push_back(controls.getId_faceVar("right pressure")); 
		
		phiC.push_back(controls.getId_cellVar("delta-pressure")); 
		phiFL.push_back(controls.getId_faceVar("left delta-pressure")); 
		phiFR.push_back(controls.getId_faceVar("right delta-pressure")); 
		
		calcHO_FaceVal.push_back(
			[&solver,phiC,phiFL,phiFR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
			
				double dpL = cellsL[phiC[1]]; double dpR = cellsR[phiC[1]];
				
				// faces[phiFL[0]] = pL; faces[phiFR[0]] = pR;
				faces[phiFL[1]] = dpL; faces[phiFR[1]] = dpR;
				
				return 0;
			}
		);
	}
}

