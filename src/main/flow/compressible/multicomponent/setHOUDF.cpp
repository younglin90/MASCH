
#include "../../../../others/solvers.h"

void MASCH_Solver::setRecoFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// 재료들
	using id_type = unsigned short;
	vector<id_type> phiC, gradCx, gradCy, gradCz, phiFL, phiFR;
	
	phiC.push_back(controls.cellVar["pressure"].id); 
	gradCx.push_back(controls.cellVar["x-gradient pressure"].id); 
	gradCy.push_back(controls.cellVar["y-gradient pressure"].id);
	gradCz.push_back(controls.cellVar["z-gradient pressure"].id);
	phiFL.push_back(controls.faceVar["left pressure"].id); 
	phiFR.push_back(controls.faceVar["right pressure"].id); 
	
	phiC.push_back(controls.cellVar["x-velocity"].id); 
	gradCx.push_back(controls.cellVar["x-gradient x-velocity"].id); 
	gradCy.push_back(controls.cellVar["y-gradient x-velocity"].id);
	gradCz.push_back(controls.cellVar["z-gradient x-velocity"].id);
	phiFL.push_back(controls.faceVar["left x-velocity"].id); 
	phiFR.push_back(controls.faceVar["right x-velocity"].id);
	
	phiC.push_back(controls.cellVar["y-velocity"].id); 
	gradCx.push_back(controls.cellVar["x-gradient y-velocity"].id); 
	gradCy.push_back(controls.cellVar["y-gradient y-velocity"].id);
	gradCz.push_back(controls.cellVar["z-gradient y-velocity"].id);
	phiFL.push_back(controls.faceVar["left y-velocity"].id); 
	phiFR.push_back(controls.faceVar["right y-velocity"].id);
	
	phiC.push_back(controls.cellVar["z-velocity"].id); 
	gradCx.push_back(controls.cellVar["x-gradient z-velocity"].id); 
	gradCy.push_back(controls.cellVar["y-gradient z-velocity"].id);
	gradCz.push_back(controls.cellVar["z-gradient z-velocity"].id);
	phiFL.push_back(controls.faceVar["left z-velocity"].id); 
	phiFR.push_back(controls.faceVar["right z-velocity"].id);
	
	phiC.push_back(controls.cellVar["temperature"].id); 
	gradCx.push_back(controls.cellVar["x-gradient temperature"].id); 
	gradCy.push_back(controls.cellVar["y-gradient temperature"].id);
	gradCz.push_back(controls.cellVar["z-gradient temperature"].id);
	phiFL.push_back(controls.faceVar["left temperature"].id); 
	phiFR.push_back(controls.faceVar["right temperature"].id); 
	
	id_type id_xLR = 
	controls.faceVar["x distance of between left and right cell"].id;
	id_type id_yLR = 
	controls.faceVar["y distance of between left and right cell"].id;
	id_type id_zLR = 
	controls.faceVar["z distance of between left and right cell"].id;
	
	// for(auto& item : phiFL){
		// cout << item << endl;
	// }
	// for(auto& item : phiFR){
		// cout << item << endl;
	// }
	
	calcHO_FaceVal.push_back(
		[&solver,phiC,gradCx,gradCy,gradCz,phiFL,phiFR,
		id_xLR,id_yLR,id_zLR](
		double* cellsL, double* cellsR, double* faces)->int{
			for(int i=0, SIZE=phiC.size(); i<SIZE; ++i){
				int phiC_id = phiC[i];
				int gradCx_id = gradCx[i];
				int gradCy_id = gradCy[i];
				int gradCz_id = gradCz[i];
				int phiFL_id = phiFL[i];
				int phiFR_id = phiFR[i];
				double phiL1 = cellsL[phiC_id];
				double phiR1 = cellsR[phiC_id];
				double phiL2 = 0.0;
				phiL2 += cellsL[gradCx_id]*faces[id_xLR];
				phiL2 += cellsL[gradCy_id]*faces[id_yLR];
				phiL2 += cellsL[gradCz_id]*faces[id_zLR];
				phiL2 = phiR1 - 2.0*phiL2;
				double phiR2 = 0.0;
				phiR2 += cellsR[gradCx_id]*faces[id_xLR];
				phiR2 += cellsR[gradCy_id]*faces[id_yLR];
				phiR2 += cellsR[gradCz_id]*faces[id_zLR];
				phiR2 = phiL1 + 2.0*phiR2;
				// faces[phiFL_id] = (solver.NVD.Minmod(phiL2,phiL1,phiR1));
				// faces[phiFR_id] = (solver.NVD.Minmod(phiR2,phiR1,phiL1));
				faces[phiFL_id] = phiL1;
				faces[phiFR_id] = phiR1;
			}
			return 0;
		}
	);
}

