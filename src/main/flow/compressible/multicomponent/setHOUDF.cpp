
#include "../../../../others/solvers.h"

void MASCH_Solver::setRecoFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// 재료들
	// using id_type = unsigned short;
	using id_type = int;
	vector<id_type> phiC, gradCx, gradCy, gradCz, phiFL, phiFR;
	
	int nSp = controls.nSp;
	
	phiC.push_back(controls.getId_cellVar("pressure")); 
	gradCx.push_back(controls.getId_cellVar("x-gradient pressure")); 
	gradCy.push_back(controls.getId_cellVar("y-gradient pressure"));
	gradCz.push_back(controls.getId_cellVar("z-gradient pressure"));
	phiFL.push_back(controls.getId_faceVar("left pressure")); 
	phiFR.push_back(controls.getId_faceVar("right pressure")); 
	
	phiC.push_back(controls.getId_cellVar("x-velocity")); 
	gradCx.push_back(controls.getId_cellVar("x-gradient x-velocity")); 
	gradCy.push_back(controls.getId_cellVar("y-gradient x-velocity"));
	gradCz.push_back(controls.getId_cellVar("z-gradient x-velocity"));
	phiFL.push_back(controls.getId_faceVar("left x-velocity")); 
	phiFR.push_back(controls.getId_faceVar("right x-velocity"));
	
	phiC.push_back(controls.getId_cellVar("y-velocity")); 
	gradCx.push_back(controls.getId_cellVar("x-gradient y-velocity")); 
	gradCy.push_back(controls.getId_cellVar("y-gradient y-velocity"));
	gradCz.push_back(controls.getId_cellVar("z-gradient y-velocity"));
	phiFL.push_back(controls.getId_faceVar("left y-velocity")); 
	phiFR.push_back(controls.getId_faceVar("right y-velocity"));
	
	phiC.push_back(controls.getId_cellVar("z-velocity")); 
	gradCx.push_back(controls.getId_cellVar("x-gradient z-velocity")); 
	gradCy.push_back(controls.getId_cellVar("y-gradient z-velocity"));
	gradCz.push_back(controls.getId_cellVar("z-gradient z-velocity"));
	phiFL.push_back(controls.getId_faceVar("left z-velocity")); 
	phiFR.push_back(controls.getId_faceVar("right z-velocity"));
	
	phiC.push_back(controls.getId_cellVar("temperature")); 
	gradCx.push_back(controls.getId_cellVar("x-gradient temperature")); 
	gradCy.push_back(controls.getId_cellVar("y-gradient temperature"));
	gradCz.push_back(controls.getId_cellVar("z-gradient temperature"));
	phiFL.push_back(controls.getId_faceVar("left temperature")); 
	phiFR.push_back(controls.getId_faceVar("right temperature")); 
	
	for(int i=0; i<nSp-1; ++i){
		string tmp_name = controls.faceVar["left mass-fraction"].sub_name[i];
		phiC.push_back(controls.getId_cellVar("mass-fraction-"+tmp_name)); 
		gradCx.push_back(controls.getId_cellVar("x-gradient mass-fraction-"+tmp_name)); 
		gradCy.push_back(controls.getId_cellVar("y-gradient mass-fraction-"+tmp_name));
		gradCz.push_back(controls.getId_cellVar("z-gradient mass-fraction-"+tmp_name));
		phiFL.push_back(controls.getId_faceVar("left mass-fraction-"+tmp_name));
		phiFR.push_back(controls.getId_faceVar("right mass-fraction-"+tmp_name));
	}
	
	id_type id_xLR = controls.getId_faceVar("x distance of between left and right cell");
	id_type id_yLR = controls.getId_faceVar("y distance of between left and right cell");
	id_type id_zLR = controls.getId_faceVar("z distance of between left and right cell");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	id_type id_area = controls.getId_faceVar("area");
	
	id_type id_dt = controls.getId_fieldVar("time-step");
	
	calcHO_FaceVal.push_back(
		[&solver,phiC,gradCx,gradCy,gradCz,phiFL,phiFR,
		id_xLR,id_yLR,id_zLR,id_area,id_dt,id_nx,id_ny,id_nz](
		double* fields, double* cellsL, double* cellsR, double* faces)->int{
		
			double pL, pR, uL, uR, vL, vR, wL, wR, TL, TR, cL, cR;
			pL = cellsL[phiC[0]]; pR = cellsR[phiC[0]];
			uL = cellsL[phiC[1]]; uR = cellsR[phiC[1]];
			vL = cellsL[phiC[2]]; vR = cellsR[phiC[2]];
			wL = cellsL[phiC[3]]; wR = cellsR[phiC[3]];
			TL = cellsL[phiC[4]]; TR = cellsR[phiC[4]];
			
			// faces[phiFL[0]] = pL; faces[phiFR[0]] = pR;
			// faces[phiFL[1]] = uL; faces[phiFR[1]] = uR;
			// faces[phiFL[2]] = vL; faces[phiFR[2]] = vR;
			// faces[phiFL[3]] = wL; faces[phiFR[3]] = wR;
			// faces[phiFL[4]] = TL; faces[phiFR[4]] = TR;
			
			double glob_dt = fields[id_dt];
			double dx = sqrt(faces[id_area]);
			double dt_tmp_L = dx/(sqrt(uL*uL+vL*vL+wL*wL)+cL);
			double dt_tmp_R = dx/(sqrt(uR*uR+vR*vR+wR*wR)+cR);			
			double coDD = max(glob_dt/dt_tmp_L,glob_dt/dt_tmp_R);
			
			for(int i=0; i<5; ++i){
				int id_phi = phiC[i];
				int id_phiL = phiFL[i]; int id_phiR = phiFR[i];
				int id_xGrad_phi = gradCx[i];
				int id_yGrad_phi = gradCy[i];
				int id_zGrad_phi = gradCz[i];
			
				double phiL1, phiR1; 
				double grad_phiL[3], grad_phiR[3]; 
				phiL1 = cellsL[id_phi]; phiR1 = cellsR[id_phi];
				grad_phiL[0] = cellsL[id_xGrad_phi]; grad_phiR[0] = cellsR[id_xGrad_phi];
				grad_phiL[1] = cellsL[id_yGrad_phi]; grad_phiR[1] = cellsR[id_yGrad_phi];
				grad_phiL[2] = cellsL[id_zGrad_phi]; grad_phiR[2] = cellsR[id_zGrad_phi];
				
				double phiL2 = 0.0;
				phiL2 += grad_phiL[0]*faces[id_xLR];
				phiL2 += grad_phiL[1]*faces[id_yLR];
				phiL2 += grad_phiL[2]*faces[id_zLR];
				phiL2 = phiR1 - 2.0*phiL2;
				double phiR2 = 0.0;
				phiR2 += grad_phiR[0]*faces[id_xLR];
				phiR2 += grad_phiR[1]*faces[id_yLR];
				phiR2 += grad_phiR[2]*faces[id_zLR];
				phiR2 = phiL1 + 2.0*phiR2;
				
				faces[id_phiL] = solver.NVD.Minmod(phiL2,phiL1,phiR1);
				faces[id_phiR] = solver.NVD.Minmod(phiL1,phiR1,phiR2);
			}
			
			// for(int i=5; i<6; ++i){
			{
				int i=5;
				
				int id_phi = phiC[i];
				int id_phiL = phiFL[i]; int id_phiR = phiFR[i];
				int id_xGrad_phi = gradCx[i];
				int id_yGrad_phi = gradCy[i];
				int id_zGrad_phi = gradCz[i];
			
				double phiL1, phiR1; 
				double grad_phiL[3], grad_phiR[3]; 
				phiL1 = cellsL[id_phi]; phiR1 = cellsR[id_phi];
				grad_phiL[0] = cellsL[id_xGrad_phi]; grad_phiR[0] = cellsR[id_xGrad_phi];
				grad_phiL[1] = cellsL[id_yGrad_phi]; grad_phiR[1] = cellsR[id_yGrad_phi];
				grad_phiL[2] = cellsL[id_zGrad_phi]; grad_phiR[2] = cellsR[id_zGrad_phi];
				
				double phiL2 = 0.0;
				phiL2 += grad_phiL[0]*faces[id_xLR];
				phiL2 += grad_phiL[1]*faces[id_yLR];
				phiL2 += grad_phiL[2]*faces[id_zLR];
				phiL2 = phiR1 - 2.0*phiL2;
				double phiR2 = 0.0;
				phiR2 += grad_phiR[0]*faces[id_xLR];
				phiR2 += grad_phiR[1]*faces[id_yLR];
				phiR2 += grad_phiR[2]*faces[id_zLR];
				phiR2 = phiL1 + 2.0*phiR2;
				
				phiL2 = max(0.0,min(1.0,phiL2));
				phiR2 = max(0.0,min(1.0,phiR2));
				
				double mfLR[3];
				mfLR[0] = 0.5*grad_phiL[0]+0.5*grad_phiR[0];
				mfLR[1] = 0.5*grad_phiL[1]+0.5*grad_phiR[1];
				mfLR[2] = 0.5*grad_phiL[2]+0.5*grad_phiR[2];
				double magMfLR = mfLR[0]*mfLR[0];
				magMfLR += mfLR[1]*mfLR[1];
				magMfLR += mfLR[2]*mfLR[2];
				magMfLR = sqrt(magMfLR);
				mfLR[0] = mfLR[0]/(magMfLR+1.e-200);
				mfLR[1] = mfLR[1]/(magMfLR+1.e-200);
				mfLR[2] = mfLR[2]/(magMfLR+1.e-200);
				double cosTheta = mfLR[0]*faces[id_nx];
				cosTheta += mfLR[1]*faces[id_ny];
				cosTheta += mfLR[2]*faces[id_nz];
				cosTheta = abs(cosTheta);
				double gamF = min(cosTheta*cosTheta*cosTheta*cosTheta,1.0);
				
				faces[id_phiL] = solver.NVD.getHO_MSTACS(phiL2,phiL1,phiR1,coDD,gamF);
				faces[id_phiR] = solver.NVD.getHO_MSTACS(phiL1,phiR1,phiR2,coDD,gamF);
				
			}
				
			return 0;
		}
	);
}

