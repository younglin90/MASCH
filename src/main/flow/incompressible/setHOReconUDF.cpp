
#include "../../../others/solvers.h"

void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pF = controls.getId_faceVar("pressure");

	int id_dp = controls.getId_cellVar("delta-pressure");
	int id_dpF = controls.getId_faceVar("delta-pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uF = controls.getId_faceVar("x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	
	{
		calcHO_FaceVal.push_back(
			[&solver,
			id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
			id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
			
				double dt = fields[id_dt];
				double rhoL = 1.0; double rhoR = 1.0;
				double dtrho = dt*0.5*(1.0/rhoL+1.0/rhoR);
				
				double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
				double pL = cellsL[id_p]; double pR = cellsR[id_p];
				double uL = cellsL[id_u]; double uR = cellsR[id_u];
				double vL = cellsL[id_v]; double vR = cellsR[id_v];
				double wL = cellsL[id_w]; double wR = cellsR[id_w];
				
				double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
				double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
				double UnF = 0.5*(UnL+UnR);
				UnF -= dtrho*(pR-pL)/dLR;
				UnF += dtrho*0.5*(dpdxL+dpdxR);
				
				faces[id_UnF] = UnF;
				
				if(UnF>=0.0){
					faces[id_uF] = uL; faces[id_vF] = vL; faces[id_wF] = wL;
				}
				else{
					faces[id_uF] = uR; faces[id_vF] = vR; faces[id_wF] = wR;
				}
				faces[id_pF] = 0.5*(pL+pR);
				
				return 0;
			}
		);
	}
	
	
	{
		calcHO_FaceVal.push_back(
			[&solver,
			id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
			id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
			
				double dt = fields[id_dt];
				double rhoL = 1.0; double rhoR = 1.0;
				double dtrho = dt*0.5*(1.0/rhoL+1.0/rhoR);
				
				double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
				double pL = cellsL[id_p]; double pR = cellsR[id_p];
				double uL = cellsL[id_u]; double uR = cellsR[id_u];
				double vL = cellsL[id_v]; double vR = cellsR[id_v];
				double wL = cellsL[id_w]; double wR = cellsR[id_w];
				
				double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
				double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
				double UnF = 0.5*(UnL+UnR);
				UnF -= dtrho*(pR-pL)/dLR;
				UnF += dtrho*0.5*(dpdxL+dpdxR);
				
				// faces[id_dtrho] = dtrho;
				faces[id_UnF] = UnF;
				
				return 0;
			}
		);
	}
	
	
	{
		calcHO_FaceVal.push_back(
			[&solver,
			id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
			id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				double dpL = cellsL[id_dp]; double dpR = cellsR[id_dp];
				faces[id_dpF] = 0.5*(dpL+dpR);
				
				return 0;
			}
		);
	}
	
}

