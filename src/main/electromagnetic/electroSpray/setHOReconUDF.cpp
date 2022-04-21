

#include "../../../../others/solvers.h"

void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	int nSp = controls.spName.size();
	
	int nSp = controls.spName.size();
	int nEq = 5+nSp-1;
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_p = controls.getId_cellVar("pressure");
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	int id_dp = controls.getId_cellVar("delta-pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	vector<int> id_alpha, id_alphaL, id_alphaR;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+controls.spName[i]));
		id_alphaL.push_back(controls.getId_cellVar("left volume-fraction-"+controls.spName[i]));
		id_alphaR.push_back(controls.getId_cellVar("right volume-fraction-"+controls.spName[i]));
	}
	int id_mu = controls.getId_cellVar("viscosity");
	int id_rho = controls.getId_cellVar("density");
	
	int id_rhoe = controls.getId_cellVar("charge-density");
	int id_phi = controls.getId_cellVar("electric-potential");
	int id_xE = controls.getId_cellVar("x-electric-field");
	int id_yE = controls.getId_cellVar("y-electric-field");
	int id_zE = controls.getId_cellVar("z-electric-field");
	int id_xFe = controls.getId_cellVar("x-electric-force");
	int id_yFe = controls.getId_cellVar("y-electric-force");
	int id_zFe = controls.getId_cellVar("z-electric-force");
	int id_k = controls.getId_cellVar("conductivity");
	int id_epsilon = controls.getId_cellVar("permittivity");
	
	// 최대 최소값
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	int id_wd = controls.getId_faceVar("distance weight"); 
	int id_alpha = controls.getId_faceVar("cosine angle of between face normal and cells"); 
	int id_nLRx = controls.getId_faceVar("x unit normal of between left and right cell");
	int id_nLRy = controls.getId_faceVar("y unit normal of between left and right cell");
	int id_nLRz = controls.getId_faceVar("z unit normal of between left and right cell");
	int id_xLF = controls.getId_faceVar("x distance of between left cell and face");
	int id_yLF = controls.getId_faceVar("y distance of between left cell and face");
	int id_zLF = controls.getId_faceVar("z distance of between left cell and face");
	int id_xRF = controls.getId_faceVar("x distance of between right cell and face");
	int id_yRF = controls.getId_faceVar("y distance of between right cell and face");
	int id_zRF = controls.getId_faceVar("z distance of between right cell and face");
	int id_xLNv = controls.getId_faceVar("left cell to face normal vector shortest x distance");
	int id_yLNv = controls.getId_faceVar("left cell to face normal vector shortest y distance");
	int id_zLNv = controls.getId_faceVar("left cell to face normal vector shortest z distance");
	int id_xRNv = controls.getId_faceVar("right cell to face normal vector shortest x distance");
	int id_yRNv = controls.getId_faceVar("right cell to face normal vector shortest y distance");
	int id_zRNv = controls.getId_faceVar("right cell to face normal vector shortest z distance");
	
    
	int id_Un = controls.getId_faceVar("contravariant-velocity");
    
	
	{
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
                
                double dt = fields[id_dt];
                
                double nvec[3];
                nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
                double area = faces[id_area];
                double dLR = faces[id_dLR]; 
                double dAlpha = faces[id_alpha];
                double nLR[3];
                nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
                
                double pL = cellsL[id_p]; double pR = cellsR[id_p];
                double uL = cellsL[id_u]; double uR = cellsR[id_u];
                double vL = cellsL[id_v]; double vR = cellsR[id_v];
                double wL = cellsL[id_w]; double wR = cellsR[id_w];
                double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
                double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
                double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
                double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
                // double rhoFL = faces[id_rhoL]; double rhoFR = faces[id_rhoR];
                
                double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
                double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
                double UnF = 0.5*(UnL+UnR);
                
                // Rhie-Chow interpolation
                UnF -= dt*0.5*(1.0/rhoL+1.0/rhoR)*(pR-pL)/dLR;
                UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
                UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
                UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
			
                faces[id_Un] = UnF;
                
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
                
                double dt = fields[id_dt];
                
                double nvec[3];
                nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
                double area = faces[id_area];
                double dLR = faces[id_dLR]; 
                double dAlpha = faces[id_alpha];
                double nLR[3];
                nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
                
                double pL = cellsL[id_p]; double pR = cellsR[id_p];
                double uL = cellsL[id_u]; double uR = cellsR[id_u];
                double vL = cellsL[id_v]; double vR = cellsR[id_v];
                double wL = cellsL[id_w]; double wR = cellsR[id_w];
                double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
                double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
                double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
                double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
                // double rhoFL = faces[id_rhoL]; double rhoFR = faces[id_rhoR];
                
                double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
                double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
                double UnF = 0.5*(UnL+UnR);
                
                // Rhie-Chow interpolation
                UnF -= dt*0.5*(1.0/rhoL+1.0/rhoR)*(pR-pL)/dLR;
                UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
                UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
                UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
			
                faces[id_Un] = UnF;
                
				return 0;
			}
		);
        
        
		calcHO_FaceVal.push_back(
			[](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
                
                double dt = fields[id_dt];
                
                double nvec[3];
                nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
                double area = faces[id_area];
                double dLR = faces[id_dLR]; 
                double dAlpha = faces[id_alpha];
                double nLR[3];
                nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
                
                double pL = cellsL[id_p]; double pR = cellsR[id_p];
                double uL = cellsL[id_u]; double uR = cellsR[id_u];
                double vL = cellsL[id_v]; double vR = cellsR[id_v];
                double wL = cellsL[id_w]; double wR = cellsR[id_w];
                double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
                double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
                double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
                double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
                // double rhoFL = faces[id_rhoL]; double rhoFR = faces[id_rhoR];
                
                double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
                double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
                double UnF = 0.5*(UnL+UnR);
                
                // Rhie-Chow interpolation
                UnF -= dt*0.5*(1.0/rhoL+1.0/rhoR)*(pR-pL)/dLR;
                UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
                UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
                UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
			
                faces[id_Un] = UnF;
                
				return 0;
			}
		);
	}
	// cout << "BBBBBBBB" << endl;
	
}


