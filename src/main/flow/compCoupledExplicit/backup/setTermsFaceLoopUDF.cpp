
#include "../../../others/solvers.h"
// convective
void MASCH_Solver::setTermsFaceLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	int nSp = controls.spName.size();
	int nEq = 5+nSp-1; 
	
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pL = controls.getId_faceVar("left pressure");
	int id_pR = controls.getId_faceVar("right pressure");
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uL = controls.getId_faceVar("left x-velocity");
	int id_uR = controls.getId_faceVar("right x-velocity");
	int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity");
	int id_vR = controls.getId_faceVar("right y-velocity");
	int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity");
	int id_wR = controls.getId_faceVar("right z-velocity");
	int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	
	int id_T = controls.getId_cellVar("temperature");
	int id_TL = controls.getId_faceVar("left temperature");
	int id_TR = controls.getId_faceVar("right temperature");
	int id_dTdx = controls.getId_cellVar("x-gradient temperature");
	int id_dTdy = controls.getId_cellVar("y-gradient temperature");
	int id_dTdz = controls.getId_cellVar("z-gradient temperature");
	
	vector<int> id_Y, id_YL, id_YR, id_dYdx, id_dYdy, id_dYdz;
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		id_Y.push_back(controls.getId_cellVar(tmp_name));
		id_YL.push_back(controls.getId_faceVar("left mass-fraction-"+controls.spName[i]));
		id_YR.push_back(controls.getId_faceVar("right mass-fraction-"+controls.spName[i]));
		id_dYdx.push_back(controls.getId_faceVar("x-gradient "+tmp_name));
		id_dYdy.push_back(controls.getId_faceVar("y-gradient "+tmp_name));
		id_dYdz.push_back(controls.getId_faceVar("z-gradient "+tmp_name));
	}
	
	int id_rho = controls.getId_cellVar("density");
	int id_rhoL = controls.getId_faceVar("left density");
	int id_rhoR = controls.getId_faceVar("right density");
	
	int id_c = controls.getId_cellVar("speed-of-sound");
	int id_cL = controls.getId_faceVar("left speed-of-sound");
	int id_cR = controls.getId_faceVar("right speed-of-sound");
	
	int id_Ht = controls.getId_cellVar("total-enthalpy");
	int id_HtL = controls.getId_faceVar("left total-enthalpy");
	int id_HtR = controls.getId_faceVar("right total-enthalpy");
	
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	int id_wd = controls.getId_faceVar("distance weight"); 
	
	
	
	// 크랭크-니콜슨 방법 계수
	double CN_coeff = 0.5;
	
	


	{
		calcConvFlux.push_back(
		[&solver,nSp,
		id_p,id_pL,id_pR,id_u,id_uL,id_uR,id_v,id_vL,id_vR,id_w,id_wL,id_wR,
		id_c,id_cL,id_cR,id_Ht,id_HtL,id_HtR,id_rho,id_rhoL,id_rhoR,
		id_dpdx,id_dpdy,id_dpdz,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,
		id_wd,nEq,CN_coeff,id_Y,id_YL,id_YR,
		id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			
			// double pL = cellsL[id_p]; double pR = cellsR[id_p];
			// double uL = cellsL[id_u]; double uR = cellsR[id_u];
			// double vL = cellsL[id_v]; double vR = cellsR[id_v];
			// double wL = cellsL[id_w]; double wR = cellsR[id_w];
			// // double TL = cellsL[id_T]; double TR = cellsR[id_T];
			// double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
			// double cL = cellsL[id_c]; double cR = cellsR[id_c];
			// double HtL = cellsL[id_Ht]; double HtR = cellsR[id_Ht];
			// double YL[nSp], YR[nSp]; 
			// for(int i=0; i<nSp-1; ++i){
				// YL[i] = cellsL[id_Y[i]]; YR[i] = cellsR[id_Y[i]];
			// }
			
			double pL = faces[id_pL]; double pR = faces[id_pR];
			double uL = faces[id_uL]; double uR = faces[id_uR];
			double vL = faces[id_vL]; double vR = faces[id_vR];
			double wL = faces[id_wL]; double wR = faces[id_wR];
			double rhoL = faces[id_rhoL]; double rhoR = faces[id_rhoR];
			double cL = faces[id_cL]; double cR = faces[id_cR];
			double HtL = faces[id_HtL]; double HtR = faces[id_HtR];
			double YL[nSp], YR[nSp]; 
			for(int i=0; i<nSp-1; ++i){
				YL[i] = faces[id_YL[i]]; YR[i] = faces[id_YR[i]];
			}
			

			double nvec[3];
			nvec[0] = faces[id_nx];
			nvec[1] = faces[id_ny];
			nvec[2] = faces[id_nz];
			double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
			double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
			
			double Unhat = 0.5*(UnL+UnR);
			double rhohat = 0.5*(rhoL+rhoR);
			double chat= 0.5*(cL+cR);
			
			double ML = UnL/chat; 
			double MR = UnR/chat;
			double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
			double MLP = 0.5*(ML+abs(ML));
			if( abs(ML) < 1.0 ) {
				MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			}
			double MRM = 0.5*(MR-abs(MR));
			if( abs(MR) < 1.0 ) {
				MRM = -0.25*(MR - 1.0)*(MR - 1.0);
			}
			double preP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			if( abs(ML) < 1.0 ) {
				preP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			} 
			double preM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			if( abs(MR) < 1.0 ) {
				preM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
			} 
			
			double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );

			// SLAU
			double Mcy = min(1.0,KLR/chat);
			double phi_c = (1.0-Mcy)*(1.0-Mcy);
			double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			double D_L = ML+(1.0-g_c)*abs(ML);
			double D_R = MR-(1.0-g_c)*abs(MR);
			double D_rho = Mbar*g_c;
			double MLPL = 0.5*(D_L+D_rho);
			double MRMR = 0.5*(D_R-D_rho);

			// double mdot = rhoFL*chat*MLPL + rhoFR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);
			double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);

			double f1L = mdot;
			double f1R = 0.0;
			if( mdot < 0.0 ) {
				f1L = 0.0; f1R = mdot;
			}

			double PLR = 0.5*(pL+pR) - 
						 0.5*Mcy*preP*preM*0.5*(pL+pR)/chat*(UnR-UnL) + 
						 Mcy*0.5*(pL+pR)*(preP+preM-1.0) - 
						 0.5*(preP-preM)*(pR-pL);
			// double area = var.faces[i][id_area];
			
			double area = faces[id_area];
			int iter=0;
			fluxB[iter++] = -(f1L+f1R)*area;
			fluxB[iter++] = -(f1L*uL + f1R*uR + PLR*nvec[0])*area;
			fluxB[iter++] = -(f1L*vL + f1R*vR + PLR*nvec[1])*area;
			fluxB[iter++] = -(f1L*wL + f1R*wR + PLR*nvec[2])*area;
			fluxB[iter++] = -(f1L*HtL+ f1R*HtR)*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = -(f1L*YL[i]+ f1R*YR[i])*area;
			}
			
			
		}); 
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	calcConvFlux_BC.resize(1);
	
	using conv_diff_Funct_BC_type = 
		function<int(
		double* fields, double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB)>;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		
		{
			int bc_id = 0;
			calcConvFlux_BC[bc_id].push_back(vector<conv_diff_Funct_BC_type>());
			
			vector<bool> bool_zeroGradient(nEq,false);
			if(controls.boundaryMap["pressure"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			if(controls.boundaryMap["velocity"][bcName+".type"]=="zeroGradient") bool_zeroGradient[1] = true;
			if(controls.boundaryMap["temperature"][bcName+".type"]=="zeroGradient") bool_zeroGradient[2] = true;
			for(int i=0; i<controls.spName.size()-1; ++i){
				string tmp_name = ("mass-fraction-"+controls.spName[i]);
				if(controls.boundaryMap[tmp_name][bcName+".type"]=="zeroGradient") bool_zeroGradient[3+i] = true;
			}
			
			
			calcConvFlux_BC[bc_id].back().push_back(
			[&solver,nSp,
			id_p,id_pL,id_pR,id_u,id_uL,id_uR,id_v,id_vL,id_vR,id_w,id_wL,id_wR,
			id_c,id_cL,id_cR,id_Ht,id_HtL,id_HtR,id_rho,id_rhoL,id_rhoR,
			id_dpdx,id_dpdy,id_dpdz,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,
			id_wd,nEq,CN_coeff,id_Y,
			id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB) ->int {
				
					
				double pL = faces[id_pL]; double pR = pL;
				double uL = faces[id_uL]; double uR = uL;
				double vL = faces[id_vL]; double vR = vL;
				double wL = faces[id_wL]; double wR = wL;
				double rhoL = faces[id_rhoL]; double rhoR = rhoL;
				double cL = faces[id_cL]; double cR = cL;
				double HtL = faces[id_HtL]; double HtR = HtL;
				double YL[nSp], YR[nSp]; 
				for(int i=0; i<nSp-1; ++i){
					YL[i] = cellsL[id_Y[i]]; YR[i] = YL[i];
				}
				
				// double YFL[nSp], YFR[nSp]; 
				// double uFL, uFR, vFL, vFR, wFL, wFR;
				// double rhoFL, rhoFR, HtFL, HtFR;
				// uFL = faces[id_uL]; uFR = faces[id_uR];
				// vFL = faces[id_vL]; vFR = faces[id_vR];
				// wFL = faces[id_wL]; wFR = faces[id_wR];
				// rhoFL = faces[id_rhoL]; rhoFR = faces[id_rhoR];
				// HtFL = faces[id_HtL]; HtFR = faces[id_HtR];
				// for(int i=0; i<nSp-1; ++i){
					// YFL[i] = faces[id_YL[i]]; YFR[i] = faces[id_YR[i]];
				// }
				

				double nvec[3];
				nvec[0] = faces[id_nx];
				nvec[1] = faces[id_ny];
				nvec[2] = faces[id_nz];
				double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
				double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
				
				double Unhat = 0.5*(UnL+UnR);
				double rhohat = 0.5*(rhoL+rhoR);
				double chat= 0.5*(cL+cR);
				
				double ML = UnL/chat; 
				double MR = UnR/chat;
				double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
				double MLP = 0.5*(ML+abs(ML));
				if( abs(ML) < 1.0 ) {
					MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				}
				double MRM = 0.5*(MR-abs(MR));
				if( abs(MR) < 1.0 ) {
					MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				}
				double preP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
				if( abs(ML) < 1.0 ) {
					preP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
				} 
				double preM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
				if( abs(MR) < 1.0 ) {
					preM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
				} 
				
				double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );

				// SLAU
				double Mcy = min(1.0,KLR/chat);
				double phi_c = (1.0-Mcy)*(1.0-Mcy);
				double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
				double D_L = ML+(1.0-g_c)*abs(ML);
				double D_R = MR-(1.0-g_c)*abs(MR);
				double D_rho = Mbar*g_c;
				double MLPL = 0.5*(D_L+D_rho);
				double MRMR = 0.5*(D_R-D_rho);

				// double mdot = rhoFL*chat*MLPL + rhoFR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);
				double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);

				double f1L = mdot;
				double f1R = 0.0;
				if( mdot < 0.0 ) {
					f1L = 0.0; f1R = mdot;
				}

				double PLR = 0.5*(pL+pR) - 
							 0.5*Mcy*preP*preM*0.5*(pL+pR)/chat*(UnR-UnL) + 
							 Mcy*0.5*(pL+pR)*(preP+preM-1.0) - 
							 0.5*(preP-preM)*(pR-pL);
				// double area = var.faces[i][id_area];
				
				double area = faces[id_area];
				int iter=0;
				// fluxB[iter++] = -(f1L+f1R)*area;
				// fluxB[iter++] = -(f1L*uFL + f1R*uFR + PLR*nvec[0])*area;
				// fluxB[iter++] = -(f1L*vFL + f1R*vFR + PLR*nvec[1])*area;
				// fluxB[iter++] = -(f1L*wFL + f1R*wFR + PLR*nvec[2])*area;
				// fluxB[iter++] = -(f1L*HtFL+ f1R*HtFR)*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[iter++] = -(f1L*YFL[i]+ f1R*YFR[i])*area;
				// }
				fluxB[iter++] = -(f1L+f1R)*area;
				fluxB[iter++] = -(f1L*uL + f1R*uR + PLR*nvec[0])*area;
				fluxB[iter++] = -(f1L*vL + f1R*vR + PLR*nvec[1])*area;
				fluxB[iter++] = -(f1L*wL + f1R*wR + PLR*nvec[2])*area;
				fluxB[iter++] = -(f1L*HtL+ f1R*HtR)*area;
				for(int i=0; i<nSp-1; ++i){
					fluxB[iter++] = -(f1L*YL[i]+ f1R*YR[i])*area;
				}
				
				
				return 0;
			});
		}
		
	}
		
	
}