
#include "../../../../others/solvers.h"
// convective
void MASCH_Solver::setConvFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// checkImplicitConvFlux.push_back(false);
	
	{
		// using US = unsigned short;
		// int nSp = controls.cellVar["mass fraction"].sub_name.size();
		int nSp = controls.nSp;
		
		// 셀값
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
			id_Y[i] = controls.getId_cellVar("mass-fraction-"+tmp_name);
		}
		int id_rho = controls.getId_cellVar("density");
		int id_c = controls.getId_cellVar("speed-of-sound");
		int id_Ht = controls.getId_cellVar("total-enthalpy");
		
		// left값
		int id_pL = controls.getId_faceVar("left pressure");
		int id_uL = controls.getId_faceVar("left x-velocity");
		int id_vL = controls.getId_faceVar("left y-velocity");
		int id_wL = controls.getId_faceVar("left z-velocity");
		int id_TL = controls.getId_faceVar("left temperature");
		vector<int> id_YL(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.faceVar["left mass-fraction"].sub_name[i];
			id_YL[i] = controls.getId_faceVar("left mass-fraction-"+tmp_name);
		}
		int id_rhoL = controls.getId_faceVar("left density");
		int id_cL = controls.getId_faceVar("left speed-of-sound");
		int id_HtL = controls.getId_faceVar("left total-enthalpy");
		
		// right값
		int id_pR = controls.getId_faceVar("right pressure");
		int id_uR = controls.getId_faceVar("right x-velocity");
		int id_vR = controls.getId_faceVar("right y-velocity");
		int id_wR = controls.getId_faceVar("right z-velocity");
		int id_TR = controls.getId_faceVar("right temperature");
		vector<int> id_YR(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.faceVar["right mass-fraction"].sub_name[i];
			id_YR[i] = controls.getId_faceVar("right mass-fraction-"+tmp_name);
		}
		int id_rhoR = controls.getId_faceVar("right density");
		int id_cR = controls.getId_faceVar("right speed-of-sound");
		int id_HtR = controls.getId_faceVar("right total-enthalpy");
		
		// 메쉬관련
		int id_nx = controls.getId_faceVar("x unit normal");
		int id_ny = controls.getId_faceVar("y unit normal");
		int id_nz = controls.getId_faceVar("z unit normal");
		int id_area = controls.getId_faceVar("area");
		
		calcConvFlux.push_back(
		[nSp,id_nx,id_ny,id_nz,id_area,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
		id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
		id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			double pL = cellsL[id_p]; double pR = 0.0;
			double uL = cellsL[id_u]; double uR = 0.0;
			double vL = cellsL[id_v]; double vR = 0.0;
			double wL = cellsL[id_w]; double wR = 0.0;
			double TL = cellsL[id_T]; double TR = 0.0;
			double rhoL = cellsL[id_rho]; double rhoR = 0.0;
			double cL = cellsL[id_c]; double cR = 0.0;
			double HtL = cellsL[id_Ht]; double HtR = 0.0;
			
			if(cellsR!=nullptr){
				pR = cellsR[id_p];
				uR = cellsR[id_u];
				vR = cellsR[id_v];
				wR = cellsR[id_w];
				TR = cellsR[id_T];
				rhoR = cellsR[id_rho];
				cR = cellsR[id_c];
				HtR = cellsR[id_Ht];
			}
			else{
				pL = faces[id_pL];
				uL = faces[id_uL];
				vL = faces[id_vL];
				wL = faces[id_wL];
				TL = faces[id_TL];
				rhoL = faces[id_rhoL];
				cL = faces[id_cL];
				HtL = faces[id_HtL];
				
				pR = faces[id_pR];
				uR = faces[id_uR];
				vR = faces[id_vR];
				wR = faces[id_wR];
				TR = faces[id_TR];
				rhoR = faces[id_rhoR];
				cR = faces[id_cR];
				HtR = faces[id_HtR];
			}
				
			double YFL[nSp], YFR[nSp]; 
			double uFL, uFR, vFL, vFR, wFL, wFR;
			double rhoFL, rhoFR, HtFL, HtFR;
			uFL = faces[id_uL]; uFR = faces[id_uR];
			vFL = faces[id_vL]; vFR = faces[id_vR];
			wFL = faces[id_wL]; wFR = faces[id_wR];
			rhoFL = faces[id_rhoL]; rhoFR = faces[id_rhoR];
			HtFL = faces[id_HtL]; HtFR = faces[id_HtR];
			for(int i=0; i<nSp-1; ++i){
				YFL[i] = faces[id_YL[i]]; YFR[i] = faces[id_YR[i]];
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

			double mdot = rhoFL*chat*MLPL + rhoFR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);

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
			fluxB[iter++] = -(f1L*uFL + f1R*uFR + PLR*nvec[0])*area;
			fluxB[iter++] = -(f1L*vFL + f1R*vFR + PLR*nvec[1])*area;
			fluxB[iter++] = -(f1L*wFL + f1R*wFR + PLR*nvec[2])*area;
			fluxB[iter++] = -(f1L*HtFL+ f1R*HtFR)*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = -(f1L*YFL[i]+ f1R*YFR[i])*area;
			}
		}); 
	}
	
	
	
	
	
}