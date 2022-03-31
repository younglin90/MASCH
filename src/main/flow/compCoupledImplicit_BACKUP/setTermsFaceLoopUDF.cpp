
#include "../../../others/solvers.h"
// convective
void MASCH_Solver::setTermsFaceLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	int nSp = controls.spName.size();
	int nEq = 5+nSp-1;
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pF = controls.getId_faceVar("pressure");

	// int id_dp = controls.getId_cellVar("delta-pressure");
	// int id_dpF = controls.getId_faceVar("delta-pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uF = controls.getId_faceVar("x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	
	vector<int> id_Y, id_YF;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_Y.push_back(controls.getId_cellVar("mass-fraction-"+controls.spName[i]));
		id_YF.push_back(controls.getId_faceVar("mass-fraction-"+controls.spName[i]));
	}
	
	int id_muF = controls.getId_faceVar("viscosity");
	
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	
	int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	
	int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	
	int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	int id_wd = controls.getId_faceVar("distance weight"); 
	
	int id_rho = controls.getId_cellVar("density");
	int id_c = controls.getId_cellVar("speed-of-sound");
	int id_Ht = controls.getId_cellVar("total-enthalpy");
	int id_rhoF = controls.getId_faceVar("density");
	int id_HtF = controls.getId_faceVar("total-enthalpy");
	// int id_cF = controls.getId_faceVar("speed-of-sound");
	
	
	// int id_rho = controls.getId_cellVar("density");
	int id_drhodp = controls.getId_cellVar("partial-density-pressure");
	int id_drhodT = controls.getId_cellVar("partial-density-temperature");
	int id_dHtdp = controls.getId_cellVar("partial-total-enthalpy-pressure");
	int id_dHtdT = controls.getId_cellVar("partial-total-enthalpy-temperature");
	vector<int> id_drhodY, id_dHtdY;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_drhodY.push_back(controls.getId_cellVar("partial-density-mass-fraction-"+controls.spName[i]));
		id_dHtdY.push_back(controls.getId_cellVar("partial-total-enthalpy-mass-fraction-"+controls.spName[i]));
	}
	
	int id_alpha = controls.getId_faceVar("cosine angle of between face normal and cells");
	
	int id_nLRx = controls.getId_faceVar("x unit normal of between left and right cell");
	int id_nLRy = controls.getId_faceVar("y unit normal of between left and right cell");
	int id_nLRz = controls.getId_faceVar("z unit normal of between left and right cell");
	
	// 곡률관련
	vector<int> id_curvature, id_alpha_VF;
	vector<double> surf_sigma;
	for(int i=0; i<controls.spName.size(); ++i){
		string name = controls.spName[i];
		string type = controls.thermophysicalProperties[name+".transport.sigma.type"];
		if(type=="constant"){
			double value = stod(controls.thermophysicalProperties[name+".transport.sigma.value"]);
			if(value<1.e-200) continue;
			id_curvature.push_back(controls.getId_cellVar("curvature-"+name));
			id_alpha_VF.push_back(controls.getId_cellVar("volume-fraction-"+name));
			surf_sigma.push_back(value);
		}
	}
	int nCurv = id_curvature.size();
	
	
	
	
	
	// 크랭크-니콜슨 방법 계수
	// double CN_coeff = 0.5;
	double CN_coeff = 1.0;
	double CN_coeff_Y = 1.0;
	
	


	{
		calcConvFlux.push_back(
		[&solver,nSp,
		id_p,id_pF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
		id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_muF,
		id_rhoF,id_HtF,id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_wd,id_rho,id_Ht,id_YF,nEq,id_c,CN_coeff,CN_coeff_Y,id_Y,
		id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz,
		id_alpha,id_nLRx,id_nLRy,id_nLRz,
		id_curvature,id_alpha_VF,nCurv,surf_sigma](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double dt = fields[id_dt];
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double cL = cellsL[id_c]; double cR = cellsR[id_c];
			double uL = cellsL[id_u]; double uR = cellsR[id_u];
			double vL = cellsL[id_v]; double vR = cellsR[id_v];
			double wL = cellsL[id_w]; double wR = cellsR[id_w];
			double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
			double dudxL = cellsL[id_dudx]; double dudxR = cellsR[id_dudx];
			double dudyL = cellsL[id_dudy]; double dudyR = cellsR[id_dudy];
			double dudzL = cellsL[id_dudz]; double dudzR = cellsR[id_dudz];
			double dvdxL = cellsL[id_dvdx]; double dvdxR = cellsR[id_dvdx];
			double dvdyL = cellsL[id_dvdy]; double dvdyR = cellsR[id_dvdy];
			double dvdzL = cellsL[id_dvdz]; double dvdzR = cellsR[id_dvdz];
			double dwdxL = cellsL[id_dwdx]; double dwdxR = cellsR[id_dwdx];
			double dwdyL = cellsL[id_dwdy]; double dwdyR = cellsR[id_dwdy];
			double dwdzL = cellsL[id_dwdz]; double dwdzR = cellsR[id_dwdz];
			double curvatureL[nCurv+1], curvatureR[nCurv+1];
			double alpha_VFL[nCurv+1], alpha_VFR[nCurv+1];
			for(int i=0; i<nCurv; ++i){
				curvatureL[i] = cellsL[id_curvature[i]];
				curvatureR[i] = cellsR[id_curvature[i]];
				alpha_VFL[i] = cellsL[id_alpha_VF[i]];
				alpha_VFR[i] = cellsR[id_alpha_VF[i]];
			}
			
			double rhoF = faces[id_rhoF];
			double uF = faces[id_uF];
			double vF = faces[id_vF];
			double wF = faces[id_wF];
			double pF = faces[id_pF];
			double UnF = faces[id_UnF];
			double HtF = faces[id_HtF];
			double YF[nSp];
			for(int i=0; i<nSp-1; ++i){
				YF[i] = faces[id_YF[i]];
			}
			
			double weiL = 1.0; double weiR = 0.0;
			if(UnF<0.0){
				weiL = 0.0; weiR = 1.0;
			}
			double wdL = faces[id_wd]; double wdR = 1.0-wdL;
			
			
			
			
			
				
				// wdL = 0.5; wdR = 0.5;
			
			
			
			double weidL = wdL; double weidR = wdR;
			double dtrho = dt*(wdL/rhoL+wdR/rhoR);
			double dp_coeff_Un = dAlpha * dtrho/dLR;
			
			double drhodpL = cellsL[id_drhodp]; double drhodpR = cellsR[id_drhodp];
			double drhodTL = cellsL[id_drhodT]; double drhodTR = cellsR[id_drhodT];
			double dHtdpL = cellsL[id_dHtdp]; double dHtdpR = cellsR[id_dHtdp];
			double dHtdTL = cellsL[id_dHtdT]; double dHtdTR = cellsR[id_dHtdT];
			double drhodYL[nSp]; double drhodYR[nSp];
			double dHtdYL[nSp]; double dHtdYR[nSp];
			for(int i=0; i<nSp-1; ++i){
				drhodYL[i] = cellsL[id_drhodY[i]]; drhodYR[i] = cellsR[id_drhodY[i]];
				dHtdYL[i] = cellsL[id_dHtdY[i]]; dHtdYR[i] = cellsR[id_dHtdY[i]];
			}
			
			double muF = faces[id_muF];
			double ubar = wdL*uL+wdR*uR;
			double vbar = wdL*vL+wdR*vR;
			double wbar = wdL*wL+wdR*wR;
			double dudxF = wdL*dudxL+wdR*dudxR;
			double dudyF = wdL*dudyL+wdR*dudyR;
			double dudzF = wdL*dudzL+wdR*dudzR;
			double dvdxF = wdL*dvdxL+wdR*dvdxR;
			double dvdyF = wdL*dvdyL+wdR*dvdyR;
			double dvdzF = wdL*dvdzL+wdR*dvdzR;
			double dwdxF = wdL*dwdxL+wdR*dwdxR;
			double dwdyF = wdL*dwdyL+wdR*dwdyR;
			double dwdzF = wdL*dwdzL+wdR*dwdzR;
			
			double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
			double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
			
			
			double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
			double w_phi = 1.0-2.0*abs(wdL-0.5);
			
			
			w_phi = 1.0;
			
			
			
			
			

			double chat= wdL*cL+wdR*cR;
			double ML = UnL/chat;
			double MR = UnR/chat;
			double Mcy = min(1.0,KLR/chat);
			double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// double Unhat = 0.5*(UnL+UnR);
			double rhohat = wdL*rhoL+wdR*rhoR;
			double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );
			// double MLP = 0.5*(ML+abs(ML));
			// if( abs(ML) < 1.0 ) {
				// MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			// }
			// double MRM = 0.5*(MR-abs(MR));
			// if( abs(MR) < 1.0 ) {
				// MRM = -0.25*(MR - 1.0)*(MR - 1.0);
			// }
			// SLAU
			double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			double D_L = ML+(1.0-g_c)*abs(ML);
			double D_R = MR-(1.0-g_c)*abs(MR);
			double D_rho = Mbar*g_c;
			double MLPL = 0.5*(D_L+D_rho);
			double MRMR = 0.5*(D_R-D_rho);
			// // 영린 개발 스킴
			// double UnF = w_phi*(MLP+MRM)*chat + (1.0-w_phi)*(wdL*UnL+wdR*UnR);
			
			// // 열역학적 보간
			// UnF -= 0.5*phi_c/rhohat/chat*(pR-pL);
			// // 셀 to 페이스 보간
			// UnF -= dAlpha * dtrho*(pR-pL)/dLR;
			// // non-orthogonal
			// UnF -= (nvec[0]-dAlpha*nLR[0]) * dtrho*(wdL*dpdxL+wdR*dpdxR);
			// UnF -= (nvec[1]-dAlpha*nLR[1]) * dtrho*(wdL*dpdyL+wdR*dpdyR);
			// UnF -= (nvec[2]-dAlpha*nLR[2]) * dtrho*(wdL*dpdzL+wdR*dpdzR);
			
			// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*nvec[0];
			// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*nvec[1];
			// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*nvec[2];
			
			
				
			// double UnF = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);
			// UnF -= dAlpha * dt*(pR-pL)/dLR;
			// UnF += dt*(wdL*dpdxL+wdR*dpdxR)*nvec[0];
			// UnF += dt*(wdL*dpdyL+wdR*dpdyR)*nvec[1];
			// UnF += dt*(wdL*dpdzL+wdR*dpdzR)*nvec[2];
			// if( UnF >= 0.0 ) {
				// UnF = UnF/rhoL;
			// }
			// else{
				// UnF = UnF/rhoR;
			// }
			
			double WUL = 0.5;
			double WUR = 0.5;
			
			
			
			
			double dp_coeff_thm = 0.5*phi_c/rhohat/chat;
			if(UnF>=0.0){
				dp_coeff_thm = 0.5*phi_c/rhoL/chat;
			}
			else{
				dp_coeff_thm = 0.5*phi_c/rhoR/chat;
			}
				

			// // double ML = UnL/chat;
			// // double MR = UnR/chat;
			// // double WUL = 0.5;
			// // double WUR = 0.5;
			// if(abs(ML) > 1.0){ 
				// WUL = 0.5*(1.0 + ((ML > 0.0) ? 1.0 : -1.0) );
			// }
			// else{
				// WUL = 0.5*(ML + 1.0);
				// // WUL += 0.125*2.0*2.0*(ML*ML-1.0);
			// }
			// if(abs(MR) > 1.0){ 
				// WUR = 0.5*(1.0 - ((MR > 0.0) ? 1.0 : -1.0) );
			// }
			// else{
				// WUR = -0.5*(MR - 1.0);
				// // WUR -= 0.125*2.0*2.0*(MR*MR-1.0);
			// }
			// // WUL = w_phi*WUL+(1.0-w_phi)*wdL;
			// // WUR = w_phi*WUR+(1.0-w_phi)*wdR;
			

			double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			if( abs(ML) < 1.0 ) {
				PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			} 
			double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			if( abs(MR) < 1.0 ) {
				PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
			} 
			double WpL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
			double WpR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
			
			// WpL = (w_phi*WpL+(1.0-w_phi)*wdL);
			// WpR = (w_phi*WpR+(1.0-w_phi)*wdR);
			
			// double WpUL = 0.0;
			// double WpUR = 0.0;
			// if( abs(ML) < 1.0 ) {
				// WpUL = 0.75*(1.0-ML*ML);
			// } 
			// if( abs(MR) < 1.0 ) {
				// WpUR = -0.75*(1.0-MR*MR);
			// }
			
			
			if(UnF>=0.0){
				double signUnL = (UnL>0.0 ? 1.0 : -1.0);
				double signUnR = (UnR>0.0 ? 1.0 : -1.0);
				WUL = (rhoL*0.5*(1.0+(1.0-g_c)*signUnL+rhoL*signUnL/(rhoL+rhoR)*g_c))/rhoL;
				WUR = (rhoR*0.5*(1.0-(1.0-g_c)*signUnR-rhoR*signUnR/(rhoL+rhoR)*g_c))/rhoL;
				// WUL = (0.5*(1.0+(1.0-g_c)*signUnL+signUnL/(rhoL+rhoR)*g_c));
				// WUR = (0.5*(1.0-(1.0-g_c)*signUnR-signUnR/(rhoL+rhoR)*g_c));
				
				dp_coeff_Un = dAlpha * dt/dLR/rhoL;
				dp_coeff_thm = 0.5*phi_c/chat/rhoL;
			}
			else{
				double signUnL = (UnL>0.0 ? 1.0 : -1.0);
				double signUnR = (UnR>0.0 ? 1.0 : -1.0);
				WUL = (rhoL*0.5*(1.0+(1.0-g_c)*signUnL+rhoL*signUnL/(rhoL+rhoR)*g_c))/rhoR;
				WUR = (rhoR*0.5*(1.0-(1.0-g_c)*signUnR-rhoR*signUnR/(rhoL+rhoR)*g_c))/rhoR;
				// WUL = (0.5*(1.0+(1.0-g_c)*signUnL+signUnL/(rhoL+rhoR)*g_c));
				// WUR = (0.5*(1.0-(1.0-g_c)*signUnR-signUnR/(rhoL+rhoR)*g_c));
				
				dp_coeff_Un = dAlpha * dt/dLR/rhoR;
				dp_coeff_thm = 0.5*phi_c/chat/rhoR;
			}
			
			
			
			// WUL = 0.5;
			// WUR = 0.5;
			// WpL = 0.5;
			// WpR = 0.5;
			// dp_coeff_Un = dAlpha * dtrho/dLR;
			// dp_coeff_thm = 0.5*phi_c/chat/rhohat;
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			int iter = 0;
			
			double drho_dUn_dpL = weiL*drhodpL*UnF+rhoF*dp_coeff_Un +rhoF*dp_coeff_thm;
			double drho_dUn_dpR = weiR*drhodpR*UnF-rhoF*dp_coeff_Un -rhoF*dp_coeff_thm;
			double rhoF_dUn_duL = rhoF*(WUL*nvec[0]);
			double rhoF_dUn_duR = rhoF*(WUR*nvec[0]);
			double rhoF_dUn_dvL = rhoF*(WUL*nvec[1]);
			double rhoF_dUn_dvR = rhoF*(WUR*nvec[1]);
			double rhoF_dUn_dwL = rhoF*(WUL*nvec[2]);
			double rhoF_dUn_dwR = rhoF*(WUR*nvec[2]);
			
			// fluxA_LL[iter] += CN_coeff * drho_dUn_dpL*area; 
			// fluxA_LR[iter] += CN_coeff * drho_dUn_dpR*area;
			// fluxA_RR[iter] -= CN_coeff * drho_dUn_dpR*area; 
			// fluxA_RL[iter] -= CN_coeff * drho_dUn_dpL*area; 
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF)*area; 
			// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF)*area;
			// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF)*area;
			// ++iter;
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF)*area; 
				// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF)*area;
				// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF)*area; 
				// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF)*area;
				// ++iter;
			// }
			
			
			
			
			// fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*uF +WpL*nvec[0])*area; 
			// fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*uF +WpR*nvec[0])*area;
			// fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*uF +WpR*nvec[0])*area; 
			// fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*uF +WpL*nvec[0])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*uF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*uF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*uF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*uF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*uF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*uF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*uF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*uF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*uF)*area; 
			// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*uF)*area;
			// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*uF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*uF)*area;
			// ++iter;
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*uF)*area; 
				// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*uF)*area;
				// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*uF)*area; 
				// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*uF)*area;
				// ++iter;
			// }
			
			
			
			
			// fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*vF +WpL*nvec[1])*area; 
			// fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*vF +WpR*nvec[1])*area;
			// fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*vF +WpR*nvec[1])*area; 
			// fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*vF +WpL*nvec[1])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*vF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*vF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*vF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*vF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*vF)*area; 
			// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*vF)*area;
			// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*vF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*vF)*area;
			// ++iter;
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*vF)*area; 
				// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*vF)*area;
				// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*vF)*area; 
				// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*vF)*area;
				// ++iter;
			// }
			
			
			
			
			// fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*wF +WpL*nvec[2])*area; 
			// fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*wF +WpR*nvec[2])*area;
			// fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*wF +WpR*nvec[2])*area; 
			// fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*wF +WpL*nvec[2])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*wF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*wF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*wF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*wF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*wF)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*wF)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*wF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*wF)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*wF)*area; 
			// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*wF)*area;
			// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*wF)*area; 
			// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*wF)*area;
			// ++iter;
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*wF)*area; 
				// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*wF)*area;
				// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*wF)*area; 
				// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*wF)*area;
				// ++iter;
			// }
			
			
			
			
			// fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*HtF +rhoF*UnF*weiL*dHtdpL)*area; 
			// fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*HtF +rhoF*UnF*weiR*dHtdpR)*area;
			// fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*HtF +rhoF*UnF*weiR*dHtdpR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*HtF +rhoF*UnF*weiL*dHtdpL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*HtF+rhoF*UnF*weiL*uL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*HtF+rhoF*UnF*weiR*uR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*HtF+rhoF*UnF*weiR*uR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*HtF+rhoF*UnF*weiL*uL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*HtF+rhoF*UnF*weiL*vL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*HtF+rhoF*UnF*weiR*vR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*HtF+rhoF*UnF*weiR*vR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*HtF+rhoF*UnF*weiL*vL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*HtF+rhoF*UnF*weiL*wL)*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*HtF+rhoF*UnF*weiR*wR)*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*HtF+rhoF*UnF*weiR*wR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*HtF+rhoF*UnF*weiL*wL)*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*HtF+rhoF*UnF*weiL*dHtdTL)*area; 
			// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*HtF+rhoF*UnF*weiR*dHtdTR)*area;
			// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*HtF+rhoF*UnF*weiR*dHtdTR)*area; 
			// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*HtF+rhoF*UnF*weiL*dHtdTL)*area;
			// ++iter;
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*HtF+rhoF*UnF*weiL*dHtdYL[i])*area; 
				// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*HtF+rhoF*UnF*weiR*dHtdYR[i])*area;
				// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*HtF+rhoF*UnF*weiR*dHtdYR[i])*area; 
				// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*HtF+rhoF*UnF*weiL*dHtdYL[i])*area;
				// ++iter;
			// }
			
			
			
			// for(int i=0; i<nSp-1; ++i){
				// fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*YF[i])*area; 
				// fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*YF[i])*area;
				// fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*YF[i])*area; 
				// fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*YF[i])*area;
				// ++iter;
				
				// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*YF[i])*area; 
				// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*YF[i])*area;
				// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*YF[i])*area; 
				// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*YF[i])*area;
				// ++iter;
				
				// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*YF[i])*area; 
				// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*YF[i])*area;
				// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*YF[i])*area; 
				// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*YF[i])*area;
				// ++iter;
				
				// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*YF[i])*area; 
				// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*YF[i])*area;
				// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*YF[i])*area; 
				// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*YF[i])*area;
				// ++iter;
				
				// fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*YF[i])*area; 
				// fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*YF[i])*area;
				// fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*YF[i])*area; 
				// fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*YF[i])*area;
				// ++iter;
				
				// for(int j=0; j<nSp-1; ++j){
					// fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[j]*UnF*YF[i])*area; 
					// fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[j]*UnF*YF[i])*area;
					// fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[j]*UnF*YF[i])*area; 
					// fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[j]*UnF*YF[i])*area;
					// if(i==j){
						// fluxA_LL[iter] += CN_coeff_Y * (rhoF*UnF*weiL)*area; 
						// fluxA_LR[iter] += CN_coeff_Y * (rhoF*UnF*weiR)*area;
						// fluxA_RR[iter] -= CN_coeff_Y * (rhoF*UnF*weiR)*area; 
						// fluxA_RL[iter] -= CN_coeff_Y * (rhoF*UnF*weiL)*area;
					// }
					// ++iter;
				// }
			// }
			
			
			
			double fluxB[nEq];
			
			// 컨벡티브 B
			fluxB[0] = -( rhoF*UnF )*area;
			fluxB[1] = -( rhoF*UnF*uF + pF*nvec[0] )*area;
			fluxB[2] = -( rhoF*UnF*vF + pF*nvec[1] )*area;
			fluxB[3] = -( rhoF*UnF*wF + pF*nvec[2] )*area;
			fluxB[4] = -( rhoF*UnF*HtF )*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[5+i] = -( rhoF*UnF*YF[i] )*area;
			}
			
			
			
			double visc_coeff = dAlpha*muF/dLR*area;
			
			// iter = nEq*1+1;
			// fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			// fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			// iter = nEq*2+2;
			// fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			// fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			// iter = nEq*3+3;
			// fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			// fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			// iter = nEq*4+1;
			// fluxA_LL[iter] += visc_coeff*ubar; fluxA_LR[iter] -= visc_coeff*ubar;
			// fluxA_RR[iter] += visc_coeff*ubar; fluxA_RL[iter] -= visc_coeff*ubar;
			
			// iter = nEq*4+2;
			// fluxA_LL[iter] += visc_coeff*vbar; fluxA_LR[iter] -= visc_coeff*vbar;
			// fluxA_RR[iter] += visc_coeff*vbar; fluxA_RL[iter] -= visc_coeff*vbar;
			
			// iter = nEq*4+3;
			// fluxA_LL[iter] += visc_coeff*wbar; fluxA_LR[iter] -= visc_coeff*wbar;
			// fluxA_RR[iter] += visc_coeff*wbar; fluxA_RL[iter] -= visc_coeff*wbar;
			
			
			
			
			// // 디퓨젼 B
			// fluxB[1] += muF*(
						// (dAlpha*(uR-uL)/dLR + dudxF*nvec[0] + dvdxF*nvec[1] + dwdxF*nvec[2]) - 
						// 2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[0] 
						// )*area;
			// // fluxB[1] -= nvec[0] * 2.0/3.0 * rhoF * tkei;
			// fluxB[2] += muF*(
						// (dAlpha*(vR-vL)/dLR + dudyF*nvec[0] + dvdyF*nvec[1] + dwdyF*nvec[2]) - 
						// 2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[1]
						// )*area;
			// // fluxB[2] -= nvec[1] * 2.0/3.0 * rhoF * tkei;
			// fluxB[3] += muF*(
						// (dAlpha*(wR-wL)/dLR + dudzF*nvec[0] + dvdzF*nvec[1] + dwdzF*nvec[2]) - 
						// 2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[2]
						// )*area;
			// // // fluxB[3] -= nvec[2] * 2.0/3.0 * rhoF * tkei;
			// fluxB[4] += muF*(
						// dAlpha*(uR-uL)/dLR*ubar + dAlpha*(vR-vL)/dLR*vbar + dAlpha*(wR-wL)/dLR*wbar + 
						// (dudxF*ubar + dudyF*vbar + dudzF*wbar)*nvec[0] +
						// (dvdxF*ubar + dvdyF*vbar + dvdzF*wbar)*nvec[1] +
						// (dwdxF*ubar + dwdyF*vbar + dwdzF*wbar)*nvec[2] -
						// 2.0/3.0*(dudxF + dvdyF + dwdzF)*(ubar*nvec[0]+vbar*nvec[1]+wbar*nvec[2])
						// )*area;
			
			// // non-orthogonal
			// fluxB[1] += (nvec[0]-dAlpha*nLR[0]) * muF*dudxF*area;
			// fluxB[1] += (nvec[1]-dAlpha*nLR[1]) * muF*dudyF*area;
			// fluxB[1] += (nvec[2]-dAlpha*nLR[2]) * muF*dudzF*area;
			
			// fluxB[2] += (nvec[0]-dAlpha*nLR[0]) * muF*dvdxF*area;
			// fluxB[2] += (nvec[1]-dAlpha*nLR[1]) * muF*dvdyF*area;
			// fluxB[2] += (nvec[2]-dAlpha*nLR[2]) * muF*dvdzF*area;
			
			// fluxB[3] += (nvec[0]-dAlpha*nLR[0]) * muF*dwdxF*area;
			// fluxB[3] += (nvec[1]-dAlpha*nLR[1]) * muF*dwdyF*area;
			// fluxB[3] += (nvec[2]-dAlpha*nLR[2]) * muF*dwdzF*area;
			
			// fluxB[4] += (nvec[0]-dAlpha*nLR[0]) * muF*(
						// dudxF*ubar+dvdxF*vbar+dwdxF*wbar)*area;
			// fluxB[4] += (nvec[1]-dAlpha*nLR[1]) * muF*(
						// dudyF*ubar+dvdyF*vbar+dwdyF*wbar)*area;
			// fluxB[4] += (nvec[2]-dAlpha*nLR[2]) * muF*(
						// dudzF*ubar+dvdzF*vbar+dwdzF*wbar)*area;
			
			

			// fluxB_LL[0] += fluxB[0]; fluxB_RR[0] -= fluxB[0];
			// fluxB_LL[1] += fluxB[1]; fluxB_RR[1] -= fluxB[1];
			// fluxB_LL[2] += fluxB[2]; fluxB_RR[2] -= fluxB[2];
			// fluxB_LL[3] += fluxB[3]; fluxB_RR[3] -= fluxB[3];
			// fluxB_LL[4] += fluxB[4]; fluxB_RR[4] -= fluxB[4];
			// for(int i=0; i<nSp-1; ++i){
				// fluxB_LL[5+i] += fluxB[5+i]; fluxB_RR[5+i] -= fluxB[5+i];
			// }
		
			
			// // 표면장력
			// for(int i=0; i<nCurv; ++i){
				// double alpha_VFF = wdL*alpha_VFL[i]+wdR*alpha_VFR[i];
				// fluxB_LL[1] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
				// fluxB_LL[2] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
				// fluxB_LL[3] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				
				// fluxB_RR[1] += surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[0] )*area;
				// fluxB_RR[2] += surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[1] )*area;
				// fluxB_RR[3] += surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[2] )*area;
				
				// // fluxB_LL[1] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
				// // fluxB_LL[2] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
				// // fluxB_LL[3] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				
				// // fluxB_RR[1] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[0] )*area;
				// // fluxB_RR[2] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[1] )*area;
				// // fluxB_RR[3] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[2] )*area;
			// }
			
			
			
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
			vector<bool> bool_inletOutlet(nEq,false);
			
			if(controls.boundaryMap["pressure"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			if(controls.boundaryMap["pressure"][bcName+".type"]=="inletOutlet") bool_inletOutlet[0] = true;
			
			if(controls.boundaryMap["velocity"][bcName+".type"]=="zeroGradient") bool_zeroGradient[1] = true;
			
			if(controls.boundaryMap["temperature"][bcName+".type"]=="zeroGradient") bool_zeroGradient[2] = true;
			if(controls.boundaryMap["temperature"][bcName+".type"]=="inletOutlet") bool_inletOutlet[2] = true;
			
			for(int i=0; i<controls.spName.size()-1; ++i){
				string tmp_name = ("mass-fraction-"+controls.spName[i]);
				if(controls.boundaryMap[tmp_name][bcName+".type"]=="zeroGradient") bool_zeroGradient[3+i] = true;
				if(controls.boundaryMap[tmp_name][bcName+".type"]=="inletOutlet") bool_inletOutlet[3+i] = true;
			}
			
			
			calcConvFlux_BC[bc_id].back().push_back(
			[&solver,nSp,
			id_p,id_pF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
			id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_muF,
			id_rhoF,id_HtF,id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
			id_wd,id_rho,id_YF,nEq,id_c,CN_coeff,CN_coeff_Y,
			id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz,
			bool_zeroGradient,bool_inletOutlet,
			id_alpha,id_nLRx,id_nLRy,id_nLRz,
			id_curvature,id_alpha_VF,nCurv,surf_sigma](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				double dt = fields[id_dt];
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				
				double cL = cellsL[id_c];
				double uL = cellsL[id_u];
				double vL = cellsL[id_v];
				double wL = cellsL[id_w];
				double rhoL = cellsL[id_rho];
				double dudxL = cellsL[id_dudx];
				double dudyL = cellsL[id_dudy];
				double dudzL = cellsL[id_dudz];
				double dvdxL = cellsL[id_dvdx];
				double dvdyL = cellsL[id_dvdy];
				double dvdzL = cellsL[id_dvdz];
				double dwdxL = cellsL[id_dwdx];
				double dwdyL = cellsL[id_dwdy];
				double dwdzL = cellsL[id_dwdz];
				double curvatureL[nCurv+1];
				double alpha_VFL[nCurv+1];
				for(int i=0; i<nCurv; ++i){
					curvatureL[i] = cellsL[id_curvature[i]];
					alpha_VFL[i] = cellsL[id_alpha_VF[i]];
				}
				
				// double cF = faces[id_cF];
				double rhoF = faces[id_rhoF];
				double uF = faces[id_uF];
				double vF = faces[id_vF];
				double wF = faces[id_wF];
				double pF = faces[id_pF];
				double HtF = faces[id_HtF];
				double YF[nSp];
				for(int i=0; i<nSp-1; ++i){
					YF[i] = faces[id_YF[i]];
				}
				
				double uR = uF;
				double vR = vF;
				double wR = wF;
				
				double UnF = uF*nvec[0]+vF*nvec[1]+wF*nvec[2];
				// double weiL = 1.0; double weiR = 0.0;
				// if(UnF<0.0){
					// weiL = 0.0; weiR = 1.0;
				// }
				// double weidL = weiL; double weidR = weiR;
				double dtrho = dt*(1.0/rhoL);
				double dp_coeff_Un = dAlpha * dtrho/dLR;
				
				double drhodpL = cellsL[id_drhodp];
				double drhodTL = cellsL[id_drhodT];
				double dHtdpL = cellsL[id_dHtdp];
				double dHtdTL = cellsL[id_dHtdT];
				double drhodYL[nSp];
				double dHtdYL[nSp];
				for(int i=0; i<nSp-1; ++i){
					drhodYL[i] = cellsL[id_drhodY[i]];
					dHtdYL[i] = cellsL[id_dHtdY[i]];
				}
				
				double muF = faces[id_muF];
				double ubar = uF;
				double vbar = vF;
				double wbar = wF;
				double dudxF = dudxL;
				double dudyF = dudyL;
				double dudzF = dudzL;
				double dvdxF = dvdxL;
				double dvdyF = dvdyL;
				double dvdzF = dvdzL;
				double dwdxF = dwdxL;
				double dwdyF = dwdyL;
				double dwdzF = dwdzL;
				
				double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uF*uF+vF*vF+wF*wF));
				double chat = cL;//0.5*(cL+cF);
				double Mdash = min(1.0,KLR/chat);
				double chi = (1.0-Mdash)*(1.0-Mdash);
				double rhohat = 0.5*rhoL+0.5*rhoF;
				double Mcy = min(1.0,KLR/chat);
				double phi_c = (1.0-Mcy)*(1.0-Mcy);
				double dp_coeff_thm = 0.5*chi/rhohat/chat;
				
				
				
			
			
				dp_coeff_Un = dAlpha * dt/dLR/rhoL;
				dp_coeff_thm = 0.5*phi_c/chat/rhoL;
				
				
				
				// if(
				// bool_zeroGradient[0]==true ||
				// (bool_inletOutlet[0]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				// ){
					// fluxA_LL[nEq*0+0] += CN_coeff * (drhodpL*UnF)*area; 
					// fluxA_LL[nEq*1+0] += CN_coeff * (drhodpL*UnF*uF + nvec[0])*area; 
					// fluxA_LL[nEq*2+0] += CN_coeff * (drhodpL*UnF*vF + nvec[1])*area; 
					// fluxA_LL[nEq*3+0] += CN_coeff * (drhodpL*UnF*wF + nvec[2])*area; 
					// fluxA_LL[nEq*4+0] += CN_coeff * (drhodpL*UnF*HtF+rhoF*UnF*dHtdpL)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+0] += CN_coeff_Y * (drhodpL*UnF*YF[i])*area; 
					// }
					
				// }
				// else{
					// // fluxA_LL[nEq*0+0] += CN_coeff * 0.5*(drhodpL*UnF)*area; 
					// // fluxA_LL[nEq*1+0] += CN_coeff * 0.5*(drhodpL*UnF*uF + nvec[0])*area; 
					// // fluxA_LL[nEq*2+0] += CN_coeff * 0.5*(drhodpL*UnF*vF + nvec[1])*area; 
					// // fluxA_LL[nEq*3+0] += CN_coeff * 0.5*(drhodpL*UnF*wF + nvec[2])*area; 
					// // fluxA_LL[nEq*4+0] += CN_coeff * 0.5*(drhodpL*UnF*HtF+rhoF*UnF*dHtdpL)*area; 
					// // for(int i=0; i<nSp-1; ++i){
						// // fluxA_LL[nEq*(5+i)+0] += CN_coeff_Y * 0.5*(drhodpL*UnF*YF[i])*area; 
					// // }
					// fluxA_LL[nEq*0+0] += CN_coeff * (rhoF*dp_coeff_Un +rhoF*dp_coeff_thm)*area;
					// fluxA_LL[nEq*1+0] += CN_coeff * (rhoF*dp_coeff_Un*uF +rhoF*dp_coeff_thm*uF)*area; 
					// fluxA_LL[nEq*2+0] += CN_coeff * (rhoF*dp_coeff_Un*vF +rhoF*dp_coeff_thm*vF)*area; 
					// fluxA_LL[nEq*3+0] += CN_coeff * (rhoF*dp_coeff_Un*wF +rhoF*dp_coeff_thm*wF)*area; 
					// fluxA_LL[nEq*4+0] += CN_coeff * (rhoF*dp_coeff_Un*HtF +rhoF*dp_coeff_thm*HtF)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+0] += CN_coeff_Y * (rhoF*dp_coeff_Un*YF[i] +rhoF*dp_coeff_thm*YF[i])*area; 
					// }
				// }
				
				// if(
				// bool_zeroGradient[1]==true ||
				// (bool_inletOutlet[1]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				// ){
					// fluxA_LL[nEq*0+1] += CN_coeff * (rhoF*(nvec[0]))*area; 
					// fluxA_LL[nEq*1+1] += CN_coeff * (rhoF*(nvec[0])*uF+rhoF*UnF)*area; 
					// fluxA_LL[nEq*2+1] += CN_coeff * (rhoF*(nvec[0])*vF)*area; 
					// fluxA_LL[nEq*3+1] += CN_coeff * (rhoF*(nvec[0])*wF)*area; 
					// fluxA_LL[nEq*4+1] += CN_coeff * (rhoF*(nvec[0])*HtF+rhoF*UnF*uL)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+1] += CN_coeff_Y * (rhoF*(nvec[0])*YF[i])*area; 
					// }
					
					// fluxA_LL[nEq*0+2] += CN_coeff * (rhoF*(nvec[1]))*area; 
					// fluxA_LL[nEq*1+2] += CN_coeff * (rhoF*(nvec[1])*uF)*area; 
					// fluxA_LL[nEq*2+2] += CN_coeff * (rhoF*(nvec[1])*vF+rhoF*UnF)*area; 
					// fluxA_LL[nEq*3+2] += CN_coeff * (rhoF*(nvec[1])*wF)*area; 
					// fluxA_LL[nEq*4+2] += CN_coeff * (rhoF*(nvec[1])*HtF+rhoF*UnF*vL)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+2] += CN_coeff_Y * (rhoF*(nvec[1])*YF[i])*area; 
					// }
					
					// fluxA_LL[nEq*0+3] += CN_coeff * (rhoF*(nvec[2]))*area; 
					// fluxA_LL[nEq*1+3] += CN_coeff * (rhoF*(nvec[2])*uF)*area; 
					// fluxA_LL[nEq*2+3] += CN_coeff * (rhoF*(nvec[2])*vF)*area; 
					// fluxA_LL[nEq*3+3] += CN_coeff * (rhoF*(nvec[2])*wF+rhoF*UnF)*area; 
					// fluxA_LL[nEq*4+3] += CN_coeff * (rhoF*(nvec[2])*HtF+rhoF*UnF*wL)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+3] += CN_coeff_Y * (rhoF*(nvec[2])*YF[i])*area; 
					// }
				// }
				// else{
					
					// // fluxA_LL[nEq*0+1] += CN_coeff * 0.5*(rhoF*(nvec[0]))*area; 
					// // fluxA_LL[nEq*1+1] += CN_coeff * 0.5*(rhoF*(nvec[0])*uF+rhoF*UnF)*area; 
					// // fluxA_LL[nEq*2+1] += CN_coeff * 0.5*(rhoF*(nvec[0])*vF)*area; 
					// // fluxA_LL[nEq*3+1] += CN_coeff * 0.5*(rhoF*(nvec[0])*wF)*area; 
					// // fluxA_LL[nEq*4+1] += CN_coeff * 0.5*(rhoF*(nvec[0])*HtF+rhoF*UnF*uL)*area; 
					// // for(int i=0; i<nSp-1; ++i){
						// // fluxA_LL[nEq*(5+i)+1] += CN_coeff_Y * 0.5*(rhoF*(nvec[0])*YF[i])*area; 
					// // }
					
					// // fluxA_LL[nEq*0+2] += CN_coeff * 0.5*(rhoF*(nvec[1]))*area; 
					// // fluxA_LL[nEq*1+2] += CN_coeff * 0.5*(rhoF*(nvec[1])*uF)*area; 
					// // fluxA_LL[nEq*2+2] += CN_coeff * 0.5*(rhoF*(nvec[1])*vF+rhoF*UnF)*area; 
					// // fluxA_LL[nEq*3+2] += CN_coeff * 0.5*(rhoF*(nvec[1])*wF)*area; 
					// // fluxA_LL[nEq*4+2] += CN_coeff * 0.5*(rhoF*(nvec[1])*HtF+rhoF*UnF*vL)*area; 
					// // for(int i=0; i<nSp-1; ++i){
						// // fluxA_LL[nEq*(5+i)+2] += CN_coeff_Y * 0.5*(rhoF*(nvec[1])*YF[i])*area; 
					// // }
					
					// // fluxA_LL[nEq*0+3] += CN_coeff * 0.5*(rhoF*(nvec[2]))*area; 
					// // fluxA_LL[nEq*1+3] += CN_coeff * 0.5*(rhoF*(nvec[2])*uF)*area; 
					// // fluxA_LL[nEq*2+3] += CN_coeff * 0.5*(rhoF*(nvec[2])*vF)*area; 
					// // fluxA_LL[nEq*3+3] += CN_coeff * 0.5*(rhoF*(nvec[2])*wF+rhoF*UnF)*area; 
					// // fluxA_LL[nEq*4+3] += CN_coeff * 0.5*(rhoF*(nvec[2])*HtF+rhoF*UnF*wL)*area; 
					// // for(int i=0; i<nSp-1; ++i){
						// // fluxA_LL[nEq*(5+i)+3] += CN_coeff_Y * 0.5*(rhoF*(nvec[2])*YF[i])*area; 
					// // }
					
					// double visc_coeff = dAlpha*muF/dLR*area;
					// fluxA_LL[nEq*1+1] += visc_coeff;
					// fluxA_LL[nEq*2+2] += visc_coeff;
					// fluxA_LL[nEq*3+3] += visc_coeff;
					// fluxA_LL[nEq*4+1] += visc_coeff*ubar;
					// fluxA_LL[nEq*4+2] += visc_coeff*vbar;
					// fluxA_LL[nEq*4+3] += visc_coeff*wbar;
					
				// }
				
				// if(
				// bool_zeroGradient[2]==true ||
				// (bool_inletOutlet[2]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				// ){
					// fluxA_LL[nEq*0+4] += CN_coeff * (drhodTL*UnF)*area; 
					// fluxA_LL[nEq*1+4] += CN_coeff * (drhodTL*UnF*uF)*area; 
					// fluxA_LL[nEq*2+4] += CN_coeff * (drhodTL*UnF*vF)*area; 
					// fluxA_LL[nEq*3+4] += CN_coeff * (drhodTL*UnF*wF)*area; 
					// fluxA_LL[nEq*4+4] += CN_coeff * (drhodTL*UnF*HtF+rhoF*UnF*dHtdTL)*area; 
					// for(int i=0; i<nSp-1; ++i){
						// fluxA_LL[nEq*(5+i)+4] += CN_coeff_Y * (drhodTL*UnF*YF[i])*area; 
					// }
				// }
				// else{
					// // fluxA_LL[nEq*0+4] += CN_coeff * 0.5*(drhodTL*UnF)*area; 
					// // fluxA_LL[nEq*1+4] += CN_coeff * 0.5*(drhodTL*UnF*uF)*area; 
					// // fluxA_LL[nEq*2+4] += CN_coeff * 0.5*(drhodTL*UnF*vF)*area; 
					// // fluxA_LL[nEq*3+4] += CN_coeff * 0.5*(drhodTL*UnF*wF)*area; 
					// // fluxA_LL[nEq*4+4] += CN_coeff * 0.5*(drhodTL*UnF*HtF+rhoF*UnF*dHtdTL)*area; 
					// // for(int i=0; i<nSp-1; ++i){
						// // fluxA_LL[nEq*(5+i)+4] += CN_coeff_Y * 0.5*(drhodTL*UnF*YF[i])*area; 
					// // }
					
				// }
				
				// for(int i=0; i<nSp-1; ++i){
					// if(
					// bool_zeroGradient[3+i]==true ||
					// (bool_inletOutlet[3+i]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
					// ){
						// fluxA_LL[nEq*0+5+i] += CN_coeff * (drhodYL[i]*UnF)*area; 
						// fluxA_LL[nEq*1+5+i] += CN_coeff * (drhodYL[i]*UnF*uF)*area; 
						// fluxA_LL[nEq*2+5+i] += CN_coeff * (drhodYL[i]*UnF*vF)*area; 
						// fluxA_LL[nEq*3+5+i] += CN_coeff * (drhodYL[i]*UnF*wF)*area; 
						// fluxA_LL[nEq*4+5+i] += CN_coeff * (drhodYL[i]*UnF*HtF+rhoF*UnF*dHtdYL[i])*area; 
						// for(int j=0; j<nSp-1; ++j){
							// fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * (drhodYL[i]*UnF*YF[j])*area; 
							// if(i==j){
								// fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * (rhoF*UnF)*area; 
							// }
						// }
					
					// }
					// else{
						// // fluxA_LL[nEq*0+5+i] += CN_coeff * 0.5*(drhodYL[i]*UnF)*area; 
						// // fluxA_LL[nEq*1+5+i] += CN_coeff * 0.5*(drhodYL[i]*UnF*uF)*area; 
						// // fluxA_LL[nEq*2+5+i] += CN_coeff * 0.5*(drhodYL[i]*UnF*vF)*area; 
						// // fluxA_LL[nEq*3+5+i] += CN_coeff * 0.5*(drhodYL[i]*UnF*wF)*area; 
						// // fluxA_LL[nEq*4+5+i] += CN_coeff * 0.5*(drhodYL[i]*UnF*HtF+rhoF*UnF*dHtdYL[i])*area; 
						// // for(int j=0; j<nSp-1; ++j){
							// // fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * 0.5*(drhodYL[i]*UnF*YF[j])*area; 
							// // if(i==j){
								// // fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * 0.5*(rhoF*UnF)*area; 
							// // }
						// // }
					// }
				// }
				
				
				fluxB_LL[0] -= ( rhoF*UnF )*area;
				fluxB_LL[1] -= ( rhoF*UnF*uF + pF*nvec[0] )*area;
				fluxB_LL[2] -= ( rhoF*UnF*vF + pF*nvec[1] )*area;
				fluxB_LL[3] -= ( rhoF*UnF*wF + pF*nvec[2] )*area;
				fluxB_LL[4] -= ( rhoF*UnF*HtF )*area;
				for(int i=0; i<nSp-1; ++i){
					fluxB_LL[5+i] -= ( rhoF*UnF*YF[i] )*area;
				}
				
				
				fluxB_LL[1] += muF*(
							(dAlpha*(uR-uL)/dLR + dudxF*nvec[0] + dvdxF*nvec[1] + dwdxF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[0] 
							)*area;
				fluxB_LL[2] += muF*(
							(dAlpha*(vR-vL)/dLR + dudyF*nvec[0] + dvdyF*nvec[1] + dwdyF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[1]
							)*area;
				fluxB_LL[3] += muF*(
							(dAlpha*(wR-wL)/dLR + dudzF*nvec[0] + dvdzF*nvec[1] + dwdzF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[2]
							)*area;
				fluxB_LL[4] += muF*(
							dAlpha*(uR-uL)/dLR*ubar + dAlpha*(vR-vL)/dLR*vbar + dAlpha*(wR-wL)/dLR*wbar + 
							(dudxF*ubar + dudyF*vbar + dudzF*wbar)*nvec[0] +
							(dvdxF*ubar + dvdyF*vbar + dvdzF*wbar)*nvec[1] +
							(dwdxF*ubar + dwdyF*vbar + dwdzF*wbar)*nvec[2] -
							2.0/3.0*(dudxF + dvdyF + dwdzF)*(ubar*nvec[0]+vbar*nvec[1]+wbar*nvec[2])
							)*area;
			
				// non-orthogonal
				fluxB_LL[1] += (nvec[0]-dAlpha*nLR[0]) * muF*dudxF*area;
				fluxB_LL[1] += (nvec[1]-dAlpha*nLR[1]) * muF*dudyF*area;
				fluxB_LL[1] += (nvec[2]-dAlpha*nLR[2]) * muF*dudzF*area;
				
				fluxB_LL[2] += (nvec[0]-dAlpha*nLR[0]) * muF*dvdxF*area;
				fluxB_LL[2] += (nvec[1]-dAlpha*nLR[1]) * muF*dvdyF*area;
				fluxB_LL[2] += (nvec[2]-dAlpha*nLR[2]) * muF*dvdzF*area;
				
				fluxB_LL[3] += (nvec[0]-dAlpha*nLR[0]) * muF*dwdxF*area;
				fluxB_LL[3] += (nvec[1]-dAlpha*nLR[1]) * muF*dwdyF*area;
				fluxB_LL[3] += (nvec[2]-dAlpha*nLR[2]) * muF*dwdzF*area;
				
				fluxB_LL[4] += (nvec[0]-dAlpha*nLR[0]) * muF*(
								dudxF*ubar+dvdxF*vbar+dwdxF*wbar)*area;
				fluxB_LL[4] += (nvec[1]-dAlpha*nLR[1]) * muF*(
								dudyF*ubar+dvdyF*vbar+dwdyF*wbar)*area;
				fluxB_LL[4] += (nvec[2]-dAlpha*nLR[2]) * muF*(
								dudzF*ubar+dvdzF*vbar+dwdzF*wbar)*area;
				
				
				

				
				// 표면장력
				for(int i=0; i<nCurv; ++i){
					double alpha_VFF = alpha_VFL[i];
					fluxB_LL[1] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
					fluxB_LL[2] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
					fluxB_LL[3] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
					
					
					// fluxB_LL[1] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
					// fluxB_LL[2] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
					// fluxB_LL[3] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				}
				
				
				return 0;
			});
		}
		
	}
		
	
}




