
#include "../../../../others/solvers.h"
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
	int id_pMax = controls.getId_cellVar("maximum pressure");
	int id_pMin = controls.getId_cellVar("minimum pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uL = controls.getId_faceVar("left x-velocity");
	int id_uR = controls.getId_faceVar("right x-velocity");
	int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	int id_uMax = controls.getId_cellVar("maximum x-velocity");
	int id_uMin = controls.getId_cellVar("minimum x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity");
	int id_vR = controls.getId_faceVar("right y-velocity");
	int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	int id_vMax = controls.getId_cellVar("maximum y-velocity");
	int id_vMin = controls.getId_cellVar("minimum y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity");
	int id_wR = controls.getId_faceVar("right z-velocity");
	int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	int id_wMax = controls.getId_cellVar("maximum z-velocity");
	int id_wMin = controls.getId_cellVar("minimum z-velocity");
	
	int id_T = controls.getId_cellVar("temperature");
	int id_TL = controls.getId_faceVar("left temperature");
	int id_TR = controls.getId_faceVar("right temperature");
	int id_dTdx = controls.getId_cellVar("x-gradient temperature");
	int id_dTdy = controls.getId_cellVar("y-gradient temperature");
	int id_dTdz = controls.getId_cellVar("z-gradient temperature");
	int id_TMax = controls.getId_cellVar("maximum temperature");
	int id_TMin = controls.getId_cellVar("minimum temperature");
	
	vector<int> id_Y, id_YL, id_YR, id_dYdx, id_dYdy, id_dYdz, id_YMax, id_YMin;
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		id_Y.push_back(controls.getId_cellVar(tmp_name));
		id_dYdx.push_back(controls.getId_cellVar("x-gradient "+tmp_name));
		id_dYdy.push_back(controls.getId_cellVar("y-gradient "+tmp_name));
		id_dYdz.push_back(controls.getId_cellVar("z-gradient "+tmp_name));
		id_YMax.push_back(controls.getId_cellVar("maximum "+tmp_name));
		id_YMin.push_back(controls.getId_cellVar("minimum "+tmp_name));
	}
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("left mass-fraction-"+controls.spName[i]);
		id_YL.push_back(controls.getId_faceVar(tmp_name));
	}
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("right mass-fraction-"+controls.spName[i]);
		id_YR.push_back(controls.getId_faceVar(tmp_name));
	}
	
	int id_rhoL = controls.getId_faceVar("left density");
	int id_rhoR = controls.getId_faceVar("right density");
	int id_HtL = controls.getId_faceVar("left total-enthalpy");
	int id_HtR = controls.getId_faceVar("right total-enthalpy");
	int id_cL = controls.getId_faceVar("left speed-of-sound");
	int id_cR = controls.getId_faceVar("right speed-of-sound");
	int id_mu = controls.getId_cellVar("viscosity");
	int id_muL = controls.getId_faceVar("left viscosity");
	int id_muR = controls.getId_faceVar("right viscosity");
	
	// int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	// int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	// int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	
	// int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	// int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	// int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	
	// int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	// int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	// int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	
	// int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	// int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	// int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	// int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
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
	
	
	// 리미터 관련
	int id_lim_p = controls.getId_cellVar("limiter-unstructured pressure");
	int id_lim_u = controls.getId_cellVar("limiter-unstructured x-velocity");
	int id_lim_v = controls.getId_cellVar("limiter-unstructured y-velocity");
	int id_lim_w = controls.getId_cellVar("limiter-unstructured z-velocity");
	
	
	// 크랭크-니콜슨 방법 계수
	// double CN_coeff = 0.5;
	double CN_coeff = 1.0;
	double CN_coeff_Y = 1.0;
	
	


	{
		calcConvFlux.push_back(
		[&solver,nSp,
		id_p,id_u,id_v,id_w,
		id_pMax,id_pMin,id_uMax,id_uMin,id_vMax,id_vMin,id_wMax,id_wMin,
		id_TMax,id_TMin,id_YMax,id_YMin,
		id_xLNv,id_yLNv,id_zLNv,id_xRNv,id_yRNv,id_zRNv,
		id_pL,id_uL,id_vL,id_wL,id_YL,id_rhoL,id_HtL,id_muL,id_cL,
		id_pR,id_uR,id_vR,id_wR,id_YR,id_rhoR,id_HtR,id_muR,id_cR,
		id_dpdx,id_dpdy,id_dpdz,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_wd,id_rho,id_Ht,nEq,id_c,CN_coeff,CN_coeff_Y,id_Y,
		id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz,
		id_alpha,id_nLRx,id_nLRy,id_nLRz,id_mu,
		id_curvature,id_alpha_VF,nCurv,surf_sigma,
		id_lim_p,id_lim_u,id_lim_v,id_lim_w,
		id_xLF,id_yLF,id_zLF,id_xRF,id_yRF,id_zRF](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double dt = fields[id_dt];
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double LNv[3];
			LNv[0] = faces[id_xLNv]; LNv[1] = faces[id_yLNv]; LNv[2] = faces[id_zLNv];
			double RNv[3];
			RNv[0] = faces[id_xRNv]; RNv[1] = faces[id_yRNv]; RNv[2] = faces[id_zRNv];
			// double xyzLR[3];
			// xyzLR[0] = faces[id_xLR]; xyzLR[1] = faces[id_yLR]; xyzLR[2] = faces[id_zLR];
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double xyzLF[3],xyzRF[3];
			xyzLF[0] = faces[id_xLF]; xyzLF[1] = faces[id_yLF]; xyzLF[2] = faces[id_zLF];
			xyzRF[0] = faces[id_xRF]; xyzRF[1] = faces[id_yRF]; xyzRF[2] = faces[id_zRF];
			
			double wdL = faces[id_wd]; double wdR = 1.0-wdL;
			
			
				// dAlpha = 1.0;
				// wdL = 1.0-wdL; wdR = 1.0-wdR; 
				wdL = 0.5; wdR = 0.5; 
				// double weiwd = 0.9;
				// wdL = weiwd*0.5 + (1.0-weiwd)*(wdL-0.5);
				// wdR = 1.0 - wdL;
			
			
			double w_phi = 1.0-2.0*abs(wdL-0.5);
			double pL = cellsL[id_p]; double pR = cellsR[id_p];
			double cL = cellsL[id_c]; double cR = cellsR[id_c];
			double uL = cellsL[id_u]; double uR = cellsR[id_u];
			double vL = cellsL[id_v]; double vR = cellsR[id_v];
			double wL = cellsL[id_w]; double wR = cellsR[id_w];
			double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
			double HtL = cellsL[id_Ht]; double HtR = cellsR[id_Ht];
			double muL = cellsL[id_mu]; double muR = cellsR[id_mu];
			double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
			double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
			double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
			double dudxL = cellsL[id_dudx]; double dudxR = cellsR[id_dudx];
			double dudyL = cellsL[id_dudy]; double dudyR = cellsR[id_dudy];
			double dudzL = cellsL[id_dudz]; double dudzR = cellsR[id_dudz];
			double dvdxL = cellsL[id_dvdx]; double dvdxR = cellsR[id_dvdx];
			double dvdyL = cellsL[id_dvdy]; double dvdyR = cellsR[id_dvdy];
			double dvdzL = cellsL[id_dvdz]; double dvdzR = cellsR[id_dvdz];
			double dwdxL = cellsL[id_dwdx]; double dwdxR = cellsR[id_dwdx];
			double dwdyL = cellsL[id_dwdy]; double dwdyR = cellsR[id_dwdy];
			double dwdzL = cellsL[id_dwdz]; double dwdzR = cellsR[id_dwdz];
			
			double pL_max = cellsL[id_pMax]; double pR_max = cellsR[id_pMax];
			double pL_min = cellsL[id_pMin]; double pR_min = cellsR[id_pMin];
			double uL_max = cellsL[id_uMax]; double uR_max = cellsR[id_uMax];
			double uL_min = cellsL[id_uMin]; double uR_min = cellsR[id_uMin];
			double vL_max = cellsL[id_vMax]; double vR_max = cellsR[id_vMax];
			double vL_min = cellsL[id_vMin]; double vR_min = cellsR[id_vMin];
			double wL_max = cellsL[id_wMax]; double wR_max = cellsR[id_wMax];
			double wL_min = cellsL[id_wMin]; double wR_min = cellsR[id_wMin];
			// double TL_max = cellsL[id_TMax]; double TR_max = cellsR[id_TMax];
			// double TL_min = cellsL[id_TMin]; double TR_min = cellsR[id_TMin];
			
			// 리미터 관련
			double gradLim_pL = cellsL[id_lim_p]; double gradLim_pR = cellsR[id_lim_p];
			double gradLim_uL = cellsL[id_lim_u]; double gradLim_uR = cellsR[id_lim_u];
			double gradLim_vL = cellsL[id_lim_v]; double gradLim_vR = cellsR[id_lim_v];
			double gradLim_wL = cellsL[id_lim_w]; double gradLim_wR = cellsR[id_lim_w];
			
			
			
			double curvatureL[nCurv+1], curvatureR[nCurv+1];
			double alpha_VFL[nCurv+1], alpha_VFR[nCurv+1];
			for(int i=0; i<nCurv; ++i){
				curvatureL[i] = cellsL[id_curvature[i]];
				curvatureR[i] = cellsR[id_curvature[i]];
				alpha_VFL[i] = cellsL[id_alpha_VF[i]];
				alpha_VFR[i] = cellsR[id_alpha_VF[i]];
			}
			double drhodpL = cellsL[id_drhodp]; double drhodpR = cellsR[id_drhodp];
			double drhodTL = cellsL[id_drhodT]; double drhodTR = cellsR[id_drhodT];
			double dHtdpL = cellsL[id_dHtdp]; double dHtdpR = cellsR[id_dHtdp];
			double dHtdTL = cellsL[id_dHtdT]; double dHtdTR = cellsR[id_dHtdT];
			double drhodYL[nSp]; double drhodYR[nSp];
			double dHtdYL[nSp]; double dHtdYR[nSp];
			// double YL_max[nSp]; double YR_max[nSp];
			// double YL_min[nSp]; double YR_min[nSp];
			double YL[nSp]; double YR[nSp];
			for(int i=0; i<nSp-1; ++i){
				drhodYL[i] = cellsL[id_drhodY[i]]; drhodYR[i] = cellsR[id_drhodY[i]];
				dHtdYL[i] = cellsL[id_dHtdY[i]]; dHtdYR[i] = cellsR[id_dHtdY[i]];
				// YL_max[i] = cellsL[id_YMax[i]]; YR_max[i] = cellsR[id_YMax[i]];
				// YL_min[i] = cellsL[id_YMin[i]]; YR_min[i] = cellsR[id_YMin[i]];
				YL[i] = cellsL[id_Y[i]]; YR[i] = cellsR[id_Y[i]];
			}
			
			
			double rhoFL = faces[id_rhoL]; double rhoFR = faces[id_rhoR];
			double uFL = faces[id_uL]; double uFR = faces[id_uR];
			double vFL = faces[id_vL]; double vFR = faces[id_vR];
			double wFL = faces[id_wL]; double wFR = faces[id_wR];
			double pFL = faces[id_pL]; double pFR = faces[id_pR];
			double HtFL = faces[id_HtL]; double HtFR = faces[id_HtR];
			double YFL[nSp]; double YFR[nSp];
			for(int i=0; i<nSp-1; ++i){
				YFL[i] = faces[id_YL[i]]; YFR[i] = faces[id_YR[i]];
			}
			double cFL = faces[id_cL]; double cFR = faces[id_cR];
			double muFL = faces[id_muL]; double muFR = faces[id_muR];
			
			double muF = wdL*muL+wdR*muR;
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
			
			double dtrho = dt*(wdL/rhoL+wdR/rhoR);
				
			
			// // double chat= wdL*cL+wdR*cR;
			// double chat= 0.5*cL+0.5*cR;
			// // double chat= wdL*cFL+wdR*cFR;
			// // double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
			// double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL)+0.5*(uR*uR+vR*vR+wR*wR));
			// double ML = UnL/chat;
			// double MR = UnR/chat;
			// double Mcy = min(1.0,KLR/chat);
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// // double Unhat = 0.5*(UnL+UnR);
			// // double rhohat = wdL*rhoL+wdR*rhoR;
			// double rhohat = 0.5*rhoL+0.5*rhoR;
			// double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );
			// double MLP = 0.5*(ML+abs(ML));
			// if( abs(ML) < 1.0 ) {
				// MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				// MLP += 0.125*(ML*ML-1.0)*(ML*ML-1.0);
			// }
			// double MRM = 0.5*(MR-abs(MR));
			// if( abs(MR) < 1.0 ) {
				// MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				// MRM -= 0.125*(MR*MR-1.0)*(MR*MR-1.0);
			// }
			// // SLAU
			// double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			// double D_L = ML+(1.0-g_c)*abs(ML);
			// double D_R = MR-(1.0-g_c)*abs(MR);
			// double D_rho = Mbar*g_c;
			// double MLPL = wdL*(D_L+D_rho);
			// double MRMR = wdR*(D_R-D_rho);



			// // 영린개발스킴0
			// double UnF = 0.5*(UnL+UnR) - dtrho*(pR-pL)/dLR;
			// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*nvec[0];
			// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*nvec[1];
			// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*nvec[2];
			// double dp_coeff_Un = dtrho/dLR;
			// double dp_coeff_thm = 0.0;
			// double WUL = 0.5; double WUR = 0.5;
			// double WpL = 0.5; double WpR = 0.5;
			// double pF = 0.5*pL + 0.5*pR;
			
			
			
			
			// 영린개발스킴1
			// double WUL = 0.5;
			// double WUR = 0.5;
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
			// double UnF = (MLP+MRM)*chat;
			
			// //----
			// double WUL = 0.25*(1.0+ML)*(1.0+ML); 
			// if(ML>0.5) WUL = (-4.0)*ML*ML*ML+(8.25)*ML*ML+(-4.5)*ML+1.25;
			// if(ML>1.0) WUL = 1.0;
			// if(ML<-1.0) WUL = 0.0;
			// double WUR = 0.25*(1.0-MR)*(1.0-MR);
			// if(MR<-0.5) WUR = (4.0)*MR*MR*MR+(8.25)*MR*MR+(4.5)*MR+1.25;
			// if(MR>1.0) WUR = 0.0;
			// if(MR<-1.0) WUR = 1.0;
			// // WUL = (w_phi*WUL+(1.0-w_phi)*wdL);
			// // WUR = (w_phi*WUR+(1.0-w_phi)*wdR);
			// double UnF = WUL*UnL+WUR*UnR;
			// //----
			
			// // ----
			// double Weps = 4.0;
			// double WUL = 1.0/(1.0+exp(-(ML*Weps)));
			// double WUR = 1.0/(1.0+exp(+(MR*Weps)));
			// double UnF = WUL*UnL+WUR*UnR;
			// // ----
			
			// //----
			// double signUnL = (UnL>0.0 ? 1.0 : -1.0);
			// double signUnR = (UnR>0.0 ? 1.0 : -1.0);
			// double WUL = 0.5*(1.0 + ( (1.0-g_c) + g_c*rhoL/(rhoL+rhoR) )*signUnL);
			// double WUR = 0.5*(1.0 - ( (1.0-g_c) - g_c*rhoR/(rhoL+rhoR) )*signUnR);
			// // WUL = (w_phi*WUL+(1.0-w_phi)*wdL);
			// // WUR = (w_phi*WUR+(1.0-w_phi)*wdR);
			// double UnF = WUL*UnL + WUR*UnR;
			// //----
			
			double chat_YYL = 0.5*(cL+cR);
			// double chat_YYL = sqrt(cL*cR);
			// double chat_YYL = 1.0/(0.5*(1.0/cL+1.0/cR));
			// double rhohat_YYL = sqrt(rhoL*rhoR);
			double rhohat_YYL = 0.5*(rhoL+rhoR);
			double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
			double Mcy = min(1.0,KLR/chat_YYL);
			double phi_c = (1.0-Mcy)*(1.0-Mcy);
			double UnL_YYL = UnL;
			double UnR_YYL = UnR;
			
			// // // high-order
			// // UnL_YYL += gradLim_uL*(dudxL*xyzLF[0] + dudyL*xyzLF[1] + dudzL*xyzLF[2])*nvec[0];
			// // UnL_YYL += gradLim_vL*(dvdxL*xyzLF[0] + dvdyL*xyzLF[1] + dvdzL*xyzLF[2])*nvec[1];
			// // UnL_YYL += gradLim_wL*(dwdxL*xyzLF[0] + dwdyL*xyzLF[1] + dwdzL*xyzLF[2])*nvec[2];
			// // UnR_YYL += gradLim_uR*(dudxR*xyzRF[0] + dudyR*xyzRF[1] + dudzR*xyzRF[2])*nvec[0];
			// // UnR_YYL += gradLim_vR*(dvdxR*xyzRF[0] + dvdyR*xyzRF[1] + dvdzR*xyzRF[2])*nvec[1];
			// // UnR_YYL += gradLim_wR*(dwdxR*xyzRF[0] + dwdyR*xyzRF[1] + dwdzR*xyzRF[2])*nvec[2];
			// // skewness
			// UnL_YYL += gradLim_uL*(dudxL*LNv[0] + dudyL*LNv[1] + dudzL*LNv[2])*nvec[0];
			// UnL_YYL += gradLim_vL*(dvdxL*LNv[0] + dvdyL*LNv[1] + dvdzL*LNv[2])*nvec[1];
			// UnL_YYL += gradLim_wL*(dwdxL*LNv[0] + dwdyL*LNv[1] + dwdzL*LNv[2])*nvec[2];
			// UnR_YYL += gradLim_uR*(dudxR*RNv[0] + dudyR*RNv[1] + dudzR*RNv[2])*nvec[0];
			// UnR_YYL += gradLim_vR*(dvdxR*RNv[0] + dvdyR*RNv[1] + dvdzR*RNv[2])*nvec[1];
			// UnR_YYL += gradLim_wR*(dwdxR*RNv[0] + dwdyR*RNv[1] + dwdzR*RNv[2])*nvec[2];
			
			double ML_YYL = UnL_YYL/chat_YYL;
			double MR_YYL = UnR_YYL/chat_YYL;
			double g_YYL = 1.0 + max( min(ML_YYL,0.0), -1.0 )*min( max(MR_YYL,0.0), 1.0 );
			double D_L_YYL = UnL_YYL+(1.0-g_YYL)*abs(UnL_YYL);
			double D_R_YYL = UnR_YYL-(1.0-g_YYL)*abs(UnR_YYL);
			double D_rho_YYL = g_YYL*( rhoL*abs(UnL_YYL)+rhoR*abs(UnR_YYL) ) / ( rhoL + rhoR );
			double UnF = wdL*(D_L_YYL+D_rho_YYL) + wdR*(D_R_YYL-D_rho_YYL);
			
			double weightThrmL = 1.0/(1.0+exp(-ML_YYL*3.5));
			double weightThrmR = 1.0/(1.0+exp(MR_YYL*3.5));
			
			UnF = weightThrmL*UnL + weightThrmR*UnR;
			
			
			// 추가적 보정텀
			double dp_coeff_thm = 0.5*phi_c/chat_YYL/rhohat_YYL;
			double dp_coeff_Un = dAlpha * dt*(wdL/rhoL+wdR/rhoR)/dLR;
			UnF -= dp_coeff_thm*(pR-pL);
			UnF -= dp_coeff_Un*(pR-pL);
			UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*dAlpha*nLR[0];
			UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*dAlpha*nLR[1];
			UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*dAlpha*nLR[2];
			
			double WUL = 0.5*(1.0+(1.0-g_YYL)*(UnL>0.0 ? 1.0 : -1.0) + 
				g_YYL*( rhoL*(UnL>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
			double WUR = 0.5*(1.0-(1.0-g_YYL)*(UnR>0.0 ? 1.0 : -1.0) - 
				g_YYL*( rhoR*(UnR>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
				
				
			WUL = weightThrmL;
			WUR = weightThrmR;
				
			
			double PLP = 0.5*(1.0 + ( ML_YYL>0.0 ? 1.0 : -1.0 ) );
			if( abs(ML_YYL) < 1.0 ) {
				PLP = 0.25*(ML_YYL+1.0)*(ML_YYL+1.0)*(2.0-ML_YYL);
				PLP += 0.1875*ML_YYL*(ML_YYL*ML_YYL-1.0)*(ML_YYL*ML_YYL-1.0);
			} 
			double PRM = 0.5*(1.0 - ( MR_YYL>0.0 ? 1.0 : -1.0 ) );
			if( abs(MR_YYL) < 1.0 ) {
				PRM = 0.25*(MR_YYL-1.0)*(MR_YYL-1.0)*(2.0+MR_YYL);
				PRM -= 0.1875*MR_YYL*(MR_YYL*MR_YYL-1.0)*(MR_YYL*MR_YYL-1.0);
			} 
			// double WpL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
			// double WpR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
			double WpL = wdL + 0.5*(PLP-PRM);
			double WpR = wdR - 0.5*(PLP-PRM);
			double pL_YYL = pL;
			double pR_YYL = pR;
			
			// // // high-order
			// // pL_YYL += gradLim_pL*(dpdxL*xyzLF[0] + dpdyL*xyzLF[1] + dpdzL*xyzLF[2]);
			// // pR_YYL += gradLim_pR*(dpdxR*xyzRF[0] + dpdyR*xyzRF[1] + dpdzR*xyzRF[2]);
			// // skewness
			// pL_YYL += gradLim_pL*(dpdxL*LNv[0] + dpdyL*LNv[1] + dpdzL*LNv[2]);
			// pR_YYL += gradLim_pR*(dpdxR*RNv[0] + dpdyR*RNv[1] + dpdzR*RNv[2]);
			
			PLP = 1.0/(1.0+exp(-ML_YYL*3.5));
			PRM = 1.0/(1.0+exp(MR_YYL*3.5));
			
			double pF = PLP*pL_YYL+PRM*pR_YYL;
			// double pF = (0.5*pL_YYL+0.5*pR_YYL) + 0.5*(PLP-PRM)*(pL_YYL-pR_YYL) + 
						// (PLP+PRM-1.0)*KLR*(wdL*rhoL+wdR*rhoR)*(wdL*cL+wdR*cR);
						// Mcy*(PLP+PRM-1.0)*(pL_YYL+pR_YYL) - 
						// 0.0;//0.5*(rhohat_YYL*(chat_YYL-abs(UnF))*(UnR-UnL));
			
			// double SL = min(UnL-cL,UnF-chat_YYL);
			// double SR = max(UnR+cR,UnF+chat_YYL);
			double SL = min(UnL-cL,0.5*(UnL+UnR)-0.5*(cL+cR));
			double SR = max(UnR+cR,0.5*(UnL+UnR)+0.5*(cL+cR));
			double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						(rhoR*(SR-UnR)-rhoL*(SL-UnL));
			double pStar = rhoL*(UnL-SL)*(UnL-SM)+pL; // rhoR*(UnR-SR)*(UnR-SM)+pR;
			if(SL<=0.0 && 0.0<SM){ pF = pStar; }
			else if(SM<=0.0 && 0.0<=SR){ pF = pStar; }
			else if(SL>0.0){ pF = pL; }
			else if(SR<0.0){ pF = pR; }
			
			// pF = 0.5*pL_YYL+0.5*pR_YYL;
			// pF = PLP*pL_YYL+PRM*pR_YYL;
			
			WpL = PLP; WpR = PRM;
			// double diffDPU = 0.5*Mcy*PLP*PRM*0.5/chat;
			double diffDPU = 0.0;
			double diffDPP = 0.0;
			double WpUL = 0.0;
			double WpUR = 0.0;
			if(
			(SL<=0.0 && 0.0<SM) || (SM<=0.0 && 0.0<=SR)
			){
				double tmp1 = (rhoR*(SR-UnR)-rhoL*(SL-UnL));
				double tmp2 = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR);
				
				WpL = 1.0 - rhoL*(UnL-SL)/(rhoR*(SR-UnR)-rhoL*(SL-UnL));
				WpR = 1.0 - rhoR*(UnR-SR)/(rhoR*(SR-UnR)-rhoL*(SL-UnL));
				
				double SMUL = (-rhoL*(SL-UnL)+rhoL*UnL)/tmp1;
				SMUL -= (+rhoL)*tmp2/tmp1/tmp1;
				double SMUR = (+rhoR*(SR-UnR)-rhoR*UnR)/tmp1;
				SMUR -= (-rhoR)*tmp2/tmp1/tmp1;
				
				WpUL = rhoL*(UnL-SM)+rhoL*(UnL-SL)-rhoL*(UnL-SL)*SMUL;
				WpUR = rhoR*(UnR-SM)+rhoR*(UnR-SR)-rhoR*(UnR-SR)*SMUR;
				
				
			}
			else if(SL>0.0){ WpL = 1.0; WpR = 0.0; }
			else if(SR<0.0){ WpL = 0.0; WpR = 1.0; }
			
			
			
			
			double weiL = 1.0; double weiR = 0.0;
			double rhoF = rhoFL;
			double uF = uFL;
			double vF = vFL;
			double wF = wFL;
			double HtF = HtFL;
			double YF[nSp];
			for(int i=0; i<nSp-1; ++i){
				YF[i] = YFL[i];
			}
			if(UnF<0.0){
				weiL = 0.0; weiR = 1.0;
				rhoF = rhoFR;
				uF = uFR;
				vF = vFR;
				wF = wFR;
				HtF = HtFR;
				for(int i=0; i<nSp-1; ++i){
					YF[i] = YFR[i];
				}
			}
			
			
			
			
			
			
			
			
			int iter = 0;
			
			double drho_dUn_dpL = weiL*drhodpL*UnF+rhoF*dp_coeff_Un +rhoF*dp_coeff_thm;
			double drho_dUn_dpR = weiR*drhodpR*UnF-rhoF*dp_coeff_Un -rhoF*dp_coeff_thm;
			double rhoF_dUn_duL = rhoF*(WUL*nvec[0]) +rhoF*diffDPU*nvec[0];
			double rhoF_dUn_duR = rhoF*(WUR*nvec[0]) -rhoF*diffDPU*nvec[0];
			double rhoF_dUn_dvL = rhoF*(WUL*nvec[1]) +rhoF*diffDPU*nvec[1];
			double rhoF_dUn_dvR = rhoF*(WUR*nvec[1]) -rhoF*diffDPU*nvec[1];
			double rhoF_dUn_dwL = rhoF*(WUL*nvec[2]) +rhoF*diffDPU*nvec[2];
			double rhoF_dUn_dwR = rhoF*(WUR*nvec[2]) -rhoF*diffDPU*nvec[2];
			
			fluxA_LL[iter] += CN_coeff * drho_dUn_dpL*area; 
			fluxA_LR[iter] += CN_coeff * drho_dUn_dpR*area;
			fluxA_RR[iter] -= CN_coeff * drho_dUn_dpR*area; 
			fluxA_RL[iter] -= CN_coeff * drho_dUn_dpL*area; 
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF)*area; 
			fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF)*area;
			fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF)*area; 
			fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF)*area; 
				fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF)*area;
				fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF)*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF)*area;
				++iter;
			}
			

			
			
			fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*uF +WpL*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*uF +WpR*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*uF +WpR*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*uF +WpL*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL +WpUL*nvec[0]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR +WpUR*nvec[0]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR +WpUR*nvec[0]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL +WpUL*nvec[0]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*uF +WpUL*nvec[1]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*uF +WpUR*nvec[1]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*uF +WpUR*nvec[1]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*uF +WpUL*nvec[1]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*uF +WpUL*nvec[2]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*uF +WpUR*nvec[2]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*uF +WpUR*nvec[2]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*uF +WpUL*nvec[2]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*uF)*area; 
			fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*uF)*area;
			fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*uF)*area; 
			fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*uF)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*uF)*area; 
				fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*uF)*area;
				fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*uF)*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*uF)*area;
				++iter;
			}
			
			
			
			
			fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*vF +WpL*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*vF +WpR*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*vF +WpR*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*vF +WpL*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*vF +WpUL*nvec[0]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*vF +WpUR*nvec[0]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*vF +WpUR*nvec[0]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*vF +WpUL*nvec[0]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL +WpUL*nvec[1]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR +WpUR*nvec[1]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR +WpUR*nvec[1]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL +WpUL*nvec[1]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF +WpUL*nvec[2]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF +WpUR*nvec[2]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF +WpUR*nvec[2]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF +WpUL*nvec[2]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*vF)*area; 
			fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*vF)*area;
			fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*vF)*area; 
			fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*vF)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*vF)*area; 
				fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*vF)*area;
				fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*vF)*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*vF)*area;
				++iter;
			}
			
			
			
			
			fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*wF +WpL*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*wF +WpR*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*wF +WpR*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*wF +WpL*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*wF +WpUL*nvec[0]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*wF +WpUR*nvec[0]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*wF +WpUR*nvec[0]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*wF +WpUL*nvec[0]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*wF +WpUL*nvec[1]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*wF +WpUR*nvec[1]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*wF +WpUR*nvec[1]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*wF +WpUL*nvec[1]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL +WpUL*nvec[2]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR +WpUR*nvec[2]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR +WpUR*nvec[2]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL +WpUL*nvec[2]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*wF)*area; 
			fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*wF)*area;
			fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*wF)*area; 
			fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*wF)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*wF)*area; 
				fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*wF)*area;
				fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*wF)*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*wF)*area;
				++iter;
			}
			
			
			
			
			fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*HtF +rhoF*UnF*weiL*dHtdpL)*area; 
			fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*HtF +rhoF*UnF*weiR*dHtdpR)*area;
			fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*HtF +rhoF*UnF*weiR*dHtdpR)*area; 
			fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*HtF +rhoF*UnF*weiL*dHtdpL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*HtF+rhoF*UnF*weiL*uL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*HtF+rhoF*UnF*weiR*uR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*HtF+rhoF*UnF*weiR*uR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*HtF+rhoF*UnF*weiL*uL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*HtF+rhoF*UnF*weiL*vL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*HtF+rhoF*UnF*weiR*vR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*HtF+rhoF*UnF*weiR*vR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*HtF+rhoF*UnF*weiL*vL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*HtF+rhoF*UnF*weiL*wL)*area; 
			fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*HtF+rhoF*UnF*weiR*wR)*area;
			fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*HtF+rhoF*UnF*weiR*wR)*area; 
			fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*HtF+rhoF*UnF*weiL*wL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*HtF+rhoF*UnF*weiL*dHtdTL)*area; 
			fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*HtF+rhoF*UnF*weiR*dHtdTR)*area;
			fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*HtF+rhoF*UnF*weiR*dHtdTR)*area; 
			fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*HtF+rhoF*UnF*weiL*dHtdTL)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[i]*UnF*HtF+rhoF*UnF*weiL*dHtdYL[i])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[i]*UnF*HtF+rhoF*UnF*weiR*dHtdYR[i])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[i]*UnF*HtF+rhoF*UnF*weiR*dHtdYR[i])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[i]*UnF*HtF+rhoF*UnF*weiL*dHtdYL[i])*area;
				++iter;
			}
			
			
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff * (drho_dUn_dpL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (drho_dUn_dpR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (drho_dUn_dpR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (drho_dUn_dpL*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (weiL*drhodTL*UnF*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (weiR*drhodTR*UnF*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (weiR*drhodTR*UnF*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (weiL*drhodTL*UnF*YF[i])*area;
				++iter;
				
				for(int j=0; j<nSp-1; ++j){
					fluxA_LL[iter] += CN_coeff_Y * (weiL*drhodYL[j]*UnF*YF[i])*area; 
					fluxA_LR[iter] += CN_coeff_Y * (weiR*drhodYR[j]*UnF*YF[i])*area;
					fluxA_RR[iter] -= CN_coeff_Y * (weiR*drhodYR[j]*UnF*YF[i])*area; 
					fluxA_RL[iter] -= CN_coeff_Y * (weiL*drhodYL[j]*UnF*YF[i])*area;
					if(i==j){
						fluxA_LL[iter] += CN_coeff_Y * (rhoF*UnF*weiL)*area; 
						fluxA_LR[iter] += CN_coeff_Y * (rhoF*UnF*weiR)*area;
						fluxA_RR[iter] -= CN_coeff_Y * (rhoF*UnF*weiR)*area; 
						fluxA_RL[iter] -= CN_coeff_Y * (rhoF*UnF*weiL)*area;
					}
					++iter;
				}
			}
			
			
			
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
			
			iter = nEq*1+1;
			fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			iter = nEq*2+2;
			fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			iter = nEq*3+3;
			fluxA_LL[iter] += visc_coeff; fluxA_LR[iter] -= visc_coeff;
			fluxA_RR[iter] += visc_coeff; fluxA_RL[iter] -= visc_coeff;
			
			iter = nEq*4+1;
			fluxA_LL[iter] += visc_coeff*ubar; fluxA_LR[iter] -= visc_coeff*ubar;
			fluxA_RR[iter] += visc_coeff*ubar; fluxA_RL[iter] -= visc_coeff*ubar;
			
			iter = nEq*4+2;
			fluxA_LL[iter] += visc_coeff*vbar; fluxA_LR[iter] -= visc_coeff*vbar;
			fluxA_RR[iter] += visc_coeff*vbar; fluxA_RL[iter] -= visc_coeff*vbar;
			
			iter = nEq*4+3;
			fluxA_LL[iter] += visc_coeff*wbar; fluxA_LR[iter] -= visc_coeff*wbar;
			fluxA_RR[iter] += visc_coeff*wbar; fluxA_RL[iter] -= visc_coeff*wbar;
			
			
			
			
			
			
			// 디퓨젼 B
			fluxB[1] += muF*(
						(dAlpha*(uR-uL)/dLR + dudxF*nvec[0] + dvdxF*nvec[1] + dwdxF*nvec[2]) - 
						2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[0] 
						)*area;
			// fluxB[1] -= nvec[0] * 2.0/3.0 * rhoF * tkei;
			fluxB[2] += muF*(
						(dAlpha*(vR-vL)/dLR + dudyF*nvec[0] + dvdyF*nvec[1] + dwdyF*nvec[2]) - 
						2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[1]
						)*area;
			// fluxB[2] -= nvec[1] * 2.0/3.0 * rhoF * tkei;
			fluxB[3] += muF*(
						(dAlpha*(wR-wL)/dLR + dudzF*nvec[0] + dvdzF*nvec[1] + dwdzF*nvec[2]) - 
						2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[2]
						)*area;
			// fluxB[3] -= nvec[2] * 2.0/3.0 * rhoF * tkei;
			fluxB[4] += muF*(
						dAlpha*(uR-uL)/dLR*ubar + dAlpha*(vR-vL)/dLR*vbar + dAlpha*(wR-wL)/dLR*wbar + 
						(dudxF*ubar + dudyF*vbar + dudzF*wbar)*nvec[0] +
						(dvdxF*ubar + dvdyF*vbar + dvdzF*wbar)*nvec[1] +
						(dwdxF*ubar + dwdyF*vbar + dwdzF*wbar)*nvec[2] -
						2.0/3.0*(dudxF + dvdyF + dwdzF)*ubar*nvec[0] -
						2.0/3.0*(dudxF + dvdyF + dwdzF)*vbar*nvec[1] -
						2.0/3.0*(dudxF + dvdyF + dwdzF)*wbar*nvec[2]
						)*area;
			
			// non-orthogonal
			fluxB[1] += (nvec[0]-dAlpha*nLR[0]) * muF*dudxF*area;
			fluxB[1] += (nvec[1]-dAlpha*nLR[1]) * muF*dudyF*area;
			fluxB[1] += (nvec[2]-dAlpha*nLR[2]) * muF*dudzF*area;
			
			fluxB[2] += (nvec[0]-dAlpha*nLR[0]) * muF*dvdxF*area;
			fluxB[2] += (nvec[1]-dAlpha*nLR[1]) * muF*dvdyF*area;
			fluxB[2] += (nvec[2]-dAlpha*nLR[2]) * muF*dvdzF*area;
			
			fluxB[3] += (nvec[0]-dAlpha*nLR[0]) * muF*dwdxF*area;
			fluxB[3] += (nvec[1]-dAlpha*nLR[1]) * muF*dwdyF*area;
			fluxB[3] += (nvec[2]-dAlpha*nLR[2]) * muF*dwdzF*area;
			
			fluxB[4] += (nvec[0]-dAlpha*nLR[0]) * muF*(
						dudxF*ubar+dvdxF*vbar+dwdxF*wbar)*area;
			fluxB[4] += (nvec[1]-dAlpha*nLR[1]) * muF*(
						dudyF*ubar+dvdyF*vbar+dwdyF*wbar)*area;
			fluxB[4] += (nvec[2]-dAlpha*nLR[2]) * muF*(
						dudzF*ubar+dvdzF*vbar+dwdzF*wbar)*area;
			
			
		
			
			// 표면장력
			for(int i=0; i<nCurv; ++i){
				double alpha_VFF = wdL*alpha_VFL[i]+wdR*alpha_VFR[i];
				fluxB_LL[1] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
				fluxB_LL[2] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
				fluxB_LL[3] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				fluxB_LL[4] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*UnF )*area;
				
				fluxB_RR[1] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[0] )*area;
				fluxB_RR[2] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[1] )*area;
				fluxB_RR[3] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[2] )*area;
				fluxB_RR[4] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*UnF )*area;
			}
			
			



			fluxB_LL[0] += fluxB[0]; fluxB_RR[0] -= fluxB[0];
			fluxB_LL[1] += fluxB[1]; fluxB_RR[1] -= fluxB[1];
			fluxB_LL[2] += fluxB[2]; fluxB_RR[2] -= fluxB[2];
			fluxB_LL[3] += fluxB[3]; fluxB_RR[3] -= fluxB[3];
			fluxB_LL[4] += fluxB[4]; fluxB_RR[4] -= fluxB[4];
			for(int i=0; i<nSp-1; ++i){
				fluxB_LL[5+i] += fluxB[5+i]; fluxB_RR[5+i] -= fluxB[5+i];
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
			id_p,id_u,id_v,id_w,
			id_pL,id_uL,id_vL,id_wL,id_YL,id_rhoL,id_HtL,id_muL,id_cL,
			id_pR,id_uR,id_vR,id_wR,id_YR,id_rhoR,id_HtR,id_muR,id_cR,
			id_dpdx,id_dpdy,id_dpdz,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,
			id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
			id_wd,id_rho,id_Ht,nEq,id_c,CN_coeff,CN_coeff_Y,id_Y,
			id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,id_dwdx,id_dwdy,id_dwdz,
			id_alpha,id_nLRx,id_nLRy,id_nLRz,
			id_curvature,id_alpha_VF,nCurv,surf_sigma,
			bool_zeroGradient,bool_inletOutlet](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				double dt = fields[id_dt];
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				
				dAlpha = 1.0;
				
				
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
				double rhoF = faces[id_rhoL];
				double uF = faces[id_uL];
				double vF = faces[id_vL];
				double wF = faces[id_wL];
				double pF = faces[id_pL];
				double HtF = faces[id_HtL];
				double YF[nSp];
				for(int i=0; i<nSp-1; ++i){
					YF[i] = faces[id_YL[i]];
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
				
				double muF = faces[id_muL];
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
				
				
				
			
				
				
				if(
				bool_zeroGradient[0]==true ||
				(bool_inletOutlet[0]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				){
					fluxA_LL[nEq*0+0] += CN_coeff * (drhodpL*UnF)*area; 
					fluxA_LL[nEq*1+0] += CN_coeff * (drhodpL*UnF*uF + nvec[0])*area; 
					fluxA_LL[nEq*2+0] += CN_coeff * (drhodpL*UnF*vF + nvec[1])*area; 
					fluxA_LL[nEq*3+0] += CN_coeff * (drhodpL*UnF*wF + nvec[2])*area; 
					fluxA_LL[nEq*4+0] += CN_coeff * (drhodpL*UnF*HtF+rhoF*UnF*dHtdpL)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+0] += CN_coeff_Y * (drhodpL*UnF*YF[i])*area; 
					}
					
				}
				else{
					fluxA_LL[nEq*0+0] += CN_coeff * (rhoF*dp_coeff_Un +rhoF*dp_coeff_thm)*area;
					fluxA_LL[nEq*1+0] += CN_coeff * (rhoF*dp_coeff_Un*uF +rhoF*dp_coeff_thm*uF)*area; 
					fluxA_LL[nEq*2+0] += CN_coeff * (rhoF*dp_coeff_Un*vF +rhoF*dp_coeff_thm*vF)*area; 
					fluxA_LL[nEq*3+0] += CN_coeff * (rhoF*dp_coeff_Un*wF +rhoF*dp_coeff_thm*wF)*area; 
					fluxA_LL[nEq*4+0] += CN_coeff * (rhoF*dp_coeff_Un*HtF +rhoF*dp_coeff_thm*HtF)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+0] += CN_coeff_Y * (rhoF*dp_coeff_Un*YF[i] +rhoF*dp_coeff_thm*YF[i])*area; 
					}
				}
				
				if(
				bool_zeroGradient[1]==true ||
				(bool_inletOutlet[1]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				){
					fluxA_LL[nEq*0+1] += CN_coeff * (rhoF*(nvec[0]))*area; 
					fluxA_LL[nEq*1+1] += CN_coeff * (rhoF*(nvec[0])*uF+rhoF*UnF)*area; 
					fluxA_LL[nEq*2+1] += CN_coeff * (rhoF*(nvec[0])*vF)*area; 
					fluxA_LL[nEq*3+1] += CN_coeff * (rhoF*(nvec[0])*wF)*area; 
					fluxA_LL[nEq*4+1] += CN_coeff * (rhoF*(nvec[0])*HtF+rhoF*UnF*uL)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+1] += CN_coeff_Y * (rhoF*(nvec[0])*YF[i])*area; 
					}
					
					fluxA_LL[nEq*0+2] += CN_coeff * (rhoF*(nvec[1]))*area; 
					fluxA_LL[nEq*1+2] += CN_coeff * (rhoF*(nvec[1])*uF)*area; 
					fluxA_LL[nEq*2+2] += CN_coeff * (rhoF*(nvec[1])*vF+rhoF*UnF)*area; 
					fluxA_LL[nEq*3+2] += CN_coeff * (rhoF*(nvec[1])*wF)*area; 
					fluxA_LL[nEq*4+2] += CN_coeff * (rhoF*(nvec[1])*HtF+rhoF*UnF*vL)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+2] += CN_coeff_Y * (rhoF*(nvec[1])*YF[i])*area; 
					}
					
					fluxA_LL[nEq*0+3] += CN_coeff * (rhoF*(nvec[2]))*area; 
					fluxA_LL[nEq*1+3] += CN_coeff * (rhoF*(nvec[2])*uF)*area; 
					fluxA_LL[nEq*2+3] += CN_coeff * (rhoF*(nvec[2])*vF)*area; 
					fluxA_LL[nEq*3+3] += CN_coeff * (rhoF*(nvec[2])*wF+rhoF*UnF)*area; 
					fluxA_LL[nEq*4+3] += CN_coeff * (rhoF*(nvec[2])*HtF+rhoF*UnF*wL)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+3] += CN_coeff_Y * (rhoF*(nvec[2])*YF[i])*area; 
					}
				}
				else{
					
					double visc_coeff = dAlpha*muF/dLR*area;
					fluxA_LL[nEq*1+1] += visc_coeff;
					fluxA_LL[nEq*2+2] += visc_coeff;
					fluxA_LL[nEq*3+3] += visc_coeff;
					fluxA_LL[nEq*4+1] += visc_coeff*ubar;
					fluxA_LL[nEq*4+2] += visc_coeff*vbar;
					fluxA_LL[nEq*4+3] += visc_coeff*wbar;
					
				}
				
				if(
				bool_zeroGradient[2]==true ||
				(bool_inletOutlet[2]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
				){
					fluxA_LL[nEq*0+4] += CN_coeff * (drhodTL*UnF)*area; 
					fluxA_LL[nEq*1+4] += CN_coeff * (drhodTL*UnF*uF)*area; 
					fluxA_LL[nEq*2+4] += CN_coeff * (drhodTL*UnF*vF)*area; 
					fluxA_LL[nEq*3+4] += CN_coeff * (drhodTL*UnF*wF)*area; 
					fluxA_LL[nEq*4+4] += CN_coeff * (drhodTL*UnF*HtF+rhoF*UnF*dHtdTL)*area; 
					for(int i=0; i<nSp-1; ++i){
						fluxA_LL[nEq*(5+i)+4] += CN_coeff_Y * (drhodTL*UnF*YF[i])*area; 
					}
				}
				
				for(int i=0; i<nSp-1; ++i){
					if(
					bool_zeroGradient[3+i]==true ||
					(bool_inletOutlet[3+i]==true && (uF*nvec[0]+vF*nvec[1]+wF*nvec[2]<0.0))
					){
						fluxA_LL[nEq*0+5+i] += CN_coeff * (drhodYL[i]*UnF)*area; 
						fluxA_LL[nEq*1+5+i] += CN_coeff * (drhodYL[i]*UnF*uF)*area; 
						fluxA_LL[nEq*2+5+i] += CN_coeff * (drhodYL[i]*UnF*vF)*area; 
						fluxA_LL[nEq*3+5+i] += CN_coeff * (drhodYL[i]*UnF*wF)*area; 
						fluxA_LL[nEq*4+5+i] += CN_coeff * (drhodYL[i]*UnF*HtF+rhoF*UnF*dHtdYL[i])*area; 
						for(int j=0; j<nSp-1; ++j){
							fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * (drhodYL[i]*UnF*YF[j])*area; 
							if(i==j){
								fluxA_LL[nEq*(5+j)+5+i] += CN_coeff_Y * (rhoF*UnF)*area; 
							}
						}
					
					}
				}
				
				
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
					fluxB_LL[4] -= surf_sigma[i]*curvatureL[i]*( alpha_VFF*UnF )*area;
				}
				
				
				return 0;
			});
		}
		
	}
		
	
}




