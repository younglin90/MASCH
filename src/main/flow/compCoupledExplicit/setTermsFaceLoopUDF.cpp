
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
			
			
				dAlpha = 1.0;
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
			
			// // skewness
			// double delphiL[5+nSp-1],delphiR[5+nSp-1];
			// double limGrad_phiL[5+nSp-1],limGrad_phiR[5+nSp-1];
			// double eta = sqrt(area);
			// eta = 1000.0*eta*eta*eta;
			// delphiL[1] = dudxL*LNv[0] + dudyL*LNv[1] + dudzL*LNv[2];
			// delphiR[1] = dudxR*RNv[0] + dudyR*RNv[1] + dudzR*RNv[2];
			// limGrad_phiL[1] = solver.limiter_MLP(uL,uL_max,uL_min,delphiL[1], eta);
			// limGrad_phiR[1] = solver.limiter_MLP(uR,uR_max,uR_min,delphiR[1], eta);
			
			// delphiL[2] = dvdxL*LNv[0] + dvdyL*LNv[1] + dvdzL*LNv[2];
			// delphiR[2] = dvdxR*RNv[0] + dvdyR*RNv[1] + dvdzR*RNv[2];
			// limGrad_phiL[2] = solver.limiter_MLP(vL,vL_max,vL_min,delphiL[2], eta);
			// limGrad_phiR[2] = solver.limiter_MLP(vR,vR_max,vR_min,delphiR[2], eta);
			
			// delphiL[3] = dwdxL*LNv[0] + dwdyL*LNv[1] + dwdzL*LNv[2];
			// delphiR[3] = dwdxR*RNv[0] + dwdyR*RNv[1] + dwdzR*RNv[2];
			// limGrad_phiL[3] = solver.limiter_MLP(wL,wL_max,wL_min,delphiL[3], eta);
			// limGrad_phiR[3] = solver.limiter_MLP(wR,wR_max,wR_min,delphiR[3], eta);
			
			
			

			double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
			double Mcy = min(1.0,KLR/(wdL*cL+wdR*cR));
			double phi_c = (1.0-Mcy)*(1.0-Mcy);
			double UnL_YYL = UnL;
			double UnR_YYL = UnR;
			double ML_YYL = UnL_YYL/(wdL*cL+wdR*cR);
			double MR_YYL = UnR_YYL/(wdL*cL+wdR*cR);
			double g_YYL = 1.0 + max( min(ML_YYL,0.0), -1.0 )*min( max(MR_YYL,0.0), 1.0 );
			double D_L_YYL = UnL_YYL+(1.0-g_YYL)*abs(UnL_YYL);
			double D_R_YYL = UnR_YYL-(1.0-g_YYL)*abs(UnR_YYL);
			double D_rho_YYL = g_YYL*( rhoL*abs(UnL_YYL)+rhoR*abs(UnR_YYL) ) / ( rhoL + rhoR );
			double UnF = wdL*(D_L_YYL+D_rho_YYL) + wdR*(D_R_YYL-D_rho_YYL);
			
			
			
			// skewness
			UnF += gradLim_uL*0.5*(dudxL*LNv[0] + dudyL*LNv[1] + dudzL*LNv[2])*nvec[0];
			UnF += gradLim_vL*0.5*(dvdxL*LNv[0] + dvdyL*LNv[1] + dvdzL*LNv[2])*nvec[1];
			UnF += gradLim_wL*0.5*(dwdxL*LNv[0] + dwdyL*LNv[1] + dwdzL*LNv[2])*nvec[2];
			UnF += gradLim_uR*0.5*(dudxR*RNv[0] + dudyR*RNv[1] + dudzR*RNv[2])*nvec[0];
			UnF += gradLim_vR*0.5*(dvdxR*RNv[0] + dvdyR*RNv[1] + dvdzR*RNv[2])*nvec[1];
			UnF += gradLim_wR*0.5*(dwdxR*RNv[0] + dwdyR*RNv[1] + dwdzR*RNv[2])*nvec[2];
			
			
			
			// 추가적 보정텀
			UnF -= 0.5*phi_c/(wdL*cL+wdR*cR)/(wdL*rhoL+wdR*rhoR)*(pR-pL);
			UnF -= dAlpha * dt*(wdL/rhoL+wdR/rhoR)*(pR-pL)/dLR;
			UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*dAlpha*nLR[0];
			UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*dAlpha*nLR[1];
			UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*dAlpha*nLR[2];
			
			double WUL = 0.5*(1.0+(1.0-g_YYL)*(UnL>0.0 ? 1.0 : -1.0) + 
				g_YYL*( rhoL*(UnL>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
			double WUR = 0.5*(1.0-(1.0-g_YYL)*(UnR>0.0 ? 1.0 : -1.0) - 
				g_YYL*( rhoR*(UnR>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
			double dp_coeff_thm = 0.5*phi_c/((wdL*cL+wdR*cR))/((wdL*rhoL+wdR*rhoR));
			double dp_coeff_Un = dAlpha * dt*(wdL/rhoL+wdR/rhoR)/dLR;
			
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
			double pF = (wdL*pL+wdR*pR) + 0.5*(PLP-PRM)*(pL-pR) + 
						// KLR*(wdL*rhoL+wdR*rhoR)*(wdL*cL+wdR*cR)*(PLP+PRM-1.0);
						Mcy*(PLP+PRM-1.0)*(pL+pR);
			


			pF = 0.5*(pL+pR);
			
			
			
			// skewness
			pF += gradLim_pL*0.5*(dpdxL*LNv[0] + dpdyL*LNv[1] + dpdzL*LNv[2]);
			pF += gradLim_pR*0.5*(dpdxR*RNv[0] + dpdyR*RNv[1] + dpdzR*RNv[2]);





			// double diffDPU = 0.5*Mcy*PLP*PRM*0.5/chat;
			double diffDPU = 0.0;
			
			
			
			// UnL_YYL += limGrad_phiL[1]*delphiL[1]*nvec[0];
			// UnL_YYL += limGrad_phiL[2]*delphiL[2]*nvec[1];
			// UnL_YYL += limGrad_phiL[3]*delphiL[3]*nvec[2];
			// UnR_YYL += limGrad_phiR[1]*delphiR[1]*nvec[0];
			// UnR_YYL += limGrad_phiR[2]*delphiR[2]*nvec[1];
			// UnR_YYL += limGrad_phiR[3]*delphiR[3]*nvec[2];
			// UnL_YYL += delphiL[1]*nvec[0];
			// UnL_YYL += delphiL[2]*nvec[1];
			// UnL_YYL += delphiL[3]*nvec[2];
			// UnR_YYL += delphiR[1]*nvec[0];
			// UnR_YYL += delphiR[2]*nvec[1];
			// UnR_YYL += delphiR[3]*nvec[2];
			
			// double ML_YYL = UnL_YYL/(0.5*(cFL+cFR));
			// double MR_YYL = UnR_YYL/(0.5*(cFL+cFR));
			// double g_YYL = 1.0 + max( min(ML_YYL,0.0), -1.0 )*min( max(MR_YYL,0.0), 1.0 );
			// double D_L_YYL = UnL_YYL+(1.0-g_YYL)*abs(UnL_YYL);
			// double D_R_YYL = UnR_YYL-(1.0-g_YYL)*abs(UnR_YYL);
			// double D_rho_YYL = g_YYL*( rhoFL*abs(UnL_YYL)+rhoFR*abs(UnR_YYL) ) / ( rhoFL + rhoFR );
			// double UnF = 0.5*(D_L_YYL+D_rho_YYL) + 0.5*(D_R_YYL-D_rho_YYL);
			// UnF -= 0.5*phi_c/(0.5*(cFL+cFR))/(0.5*(rhoFL+rhoFR))*(pR-pL);
			
			
			// UnF -= 0.5*phi_c/(0.5*(rhoL+rhoR))/(0.5*(cL+cR))*(pR-pL);
			// UnF -= dAlpha * dtrho*(pR-pL)/dLR;
			// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*dAlpha*nLR[0];
			// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*dAlpha*nLR[1];
			// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*dAlpha*nLR[2];
			
			// double dp_coeff_Un = dAlpha * dtrho/dLR;
			// double dp_coeff_thm = 0.5*phi_c/rhohat/chat;
			
			// double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML) < 1.0 ) {
				// PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
				// PLP += 0.1875*ML*(ML*ML-1.0)*(ML*ML-1.0);
			// } 
			// double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR) < 1.0 ) {
				// PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
				// PRM -= 0.1875*MR*(MR*MR-1.0)*(MR*MR-1.0);
			// } 
			// // double WpL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
			// // double WpR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
			// double WpL = PLP;
			// double WpR = PRM;
			// // WpL = (w_phi*WpL+(1.0-w_phi)*wdL);
			// // WpR = (w_phi*WpR+(1.0-w_phi)*wdR);
			// // double tmp_sum2 = (WpL+WpR);
			// // WpL = 0.5;//WpL/tmp_sum2;
			// // WpR = 0.5;//WpR/tmp_sum2;
			// // double pF = WpL*pL + WpR*pR;
			// // double pF = 0.5*(pL+pR) + 0.5*(PLP-PRM)*(pL-pR) + (1.0-phi_c)*(PLP+PRM-1.0)*0.5*(pL+pR);
			
			// double pL_YYL = pL;
			// double pR_YYL = pR;
			// // // skewness
			// // delphiL[0] = dpdxL*LNv[0] + dpdyL*LNv[1] + dpdzL*LNv[2];
			// // delphiR[0] = dpdxR*RNv[0] + dpdyR*RNv[1] + dpdzR*RNv[2];
			// // limGrad_phiL[0] = solver.limiter_MLP(pL,pL_max,pL_min,delphiL[0], eta);
			// // limGrad_phiR[0] = solver.limiter_MLP(pR,pR_max,pR_min,delphiR[0], eta);
			// // // pL_YYL += limGrad_phiL[0]*delphiL[0];
			// // // pR_YYL += limGrad_phiR[0]*delphiR[0];
			// // pL_YYL += delphiL[0];
			// // pR_YYL += delphiR[0];
			// double pF = 0.5*(pL_YYL+pR_YYL) + 0.5*(PLP-PRM)*(pL-pR) + 
						// KLR*(0.5*(rhoFL+rhoFR))*(0.5*(cFL+cFR))*(PLP+PRM-1.0);
			
			
			
			// pF = 0.5*pL + 0.5*pR;
			
			
			// double weiL = 1.0; double weiR = 0.0;
			// double rhoF = rhoFL;
			// double uF = uFL;
			// double vF = vFL;
			// double wF = wFL;
			// double HtF = HtFL;
			// double YF[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// YF[i] = YFL[i];
			// }
			// if(UnF<0.0){
				// weiL = 0.0; weiR = 1.0;
				// rhoF = rhoFR;
				// uF = uFR;
				// vF = vFR;
				// wF = wFR;
				// HtF = HtFR;
				// for(int i=0; i<nSp-1; ++i){
					// YF[i] = YFR[i];
				// }
			// }
			
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
			
			
			
			
			int iter;
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
			
			
			
			
			
			
			// // 영린개발스킴2
			// double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
			// double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
			// double w = 1.0 - pow( min(preLs/preRs,preRs/preLs),2.0);
			// double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML) < 1.0 ) {
				// PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
				// // PLP += 0.1875*( (ML*ML-1.0)*(ML*ML-1.0) + ML*2.0*(ML*ML-1.0) );
			// } 
			// double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR) < 1.0 ) {
				// PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
				// // PRM -= 0.1875*( (MR*MR-1.0)*(MR*MR-1.0) + MR*2.0*(MR*MR-1.0) );
			// } 

			// double MPP = MLP+MRM;
			// double MLP_AUSM = 0.5*(MPP+abs(MPP));
			// double MRM_AUSM = 0.5*(MPP-abs(MPP));

			// double MLP_SLAU = 0.5*(D_L+D_rho);
			// double MRM_SLAU = 0.5*(D_R-D_rho);

			// double fa1 = w;

			// double MLPL_HAUS = fa1*MLP_AUSM + (1.0-fa1)*MLP_SLAU;
			// double MRMR_HAUS = fa1*MRM_AUSM + (1.0-fa1)*MRM_SLAU;

			// double mdot = rhoFL*chat*MLPL_HAUS + rhoFR*chat*MRMR_HAUS - 0.5*phi_c/chat*(pR-pL);
			
			

			// double f1L;
			// double f1R;
			// if( mdot >= 0.0 ) {
				// f1L = mdot;
				// f1R = 0.0;
			// }
			// else{
				// f1L = 0.0;
				// f1R = mdot;
			// }
	
			// // double UnF = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);
			// // // UnF -= dAlpha * dt*(pR-pL)/dLR;
			// // // UnF += dt*(wdL*dpdxL+wdR*dpdxR)*nvec[0];
			// // // UnF += dt*(wdL*dpdyL+wdR*dpdyR)*nvec[1];
			// // // UnF += dt*(wdL*dpdzL+wdR*dpdzR)*nvec[2];
			// // // // non-orthogonal
			// // // UnF -= (nvec[0]-dAlpha*nLR[0]) * dt*(wdL*dpdxL+wdR*dpdxR);
			// // // UnF -= (nvec[1]-dAlpha*nLR[1]) * dt*(wdL*dpdyL+wdR*dpdyR);
			// // // UnF -= (nvec[2]-dAlpha*nLR[2]) * dt*(wdL*dpdzL+wdR*dpdzR);
			// // // double dp_coeff_Un;
			// // // double dp_coeff_thm;
			// // double WUL = 0.5;
			// // double WUR = 0.5;
			// // if(UnF>=0.0){
				// // UnF = UnF/rhoL;
				// // // double signUnL = (UnL>0.0 ? 1.0 : -1.0);
				// // // double signUnR = (UnR>0.0 ? 1.0 : -1.0);
				// // // WUL = (rhoL*0.5*(1.0+(1.0-g_c)*signUnL+rhoL*signUnL/(rhoL+rhoR)*g_c))/rhoL;
				// // // WUR = (rhoR*0.5*(1.0-(1.0-g_c)*signUnR-rhoR*signUnR/(rhoL+rhoR)*g_c))/rhoL;
				
				// // // dp_coeff_Un = dAlpha * dt/dLR/rhoL;
				// // // dp_coeff_thm = 0.5*phi_c/chat/rhoL;
			// // }
			// // else{
				// // UnF = UnF/rhoR;
				// // // double signUnL = (UnL>0.0 ? 1.0 : -1.0);
				// // // double signUnR = (UnR>0.0 ? 1.0 : -1.0);
				// // // WUL = (rhoL*0.5*(1.0+(1.0-g_c)*signUnL+rhoL*signUnL/(rhoL+rhoR)*g_c))/rhoR;
				// // // WUR = (rhoR*0.5*(1.0-(1.0-g_c)*signUnR-rhoR*signUnR/(rhoL+rhoR)*g_c))/rhoR;
				
				// // // dp_coeff_Un = dAlpha * dt/dLR/rhoR;
				// // // dp_coeff_thm = 0.5*phi_c/chat/rhoR;
			// // }
			// double WpL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
			// double WpR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
			// // double WpL = 0.5 - 0.5*(KLR/chat)*PLP*PRM*0.5/chat*(UnR-UnL) + 
						// // 0.6*(KLR/chat)*0.5*(PLP+PRM-1.0) + 0.5*(PLP-PRM);
			// // double WpR = 0.5 - 0.5*(KLR/chat)*PLP*PRM*0.5/chat*(UnR-UnL) + 
						// // 0.6*(KLR/chat)*0.5*(PLP+PRM-1.0) - 0.5*(PLP-PRM);
			// // WpL = (w_phi*WpL+(1.0-w_phi)*wdL);
			// // WpR = (w_phi*WpR+(1.0-w_phi)*wdR);
			// // double tmp_sum2 = (WpL+WpR);
			// // WpL = 0.5;//WpL/tmp_sum2;
			// // WpR = 0.5;//WpR/tmp_sum2;
			// double pF = WpL*pL + WpR*pR;
			
			
			
			// // pF = 0.5*pL + 0.5*pR;
			
			
			// int iter;
			// double fluxB[nEq];
			
			// // 컨벡티브 B
			// fluxB[0] = -( f1L + f1R )*area;
			// fluxB[1] = -( f1L*uFL + f1R*uFR + pF*nvec[0] )*area;
			// fluxB[2] = -( f1L*vFL + f1R*vFR + pF*nvec[1] )*area;
			// fluxB[3] = -( f1L*wFL + f1R*wFR + pF*nvec[2] )*area;
			// fluxB[4] = -( f1L*HtFL + f1R*HtFR )*area;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[5+i] = -( f1L*YFL[i] + f1R*YFR[i] )*area;
			// }
			
			
				
				
				
				
			// // double rhoUnF = rhoFL*chat*MLPL + rhoFR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);
				
			
			
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
				
				
				dAlpha = 1.0;
				
				
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
				double cbar = cL;//0.5*(cL+cF);
				double Mdash = min(1.0,KLR/cbar);
				double chi = (1.0-Mdash)*(1.0-Mdash);
				double rhohat = 0.5*rhoL+0.5*rhoF;
				double dp_coeff_thm = 0.5*chi/rhohat/cbar;
				
				
				
				fluxB_LL[0] -= ( rhoF*UnF )*area;
				fluxB_LL[1] -= ( rhoF*UnF*uF + pF*nvec[0] )*area;
				fluxB_LL[2] -= ( rhoF*UnF*vF + pF*nvec[1] )*area;
				fluxB_LL[3] -= ( rhoF*UnF*wF + pF*nvec[2] )*area;
				fluxB_LL[4] -= ( rhoF*UnF*HtF )*area;
				for(int i=0; i<nSp-1; ++i){
					fluxB_LL[5+i] -= ( rhoF*UnF*YF[i] )*area;
				}
				
					
				
				
				
				
				fluxB_LL[1] += muF*(
							(dAlpha*(uF-uL)/dLR + dudxF*nvec[0] + dvdxF*nvec[1] + dwdxF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[0] 
							)*area;
				fluxB_LL[2] += muF*(
							(dAlpha*(vF-vL)/dLR + dudyF*nvec[0] + dvdyF*nvec[1] + dwdyF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[1]
							)*area;
				fluxB_LL[3] += muF*(
							(dAlpha*(wF-wL)/dLR + dudzF*nvec[0] + dvdzF*nvec[1] + dwdzF*nvec[2]) - 
							2.0/3.0*(dudxF + dvdyF + dwdzF)*nvec[2]
							)*area;
				fluxB_LL[4] += muF*(
							dAlpha*(uF-uL)/dLR*ubar + dAlpha*(vF-vL)/dLR*vbar + dAlpha*(wF-wL)/dLR*wbar + 
							(dudxF*ubar + dudyF*vbar + dudzF*wbar)*nvec[0] +
							(dvdxF*ubar + dvdyF*vbar + dvdzF*wbar)*nvec[1] +
							(dwdxF*ubar + dwdyF*vbar + dwdzF*wbar)*nvec[2] -
							2.0/3.0*(dudxF + dvdyF + dwdzF)*ubar*nvec[0] -
							2.0/3.0*(dudxF + dvdyF + dwdzF)*vbar*nvec[1] -
							2.0/3.0*(dudxF + dvdyF + dwdzF)*wbar*nvec[2]
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
					fluxB_LL[1] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
					fluxB_LL[2] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
					fluxB_LL[3] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				}
				
				
				return 0;
			});
		}
		
	}
		
	
}


