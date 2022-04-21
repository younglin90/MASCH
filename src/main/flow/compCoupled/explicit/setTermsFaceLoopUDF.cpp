
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
    
	int id_vol = controls.getId_cellVar("volume");
	
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
		id_xLF,id_yLF,id_zLF,id_xRF,id_yRF,id_zRF,
        id_dYdx,id_dYdy,id_dYdz,id_vol](
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
				
			
            
            // =================================================================
            // turbulent
            // WALE
            double C_I = 0.007;
            double C_eta = 0.92;
            
            double avgLij[3][3];
            double avgMij[3][3];
            double avgSij[3][3];
            double sumRho = 0.0;
            double sumV = 0.0;
            double sumU = 0.0;
            double sumW = 0.0;
            double sumSabs = 0.0;
            double sumVol = 0.0;
            double volL = cellsL[id_vol];
            double volR = cellsR[id_vol];
            double SabsL = 0.0;
            double SabsR = 0.0;
            {
                // calculate strain rate tensor, Sij
                double SijL[3][3];
                SijL[0][0] = dudxL;
                SijL[0][1] = 0.5*(dudyL+dvdxL);
                SijL[0][2] = 0.5*(dudzL+dwdxL);
                SijL[1][0] = SijL[0][1];
                SijL[1][1] = dvdyL;
                SijL[1][2] = 0.5*(dvdzL+dwdyL);
                SijL[2][0] = SijL[0][2];
                SijL[2][1] = SijL[1][2];
                SijL[2][2] = dwdzL;

                // calculate S_ij * S_ij
                double SSL = SijL[0][0]*SijL[0][0] + SijL[1][1]*SijL[1][1] + SijL[2][2]*SijL[2][2] +
                    2.0*(SijL[0][1]*SijL[0][1] + SijL[0][2]*SijL[0][2] + SijL[1][2]*SijL[1][2]);
                
                SabsL = sqrt(2.0*SSL);
                
                // filtering process f(rho*ui*uj), f(rho), f(ui)
                // double LijL[3][3];
                avgLij[0][0] = rhoL * uL * uL * volL;
                avgLij[0][1] = rhoL * uL * vL * volL;
                avgLij[0][2] = rhoL * uL * wL * volL;
                avgLij[1][1] = rhoL * vL * vL * volL;
                avgLij[1][2] = rhoL * vL * wL * volL;
                avgLij[2][2] = rhoL * wL * wL * volL;
                    
                // double MijL[3][3];
                double tmpVol062L = pow(volL,0.66666666);
                for(int j=0; j<3; ++j){
                    for(int k=0; k<3; ++k){
                        avgMij[j][k] = -tmpVol062L*SabsL*SijL[j][k]*volL;
                        avgSij[j][k] = SijL[j][k]*volL;
                    }
                }
                sumRho += rhoL * volL;
                sumV += uL * volL;
                sumU += vL * volL;
                sumW += wL * volL;
                sumSabs += SabsL * volL;
                sumVol += volL;
            }
            {
                // calculate strain rate tensor, Sij
                double SijR[3][3];
                SijR[0][0] = dudxR;
                SijR[0][1] = 0.5*(dudyR+dvdxR);
                SijR[0][2] = 0.5*(dudzR+dwdxR);
                SijR[1][0] = SijR[0][1];
                SijR[1][1] = dvdyR;
                SijR[1][2] = 0.5*(dvdzR+dwdyR);
                SijR[2][0] = SijR[0][2];
                SijR[2][1] = SijR[1][2];
                SijR[2][2] = dwdzR;

                // calculate S_ij * S_ij
                double SSR = SijR[0][0]*SijR[0][0] + SijR[1][1]*SijR[1][1] + SijR[2][2]*SijR[2][2] +
                    2.0*(SijR[0][1]*SijR[0][1] + SijR[0][2]*SijR[0][2] + SijR[1][2]*SijR[1][2]);
                
                SabsR = sqrt(2.0*SSR);
                
                // filtering process f(rho*ui*uj), f(rho), f(ui)
                // double LijR[3][3];
                avgLij[0][0] += rhoR * uR * uR * volR;
                avgLij[0][1] += rhoR * uR * vR * volR;
                avgLij[0][2] += rhoR * uR * wR * volR;
                avgLij[1][1] += rhoR * vR * vR * volR;
                avgLij[1][2] += rhoR * vR * wR * volR;
                avgLij[2][2] += rhoR * wR * wR * volR;
                    
                // double MijR[3][3];
                double tmpVol062R = pow(volR,0.66666666);
                for(int j=0; j<3; ++j){
                    for(int k=0; k<3; ++k){
                        avgMij[j][k] += (-tmpVol062R*SabsR*SijR[j][k]*volR);
                        avgSij[j][k] += SijR[j][k]*volR;
                    }
                }
                sumRho += rhoR * volR;
                sumV += uR * volR;
                sumU += vR * volR;
                sumW += wR * volR;
                sumSabs += SabsR * volR;
                sumVol += volR;
            }
                    

            for(int j=0; j<3; ++j){
                for(int k=0; k<3; ++k){
                    avgSij[j][k] /= sumVol;
                    avgLij[j][k] /= sumVol;
                    avgMij[j][k] /= sumVol;
                }
            }
            sumRho /= sumVol;
            sumV /= sumVol;
            sumU /= sumVol;
            sumW /= sumVol;
            sumSabs /= sumVol;
		
            
            avgLij[0][0] -= sumRho * sumU * sumU;
            avgLij[0][1] -= sumRho * sumU * sumV;
            avgLij[0][2] -= sumRho * sumU * sumW;
            avgLij[1][1] -= sumRho * sumV * sumV;
            avgLij[1][2] -= sumRho * sumV * sumW;
            avgLij[2][2] -= sumRho * sumW * sumW;
            avgLij[1][0] = avgLij[0][1];
            avgLij[2][0] = avgLij[0][2];
            avgLij[2][1] = avgLij[1][2];

            double sumDelta = pow(sumVol,0.66666666666);
            for(int j=0; j<3; ++j){
                for(int k=0; k<3; ++k){
                    avgMij[j][k] = sumDelta * sumSabs * avgSij[j][k] - avgMij[j][k];
                }
            }
            
            // ensemble average for numerator and denominator
            double sum1 = 0.0;
            double sum2 = 0.0;
            for(int j=0; j<3; ++j){
                for(int k=0; k<3; ++k){
                    sum1 += avgLij[j][k] * avgMij[j][k];
                    sum2 += avgMij[j][k] * avgMij[j][k];
                }
            }
            sum1 /= 9.0;
            sum2 /= 9.0;
            
            // calculate parameter C_R
            double C_R = 0.0;
            if(sum2 != 0.0){
                C_R = -0.5 * sum1 / sum2;
            }
            
            // calculate Delta_bar
            double DeltaBar2 = pow(0.5*(volL+volR),0.666666666);
            
            // calculate muT & kSGS
            double SabsF = 0.5*(SabsL+SabsR);
            double muT_F = C_R * 0.5*(rhoL+rhoR)*DeltaBar2*SabsF;
            double kSGS_F = C_I * DeltaBar2 * (0.5*SabsF*SabsF);
            
            // muT_F = max(-4.0*muF,min(2.0*muF,muT_F));
            muT_F = max(-muF,min(muF,muT_F));
            
            muF += muT_F;
            // =================================================================
            
            
            
            
            
            
            
            
            
            // =================================================================
			// AUSM-like expression of HLLC and its all-speed extension
            // Keiichi Kitamura
            double RT = sqrt(rhoFR/rhoFL); 
			double chat = 0.5*(cFL+cFR);
			double rhohat = sqrt(rhoFL*rhoFR);
			double Unhat = (UnL+RT*UnR)/(1.0+RT); 
			// double uhat = (uFL+RT*uFR)/(1.0+RT); 
			// double vhat = (vFL+RT*vFR)/(1.0+RT); 
			// double what = (wFL+RT*wFR)/(1.0+RT); 
			// double Yhat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yhat[i] = (YFL[i]+RT*YFR[i])/(1.0+RT);
            // }
            
            double preLs = abs(pL) + 0.1 * rhoFL*cFL*cFL;
            double preRs = abs(pR) + 0.1 * rhoFR*cFR*cFR;
            double fp = 0.0;
            fp = min(preLs/preRs,preRs/preLs);
            fp = fp*fp*fp;
            
            double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
            double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
            double KLR = sqrt(0.5*(U2L+U2R));
            double Mcy = min(1.0,KLR/chat);
            double phi_c = (1.0-Mcy)*(1.0-Mcy);
            double ML = UnL/chat;
            double MR = UnR/chat;
            
            double SL = min(0.0,min(UnL-cFL,Unhat-chat));
            double SR = max(0.0,max(UnR+cFR,Unhat+chat));
            
			double fluxB[nEq];
            { 
                // HLLCL
                double SL0 = min(UnL-cFL,UnR-cFR);
                double SR0 = max(UnR+cFR,UnL+cFL);
                // double theta = (rhoFL*(UnL-SL)+rhoFR*(SR-UnR))/ (rhoFL*(UnL-SL0)+rhoFR*(SR0-UnR));
                double theta = 1.0;
                double SM = (rhoFL*(UnL-SL)*UnL + rhoFR*(SR-UnR)*UnR - theta*(pR-pL))/(rhoFL*(UnL-SL)+rhoFR*(SR-UnR));
                double pStar = 0.5*(pL+pR+rhoFL*(UnL-SL)*(UnL-SM)+rhoFR*(UnR-SR)*(UnR-SM));
                
                double mdot = 0.0;
                double pF = pStar;
                if(SM>0.0){
                    mdot = rhoFL*(UnL + SL*((SL-UnL)/(SL-SM)-1.0));
                    
                    pF += SM/(SM-SL)*(pL-pStar);
                }
                else{
                    mdot = rhoFR*(UnR + SR*((SR-UnR)/(SR-SM)-1.0));
                    
                    pF += SM/(SM-SR)*(pR-pStar);
                }
                
                // // double absMhat = abs(Unhat)/chat;
                // // mdot += (fp-1.0)*SL*SR/(SR-SL)/(1.0+absMhat)*(pR-pL)/chat/chat;
                // // mdot -= (1.0-fp)*phi_c*SL*SR/(SR-SL)*(pR-pL)/chat/chat;
                // // mdot -= 0.5*(1.0-fp)*phi_c/chat*(pR-pL);
                // mdot -= 0.5*phi_c/chat*(pR-pL);
                
                // in pressure-based
                mdot -= dAlpha * dt/dLR*(pR-pL);
                mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
                mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
                mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
                
                double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
                if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
                double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
                if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
                // pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
                pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
                
                if(mdot>=0.0){
                    fluxB[0] = -( mdot )*area;
                    fluxB[1] = -( mdot*uFL + pF*nvec[0] )*area;
                    fluxB[2] = -( mdot*vFL + pF*nvec[1] )*area;
                    fluxB[3] = -( mdot*wFL + pF*nvec[2] )*area;
                    fluxB[4] = -( mdot*HtFL  + SL*(pStar-pL)/(SL-UnL) )*area;
                    // fluxB[4] = -( mdot*HtFL )*area;
                    for(int i=0; i<nSp-1; ++i){
                        fluxB[5+i] = -( mdot*YFL[i] )*area;
                    }
                }
                else{
                    fluxB[0] = -( mdot )*area;
                    fluxB[1] = -( mdot*uFR + pF*nvec[0] )*area;
                    fluxB[2] = -( mdot*vFR + pF*nvec[1] )*area;
                    fluxB[3] = -( mdot*wFR + pF*nvec[2] )*area;
                    fluxB[4] = -( mdot*HtFR  + SR*(pStar-pR)/(SR-UnR) )*area;
                    // fluxB[4] = -( mdot*HtFR )*area;
                    for(int i=0; i<nSp-1; ++i){
                        fluxB[5+i] = -( mdot*YFR[i] )*area;
                    }
                }
                
                
            }
            // =================================================================
            
            
            
            
            
				
            
            
            
            
            // // =================================================================
			// // SLAU2
			// double chat = 0.5*(cFL+cFR);
			// double rhohat = 0.5*(rhoFL+rhoFR);
            // double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
            // double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
			// double Ubar = ( rhoFL*abs(UnL)+rhoFR*abs(UnR) ) / ( rhoFL + rhoFR );
			// // double chat = 0.5*(cL+cR);
			// // double rhohat = 0.5*(rhoL+rhoR);
            // // double U2L = uL*uL+vL*vL+wL*wL;
            // // double U2R = uR*uR+vR*vR+wR*wR;
			// // double Ubar = ( rhoL*abs(UnL)+rhoR*abs(UnR) ) / ( rhoL + rhoR );
            
            // double Unhat = 0.5*(UnL+UnR);   
			// double KLR = sqrt(0.5*(U2L+U2R));
			// double Mcy = min(1.0,KLR/chat);
            
			// double ML = UnL/chat;
			// double MR = UnR/chat;
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// double g_c = -max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
            // double mdot = 0.5*(
                // rhoFL*(UnL+(1.0-g_c)*Ubar+g_c*abs(UnL)) +
                // rhoFR*(UnR-(1.0-g_c)*Ubar-g_c*abs(UnR)) -
                // phi_c/chat*(pR-pL));
                
			// // double MLP = 0.5*(ML+abs(ML));
			// // if( abs(ML) < 1.0 ) MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			// // double MRM = 0.5*(MR-abs(MR));
            // // if(MLP+MRM>=0.0){
                // // mdot = rhoFL*(MLP+MRM)*chat - 0.5*phi_c/chat*(pR-pL);
            // // }
            // // else{
                // // mdot = rhoFR*(MLP+MRM)*chat - 0.5*phi_c/chat*(pR-pL);
            // // }
                
			// // in pressure-based
			// mdot -= dAlpha * dt/dLR*(pR-pL);
			// mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
			// mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
			// mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
            
			// double PLP = ( ML>0.0 ? 1.0 : 0.0 );
			// double PRM = ( MR<0.0 ? 1.0 : 0.0 );
			// if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            
			// // // K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
			// // double UnF = 0.5*(mdot+abs(mdot))/rhoL + 0.5*(mdot-abs(mdot))/rhoR;
            // // mdot = 0.5*(UnF+abs(UnF))*rhoFL + 0.5*(UnF-abs(UnF))*rhoFR;
			
			// // // in pressure-based
			// // UnF -= dAlpha * dt*0.5*(1.0/rhoL+1.0/rhoR)/dLR*(pR-pL);
			// // UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
			// // UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
			// // UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
            
			// // double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
			// double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
            
            // double f1L = 0.5*(mdot+abs(mdot));
            // double f1R = 0.5*(mdot-abs(mdot));
			
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
            // // =================================================================
            
            
            








            // // =================================================================
            // // anti-diffusion
            // double anti_diffusion_coeff = 0.01;
            
            // for(int i=0; i<nSp-1; ++i){
                
                // double L_flux;
                // double R_flux;
                // double nUnL = 0.0;
                // double nUnR = 0.0;
                // {
                    // double tmp_gradx = (cellsL[id_dYdx[i]]);
                    // double tmp_grady = (cellsL[id_dYdy[i]]);
                    // double tmp_gradz = (cellsL[id_dYdz[i]]);

                    // double magGrad = tmp_gradx*tmp_gradx+tmp_grady*tmp_grady+tmp_gradz*tmp_gradz;
                    // magGrad = sqrt(magGrad);
                    
                    // double normalVec[3];
                    // normalVec[0] = 0.0; normalVec[1] = 0.0; normalVec[2] = 0.0;
                    // if(magGrad!=0.0){
                        // normalVec[0] = tmp_gradx/magGrad; normalVec[1] = tmp_grady/magGrad; normalVec[2] = tmp_gradz/magGrad;
                    // }
                    
                    // double absU = anti_diffusion_coeff*sqrt(uL*uL+vL*vL+wL*wL);
                    
                    // nUnL += absU*normalVec[0]*nvec[0];
                    // nUnL += absU*normalVec[1]*nvec[1];
                    // nUnL += absU*normalVec[2]*nvec[2];
                
                    // L_flux = ( rhoL*nUnL*YL[i]*(1.0-YL[i]) )*area;
                // }
                // {
                    // double tmp_gradx = (cellsR[id_dYdx[i]]);
                    // double tmp_grady = (cellsR[id_dYdy[i]]);
                    // double tmp_gradz = (cellsR[id_dYdz[i]]);

                    // double magGrad = tmp_gradx*tmp_gradx+tmp_grady*tmp_grady+tmp_gradz*tmp_gradz;
                    // magGrad = sqrt(magGrad);
                    
                    // double normalVec[3];
                    // normalVec[0] = 0.0; normalVec[1] = 0.0; normalVec[2] = 0.0;
                    // if(magGrad!=0.0){
                        // normalVec[0] = tmp_gradx/magGrad; normalVec[1] = tmp_grady/magGrad; normalVec[2] = tmp_gradz/magGrad;
                    // }
                    
                    // double absU = anti_diffusion_coeff*sqrt(uR*uR+vR*vR+wR*wR);
                    
                    // nUnR += absU*normalVec[0]*nvec[0];
                    // nUnR += absU*normalVec[1]*nvec[1];
                    // nUnR += absU*normalVec[2]*nvec[2];
                    
                    // R_flux = ( rhoR*nUnR*YR[i]*(1.0-YR[i]) )*area;
                // }
                
                // // fluxB[5+i] -= 0.5*(L_flux+R_flux);
                
                // double nUn = 0.5*(nUnL+nUnR);
                // if(nUn>=0.0){fluxB[5+i] -= ( nUn*rhoL*YL[i]*(1.0-YL[i]) )*area;}
                // else{fluxB[5+i] -= ( nUn*rhoR*YR[i]*(1.0-YR[i]) )*area;}
                
                
            // }
            // // =================================================================
            
            
            
            
            
			
			// int iter = 0;
			
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
                
                // MY Oomar, 2021
                double st_SL = min(UnL-chat,0.5*(UnL+UnR)-0.5*(cL+cR));
                double st_SR = max(UnL+chat,0.5*(UnL+UnR)+0.5*(cL+cR));
                double st_SM = (pL-pR-surf_sigma[i]*0.5*(curvatureL[i]+curvatureR[i])*(alpha_VFL[i]-alpha_VFR[i]))/
                                (rhoR*st_SR-rhoL*st_SL);
                double st_weiL = (st_SM>=0.0 ? 1.0 : 0.0); double st_weiR = 1.0-st_weiL;
                alpha_VFF = st_weiL*alpha_VFL[i] + st_weiR*alpha_VFR[i];
                
                
				fluxB_LL[1] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[0] )*area;
				fluxB_LL[2] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[1] )*area;
				fluxB_LL[3] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*nvec[2] )*area;
				fluxB_LL[4] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*UnL )*area;
				
				fluxB_RR[1] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[0] )*area;
				fluxB_RR[2] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[1] )*area;
				fluxB_RR[3] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*nvec[2] )*area;
				fluxB_RR[4] -= surf_sigma[i]*curvatureR[i]*( alpha_VFF*UnR )*area;
				
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
            id_dYdx,id_dYdy,id_dYdz](
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
				
				// double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uF*uF+vF*vF+wF*wF));
				// double cbar = cL;//0.5*(cL+cF);
				// double Mdash = min(1.0,KLR/cbar);
				// double chi = (1.0-Mdash)*(1.0-Mdash);
				// double rhohat = 0.5*rhoL+0.5*rhoF;
				// double dp_coeff_thm = 0.5*chi/rhohat/cbar;
				
				
				
				fluxB_LL[0] -= ( rhoF*UnF )*area;
				fluxB_LL[1] -= ( rhoF*UnF*uF + pF*nvec[0] )*area;
				fluxB_LL[2] -= ( rhoF*UnF*vF + pF*nvec[1] )*area;
				fluxB_LL[3] -= ( rhoF*UnF*wF + pF*nvec[2] )*area;
				fluxB_LL[4] -= ( rhoF*UnF*HtF )*area;
				for(int i=0; i<nSp-1; ++i){
					fluxB_LL[5+i] -= ( rhoF*UnF*YF[i] )*area;
				}
				
					
                    
                    
                    
                // // anti-diffusion
                // double anti_diffusion_coeff = 0.01;
                // for(int i=0; i<nSp-1; ++i){
                    // double magVelL = sqrt(uL*uL+vL*vL+wL*wL);
                    // double magVelF = (magVelL);
                    
                    // double tmp_gradx = (cellsL[id_dYdx[i]]);
                    // double tmp_grady = (cellsL[id_dYdy[i]]);
                    // double tmp_gradz = (cellsL[id_dYdz[i]]);

                    // double magGrad = tmp_gradx*tmp_gradx+tmp_grady*tmp_grady+tmp_gradz*tmp_gradz;
                    // magGrad = sqrt(magGrad);
                    
                    // double normalVec[3];
                    // normalVec[0] = 0.0; normalVec[1] = 0.0; normalVec[2] = 0.0;
                    // if(magGrad!=0.0){
                        // normalVec[0] = tmp_gradx/magGrad; normalVec[1] = tmp_grady/magGrad; normalVec[2] = tmp_gradz/magGrad;
                    // }
                
                    // double absU = anti_diffusion_coeff*abs(UnF);
                    // double normS = 0.0;
                    // normS += nvec[0]*normalVec[0];
                    // normS += nvec[1]*normalVec[1];
                    // normS += nvec[2]*normalVec[2];
                        
                    // double tmp_YF = (YF[i]);
                    // fluxB_LL[5+i] -= ( rhoF*absU*normS*tmp_YF*(1.0-tmp_YF) )*area;
                    // // fluxB_LL[5+i] += ( rhoF*absU*normS*tmp_YF*(1.0-tmp_YF) )*area;
                // }
            
				
				
				
				
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
					fluxB_LL[4] += surf_sigma[i]*curvatureL[i]*( alpha_VFF*UnF )*area;
					
				}
				
				
				return 0;
			});
		}
		
	}
		
	
}



                
                
				
				
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
			
			
			// double chat_YYL = 0.5*(cL+cR);
			// // double chat_YYL = sqrt(cL*cR);
			// // double chat_YYL = 1.0/(0.5*(1.0/cL+1.0/cR));
			// // double rhohat_YYL = sqrt(rhoL*rhoR);
			// double rhohat_YYL = 0.5*(rhoL+rhoR);
			// double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
			// double Mcy = min(1.0,KLR/chat_YYL);
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// double UnL_YYL = UnL;
			// double UnR_YYL = UnR;
			
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
			
			// double ML_YYL = UnL_YYL/chat_YYL;
			// double MR_YYL = UnR_YYL/chat_YYL;
			// double g_YYL = 1.0 + max( min(ML_YYL,0.0), -1.0 )*min( max(MR_YYL,0.0), 1.0 );
			// double D_L_YYL = UnL_YYL+(1.0-g_YYL)*abs(UnL_YYL);
			// double D_R_YYL = UnR_YYL-(1.0-g_YYL)*abs(UnR_YYL);
			// double D_rho_YYL = g_YYL*( rhoL*abs(UnL_YYL)+rhoR*abs(UnR_YYL) ) / ( rhoL + rhoR );
			// double MLP = 0.5*(ML_YYL+abs(ML_YYL));
			// if( abs(ML_YYL) < 1.0 ) {
				// MLP = 0.25*(ML_YYL + 1.0)*(ML_YYL + 1.0);
				// // MLP += 0.125*(ML_YYL*ML_YYL-1.0)*(ML_YYL*ML_YYL-1.0);
			// }
			// double MRM = 0.5*(MR_YYL-abs(MR_YYL));
			// if( abs(MR_YYL) < 1.0 ) {
				// MRM = -0.25*(MR_YYL - 1.0)*(MR_YYL - 1.0);
				// // MRM -= 0.125*(MR_YYL*MR_YYL-1.0)*(MR_YYL*MR_YYL-1.0);
			// }
			// double UnF = wdL*(D_L_YYL+D_rho_YYL) + wdR*(D_R_YYL-D_rho_YYL);
			
			// double weightThrmL = 1.0/(1.0+exp(-ML_YYL*3.5));
			// double weightThrmR = 1.0/(1.0+exp(MR_YYL*3.5));
			// double total_weightThrm = weightThrmL + weightThrmR;
			// weightThrmL /= total_weightThrm;
			// weightThrmR /= total_weightThrm;
			
			// UnF = weightThrmL*UnL + weightThrmR*UnR;
			// UnF = 0.5*UnL + 0.5*UnR;
			// UnF = MLP*UnL + MRM*UnR;
			
			
			// // 추가적 보정텀
			// double dp_coeff_thm = 0.5*phi_c/chat_YYL/rhohat_YYL;
			// double dp_coeff_Un = dAlpha * dt*(wdL/rhoL+wdR/rhoR)/dLR;
			// UnF -= dp_coeff_thm*(pR-pL);
			// UnF -= dp_coeff_Un*(pR-pL);
			// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*dAlpha*nLR[0];
			// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*dAlpha*nLR[1];
			// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*dAlpha*nLR[2];
			
			// double WUL = 0.5*(1.0+(1.0-g_YYL)*(UnL>0.0 ? 1.0 : -1.0) + 
				// g_YYL*( rhoL*(UnL>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
			// double WUR = 0.5*(1.0-(1.0-g_YYL)*(UnR>0.0 ? 1.0 : -1.0) - 
				// g_YYL*( rhoR*(UnR>0.0 ? 1.0 : -1.0) ) / ( rhoL + rhoR ));
				
				
			// WUL = weightThrmL;
			// WUR = weightThrmR;
				
			
			// double PLP = 0.5*(1.0 + ( ML_YYL>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML_YYL) < 1.0 ) {
				// PLP = 0.25*(ML_YYL+1.0)*(ML_YYL+1.0)*(2.0-ML_YYL);
				// // PLP += 0.1875*ML_YYL*(ML_YYL*ML_YYL-1.0)*(ML_YYL*ML_YYL-1.0);
			// } 
			// double PRM = 0.5*(1.0 - ( MR_YYL>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR_YYL) < 1.0 ) {
				// PRM = 0.25*(MR_YYL-1.0)*(MR_YYL-1.0)*(2.0+MR_YYL);
				// // PRM -= 0.1875*MR_YYL*(MR_YYL*MR_YYL-1.0)*(MR_YYL*MR_YYL-1.0);
			// } 
			// // double WpL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
			// // double WpR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
			// double WpL = wdL + 0.5*(PLP-PRM);
			// double WpR = wdR - 0.5*(PLP-PRM);
			// double pL_YYL = pL;
			// double pR_YYL = pR;
			
			// // // high-order
			// // pL_YYL += gradLim_pL*(dpdxL*xyzLF[0] + dpdyL*xyzLF[1] + dpdzL*xyzLF[2]);
			// // pR_YYL += gradLim_pR*(dpdxR*xyzRF[0] + dpdyR*xyzRF[1] + dpdzR*xyzRF[2]);
			// // skewness
			// pL_YYL += gradLim_pL*(dpdxL*LNv[0] + dpdyL*LNv[1] + dpdzL*LNv[2]);
			// pR_YYL += gradLim_pR*(dpdxR*RNv[0] + dpdyR*RNv[1] + dpdzR*RNv[2]);
			
			// PLP = weightThrmL;
			// PRM = weightThrmR;
			
			// double pF = PLP*pL_YYL+PRM*pR_YYL;
			// double pF = (0.5*pL_YYL+0.5*pR_YYL) + 0.5*(PLP-PRM)*(pL_YYL-pR_YYL) + 
						// (PLP+PRM-1.0)*KLR*(wdL*rhoL+wdR*rhoR)*(wdL*cL+wdR*cR);
						// Mcy*(PLP+PRM-1.0)*(pL_YYL+pR_YYL) - 
						// 0.0;//0.5*(rhohat_YYL*(chat_YYL-abs(UnF))*(UnR-UnL));
			
			// double SL = min(UnL-cL,UnF-chat_YYL);
			// double SR = max(UnR+cR,UnF+chat_YYL);
			// double SL = min(UnL-cL,0.5*(UnL+UnR)-0.5*(cL+cR));
			// double SR = max(UnR+cR,0.5*(UnL+UnR)+0.5*(cL+cR));
			// double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						// (rhoR*(SR-UnR)-rhoL*(SL-UnL));
			// double pStar = rhoL*(UnL-SL)*(UnL-SM)+pL; // rhoR*(UnR-SR)*(UnR-SM)+pR;
			// if(SL<=0.0 && 0.0<SM){ 
				// pF = pStar;
				// pF += SM/(SL-SM)*(pStar-pL)*SM;
			// }
			// else if(SM<=0.0 && 0.0<=SR){ 
				// pF = pStar; 
				// pF += SM/(SR-SM)*(pStar-pR)*SM;
			// }
			// else if(SL>0.0){ pF = pL; }
			// else if(SR<0.0){ pF = pR; }
			
			
			
			
			// double Un_roe = (sqrt(rhoL)*UnL+sqrt(rhoR)*UnR)/(sqrt(rhoL)+sqrt(rhoR));
			// double c_roe = (sqrt(rhoL)*cL+sqrt(rhoR)*cR)/(sqrt(rhoL)+sqrt(rhoR));
			// double Un_roe = 0.5*(UnL+UnR);
			// double c_roe = 0.5*(cL+cR);
			// dp_coeff_Un = dAlpha * dt*(wdL/rhoL+wdR/rhoR)/dLR;
			// dp_coeff_thm = 0.5*phi_c/c_roe/(0.5*(rhoL+rhoR));
			// // diffDPU = 0.0;
			// WUL = MLP;
			// WUR = MRM;
			
			// double chat = 0.5*(cL+cR);
			// double rhohat = 0.5*(rhoL+rhoR);
			// double Unhat = 0.5*(UnL+UnR);
			// double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
			// double Mcy = min(1.0,KLR/chat);
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);
			
			// double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
			// double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
			// double weiSup = 1.0 - pow( min(preLs/preRs,preRs/preLs),2.0);
			
			// double ML = UnL/chat;
			// double MR = UnR/chat;
			// double g = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			// double D_L = UnL+(1.0-g)*abs(UnL);
			// double D_R = UnR-(1.0-g)*abs(UnR);
			// double D_rho = g*( rhoL*abs(UnL)+rhoR*abs(UnR) ) / ( rhoL + rhoR );
			
			// double UnF = ( 0.5*(D_L+D_rho) + 0.5*(D_R-D_rho) );
			// double dp_coeff_thm = 0.5*phi_c/chat/rhohat;
			// UnF -= dp_coeff_thm*(pR-pL);
					
			// double PLP = 0.5*(1.0 +1.0 * ( ML>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML) < 1.0 ) {
				// PLP = 0.25*pow((ML+1.0),2.0)*(2.0-ML);
			// }
			// double PRM = 0.5*(1.0 -1.0 * ( MR>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR) < 1.0 ) {
				// PRM = 0.25*pow((MR-1.0),2.0)*(2.0+MR);
			// }

			// double pF = PLP*pL+PRM*pR;
			// double UnF = ( MLP + MRM )*chat_YYL;
			// UnF = weiSup*( MLP + MRM )*chat_YYL;
			// UnF += (1.0-weiSup)*( 0.5*(D_L_YYL+D_rho_YYL) + 0.5*(D_R_YYL-D_rho_YYL) );
			
			// UnF = phi_c*0.5*(UnL+UnR) + (1.0-phi_c)*UnF;
			
			// UnF -= dp_coeff_thm*(pR-pL);
			// UnF -= dp_coeff_Un*(pR-pL);
			// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*dAlpha*nLR[0];
			// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*dAlpha*nLR[1];
			// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*dAlpha*nLR[2];
			
			// double pF = (PLP*pL+PRM*pR) + 0.5*(PLP-PRM)*(pL-pR) + 
				// (PLP+PRM-1.0)*KLR*rhohat*chat;
				
				
				
				
				
			// // HLLC-AP
			// double chat = 0.5*(cL+cR);
			// double rhohat = 0.5*(rhoL+rhoR);
			// double Unhat = 0.5*(UnL+UnR);
			// double SL = min(0.0,min(UnL-cL,Unhat-chat));
			// double SR = max(0.0,max(UnR+cR,Unhat+chat));
			// double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						// (rhoR*(SR-UnR)-rhoL*(SL-UnL));
			// double pStar = 0.5*(pL+pR+rhoL*(UnL-SL)*(UnL-SM)+rhoR*(UnR-SR)*(UnR-SM));
            
            // double U2L = uL*uL+vL*vL+wL*wL;
            // double U2R = uR*uR+vR*vR+wR*wR;
            // double Mbar0 = max(sqrt(U2L)/chat,sqrt(U2R)/chat);
            // double kL=0.0; double kR=0.0;
            // if( (UnL-cL)>0.0 && (UnR-cR)<0.0 ) kL=1.0;
            // if( (UnL+cL)>0.0 && (UnR+cR)<0.0 ) kL=1.0;
            // if( (UnR-cR)>0.0 && (UnL-cL)<0.0 ) kR=1.0;
            // if( (UnR+cR)>0.0 && (UnL+cL)<0.0 ) kR=1.0;
            // double fM = 1.0;
            // if(kL<1.e-200 && kR<1.e-200) fM=min(1.0,Mbar0);
            // double apc_term = rhoL*rhoR*(SL-UnL)*(SR-UnR)/(rhoR*(SR-UnR)-rhoL*(SL-UnL));
            
			// double fluxB[nEq];
            // if(0.0<=SL){
                // fluxB[0] = -( rhoFL*UnL )*area;
                // fluxB[1] = -( rhoFL*UnL*uFL + pL*nvec[0] )*area;
                // fluxB[2] = -( rhoFL*UnL*vFL + pL*nvec[1] )*area;
                // fluxB[3] = -( rhoFL*UnL*wFL + pL*nvec[2] )*area;
                // fluxB[4] = -( rhoFL*UnL*HtFL )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( rhoFL*UnL*YFL[i] )*area;
                // }
            // }
            // else if(SL<0.0 && 0.0<=SM){
                // double Coef0 = SM/(SL-SM)*SL;
                // double Coef1 = SM/(SL-SM);
                // double Coef2 = SL/(SL-SM);
                // fluxB[0] = -( Coef0*rhoFL - Coef1*rhoFL*UnL )*area;
                // fluxB[1] = -( Coef0*rhoFL*uFL - Coef1*(rhoFL*uFL*UnL+pL*nvec[0]) + Coef2*pStar*nvec[0] -
                              // (1.0-fM)*( apc_term*nvec[0] ))*area;
                // fluxB[2] = -( Coef0*rhoFL*vFL - Coef1*(rhoFL*vFL*UnL+pL*nvec[1]) + Coef2*pStar*nvec[1] -
                              // (1.0-fM)*( apc_term*nvec[1] ))*area;
                // fluxB[3] = -( Coef0*rhoFL*wFL - Coef1*(rhoFL*wFL*UnL+pL*nvec[2]) + Coef2*pStar*nvec[2] -
                              // (1.0-fM)*( apc_term*nvec[2] ))*area;
                // fluxB[4] = -( Coef0*(rhoFL*HtFL-pL) - Coef1*(rhoFL*HtFL*UnL) + Coef2*pStar*SM -
                              // (1.0-fM)*( apc_term*SM ) )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( Coef0*(rhoFL*YFL[i]) - Coef1*(rhoFL*YFL[i]*UnL) )*area;
                // }
            // }
            // else if(SM<0.0 && 0.0<=SR){
                // double Coef0 = SM/(SR-SM)*SR;
                // double Coef1 = SM/(SR-SM);
                // double Coef2 = SR/(SR-SM);
                // fluxB[0] = -( Coef0*rhoFR - Coef1*rhoFR*UnR )*area;
                // fluxB[1] = -( Coef0*rhoFR*uFR - Coef1*(rhoFR*uFR*UnR+pR*nvec[0]) + Coef2*pStar*nvec[0] -
                              // (1.0-fM)*( apc_term*nvec[0] ) )*area;
                // fluxB[2] = -( Coef0*rhoFR*vFR - Coef1*(rhoFR*vFR*UnR+pR*nvec[1]) + Coef2*pStar*nvec[1] -
                              // (1.0-fM)*( apc_term*nvec[1] ) )*area;
                // fluxB[3] = -( Coef0*rhoFR*wFR - Coef1*(rhoFR*wFR*UnR+pR*nvec[2]) + Coef2*pStar*nvec[2] -
                              // (1.0-fM)*( apc_term*nvec[2] ) )*area;
                // fluxB[4] = -( Coef0*(rhoFR*HtFR-pR) - Coef1*(rhoFR*HtFR*UnR) + Coef2*pStar*SM -
                              // (1.0-fM)*( apc_term*SM ) )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( Coef0*(rhoFR*YFR[i]) - Coef1*(rhoFR*YFR[i]*UnR) )*area;
                // }
            // }
            // else{
                // fluxB[0] = -( rhoFR*UnR )*area;
                // fluxB[1] = -( rhoFR*UnR*uFR + pR*nvec[0] )*area;
                // fluxB[2] = -( rhoFR*UnR*vFR + pR*nvec[1] )*area;
                // fluxB[3] = -( rhoFR*UnR*wFR + pR*nvec[2] )*area;
                // fluxB[4] = -( rhoFR*UnR*HtFR )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( rhoFR*UnR*YFR[i] )*area;
                // }
            // }
            
            
            
            
            
            
            
            
            
            
            
            
            
            
			// // HLLC - AP _ YYL
            // double RT = sqrt(rhoFR/rhoFL); 
			// double chat = 0.5*(cFL+cFR);
			// // double chat = sqrt(cL*cR);
			// double rhohat = sqrt(rhoFL*rhoFR);//(rhoL+RT*rhoR)/(1.0+RT); 
            
            // // double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
            // // double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
            // // double fp = 0.0;
            // // fp = min(preLs/preRs,preRs/preLs);
            // // fp = fp*fp*fp;
            
            // // double zF = min(max(sqrt(uL*uL+vL*vL+wL*wL)/cL,sqrt(uR*uR+vR*vR+wR*wR)/cR),1.0);
            // // double UnLs = 0.5*( (1.0+zF)*UnL + (1.0-zF)*UnR );
            // // double UnRs = 0.5*( (1.0+zF)*UnR + (1.0-zF)*UnL );
            // // double UnLas = fp*UnLs + (1.0-fp)*UnL;
            // // double UnRas = fp*UnRs + (1.0-fp)*UnR;
            
			// double Unhat = (UnL+RT*UnR)/(1.0+RT); 
			// double uhat = (uFL+RT*uFR)/(1.0+RT); 
			// double vhat = (vFL+RT*vFR)/(1.0+RT); 
			// double what = (wFL+RT*wFR)/(1.0+RT); 
			// double Yhat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yhat[i] = (YFL[i]+RT*YFR[i])/(1.0+RT); 
                // // Yhat[i] = 0.5*(YFL[i]+YFR[i]);
            // }
			// // double SL = min(UnL-cL,UnR-cR);
			// // double SR = max(UnR+cR,UnR+cR);
            // // double SL = min(UnL-cL,Unhat-chat);
            // // double SR = max(UnR+cR,Unhat+chat);
            // double SL = min(0.0,min(UnL-cFL,Unhat-chat));
            // double SR = max(0.0,max(UnR+cFR,Unhat+chat));
            
			// double fluxB[nEq];
			// // 컨벡티브 B
            // if(SL>=0.0){
                // fluxB[0] = -( rhoFL*UnL )*area;
                // fluxB[1] = -( rhoFL*UnL*uFL + pL*nvec[0] )*area;
                // fluxB[2] = -( rhoFL*UnL*vFL + pL*nvec[1] )*area;
                // fluxB[3] = -( rhoFL*UnL*wFL + pL*nvec[2] )*area;
                // fluxB[4] = -( rhoFL*UnL*HtFL )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( rhoFL*UnL*YFL[i] )*area;
                // }
            // }
            // else if(0.0>=SR){
                // fluxB[0] = -( rhoFR*UnR )*area;
                // fluxB[1] = -( rhoFR*UnR*uFR + pR*nvec[0] )*area;
                // fluxB[2] = -( rhoFR*UnR*vFR + pR*nvec[1] )*area;
                // fluxB[3] = -( rhoFR*UnR*wFR + pR*nvec[2] )*area;
                // fluxB[4] = -( rhoFR*UnR*HtFR )*area;
                // for(int i=0; i<nSp-1; ++i){
                    // fluxB[5+i] = -( rhoFR*UnR*YFR[i] )*area;
                // }
            // }
            // else{
                // // UnL = UnLas;
                // // UnR = UnRas;

                // // HLLC
                // // double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
                            // // (rhoR*(SR-UnR)-rhoL*(SL-UnL));
                // // // double pStar = 0.5*(pL+pR+rhoL*(UnL-SL)*(UnL-SM)+rhoR*(UnR-SR)*(UnR-SM));
                // // double pStar = rhoL*(UnL-SL)*(UnL-SM)+pL; // rhoR*(UnR-SR)*(UnR-SM)+pR;
                // // if(SL<=0.0 && 0.0<SM){ 
                    // // double UnStar = (SL-UnL)/(SL-SM)*SM;
                    // // double dpStar = (pStar-pL)/(SL-SM)*SM;
                    // // double pF = pStar;
                    
                    // // fluxB[0] = -( rhoL*UnStar )*area;
                    // // fluxB[1] = -( rhoL*uL*UnStar + dpStar*nvec[0] + pF*nvec[0])*area;
                    // // fluxB[2] = -( rhoL*vL*UnStar + dpStar*nvec[1] + pF*nvec[1])*area;
                    // // fluxB[3] = -( rhoL*wL*UnStar + dpStar*nvec[2] + pF*nvec[2])*area;
                    // // fluxB[4] = -( rhoL*HtL*UnStar )*area;
                    // // for(int i=0; i<nSp-1; ++i){
                        // // fluxB[5+i] = -( rhoL*YL[i]*UnStar )*area;
                    // // }
                    
                // // }
                // // else { 
                    // // double UnStar = (SR-UnR)/(SR-SM)*SM;
                    // // double dpStar = (pStar-pR)/(SR-SM)*SM;
                    // // double pF = pStar;
                    
                    // // fluxB[0] = -( rhoR*UnStar )*area;
                    // // fluxB[1] = -( rhoR*uR*UnStar + dpStar*nvec[0] + pF*nvec[0])*area;
                    // // fluxB[2] = -( rhoR*vR*UnStar + dpStar*nvec[1] + pF*nvec[1])*area;
                    // // fluxB[3] = -( rhoR*wR*UnStar + dpStar*nvec[2] + pF*nvec[2])*area;
                    // // fluxB[4] = -( rhoR*HtR*UnStar )*area;
                    // // for(int i=0; i<nSp-1; ++i){
                        // // fluxB[5+i] = -( rhoR*YR[i]*UnStar )*area;
                    // // }
                // // }
                
                // // HLLC
                // double SM = (rhoFR*UnR*(SR-UnR)-rhoFL*UnL*(SL-UnL)+pL-pR)/
                            // (rhoFR*(SR-UnR)-rhoFL*(SL-UnL));
                // double pStar = 0.5*(pL+pR+rhoFL*(UnL-SL)*(UnL-SM)+rhoFR*(UnR-SR)*(UnR-SM));
                
                // double mdot = 0.0;
                // double pF = pStar;
                // if(SM>0.0){
                    // // double rhoLStar = (SL-UnL)/(SL-SM)*rhoL;
                    // // mdot = rhoL*UnL + SL*(rhoLStar-rhoL);
                    
                    // double rhoLStar = (SL-UnL)/(SL-SM)*rhoFL;
                    // mdot = rhoFL*UnL + SL*(rhoLStar-rhoFL);
                    
                    // pF += SM/(SM-SL)*(pL-pStar);
                // }
                // else{
                    // // double rhoRStar = (SR-UnR)/(SR-SM)*rhoR;
                    // // mdot = rhoR*UnR + SR*(rhoRStar-rhoR);
                    
                    // double rhoRStar = (SR-UnR)/(SR-SM)*rhoFR;
                    // mdot = rhoFR*UnR + SR*(rhoRStar-rhoFR);
                    
                    // pF += SM/(SM-SR)*(pR-pStar);
                // }
                
                
                // // double SM = (rhoFR*UnR*(SR-UnR)-rhoFL*UnL*(SL-UnL)+pL-pR)/
                            // // (rhoFR*(SR-UnR)-rhoFL*(SL-UnL));
                // // double pStar = 0.5*(pL+pR+rhoFL*(UnL-SL)*(UnL-SM)+rhoFR*(UnR-SR)*(UnR-SM));
                
                // // double mdot = 0.0;
                // // double pF = pStar;
                // // if(SM>0.0){
                    // // double rhoLStar = (SL-UnL)/(SL-SM)*rhoFL;
                    // // mdot = rhoFL*UnL + SL*(rhoLStar-rhoFL);
                    // // pF += SM/(SM-SL)*(pL-pStar);
                // // }
                // // else{
                    // // double rhoRStar = (SR-UnR)/(SR-SM)*rhoFR;
                    // // mdot = rhoFR*UnR + SR*(rhoRStar-rhoFR);
                    // // pF += SM/(SM-SR)*(pR-pStar);
                // // }
                
                   
                // double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
                // double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
                // double KLR = sqrt(0.5*(U2L+U2R));
                // double Mcy = min(1.0,KLR/chat);
                // double ML = UnL/chat;
                // double MR = UnR/chat;
                // double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
                // if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
                // double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
                // if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
                // // pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
                // pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
                
                
                // // // in pressure-based
                // // mdot -= dAlpha * dt/dLR*(pR-pL);
                // // mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
                // // mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
                // // mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
                
                
                // if(mdot>=0.0){
                    
                    // // mdot = mdot/rhoL*rhoFL;
                    
                    // fluxB[0] = -( mdot )*area;
                    // fluxB[1] = -( mdot*uFL + pF*nvec[0] )*area;
                    // fluxB[2] = -( mdot*vFL + pF*nvec[1] )*area;
                    // fluxB[3] = -( mdot*wFL + pF*nvec[2] )*area;
                    // fluxB[4] = -( mdot*HtFL  + SL*(pStar-pL)/(SL-UnL) )*area;
                    // // fluxB[4] = -( mdot*HtFL )*area;
                    // for(int i=0; i<nSp-1; ++i){
                        // fluxB[5+i] = -( mdot*YFL[i] )*area;
                    // }
                // }
                // else{
                    
                    // // mdot = mdot/rhoR*rhoFR;
                    
                    // fluxB[0] = -( mdot )*area;
                    // fluxB[1] = -( mdot*uFR + pF*nvec[0] )*area;
                    // fluxB[2] = -( mdot*vFR + pF*nvec[1] )*area;
                    // fluxB[3] = -( mdot*wFR + pF*nvec[2] )*area;
                    // fluxB[4] = -( mdot*HtFR  + SR*(pStar-pR)/(SR-UnR) )*area;
                    // // fluxB[4] = -( mdot*HtFR )*area;
                    // for(int i=0; i<nSp-1; ++i){
                        // fluxB[5+i] = -( mdot*YFR[i] )*area;
                    // }
                // }
                
                
            // }
            
            
            
            
            
            
			
			// // HLLC, AUSM version -> kitamura, shima, 2018
			// double chat = 0.5*(cL+cR);
			// double rhohat = 0.5*(rhoL+rhoR);
			// double Unhat = 0.5*(UnL+UnR);
			// double SL = min(0.0,min(UnL-cL,Unhat-chat));
			// double SR = max(0.0,max(UnR+cR,Unhat+chat));
			// double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						// (rhoR*(SR-UnR)-rhoL*(SL-UnL));
			// // double pStar = rhoL*(UnL-SL)*(UnL-SM)+pL; // rhoR*(UnR-SR)*(UnR-SM)+pR;
			// double pStar = 0.5*(pL+pR+rhoL*(UnL-SL)*(UnL-SM)+rhoR*(UnR-SR)*(UnR-SM));
			
			// double mdot = 0.0;
			// double pF = pStar;
			// if(SM>0.0){
				// double rhoLStar = (SL-UnL)/(SL-SM)*rhoL;
				// mdot = rhoL*UnL + SL*(rhoLStar-rhoL);
				// pF += SM/(SM-SL)*(pL-pStar);
			// }
			// else{
				// double rhoRStar = (SR-UnR)/(SR-SM)*rhoR;
				// mdot = rhoR*UnR + SR*(rhoRStar-rhoR);
				// pF += SM/(SM-SR)*(pR-pStar);
			// }
			
			// // // in pressure-based
			// // mdot -= dAlpha * dt/dLR*(pR-pL);
			// // mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
			// // mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
			// // mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
			
			// // K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
			// double UnF = 0.5*(mdot+abs(mdot))/rhoL + 0.5*(mdot-abs(mdot))/rhoR;
			// // double UnF = 0.5*(mdot+abs(mdot))/rhoFL + 0.5*(mdot-abs(mdot))/rhoFR;
			
			// double weiL = (UnF>=0.0 ? 1.0 : 0.0); double weiR = 1.0-weiL;
			
			// double rhoF = weiL*rhoFL + weiR*rhoFR;
			// double uF = weiL*uFL + weiR*uFR;
			// double vF = weiL*vFL + weiR*vFR;
			// double wF = weiL*wFL + weiR*wFR;
			// double HtF = weiL*HtFL + weiR*HtFR;
			// double YF[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// YF[i] = weiL*YFL[i] + weiR*YFR[i];
			// }
			// double dissEg = weiL*SL*(pStar-pL)/(SL-UnL) + weiR*SR*(pStar-pR)/(SR-UnR);
			
			
			// double fluxB[nEq];
			// // 컨벡티브 B
			// fluxB[0] = -( rhoF*UnF )*area;
			// fluxB[1] = -( rhoF*UnF*uF + pF*nvec[0] )*area;
			// fluxB[2] = -( rhoF*UnF*vF + pF*nvec[1] )*area;
			// fluxB[3] = -( rhoF*UnF*wF + pF*nvec[2] )*area;
			// fluxB[4] = -( rhoF*UnF*HtF  + dissEg )*area;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[5+i] = -( rhoF*UnF*YF[i] )*area;
			// }
			
			
			
            
            
            
            
            
			// // SLAU2
			// double chat = 0.5*(cL+cR);
			// double rhohat = 0.5*(rhoL+rhoR);
            // double Unhat = 0.5*(UnL+UnR);   
            // double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
            // double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
			// double KLR = sqrt(0.5*(U2L+U2R));
			// double Mcy = min(1.0,KLR/chat);
            
			// double ML = UnL/chat;
			// double MR = UnR/chat;
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// double Ubar = ( rhoFL*abs(UnL)+rhoFR*abs(UnR) ) / ( rhoFL + rhoFR );
			// double g_c = -max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
            // double mdot = 0.5*(
                // rhoFL*(UnL+(1.0-g_c)*Ubar+g_c*abs(UnL)) +
                // rhoFR*(UnR-(1.0-g_c)*Ubar-g_c*abs(UnR)) -
                // phi_c/chat*(pR-pL));
                
			// // double MLP = 0.5*(ML+abs(ML));
			// // if( abs(ML) < 1.0 ) MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			// // double MRM = 0.5*(MR-abs(MR));
            // // if(MLP+MRM>=0.0){
                // // mdot = rhoFL*(MLP+MRM)*chat - 0.5*phi_c/chat*(pR-pL);
            // // }
            // // else{
                // // mdot = rhoFR*(MLP+MRM)*chat - 0.5*phi_c/chat*(pR-pL);
            // // }
                
			// // in pressure-based
			// mdot -= dAlpha * dt/dLR*(pR-pL);
			// mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
			// mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
			// mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
            
			// double PLP = ( ML>0.0 ? 1.0 : 0.0 );
			// double PRM = ( MR<0.0 ? 1.0 : 0.0 );
			// if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            
			// // // K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
			// // double UnF = 0.5*(mdot+abs(mdot))/rhoL + 0.5*(mdot-abs(mdot))/rhoR;
            // // mdot = 0.5*(UnF+abs(UnF))*rhoFL + 0.5*(UnF-abs(UnF))*rhoFR;
			
			// // // in pressure-based
			// // UnF -= dAlpha * dt*0.5*(1.0/rhoL+1.0/rhoR)/dLR*(pR-pL);
			// // UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
			// // UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
			// // UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
            
			// // double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
			// double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
            
            // double f1L = 0.5*(mdot+abs(mdot));
            // double f1R = 0.5*(mdot-abs(mdot));
			
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
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            // // RoeM_N
            // double RT = sqrt(rhoR/rhoL); 
			// double chat = 0.5*(cL+cR);
			// // double chat = sqrt(cL*cR);
			// double rhohat = sqrt(rhoL*rhoR);//(rhoL+RT*rhoR)/(1.0+RT); 
			// double Unhat = (UnL+RT*UnR)/(1.0+RT); 
			// double uhat = (uL+RT*uR)/(1.0+RT); 
			// double vhat = (vL+RT*vR)/(1.0+RT); 
			// double what = (wL+RT*wR)/(1.0+RT); 
			// double Hthat = (HtL+RT*HtR)/(1.0+RT); 
			// double Yhat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yhat[i] = (YL[i]+RT*YR[i])/(1.0+RT); 
            // }
            // double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
            // double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
            // double w3 = min(preLs/preRs,preRs/preLs);
            // double Mhat = Unhat/chat;
            // double h_w = 1.0 - w3;
            // double f_w = 1.0;
            // if( uhat*uhat+vhat*vhat+what*what != 0.0 ) f_w=pow(abs(Mhat),h_w);
            // double g_w = 1.0;
            // if( Mhat != 0.0 ) g_w = pow(abs(Mhat),1.0-w3);

            // double fluxL[nEq], fluxR[nEq], delWdash[nEq], BdelW1dash[nEq], BdelW2dash[nEq];

			// fluxL[0] = ( rhoL*UnL );
			// fluxL[1] = ( rhoL*UnL*uL + pL*nvec[0] );
			// fluxL[2] = ( rhoL*UnL*vL + pL*nvec[1] );
			// fluxL[3] = ( rhoL*UnL*wL + pL*nvec[2] );
			// fluxL[4] = ( rhoL*UnL*HtL );
			// for(int i=0; i<nSp-1; ++i){
				// fluxL[5+i] = ( rhoL*UnL*YL[i] );
			// }
			// fluxR[0] = ( rhoR*UnR );
			// fluxR[1] = ( rhoR*UnR*uR + pR*nvec[0] );
			// fluxR[2] = ( rhoR*UnR*vR + pR*nvec[1] );
			// fluxR[3] = ( rhoR*UnR*wR + pR*nvec[2] );
			// fluxR[4] = ( rhoR*UnR*HtR );
			// for(int i=0; i<nSp-1; ++i){
				// fluxR[5+i] = ( rhoR*UnR*YR[i] );
			// }
			// delWdash[0] = ( rhoR - rhoL );
			// delWdash[1] = ( rhoR*uR - rhoL*uL );
			// delWdash[2] = ( rhoR*vR - rhoL*vL );
			// delWdash[3] = ( rhoR*wR - rhoL*wL );
			// delWdash[4] = ( rhoR*HtR - rhoL*HtL );
			// for(int i=0; i<nSp-1; ++i){
				// delWdash[5+i] = ( rhoR*YR[i] - rhoL*YL[i] );
			// }
			// BdelW1dash[0] = ( rhoR - rhoL - f_w*(pR-pL)/chat/chat );
			// BdelW1dash[1] = ( BdelW1dash[0]*uhat + rhohat*(uR-uL) );
			// BdelW1dash[2] = ( BdelW1dash[0]*vhat + rhohat*(vR-vL) );
			// BdelW1dash[3] = ( BdelW1dash[0]*what + rhohat*(wR-wL) );
			// BdelW1dash[4] = ( BdelW1dash[0]*Hthat + rhohat*(HtR-HtL) );
			// for(int i=0; i<nSp-1; ++i){
				// BdelW1dash[5+i] = ( BdelW1dash[0]*Yhat[i] + rhohat*(YR[i]-YL[i]) );
			// }
			// BdelW2dash[0] = ( 0.0 );
			// BdelW2dash[1] = ( rhohat*(UnR-UnL)*(-nvec[0]) );
			// BdelW2dash[2] = ( rhohat*(UnR-UnL)*(-nvec[1]) );
			// BdelW2dash[3] = ( rhohat*(UnR-UnL)*(-nvec[2]) );
			// BdelW2dash[4] = ( 0.0 );
			// for(int i=0; i<nSp-1; ++i){
				// BdelW2dash[5+i] = ( 0.0 );
			// }

            // double b1 = max(Unhat+chat,UnR+chat); b1 = max(b1,0.0);
            // double b2 = min(Unhat-chat,UnL-chat); b2 = min(b2,0.0);
            // double Mtilde = ( Mhat>0.0 ? 1.0 : -1.0 ) * min(1.0,abs(Mhat));
            // double b1v = max(Unhat+chat,UnR+chat); b1v = max(b1v,0.0);
            // double b2v = min(Unhat-chat,UnL-chat); b2v = min(b2v,0.0);
            // double Mvtilde = ( Unhat/chat>0.0 ? 1.0 : -1.0 )*min(1.0,abs(Unhat/chat));

            // //> comp. convective flux
            // double fluxB[nEq];
            // // cout << b1 << " " << b2 << endl;
            // for(int i=0; i<5+nSp-1; ++i){
                // double fluxTemp = (b1*fluxL[i] - b2*fluxR[i])/(b1-b2)
                    // + b1*b2/(b1-b2)*(delWdash[i]-g_w/(1.0+abs(Mtilde))*BdelW1dash[i])
                    // - b1v*b2v/(b1v-b2v)*g_w/(1.0+abs(Mvtilde))*BdelW2dash[i];
                // fluxB[i] = -fluxTemp*area;
            // }
            
            
            
            
                    
            // double RT = sqrt(rhoR/rhoL); 
			// double chat = 0.5*(cL+cR);
			// // double chat = sqrt(cL*cR);
			// double rhohat = sqrt(rhoL*rhoR);//(rhoL+RT*rhoR)/(1.0+RT); 
			// double Unhat = (UnL+RT*UnR)/(1.0+RT); 
			// double uhat = (uL+RT*uR)/(1.0+RT); 
			// double vhat = (vL+RT*vR)/(1.0+RT); 
			// double what = (wL+RT*wR)/(1.0+RT); 
			// double Hthat = (HtL+RT*HtR)/(1.0+RT); 
			// double Yhat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yhat[i] = (YL[i]+RT*YR[i])/(1.0+RT); 
            // }
            // double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
            // double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
            // double w3 = min(preLs/preRs,preRs/preLs);
            // double rhoLR = rhoL;
            // double ML = UnL/chat;
            // double MR = UnR/chat;
			// double MLP = 0.5*(ML+abs(ML));
			// if( abs(ML) < 1.0 ) MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			// double MRM = 0.5*(MR-abs(MR));
			// if( abs(MR) < 1.0 ) MRM = -0.25*(MR - 1.0)*(MR - 1.0);
			// double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            // double MLR = MLP + MRM;
            // if( MLR < 0.0 ) rhoLR = rhoR;
            // double w4 = 1.0-w3;
            // double fL = pL/(rhohat*chat*chat)*(1.0-w4)*rhohat/rhoLR;
            // double fR = pR/(rhohat*chat*chat)*(1.0-w4)*rhohat/rhoLR;
            // // =======================================
            // double MBLP, MBLM;
            // if(MLR >= 0.0) {
              // MBLP = MLP + MRM*((1.0-w4)*(1.0+fR)-fL); 
              // MBLM = MRM * w4 * (1.0+fR);
            // }
            // else{
              // MBLP = MLP * w4 * (1.0+fL); 
              // MBLM = MRM + MLP*((1.0-w4)*(1.0+fL)-fR);
            // }
            // double f1L = chat*rhoL*MBLP; 
            // double f1R = chat*rhoR*MBLM;
            // double Ku = 0.5;
            // double pu = -2.0*Ku*PLP*PRM*rhohat*chat*(UnR-UnL);
            // double PLR = PLP*pL + PRM*pR +  pu;

            // //> comp. convective flux
            // double fluxB[nEq];
            // fluxB[0] = -( f1L+f1R )*area;
            // fluxB[1] = -( f1L*uL + f1R*uR + PLR*nvec[0] )*area;
            // fluxB[2] = -( f1L*vL + f1R*vR + PLR*nvec[1] )*area;
            // fluxB[3] = -( f1L*wL + f1R*wR + PLR*nvec[2] )*area;
            // fluxB[4] = -( f1L*HtL+ f1R*HtR )*area;
            // for(int i=0; i<nSp-1; ++i){
                // fluxB[5+i] = -( f1L*YL[i]+ f1R*YR[i] )*area;
            // }
            
            
			
			
			
			// // SLAU2
			// double chat = 0.5*(cL+cR);
			// double Unhat = 0.5*(UnL+UnR);
			// double SL = min(UnL-cL,Unhat-chat);
			// double SR = max(UnR+cR,Unhat+chat);
			// double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						// (rhoR*(SR-UnR)-rhoL*(SL-UnL));
            // // chat = (rhoR*(SR-UnR)-rhoL*(SL-UnL))/(rhoL+rhoR);
			// double rhohat = 0.5*(rhoL+rhoR);
            // double U2L = uL*uL+vL*vL+wL*wL;
            // double U2R = uR*uR+vR*vR+wR*wR;
			// double KLR = sqrt(0.5*(U2L+U2R));
			// double ML = UnL/chat;
			// double MR = UnR/chat;
            
            
			// double Mcy = min(1.0,KLR/chat);
			// double phi_c = (1.0-Mcy)*(1.0-Mcy);

            // // double Mbar0 = max(sqrt(U2L)/chat,sqrt(U2R)/chat);
            // // double kL=0.0; double kR=0.0;
            // // if( (UnL-cL)>0.0 && (UnR-cR)<0.0 ) kL=1.0;
            // // if( (UnL+cL)>0.0 && (UnR+cR)<0.0 ) kL=1.0;
            // // if( (UnR-cR)>0.0 && (UnL-cL)<0.0 ) kR=1.0;
            // // if( (UnR+cR)>0.0 && (UnL+cL)<0.0 ) kR=1.0;
            // // double fM = 1.0;
            // // if(kL<1.e-200 && kR<1.e-200) fM=min(1.0,Mbar0);
            
			// // double Mcy = fM;
			// // double phi_c = (1.0-Mcy)*(1.0-Mcy);
            
			// double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );
			// double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			// double D_L = ML+(1.0-g_c)*abs(ML);
			// double D_R = MR-(1.0-g_c)*abs(MR);
			// double D_rho = Mbar*g_c;
			// double MUPL = D_L+D_rho;
			// double MUMR = D_R-D_rho;
            // double dpm = - 0.5*phi_c/chat*(pR-pL);
			// double mdot = 0.5*rhoL*MUPL*chat + 0.5*rhoR*MUMR*chat + dpm;
            
			// double MLP = 0.5*(ML+abs(ML));
			// if( abs(ML) < 1.0 ) { MLP = 0.25*(ML + 1.0)*(ML + 1.0); }
			// double MRM = 0.5*(MR-abs(MR));
			// if( abs(MR) < 1.0 ) { MRM = -0.25*(MR - 1.0)*(MR - 1.0); }
			// double PLP = ( ML>0.0 ? 1.0 : 0.0 );
			// double PRM = ( MR<0.0 ? 1.0 : 0.0 );
			// if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            
                    
            // // double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
            // // double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
            // // double fa1 = 1.0 - pow( min(preLs/preRs,preRs/preLs),2.0);

            // // double ps = PLP*pL+PRM*pR;
            // // double pll = 0.0;
            // // if( 3.0/4.0 <= min(pL/pR,pR/pL) && 1.0 > min(pL/pR,pR/pL) ){ pll=4.0*min(pL/pR,pR/pL)-3.0; }
            // // double fL = 0.0;
            // // if( abs(ML) <= 1.0 ) { fL = (pL/ps-1.0)*pll*abs(MLP)*min(1.0,pow( (abs(UnL)/chat),0.25 )); }
            // // double fR = 0.0;
            // // if( abs(MR) <= 1.0 ) { fR = (pL/ps-1.0)*pll*abs(MRM)*min(1.0,pow( (abs(UnR)/chat),0.25 )); }

            // // double MPP = MLP+MRM;
            // // double MLP_AUSM = 0.5*(MPP+abs(MPP));
            // // double MRM_AUSM = 0.5*(MPP-abs(MPP));

            // // double MLP_SLAU = 0.5*(D_L+D_rho);
            // // double MRM_SLAU = 0.5*(D_R-D_rho);
            
            // // double MLPL = fa1*MLP_AUSM + (1.0-fa1)*MLP_SLAU;
            // // double MRMR = fa1*MRM_AUSM + (1.0-fa1)*MRM_SLAU;

            // // double mdot = rhoFL*chat*MLPL + rhoFR*chat*MRMR + dpm;

            // // MLP = 0.5*(ML+abs(ML));
            // // MRM = 0.5*(MR-abs(MR));
            // // double f1L;
            // // double f1R;
            // // if( mdot >= 0.0 ) {
                // // f1L = mdot - (rhoFL*chat*MRM)*( fa1*(1.0+fR)-fR+fL );
                // // f1R = (rhoFR*chat*MRM)*( fa1*(1.0+fR) );
            // // }
            // // else{
                // // f1L = (rhoFL*chat*MLP)*( fa1*(1.0+fL) );
                // // f1R = mdot - (rhoFR*chat*MLP)*( fa1*(1.0+fL)-fL+fR );
            // // }
            
            
			// // double MLP = 0.5*(ML+abs(ML));
			// // if( abs(ML) < 1.0 ) {
				// // MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				// // // MLP += 0.125*(ML*ML-1.0)*(ML*ML-1.0);
			// // }
			// // double MRM = 0.5*(MR-abs(MR));
			// // if( abs(MR) < 1.0 ) {
				// // MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				// // // MRM -= 0.125*(MR*MR-1.0)*(MR*MR-1.0);
			// // }
			
			// // K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
			// double UnF = 0.5*(mdot+abs(mdot))/rhoL + 0.5*(mdot-abs(mdot))/rhoR;
			
			// // in pressure-based
			// UnF -= dAlpha * dt*0.5*(1.0/rhoL+1.0/rhoR)/dLR*(pR-pL);
			// UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
			// UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
			// UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
            
            // // double RT = sqrt(rhoR/rhoL); 
            // // double croe = (cL+RT*cR)/(1.0+RT); 
            // // double Unroe = (UnL+RT*UnR)/(1.0+RT); 
			// // double SL = min(UnL-cL,Unroe-croe);
			// // double SR = max(UnR+cR,Unroe+croe);
            // double SLM = min(0.0,SL);
            // double SRP = max(0.0,SR);
            // double Ustar = (SRP*UnL-SLM*UnR)/(SRP-SLM) + 2.0*dpm/(rhoL+rhoR);
            
            // // UnF = Mcy*(MLP+MRM)*chat + (1.0-Mcy)*UnF;
			
            // // UnF = UnF*(1.0-phi_c) + 0.5*(UnL+UnR)*phi_c;
				
			// // double PLP = (ML>=0.0 ? 1.0 : -1.0);
            // // double PRM = (MR>=0.0 ? 1.0 : -1.0);
			// // if( abs(ML) < 1.0 ) PLP = sin(3.141592*0.5*ML);
			// // if( abs(MR) < 1.0 ) PRM = sin(3.141592*0.5*MR);
            // // double pF = 0.5*(pL+pR) - 0.5*(0.5*(PLP+PRM)*(pR-pL)+Mcy*0.5*(pL+pR)*(PRM-PLP));
            
            // // double preLs = pL + 0.1 * min(rhoL*cL*cL,rhoR*cR*cR);
            // // double preRs = pR + 0.1 * min(rhoL*cL*cL,rhoR*cR*cR);
            // // double cpi = min(preLs/preRs,preRs/preLs);
            // // double gammaw = 0.5*(1.0-tanh(5.0*3.141592*cpi));
            // // double gamma2 = 1.0;
            // // if( dAlpha > 1.5 ) gamma2 = 1.0 - (2.0/3.0/dAlpha + 1.0/3.0);
            // // double gamma = max(0.1, max(gamma2, gammaw));
			// // double PLP = ( ML>0.0 ? 1.0 : 0.0 );
			// // double PRM = ( MR<0.0 ? 1.0 : 0.0 );
			// // if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// // if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            
			// // double pF = PLP*pL+PRM*pR;
            // // pF = pF - (1.0-fM)*(PLP+PRM-1.0)*0.5*(pL+pR);
            
			// double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL);
            // // Journal of Computational Physics 456 (2022) 11102
            // double nPd = 1.0;
            // pF += 1.0/nPd * KLR*(PLP+PRM-1.0)*rhohat*chat;
            // // pF += 1.0/nPd * Mcy*(-0.25*(abs(Unhat-chat)+abs(Unhat+chat))*rhohat*(UnR-UnL));
            // // pF += 1.0/nPd * Mcy*(rhoFL*rhoFR*(SL-UnL)*(SR-UnR)/(rhoFR*(SR-UnR)-rhoFL*(SL-UnL))*(UnR-UnL));
            // // // pF += 1.0/nPd * Mcy*(-0.25*(SRP-SLM)*0.5*(rhoL+rhoR)*(UnR-UnL));
            // // // pF += 1.0/nPd * Mcy*(-0.5*chat*0.5*(rhoL+rhoR)*(UnR-UnL));
            
            // // double dissEng = 0.0;
            // // double nEd = 2.0;
            // // dissEng += 1.0/nEd * Mcy*(-0.25*(abs(Unhat-chat)+abs(Unhat+chat))*rhohat*(UnR-UnL)*Unhat);
            // // dissEng += 1.0/nEd * Mcy*(rhoL*rhoR*(SL-UnL)*(SR-UnR)/(rhoR*(SR-UnR)-rhoL*(SL-UnL))*(UnR-UnL)*SM);
            // // // dissEng += 1.0/nEd * Mcy*(-0.25*(SRP-SLM)*0.5*(rhoL+rhoR)*0.5*(U2R-U2L));
            // // // dissEng += 1.0/nEd * Mcy*(-0.5*chat*0.5*(rhoL+rhoR)*0.5*(U2R-U2L));
            
            
			// // double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) +
						// // gamma*(PLP+PRM-1.0)*(pL+pR);
			// // double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) +
						// // Mcy*(PLP+PRM-1.0)*0.5*(pL+pR) - 
                        // // 0.3*0.5*Mcy*PLP*PRM*0.5*(pL+pR)/chat*(UnR-UnL);
            // // pF = (PLP*pL+PRM*pR);
            // // pF = pF*(1.0-phi_c) + 0.5*(pL+pR)*phi_c;
			
			// // double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
						// // (rhoR*(SR-UnR)-rhoL*(SL-UnL));
			// // double pF = 0.5*(pL+pR+rhoL*(UnL-SL)*(UnL-SM)+rhoR*(UnR-SR)*(UnR-SM));
            // // if(0.0<=SL) pF = pL; 
            // // if(SL<0.0 && 0.0<=SM) pF = (pF*SL-pL*SM)/(SL-SM);
            // // if(SM<0.0 && 0.0<=SR) pF = (pF*SR-pR*SM)/(SR-SM);
            // // if(0.0>SR) pF = pR; 
            
			
			// double weiL = (UnF>=0.0 ? 1.0 : 0.0); double weiR = 1.0-weiL;
			
			// double rhoF = weiL*rhoFL + weiR*rhoFR;
			// double uF = weiL*uFL + weiR*uFR;
			// double vF = weiL*vFL + weiR*vFR;
			// double wF = weiL*wFL + weiR*wFR;
			// double HtF = weiL*HtFL + weiR*HtFR;
			// double YF[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// YF[i] = weiL*YFL[i] + weiR*YFR[i];
			// }
			
            // double U2F = weiL*U2L + weiR*U2R;
            // double HL = HtFL-0.5*U2L;
            // double HR = HtFR-0.5*U2R;
            // double HF = weiL*HL + weiR*HR;
			
			// double fluxB[nEq];
			// // 컨벡티브 B
			// fluxB[0] = -( rhoF*UnF )*area;
			// fluxB[1] = -( rhoF*UnF*uF + pF*nvec[0] )*area;
			// fluxB[2] = -( rhoF*UnF*vF + pF*nvec[1] )*area;
			// fluxB[3] = -( rhoF*UnF*wF + pF*nvec[2] )*area;
			// fluxB[4] = -( rhoF*UnF*HtF )*area;
			// // fluxB[4] = -( rhoF*UnF*HtF + dissEng )*area;
			// // fluxB[4] = -( rhoF*( HF*Ustar + (0.5*U2F)*UnF ) )*area;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[5+i] = -( rhoF*UnF*YF[i] )*area;
			// }
			// // fluxB[0] = -( f1L + f1R )*area;
			// // fluxB[1] = -( f1L*uFL + f1R*uFR + pF*nvec[0] )*area;
			// // fluxB[2] = -( f1L*vFL + f1R*vFR + pF*nvec[1] )*area;
			// // fluxB[3] = -( f1L*wFL + f1R*wFR + pF*nvec[2] )*area;
			// // fluxB[4] = -( f1L*HtFL + f1R*HtFR )*area;
			// // // fluxB[4] = -( f1L*HtF + dissEng )*area;
			// // // fluxB[4] = -( rhoF*( HF*Ustar + (0.5*U2F)*UnF ) )*area;
			// // for(int i=0; i<nSp-1; ++i){
				// // fluxB[5+i] = -( f1L*YFL[i] + f1R*YFR[i] )*area;
			// // }
            
            
            
			
			
			
            
            
			// // TVAP, S. CHEN et al., 2021
            // double RT = sqrt(rhoR/rhoL); 
            // double croe = (cL+RT*cR)/(1.0+RT); 
            // double Unroe = (UnL+RT*UnR)/(1.0+RT); 
			// double SL = min(UnL-cL,Unroe-croe);
			// double SR = max(UnR+cR,Unroe+croe);
            // double chat = (rhoR*(SR-UnR)-rhoL*(SL-UnL))/(rhoL+rhoR);
            // double ML = UnL/chat;
            // double MR = UnR/chat;
			// double MLP = 0.5*(ML+abs(ML));
			// if( abs(ML) < 1.0 ) {
				// MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				// // MLP += 0.125*(ML*ML-1.0)*(ML*ML-1.0);
			// }
			// double MRM = 0.5*(MR-abs(MR));
			// if( abs(MR) < 1.0 ) {
				// MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				// // MRM -= 0.125*(MR*MR-1.0)*(MR*MR-1.0);
			// }
            // double UnF = (MLP+MRM)*chat;
            // double U2L = uL*uL+vL*vL+wL*wL;
            // double U2R = uR*uR+vR*vR+wR*wR;
			// double KLR = sqrt(0.5*(U2L+U2R));
			// // double Mcy = min(1.0,KLR/chat);
            // double Mbar = max(sqrt(U2L)/chat,sqrt(U2R)/chat);
            // double kL=0.0; double kR=0.0;
            // if( (UnL-cL)>0.0 && (UnR-cR)<0.0 ) kL=1.0;
            // if( (UnL+cL)>0.0 && (UnR+cR)<0.0 ) kL=1.0;
            // if( (UnR-cR)>0.0 && (UnL-cL)<0.0 ) kR=1.0;
            // if( (UnR+cR)>0.0 && (UnL+cL)<0.0 ) kR=1.0;
            // double fM = 1.0;
            // if(kL<1.e-200 && kR<1.e-200) fM=min(1.0,Mbar);
            
            // double dpm = -0.5*fM*(pR-pL)/chat;
            // double SLM = min(0.0,SL);
            // double SRP = max(0.0,SR);
            // double Ustar = (SRP*UnL-SLM*UnR)/(SRP-SLM) + 2.0*dpm/(rhoL+rhoR);
            // double HL = HtFL-0.5*U2L;
            // double HR = HtFR-0.5*U2R;
			// double PLP = (ML>=0.0 ? 1.0 : -1.0);
            // double PRM = (MR>=0.0 ? 1.0 : -1.0);
			// if( abs(ML) < 1.0 ) PLP = sin(3.141592*0.5*ML);
			// if( abs(MR) < 1.0 ) PRM = sin(3.141592*0.5*MR);
			// // double PLP = 0.5*(1.0 +1.0 * ( ML>0.0 ? 1.0 : -1.0 ) );
			// // if( abs(ML) < 1.0 ) {
				// // PLP = 0.25*pow((ML+1.0),2.0)*(2.0-ML);
			// // }
			// // double PRM = 0.5*(1.0 -1.0 * ( MR>0.0 ? 1.0 : -1.0 ) );
			// // if( abs(MR) < 1.0 ) {
				// // PRM = 0.25*pow((MR-1.0),2.0)*(2.0+MR);
			// // }
            // double pF = 0.5*(pL+pR) - 0.5*(0.5*(PLP+PRM)*(pR-pL)+fM*0.5*(pL+pR)*(PRM-PLP));
			// // double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) +
						// // KLR*(PLP+PRM-1.0)*rhohat*chat;
            // // double pF = 0.5*(pL+pR) - 0.5*(SR+SL)/(SR-SL)*(pR-pL) + SR*SL/(SR-SL)*(pR*UnR-pL*UnL)/chat/chat;           
            
            
			// double weiL = (UnF>=0.0 ? 1.0 : 0.0); double weiR = 1.0-weiL;
			// double rhoF = weiL*rhoFL + weiR*rhoFR;
			// double uF = weiL*uFL + weiR*uFR;
			// double vF = weiL*vFL + weiR*vFR;
			// double wF = weiL*wFL + weiR*wFR;
			// double HtF = weiL*HtFL + weiR*HtFR;
			// double YF[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// YF[i] = weiL*YFL[i] + weiR*YFR[i];
			// }
			// double U2F = weiL*U2L + weiR*U2R;
			// double HF = weiL*HL + weiR*HR;
            
            // double mdot = rhoF*UnF + dpm;
			
			// double fluxB[nEq];
			// // 컨벡티브 B
			// fluxB[0] = -( mdot )*area;
			// fluxB[1] = -( mdot*uF + pF*nvec[0] )*area;
			// fluxB[2] = -( mdot*vF + pF*nvec[1] )*area;
			// fluxB[3] = -( mdot*wF + pF*nvec[2] )*area;
			// fluxB[4] = -( mdot*(0.5*U2F) + rhoF*HF*Ustar )*area;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[5+i] = -( mdot*YF[i] )*area;
			// }
			
			
			
            
            
            
			// // Roe
            // double rhoL_t = rhoFL;
            // double rhoR_t = rhoFR;
            // double RT = sqrt(rhoR_t/rhoL_t); 
            // double uhat = (uL+RT*uR)/(1.0+RT); 
            // double vhat = (vL+RT*vR)/(1.0+RT);
            // double what = (wL+RT*wR)/(1.0+RT); 
            // double Hthat = (HtL+RT*HtR)/(1.0+RT); 
            // double rhohat = RT*rhoL_t;
            // // double chat= (cL+RT*cR)/(1.0+RT); 
            // // double Unhat = (UnL+RT*UnR)/(1.0+RT); 
            // double chat= 0.5*(cL+cR); 
            // double Unhat = 0.5*(UnL+UnR); 
            // double Yihat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yihat[i] = (YL[i]+RT*YR[i])/(1.0+RT);
            // }
            // //======= carefully, sensitive ===========
            // double preLs = pL + 0.1 * min(rhoL_t*cL*cL,rhoR_t*cR*cR);
            // double preRs = pR + 0.1 * min(rhoL_t*cL*cL,rhoR_t*cR*cR);
            // //========================================
            // double cpi = min(preLs/preRs,preRs/preLs);
            // double cpi3 = cpi;
            
            // // left advection flux
            // double fluxAdvL[5+nSp-1];
            // fluxAdvL[0] = (rhoL_t*UnL);
            // fluxAdvL[1] = (rhoL_t*UnL*uL + pL*nvec[0]);
            // fluxAdvL[2] = (rhoL_t*UnL*vL + pL*nvec[1]);
            // fluxAdvL[3] = (rhoL_t*UnL*wL + pL*nvec[2]);
            // fluxAdvL[4] = (rhoL_t*UnL*HtL);
			// for(int i=0; i<nSp-1; ++i){
                // fluxAdvL[5+i] = (rhoL_t*UnL*YL[i]);
            // }
            
            // // right advection flux
            // double fluxAdvR[5+nSp-1];
            // fluxAdvR[0] = (rhoR_t*UnR);
            // fluxAdvR[1] = (rhoR_t*UnR*uR + pR*nvec[0]);
            // fluxAdvR[2] = (rhoR_t*UnR*vR + pR*nvec[1]);
            // fluxAdvR[3] = (rhoR_t*UnR*wR + pR*nvec[2]);
            // fluxAdvR[4] = (rhoR_t*UnR*HtR);
			// for(int i=0; i<nSp-1; ++i){
                // fluxAdvR[5+i] = (rhoR_t*UnR*YR[i]);
            // }
                    
            // // dissipation flux
            // double DW[5+nSp-1];
            // DW[0] = ( (rhoR_t-rhoL_t) );
            // DW[1] = ( (rhoR_t*uR-rhoL_t*uL) );
            // DW[2] = ( (rhoR_t*vR-rhoL_t*vL) );
            // DW[3] = ( (rhoR_t*wR-rhoL_t*wL) );
            // DW[4] = ( (rhoR_t*(HtR-pR/rhoR_t)-rhoL_t*(HtL-pL/rhoL_t)) );
			// for(int i=0; i<nSp-1; ++i){
                // DW[5+i] = ( (rhoR_t*YR[i]-rhoL_t*YL[i]) );
            // }
            
			// double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL)+0.5*(uR*uR+vR*vR+wR*wR));
            // double Mcy = min(1.0,KLR/chat);
            // double theta = min(max(Mcy*Mcy,0.01*0.01),1.0);
            // double Um = Unhat*min(max(Mcy,0.01),1.0);
            // double Cm = chat*min(max(Mcy,0.01),1.0);
            // double thetad = min(Mcy*Mcy,1.0);
            // double Ud = 0.5*(1.0+thetad);
            // double Cd = 0.5*sqrt(4.0*chat*chat*thetad+(1.0-thetad)*(1.0-thetad)*Unhat*Unhat);
            // double Und = abs(Unhat);				
                    
            // // double xi = abs(Unhat);	double Mbar = abs(unhat)/chat;
            // double ksqrt = sqrt((uR-uL)*(uR-uL)+(vR-vL)*(vR-vL)+(wR-wL)*(wR-wL));
            // double n1x, n1y, n1z;
            // if( ksqrt < 1.e-5 ) {n1x = nvec[0]; n1y = nvec[1]; n1z = nvec[2];}
            // else{n1x = (uR-uL)/ksqrt; n1y = (vR-vL)/ksqrt; n1z = (wR-wL)/ksqrt;}
            // double n2x = - n1y*(n1x*nvec[1] - n1y*nvec[0]) - n1z*(n1x*nvec[2] - n1z*nvec[0]);
            // double n2y = n1x*(n1x*nvec[1] - n1y*nvec[0]) - n1z*(n1y*nvec[2] - n1z*nvec[1]);
            // double n2z = n1x*(n1x*nvec[2] - n1z*nvec[0]) + n1y*(n1y*nvec[2] - n1z*nvec[1]);
            // double a1=n1x*nvec[0]+n1y*nvec[1]+n1z*nvec[2];	
            // double a2=n2x*nvec[0]+n2y*nvec[1]+n2z*nvec[2];
            // double U1=n1x*uhat+n1y*vhat+n1z*what;	
            // double U2=n2x*uhat+n2y*vhat+n2z*what;
            // double frr = abs(a1*U1)+abs(a2*U2);
                     
            // double xi = max( abs(0.5*(UnR+UnL))+0.5*(UnR-UnL) , (1.0-Mcy)*Mcy*min(0.05*chat,frr) );
                   
                   
            // double DU[5+nSp-1];
            // // DU[0] = (Cm - 0.5*(1.0-theta)*Unhat*Um/Cm - theta*abs(Unhat))*(pR-pL)/rhohat/theta/chat/chat
                    // // + Ud/Cm*(UnR-UnL);
            // // DU[0] = (chat-abs(Unhat))*(pR-pL)/rhohat/chat/chat + Unhat/chat*(UnR-UnL);
            // // DU[0] = (chat-abs(Unhat))*(pR-pL)/rhohat/chat/chat;
            // DU[0] = cpi*(1.0-Mcy*Mcy)*max(0.0,chat-Und)*(pR-pL)/(chat*chat);
            // DU[0] = DU[0]*rhohat;
            // DU[1] = DU[0]*uhat; 
            // DU[2] = DU[0]*vhat;	
            // DU[3] = DU[0]*what;
            // DU[4] = DU[0]*Hthat;
			// for(int i=0; i<nSp-1; ++i){
                // DU[5+i] = DU[0]*Yihat[i];
            // }
            
            // //> X. Li et al., 2017
            // double DP[5+nSp-1];	
            // // DP[0] = (Cd - abs(Unhat) + 0.5*(1.0-thetad)*Unhat*Ud/Cm)*rhohat*(UnR-UnL) + Ud*(pR-pL)/Cm;
            // // DP[0] = Unhat/chat*(pR-pL) + (chat - abs(Unhat))*rhohat*(UnR-UnL);
            // // DP[0] = Unhat/chat*(pR-pL);
            // DP[0] = (1.0-Mcy)*max(0.0,chat - Und)*rhohat*(UnR-UnL) + 
                    // ( (Unhat > 0.0) ? 1.0 : -1.0 )*min(Und,chat)*(pR-pL)/chat;
            // DP[1] = DP[0]*nvec[0]; 
            // DP[2] = DP[0]*nvec[1];
            // DP[3] = DP[0]*nvec[2]; 
            // DP[4] = 0.0;//DP[0]*Unhat; 
            // DP[0] = 0.0;
			// for(int i=0; i<nSp-1; ++i){
                // DP[5+i] = 0.0;
            // }
            
            
            // double dUuLR[5+nSp-1];	
            // double xi_LR = 0.5*( (Unhat > 0.0) ? 1.0 : -1.0 )*min(Und,chat)*(UnR-UnL)/chat;
            // dUuLR[0] = xi_LR*(rhoL_t+rhoR_t);
            // dUuLR[1] = xi_LR*(rhoL_t*uL+rhoR_t*uR);
            // dUuLR[2] = xi_LR*(rhoL_t*vL+rhoR_t*vR);
            // dUuLR[3] = xi_LR*(rhoL_t*wL+rhoR_t*wR);
            // dUuLR[4] = xi_LR*(rhoL_t*HtL+rhoR_t*HtR);
			// for(int i=0; i<nSp-1; ++i){
                // dUuLR[5+i] = xi_LR*(rhoL_t*YL[i]+rhoR_t*YR[i]);
            // }
            
			// double fluxB[nEq];
            // for(int i=0; i<5+nSp-1; ++i){
                // fluxB[i] = -( 0.5*(fluxAdvL[i]+fluxAdvR[i]) - 0.5*(xi*DW[i] + DP[i] + DU[i] + dUuLR[i]) );
            // }
			
            
			
            
            
            
            
			// // Roe+
            
            // double gamma1 = 0.5*((pR-pL)*chat/chat-rhohat/chat*(UnR-UnL));
            // double gamma2 = ((rhoR-rhoL)-(pR-pL)/chat/chat)*nvec[0] + rhohat*(nvec[2]*(vR-vL)-nvec[1]*(wR-wL));
            // double gamma3 = ((rhoR-rhoL)-(pR-pL)/chat/chat)*nvec[1] + rhohat*(nvec[0]*(wR-wL)-nvec[2]*(uR-uL));
            // double gamma4 = ((rhoR-rhoL)-(pR-pL)/chat/chat)*nvec[2] + rhohat*(nvec[1]*(uR-uL)-nvec[0]*(vR-vL));
            // double gamma5 = 0.5*((pR-pL)*chat/chat+rhohat/chat*(UnR-UnL));
            
            
            
            
			// // Kurganov and Tadmor central scheme
            // double RT = sqrt(rhoFR/rhoFL); 
            // // double uhat = (uL+RT*uR)/(1.0+RT); 
            // // double vhat = (vL+RT*vR)/(1.0+RT);
            // // double what = (wL+RT*wR)/(1.0+RT); 
            // // double Hthat = (HtFL+RT*HtFR)/(1.0+RT); 
            // // double rhohat = RT*rhoFL;
            // // double chat= (cL+RT*cR)/(1.0+RT); 
            // // double Unhat = (UnL+RT*UnR)/(1.0+RT); 
            // double uhat = 0.5*(uL+uR);
            // double vhat = 0.5*(vL+vR);
            // double what = 0.5*(wL+wR);
            // double Hthat = 0.5*(HtFL+HtFR);
            // double rhohat = 0.5*(rhoFL+rhoFR);
            // double chat= 0.5*(cFL+cFR); 
            // double Unhat = 0.5*(UnL+UnR); 
            // double Yihat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yihat[i] = (YFL[i]+RT*YFR[i])/(1.0+RT);
            // }
            
            // // left advection flux
            // double fluxAdvL[5+nSp-1];
            // fluxAdvL[0] = (rhoFL*UnL);
            // fluxAdvL[1] = (rhoFL*UnL*uFL + pL*nvec[0]);
            // fluxAdvL[2] = (rhoFL*UnL*vFL + pL*nvec[1]);
            // fluxAdvL[3] = (rhoFL*UnL*wFL + pL*nvec[2]);
            // fluxAdvL[4] = (rhoFL*UnL*HtFL);
			// for(int i=0; i<nSp-1; ++i){
                // fluxAdvL[5+i] = (rhoFL*UnL*YFL[i]);
            // }
            
            // // right advection flux
            // double fluxAdvR[5+nSp-1];
            // fluxAdvR[0] = (rhoFR*UnR);
            // fluxAdvR[1] = (rhoFR*UnR*uFR + pR*nvec[0]);
            // fluxAdvR[2] = (rhoFR*UnR*vFR + pR*nvec[1]);
            // fluxAdvR[3] = (rhoFR*UnR*wFR + pR*nvec[2]);
            // fluxAdvR[4] = (rhoFR*UnR*HtFR);
			// for(int i=0; i<nSp-1; ++i){
                // fluxAdvR[5+i] = (rhoFR*UnR*YFR[i]);
            // }
                    
            // // dissipation flux
            // double DW[5+nSp-1];
            // DW[0] = ( (rhoFR-rhoFL) );
            // DW[1] = ( (rhoFR*uFR-rhoFL*uFL) );
            // DW[2] = ( (rhoFR*vFR-rhoFL*vFL) );
            // DW[3] = ( (rhoFR*wFR-rhoFL*wFL) );
            // DW[4] = ( (rhoFR*(HtFR-pR/rhoFR)-rhoFL*(HtFL-pL/rhoFL)) );
			// for(int i=0; i<nSp-1; ++i){
                // DW[5+i] = ( (rhoFR*YFR[i]-rhoFL*YFL[i]) );
            // }
            
            // double lambdaL = max(abs(UnL+cFL),abs(UnL-cFL));
            // double lambdaR = max(abs(UnR+cFR),abs(UnR-cFR));
            // double xi = max(lambdaL,lambdaR);
            
			// double fluxB[nEq];
            // for(int i=0; i<5+nSp-1; ++i){
                // fluxB[i] = -( 0.5*(fluxAdvL[i]+fluxAdvR[i]) - 0.5*(xi*DW[i]) );
            // }
            
            
            
			
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
			
			
			
			
			int iter;
			// double fluxB[nEq];
			
			// // 컨벡티브 B
			// fluxB[0] = -( rhoF*UnF )*area;
			// fluxB[1] = -( rhoF*UnF*uF + pF*nvec[0] )*area;
			// fluxB[2] = -( rhoF*UnF*vF + pF*nvec[1] )*area;
			// fluxB[3] = -( rhoF*UnF*wF + pF*nvec[2] )*area;
			// fluxB[4] = -( rhoF*UnF*HtF )*area;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[5+i] = -( rhoF*UnF*YF[i] )*area;
			// }
			
			
			// // HLLC
			// double Un_roe = (sqrt(rhoL)*UnL+sqrt(rhoR)*UnR)/(sqrt(rhoL)+sqrt(rhoR));
			// double c_roe = (sqrt(rhoL)*cFL+sqrt(rhoR)*cFR)/(sqrt(rhoL)+sqrt(rhoR));
			// double SL = min(UnL-cFL,Un_roe-c_roe);
			// double SR = max(UnR+cFR,Un_roe+c_roe);
			// double SM = (rhoFR*UnR*(SR-UnR)-rhoFL*UnL*(SL-UnL)+pL-pR)/
						// (rhoFR*(SR-UnR)-rhoFL*(SL-UnL));
			// double pStar = rhoL*(UnL-SL)*(UnL-SM)+pL; // rhoR*(UnR-SR)*(UnR-SM)+pR;
			// if(SL<=0.0 && 0.0<SM){ 
				// double UnStar = (SL-UnL)/(SL-SM)*SM;
				// double dpStar = (pStar-pL)/(SL-SM)*SM;
				// pF = pStar;
				
				// fluxB[0] = -( rhoFL*UnStar )*area;
				// fluxB[1] = -( rhoFL*uFL*UnStar + dpStar*nvec[0] + pF*nvec[0])*area;
				// fluxB[2] = -( rhoFL*vFL*UnStar + dpStar*nvec[1] + pF*nvec[1])*area;
				// fluxB[3] = -( rhoFL*wFL*UnStar + dpStar*nvec[2] + pF*nvec[2])*area;
				// fluxB[4] = -( rhoFL*HtFL*UnStar )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoFL*YFL[i]*UnStar )*area;
				// }
				
			// }
			// else if(SM<=0.0 && 0.0<=SR){ 
				// double UnStar = (SR-UnR)/(SR-SM)*SM;
				// double dpStar = (pStar-pR)/(SR-SM)*SM;
				// pF = pStar;
				
				// fluxB[0] = -( rhoFR*UnStar )*area;
				// fluxB[1] = -( rhoFR*uFR*UnStar + dpStar*nvec[0] + pF*nvec[0])*area;
				// fluxB[2] = -( rhoFR*vFR*UnStar + dpStar*nvec[1] + pF*nvec[1])*area;
				// fluxB[3] = -( rhoFR*wFR*UnStar + dpStar*nvec[2] + pF*nvec[2])*area;
				// fluxB[4] = -( rhoFR*HtFR*UnStar )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoFR*YFR[i]*UnStar )*area;
				// }
			// }
			// else if(SL>0.0){ 
				// pF = pL; 
				
				// fluxB[0] = -( rhoFL*UnL )*area;
				// fluxB[1] = -( rhoFL*uFL*UnL + pF*nvec[0])*area;
				// fluxB[2] = -( rhoFL*vFL*UnL + pF*nvec[1])*area;
				// fluxB[3] = -( rhoFL*wFL*UnL + pF*nvec[2])*area;
				// fluxB[4] = -( rhoFL*HtFL*UnL )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoFL*YFL[i]*UnL )*area;
				// }
			// }
			// else if(SR<0.0){ 
				// pF = pR; 
				
				// fluxB[0] = -( rhoFR*UnR )*area;
				// fluxB[1] = -( rhoFR*uFR*UnR + pF*nvec[0])*area;
				// fluxB[2] = -( rhoFR*vFR*UnR + pF*nvec[1])*area;
				// fluxB[3] = -( rhoFR*wFR*UnR + pF*nvec[2])*area;
				// fluxB[4] = -( rhoFR*HtFR*UnR )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoFR*YFR[i]*UnR )*area;
				// }
			// }
			// else{
				// cout << "#WARNING : flux error" << endl;
			// }
			
			
			
			
			
			// pF = wdL*pL + wdR*pR;
			// pF = wdL*pL_YYL+wdR*pR_YYL;
			// pF = (0.5*pL_YYL+0.5*pR_YYL) + 0.5*(PLP-PRM)*(pL_YYL-pR_YYL) + 
						// (PLP+PRM-1.0)*KLR*(wdL*rhoL+wdR*rhoR)*(wdL*cL+wdR*cR);
						// Mcy*(PLP+PRM-1.0)*(pL_YYL+pR_YYL) - 
						// 0.0;//0.5*(rhohat_YYL*(chat_YYL-abs(UnF))*(UnR-UnL));
			
			
			
			
			
			// // TV-MAS
			// double Un_roe = (sqrt(rhoL)*UnL+sqrt(rhoR)*UnR)/(sqrt(rhoL)+sqrt(rhoR));
			// double c_roe = (sqrt(rhoL)*cL+sqrt(rhoR)*cR)/(sqrt(rhoL)+sqrt(rhoR));
			// double SL = min(UnL-cL,Un_roe-c_roe);
			// double SR = max(UnR+cR,Un_roe+c_roe);
			// if(0.5*(UnL+UnR) >= 0.0){ 
				// double ck = UnL-SL;
				// double Mk = Un_roe/(Un_roe-SL);
				// double UnStar = Mk*ck;
				// double diss_p0 = SL*SR/(SR-SL)*c_roe/c_roe*(pL-pR);
				// double diss_p1 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[0] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*uL-pR*uR);
				// double diss_p2 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[1] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*vL-pR*vR);
				// double diss_p3 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[2] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*wL-pR*wR);
				// pF = 0.5*(pL+pR);
				
				
				// fluxB[0] = -( rhoL*UnStar )*area;
				// fluxB[1] = -( rhoL*uL*UnStar + pF*nvec[0] - diss_p1)*area;
				// fluxB[2] = -( rhoL*vL*UnStar + pF*nvec[1] - diss_p2)*area;
				// fluxB[3] = -( rhoL*wL*UnStar + pF*nvec[2] - diss_p3)*area;
				// fluxB[4] = -( rhoL*HtL*UnStar )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoL*YFL[i]*UnStar )*area;
				// }
				
			// }
			// else{
				// double ck = UnR-SR;
				// double Mk = Un_roe/(Un_roe-SR);
				// double UnStar = Mk*ck;
				// double diss_p0 = SL*SR/(SR-SL)*c_roe/c_roe*(pL-pR);
				// double diss_p1 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[0] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*uL-pR*uR);
				// double diss_p2 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[1] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*vL-pR*vR);
				// double diss_p3 = (SR+SL)/(2.0*(SR-SL))*(pR-pL)*nvec[2] + Mcy*SL*SR/(SR-SL)*c_roe/c_roe*(pL*wL-pR*wR);
				// pF = 0.5*(pL+pR);
				
				
				// fluxB[0] = -( rhoR*UnStar )*area;
				// fluxB[1] = -( rhoR*uR*UnStar + pF*nvec[0] - diss_p1)*area;
				// fluxB[2] = -( rhoR*vR*UnStar + pF*nvec[1] - diss_p2)*area;
				// fluxB[3] = -( rhoR*wR*UnStar + pF*nvec[2] - diss_p3)*area;
				// fluxB[4] = -( rhoR*HtR*UnStar )*area;
				// for(int i=0; i<nSp-1; ++i){
					// fluxB[5+i] = -( rhoR*YFR[i]*UnStar )*area;
				// }
				
			// }
			
			
			
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
				
			
			
			
			