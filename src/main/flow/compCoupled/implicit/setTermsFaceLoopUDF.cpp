
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
	
    
	vector<int> id_gammaYL, id_gammaYR;
    id_gammaYL.push_back(controls.getId_faceVar("left limiter-pressure"));
    id_gammaYR.push_back(controls.getId_faceVar("right limiter-pressure"));
    id_gammaYL.push_back(controls.getId_faceVar("left limiter-x-velocity"));
    id_gammaYR.push_back(controls.getId_faceVar("right limiter-x-velocity"));
    id_gammaYL.push_back(controls.getId_faceVar("left limiter-y-velocity"));
    id_gammaYR.push_back(controls.getId_faceVar("right limiter-y-velocity"));
    id_gammaYL.push_back(controls.getId_faceVar("left limiter-z-velocity"));
    id_gammaYR.push_back(controls.getId_faceVar("right limiter-z-velocity"));
    id_gammaYL.push_back(controls.getId_faceVar("left limiter-temperature"));
    id_gammaYR.push_back(controls.getId_faceVar("right limiter-temperature"));
	for(int i=0; i<controls.spName.size()-1; ++i){
        id_gammaYL.push_back(controls.getId_faceVar("left limiter-mass-fraction-"+controls.spName[i]));
        id_gammaYR.push_back(controls.getId_faceVar("right limiter-mass-fraction-"+controls.spName[i]));
    }
    
    
	
	// 크랭크-니콜슨 방법 계수
	// double CN_coeff = 0.5;
	double CN_coeff = 1.0;
	// double CN_coeff = 0.0;
	double CN_coeff_Y = 1.0;
	// double CN_coeff_Y = 0.5;
	// double CN_coeff_Y = 0.0;
	
	


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
        id_gammaYL,id_gammaYR,
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
				// double weiwd = 0.7;
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
				
			
            double gamPhiL[5+nSp], gamPhiR[5+nSp];
            for(int i=0; i<5+nSp-1; ++i){
                gamPhiL[i] = faces[id_gammaYL[i]];
                gamPhiR[i] = faces[id_gammaYR[i]];
            }
			
			

            
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
            
            if(muF>1.e-200) muF += muT_F;
            // =================================================================
            
            
            
            
            
            
            
            
            // // =================================================================
			// // AUSM-like expression of HLLC and its all-speed extension
            // // Keiichi Kitamura
			// double UnFL = uFL*nvec[0] + vFL*nvec[1] + wFL*nvec[2];
			// double UnFR = uFR*nvec[0] + vFR*nvec[1] + wFR*nvec[2];
            // double RT = sqrt(rhoFR/rhoFL); 
			// double chat = 0.5*(cFL+cFR);
			// double rhohat = sqrt(rhoFL*rhoFR);
			// double Unhat = (UnFL+RT*UnFR)/(1.0+RT); 
			// double uhat = (uFL+RT*uFR)/(1.0+RT); 
			// double vhat = (vFL+RT*vFR)/(1.0+RT); 
			// double what = (wFL+RT*wFR)/(1.0+RT); 
			// double Yhat[nSp];
			// for(int i=0; i<nSp-1; ++i){
                // Yhat[i] = (YFL[i]+RT*YFR[i])/(1.0+RT);
            // }
            
            // double preLs = abs(pL) + 0.1 * rhoFL*cFL*cFL;
            // double preRs = abs(pR) + 0.1 * rhoFR*cFR*cFR;
            // double fp = 0.0;
            // fp = min(preLs/preRs,preRs/preLs);
            // fp = fp*fp*fp;
            
            // double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
            // double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
            // double KLR = sqrt(0.5*(U2L+U2R));
            // double Mcy = min(1.0,KLR/chat);
            // double phi_c = (1.0-Mcy)*(1.0-Mcy);
            // double ML = UnFL/chat;
            // double MR = UnFR/chat;
            
            // double SL = min(0.0,min(UnFL-cFL,Unhat-chat));
            // double SR = max(0.0,max(UnFR+cFR,Unhat+chat));
            
			// double weiL; double weiR;
			
			// double rhoF;
			// double uF;
			// double vF;
			// double wF;
			// double HtF;
			// double YF[nSp];
            
			// double mdot_pL,mdot_pR,mdot_UnL,mdot_UnR,mdot_TL,mdot_TR;
			// double mdot_YL[nSp];
			// double mdot_YR[nSp];
			// double pF_pL,pF_pR,pF_UnL,pF_UnR,pF_TL,pF_TR;
			// double pF_YL[nSp];
			// double pF_YR[nSp];
			// double dissEg_pL,dissEg_pR,dissEg_UnL,dissEg_UnR,dissEg_TL,dissEg_TR;
			// double dissEg_YL[nSp];
			// double dissEg_YR[nSp];

			// double fluxB[nEq];
            // double mdot = 0.0;
            // {
                // // HLLCL
                // // double SL0 = min(UnL-cFL,UnR-cFR);
                // // double SR0 = max(UnR+cFR,UnL+cFL);
                // // double theta = (rhoFL*(UnL-SL)+rhoFR*(SR-UnR))/ (rhoFL*(UnL-SL0)+rhoFR*(SR0-UnR));
                // double theta = 1.0;
                // double SM = (rhoFL*(UnFL-SL)*UnFL + rhoFR*(SR-UnR)*UnR - theta*(pR-pL))/(rhoFL*(UnFL-SL)+rhoFR*(SR-UnR));
                // double pStar = 0.5*(pL+pR+rhoFL*(UnFL-SL)*(UnFL-SM)+rhoFR*(UnR-SR)*(UnR-SM));
                
                // double mdot = 0.0;
                // double pF = pStar;
                // if(SM>0.0){
                    // // mdot = rhoFL*(UnL + SL*((SL-UnL)/(SL-SM)-1.0));
                    // mdot = rhoFL*(UnFL + SL*(SM-UnFL)/(SL-SM));
                    
                    // pF += SM/(SM-SL)*(pL-pStar);
                // }
                // else{
                    // // mdot = rhoFR*(UnFR + SR*((SR-UnFR)/(SR-SM)-1.0));
                    // mdot = rhoFR*(UnFR + SR*(SM-UnFR)/(SR-SM));
                    
                    // pF += SM/(SM-SR)*(pR-pStar);
                // }
                
                // // // double absMhat = abs(Unhat)/chat;
                // // // mdot += (fp-1.0)*SL*SR/(SR-SL)/(1.0+absMhat)*(pR-pL)/chat/chat;
                // // // mdot -= (1.0-fp)*phi_c*SL*SR/(SR-SL)*(pR-pL)/chat/chat;
                // // mdot -= 0.5*(1.0-fp)*phi_c/chat*(pR-pL);
                // // // mdot -= 0.5*phi_c/chat*(pR-pL);
                
                
                
                // // // in pressure-based
                // // mdot -= dAlpha * dt/dLR*(pR-pL);
                // // mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
                // // mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
                // // mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
                
                
                
                // double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
                // if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
                // double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
                // if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
                // // pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
                // pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
                
                // if(mdot>=0.0){
                    // fluxB[0] = -( mdot )*area;
                    // fluxB[1] = -( mdot*uFL + pF*nvec[0] )*area;
                    // fluxB[2] = -( mdot*vFL + pF*nvec[1] )*area;
                    // fluxB[3] = -( mdot*wFL + pF*nvec[2] )*area;
                    // fluxB[4] = -( mdot*HtFL  + SL*(pStar-pL)/(SL-UnFL) )*area;
                    // // fluxB[4] = -( mdot*HtFL )*area;
                    // for(int i=0; i<nSp-1; ++i){
                        // fluxB[5+i] = -( mdot*YFL[i] )*area;
                    // }
                // }
                // else{
                    // fluxB[0] = -( mdot )*area;
                    // fluxB[1] = -( mdot*uFR + pF*nvec[0] )*area;
                    // fluxB[2] = -( mdot*vFR + pF*nvec[1] )*area;
                    // fluxB[3] = -( mdot*wFR + pF*nvec[2] )*area;
                    // fluxB[4] = -( mdot*HtFR  + SR*(pStar-pR)/(SR-UnFR) )*area;
                    // // fluxB[4] = -( mdot*HtFR )*area;
                    // for(int i=0; i<nSp-1; ++i){
                        // fluxB[5+i] = -( mdot*YFR[i] )*area;
                    // }
                // }
                
                // // // K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
                // // double UnF = 0.5*(mdot+abs(mdot))/rhoFL + 0.5*(mdot-abs(mdot))/rhoFR;
                // // // in pressure-based
                // // UnF -= dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*(pR-pL);
                // // UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
                // // UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
                // // UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
                // // if(UnF>=0.0) { mdot = rhoFL*UnF; }
                // // else { mdot = rhoFR*UnF; }
                
                
                
                
                
                
                // weiL = (mdot>=0.0 ? 1.0 : 0.0); weiR = 1.0-weiL;
                
                // rhoF = weiL*rhoFL + weiR*rhoFR;
                // uF = weiL*uFL + weiR*uFR;
                // vF = weiL*vFL + weiR*vFR;
                // wF = weiL*wFL + weiR*wFR;
                // HtF = weiL*HtFL + weiR*HtFR;
                // for(int i=0; i<nSp-1; ++i){
                    // YF[i] = weiL*YFL[i] + weiR*YFR[i];
                // }
                
                
                
                // // A matrix
                
                // // double weiYLL = gamYL[0];
                // // double weiYLR = 1.0-weiYLL;
                // // double weiYRR = gamYR[0];
                // // double weiYRL = 1.0-weiYRR;
                // // double new_drhodpL = weiYLL*drhodpL + weiYLR*drhodpR;
                // // double new_drhodpR = weiYRR*drhodpR + weiYRL*drhodpL;
                // // double new_drhodTL = weiYLL*drhodTL + weiYLR*drhodTR;
                // // double new_drhodTR = weiYRR*drhodTR + weiYRL*drhodTL;
                
				// double SM_up = (rhoFL*(UnFL-SL)*UnFL + rhoFR*(SR-UnFR)*UnFR - theta*(pR-pL));
				// double SM_down = (rhoFL*(UnFL-SL)+rhoFR*(SR-UnFR));
                
				// double dSMdpL = (drhodpL*(UnFL-SL)*UnFL+theta)/SM_down - SM_up/SM_down/SM_down*(drhodpL*(UnFL-SL));
				// double dSMdpR = (drhodpR*(SR-UnFR)*UnFR-theta)/SM_down - SM_up/SM_down/SM_down*(drhodpR*(SR-UnFR));
				// // double dSMdpL = (new_drhodpL*(UnFL-SL)*UnFL+theta)/SM_down - SM_up/SM_down/SM_down*(new_drhodpL*(UnFL-SL));
				// // double dSMdpR = (new_drhodpR*(SR-UnFR)*UnFR-theta)/SM_down - SM_up/SM_down/SM_down*(new_drhodpR*(SR-UnFR));
                
				// double dSMdUnL = (rhoFL*(UnFL-SL) + rhoFL*UnFL*(+1.0))/SM_down - SM_up/SM_down/SM_down*(rhoFL*(+1.0));
				// double dSMdUnR = (rhoFR*(SR-UnFR) + rhoFR*UnFR*(-1.0))/SM_down - SM_up/SM_down/SM_down*(rhoFR*(-1.0));
                
				// double dSMdTL = (drhodTL*UnFL*(UnFL-SL))/SM_down - SM_up/SM_down/SM_down*(drhodTL*(UnFL-SL));
				// double dSMdTR = (drhodTR*UnFR*(SR-UnFR))/SM_down - SM_up/SM_down/SM_down*(drhodTR*(SR-UnFR));
				// // double dSMdTL = (new_drhodTL*UnFL*(UnFL-SL))/SM_down - SM_up/SM_down/SM_down*(new_drhodTL*(UnFL-SL));
				// // double dSMdTR = (new_drhodTR*UnFR*(SR-UnFR))/SM_down - SM_up/SM_down/SM_down*(new_drhodTR*(SR-UnFR));
                
                // double dSMdYL[nSp], dSMdYR[nSp];
                // for(int i=0; i<nSp-1; ++i){
                    // // double new_drhodYL = weiYLL*drhodYL[i] + weiYLR*drhodYR[i];
                    // // double new_drhodYR = weiYRR*drhodYR[i] + weiYRL*drhodYL[i];
                    
                    // dSMdYL[i] = (drhodYL[i]*UnFL*(UnFL-SL))/SM_down - SM_up/SM_down/SM_down*(drhodYL[i]*(UnFL-SL));
                    // dSMdYR[i] = (drhodYR[i]*UnFR*(SR-UnFR))/SM_down - SM_up/SM_down/SM_down*(drhodYR[i]*(SR-UnFR));
                    // // dSMdYL[i] = (new_drhodYL*UnFL*(UnFL-SL))/SM_down - SM_up/SM_down/SM_down*(new_drhodYL*(UnFL-SL));
                    // // dSMdYR[i] = (new_drhodYR*UnFR*(SR-UnFR))/SM_down - SM_up/SM_down/SM_down*(new_drhodYR*(SR-UnFR));
                // }
                
                // // double tnahCoeff = 2.0;
                // // double tanhSM = tanh(tnahCoeff*SM);
                // // double weiSML = 0.5*(1.0+tanhSM);
                // // double weiSMR = 1.0-weiSML;
                // // double weidSML = 0.5*tnahCoeff*(1.0-tanhSM*tanhSM);
                // // double weidSMR = -weidSML;
                
                     
                // // { 
                    // // double mdot_L = rhoFL*(UnFL + SL*(SM-UnFL)/(SL-SM));
                    // // double mdot_R = rhoFR*(UnFR + SR*(SM-UnFR)/(SR-SM));
                    // // double mdot_pL_L = drhodpL*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdpL/(SL-SM)/(SL-SM);
                    // // double mdot_pL_R = rhoFR*SR*(SR-UnFR)*dSMdpL/(SR-SM)/(SR-SM);
                    // // double mdot_pR_L = rhoFL*SL*(SL-UnFL)*dSMdpR/(SL-SM)/(SL-SM);
                    // // double mdot_pR_R = drhodpR*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdpR/(SR-SM)/(SR-SM);
                    // // mdot_pL = weiSML*mdot_pL_L + weiSMR*mdot_pL_R + dSMdpL*(weidSML*mdot_L + weidSMR*mdot_R);
                    // // mdot_pR = weiSML*mdot_pR_L + weiSMR*mdot_pR_R + dSMdpR*(weidSML*mdot_L + weidSMR*mdot_R);
                    
                    // // double mdot_UnL_L = rhoFL*(1.0 + SL*(dSMdUnL*(SL-UnFL)+(-1.0)*(SL-SM))/(SL-SM)/(SL-SM));
                    // // double mdot_UnL_R = rhoFR*SL*dSMdUnL*(SR-UnFR)/(SR-SM)/(SR-SM);
                    // // double mdot_UnR_L = rhoFL*SL*dSMdUnR*(SL-UnFL)/(SL-SM)/(SL-SM);
                    // // double mdot_UnR_R = rhoFR*(1.0 + SR*(dSMdUnR*(SR-UnFR)+(-1.0)*(SR-SM))/(SR-SM)/(SR-SM));
                    // // mdot_UnL = weiSML*mdot_UnL_L + weiSMR*mdot_UnL_R + dSMdUnL*(weidSML*mdot_L + weidSMR*mdot_R);
                    // // mdot_UnR = weiSML*mdot_UnR_L + weiSMR*mdot_UnR_R + dSMdUnR*(weidSML*mdot_L + weidSMR*mdot_R);
                    
                    // // double mdot_TL_L = drhodTL*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdTL/(SL-SM)/(SL-SM);
                    // // double mdot_TL_R = rhoFR*SR*(SR-UnFR)*dSMdTL/(SR-SM)/(SR-SM);
                    // // double mdot_TR_L = rhoFL*SL*(SL-UnFL)*dSMdTR/(SL-SM)/(SL-SM);
                    // // double mdot_TR_R = drhodTR*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdTR/(SR-SM)/(SR-SM);
                    // // mdot_TL = weiSML*mdot_TL_L + weiSMR*mdot_TL_R + dSMdTL*(weidSML*mdot_L + weidSMR*mdot_R);
                    // // mdot_TR = weiSML*mdot_TR_L + weiSMR*mdot_TR_R + dSMdTR*(weidSML*mdot_L + weidSMR*mdot_R);
                    
                    // // for(int i=0; i<nSp-1; ++i){
                        // // double mdot_YL_L = drhodYL[i]*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdYL[i]/(SL-SM)/(SL-SM);
                        // // double mdot_YL_R = rhoFR*SR*(SR-UnFR)*dSMdYL[i]/(SR-SM)/(SR-SM);
                        // // double mdot_YR_L = rhoFL*SL*(SL-UnFL)*dSMdYR[i]/(SL-SM)/(SL-SM);
                        // // double mdot_YR_R = drhodYR[i]*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdYR[i]/(SR-SM)/(SR-SM);
                        // // mdot_YL[i] = weiSML*mdot_YL_L + weiSMR*mdot_YL_R + dSMdYL[i]*(weidSML*mdot_L + weidSMR*mdot_R);
                        // // mdot_YR[i] = weiSML*mdot_YR_L + weiSMR*mdot_YR_R + dSMdYR[i]*(weidSML*mdot_L + weidSMR*mdot_R);
                    // // }
                // // }
                            
                // if(SM>0.0){ 
                    // // mdot_pL = drhodpL*UnFL + SL*drhodpL*((SL-UnFL)/(SL-SM)-1.0) - SL*rhoFL*((SL-UnFL)/(SL-SM)/(SL-SM)*(-dSMdpL));
                    // // mdot_pR = rhoFL*SL*(SL-UnFL)/(SL-SM)/(SL-SM)*(dSMdpR);
                    // mdot_pL = drhodpL*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdpL/(SL-SM)/(SL-SM);
                    // mdot_pR = rhoFL*SL*(SL-UnFL)*dSMdpR/(SL-SM)/(SL-SM);
                    
                    // // mdot_UnFL = rhoFL*(1.0 + SL*(-1.0/(SL-SM) - (SL-UnFL)/(SL-SM)/(SL-SM)*(-dSMdUnFL)));
                    // // mdot_UnFR = rhoFL*SL*(- (SL-UnFL)/(SL-SM)/(SL-SM)*(-dSMdUnFR));
                    // mdot_UnL = rhoFL*(1.0 + SL*(dSMdUnL*(SL-UnFL)+(-1.0)*(SL-SM))/(SL-SM)/(SL-SM));
                    // mdot_UnR = rhoFL*SL*dSMdUnR*(SL-UnFL)/(SL-SM)/(SL-SM);
                    
                    // // mdot_TL = drhodTL*UnFL + SL*drhodTL*((SL-UnFL)/(SL-SM)-1.0) - SL*rhoFL*((SL-UnFL)/(SL-SM)/(SL-SM)*(-dSMdTL));
                    // // mdot_TR = rhoFL*SL*(SL-UnFL)/(SL-SM)/(SL-SM)*(dSMdTR);
                    // mdot_TL = drhodTL*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdTL/(SL-SM)/(SL-SM);
                    // mdot_TR = rhoFL*SL*(SL-UnFL)*dSMdTR/(SL-SM)/(SL-SM);
                    
                    // for(int i=0; i<nSp-1; ++i){
                        // // mdot_YL[i] = drhodYL[i]*UnFL + SL*drhodYL[i]*((SL-UnFL)/(SL-SM)-1.0) - 
                                     // // SL*rhoFL*((SL-UnFL)/(SL-SM)/(SL-SM)*(-dSMdYL[i]));
                        // // mdot_YR[i] = rhoFL*SL*(SL-UnFL)/(SL-SM)/(SL-SM)*(dSMdYR[i]);
                        // mdot_YL[i] = drhodYL[i]*(UnFL+SL*(SM-UnFL)/(SL-SM)) + rhoFL*SL*(SL-UnFL)*dSMdYL[i]/(SL-SM)/(SL-SM);
                        // mdot_YR[i] = rhoFL*SL*(SL-UnFL)*dSMdYR[i]/(SL-SM)/(SL-SM);
                    // }
                    
                    
                    // // mdot_pL += dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*rhoFL;
                    // // mdot_pR -= dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*rhoFL;
                    
                    
                // }
                // else{
                    // // mdot_pL = rhoFR*SR*(SR-UnFR)/(SR-SM)/(SR-SM)*(dSMdpL);
                    // // mdot_pR = drhodpR*UnFR + SR*drhodpR*((SR-UnFR)/(SR-SM)-1.0) - SR*rhoFR*((SR-UnFR)/(SR-SM)/(SR-SM)*(-dSMdpR));
                    // mdot_pL = rhoFR*SR*(SR-UnFR)*dSMdpL/(SR-SM)/(SR-SM);
                    // mdot_pR = drhodpR*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdpR/(SR-SM)/(SR-SM);
                    
                    // // mdot_UnFL = rhoFR*SR*(- (SR-UnFR)/(SR-SM)/(SR-SM)*(-dSMdUnFL));
                    // // mdot_UnFR = rhoFR*(1.0 + SR*(-1.0/(SR-SM) - (SR-UnFR)/(SR-SM)/(SR-SM)*(-dSMdUnFR)));
                    // mdot_UnL = rhoFR*SL*dSMdUnL*(SR-UnFR)/(SR-SM)/(SR-SM);
                    // mdot_UnR = rhoFR*(1.0 + SR*(dSMdUnR*(SR-UnFR)+(-1.0)*(SR-SM))/(SR-SM)/(SR-SM));
                    
                    // // mdot_TL = rhoFR*SR*(SR-UnFR)/(SR-SM)/(SR-SM)*(dSMdTL);
                    // // mdot_TR = drhodTR*UnFR + SR*drhodTR*((SR-UnFR)/(SR-SM)-1.0) - SR*rhoFR*((SR-UnFR)/(SR-SM)/(SR-SM)*(-dSMdTR));
                    // mdot_TL = rhoFR*SR*(SR-UnFR)*dSMdTL/(SR-SM)/(SR-SM);
                    // mdot_TR = drhodTR*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdTR/(SR-SM)/(SR-SM);
                    
                    // for(int i=0; i<nSp-1; ++i){
                        // // mdot_YL[i] = rhoFR*SR*(SR-UnFR)/(SR-SM)/(SR-SM)*(dSMdYL[i]);
                        // // mdot_YR[i] = drhodYR[i]*UnFR + SR*drhodYR[i]*((SR-UnFR)/(SR-SM)-1.0) - 
                                     // // SR*rhoFR*((SR-UnFR)/(SR-SM)/(SR-SM)*(-dSMdYR[i]));
                        // mdot_YL[i] = rhoFR*SR*(SR-UnFR)*dSMdYL[i]/(SR-SM)/(SR-SM);
                        // mdot_YR[i] = drhodYR[i]*(UnFR+SR*(SM-UnFR)/(SR-SM)) + rhoFR*SR*(SR-UnFR)*dSMdYR[i]/(SR-SM)/(SR-SM);
                    // }
                    
                    
                    
                    // // mdot_pL += dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*rhoFR;
                    // // mdot_pR -= dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*rhoFR;
                    
                    
                // }
                // // // mdot_pL += 0.5*(1.0-fp)*phi_c/chat;
                // // // mdot_pR -= 0.5*(1.0-fp)*phi_c/chat;
                
                // mdot_pL += dAlpha * dt/dLR;
                // mdot_pR -= dAlpha * dt/dLR;
                
                
                
                
                
                
                // // dissEg_pL = 0.0; dissEg_pR = 0.0;
                // // dissEg_TL = 0.0; dissEg_TR = 0.0;
                // // for(int i=0; i<nSp-1; ++i){
                    // // dissEg_YL[i] = 0.0; dissEg_YR[i] = 0.0;
                    
                // // }
                // // dissEg_UnFL = 0.0; dissEg_UnFR = 0.0;
            
            
				// // double pStar_pL = 0.5*( 1.0 + drhodpL*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdpL) +
								  // // rhoFR*(UnFR-SR)*(-dSMdpL) );
				// // double pStar_pR = 0.5*( 1.0 + rhoFL*(UnFL-SL)*(-dSMdpR) + drhodpR*(UnFR-SR)*(UnFR-SM) +
								  // // rhoFR*(UnFR-SR)*(-dSMdpR) );
				// // double pStar_UnFL = 0.5*(rhoFL*(1.0)*(UnFL-SM)+rhoFL*(UnFL-SL)*(1.0-dSMdUnFL)+rhoFR*(UnFR-SR)*(-dSMdUnFL));
				// // double pStar_UnFR = 0.5*(rhoFR*(1.0)*(UnFR-SM)+rhoFR*(UnFR-SR)*(1.0-dSMdUnFR)+rhoFL*(UnFL-SL)*(-dSMdUnFR));
				// // double pStar_TL = 0.5*(drhodTL*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdTL) +
								  // // rhoFR*(UnFR-SR)*(-dSMdTL) );
				// // double pStar_TR = 0.5*(rhoFL*(UnFL-SL)*(-dSMdTR) + drhodTR*(UnFR-SR)*(UnFR-SM) +
								  // // rhoFR*(UnFR-SR)*(-dSMdTR) );
				// // double pStar_YL[nSp];
				// // double pStar_YR[nSp];
				// // for(int i=0; i<nSp-1; ++i){
					// // pStar_YL[i] = 0.5*(drhodYL[i]*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdYL[i]) +
									  // // rhoFR*(UnFR-SR)*(-dSMdYL[i]) );
					// // pStar_YR[i] = 0.5*(rhoFL*(UnFL-SL)*(-dSMdYR[i]) + drhodYR[i]*(UnFR-SR)*(UnFR-SM) +
									  // // rhoFR*(UnFR-SR)*(-dSMdYR[i]) );
				// // }
            
				// double pStar_pL = 0.5*( 1.0 + drhodpL*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdpL) + rhoFR*(UnFR-SR)*(-dSMdpL) );
				// double pStar_pR = 0.5*( 1.0 + rhoFL*(UnFL-SL)*(-dSMdpR) + drhodpR*(UnFR-SR)*(UnFR-SM) + rhoFR*(UnFR-SR)*(-dSMdpR) );
				// double pStar_UnL = 0.5*(rhoFL*(1.0)*(UnFL-SM)+rhoFL*(UnFL-SL)*(1.0-dSMdUnL)+rhoFR*(UnFR-SR)*(-dSMdUnL));
				// double pStar_UnR = 0.5*(rhoFR*(1.0)*(UnFR-SM)+rhoFR*(UnFR-SR)*(1.0-dSMdUnR)+rhoFL*(UnFL-SL)*(-dSMdUnR));
				// double pStar_TL = 0.5*(drhodTL*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdTL) + rhoFR*(UnFR-SR)*(-dSMdTL) );
				// double pStar_TR = 0.5*(rhoFL*(UnFL-SL)*(-dSMdTR) + drhodTR*(UnFR-SR)*(UnFR-SM) + rhoFR*(UnFR-SR)*(-dSMdTR) );
				// double pStar_YL[nSp];
				// double pStar_YR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// pStar_YL[i] = 0.5*(drhodYL[i]*(UnFL-SL)*(UnFL-SM) + rhoFL*(UnFL-SL)*(-dSMdYL[i]) + rhoFR*(UnFR-SR)*(-dSMdYL[i]) );
					// pStar_YR[i] = 0.5*(rhoFL*(UnFL-SL)*(-dSMdYR[i]) + drhodYR[i]*(UnFR-SR)*(UnFR-SM) + rhoFR*(UnFR-SR)*(-dSMdYR[i]) );
				// }
				
								  
                // dissEg_pL = weiL*SL*(pStar_pL-1.0)/(SL-UnFL) + weiR*SR*(pStar_pL)/(SR-UnFR);
                // dissEg_pR = weiL*SL*(pStar_pR)/(SL-UnFL) + weiR*SR*(pStar_pR-1.0)/(SR-UnFR);
                // dissEg_TL = weiL*SL*(pStar_TL)/(SL-UnFL) + weiR*SR*(pStar_TL)/(SR-UnFR);
                // dissEg_TR = weiL*SL*(pStar_TR)/(SL-UnFL) + weiR*SR*(pStar_TR)/(SR-UnFR);
                // for(int i=0; i<nSp-1; ++i){
                    // dissEg_YL[i] = weiL*SL*(pStar_YL[i])/(SL-UnFL) + weiR*SR*(pStar_YL[i])/(SR-UnFR);
                    // dissEg_YR[i] = weiL*SL*(pStar_YR[i])/(SL-UnFL) + weiR*SR*(pStar_YR[i])/(SR-UnFR);
                // }
                // dissEg_UnL = weiL*SL*(pStar_UnL)/(SL-UnFL) - weiL*SL*(pStar-pL)/(SL-UnFL)/(SL-UnFL)*(-1.0) + 
                             // weiR*SR*(pStar_UnL)/(SR-UnFR);
                // dissEg_UnR = weiL*SL*(pStar_UnR)/(SL-UnFL) - weiR*SR*(pStar-pR)/(SR-UnFR)/(SR-UnFR)*(-1.0) + 
                             // weiR*SR*(pStar_UnR)/(SR-UnFR);
                
                
                
                
                
                // // pressure term
                // pF_pL = 0.5 - 0.5*(PLP-PRM)*(-1.0) + Mcy*(PLP+PRM-1.0)*0.5;
                // pF_pR = 0.5 - 0.5*(PLP-PRM)*(+1.0) + Mcy*(PLP+PRM-1.0)*0.5;
                
                // // double PLP_UnFL = 0.0;
                // // if( abs(ML) < 1.0 ) {
                    // // PLP_UnFL = 0.25*(1.0/chat)*(ML+1.0)*(2.0-ML) + 
                                // // 0.25*(ML+1.0)*(1.0/chat)*(2.0-ML) + 
                                // // 0.25*(ML+1.0)*(ML+1.0)*(-1.0/chat);
                // // }
                // // double PRM_UnFR = 0.0;
                // // if( abs(MR) < 1.0 ) {
                    // // PRM_UnFR = 0.25*(1.0/chat)*(MR-1.0)*(2.0+MR) + 
                                // // 0.25*(MR-1.0)*(1.0/chat)*(2.0+MR) + 
                                // // 0.25*(MR-1.0)*(MR-1.0)*(1.0/chat);
                // // }
                // // pF_UnL = -0.5*(+PLP_UnFL)*(pR-pL) + Mcy*(PLP_UnFL)*0.5*(pL+pR);
                // // pF_UnR = -0.5*(-PRM_UnFR)*(pR-pL) + Mcy*(PRM_UnFR)*0.5*(pL+pR);
                // pF_UnL = 0.0;
                // pF_UnR = 0.0;
                
                // pF_TL = 0.0;
                // pF_TR = 0.0;
                
                // for(int i=0; i<nSp-1; ++i){
                    // pF_YL[i] = 0.0;
                    // pF_YR[i] = 0.0;
                // }
                
                
            // }
            // // =================================================================
            
            
            
            
            
            
            
			
			
			
			
			
			
			
			
			
			
            // =================================================================
			// SLAU2
			double UnFL = uFL*nvec[0] + vFL*nvec[1] + wFL*nvec[2];
			double UnFR = uFR*nvec[0] + vFR*nvec[1] + wFR*nvec[2];
			double chat = 0.5*(cFL+cFR);
			double rhohat = 0.5*(rhoFL+rhoFR);    
            double U2L = uFL*uFL+vFL*vFL+wFL*wFL;
            double U2R = uFR*uFR+vFR*vFR+wFR*wFR;
            double Ubar = ( rhoFL*abs(UnFL)+rhoFR*abs(UnFR) ) / ( rhoFL + rhoFR );
            
			double Unhat = 0.5*(UnFL+UnFR);        
            double KLR = sqrt(0.5*(U2L+U2R));
			double Mcy = min(1.0,KLR/chat);
			double phi_c = (1.0-Mcy)*(1.0-Mcy);
			double ML = UnFL/chat;
			double MR = UnFR/chat;			
			double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			double D_L = UnFL+(1.0-g_c)*abs(UnFL);
			double D_R = UnFR-(1.0-g_c)*abs(UnFR);
			double D_rho = Ubar*g_c;
			double UPL = D_L+D_rho;
			double UMR = D_R-D_rho;
			double mdot = 0.5*rhoFL*UPL + 0.5*rhoFR*UMR - 0.5*phi_c/chat*(pR-pL);
            // UPL = (0.5*(UnL+UnR)+abs(0.5*(UnL+UnR)));
            // UMR = (0.5*(UnL+UnR)-abs(0.5*(UnL+UnR)));
			// double mdot = 0.5*rhoFL*UPL + 0.5*rhoFR*UMR;
			
			// in pressure-based
			mdot -= dAlpha * dt/dLR*(pR-pL);
			mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
			mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
			mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
			
			// K. Kitamura, E. Shima / Computers and Fluids 163 (2018) 86–96
			// double UnF = 0.5*(mdot+abs(mdot))/rhoL + 0.5*(mdot-abs(mdot))/rhoR;
			// double UnF = 0.5*(UnL+UnR);
			// // in pressure-based
			// UnF -= dAlpha * dt*(1.0/rhoL+1.0/rhoR)*(pR-pL)/dLR;
			// UnF += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
			// UnF += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
			// UnF += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
            // mdot = 0.5*(UnF+abs(UnF))*rhoFL + 0.5*(UnF-abs(UnF))*rhoFR;
            
            double PLP = ( ML>0.0 ? 1.0 : 0.0 );
			double PRM = ( MR<0.0 ? 1.0 : 0.0 );
			if( abs(ML) < 1.0 ) PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			if( abs(MR) < 1.0 ) PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
            
			// double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + KLR*(PLP+PRM-1.0)*rhohat*chat;
			double pF = 0.5*(pL+pR) - 0.5*(PLP-PRM)*(pR-pL) + Mcy*Mcy*(PLP+PRM-1.0)*0.5*(pL+pR);
			// double pF = 0.5*(pL+pR);
			// double pF = PLP*pL+PRM*pR;
			
			
			double weiL = (mdot>=0.0 ? 1.0 : 0.0); double weiR = 1.0-weiL;
			
			double rhoF = weiL*rhoFL + weiR*rhoFR;
			double uF = weiL*uFL + weiR*uFR;
			double vF = weiL*vFL + weiR*vFR;
			double wF = weiL*wFL + weiR*wFR;
			double HtF = weiL*HtFL + weiR*HtFR;
			double YF[nSp];
			for(int i=0; i<nSp-1; ++i){
				YF[i] = weiL*YFL[i] + weiR*YFR[i];
			}
			
			
			double fluxB[nEq];
			// 컨벡티브 B
			fluxB[0] = -( mdot )*area;
			fluxB[1] = -( mdot*uF + pF*nvec[0] )*area;
			fluxB[2] = -( mdot*vF + pF*nvec[1] )*area;
			fluxB[3] = -( mdot*wF + pF*nvec[2] )*area;
			fluxB[4] = -( mdot*HtF )*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[5+i] = -( mdot*YF[i] )*area;
			}
			
			
			
			double sign_UnL = (UnFL>=0.0 ? 1.0 : -1.0);
			double sign_UnR = (UnFR>=0.0 ? 1.0 : -1.0);
			double D_rho_dL = g_c*(abs(UnFL)-Ubar)/(rhoFL+rhoFR);
			double D_rho_dR = g_c*(abs(UnFR)-Ubar)/(rhoFL+rhoFR);
			// double D_rho_dL = g_c*(abs(UnL)-Ubar)/(rhoL+rhoR);
			// double D_rho_dR = g_c*(abs(UnR)-Ubar)/(rhoL+rhoR);
			
			// double mdot_pL = 0.5*( drhodpL*UPL + rhoFL*(+drhodpL*D_rho_dL) + rhoFR*(-drhodpL*D_rho_dL) + phi_c/chat );
			// double mdot_pR = 0.5*( drhodpR*UMR + rhoFR*(-drhodpR*D_rho_dR) + rhoFL*(+drhodpR*D_rho_dR) - phi_c/chat );
			// double mdot_pL = 0.5*( drhodpL*UPL + rhoFL*(+drhodpL*D_rho_dL) + phi_c/chat );
			// double mdot_pR = 0.5*( drhodpR*UMR + rhoFR*(-drhodpR*D_rho_dR) - phi_c/chat );
			// double mdot_pL = 0.5*( drhodpL*UPL + phi_c/chat );
			// double mdot_pR = 0.5*( drhodpR*UMR - phi_c/chat );
			double mdot_pL = 0.5*( drhodpL*UPL );
			double mdot_pR = 0.5*( drhodpR*UMR );
            mdot_pL += phi_c/chat; mdot_pR -= phi_c/chat;
            mdot_pL += dAlpha * dt/dLR; mdot_pR -= dAlpha * dt/dLR;
            // if(mdot>=0.0){
                // mdot_pL += dAlpha * dt*rhoFL*(1.0/rhoL+1.0/rhoR)/dLR;
                // mdot_pR -= dAlpha * dt*rhoFL*(1.0/rhoL+1.0/rhoR)/dLR;
            // }
            // else{
                // mdot_pL += dAlpha * dt*rhoFR*(1.0/rhoL+1.0/rhoR)/dLR;
                // mdot_pR -= dAlpha * dt*rhoFR*(1.0/rhoL+1.0/rhoR)/dLR;
            // }
            
			// double UnF_pL = mdot_pL*(weiL/rhoL+weiR/rhoR) + mdot*drhodpL*(-weiL/rhoL/rhoL);
			// double UnF_pR = mdot_pR*(weiL/rhoL+weiR/rhoR) + mdot*drhodpR*(-weiR/rhoR/rhoR);
			// mdot_pL = weiL*drhodpL*UnF + rhoF*UnF_pL;
			// mdot_pR = weiR*drhodpR*UnF + rhoF*UnF_pR;
            
            
			
			// double mdot_UnL = 0.5*( rhoFL*(1.0+(1.0-g_c)*sign_UnL+g_c*rhoFL/(rhoL+rhoR)*sign_UnL) + 
                                    // rhoFR*(-rhoFR*sign_UnR/(rhoL+rhoR)*g_c) );
			// double mdot_UnR = 0.5*( rhoFR*(1.0-(1.0-g_c)*sign_UnR-g_c*rhoFR/(rhoL+rhoR)*sign_UnR) +
                                    // rhoFL*(+rhoFL*sign_UnL/(rhoL+rhoR)*g_c) );
			// double mdot_UnL = 0.5*( rhoFL*(1.0+(1.0-g_c)*sign_UnL+g_c*rhoFL/(rhoFL+rhoFR)*sign_UnL) );
			// double mdot_UnR = 0.5*( rhoFR*(1.0-(1.0-g_c)*sign_UnR-g_c*rhoFR/(rhoFL+rhoFR)*sign_UnR) );
			// double mdot_UnL = 0.5*( rhoFL*(1.0+(1.0-g_c)*sign_UnL+g_c*rhoL/(rhoL+rhoR)*sign_UnL) );
			// double mdot_UnR = 0.5*( rhoFR*(1.0-(1.0-g_c)*sign_UnR-g_c*rhoR/(rhoL+rhoR)*sign_UnR) );
			double mdot_UnL = 0.5*( rhoFL );
			double mdot_UnR = 0.5*( rhoFR );
			// double UnF_UnL = mdot_UnL*(weiL/rhoL+weiR/rhoR);
			// double UnF_UnR = mdot_UnR*(weiL/rhoL+weiR/rhoR);
			// mdot_UnL = rhoF*UnF_UnL;
			// mdot_UnR = rhoF*UnF_UnR;
            
			
			// double mdot_TL = 0.5*( drhodTL*UPL + rhoFL*(+drhodTL*D_rho_dL) + rhoFR*(-drhodTL*D_rho_dL) );
			// double mdot_TR = 0.5*( drhodTR*UMR + rhoFR*(-drhodTR*D_rho_dR) + rhoFL*(+drhodTR*D_rho_dR) );
			// double mdot_TL = 0.5*( drhodTL*UPL + rhoFL*(+drhodTL*D_rho_dL) );
			// double mdot_TR = 0.5*( drhodTR*UMR + rhoFR*(-drhodTR*D_rho_dR) );
			double mdot_TL = 0.5*( drhodTL*UPL );
			double mdot_TR = 0.5*( drhodTR*UMR );
			// double UnF_TL = mdot_TL*(weiL/rhoL+weiR/rhoR) + mdot*drhodTL*(-weiL/rhoL/rhoL);
			// double UnF_TR = mdot_TR*(weiL/rhoL+weiR/rhoR) + mdot*drhodTR*(-weiR/rhoR/rhoR);
			// mdot_TL = weiL*drhodTL*UnF + rhoF*UnF_TL;
			// mdot_TR = weiR*drhodTR*UnF + rhoF*UnF_TR;
			
			double mdot_YL[nSp];
			double mdot_YR[nSp];
			// double UnF_YL[nSp];
			// double UnF_YR[nSp];
			for(int i=0; i<nSp-1; ++i){
                
                // mdot_YL[i] = 0.5*( drhodYL[i]*UPL + rhoFL*(+drhodYL[i]*D_rho_dL) + rhoFR*(-drhodYL[i]*D_rho_dL) );
                // mdot_YR[i] = 0.5*( drhodYR[i]*UMR + rhoFR*(-drhodYR[i]*D_rho_dR) + rhoFL*(+drhodYR[i]*D_rho_dR) );
                // mdot_YL[i] = 0.5*( drhodYL[i]*UPL + rhoFL*(+drhodYL[i]*D_rho_dL) );
                // mdot_YR[i] = 0.5*( drhodYR[i]*UMR + rhoFR*(-drhodYR[i]*D_rho_dR) );
                mdot_YL[i] = 0.5*( drhodYL[i]*UPL );
                mdot_YR[i] = 0.5*( drhodYR[i]*UMR );
				// mdot_YL[i] = 0.5*drhodYL[i]*(UPL + D_rho_dL);
				// mdot_YR[i] = 0.5*drhodYR[i]*(UMR - D_rho_dR);
				// UnF_YL[i] = mdot_YL[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYL[i]*(-weiL/rhoL/rhoL);
				// UnF_YR[i] = mdot_YR[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYR[i]*(-weiR/rhoR/rhoR);
				// mdot_YL[i] = weiL*drhodYL[i]*UnF + rhoF*UnF_YL[i];
				// mdot_YR[i] = weiR*drhodYR[i]*UnF + rhoF*UnF_YR[i];
				
				
				// //=========================
				// double tmp_weiL = 1.0; double tmp_weiR = 0.0;
				// if(YF[i] < YL[i]+1.e-50 && YF[i] > YL[i]-1.e-50){
					// tmp_weiL = 1.0; tmp_weiR = 0.0;
				// }
				// else{
					// tmp_weiL = 0.0; tmp_weiR = 1.0;
				// }
				// mdot_YL[i] = 0.5*drhodYL[i]*(UPL + D_rho_dL);
				// mdot_YR[i] = 0.5*drhodYR[i]*(UMR - D_rho_dR);
				// UnF_YL[i] = mdot_YL[i]*(tmp_weiL/rhoL+tmp_weiR/rhoR) + mdot*drhodYL[i]*(-tmp_weiL/rhoL/rhoL);
				// UnF_YR[i] = mdot_YR[i]*(tmp_weiL/rhoL+tmp_weiR/rhoR) + mdot*drhodYR[i]*(-tmp_weiR/rhoR/rhoR);
				// mdot_YL[i] = tmp_weiL*drhodYL[i]*UnF + rhoF*UnF_YL[i];
				// mdot_YR[i] = tmp_weiR*drhodYR[i]*UnF + rhoF*UnF_YR[i];
				// //=========================
				
				
				
			}
			
			// double mdot_UnL = 0.5*rhoL*(1.0+(1.0-g_c)*sign_UnL+g_c*rhoL/(rhoL+rhoR)*sign_UnL);
			// double mdot_UnR = 0.5*rhoR*(1.0+(1.0-g_c)*sign_UnR-g_c*rhoR/(rhoL+rhoR)*sign_UnR);
			// double UnF_UnL = mdot_UnL*(weiL/rhoL+weiR/rhoR);
			// double UnF_UnR = mdot_UnR*(weiL/rhoL+weiR/rhoR);
			// mdot_UnL = rhoF*UnF_UnL;
			// mdot_UnR = rhoF*UnF_UnR;
			
			// double mdot_TL = 0.5*drhodTL*(UPL + D_rho_dL);
			// double mdot_TR = 0.5*drhodTR*(UMR - D_rho_dR);
			// double UnF_TL = mdot_TL*(weiL/rhoL+weiR/rhoR) + mdot*drhodTL*(-weiL/rhoL/rhoL);
			// double UnF_TR = mdot_TR*(weiL/rhoL+weiR/rhoR) + mdot*drhodTR*(-weiR/rhoR/rhoR);
			// mdot_TL = weiL*drhodTL*UnF + rhoF*UnF_TL;
			// mdot_TR = weiR*drhodTR*UnF + rhoF*UnF_TR;
			
			// double mdot_YL[nSp];
			// double mdot_YR[nSp];
			// double UnF_YL[nSp];
			// double UnF_YR[nSp];
			// for(int i=0; i<nSp-1; ++i){
				
				// mdot_YL[i] = 0.5*drhodYL[i]*(UPL + D_rho_dL);
				// mdot_YR[i] = 0.5*drhodYR[i]*(UMR - D_rho_dR);
				// UnF_YL[i] = mdot_YL[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYL[i]*(-weiL/rhoL/rhoL);
				// UnF_YR[i] = mdot_YR[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYR[i]*(-weiR/rhoR/rhoR);
				// mdot_YL[i] = weiL*drhodYL[i]*UnF + rhoF*UnF_YL[i];
				// mdot_YR[i] = weiR*drhodYR[i]*UnF + rhoF*UnF_YR[i];
			// }
            
            
            
            
            // pressure term
			// double dcdphi_coeffL = (cL/2.0/rhoL*(1.0-cL*cL*(drhodpL-drhodTL*dHtdpL/dHtdTL)));
			// double dcdphi_coeffR = (cR/2.0/rhoR*(1.0-cR*cR*(drhodpR-drhodTR*dHtdpR/dHtdTR)));
			
			// double pF_pL = 0.5 - 0.5*(PLP-PRM)*(-1.0) + KLR*(PLP+PRM-1.0)*0.5*drhodpL*chat;
			// double pF_pR = 0.5 - 0.5*(PLP-PRM)*(+1.0) + KLR*(PLP+PRM-1.0)*0.5*drhodpR*chat;
			// double pF_pL = 0.5 - 0.5*(PLP-PRM)*(-1.0) + Mcy*(PLP+PRM-1.0)*0.5;
			// double pF_pR = 0.5 - 0.5*(PLP-PRM)*(+1.0) + Mcy*(PLP+PRM-1.0)*0.5;
			double pF_pL = 0.5 - 0.5*(PLP-PRM)*(-1.0) + Mcy*Mcy*(PLP+PRM-1.0)*0.5;
			double pF_pR = 0.5 - 0.5*(PLP-PRM)*(+1.0) + Mcy*Mcy*(PLP+PRM-1.0)*0.5;
			// double pF_pL = 0.5;
			// double pF_pR = 0.5;
			// double pF_pL = PLP;
			// double pF_pR = PRM;
			// double dcdpL = drhodpL*dcdphi_coeffL;
			// double PLP_pL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdpL);
			// pF_pL -= 0.5*(pR-pL)*PLP_pL;
			// pF_pL += KLR*rhohat*chat*PLP_pL;
			// pF_pL += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdpL;
			// double dcdpR = drhodpR*dcdphi_coeffR;
			// double PRM_pR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdpR);
			// pF_pR -= 0.5*(pR-pL)*PRM_pR;
			// pF_pR += KLR*rhohat*chat*PRM_pR;
			// pF_pR += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdpR;
			
			// double PLP_UnL = 0.0;
            // if( abs(ML) < 1.0 ) {
                // PLP_UnL = 0.25*(1.0/chat)*(ML+1.0)*(2.0-ML) + 
                            // 0.25*(ML+1.0)*(1.0/chat)*(2.0-ML) + 
                            // 0.25*(ML+1.0)*(ML+1.0)*(-1.0/chat);
            // }
			// double PRM_UnR = 0.0;
            // if( abs(MR) < 1.0 ) {
                // PRM_UnR = 0.25*(1.0/chat)*(MR-1.0)*(2.0+MR) + 
                            // 0.25*(MR-1.0)*(1.0/chat)*(2.0+MR) + 
                            // 0.25*(MR-1.0)*(MR-1.0)*(1.0/chat);
            // }
			// // double pF_UnL = -0.5*(PLP_UnL)*(pR-pL) + KLR*(PLP_UnL)*rhohat*chat;
			// // double pF_UnR = -0.5*(PLP_UnR)*(pR-pL) + KLR*(PLP_UnR)*rhohat*chat;
			// double pF_UnL = -0.5*(+PLP_UnL)*(pR-pL) + Mcy*(PLP_UnL)*0.5*(pL+pR);
			// double pF_UnR = -0.5*(-PRM_UnR)*(pR-pL) + Mcy*(PRM_UnR)*0.5*(pL+pR);
			double pF_UnL = 0.0;
			double pF_UnR = 0.0;
			
			// double pF_TL = KLR*(PLP+PRM-1.0)*0.5*drhodTL*chat;
			// double pF_TR = KLR*(PLP+PRM-1.0)*0.5*drhodTR*chat; 
			double pF_TL = 0.0;
			double pF_TR = 0.0;
			// double dcdTL = drhodTL*dcdphi_coeffL;
			// double PLP_TL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdTL);
			// pF_TL -= 0.5*(pR-pL)*PLP_TL;
			// pF_TL += KLR*rhohat*chat*PLP_TL;
			// pF_TL += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdTL;
			// double dcdTR = drhodTR*dcdphi_coeffR;
			// double PRM_TR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdTR);
			// pF_TR -= 0.5*(pR-pL)*PRM_TR;
			// pF_TR += KLR*rhohat*chat*PRM_TR;
			// pF_TR += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdTR;
			
			double pF_YL[nSp];
			double pF_YR[nSp];
			for(int i=0; i<nSp-1; ++i){
				// pF_YL[i] = KLR*(PLP+PRM-1.0)*0.5*drhodYL[i]*chat;
				// pF_YR[i] = KLR*(PLP+PRM-1.0)*0.5*drhodYR[i]*chat;
				pF_YL[i] = 0.0;
				pF_YR[i] = 0.0;
				// double dcdYL = drhodYL[i]*dcdphi_coeffL;
				// double PLP_YL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdYL);
				// pF_YL[i] -= 0.5*(pR-pL)*PLP_YL;
				// pF_YL[i] += KLR*rhohat*chat*PLP_YL;
				// pF_YL[i] += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdYL;
				// double dcdYR = drhodYR[i]*dcdphi_coeffR;
				// double PRM_YR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdYR);
				// pF_YR[i] -= 0.5*(pR-pL)*PRM_YR;
				// pF_YR[i] += KLR*rhohat*chat*PRM_YR;
				// pF_YR[i] += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdYR;
			}
			
            
			double dissEg_pL,dissEg_pR,dissEg_UnL,dissEg_UnR,dissEg_TL,dissEg_TR;
			double dissEg_YL[nSp];
			double dissEg_YR[nSp];
			dissEg_pL = 0.0; dissEg_pR = 0.0;
			dissEg_UnL = 0.0; dissEg_UnR = 0.0;
			dissEg_TL = 0.0; dissEg_TR = 0.0;
			for(int i=0; i<nSp-1; ++i){
				dissEg_YL[i] = 0.0;
				dissEg_YR[i] = 0.0;
			}
            // =================================================================
			
			
            
            
            
            

            // // // =================================================================
			// // // Upwind
            // // double Unhat = 0.5*(UnL+UnR);
            
            // // Unhat -= dAlpha * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR)*(pR-pL);
			// // Unhat += dt*0.5*(dpdxL/rhoL+dpdxR/rhoR)*dAlpha*nLR[0];
			// // Unhat += dt*0.5*(dpdyL/rhoL+dpdyR/rhoR)*dAlpha*nLR[1];
			// // Unhat += dt*0.5*(dpdzL/rhoL+dpdzR/rhoR)*dAlpha*nLR[2];
            
			// // double mdot = 0.0;
            // // if(Unhat>=0.0){
                // // mdot = rhoFL*Unhat;
            // // }
            // // else{
                // // mdot = rhoFR*Unhat;
            // // }
            
            
			
			// // // // in pressure-based
			// // // mdot -= dAlpha * dt/dLR*(pR-pL);
			// // // mdot += dt*0.5*(dpdxL+dpdxR)*dAlpha*nLR[0];
			// // // mdot += dt*0.5*(dpdyL+dpdyR)*dAlpha*nLR[1];
			// // // mdot += dt*0.5*(dpdzL+dpdzR)*dAlpha*nLR[2];
            
			// // double pF = 0.5*(pL+pR);
			
			
			// // double weiL = (mdot>=0.0 ? 1.0 : 0.0); double weiR = 1.0-weiL;
			
			// // double rhoF = weiL*rhoFL + weiR*rhoFR;
			// // double uF = weiL*uFL + weiR*uFR;
			// // double vF = weiL*vFL + weiR*vFR;
			// // double wF = weiL*wFL + weiR*wFR;
			// // double HtF = weiL*HtFL + weiR*HtFR;
			// // double YF[nSp];
			// // for(int i=0; i<nSp-1; ++i){
				// // YF[i] = weiL*YFL[i] + weiR*YFR[i];
			// // }
			
			
			// // double fluxB[nEq];
			// // // 컨벡티브 B
			// // fluxB[0] = -( mdot )*area;
			// // fluxB[1] = -( mdot*uF + pF*nvec[0] )*area;
			// // fluxB[2] = -( mdot*vF + pF*nvec[1] )*area;
			// // fluxB[3] = -( mdot*wF + pF*nvec[2] )*area;
			// // fluxB[4] = -( mdot*HtF )*area;
			// // for(int i=0; i<nSp-1; ++i){
				// // fluxB[5+i] = -( mdot*YF[i] )*area;
			// // }
			
			
			// // double mdot_pL = drhodpL*Unhat;
			// // double mdot_pR = drhodpR*Unhat;
            
            // // // mdot_pL += dAlpha * dt/dLR;
            // // // mdot_pR -= dAlpha * dt/dLR;
            // // if(mdot>=0.0){
                // // mdot_pL += rhoFL * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR);
                // // mdot_pR -= rhoFL * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR);
            // // }
            // // else{
                // // mdot_pL += rhoFR * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR);
                // // mdot_pR -= rhoFR * dt/dLR*0.5*(1.0/rhoL+1.0/rhoR);
            // // }
            
			// // double mdot_UnL = rhoF*0.5;
			// // double mdot_UnR = rhoF*0.5;
            
			// // double mdot_TL = drhodTL*Unhat;
			// // double mdot_TR = drhodTR*Unhat;
			
			// // double mdot_YL[nSp];
			// // double mdot_YR[nSp];
			// // for(int i=0; i<nSp-1; ++i){
                // // mdot_YL[i] = drhodYL[i]*Unhat;
                // // mdot_YR[i] = drhodYR[i]*Unhat;
			// // }
			
			// // double pF_pL = 0.5;
			// // double pF_pR = 0.5;
			
			// // double pF_UnL = 0.0;
			// // double pF_UnR = 0.0;
			// // double pF_TL = 0.0;
			// // double pF_TR = 0.0;
			// // double pF_YL[nSp];
			// // double pF_YR[nSp];
			// // for(int i=0; i<nSp-1; ++i){
				// // pF_YL[i] = 0.0;
				// // pF_YR[i] = 0.0;
			// // }
			
            
			// // double dissEg_pL,dissEg_pR,dissEg_UnL,dissEg_UnR,dissEg_TL,dissEg_TR;
			// // double dissEg_YL[nSp];
			// // double dissEg_YR[nSp];
			// // dissEg_pL = 0.0; dissEg_pR = 0.0;
			// // dissEg_UnL = 0.0; dissEg_UnR = 0.0;
			// // dissEg_TL = 0.0; dissEg_TR = 0.0;
			// // for(int i=0; i<nSp-1; ++i){
				// // dissEg_YL[i] = 0.0;
				// // dissEg_YR[i] = 0.0;
			// // }
            // // // =================================================================
			
			
            
            
            
            
            
            
            
            
            

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
            
            
            
            
            
            
            
            double weipL = weiL*gamPhiL[0] + weiR*(1.0-gamPhiR[0]);
            double weipR = weiL*(1.0-gamPhiL[0]) + weiR*gamPhiR[0];
            double weiuL = weiL*gamPhiL[1] + weiR*(1.0-gamPhiR[1]);
            double weiuR = weiL*(1.0-gamPhiL[1]) + weiR*gamPhiR[1];
            double weivL = weiL*gamPhiL[2] + weiR*(1.0-gamPhiR[2]);
            double weivR = weiL*(1.0-gamPhiL[2]) + weiR*gamPhiR[2];
            double weiwL = weiL*gamPhiL[3] + weiR*(1.0-gamPhiR[3]);
            double weiwR = weiL*(1.0-gamPhiL[3]) + weiR*gamPhiR[3];
            double weiTL = weiL*gamPhiL[4] + weiR*(1.0-gamPhiR[4]);
            double weiTR = weiL*(1.0-gamPhiL[4]) + weiR*gamPhiR[4];
            
            
			
			int iter = 0;
			
			fluxA_LL[iter] += CN_coeff * (mdot_pL)*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_pR)*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_pR)*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_pL)*area; 
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_TL)*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_TR)*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_TR)*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_TL)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i])*area;
				++iter;
			}
			

			
			
			fluxA_LL[iter] += CN_coeff * (mdot_pL*uF +pF_pL*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_pR*uF +pF_pR*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_pR*uF +pF_pR*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_pL*uF +pF_pL*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0]*uF+mdot*weiuL +pF_UnL*nvec[0]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0]*uF+mdot*weiuR +pF_UnR*nvec[0]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0]*uF+mdot*weiuR +pF_UnR*nvec[0]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0]*uF+mdot*weiuL +pF_UnL*nvec[0]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1]*uF +pF_UnL*nvec[1]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1]*uF +pF_UnR*nvec[1]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1]*uF +pF_UnR*nvec[1]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1]*uF +pF_UnL*nvec[1]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2]*uF +pF_UnL*nvec[2]*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2]*uF +pF_UnR*nvec[2]*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2]*uF +pF_UnR*nvec[2]*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2]*uF +pF_UnL*nvec[2]*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_TL*uF + pF_TL*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_TR*uF + pF_TR*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_TR*uF + pF_TR*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_TL*uF + pF_TL*nvec[0])*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*uF + pF_YL[i]*nvec[0])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*uF + pF_YR[i]*nvec[0])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*uF + pF_YR[i]*nvec[0])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*uF + pF_YL[i]*nvec[0])*area;
				++iter;
			}
			
			
			
			
			
			
			fluxA_LL[iter] += CN_coeff * (mdot_pL*vF +pF_pL*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_pR*vF +pF_pR*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_pR*vF +pF_pR*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_pL*vF +pF_pL*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0]*vF +pF_UnL*nvec[0]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0]*vF +pF_UnR*nvec[0]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0]*vF +pF_UnR*nvec[0]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0]*vF +pF_UnL*nvec[0]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1]*vF+mdot*weivL +pF_UnL*nvec[1]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1]*vF+mdot*weivR +pF_UnR*nvec[1]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1]*vF+mdot*weivR +pF_UnR*nvec[1]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1]*vF+mdot*weivL +pF_UnL*nvec[1]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2]*vF +pF_UnL*nvec[2]*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2]*vF +pF_UnR*nvec[2]*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2]*vF +pF_UnR*nvec[2]*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2]*vF +pF_UnL*nvec[2]*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_TL*vF + pF_TL*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_TR*vF + pF_TR*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_TR*vF + pF_TR*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_TL*vF + pF_TL*nvec[1])*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*vF + pF_YL[i]*nvec[1])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*vF + pF_YR[i]*nvec[1])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*vF + pF_YR[i]*nvec[1])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*vF + pF_YL[i]*nvec[1])*area;
				++iter;
			}
			
			
			
			
			
			
			fluxA_LL[iter] += CN_coeff * (mdot_pL*wF +pF_pL*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_pR*wF +pF_pR*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_pR*wF +pF_pR*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_pL*wF +pF_pL*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0]*wF +pF_UnL*nvec[0]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0]*wF +pF_UnR*nvec[0]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0]*wF +pF_UnR*nvec[0]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0]*wF +pF_UnL*nvec[0]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1]*wF +pF_UnL*nvec[1]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1]*wF +pF_UnR*nvec[1]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1]*wF +pF_UnR*nvec[1]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1]*wF +pF_UnL*nvec[1]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2]*wF+mdot*weiwL +pF_UnL*nvec[2]*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2]*wF+mdot*weiwR +pF_UnR*nvec[2]*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2]*wF+mdot*weiwR +pF_UnR*nvec[2]*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2]*wF+mdot*weiwL +pF_UnL*nvec[2]*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_TL*wF + pF_TL*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_TR*wF + pF_TR*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_TR*wF + pF_TR*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_TL*wF + pF_TL*nvec[2])*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*wF + pF_YL[i]*nvec[2])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*wF + pF_YR[i]*nvec[2])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*wF + pF_YR[i]*nvec[2])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*wF + pF_YL[i]*nvec[2])*area;
				++iter;
			}
			
			
			
			fluxA_LL[iter] += CN_coeff * (mdot_pL*HtF + mdot*weiL*dHtdpL  + dissEg_pL)*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_pR*HtF + mdot*weiR*dHtdpR  + dissEg_pR)*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_pR*HtF + mdot*weiR*dHtdpR  + dissEg_pR)*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_pL*HtF + mdot*weiL*dHtdpL  + dissEg_pL)*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0]*HtF + mdot*weiuL*uL  + dissEg_UnL*nvec[0])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0]*HtF + mdot*weiuR*uR  + dissEg_UnR*nvec[0])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0]*HtF + mdot*weiuR*uR  + dissEg_UnR*nvec[0])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0]*HtF + mdot*weiuL*uL  + dissEg_UnL*nvec[0])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1]*HtF + mdot*weivL*vL  + dissEg_UnL*nvec[1])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1]*HtF + mdot*weivR*vR  + dissEg_UnR*nvec[1])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1]*HtF + mdot*weivR*vR  + dissEg_UnR*nvec[1])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1]*HtF + mdot*weivL*vL  + dissEg_UnL*nvec[1])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2]*HtF + mdot*weiwL*wL  + dissEg_UnL*nvec[2])*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2]*HtF + mdot*weiwR*wR  + dissEg_UnR*nvec[2])*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2]*HtF + mdot*weiwR*wR  + dissEg_UnR*nvec[2])*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2]*HtF + mdot*weiwL*wL  + dissEg_UnL*nvec[2])*area;
			++iter;
			
			fluxA_LL[iter] += CN_coeff * (mdot_TL*HtF + mdot*weiL*dHtdTL  + dissEg_TL)*area; 
			fluxA_LR[iter] += CN_coeff * (mdot_TR*HtF + mdot*weiR*dHtdTR  + dissEg_TR)*area;
			fluxA_RR[iter] -= CN_coeff * (mdot_TR*HtF + mdot*weiR*dHtdTR  + dissEg_TR)*area; 
			fluxA_RL[iter] -= CN_coeff * (mdot_TL*HtF + mdot*weiL*dHtdTL  + dissEg_TL)*area;
			++iter;
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*HtF + mdot*weiL*dHtdYL[i]  + dissEg_YL[i])*area; 
				fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*HtF + mdot*weiR*dHtdYR[i]  + dissEg_YR[i])*area;
				fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*HtF + mdot*weiR*dHtdYR[i]  + dissEg_YR[i])*area; 
				fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*HtF + mdot*weiL*dHtdYL[i]  + dissEg_YL[i])*area;
				++iter;
			}
			
			
			
			for(int i=0; i<nSp-1; ++i){
				fluxA_LL[iter] += CN_coeff * (mdot_pL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (mdot_pR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (mdot_pR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (mdot_pL*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[0]*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[0]*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[0]*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[0]*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[1]*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[1]*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[1]*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[1]*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (mdot_UnL*nvec[2]*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (mdot_UnR*nvec[2]*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (mdot_UnR*nvec[2]*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (mdot_UnL*nvec[2]*YF[i])*area;
				++iter;
				
				fluxA_LL[iter] += CN_coeff * (mdot_TL*YF[i])*area; 
				fluxA_LR[iter] += CN_coeff * (mdot_TR*YF[i])*area;
				fluxA_RR[iter] -= CN_coeff * (mdot_TR*YF[i])*area; 
				fluxA_RL[iter] -= CN_coeff * (mdot_TL*YF[i])*area;
				++iter;
				
				for(int j=0; j<nSp-1; ++j){
					fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*YF[i])*area; 
					fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*YF[i])*area;
					fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*YF[i])*area; 
					fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*YF[i])*area;
					if(i==j){
						// fluxA_LL[iter] += CN_coeff_Y * (mdot*weiL)*area; 
						// fluxA_LR[iter] += CN_coeff_Y * (mdot*weiR)*area;
						// fluxA_RR[iter] -= CN_coeff_Y * (mdot*weiR)*area; 
						// fluxA_RL[iter] -= CN_coeff_Y * (mdot*weiL)*area;
                        
                        // double gamYL = faces[id_gammaYL[j]];
                        // double gamYR = faces[id_gammaYR[j]];
                        
                        double weiYL = weiL*gamPhiL[5+i] + weiR*(1.0-gamPhiR[5+i]);
                        double weiYR = weiL*(1.0-gamPhiL[5+i]) + weiR*gamPhiR[5+i];
                        
						fluxA_LL[iter] += CN_coeff_Y * (mdot*weiYL)*area; 
						fluxA_LR[iter] += CN_coeff_Y * (mdot*weiYR)*area;
						fluxA_RR[iter] -= CN_coeff_Y * (mdot*weiYR)*area; 
						fluxA_RL[iter] -= CN_coeff_Y * (mdot*weiYL)*area;
					}
					
					
					// //=========================
					// double tmp_weiL = 1.0; double tmp_weiR = 0.0;
					// if(YF[i] < YL[i]+1.e-50 && YF[i] > YL[i]-1.e-50){
						// tmp_weiL = 1.0; tmp_weiR = 0.0;
					// }
					// else{
						// tmp_weiL = 0.0; tmp_weiR = 1.0;
					// }
					// fluxA_LL[iter] += CN_coeff_Y * (mdot_YL[i]*YF[i])*area; 
					// fluxA_LR[iter] += CN_coeff_Y * (mdot_YR[i]*YF[i])*area;
					// fluxA_RR[iter] -= CN_coeff_Y * (mdot_YR[i]*YF[i])*area; 
					// fluxA_RL[iter] -= CN_coeff_Y * (mdot_YL[i]*YF[i])*area;
					// if(i==j){
						// fluxA_LL[iter] += CN_coeff_Y * (mdot*tmp_weiL)*area; 
						// fluxA_LR[iter] += CN_coeff_Y * (mdot*tmp_weiR)*area;
						// fluxA_RR[iter] -= CN_coeff_Y * (mdot*tmp_weiR)*area; 
						// fluxA_RL[iter] -= CN_coeff_Y * (mdot*tmp_weiL)*area;
					// }
					// //=========================
					
					
					
					++iter;
				}
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
			bool_zeroGradient,bool_inletOutlet,
            id_dYdx,id_dYdy,id_dYdz](
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
                
                
                
				double dp_coeff_Un = dAlpha * dtrho/dLR;
                
				// double dp_coeff_thm = 0.5*chi/rhohat/chat;
				double dp_coeff_thm = 0.0;
				
				
				
			
				
				
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









			// double drho_dUn_dpL = weiL*drhodpL*UnF+rhoF*dp_coeff_Un +rhoF*dp_coeff_thm;
			// double drho_dUn_dpR = weiR*drhodpR*UnF-rhoF*dp_coeff_Un -rhoF*dp_coeff_thm;
			// double rhoF_dUn_duL = rhoF*(WUL*nvec[0]) +rhoF*diffDPU*nvec[0];
			// double rhoF_dUn_duR = rhoF*(WUR*nvec[0]) -rhoF*diffDPU*nvec[0];
			// double rhoF_dUn_dvL = rhoF*(WUL*nvec[1]) +rhoF*diffDPU*nvec[1];
			// double rhoF_dUn_dvR = rhoF*(WUR*nvec[1]) -rhoF*diffDPU*nvec[1];
			// double rhoF_dUn_dwL = rhoF*(WUL*nvec[2]) +rhoF*diffDPU*nvec[2];
			// double rhoF_dUn_dwR = rhoF*(WUR*nvec[2]) -rhoF*diffDPU*nvec[2];
			
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
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL +WpUL*nvec[0]*nvec[0])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR +WpUR*nvec[0]*nvec[0])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*uF+rhoF*UnF*weiR +WpUR*nvec[0]*nvec[0])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*uF+rhoF*UnF*weiL +WpUL*nvec[0]*nvec[0])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*uF +WpUL*nvec[1]*nvec[0])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*uF +WpUR*nvec[1]*nvec[0])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*uF +WpUR*nvec[1]*nvec[0])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*uF +WpUL*nvec[1]*nvec[0])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*uF +WpUL*nvec[2]*nvec[0])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*uF +WpUR*nvec[2]*nvec[0])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*uF +WpUR*nvec[2]*nvec[0])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*uF +WpUL*nvec[2]*nvec[0])*area;
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
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*vF +WpUL*nvec[0]*nvec[1])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*vF +WpUR*nvec[0]*nvec[1])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*vF +WpUR*nvec[0]*nvec[1])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*vF +WpUL*nvec[0]*nvec[1])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL +WpUL*nvec[1]*nvec[1])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR +WpUR*nvec[1]*nvec[1])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*vF+rhoF*UnF*weiR +WpUR*nvec[1]*nvec[1])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*vF+rhoF*UnF*weiL +WpUL*nvec[1]*nvec[1])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF +WpUL*nvec[2]*nvec[1])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF +WpUR*nvec[2]*nvec[1])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF +WpUR*nvec[2]*nvec[1])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF +WpUL*nvec[2]*nvec[1])*area;
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
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_duL*wF +WpUL*nvec[0]*nvec[2])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_duR*wF +WpUR*nvec[0]*nvec[2])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_duR*wF +WpUR*nvec[0]*nvec[2])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_duL*wF +WpUL*nvec[0]*nvec[2])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dvL*wF +WpUL*nvec[1]*nvec[2])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dvR*wF +WpUR*nvec[1]*nvec[2])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dvR*wF +WpUR*nvec[1]*nvec[2])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dvL*wF +WpUL*nvec[1]*nvec[2])*area;
			// ++iter;
			
			// fluxA_LL[iter] += CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL +WpUL*nvec[2]*nvec[2])*area; 
			// fluxA_LR[iter] += CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR +WpUR*nvec[2]*nvec[2])*area;
			// fluxA_RR[iter] -= CN_coeff * (rhoF_dUn_dwR*vF+rhoF*UnF*weiR +WpUR*nvec[2]*nvec[2])*area; 
			// fluxA_RL[iter] -= CN_coeff * (rhoF_dUn_dwL*vF+rhoF*UnF*weiL +WpUL*nvec[2]*nvec[2])*area;
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
			
			
			
			// double visc_coeff = dAlpha*muF/dLR*area;
			
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
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			// double sign_UnL = (UnL>=0.0 ? 1.0 : -1.0);
			// double sign_UnR = (UnR>=0.0 ? 1.0 : -1.0);
			// double D_rho_dL = g_c*(abs(UnL)-Mbar)/(rhoL+rhoR);
			// double D_rho_dR = g_c*(abs(UnR)-Mbar)/(rhoL+rhoR);
			
			// double mdot_pL = 0.5*drhodpL*(UPL + D_rho_dL) +
								// 0.5*phi_c/chat + dAlpha*dt/dLR;
			// double mdot_pR = 0.5*drhodpR*(UMR - D_rho_dR) -
								// 0.5*phi_c/chat - dAlpha*dt/dLR;
			// mdot -= dAlpha * dt/dLR*(pR-pL);
			
			// double mdot_UnL = 0.5*rhoL*(1.0+(1.0-g_c)*sign_UnL+g_c*rhoL/(rhoL+rhoR)*sign_UnL);
			// double mdot_UnR = 0.5*rhoR*(1.0+(1.0-g_c)*sign_UnR-g_c*rhoR/(rhoL+rhoR)*sign_UnR);
			
			// double mdot_TL = 0.5*drhodTL*(UPL + D_rho_dL);
			// double mdot_TR = 0.5*drhodTR*(UMR - D_rho_dR);
			
			// double mdot_YL[nSp];
			// double mdot_YR[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// mdot_YL[i] = 0.5*drhodYL[i]*(UPL + D_rho_dL);
				// mdot_YR[i] = 0.5*drhodYR[i]*(UMR - D_rho_dR);
			// }
			
			
			
			// double dcdphi_coeffL = (cL/2.0/rhoL*(1.0-cL*cL*(drhodpL-drhodTL*dHtdpL/dHtdTL)));
			// double dcdphi_coeffR = (cR/2.0/rhoR*(1.0-cR*cR*(drhodpR-drhodTR*dHtdpR/dHtdTR)));
			
			// double pF_pL = 0.5 - 0.5*(PLP-PRM)*(-1.0) + KLR*(PLP+PRM-1.0)*0.5*drhodpL*chat;
			// double pF_pR = 0.5 - 0.5*(PLP-PRM)*(+1.0) + KLR*(PLP+PRM-1.0)*0.5*drhodpR*chat;
			// double dcdpL = drhodpL*dcdphi_coeffL;
			// double PLP_pL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdpL);
			// pF_pL -= 0.5*(pR-pL)*PLP_pL;
			// pF_pL += KLR*rhohat*chat*PLP_pL;
			// pF_pL += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdpL;
			// double dcdpR = drhodpR*dcdphi_coeffR;
			// double PRM_pR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdpR);
			// pF_pR -= 0.5*(pR-pL)*PRM_pR;
			// pF_pR += KLR*rhohat*chat*PRM_pR;
			// pF_pR += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdpR;
			
			// double PLP_UnL = (-0.75*ML*ML+0.75)/chat;
			// double pF_UnL = -0.5*(PLP_UnL)*(pR-pL) + KLR*(PLP_UnL)*rhohat*chat;
			// double PLP_UnR = (0.75*MR*MR-0.75)/chat;
			// double pF_UnR = -0.5*(PLP_UnR)*(pR-pL) + KLR*(PLP_UnR)*rhohat*chat;
			
			// double pF_TL = KLR*(PLP+PRM-1.0)*0.5*drhodTL*chat;
			// double pF_TR = KLR*(PLP+PRM-1.0)*0.5*drhodTR*chat;
			// double dcdTL = drhodTL*dcdphi_coeffL;
			// double PLP_TL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdTL);
			// pF_TL -= 0.5*(pR-pL)*PLP_TL;
			// pF_TL += KLR*rhohat*chat*PLP_TL;
			// pF_TL += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdTL;
			// double dcdTR = drhodTR*dcdphi_coeffR;
			// double PRM_TR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdTR);
			// pF_TR -= 0.5*(pR-pL)*PRM_TR;
			// pF_TR += KLR*rhohat*chat*PRM_TR;
			// pF_TR += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdTR;
			
			// double pF_YL[nSp];
			// double pF_YR[nSp];
			// for(int i=0; i<nSp-1; ++i){
				// pF_YL[i] = KLR*(PLP+PRM-1.0)*0.5*drhodYL[i]*chat;
				// pF_YR[i] = KLR*(PLP+PRM-1.0)*0.5*drhodYR[i]*chat;
				// double dcdYL = drhodYL[i]*dcdphi_coeffL;
				// double PLP_YL = (-0.75*ML*ML+0.75)*(-UnL/chat/chat*0.5*dcdYL);
				// pF_YL[i] -= 0.5*(pR-pL)*PLP_YL;
				// pF_YL[i] += KLR*rhohat*chat*PLP_YL;
				// pF_YL[i] += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdYL;
				// double dcdYR = drhodYR[i]*dcdphi_coeffR;
				// double PRM_YR = (-0.75*MR*MR+0.75)*(-UnR/chat/chat*0.5*dcdYR);
				// pF_YR[i] -= 0.5*(pR-pL)*PRM_YR;
				// pF_YR[i] += KLR*rhohat*chat*PRM_YR;
				// pF_YR[i] += KLR*(PLP+PRM-1.0)*rhohat*0.5*dcdYR;
			// }
			
			
            
            
            
            
            
            
            
            
            
            
            
            
			// double dissEg = weiL*SL*(pStar-pL)/(SL-UnL) + weiR*SR*(pStar-pR)/(SR-UnR);
			
			// double mdot_pL,mdot_pR,mdot_UnL,mdot_UnR,mdot_TL,mdot_TR;
			// double mdot_YL[nSp];
			// double mdot_YR[nSp];
			// double pF_pL,pF_pR,pF_UnL,pF_UnR,pF_TL,pF_TR;
			// double pF_YL[nSp];
			// double pF_YR[nSp];
			// double dissEg_pL,dissEg_pR,dissEg_UnL,dissEg_UnR,dissEg_TL,dissEg_TR;
			// double dissEg_YL[nSp];
			// double dissEg_YR[nSp];
            
            // if(SL>=0.0){
            // }
            // else if(0.0>=SR){
            // }
            // else{
                
                // // UnL = UnLas;
                // // UnR = UnRas;

                // // // HLLC
                // // double SM = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR)/
                            // // (rhoR*(SR-UnR)-rhoL*(SL-UnL));
                            
            // }
            

			
            
            
            
			
			// if(SM>0.0){
				
				// double SM_up = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR);
				// double SM_down = (rhoR*(SR-UnR)-rhoL*(SL-UnL));
				// double SL_pL = 0.0;
				// double SM_pL = (-drhodpL*UnL*(SL-UnL)+1.0)/SM_down - SM_up/SM_down/SM_down*(-drhodpL*(SL-UnL));
				// double SM_pR = (drhodpR*UnR*(SR-UnR)-1.0)/SM_down - SM_up/SM_down/SM_down*(drhodpR*(SR-UnR));
				// double rhoLStar_pL = - (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(-SM_pL) +
									 // (SL-UnL)/(SL-SM)*drhodpL;
				// double rhoLStar_pR = (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(SM_pR);
				// double SM_UnL = (-rhoL*(SL-UnL) - rhoL*UnL*(-1.0))/SM_down - SM_up/SM_down/SM_down*(-rhoL*(-1.0));
				// double SM_UnR = (rhoR*(SR-UnR) + rhoR*UnR*(-1.0))/SM_down - SM_up/SM_down/SM_down*(rhoR*(-1.0));
				// double rhoLStar_UnL = -1.0/(SL-SM)*rhoL - (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(-SM_UnL);
				// double rhoLStar_UnR = - (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(-SM_UnR);
				// double SM_TL = (-drhodTL*UnL*(SL-UnL))/SM_down - SM_up/SM_down/SM_down*(-drhodTL*(SL-UnL));
				// double SM_TR = (drhodTR*UnR*(SR-UnR))/SM_down - SM_up/SM_down/SM_down*(drhodTR*(SR-UnR));
				// double rhoLStar_TL = - (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(-SM_TL) +
									 // (SL-UnL)/(SL-SM)*drhodTL;
				// double rhoLStar_TR = (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(SM_TR);
				
				// mdot_pL = drhodpL*UnL + SL*(rhoLStar_pL-drhodpL);
				// mdot_pR = SL*rhoLStar_pR;
				// double UnF_pL = mdot_pL*(weiL/rhoL+weiR/rhoR) + mdot*drhodpL*(-weiL/rhoL/rhoL);
				// double UnF_pR = mdot_pR*(weiL/rhoL+weiR/rhoR) + mdot*drhodpR*(-weiR/rhoR/rhoR);
				// mdot_pL = weiL*drhodpL*UnF + rhoF*UnF_pL;
				// mdot_pR = weiR*drhodpR*UnF + rhoF*UnF_pR;
				
				// mdot_UnL = rhoL*1.0 + SL*(rhoLStar_UnL);
				// mdot_UnR = SL*(rhoLStar_UnR);
				// double UnF_UnL = mdot_UnL*(weiL/rhoL+weiR/rhoR);
				// double UnF_UnR = mdot_UnR*(weiL/rhoL+weiR/rhoR);
				// mdot_UnL = rhoF*UnF_UnL;
				// mdot_UnR = rhoF*UnF_UnR;
				
				// mdot_TL = drhodTL*UnL + SL*(rhoLStar_TL-drhodTL);
				// mdot_TR = SL*(rhoLStar_TR);
				// double UnF_TL = mdot_TL*(weiL/rhoL+weiR/rhoR) + mdot*drhodTL*(-weiL/rhoL/rhoL);
				// double UnF_TR = mdot_TR*(weiL/rhoL+weiR/rhoR) + mdot*drhodTR*(-weiR/rhoR/rhoR);
				// mdot_TL = weiL*drhodTL*UnF + rhoF*UnF_TL;
				// mdot_TR = weiR*drhodTR*UnF + rhoF*UnF_TR;
				
				// double SM_YL[nSp];
				// double SM_YR[nSp];
				// double UnF_YL[nSp];
				// double UnF_YR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// SM_YL[i] = (-drhodYL[i]*UnL*(SL-UnL))/SM_down - SM_up/SM_down/SM_down*(-drhodYL[i]*(SL-UnL));
					// SM_YR[i] = (drhodYR[i]*UnR*(SR-UnR))/SM_down - SM_up/SM_down/SM_down*(drhodYR[i]*(SR-UnR));
					// double rhoLStar_YL = - (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(-SM_YL[i]) +
										 // (SL-UnL)/(SL-SM)*drhodYL[i];
					// double rhoLStar_YR = (SL-UnL)/(SL-SM)/(SL-SM)*rhoL*(SM_YR[i]);
					// mdot_YL[i] = drhodYL[i]*UnL + SL*(rhoLStar_YL-drhodYL[i]);
					// mdot_YR[i] = SL*(rhoLStar_YR);
					// UnF_YL[i] = mdot_YL[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYL[i]*(-weiL/rhoL/rhoL);
					// UnF_YR[i] = mdot_YR[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYR[i]*(-weiR/rhoR/rhoR);
					// mdot_YL[i] = weiL*drhodYL[i]*UnF + rhoF*UnF_YL[i];
					// mdot_YR[i] = weiR*drhodYR[i]*UnF + rhoF*UnF_YR[i];
				// }
				
			
				// double pStar_pL = 0.5*( 1.0 + drhodpL*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_pL) +
								  // rhoR*(UnR-SR)*(-SM_pL) );
				// double pStar_pR = 0.5*( 1.0 + rhoL*(UnL-SL)*(-SM_pR) + drhodpR*(UnR-SR)*(UnR-SM) +
								  // rhoR*(UnR-SR)*(-SM_pR) );
				// double pStar_UnL = 0.5*(rhoL*(1.0)*(UnL-SM)+rhoL*(UnL-SL)*(1.0-SM_UnL)+rhoR*(UnR-SR)*(-SM_UnL));
				// double pStar_UnR = 0.5*(rhoR*(1.0)*(UnR-SM)+rhoR*(UnR-SR)*(1.0-SM_UnR)+rhoL*(UnL-SL)*(-SM_UnR));
				// double pStar_TL = 0.5*(drhodTL*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_TL) +
								  // rhoR*(UnR-SR)*(-SM_TL) );
				// double pStar_TR = 0.5*(rhoL*(UnL-SL)*(-SM_TR) + drhodTR*(UnR-SR)*(UnR-SM) +
								  // rhoR*(UnR-SR)*(-SM_TR) );
								  
				// pF_pL = pStar_pL + SM_pL/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_pL +
						// SM/(SM-SL)*(1.0-pStar_pL);
				// pF_pR = pStar_pR + SM_pR/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_pR +
						// SM/(SM-SL)*(-pStar_pR);
				
				// pF_UnL = pStar_UnL + SM_UnL/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(SM_UnL) + SM/(SM-SL)*(-pStar_UnL);
				// pF_UnR = pStar_UnR + SM_UnR/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(SM_UnR) + SM/(SM-SL)*(-pStar_UnR);
				
				// pF_TL = pStar_TL + SM_TL/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_TL +
						// SM/(SM-SL)*(-pStar_TL);
				// pF_TR = pStar_TR + SM_TR/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_TR +
						// SM/(SM-SL)*(-pStar_TR);
				
				// double pStar_YL[nSp];
				// double pStar_YR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// pStar_YL[i] = 0.5*(drhodYL[i]*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_YL[i]) +
									  // rhoR*(UnR-SR)*(-SM_YL[i]) );
					// pStar_YR[i] = 0.5*(rhoL*(UnL-SL)*(-SM_YR[i]) + drhodYR[i]*(UnR-SR)*(UnR-SM) +
									  // rhoR*(UnR-SR)*(-SM_YR[i]) );
					// pF_YL[i] = pStar_YL[i] + SM_YL[i]/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_YL[i] +
								// SM/(SM-SL)*(-pStar_YL[i]);
					// pF_YR[i] = pStar_YR[i] + SM_YR[i]/(SM-SL)*(pL-pStar) - SM/(SM-SL)/(SM-SL)*(pL-pStar)*SM_YR[i] +
								// SM/(SM-SL)*(-pStar_YR[i]);
				// }
				
				
				// dissEg_pL = weiL*SL*(pStar_pL-1.0)/(SL-UnL) + weiR*SR*(pStar_pL)/(SR-UnR);
				// dissEg_pR = weiL*SL*(pStar_pR)/(SL-UnL) + weiR*SR*(pStar_pR-1.0)/(SR-UnR);
				// dissEg_TL = weiL*SL*(pStar_TL)/(SL-UnL) + weiR*SR*(pStar_TL)/(SR-UnR);
				// dissEg_TR = weiL*SL*(pStar_TR)/(SL-UnL) + weiR*SR*(pStar_TR)/(SR-UnR);
				// for(int i=0; i<nSp-1; ++i){
					// dissEg_YL[i] = weiL*SL*(pStar_YL[i])/(SL-UnL) + weiR*SR*(pStar_YL[i])/(SR-UnR);
					// dissEg_YR[i] = weiL*SL*(pStar_YR[i])/(SL-UnL) + weiR*SR*(pStar_YR[i])/(SR-UnR);
				// }
				// dissEg_UnL = weiL*SL*(pStar_UnL)/(SL-UnL) - weiL*SL*(pStar-pL)/(SL-UnL)/(SL-UnL)*(-1.0) + 
							 // weiR*SR*(pStar_UnL)/(SR-UnR);
				// dissEg_UnR = weiL*SL*(pStar_UnR)/(SL-UnL) - weiR*SR*(pStar-pR)/(SR-UnR)/(SR-UnR)*(-1.0) + 
							 // weiR*SR*(pStar_UnR)/(SR-UnR);
				
			// }
			// else{
				
				// double SM_up = (rhoR*UnR*(SR-UnR)-rhoL*UnL*(SL-UnL)+pL-pR);
				// double SM_down = (rhoR*(SR-UnR)-rhoL*(SL-UnL));
				// double SM_pL = (-drhodpL*UnL*(SL-UnL)+1.0)/SM_down - SM_up/SM_down/SM_down*(-drhodpL*(SL-UnL));
				// double SM_pR = (drhodpR*UnR*(SR-UnR)-1.0)/SM_down - SM_up/SM_down/SM_down*(drhodpR*(SR-UnR));
				// double SM_UnL = (-rhoL*(SL-UnL) - rhoL*UnL*(-1.0))/SM_down - SM_up/SM_down/SM_down*(-rhoL*(-1.0));
				// double SM_UnR = (rhoR*(SR-UnR) + rhoR*UnR*(-1.0))/SM_down - SM_up/SM_down/SM_down*(rhoR*(-1.0));
				// double SM_TL = (-drhodTL*UnL*(SL-UnL))/SM_down - SM_up/SM_down/SM_down*(-drhodTL*(SL-UnL));
				// double SM_TR = (drhodTR*UnR*(SR-UnR))/SM_down - SM_up/SM_down/SM_down*(drhodTR*(SR-UnR));
				
				// double rhoRStar_pL = - (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(-SM_pL);
				// double rhoRStar_pR = (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(SM_pR) +
									 // (SR-UnR)/(SR-SM)*drhodpR;
				// double rhoRStar_UnL = - (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(-SM_UnR);
				// double rhoRStar_UnR = -1.0/(SR-SM)*rhoR - (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(-SM_UnR);
				// double rhoRStar_TL = - (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(-SM_TL);
				// double rhoRStar_TR = (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(SM_TR) +
									 // (SR-UnR)/(SR-SM)*drhodTR;
				
				// mdot_pL = SR*rhoRStar_pL;
				// mdot_pR = drhodpR*UnR + SR*(rhoRStar_pL-drhodpR);
				// double UnF_pL = mdot_pL*(weiL/rhoL+weiR/rhoR) + mdot*drhodpL*(-weiL/rhoL/rhoL);
				// double UnF_pR = mdot_pR*(weiL/rhoL+weiR/rhoR) + mdot*drhodpR*(-weiR/rhoR/rhoR);
				// mdot_pL = weiL*drhodpL*UnF + rhoF*UnF_pL;
				// mdot_pR = weiR*drhodpR*UnF + rhoF*UnF_pR;
				
				// mdot_UnL = SR*(rhoRStar_UnL);
				// mdot_UnR = rhoR*1.0 + SR*(rhoRStar_UnR);
				// double UnF_UnL = mdot_UnL*(weiL/rhoL+weiR/rhoR);
				// double UnF_UnR = mdot_UnR*(weiL/rhoL+weiR/rhoR);
				// mdot_UnL = rhoF*UnF_UnL;
				// mdot_UnR = rhoF*UnF_UnR;
				
				// mdot_TL = SR*(rhoRStar_TL);
				// mdot_TR = drhodTR*UnR + SR*(rhoRStar_TR-drhodTR);
				// double UnF_TL = mdot_TL*(weiL/rhoL+weiR/rhoR) + mdot*drhodTL*(-weiL/rhoL/rhoL);
				// double UnF_TR = mdot_TR*(weiL/rhoL+weiR/rhoR) + mdot*drhodTR*(-weiR/rhoR/rhoR);
				// mdot_TL = weiL*drhodTL*UnF + rhoF*UnF_TL;
				// mdot_TR = weiR*drhodTR*UnF + rhoF*UnF_TR;
				
				
				// double SM_YL[nSp];
				// double SM_YR[nSp];
				// double UnF_YL[nSp];
				// double UnF_YR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// SM_YL[i] = (-drhodYL[i]*UnL*(SL-UnL))/SM_down - SM_up/SM_down/SM_down*(-drhodYL[i]*(SL-UnL));
					// SM_YR[i] = (drhodYR[i]*UnR*(SR-UnR))/SM_down - SM_up/SM_down/SM_down*(drhodYR[i]*(SR-UnR));
					// double rhoRStar_YL = - (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(-SM_YL[i]);
					// double rhoRStar_YR = (SR-UnR)/(SR-SM)/(SR-SM)*rhoR*(SM_YR[i]) +
										 // (SR-UnR)/(SR-SM)*drhodYR[i];
					// mdot_YL[i] = SR*(rhoRStar_YL);
					// mdot_YR[i] = drhodYR[i]*UnR + SR*(rhoRStar_YR-drhodYR[i]);
					// UnF_YL[i] = mdot_YL[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYL[i]*(-weiL/rhoL/rhoL);
					// UnF_YR[i] = mdot_YR[i]*(weiL/rhoL+weiR/rhoR) + mdot*drhodYR[i]*(-weiR/rhoR/rhoR);
					// mdot_YL[i] = weiL*drhodYL[i]*UnF + rhoF*UnF_YL[i];
					// mdot_YR[i] = weiR*drhodYR[i]*UnF + rhoF*UnF_YR[i];
				// }
				
			
				// double pStar_pL = 0.5*( 1.0 + drhodpL*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_pL) +
								  // rhoR*(UnR-SR)*(-SM_pL) );
				// double pStar_pR = 0.5*( 1.0 + rhoL*(UnL-SL)*(-SM_pR) + drhodpR*(UnR-SR)*(UnR-SM) +
								  // rhoR*(UnR-SR)*(-SM_pR) );
				// double pStar_UnL = 0.5*(rhoL*(1.0)*(UnL-SM)+rhoL*(UnL-SL)*(1.0-SM_UnL)+rhoR*(UnR-SR)*(-SM_UnL));
				// double pStar_UnR = 0.5*(rhoR*(1.0)*(UnR-SM)+rhoR*(UnR-SR)*(1.0-SM_UnR)+rhoL*(UnL-SL)*(-SM_UnR));
				// double pStar_TL = 0.5*(drhodTL*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_TL) +
								  // rhoR*(UnR-SR)*(-SM_TL) );
				// double pStar_TR = 0.5*(rhoL*(UnL-SL)*(-SM_TR) + drhodTR*(UnR-SR)*(UnR-SM) +
								  // rhoR*(UnR-SR)*(-SM_TR) );
								  
				// pF_pL = pStar_pL + SM_pL/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_pL +
						// SM/(SM-SR)*(-pStar_pL);
				// pF_pR = pStar_pR + SM_pR/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_pR +
						// SM/(SM-SR)*(1.0-pStar_pR);
				
				// pF_UnL = pStar_UnL + SM_UnL/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(SM_UnL) + SM/(SM-SR)*(-pStar_UnL);
				// pF_UnR = pStar_UnR + SM_UnR/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(SM_UnR) + SM/(SM-SR)*(-pStar_UnR);
				
				// pF_TL = pStar_TL + SM_TL/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_TL +
						// SM/(SM-SR)*(-pStar_TL);
				// pF_TR = pStar_TR + SM_TR/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_TR +
						// SM/(SM-SR)*(-pStar_TR);
				
				// double pStar_YL[nSp];
				// double pStar_YR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// pStar_YL[i] = 0.5*(drhodYL[i]*(UnL-SL)*(UnL-SM) + rhoL*(UnL-SL)*(-SM_YL[i]) +
									  // rhoR*(UnR-SR)*(-SM_YL[i]) );
					// pStar_YR[i] = 0.5*(rhoL*(UnL-SL)*(-SM_YR[i]) + drhodYR[i]*(UnR-SR)*(UnR-SM) +
									  // rhoR*(UnR-SR)*(-SM_YR[i]) );
					// pF_YL[i] = pStar_YL[i] + SM_YL[i]/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_YL[i] +
								// SM/(SM-SR)*(-pStar_YL[i]);
					// pF_YR[i] = pStar_YR[i] + SM_YR[i]/(SM-SR)*(pR-pStar) - SM/(SM-SR)/(SM-SR)*(pR-pStar)*SM_YR[i] +
								// SM/(SM-SR)*(-pStar_YR[i]);
				// }
				
				
				
				// dissEg_pL = weiL*SL*(pStar_pL-1.0)/(SL-UnL) + weiR*SR*(pStar_pL)/(SR-UnR);
				// dissEg_pR = weiL*SL*(pStar_pR)/(SL-UnL) + weiR*SR*(pStar_pR-1.0)/(SR-UnR);
				// dissEg_TL = weiL*SL*(pStar_TL)/(SL-UnL) + weiR*SR*(pStar_TL)/(SR-UnR);
				// dissEg_TR = weiL*SL*(pStar_TR)/(SL-UnL) + weiR*SR*(pStar_TR)/(SR-UnR);
				// for(int i=0; i<nSp-1; ++i){
					// dissEg_YL[i] = weiL*SL*(pStar_YL[i])/(SL-UnL) + weiR*SR*(pStar_YL[i])/(SR-UnR);
					// dissEg_YR[i] = weiL*SL*(pStar_YR[i])/(SL-UnL) + weiR*SR*(pStar_YR[i])/(SR-UnR);
				// }
				// dissEg_UnL = weiL*SL*(pStar_UnL)/(SL-UnL) - weiL*SL*(pStar-pL)/(SL-UnL)/(SL-UnL)*(-1.0) + 
							 // weiR*SR*(pStar_UnL)/(SR-UnR);
				// dissEg_UnR = weiL*SL*(pStar_UnR)/(SL-UnL) - weiR*SR*(pStar-pR)/(SR-UnR)/(SR-UnR)*(-1.0) + 
							 // weiR*SR*(pStar_UnR)/(SR-UnR);
				
			// }
			
			
			
			
			
			
			