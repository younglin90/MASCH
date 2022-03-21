
#include "../../../../others/solvers.h"

void MASCH_Solver::setDiffFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// checkImplicitLaplFlux.push_back(false);
	
	{
		using US = unsigned short;
		int nSp = controls.cellVar["mass fraction"].sub_name.size();
		
		// 셀값
		int id_p = controls.cellVar["pressure"].id;
		int id_u = controls.cellVar["x-velocity"].id;
		int id_v = controls.cellVar["y-velocity"].id;
		int id_w = controls.cellVar["z-velocity"].id;
		int id_T = controls.cellVar["temperature"].id;
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass fraction"].sub_name[i];
			id_Y[i] = controls.cellVar[tmp_name].id;
		}
		int id_rho = controls.cellVar["density"].id;
		int id_c = controls.cellVar["speed of sound"].id;
		int id_Ht = controls.cellVar["total enthalpy"].id;
		
		// left값
		int id_pL = controls.faceVar["left pressure"].id;
		int id_uL = controls.faceVar["left x-velocity"].id;
		int id_vL = controls.faceVar["left y-velocity"].id;
		int id_wL = controls.faceVar["left z-velocity"].id;
		int id_TL = controls.faceVar["left temperature"].id;
		vector<int> id_YL(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.faceVar["left mass fraction"].sub_name[i];
			id_YL[i] = controls.faceVar[tmp_name].id;
		}
		int id_rhoL = controls.faceVar["left pressure"].id;
		int id_cL = controls.faceVar["left x-velocity"].id;
		int id_HtL = controls.faceVar["left y-velocity"].id;
		
		// right값
		int id_pR = controls.faceVar["right pressure"].id;
		int id_uR = controls.faceVar["right x-velocity"].id;
		int id_vR = controls.faceVar["right y-velocity"].id;
		int id_wR = controls.faceVar["right z-velocity"].id;
		int id_TR = controls.faceVar["right temperature"].id;
		vector<int> id_YR(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.faceVar["right mass fraction"].sub_name[i];
			id_YR[i] = controls.faceVar[tmp_name].id;
		}
		int id_rhoR = controls.faceVar["right pressure"].id;
		int id_cR = controls.faceVar["right x-velocity"].id;
		int id_HtR = controls.faceVar["right y-velocity"].id;
		
		// 디퓨젼 계수 관련
		int id_mu = controls.faceVar["viscosity"].id;
		
		// 그레디언트 텀 관련
		int id_dudx = controls.faceVar["x-gradient x-velocity"].id;
		int id_dudy = controls.faceVar["y-gradient x-velocity"].id;
		int id_dudz = controls.faceVar["z-gradient x-velocity"].id;
		int id_dvdx = controls.faceVar["x-gradient y-velocity"].id;
		int id_dvdy = controls.faceVar["y-gradient y-velocity"].id;
		int id_dvdz = controls.faceVar["z-gradient y-velocity"].id;
		int id_dwdx = controls.faceVar["x-gradient z-velocity"].id;
		int id_dwdy = controls.faceVar["y-gradient z-velocity"].id;
		int id_dwdz = controls.faceVar["z-gradient z-velocity"].id;
		
		// 메쉬관련
		int id_nx = controls.faceVar["x unit normal"].id;
		int id_ny = controls.faceVar["y unit normal"].id;
		int id_nz = controls.faceVar["z unit normal"].id;
		int id_area = controls.faceVar["area"].id;
		int id_Wc = controls.faceVar["distance weight"].id;
		int id_dLR = controls.faceVar["distance of between left and right cell"].id;
		
		calcLaplFlux.push_back(
		[nSp,id_nx,id_ny,id_nz,id_area,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
		id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
		id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR,
		id_mu,id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,
		id_dwdx,id_dwdy,id_dwdz,id_Wc,id_dLR](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			double muL, muR;
			double uL, uR, vL, vR, wL, wR, rhoL, rhoR;
			double dudxL,dudyL,dudzL,dvdxL,dvdyL,dvdzL,dwdxL,dwdyL,dwdzL;
			double dudxR,dudyR,dudzR,dvdxR,dvdyR,dvdzR,dwdxR,dwdyR,dwdzR;
			if(cellsR!=nullptr){
				muL = cellsL[id_mu]; muR = cellsR[id_mu];
				uL = cellsL[id_u]; uR = cellsR[id_u];
				vL = cellsL[id_u]; vR = cellsR[id_u];
				wL = cellsL[id_u]; wR = cellsR[id_u];
				rhoL = cellsL[id_u]; rhoR = cellsR[id_u];
				
				dudxL = cellsL[id_dudx]; dudxR = cellsR[id_dudx];
				dudyL = cellsL[id_dudy]; dudyR = cellsR[id_dudy];
				dudzL = cellsL[id_dudz]; dudzR = cellsR[id_dudz];
				dvdxL = cellsL[id_dvdx]; dvdxR = cellsR[id_dvdx];
				dvdyL = cellsL[id_dvdy]; dvdyR = cellsR[id_dvdy];
				dvdzL = cellsL[id_dvdz]; dvdzR = cellsR[id_dvdz];
				dwdxL = cellsL[id_dwdx]; dwdxR = cellsR[id_dwdx];
				dwdyL = cellsL[id_dwdy]; dwdyR = cellsR[id_dwdy];
				dwdzL = cellsL[id_dwdz]; dwdzR = cellsR[id_dwdz];
			}
			else{
				muL = cellsL[id_mu]; muR = cellsL[id_mu];
				uL = cellsL[id_u]; uR = faces[id_uR];
				vL = cellsL[id_u]; vR = faces[id_vR];
				wL = cellsL[id_u]; wR = faces[id_wR];
				rhoL = cellsL[id_u]; rhoR = faces[id_rhoR];
				
				dudxL = cellsL[id_dudx]; dudxR = cellsL[id_dudx];
				dudyL = cellsL[id_dudy]; dudyR = cellsL[id_dudy];
				dudzL = cellsL[id_dudz]; dudzR = cellsL[id_dudz];
				dvdxL = cellsL[id_dvdx]; dvdxR = cellsL[id_dvdx];
				dvdyL = cellsL[id_dvdy]; dvdyR = cellsL[id_dvdy];
				dvdzL = cellsL[id_dvdz]; dvdzR = cellsL[id_dvdz];
				dwdxL = cellsL[id_dwdx]; dwdxR = cellsL[id_dwdx];
				dwdyL = cellsL[id_dwdy]; dwdyR = cellsL[id_dwdy];
				dwdzL = cellsL[id_dwdz]; dwdzR = cellsL[id_dwdz];
			}
			
			double nvec[3];
			nvec[0] = faces[id_nx];
			nvec[1] = faces[id_ny];
			nvec[2] = faces[id_nz];
			
			double gc = faces[id_Wc];
			
			double muF = gc*muL + (1.0-gc)*muR;
			double rhoF = gc*rhoL + (1.0-gc)*rhoR;
			double uF = gc*uL + (1.0-gc)*uR;
			double vF = gc*vL + (1.0-gc)*vR;
			double wF = gc*wL + (1.0-gc)*wR;
			
			double dLR = faces[id_dLR];
			double dudn = (uR-uL)/dLR;
			double dvdn = (vR-vL)/dLR;
			double dwdn = (wR-wL)/dLR;
			
			double area = faces[id_area];
			int iter=0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = (muF * dudn)*area;
			fluxB[iter++] = (muF * dvdn)*area;
			fluxB[iter++] = (muF * dwdn)*area;
			fluxB[iter++] = (muF * (uF*dudn + vF*dvdn + wF*dwdn))*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = 0.0;
			}
			return 0;
		}); 
		
		
		
		
		
		calcNLaplFlux.push_back(
		[nSp,id_nx,id_ny,id_nz,id_area,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
		id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
		id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR,
		id_mu,id_dudx,id_dudy,id_dudz,id_dvdx,id_dvdy,id_dvdz,
		id_dwdx,id_dwdy,id_dwdz,id_Wc](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			double muL, muR;
			double uL, uR, vL, vR, wL, wR, rhoL, rhoR;
			double dudxL,dudyL,dudzL,dvdxL,dvdyL,dvdzL,dwdxL,dwdyL,dwdzL;
			double dudxR,dudyR,dudzR,dvdxR,dvdyR,dvdzR,dwdxR,dwdyR,dwdzR;
			if(cellsR!=nullptr){
				muL = cellsL[id_mu]; muR = cellsR[id_mu];
				uL = cellsL[id_u]; uR = cellsR[id_u];
				vL = cellsL[id_u]; vR = cellsR[id_u];
				wL = cellsL[id_u]; wR = cellsR[id_u];
				rhoL = cellsL[id_u]; rhoR = cellsR[id_u];
				
				dudxL = cellsL[id_dudx]; dudxR = cellsR[id_dudx];
				dudyL = cellsL[id_dudy]; dudyR = cellsR[id_dudy];
				dudzL = cellsL[id_dudz]; dudzR = cellsR[id_dudz];
				dvdxL = cellsL[id_dvdx]; dvdxR = cellsR[id_dvdx];
				dvdyL = cellsL[id_dvdy]; dvdyR = cellsR[id_dvdy];
				dvdzL = cellsL[id_dvdz]; dvdzR = cellsR[id_dvdz];
				dwdxL = cellsL[id_dwdx]; dwdxR = cellsR[id_dwdx];
				dwdyL = cellsL[id_dwdy]; dwdyR = cellsR[id_dwdy];
				dwdzL = cellsL[id_dwdz]; dwdzR = cellsR[id_dwdz];
			}
			else{
				muL = cellsL[id_mu]; muR = cellsL[id_mu];
				uL = cellsL[id_u]; uR = faces[id_uR];
				vL = cellsL[id_u]; vR = faces[id_vR];
				wL = cellsL[id_u]; wR = faces[id_wR];
				rhoL = cellsL[id_u]; rhoR = faces[id_rhoR];
				
				dudxL = cellsL[id_dudx]; dudxR = cellsL[id_dudx];
				dudyL = cellsL[id_dudy]; dudyR = cellsL[id_dudy];
				dudzL = cellsL[id_dudz]; dudzR = cellsL[id_dudz];
				dvdxL = cellsL[id_dvdx]; dvdxR = cellsL[id_dvdx];
				dvdyL = cellsL[id_dvdy]; dvdyR = cellsL[id_dvdy];
				dvdzL = cellsL[id_dvdz]; dvdzR = cellsL[id_dvdz];
				dwdxL = cellsL[id_dwdx]; dwdxR = cellsL[id_dwdx];
				dwdyL = cellsL[id_dwdy]; dwdyR = cellsL[id_dwdy];
				dwdzL = cellsL[id_dwdz]; dwdzR = cellsL[id_dwdz];
			}
			
			double nvec[3];
			nvec[0] = faces[id_nx];
			nvec[1] = faces[id_ny];
			nvec[2] = faces[id_nz];
			
			double gc = faces[id_Wc];
			
			double muF = gc*muL + (1.0-gc)*muR;
			double rhoF = gc*rhoL + (1.0-gc)*rhoR;
			double uF = gc*uL + (1.0-gc)*uR;
			double vF = gc*vL + (1.0-gc)*vR;
			double wF = gc*wL + (1.0-gc)*wR;
			
			double dudxF = gc*dudxL + (1.0-gc)*dudxR;
			double dudyF = gc*dudyL + (1.0-gc)*dudyR;
			double dudzF = gc*dudzL + (1.0-gc)*dudzR;
			double dvdxF = gc*dvdxL + (1.0-gc)*dvdxR;
			double dvdyF = gc*dvdyL + (1.0-gc)*dvdyR;
			double dvdzF = gc*dvdzL + (1.0-gc)*dvdzR;
			double dwdxF = gc*dwdxL + (1.0-gc)*dwdxR;
			double dwdyF = gc*dwdyL + (1.0-gc)*dwdyR;
			double dwdzF = gc*dwdzL + (1.0-gc)*dwdzR;
			
			// Ref : Blazek's book, pp.20-21
			double lambda = -2.0/3.0;
			double div_dot_U = dudxF + dvdyF + dwdzF;
			
			double x_mom = dudxF*nvec[0] + dvdxF*nvec[1] + dwdxF*nvec[2] + 
							lambda*div_dot_U*nvec[0];
			double y_mom = dudyF*nvec[0] + dvdyF*nvec[1] + dwdyF*nvec[2] + 
							lambda*div_dot_U*nvec[1];
			double z_mom = dudzF*nvec[0] + dvdzF*nvec[1] + dwdzF*nvec[2] + 
							lambda*div_dot_U*nvec[2];
			double ener = 
			(dudxF*uF + dudyF*vF + dudzF*wF + lambda*div_dot_U*uF)*nvec[0] +
			(dvdxF*uF + dvdyF*vF + dvdzF*wF + lambda*div_dot_U*vF)*nvec[1] +
			(dwdxF*uF + dwdyF*vF + dwdzF*wF + lambda*div_dot_U*wF)*nvec[2];
			
			double area = faces[id_area];
			int iter=0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = (muF * x_mom)*area;
			fluxB[iter++] = (muF * y_mom)*area;
			fluxB[iter++] = (muF * z_mom)*area;
			fluxB[iter++] = (muF * ener)*area;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = 0.0;
			}
			return 0;
			
		}); 
		
		
		
	}
	
	
	
	
	
	
}