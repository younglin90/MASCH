
#include "../../../others/solvers.h"

void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int nSp = controls.spName.size();
	
	
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
	
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	// int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell"); 
	int id_wd = controls.getId_faceVar("distance weight"); 
	int id_xLR = controls.getId_faceVar("x distance of between left and right cell");
	int id_yLR = controls.getId_faceVar("y distance of between left and right cell");
	int id_zLR = controls.getId_faceVar("z distance of between left and right cell");
	
	{
		calcHO_FaceVal.push_back(
			[&solver,
			id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_xLR,id_yLR,id_zLR,
			id_p,id_pL,id_pR,id_dpdx,id_dpdy,id_dpdz,
			id_u,id_uL,id_uR,id_dudx,id_dudy,id_dudz,
			id_v,id_vL,id_vR,id_dvdx,id_dvdy,id_dvdz,
			id_w,id_wL,id_wR,id_dwdx,id_dwdy,id_dwdz,
			id_T,id_TL,id_TR,id_dTdx,id_dTdy,id_dTdz,
			id_Y,id_YL,id_YR,id_dYdx,id_dYdy,id_dYdz,
			id_wd,id_rho,id_rhoL,id_rhoR,id_Ht,id_HtL,id_HtR,
			id_c,id_cL,id_cR,nSp](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{

				faces[id_pL] = cellsL[id_p];
				faces[id_uL] = cellsL[id_u];
				faces[id_vL] = cellsL[id_v];
				faces[id_wL] = cellsL[id_w];
				faces[id_TL] = cellsL[id_T];
				faces[id_rhoL] = cellsL[id_rho];
				faces[id_cL] = cellsL[id_c];
				faces[id_HtL] = cellsL[id_Ht];
				
				faces[id_pR] = cellsR[id_p];
				faces[id_uR] = cellsR[id_u];
				faces[id_vR] = cellsR[id_v];
				faces[id_wR] = cellsR[id_w];
				faces[id_TR] = cellsR[id_T];
				faces[id_rhoR] = cellsR[id_rho];
				faces[id_cR] = cellsR[id_c];
				faces[id_HtR] = cellsR[id_Ht];
				
				// double nvec[3];
				// nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				// double area = faces[id_area];
				// double dLR = faces[id_dLR]; 
				// double wdL = faces[id_wd]; double wdR = 1.0-wdL;
				// double dt = fields[id_dt];
				// double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
				// double dtrho = dt*(wdL/rhoL+wdR/rhoR);
				
				// double cL = cellsL[id_c]; double cR = cellsR[id_c];
				// double pL = cellsL[id_p]; double pR = cellsR[id_p];
				// double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
				// double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
				// double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
				// double uL = cellsL[id_u]; double uR = cellsR[id_u];
				// double dudxL = cellsL[id_dudx]; double dudxR = cellsR[id_dudx];
				// double dudyL = cellsL[id_dudy]; double dudyR = cellsR[id_dudy];
				// double dudzL = cellsL[id_dudz]; double dudzR = cellsR[id_dudz];
				// double vL = cellsL[id_v]; double vR = cellsR[id_v];
				// double dvdxL = cellsL[id_dvdx]; double dvdxR = cellsR[id_dvdx];
				// double dvdyL = cellsL[id_dvdy]; double dvdyR = cellsR[id_dvdy];
				// double dvdzL = cellsL[id_dvdz]; double dvdzR = cellsR[id_dvdz];
				// double wL = cellsL[id_w]; double wR = cellsR[id_w];
				// double dwdxL = cellsL[id_dwdx]; double dwdxR = cellsR[id_dwdx];
				// double dwdyL = cellsL[id_dwdy]; double dwdyR = cellsR[id_dwdy];
				// double dwdzL = cellsL[id_dwdz]; double dwdzR = cellsR[id_dwdz];
				// double TL = cellsL[id_T]; double TR = cellsR[id_T];
				// double dTdxL = cellsL[id_dTdx]; double dTdxR = cellsR[id_dTdx];
				// double dTdyL = cellsL[id_dTdy]; double dTdyR = cellsR[id_dTdy];
				// double dTdzL = cellsL[id_dTdz]; double dTdzR = cellsR[id_dTdz];
				// double YL[nSp]; double YR[nSp];
				// double dYdxL[nSp]; double dYdxR[nSp];
				// double dYdyL[nSp]; double dYdyR[nSp];
				// double dYdzL[nSp]; double dYdzR[nSp];
				// for(int i=0; i<nSp-1; ++i){
					// YL[i] = cellsL[id_Y[i]]; YR[i] = cellsR[id_Y[i]];
					// dYdxL[i] = cellsL[id_dYdx[i]]; dYdxR[i] = cellsR[id_dYdx[i]];
					// dYdyL[i] = cellsL[id_dYdy[i]]; dYdyR[i] = cellsR[id_dYdy[i]];
					// dYdzL[i] = cellsL[id_dYdz[i]]; dYdzR[i] = cellsR[id_dYdz[i]];
				// }
				
				// double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
				// double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
				
				// double UnF = wdL*UnL+wdR*UnR;
				// UnF -= dtrho*(pR-pL)/dLR;
				// UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*nvec[0];
				// UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*nvec[1];
				// UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*nvec[2];
				// double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
				// double cbar = wdL*cL+wdR*cR;
				// double Mdash = min(1.0,KLR/cbar);
				// double chi = (1.0-Mdash)*(1.0-Mdash);
				// double rhohat = wdL*rhoL+wdR*rhoR;
				// UnF -= 0.5*chi/rhohat/cbar*(pR-pL);
				
				// faces[id_UnF] = UnF;
				
				// double phiL2[5+nSp], phiL1[5+nSp], phiR1[5+nSp];
				// if(UnF>=0.0){
					// phiL1[0] = uL; phiR1[0] = uR;
					// phiL2[0] = 0.0;
					// phiL2[0] += dudxL*faces[id_xLR];
					// phiL2[0] += dudyL*faces[id_yLR];
					// phiL2[0] += dudzL*faces[id_zLR];
					// phiL2[0] = phiR1[0] - 2.0*phiL2[0];
					
					// phiL1[1] = vL; phiR1[1] = vR;
					// phiL2[1] = 0.0;
					// phiL2[1] += dvdxL*faces[id_xLR];
					// phiL2[1] += dvdyL*faces[id_yLR];
					// phiL2[1] += dvdzL*faces[id_zLR];
					// phiL2[1] = phiR1[1] - 2.0*phiL2[1];
					
					// phiL1[2] = wL; phiR1[2] = wR;
					// phiL2[2] = 0.0;
					// phiL2[2] += dwdxL*faces[id_xLR];
					// phiL2[2] += dwdyL*faces[id_yLR];
					// phiL2[2] += dwdzL*faces[id_zLR];
					// phiL2[2] = phiR1[2] - 2.0*phiL2[2];
					
					// phiL1[3] = TL; phiR1[3] = TR;
					// phiL2[3] = 0.0;
					// phiL2[3] += dTdxL*faces[id_xLR];
					// phiL2[3] += dTdyL*faces[id_yLR];
					// phiL2[3] += dTdzL*faces[id_zLR];
					// phiL2[3] = phiR1[3] - 2.0*phiL2[3];
					
					// for(int i=0; i<nSp-1; ++i){
						// int tmp_i = 4+i;
						// phiL1[tmp_i] = YL[i]; phiR1[tmp_i] = YR[i];
						// phiL2[tmp_i] = 0.0;
						// phiL2[tmp_i] += dYdxL[i]*faces[id_xLR];
						// phiL2[tmp_i] += dYdyL[i]*faces[id_yLR];
						// phiL2[tmp_i] += dYdzL[i]*faces[id_zLR];
						// phiL2[tmp_i] = phiR1[tmp_i] - 2.0*phiL2[tmp_i];
					// }
					
					
				// }
				// else{
					// phiL1[0] = uR; phiR1[0] = uL;
					// phiL2[0] = 0.0;
					// phiL2[0] += dudxR*faces[id_xLR];
					// phiL2[0] += dudyR*faces[id_yLR];
					// phiL2[0] += dudzR*faces[id_zLR];
					// phiL2[0] = phiR1[0] + 2.0*phiL2[0];
					
					// phiL1[1] = vR; phiR1[1] = vL;
					// phiL2[1] = 0.0;
					// phiL2[1] += dvdxR*faces[id_xLR];
					// phiL2[1] += dvdyR*faces[id_yLR];
					// phiL2[1] += dvdzR*faces[id_zLR];
					// phiL2[1] = phiR1[1] + 2.0*phiL2[1];
					
					// phiL1[2] = wR; phiR1[2] = wL;
					// phiL2[2] = 0.0;
					// phiL2[2] += dwdxR*faces[id_xLR];
					// phiL2[2] += dwdyR*faces[id_yLR];
					// phiL2[2] += dwdzR*faces[id_zLR];
					// phiL2[2] = phiR1[2] + 2.0*phiL2[2];
					
					// phiL1[3] = TR; phiR1[3] = TL;
					// phiL2[3] = 0.0;
					// phiL2[3] += dTdxR*faces[id_xLR];
					// phiL2[3] += dTdyR*faces[id_yLR];
					// phiL2[3] += dTdzR*faces[id_zLR];
					// phiL2[3] = phiR1[3] + 2.0*phiL2[3];
					
					// for(int i=0; i<nSp-1; ++i){
						// int tmp_i = 4+i;
						// phiL1[tmp_i] = YR[i]; phiR1[tmp_i] = YL[i];
						// phiL2[tmp_i] = 0.0;
						// phiL2[tmp_i] += dYdxR[i]*faces[id_xLR];
						// phiL2[tmp_i] += dYdyR[i]*faces[id_yLR];
						// phiL2[tmp_i] += dYdzR[i]*faces[id_zLR];
						// phiL2[tmp_i] = phiR1[tmp_i] + 2.0*phiL2[tmp_i];
					// }
				// }
				
				// faces[id_pF] = wdL*pL+wdR*pR;
				
				// faces[id_uF] = phiL1[0];
				// faces[id_vF] = phiL1[1];
				// faces[id_wF] = phiL1[2];
				// faces[id_TF] = phiL1[3];
				
				// // faces[id_uF] = solver.NVD.Minmod(phiL2[0],phiL1[0],phiR1[0]);
				// // faces[id_vF] = solver.NVD.Minmod(phiL2[1],phiL1[1],phiR1[1]);
				// // faces[id_wF] = solver.NVD.Minmod(phiL2[2],phiL1[2],phiR1[2]);
				// // faces[id_TF] = solver.NVD.Minmod(phiL2[3],phiL1[3],phiR1[3]);
				
				// {
					// double dx = sqrt(faces[id_area]);
					// double dt_tmp_L = dx/(sqrt(uL*uL+vL*vL+wL*wL)+1.e-100);
					// double dt_tmp_R = dx/(sqrt(uR*uR+vR*vR+wR*wR)+1.e-100);			
					// double coDD = max(dt/dt_tmp_L,dt/dt_tmp_R);
					
					// for(int i=0; i<nSp-1; ++i){
						// double mfLR[3];
						// mfLR[0] = 0.5*dYdxL[i]+0.5*dYdxR[i];
						// mfLR[1] = 0.5*dYdyL[i]+0.5*dYdyR[i];
						// mfLR[2] = 0.5*dYdzL[i]+0.5*dYdzR[i];
						// double magMfLR = mfLR[0]*mfLR[0];
						// magMfLR += mfLR[1]*mfLR[1];
						// magMfLR += mfLR[2]*mfLR[2];
						// magMfLR = sqrt(magMfLR);
						// mfLR[0] = mfLR[0]/(magMfLR+1.e-200);
						// mfLR[1] = mfLR[1]/(magMfLR+1.e-200);
						// mfLR[2] = mfLR[2]/(magMfLR+1.e-200);
						// double cosTheta = mfLR[0]*faces[id_nx];
						// cosTheta += mfLR[1]*faces[id_ny];
						// cosTheta += mfLR[2]*faces[id_nz];
						// cosTheta = abs(cosTheta);
						// double gamF = min(cosTheta*cosTheta*cosTheta*cosTheta,1.0);
						
						// // faces[id_YF[i]] = solver.NVD.Minmod(phiL2[4+i],phiL1[4+i],phiR1[4+i]);
						// faces[id_YF[i]] = solver.NVD.getHO_MSTACS(phiL2[4+i],phiL1[4+i],phiR1[4+i],coDD,gamF);
					// }
				// }
				
				
				return 0;
			}
		);
	}
	// cout << "BBBBBBBB" << endl;
	
}

