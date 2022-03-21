
#include "../../../others/solvers.h"

void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int nSp = controls.spName.size();
	
	int id_c = controls.getId_cellVar("speed-of-sound");
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pF = controls.getId_faceVar("pressure");
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uF = controls.getId_faceVar("x-velocity");
	int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	
	int id_T = controls.getId_cellVar("temperature");
	int id_TF = controls.getId_faceVar("temperature");
	int id_dTdx = controls.getId_cellVar("x-gradient temperature");
	int id_dTdy = controls.getId_cellVar("y-gradient temperature");
	int id_dTdz = controls.getId_cellVar("z-gradient temperature");
	
	vector<int> id_Y, id_YF, id_dYdx, id_dYdy, id_dYdz;
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		id_Y.push_back(controls.getId_cellVar(tmp_name));
		id_YF.push_back(controls.getId_faceVar(tmp_name));
		id_dYdx.push_back(controls.getId_cellVar("x-gradient "+tmp_name));
		id_dYdy.push_back(controls.getId_cellVar("y-gradient "+tmp_name));
		id_dYdz.push_back(controls.getId_cellVar("z-gradient "+tmp_name));
	}
	
	int id_rho = controls.getId_cellVar("density");
	int id_rhoF = controls.getId_faceVar("density");
	
	int id_Ht = controls.getId_cellVar("total-enthalpy");
	int id_HtF = controls.getId_faceVar("total-enthalpy");
	
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
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
	
	int id_xLNv = controls.getId_faceVar("left cell to face normal vector shortest x distance");
	int id_yLNv = controls.getId_faceVar("left cell to face normal vector shortest y distance");
	int id_zLNv = controls.getId_faceVar("left cell to face normal vector shortest z distance");
	int id_xRNv = controls.getId_faceVar("right cell to face normal vector shortest x distance");
	int id_yRNv = controls.getId_faceVar("right cell to face normal vector shortest y distance");
	int id_zRNv = controls.getId_faceVar("right cell to face normal vector shortest z distance");
	
	int id_alpha = controls.getId_faceVar("cosine angle of between face normal and cells");
	
	int id_nLRx = controls.getId_faceVar("x unit normal of between left and right cell");
	int id_nLRy = controls.getId_faceVar("y unit normal of between left and right cell");
	int id_nLRz = controls.getId_faceVar("z unit normal of between left and right cell");
	
	{
		calcHO_FaceVal.push_back(
			[&solver,
			id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_xLR,id_yLR,id_zLR,
			id_p,id_pF,id_dpdx,id_dpdy,id_dpdz,
			id_u,id_uF,id_dudx,id_dudy,id_dudz,
			id_v,id_vF,id_dvdx,id_dvdy,id_dvdz,
			id_w,id_wF,id_dwdx,id_dwdy,id_dwdz,
			id_T,id_TF,id_dTdx,id_dTdy,id_dTdz,
			id_Y,id_YF,id_dYdx,id_dYdy,id_dYdz,
			id_wd,id_rho,id_rhoF,id_Ht,id_HtF,
			id_c,nSp,id_alpha,id_nLRx,id_nLRy,id_nLRz,
			id_xLNv,id_yLNv,id_zLNv,id_xRNv,id_yRNv,id_zRNv](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				double LNv[3];
				LNv[0] = faces[id_xLNv]; LNv[1] = faces[id_yLNv]; LNv[2] = faces[id_zLNv];
				double RNv[3];
				RNv[0] = faces[id_xRNv]; RNv[1] = faces[id_yRNv]; RNv[2] = faces[id_zRNv];
				double wdL = faces[id_wd]; double wdR = 1.0-wdL;
				
				
				
				wdL = 0.5; wdR = 0.5;
				
				
				
				
				double dt = fields[id_dt];
				double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
				double dtrho = dt*(wdL/rhoL+wdR/rhoR);
				
				double cL = cellsL[id_c]; double cR = cellsR[id_c];
				double pL = cellsL[id_p]; double pR = cellsR[id_p];
				double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
				double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
				double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
				double uL = cellsL[id_u]; double uR = cellsR[id_u];
				double dudxL = cellsL[id_dudx]; double dudxR = cellsR[id_dudx];
				double dudyL = cellsL[id_dudy]; double dudyR = cellsR[id_dudy];
				double dudzL = cellsL[id_dudz]; double dudzR = cellsR[id_dudz];
				double vL = cellsL[id_v]; double vR = cellsR[id_v];
				double dvdxL = cellsL[id_dvdx]; double dvdxR = cellsR[id_dvdx];
				double dvdyL = cellsL[id_dvdy]; double dvdyR = cellsR[id_dvdy];
				double dvdzL = cellsL[id_dvdz]; double dvdzR = cellsR[id_dvdz];
				double wL = cellsL[id_w]; double wR = cellsR[id_w];
				double dwdxL = cellsL[id_dwdx]; double dwdxR = cellsR[id_dwdx];
				double dwdyL = cellsL[id_dwdy]; double dwdyR = cellsR[id_dwdy];
				double dwdzL = cellsL[id_dwdz]; double dwdzR = cellsR[id_dwdz];
				double TL = cellsL[id_T]; double TR = cellsR[id_T];
				double dTdxL = cellsL[id_dTdx]; double dTdxR = cellsR[id_dTdx];
				double dTdyL = cellsL[id_dTdy]; double dTdyR = cellsR[id_dTdy];
				double dTdzL = cellsL[id_dTdz]; double dTdzR = cellsR[id_dTdz];
				double YL[nSp]; double YR[nSp];
				double dYdxL[nSp]; double dYdxR[nSp];
				double dYdyL[nSp]; double dYdyR[nSp];
				double dYdzL[nSp]; double dYdzR[nSp];
				for(int i=0; i<nSp-1; ++i){
					YL[i] = cellsL[id_Y[i]]; YR[i] = cellsR[id_Y[i]];
					dYdxL[i] = cellsL[id_dYdx[i]]; dYdxR[i] = cellsR[id_dYdx[i]];
					dYdyL[i] = cellsL[id_dYdy[i]]; dYdyR[i] = cellsR[id_dYdy[i]];
					dYdzL[i] = cellsL[id_dYdz[i]]; dYdzR[i] = cellsR[id_dYdz[i]];
				}
				
				double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
				double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
				
				
				
				// // skewness
				// UnL += (dudxL*LNv[0]+dudyL*LNv[1]+dudzL*LNv[2])*nvec[0];
				// UnL += (dvdxL*LNv[0]+dvdyL*LNv[1]+dvdzL*LNv[2])*nvec[1];
				// UnL += (dwdxL*LNv[0]+dwdyL*LNv[1]+dwdzL*LNv[2])*nvec[2];
				// UnR += (dudxR*RNv[0]+dudyR*RNv[1]+dudzR*RNv[2])*nvec[0];
				// UnR += (dvdxR*RNv[0]+dvdyR*RNv[1]+dvdzR*RNv[2])*nvec[1];
				// UnR += (dwdxR*RNv[0]+dwdyR*RNv[1]+dwdzR*RNv[2])*nvec[2];
				
				
				
				

				double KLR = sqrt(wdL*(uL*uL+vL*vL+wL*wL)+wdR*(uR*uR+vR*vR+wR*wR));
				// double wdiL = faces[id_wd]; double wdiR = 1.0-wdiL;
				double w_phi = 1.0-2.0*abs(wdL-0.5);
				
				
				
				
				
				w_phi = 1.0;
				
				
				
				
				double chat= wdL*cL+wdR*cR;
				double Mcy = min(1.0,KLR/chat);
				double phi_c = (1.0-Mcy)*(1.0-Mcy);
				// double Unhat = 0.5*(UnL+UnR);
				double rhohat = wdL*rhoL+wdR*rhoR;
				double ML = UnL/chat; 
				double MR = UnR/chat;
				double MLP = 0.5*(ML+abs(ML));
				if( abs(ML) < 1.0 ) {
					MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				}
				double MRM = 0.5*(MR-abs(MR));
				if( abs(MR) < 1.0 ) {
					MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				}
				// 영린 개발 스킴
				double UnF = w_phi*(MLP+MRM)*chat + (1.0-w_phi)*(wdL*UnL+wdR*UnR);
				
				// 열역학적 보간
				UnF -= 0.5*phi_c/rhohat/chat*(pR-pL);
				// 셀 to 페이스 보간
				UnF -= dAlpha * dtrho*(pR-pL)/dLR;
				// non-orthogonal
				UnF -= (nvec[0]-dAlpha*nLR[0]) * dtrho*(wdL*dpdxL+wdR*dpdxR);
				UnF -= (nvec[1]-dAlpha*nLR[1]) * dtrho*(wdL*dpdyL+wdR*dpdyR);
				UnF -= (nvec[2]-dAlpha*nLR[2]) * dtrho*(wdL*dpdzL+wdR*dpdzR);
				
				UnF += dt*(wdL*dpdxL/rhoL+wdR*dpdxR/rhoR)*nvec[0];
				UnF += dt*(wdL*dpdyL/rhoL+wdR*dpdyR/rhoR)*nvec[1];
				UnF += dt*(wdL*dpdzL/rhoL+wdR*dpdzR/rhoR)*nvec[2];
				
				
				
				
				// UnF = 0.5*(UnL+UnR) - dtrho*(pR-pL)/dLR + dtrho*(
					// 0.5*(dpdxL+dpdxR)*nvec[0] + 
					// 0.5*(dpdyL+dpdyR)*nvec[1] + 
					// 0.5*(dpdzL+dpdzR)*nvec[2]);
				
				// UnF = 0.5*(UnL+UnR) - 0.5*phi_c/rhohat/chat*(pR-pL);
				// UnF = (MLP+MRM)*chat  - 0.5*phi_c/rhohat/chat*(pR-pL);
				
				
				
				faces[id_UnF] = UnF;
			
			
				
				double PLP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
				if( abs(ML) < 1.0 ) {
					PLP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
				} 
				double PRM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
				if( abs(MR) < 1.0 ) {
					PRM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
				} 
				double WPL = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
				double WPR = 0.5 - 0.5*Mcy*PLP*PRM*0.5/chat*(UnR-UnL) + 0.5*Mcy*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
				// 영린 개발 스킴
				double pF = (w_phi*WPL+(1.0-w_phi)*wdL)*pL + (w_phi*WPR+(1.0-w_phi)*wdR)*pR;
				
				
				
				
				// pF = 0.5*(pL+pR);
				
				
				
				
				faces[id_pF] = pF;
				
				
				
				
				
				
				
				
				
				


			// double Unhat = 0.5*(UnL+UnR);
			// // double rhohat = 0.5*(rhoL+rhoR);
			// // double chat= 0.5*(cL+cR);
			
			// // double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
			// // double MLP = 0.5*(ML+abs(ML));
			// // if( abs(ML) < 1.0 ) {
				// // MLP = 0.25*(ML + 1.0)*(ML + 1.0);
			// // }
			// // double MRM = 0.5*(MR-abs(MR));
			// // if( abs(MR) < 1.0 ) {
				// // MRM = -0.25*(MR - 1.0)*(MR - 1.0);
			// // }
			// double preP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
			// if( abs(ML) < 1.0 ) {
				// preP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
			// } 
			// double preM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
			// if( abs(MR) < 1.0 ) {
				// preM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
			// } 
			
			// double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );

			// // SLAU
			// // double Mcy = min(1.0,KLR/chat);
			// // double phi_c = (1.0-Mcy)*(1.0-Mcy);
			// double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
			// double D_L = ML+(1.0-g_c)*abs(ML);
			// double D_R = MR-(1.0-g_c)*abs(MR);
			// double D_rho = Mbar*g_c;
			// double MLPL = 0.5*(D_L+D_rho);
			// double MRMR = 0.5*(D_R-D_rho);

			// double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);

			// UnF = mdot/rhoL;
			// if( mdot < 0.0 ) {
				// UnF = mdot/rhoR;
			// }
			
			// // UnF = (MRM+MLP)*chat;

			// double PLR = 0.5*(pL+pR) - 
						 // 0.5*Mcy*preP*preM*0.5*(pL+pR)/chat*(UnR-UnL) + 
						 // Mcy*0.5*(pL+pR)*(preP+preM-1.0) - 
						 // 0.5*(preP-preM)*(pR-pL);
						 
						 
			// faces[id_UnF] = UnF;
			// faces[id_pF] = PLR;
				
				
				
				
		
		
				
				
				double phiL2[5+nSp], phiL1[5+nSp], phiR1[5+nSp];
				if(UnF>=0.0){
					phiL1[0] = uL; phiR1[0] = uR;
					// skewness
					phiL1[0] += (dudxL*LNv[0]+dudyL*LNv[1]+dudzL*LNv[2]);
					phiR1[0] += (dudxR*RNv[0]+dudyR*RNv[1]+dudzR*RNv[2]);
					phiL2[0] = 0.0;
					phiL2[0] += dudxL*dLR*nvec[0]; // phiL2[0] += dudxL*faces[id_xLR];
					phiL2[0] += dudyL*dLR*nvec[1]; // phiL2[0] += dudyL*faces[id_yLR];
					phiL2[0] += dudzL*dLR*nvec[2]; // phiL2[0] += dudzL*faces[id_zLR];
					phiL2[0] = phiR1[0] - 2.0*phiL2[0];
					// skewness
					phiL2[0] += 2.0*(dudxL*LNv[0]+dudyL*LNv[1]+dudzL*LNv[2]);
					
					phiL1[1] = vL; phiR1[1] = vR;
					// skewness
					phiL1[1] += (dvdxL*LNv[0]+dvdyL*LNv[1]+dvdzL*LNv[2]);
					phiR1[1] += (dvdxR*RNv[0]+dvdyR*RNv[1]+dvdzR*RNv[2]);
					phiL2[1] = 0.0;
					phiL2[1] += dvdxL*dLR*nvec[0];
					phiL2[1] += dvdyL*dLR*nvec[1];
					phiL2[1] += dvdzL*dLR*nvec[2];
					phiL2[1] = phiR1[1] - 2.0*phiL2[1];
					// skewness
					phiL2[1] += 2.0*(dvdxL*LNv[0]+dvdyL*LNv[1]+dvdzL*LNv[2]);
					
					phiL1[2] = wL; phiR1[2] = wR;
					// skewness
					phiL1[2] += (dwdxL*LNv[0]+dwdyL*LNv[1]+dwdzL*LNv[2]);
					phiR1[2] += (dwdxR*RNv[0]+dwdyR*RNv[1]+dwdzR*RNv[2]);
					phiL2[2] = 0.0;
					phiL2[2] += dwdxL*dLR*nvec[0];
					phiL2[2] += dwdyL*dLR*nvec[1];
					phiL2[2] += dwdzL*dLR*nvec[2];
					phiL2[2] = phiR1[2] - 2.0*phiL2[2];
					// skewness
					phiL2[2] += 2.0*(dwdxL*LNv[0]+dwdyL*LNv[1]+dwdzL*LNv[2]);
					
					phiL1[3] = TL; phiR1[3] = TR;
					// skewness
					phiL1[3] += (dTdxL*LNv[0]+dTdyL*LNv[1]+dTdzL*LNv[2]);
					phiR1[3] += (dTdxR*RNv[0]+dTdyR*RNv[1]+dTdzR*RNv[2]);
					phiL2[3] = 0.0;
					phiL2[3] += dTdxL*dLR*nvec[0];
					phiL2[3] += dTdyL*dLR*nvec[1];
					phiL2[3] += dTdzL*dLR*nvec[2];
					phiL2[3] = phiR1[3] - 2.0*phiL2[3];
					// skewness
					phiL2[3] += 2.0*(dTdxL*LNv[0]+dTdyL*LNv[1]+dTdzL*LNv[2]);
					
					for(int i=0; i<nSp-1; ++i){
						int tmp_i = 4+i;
						phiL1[tmp_i] = YL[i]; phiR1[tmp_i] = YR[i];
						// skewness
						phiL1[tmp_i] += (dYdxL[i]*LNv[0]+dYdyL[i]*LNv[1]+dYdzL[i]*LNv[2]);
						phiR1[tmp_i] += (dYdxR[i]*RNv[0]+dYdyR[i]*RNv[1]+dYdzR[i]*RNv[2]);
						phiL2[tmp_i] = 0.0;
						phiL2[tmp_i] += dYdxL[i]*dLR*nvec[0];
						phiL2[tmp_i] += dYdyL[i]*dLR*nvec[1];
						phiL2[tmp_i] += dYdzL[i]*dLR*nvec[2];
						phiL2[tmp_i] = phiR1[tmp_i] - 2.0*phiL2[tmp_i];
						// skewness
						phiL2[tmp_i] += 2.0*(dYdxL[i]*LNv[0]+dYdyL[i]*LNv[1]+dYdzL[i]*LNv[2]);
					}
					
					
				}
				else{
					phiL1[0] = uR; phiR1[0] = uL;
					// skewness
					phiL1[0] += (dudxR*RNv[0]+dudyR*RNv[1]+dudzR*RNv[2]);
					phiR1[0] += (dudxL*LNv[0]+dudyL*LNv[1]+dudzL*LNv[2]);
					phiL2[0] = 0.0;
					phiL2[0] += dudxR*dLR*nvec[0];
					phiL2[0] += dudyR*dLR*nvec[1];
					phiL2[0] += dudzR*dLR*nvec[2];
					phiL2[0] = phiR1[0] + 2.0*phiL2[0];
					// skewness
					phiL2[0] += 2.0*(dudxL*RNv[0]+dudyL*RNv[1]+dudzL*RNv[2]);
					
					phiL1[1] = vR; phiR1[1] = vL;
					// skewness
					phiL1[1] += (dvdxR*RNv[0]+dvdyR*RNv[1]+dvdzR*RNv[2]);
					phiR1[1] += (dvdxL*LNv[0]+dvdyL*LNv[1]+dvdzL*LNv[2]);
					phiL2[1] = 0.0;
					phiL2[1] += dvdxR*dLR*nvec[0];
					phiL2[1] += dvdyR*dLR*nvec[1];
					phiL2[1] += dvdzR*dLR*nvec[2];
					phiL2[1] = phiR1[1] + 2.0*phiL2[1];
					// skewness
					phiL2[1] += 2.0*(dvdxL*RNv[0]+dvdyL*RNv[1]+dvdzL*RNv[2]);
					
					phiL1[2] = wR; phiR1[2] = wL;
					// skewness
					phiL1[2] += (dwdxR*RNv[0]+dwdyR*RNv[1]+dwdzR*RNv[2]);
					phiR1[2] += (dwdxL*LNv[0]+dwdyL*LNv[1]+dwdzL*LNv[2]);
					phiL2[2] = 0.0;
					phiL2[2] += dwdxR*dLR*nvec[0];
					phiL2[2] += dwdyR*dLR*nvec[1];
					phiL2[2] += dwdzR*dLR*nvec[2];
					phiL2[2] = phiR1[2] + 2.0*phiL2[2];
					// skewness
					phiL2[2] += 2.0*(dwdxL*RNv[0]+dwdyL*RNv[1]+dwdzL*RNv[2]);
					
					phiL1[3] = TR; phiR1[3] = TL;
					// skewness
					phiL1[3] += (dTdxR*RNv[0]+dTdyR*RNv[1]+dTdzR*RNv[2]);
					phiR1[3] += (dTdxL*LNv[0]+dTdyL*LNv[1]+dTdzL*LNv[2]);
					phiL2[3] = 0.0;
					phiL2[3] += dTdxR*dLR*nvec[0];
					phiL2[3] += dTdyR*dLR*nvec[1];
					phiL2[3] += dTdzR*dLR*nvec[2];
					phiL2[3] = phiR1[3] + 2.0*phiL2[3];
					// skewness
					phiL2[3] += 2.0*(dTdxL*RNv[0]+dTdyL*RNv[1]+dTdzL*RNv[2]);
					
					for(int i=0; i<nSp-1; ++i){
						int tmp_i = 4+i;
						phiL1[tmp_i] = YR[i]; phiR1[tmp_i] = YL[i];
						// skewness
						phiL1[tmp_i] += (dYdxR[i]*RNv[0]+dYdyR[i]*RNv[1]+dYdzR[i]*RNv[2]);
						phiR1[tmp_i] += (dYdxL[i]*LNv[0]+dYdyL[i]*LNv[1]+dYdzL[i]*LNv[2]);
						phiL2[tmp_i] = 0.0;
						phiL2[tmp_i] += dYdxR[i]*dLR*nvec[0];
						phiL2[tmp_i] += dYdyR[i]*dLR*nvec[1];
						phiL2[tmp_i] += dYdzR[i]*dLR*nvec[2];
						phiL2[tmp_i] = phiR1[tmp_i] + 2.0*phiL2[tmp_i];
						// skewness
						phiL2[tmp_i] += 2.0*(dYdxL[i]*RNv[0]+dYdyL[i]*RNv[1]+dYdzL[i]*RNv[2]);
					}
				}
				
				// // faces[id_pF] = wdL*pL+wdR*pR;
				
				
				faces[id_uF] = solver.NVD.Minmod(phiL2[0],phiL1[0],phiR1[0]);
				faces[id_vF] = solver.NVD.Minmod(phiL2[1],phiL1[1],phiR1[1]);
				faces[id_wF] = solver.NVD.Minmod(phiL2[2],phiL1[2],phiR1[2]);
				faces[id_TF] = solver.NVD.Minmod(phiL2[3],phiL1[3],phiR1[3]);
				
				{
					double dx = sqrt(faces[id_area]);
					double dt_tmp_L = dx/(sqrt(uL*uL+vL*vL+wL*wL)+cL);
					double dt_tmp_R = dx/(sqrt(uR*uR+vR*vR+wR*wR)+cR);			
					double coDD = 1.0*max(dt/dt_tmp_L,dt/dt_tmp_R);
					
					for(int i=0; i<nSp-1; ++i){
						double mfLR[3];
						mfLR[0] = 0.5*dYdxL[i]+0.5*dYdxR[i];
						mfLR[1] = 0.5*dYdyL[i]+0.5*dYdyR[i];
						mfLR[2] = 0.5*dYdzL[i]+0.5*dYdzR[i];
						double magMfLR = mfLR[0]*mfLR[0];
						magMfLR += mfLR[1]*mfLR[1];
						magMfLR += mfLR[2]*mfLR[2];
						magMfLR = sqrt(magMfLR);
						mfLR[0] = mfLR[0]/(magMfLR+1.e-200);
						mfLR[1] = mfLR[1]/(magMfLR+1.e-200);
						mfLR[2] = mfLR[2]/(magMfLR+1.e-200);
						double cosTheta = mfLR[0]*faces[id_nx];
						cosTheta += mfLR[1]*faces[id_ny];
						cosTheta += mfLR[2]*faces[id_nz];
						cosTheta = abs(cosTheta);
						double gamF = min(cosTheta*cosTheta*cosTheta*cosTheta*cosTheta*cosTheta,1.0);
						
						// faces[id_YF[i]] = solver.NVD.Minmod(phiL2[4+i],phiL1[4+i],phiR1[4+i]);
						faces[id_YF[i]] = solver.NVD.getHO_MSTACS(phiL2[4+i],phiL1[4+i],phiR1[4+i],coDD,gamF);
						
						faces[id_YF[i]] = max(0.0,min(1.0,faces[id_YF[i]]));
					}
				}
				
				
				
				if(UnF>=0.0){
					// faces[id_poF] = pL;
					faces[id_uF] = uL; 
					faces[id_vF] = vL; 
					faces[id_wF] = wL; 
					faces[id_TF] = TL;
					// for(int i=0; i<nSp-1; ++i){
						// faces[id_YF[i]] = YL[i];
					// }
				}
				else{
					// faces[id_poF] = pR;
					faces[id_uF] = uR; 
					faces[id_vF] = vR; 
					faces[id_wF] = wR; 
					faces[id_TF] = TR;
					// for(int i=0; i<nSp-1; ++i){
						// faces[id_YF[i]] = YR[i];
					// }
				}
				
				
				
				return 0;
			}
		);
	}
	// cout << "BBBBBBBB" << endl;
	
}

