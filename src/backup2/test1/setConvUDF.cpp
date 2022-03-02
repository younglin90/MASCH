
#include "./solvers.h"
// convective
void MASCH_Solver::setConvFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	checkImplicitConvFlux.push_back(false);
	
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
		
		// 메쉬관련
		int id_nx = controls.faceVar["x unit normal"].id;
		int id_ny = controls.faceVar["y unit normal"].id;
		int id_nz = controls.faceVar["z unit normal"].id;
		int id_area = controls.faceVar["area"].id;
		
		calcConvFlux.push_back(
		[nSp,id_nx,id_ny,id_nz,id_area,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
		id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
		id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR](
		double* cellsL, double* cellsR, 
		double* faces, double* fluxA, double* fluxB) ->int {
			// double pL = faces[id_pL]; double pR = faces[id_pR];
			double uL = faces[id_uL]; double uR = faces[id_uR];
			double vL = faces[id_vL]; double vR = faces[id_vR];
			double wL = faces[id_wL]; double wR = faces[id_wR];
			double TL = faces[id_TL]; double TR = faces[id_TR];
			double rhoL = faces[id_rhoL]; double rhoR = faces[id_rhoR];
			double cL = faces[id_cL]; double cR = faces[id_cR];
			double HtL = faces[id_HtL]; double HtR = faces[id_HtR];
			double YL[nSp]; double YR[nSp];
			for(int i=0; i<nSp-1; ++i){
				YL[i] = faces[id_YL[i]]; YR[i] = faces[id_YR[i]];
			}
			
			double pL, pR, uCL, uCR, vCL, vCR, wCL, wCR;
			if(cellsR!=nullptr){
				pL = cellsL[id_p]; pR = cellsR[id_p];
				uCL = cellsL[id_u]; uCR = cellsR[id_u];
				vCL = cellsL[id_v]; vCR = cellsR[id_v];
				wCL = cellsL[id_w]; wCR = cellsR[id_w];
			}
			else{
				pL = faces[id_pL]; pR = faces[id_pR];
				uCL = faces[id_uL]; uCR = faces[id_uR];
				vCL = faces[id_vL]; vCR = faces[id_vR];
				wCL = faces[id_wL]; wCR = faces[id_wR];
			}
			
			double nvec[3];
			nvec[0] = faces[id_nx];
			nvec[1] = faces[id_ny];
			nvec[2] = faces[id_nz];
			// properties of Left
			double UnL = uCL*nvec[0] + vCL*nvec[1] + wCL*nvec[2];
			// properties of Right
			double UnR = uCR*nvec[0] + vCR*nvec[1] + wCR*nvec[2];
			
			double unhat = (UnL+UnR)/2.0;
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

			double preLs = abs(pL) + 0.1 * rhoL*cL*cL;
			double preRs = abs(pR) + 0.1 * rhoR*cR*cR;
			double w = min(preLs/preRs,preRs/preLs);
			w = 1.0 - w*w;

			double MLPL = 0.5*(D_L+D_rho);
			double MRMR = 0.5*(D_R-D_rho);

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
	
	
	
	
	
}