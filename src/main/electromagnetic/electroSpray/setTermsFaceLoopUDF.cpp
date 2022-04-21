
#include "../../../../others/solvers.h"
// convective
void MASCH_Solver::setTermsFaceLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	int nSp = controls.spName.size();
	int nEq = 5+nSp-1;
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_p = controls.getId_cellVar("pressure");
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	int id_dp = controls.getId_cellVar("delta-pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	vector<int> id_alpha, id_alphaL, id_alphaR;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+controls.spName[i]));
		id_alphaL.push_back(controls.getId_cellVar("left volume-fraction-"+controls.spName[i]));
		id_alphaR.push_back(controls.getId_cellVar("right volume-fraction-"+controls.spName[i]));
	}
	int id_mu = controls.getId_cellVar("viscosity");
	int id_rho = controls.getId_cellVar("density");
	
	int id_rhoe = controls.getId_cellVar("charge-density");
	int id_phi = controls.getId_cellVar("electric-potential");
	int id_xE = controls.getId_cellVar("x-electric-field");
	int id_yE = controls.getId_cellVar("y-electric-field");
	int id_zE = controls.getId_cellVar("z-electric-field");
	int id_xFe = controls.getId_cellVar("x-electric-force");
	int id_yFe = controls.getId_cellVar("y-electric-force");
	int id_zFe = controls.getId_cellVar("z-electric-force");
	int id_k = controls.getId_cellVar("conductivity");
	int id_epsilon = controls.getId_cellVar("permittivity");
    
	int id_Un = controls.getId_faceVar("contravariant-velocity");
	
	// 최대 최소값
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	int id_wd = controls.getId_faceVar("distance weight"); 
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
	
	
	
	


	{
		
		//=================================================
		// 1번째
		calcConvFlux.push_back(
		[nSp,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
		id_alphaL,id_alphaR](
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
			
			double UnF = faces[id_Un];
			
			// upwind weighting
			double weiL = 1.0; double weiR = 0.0;
			if(UnF<0.0){
				weiL = 0.0; weiR = 1.0;
			}
			
			double alphaF[nSp];
			for(int i=0; i<nSp-1; ++i){
				double alphaL = faces[id_alphaL[i]];
				double alphaR = faces[id_alphaR[i]];
				alphaF[i] = weiL*alphaL + weiR*alphaR;
			}
			
			int iter = 0;
			for(int i=0; i<nSp-1; ++i){
				for(int j=0; j<nSp-1; ++j){
					if(i==j){
						fluxA_LL[iter] += (UnF*weiL)*area; fluxA_LR[iter] += (UnF*weiR)*area;
						fluxA_RR[iter] -= (UnF*weiR)*area; fluxA_RL[iter] -= (UnF*weiL)*area;
					}
					++iter;
				}
			}
			
			for(int i=0; i<nSp-1; ++i){
				double flux = -alphaF[i]*UnF*area;
				fluxB_LL[i] += flux; fluxB_RR[i] -= flux;
			}
			
		}); 
		
		
		
		
		//=================================================
		// 2번째
		calcConvFlux.push_back(
		[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
		id_rhoe,id_k,id_xE,id_yE,id_zE](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double UnF = faces[id_Un];
            
			// upwind weighting
			double weiL = 1.0; double weiR = 0.0;
			if(UnF<0.0){
				weiL = 0.0; weiR = 1.0;
			}
			
			fluxA_LL[0] += (UnF*weiL)*area; fluxA_LR[0] += (UnF*weiR)*area;
			fluxA_RR[0] -= (UnF*weiR)*area; fluxA_RL[0] -= (UnF*weiL)*area;
			
			
			double rhoeL = cellsL[id_rhoe]; double rhoeR = cellsR[id_rhoe];
			double kL = cellsL[id_k]; double kR = cellsR[id_k];
			double xEL = cellsL[id_xE]; double xER = cellsR[id_xE];
			double yEL = cellsL[id_yE]; double yER = cellsR[id_yE];
			double zEL = cellsL[id_zE]; double zER = cellsR[id_zE];
			
			double rhoeF = weiL*rhoeL + weiR*rhoeR;
			
			double EL = xEL*nvec[0] + yEL*nvec[1] + zEL*nvec[2];
			double ER = xER*nvec[0] + yER*nvec[1] + zER*nvec[2];
			double EF = 0.5*(EL+ER);
			
			// upwind 기법
			double kF = kL;
			if(EF<0.0) kF = kR;
			
			double flux = -(rhoeF*UnF + kF*EF)*area;
			fluxB_LL[0] += flux; fluxB_RR[0] -= flux;
			
		}); 
		
		
		
		
		//=================================================
		// 3번째
		calcConvFlux.push_back(
		[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_phi,id_epsilon](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double phiL = cellsL[id_phi]; double phiR = cellsR[id_phi];
			double epsilonL = cellsL[id_epsilon]; double epsilonR = cellsR[id_epsilon];
			
			double epsilonF = 0.5*(epsilonL + epsilonR);
			
			double vflux_coeff = (epsilonF/dLR)*area;
			
			fluxA_LL[0] += (-vflux_coeff); fluxA_LR[0] += (+vflux_coeff);
			fluxA_RR[0] -= (+vflux_coeff); fluxA_RL[0] -= (-vflux_coeff);
			
			double vflux = -(epsilonF*(phiR-phiL)/dLR)*area;
			
			fluxB_LL[0] += vflux; fluxB_RR[0] -= vflux;
			
		}); 
		
		
		
		
		
		
		//=================================================
		// 4번째
		calcConvFlux.push_back(
		[](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double phiL = cellsL[id_phi]; double phiR = cellsR[id_phi];
			
			double phiF = 0.5*(phiL+phiR);
			
			double flux = -phiF*area;
			
			fluxB_LL[0] += flux; fluxB_RR[0] -= flux;
			
		}); 
		
		
		
		
		
		
		//=================================================
		// 5번째
		calcConvFlux.push_back(
		[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_phi,id_epsilon,id_xE,id_yE,id_zE](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double phiL = cellsL[id_phi]; double phiR = cellsR[id_phi];
			double epsilonL = cellsL[id_epsilon]; double epsilonR = cellsR[id_epsilon];
			double xEL = cellsL[id_xE]; double xER = cellsR[id_xE];
			double yEL = cellsL[id_yE]; double yER = cellsR[id_yE];
			double zEL = cellsL[id_zE]; double zER = cellsR[id_zE];
			
			double EL = xEL*nvec[0] + yEL*nvec[1] + zEL*nvec[2];
			double ER = xER*nvec[0] + yER*nvec[1] + zER*nvec[2];
			
			double epsilonF = 0.5*(epsilonL + epsilonR);
			double EF = 0.5*(EL+ER);
			double E2L = xEL*xEL + yEL*yEL + zEL*zEL;
			double E2R = xER*xER + yER*yER + zER*zER;
			double E2 = 0.5*(E2L+E2R);
			
			double xEF = xEL; double yEF = yEL; double zEF = zEL;
			if(EF<0.0){
				xEF = xER; yEF = yER; zEF = zER;
			}
			
			double xflux = epsilonF*(xEF*EF-0.5*E2*nvec[0])*area;
			double yflux = epsilonF*(yEF*EF-0.5*E2*nvec[1])*area;
			double zflux = epsilonF*(zEF*EF-0.5*E2*nvec[2])*area;
			
			fluxB_LL[0] += xflux; fluxB_RR[0] -= xflux;
			fluxB_LL[1] += yflux; fluxB_RR[1] -= yflux;
			fluxB_LL[2] += zflux; fluxB_RR[2] -= zflux;
			
		}); 
		
		
		
		
		//=================================================
		// 6번째
		calcConvFlux.push_back(
		[nSp,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
		id_mu,surf_sigma,id_kappa,nCurv,id_alpha_VF](
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
			
			double pL = cellsL[id_p]; double pR = cellsR[id_p];
			double uL = cellsL[id_u]; double uR = cellsR[id_u];
			double vL = cellsL[id_v]; double vR = cellsR[id_v];
			double wL = cellsL[id_w]; double wR = cellsR[id_w];
			double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
			double muL = cellsL[id_mu]; double muR = cellsR[id_mu];
			// double rhoFL = faces[id_rhoL]; double rhoFR = faces[id_rhoR];
			double alphaL[nSp], alphaR[nSp];
			for(int i=0; i<nSp-1; ++i){
				alphaL[i] = cellsL[id_alpha[i]];
				alphaR[i] = cellsR[id_alpha[i]];
			}
			
			double UnF = faces[id_Un];
			
			// upwind weighting
			double weiL = 1.0; double weiR = 0.0;
			double uF = uL; double vF = vL; double wF = wL;
			if(UnF<0.0){
				uF = uR; vF = vR; wF = wR;
				weiL = 0.0; weiR = 1.0;
			}
			
			double pF = 0.5*(pL+pR);
			double muF = 0.5*(muL+muR);
			
			
			int iter=0;
			
			double xfluxL = -(rhoL*uF*UnF + pF*nvec[0])*area;
			double yfluxL = -(rhoL*vF*UnF + pF*nvec[1])*area;
			double zfluxL = -(rhoL*wF*UnF + pF*nvec[2])*area;
			double xfluxR = -(rhoR*uF*UnF + pF*nvec[0])*area;
			double yfluxR = -(rhoR*vF*UnF + pF*nvec[1])*area;
			double zfluxR = -(rhoR*wF*UnF + pF*nvec[2])*area;
			
			
			iter = nEq*0+0;
			fluxA_LL[iter] += (rhoL*weiL*UnF*area); fluxA_LR[iter] += (rhoL*weiR*UnF*area);
			fluxA_RR[iter] -= (rhoR*weiR*UnF*area); fluxA_RL[iter] -= (rhoR*weiL*UnF*area);
			
			iter = nEq*1+1;
			fluxA_LL[iter] += (rhoL*weiL*UnF*area); fluxA_LR[iter] += (rhoL*weiR*UnF*area);
			fluxA_RR[iter] -= (rhoR*weiR*UnF*area); fluxA_RL[iter] -= (rhoR*weiL*UnF*area);
			
			iter = nEq*2+2;
			fluxA_LL[iter] += (rhoL*weiL*UnF*area); fluxA_LR[iter] += (rhoL*weiR*UnF*area);
			fluxA_RR[iter] -= (rhoR*weiR*UnF*area); fluxA_RL[iter] -= (rhoR*weiL*UnF*area);
			
			fluxB_LL[0] += (xfluxL); fluxB_RR[0] -= (xfluxR);
			fluxB_LL[1] += (yfluxL); fluxB_RR[1] -= (yfluxR);
			fluxB_LL[2] += (zfluxL); fluxB_RR[2] -= (zfluxR);
			
			
			double vflux_coeff = (muF/dLR)*area;
			double vxflux = (muF*(uR-uL)/dLR)*area;
			double vyflux = (muF*(vR-vL)/dLR)*area;
			double vzflux = (muF*(wR-wL)/dLR)*area;
			
			iter = nEq*0+0;
			fluxA_LL[iter] += (+vflux_coeff); fluxA_LR[iter] += (-vflux_coeff);
			fluxA_RR[iter] -= (-vflux_coeff); fluxA_RL[iter] -= (+vflux_coeff);
			
			iter = nEq*1+1;
			fluxA_LL[iter] += (+vflux_coeff); fluxA_LR[iter] += (-vflux_coeff);
			fluxA_RR[iter] -= (-vflux_coeff); fluxA_RL[iter] -= (+vflux_coeff);
			
			iter = nEq*2+2;
			fluxA_LL[iter] += (+vflux_coeff); fluxA_LR[iter] += (-vflux_coeff);
			fluxA_RR[iter] -= (-vflux_coeff); fluxA_RL[iter] -= (+vflux_coeff);
			
			fluxB_LL[0] += (vxflux); fluxB_RR[0] -= (vxflux);
			fluxB_LL[1] += (vyflux); fluxB_RR[1] -= (vyflux);
			fluxB_LL[2] += (vzflux); fluxB_RR[2] -= (vzflux);
			
			
			// 표면장력
			for(int i=0; i<nCurv; ++i){
				double alphaF = 0.5*(cellsL[id_alpha_VF[i]]+cellsR[id_alpha_VF[i]]);
				double curvatureL = cellsL[id_kappa];
				double curvatureR = cellsR[id_kappa];
				fluxB_LL[1] += surf_sigma[i]*curvatureL*( alphaF*nvec[0] )*area;
				fluxB_LL[2] += surf_sigma[i]*curvatureL*( alphaF*nvec[1] )*area;
				fluxB_LL[3] += surf_sigma[i]*curvatureL*( alphaF*nvec[2] )*area;
				
				fluxB_RR[1] -= surf_sigma[i]*curvatureR*( alphaF*nvec[0] )*area;
				fluxB_RR[2] -= surf_sigma[i]*curvatureR*( alphaF*nvec[1] )*area;
				fluxB_RR[3] -= surf_sigma[i]*curvatureR*( alphaF*nvec[2] )*area;
			}
			
			
		}); 
		
		
		
		
		
		
		//=================================================
		// 7번째
		calcConvFlux.push_back(
		[id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho](
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
            
			double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
			
			double UnF = faces[id_Un];
			
			// upwind weighting
			double weiL = 1.0; double weiR = 0.0;
			double uF = uL; double vF = vL; double wF = wL;
			if(UnF<0.0){
				uF = uR; vF = vR; wF = wR;
				weiL = 0.0; weiR = 1.0;
			}
			
			double vflux_coeff = (dt*0.5*(1.0/rhoL+1.0/rhoR)/dLR)*area;
			
			fluxA_LL[0] += (-vflux_coeff); fluxA_LR[0] += (+vflux_coeff);
			fluxA_RR[0] -= (+vflux_coeff); fluxA_RL[0] -= (-vflux_coeff);
			
			
			double flux = UnF*area;
			
			fluxB_LL[0] += flux; fluxB_RR[0] -= flux;
			
			
		}); 
		
		
		
		
		//=================================================
		// 8번째
		calcConvFlux.push_back(
		[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
		id_nLRx,id_nLRy,id_nLRz,
		id_dp](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, 
		double* fluxB_LL, double* fluxB_RR) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			double dAlpha = faces[id_alpha];
			double nLR[3];
			nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
			
			double dpL = cellsL[id_dp]; double dpR = cellsR[id_dp];
			
			double dpF = 0.5*(dpL+dpR);
			
			double xflux = -dpF*nvec[0]*area;
			double yflux = -dpF*nvec[1]*area;
			double zflux = -dpF*nvec[2]*area;
			
			fluxB_LL[0] += xflux; fluxB_RR[0] -= xflux;
			fluxB_LL[1] += yflux; fluxB_RR[1] -= yflux;
			fluxB_LL[2] += zflux; fluxB_RR[2] -= zflux;
			
			
		}); 
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	calcConvFlux_BC.resize(8);
	
	using conv_diff_Funct_BC_type = 
		function<int(
		double* fields, double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB)>;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		
		{
			//=================================================
			// 1번째
			calcConvFlux_BC[0].push_back(vector<conv_diff_Funct_BC_type>());
			
			int nEq = controls.spName.size()-1;
			vector<bool> bool_zeroGradient(nEq,false);
			for(int i=0; i<controls.spName.size()-1; ++i){
				string tmp_name = ("volume-fraction-"+controls.spName[i]);
				if(controls.boundaryMap[tmp_name][bcName+".type"]=="zeroGradient") bool_zeroGradient[i] = true;
			}
				
			calcConvFlux_BC[0].back().push_back(
			[nSp,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
			id_alphaF,
			bool_zeroGradient](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				
				double uF = cellsL[id_uF];// double uR = cellsR[id_u];
				double vF = cellsL[id_vF];// double vR = cellsR[id_v];
				double wF = cellsL[id_wF];// double wR = cellsR[id_w];
				
				double UnF = uF*nvec[0] + vF*nvec[1] + wF*nvec[2];
				
				double alphaF[nSp];
				for(int i=0; i<nSp-1; ++i){
					alphaF[i] = faces[id_alphaF[i]];
				}
				
				int iter = 0;
				for(int i=0; i<nSp-1; ++i){
					for(int j=0; j<nSp-1; ++j){
						if(i==j){
							if(bool_zeroGradient[i]==true) 
								fluxA_LL[iter] += (UnF)*area;
						}
						++iter;
					}
				}
				
				
				for(int i=0; i<nSp-1; ++i){
					double flux = -alphaF[i]*UnF*area;
					fluxB_LL[i] += flux;
				}
				
			}); 
		}
			
			
		{
			//=================================================
			// 2번째
			calcConvFlux_BC[1].push_back(vector<conv_diff_Funct_BC_type>());
			int nEq = 1;
			vector<bool> bool_zeroGradient(nEq,false);
			if(controls.boundaryMap["charge-density"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			
			calcConvFlux_BC[1].back().push_back(
			[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
			id_rhoe,id_k,id_xE,id_yE,id_zE,
			bool_zeroGradient](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				
				double uF = cellsL[id_uF];// double uR = cellsR[id_u];
				double vF = cellsL[id_vF];// double vR = cellsR[id_v];
				double wF = cellsL[id_wF];// double wR = cellsR[id_w];
				
				double UnF = uF*nvec[0] + vF*nvec[1] + wF*nvec[2];
				
				if(bool_zeroGradient[0]==true)
					fluxA_LL[0] += (UnF)*area;
				
				
				double rhoeF = faces[id_rhoeF];// double rhoeR = cellsR[id_rhoe];
				double kF = faces[id_kF];// double kR = cellsR[id_k];
				// double xEF = faces[id_xEF];// double xER = cellsR[id_xE];
				// double yEF = faces[id_yEF];// double yER = cellsR[id_yE];
				// double zEF = faces[id_zEF];// double zER = cellsR[id_zE];
				
				double rhoeF = rhoeF;
				
				double EF = -(faces[id_phiF]-cellsL[id_phi])/dLR;
				
				// upwind 기법
				double kF = faces[id_kF];
				
				double flux = -(rhoeF*UnF + kF*EF)*area;
				fluxB_LL[0] += flux;
				
			}); 
		}
			
			
		{
			//=================================================
			// 3번째
			calcConvFlux_BC[2].push_back(vector<conv_diff_Funct_BC_type>());
			int nEq = 1;
			vector<bool> bool_zeroGradient(nEq,false);
			if(controls.boundaryMap["charge-potential"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			
			calcConvFlux_BC[2].back().push_back(
			[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_phi,id_epsilon,
			bool_zeroGradient](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				
				double phiL = cellsL[id_phi]; double phiF = faces[id_phiF];
				
				double epsilonF = faces[id_epsilonF];
				
				double vflux_coeff = (epsilonF/dLR)*area;
				
				if(bool_zeroGradient[0]==false)
					fluxA_LL[0] += (-vflux_coeff);
				
				double vflux = -(epsilonF*(phiF-phiL)/dLR)*area;
				
				fluxB_LL[0] += vflux;
				
			}); 
			
		}
			
			
			
		{
			//=================================================
			// 4번째
			calcConvFlux_BC[3].push_back(vector<conv_diff_Funct_BC_type>());
			
			calcConvFlux_BC[3].back().push_back(
			[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_phi,id_epsilon](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				
				double phiF = faces[id_phiF];
				
				double flux = -phiF*area;
				
				fluxB_LL[0] += flux; fluxB_RR[0] -= flux;
				
			}); 
			
		}
			
			
			
		{
			//=================================================
			// 5번째
			calcConvFlux_BC[4].push_back(vector<conv_diff_Funct_BC_type>());
			
			calcConvFlux_BC[4].back().push_back(
			[id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_phi,id_epsilon,id_xE,id_yE,id_zE](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				
				// double xEF = faces[id_xEF];
				// double yEF = faces[id_yEF];
				// double zEF = faces[id_zEF];
				
				double epsilonF = faces[id_epsilonF];
				double EF = -(faces[id_phiF]-cellsL[id_phi])/dLR;
				double xEF = EF*nvec[0];
				double yEF = EF*nvec[1];
				double zEF = EF*nvec[2];
				double E2 = xEF*xEF + yEF*yEF + zEF*zEF;
				
				double xflux = epsilonF*(xEF*EF-0.5*E2*nvec[0])*area;
				double yflux = epsilonF*(yEF*EF-0.5*E2*nvec[1])*area;
				double zflux = epsilonF*(zEF*EF-0.5*E2*nvec[2])*area;
				
				fluxB_LL[0] += xflux;
				fluxB_LL[1] += yflux;
				fluxB_LL[2] += zflux;
				
			}); 
		}
			
			
		{
			//=================================================
			// 6번째
			calcConvFlux_BC[5].push_back(vector<conv_diff_Funct_BC_type>());
			int nEq = 1;
			vector<bool> bool_zeroGradient(nEq,false);
			if(controls.boundaryMap["velocity"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			
			calcConvFlux_BC[5].back().push_back(
			[nSp,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
			id_mu,surf_sigma,id_kappa,nCurv,id_alpha_VF,
			bool_zeroGradient](
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
				
				double uL = cellsL[id_u]; double uF = faces[id_uF];
				double vL = cellsL[id_v]; double vF = faces[id_vF];
				double wL = cellsL[id_w]; double wF = faces[id_wF];
				double rhoL = cellsL[id_rho];
				
				double UnF = uF*nvec[0] + vF*nvec[1] + wF*nvec[2];
				
				// upwind
				double pF = faces[id_pF];
				double muF = faces[id_muF];
				
				
				int iter=0;
				
				double xfluxL = -(rhoL*uF*UnF + pF*nvec[0])*area;
				double yfluxL = -(rhoL*vF*UnF + pF*nvec[1])*area;
				double zfluxL = -(rhoL*wF*UnF + pF*nvec[2])*area;
				
				if(bool_zeroGradient[0]==true){
					iter = nEq*0+0;
					fluxA_LL[iter] += (rhoL*UnF*area);
					
					iter = nEq*1+1;
					fluxA_LL[iter] += (rhoL*UnF*area);
					
					iter = nEq*2+2;
					fluxA_LL[iter] += (rhoL*UnF*area);
				}
				
				fluxB_LL[0] += (xfluxL);
				fluxB_LL[1] += (yfluxL);
				fluxB_LL[2] += (zfluxL);
				
				
				double vflux_coeff = (muF/dLR)*area;
				double vxflux = (muF*(uF-uL)/dLR)*area;
				double vyflux = (muF*(vF-vL)/dLR)*area;
				double vzflux = (muF*(wF-wL)/dLR)*area;
				
				
				if(bool_zeroGradient[0]==false){
					iter = nEq*0+0;
					fluxA_LL[iter] += (+vflux_coeff);
					
					iter = nEq*1+1;
					fluxA_LL[iter] += (+vflux_coeff);
					
					iter = nEq*2+2;
					fluxA_LL[iter] += (+vflux_coeff);
				}
				
				
				fluxB_LL[0] += (vxflux);
				fluxB_LL[1] += (vyflux);
				fluxB_LL[2] += (vzflux);
				
				
				// 표면장력
				for(int i=0; i<nCurv; ++i){
					double alphaF = faces[id_alpha_VFF[i]];
					double curvatureL = cellsL[id_kappa];
					fluxB_LL[1] += surf_sigma[i]*curvatureL*( alphaF*nvec[0] )*area;
					fluxB_LL[2] += surf_sigma[i]*curvatureL*( alphaF*nvec[1] )*area;
					fluxB_LL[3] += surf_sigma[i]*curvatureL*( alphaF*nvec[2] )*area;
				}
				
				
			}); 
		}
			
			
			
			
		{
			//=================================================
			// 7번째
			calcConvFlux_BC[6].push_back(vector<conv_diff_Funct_BC_type>());
			int nEq = 1;
			vector<bool> bool_zeroGradient(nEq,false);
			if(controls.boundaryMap["delta-pressure"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;
			
			calcConvFlux_BC[6].back().push_back(
			[id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_p,id_u,id_v,id_w,id_dpdx,id_dpdy,id_dpdz,id_rho,
			bool_zeroGradient](
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
				
				double uL = cellsL[id_u]; double uF = faces[id_uF];
				double vL = cellsL[id_v]; double vF = faces[id_vF];
				double wL = cellsL[id_w]; double wF = faces[id_wF];
				double rhoL = cellsL[id_rho];
				
				double UnF = uF*nvec[0] + vF*nvec[1] + wF*nvec[2];
				
				double vflux_coeff = (dt*(1.0/rhoL)/dLR)*area;
				
				if(bool_zeroGradient[0]==false){
					fluxA_LL[0] += (-vflux_coeff);
				}
				
				double flux = UnF*area;
				
				fluxB_LL[0] += flux;
				
				
			}); 
		}
			
			
		{
			//=================================================
			// 8번째
            calcConvFlux_BC[7].push_back(vector<conv_diff_Funct_BC_type>());
            int nEq = 1;
            vector<bool> bool_zeroGradient(nEq,false);
            if(controls.boundaryMap["delta-pressure"][bcName+".type"]=="zeroGradient") bool_zeroGradient[0] = true;

			calcConvFlux_BC[7].back().push_back(
            [id_nx,id_ny,id_nz,id_area,id_dLR,id_alpha,
			id_nLRx,id_nLRy,id_nLRz,
			id_dp,
			bool_zeroGradient](
			double* fields, double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB_LL) ->int {
				
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				
				double dpL = cellsL[id_dp];
				
				double dpF = 0.0;
				if(bool_zeroGradient[0]==true){
					dpF = dpL;
				}
				else{
					dpF = 0.5*dpL;
				}
				
				double xflux = -dpF*nvec[0]*area;
				double yflux = -dpF*nvec[1]*area;
				double zflux = -dpF*nvec[2]*area;
				
				fluxB_LL[0] += xflux;
				fluxB_LL[1] += yflux;
				fluxB_LL[2] += zflux;
				
				
			}); 
		}
		
	}
		
	
}





