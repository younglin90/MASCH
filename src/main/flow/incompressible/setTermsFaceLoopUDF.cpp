
#include "../../../others/solvers.h"
// convective
void MASCH_Solver::setTermsCellLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pF = controls.getId_faceVar("pressure");

	int id_dp = controls.getId_cellVar("delta-pressure");
	int id_dpF = controls.getId_faceVar("delta-pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uF = controls.getId_faceVar("x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	

	double muF = 0.001;

	{
		calcConvFlux.push_back(
		[&solver,
		id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
		id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,muF](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			
			double uL = cellsL[id_u]; double uR = cellsR[id_u];
			double vL = cellsL[id_v]; double vR = cellsR[id_v];
			double wL = cellsL[id_w]; double wR = cellsR[id_w];
			
			double uF = faces[id_uF];
			double vF = faces[id_vF];
			double wF = faces[id_wF];
			double pF = faces[id_pF];
			double UnF = faces[id_UnF];
			double weiL = 1.0; double weiR = 0.0;
			if(UnF<0.0){
				weiL = 0.0; weiR = 1.0;
			}
			
			
			fluxA_LL[0] += weiL*UnF*area;    fluxA_LR[0] += weiR*UnF*area;
			fluxA_RR[0] += weiR*(-UnF)*area; fluxA_RL[0] += weiL*(-UnF)*area;
			
			fluxA_LL[4] += weiL*UnF*area;    fluxA_LR[4] += weiR*UnF*area;
			fluxA_RR[4] += weiR*(-UnF)*area; fluxA_RL[4] += weiL*(-UnF)*area;
			
			fluxA_LL[8] += weiL*UnF*area;    fluxA_LR[8] += weiR*UnF*area;
			fluxA_RR[8] += weiR*(-UnF)*area; fluxA_RL[8] += weiL*(-UnF)*area;
			
			fluxA_LL[0] += muF/dLR*area; fluxA_LR[0] -= muF/dLR*area;
			fluxA_RR[0] += muF/dLR*area; fluxA_RL[0] -= muF/dLR*area;
			
			fluxA_LL[4] += muF/dLR*area; fluxA_LR[4] -= muF/dLR*area;
			fluxA_RR[4] += muF/dLR*area; fluxA_RL[4] -= muF/dLR*area;
			
			fluxA_LL[8] += muF/dLR*area; fluxA_LR[8] -= muF/dLR*area;
			fluxA_RR[8] += muF/dLR*area; fluxA_RL[8] -= muF/dLR*area;
			
			fluxB[0] -= ( uF*UnF + pF*nvec[0] )*area;
			fluxB[1] -= ( vF*UnF + pF*nvec[1] )*area;
			fluxB[2] -= ( wF*UnF + pF*nvec[2] )*area;
			
			fluxB[0] += ( muF*(uR-uL)/dLR )*area;
			fluxB[1] += ( muF*(vR-vL)/dLR )*area;
			fluxB[2] += ( muF*(wR-wL)/dLR )*area;
		}); 
	}
	

	{
		calcConvFlux.push_back(
		[&solver,
		id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
		id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			
			double uF = faces[id_uF];
			double vF = faces[id_vF];
			double wF = faces[id_wF];
			double pF = faces[id_pF];
			double UnF = faces[id_UnF];
			double dtrho = fields[id_dt]/1.0;
			
			fluxA_LL[0] -= dtrho/dLR*area; fluxA_RR[0] -= dtrho/dLR*area;
			fluxA_LR[0] += dtrho/dLR*area; fluxA_RL[0] += dtrho/dLR*area;
			
			fluxB[0] += UnF*area;
		}); 
	}
	
	

	{
		calcConvFlux.push_back(
		[&solver,
		id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
		id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
		double* fields, double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			double area = faces[id_area];
			double dpF = faces[id_dpF];
			
			fluxB[0] += ( dpF*nvec[0] )*area;
			fluxB[1] += ( dpF*nvec[1] )*area;
			fluxB[2] += ( dpF*nvec[2] )*area;
			
			// fluxB[0] += ( 0.5*(cellsL[id_dp]+cellsR[id_dp])*nvec[0] )*area;
			// fluxB[1] += ( 0.5*(cellsL[id_dp]+cellsR[id_dp])*nvec[1] )*area;
			// fluxB[2] += ( 0.5*(cellsL[id_dp]+cellsR[id_dp])*nvec[2] )*area;
		}); 
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	calcConvFlux_BC.resize(3);
	
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
			string type = controls.boundaryMap["velocity"][bcName+".type"];
			if(type=="zeroGradient"){
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,muF](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {

					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
							
					double uL = cellsL[id_u];
					double vL = cellsL[id_v];
					double wL = cellsL[id_w];
			
					double pF = faces[id_pF];
					double uF = faces[id_uF];
					double vF = faces[id_vF];
					double wF = faces[id_wF];
					
					double UnF = uF*nvec[0]+vF*nvec[1]+wF*nvec[2];
					
					// double muF = 0.001;
					
					fluxB[0] -= ( uF*UnF + pF*nvec[0] )*area;
					fluxB[1] -= ( vF*UnF + pF*nvec[1] )*area;
					fluxB[2] -= ( wF*UnF + pF*nvec[2] )*area;
					
					fluxB[0] += ( muF*(uF-uL)/dLR )*area;
					fluxB[1] += ( muF*(vF-vL)/dLR )*area;
					fluxB[2] += ( muF*(wF-wL)/dLR )*area;
					
					return 0;
				});
			}
			else {
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,muF](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {

					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
					
					double uL = cellsL[id_u];
					double vL = cellsL[id_v];
					double wL = cellsL[id_w];
					
					double pF = faces[id_pF];
					double uF = faces[id_uF];
					double vF = faces[id_vF];
					double wF = faces[id_wF];
					// cout << uF << " " << vF << " " << wF << endl;
					
					double UnF = uF*nvec[0]+vF*nvec[1]+wF*nvec[2];
					
					// double muF = 1.0;
					
					fluxA_LL[0] += UnF*area; 
					fluxA_LL[4] += UnF*area; 
					fluxA_LL[8] += UnF*area; 
					
					fluxA_LL[0] += muF/dLR*area; 
					fluxA_LL[4] += muF/dLR*area; 
					fluxA_LL[8] += muF/dLR*area; 
					
					fluxB[0] -= ( uF*UnF + pF*nvec[0] )*area;
					fluxB[1] -= ( vF*UnF + pF*nvec[1] )*area;
					fluxB[2] -= ( wF*UnF + pF*nvec[2] )*area;
					
					fluxB[0] += ( muF*(uF-uL)/dLR )*area;
					fluxB[1] += ( muF*(vF-vL)/dLR )*area;
					fluxB[2] += ( muF*(wF-wL)/dLR )*area;
					
					return 0;
				});
			}
		}
		
		
		
		

		{
			int bc_id = 1;
			calcConvFlux_BC[bc_id].push_back(vector<conv_diff_Funct_BC_type>());
			string type = controls.boundaryMap["pressure"][bcName+".type"];
			if(type=="zeroGradient"){
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {

					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
					
					double pF = faces[id_pF];
					double uF = faces[id_uF];
					double vF = faces[id_vF];
					double wF = faces[id_wF];
					
					double UnF = uF*nvec[0]+vF*nvec[1]+wF*nvec[2];
					
					
					fluxB[0] += UnF*area;
					
					return 0;
				});
			}
			else {
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {

					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
					
					double pF = faces[id_pF];
					double uF = faces[id_uF];
					double vF = faces[id_vF];
					double wF = faces[id_wF];
					
					double UnF = uF*nvec[0]+vF*nvec[1]+wF*nvec[2];
					double dtrho = fields[id_dt]/1.0;
					
					fluxA_LL[0] -= dtrho/dLR*area;
					
					fluxB[0] += UnF*area;
					
					return 0;
				});
			}
		}
		
		
		
		
		

		{
			int bc_id = 2;
			calcConvFlux_BC[bc_id].push_back(vector<conv_diff_Funct_BC_type>());
			string type = controls.boundaryMap["pressure"][bcName+".type"];
			if(type=="zeroGradient"){
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {
					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
					
					// double dpF = cellsL[id_dp];
					double dpF = faces[id_dpF];
					
					fluxB[0] += dpF*nvec[0]*area;
					fluxB[1] += dpF*nvec[1]*area;
					fluxB[2] += dpF*nvec[2]*area;
					return 0;
				});
			}
			else {
				calcConvFlux_BC[bc_id].back().push_back(
				[&solver,
				id_p,id_pF,id_dp,id_dpF,id_u,id_uF,id_v,id_vF,id_w,id_wF,
				id_dpdx,id_dpdy,id_dpdz,id_UnF,id_dt,id_nx,id_ny,id_nz,id_area,id_dLR](
				double* fields, double* cellsL, double* faces, 
				double* fluxA_LL, double* fluxB) ->int {

					double nvec[3];
					nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
					double area = faces[id_area];
					double dLR = faces[id_dLR]; 
					
					double dpF = cellsL[id_dp];
					
					// fluxB[0] += 0.5*dpF*nvec[0]*area;
					// fluxB[1] += 0.5*dpF*nvec[1]*area;
					// fluxB[2] += 0.5*dpF*nvec[2]*area;
					return 0;
				});
			}
		}
	}
	
	
}