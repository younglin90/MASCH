
#include "../../../others/solvers.h"
// convective
void MASCH_Solver::setConvFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	// 셀값
	int id_dp = controls.getId_cellVar("delta-pressure");
	int id_p = controls.getId_cellVar("pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_rho = controls.getId_cellVar("density");
	
	// left값
	int id_dpL = controls.getId_faceVar("left delta-pressure");
	int id_pL = controls.getId_faceVar("left pressure");
	int id_uL = controls.getId_faceVar("left x-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity");
	int id_rhoL = controls.getId_faceVar("left density");
	
	// right값
	int id_dpR = controls.getId_faceVar("right delta-pressure");
	int id_pR = controls.getId_faceVar("right pressure");
	int id_uR = controls.getId_faceVar("right x-velocity");
	int id_vR = controls.getId_faceVar("right y-velocity");
	int id_wR = controls.getId_faceVar("right z-velocity");
	int id_rhoR = controls.getId_faceVar("right density");
	
	// 메쉬관련
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");

	
	
	// checkImplicitConvFlux.push_back(true);
	{
		
		calcConvFlux.push_back(
		[id_dp,id_dpL,id_dpR,
		id_nx,id_ny,id_nz,id_area,id_dLR,
		id_p,id_u,id_v,id_w,id_rho,
		id_pL,id_uL,id_vL,id_wL,id_rhoL,
		id_pR,id_uR,id_vR,id_wR,id_rhoR](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			
			double dpL = cellsL[id_dp]; double dpR = 0.0;
			double uL = cellsL[id_u]; double uR = 0.0;
			double vL = cellsL[id_v]; double vR = 0.0;
			double wL = cellsL[id_w]; double wR = 0.0;
			double rhoL = cellsL[id_rho]; double rhoR = 0.0;
			if(cellsR!=nullptr){
				uR = cellsR[id_u];
				vR = cellsR[id_v];
				wR = cellsR[id_w];
				rhoR = cellsR[id_rho];
			}
			else{
				dpL = faces[id_dpL]; dpR = faces[id_dpR];
				uL = faces[id_uL]; uR = faces[id_uR];
				vL = faces[id_vL]; vR = faces[id_vR];
				wL = faces[id_wL]; wR = faces[id_wR];
				rhoL = faces[id_rhoL]; rhoR = faces[id_rhoR];
			}
			
			double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
			double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
			
			double area = faces[id_area];
			double avg_rho = 0.5*(1.0/rhoL+1.0/rhoR);
			double dLR = faces[id_dLR]; 
			
			// // zeroGradient
			// if(cellsR!=nullptr){
				// if(abs(dpL-dpR) < 1.e-200) avg_rho = 0.0;
				// // avg_rho = 0.0;
			// }
			fluxA_LL[0] -= avg_rho/dLR*area;
			fluxA_RR[0] -= avg_rho/dLR*area;
			fluxA_LR[0] += avg_rho/dLR*area;
			fluxA_RL[0] += avg_rho/dLR*area;
			
			fluxB[0] += 0.5*(UnL+UnR)*area;
			
		}); 
	}
	
	
	
	
	// checkImplicitConvFlux.push_back(false);
	{
		calcConvFlux.push_back(
		[id_nx,id_ny,id_nz,id_area,id_dLR,
		id_dp,id_dpL,id_dpR](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double dpL = cellsL[id_dp]; double dpR = 0.0;
			if(cellsR!=nullptr){
				dpR = cellsR[id_dp];
			}
			else{
				dpL = faces[id_dpL]; dpR = faces[id_dpR];
			}
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			
			double area = faces[id_area];
			fluxB[0] += 0.5*(dpL+dpR)*nvec[0]*area;
			fluxB[1] += 0.5*(dpL+dpR)*nvec[1]*area;
			fluxB[2] += 0.5*(dpL+dpR)*nvec[2]*area;
			
		
		}); 
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	using conv_diff_Funct_BC_type = 
		function<int(
		double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB)>;
	
	calcConvFlux_BC.resize(2);

	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcConvFlux_BC[0].push_back(vector<conv_diff_Funct_BC_type>());
		
		string name = "pressure";
		string type = controls.boundaryMap[name][bcName+".type"];
		
		if(type=="fixedValue"){
			double value = stod(controls.boundaryMap[name][bcName+".value"]);
			calcConvFlux_BC[0].back().push_back(
			[id_dpL,id_dpR,
			id_nx,id_ny,id_nz,id_area,id_dLR,
			id_pL,id_uL,id_vL,id_wL,id_rhoL,
			id_pR,id_uR,id_vR,id_wR,id_rhoR](
			double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB) ->int {

				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				
				double dpL = faces[id_dpL]; double dpR = faces[id_dpR];
				double uL = faces[id_uL]; double uR = faces[id_uR];
				double vL = faces[id_vL]; double vR = faces[id_vR];
				double wL = faces[id_wL]; double wR = faces[id_wR];
				double rhoL = faces[id_rhoL]; double rhoR = faces[id_rhoR];
				
				double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
				double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
				
				double area = faces[id_area];
				double avg_rho = 0.5*(1.0/rhoL+1.0/rhoR);
				double dLR = faces[id_dLR]; 
				
				fluxA_LL[0] -= avg_rho/dLR*area;
				
				fluxB[0] += 0.5*(UnL+UnR)*area;
				
				return 0;
			});
		}
		else if(type=="zeroGradient"){
			calcConvFlux_BC[0].back().push_back(
			[id_dpL,id_dpR,
			id_nx,id_ny,id_nz,id_area,id_dLR,
			id_pL,id_uL,id_vL,id_wL,id_rhoL,
			id_pR,id_uR,id_vR,id_wR,id_rhoR](
			double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB) ->int {

				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				
				double dpL = faces[id_dpL]; double dpR = faces[id_dpR];
				double uL = faces[id_uL]; double uR = faces[id_uR];
				double vL = faces[id_vL]; double vR = faces[id_vR];
				double wL = faces[id_wL]; double wR = faces[id_wR];
				double rhoL = faces[id_rhoL]; double rhoR = faces[id_rhoR];
				
				double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
				double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
				
				double area = faces[id_area];
				double avg_rho = 0.5*(1.0/rhoL+1.0/rhoR);
				double dLR = faces[id_dLR]; 
				
				// fluxA_LL[0] -= avg_rho/dLR*area;
				
				fluxB[0] += 0.5*(UnL+UnR)*area;
				
				return 0;
			});
		}
	}
	
	
	
	
	
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcConvFlux_BC[1].push_back(vector<conv_diff_Funct_BC_type>());
		
		// string type = controls.boundaryMap["velocity"][bcName+".type"];
			
		calcConvFlux_BC[1].back().push_back(
		[id_dpL,id_dpR,
		id_nx,id_ny,id_nz,id_area,id_dLR,
		id_pL,id_uL,id_vL,id_wL,id_rhoL,
		id_pR,id_uR,id_vR,id_wR,id_rhoR](
		double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB) ->int {

			double dpL = faces[id_dpL]; double dpR = faces[id_dpR];
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			
			double area = faces[id_area];
			fluxB[0] += 0.5*(dpL+dpR)*nvec[0]*area;
			fluxB[1] += 0.5*(dpL+dpR)*nvec[1]*area;
			fluxB[2] += 0.5*(dpL+dpR)*nvec[2]*area;
			
		}); 
		
	}
	
	
	
}