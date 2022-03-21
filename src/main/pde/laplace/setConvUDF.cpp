
#include "../../../others/solvers.h"
// convective
void MASCH_Solver::setConvFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	
	// 셀값
	int id_phi = controls.getId_cellVar("unknown");
	// int id_dphi = controls.getId_cellVar("delta-unknown");
	
	// left값
	int id_phiL = controls.getId_faceVar("left unknown");
	// int id_dphiL = controls.getId_faceVar("left delta-unknown");
	
	// right값
	int id_phiR = controls.getId_faceVar("left unknown");
	// int id_dphiR = controls.getId_faceVar("left delta-unknown");
	
	// 메쉬관련
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");


	{
		
		calcConvFlux.push_back(
		[id_phi,
		id_phiL,id_phiR,
		id_nx,id_ny,id_nz,id_area,id_dLR](
		double* cellsL, double* cellsR, double* faces, 
		double* fluxA_LL, double* fluxA_RR, 
		double* fluxA_LR, double* fluxA_RL, double* fluxB) ->int {
			
			double nvec[3];
			nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
			
			double phiL = cellsL[id_phi]; double phiR = cellsR[id_phi];
			
			double area = faces[id_area];
			double dLR = faces[id_dLR]; 
			
			fluxA_LL[0] -= 1.0/dLR*area;
			fluxA_RR[0] -= 1.0/dLR*area;
			fluxA_LR[0] += 1.0/dLR*area;
			fluxA_RL[0] += 1.0/dLR*area;
			
			fluxB[0] -= (phiR-phiL)/dLR*area;
			
		}); 
	}
	
	
	
	
	
	
	
	
	
	
	
	using conv_diff_Funct_BC_type = 
		function<int(
		double* cellsL, double* faces, 
		double* fluxA_LL, double* fluxB)>;
	
	calcConvFlux_BC.resize(1);

	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcConvFlux_BC[0].push_back(vector<conv_diff_Funct_BC_type>());
		
		string name = "unknown";
		string type = controls.boundaryMap[name][bcName+".type"];
		
		if(type=="fixedValue"){
			double value = stod(controls.boundaryMap[name][bcName+".value"]);
			calcConvFlux_BC[0].back().push_back(
			[id_phi,
			id_phiL,id_phiR,
			id_nx,id_ny,id_nz,id_area,id_dLR](
			double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB) ->int {

				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				
				double phiL = cellsL[id_phi]; double phiR = faces[id_phiR];
				
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				
				fluxA_LL[0] -= 1.0/dLR*area;
				
				fluxB[0] -= (phiR-phiL)/dLR*area;
				
				return 0;
			});
		}
		else if(type=="zeroGradient"){
			calcConvFlux_BC[0].back().push_back(
			[id_phi,
			id_phiL,id_phiR,
			id_nx,id_ny,id_nz,id_area,id_dLR](
			double* cellsL, double* faces, 
			double* fluxA_LL, double* fluxB) ->int {

				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				
				double phiL = cellsL[id_phi]; double phiR = faces[id_phiR];
				
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				
				// fluxA_LL[0] -= 1.0/dLR*area;
				
				fluxB[0] -= (phiR-phiL)/dLR*area;
				
				return 0;
			});
		}
	}
	
	
	
	
	
	
}