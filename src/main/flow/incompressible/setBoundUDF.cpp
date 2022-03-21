
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcBoundFacePrimVal.push_back(vector<Bound_Funct_type>());
		
		{
			string type = controls.boundaryMap["pressure"][bcName+".type"];
			int id_p = controls.getId_cellVar("pressure");
			int id_pF = controls.getId_faceVar("pressure");
			int id_dp = controls.getId_cellVar("delta-pressure");
			int id_dpF = controls.getId_faceVar("delta-pressure");
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pF,id_dp,id_dpF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_pF] = value;
					faces[id_dpF] = 0.0;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_p,id_pF,id_dp,id_dpF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					faces[id_pF] = tmp_val;
					faces[id_dpF] = cells[id_dp];
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << "pressure" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
		{
			string type = controls.boundaryMap["velocity"][bcName+".type"];
			
			int id_u = controls.getId_cellVar("x-velocity");
			int id_uF = controls.getId_faceVar("x-velocity");
			
			int id_v = controls.getId_cellVar("y-velocity");
			int id_vF = controls.getId_faceVar("y-velocity");
			
			int id_w = controls.getId_cellVar("z-velocity");
			int id_wF = controls.getId_faceVar("z-velocity");
			
			if(type=="fixedValue"){
				vector<string> s_value = load.extractVector(controls.boundaryMap["velocity"][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				calcBoundFacePrimVal.back().push_back(
				[value,id_u, id_v, id_w,id_uF, id_vF, id_wF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uF] = value[0]; faces[id_vF] = value[1]; faces[id_wF] = value[2];
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_u, id_v, id_w,id_uF, id_vF, id_wF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uF] = cells[id_u]; faces[id_vF] = cells[id_v]; faces[id_wF] = cells[id_w];
					return 0;
				});
			}
			else if(type=="slip"){
				vector<string> s_value = load.extractVector(controls.boundaryMap["velocity"][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				
				int id_nx = controls.getId_faceVar("x unit normal");
				int id_ny = controls.getId_faceVar("y unit normal");
				int id_nz = controls.getId_faceVar("z unit normal");
					
				calcBoundFacePrimVal.back().push_back(
				[value, id_u, id_v, id_w,id_uF, id_vF, id_wF,
				id_nx, id_ny, id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double u_vel = cells[id_u];
					double v_vel = cells[id_v];
					double w_vel = cells[id_w];
					double norVel = u_vel * faces[id_nx] + 
									v_vel * faces[id_ny] + 
									w_vel * faces[id_nz];
					double invU = u_vel - norVel * faces[id_nx];
					double invV = v_vel - norVel * faces[id_ny];
					double invW = w_vel - norVel * faces[id_nz];
					faces[id_uF] = invU; faces[id_vF] = invV; faces[id_wF] = invW;
					return 0;
				});
			}
			else if(type=="noSlip"){
				calcBoundFacePrimVal.back().push_back(
				[id_uF,id_vF,id_wF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uF] = 0.0; faces[id_vF] = 0.0; faces[id_wF] = 0.0;
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << "velocity" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



