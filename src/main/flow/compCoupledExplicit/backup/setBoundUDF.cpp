
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
		

	int id_p = controls.getId_cellVar("pressure");
	// int id_pF = controls.getId_faceVar("pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	// int id_uF = controls.getId_faceVar("x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	// int id_vF = controls.getId_faceVar("y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	// int id_wF = controls.getId_faceVar("z-velocity");

	int id_T = controls.getId_cellVar("temperature");
	// int id_TF = controls.getId_faceVar("temperature");
			
	int id_pL = controls.getId_faceVar("left pressure"); int id_pR = controls.getId_faceVar("right pressure");
	int id_uL = controls.getId_faceVar("left x-velocity"); int id_uR = controls.getId_faceVar("right x-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity"); int id_vR = controls.getId_faceVar("right y-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity"); int id_wR = controls.getId_faceVar("right z-velocity");
	int id_TL = controls.getId_faceVar("left temperature"); int id_TR = controls.getId_faceVar("right temperature");
	vector<int> id_YL, id_YR;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_YL.push_back(controls.getId_faceVar("left mass-fraction-"+controls.spName[i]));
		id_YR.push_back(controls.getId_faceVar("right mass-fraction-"+controls.spName[i]));
	}
		
	
	for(int ibc=0; ibc<mesh.boundaries.size(); ++ibc){
		auto& boundary = mesh.boundaries[ibc];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcBoundFacePrimVal.push_back(vector<Bound_Funct_type>());
		
		{
			string type = controls.boundaryMap["pressure"][bcName+".type"];
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pL,id_pR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_pL] = value;
					faces[id_pR] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_p,id_pL,id_pR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					faces[id_pL] = tmp_val;
					faces[id_pR] = tmp_val;
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << "pressure" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
		{
			string type = controls.boundaryMap["velocity"][bcName+".type"];
			
			if(type=="fixedValue"){
				vector<string> s_value = load.extractVector(controls.boundaryMap["velocity"][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				calcBoundFacePrimVal.back().push_back(
				[value,id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = value[0]; faces[id_vL] = value[1]; faces[id_wL] = value[2];
					faces[id_uR] = value[0]; faces[id_vR] = value[1]; faces[id_wR] = value[2];
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = cells[id_u]; faces[id_vL] = cells[id_v]; faces[id_wL] = cells[id_w];
					faces[id_uR] = cells[id_u]; faces[id_vR] = cells[id_v]; faces[id_wR] = cells[id_w];
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
				[value, id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,
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
					faces[id_uL] = invU; faces[id_vL] = invV; faces[id_wL] = invW;
					faces[id_uR] = invU; faces[id_vR] = invV; faces[id_wR] = invW;
					return 0;
				});
			}
			else if(type=="noSlip"){
				calcBoundFacePrimVal.back().push_back(
				[id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = 0.0; faces[id_vL] = 0.0; faces[id_wL] = 0.0;
					faces[id_uR] = 0.0; faces[id_vR] = 0.0; faces[id_wR] = 0.0;
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << "velocity" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		{
			string type = controls.boundaryMap["temperature"][bcName+".type"];
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["temperature"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_T,id_TL,id_TR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_TL] = value;
					faces[id_TR] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_T,id_TL,id_TR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_T];
					faces[id_TL] = tmp_val;
					faces[id_TR] = tmp_val;
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << "temperature" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
		
		{
	
			for(int i=0; i<controls.spName.size()-1; ++i){
				string spName = controls.spName[i];
				string type_name = (bcName+spName+".type");
				string value_name = (bcName+spName+".value");
				string type = controls.boundaryMap["mass-fraction"][type_name];
	
				int id_Y = controls.getId_cellVar("mass-fraction-"+spName);
				int id_YL = controls.getId_faceVar("left mass-fraction-"+spName);
				int id_YR = controls.getId_faceVar("right mass-fraction-"+spName);
				
				if(type=="fixedValue"){
					double value = stod(controls.boundaryMap["mass-fraction"][value_name]);
					calcBoundFacePrimVal.back().push_back(
					[value,id_Y,id_YL,id_YR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						faces[id_YL] = value;
						faces[id_YR] = value;
						return 0;
					});
				}
				else if(type=="zeroGradient"){
					calcBoundFacePrimVal.back().push_back(
					[id_Y,id_YL,id_YR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[id_Y];
						faces[id_YL] = tmp_val;
						faces[id_YR] = tmp_val;
						return 0;
					});
				}
				else{
					cout << "#WARNING : not defiend " << "temperature" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
			}
		}
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



