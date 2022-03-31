
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
		

	int id_p = controls.getId_cellVar("pressure");
	int id_pL = controls.getId_faceVar("left pressure");
	int id_pR = controls.getId_faceVar("right pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uL = controls.getId_faceVar("left x-velocity");
	int id_uR = controls.getId_faceVar("right x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity");
	int id_vR = controls.getId_faceVar("right y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity");
	int id_wR = controls.getId_faceVar("right z-velocity");

	int id_T = controls.getId_cellVar("temperature");
	int id_TL = controls.getId_faceVar("left temperature");
	int id_TR = controls.getId_faceVar("right temperature");
				
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
		
	
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
					faces[id_pR] = faces[id_pL];
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
					faces[id_pR] = faces[id_pL];
					return 0;
				});
			}
			else if(type=="inletOutlet"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".inletValue"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pL,id_pR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
					faces[id_pL] = tmp_val;
					if(Un<0.0) faces[id_pL] = value;
					faces[id_pR] = faces[id_pL];
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap["pressure"][bcName+".file"];
				string inp_funct_name = controls.boundaryMap["pressure"][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction,id_p,id_pL,id_pR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_p,id_pL,cells,faces);
						faces[id_pR] = faces[id_pL];
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
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
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = cells[id_u]; faces[id_vL] = cells[id_v]; faces[id_wL] = cells[id_w];
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="slip"){
				vector<string> s_value = load.extractVector(controls.boundaryMap["velocity"][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
					
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
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="noSlip"){
				calcBoundFacePrimVal.back().push_back(
				[id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = 0.0; faces[id_vL] = 0.0; faces[id_wL] = 0.0;
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="surfaceNormalFixedValue"){
				double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
				
				calcBoundFacePrimVal.back().push_back(
				[value, id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,
				id_nx, id_ny, id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uL] = value*faces[id_nx]; 
					faces[id_vL] = value*faces[id_ny]; 
					faces[id_wL] = value*faces[id_nz]; 
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="massFlowRatePerAreaIdealGas"){
				double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
				double R = stod(controls.boundaryMap["velocity"][bcName+".R"]);
				
				calcBoundFacePrimVal.back().push_back(
				[value, R, id_u, id_v, id_w,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,
				id_nx, id_ny, id_nz, id_T, id_p](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double mag_U = value/(cells[id_p]/R/cells[id_T]);
					faces[id_uL] = -mag_U*faces[id_nx];
					faces[id_vL] = -mag_U*faces[id_ny];
					faces[id_wL] = -mag_U*faces[id_nz];
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap["velocity"][bcName+".file"];
				string inp_funct_name = controls.boundaryMap["velocity"][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, int, int, int, int, int, int, int, int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction, id_u,id_v,id_w, id_p, id_T, 
					id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_u,id_v,id_w,id_p,id_T,id_uL,id_vL,id_wL,
						id_nx,id_ny,id_nz,
						cells,faces);
						faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << "velocity" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		{
			string type = controls.boundaryMap["temperature"][bcName+".type"];
			
			// for(auto& [key, value] : controls.boundaryMap["temperature"]){
				// cout << key << " " << value << endl;
				// cout << bcName+".type" << "AAA" << endl;
				// cout << "suboutlet.type" << " " << controls.boundaryMap["temperature"]["suboutlet.type"] << endl;
				// cout << "invwall.type" << " " << controls.boundaryMap["temperature"]["invwall.type"] << endl;
				// cout << "viswall.type" << " " << controls.boundaryMap["temperature"]["viswall.type"] << endl;
				// cout << "liqinlet.type" << " " << controls.boundaryMap["temperature"]["liqinlet.type"] << endl;
				// cout << "gasinlet.type" << " " << controls.boundaryMap["temperature"]["gasinlet.type"] << endl;
			// }
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["temperature"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_T,id_TL,id_TR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_TL] = value;
					faces[id_TR] = faces[id_TL];
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
					faces[id_TR] = faces[id_TL];
					return 0;
				});
			}
			else if(type=="inletOutlet"){
				double value = stod(controls.boundaryMap["temperature"][bcName+".inletValue"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_T,id_TL,id_TR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_T];
					double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
					faces[id_TL] = tmp_val;
					if(Un<0.0) faces[id_TL] = value;
					faces[id_TR] = faces[id_TL];
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap["temperature"][bcName+".file"];
				string inp_funct_name = controls.boundaryMap["temperature"][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction,id_T,id_TL,id_TR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_T,id_TL,cells,faces);
						faces[id_TR] = faces[id_TL];
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << "temperature" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
		
		{
	
			for(int i=0; i<controls.spName.size()-1; ++i){
				string spName = controls.spName[i];
				string type_name = (bcName+"."+spName+".type");
				string value_name = (bcName+"."+spName+".value");
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
						faces[id_YR] = faces[id_YL];
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
						faces[id_YR] = faces[id_YL];
						return 0;
					});
				}
				else if(type=="inletOutlet"){
					double value = stod(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".inletValue"]);
					calcBoundFacePrimVal.back().push_back(
					[value,id_Y,id_YL,id_YR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[id_Y];
						double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
						faces[id_YL] = tmp_val;
						if(Un<0.0) faces[id_YL] = value;
						faces[id_YR] = faces[id_YL];
						return 0;
					});
				}
				else if(type=="function"){
					string inp_file = controls.boundaryMap["mass-fraction"][bcName+"."+spName+".file"];
					string inp_funct_name = controls.boundaryMap["mass-fraction"][bcName+"."+spName+".name"];
					void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
					char *error = nullptr;
					if (handle) {
						using setFunc_t = int(*)(double, double, double, double, 
						int, int, double*, double*);
						setFunc_t setFunction;
						*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
						calcBoundFacePrimVal.back().push_back(
						[setFunction,id_Y,id_YL,id_YR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
						double time, double x, double y, double z, 
						double* cells, double* faces) ->int {
							(*setFunction)(time,x,y,z,
							id_Y,id_YL,cells,faces);
							faces[id_YR] = faces[id_YL];
							return 0;
						});
					}
					else{
						cout << "#WARNING : file not there, " << inp_file << endl;
					}
				}
				else{
					cout << "#WARNING : not defiend " << "mass-fraction" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
			}
		}
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



