
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
		

	int id_p = controls.getId_cellVar("pressure");
	int id_pF = controls.getId_faceVar("pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uF = controls.getId_faceVar("x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");

	int id_T = controls.getId_cellVar("temperature");
	int id_TF = controls.getId_faceVar("temperature");
				
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
				[value,id_p,id_pF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_pF] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_p,id_pF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					faces[id_pF] = tmp_val;
					return 0;
				});
			}
			else if(type=="inletOutlet"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".inletValue"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
					faces[id_pF] = tmp_val;
					if(Un<0.0) faces[id_pF] = value;
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
					[setFunction,id_p,id_pF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_p,id_pF,cells,faces);
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
			else if(type=="surfaceNormalFixedValue"){
				double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
				
				calcBoundFacePrimVal.back().push_back(
				[value, id_u, id_v, id_w,id_uF, id_vF, id_wF,
				id_nx, id_ny, id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_uF] = value*faces[id_nx]; 
					faces[id_vF] = value*faces[id_ny]; 
					faces[id_wF] = value*faces[id_nz]; 
					return 0;
				});
			}
			else if(type=="massFlowRatePerAreaIdealGas"){
				double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
				double R = stod(controls.boundaryMap["velocity"][bcName+".R"]);
				
				calcBoundFacePrimVal.back().push_back(
				[value, R, id_u, id_v, id_w,id_uF, id_vF, id_wF,
				id_nx, id_ny, id_nz, id_T, id_p](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double mag_U = value/(cells[id_p]/R/cells[id_T]);
					faces[id_uF] = -mag_U*cells[id_nx];
					faces[id_vF] = -mag_U*cells[id_ny];
					faces[id_wF] = -mag_U*cells[id_nz];
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
					int, int, int, int, int, int, int, int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction, id_u,id_v,id_w,id_uF,id_vF,id_wF,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_u,id_v,id_w,id_uF,id_vF,id_wF,
						id_nx,id_ny,id_nz,
						cells,faces);
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
				[value,id_T,id_TF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_TF] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_T,id_TF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_T];
					faces[id_TF] = tmp_val;
					return 0;
				});
			}
			else if(type=="inletOutlet"){
				double value = stod(controls.boundaryMap["temperature"][bcName+".inletValue"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_T,id_TF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_T];
					double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
					faces[id_TF] = tmp_val;
					if(Un<0.0) faces[id_TF] = value;
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
					[setFunction,id_T,id_TF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_T,id_TF,cells,faces);
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
				int id_YF = controls.getId_faceVar("mass-fraction-"+spName);
				
				if(type=="fixedValue"){
					double value = stod(controls.boundaryMap["mass-fraction"][value_name]);
					calcBoundFacePrimVal.back().push_back(
					[value,id_Y,id_YF](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						faces[id_YF] = value;
						return 0;
					});
				}
				else if(type=="zeroGradient"){
					calcBoundFacePrimVal.back().push_back(
					[id_Y,id_YF](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[id_Y];
						faces[id_YF] = tmp_val;
						return 0;
					});
				}
				else if(type=="inletOutlet"){
					double value = stod(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".inletValue"]);
					calcBoundFacePrimVal.back().push_back(
					[value,id_Y,id_YF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[id_Y];
						double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
						faces[id_YF] = tmp_val;
						if(Un<0.0) faces[id_YF] = value;
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
						[setFunction,id_Y,id_YF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
						double time, double x, double y, double z, 
						double* cells, double* faces) ->int {
							(*setFunction)(time,x,y,z,
							id_Y,id_YF,cells,faces);
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



