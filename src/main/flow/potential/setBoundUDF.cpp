
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	// using Bound_Funct_type = function<int(double* cells, double* faces)>;
	using funct_type = 
	function<int(
	double time, double x, double y, double z, 
	double* cells, double* faces)>;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcBoundFacePrimVal.push_back(vector<funct_type>());
		
		// string name = "delta-pressure";
		for(auto& name : controls.primScalarNames)
		{
			string type = controls.boundaryMap[name][bcName+".type"];
			int id_C = controls.getId_cellVar(name);
			int id_FL = controls.getId_faceVar("left "+name);
			int id_FR = controls.getId_faceVar("right "+name);
		
			
			// vec_inpId.push_back(vector<int>());
			// vec_inpId.back().push_back(id);
			
			funct_type setFunct;
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap[name][bcName+".value"]);
				setFunct = (
				[value, id_C, id_FL, id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_FL] = value;
					faces[id_FR] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				setFunct = (
				[id_C, id_FL, id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_C];
					faces[id_FL] = tmp_val;
					faces[id_FR] = tmp_val;
					return 0;
				});
			}
			else if(type=="switch"){
				double value = stod(controls.boundaryMap[name][bcName+".value"]);
				int id_u = controls.getId_cellVar("x-velocity");
				int id_v = controls.getId_cellVar("y-velocity");
				int id_w = controls.getId_cellVar("z-velocity");
				int id_c = controls.getId_cellVar("speed-of-sound");
				setFunct = (
				[value,id_C, id_FL, id_FR, id_u,id_v,id_w,id_c](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double M = sqrt(
							cells[id_u]*cells[id_u]+cells[id_v]*cells[id_v]+cells[id_w]*cells[id_w])/
							cells[id_c];
					if(M>1.0){
						double tmp_val = cells[id_C];
						faces[id_FL] = tmp_val;
						faces[id_FR] = tmp_val;
					}
					else{
						faces[id_FL] = value;
						faces[id_FR] = value;
					}
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap[name][bcName+".file"];
				string inp_funct_name = controls.boundaryMap[name][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, int,
					double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					// setFunct = (*setFunction);
					setFunct = (
					[setFunction, id_C, id_FL, id_FR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_C,id_FL,id_FR,
						cells,faces);
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << name << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
			
			calcBoundFacePrimVal.back().push_back(setFunct);
		}
		
		
		
		for(auto& name : controls.primVector3Names){
			string type = controls.boundaryMap[name][bcName+".type"];
			vector<string> sub_names = controls.cellVar[name].sub_name;
			vector<int> sub_id_C;
			vector<int> sub_id_FL, sub_id_FR;
			for(auto& item : sub_names){
				sub_id_C.push_back(controls.getId_cellVar(item));
				sub_id_FL.push_back(controls.getId_faceVar("left "+item));
				sub_id_FR.push_back(controls.getId_faceVar("right "+item));
			}
			if(sub_id_C.size()!=3) cout << "#WARNING, sub_id_C.size()!=3" << endl;
			if(sub_id_FL.size()!=3) cout << "#WARNING, sub_id_FL.size()!=3" << endl;
			if(sub_id_FR.size()!=3) cout << "#WARNING, sub_id_FR.size()!=3" << endl;
			
			funct_type setFunct;
			
			// vec_inpId.push_back(sub_id);
			
			if(type=="fixedValue"){
				vector<string> s_value = load.extractVector(controls.boundaryMap[name][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				setFunct = (
				[value, sub_id_C, sub_id_FL, sub_id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
				// cout << value[0] << endl;
					faces[sub_id_FL[0]] = value[0];
					faces[sub_id_FL[1]] = value[1];
					faces[sub_id_FL[2]] = value[2];
					faces[sub_id_FR[0]] = value[0];
					faces[sub_id_FR[1]] = value[1];
					faces[sub_id_FR[2]] = value[2];
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				vector<string> s_value = load.extractVector(controls.boundaryMap[name][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				setFunct = (
				[value, sub_id_C, sub_id_FL, sub_id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val0 = cells[sub_id_C[0]];
					double tmp_val1 = cells[sub_id_C[1]];
					double tmp_val2 = cells[sub_id_C[2]];
					faces[sub_id_FL[0]] = tmp_val0;
					faces[sub_id_FL[1]] = tmp_val1;
					faces[sub_id_FL[2]] = tmp_val2;
					faces[sub_id_FR[0]] = tmp_val0;
					faces[sub_id_FR[1]] = tmp_val1;
					faces[sub_id_FR[2]] = tmp_val2;
					return 0;
				});
			}
			else if(type=="slip"){
				vector<string> s_value = load.extractVector(controls.boundaryMap[name][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				
				int id_nx = controls.getId_faceVar("x unit normal");
				int id_ny = controls.getId_faceVar("y unit normal");
				int id_nz = controls.getId_faceVar("z unit normal");
					
				setFunct = (
				[value, sub_id_C, sub_id_FL, sub_id_FR, id_nx, id_ny, id_nz](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double u_vel = cells[sub_id_C[0]];
					double v_vel = cells[sub_id_C[1]];
					double w_vel = cells[sub_id_C[2]];
					double norVel = u_vel * faces[id_nx] + 
									v_vel * faces[id_ny] + 
									w_vel * faces[id_nz];
					double invU = u_vel - norVel * faces[id_nx];
					double invV = v_vel - norVel * faces[id_ny];
					double invW = w_vel - norVel * faces[id_nz];
					faces[sub_id_FL[0]] = invU;
					faces[sub_id_FL[1]] = invV;
					faces[sub_id_FL[2]] = invW;
					faces[sub_id_FR[0]] = invU;
					faces[sub_id_FR[1]] = invV;
					faces[sub_id_FR[2]] = invW;
					return 0;
				});
			}
			else if(type=="noSlip"){
				vector<string> s_value = load.extractVector(controls.boundaryMap[name][bcName+".value"]);
				vector<double> value;
				for(auto& item : s_value) value.push_back(stod(item));
				setFunct = (
				[value, sub_id_C, sub_id_FL, sub_id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[sub_id_FL[0]] = 0.0;
					faces[sub_id_FL[1]] = 0.0;
					faces[sub_id_FL[2]] = 0.0;
					faces[sub_id_FR[0]] = 0.0;
					faces[sub_id_FR[1]] = 0.0;
					faces[sub_id_FR[2]] = 0.0;
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.initialMap[name]["file"];
				string inp_funct_name = controls.initialMap[name]["name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, int, int, int, int, int, int, int, 
					double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					setFunct = (
					[setFunction, sub_id_C, sub_id_FL, sub_id_FR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						sub_id_C[0],sub_id_C[1],sub_id_C[2],
						sub_id_FL[0],sub_id_FL[1],sub_id_FL[2],
						sub_id_FR[0],sub_id_FR[1],sub_id_FR[2],
						cells,faces);
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << name << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
			calcBoundFacePrimVal.back().push_back(setFunct);
		}
		
		
		
		for(auto& sup_name : controls.primVectorNames){
			vector<string> sub_names = controls.cellVar[sup_name].sub_name;
			vector<string> sub_roles = controls.cellVar[sup_name].sub_role;
			int iter=0;
			for(auto& name : sub_names){
				if(sub_roles[iter++]!="primitive") continue;
				string type = controls.boundaryMap[sup_name][bcName+"."+name+".type"]; 
				int id_C = controls.getId_cellVar(sup_name+"-"+name);
				int id_FL = controls.getId_faceVar("left "+sup_name+"-"+name);
				int id_FR = controls.getId_faceVar("right "+sup_name+"-"+name);
				
				// vec_inpId.push_back(vector<int>());
				// vec_inpId.back().push_back(id);
			

				funct_type setFunct;
				
				if(type=="fixedValue"){
					double value = stod(controls.boundaryMap[sup_name][bcName+"."+name+".value"]);
					setFunct = (
					[value, id_C, id_FL, id_FR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						faces[id_FL] = value;
						faces[id_FR] = value;
						return 0;
					});
				}
				else if(type=="zeroGradient"){
					setFunct = (
					[id_C, id_FL, id_FR](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[id_C];
						faces[id_FL] = tmp_val;
						faces[id_FR] = tmp_val;
						return 0;
					});
				}
				else if(type=="function"){
					string inp_file = controls.boundaryMap[sup_name][bcName+"."+name+".file"];
					string inp_funct_name = controls.boundaryMap[sup_name][bcName+"."+name+".name"];
					void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
					char *error = nullptr;
					if (handle) {
						using setFunc_t = int(*)(double, double, double, double, 
						int, int, int,
						double*, double*);
						setFunc_t setFunction;
						*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
						// setFunct = (*setFunction);
						setFunct = (
						[setFunction, id_C, id_FL, id_FR](
						double time, double x, double y, double z, 
						double* cells, double* faces) ->int {
							(*setFunction)(time,x,y,z,
							id_C,id_FL,id_FR,
							cells,faces);
							return 0;
						});
					}
					else{
						cout << "#WARNING : file not there, " << inp_file << endl;
					}
				}
				else{
					cout << "#WARNING : not defiend " << name << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
				calcBoundFacePrimVal.back().push_back(setFunct);
			}
		}
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



