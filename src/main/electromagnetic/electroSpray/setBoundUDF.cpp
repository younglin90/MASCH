
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using Bound_Funct_type = function<int(
		double time, double x, double y, double z, 
		double* cells, double* faces)>;
		
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");

	int id_p = controls.getId_cellVar("pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_rho = controls.getId_cellVar("density");
	int id_mu = controls.getId_cellVar("viscosity");
	int id_rhoe = controls.getId_cellVar("charge-density");
	int id_xE = controls.getId_cellVar("x-electric-field");
	int id_yE = controls.getId_cellVar("y-electric-field");
	int id_zE = controls.getId_cellVar("z-electric-field");
	int id_phi = controls.getId_cellVar("electric-potential");
	int id_k = controls.getId_cellVar("conductivity");
	int id_epsilon = controls.getId_cellVar("permittivity");
    vector<int> id_alpha;
	for(int i=0; i<controls.spName.size()-1; ++i){
        id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+controls.spName[i]));
	}

	int id_pF = controls.getId_faceVar("pressure");
	int id_uF = controls.getId_faceVar("x-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	int id_rhoF = controls.getId_faceVar("density");
	int id_muF = controls.getId_faceVar("viscosity");
	int id_rhoeF = controls.getId_faceVar("charge-density");
	int id_xEF = controls.getId_faceVar("x-electric-field");
	int id_yEF = controls.getId_faceVar("y-electric-field");
	int id_zEF = controls.getId_faceVar("z-electric-field");
	int id_phiF = controls.getId_faceVar("electric-potential");
	int id_kF = controls.getId_faceVar("conductivity");
	int id_epsilonF = controls.getId_faceVar("permittivity");
    vector<int> id_alphaF, id_alphaL, id_alphaR;
	for(int i=0; i<controls.spName.size()-1; ++i){
        id_alphaF.push_back(controls.getId_faceVar("volume-fraction-"+controls.spName[i]));
        // id_alphaL.push_back(controls.getId_faceVar("left volume-fraction-"+controls.spName[i]));
        // id_alphaR.push_back(controls.getId_faceVar("right volume-fraction-"+controls.spName[i]));
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
			// else if(type=="inletOutlet"){
				// double value = stod(controls.boundaryMap["pressure"][bcName+".inletValue"]);
				// calcBoundFacePrimVal.back().push_back(
				// [value,id_p,id_pL,id_pR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
				// double time, double x, double y, double z, 
				// double* cells, double* faces) ->int {
					// double tmp_val = cells[id_p];
					// double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
					// faces[id_pL] = tmp_val;
					// if(Un<0.0) faces[id_pL] = value;
					// faces[id_pR] = faces[id_pL];
					// return 0;
				// });
			// }
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
				[id_uF, id_vF, id_wF](
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
			// else if(type=="massFlowRatePerAreaIdealGas"){
				// double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
				// double R = stod(controls.boundaryMap["velocity"][bcName+".R"]);
				
				// calcBoundFacePrimVal.back().push_back(
				// [value, R, id_u, id_v, id_w,id_uF, id_vF, id_wF,
				// id_nx, id_ny, id_nz, id_T, id_p](
				// double time, double x, double y, double z, 
				// double* cells, double* faces) ->int {
					// double mag_U = value/(cells[id_p]/R/cells[id_T]);
					// faces[id_uF] = -mag_U*faces[id_nx];
					// faces[id_vF] = -mag_U*faces[id_ny];
					// faces[id_wF] = -mag_U*faces[id_nz];
					// return 0;
				// });
			// }
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
					id_uF, id_vF, id_wF,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_u,id_v,id_w,id_p,id_T,id_uF,id_vF,id_wF,
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
	
			for(int i=0; i<controls.spName.size()-1; ++i){
				string spName = controls.spName[i];
				string type_name = (bcName+"."+spName+".type");
				string value_name = (bcName+"."+spName+".value");
				string type = controls.boundaryMap["volume-fraction"][type_name];
                
                int tmp_id_alpha = id_alpha[i];
                int tmp_id_alphaF = id_alphaF[i];
                // int tmp_id_alphaL = id_alphaL[i];
                // int tmp_id_alphaR = id_alphaR[i];
				
				if(type=="fixedValue"){
					double value = stod(controls.boundaryMap["volume-fraction"][value_name]);
					calcBoundFacePrimVal.back().push_back(
					// [value,tmp_id_alpha,tmp_id_alphaL,tmp_id_alphaR](
					[value,tmp_id_alpha,tmp_id_alphaF](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						faces[tmp_id_alphaF] = value;
						// faces[tmp_id_alphaR] = faces[tmp_id_alpha];
						return 0;
					});
				}
				else if(type=="zeroGradient"){
					calcBoundFacePrimVal.back().push_back(
					// [tmp_id_alpha,tmp_id_alphaL,tmp_id_alphaR](
					[tmp_id_alpha,tmp_id_alphaF](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						double tmp_val = cells[tmp_id_alpha];
						faces[tmp_id_alphaF] = tmp_val;
						// faces[tmp_id_alphaR] = faces[tmp_id_alpha];
						return 0;
					});
				}
				// else if(type=="inletOutlet"){
					// double value = stod(controls.boundaryMap["volume-fraction"][bcName+"."+spName+".inletValue"]);
					// calcBoundFacePrimVal.back().push_back(
					// [value,id_Y,id_YL,id_YR,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					// double time, double x, double y, double z, 
					// double* cells, double* faces) ->int {
						// double tmp_val = cells[id_Y];
						// double Un = cells[id_u]*faces[id_nx] + cells[id_v]*faces[id_ny] + cells[id_w]*faces[id_nz];
						// faces[id_YL] = tmp_val;
						// if(Un<0.0) faces[id_YL] = value;
						// faces[id_YR] = faces[id_YL];
						// return 0;
					// });
				// }
				else if(type=="function"){
					string inp_file = controls.boundaryMap["volume-fraction"][bcName+"."+spName+".file"];
					string inp_funct_name = controls.boundaryMap["volume-fraction"][bcName+"."+spName+".name"];
					void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
					char *error = nullptr;
					if (handle) {
						using setFunc_t = int(*)(double, double, double, double, 
						int, int, double*, double*);
						setFunc_t setFunction;
						*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
						calcBoundFacePrimVal.back().push_back(
						// [setFunction,tmp_id_alpha,tmp_id_alphaL,tmp_id_alphaR,
                        // id_u,id_v,id_w,id_nx,id_ny,id_nz](
						[setFunction,tmp_id_alpha,tmp_id_alphaF,
                        id_u,id_v,id_w,id_nx,id_ny,id_nz](
						double time, double x, double y, double z, 
						double* cells, double* faces) ->int {
							(*setFunction)(time,x,y,z,
							tmp_id_alpha,tmp_id_alphaF,cells,faces);
							// faces[tmp_id_alphaR] = faces[tmp_id_alphaL];
							return 0;
						});
					}
					else{
						cout << "#WARNING : file not there, " << inp_file << endl;
					}
				}
				else{
					cout << "#WARNING : not defiend " << "volume-fraction" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
			}
		}
        
        
        
        
        
		{
			string type = controls.boundaryMap["charge-density"][bcName+".type"];
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["charge-density"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_rhoe,id_rhoeF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_rhoeF] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_rhoe,id_rhoeF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_rhoe];
					faces[id_rhoeF] = tmp_val;
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap["charge-density"][bcName+".file"];
				string inp_funct_name = controls.boundaryMap["charge-density"][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction,id_rhoe,id_rhoeF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_rhoe,id_rhoeF,cells,faces);
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << "charge-density" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
        
        
        
		{
			string type = controls.boundaryMap["electric-potential"][bcName+".type"];
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap["electric-potential"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_phi,id_phiF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_phiF] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				calcBoundFacePrimVal.back().push_back(
				[id_phi,id_phiF](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_phi];
					faces[id_phiF] = tmp_val;
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.boundaryMap["electric-potential"][bcName+".file"];
				string inp_funct_name = controls.boundaryMap["electric-potential"][bcName+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, 
					int, int, double*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcBoundFacePrimVal.back().push_back(
					[setFunction,id_phi,id_phiF,id_u,id_v,id_w,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_phi,id_phiF,cells,faces);
						return 0;
					});
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			else{
				cout << "#WARNING : not defiend " << "electric-potential" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
		}
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



