
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
		
	
    
    
    
    
    int nSp = controls.spName.size();
    controls.timeVaryingMappedFixedValueNCycle.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTimeCycle.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTime1.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTime2.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTimeOrder.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTime.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueValue1.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueValue2.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueFileName.resize(5+nSp-1);
    
    
    double input_time = stod(controls.startFrom);
    
    
    
    
    
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
			else if(type=="timeVaryingMappedFixedValue"){
                int id_time_phi = 0;
                
				controls.timeVaryingMappedFixedValueFileName[id_time_phi] = 
                load.extractVector(controls.boundaryMap["pressure"][bcName+".file"]);
                
				vector<string> s_value = load.extractVector(controls.boundaryMap["pressure"][bcName+".time"]);
                for(auto& item : s_value){
                    controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                }
                controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                stod(controls.boundaryMap["pressure"][bcName+".timeCycle"]);
                
                int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                double namuji = input_time - (double)jung * 
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
                controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                
                for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                    auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                    if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                        controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                        controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                        controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                        break;
                    }
                }
                
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls, id_pL,id_pR, id_time_phi](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                    if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                        controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                    }
                    double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                    double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                    double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                    double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                    double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                    double dTime = time-timeCycle*nCycle;
                    double tmp_val = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                    
					faces[id_pL] = tmp_val;
					faces[id_pR] = faces[id_pL];
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
			else if(type=="timeVaryingMappedFixedValue"){
                
                vector<string> s_value = load.extractVector(controls.boundaryMap["velocity"][bcName+".time"]);
                
                for(int id_time_phi=1; id_time_phi<4; ++id_time_phi){
                
                    controls.timeVaryingMappedFixedValueFileName[id_time_phi] = 
                    load.extractVector(controls.boundaryMap["velocity"][bcName+".file"]);
                
                    for(auto& item : s_value){
                        controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                    }
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                    stod(controls.boundaryMap["velocity"][bcName+".timeCycle"]);
                    
                    int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    double namuji = input_time - (double)jung * 
                        controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    
                    controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                    
                    for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                        auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                        if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                            controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                            controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                            controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                            break;
                        }
                    }
                }
                
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    double tmp_val1 = 0.0;
                    {
                        int id_time_phi = 1;
                        int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                        if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                            controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                        }
                        double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                        double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                        double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                        double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                        double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                        double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                        double dTime = time-timeCycle*nCycle;
                        tmp_val1 = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                    }
                    double tmp_val2 = 0.0;
                    {
                        int id_time_phi = 2;
                        int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                        if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                            controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                        }
                        double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                        double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                        double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                        double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                        double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                        double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                        double dTime = time-timeCycle*nCycle;
                        tmp_val2 = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                    }
                    double tmp_val3 = 0.0;
                    {
                        int id_time_phi = 3;
                        int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                        if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                            controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                        }
                        double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                        double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                        double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                        double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                        double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                        double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                        double dTime = time-timeCycle*nCycle;
                        tmp_val3 = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                    }
                    
					faces[id_uL] = tmp_val1; faces[id_vL] = tmp_val2; faces[id_wL] = tmp_val3;
					faces[id_uR] = tmp_val1; faces[id_vR] = tmp_val2; faces[id_wR] = tmp_val3;
					return 0;
				});
                
                
                
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
			else if(type=="timeVaryingMappedFixedValue"){
                int id_time_phi = 4;
                
                controls.timeVaryingMappedFixedValueFileName[id_time_phi] = 
                load.extractVector(controls.boundaryMap["temperature"][bcName+".file"]);
                
                vector<string> s_value = load.extractVector(controls.boundaryMap["temperature"][bcName+".time"]);
                
                for(auto& item : s_value){
                    controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                }
                controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                stod(controls.boundaryMap["temperature"][bcName+".timeCycle"]);
                
                int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                double namuji = input_time - (double)jung * 
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
                controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                
                for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                    auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                    if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                        controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                        controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                        controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                        break;
                    }
                }
            
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls, id_TL, id_TR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    double tmp_val1 = 0.0;
                    {
                        int id_time_phi = 4;
                        int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                        if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                            controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                        }
                        double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                        double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                        double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                        double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                        double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                        double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                        double dTime = time-timeCycle*nCycle;
                        tmp_val1 = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                    }
                    
					faces[id_TL] = tmp_val1;
					faces[id_TR] = faces[id_TL];
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
                else if(type=="timeVaryingMappedFixedValue"){
                    
                    int id_time_phi = 5+i;
                    
                    controls.timeVaryingMappedFixedValueFileName[id_time_phi] = 
                    load.extractVector(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".file"]);
                        
                    vector<string> s_value = load.extractVector(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".time"]);
                    
                    for(auto& item : s_value){
                        controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                    }
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                    stod(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".timeCycle"]);
                    
                    int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    double namuji = input_time - (double)jung * 
                        controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    
                    controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                    
                    for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                        auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                        if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                            controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                            controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                            controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                            break;
                        }
                    }
                
                    calcBoundFacePrimVal.back().push_back(
                    [&solver, &controls,id_YL,id_YR, i](
                    double time, double x, double y, double z, 
                    double* cells, double* faces) ->int {
                        double tmp_val1 = 0.0;
                        {
                            int id_time_phi = 5+i;
                            int iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi]++;
                            if(iter+1 == controls.timeVaryingMappedFixedValueValue1[id_time_phi].size()) {
                                controls.timeVaryingMappedFixedValueValueIter[id_time_phi] = 0;
                            }
                            double value1 = controls.timeVaryingMappedFixedValueValue1[id_time_phi][iter];
                            double value2 = controls.timeVaryingMappedFixedValueValue2[id_time_phi][iter];
                            double time1 = controls.timeVaryingMappedFixedValueTime1[id_time_phi];
                            double time2 = controls.timeVaryingMappedFixedValueTime2[id_time_phi];
                            double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                            double nCycle = controls.timeVaryingMappedFixedValueNCycle[id_time_phi];
                            double dTime = time-timeCycle*nCycle;
                            tmp_val1 = value1 + (value2-value1)/(time2 - time1)*(dTime - time1);
                        }
                        
                        faces[id_YL] = tmp_val1;
                        faces[id_YR] = faces[id_YL];
                        return 0;
                    });
                    
                    
                    
                    
                    
                }
				else{
					cout << "#WARNING : not defiend " << "mass-fraction" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
			}
		}
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}













// update time varying boundary values
void MASCH_Solver::updateTimeVaryingMappedFixedValue(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
    
	auto& solver = (*this);
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
    
    
	for(int ibc=0; ibc<mesh.boundaries.size(); ++ibc){
		auto& boundary = mesh.boundaries[ibc];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
        int str = boundary.startFace;
        int end = str + boundary.nFaces;
		
		string bcName = boundary.name;
        {
            string type = controls.boundaryMap["pressure"][bcName+".type"];
            if(type == "timeVaryingMappedFixedValue"){
                
                solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "pressure", 0);
                
            }
        }
        {
            string type = controls.boundaryMap["velocity"][bcName+".type"];
            if(type == "timeVaryingMappedFixedValue"){
                
                solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "x-velocity", 1);
                solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "y-velocity", 2);
                solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "z-velocity", 3);
                
            }
        }
        {
            string type = controls.boundaryMap["temperature"][bcName+".type"];
            if(type == "timeVaryingMappedFixedValue"){
                
                solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "temperature", 4);
                
            }
        }
        {
			for(int i=0; i<controls.spName.size()-1; ++i){
				string spName = controls.spName[i];
                string type = controls.boundaryMap["mass-fraction-"+spName][bcName+".type"];
                if(type == "timeVaryingMappedFixedValue"){
                    
                    solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, 
                        "mass-fraction-"+spName, 5+i);
                    
                }
            }
        }
        
    }
    
	
}



void MASCH_Solver::updateTimeVaryingMappedFixedValue12(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
int str, int end, string phi_name, int id_phi){
    
	auto fieldVar = var.fields.data();
	int id_t = controls.getId_fieldVar("time");
    
    double time = fieldVar[id_t];

    double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
    double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
    double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
    double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
    int& order = controls.timeVaryingMappedFixedValueTimeOrder[id_phi];
    
    double dTime = time-timeCycle*nCycle;
    bool boolSwitchTime = false;
    if(dTime>=timeCycle){
        nCycle += 1.0;
        order = 0;
        boolSwitchTime = true;
    }
    else{
        if(dTime>=time2){
            ++order;
            boolSwitchTime = true;
            
        }
    }
    string fileName1, fileName2;
    if(controls.timeVaryingMappedFixedValueTime[id_phi].size() == order+1){
        time1 = controls.timeVaryingMappedFixedValueTime[id_phi][order];
        time2 = timeCycle;
        
        fileName1 = controls.timeVaryingMappedFixedValueFileName[id_phi][order];
        fileName2 = controls.timeVaryingMappedFixedValueFileName[id_phi][0];
        
    }
    else{
        time1 = controls.timeVaryingMappedFixedValueTime[id_phi][order];
        time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order+1];
        
        fileName1 = controls.timeVaryingMappedFixedValueFileName[id_phi][order];
        fileName2 = controls.timeVaryingMappedFixedValueFileName[id_phi][order+1];
    }
    
    if(boolSwitchTime==true){
    
        // 파일 열기
        vector<vector<string>> csv_contents1;
        vector<vector<string>> csv_contents2;
        int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1;
        int id_x2=-1, id_y2=-1, id_z2=-1, id_phi2=-1;
        for(int ii=0; ii<2; ++ii)
        {
            string timeName;
            if(ii==0) timeName = fileName1;
            if(ii==1) timeName = fileName2;
            string filename("./setting/boundary/"+timeName+".csv");
            string file_contents;
            char delimiter = ',';

            {
                auto ss = ostringstream{};
                ifstream input_file(filename);
                if (!input_file.is_open()) {
                    cerr << "Could not open the file - '"
                        << filename << "'" << endl;
                    exit(EXIT_FAILURE);
                }
                ss << input_file.rdbuf();
                file_contents = ss.str();
            }

            istringstream sstream(file_contents);
            string record;

            int counter = 0;
            while (std::getline(sstream, record)) {
                vector<string> items;
                istringstream line(record);
                while (std::getline(line, record, delimiter)) {
                    record.erase(record.begin(), std::find_if(record.begin(), record.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
                    record.erase(std::find_if(record.rbegin(), record.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), record.end()); 
                    // record.erase(std::remove_if(record.begin(), record.end(), std::isspace), record.end());
                    items.push_back(record);
                }

                if(ii==0) csv_contents1.push_back(items);
                if(ii==1) csv_contents2.push_back(items);
                counter += 1;
            }

            int tmp_i = 0;
            if(ii==0){
                for (auto& name : csv_contents1[0]) {
                    if(name==phi_name) id_phi1 = tmp_i;
                    if(name=="x") id_x1 = tmp_i;
                    if(name=="y") id_y1 = tmp_i;
                    if(name=="z") id_z1 = tmp_i;
                    ++tmp_i;
                }
                if(id_phi1==-1) cout << "#WARNING, id_phi == -1" << endl;
                if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;
            }
            else{
                for (auto& name : csv_contents2[0]) {
                    if(name==phi_name) id_phi2 = tmp_i;
                    if(name=="x") id_x2 = tmp_i;
                    if(name=="y") id_y2 = tmp_i;
                    if(name=="z") id_z2 = tmp_i;
                    ++tmp_i;
                }
                if(id_phi2==-1) cout << "#WARNING, id_phi == -1" << endl;
                if(id_x2==-1) cout << "#WARNING, id_x == -1" << endl;
                if(id_y2==-1) cout << "#WARNING, id_y == -1" << endl;
                if(id_z2==-1) cout << "#WARNING, id_z == -1" << endl;
            }
        }
        
        auto& values1 = controls.timeVaryingMappedFixedValueValue1[id_phi];
        values1.clear();
        auto& values2 = controls.timeVaryingMappedFixedValueValue2[id_phi];
        values2.clear();
        for(int i=str; i<end; ++i){
            auto& face = mesh.faces[i];
            {
                double min_dist = 1.e15;
                double inp_phi = -1.e15;
                for (int i = 1; i < csv_contents1.size(); ++i) {
                    double tmp_x = stod(csv_contents1[i][id_x1]);
                    double tmp_y = stod(csv_contents1[i][id_y1]);
                    double tmp_z = stod(csv_contents1[i][id_z1]);
                    double tmp_phi = stod(csv_contents1[i][id_phi1]);
                
                    double dx = face.x - tmp_x;
                    double dy = face.y - tmp_y;
                    double dz = face.z - tmp_z;
                    double dist = sqrt(dx*dx+dy*dy+dz*dz);
                    
                    if(dist<min_dist){
                        inp_phi = tmp_phi;
                        min_dist = dist;
                    }
                }
                values1.push_back(inp_phi);
            }
            {
                double min_dist = 1.e15;
                double inp_phi = -1.e15;
                for (int i = 1; i < csv_contents2.size(); ++i) {
                    double tmp_x = stod(csv_contents2[i][id_x2]);
                    double tmp_y = stod(csv_contents2[i][id_y2]);
                    double tmp_z = stod(csv_contents2[i][id_z2]);
                    double tmp_phi = stod(csv_contents2[i][id_phi2]);
                
                    double dx = face.x - tmp_x;
                    double dy = face.y - tmp_y;
                    double dz = face.z - tmp_z;
                    double dist = sqrt(dx*dx+dy*dy+dz*dz);
                    
                    if(dist<min_dist){
                        inp_phi = tmp_phi;
                        min_dist = dist;
                    }
                }
                values2.push_back(inp_phi);
            }
        }
        
        
    }
    
}