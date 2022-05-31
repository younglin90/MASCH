
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){


    // double eps = 1.e-50;
    double eps = -10.0;


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
		
	int id_c = controls.getId_cellVar("speed-of-sound");
	
	int id_rho = controls.getId_cellVar("density");
    
    
    
    
    int nSp = controls.spName.size();
    // controls.timeVaryingMappedFixedValueNCycle.resize(5+nSp-1,0.0);
    controls.timeVaryingMappedFixedValueTimeCycle.resize(5+nSp-1);
    // controls.timeVaryingMappedFixedValueTime1.resize(5+nSp-1);
    // controls.timeVaryingMappedFixedValueTime2.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTimeOrder1.resize(5+nSp-1);
    // controls.timeVaryingMappedFixedValueTimeOrder2.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueTime.resize(5+nSp-1);
    // controls.timeVaryingMappedFixedValueValue1.resize(5+nSp-1);
    // controls.timeVaryingMappedFixedValueValue2.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueFileName.resize(5+nSp-1);
    controls.timeVaryingMappedFixedValueValueIter.resize(5+nSp-1,0);
    
    
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
			else if(type=="fixedMean"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pL,id_pR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_p];
					faces[id_pL] = 0.5*(tmp_val+value);
					faces[id_pR] = faces[id_pL];
					return 0;
				});
			}
			else if(type=="switchMach"){
				double value = stod(controls.boundaryMap["pressure"][bcName+".value"]);
				calcBoundFacePrimVal.back().push_back(
				[value,id_p,id_pL,id_pR, id_c,id_u,id_v,id_w](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    double tmp_c = cells[id_c];
					double tmp_u = cells[id_u];
					double tmp_v = cells[id_v];
					double tmp_w = cells[id_w];
                    double tmp_U = sqrt(tmp_u*tmp_u+tmp_v*tmp_v+tmp_w*tmp_w);
                    double machNumber = tmp_U/(tmp_c+1.e-200);
					double tmp_val = cells[id_p];
                    if(machNumber<1.0){
                        tmp_val = value;
                    }
					faces[id_pL] = tmp_val;
					faces[id_pR] = faces[id_pL];
					return 0;
				});
			}
			else if(type=="totalPressure"){
				// double p0 = stod(controls.boundaryMap["pressure"][bcName+".p0"]);
				// calcBoundFacePrimVal.back().push_back(
				// [id_p,id_pL,id_pR, id_c,id_u,id_v,id_w,id_rho,p0](
				// double time, double x, double y, double z, 
				// double* cells, double* faces) ->int {
                    // double tmp_rho = cells[id_rho];
					// double tmp_u = cells[id_u];
					// double tmp_v = cells[id_v];
					// double tmp_w = cells[id_w];
                    // double tmp_U = sqrt(tmp_u*tmp_u+tmp_v*tmp_v+tmp_w*tmp_w);
					// faces[id_pL] = p0 - 0.5*tmp_rho*tmp_U*tmp_U;
					// faces[id_pR] = faces[id_pL];
					// return 0;
				// });
                
				double U0 = stod(controls.boundaryMap["pressure"][bcName+".U"]);
				calcBoundFacePrimVal.back().push_back(
				[id_p,id_pL,id_pR, id_c,id_u,id_v,id_w,id_rho,U0](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    double tmp_rho = cells[id_rho];
                    double tmp_p = cells[id_p];
					double tmp_u = cells[id_u];
					double tmp_v = cells[id_v];
					double tmp_w = cells[id_w];
                    double tmp_U = sqrt(tmp_u*tmp_u+tmp_v*tmp_v+tmp_w*tmp_w);
					faces[id_pL] = tmp_p + 0.5*tmp_rho*(tmp_U-U0)*(tmp_U-U0);
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
                
                int tmp_size = controls.timeVaryingMappedFixedValueFileName[id_time_phi].size();
                
                vector<vector<double>> inp_values(tmp_size);
                int str = boundary.startFace;
                int nFaces = boundary.nFaces;
                int end = str + nFaces;
                
                // 파일 열기
                for(int idt=0; idt<tmp_size; ++idt){
                        
                    vector<double> tmp_x;
                    vector<double> tmp_y;
                    vector<double> tmp_z;
                    vector<double> tmp_values;
                    tmp_x.reserve(end-str);
                    tmp_y.reserve(end-str);
                    tmp_z.reserve(end-str);
                    tmp_values.reserve(end-str);
                    
                    vector<vector<string>> csv_contents;
                    int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1;
                    string timeName = controls.timeVaryingMappedFixedValueFileName[id_time_phi][idt];
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

                        csv_contents.push_back(items);
                        counter += 1;
                    }

                    int tmp_i = 0;
                    for (auto& name : csv_contents[0]) {
                        if(name=="pressure") id_phi1 = tmp_i;
                        if(name=="x") id_x1 = tmp_i;
                        if(name=="y") id_y1 = tmp_i;
                        if(name=="z") id_z1 = tmp_i;
                        ++tmp_i;
                    }
                    if(id_phi1==-1) cout << "#WARNING, id_phi == -1" << endl;
                    if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                    if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                    if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;     
                    
                    for(int j=1; j<csv_contents.size(); ++j){
                        tmp_x.push_back(stod(csv_contents[j][id_x1]));
                        tmp_y.push_back(stod(csv_contents[j][id_y1]));
                        tmp_z.push_back(stod(csv_contents[j][id_z1]));
                        tmp_values.push_back(stod(csv_contents[j][id_phi1]));
                    }
                    
                    int value_size = tmp_x.size();
                    
                    for(int k=str; k<end; ++k){
                        auto& face = mesh.faces[k];
                        
                        // //=======================
                        // mesh.cells[face.iL].level = -1;
                        // //=======================
                        
                        double face_x = face.x;
                        double face_y = face.y;
                        double face_z = face.z;
                        double tmp_v = 0.0;
                        double max_dist = 1.e50;
                        for(int j=0; j<value_size; ++j){
                            double dx = face_x - tmp_x[j];
                            double dy = face_y - tmp_y[j];
                            double dz = face_z - tmp_z[j];
                            
                            double dist = dx*dx + dy*dy + dz*dz;
                            
                            if(dist<max_dist){
                                tmp_v = tmp_values[j];
                                max_dist = dist;
                            }
                            if(dist<eps) break;
                            
                        }
                        inp_values[idt].push_back(tmp_v);
                    }
                }
                
                
                
                
				vector<string> s_value = load.extractVector(controls.boundaryMap["pressure"][bcName+".time"]);
                for(auto& item : s_value){
                    controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                }
                controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                stod(controls.boundaryMap["pressure"][bcName+".timeCycle"]);
                
                // int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                // double namuji = input_time - (double)jung * 
                    // controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
                // controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = (double)jung;
                
                // for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                    // auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                    // if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                        // controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = ii;
                        // controls.timeVaryingMappedFixedValueTimeOrder2[id_time_phi] = (ii+1 >= tmp_size ? 0 : ii+1);
                        // // controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                        // // controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                        // break;
                    // }
                // }
                
                
                controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = 0;
                
                double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls, id_pL,id_pR, id_time_phi, inp_values, nFaces,
                timeCycle](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    int& iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi];
                    double tmp_val = -10.0;
                    {
                        int nCycle = time / timeCycle;
                        double new_time = time - timeCycle*((double)nCycle);
                        int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                        int order2 = order1 + 1;
                        double time1, time2;
                        int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                        int tmp_iter = 0;
                        while(1){
                            time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                            if(order2>=order_size) {
                                order2 = 0;
                                time2 = timeCycle;
                            }
                            else{
                                time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                            }
                            if(time1<=new_time && new_time<=time2){
                                
                                break;
                            }
                            else{
                                order1 += 1;
                                if(order1>=order_size) order1 = 0;
                                order2 = order1 + 1;
                            }
                            
                            if(tmp_iter==order_size) new_time -= timeCycle;
                            if(new_time<=0.0) new_time = 1.e-200;
                            if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                            
                            ++tmp_iter;
                        }
                        
                        double value1 = inp_values.at(order1).at(iter);
                        double value2 = inp_values.at(order2).at(iter);
                        
                        tmp_val = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                    }
                    
                    
					faces[id_pL] = tmp_val;
					faces[id_pR] = faces[id_pL];
                    
                    ++iter;
                    if(inp_values[0].size() <= iter) iter = 0;
                    
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
			else if(type=="inletOutlet"){
				double value = stod(controls.boundaryMap["velocity"][bcName+".value"]);
					
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
                    if(norVel>=0.0){
                        // cout << value << " " << norVel << " " << u_vel << " " << v_vel << " " << w_vel << endl;
                        faces[id_uL] = u_vel; 
                        faces[id_vL] = v_vel; 
                        faces[id_wL] = w_vel; 
                    }
                    else{
                        faces[id_uL] = value * faces[id_nx]; 
                        faces[id_vL] = value * faces[id_ny]; 
                        faces[id_wL] = value * faces[id_nz];
                    }
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
            // synthetic eddy method
			else if(type=="SEM"){
                
                solver.sem.boolExcute = true;
                
                
				int nEddies = stoi(controls.boundaryMap["velocity"][bcName+".nEddies"]);
                double sigma = stod(controls.boundaryMap["velocity"][bcName+".eddySize"]);
                double dR = stod(controls.boundaryMap["velocity"][bcName+".domainSize"]);
                
				vector<string> s_str = load.extractVector(controls.boundaryMap["velocity"][bcName+".startBoundaryFacePoints"]);
				vector<double> str;
				for(auto& item : s_str) str.push_back(stod(item));
                
				vector<string> s_end = load.extractVector(controls.boundaryMap["velocity"][bcName+".endBoundaryFacePoints"]);
				vector<double> end;
				for(auto& item : s_end) end.push_back(stod(item));
                
				vector<string> s_setU = load.extractVector(controls.boundaryMap["velocity"][bcName+".velocities"]);
				vector<double> setU;
				for(auto& item : s_setU) setU.push_back(stod(item));
                
				vector<string> s_Rij = load.extractVector(controls.boundaryMap["velocity"][bcName+".reynoldsTensor"]);
				vector<double> Rij;
				for(auto& item : s_Rij) Rij.push_back(stod(item));
                
                
                

                // nEddies, sigma, dR, str, end, setU, Rij
                solver.sem.initialSetting(
                    nEddies, sigma, dR, str, end, setU, Rij);
                      
                // int id_dt = controls.getId_fieldVar("time-step");
                // double& tmp_dt = var.fields[id_dt];

				calcBoundFacePrimVal.back().push_back(
				[&solver, id_uL, id_vL, id_wL,id_uR, id_vR, id_wR](
				double time, double x, double y, double z, 
				double* cells, double* faces) mutable ->int {
                    
                    vector<double> u_fluct = solver.sem.calcFluctuationVelocities(x, y, z);
					faces[id_uL] = solver.sem.setU[0] + u_fluct[0];
					faces[id_vL] = solver.sem.setU[1] + u_fluct[1];
					faces[id_wL] = solver.sem.setU[2] + u_fluct[2];
					faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];


					return 0;
				});

                
            }
            // synthetic eddy method
			else if(type=="function+SEM"){
                
                solver.sem.boolExcute = true;
                
				int nEddies = stoi(controls.boundaryMap["velocity"][bcName+".nEddies"]);
                double sigma = stod(controls.boundaryMap["velocity"][bcName+".eddySize"]);
                double dR = stod(controls.boundaryMap["velocity"][bcName+".domainSize"]);
                
				vector<string> s_str = load.extractVector(controls.boundaryMap["velocity"][bcName+".startBoundaryFacePoints"]);
				vector<double> str;
				for(auto& item : s_str) str.push_back(stod(item));
                
				vector<string> s_end = load.extractVector(controls.boundaryMap["velocity"][bcName+".endBoundaryFacePoints"]);
				vector<double> end;
				for(auto& item : s_end) end.push_back(stod(item));
                
				vector<string> s_setU = load.extractVector(controls.boundaryMap["velocity"][bcName+".velocities"]);
				vector<double> setU;
				for(auto& item : s_setU) setU.push_back(stod(item));
                
				vector<string> s_Rij = load.extractVector(controls.boundaryMap["velocity"][bcName+".reynoldsTensor"]);
				vector<double> Rij;
				for(auto& item : s_Rij) Rij.push_back(stod(item));
                
                
                

                // nEddies, sigma, dR, str, end, setU, Rij
                solver.sem.initialSetting(
                    nEddies, sigma, dR, str, end, setU, Rij);
                      
                      
                      
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
					[&solver, setFunction, id_u,id_v,id_w, id_p, id_T, 
					id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,id_nx,id_ny,id_nz](
					double time, double x, double y, double z, 
					double* cells, double* faces) ->int {
						(*setFunction)(time,x,y,z,
						id_u,id_v,id_w,id_p,id_T,id_uL,id_vL,id_wL,
						id_nx,id_ny,id_nz,
						cells,faces);
						faces[id_uR] = faces[id_uL]; faces[id_vR] = faces[id_vL]; faces[id_wR] = faces[id_wL];
                        
                        vector<double> u_fluct = solver.sem.calcFluctuationVelocities(x, y, z);
                        faces[id_uL] += u_fluct[0];
                        faces[id_vL] += u_fluct[1];
                        faces[id_wL] += u_fluct[2];
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
                    
                    int tmp_size = controls.timeVaryingMappedFixedValueFileName[id_time_phi].size();
                    
                
                    for(auto& item : s_value){
                        controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                    }
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                    stod(controls.boundaryMap["velocity"][bcName+".timeCycle"]);
                    
                    // int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    // double namuji = input_time - (double)jung * 
                        // controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    
                    // controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                    
                    // for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                        // auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                        // if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                            // controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = ii;
                            // controls.timeVaryingMappedFixedValueTimeOrder2[id_time_phi] = (ii+1 >= tmp_size ? 0 : ii+1);
                            // // controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                            // // controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                            // break;
                        // }
                    // }
                }
                
                
                    
                int tmp_size = controls.timeVaryingMappedFixedValueFileName[1].size();
                
                vector<vector<double>> inp_values1(tmp_size);
                vector<vector<double>> inp_values2(tmp_size);
                vector<vector<double>> inp_values3(tmp_size);
                int str = boundary.startFace;
                int end = str + boundary.nFaces;
                
                // 파일 열기
                for(int idt=0; idt<tmp_size; ++idt){
                        
                    vector<double> tmp_x;
                    vector<double> tmp_y;
                    vector<double> tmp_z;
                    vector<double> tmp_values1;
                    vector<double> tmp_values2;
                    vector<double> tmp_values3;
                    tmp_x.reserve(end-str);
                    tmp_y.reserve(end-str);
                    tmp_z.reserve(end-str);
                    tmp_values1.reserve(end-str);
                    tmp_values2.reserve(end-str);
                    tmp_values3.reserve(end-str);
                    
                    vector<vector<string>> csv_contents;
                    int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1, id_phi2=-1, id_phi3=-1;
                    string timeName = controls.timeVaryingMappedFixedValueFileName[1][idt];
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

                        csv_contents.push_back(items);
                        counter += 1;
                    }

                    int tmp_i = 0;
                    for (auto& name : csv_contents[0]) {
                        if(name=="x-velocity") id_phi1 = tmp_i;
                        if(name=="y-velocity") id_phi2 = tmp_i;
                        if(name=="z-velocity") id_phi3 = tmp_i;
                        if(name=="x") id_x1 = tmp_i;
                        if(name=="y") id_y1 = tmp_i;
                        if(name=="z") id_z1 = tmp_i;
                        ++tmp_i;
                    }
                    if(id_phi1==-1) cout << "#WARNING, id_phi1 == -1" << endl;
                    if(id_phi2==-1) cout << "#WARNING, id_phi2 == -1" << endl;
                    if(id_phi3==-1) cout << "#WARNING, id_phi3 == -1" << endl;
                    if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                    if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                    if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;     
                    
                    for(int j=1; j<csv_contents.size(); ++j){
                        tmp_x.push_back(stod(csv_contents[j][id_x1]));
                        tmp_y.push_back(stod(csv_contents[j][id_y1]));
                        tmp_z.push_back(stod(csv_contents[j][id_z1]));
                        tmp_values1.push_back(stod(csv_contents[j][id_phi1]));
                        tmp_values2.push_back(stod(csv_contents[j][id_phi2]));
                        tmp_values3.push_back(stod(csv_contents[j][id_phi3]));
                    }
                    
                    int value_size = tmp_x.size();
                    
                    for(int k=str; k<end; ++k){
                        auto& face = mesh.faces[k];
                        
                        // //=======================
                        // mesh.cells[face.iL].level = -1;
                        // //=======================
                        
                        double face_x = face.x;
                        double face_y = face.y;
                        double face_z = face.z;
                        double tmp_v1 = 0.0;
                        double tmp_v2 = 0.0;
                        double tmp_v3 = 0.0;
                        double max_dist = 1.e50;
                        for(int j=0; j<value_size; ++j){
                            double dx = face_x - tmp_x[j];
                            double dy = face_y - tmp_y[j];
                            double dz = face_z - tmp_z[j];
                            
                            double dist = dx*dx + dy*dy + dz*dz;
                            
                            if(dist<max_dist){
                                tmp_v1 = tmp_values1[j];
                                tmp_v2 = tmp_values2[j];
                                tmp_v3 = tmp_values3[j];
                                max_dist = dist;
                            }
                            if(dist<eps) break;
                            
                        }
                        inp_values1[idt].push_back(tmp_v1);
                        inp_values2[idt].push_back(tmp_v2);
                        inp_values3[idt].push_back(tmp_v3);
                    }
                }
                
                controls.timeVaryingMappedFixedValueTimeOrder1[1] = 0;
                controls.timeVaryingMappedFixedValueTimeOrder1[2] = 0;
                controls.timeVaryingMappedFixedValueTimeOrder1[3] = 0;
                
                double timeCycle1 = controls.timeVaryingMappedFixedValueTimeCycle[1];
                double timeCycle2 = controls.timeVaryingMappedFixedValueTimeCycle[2];
                double timeCycle3 = controls.timeVaryingMappedFixedValueTimeCycle[3];
                
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls,id_uL, id_vL, id_wL,id_uR, id_vR, id_wR,
                inp_values1,inp_values2,inp_values3,timeCycle1,timeCycle2,timeCycle3](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    
                    int& iter = controls.timeVaryingMappedFixedValueValueIter[1];
                    
                    double tmp_val1 = 0.0;
                    {
                        int id_time_phi = 1;
                        double timeCycle = timeCycle1;
                        
                        int nCycle = time / timeCycle;
                        double new_time = time - timeCycle*((double)nCycle);
                        int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                        int order2 = order1 + 1;
                        double time1, time2;
                        int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                        int tmp_iter = 0;
                        while(1){
                            time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                            if(order2>=order_size) {
                                order2 = 0;
                                time2 = timeCycle;
                            }
                            else{
                                time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                            }
                            if(time1<=new_time && new_time<=time2){
                                
                                break;
                            }
                            else{
                                order1 += 1;
                                if(order1>=order_size) order1 = 0;
                                order2 = order1 + 1;
                            }
                            
                            if(tmp_iter==order_size) new_time -= timeCycle;
                            if(new_time<=0.0) new_time = 1.e-200;
                            if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                            
                            ++tmp_iter;
                        }
                        
                        double value1 = inp_values1.at(order1).at(iter);
                        double value2 = inp_values1.at(order2).at(iter);
                        
                        tmp_val1 = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                        
                        // int order1 = 0;
                        // double value1 = inp_values1.at(order1).at(iter);
                        // tmp_val1 = value1;
                    }
                    
                    double tmp_val2 = 0.0;
                    {
                        int id_time_phi = 2;
                        double timeCycle = timeCycle2;
                        
                        int nCycle = time / timeCycle;
                        double new_time = time - timeCycle*((double)nCycle);
                        int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                        int order2 = order1 + 1;
                        double time1, time2;
                        int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                        int tmp_iter = 0;
                        while(1){
                            time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                            if(order2>=order_size) {
                                order2 = 0;
                                time2 = timeCycle;
                            }
                            else{
                                time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                            }
                            if(time1<=new_time && new_time<=time2){
                                
                                break;
                            }
                            else{
                                order1 += 1;
                                if(order1>=order_size) order1 = 0;
                                order2 = order1 + 1;
                            }
                            
                            if(tmp_iter==order_size) new_time -= timeCycle;
                            if(new_time<=0.0) new_time = 1.e-200;
                            if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                            
                            ++tmp_iter;
                        }
                        
                        double value1 = inp_values2.at(order1).at(iter);
                        double value2 = inp_values2.at(order2).at(iter);
                        
                        tmp_val2 = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                        
                        // int order1 = 0;
                        // double value1 = inp_values2.at(order1).at(iter);
                        // tmp_val2 = value1;
                    }
                    
                    double tmp_val3 = 0.0;
                    {
                        int id_time_phi = 3;
                        double timeCycle = timeCycle3;
                        
                        int nCycle = time / timeCycle;
                        double new_time = time - timeCycle*((double)nCycle);
                        int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                        int order2 = order1 + 1;
                        double time1, time2;
                        int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                        int tmp_iter = 0;
                        while(1){
                            time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                            if(order2>=order_size) {
                                order2 = 0;
                                time2 = timeCycle;
                            }
                            else{
                                time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                            }
                            if(time1<=new_time && new_time<=time2){
                                
                                break;
                            }
                            else{
                                order1 += 1;
                                if(order1>=order_size) order1 = 0;
                                order2 = order1 + 1;
                            }
                            
                            if(tmp_iter==order_size) new_time -= timeCycle;
                            if(new_time<=0.0) new_time = 1.e-200;
                            if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                            
                            ++tmp_iter;
                        }
                        
                        double value1 = inp_values3.at(order1).at(iter);
                        double value2 = inp_values3.at(order2).at(iter);
                        
                        tmp_val3 = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                        
                        // int order1 = 0;
                        // double value1 = inp_values3.at(order1).at(iter);
                        // tmp_val3 = value1;
                    }
                    
                    
					faces[id_uL] = tmp_val1; faces[id_vL] = tmp_val2; faces[id_wL] = tmp_val3;
					faces[id_uR] = tmp_val1; faces[id_vR] = tmp_val2; faces[id_wR] = tmp_val3;
                    
                    ++iter;
                    if(inp_values1[0].size() <= iter) iter = 0;
                    
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
                
                
                int tmp_size = controls.timeVaryingMappedFixedValueFileName[id_time_phi].size();
                
                vector<vector<double>> inp_values(tmp_size);
                int str = boundary.startFace;
                int end = str + boundary.nFaces;
                
                // 파일 열기
                for(int idt=0; idt<tmp_size; ++idt){
                        
                    vector<double> tmp_x;
                    vector<double> tmp_y;
                    vector<double> tmp_z;
                    vector<double> tmp_values;
                    tmp_x.reserve(end-str);
                    tmp_y.reserve(end-str);
                    tmp_z.reserve(end-str);
                    tmp_values.reserve(end-str);
                    
                    vector<vector<string>> csv_contents;
                    int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1;
                    string timeName = controls.timeVaryingMappedFixedValueFileName[id_time_phi][idt];
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

                        csv_contents.push_back(items);
                        counter += 1;
                    }

                    int tmp_i = 0;
                    for (auto& name : csv_contents[0]) {
                        if(name=="temperature") id_phi1 = tmp_i;
                        if(name=="x") id_x1 = tmp_i;
                        if(name=="y") id_y1 = tmp_i;
                        if(name=="z") id_z1 = tmp_i;
                        ++tmp_i;
                    }
                    if(id_phi1==-1) cout << "#WARNING, id_phi == -1" << endl;
                    if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                    if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                    if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;     
                    
                    for(int j=1; j<csv_contents.size(); ++j){
                        tmp_x.push_back(stod(csv_contents[j][id_x1]));
                        tmp_y.push_back(stod(csv_contents[j][id_y1]));
                        tmp_z.push_back(stod(csv_contents[j][id_z1]));
                        tmp_values.push_back(stod(csv_contents[j][id_phi1]));
                    }
                    
                    int value_size = tmp_x.size();
                    
                    for(int k=str; k<end; ++k){
                        auto& face = mesh.faces[k];
                        
                        // //=======================
                        // mesh.cells[face.iL].level = -1;
                        // //=======================
                        
                        double face_x = face.x;
                        double face_y = face.y;
                        double face_z = face.z;
                        double tmp_v = 0.0;
                        double max_dist = 1.e50;
                        for(int j=0; j<value_size; ++j){
                            double dx = face_x - tmp_x[j];
                            double dy = face_y - tmp_y[j];
                            double dz = face_z - tmp_z[j];
                            
                            double dist = dx*dx + dy*dy + dz*dz;
                            
                            if(dist<max_dist){
                                tmp_v = tmp_values[j];
                                max_dist = dist;
                            }
                            if(dist<eps) break;
                            
                        }
                        inp_values[idt].push_back(tmp_v);
                    }
                }
                
                
                
                vector<string> s_value = load.extractVector(controls.boundaryMap["temperature"][bcName+".time"]);
                
                for(auto& item : s_value){
                    controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                }
                controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                stod(controls.boundaryMap["temperature"][bcName+".timeCycle"]);
                
                // int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                // double namuji = input_time - (double)jung * 
                    // controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
                // controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                
                // for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                    // auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                    // if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                        // // controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                        // controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = ii;
                        // controls.timeVaryingMappedFixedValueTimeOrder2[id_time_phi] = (ii+1 >= tmp_size ? 0 : ii+1);
                        // controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                        // controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                        // break;
                    // }
                // }
                
                controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = 0;
                
                double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
            
				calcBoundFacePrimVal.back().push_back(
				[&solver, &controls, id_TL, id_TR, id_time_phi, inp_values,
                timeCycle](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
                    int& iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi];
                    double tmp_val = -10.0;
                    {
                        int nCycle = time / timeCycle;
                        double new_time = time - timeCycle*((double)nCycle);
                        int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                        int order2 = order1 + 1;
                        double time1, time2;
                        int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                        int tmp_iter = 0;
                        while(1){
                            time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                            if(order2>=order_size) {
                                order2 = 0;
                                time2 = timeCycle;
                            }
                            else{
                                time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                            }
                            if(time1<=new_time && new_time<=time2){
                                
                                break;
                            }
                            else{
                                order1 += 1;
                                if(order1>=order_size) order1 = 0;
                                order2 = order1 + 1;
                            }
                            
                            if(tmp_iter==order_size) new_time -= timeCycle;
                            if(new_time<=0.0) new_time = 1.e-200;
                            if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                            
                            ++tmp_iter;
                        }
                        
                        double value1 = inp_values.at(order1).at(iter);
                        double value2 = inp_values.at(order2).at(iter);
                        
                        tmp_val = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                        
                        // int order1 = 0;
                        // double value1 = inp_values.at(order1).at(iter);
                        // tmp_val = value1;
                    }
                    
					faces[id_TL] = tmp_val;
					faces[id_TR] = faces[id_TL];
                    
                    ++iter;
                    if(inp_values[0].size() <= iter) iter = 0;
                    
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
                    
                
                    int tmp_size = controls.timeVaryingMappedFixedValueFileName[id_time_phi].size();
                    
                    vector<vector<double>> inp_values(tmp_size);
                    int str = boundary.startFace;
                    int end = str + boundary.nFaces;
                    
                    // 파일 열기
                    for(int idt=0; idt<tmp_size; ++idt){
                            
                        vector<double> tmp_x;
                        vector<double> tmp_y;
                        vector<double> tmp_z;
                        vector<double> tmp_values;
                        tmp_x.reserve(end-str);
                        tmp_y.reserve(end-str);
                        tmp_z.reserve(end-str);
                        tmp_values.reserve(end-str);
                        
                        vector<vector<string>> csv_contents;
                        int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1;
                        string timeName = controls.timeVaryingMappedFixedValueFileName[id_time_phi][idt];
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

                            csv_contents.push_back(items);
                            counter += 1;
                        }

                        int tmp_i = 0;
                        for (auto& name : csv_contents[0]) {
                            if(name=="mass-fraction-"+spName) id_phi1 = tmp_i;
                            if(name=="x") id_x1 = tmp_i;
                            if(name=="y") id_y1 = tmp_i;
                            if(name=="z") id_z1 = tmp_i;
                            ++tmp_i;
                        }
                        if(id_phi1==-1) cout << "#WARNING, id_phi == -1" << endl;
                        if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                        if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                        if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;     
                        
                        for(int j=1; j<csv_contents.size(); ++j){
                            tmp_x.push_back(stod(csv_contents[j][id_x1]));
                            tmp_y.push_back(stod(csv_contents[j][id_y1]));
                            tmp_z.push_back(stod(csv_contents[j][id_z1]));
                            tmp_values.push_back(stod(csv_contents[j][id_phi1]));
                        }
                        
                        int value_size = tmp_x.size();
                        
                        for(int k=str; k<end; ++k){
                            auto& face = mesh.faces[k];
                        
                            // //=======================
                            // if(mesh.cells[face.iL].level==-1) mesh.cells[face.iL].level=0;
                            // mesh.cells[face.iL].level = -1;
                            // //=======================
                        
                            double face_x = face.x;
                            double face_y = face.y;
                            double face_z = face.z;
                            double tmp_v = 0.0;
                            double max_dist = 1.e50;
                            for(int j=0; j<value_size; ++j){
                                double dx = face_x - tmp_x[j];
                                double dy = face_y - tmp_y[j];
                                double dz = face_z - tmp_z[j];
                                
                                double dist = dx*dx + dy*dy + dz*dz;
                                
                                if(dist<max_dist){
                                    tmp_v = tmp_values[j];
                                    max_dist = dist;
                                }
                                if(dist<eps) break;
                                
                            }
                            inp_values[idt].push_back(tmp_v);
                        }
                    }
                    
                    
                        
                    vector<string> s_value = load.extractVector(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".time"]);
                    
                    for(auto& item : s_value){
                        controls.timeVaryingMappedFixedValueTime[id_time_phi].push_back(stod(item));
                    }
                    controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi] = 
                    stod(controls.boundaryMap["mass-fraction"][bcName+"."+spName+".timeCycle"]);
                    
                    // int jung = input_time / controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    // double namuji = input_time - (double)jung * 
                        // controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                    
                    // controls.timeVaryingMappedFixedValueNCycle[id_time_phi] = jung;
                    
                    // for(int ii=0; ii<controls.timeVaryingMappedFixedValueTime[id_time_phi].size()-1; ++ii){
                        // auto& tmp = controls.timeVaryingMappedFixedValueTime[id_time_phi];
                        // if(tmp[ii]<=namuji && namuji<tmp[ii+1]){
                            // // controls.timeVaryingMappedFixedValueTimeOrder[id_time_phi] = ii;
                            // controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = ii;
                            // controls.timeVaryingMappedFixedValueTimeOrder2[id_time_phi] = (ii+1 >= tmp_size ? 0 : ii+1);
                            // controls.timeVaryingMappedFixedValueTime1[id_time_phi] = tmp[ii];
                            // controls.timeVaryingMappedFixedValueTime2[id_time_phi] = tmp[ii+1];
                            // break;
                        // }
                    // }
                    
                    controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi] = 0;
                    
                    double timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_time_phi];
                
                    calcBoundFacePrimVal.back().push_back(
                    [&solver, &controls,id_YL,id_YR, id_time_phi, inp_values,
                    timeCycle](
                    double time, double x, double y, double z, 
                    double* cells, double* faces) ->int {
                        int& iter = controls.timeVaryingMappedFixedValueValueIter[id_time_phi];
                        double tmp_val = -10.0;
                        {
                            int nCycle = time / timeCycle;
                            double new_time = time - timeCycle*((double)nCycle);
                            int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_time_phi];
                            int order2 = order1 + 1;
                            double time1, time2;
                            int order_size = controls.timeVaryingMappedFixedValueTime[id_time_phi].size();
                            int tmp_iter = 0;
                            while(1){
                                time1 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order1);
                                if(order2>=order_size) {
                                    order2 = 0;
                                    time2 = timeCycle;
                                }
                                else{
                                    time2 = controls.timeVaryingMappedFixedValueTime[id_time_phi].at(order2);
                                }
                                if(time1<=new_time && new_time<=time2){
                                    
                                    break;
                                }
                                else{
                                    order1 += 1;
                                    if(order1>=order_size) order1 = 0;
                                    order2 = order1 + 1;
                                }
                                
                                if(tmp_iter==order_size) new_time -= timeCycle;
                                if(new_time<=0.0) new_time = 1.e-200;
                                if(tmp_iter>1000) cout << "#WARNING, timeVaryingMapped boundary error" << endl;
                                
                                ++tmp_iter;
                            }
                            
                            double value1 = inp_values.at(order1).at(iter);
                            double value2 = inp_values.at(order2).at(iter);
                            
                            tmp_val = value1 + (value2-value1)/(time2 - time1)*(new_time - time1);
                            
                            // int order1 = 0;
                            // double value1 = inp_values[order1].at(iter);
                            // tmp_val = value1;
                        }
                        
                        faces[id_YL] = tmp_val;
                        faces[id_YR] = faces[id_YL];
                    
                        ++iter;
                        if(inp_values[0].size() <= iter) iter = 0;
                    
                        return 0;
                    });
                    
                    
                    
                    
                    
                }
				else{
					cout << "#WARNING : not defiend " << "mass-fraction" << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
				}
			}
		}
		
		
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}













// update time varying boundary values
void MASCH_Solver::updateTimeVaryingMappedFixedValue(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
    
	// auto& solver = (*this);
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size();
	
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
    // auto fieldVar = var.fields.data();
    // int id_t = controls.getId_fieldVar("time");
    // double time = fieldVar[id_t];
    
    
	// for(int ibc=0; ibc<mesh.boundaries.size(); ++ibc){
		// auto& boundary = mesh.boundaries[ibc];
		// if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
        // int str = boundary.startFace;
        // int end = str + boundary.nFaces;
		
		// string bcName = boundary.name;
        // {
            // string type = controls.boundaryMap["pressure"][bcName+".type"];
            // if(type == "timeVaryingMappedFixedValue"){
                
                // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "pressure", 0);
                
                // int id_phi = 0;
                
                // double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
                // double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
                // // double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
                // // double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
                // int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_phi];
                // int& order2 = controls.timeVaryingMappedFixedValueTimeOrder2[id_phi];
                // int tmp_size = controls.timeVaryingMappedFixedValueTime[id_phi].size();
                
                // double time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order2];
                
                // double dTime = time-timeCycle*nCycle;
                // if(dTime>=timeCycle*nCycle){
                    // nCycle += 1.0;
                    // order1 = 0;
                    // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                // }
                // else{
                    // if(dTime>=time2){
                        // ++order1;
                        // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                        
                    // }
                // }
                
            // }
        // }
        // {
            // string type = controls.boundaryMap["velocity"][bcName+".type"];
            // if(type == "timeVaryingMappedFixedValue"){
                
                // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "x-velocity", 1);
                // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "y-velocity", 2);
                // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "z-velocity", 3);
                
                // for(int id_phi=1; id_phi<4; ++id_phi){
                
                    // double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
                    // double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
                    // // double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
                    // // double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
                    // int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_phi];
                    // int& order2 = controls.timeVaryingMappedFixedValueTimeOrder2[id_phi];
                    // int tmp_size = controls.timeVaryingMappedFixedValueTime[id_phi].size();
                    
                    // double time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order2];
                    
                    // double dTime = time-timeCycle*nCycle;
                    // if(dTime>=timeCycle*nCycle){
                        // nCycle += 1.0;
                        // order1 = 0;
                        // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                    // }
                    // else{
                        // if(dTime>=time2){
                            // ++order1;
                            // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                            
                        // }
                    // }
                // }
                
            // }
        // }
        // {
            // string type = controls.boundaryMap["temperature"][bcName+".type"];
            // if(type == "timeVaryingMappedFixedValue"){
                
                // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, "temperature", 4);
                
                // int id_phi = 4;
                
                // double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
                // double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
                // // double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
                // // double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
                // int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_phi];
                // int& order2 = controls.timeVaryingMappedFixedValueTimeOrder2[id_phi];
                // int tmp_size = controls.timeVaryingMappedFixedValueTime[id_phi].size();
                
                // double time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order2];
                
                // double dTime = time-timeCycle*nCycle;
                // if(dTime>=timeCycle*nCycle){
                    // nCycle += 1.0;
                    // order1 = 0;
                    // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                // }
                // else{
                    // if(dTime>=time2){
                        // ++order1;
                        // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                        
                    // }
                // }
                
            // }
        // }
        // {
			// for(int i=0; i<controls.spName.size()-1; ++i){
				// string spName = controls.spName[i];
                // string type = controls.boundaryMap["mass-fraction-"+spName][bcName+".type"];
                // if(type == "timeVaryingMappedFixedValue"){
                    
                    // // solver.updateTimeVaryingMappedFixedValue12(mesh, controls, var, str, end, 
                        // // "mass-fraction-"+spName, 5+i);
                        
                    // int id_phi = 5+i;
                    
                    // double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
                    // double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
                    // // double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
                    // // double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
                    // int& order1 = controls.timeVaryingMappedFixedValueTimeOrder1[id_phi];
                    // int& order2 = controls.timeVaryingMappedFixedValueTimeOrder2[id_phi];
                    // int tmp_size = controls.timeVaryingMappedFixedValueTime[id_phi].size();
                
                    // double time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order2];
                    
                    // double dTime = time-timeCycle*nCycle;
                    // if(dTime>=timeCycle*nCycle){
                        // nCycle += 1.0;
                        // order1 = 0;
                        // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                    // }
                    // else{
                        // if(dTime>=time2){
                            // ++order1;
                            // order2 = (order1+1 >= tmp_size ? 0 : order1+1);
                            
                        // }
                    // }
                    
                // }
            // }
        // }
        
    // }
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    
	
}



void MASCH_Solver::updateTimeVaryingMappedFixedValue12(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
int str, int end, string phi_name, int id_phi){
    
	// auto fieldVar = var.fields.data();
	// int id_t = controls.getId_fieldVar("time");
    
    // double time = fieldVar[id_t];

    // double& nCycle = controls.timeVaryingMappedFixedValueNCycle[id_phi];
    // double& timeCycle = controls.timeVaryingMappedFixedValueTimeCycle[id_phi];
    // double& time1 = controls.timeVaryingMappedFixedValueTime1[id_phi];
    // double& time2 = controls.timeVaryingMappedFixedValueTime2[id_phi];
    // int& order = controls.timeVaryingMappedFixedValueTimeOrder[id_phi];
    
    // double dTime = time-timeCycle*nCycle;
    // bool boolSwitchTime = false;
    
    // // cout << dTime << " " << timeCycle << " " << nCycle << " " << time2 << endl;
    
    
    // if(dTime>=timeCycle*nCycle){
        // nCycle += 1.0;
        // order = 0;
        // boolSwitchTime = true;
    // }
    // else{
        // if(dTime>=time2){
            // ++order;
            // boolSwitchTime = true;
            
        // }
    // }
    // string fileName1, fileName2;
    // if(controls.timeVaryingMappedFixedValueTime[id_phi].size() == order+1){
        // time1 = controls.timeVaryingMappedFixedValueTime[id_phi][order];
        // time2 = timeCycle;
        
        // fileName1 = controls.timeVaryingMappedFixedValueFileName[id_phi][order];
        // fileName2 = controls.timeVaryingMappedFixedValueFileName[id_phi][0];
        
    // }
    // else{
        // time1 = controls.timeVaryingMappedFixedValueTime[id_phi][order];
        // time2 = controls.timeVaryingMappedFixedValueTime[id_phi][order+1];
        
        // fileName1 = controls.timeVaryingMappedFixedValueFileName[id_phi][order];
        // fileName2 = controls.timeVaryingMappedFixedValueFileName[id_phi][order+1];
    // }
    
    // if(boolSwitchTime==true){
    
        // // 파일 열기
        // vector<vector<string>> csv_contents1;
        // vector<vector<string>> csv_contents2;
        // int id_x1=-1, id_y1=-1, id_z1=-1, id_phi1=-1;
        // int id_x2=-1, id_y2=-1, id_z2=-1, id_phi2=-1;
        // for(int ii=0; ii<2; ++ii)
        // {
            // string timeName;
            // if(ii==0) timeName = fileName1;
            // if(ii==1) timeName = fileName2;
            // string filename("./setting/boundary/"+timeName+".csv");
            // string file_contents;
            // char delimiter = ',';

            // {
                // auto ss = ostringstream{};
                // ifstream input_file(filename);
                // if (!input_file.is_open()) {
                    // cerr << "Could not open the file - '"
                        // << filename << "'" << endl;
                    // exit(EXIT_FAILURE);
                // }
                // ss << input_file.rdbuf();
                // file_contents = ss.str();
            // }

            // istringstream sstream(file_contents);
            // string record;

            // int counter = 0;
            // while (std::getline(sstream, record)) {
                // vector<string> items;
                // istringstream line(record);
                // while (std::getline(line, record, delimiter)) {
                    // record.erase(record.begin(), std::find_if(record.begin(), record.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
                    // record.erase(std::find_if(record.rbegin(), record.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), record.end()); 
                    // // record.erase(std::remove_if(record.begin(), record.end(), std::isspace), record.end());
                    // items.push_back(record);
                // }

                // if(ii==0) csv_contents1.push_back(items);
                // if(ii==1) csv_contents2.push_back(items);
                // counter += 1;
            // }

            // int tmp_i = 0;
            // if(ii==0){
                // for (auto& name : csv_contents1[0]) {
                    // if(name==phi_name) id_phi1 = tmp_i;
                    // if(name=="x") id_x1 = tmp_i;
                    // if(name=="y") id_y1 = tmp_i;
                    // if(name=="z") id_z1 = tmp_i;
                    // ++tmp_i;
                // }
                // if(id_phi1==-1) cout << "#WARNING, id_phi == -1" << endl;
                // if(id_x1==-1) cout << "#WARNING, id_x == -1" << endl;
                // if(id_y1==-1) cout << "#WARNING, id_y == -1" << endl;
                // if(id_z1==-1) cout << "#WARNING, id_z == -1" << endl;
            // }
            // else{
                // for (auto& name : csv_contents2[0]) {
                    // if(name==phi_name) id_phi2 = tmp_i;
                    // if(name=="x") id_x2 = tmp_i;
                    // if(name=="y") id_y2 = tmp_i;
                    // if(name=="z") id_z2 = tmp_i;
                    // ++tmp_i;
                // }
                // if(id_phi2==-1) cout << "#WARNING, id_phi == -1" << endl;
                // if(id_x2==-1) cout << "#WARNING, id_x == -1" << endl;
                // if(id_y2==-1) cout << "#WARNING, id_y == -1" << endl;
                // if(id_z2==-1) cout << "#WARNING, id_z == -1" << endl;
            // }
        // }
        
        // auto& values1 = controls.timeVaryingMappedFixedValueValue1[id_phi];
        // values1.clear();
        // auto& values2 = controls.timeVaryingMappedFixedValueValue2[id_phi];
        // values2.clear();
        // double eps = 1.e-12;
        // // cout << csv_contents1.size() << " " << csv_contents2.size() << endl;
        // for(int i=str; i<end; ++i){
            // // cout << i << endl;
            // auto& face = mesh.faces[i];
            // {
                // double min_dist = 1.e15;
                // double inp_phi = -1.e15;
                // auto csv_contents1_ptr = csv_contents1.data();
                // for (int i = 1; i < csv_contents1.size(); ++i) {
                    // auto csv_contents1_ptr_i = csv_contents1_ptr[i].data();
                    // double tmp_x = stod(csv_contents1_ptr_i[id_x1]);
                    // double tmp_y = stod(csv_contents1_ptr_i[id_y1]);
                    // double tmp_z = stod(csv_contents1_ptr_i[id_z1]);
                
                    // double dx = face.x - tmp_x;
                    // double dy = face.y - tmp_y;
                    // double dz = face.z - tmp_z;
                    // double dist = (dx*dx+dy*dy+dz*dz);
                    
                    // if(dist<min_dist){
                        // inp_phi = stod(csv_contents1_ptr_i[id_phi1]);
                        // min_dist = dist;
                    // }
                    
                    // if(dist<eps) break;
                // }
                // values1.push_back(inp_phi);
            // }
            // {
                // double min_dist = 1.e15;
                // double inp_phi = -1.e15;
                // auto csv_contents2_ptr = csv_contents2.data();
                // for (int i = 1; i < csv_contents2.size(); ++i) {
                    // auto csv_contents2_ptr_i = csv_contents2_ptr[i].data();
                    // double tmp_x = stod(csv_contents2_ptr_i[id_x2]);
                    // double tmp_y = stod(csv_contents2_ptr_i[id_y2]);
                    // double tmp_z = stod(csv_contents2_ptr_i[id_z2]);
                
                    // double dx = face.x - tmp_x;
                    // double dy = face.y - tmp_y;
                    // double dz = face.z - tmp_z;
                    // double dist = (dx*dx+dy*dy+dz*dz);
                    
                    // if(dist<min_dist){
                        // inp_phi = stod(csv_contents2_ptr_i[id_phi2]);
                        // min_dist = dist;
                    // }
                    
                    // if(dist<eps) break;
                // }
                // values2.push_back(inp_phi);
            // }
        // }
        
        
        
    // }
    
}