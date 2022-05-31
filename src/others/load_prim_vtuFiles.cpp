
#include "./load.h" 


void MASCH_Load::fvmFiles(string folderName, int rank,
	MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
		
	// MASCH_Mesh_Load mesh_load;
	(*this).vtuPrimitive(folderName, rank, controls, mesh, var);
}

void MASCH_Mesh_Load::vtuPrimitive(string folderName, int rank,
	MASCH_Control& controls, MASCH_Mesh& mesh, MASCH_Variables& var){
		
		
	string saveFolderName = folderName;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	string nextToken;
	// boolBinary = false;
	boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				boolCompress = true;
			}
			break;
		}
	}
	
	
	vector<string> scal_cell_save_name;
	for(auto& item : controls.primScalarNames) {
		scal_cell_save_name.push_back(item);
	}
	vector<string> vec_cell_save_sup_name;
	vector<vector<string>> vec_cell_save_name;
	for(auto& item : controls.primVector3Names) {
		vec_cell_save_sup_name.push_back(item);
		vec_cell_save_name.push_back(vector<string>());
		vector<string> sub_names = controls.cellVar[item].sub_name;
		for(int i=0; i<sub_names.size(); ++i){
			vec_cell_save_name.back().push_back(sub_names[i]);
		}
	}
	for(auto& item : controls.primVectorNames) {
		// vec_cell_save_name.push_back(vector<string>());
		vector<string> sub_roles = controls.cellVar[item].sub_role;
		vector<string> sub_names = controls.cellVar[item].sub_name;
		for(int i=0; i<sub_names.size(); ++i){
			if(sub_roles[i]!="primitive") continue;
			string name = item;
			name += "-"+sub_names[i];
			scal_cell_save_name.push_back(name);
		}
	}
    
    
    // 평균값 로드 
    for(auto& item : controls.saveMeanCellValues){
        // cout << "AA " << item << endl;
        scal_cell_save_name.push_back("mean-"+item);
    }
    for(auto& item : controls.saveSMDValues){
        // cout << "AA " << item << endl;
        scal_cell_save_name.push_back("fvm-mean-surface-area-"+item);
        scal_cell_save_name.push_back("fvm-mean-volume-"+item);
        
        scal_cell_save_name.push_back("parcel-mean-surface-area-"+item);
        scal_cell_save_name.push_back("parcel-mean-volume-"+item);
        
        scal_cell_save_name.push_back("fvm-sauter-mean-diameter-"+item);
        scal_cell_save_name.push_back("parcel-sauter-mean-diameter-"+item);
        scal_cell_save_name.push_back("sauter-mean-diameter-"+item);
    }
    // 평균값의 시간
    vector<string> total_times;
    for(auto& item : controls.saveMeanCellValues){
        total_times.push_back("total-time-of-mean-"+item);
    }
    for(auto& item : controls.saveSMDValues){
        total_times.push_back("total-time-of-fvm-mean-surface-area-"+item);
        total_times.push_back("total-time-of-fvm-mean-volume-"+item);
        total_times.push_back("total-time-of-parcel-mean-surface-area-"+item);
        total_times.push_back("total-time-of-parcel-mean-volume-"+item);
    }
	// 필드값 로드
	for(auto& name : total_times){
		vector<double> tmp_vars;
		loadDatasAtVTU(inputFile, name, tmp_vars);
		int id = controls.getId_fieldVar(name);
        if(tmp_vars.size()>0){
            var.fields[id] = tmp_vars[0];
        }
        else{
            var.fields[id] = 0.0;
        }
	}
    
	// 셀 값 로드
	for(auto& name : scal_cell_save_name){
		vector<double> tmp_cellVars;
		loadDatasAtVTU(inputFile, name, tmp_cellVars);
		// cout << name << " " << tmp_cellVars.size() << endl;
		int id = controls.getId_cellVar(name);
		int iter = 0;
        if(tmp_cellVars.size()>0){
            for(auto& cell : var.cells){
                cell[id] = tmp_cellVars.at(iter++);
            }
        }
        else{
            for(auto& cell : var.cells){ cell[id] = 0.0; }
        }
	}
	{
		int iter_main = 0;
		for(auto& name : controls.primVector3Names){
			vector<double> tmp_cellVars;
			loadDatasAtVTU(inputFile, name, tmp_cellVars);
			// cout << name << " " << tmp_cellVars.size() << endl;
			int id0 = controls.getId_cellVar(vec_cell_save_name[iter_main][0]);
			int id1 = controls.getId_cellVar(vec_cell_save_name[iter_main][1]);
			int id2 = controls.getId_cellVar(vec_cell_save_name[iter_main][2]);
			int iter = 0;
			for(auto& cell : var.cells){
				cell[id0] = tmp_cellVars[iter++];
				cell[id1] = tmp_cellVars[iter++];
				cell[id2] = tmp_cellVars[iter++];
			}
			++iter_main;
		}
	}
	
	
	inputFile.close();
	
	
}





void MASCH_Load::dpmFiles(string folderName, int rank,
	MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
		
		
	string saveFolderName = folderName;
	string saveFileName = "parcels";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		return;
		cerr << "Unable to open file for reading : " << openFileName << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	string nextToken;
	// boolBinary = false;
	boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				boolCompress = true;
			}
			break;
		}
	}
	
	
	vector<string> scal_cell_save_name;
	// for(auto& item : controls.primScalarNames) {
		// scal_cell_save_name.push_back(item);
	// }
	vector<string> vec_cell_save_sup_name;
	vector<string> parcelVector3Names;
	vector<vector<string>> vec_cell_save_name;
	// for(auto& item : controls.primVector3Names) {
		// vec_cell_save_sup_name.push_back(item);
		// vec_cell_save_name.push_back(vector<string>());
		// vector<string> sub_names = controls.cellVar[item].sub_name;
		// for(int i=0; i<sub_names.size(); ++i){
			// vec_cell_save_name.back().push_back(sub_names[i]);
		// }
	// }
	// for(auto& item : controls.primVectorNames) {
		// // vec_cell_save_name.push_back(vector<string>());
		// vector<string> sub_roles = controls.cellVar[item].sub_role;
		// vector<string> sub_names = controls.cellVar[item].sub_name;
		// for(int i=0; i<sub_names.size(); ++i){
			// if(sub_roles[i]!="primitive") continue;
			// string name = item;
			// name += "-"+sub_names[i];
			// scal_cell_save_name.push_back(name);
		// }
	// }
	
	
	// for(auto& [key, value] : controls.parcelVar){
		// if(value.shape == "scalar"){
			// scal_cell_save_name.push_back(key);
		// }
		// else if(value.shape == "vector"){
			// for(auto& name : value.sub_name){
				// scal_cell_save_name.push_back(key+"-"+name);
			// }
		// }
		// else if(value.shape == "vector3"){
			// parcelVector3Names.push_back(key);
			// vec_cell_save_name.push_back(vector<string>());
			// for(auto& name : value.sub_name){
				// vec_cell_save_name.back().push_back(name);
			// }
		// }
		// else {
			// cout << "#WARNING : not defined parcelVar shape" << endl;
		// }
	// }
	vector<string> tmp_dummy;
	tmp_dummy.push_back("x-position");
	tmp_dummy.push_back("y-position");
	tmp_dummy.push_back("z-position");
	for(auto& [key, value] : controls.parcelVar){
		if(value.shape == "vector3"){
			if(key=="position") continue;
			parcelVector3Names.push_back(key);
			vec_cell_save_name.push_back(vector<string>());
			for(auto& name : value.sub_name){
				vec_cell_save_name.back().push_back(name);
				tmp_dummy.push_back(name);
			}
		}
	}
	// vec_save_sup_name.erase(remove(vec_save_sup_name.begin(), vec_save_sup_name.end(), "position"), vec_save_sup_name.end());
	for(auto& [key, value] : controls.parcelVar){
		if(value.shape == "scalar") scal_cell_save_name.push_back(key);
	}
	for(auto& name : tmp_dummy){
		// scal_cell_save_name.erase(std::find(scal_cell_save_name.begin(), scal_cell_save_name.end(), name));
		scal_cell_save_name.erase(remove(scal_cell_save_name.begin(), scal_cell_save_name.end(), name), scal_cell_save_name.end());
	}
	
	
	for(auto& name : scal_cell_save_name){
		vector<double> tmp_cellVars;
		loadDatasAtVTU(inputFile, name, tmp_cellVars);
		int id = controls.getId_parcelVar(name);
		int iter = 0;
        if(tmp_cellVars.size()!=0){
            for(auto& parcel : var.parcels){
                parcel[id] = tmp_cellVars[iter++];
            }
        }
        else{
            if(name=="number-of-parcel"){
                for(auto& parcel : var.parcels){
                    parcel[id] = 1.0;
                }
            }
            else{
                for(auto& parcel : var.parcels){
                    parcel[id] = 0.0;
                }
            }
        }
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	{
		int iter_main = 0;
		for(auto& name : parcelVector3Names){
			vector<double> tmp_cellVars;
			loadDatasAtVTU(inputFile, name, tmp_cellVars);
			// cout << name << " " << tmp_cellVars.size() << endl;
			int id0 = controls.getId_parcelVar(vec_cell_save_name[iter_main][0]);
			int id1 = controls.getId_parcelVar(vec_cell_save_name[iter_main][1]);
			int id2 = controls.getId_parcelVar(vec_cell_save_name[iter_main][2]);
			int iter = 0;
			for(auto& parcel : var.parcels){
				parcel[id0] = tmp_cellVars[iter++];
				parcel[id1] = tmp_cellVars[iter++];
				parcel[id2] = tmp_cellVars[iter++];
			}
			++iter_main;
		}
	}
	{
		vector<double> tmp_cellVars;
		loadDatasAtVTU(inputFile, "Position", tmp_cellVars);
		int id0 = controls.getId_parcelVar("x-position");
		int id1 = controls.getId_parcelVar("y-position");
		int id2 = controls.getId_parcelVar("z-position");
		int iter = 0;
		for(auto& parcel : var.parcels){
			parcel[id0] = tmp_cellVars[iter++];
			parcel[id1] = tmp_cellVars[iter++];
			parcel[id2] = tmp_cellVars[iter++];
		}
	}
	
	
	
	inputFile.close();
	
	
}




