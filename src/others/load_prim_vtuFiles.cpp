
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
	
	for(auto& name : scal_cell_save_name){
		vector<double> tmp_cellVars;
		loadDatasAtVTU(inputFile, name, tmp_cellVars);
		// cout << name << " " << tmp_cellVars.size() << endl;
		int id = controls.getId_cellVar(name);
		int iter = 0;
		for(auto& cell : var.cells){
			cell[id] = tmp_cellVars[iter++];
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




