
#include "./load.h" 

void MASCH_Mesh_Load::OpenFoam(string folder, MASCH_Mesh &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute load OpenFoam files ... ";
	
	// if(rank == 0)
	{
		
		string gridFolderName = folder;
		string pointsName = "points";
		string facesName = "faces";
		string ownerName = "owner";
		string neighbourName = "neighbour";
		string boundaryName = "boundary";
		
		ifstream inputFile;
		string openFileName;
		
		// points 읽기
		openFileName = gridFolderName + "/" + pointsName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		string nextToken;
		bool startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					vector<double> xyz(3,0.0);
					nextToken.erase(nextToken.find("("),1); 
					nextToken.erase(nextToken.find(")"),1); 
					stringstream sstream(nextToken);
					string word;
					char del = ' ';
					int num=0;
					while (getline(sstream, word, del)){
						xyz[num] = stold(word);
						++num;
					}
					mesh.addPoint();
					mesh.points.back().x = xyz[0];
					mesh.points.back().y = xyz[1];
					mesh.points.back().z = xyz[2];
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		

		// faces 읽기
		openFileName = gridFolderName + "/" + facesName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		startInput=false;
		bool continueInput=false;
		string saveToken;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" && !continueInput ){
					break;
				}
				else{
					if(nextToken.size()==1) continue;
					
					if(nextToken.find(")") == string::npos){
						saveToken.append(" ");
						rtrim(nextToken);
						saveToken.append(nextToken);
						continueInput = true;
						continue;
					}
					saveToken.append(" ");
					rtrim(nextToken);
					saveToken.append(nextToken);
					
					saveToken.replace(saveToken.find("("),1," ");
					saveToken.replace(saveToken.find(")"),1," ");
					istringstream iss(saveToken);
					int tempint;
					iss >> tempint;
					
					mesh.addFace();
					
					while(iss >> tempint){
						mesh.faces.back().ipoints.push_back(
						static_cast<int>(tempint));
					}
					saveToken.clear();
					continueInput = false;
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		// owner
		openFileName = gridFolderName + "/" + ownerName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		int temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					while(iss >> tempint){
						mesh.faces[temp_num].iL = static_cast<int>(tempint);
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		
		// neighbour
		openFileName = gridFolderName + "/" + neighbourName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					while(iss >> tempint){
						if(tempint<0) break;
						mesh.faces[temp_num].iR = static_cast<int>(tempint);
						// mesh.faces[temp_num].thereR = true;
						mesh.faces[temp_num].setType(MASCH_Face_Types::INTERNAL);
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		// boundary
		openFileName = gridFolderName + "/" + boundaryName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		vector<string> boundary_name;
		vector<char> boundary_type;
		vector<int> boundary_nFaces;
		vector<int> boundary_startFace;
		vector<int> boundary_myProcNo;
		vector<int> boundary_neighbProcNo;
		
		string backToken;
		startInput=false;
		vector<string> setToken;
		while(getline(inputFile, nextToken)){
			string asignToken;
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				setToken.push_back(nextToken.c_str());
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}
			backToken = nextToken;
		}
		
		string names;
		vector<string> setToken2;
		startInput=false;
		for (auto item : setToken) {
			string asignToken;
			if(startInput){
				if( item.find("}") != string::npos ){
					trim(names);
					
					boundary_name.push_back(names);
					int l=0;
					for (auto item2 : setToken2) {
						if( item2.find("nFaces") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_nFaces.push_back(temptempint);
						}
						if( item2.find("startFace") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_startFace.push_back(temptempint);
						}
						if( item2.find("myProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_myProcNo.push_back(temptempint);
							++l;
						}
						if( item2.find("neighbProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_neighbProcNo.push_back(temptempint);
							++l;
						}
					}
					if(l==0) {
						boundary_myProcNo.push_back(-1);
						boundary_neighbProcNo.push_back(-1);
					}
					
					startInput=false;
					setToken2.clear();
				}
				setToken2.push_back(item.c_str());
			}
			else{ 
				if( item.find("{") != string::npos ){
					names.clear();
					names = backToken;
					startInput=true;
				}
			}
			backToken = item;
		}
		
		
		inputFile.close();
		mesh.boundaries.clear();
		int nbcs=boundary_name.size();
		for (int i=0; i<boundary_name.size(); ++i) {
			mesh.addBoundary();
		}
		for (int i=0; i<mesh.boundaries.size(); ++i) {
			mesh.boundaries[i].name = trim(boundary_name[i]);
			mesh.boundaries[i].nFaces = static_cast<int>(boundary_nFaces[i]);
			mesh.boundaries[i].startFace = static_cast<int>(boundary_startFace[i]);
			mesh.boundaries[i].myProcNo = static_cast<int>(boundary_myProcNo[i]);
			mesh.boundaries[i].rightProcNo = static_cast<int>(boundary_neighbProcNo[i]);
			// mesh.boundaries[i].thereR = false;
			mesh.boundaries[i].setType(MASCH_Face_Types::BOUNDARY);
		}
	}
	
	if(rank!=0){
		mesh.points.clear();
		mesh.faces.clear();
		mesh.cells.clear();
		for (int i=0; i<mesh.boundaries.size(); ++i) {
			mesh.boundaries[i].nFaces = 0;
			mesh.boundaries[i].startFace = 0;
		}
	}
	
	
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	

	mesh.check();
	mesh.setFaceTypes();
	mesh.buildCells();
	mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	mesh.informations();
	
	
}

