
#include "./load.h" 

void MASCH_Load::meshFiles(string folderName, MASCH_Control& controls, MASCH_Mesh& mesh){
	
	// MASCH_Mesh_Load mesh_load;
	
	(*this).vtu(folderName, controls, mesh);
	
}


void MASCH_Mesh_Load::vtu(
	string folder, 
	MASCH_Control& controls, 
	MASCH_Mesh& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load vtu data files ... ";
	}
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	// points
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
	
	vector<double> timevalue;
	loadDatasAtVTU(inputFile, "TimeValue", timevalue);
	if(timevalue.size()==0) timevalue.push_back(0.0);
	
	vector<int> pointLevels;
	loadDatasAtVTU(inputFile, "pointLevels", pointLevels);
	
	vector<int> cellLevels;
	loadDatasAtVTU(inputFile, "cellLevels", cellLevels);
	
	vector<int> cellGroups;
	loadDatasAtVTU(inputFile, "cellGroups", cellGroups);
	
	vector<double> NodeCoordinates;
	loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	
	vector<int> connectivity;
	loadDatasAtVTU(inputFile, "connectivity", connectivity);
	
	vector<int> offsets;
	loadDatasAtVTU(inputFile, "offsets", offsets);
	
	vector<int> read_ifaces;
	loadDatasAtVTU(inputFile, "faces", read_ifaces);
	
	vector<int> faceoffsets;
	loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
	
	vector<int> read_iL;
	loadDatasAtVTU(inputFile, "owner", read_iL);
	
	vector<int> read_iR;
	loadDatasAtVTU(inputFile, "neighbour", read_iR);
	
	vector<string> bcName;
	loadDatasAtVTU(inputFile, "bcName", bcName);
	
	vector<int> bcStartFace;
	loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
	
	vector<int> bcNFaces;
	loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
	
	vector<int> bcNeighbProcNo;
	loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);

	vector<int> connPoints;
	loadDatasAtVTU(inputFile, "connPoints", connPoints);

	inputFile.close();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "B1" << endl;
	
	
	for(int i=0; i<NodeCoordinates.size()/3; ++i){
		mesh.addPoint();
		mesh.points.back().x = NodeCoordinates.at(i*3+0);
		mesh.points.back().y = NodeCoordinates.at(i*3+1);
		mesh.points.back().z = NodeCoordinates.at(i*3+2);
		mesh.points.back().level = 0;
	}
	NodeCoordinates.clear();
	
	// 포인트 레벨
	for(int i=0; i<pointLevels.size(); ++i){
		mesh.points.at(i).level = pointLevels[i];
	}
	
	
	int ncells=-1;
	mesh.faces.clear();
	for(int i=0; i<read_iL.size(); ++i){
		int tmp_iL = read_iL[i];
		mesh.addFace();
		mesh.faces.back().iL = tmp_iL;
		ncells = max(ncells, tmp_iL);
	}
	
	
	// cout << ncells << endl;
	
	vector<int> str_icells(size+1,0);
	if(size>1){
		int send_ncells = ncells+1;
		vector<int> recv_ncells(size);
		MPI_Allgather(&send_ncells, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
		for(int ip=0; ip<size; ++ip){
			str_icells[ip+1] = str_icells[ip] + recv_ncells[ip];
		}
	}
	
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();
		mesh.cells.back().level = 0;
		mesh.cells.back().group = str_icells[rank] + i;
	}
	
	// 셀 레벨
	for(int i=0; i<cellLevels.size(); ++i){
		mesh.cells.at(i).level = cellLevels.at(i);
	}
	// 셀 그룹
	for(int i=0; i<cellGroups.size(); ++i){
		mesh.cells.at(i).group = cellGroups.at(i);
	}
	
	
	for(int i=0; i<read_iR.size(); ++i){
		int tmp_iR = read_iR[i];
		mesh.faces.at(i).iR = tmp_iR;
		mesh.faces.at(i).setType(MASCH_Face_Types::INTERNAL);
	}
	for(int i=read_iR.size(); i<mesh.faces.size(); ++i){
		mesh.faces.at(i).setType(MASCH_Face_Types::PROCESSOR);
	}
	
	
	for(int i=0, str=0, iter=0; i<offsets.size(); ++i){
		int tmp_offset = offsets[i];
		for(int j=str; j<tmp_offset; ++j){
			int ipoint = connectivity[j];
			mesh.cells.at(iter).ipoints.push_back( ipoint );
		}
		str=tmp_offset;
		++iter;
	}
	
	
	
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			mesh.cells[ face.iL ].ifaces.push_back( i );
			mesh.cells[ face.iR ].ifaces.push_back( i );
		}
		else{
			mesh.cells[ face.iL ].ifaces.push_back( i );
		}
	}
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// cout << mesh.cells.size() << " " << faceoffsets.size() << endl;
	
	for(int i=0, str=0; i<faceoffsets.size(); ++i){
		int face_size = read_ifaces.at(str++);
		for(auto& iface : mesh.cells[i].ifaces){
			int point_size = read_ifaces.at(str++);
			for(int k=0; k<point_size; ++k){
				int ipoint = read_ifaces.at(str++);
				if(mesh.faces[ iface ].ipoints.size() == point_size) continue;
				mesh.faces[ iface ].ipoints.push_back( ipoint );
			}
		}
		if(str!=faceoffsets[i]){
		cout << i << " " << face_size << " " << str << " " << faceoffsets[i] << endl;
		}
	}
	
	for(int i=0; i<bcStartFace.size(); ++i){
		int startFace = bcStartFace[i];
		
		mesh.addBoundary();
		
		string tmp_bcName = bcName[i];
		tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
		tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
		tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
		tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
		mesh.boundaries.back().name = tmp_bcName;
		mesh.boundaries.back().startFace = bcStartFace[i];
		mesh.boundaries.back().nFaces = bcNFaces[i];
		if(bcNeighbProcNo[i]<0){
			mesh.boundaries.back().rightProcNo = -1;
			mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
		}
		else{
			mesh.boundaries.back().rightProcNo = bcNeighbProcNo[i];
			mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
		}
		mesh.boundaries.back().myProcNo = rank;
		
	}
	
	int maxBCNum = mesh.boundaries.size()-1;
	mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
	for(int i=maxBCNum-1; i>=0; --i){
		mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
	}
	
	for(int i=0; i<connPoints.size(); ++i){
		mesh.points[connPoints[i]].connPoints.push_back(
		make_pair(connPoints[i+1],connPoints[i+2])
		);
		++i; ++i;
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "B2" << endl;
	
	
	// // mesh.set(NodeCoordinates, connectivity, offsets, faces, faceoffsets,
		// // iL, iR, bcName, bcStartFace, bcNFaces, bcNeighbProcNo, connPoints,
		// // pointLevels, cellLevels, cellGroups);
	
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	mesh.setFaceLevels();
	mesh.informations();

	MPI_Barrier(MPI_COMM_WORLD);

		
}






