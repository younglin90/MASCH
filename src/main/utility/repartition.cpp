#include <iostream>
// #include <algorithm>
// #include <cmath>
// #include <mpi.h>
// #include <iomanip>
// #include "parmetis.h" 
// #include "scotch.h" 

#include "../../others/mesh.h"
#include "../../others/controls.h"
#include "../../others/save.h"
#include "../../others/load.h"
#include "../../others/variables.h"

void print_help();

// 필요한 객체들 생성
MASCH_Mesh_Save save;
MASCH_Control controls;
MASCH_Load load;
MASCH_Mesh mesh;
MASCH_Variables var;
vector<string> primScalarNames;
vector<string> primVector3Names;
vector<string> parcel_primScalarNames;
vector<string> parcel_primVector3Names;

vector<int> cell_ip_g;
vector<int> parcel_ip_g;

int varSize;
int varScalarSize;
int varVector3Size;
// int parcelVarSize;

void loadVTUfiles(string folderName){

	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load vtu data files ... ";
	}
		
	string saveFolderName = folderName;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	
	bool presentFiles = true;
	
	if(inputFile.fail()){
		presentFiles = false;
		// cerr << "Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	inputFile.close();
	
	
	string nextToken;
	load.boolCompress = false;
	vector<double> timevalue;
	vector<int> pointLevels;
	vector<int> vtkGhostType;
	vector<int> cellLevels;
	vector<int> cellGroups;
	vector<double> NodeCoordinates;
	vector<int> connectivity;
	vector<int> offsets;
	vector<int> read_ifaces;
	vector<int> faceoffsets;
	vector<int> read_iL;
	vector<int> read_iR;
	vector<string> bcName;
	vector<int> bcStartFace;
	vector<int> bcNFaces;
	vector<int> connPoints;
	vector<int> bcNeighbProcNo;
	bool ghostTypePres = false;
	bool pointLevelsPres = false;
	int ncells=-1;
	// if(presentFiles)
	{
		string tmp_openFileName = saveFolderName + "/" + saveFileName + ".0.vtu";
		inputFile.open(tmp_openFileName);
	
		while(getline(inputFile, nextToken)){
			if( nextToken.find("VTKFile") != string::npos ){
				if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
					load.boolCompress = true;
				}
				break;
			}
		}
		
		while(getline(inputFile, nextToken)){
			if( nextToken.find("DataArray") != string::npos ){
				if( nextToken.find("Name") != string::npos ){
					if( nextToken.find("NumberOfTuples") != string::npos ) continue;
					
					if( nextToken.find("NumberOfComponents") != string::npos ){
						// string tmp_nextToken = nextToken;
						int str = nextToken.find("Name")+6;
						int end = nextToken.find("NumberOfComponents")-2;
						string tmp_name = nextToken.substr(str,end-str);
						primVector3Names.push_back(tmp_name);
						// if(rank==0) cout << tmp_name << endl;
						// if(rank==0) cout << endl;
						// for(int i=str; i<end; ++i){
							// if(rank==0) cout << tmp_nextToken[i];
						// }
						// if(rank==0) cout << endl;
					}
					else{
						// string tmp_nextToken = nextToken;
						int str = nextToken.find("Name")+6;
						int end = nextToken.find("format")-2;
						string tmp_name = nextToken.substr(str,end-str);
						primScalarNames.push_back(tmp_name);
						// if(rank==0) cout << tmp_name << endl;
						// if(rank==0) cout << endl;
						// for(int i=str; i<end; ++i){
							// if(rank==0) cout << tmp_nextToken[i];
						// }
						// if(rank==0) cout << endl;
					}
				}
			}
		}
		
		
		
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "pointLevels"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "vtkGhostType"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "cellLevels"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "cellGroups"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "connectivity"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "offsets"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "types"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "faces"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "faceoffsets"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "owner"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "neighbour"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "bcName"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "bcStartFace"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "bcNFaces"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "bcNeighbProcNo"), primScalarNames.end());
		primScalarNames.erase(remove(primScalarNames.begin(), primScalarNames.end(), "connPoints"), primScalarNames.end());
		
		primVector3Names.erase(remove(primVector3Names.begin(), primVector3Names.end(), "NodeCoordinates"), primVector3Names.end());
		
		// for(auto& item : primScalarNames){
			// if(rank==0) cout << item << endl;
		// }
		// for(auto& item : primVector3Names){
			// if(rank==0) cout << item << endl;
		// }
		
		
		inputFile.close();
	}
		
		
	if(presentFiles){
		
		


		// string saveFolderName = folder;
		// string saveFileName = "plot";
		// string saveRankName = to_string(rank);
		
		// ifstream inputFile;
		// string openFileName;
		
		// points
		// openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
		inputFile.open(openFileName);
		if(inputFile.fail()){
			// cerr << "Unable to open file for reading : " << openFileName << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		load.loadDatasAtVTU(inputFile, "TimeValue", timevalue);
		if(timevalue.size()==0) timevalue.push_back(0.0);
		
		load.loadDatasAtVTU(inputFile, "pointLevels", pointLevels);
		
		load.loadDatasAtVTU(inputFile, "vtkGhostType", vtkGhostType);
		
		load.loadDatasAtVTU(inputFile, "cellLevels", cellLevels);
		
		load.loadDatasAtVTU(inputFile, "cellGroups", cellGroups);
		
		load.loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
		
		load.loadDatasAtVTU(inputFile, "connectivity", connectivity);
		
		load.loadDatasAtVTU(inputFile, "offsets", offsets);
		
		load.loadDatasAtVTU(inputFile, "faces", read_ifaces);
		
		load.loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
		
		load.loadDatasAtVTU(inputFile, "owner", read_iL);
		
		load.loadDatasAtVTU(inputFile, "neighbour", read_iR);
		
		load.loadDatasAtVTU(inputFile, "bcName", bcName);
		
		load.loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
		
		load.loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
		
		load.loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);

		load.loadDatasAtVTU(inputFile, "connPoints", connPoints);

		inputFile.close();
		
		
		if(vtkGhostType.size()>0) ghostTypePres = true;
		if(pointLevels.size()>0) pointLevelsPres = true;
		
		// cout << NodeCoordinates.size() << endl;
		for(int i=0; i<NodeCoordinates.size()/3; ++i){
			if(pointLevelsPres==true) { if(pointLevels[i]==-100) break; }
			mesh.addPoint();
			mesh.points.back().x = NodeCoordinates.at(i*3+0);
			mesh.points.back().y = NodeCoordinates.at(i*3+1);
			mesh.points.back().z = NodeCoordinates.at(i*3+2);
			mesh.points.back().level = 0;
			// mesh.points.back().level = pointLevels[i];
		}
		NodeCoordinates.clear();
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// 포인트 레벨
		for(int i=0; i<pointLevels.size(); ++i){
			if(pointLevelsPres==true) { if(pointLevels[i]==-100) break; }
			mesh.points.at(i).level = pointLevels[i];
		}
		
		
		mesh.faces.clear();
		for(int i=0; i<read_iL.size(); ++i){
			// if(ghostTypePres==true) { 
				// if(vtkGhostType[i]>0) {
					// cout << "AAAAAA" << endl;
					// break;
				// }
			// }
			int tmp_iL = read_iL[i];
			mesh.addFace();
			mesh.faces.back().iL = tmp_iL;
			ncells = max(ncells, tmp_iL);
		}
		
	}
	
	
	
	vector<int> str_icells(size+1,0);
	{
		if(size>1){
			int send_ncells = ncells+1;
			vector<int> recv_ncells(size);
			MPI_Allgather(&send_ncells, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
			for(int ip=0; ip<size; ++ip){
				str_icells[ip+1] = str_icells[ip] + recv_ncells[ip];
			}
		}
	}
		
	if(presentFiles){
		mesh.cells.clear();
		for(int i=0; i<ncells+1; ++i){
			mesh.addCell();
			mesh.cells.back().level = 0;
			mesh.cells.back().group = str_icells[rank] + i;
		}
		
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// 셀 레벨
		for(int i=0; i<cellLevels.size(); ++i){
			if(ghostTypePres==true) { if(vtkGhostType[i]>0) break; }
			mesh.cells.at(i).level = cellLevels.at(i);
		}
		// 셀 그룹
		for(int i=0; i<cellGroups.size(); ++i){
			if(ghostTypePres==true) { if(vtkGhostType[i]>0) break; }
			mesh.cells.at(i).group = cellGroups.at(i);
		}
		
		// cout << mesh.faces.size() << " " << read_iR.size() << endl;
		
		for(int i=0; i<read_iR.size(); ++i){
			int tmp_iR = read_iR[i];
			mesh.faces.at(i).iR = tmp_iR;
			mesh.faces.at(i).setType(MASCH_Face_Types::INTERNAL);
		}
		for(int i=read_iR.size(); i<mesh.faces.size(); ++i){
			mesh.faces.at(i).setType(MASCH_Face_Types::PROCESSOR);
		}
		
		
		for(int i=0, str=0, iter=0; i<offsets.size(); ++i){
			if(ghostTypePres==true) { if(vtkGhostType[i]>0) break; }
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
		
		
		
		
		// cout << mesh.cells.size() << " " << faceoffsets.size() << endl;
		
		for(int i=0, str=0; i<faceoffsets.size(); ++i){
			if(ghostTypePres==true) { if(vtkGhostType[i]>0) break; }
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
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// if(rank==0) cout << "B2" << endl;
		
		
		// // mesh.set(NodeCoordinates, connectivity, offsets, faces, faceoffsets,
			// // iL, iR, bcName, bcStartFace, bcNFaces, bcNeighbProcNo, connPoints,
			// // pointLevels, cellLevels, cellGroups);
		
			
		// MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			cout << "-> completed" << endl;
			cout << "└────────────────────────────────────────────────────" << endl;
		}

		// MPI_Barrier(MPI_COMM_WORLD);

		
		
		inputFile.close();
		
		
	}
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	mesh.setFaceLevels();
	mesh.informations();
	
	
	{
		int tmp_varSize = primScalarNames.size() + 3*primVector3Names.size();
		
		var.cells.resize(mesh.cells.size());
		for(auto& cell : var.cells){
			cell.resize(tmp_varSize);
		}
	}
	
		
	if(presentFiles){
		inputFile.open(openFileName);
		
		while(getline(inputFile, nextToken)){
			if( nextToken.find("VTKFile") != string::npos ){
				if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
					load.boolCompress = true;
				}
				break;
			}
		}
		
		int id = 0;
		for(auto& name : primScalarNames){
			vector<double> tmp_cellVars;
			load.loadDatasAtVTU(inputFile, name, tmp_cellVars);
			int iter = 0;
			// if(rank==0) cout << "AAAAAAA " << mesh.cells.size() << " " << tmp_cellVars.size() << endl;
			for(auto& cell : var.cells){
				cell[id] = tmp_cellVars[iter++];
			}
			++id;
		}
		{
			int iter_main = 0;
			for(auto& name : primVector3Names){
				vector<double> tmp_cellVars;
				load.loadDatasAtVTU(inputFile, name, tmp_cellVars);
				int iter = 0;
				for(auto& cell : var.cells){
					cell[id] = tmp_cellVars[iter++];
					cell[id+1] = tmp_cellVars[iter++];
					cell[id+2] = tmp_cellVars[iter++];
				}
				++iter_main;
				++id; ++id; ++id;
			}
		}
		
		
		
		inputFile.close();
		
	}
	
}





void loadDPMFiles(string folderName){
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	{
		string saveFolderName = folderName;
		string saveFileName = "parcels";
		string saveRankName = to_string(0);
		
		ifstream inputFile;
		string openFileName;
		
		openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
		inputFile.open(openFileName);
		if(inputFile.fail()){
			return;
			// cerr << "Unable to open file for reading : " << openFileName << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		string nextToken;
		// boolBinary = false;
		load.boolCompress = false;
		while(getline(inputFile, nextToken)){
			if( nextToken.find("VTKFile") != string::npos ){
				if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
					load.boolCompress = true;
				}
				break;
			}
		}
		

		while(getline(inputFile, nextToken)){
			if( nextToken.find("DataArray") != string::npos ){
				if( nextToken.find("Name") != string::npos ){
					if( nextToken.find("NumberOfTuples") != string::npos ) continue;
					
					if( nextToken.find("NumberOfComponents") != string::npos ){
						int str = nextToken.find("Name")+6;
						int end = nextToken.find("NumberOfComponents")-2;
						string tmp_name = nextToken.substr(str,end-str);
						parcel_primVector3Names.push_back(tmp_name);
					}
					else{
						int str = nextToken.find("Name")+6;
						int end = nextToken.find("format")-2;
						string tmp_name = nextToken.substr(str,end-str);
						parcel_primScalarNames.push_back(tmp_name);
					}
				}
			}
		}
		
		inputFile.close();
		
		
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "id"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "icell"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "connectivity"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "offsets"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "types"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "faces"), parcel_primScalarNames.end());
		parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "faceoffsets"), parcel_primScalarNames.end());
		
		parcel_primVector3Names.erase(remove(parcel_primVector3Names.begin(), parcel_primVector3Names.end(), "Position"), parcel_primVector3Names.end());
			
	}
	
	
	
		
	string saveFolderName = folderName;
	string saveFileName = "parcels";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		return;
		// cerr << "Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	string nextToken;
	// boolBinary = false;
	load.boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				load.boolCompress = true;
			}
			break;
		}
	}
	

	// while(getline(inputFile, nextToken)){
		// if( nextToken.find("DataArray") != string::npos ){
			// if( nextToken.find("Name") != string::npos ){
				// if( nextToken.find("NumberOfTuples") != string::npos ) continue;
				
				// if( nextToken.find("NumberOfComponents") != string::npos ){
					// int str = nextToken.find("Name")+6;
					// int end = nextToken.find("NumberOfComponents")-2;
					// string tmp_name = nextToken.substr(str,end-str);
					// parcel_primVector3Names.push_back(tmp_name);
				// }
				// else{
					// int str = nextToken.find("Name")+6;
					// int end = nextToken.find("format")-2;
					// string tmp_name = nextToken.substr(str,end-str);
					// parcel_primScalarNames.push_back(tmp_name);
				// }
			// }
		// }
	// }
	
	
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "id"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "icell"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "connectivity"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "offsets"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "types"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "faces"), parcel_primScalarNames.end());
	// parcel_primScalarNames.erase(remove(parcel_primScalarNames.begin(), parcel_primScalarNames.end(), "faceoffsets"), parcel_primScalarNames.end());
	
	// parcel_primVector3Names.erase(remove(parcel_primVector3Names.begin(), parcel_primVector3Names.end(), "Position"), parcel_primVector3Names.end());
		
	
	
	vector<int> parcel_id;
	load.loadDatasAtVTU(inputFile, "id", parcel_id);
	vector<int> parcel_icell;
	load.loadDatasAtVTU(inputFile, "icell", parcel_icell);
	
	int NumberOfPoints = parcel_id.size();
	
	mesh.parcels.resize(NumberOfPoints);
	{
		int iter = 0;
		for(auto& parcel : mesh.parcels){
			parcel.id = parcel_id[iter];
			parcel.icell = parcel_icell[iter];
			parcel.setType(MASCH_Parcel_Types::INSIDE);
			++iter;
		}
	}
	
	
	var.parcels.resize(mesh.parcels.size());
	
	// int parcelVarSize = controls.nParcelVar;
	int parcelVarSize = parcel_primScalarNames.size()+3*parcel_primVector3Names.size()+3;

	
	for(auto& parcel : var.parcels){
		parcel.resize(parcelVarSize);
	}
	
	
	int id = 0;
	for(auto& name : parcel_primScalarNames){
		vector<double> tmp_cellVars;
		load.loadDatasAtVTU(inputFile, name, tmp_cellVars);
		int iter = 0;
		for(auto& parcel : var.parcels){
			parcel[id] = tmp_cellVars[iter++];
		}
		++id;
	}
	{
		int iter_main = 0;
		for(auto& name : parcel_primVector3Names){
			vector<double> tmp_cellVars;
			load.loadDatasAtVTU(inputFile, name, tmp_cellVars);
			int id0 = parcel_primScalarNames.size() + 3*iter_main;
			int id1 = id0+1;
			int id2 = id0+2;
			int iter = 0;
			for(auto& parcel : var.parcels){
				parcel[id0] = tmp_cellVars[iter++];
				parcel[id1] = tmp_cellVars[iter++];
				parcel[id2] = tmp_cellVars[iter++];
			}
			++iter_main;
			++id; ++id; ++id;
		}
	}
	{
		vector<double> tmp_cellVars;
		load.loadDatasAtVTU(inputFile, "Position", tmp_cellVars);
		int id0 = id;
		int id1 = id+1;
		int id2 = id+2;
		int iter = 0;
		for(auto& parcel : var.parcels){
			parcel[id0] = tmp_cellVars[iter++];
			parcel[id1] = tmp_cellVars[iter++];
			parcel[id2] = tmp_cellVars[iter++];
		}
	}
	
	
	
	inputFile.close();
	
}










void repartParMETIS(int nSizeOrg, int nSizeTar){

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute ParMETIS ... ";
	}
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	options[0] = 0;
	options[1] = 0;
	options[2] = 0;
	
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	int wgtflag=0;
	int numflag=0;
	vector<real_t> tpwgts(size*ncon);
	for(int i=0;i<nSizeTar*ncon;++i) tpwgts[i]=1.0/(double)nSizeTar;
	for(int i=nSizeTar*ncon;i<size*ncon;++i) tpwgts[i]=0.0;

	// tpwgts[60] = 0.0;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.02);
	
	real_t ubvec = 1.02;
	
	
	bool boolZeroCell = false;
	if(ncells==0) {
		ncells = 1;
		boolZeroCell=true;
	}
	
	
	vector<int> temp_vtxdist(size,0);
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	// vector<int> vtxdist(nSize+1);
	vector<int> vtxdist(size+1,0);
	for(int i=1; i<size+1; ++i){
		vtxdist.at(i) = vtxdist.at(i-1) + (temp_vtxdist.at(i-1));
	}
	
	
	vector<int> xadj(ncells+1,0);
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		// for(auto& face : cell.faces()){
			// if(
			// (*face).getType() == MASCH_Face_Types::INTERNAL ||
			// (*face).getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		// }
		for(auto& iface : cell.ifaces){
			auto& face = mesh.faces[iface];
			if(
			face.getType() == MASCH_Face_Types::INTERNAL ||
			face.getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		}
		++numt;
		xadj.at(numt) = xadj.at(numt-1) + (numt2);
	}
	
	if(boolZeroCell){
		xadj[1] = 1;
	}
	int saveOrg = 0;
	if(rank==nSizeOrg-1){
		saveOrg = xadj[ncells];
		xadj[ncells] += size-nSizeOrg;
	}
	
	
	vector<int> adjncy(xadj.at(ncells),0);
	// if(ncells==0) adjncy.resize(1,0);
	
	
	vector<int> recv_right_rank;
	vector<int> recv_right_icell;
	if(size>1){
		vector<int> send_right_rank;
		vector<int> send_right_icell;
		for(auto& face : mesh.faces){
			if(face.getType() == MASCH_Face_Types::PROCESSOR) {
				send_right_rank.push_back(rank);
				send_right_icell.push_back(face.iL);
			}
		}
		recv_right_rank.resize(send_right_rank.size());
		recv_right_icell.resize(send_right_icell.size());
		MPI_Alltoallv( send_right_rank.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_right_rank.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( send_right_icell.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_right_icell.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	vector<int> num_ncell(ncells,0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL) {
			
			int tnum = xadj.at(face.iL) + num_ncell.at(face.iL);
			adjncy.at(tnum) = vtxdist.at(rank) + face.iR;
			++num_ncell.at(face.iL);
			
			tnum = xadj.at(face.iR) + num_ncell.at(face.iR);
			adjncy.at(tnum) = vtxdist.at(rank) + face.iL;
			++num_ncell.at(face.iR);
			
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR) {
			
			int tnum = xadj.at(face.iL) + num_ncell.at(face.iL);
			int right_rank = recv_right_rank.at(proc_num);
			adjncy.at(tnum) = vtxdist.at(right_rank) + recv_right_icell.at(proc_num);
			++num_ncell.at(face.iL);
			
			++proc_num;
		}
	}
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// for(auto& item : adjncy){
		// cout << item << endl;
	// }

	cell_ip_g.resize(ncells);
	
	
	if(boolZeroCell){
		adjncy[0] = vtxdist[nSizeOrg]-1;
	}
	if(rank==nSizeOrg-1){
		int iter2 = 0;
		for(int i=saveOrg; i<adjncy.size(); ++i){
			adjncy[i] = vtxdist[nSizeOrg+iter2];
			++iter2;
		}
	}
	
	
	// if(rank==60){
		// for(auto& item : adjncy){
			// cout << item << endl;
		// }
		// adjncy[0] = 100;
	// }


	ParMETIS_V3_PartKway(
		vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, &wgtflag, &numflag,
		&ncon, &nSizeTar, tpwgts.data(), &ubvec,
		options, &objval, cell_ip_g.data(), &comm);
		
		
		
	// for(auto& item : cell_ip_g){
		// if(rank==60) cout << item << endl;
	// }
		
	// idx_t vsize = NULL;
	// real_t itr = 1000.0;
	// ParMETIS_V3_AdaptiveRepart(
		// vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, &wgtflag, &numflag,
		// &ncon, &nSizeTar, tpwgts.data(), &ubvec, &itr, 
		// options, &objval, cell_ip_g.data(), &comm);
		
	
	// 그룹 재정립
	vector<int> nIps(size,0);
	for(int i=0; i<mesh.cells.size(); ++i){
		int group0 = mesh.cells[i].group;
		int cell_ip_g0 = cell_ip_g[i];
		int group = group0;
		vector<int> tmp_icell;
		tmp_icell.push_back(i);
		vector<int> tmp_ip;
		tmp_ip.push_back(cell_ip_g[i]); ++nIps[cell_ip_g[i]];
		int ip_max = cell_ip_g[i];
		while(1){
			++i;
			if(i==mesh.cells.size()) break;
			group = mesh.cells[i].group;
			if(group0!=group) break;
			tmp_icell.push_back(i);
			tmp_ip.push_back(cell_ip_g[i]); ++nIps[cell_ip_g[i]];
			if(nIps[cell_ip_g[i-1]]<nIps[cell_ip_g[i]]) ip_max = cell_ip_g[i];
		}
		--i;
		
		for(auto& icell : tmp_icell){
			cell_ip_g[icell] = ip_max;
		}
	}
	
	
	varSize = primScalarNames.size() + 3*primVector3Names.size();
	varScalarSize = primScalarNames.size();
	varVector3Size = primVector3Names.size();

	// if(size>1){
		// int nbc_glo;
		// MPI_Allreduce(&varSize, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		// varSize = nbc_glo;
		// MPI_Allreduce(&varScalarSize, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		// varScalarSize = nbc_glo;
		// MPI_Allreduce(&varVector3Size, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		// varVector3Size = nbc_glo;
	// }

	
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}







void resetVariableArray(vector<int>& cellConn){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();

	
	// int org_nCells = cellConn.size();
	int org_nCells = var.cells.size();
	int org_nFaces = var.faces.size();
	int org_nProcFaces = var.cells.size();
	int org_nPoints = var.points.size();
	int org_nBoundaries = var.boundaries.size();
	int org_nParcels = var.parcels.size();
	
	
	// var.cells.resize(org_nCells);
	// for(auto& cell : var.cells){
		// cell.resize(varSize);
	// }
	
	
	// vector<vector<double>> send_val(varSize);
	// vector<vector<double>> send_parcel_val(parcelVarSize);
	// vector<int> send_parcel_size;
	// for(int i=0; i<org_nCells; ++i){
		// for(int iprim=0; iprim<varSize; ++iprim){
			// double value = var.cells.at(i).at(iprim);
			// send_val[iprim].push_back(value);
		// }
		// send_parcel_size.push_back(var.parcels.at(i).size());
		// // if(var.parcels.at(i).size()>0){
			// // for(int iprim=0; iprim<parcelVarSize; ++iprim){
				// // double value = var.parcels.at(i).at(iprim);
				// // send_parcel_val[iprim].push_back(value);
			// // }
		// // }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	int nCells = mesh.cells.size();
	int nFaces = mesh.faces.size();
	int nProcFaces = mesh.nProcessorFaces;
	int nPoints = mesh.points.size();
	int nBoundaries = mesh.boundaries.size();
	int nParcels = mesh.parcels.size();
	
	
	
	
	int nCellsMax = max(nCells,org_nCells);
	int nFacesMax = max(nFaces,org_nFaces);
	int nProcFacesMax = max(nProcFaces,org_nProcFaces);
	int nPointsMax = max(nPoints,org_nPoints);
	int nBoundariesMax = max(nBoundaries,org_nBoundaries);
	int nParcelsMax = max(nParcels,org_nParcels);
	
	
	
	var.cells.resize(nCellsMax);
	var.faces.resize(nFacesMax);
	var.procRightCells.resize(nProcFacesMax);
	var.points.resize(nPointsMax);
	var.boundaries.resize(nBoundariesMax);
	var.parcels.resize(nParcelsMax);
	
	
	// if(size>1){
		// int nbc_glo;
		// MPI_Allreduce(&parcelVarSize, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		// parcelVarSize = nbc_glo;
	// }
	int parcelVarSize = parcel_primScalarNames.size()+3*parcel_primVector3Names.size()+3;
	
	
	// int nVarBoundary = 0;
	// for(auto& [name, tmp_var] : controls.cellVar){
		// if(tmp_var.id>=0 && tmp_var.role=="primitive"){
			// ++nVarBoundary;
		// }
	// }
	
	// var.fields.resize(controls.nFieldVar);
	for(auto& item : var.cells){
		item.resize(varSize);
	}
	// if(controls.nCellVar>0){
		// for(int i=0; i<nCells; ++i) var.cells[i].resize(varSize);
		// // for(int i=0; i<nProcFaces; ++i) var.procRightCells[i].resize(controls.nCellVar);
	// }
	// if(controls.nFaceVar>0){
		// for(int i=0; i<nFaces; ++i) var.faces[i].resize(controls.nFaceVar);
	// }
	// if(controls.nPointVar>0){
		// for(int i=0; i<nPoints; ++i) var.points[i].resize(controls.nPointVar);
	// }
	// for(int i=0; i<nBoundaries; ++i) var.boundaries[i].resize(nVarBoundary);
	for(auto& item : var.parcels){
		item.resize(parcelVarSize);
	}
	// for(int i=0; i<nParcels; ++i) var.parcels[i].resize(parcelVarSize);

	
	MASCH_MPI mpi;
	for(int iprim=0; iprim<varSize; ++iprim){
		vector<vector<double>> send_value(size);
		for(int i=0; i<org_nCells; ++i){
			double org_value = var.cells.at(i).at(iprim);
			send_value[cell_ip_g[i]].push_back(org_value);
		}
		vector<vector<double>> recv_value(size);
		mpi.Alltoallv(send_value, recv_value);
		
		for(int ip=0, iter=0; ip<size; ++ip){
			for(auto& item : recv_value[ip]){
				var.cells.at(iter).at(iprim) = item;
				++iter;
			}
		}
	}
	
	// if(rank==0) cout << parcelVarSize << endl;
	// cout << org_nParcels << endl;
	
	for(int iprim=0; iprim<parcelVarSize; ++iprim){
		vector<vector<double>> send_value(size);
		for(int i=0; i<org_nParcels; ++i){
			double org_value = var.parcels.at(i).at(iprim);
			// cout << org_value << endl;
			send_value[parcel_ip_g[i]].push_back(org_value);
		}
		vector<vector<double>> recv_value(size);
		mpi.Alltoallv(send_value, recv_value);
		
		for(int ip=0, iter=0; ip<size; ++ip){
			for(auto& item : recv_value[ip]){
				var.parcels.at(iter).at(iprim) = item;
				++iter;
			}
		}
	}
	
	// 재정립
	var.cells.resize(nCells);
	var.faces.resize(nFaces);
	var.procRightCells.resize(nProcFaces);
	var.points.resize(nPoints);
	var.boundaries.resize(nBoundaries);
	var.parcels.resize(nParcels);
	
	// for(auto& item : var.cells){
		// item.resize(varSize);
	// }
	// for(auto& item : var.parcels){
		// item.resize(parcelVarSize);
	// }
	
	
}












void saveVTUFiles(string folder, int rank){
	
	// int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// 폴더 만들기
    // auto ret = filesystem::create_directories(folder);
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	save.mkdirs(folder_name);
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	if(rank==0){
		string printPlot = folder + "plot.proc.vtu";
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << printPlot << ") ... ";
	}
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	string saveFormat = controls.saveFormat;
	
	
	
	
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;

	if(controls.saveFormat == "ascii"){
		outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		}
		else{
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
	
	outputFile << "  <UnstructuredGrid>" << endl;
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);


	// vtk 파일 만들기
	// Field data
	outputFile << "    <FieldData>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		// values.push_back(var.fields[controls.fieldVar["time"].id]);
		values.push_back(stod(controls.startFrom));
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << 
	mesh.points.size()  << 
	"\" NumberOfCells=\"" << 
	mesh.cells.size()  << 
	// mesh.cells.size() + 1 << 
	"\">" << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// Points data
	outputFile << "    <PointData>" << endl;
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"pointLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.points.size());
		for(auto& point : mesh.points) values.push_back(point.level);
		// for(auto& point : ro_proc_points_x) values.push_back(-100);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </PointData>" << endl;
	
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	
	
	// {
		// outputFile << "     <DataArray type=\"UInt8\" Name=\"vtkGhostType\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// values.reserve(mesh.cells.size());
		// for(auto& cell : mesh.cells) values.push_back(0);
		// for(auto& cell : ro_proc_cell_ip_goints) values.push_back(1);
		// save.writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
		
	// }
	
	
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.cells.size());
		for(auto& cell : mesh.cells) values.push_back(cell.level);
		// for(auto& cell : ro_proc_cell_ip_goints) values.push_back(-100);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellGroups\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.cells.size());
		for(auto& cell : mesh.cells) values.push_back(cell.group);
		// for(auto& cell : ro_proc_cell_ip_goints) values.push_back(-100);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	
	
	
	
	// 스칼라 형식 데이터 저장
	for(int id=0; id<varScalarSize; ++id)
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
		primScalarNames[id] << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.reserve(mesh.cells.size());
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i) {
			auto cellVar_i = cellVar[i];
			values.push_back(cellVar_i[id]);
		}
	
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	// 벡터형식 데이터 저장
	for(int id=0; id<varVector3Size; ++id)
	{
		{
			outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
			primVector3Names[id] << "\" NumberOfComponents=\"" << 3 << 
			"\" format=\"" << saveFormat << "\">" << endl;
			vector<double> values;
			values.reserve(3*mesh.cells.size());
			vector<int> id_phi;
			auto cellVar = var.cells.data();
			int real_id0 = varScalarSize + 3*id;
			int real_id1 = varScalarSize + 3*id + 1;
			int real_id2 = varScalarSize + 3*id + 2;
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i) {
				auto cellVar_i = cellVar[i];
				values.push_back(cellVar_i[real_id0]);
				values.push_back(cellVar_i[real_id1]);
				values.push_back(cellVar_i[real_id2]);
			}
			save.writeDatasAtVTU(controls, outputFile, values);
			outputFile << "     </DataArray>" << endl;
			
		}
	}
	
	outputFile << "    </CellData>" << endl;
	
	
	
	
	
	
	// Points
	outputFile << "    <Points>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& point : mesh.points) {
			values.push_back(point.x);
			values.push_back(point.y);
			values.push_back(point.z);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "   </Points>" << endl;
	
	
	
	
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			for(auto i : cell.ipoints){
				values.push_back(i);
			}
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
			// if(rank==0) cout << endl;
	// offsets (cell's points offset)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFaceOffset = 0;
		for(auto& cell : mesh.cells){
			cellFaceOffset += cell.ipoints.size();
			values.push_back(cellFaceOffset);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
			// if(rank==0) cout << endl;
	
	// types (cell's type, 42 = polyhedron)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"types\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values(mesh.cells.size(),42);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faces\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			values.push_back(cell.ifaces.size());
			for(auto& i : cell.ifaces){
				values.push_back(mesh.faces[i].ipoints.size());
				for(auto& j : mesh.faces[i].ipoints){
					values.push_back(j);
				}
			}
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// if(rank==0) cout << endl;
	
	
	// faceoffsets (cell's face offset)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faceoffsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFacePointOffset = 0;
		for(auto& cell : mesh.cells){
			int numbering = 1 + cell.ifaces.size();
			for(auto& i : cell.ifaces){
				numbering += mesh.faces[i].ipoints.size();
			}
			cellFacePointOffset += numbering;
			values.push_back(cellFacePointOffset);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	// additional informations
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			values.push_back(face.iL);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << " </DataArray>" << endl;
	}
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				values.push_back(face.iR);
			}
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	
	
	
	// boundary informations
	{
		outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.name << endl;
			// trim;
			string bcName = boundary.name;
			
			bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// outputFile << boundary.name << " ";
			outputFile << bcName << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.startFace << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.nFaces << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.rightProcNo << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"connPoints\" format=\"" << "ascii" << "\">" << endl;
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [item, item2] : point.connPoints){
				outputFile << i << " " << item << " " << item2 << " ";
			}
			// if(!point.connPoints.empty()) outputFile << endl;
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		
	}
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}






void saveDPMFiles(string folder, int rank){
	
	// if(controls.nameParcels.size()==0) return;
	

	// int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// 폴더 만들기
    // auto ret = filesystem::create_directories(folder);
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	save.mkdirs(folder_name);
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save parcel files ... ";
	}
	
	ofstream outputFile;
	string filenamePlot = folder + "parcels." + to_string(rank) + ".vtu";
	
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	string saveFormat = controls.saveFormat;
	
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	if(controls.saveFormat == "ascii"){
		outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		}
		else{
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
	
	outputFile << "  <UnstructuredGrid>" << endl;

	outputFile << "    <FieldData>" << endl;
	outputFile << "    </FieldData>" << endl;

	outputFile << "   <Piece NumberOfPoints=\"" << mesh.parcels.size() << "\" NumberOfCells=\"" << 0 << "\">" << endl;
	
	outputFile << "    <PointData>" << endl;
	
	
	// 필수 데이터 저장
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"" <<
		"id" << "\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.parcels.size());
		for(auto& parcel : mesh.parcels){
			values.push_back(parcel.id);
		}
		if(values.size()!=0) save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"" <<
		"icell" << "\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.parcels.size());
		for(auto& parcel : mesh.parcels){
			values.push_back(parcel.icell);
		}
		if(values.size()!=0) save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	// 스칼라 형식 데이터 저장
	for(int id=0; id<parcel_primScalarNames.size(); ++id)
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
		parcel_primScalarNames[id] << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.reserve(mesh.parcels.size());
		auto parcelVar = var.parcels.data();
		for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
			auto parcelVar_i = parcelVar[i];
			values.push_back(parcelVar_i[id]);
		}
	
		if(values.size()!=0) save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	// 벡터형식 데이터 저장
	for(int id=0; id<parcel_primVector3Names.size(); ++id)
	{
		{
			outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
			parcel_primVector3Names[id] << "\" NumberOfComponents=\"" << 3 << 
			"\" format=\"" << saveFormat << "\">" << endl;
			vector<double> values;
			values.reserve(3*mesh.parcels.size());
			vector<int> id_phi;
			auto parcelVar = var.parcels.data();
			int real_id0 = parcel_primScalarNames.size() + 3*id;
			int real_id1 = parcel_primScalarNames.size() + 3*id + 1;
			int real_id2 = parcel_primScalarNames.size() + 3*id + 2;
			for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
				auto parcelVar_i = parcelVar[i];
				values.push_back(parcelVar_i[real_id0]);
				values.push_back(parcelVar_i[real_id1]);
				values.push_back(parcelVar_i[real_id2]);
			}
			if(values.size()!=0) save.writeDatasAtVTU(controls, outputFile, values);
			outputFile << "     </DataArray>" << endl;
			
		}
	}
	
	
	
	
	outputFile << "    </PointData>" << endl;

	
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	
	// Points
	outputFile << "    <Points>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		auto parcelVar = var.parcels.data();
		int str = parcel_primScalarNames.size()+3*parcel_primVector3Names.size();
		int id_x = str;
		int id_y = str+1;
		int id_z = str+2;
		for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
			auto parcelVar_i = parcelVar[i].data();
			values.push_back(parcelVar_i[id_x]);
			values.push_back(parcelVar_i[id_y]);
			values.push_back(parcelVar_i[id_z]);
		}
		if(values.size()!=0) save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "   </Points>" << endl;
	
	
	outputFile << "   <Cells>" << endl; 
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	outputFile << "</VTKFile>" << endl;
	outputFile.close();
	// MPI_Barrier(MPI_COMM_WORLD);
	
	
	


	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
}



//============================================


int main(int argc, char* argv[]) {
	
	MPI::Init(); 
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	map<string,string> mapArgv;
	for(int i=1; i<argc; i+=2){
		string first = argv[i]; string second;
		if(i+1==argc){ second = "nan"; }
		else{ second = argv[i+1]; }
		mapArgv.insert(make_pair(first, second));
	}
	
	if(
	(mapArgv.find("-help") != mapArgv.end() ||
	 mapArgv.find("-h") != mapArgv.end() ) ||
	(mapArgv.find("-from") == mapArgv.end() ||
	mapArgv.find("-to") == mapArgv.end())
	){
		if(rank==0) print_help();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	// set_mapArgv(argc, argv, mapArgv);
	int nSizeOrg = stoi(mapArgv["-from"]);
	int nSizeTar = stoi(mapArgv["-to"]);
	
	
	// 기본 variables 셋팅
	controls.setVariablesBasic();
	// 셋팅 파일 로드
	load.settingFiles("./setting/", controls);
	// 파일 로드
	loadVTUfiles(controls.getLoadFolderName());
	loadDPMFiles(controls.getLoadFolderName());
	
    // mesh.repartParMETIS(size, cell_ip, mesh);
	repartParMETIS(nSizeOrg, nSizeTar); 
	if(rank>=nSizeOrg) cell_ip_g.clear();
	{
		parcel_ip_g.resize(mesh.parcels.size());
		int iter=0;
		for(auto& parcel : mesh.parcels){
			parcel_ip_g[iter] = cell_ip_g[parcel.icell];
			++iter;
		}
	}
	
	vector<int> dummy;
	vector<int> to_new_cell_id(mesh.cells.size());
	
	int maxLevel = -1;
	for(auto& cell : mesh.cells){
		maxLevel = max(cell.level,maxLevel);
	}
	if(size>1){
		int nbc_glo;
		MPI_Allreduce(&maxLevel, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		maxLevel = nbc_glo;
	}
	
	mesh.repartitioning(cell_ip_g, maxLevel, to_new_cell_id, dummy);
    
	// MPI_Barrier(MPI_COMM_WORLD);
    // if(rank==0) cout << "AA1" << endl;
	
	
	resetVariableArray(to_new_cell_id);
    
	// MPI_Barrier(MPI_COMM_WORLD);
    // if(rank==0) cout << "AA2" << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	if(rank<nSizeTar) saveVTUFiles("./save/repartition/",rank);
	if(rank<nSizeTar) saveDPMFiles("./save/repartition/",rank);
	
	// if(rank==61) cout << mesh.cells.size() << " " << nSizeTar<< endl;
	
	// ==========================================
	// pvtu file
	if(rank==0){
		ofstream outputFile;
		string filenamePvtu = "./save/plot.";
		string stime = "./save/repartition/";
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		filenamePvtu += stime;
		filenamePvtu += ".pvtu";
		
		outputFile.open(filenamePvtu);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// string out_line;
		outputFile << "<?xml version=\"1.0\"?>" << endl;
		outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		outputFile << "  <PUnstructuredGrid>" << endl;

		outputFile << "   <PFieldData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"TimeValue\"/>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0; ip<nSizeTar; ++ip){
			string filenamevtus = "./" + stime;
			filenamevtus += "/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;

		outputFile << "    <PDataArray type=\"Int32\" Name=\"pointLevels\"/>" << endl;
		
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		// outputFile << "     <PDataArray type=\"UInt8\" Name=\"vtkGhostType\"/>" << endl;
		outputFile << "     <PDataArray type=\"Int32\" Name=\"cellLevels\"/>" << endl;
		outputFile << "     <PDataArray type=\"Int32\" Name=\"cellGroups\"/>" << endl;
		for(auto& name : primScalarNames){
			outputFile << "    <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << endl;
		}
		{
			int tmp_iter=0;
			for(auto& sup_name : primVector3Names){
				outputFile << "    <PDataArray type=\"Float64\" Name=\"" <<
				sup_name << "\" NumberOfComponents=\"" << 3 << "\"/>" << endl;
				// outputFile << "    <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << endl;
				++tmp_iter;
			}
		}
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	// ==========================================
	// pvtu file
	if(rank==0){
		ofstream outputFile;
		string filenamePvtu = "./save/parcels.";
		string stime = "./save/repartition/";
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		filenamePvtu += stime;
		filenamePvtu += ".pvtu";
		
		outputFile.open(filenamePvtu);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// string out_line;
		outputFile << "<?xml version=\"1.0\"?>" << endl;
		outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		outputFile << "  <PUnstructuredGrid>" << endl;

		outputFile << "   <PFieldData>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0; ip<nSizeTar; ++ip){
			string filenamevtus = "./" + stime;
			filenamevtus += "/parcels.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;

		outputFile << "    <PDataArray type=\"Int32\" Name=\"" << "id" << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int32\" Name=\"" << "icell" << "\"/>" << endl;
		for(auto& name : parcel_primScalarNames)
		{
			outputFile << "    <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << endl;
		}
		for(auto& sup_name : parcel_primVector3Names){
			outputFile << "    <PDataArray type=\"Float64\" Name=\"" <<
			sup_name << "\" NumberOfComponents=\"" << 3 << "\"/>" << endl;
		}
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	if(rank==0) cout << endl << "| COMPLETED Repartition ^_^ |" << endl << endl;
	
	MPI::Finalize();
	return EXIT_SUCCESS;
}




void print_help(){

	cout << endl;
	cout << "┌─────── Repartitioning helper ─────────────────────────────── " << endl;
	cout << "| # CAUTION # : mpi np >= max(from N, to N)" << endl;
	cout << "| -from \"int num.\"  : # of original mesh block" << endl;
	cout << "| -to   \"int num.\"  : # of target mesh block" << endl;
	cout << "└───────────────────────────────────────────────────────────────── " << endl;
	cout << endl;
}





