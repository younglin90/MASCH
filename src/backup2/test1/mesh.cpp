
#include "./mesh.h"

void MASCH_Mesh::check(){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	// MASCH_Log log;
	
	// 셀의 최대 id 찾기
	int maxL=0, maxR=0;
	int check_LR_Reverse=0;
	for(auto& face : mesh.faces){
		maxL = max(maxL , face.iL);
		maxR = max(maxR , face.iR);
	}
	int icell_max = max(maxL,maxR);
	
	// 바운더리 벡터 방향 제대로 해주기 = 포인트 순서 바꾸기
	for(auto& face : mesh.faces){
		if(face.iR > maxL){
			int tempnum = face.iL;
			face.iL = face.iR;
			face.iR = tempnum;
			std::reverse(face.ipoints.begin()+1,face.ipoints.end());
			++check_LR_Reverse;
		}
	}
	
	// 포인터 순서 바꾼 페이스들 mpi로 옮기고 경고창 띄우기
	if(size>1){
		if(rank == 0){
			vector<int> buffer(size,0);
			int gatherValue = check_LR_Reverse;
			MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
			int sumBuffer = 0;
			for(auto& i : buffer) {
				sumBuffer += i;
			}
			if(sumBuffer>0){
				std::stringstream ss;
				ss << "check face iL < iR, executed reverse : ";
				for(auto& i : buffer) {
					ss << i << " | ";
				}
				// log.warning.push(ss.str(),MASCH_FFL);
			}
		}
		else{
			int gatherValue = check_LR_Reverse;
			MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	else{		
		if(check_LR_Reverse>0){	
			std::stringstream ss;
			ss << "check face iL < iR, executed reverse : ";
			ss << check_LR_Reverse;
			// log.warning.push(ss.str(),MASCH_FFL);
		}
	}
	
	// 내부면인데 right id 가 없는 경우
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		// int id = static_cast<int>(i);
		
		// if(face.thereR == false){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			if(i < mesh.boundaries[0].startFace){
				std::stringstream ss;
				ss << "face i(" << i << ") < boundaries start face(" << 
				mesh.boundaries[0].startFace << ")";
				// log.warning.push(ss.str(),MASCH_FFL);
			}
		}
		
	}
	
	
}


// 페이스 타입 셋팅
void MASCH_Mesh::setFaceTypes(){
	
	auto& mesh = (*this);

	// set types
	for(auto& face : mesh.faces){
		face.setType(MASCH_Face_Types::INTERNAL);
	}
	
	for(int ibc=0; ibc<mesh.boundaries.size(); ++ibc){
		int startFace = mesh.boundaries[ibc].startFace;
		for(int i=startFace; i<startFace+mesh.boundaries[ibc].nFaces; ++i){
			// if( mesh.boundaries[ibc].thereR == false ){
			if(mesh.boundaries[ibc].getType() == MASCH_Face_Types::BOUNDARY){
				mesh.faces[i].setType(MASCH_Face_Types::BOUNDARY);
				// mesh.faces[i].setTypeBC(ibc);
			}
			else if(mesh.boundaries[ibc].getType() == MASCH_Face_Types::PROCESSOR){
				mesh.faces[i].setType(MASCH_Face_Types::PROCESSOR);
			}
			else{
				cout << "ERROR : 1654" << endl;
			}
		}
	}
	
}



// cell 생성
void MASCH_Mesh::buildCells(){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	// 원래 셀 삭제
	for(auto& cell : mesh.cells){
		cell.ipoints.clear();
		cell.ifaces.clear();
	}
	mesh.cells.clear();
	// 셀 전체 갯수 구하기
	int cell_num=0;
	for(auto& face : mesh.faces){
		cell_num = max(cell_num , face.iL);
		cell_num = max(cell_num , face.iR);
	}
	// add Cells
	if(mesh.faces.size() > 3) {
		for(int i=0; i<cell_num+1; ++i){
			mesh.addCell();
		}
	}
	
}




// 페이스 connection 
void MASCH_Mesh::connectFacetoPointsCells(){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	// 리셋
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		face.setPoints().clear();
	}
	
	// 페이스 포인트 및 셀 포인터 셋팅
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		face.setL(&mesh.cells[face.iL]);
		face.setR(nullptr);
		// if(face.thereR == true) face.setR(&mesh.cells[face.iR]);
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			face.setR(&mesh.cells[face.iR]);
		}
		for(auto& ipoint : face.ipoints){
			face.setPoints().push_back(&mesh.points[ipoint]);
		}
	}
}


// cell connection (cell's face)
void MASCH_Mesh::connectCelltoFaces(){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	// 리셋
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		cell.setFaces().clear();
		cell.ifaces.clear();
	}
	
	// 셀 페이스 넣기
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			mesh.cells[face.iL].ifaces.push_back(i);
			mesh.cells[face.iR].ifaces.push_back(i);
		}
		else if(
		face.getType() == MASCH_Face_Types::BOUNDARY ||
		face.getType() == MASCH_Face_Types::PROCESSOR ){
			mesh.cells[face.iL].ifaces.push_back(i);
		}
	}
	
	// 셀 페이스 포인터 셋팅
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		for(auto& iface : cell.ifaces){
			cell.setFaces().push_back(&mesh.faces[iface]);
		}
	}
}



// cell connection (cell's points)
void MASCH_Mesh::connectCelltoPoints(){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	// 리셋
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		cell.setPoints().clear();
		cell.ipoints.clear();
	}
	
	// 셀 포인트 셋팅
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			for(int i=0; i<face.ipoints.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face.iL].ipoints.size(); ++j){
					if(
					mesh.cells[face.iL].ipoints[j] == 
					face.ipoints[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face.iL].ipoints.push_back(face.ipoints[i]);
				
				l=0;
				for(int j=0; j<mesh.cells[face.iR].ipoints.size(); ++j){
					if(
					mesh.cells[face.iR].ipoints[j] == 
					face.ipoints[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face.iR].ipoints.push_back(face.ipoints[i]);
			}
			
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
			for(int i=0; i<face.ipoints.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face.iL].ipoints.size(); ++j){
					if(
					mesh.cells[face.iL].ipoints[j] == 
					face.ipoints[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face.iL].ipoints.push_back(face.ipoints[i]);
			}
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=0; i<face.ipoints.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face.iL].ipoints.size(); ++j){
					if(
					mesh.cells[face.iL].ipoints[j] == 
					face.ipoints[i] ) ++l;
				}
				if(l==0) mesh.cells[face.iL].ipoints.push_back(face.ipoints[i]);
			}
		}
	}
	
	// 셀 포인트 포인터 셋팅
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		for(auto& ipoint : cell.ipoints){
			cell.setPoints().push_back(&mesh.points[ipoint]);
		}
	}
}


// 프로세스 페이스 갯수 셋팅
void MASCH_Mesh::setCountsProcFaces(){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	mesh.countsProcFaces.clear();
	mesh.countsProcFaces.resize(size,0);
	for(int ip=0; ip<size; ++ip){
		mesh.countsProcFaces[ip] = 0;
		for(auto& bc : mesh.boundaries){
			// if(ip == bc.rightProcNo && bc.thereR==true){
			if(ip == bc.rightProcNo && bc.getType() == MASCH_Face_Types::PROCESSOR){
				mesh.countsProcFaces[ip] = bc.nFaces;
				break;
			}
		}
	}
}


// 프로세스 페이스 포인터 위치 셋팅
void MASCH_Mesh::setDisplsProcFaces(){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	mesh.displsProcFaces.clear();
	mesh.displsProcFaces.resize(size,0);
	mesh.displsProcFaces[0] = 0;
	for(int ip=1; ip<size; ++ip){
		mesh.displsProcFaces[ip] = mesh.displsProcFaces[ip-1] + mesh.countsProcFaces[ip-1];
		
	}
}



// 전체 셀 셋팅
void MASCH_Mesh::cellsGlobal(){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	int ncells = mesh.cells.size();
	vector<int> procNcells(size,0);
	MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	int myStart=0;
	for(int i=0; i<rank; ++i){
		myStart += procNcells[i];
	}
	
	int ncellTot;
	MPI_Allreduce(&ncells, &ncellTot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	vector<int> procStart(size,0);
	MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	mesh.startCellGlobal = myStart;
	mesh.startProcCellGlobal.resize(size+1,0);
	for(int i=0; i<size; ++i){
		mesh.startProcCellGlobal[i] = procStart[i];
	}
	mesh.ncellsTotal = ncellTot;
	mesh.startProcCellGlobal[size] = ncellTot;
	
	// vector<int> rightProcNo;
	// vector<int> ownerNo;
	// vector<int> neighbNo;
	
	// for(auto& boundaries : mesh.boundaries){
		// if(boundaries.thereR == true){
			// int str=boundaries.startFace;
			// int end=str+boundaries.nFaces;
			// for(int i=str; i<end; ++i){
				// rightProcNo.push_back(boundaries.rightProcNo);
				// ownerNo.push_back(mesh.faces[i].owner);
			// }
		// }
	// }
	
	// if(size>1){
		// neighbNo.clear();
		// neighbNo.resize(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		
		// MPI_Alltoallv( ownerNo.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_UNSIGNED, 
					   // neighbNo.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_UNSIGNED, 
					   // MPI_COMM_WORLD);
	
	// }
	
	// mesh.rightProcNo.resize(rightProcNo.size(),0);
	// for(int i=0; i<rightProcNo.size(); ++i){
		// mesh.rightProcNo[i] = rightProcNo[i];
	// }
	
	// mesh.procNeighbCellNo.resize(neighbNo.size(),0);
	// for(int i=0; i<neighbNo.size(); ++i){
		// mesh.procNeighbCellNo[i] = neighbNo[i];
	// }
	
	
	
	// // temporary COO format
	// int strRow = mesh.startCellGlobal;
	// vector<int> COO_row(mesh.cells.size(),0);
	// vector<int> COO_col(mesh.cells.size(),0);
	// vector<int> COO_iface;
	// vector<string> COO_iface_LR;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// COO_row[i] = strRow + i;
		// COO_col[i] = strRow + i;
	// }
	// for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// COO_row.push_back(strRow + face.owner);
			// COO_col.push_back(strRow + face.neighbour);
			
			// COO_iface_LR.push_back("LR");
			// COO_iface.push_back(i);
			
			// COO_row.push_back(strRow + face.neighbour);
			// COO_col.push_back(strRow + face.owner);
			
			// COO_iface_LR.push_back("RL");
			// COO_iface.push_back(i);
		// }
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// int strNeigbRow_Glo = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
			// int idNeigbCell_Loc = mesh.procNeighbCellNo[proc_num];
			
			// COO_row.push_back(strRow + face.owner);
			// COO_col.push_back(strNeigbRow_Glo + idNeigbCell_Loc);
			
			// COO_iface_LR.push_back("LR");
			// COO_iface.push_back(i);
			
			// ++proc_num;
		// }
	// }
	
	// // save CSR format
    // //compute number of non-zero entries per row of A 
	// int n_row = mesh.cells.size();
	// int nnz = COO_row.size();
	
	// mesh.non_zeros = nnz;
	
	// mesh.CRS_ptr.clear();
	// mesh.CRS_col.clear();
	// mesh.CRS_col_ptr_dig.clear();
	// mesh.CRS_col_ptr_LR.clear();
	// mesh.CRS_col_ptr_RL.clear();
	
	// mesh.CRS_col.resize(nnz,0);
	// mesh.CRS_ptr.resize(n_row+1,0);
	// mesh.CRS_col_ptr_dig.resize(mesh.cells.size(),-1);
	// mesh.CRS_col_ptr_LR.resize(mesh.faces.size(),-1);
	// mesh.CRS_col_ptr_RL.resize(mesh.faces.size(),-1);

    // for (int n = 0; n < nnz; n++){            
        // mesh.CRS_ptr[COO_row[n]-strRow]++;
    // }

    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = mesh.CRS_ptr[i];
        // mesh.CRS_ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // mesh.CRS_ptr[n_row] = nnz; 

    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = COO_row[n] - strRow;
        // int dest = mesh.CRS_ptr[row];

        // mesh.CRS_col[dest] = COO_col[n];
        // // val[dest] = A_vals[n];
		
		// if(n<n_row){
			// mesh.CRS_col_ptr_dig[n] = dest;
		// }
		// else{
			// int iface = COO_iface[n-n_row];
			
			// if(COO_iface_LR[n-n_row] == "LR"){
				// mesh.CRS_col_ptr_LR[iface] = dest;
			// }
			// else{
				// mesh.CRS_col_ptr_RL[iface] = dest;
			// }
			
		// }

        // mesh.CRS_ptr[row]++;
    // }

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = mesh.CRS_ptr[i];
        // mesh.CRS_ptr[i]  = last;
        // last = temp;
    // }
	
	
}



void MASCH_Mesh::informations(){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	nBoundFaces = 0;
	nProcFaces = 0;
	nBoundTypes = 0;
	for(auto& boundary : mesh.boundaries){
		// if(boundary.thereR == false){
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			nBoundFaces += boundary.nFaces;
			++nBoundTypes;
		}
		else{
			nProcFaces += boundary.nFaces;
		}
	}
	
	// faces type
	nInterFaces = 0;
	nTriangle = 0;
	nQuadrangle = 0;
	nPolygon = 0;
	for(auto& face : mesh.faces){
		// if(face.thereR == true){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			++nInterFaces;
		}
		
		if(face.ipoints.size() == 3){
			++nTriangle;
		}
		else if(face.ipoints.size() == 4){
			++nQuadrangle;
		}
		else{
			++nPolygon;
		}
	}
	
	
	// cells type
	nTetrahedron = 0;
	nHexahedron = 0;
	nPrism = 0;
	nPyramid = 0;
	nPolyhedron = 0;
	for(auto& cell : mesh.cells){
		if( cell.ipoints.size() == 4 && cell.ifaces.size() == 4 ){
			++nTetrahedron;
		}
		else if( cell.ipoints.size() == 5 && cell.ifaces.size() == 5 ){
			++nPyramid;
		}
		else if( cell.ipoints.size() == 6 && cell.ifaces.size() == 5 ){
			++nPrism;
		}
		else if( cell.ipoints.size() == 8 && cell.ifaces.size() == 6 ){
			++nHexahedron;
		}
		else {
			++nPolyhedron;
		}
	}
	
	
	
    if(rank == 0){
        vector<int> buffer(size,0);
		int gatherValue = mesh.points.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| present MPI size : " << size;
		cout << "| load data MPI size : " << size << endl;
		cout << "| points size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = mesh.faces.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = mesh.cells.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| cells size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nInterFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| internal faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nBoundFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| boundary faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nProcFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| processor faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "| boundary types : " << nBoundTypes << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nTriangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Triangle faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nQuadrangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Quadrangle faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPolygon;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Polygon faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nTetrahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Tetrahedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nHexahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Hexahedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPrism;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Prism cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPyramid;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Pyramid cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPolyhedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Polyhedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
    }
    else{
		int gatherValue = mesh.points.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = mesh.faces.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = mesh.cells.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nInterFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nBoundFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nProcFaces;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		gatherValue = nTriangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nQuadrangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPolygon;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		gatherValue = nTetrahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nHexahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPrism;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPyramid;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPolyhedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}



void MASCH_Mesh::setFaceLevels(){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	// face level
	vector<int> cLevel_recv;
	if(size>1){
		vector<int> cLevel_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				cLevel_send.push_back(mesh.cells[face.iL].level);
			}
		}
		
		cLevel_recv.clear();
		cLevel_recv.resize(cLevel_send.size(),0);

		MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	int proc_num=0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			face.level = max(
				mesh.cells[face.iL].level, mesh.cells[face.iR].level);
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			face.level = max(
				mesh.cells[face.iL].level, cLevel_recv[proc_num]);
			++proc_num;
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
			face.level = mesh.cells[face.iL].level;
		}
	}
}




void MASCH_TMP_MPI_Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<int> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	vector<int> tmp_recv_value(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_INT, 
				   tmp_recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
	recv_value.clear();
	recv_value.resize(size);
	for(int ip=0; ip<size; ++ip){
		int str = recv_displs[ip];
		int end = recv_displs[ip+1];
		for(int i=str; i<end; ++i){
			recv_value[ip].push_back(tmp_recv_value[i]);
		}
	}


}




void MASCH_Mesh::setCellStencils(){
	
	auto& mesh = *this;
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	
	vector<vector<int>> pointStencils(mesh.points.size()); 
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		
		cell.iStencils.clear();
		for(auto& ipoint : cell.ipoints){
			auto& point = mesh.points[ipoint];
			pointStencils[ipoint].push_back(i);
		}
	}
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		
		for(auto& ipoint : cell.ipoints){
			for(auto& icellSten : pointStencils[ipoint]){
				if(icellSten==i) continue;
				if(find(cell.iStencils.begin(),
				cell.iStencils.end(),
				icellSten)==cell.iStencils.end()){
					cell.iStencils.push_back(icellSten);
				}
			}
		}
	}
	
	
	vector<vector<int>> send_cells(size); 
	vector<vector<int>> send_point_id(size); 
	for(int i=0, iter=0, SIZE=mesh.points.size(); i<SIZE; ++i){
		auto& point = mesh.points[i];
		for(auto& [ip, id] : point.connPoints){
			for(auto& icell : pointStencils[i]){
				send_cells[ip].push_back(icell);
				send_point_id[ip].push_back(id);
			}
		}
	}
	
	vector<vector<int>> send_real_cells(size); 
	vector<vector<int>> send_real_cells_idLocal(size); 
	vector<vector<int>> send_real_points_id(size); 
	vector<vector<int>> send_n_cells(size,vector<int>(1,0));
	for(int ip=0; ip<size; ++ip){
		vector<int> cells_idLocal(mesh.cells.size(),-1);
		int tmp_size = send_cells[ip].size();
		for(int i=0; i<tmp_size; ++i){
			int id_cell = send_cells[ip][i];
			int id_point = send_point_id[ip][i];
			if(cells_idLocal[id_cell]==-1){
				int send_id = send_n_cells[ip][0]++;
				send_real_cells[ip].push_back(id_cell);
				cells_idLocal[id_cell] = send_id;
				send_real_cells_idLocal[ip].push_back(send_id);
				send_real_points_id[ip].push_back(id_point);
			}
			else{
				int send_id = cells_idLocal[id_cell];
				cells_idLocal[id_cell] = send_id;
				send_real_cells_idLocal[ip].push_back(send_id);
				send_real_points_id[ip].push_back(id_point);
			}
		}
	}
	
	vector<vector<int>> save_recv_cellStencils(mesh.cells.size()); 
	if(size>1){ 
		vector<vector<int>> recv_real_cells_idLocal(size); 
		vector<vector<int>> recv_real_points_id(size);
		vector<vector<int>> recv_n_cells(size);
		MASCH_TMP_MPI_Alltoallv(send_real_cells_idLocal, recv_real_cells_idLocal);
		MASCH_TMP_MPI_Alltoallv(send_real_points_id, recv_real_points_id);
		MASCH_TMP_MPI_Alltoallv(send_n_cells, recv_n_cells);
		vector<int> str_n_cells(size,0);
		for(int ip=0; ip<size; ++ip){
			str_n_cells[ip+1] = str_n_cells[ip] + recv_n_cells[ip][0];
		}
		
		for(int ip=0; ip<size; ++ip){
			int str = str_n_cells[ip];
			int iter=0;
			for(auto& ipoint : recv_real_points_id[ip]){
				int recv_icell = str + recv_real_cells_idLocal[ip][iter];
				for(auto& icell : pointStencils[ipoint]){
					if(find(save_recv_cellStencils[icell].begin(),
					save_recv_cellStencils[icell].end(),
					recv_icell)==save_recv_cellStencils[icell].end()){
						save_recv_cellStencils[icell].push_back(recv_icell);
					}
				}
				++iter;
			}
		}
		
		
	}
	
	
	send_StencilCellsId.clear();
	send_countsStencilCells.clear();
	send_displsStencilCells.clear();
	recv_countsStencilCells.clear();
	recv_displsStencilCells.clear();
	
	send_countsStencilCells.resize(size,0);
	for(int ip=0, iter=0; ip<size; ++ip){
		for(auto& icell : send_real_cells[ip]){
			send_StencilCellsId.push_back(icell);
			++iter;
		}
		send_countsStencilCells[ip] = iter;
	}
	send_displsStencilCells.resize(size+1,0);
	for(int ip=0, iter=0; ip<size; ++ip){
		send_displsStencilCells[ip+1] = send_displsStencilCells[ip] + send_countsStencilCells[ip];
	}
	if(size>1){ 
		recv_countsStencilCells.resize(size,0);
		recv_displsStencilCells.resize(size+1,0);
				
		MPI_Alltoall(
		send_countsStencilCells.data(), 1, MPI_INT, 
		recv_countsStencilCells.data(), 1, MPI_INT, MPI_COMM_WORLD);
		
		for(int ip=0; ip<size; ++ip) recv_displsStencilCells[ip+1] = recv_displsStencilCells[ip] + recv_countsStencilCells[ip];
	}
	
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		mesh.cells[i].recv_iStencils.clear();
		for(auto& icell : save_recv_cellStencils[i]){
			mesh.cells[i].recv_iStencils.push_back(icell);
		}
	}
	
	
}


void MASCH_Mesh::setNumberOfFaces(){
	
	auto& mesh = *this;
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	nInternalFaces = 0;
	nBoundaryFaces = 0;
	nProcessorFaces = 0;
	for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		if(mesh.faces[i].getType()==MASCH_Face_Types::INTERNAL){
			++nInternalFaces;
		}
		else if(mesh.faces[i].getType()==MASCH_Face_Types::BOUNDARY){
			++nBoundaryFaces;
		}
		else if(mesh.faces[i].getType()==MASCH_Face_Types::PROCESSOR){
			++nProcessorFaces;
		}
	}
}


void MASCH_Mesh::set(
vector<double>& NodeCoordinates, vector<int>& connectivity, vector<int>& offsets, 
vector<int>& inp_faces, vector<int>& faceoffsets,
vector<int>& inp_owner, vector<int>& inp_neighbour, vector<string>& bcName, vector<int>& bcStartFace, 
vector<int>& bcNFaces, vector<int>& bcNeighbProcNo, vector<int>& connPoints,
vector<int>& pointLevels, vector<int>& cellLevels, vector<int>& cellGroups){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
	
	
	{
		int n=0;
		for(auto& i : inp_neighbour){
			mesh.faces[n].iR = i;
			mesh.faces[n].setType(MASCH_Face_Types::INTERNAL);
			++n;
		}
		inp_neighbour.clear();
	}
	

	
	{
		int m=0;
		int n=0;
		for(auto& i : offsets){
			for(int j=n; j<i; ++j){
				int inp_point = connectivity[j];
				// cout << m << endl;
				mesh.cells[m].ipoints.push_back( static_cast<int>(inp_point) );
				// if(rank==1) cout << point << endl;
			}
			
			// if(rank==0 && mesh.cells[m].ipoints.size()!=8) cout << mesh.cells[m].ipoints.size() << endl;
			
			n=i;
			++m;
		}
	}
	
	{
		int n=0;
		int nFacesInt=0;
		for(auto& face : mesh.faces){
			// if(face.iR != -1){
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
				mesh.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
				++nFacesInt;
			}
			else{
				mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			}
			++n;
		}
	}
	
	{
		int m=0;
		int n=0;
		for(auto& i : faceoffsets){
			// if(faces[n]>5) cout << faces[n] << endl;
			int N=0;
			int face_size = inp_faces[m+N];
			for(int j=0; j<face_size; ++j){
				int inp_face = mesh.cells[n].ifaces[j];
				++N;
				int point_size = inp_faces[m+N];
				for(int k=0; k<point_size; ++k){
					++N;
					int inp_point = inp_faces[m+N];
					if(mesh.faces[ inp_face ].ipoints.size() == point_size) continue;
					mesh.faces[ inp_face ].ipoints.push_back( static_cast<int>(inp_point) );
					// if(rank==1) cout << point << endl;
				}
			}
			m=i;
			++n;
		}
		inp_faces.clear();
		faceoffsets.clear();
	}
	
	
	{
		int n=0;
		for(auto& startFace : bcStartFace){
			
			mesh.addBoundary();
			
			string tmp_bcName = bcName[n];
			tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
			tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
			tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
			tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
				
			mesh.boundaries.back().name = tmp_bcName;
			mesh.boundaries.back().startFace = bcStartFace[n];
			mesh.boundaries.back().nFaces = bcNFaces[n];
			if(bcNeighbProcNo[n]<0){
				mesh.boundaries.back().rightProcNo = 0;
				mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			}
			else{
				mesh.boundaries.back().rightProcNo = bcNeighbProcNo[n];
				mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
			}
			mesh.boundaries.back().myProcNo = rank;
			// if(bcNeighbProcNo[n] < 0){
				// mesh.boundaries.back().myProcNo = -1;
			// }
			
			++n;
			
		}
	}
	
	{
		int maxBCNum = mesh.boundaries.size()-1;
		mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
		for(int i=maxBCNum-1; i>=0; --i){
			mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
		}
		
		bcName.clear();
		bcStartFace.clear();
		bcNFaces.clear();
		bcNeighbProcNo.clear();
	}
	
		
	for(int i=0; i<connPoints.size(); ++i){
		// cout << mesh.points.size() << " " << connPoints[i] << " " << connPoints[i+1] << " " << connPoints[i+2] << " " << endl;
		mesh.points[connPoints[i]].connPoints.push_back(
		make_pair(connPoints[i+1],connPoints[i+2])
		);
		++i; ++i;
	}
	
	
	
	
	
	
	
	
	
	// for(int i=0; i<NodeCoordinates.size()/3; ++i){
		// mesh.addPoint();
		// mesh.points.back().x = NodeCoordinates[i*3+0];
		// mesh.points.back().y = NodeCoordinates[i*3+1];
		// mesh.points.back().z = NodeCoordinates[i*3+2];
	// }
	// NodeCoordinates.clear();
	
	
	// int n=0;
	// int ncells=0;
	// mesh.faces.clear();
	// for(auto& i : owner){
		// mesh.addFace();
		// mesh.faces.back().iL = i;
		// ncells = max(ncells, mesh.faces.back().iL);
	// }
	// owner.clear();
	
	
	// for(int i=0; i<pointLevels.size(); ++i){
		// mesh.points[i].level = pointLevels[i];
	// }
	
	// mesh.cells.clear();
	// for(int i=0; i<ncells+1; ++i){
		// mesh.addCell();
	// }
	// for(int i=0; i<cellLevels.size(); ++i){
		// mesh.cells[i].level = cellLevels[i];
	// }
	// for(int i=0; i<cellGroups.size(); ++i){
		// mesh.cells[i].group = cellGroups[i];
	// }
	
	
	// n=0;
	// for(auto& i : neighbour){
		// mesh.faces[n].iR = i;
		// mesh.faces[n].setType(MASCH_Face_Types::INTERNAL);
		// ++n;
	// }
	// neighbour.clear();
	
	
	
	// int m=0;
	// n=0;
	// for(auto& i : offsets){
		// for(int j=n; j<i; ++j){
			// int point = connectivity[j];
			// // cout << m << endl;
			// mesh.cells[m].ipoints.push_back( static_cast<int>(point) );
			// // if(rank==1) cout << point << endl;
		// }
		
		// // if(rank==0 && mesh.cells[m].ipoints.size()!=8) cout << mesh.cells[m].ipoints.size() << endl;
		
		// n=i;
		// ++m;
	// }
	
	
	
	// n=0;
	// int nFacesInt=0;
	// for(auto& face : mesh.faces){
		// // if(face.iR != -1){
		// if(face.getType() == MASCH_Face_Types::INTERNAL){
			// mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			// mesh.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
			// ++nFacesInt;
		// }
		// else{
			// mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
		// }
		// ++n;
	// }
	
	
	// m=0;
	// n=0;
	// for(auto& i : faceoffsets){
		// // if(faces[n]>5) cout << faces[n] << endl;
		// int N=0;
		// int face_size = faces[m+N];
		// for(int j=0; j<face_size; ++j){
			// int face = mesh.cells[n].ifaces[j];
			// ++N;
			// int point_size = faces[m+N];
			// for(int k=0; k<point_size; ++k){
				// ++N;
				// int point = faces[m+N];
				// if(mesh.faces[ face ].ipoints.size() == point_size) continue;
				// mesh.faces[ face ].ipoints.push_back( static_cast<int>(point) );
				// // if(rank==1) cout << point << endl;
			// }
		// }
		// m=i;
		// ++n;
	// }
	// faces.clear();
	// faceoffsets.clear();
	
	
	
	// n=0;
	// for(auto& startFace : bcStartFace){
		
		
		// mesh.addBoundary();
		
		// string tmp_bcName = bcName[n];
		// tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
		// tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
		// tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
		// tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
		// mesh.boundaries.back().name = tmp_bcName;
		// mesh.boundaries.back().startFace = bcStartFace[n];
		// mesh.boundaries.back().nFaces = bcNFaces[n];
		// if(bcNeighbProcNo[n]<0){
			// mesh.boundaries.back().rightProcNo = 0;
			// mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
		// }
		// else{
			// mesh.boundaries.back().rightProcNo = bcNeighbProcNo[n];
			// mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
		// }
		// mesh.boundaries.back().myProcNo = rank;
		// // if(bcNeighbProcNo[n] < 0){
			// // mesh.boundaries.back().myProcNo = -1;
		// // }
		
		// ++n;
		
	// }
	
	// // cout << bcStartFace.size() << endl;
	// int maxBCNum = mesh.boundaries.size()-1;
	// // cout << maxBCNum << endl;
	// mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
	// for(int i=maxBCNum-1; i>=0; --i){
		// mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
	// }
	
	// bcName.clear();
	// bcStartFace.clear();
	// bcNFaces.clear();
	// bcNeighbProcNo.clear();
	
		
	// for(int i=0; i<connPoints.size(); ++i){
		// // cout << mesh.points.size() << " " << connPoints[i] << " " << connPoints[i+1] << " " << connPoints[i+2] << " " << endl;
		// mesh.points[connPoints[i]].connPoints.push_back(
		// make_pair(connPoints[i+1],connPoints[i+2])
		// );
		// ++i; ++i;
	// }
	
	
	
	
}





// void SEMO_Mesh_Builder::setFaceLevels(vector<double>& output){
	
	// // face level
	// vector<int> cLevel_recv;
	// if(size>1){
		// vector<int> cLevel_send;
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// cLevel_send.push_back(mesh.cells[face.owner].level);
			// }
		// }
		
		// cLevel_recv.clear();
		// cLevel_recv.resize(cLevel_send.size(),0);

		// MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // MPI_COMM_WORLD);
	// }
	
	// int proc_num=0;
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// face.level = max(
				// mesh.cells[face.owner].level, mesh.cells[face.neighbour].level);
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// face.level = max(
				// mesh.cells[face.owner].level, cLevel_recv[proc_num]);
			// ++proc_num;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// face.level = mesh.cells[face.owner].level;
		// }
	// }
	


		
// }


// // void SEMO_Mesh_Builder::calcSkewness(vector<double>& output){
	
	// // SEMO_Mesh_Builder& mesh = *this;
	
	// // output.clear();
	// // output.resize(mesh.cells.size(),0.0);

	// // for(auto& face : mesh.faces){
		// // double value1=0.0;
		// // double value2=0.0;
		// // value1 += face.vecSkewness[0]*face.vecSkewness[0];
		// // value1 += face.vecSkewness[1]*face.vecSkewness[1];
		// // value1 += face.vecSkewness[2]*face.vecSkewness[2];
		// // value2 += face.vecPN[0]*face.vecPN[0];
		// // value2 += face.vecPN[1]*face.vecPN[1];
		// // value2 += face.vecPN[2]*face.vecPN[2];
		// // double value = sqrt(value1)/sqrt(value2);
		// // output[face.owner] = max(output[face.owner],value);
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // output[face.neighbour] = max(output[face.neighbour],value);
		// // }
	// // }
		
// // }

// // void SEMO_Mesh_Builder::calcNonOrthogonality(vector<double>& output){
	
	// // SEMO_Mesh_Builder& mesh = *this;
	
	// // output.clear();
	// // output.resize(mesh.cells.size(),0.0);
	// // for(auto& face : mesh.faces){
		// // double value1=0.0;
		// // double value2=0.0;
		// // value1 += face.vecPN[0]*face.unitNormals[0];
		// // value1 += face.vecPN[1]*face.unitNormals[1];
		// // value1 += face.vecPN[2]*face.unitNormals[2];
		// // value2 += face.vecPN[0]*face.vecPN[0];
		// // value2 += face.vecPN[1]*face.vecPN[1];
		// // value2 += face.vecPN[2]*face.vecPN[2];
		// // double value = acos(abs(value1)/sqrt(value2));
		// // // cout << face.owner << endl;
		// // output[face.owner] = max(output[face.owner],value);
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // output[face.neighbour] = max(output[face.neighbour],value);
		// // }
	// // }
	
// // }




// // void SEMO_Mesh_Builder::calcUniformity(vector<double>& output){
	
	// // SEMO_Mesh_Builder& mesh = *this;
	
	// // output.clear();
	// // output.resize(mesh.cells.size(),1.e8);
	// // for(auto& face : mesh.faces){
		// // double value1=0.0;
		// // double value2=0.0;
		// // for(int ii=0; ii<3; ++ii){
			// // value1 += face.vecPF[ii]*face.vecPF[ii];
			// // value2 += face.vecNF[ii]*face.vecNF[ii];
		// // }
		// // value1 = sqrt(value1);
		// // value2 = sqrt(value2);
		// // value2 = value1 + value2;
		// // double value = value1/value2;
		// // value = 1.0-2.0*abs(value-0.5);
		// // output[face.owner] = min(output[face.owner],value);
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // output[face.neighbour] = min(output[face.neighbour],value);
		// // }
	// // }
		
// // }