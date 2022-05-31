
#include "./mesh.h"
#include "./polyAMR.h"
// #include "geometric.h"
#include "./mpi.h"


void MASCH_Poly_AMR_Builder::polyRefine(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	int maxLevel_AMR, int maxCells_AMR, double minVolume_AMR, 
	vector<vector<double>> indicatorCriterion,
	vector<vector<double>>& indicatorValues,
	vector<vector<int>>& child_new_cell_id_of_org,
	vector<bool>& boolCellPreserved,
	vector<bool>& boolCellRefine,
	vector<bool>& boolCellUnrefine,
	int iter){
		
		

	// int nBuffers = 3;
		
		
		
		
	// int maxLevel_AMR = controls.maxLevelRefine;
	// double minVolume_AMR = controls.minVolumeRefine;
	// int maxCells_AMR = controls.maxCellsRefine;
	// // ===========================================
	// // connPoints 디버깅
	// mesh.debug_connPoints(1.e-36);
	// // proc face points 디버깅
	// mesh.debug_procFacePoints(1.e-36);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // ===========================================

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_MPI mpi;
	
	int beforeCellSize = mesh.cells.size();
	int beforeFaceSize = mesh.faces.size();
	int beforePointSize = mesh.points.size();
	int afterCellSize = 0;
	int afterFaceSize = 0;
	int afterPointSize = 0;
	
	
	child_new_cell_id_of_org.clear();
	child_new_cell_id_of_org.resize(beforeCellSize);
	
	

	// //=======================================
	// // 디버그
	// {
		
		// vector<vector<int>> tmp_procFace_id(size);
		// vector<vector<int>> tmp_procFace_group(size);
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			// int rightProcNo = boundary.rightProcNo;
			// for(int i=str, iter=0; i<end; ++i){
				// auto& face = mesh.faces[i];
				// int iL = face.iL;
				// // int igroup = mesh.cells[iL].group-min_group_id;
				// int igroup = mesh.cells[iL].group;
				// tmp_procFace_id[rightProcNo].push_back(i);
				// tmp_procFace_group[rightProcNo].push_back(igroup);
			// }
		// }
		
		
		// if(rank==0){
			// int iter=0;
			// for(auto& item : tmp_procFace_id[1]){
				// cout << rank << ", " << item << ", " << tmp_procFace_group[1][iter] << endl;
				// ++iter;
			// }
		// }
		
		// MPI_Barrier(MPI_COMM_WORLD);
		
		// if(rank==1){
			// int iter=0;
			// for(auto& item : tmp_procFace_id[0]){
				// cout << rank << ", " << item << ", " << tmp_procFace_group[0][iter] << endl;
				// ++iter;
			// }
		// }
	// }
			
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "START0" << endl;
	
    bool boolDebug = false;
	
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute AMR - Refine ";
	}
	
	// SEMO_Mesh_Geometric geometric;
	
	// int nTotalFaceLRVar = mesh.faces[0].varL.size();
	// int nTotalFaceVar = mesh.faces[0].var.size();
	
	//====================================================
	// Refine 되는 셀 & 면 조사
	
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	
	MASCH_Poly_Mesh_Refine refineMesh;
	
	// int maxLevel_AMR = controls.maxLevelRefine;
	// double minVolume_AMR = controls.minVolumeRefine;
	// int maxCells_AMR = controls.maxCellsRefine;
	
	// refine 되는 셀 조사
	// vector<bool> boolCellRefine(mesh.cells.size(),false);
	// cout << maxCells_AMR << endl;
	// if(mesh.cells.size()<maxCells_AMR){
		bool boolMaxCells = false;
		if(mesh.cells.size()>maxCells_AMR) boolMaxCells = true;
		// cout << maxCells_AMR << endl;
		// cout << mesh.cells.size() << " " << minVolume_AMR << endl;
		for(int i=0; i<mesh.cells.size(); ++i){
            auto& cell = mesh.cells[i];
            // boolCellRefine[i] = false;
			// for(int indi=0; indi<indicatorCriterion.size(); ++indi)
			// {
				// for(int level=0; level<indicatorCriterion.at(indi).size(); ++level)
				// {
					// double indicatorRefine_AMR = indicatorCriterion.at(indi).at(level);
					// if( mesh.cells.at(i).level < level+1 ){
						// if( indicatorValues.at(indi).at(i) > indicatorRefine_AMR ){
							// boolCellRefine[i] = true;
						// }
					// }				
				// }
			// }
							// // boolCellRefine[i] = true;
			
				// if( 
				// (rank==0 && distr(eng) > 0.8) ||
				// (rank==1 && distr(eng) > 100.0) ||
				// (rank==2 && distr(eng) > 100.0) ||
				// (rank==3 && distr(eng) > 100.0) 
				// ){
					// cout << mesh.cells[i].volume << endl;
					// boolCellRefine[i] = true;
				// }

					// boolCellRefine[i] = false;
				// if(distr(eng) > 0.8){
					// boolCellRefine[i] = true;
                // }                
			
			// if(mesh.cells[i].volume < minVolume_AMR) {
				// cout << mesh.cells[i].volume << " " << minVolume_AMR << endl;
				// boolCellRefine[i] = false;
			// }
            
            // boolCellRefine[i] = true;
            
			if(mesh.cells[i].level >= maxLevel_AMR) boolCellRefine[i] = false;
			if(mesh.cells[i].level < 0) boolCellRefine[i] = false;
			if(boolMaxCells==true) boolCellRefine[i] = false;
            
            
            // if(cell.ipoints.size()==5 && cell.ifaces.size()==5 && boolCellRefine[i]==true) cout << "AAAAAAAAA" << endl;
			
			
			// if(boolCellPreserved[i] == true) boolCellRefine[i] = false;
		} 
        
        
		// for(int i=0; i<mesh.nInternalFaces; ++i){
            // auto& face = mesh.faces[i];
            // if(mesh.cells[face.iL].level < 0) boolCellRefine[face.iR] = false;
            // if(mesh.cells[face.iR].level < 0) boolCellRefine[face.iL] = false; 
		// } 
        
        
		// for(int i=0; i<mesh.faces.size(); ++i){
            // auto& face = mesh.faces[i];
            // if(face.getType()==MASCH_Face_Types::PROCESSOR){
                // boolCellRefine[face.iL] = false;
                
            // }
            // // if(mesh.cells[face.iL].level < 0) boolCellRefine[face.iR] = false;
            // // if(mesh.cells[face.iR].level < 0) boolCellRefine[face.iL] = false; 
		// } 
        
	// // processor faces
    // {
        // vector<int> recv_value2;
        // if(size>1){

            // vector<int> send_value2;
            // send_value2.reserve(mesh.send_StencilCellsId.size());
            // for(auto& icell : mesh.send_StencilCellsId){
                // if(mesh.cells[icell].level<0){
                    // send_value2.push_back(1);
                // }
                // else{
                    // send_value2.push_back(0);
                // }
            // }
            // recv_value2.resize(mesh.recv_displsStencilCells[size]);
            // MPI_Alltoallv( send_value2.data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_INT, 
                           // recv_value2.data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_INT, 
                           // MPI_COMM_WORLD);
            
        // }
        // for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
            // auto& cell = mesh.cells[i];
            
            // bool tmp_bool = false;
            // if(cell.level<0) tmp_bool = true;
            // for(auto& icell : cell.iStencils){
                // auto& cellSten = mesh.cells[icell];
                // if(cellSten.level<0) tmp_bool = true;
            // }
            
        // // controls.log.pop();
            // for(auto& icell : cell.recv_iStencils){
                // if(recv_value2[icell]==1) tmp_bool = true;
            // }
            
            // if(tmp_bool==true) boolCellRefine[i] = false;
        // }
    // }
    
		
	// }
	
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// boolCellPreserved[i] = true;
		// }
	
	// this->bufferLayerRefine(mesh, boolCellRefine, maxLevel_AMR, nBuffers);
	
	
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START1" << endl;
	
	
	// amr 정보 mpi 교환
	vector<int> cLevel_recv;
	vector<int> cRefine_recv;
	this->mpiLevelRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	this->restrictCellRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	this->mpiRefines(mesh, boolCellRefine, cRefine_recv);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0 && iter==1) cout << "REGION1" << endl;
	
	// 엣지 만들기
	vector<int> edgesPoint0;
	vector<int> edgesPoint1;
	vector<vector<int>> facesEdges;
	vector<vector<int>> edgesFaces;
	vector<int> edgesLevel;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0 && iter==1) cout << "REGION2" << endl;
	
	this->createEdges(mesh, edgesPoint0, edgesPoint1, facesEdges, edgesFaces, edgesLevel);
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0 && iter==1) cout << "REGION3" << endl;
	
	// ==========================================================
	// proc 엣지 만들기 (proc, id)
	vector<vector<pair<int,int>>> procEdges(edgesPoint0.size());
	// if(iter!=1)
	{
		vector<vector<int>> send_procEdges(size);
		vector<vector<int>> send_procEdgesPoint0(size);
		vector<vector<int>> send_procEdgesPoint1(size);
		vector<vector<int>> pointsEdges(mesh.points.size());
		for(int i=0, SIZE=edgesPoint0.size(); i<SIZE; ++i){
			int ipoint0 = edgesPoint0[i];
			int ipoint1 = edgesPoint1[i];
			for(auto& [proc0, id0] : mesh.points[ipoint0].connPoints){
				for(auto& [proc1, id1] : mesh.points[ipoint1].connPoints){
					if(proc0==proc1){
						send_procEdges[proc0].push_back(i);
						send_procEdgesPoint0[proc0].push_back(id0);
						send_procEdgesPoint1[proc1].push_back(id1);
					}
				}
			}
			if(
			find(pointsEdges[ipoint0].begin(),
			pointsEdges[ipoint0].end(),
			i) == pointsEdges[ipoint0].end()){
				pointsEdges[ipoint0].push_back(i);
			}
			if(
			find(pointsEdges[ipoint1].begin(),
			pointsEdges[ipoint1].end(),
			i) == pointsEdges[ipoint1].end()){
				pointsEdges[ipoint1].push_back(i);
			}
		}
		
		vector<vector<int>> recv_procEdges;
		mpi.Alltoallv(send_procEdges, recv_procEdges);
		vector<vector<int>> recv_procEdgesPoint0;
		mpi.Alltoallv(send_procEdgesPoint0, recv_procEdgesPoint0);
		vector<vector<int>> recv_procEdgesPoint1;
		mpi.Alltoallv(send_procEdgesPoint1, recv_procEdgesPoint1);
		
		for(int ip=0; ip<size; ++ip){
			int tmp_size = recv_procEdges[ip].size();
			for(int i=0; i<tmp_size; ++i){
				int recv_iedge = recv_procEdges[ip][i];
				int ipoint0 = recv_procEdgesPoint0[ip][i];
				int ipoint1 = recv_procEdgesPoint1[ip][i];
				
				for(auto& iedge : pointsEdges.at(ipoint0)){
					if(
					(edgesPoint0.at(iedge)==ipoint0 && 
					 edgesPoint1.at(iedge)==ipoint1) ||
					(edgesPoint0.at(iedge)==ipoint1 && 
					 edgesPoint1.at(iedge)==ipoint0)
					){
						procEdges.at(iedge).push_back(make_pair(ip,recv_iedge));
						break;
					}
				}
			}
		}
	}
	// ==========================================================
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START2" << endl;
	
	
	
	// 엣지 리파인 
	vector<bool> boolEdgeRefine(edgesPoint0.size(),false);
	
	// 그룹 차일드 페이스들 만들기
	vector<int> groupChildFaces_id(mesh.faces.size(),-1);
	vector<bool> groupChildFaces_HighLevel(mesh.faces.size(),false);
	vector<groupMesh_Refine> groupChildFaces;
	//====
	vector<int> proc_edgeCenterPoints_id(edgesPoint0.size(),-1);
	//====
	for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];
		int faceLevel = face.level;
		
		bool ownRefine = boolCellRefine[face.iL];
		bool ngbRefine;
		int ownLevel = mesh.cells[face.iL].level;
		int ngbLevel;
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			ngbRefine = boolCellRefine[face.iR];
			ngbLevel = mesh.cells[face.iR].level;
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			ngbRefine = false;
			if(cRefine_recv[proc_num]==1) ngbRefine = true;
			ngbLevel = cLevel_recv[proc_num];
			++proc_num;
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
			ngbRefine = false;
			ngbLevel = faceLevel;
		}
			
			
		if(
		(ownRefine==true && ownLevel==faceLevel) ||
		(ngbRefine==true && ngbLevel==faceLevel) 
		){
			
			groupChildFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			
			groupChildFace.parent = i;
			groupChildFace.type = 0;
			
			vector<int> tmpFacePoints;
			vector<int> tmpFacePointLevels;
			for(auto& j : face.ipoints){
				// cout << i << " " << j << endl;
				tmpFacePoints.push_back(j);
				tmpFacePointLevels.push_back(mesh.points[j].level);
			}
			
			extractVertexPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				groupChildFace.vertexPoints);

			vector<bool> addVertexCeterPoints;
			extractVertexCenterPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints,
				addVertexCeterPoints);
				
			vector<int> canRefineEdgeOrders;
			for(int j=0; j<facesEdges[i].size(); ++j){
				int iEdge = facesEdges[i][j];
				if( 
				edgesLevel[iEdge]==faceLevel && 
				boolEdgeRefine[iEdge]==false
				){
					canRefineEdgeOrders.push_back(iEdge);
				}
			}
			
			if( groupChildFace.vertexPoints.size() != groupChildFace.vertexCenterPoints.size() ){
				cout << "| WARNING : vertexPoints.size != vertexCenterPoints.size " << groupChildFace.vertexPoints.size() << " " << groupChildFace.vertexCenterPoints.size() << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			// //==============
			// // if(iter!=1)
			// vector<int> test0;
			// {
				// int tmp_iter = 0;
				// int tmp_iter2 = 0;
				// for(auto& vCP : groupChildFace.vertexCenterPoints){
					// if(vCP == -1){
						// proc_edgeCenterPoints_id[facesEdges[i].at(tmp_iter)] =
						// mesh.points.size()+tmp_iter2;
						// // proc_point_id_iedge.push_back(make_pair(
						// // mesh.points.size()+tmp_iter2,facesEdges[tmp_iter]));
						// test0.push_back(tmp_iter);
						// ++tmp_iter2;
					// }
					// ++tmp_iter;
				// }
			// }
			// //===========
				
			addVertexCenterPoint(mesh, faceLevel, 
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints);
				
			addEdgeFacesPointsToVertexCenterPoint(
				mesh,
				i,
				facesEdges,
				edgesFaces,
				boolEdgeRefine,
				addVertexCeterPoints,
				canRefineEdgeOrders,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints);
				
			//==============
			{
				int vertexSize = groupChildFace.vertexPoints.size();
				int tmpNewNum = 0;
				for(int j=0; j<vertexSize; ++j){
					if(addVertexCeterPoints[j]==true){
						int iEdge = canRefineEdgeOrders[tmpNewNum];
						++tmpNewNum;
						int isertP = groupChildFace.vertexCenterPoints[j];
						
						proc_edgeCenterPoints_id[iEdge] = isertP;
                    // if(rank==419 && iEdge==9321) cout << "AAAAAAAAAA" << endl;
					}
				}
			}
			//==============
				
			// for(auto& tes : test0){
				// if(boolEdgeRefine.at(facesEdges[i].at(tes))==false){
					// cout << "#ERROR 23" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
			// }
			
		}
		else if(
		(ownRefine==true && ownLevel+1==faceLevel) ||
		(ngbRefine==true && ngbLevel+1==faceLevel) 
		){
			groupChildFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces_HighLevel[i] = true;
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			
			groupChildFace.parent = i;
			groupChildFace.type = 0;
			
			vector<int> tmpFacePoints;
			vector<int> tmpFacePointLevels;
			for(auto& j : face.ipoints){
				tmpFacePoints.push_back(j);
				tmpFacePointLevels.push_back(mesh.points[j].level);
			}
			
			vector<int> vertexPoints;
			extractVertexPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				vertexPoints);
			
			groupChildFace.vertexPoints.push_back(vertexPoints[0]);
			groupChildFace.vertexPoints.push_back(vertexPoints[1]);
			groupChildFace.vertexPoints.push_back(vertexPoints[2]);
			groupChildFace.vertexPoints.push_back(vertexPoints[3]);
			
			groupChildFace.vertexCenterPoints.push_back(vertexPoints[1]);
			groupChildFace.centerPoint = vertexPoints[2];
			groupChildFace.vertexCenterPoints.push_back(vertexPoints[3]);
			
			
		}
	}
	
	
	
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START3" << endl;
	
	
	
	
	
	
	
	// ==========================================================
	// proc 엣지 판단해서 refine 하기
	// proc 엣지가 하나라도 refine 이면 모두 refine
	// vector<faces_proc_refind> addiFaces(mesh.faces.size());
	// if(iter!=0)
	{
		vector<vector<int>> send_edgeId(size);
		vector<vector<int>> send_boolEdgeRefine(size);
		int tmp_iter = 0;
		for(auto& procEdge : procEdges){
			for(auto& [proc, id] : procEdge){
				send_edgeId[proc].push_back(id);
				if(boolEdgeRefine[tmp_iter]==true){
					send_boolEdgeRefine[proc].push_back(1);
				}
				else{
					send_boolEdgeRefine[proc].push_back(0);
				}
			}
			++tmp_iter;
		}
		
		vector<vector<int>> recv_edgeId;
		mpi.Alltoallv(send_edgeId, recv_edgeId);
		vector<vector<int>> recv_boolEdgeRefine;
		mpi.Alltoallv(send_boolEdgeRefine, recv_boolEdgeRefine);
		
		
		vector<bool> bool_addiFaces(mesh.faces.size(),false);
		vector<bool> bool_edges(edgesPoint0.size(),false);
		vector<vector<int>> send_edgeId2(size);
		vector<vector<int>> send_pointId2(size);
		// vector<int> new_edge_center_point_id(edgesPoint0.size(),-1);
		for(int ip=0, tmp_point_size=0; ip<size; ++ip){
			int tmp_size = recv_edgeId[ip].size();
			for(int i=0; i<tmp_size; ++i){
				int id = recv_edgeId[ip][i];
                
                // if(bool_edges[id]==true) continue;
                // bool_edges[id] = true;
					
				int ipoint0 = edgesPoint0[id];
				int ipoint1 = edgesPoint1[id];
				
				int ibool = recv_boolEdgeRefine[ip][i];
				if(ibool==1 && boolEdgeRefine[id]==false){
					// F2T_proc_boolEdgeRefine[id] = true;
                    
                    if(bool_edges[id]==true) continue;
                    bool_edges[id] = true;
					
					// 포인트 생성
					int tmp_point_id = mesh.points.size();
					// new_edge_center_point_id[id] = tmp_point_id;
					proc_edgeCenterPoints_id[id] = tmp_point_id;
                    
                    // if(rank==419 && id==9321) cout << "AAAAAAAAAA" << endl;
					
					mesh.addPoint();
					mesh.points.back().x = 0.5*(mesh.points[ipoint0].x + mesh.points[ipoint1].x);
					mesh.points.back().y = 0.5*(mesh.points[ipoint0].y + mesh.points[ipoint1].y);
					mesh.points.back().z = 0.5*(mesh.points[ipoint0].z + mesh.points[ipoint1].z);
					mesh.points.back().level = max(
					mesh.points[ipoint0].level,mesh.points[ipoint1].level)+1;
					
					// connPoints 만들기 위한 재료
					for(auto& [proc, proc_id] : procEdges[id]){
						send_edgeId2[proc].push_back(proc_id);
						send_pointId2[proc].push_back(tmp_point_id);
					}
					
					// 페이스 포인트 추가
					for(auto& iface : edgesFaces[id]){
						// if(bool_addiFaces[iface]==true) continue;
						// bool_addiFaces[iface]=true;
						auto& face = mesh.faces[iface];
						int str = find(
							face.ipoints.begin(),face.ipoints.end(),
							ipoint0) - face.ipoints.begin();
						int end = find(
							face.ipoints.begin(),face.ipoints.end(),
							ipoint1) - face.ipoints.begin();
						if(str>end){
							int dummy = str;
							str = end;
							end = dummy;
						}
						
						// if(groupChildFaces_id[iface]!=-1){
							// cout << "ERROR1" << endl;
						// }
						// if(addiFaces[iface].ipoints.size()!=0){
							// cout << "ERROR2" << endl;
						// }
						// if(rank==2) cout << str << " " << end << " " 
						// << face.ipoints.size() << endl;
						if(str==0 && end==face.ipoints.size()-1){
							face.ipoints.insert(face.ipoints.end(),tmp_point_id);
						}
						else{
							if(str==0 && end!=1) cout << "#ERROR 22" << endl;
							if(end-str!=1) cout << "#ERROR 22" << endl;
							face.ipoints.insert(face.ipoints.begin()+end,tmp_point_id);
						}
						
					}
				}
				else if(boolEdgeRefine[id]==true){
					// T_proc_boolEdges[id] = true;
                    
                    if(bool_edges[id]==true) continue;
                    bool_edges[id] = true;
					
					int ipoint = proc_edgeCenterPoints_id[id];
					if(ipoint==-1) cout << "ERROR3" << endl;
					
					for(auto& [proc, proc_id] : procEdges[id]){
						send_edgeId2[proc].push_back(proc_id);
						send_pointId2[proc].push_back(ipoint);
					}
					
				}
			}
		}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// if(rank==0 && iter==1) cout << "GOOD2" << endl;
		vector<vector<int>> recv_edgeId2;
		mpi.Alltoallv(send_edgeId2, recv_edgeId2);
		vector<vector<int>> recv_pointId2;
		mpi.Alltoallv(send_pointId2, recv_pointId2);
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0 && iter==1) cout << "GOOD3" << endl;
		// connPoints 추가
		for(int ip=0, tmp_point_size=0; ip<size; ++ip){
			int tmp_size = recv_edgeId2[ip].size();
			for(int i=0; i<tmp_size; ++i){
				int iedge = recv_edgeId2[ip][i];
				int recv_ipoint = recv_pointId2[ip][i];
				
				int ipoint = proc_edgeCenterPoints_id[iedge];
				if(ipoint==-1) {
                    cout << "ERROR4" << endl;
                    cout << rank << " " << iedge << " " << endl;
                    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
                }
				
				int tmp_iter = 0;
				for(auto& [proc, id] : mesh.points[ipoint].connPoints){
					if(proc==ip && recv_ipoint==id) ++tmp_iter;
				}
				if(tmp_iter==0){
					mesh.points[ipoint].connPoints.push_back(
					make_pair(ip,recv_ipoint));
				}
			}
		}
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0 && iter==1) cout << "GOOD4" << endl;
	}
	
	// ==========================================================
	
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START4" << endl;
	
	
	
	
	
	// 그룹 페이스 제작
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];
		int faceLevel = face.level;
		
		if(
		groupChildFaces_id[i] != -1 &&
		groupChildFaces_HighLevel[i] == false){
			
			auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
			addFaceCenterPoint(mesh, faceLevel, 
				groupChildFace.vertexPoints,
				groupChildFace.centerPoint);
				
			// double residual = 0.0;
			// residual += abs(face.x-mesh.points[groupChildFace.centerPoint].x);
			// residual += abs(face.y-mesh.points[groupChildFace.centerPoint].y);
			// residual += abs(face.z-mesh.points[groupChildFace.centerPoint].z);
			// if(residual>1.e-16 && rank==0){
				// cout << residual << endl;
				// cout << face.x << " " << mesh.points[groupChildFace.centerPoint].x << endl;
				// cout << face.y << " " << mesh.points[groupChildFace.centerPoint].y << endl;
				// cout << face.z << " " << mesh.points[groupChildFace.centerPoint].z << endl;
			// }
				
			vector<int> tmpFacePoints;
			for(auto& j : face.ipoints){
				tmpFacePoints.push_back(j);
			}
			
			extractSubEdgePoints(
				tmpFacePoints,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints,
				groupChildFace.subOutEdgePoints);
				
			addSubOuterFacesPoints(
				groupChildFace.vertexPoints,
				groupChildFace.centerPoint,
				groupChildFace.subOutEdgePoints,
				groupChildFace.faces);
				
			extractSubInternalEdgesPoints(
				groupChildFace.vertexCenterPoints,
				groupChildFace.centerPoint,
				groupChildFace.subIntEdgePoints);
				
			for(auto& aaaaaaaaaaaaa : groupChildFace.faces){
				if(aaaaaaaaaaaaa.ipoints.size()>100){
					cout << "| WARNING 1 : face.ipoints.size > 10" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
		}
		else if(
		groupChildFaces_id[i] != -1 &&
		groupChildFaces_HighLevel[i] == true){
			
			auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
			
			vector<int> tmpFacePoints;
			for(auto& j : face.ipoints){
				tmpFacePoints.push_back(j);
			}
			
			
			groupChildFace.faces.push_back(faces_refind());
			for(auto& j : face.ipoints){
				groupChildFace.faces.back().ipoints.push_back(j);
			}
			
			auto iter0 = std::find(face.ipoints.begin(), face.ipoints.end(), groupChildFace.vertexPoints[1]);
			auto iter1 = std::find(face.ipoints.begin(), face.ipoints.end(), groupChildFace.vertexPoints[2]);
			auto iter2 = std::find(face.ipoints.begin(), face.ipoints.end(), groupChildFace.vertexPoints[3]);
				
			groupChildFace.subIntEdgePoints.resize(2);
			std::copy(iter0, iter1+1, std::back_inserter(groupChildFace.subIntEdgePoints[0]));
			std::copy(iter1, iter2+1, std::back_inserter(groupChildFace.subIntEdgePoints[1]));
			std::reverse(
				groupChildFace.subIntEdgePoints[1].begin(),
				groupChildFace.subIntEdgePoints[1].end());
			
			if(face.ipoints.size()>100){
				cout << "| WARNING 2 : face.ipoints.size > 10" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}
		
			
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START5" << endl;
	
	
	// =====================================
	// 페이스 센터 포인트 conn 연결
	// if(iter!=0)
	{
		vector<vector<int>> send_ipoint(size);
		vector<vector<double>> send_debug_ipoint_x(size);
		vector<vector<double>> send_debug_ipoint_y(size);
		vector<vector<double>> send_debug_ipoint_z(size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			int proc = boundary.rightProcNo;
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				if(
				groupChildFaces_id[i] != -1 &&
				groupChildFaces_HighLevel[i] == false){
					auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
					int ipoint = groupChildFace.centerPoint;
					send_ipoint[proc].push_back(ipoint);
					send_debug_ipoint_x[proc].push_back(mesh.points[ipoint].x);
					send_debug_ipoint_y[proc].push_back(mesh.points[ipoint].y);
					send_debug_ipoint_z[proc].push_back(mesh.points[ipoint].z);
				}
				else{
					send_ipoint[proc].push_back(-100);
					send_debug_ipoint_x[proc].push_back(face.x);
					send_debug_ipoint_y[proc].push_back(face.y);
					send_debug_ipoint_z[proc].push_back(face.z);
				}
			}
		}
		vector<vector<int>> recv_ipoint;
		mpi.Alltoallv(send_ipoint, recv_ipoint);
		vector<vector<double>> recv_debug_ipoint_x;
		mpi.Alltoallv(send_debug_ipoint_x, recv_debug_ipoint_x);
		vector<vector<double>> recv_debug_ipoint_y;
		mpi.Alltoallv(send_debug_ipoint_y, recv_debug_ipoint_y);
		vector<vector<double>> recv_debug_ipoint_z;
		mpi.Alltoallv(send_debug_ipoint_z, recv_debug_ipoint_z);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			int proc = boundary.rightProcNo;
			for(int i=str, tmpSize=0; i<end; ++i){
				auto& face = mesh.faces[i];
				if(
				groupChildFaces_id[i] != -1 &&
				groupChildFaces_HighLevel[i] == false){
					auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
					int ipoint = groupChildFace.centerPoint;
					int conn_ipoint = recv_ipoint[proc].at(tmpSize);
					mesh.points.at(ipoint).connPoints.push_back(
					make_pair(proc,conn_ipoint));
					
					double residual = 0.0;
					residual += abs(recv_debug_ipoint_x[proc][tmpSize]-mesh.points.at(ipoint).x);
					residual += abs(recv_debug_ipoint_y[proc][tmpSize]-mesh.points.at(ipoint).y);
					residual += abs(recv_debug_ipoint_z[proc][tmpSize]-mesh.points.at(ipoint).z);
					
					if(residual>1.e-7){
						cout << "#ERROR residual = " << residual << endl;
					}
					
					double x_avg = 0.5*(recv_debug_ipoint_x[proc][tmpSize]+mesh.points.at(ipoint).x);
					double y_avg = 0.5*(recv_debug_ipoint_y[proc][tmpSize]+mesh.points.at(ipoint).y);
					double z_avg = 0.5*(recv_debug_ipoint_z[proc][tmpSize]+mesh.points.at(ipoint).z);
					
					// mesh.points.at(ipoint).x = x_avg;
					// mesh.points.at(ipoint).y = y_avg;
					// mesh.points.at(ipoint).z = z_avg;
					
				}
				else{
					if(recv_ipoint[proc].at(tmpSize)!=-100){
						cout << "#ERROR send_ipoint[proc].at(tmpSize++)!=-100" << endl;
					}
					
					double residual = 0.0;
					residual += abs(recv_debug_ipoint_x[proc][tmpSize]-face.x);
					residual += abs(recv_debug_ipoint_y[proc][tmpSize]-face.y);
					residual += abs(recv_debug_ipoint_z[proc][tmpSize]-face.z);
					if(residual>1.e-7){
						cout << "#ERROR residual = " << residual << endl;
					}
					
					double x_avg = 0.5*(recv_debug_ipoint_x[proc][tmpSize]+face.x);
					double y_avg = 0.5*(recv_debug_ipoint_y[proc][tmpSize]+face.y);
					double z_avg = 0.5*(recv_debug_ipoint_z[proc][tmpSize]+face.z);
					
					// face.x = x_avg;
					// face.y = y_avg;
					// face.z = z_avg;
				}
				++tmpSize;
			}
		}
		
	}
	// =====================================
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START6" << endl;
	
	
	
	// 내부페이스 만들기
	vector<int> groupCellInternalFaces_id(mesh.cells.size(),-1);
	
	vector<bool> boolInputFacesiL(mesh.faces.size(),false);
	vector<bool> boolInputFacesiR(mesh.faces.size(),false);
	
	vector<vector<int>> groupCell_Levels;
	vector<vector<int>> groupCell_Groups;
	int totalCellNum = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		int cellLevel = cell.level;
		int cellGroup = cell.group;
		
		if(boolCellRefine[i]==true){
			
			groupCellInternalFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			

			vector<int> cellVertexPoints;
			for(auto& j : cell.ipoints){
				if(mesh.points[j].level > cellLevel) continue;
				cellVertexPoints.push_back(j);
			}
			
			int cellCenterPoint;
			addCellCenterPoint(mesh, cellLevel, 
				cellVertexPoints, cellCenterPoint);
			
			
			map<int,int> cellInternalFaces;
			for(auto& j : cell.ifaces){
				auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
				for(auto& k : groupOuterFace.vertexCenterPoints){
					if(cellInternalFaces.find(k) == cellInternalFaces.end()){
						int tmpSize = cellInternalFaces.size();
						cellInternalFaces.insert(make_pair(k,tmpSize));
						groupChildFace.faces.push_back(faces_refind());
					}
				}
			}
			
			
			map<int,int> cellVertexOrder;
			extractCellVertexOrder(cellVertexPoints, cellVertexOrder);
			
			
			vector<vector<int>> intFacesVertexs;
			extractInternalFacesVertexs(cell, cellInternalFaces,
				groupChildFaces_id, groupChildFaces, groupChildFace,
				cellCenterPoint, intFacesVertexs);
				
			for(auto& test : intFacesVertexs){
				if( test.size() != 2 ){
					cout << "| WARNING : no own ngb" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
				
			addInternalFacesiLiR(mesh, cellLevel, totalCellNum,
				cellVertexOrder, intFacesVertexs, groupChildFace.faces);
				
			for(auto& j : cell.ifaces){
				auto& face = mesh.faces[j];
				auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
				
				int verPSize = groupOuterFace.faces.size();
				
				for(int k=0; k<verPSize; ++k){
					int verP = groupOuterFace.vertexPoints[k];
					auto& subFace = groupOuterFace.faces[k];
					// cout << subface.ipoints.size() << endl;
					
					if(cellVertexOrder.find(verP) == cellVertexOrder.end()){
						cout << "| WARNING : no find cellVertexOrder.find(verP)" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					
					if(face.iL == i && subFace.booliL==false){
						subFace.iL = totalCellNum + cellVertexOrder[verP];
						subFace.booliL = true;
					}
					else if(face.iR == i && subFace.booliR==false){
						subFace.iR = totalCellNum + cellVertexOrder[verP];
						subFace.booliR = true;
					}
					else{
						cout << "NONONONONONONONONO 33" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
				}
			}
			
			groupCell_Levels.push_back(vector<int>());
			for(int jj=0; jj<cellVertexOrder.size(); ++jj){
				groupCell_Levels.back().push_back(cellLevel+1);
			}
			
			groupCell_Groups.push_back(vector<int>());
			for(int jj=0; jj<cellVertexOrder.size(); ++jj){
				groupCell_Groups.back().push_back(cellGroup);
			}
			
			totalCellNum += cellVertexOrder.size();
			
		}
		else{
			reorderOuterFacesiLiR(mesh, cell, i, totalCellNum,
				groupChildFaces_id, boolInputFacesiL, boolInputFacesiR, 
				groupChildFaces);
				
			
			groupCell_Levels.push_back(vector<int>());
			groupCell_Levels.back().push_back(cellLevel);
			
			groupCell_Groups.push_back(vector<int>());
			groupCell_Groups.back().push_back(cellGroup);
			
			++totalCellNum;
			
		}
		
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7" << endl;
	
	
	// 페이스 넘버링
	int faceResizeNum = 0;
	vector<int> startFaces(mesh.faces.size(),-1);
	
	// 셀 겉면 페이스 넘버링
	for(int i=0; i<mesh.faces.size(); ++i){
		if(mesh.faces.at(i).getType() == MASCH_Face_Types::INTERNAL){
			if(groupChildFaces_id.at(i) != -1){
				auto& groupChildFace = groupChildFaces.at(groupChildFaces_id.at(i));
				startFaces.at(i) = faceResizeNum;
				faceResizeNum += groupChildFace.faces.size();
			}
			else{
				startFaces.at(i) = faceResizeNum;
				++faceResizeNum;
			}
		}
	}
    
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.1" << endl;
	
	// 셀 안쪽 페이스 넘버링
	vector<int> startIntFaces(mesh.cells.size(),-1);
	for(int i=0; i<mesh.cells.size(); ++i){
		if(groupCellInternalFaces_id.at(i) != -1){
			auto& groupChildFace = groupChildFaces.at(groupCellInternalFaces_id.at(i));
			startIntFaces.at(i) = faceResizeNum;
			faceResizeNum += groupChildFace.faces.size();
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.2" << endl;
	
	// proc & 바운더리 페이스 넘버링
	for(int i=0; i<mesh.faces.size(); ++i){
		if(mesh.faces.at(i).getType() != MASCH_Face_Types::INTERNAL){
			if(groupChildFaces_id.at(i) != -1){
				auto& groupChildFace = groupChildFaces.at(groupChildFaces_id.at(i));
				startFaces.at(i) = faceResizeNum;
				faceResizeNum += groupChildFace.faces.size();
			}
			else{
				startFaces.at(i) = faceResizeNum;
				++faceResizeNum;
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.3" << endl;
	
	// 원래 페이스 다시 제작
	int orgFaceSize = mesh.faces.size();
	int orgCellSize = mesh.cells.size();
	
	MPI_Barrier(MPI_COMM_WORLD); 
	if(rank==0 && boolDebug) cout << "START7.3.1" << endl;
    
    if(faceResizeNum>100000000) cout << "rank = " << rank << ", face number > 1 million ";
	mesh.faces.resize(faceResizeNum);
    
    
    if(faceResizeNum==0){
        cout << "NONONONOONONON " << beforeCellSize << " " << beforeFaceSize << endl;
    }
    
    
    
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.3.2" << endl;
	
	// Proc & B.C. faces	
	for(auto& boundary : mesh.boundaries){
		// if(boundary.neighbProcNo == -1) continue;
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY) continue;
		if(boundary.rightProcNo <= rank) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
        // cout << "AA" << endl;
		for(int i=str; i<end; ++i){
			if(groupChildFaces_id.at(i) != -1){
				auto& groupChildFace = groupChildFaces.at(groupChildFaces_id.at(i));
				std::reverse(groupChildFace.faces.begin()+1,groupChildFace.faces.end());
			}
		}
        // cout << "BB" << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.4" << endl;
	
	int saveI = 0;
	for(int i=orgFaceSize-1; i>=0; --i){
		if(mesh.faces.at(i).getType() != MASCH_Face_Types::INTERNAL){
			int str = startFaces.at(i);
			if(groupChildFaces_id.at(i) != -1){
				auto& groupChildFace = groupChildFaces.at(groupChildFaces_id.at(i));

				int tmpNum = 0;
				for(auto& face : groupChildFace.faces){
					mesh.faces[str+tmpNum].ipoints.clear();
					for(auto& j : face.ipoints){
						mesh.faces.at(str+tmpNum).ipoints.push_back(j);
					}
					mesh.faces.at(str+tmpNum).iL = face.iL;
					mesh.faces.at(str+tmpNum).iR = -1;
					mesh.faces.at(str+tmpNum).setType(mesh.faces.at(i).getType());
					
					
					++tmpNum;
				}
			}
			else{
				mesh.faces[str] = mesh.faces[i];
				
				// // ========================
				// if(addiFaces[i].ipoints.size()>0){
					// mesh.faces[str].ipoints.clear();
					// for(auto& ipoint : addiFaces[i].ipoints){
						// mesh.faces[str].ipoints.push_back(ipoint);
					// }
				// }
				// // ========================
				
				
			}
			saveI = i;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.5" << endl;
	
	
	// cell internal faces
	for(int i=orgCellSize-1; i>=0; --i){
		if(groupCellInternalFaces_id.at(i) != -1){
			int str = startIntFaces.at(i);
			auto& groupChildFace = groupChildFaces.at(groupCellInternalFaces_id.at(i));
			int tmpNum = 0;
			for(auto& face : groupChildFace.faces){
				mesh.faces.at(str+tmpNum).ipoints.clear();
				for(auto& j : face.ipoints){
					mesh.faces.at(str+tmpNum).ipoints.push_back(j);
				}
				mesh.faces.at(str+tmpNum).iL = face.iL;
				mesh.faces.at(str+tmpNum).iR = face.iR;
				mesh.faces.at(str+tmpNum).setType(MASCH_Face_Types::INTERNAL);
				++tmpNum;
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.6" << endl;
	
	// cell outer faces
	for(int i=saveI-1; i>=0; --i){
		if(mesh.faces.at(i).getType() == MASCH_Face_Types::INTERNAL){
			int str = startFaces.at(i);
			if(groupChildFaces_id.at(i) != -1){
				auto& groupChildFace = groupChildFaces.at(groupChildFaces_id.at(i));
				int tmpNum = 0;
				for(auto& face : groupChildFace.faces){
					mesh.faces.at(str+tmpNum).ipoints.clear();
					for(auto& j : face.ipoints){
						mesh.faces.at(str+tmpNum).ipoints.push_back(j);
					}
					mesh.faces.at(str+tmpNum).iL = face.iL;
					mesh.faces.at(str+tmpNum).iR = face.iR;
					mesh.faces.at(str+tmpNum).setType(MASCH_Face_Types::INTERNAL);
					++tmpNum;
				}
			}
			else{
				mesh.faces.at(str) = mesh.faces.at(i);
				
				// // ========================
				// if(addiFaces[i].ipoints.size()>0){
					// mesh.faces[str].ipoints.clear();
					// for(auto& ipoint : addiFaces[i].ipoints){
						// mesh.faces[str].ipoints.push_back(ipoint);
					// }
				// }
				// // ========================
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.7" << endl;
	

	// boundary setting
	for (int i=0; i<mesh.boundaries.size(); ++i) {
        // if(mesh.boundaries.at(i).startFace>=startFaces.size()) cout << i << " " << mesh.boundaries.at(i).startFace << endl;
        
        // if(mesh.boundaries.at(i).startFace==0){
            // cout << startFaces.size() << " " << mesh.faces.size() << endl;
        // }
        
		mesh.boundaries.at(i).startFace = startFaces.at( mesh.boundaries.at(i).startFace );
	}
    
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.7.1" << endl;
	
	
	for (int i=0; i<mesh.boundaries.size()-1; ++i) {;
		mesh.boundaries.at(i).nFaces = mesh.boundaries.at(i+1).startFace-mesh.boundaries.at(i).startFace;
	}
    
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.7.2" << endl;
	
	int maxBDsize = mesh.boundaries.size()-1;
	mesh.boundaries.at(maxBDsize).nFaces = mesh.faces.size()-mesh.boundaries.at(maxBDsize).startFace;
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7.8" << endl;
	
	int tmpCellNum = totalCellNum-1;
	mesh.cells.resize(totalCellNum,MASCH_Cell());
	// for(int i=orgCellSize-1; i>=0; --i){
	boolCellPreserved.resize(totalCellNum);
	boolCellRefine.resize(totalCellNum);
	boolCellUnrefine.resize(totalCellNum);
	for(int i=orgCellSize-1; i>=0; --i){
		int subCellSize = groupCell_Levels.at(i).size();
		
		child_new_cell_id_of_org.at(i).resize(subCellSize);
		
		for(int j=0; j<subCellSize; ++j){
			// mesh.cells[tmpCellNum].var.assign(mesh.cells[i].var.begin(),mesh.cells[i].var.end());
			
			child_new_cell_id_of_org.at(i).at(j) = tmpCellNum;
			
			boolCellRefine.at(tmpCellNum) = boolCellRefine.at(i);
			boolCellUnrefine.at(tmpCellNum) = boolCellUnrefine.at(i);
			boolCellPreserved.at(tmpCellNum) = boolCellPreserved.at(i);
			// boolCellPreserved[tmpCellNum] = true;
			
			
			
			mesh.cells.at(tmpCellNum).ipoints.clear();
			mesh.cells.at(tmpCellNum).ifaces.clear();
			--tmpCellNum;
		}
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START8" << endl;
	
	
	// mesh.check();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces(); 
	mesh.cellsGlobal();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	// mesh.informations();
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START9" << endl;
	
	
	
	
	// level setting
	tmpCellNum = 0;
	for(auto& groupCell : groupCell_Levels){
		for(auto& level : groupCell){
			mesh.cells[tmpCellNum].level = level;
			++tmpCellNum;
		}
	}
	tmpCellNum = 0;
	for(auto& groupCell : groupCell_Groups){
		for(auto& group : groupCell){
			mesh.cells[tmpCellNum].group = group;
			++tmpCellNum;
		}
	}
	{
		int total_nGroup = 0;
		int tmp_cell_size = mesh.cells.size();
		vector<int> str_icells(size+1,0);
		for(int i=0; i<tmp_cell_size; ++i){
			auto& cell0 = mesh.cells[i];
			int group0 = cell0.group;
			cell0.group = total_nGroup;
			while(1){
				++i;
				if(i==tmp_cell_size) break;
				auto& cell = mesh.cells[i];
				int group = cell.group;
				if(group0!=group) break;
				cell.group = total_nGroup;
			}
			--i;
			++total_nGroup;
		}
		
		if(size>1){
			vector<int> recv_ncells(size);
			MPI_Allgather(&total_nGroup, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
			for(int ip=0; ip<size; ++ip){
				str_icells[ip+1] = str_icells[ip] + recv_ncells[ip];
			}
		}
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			mesh.cells[i].group += str_icells[rank];
		}
	}
	

	// mesh.debug_group_procFaces(0.99);
	
	
	
	

	// //=======================================
	// // 디버그
	// {
		
		// vector<vector<int>> tmp_procFace_id(size);
		// vector<vector<int>> tmp_procFace_group(size);
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			// int rightProcNo = boundary.rightProcNo;
			// for(int i=str, iter=0; i<end; ++i){
				// auto& face = mesh.faces[i];
				// int iL = face.iL;
				// // int igroup = mesh.cells[iL].group-min_group_id;
				// int igroup = mesh.cells[iL].group;
				// tmp_procFace_id[rightProcNo].push_back(i);
				// tmp_procFace_group[rightProcNo].push_back(igroup);
			// }
		// }
	
		// MASCH_Math math;
	
		
		// vector<vector<vector<int>>> reorder_procFace_id(size);
		// vector<vector<vector<double>>> reorder_procFace_x(size);
		// vector<vector<vector<double>>> reorder_procFace_y(size);
		// vector<vector<vector<double>>> reorder_procFace_z(size);
		// for(int ip=0; ip<size; ++ip){
			// int size_procFace = tmp_procFace_group[ip].size();
				
			// vector<vector<int>> tmp_face_group;
			// vector<vector<double>> tmp_face_group_x;
			// vector<vector<double>> tmp_face_group_y;
			// vector<vector<double>> tmp_face_group_z;
				
			// // sort(tmp_procFace_group[ip].begin(),tmp_procFace_group[ip].end());
			// for(int i=0; i<size_procFace; ++i){
				// int id0 = tmp_procFace_id[ip][i];
				// int igroup0 = tmp_procFace_group[ip][i];
				// vector<int> tmp_id;
				// tmp_id.push_back(id0);
				// while(1){
					// ++i;
					// if(i==size_procFace) break;
					// int id = tmp_procFace_id[ip][i];
					// int igroup = tmp_procFace_group[ip][i];
					// if(igroup0!=igroup) break;
					// tmp_id.push_back(id);
				// }
				// --i;
				// {
					// vector<vector<double>> procFace_unitNormals;
					// for(auto& iface : tmp_id){
						// auto& face = mesh.faces[iface];
						// vector<double> Vx, Vy, Vz;
						// for(auto& ipoint : face.ipoints){
							// Vx.push_back(mesh.points[ipoint].x);
							// Vy.push_back(mesh.points[ipoint].y);
							// Vz.push_back(mesh.points[ipoint].z);
						// }
						
						// double VSn=0.0;
						// vector<double> cellCentroid;
						// vector<double> tmp_unitNormals(3,0.0);
						// double area;
						// math.calcUnitNormals_Area3dPolygon(
							// face.ipoints.size(), Vx,Vy,Vz,
							// tmp_unitNormals, area,
							// face.x, face.y, face.z,
							// VSn, cellCentroid);
							
						// procFace_unitNormals.push_back(tmp_unitNormals);
					// }
					
					// int tmp_faceSize = procFace_unitNormals.size();
					// for(int j=0; j<tmp_faceSize; ++j){
						// double normal0x = procFace_unitNormals[j][0];
						// double normal0y = procFace_unitNormals[j][1];
						// double normal0z = procFace_unitNormals[j][2];
						
						// vector<int> tmp_face;
						// vector<double> tmp_x;
						// vector<double> tmp_y;
						// vector<double> tmp_z;
						// tmp_face.push_back(tmp_id[j]);
						// tmp_x.push_back(normal0x);
						// tmp_y.push_back(normal0y);
						// tmp_z.push_back(normal0z);
						
						// while(1){
							// ++j;
							// if(j==tmp_faceSize) break;
							// double normalx = procFace_unitNormals[j][0];
							// double normaly = procFace_unitNormals[j][1];
							// double normalz = procFace_unitNormals[j][2];
							
							// double resi = 0.0;
							// resi += normalx*normal0x;
							// resi += normaly*normal0y;
							// resi += normalz*normal0z;
							// if(resi>0.5){
								// tmp_face.push_back(tmp_id[j]);
								// tmp_x.push_back(normalx);
								// tmp_y.push_back(normaly);
								// tmp_z.push_back(normalz);
							// }
							// else{
								// break;
							// }
							
						// }
						// --j;
						
						// tmp_face_group.push_back(tmp_face);
						// tmp_face_group_x.push_back(tmp_x);
						// tmp_face_group_y.push_back(tmp_y);
						// tmp_face_group_z.push_back(tmp_z);
					// }
				// }
			// }
				
			// reorder_procFace_id[ip] = (tmp_face_group);
			// reorder_procFace_x[ip] = (tmp_face_group_x);
			// reorder_procFace_y[ip] = (tmp_face_group_y);
			// reorder_procFace_z[ip] = (tmp_face_group_z);
		// }
		
		
		// if(rank==0){
			// int iter=0;
			// for(auto& item : reorder_procFace_id[1]){
				// int iter2=0;
				// for(auto& item2 : item){
					// cout << rank << ", " << iter << ", " << 
					// mesh.cells[mesh.faces[item2].iL].group << ", " << 
					// iter2 << ", " << item2 << 
					// ", " << mesh.faces[item2].x <<
					// ", " << mesh.faces[item2].y <<
					// ", " << mesh.faces[item2].z <<
					// ", " << reorder_procFace_x[1][iter][iter2] <<
					// ", " << reorder_procFace_y[1][iter][iter2] <<
					// ", " << reorder_procFace_z[1][iter][iter2] <<
					// endl;
					// for(auto& ipoint : mesh.faces[item2].ipoints){
						// cout << 
						// mesh.points[ipoint].x <<
						// ", " << mesh.points[ipoint].y <<
						// ", " << mesh.points[ipoint].z << endl;
					// }
					// ++iter2;
				// }
				// ++iter;
			// }
		// }
		
		// MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==1 || rank==0){
			// cout << endl;
			// cout << endl;
			// cout << endl;
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
		
		// if(rank==1){
			// int iter=0;
			// for(auto& item : reorder_procFace_id[0]){
				// int iter2=0;
				// for(auto& item2 : item){
					// cout << rank << ", " << iter << ", " << 
					// mesh.cells[mesh.faces[item2].iL].group << ", " << 
					// iter2 << ", " << item2 << 
					// ", " << mesh.faces[item2].x <<
					// ", " << mesh.faces[item2].y <<
					// ", " << mesh.faces[item2].z <<
					// ", " << reorder_procFace_x[0][iter][iter2] <<
					// ", " << reorder_procFace_y[0][iter][iter2] <<
					// ", " << reorder_procFace_z[0][iter][iter2] <<
					// endl;
					// for(auto& ipoint : mesh.faces[item2].ipoints){
						// cout << 
						// mesh.points[ipoint].x <<
						// ", " << mesh.points[ipoint].y <<
						// ", " << mesh.points[ipoint].z << endl;
					// }
					// ++iter2;
				// }
				// ++iter;
			// }
		// }
	// }
			


		// //=======================================
		// // 디버그
		// {
			// int total_nGroup2 = 0;
			// int tmp_cell_size = mesh.cells.size();
			// vector<int> str_icells2(size+1,0);
			// for(int i=0; i<tmp_cell_size; ++i){
				// auto& cell0 = mesh.cells[i];
				// int group0 = cell0.group;
				// // cell0.group = total_nGroup2;
				// while(1){
					// ++i;
					// if(i==tmp_cell_size) break;
					// auto& cell = mesh.cells[i];
					// int group = cell.group;
					// if(group0!=group) break;
					// // cell.group = total_nGroup2;
				// }
				// --i;
				// ++total_nGroup2;
			// }
			// if(size>1){
				// vector<int> recv_ncells(size);
				// MPI_Allgather(&total_nGroup2, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
				// for(int ip=0; ip<size; ++ip){
					// str_icells2[ip+1] = str_icells2[ip] + recv_ncells[ip];
				// }
			// }
			// vector<vector<int>> groupCellsProcFaces(total_nGroup2);
			// vector<vector<int>> groupCellsRightProcNo(total_nGroup2);
			// {
				// for(auto& boundary : mesh.boundaries){
					// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
					// int str = boundary.startFace;
					// int end = str + boundary.nFaces;
					// int rightProcNo = boundary.rightProcNo;
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
						// int igroup = mesh.cells[iL].group-str_icells2[rank];
						// groupCellsProcFaces[igroup].push_back(i);
						// groupCellsRightProcNo[igroup].push_back(rightProcNo);
					// }
				// }
			// }
			// for(int i=0; i<groupCellsProcFaces.size(); ++i){
				// int org_Rip0=-1;
				// if(groupCellsRightProcNo[i].size()!=0) 
					// org_Rip0 = groupCellsRightProcNo[i][0];
				// for(int j=0; j<groupCellsProcFaces[i].size(); ++j){
					// int org_id = groupCellsProcFaces[i][j];
					// int org_Rip = groupCellsRightProcNo[i][j];
					// if(org_Rip0==org_Rip) cout << i << " " << org_Rip << " " << org_id << endl;
				// }
			// }
		// }
		// //=======================================

	
	
	
	
	
	this->mpiLevels(mesh, cLevel_recv);
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// face.var.clear();
		// face.varL.clear();
		// face.varR.clear();
	// }
	
	for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		// face.var.resize(controls.nTotalFaceVar,0.0);
		// face.varL.resize(controls.nTotalFaceLRVar,0.0);
		// face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			
			if(face.iL > mesh.cells.size()-1){
				cout << " face.iL > mesh.cells.size()-1 " << face.iL << endl;
			}
			if(face.iR < 0 || face.iR > mesh.cells.size()-1){
				cout << " face.iR < 0 || face.iR > mesh.cells.size()-1 " << face.iR << endl;
			}
			
			int maxLevel = 
				max(mesh.cells[face.iL].level,
					mesh.cells[face.iR].level);
			face.level = maxLevel;
			
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			
			if(proc_num > cLevel_recv.size()-1){
				cout << " proc_num > cLevel_recv.size()-1 " << " " << cLevel_recv.size() << " " << proc_num << endl;
			}
			
			int maxLevel = 
				max(mesh.cells[face.iL].level,
					cLevel_recv[proc_num]);
			face.level = maxLevel;
			
			++proc_num;
			
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
			
			if(face.iL > mesh.cells.size()-1){
				cout << " face.iL > mesh.cells.size()-1 " << face.iL << endl;
			}
			
			face.level = mesh.cells[face.iL].level;
		}
	}
	
	
	
	// ===========================================
	// connPoints 디버깅
	mesh.debug_connPoints(1.e-12);
	// proc face points 디버깅
	mesh.debug_procFacePoints(1.e-12);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// ===========================================
	
	
	
	
	// AMR 결과 프린트
	afterCellSize = mesh.cells.size();
	afterFaceSize = mesh.faces.size();
	afterPointSize = mesh.points.size();
	
	int addedCellSize = afterCellSize - beforeCellSize;
	int addedFaceSize = afterFaceSize - beforeFaceSize;
	int addedPointSize = afterPointSize - beforePointSize;
	
	int addedCellSize_glb = addedCellSize, 
    addedFaceSize_glb = addedFaceSize, 
    addedPointSize_glb = addedPointSize;
	if(size>1){
		MPI_Allreduce(&addedCellSize, &addedCellSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&addedFaceSize, &addedFaceSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&addedPointSize, &addedPointSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		cout << "-> completed" << endl;
		cout << 
		"| cell = +" << addedCellSize_glb <<
		" | face = +" << addedFaceSize_glb <<
		" | point = +" << addedPointSize_glb << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}




void MASCH_Poly_Mesh_Refine::separateEdgesPoints(
	MASCH_Mesh& mesh, 
	int faceLevel,
	vector<int>& points, 
	vector<int>& edges,
	vector<bool>& boolEdgeRefine,
	vector<int>& edgesCenterPointNumber,
	vector<int>& edgesLevel,
	vector<int>& faceVertex,
	vector<int>& edgeCenterPoints,
	vector<vector<int>>& subFaceEdgePoints,
	int& newPointNumber){


	faceVertex.clear();
	edgeCenterPoints.clear();
	subFaceEdgePoints.clear();

	vector<int> tmpVector;
	int iEdge;
	for(int k=0; k<points.size(); ++k){
		
		
		
		int iPoint = points[k];
		int pointLevel = mesh.points[iPoint].level;
		iEdge = edges[k];
		bool saveEdgeRefine = boolEdgeRefine[iEdge];
		int edgeLevel = edgesLevel[iEdge];
		
		if( pointLevel <= faceLevel ){
			faceVertex.push_back(iPoint);
			
			if( tmpVector.size() > 0 ){
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
			
	// cout << k << " " << points.size() << endl;
			if(
			faceVertex.size() != 1 &&
			edgeCenterPoints.size() != faceVertex.size()-1 &&
			saveEdgeRefine == false
			) {
				
				int backVertexNum = faceVertex.back();
				
				boolEdgeRefine[iEdge] = true;
				
				// edgesCenterPointNumber[iEdge] = newPointNumber++;
				edgesCenterPointNumber[iEdge] = mesh.points.size();
				// this->addCenterPoint(mesh, cellVertex, cellLevel);
				
				vector<int> vertex;
				vertex.push_back(backVertexNum);
				vertex.push_back(iPoint);
				this->addCenterPoint(mesh, vertex, faceLevel);
				
				edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
				
				tmpVector.push_back(edgesCenterPointNumber[iEdge]);
			}
			
			
		}
		else if( pointLevel == faceLevel+1 ){
			
			edgeCenterPoints.push_back(iPoint);
			
			if( tmpVector.size() > 0 ){
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
		}
		
		tmpVector.push_back(iPoint);
		
		if( saveEdgeRefine == true ){
			if( faceLevel == edgeLevel ){
				edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
			
			tmpVector.push_back(edgesCenterPointNumber[iEdge]);
			
		}
	}
	
	subFaceEdgePoints.push_back(tmpVector);
	
	if( edgeCenterPoints.size() != faceVertex.size() ){
		boolEdgeRefine[iEdge] = true;
		// edgesCenterPointNumber[iEdge] = newPointNumber++;
		edgesCenterPointNumber[iEdge] = mesh.points.size();
		
		vector<int> vertex;
		vertex.push_back(faceVertex[0]);
		vertex.push_back(faceVertex.back());
		this->addCenterPoint(mesh, vertex, faceLevel);
		
		edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
	}
	
	if( edgeCenterPoints.size() != faceVertex.size() ){
		cout << "| #Error 1" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
}






void MASCH_Poly_AMR_Builder::extractVertexPoints(
	int faceLevel,
	vector<int>& facePoints,
	vector<int>& facePointLevels,
	vector<int>& vertexPoints
) {
	vertexPoints.clear();

	int facePointsSize = (int)facePoints.size();
	for (int i = 0; i < facePointsSize; ++i) {
		if(facePointLevels[i] <= faceLevel){
			vertexPoints.push_back(facePoints[i]);
		}
	}
};


void MASCH_Poly_AMR_Builder::extractVertexCenterPoints(
	int faceLevel,
	vector<int>& facePoints,
	vector<int>& facePointLevels,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints,
	vector<bool>& addVertexCeterPoints
) {
	vertexCenterPoints.clear();
	addVertexCeterPoints.clear();
	// edgeOrders.clear();

	int vertexPointsSize = (int)vertexPoints.size();
	int facePointsSize = (int)facePoints.size();
	
	// vector<int> facePoints_copy = facePoints;

	int num = 0;
	// auto it = facePoints_copy.begin();
	for (int i = 0; i < facePointsSize; ++i) {
		
		// cout << facePoints[i] << " " << facePointLevels[i] << endl;
		
		
		// cout << i << " " << facePoints.size() << " " << num << " " << vertexPoints.size() << endl;
		// ++it;
		if (facePoints[i] == vertexPoints[num]) {
			int iNext = (i + 1 < facePointsSize) ? i + 1 : 0;
			int numNext = (num + 1 < vertexPointsSize) ? num + 1 : 0;
			
			// cout << iNext << " " << facePoints.size() << " " << numNext << " " << vertexPoints.size() << endl;
			
			
			if (facePoints[iNext] == vertexPoints[numNext]) {
				// // create vertexCenterPoints
				// int newPoint = (int)facePoints_copy.size();
				// it = ++facePoints_copy.insert(it, newPoint);
				int newPoint = -1;
				vertexCenterPoints.push_back(newPoint);
				addVertexCeterPoints.push_back(true);
				// edgeOrders.push_back(i);
			}
			++num;
			if(num >= vertexPointsSize-1) num = vertexPointsSize-1;
		}
		else {
			if (facePointLevels[i] == faceLevel + 1) {
				vertexCenterPoints.push_back(facePoints[i]);
				addVertexCeterPoints.push_back(false);
				// edgeOrders.push_back(i);
			}
		}
	}

	// facePoints = facePoints_copy;
};



void MASCH_Poly_AMR_Builder::extractSubEdgePoints(
	vector<int>& facePoints,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints,
	vector<vector<int>>& subEdgePoints
) {

	int vertexCeterPointsSize = (int)vertexCenterPoints.size();
	int subEdgeSize = vertexCeterPointsSize * 2;
	subEdgePoints.resize(subEdgeSize,vector<int>());

	int i = 0;
	auto iter0 = find(facePoints.begin(), facePoints.end(), vertexPoints[0]);
	for (i = 0; i < vertexCeterPointsSize-1; ++i) {
		auto iter1 = find(iter0, facePoints.end(), vertexCenterPoints[i]);
		auto iter2 = find(iter1, facePoints.end(), vertexPoints[i+1]);

		int dist = (int)distance(iter0, iter1);
		subEdgePoints[i * 2].insert(subEdgePoints[i * 2].begin(), iter0, iter1);
		//cout << dist << endl;

		dist = (int)distance(iter1, iter2);
		subEdgePoints[i * 2+1].insert(subEdgePoints[i * 2+1].begin(), iter1, iter2);
		//cout << dist << endl;

		iter0 = iter2;
	}
	//cout << "end" << endl;
	auto iter1 = find(iter0, facePoints.end(), vertexCenterPoints[i]);
	auto iter2 = facePoints.end();
	int dist = (int)distance(iter0, iter1);
	subEdgePoints[i * 2].insert(subEdgePoints[i * 2].begin(), iter0, iter1);
	//cout << dist << endl;

	dist = (int)distance(iter1, iter2);
	subEdgePoints[i * 2+1].insert(subEdgePoints[i * 2+1].begin(), iter1, iter2);
	//cout << dist << endl;

};



void MASCH_Poly_AMR_Builder::addVertexCenterPoint(
	MASCH_Mesh& mesh, 
	int faceLevel,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints
) {

	int verNum = vertexCenterPoints.size();
	for(int j=0; j<verNum; ++j){
		// if(addVertexCeterPoints[j]==true){
		if(vertexCenterPoints[j] == -1){
			vertexCenterPoints[j] = mesh.points.size();
			mesh.addPoint();
			int pointNum0 = j;
			int pointNum1 = (j+1<verNum) ? j+1 : 0;
			int point0 = vertexPoints[pointNum0];
			int point1 = vertexPoints[pointNum1];
			
			if( point0 >= mesh.points.size() ) cout << "| WARNING : point0 >= point.size" << endl;
			if( point1 >= mesh.points.size() ) cout << "| WARNING : point0 >= point.size" << endl;
			// cout << mesh.points.size() << endl;
			// cout << pointNum0 << " " << pointNum1 << " " << point0 << " " << point1 << endl;
			
			mesh.points.back().x = 0.5*(mesh.points[point0].x + mesh.points[point1].x);
			mesh.points.back().y = 0.5*(mesh.points[point0].y + mesh.points[point1].y);
			mesh.points.back().z = 0.5*(mesh.points[point0].z + mesh.points[point1].z);
			mesh.points.back().level = faceLevel+1;
		}
	}
				

};



void MASCH_Poly_AMR_Builder::addFaceCenterPoint(
	MASCH_Mesh& mesh, 
	int faceLevel,
	vector<int>& vertexPoints,
	int& centerPoint
) {

	centerPoint = mesh.points.size();
	mesh.addPoint();
	mesh.points.back().x = 0.0;
	mesh.points.back().y = 0.0;
	mesh.points.back().z = 0.0;
	int verNum = vertexPoints.size();
	double dVerNum = static_cast<double>(verNum);
	for(auto& j : vertexPoints){
		mesh.points.back().x += mesh.points[j].x;
		mesh.points.back().y += mesh.points[j].y;
		mesh.points.back().z += mesh.points[j].z;
	}
	mesh.points.back().x /= dVerNum;
	mesh.points.back().y /= dVerNum;
	mesh.points.back().z /= dVerNum;
	mesh.points.back().level = faceLevel+1;
	
};


void MASCH_Poly_AMR_Builder::addEdgeFacesPointsToVertexCenterPoint(
	MASCH_Mesh& mesh,
	int i,
	vector<vector<int>>& facesEdges,
	vector<vector<int>>& edgesFaces,
	vector<bool>& boolEdgeRefine,
	vector<bool>& addVertexCeterPoints,
	vector<int>& canRefineEdgeOrders,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints
) {
	int vertexSize = vertexPoints.size();
	int tmpNewNum = 0;
	for(int j=0; j<vertexSize; ++j){
		if(addVertexCeterPoints[j]==true){
			// int ordEdge = canRefineEdgeOrders[tmpNewNum];
			// ++tmpNewNum;
			// cout << ordEdge << " " << j << endl;
			// int iEdge = facesEdges[i][ordEdge];
			
			int iEdge = canRefineEdgeOrders[tmpNewNum];
			++tmpNewNum;
			
			// if(boolEdgeRefine[iEdge]==false){
				boolEdgeRefine[iEdge] = true;
				int isertP = vertexCenterPoints[j];
				for(auto& k : edgesFaces[iEdge]){
					
					int dist = 0;
					for(int l=0; l<facesEdges[k].size(); ++l){
						int iEdge2 = facesEdges[k][l];
						++dist;
						if(iEdge2 == iEdge) break;
						if(boolEdgeRefine[iEdge2]==true){
							++dist;
						}
					}
					mesh.faces[k].ipoints.insert(
						mesh.faces[k].ipoints.begin()+dist,
						isertP);
					
				}
			// }
		}
	}
			
};



void MASCH_Poly_AMR_Builder::addCellCenterPoint(
	MASCH_Mesh& mesh, 
	int cellLevel,
	vector<int>& vertexPoints,
	int& centerPoint
) {

	centerPoint = mesh.points.size();
	mesh.addPoint();
	mesh.points.back().x = 0.0;
	mesh.points.back().y = 0.0;
	mesh.points.back().z = 0.0;
	int verNum = vertexPoints.size();
	double dVerNum = 1.0/(double)verNum;
	for(auto& j : vertexPoints){
		mesh.points.back().x += mesh.points[j].x*dVerNum;
		mesh.points.back().y += mesh.points[j].y*dVerNum;
		mesh.points.back().z += mesh.points[j].z*dVerNum;
	}
	mesh.points.back().level = cellLevel+1;
	
};




void MASCH_Poly_AMR_Builder::addSubOuterFacesPoints(
	vector<int>& vertexPoints,
	int& centerPoint,
	vector<vector<int>>& subOutEdgePoints,
	vector<faces_refind>& faces
) {

	int vertexSize = (int)vertexPoints.size();
	for(int j=0; j<vertexSize; ++j){
		int leftEdge = ( j==0 ? vertexSize*2-1 : j*2-1 );
		int centerEdge = j*2; 
		int rightEdge = j*2+1;
		
		faces.push_back(faces_refind());
		auto& tmp_Face = faces.back();
		for(auto& k : subOutEdgePoints[centerEdge]){
			tmp_Face.ipoints.push_back(k);
		}
		tmp_Face.ipoints.push_back(subOutEdgePoints[rightEdge][0]);
		tmp_Face.ipoints.push_back(centerPoint);
		for(auto& k : subOutEdgePoints[leftEdge]){
			tmp_Face.ipoints.push_back(k);
		}
		
	}
	
};



void MASCH_Poly_AMR_Builder::extractSubInternalEdgesPoints(
	vector<int>& vertexCenterPoints,
	int& centerPoint,
	vector<vector<int>>& subIntEdgePoints
) {

	int verCentPSize = vertexCenterPoints.size();
	// cout << verCentPSize << endl;
	subIntEdgePoints.resize(verCentPSize);
	for(int j=0; j<verCentPSize; ++j){
		subIntEdgePoints[j].push_back(vertexCenterPoints[j]);
		subIntEdgePoints[j].push_back(centerPoint);
	}

};




void MASCH_Poly_AMR_Builder::extractCellVertexOrder(
	vector<int>& cellVertexPoints,
	map<int,int>& cellVertexOrder
) {

	for(int j=0; j<cellVertexPoints.size(); ++j){
		cellVertexOrder.insert(make_pair(cellVertexPoints[j],j));
	}
	// int tmpNum = 0;
	// for(int j=0; j<cell.ipoints.size(); ++j){
		// int vertex = cell.ipoints[j];
		// if(mesh.points[vertex].level > cellLevel) continue;
		// cellVertexOrder.insert(make_pair(vertex,tmpNum));
		// ++tmpNum;
	// }
};


void MASCH_Poly_AMR_Builder::extractInternalFacesVertexs(
	MASCH_Cell& cell,
	map<int,int>& cellInternalFaces,
	vector<int>& groupChildFaces_id,
	vector<groupMesh_Refine>& groupChildFaces,
	groupMesh_Refine& groupChildFace,
	int& cellCenterPoint,
	vector<vector<int>>& intFacesVertexs
) {

	intFacesVertexs.resize(cellInternalFaces.size());
	
	for(auto& j : cell.ifaces){
		auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
		
		int faceCenterPoint = groupOuterFace.centerPoint;
		
		int verCentPSize = groupOuterFace.vertexCenterPoints.size();
		for(int k=0; k<verCentPSize; ++k){
			int verCentP = groupOuterFace.vertexCenterPoints[k];
			int iFace = cellInternalFaces[verCentP];
			auto& face = groupChildFace.faces[iFace];
			// int facePointsSize = face.ipoints.size();
			
			if(
			(std::find(face.ipoints.begin(),face.ipoints.end(),cellCenterPoint) 
				== face.ipoints.end()) 
			){
				
				// cout << "POINTS0 : " << face.ipoints.size() << endl;
				// cout << groupOuterFace.subIntEdgePoints.size() << endl;
				for(auto& l : groupOuterFace.subIntEdgePoints[k]){
					face.ipoints.push_back(l);
				}
				face.ipoints.push_back(cellCenterPoint);
				// cout << "POINTS1 : " << face.ipoints.size() << endl;
			}
			else if(
			(std::find(face.ipoints.begin(),face.ipoints.end(),faceCenterPoint) 
				== face.ipoints.end())
			){
				// cout << "POINTS0 : " << face.ipoints.size() << endl;
				int subEdgeSize = groupOuterFace.subIntEdgePoints[k].size();
				for(int l=subEdgeSize-1; l>=1; --l){
					face.ipoints.push_back(groupOuterFace.subIntEdgePoints[k][l]);
				}
				
				// cout << "POINTS1 : " << face.ipoints.size() << endl;
				
			}
			
			
			if(verCentPSize==2){
				int verPoint = groupOuterFace.vertexPoints[0];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint);
				}
				
			}
			else{
				int verPoint0 = groupOuterFace.vertexPoints[k];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint0) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint0);
				}
				
				int verPoint1 = groupOuterFace.vertexPoints[ ( (k==verCentPSize-1) ? 0 : k+1 ) ];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint1) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint1);
				}
				
			}
		}
	}
};



void MASCH_Poly_AMR_Builder::addInternalFacesiLiR(
	MASCH_Mesh& mesh, 
	int& cellLevel,
	int& totalCellNum,
	map<int,int>& cellVertexOrder,
	vector<vector<int>>& intFacesVertexs,
	vector<faces_refind>& faces
) {

	
	for(int j=0; j<faces.size(); ++j){
		// cout << intFacesVertexs[j].size() << endl;
		int iiL = intFacesVertexs[j][0];
		int iiR = intFacesVertexs[j][1];
		
		vector<int> tmpVertexPoints;
		for(auto& j : faces[j].ipoints){
			if(mesh.points[j].level > cellLevel+1) continue;
			tmpVertexPoints.push_back(j);
		}
		
		vector<double> unitNormals(3,0.0);
		double x1 = mesh.points[tmpVertexPoints[0]].x;
		double x2 = mesh.points[tmpVertexPoints[1]].x;
		double x3 = mesh.points[tmpVertexPoints[2]].x;
		double y1 = mesh.points[tmpVertexPoints[0]].y;
		double y2 = mesh.points[tmpVertexPoints[1]].y;
		double y3 = mesh.points[tmpVertexPoints[2]].y;
		double z1 = mesh.points[tmpVertexPoints[0]].z;
		double z2 = mesh.points[tmpVertexPoints[1]].z;
		double z3 = mesh.points[tmpVertexPoints[2]].z;
		double v1[3] = {x2-x1,y2-y1,z2-z1};
		double v2[3] = {x3-x1,y3-y1,z3-z1};
		
		unitNormals[0] = v1[1] * v2[2] - v1[2] * v2[1];
		unitNormals[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
		unitNormals[2] = v1[0] * v2[1] - v1[1] * v2[0];
		
		double ownVecX = mesh.points[iiL].x-mesh.points[iiR].x;
		double ownVecY = mesh.points[iiL].y-mesh.points[iiR].y;
		double ownVecZ = mesh.points[iiL].z-mesh.points[iiR].z;
		
		double cosTheta = 
			unitNormals[0]*ownVecX + 
			unitNormals[1]*ownVecY + 
			unitNormals[2]*ownVecZ;
		
		// cout << cosTheta << endl;
		
		if(cosTheta>=0.0){
			faces[j].iR = totalCellNum + cellVertexOrder[iiL];
			faces[j].iL = totalCellNum + cellVertexOrder[iiR];
		}
		else{
			faces[j].iL = totalCellNum + cellVertexOrder[iiL];
			faces[j].iR = totalCellNum + cellVertexOrder[iiR];
		}
		
		
	}
	
			
};


void MASCH_Poly_AMR_Builder::addOuterFacesiLiR(
	MASCH_Mesh& mesh, 
	MASCH_Cell& cell,
	int& i,
	int& totalCellNum,
	map<int,int>& cellVertexOrder,
	vector<int>& groupChildFaces_id,
	vector<groupMesh_Refine>& groupChildFaces
) {
	for(auto& j : cell.ifaces){
		auto& face = mesh.faces[j];
		if(groupChildFaces_id[j]==-1) cout << "NONONONONON" << endl;
		auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
		
		int verPSize = groupOuterFace.vertexPoints.size();
		for(int k=0; k<verPSize; ++k){
			int verP = groupOuterFace.vertexPoints[k];
			auto& subFace = groupOuterFace.faces[k];
			if(face.iL == i){
				subFace.iL = totalCellNum + cellVertexOrder[verP];
				cout << " iL = " << subFace.iL << endl;
			}
			else{
				subFace.iR = totalCellNum + cellVertexOrder[verP];
				cout << " iR = " << subFace.iR << endl;
			}
		}
	}

};


void MASCH_Poly_AMR_Builder::reorderOuterFacesiLiR(
	MASCH_Mesh& mesh, 
	MASCH_Cell& cell,
	int& i,
	int& totalCellNum,
	vector<int>& groupChildFaces_id,
	vector<bool>& boolInputFacesiL,
	vector<bool>& boolInputFacesiR, 
	vector<groupMesh_Refine>& groupChildFaces
) {

	for(auto& j : cell.ifaces){
		auto& face = mesh.faces[j];
		
		if(groupChildFaces_id[j]==-1){
			if(face.iL == i && boolInputFacesiL[j]==false){
				face.iL = totalCellNum;
				boolInputFacesiL[j] = true;
				// cout << " iL = " << face.iL << " " << face.level << endl;
			}
			else if(face.iR == i && boolInputFacesiR[j]==false){
				face.iR = totalCellNum;
				boolInputFacesiR[j] = true;
				// cout << " iR = " << face.iR << " " << face.level << endl;
			}
			else{
				cout << "ERROR 25" << endl;
				cout << face.iL << " " <<  face.iR << " " << i << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}
		else{
			auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
			for(auto& subFace : groupOuterFace.faces){
				if(face.iL == i && subFace.booliL==false){
					subFace.iL = totalCellNum;
					subFace.booliL = true;
					// cout << " iL = " << subFace.iL << endl;
				}
				else if(face.iR == i && subFace.booliR==false){
					subFace.iR = totalCellNum;
					subFace.booliR = true;
					// cout << " iR = " << subFace.iR << endl;
				}
				else{
					cout << "ERROR 35" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
		}
	}
};

