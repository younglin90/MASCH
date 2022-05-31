
#include "./mesh.h"
#include "./polyAMR.h"
// #include "geometric.h"
#include "./mpi.h"




void MASCH_Poly_AMR_Builder::polyUnrefine(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	int maxLevel_AMR, 
	vector<vector<double>> indicatorCriterion,
	vector<vector<double>>& indicatorValues,
	vector<vector<int>>& child_org_cell_id_of_new,
	vector<bool>& boolCellPreserved,
	vector<bool>& boolCellRefine,
	vector<bool>& boolCellUnrefine,
	int iter){
		
		
		
	int nBuffers = 3;
	
	
		

	MASCH_MPI mpi;

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int beforeCellSize = mesh.cells.size();
	int beforeFaceSize = mesh.faces.size();
	int beforePointSize = mesh.points.size();
	int afterCellSize = 0;
	int afterFaceSize = 0;
	int afterPointSize = 0;
	
	
	
	
	int proc_num = 0;
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute AMR - Unrefine ";
	}
	

	// SEMO_Mesh_Geometric geometric;
	
	//====================================================
	// 셀 Unrefine : Unrefine 되는 셀 & 면 조사
	
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	// cout << indicatorValues.size() << " " <<
	// indicatorCriterion.size() << " " <<
	// indicatorValues[0].size() << " " <<
	// indicatorCriterion[0].size() << " " <<
	// mesh.cells.size() <<
	// endl;
	
	
	
	
	
	
	
	
	
	
	
	

	// vector<bool> boolBufferPreserved(mesh.cells.size(),false);
	// {
		// vector<int> cLevel_recv;
		// this->mpiLevels(mesh, cLevel_recv);
		
		// for(int iLevel=maxLevel_AMR-1; iLevel>=1; --iLevel){
			// for(int i=0; i<mesh.nInternalFaces; ++i){
				// auto& face = mesh.faces[i];
				// int iL = face.iL;
				// int iR = face.iR;
				// auto& cellL = mesh.cells[iL];
				// auto& cellR = mesh.cells[iR];
				// if(cellL.level==iLevel && cellL.level < cellR.level){
					// boolBufferPreserved[iL]=true;
				// }
				// if(cellR.level==iLevel && cellL.level > cellR.level){
					// boolBufferPreserved[iR]=true;
				// }
			// }
			// int ip=0;
			// for(auto& boundary : mesh.boundaries){
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				// if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
						// auto& cellL = mesh.cells[iL];
						// if(cellL.level==iLevel && cellL.level < cLevel_recv[ip]){
							// boolBufferPreserved[iL]=true;
						// }
						// ++ip;
					// }
				// }
			// }	
			// for(int iBuffer=0; iBuffer<nBuffers; ++iBuffer){
				// vector<int> boolBufferPreserved_recv;
				// if(size>1){
					// vector<int> tmp_send;
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){
							// if(boolBufferPreserved[face.iL]==true){
								// tmp_send.push_back(1);
							// }
							// else{
								// tmp_send.push_back(0);
							// }
						// }
					// }
					// boolBufferPreserved_recv.resize(tmp_send.size(),0);
					// MPI_Alltoallv( tmp_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
								   // boolBufferPreserved_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
								   // MPI_COMM_WORLD);
							   
				// }
				// for(int i=0; i<mesh.nInternalFaces; ++i){
					// auto& face = mesh.faces[i];
					// int iL = face.iL;
					// int iR = face.iR;
					// auto& cellL = mesh.cells[iL];
					// auto& cellR = mesh.cells[iR];
					// if(cellL.level==iLevel && 
					// cellL.level <= cellR.level &&
					// boolBufferPreserved[iR]==true){
						// boolBufferPreserved[iL]=true;
					// }
					// if(cellR.level==iLevel && 
					// cellR.level <= cellL.level &&
					// boolBufferPreserved[iL]==true){
						// boolBufferPreserved[iR]=true;
					// }
				// }
				// int ip=0;
				// for(auto& boundary : mesh.boundaries){
					// int str = boundary.startFace;
					// int end = str + boundary.nFaces;
					// if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
						// for(int i=str; i<end; ++i){
							// auto& face = mesh.faces[i];
							// int iL = face.iL;
							// auto& cellL = mesh.cells[iL];
							// if(cellL.level==iLevel && 
							// cellL.level <= cLevel_recv[ip] &&
							// boolBufferPreserved_recv[ip]==1){
								// boolBufferPreserved[iL]=true;
							// }
							// ++ip;
						// }
					// }
				// }	
			// }
		// }
	// }		
					
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	// vector<bool> boolCellUnrefine(mesh.cells.size(),true);
	for(int i=0; i<mesh.cells.size(); ++i){
		
		// for(int indi=0; indi<indicatorCriterion.size(); ++indi)
		// {
			// for(int level=0; level<indicatorCriterion.at(indi).size(); ++level)
			// {
				// double indicatorRefine_AMR = indicatorCriterion.at(indi).at(level);
				// if( mesh.cells[i].level > level ){
					// if( indicatorValues.at(indi).at(i) <= indicatorRefine_AMR ){
						
					// }
					// else{
						// boolCellUnrefine[i] = false;
					// }
				// }				
			// }
		// }
		// // for(int level=0; level<controls.indicatorCriterion.size(); ++level){
			// // double indicatorRefine_AMR = controls.indicatorCriterion[level];
			// // if( mesh.cells[i].level > level ){
				// // if( mesh.cells[i].var[controls.indicatorAMR[0]] <= indicatorRefine_AMR ){
					// // boolCellUnrefine[i] = true;
				// // }
			// // }					
		// // }
		// // if( 
		// // (rank==0 && distr(eng) > 100.0) ||
		// // (rank==1 && distr(eng) > 0.0) ||
		// // (rank==2 && distr(eng) > 0.0) ||
		// // (rank==3 && distr(eng) > 0.0) 
		// // ){
			// // boolCellUnrefine[i] = true;
		// // }
        
                // boolCellUnrefine[i] = false;
        // if(distr(eng) > 0.4){
            // boolCellUnrefine[i] = true;
        // }                
		
		// 만약 셀의 레벨이 0 이면, false
		if(mesh.cells[i].level <= 0) boolCellUnrefine[i] = false;
		if(mesh.cells[i].level > maxLevel_AMR) boolCellUnrefine[i] = true;
		// if(boolCellPreserved[i] == true) boolCellUnrefine[i] = false;
		
		
		// if(boolBufferPreserved[i] == true) boolCellUnrefine[i] = false;
		
		
	}
	
	
	
	
	
	
	
	
	
	
	

	// // 대각선 방향에서, 레벨 차이 2 이상이면 리파인
	// {
		// // processor faces
		// vector<int> recv_value;
		// if(size>1){

			// vector<int> send_value;
			// send_value.reserve(mesh.send_StencilCellsId.size());
			// for(auto& icell : mesh.send_StencilCellsId){
				// send_value.push_back(mesh.cells[icell].level);
			// }
			// recv_value.resize(mesh.recv_displsStencilCells[size]);
			// MPI_Alltoallv( send_value.data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_INT, 
						   // recv_value.data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_INT, 
						   // MPI_COMM_WORLD);
			
		// }
		// vector<int> tmp_maxLevel(mesh.cells.size());
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = mesh.cells[i];
			// int maxInd = -100;
			// for(auto& icell : cell.iStencils){
				// maxInd = max(mesh.cells[icell].level,maxInd);
			// }
			// for(auto& icell : cell.recv_iStencils){
				// maxInd = max(recv_value[icell],maxInd);
			// }
			// tmp_maxLevel[i] = maxInd;
		// }
		// int inp_size = indicatorValues.size();
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			// int my_level = cell.level;
			// if(my_level <= 0) continue;
			// if(boolCellPreserved[i] == true) continue;
			// if(boolBufferPreserved[i] == true) continue;
			// if(my_level<tmp_maxLevel[i]-2){
				// boolCellUnrefine[i] = false;
			// }
		// }
	// }
	
	
	
	
	
	
	
	
	
	
	// this->bufferLayerUnrefine(mesh, boolCellUnrefine, nBuffers);
	
	
	
	
	
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	vector<int> cLevel_recv;
	vector<int> cUnrefine_recv;
	this->mpiLevelRefine(mesh, boolCellUnrefine, cLevel_recv, cUnrefine_recv);
	
	this->restrictCellUnrefine(mesh, boolCellUnrefine, cLevel_recv, cUnrefine_recv);
	
	this->mpiRefines(mesh, boolCellUnrefine, cUnrefine_recv);
	
	
	
	//====================================================
	// 오리지널 그룹 
	vector<vector<int>> groupCellListsCanUnrefine;
	
	extractGroupCellListsCanUnrefine(mesh, groupCellListsCanUnrefine);
	
	//====================================================
	
	vector<groupCells_Unrefine> groupCellsUnrefine;
	extractGroupUnrefineCells(mesh, boolCellUnrefine, groupCellListsCanUnrefine, 
		groupCellsUnrefine);
	
	
	std::fill(boolCellUnrefine.begin(), boolCellUnrefine.end(), false);
	vector<int> groupCells_id(mesh.cells.size(),-1);
	for(int i=0; i<groupCellsUnrefine.size(); ++i){
		for(auto& j : groupCellsUnrefine[i].ichild){
			boolCellUnrefine[j] = true;
			groupCells_id[j] = i;
		}
	}
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(boolCellUnrefine[i]!=true && mesh.cells[i].level > 0){
			// cout << "EEEE" << endl;
		// }
	// }
	
	this->mpiRefines(mesh, boolCellUnrefine, cUnrefine_recv);
	
	//====================================================
	// 셀 넘버링
	
	
	vector<int> newCellsNumber(mesh.cells.size(),-1);
	vector<int> cellsLevel(mesh.cells.size(),-1);
	vector<int> cellsGroup(mesh.cells.size(),-1);
	vector<vector<int>> cells_iparcels(mesh.cells.size());
	int newCellNum = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		if(groupCells_id[i] == -1){
			// child_org_cell_id_of_new.push_back(vector<int>());
			// child_org_cell_id_of_new.back().push_back(i);
			
			cellsLevel[newCellNum] = mesh.cells[i].level;
			// if(boolCanNotUnrefineCells[i]==true) cellsLevel[newCellNum] = -1;
			cellsGroup[newCellNum] = mesh.cells[i].group;
			
			cells_iparcels[newCellNum] = mesh.cells[i].iparcels;
			
			newCellsNumber[i] = newCellNum;
			
			++newCellNum;
		}
		else{
			// child_org_cell_id_of_new.push_back(vector<int>());
			
			cellsLevel[newCellNum] = mesh.cells[i].level-1;
			// if(boolCanNotUnrefineCells[i]==true) cellsLevel[newCellNum] = -1;
			cellsGroup[newCellNum] = mesh.cells[i].group;
			
			vector<int> tmp_iparcels;
			int tmp_id = groupCells_id[i];
			for(auto& j : groupCellsUnrefine[tmp_id].ichild){
				// child_org_cell_id_of_new.back().push_back(i);
				
				newCellsNumber[j] = newCellNum;
				
				for(auto& item : mesh.cells[i].iparcels){
					tmp_iparcels.push_back(item);
				}
				
				++i;
			}
			--i;
			
			cells_iparcels[newCellNum] = tmp_iparcels;
			
			++newCellNum;
		}
	}
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	int totalCellNum = newCellNum;
	
	
	child_org_cell_id_of_new.clear();
	child_org_cell_id_of_new.resize(totalCellNum);
	boolCellPreserved.resize(totalCellNum);
	// boolCellUnrefine.resize(totalCellNum);
	for(int i=0, iter=0; i<mesh.cells.size(); ++i){
		if(groupCells_id[i] == -1){
			child_org_cell_id_of_new[iter].push_back(i);
			boolCellPreserved[iter]=false;
			// boolCellUnrefine[iter]=false;
			++iter;
		}
		else{
			int tmp_id = groupCells_id[i];
			for(auto& j : groupCellsUnrefine[tmp_id].ichild){
				child_org_cell_id_of_new[iter].push_back(i);
				boolCellPreserved[iter]=true;
				// boolCellUnrefine[iter]=true;
				++i;
			}
			--i;
			++iter;
		}
	}
	
	// // 그룹 리넘버링
	// {
		// int total_nGroup = 0;
		// int tmp_cell_size = mesh.cells.size();
		// vector<int> str_icells(size+1,0);
		// for(int i=0; i<tmp_cell_size; ++i){
			// auto& cell0 = mesh.cells[i];
			// int group0 = cell0.group;
			// cell0.group = total_nGroup;
			// while(1){
				// ++i;
				// if(i==tmp_cell_size) break;
				// auto& cell = mesh.cells[i];
				// int group = cell.group;
				// if(group0!=group) break;
				// cell.group = total_nGroup;
			// }
			// --i;
			// ++total_nGroup;
		// }
		
		// if(size>1){
			// vector<int> recv_ncells(size);
			// MPI_Allgather(&total_nGroup, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
			// for(int ip=0; ip<size; ++ip){
				// str_icells[ip+1] = str_icells[ip] + recv_ncells[ip];
			// }
		// }
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// mesh.cells[i].group += str_icells[rank];
		// }
	// }
	
	
	
	
	//====================================================
	// 포인트 삭제 및 넘버링

	vector<bool> boolDeletePoints(mesh.points.size(),true);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		int cellLevel = cell.level;
		
		if(cellLevel==-1) cellLevel=0;
		
		if(boolCellUnrefine[i] == true){
			for(auto& j : cell.ipoints){
				int pointLevel = mesh.points[j].level;
				if(pointLevel <= cellLevel-1){
					boolDeletePoints[j] = false;
				}
			}
		}
		else{
			for(auto& j : cell.ipoints){
				int pointLevel = mesh.points[j].level;
				if(pointLevel <= cellLevel)
					boolDeletePoints[j] = false;
			}
		}
	}
	
	// proc
	proc_num = 0;
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType() != MASCH_Face_Types::PROCESSOR) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			auto& face = mesh.faces[i];
			int faceLevel = face.level;
			int ownLevel = mesh.cells[face.iL].level;
			int ngbLevel = cLevel_recv[proc_num];
			bool ownUnrefine = boolCellUnrefine[face.iL];
			bool ngbUnrefine = cUnrefine_recv[proc_num];
			
			if(faceLevel==-1) faceLevel=0;
			if(ownLevel==-1) ownLevel=0;
			if(ngbLevel==-1) ngbLevel=0;
			
			if(
			(ownLevel<ngbLevel && ngbUnrefine==false) ||
			(ownLevel==ngbLevel && ownUnrefine==true && ngbUnrefine==false) 
			){
				for(auto& j : face.ipoints){
					int pointLevel = mesh.points[j].level;
					if(pointLevel <= faceLevel){
						boolDeletePoints[j] = false;
					}
				}
			}
			++proc_num;
		}
	}

	
	// *******************************
	// connPoints
	{
		vector<vector<int>> send_connPoint_id(size);
		vector<vector<int>> send_boolDeletePoints(size);
		for(int i=0, SIZE=mesh.points.size(); i<SIZE; ++i){
			auto& point = mesh.points[i];
			if(boolDeletePoints[i]==true){
				for(auto& [proc, id] : point.connPoints){
					send_connPoint_id[proc].push_back(id);
					send_boolDeletePoints[proc].push_back(1);
				}
			}
			else{
				for(auto& [proc, id] : point.connPoints){
					send_connPoint_id[proc].push_back(id);
					send_boolDeletePoints[proc].push_back(0);
				}
			}
		}
		vector<vector<int>> recv_connPoint_id;
		mpi.Alltoallv(send_connPoint_id, recv_connPoint_id);
		vector<vector<int>> recv_boolDeletePoints;
		mpi.Alltoallv(send_boolDeletePoints, recv_boolDeletePoints);
		for(int ip=0; ip<size; ++ip){
			int tmp_size = recv_connPoint_id[ip].size();
			for(int i=0; i<tmp_size; ++i){
				int id = recv_connPoint_id[ip][i];
				int recv_boolDelPoint = recv_boolDeletePoints[ip][i];
				if(recv_boolDelPoint==0){
					boolDeletePoints[id] = false;
				}
			}
		}
	}
	// *******************************
	
	
	
	vector<int> newPointsNumber(mesh.points.size(),-1);
	int newPointNum = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		if(boolDeletePoints[i] == false){
			newPointsNumber[i] = newPointNum;
			++newPointNum;
		}
	}
	
	
	
	// *******************************
	// connPoints 번호 재설정
	{
		vector<vector<int>> send_connPoint_id(size);
		vector<vector<int>> send_connPoint_new_id(size);
		for(int i=0, SIZE=mesh.points.size(); i<SIZE; ++i){
			auto& point = mesh.points[i];
			for(auto& [proc, id] : point.connPoints){
				send_connPoint_id[proc].push_back(id);
				send_connPoint_new_id[proc].push_back(newPointsNumber[i]);
			}
			point.connPoints.clear();
		}
		vector<vector<int>> recv_connPoint_id;
		mpi.Alltoallv(send_connPoint_id, recv_connPoint_id);
		vector<vector<int>> recv_connPoint_new_id;
		mpi.Alltoallv(send_connPoint_new_id, recv_connPoint_new_id);
		for(int ip=0; ip<size; ++ip){
			int tmp_size = recv_connPoint_id[ip].size();
			for(int i=0; i<tmp_size; ++i){
				int id = recv_connPoint_id[ip][i];
				int new_id = recv_connPoint_new_id[ip][i];
				auto& point = mesh.points[id];
				int lll = 0;
				for(auto& [proc, recvId] : point.connPoints){
					if(proc==ip && recvId==new_id) ++lll;
				}
				if(lll==0){
					point.connPoints.push_back(make_pair(ip,new_id));
				}
			}
		}
		
	}
	// *******************************
	
	
	//====================================================
	// Unrefine : 면 합침 , 면 - 포인트 연결
	
	int proc_total_num = 0;
	
	vector<bool> boolIntFaces_Deleted(mesh.faces.size(),false);
	vector<int> groupFaces_id(mesh.faces.size(),-1);
	vector<faces_Unrefind> groupOutFaces;
	for(int i=0; i<groupCellsUnrefine.size(); ++i){
		auto& group = groupCellsUnrefine[i];
		
		auto& target_Cell = mesh.cells[group.ichild[0]];
		group.cellCenterPoint = -1;
		for(auto& targetPoint : target_Cell.ipoints){
			bool boolCellCenterPoint = true;
			for(int j=1; j<group.ichild.size(); ++j){
				auto& cell = mesh.cells[group.ichild[j]];
				if(
				std::find(cell.ipoints.begin(),cell.ipoints.end(),targetPoint)
				== cell.ipoints.end()
				){
					boolCellCenterPoint = false;
					break;
				}
			}
			if(boolCellCenterPoint==true){
				group.cellCenterPoint = targetPoint;
				break;
			}
		}
		
		// search cell internal & outer faces
		int cellCenterPoint = group.cellCenterPoint;
		vector<int> outChildFaces;
		set<int> intChildFaces;
		for(auto& j : group.ichild){
			auto& cell = mesh.cells[j];
			for(auto& k : cell.ifaces){
				auto& face = mesh.faces[k];
				if(
				std::find(face.ipoints.begin(),face.ipoints.end(),cellCenterPoint)
				== face.ipoints.end()
				){
					outChildFaces.push_back(k);
				}
				else{
					boolIntFaces_Deleted[k] = true;
					groupFaces_id[k] = i;
					intChildFaces.insert(k);
				}
			}
		}
		
		group.intChildFaces.assign(intChildFaces.begin(),intChildFaces.end());
		
		// search face center points
		set<int> tmpFaceCenterPoints;
		for(int j=0; j<group.intChildFaces.size(); ++j){
			auto& face0 = mesh.faces[group.intChildFaces[j]];
			for(int k=j+1; k<group.intChildFaces.size(); ++k){
				auto& face1 = mesh.faces[group.intChildFaces[k]];
				for(auto& p0 : face0.ipoints){
					for(auto& p1 : face1.ipoints){
						if(p0==p1){
							tmpFaceCenterPoints.insert(p0);
						}
					}
				}
			}
		}
		
		if(tmpFaceCenterPoints.find(cellCenterPoint)==tmpFaceCenterPoints.end()){
			cout << "| #WARNING : NO searching cellCenterPoint " << rank << " " << group.cellCenterPoint << endl;
		}
		
		tmpFaceCenterPoints.erase(tmpFaceCenterPoints.find(cellCenterPoint));
		group.faceCenterPoints.assign(tmpFaceCenterPoints.begin(),tmpFaceCenterPoints.end());
		
		// groupping outer faces, save points, own, ngb
		for(auto& centP : group.faceCenterPoints){
			vector<int> tmpChildFaces;
			for(auto& j : outChildFaces){
				auto& face = mesh.faces[j];
				if(
				std::find(face.ipoints.begin(),face.ipoints.end(),centP)
				!= face.ipoints.end()
				){
					tmpChildFaces.push_back(j);
				}
			}
			
			int face0 = tmpChildFaces[0];
			if(groupFaces_id[face0]==-1){
				for(auto& j : tmpChildFaces){
					groupFaces_id[j] = groupOutFaces.size();
				}
				group.groupOutFaces_id.push_back(groupOutFaces.size());
				groupOutFaces.push_back(faces_Unrefind());
				groupOutFaces.back().ichild = tmpChildFaces;
				
				// points
				std::sort(groupOutFaces.back().ichild.begin(), 
						  groupOutFaces.back().ichild.end());
						  
				vector<vector<int>> vertexLists;
				for(auto& j : groupOutFaces.back().ichild){
					auto& face = mesh.faces[j];
					int faceLevel = face.level;
					vertexLists.push_back(vector<int>());
					for(auto& p0 : face.ipoints){
						auto& point = mesh.points[p0];
						int pointLevel = point.level;
						if(pointLevel<=faceLevel){
							vertexLists.back().push_back(p0);
						}
					}
				}
				
				
				vector<int> groupFacePoints0;
				int tmptmp = 0;
				for(auto& j : groupOutFaces.back().ichild){
					auto& face = mesh.faces[j];
					auto& points = face.ipoints;
					
					if(tmptmp!=0){
						
						auto iStr = std::find(points.begin(),points.end(),vertexLists[tmptmp][3]);
						
						groupFacePoints0.insert( groupFacePoints0.end(), iStr, points.end() );
						
					}
					
					auto iEnd = std::find(points.begin(),points.end(),vertexLists[tmptmp][1]);
					groupFacePoints0.insert( groupFacePoints0.end(), points.begin(), iEnd );
					
					
					++tmptmp;
				}
				auto& faceZero = mesh.faces[groupOutFaces.back().ichild[0]];
				auto& pointsZero = faceZero.ipoints;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists[0][3]);
				groupFacePoints0.insert( groupFacePoints0.end(), iStr0, pointsZero.end() );
				
				vector<int> groupFacePoints;
				for(auto& j : groupFacePoints0){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						groupFacePoints.push_back(newPN);
					}
				}
							
				groupOutFaces.back().ipoints = groupFacePoints;
				
				// own, ngb
				groupOutFaces.back().iL = newCellsNumber[mesh.faces[face0].iL];
				if(mesh.faces[face0].iR==-1){
					groupOutFaces.back().iR = -1;
				}
				else{
					groupOutFaces.back().iR = newCellsNumber[mesh.faces[face0].iR];
				}
				
				
			}
			else{
				int group_id = groupFaces_id[face0];
				group.groupOutFaces_id.push_back(group_id);
			}
			
		}
	}
	
	

	// proc
	vector<bool> boolProcFaces_Combine(proc_num,false);
	vector<int> groupProcFaces_id(proc_num,-1);
	vector<faces_Unrefind> groupProcFaces;
	
	
	proc_num = 0;
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType() != MASCH_Face_Types::PROCESSOR) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			auto& face = mesh.faces[i];
			int faceLevel = face.level;
			int ownLevel = mesh.cells[face.iL].level;
			int ngbLevel = cLevel_recv[proc_num];
			bool ownUnrefine = boolCellUnrefine[face.iL];
			bool ngbUnrefine = cUnrefine_recv[proc_num];
			if(
			(ownLevel<ngbLevel && ngbUnrefine==true)
			){
				
				if(boolCellUnrefine[face.iL]==true){
					cout << "NONONONONONONONOON1" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				// face center point
				auto& face1 = mesh.faces[i+1];
				auto& face2 = mesh.faces[i+2];
				int centerP = -1;
				for(auto& p0 : face.ipoints){
					if(std::find(face1.ipoints.begin(),face1.ipoints.end(),p0)==face1.ipoints.end()){
						continue;
					}
					if(std::find(face2.ipoints.begin(),face2.ipoints.end(),p0)==face2.ipoints.end()){
						continue;
					}
					centerP = p0;
					break;
				}
				
				if(centerP==-1){
					cout << "NONONONONONONONOON2" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				vector<int> tmpProcFace;
				int ownTarget = face.iL;
				for(int j=i; j<end; ++j){
					auto& face_tmp = mesh.faces[j];
					if(std::find(face_tmp.ipoints.begin(),face_tmp.ipoints.end(),centerP)!=face_tmp.ipoints.end()){
						tmpProcFace.push_back(j);
					}
					else{
						break;
					}
				}
				
				int groupNum = groupProcFaces.size();
				groupProcFaces.push_back(faces_Unrefind());
				for(auto& j : tmpProcFace){
					boolProcFaces_Combine[proc_num] = true;
					groupProcFaces_id[proc_num] = groupNum;
					groupProcFaces.back().ichild.push_back(j);
					++i;
					++proc_num;
				}
				--i;
				--proc_num;
				
				
				// proc 면 포인트 순서 바꾸기
				vector<int> tmpPoints;
				if(boundary.rightProcNo > rank){
					
					std::reverse(tmpProcFace.begin()+1, tmpProcFace.end());
						
				}
				

				vector<int> vertexLists0;
				int tmptmp = 0;
				for(auto& j : tmpProcFace){
					auto& face_tmp = mesh.faces[j];
					int faceLevel_tmp = face_tmp.level;
					auto& points_tmp = face_tmp.ipoints;
					
					vector<int> vertexLists;
					for(auto& p0 : points_tmp){
						if(mesh.points[p0].level<=faceLevel_tmp){
							vertexLists.push_back(p0);
						}
					}
					if(tmptmp==0){
						vertexLists0 = vertexLists;
					}
			
					if(tmptmp!=0){
						auto iStr = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[3]);
						tmpPoints.insert( tmpPoints.end(), iStr, points_tmp.end() );
					}
					auto iEnd = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[1]);
					tmpPoints.insert( tmpPoints.end(), points_tmp.begin(), iEnd );
					++tmptmp;
				}
				auto& faceZero = mesh.faces[tmpProcFace[0]];
				auto& pointsZero = faceZero.ipoints;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists0[3]);
				tmpPoints.insert( tmpPoints.end(), iStr0, pointsZero.end() );
						
				vector<int> groupFacePoints;
				for(auto& j : tmpPoints){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						groupFacePoints.push_back(newPN);
					}
				}
				
				groupProcFaces.back().ipoints = groupFacePoints;
				
				
			}
			++proc_num;
		}
	}
	
	
	//====================================================
	// proc 면 포인트 순서 바꾸기 & 재부여
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType() != MASCH_Face_Types::PROCESSOR) continue;
		if(boundary.rightProcNo <= rank) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			if(groupFaces_id[i] != -1){
				
				auto& group_face = groupOutFaces[groupFaces_id[i]];
				
				
				vector<int> tmpProcFace = group_face.ichild;
				
				std::reverse(tmpProcFace.begin()+1, tmpProcFace.end());
				
				
				vector<int> tmpPoints;
				vector<int> vertexLists0;
				int tmptmp = 0;
				for(auto& j : tmpProcFace){
					auto& face_tmp = mesh.faces[j];
					int faceLevel_tmp = face_tmp.level;
					auto& points_tmp = face_tmp.ipoints;
					
					vector<int> vertexLists;
					for(auto& p0 : points_tmp){
						if(mesh.points[p0].level<=faceLevel_tmp){
							vertexLists.push_back(p0);
						}
					}
					if(tmptmp==0){
						vertexLists0 = vertexLists;
					}
			
					if(tmptmp!=0){
						auto iStr = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[3]);
						tmpPoints.insert( tmpPoints.end(), iStr, points_tmp.end() );
					}
					auto iEnd = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[1]);
					tmpPoints.insert( tmpPoints.end(), points_tmp.begin(), iEnd );
					++tmptmp;
				}
				auto& faceZero = mesh.faces[tmpProcFace[0]];
				auto& pointsZero = faceZero.ipoints;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists0[3]);
				tmpPoints.insert( tmpPoints.end(), iStr0, pointsZero.end() );
						
						
				auto& points = group_face.ipoints;
				points.clear();
				for(auto& j : tmpPoints){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						points.push_back(newPN);
					}
				}
				
				
				i += group_face.ichild.size();
				--i;
			}
		}
	}
	
	proc_total_num = 0;
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == MASCH_Face_Types::PROCESSOR){
			int faceLevel = face.level;
			
			bool ownUnrefine = boolCellUnrefine[face.iL];
			bool ngbUnrefine = false;
			if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
			int ownLevel = mesh.cells[face.iL].level;
			int ngbLevel = cLevel_recv[proc_num];
			
			if(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel){
				++proc_total_num;
			}
			else if(ownUnrefine==true && ownLevel>ngbLevel){
				++proc_total_num;
			}
			else if(ngbUnrefine==true && ownLevel<ngbLevel){
				
				if(boolProcFaces_Combine[proc_num]==false){
					cout << "AAAAAAAAAAAAAA" << endl;
				}
				++proc_total_num;
			}
			else{
				
			}
			++proc_num;
			
		}
		
	}
	
	
	//====================================================
	// 포인트 삭제
	int numN = 0;
	mesh.points.erase( std::remove_if( mesh.points.begin(), mesh.points.end(), 
		[&boolDeletePoints, &numN](MASCH_Point const& v) { 
		return boolDeletePoints[numN++]; 
		}), mesh.points.end());
		
		
	//====================================================
	// 면 재정립
	
	proc_total_num=0;
	int nBC = 0;
	int saveI;
	numN = 0;
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		int orgnBC = nBC;
		for(int j=orgnBC; j<mesh.boundaries.size(); ++j){
			if(
			numN==mesh.boundaries[j].startFace ||
			mesh.boundaries[j].startFace == 0 ||
			mesh.boundaries[j].nFaces == 0){
				mesh.boundaries[j].startFace = i;
				++nBC;
			}
			else{
				break;
			}
		}
		// internal faces
		if(boolIntFaces_Deleted[numN]==true){
			auto& group = groupCellsUnrefine[ groupFaces_id[numN] ];
			numN += group.intChildFaces.size();
			--i;
		}
		// outer faces
		else{
			if(mesh.faces[numN].getType() != MASCH_Face_Types::PROCESSOR){
				// let it be
				if(groupFaces_id[numN]==-1){
					auto& face_copy = mesh.faces[numN];
					
					// proc faces
					if(face_copy.getType() == MASCH_Face_Types::PROCESSOR){
						
						if(boolProcFaces_Combine[proc_num]==true){
							
							++proc_total_num;
							
							auto& groupPrcoc_face = groupProcFaces[ groupProcFaces_id[proc_num] ];
							
							face.ipoints = groupPrcoc_face.ipoints;
							face.iL = newCellsNumber[ face_copy.iL ];
							face.iR = -1;
							face.setType(face_copy.getType());
							
							for(auto& j : groupPrcoc_face.ichild){
								++numN;
								++proc_num;
							}
							--numN;
							--proc_num;
						}
						else{
							vector<int> tmpPoints;
							for(auto& j : face_copy.ipoints){
								if(boolDeletePoints[j]==false){
									tmpPoints.push_back( newPointsNumber[j] );
								}
							}
							face.ipoints.clear();
							face.ipoints = tmpPoints;
							face.iL = newCellsNumber[ face_copy.iL ];
							face.iR = -1;
							face.setType(face_copy.getType());
							
						}
						++proc_num;
					}
					else{
						vector<int> tmpPoints;
						for(auto& j : face_copy.ipoints){
							if(boolDeletePoints[j]==false){
								tmpPoints.push_back( newPointsNumber[j] );
							}
						}
						face.ipoints.clear();
						face.ipoints = tmpPoints;
						face.iL = newCellsNumber[ face_copy.iL ];
						if(face_copy.iR != -1){
							face.iR = newCellsNumber[ face_copy.iR ];
						}
						else{
							face.iR = -1;
						}
						face.setType(face_copy.getType());
					}
					++numN;
				}
				// Unrefine
				else{
					
					auto& group_face = groupOutFaces[ groupFaces_id[numN] ];
					
					auto& face_copy = mesh.faces[numN];
					
					int faceLevel = face_copy.level;
					
					bool ownUnrefine = boolCellUnrefine[face_copy.iL];
					bool ngbUnrefine;
					int ownLevel = mesh.cells[face_copy.iL].level;
					int ngbLevel;
					if(face_copy.getType() == MASCH_Face_Types::INTERNAL){
						ngbUnrefine = boolCellUnrefine[face_copy.iR];
						ngbLevel = mesh.cells[face_copy.iR].level;
					}
					else if(face_copy.getType() == MASCH_Face_Types::PROCESSOR){
						ngbUnrefine = false;
						if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
						ngbLevel = cLevel_recv[proc_num];
						proc_num += group_face.ichild.size();
					}
					else if(face_copy.getType() == MASCH_Face_Types::BOUNDARY){
						ngbUnrefine = true;
						ngbLevel = 0;
					}
					
					
					if( 
					(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel) ||
					(ownUnrefine==true && ownLevel>ngbLevel) ||
					(ngbUnrefine==true && ownLevel<ngbLevel)
					){
						
						if(face_copy.getType() == MASCH_Face_Types::PROCESSOR){
							++proc_total_num;
						}
						
						face.ipoints.clear();
						for(auto& j : group_face.ipoints){
							face.ipoints.push_back( j );
						}
						face.iL = group_face.iL;
						face.iR = group_face.iR;
						face.setType(face_copy.getType());
						
						numN += group_face.ichild.size();
						
					}
					else{
						for(auto& j : group_face.ichild){
							auto& face_target = mesh.faces[i];
							auto& face_child = mesh.faces[j];
							
							vector<int> tmpPoints;
							for(auto& k : face_child.ipoints){
								if(boolDeletePoints[k]==false){
									tmpPoints.push_back( newPointsNumber[k] );
								}
							}
							
							face_target.ipoints.clear();
							face_target.ipoints = tmpPoints;
							face_target.iL = newCellsNumber[ face_child.iL ];
							if(face_child.iR != -1){
								face_target.iR = newCellsNumber[ face_child.iR ];
							}
							else{
								face_target.iR = -1;
							}
							face_target.setType(face_child.getType());
					
							++numN;
							++i;
						}
						--i;
					}
				}	
			}
			else{
				auto& face_copy = mesh.faces[numN];
				
				int faceLevel = face_copy.level;
				
				bool ownUnrefine = boolCellUnrefine[face_copy.iL];
				bool ngbUnrefine = false;
				if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
				int ownLevel = mesh.cells[face_copy.iL].level;
				int ngbLevel = cLevel_recv[proc_num];
				
				if(
				(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel) ||
				(ownUnrefine==true && ownLevel>ngbLevel)
				){
					
					
					auto& group_face = groupOutFaces[ groupFaces_id[numN] ];
					
					face.ipoints.clear();
					for(auto& j : group_face.ipoints){
						face.ipoints.push_back( j );
					}
					face.iL = group_face.iL;
					face.iR = group_face.iR;
					face.setType(face_copy.getType());
					
					numN += group_face.ichild.size();
					proc_num += group_face.ichild.size();
					
					proc_total_num += group_face.ichild.size();
					
					
				}
				else if(ngbUnrefine==true && ownLevel<ngbLevel){
					
					auto& groupProc_face = groupProcFaces[ groupProcFaces_id[proc_num] ];
					
					face.ipoints = groupProc_face.ipoints;
					face.iL = newCellsNumber[ face_copy.iL ];
					face.iR = -1;
					face.setType(face_copy.getType());
					
					
					numN += groupProc_face.ichild.size();
					proc_num += groupProc_face.ichild.size();
						
					proc_total_num += groupProc_face.ichild.size();
					
				}
				else{
					if(groupFaces_id[numN]==-1){
					
						vector<int> tmpPoints;
						for(auto& j : face_copy.ipoints){
							if(boolDeletePoints[j]==false){
								tmpPoints.push_back( newPointsNumber[j] );
							}
						}
						face.ipoints.clear();
						face.ipoints = tmpPoints;
						face.iL = newCellsNumber[ face_copy.iL ];
						if(face_copy.iR != -1){
							face.iR = newCellsNumber[ face_copy.iR ];
						}
						else{
							face.iR = -1;
						}
						
						face.setType(face_copy.getType());
					
						++numN;
						++proc_num;
						
					
					}
					else{
						
						auto& groupFace = groupOutFaces[ groupFaces_id[numN] ];
						
						for(auto& j : groupFace.ichild){
							auto& face_target = mesh.faces[i];
							auto& face_child = mesh.faces[j];
							
							vector<int> tmpPoints;
							for(auto& k : face_child.ipoints){
								if(boolDeletePoints[k]==false){
									tmpPoints.push_back( newPointsNumber[k] );
								}
							}
							
							face_target.ipoints.clear();
							face_target.ipoints = tmpPoints;
							face_target.iL = newCellsNumber[ face_child.iL ];
							if(face_child.iR != -1){
								face_target.iR = newCellsNumber[ face_child.iR ];
							}
							else{
								face_target.iR = -1;
							}
							face_target.setType(face_child.getType());
					
							++numN;
							++proc_num;
							++i;
						}
						--i;
					}
				}
			}
		}
		if(numN > mesh.faces.size()-1){
			saveI = i;
			break;
		}
	}
	
	
	if(nBC!=mesh.boundaries.size()){
		cout << rank << " NO BOUNDARY MATCHING" << endl;
		for(auto& boundary : mesh.boundaries){
			cout << rank << " : " << boundary.rightProcNo << " " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
		}
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	int orgFacesSize = mesh.faces.size();
	for(int i=orgFacesSize-1; i>saveI; --i){
		mesh.faces.pop_back();
	}
	
	//====================================================
	// boundary setting

	int maxBCnum = mesh.boundaries.size()-1;
	if(mesh.boundaries[maxBCnum].nFaces == 0){
		mesh.boundaries[maxBCnum].startFace = mesh.faces.size();
	}
	for(int i=maxBCnum-1; i>=0; --i){
		if(mesh.boundaries[i].nFaces == 0){
			mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace;
		}
	} 
	
	for (int i=0; i<mesh.boundaries.size()-1; ++i) {
		mesh.boundaries[i].nFaces = mesh.boundaries[i+1].startFace-mesh.boundaries[i].startFace;
	}
	int maxBDsize = mesh.boundaries.size()-1;
	mesh.boundaries[maxBDsize].nFaces = mesh.faces.size()-mesh.boundaries[maxBDsize].startFace;
	
	// mesh.buildCells();
	int orgCellSize = mesh.cells.size();
	int tmpCellNum = 0;
	for(int i=0; i<orgCellSize; ++i){
		// if(groupCells_id[i] == -1){
			// mesh.cells[tmpCellNum].var.assign(
				// mesh.cells[i].var.begin(),mesh.cells[i].var.end());
		// }
		// else{
			// int varSize = mesh.cells[tmpCellNum].var.size();
			// int subCellSize = groupCellsUnrefine[groupCells_id[i]].ichild.size();
			// double dSubCellSize = (double)subCellSize;
			// vector<double> tmpVars(varSize,0.0);
			// for(auto& j : groupCellsUnrefine[groupCells_id[i]].ichild){
				// for(int k=0; k<varSize; ++k){
					// tmpVars[k] += mesh.cells[i].var[k]/dSubCellSize;
				// }
				// ++i;
			// }
			// --i;
			// mesh.cells[tmpCellNum].var.assign(tmpVars.begin(),tmpVars.end());
		// }
		
		mesh.cells[tmpCellNum].ipoints.clear();
		mesh.cells[tmpCellNum].ifaces.clear();
		
		++tmpCellNum;
		
	}
	
	
	for(int i=totalCellNum; i<orgCellSize; ++i){
		mesh.cells.pop_back();
	}
	
	// mesh.check();
	// mesh.setFaceTypes();
	// mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	// mesh.setFaceLevels();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	
	
	
	
	// level setting
	for(int i=0; i<mesh.cells.size(); ++i){
		mesh.cells[i].level = cellsLevel[i];
		mesh.cells[i].group = cellsGroup[i];
		mesh.cells[i].iparcels = cells_iparcels[i];
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
	
	this->mpiLevels(mesh, cLevel_recv);
	
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		// face.var.resize(controls.nTotalFaceVar,0.0);
		// face.varL.resize(controls.nTotalFaceLRVar,0.0);
		// face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			
			int maxLevel = 
				max(mesh.cells[face.iL].level,
					mesh.cells[face.iR].level);
			face.level = maxLevel;
			
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			
			int maxLevel = 
				max(mesh.cells[face.iL].level,
					cLevel_recv[proc_num]);
			face.level = maxLevel;
			
			++proc_num;
			
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
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
	
	
	
	
	afterCellSize = mesh.cells.size();
	afterFaceSize = mesh.faces.size();
	afterPointSize = mesh.points.size();
	
	int deletedCellSize = beforeCellSize - afterCellSize;
	int deletedFaceSize = beforeFaceSize - afterFaceSize;
	int deletedPointSize = beforePointSize - afterPointSize;
	
	int deletedCellSize_glb, deletedFaceSize_glb, deletedPointSize_glb;
	MPI_Allreduce(&deletedCellSize, &deletedCellSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&deletedFaceSize, &deletedFaceSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&deletedPointSize, &deletedPointSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << 
		"| cell = -" << deletedCellSize_glb <<
		" | face = -" << deletedFaceSize_glb <<
		" | point = -" << deletedPointSize_glb << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	if(rank==0){
		// cout << "| AMR - Unrefine completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}

























void MASCH_Poly_AMR_Builder::sortCellCanUnrefine(
	vector<int>& vLevel, 
	int saveLevel, 
	int& num, 
	vector<vector<int>>& vGroupLevel, 
	vector<vector<int>>& vGroupNumber
	) {
		
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
		
	vector<int> tmpLevels;
	vector<int> tmpNumbers;
	int limitNum = 7;
	if (saveLevel <= 1) limitNum = 2147483647;
	while (1) {
		if (vLevel.size() - 1 >= num) {
			int level = vLevel.at(num);
			// cout << num << " " << tmpLevels.size()<< " " << saveLevel << " " << level << endl;

			if(
			saveLevel == level &&
			tmpLevels.size() < limitNum
			){
				tmpLevels.push_back(saveLevel);
				tmpNumbers.push_back(num);
			}
			else if(
			saveLevel == level &&
			tmpLevels.size() == limitNum
			){
				if (std::find(tmpNumbers.begin(), tmpNumbers.end(), -1) == tmpNumbers.end()) {
					tmpLevels.push_back(saveLevel);
					tmpNumbers.push_back(num);
					vGroupLevel.push_back(tmpLevels);
					vGroupNumber.push_back(tmpNumbers);
				}
				break;
			}
			else if(
			saveLevel < level
			){
				// cout << "start" << endl;

				if (tmpLevels.size() == limitNum+1) {
					tmpLevels.clear();
					tmpNumbers.clear();
				}

				sortCellCanUnrefine(vLevel, saveLevel + 1, num, vGroupLevel, vGroupNumber);
				tmpLevels.push_back(-1);
				tmpNumbers.push_back(-1);
				
				if (tmpLevels.size() == limitNum + 1) {
					break;
				}
				// cout << tmpLevels.size() << endl;
			}
			else {
				--num;
				break;
			}
			++num;
		}
		else {
			
			if(saveLevel>1 && tmpLevels.size()!=8){
				cout << "WWWWWWWWWWWWWWW" << " " << rank << " " << saveLevel << " " << tmpLevels.size() << " " << vLevel.size() << " " << num << endl;
				
				cout << rank << " ";
				for(auto& k : vLevel){
					cout << k << " ";
				}
				cout << endl;
				
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			break;
		}
	}
   
}






void MASCH_Poly_AMR_Builder::restrictCellUnrefine(
	MASCH_Mesh& mesh, 
	vector<bool>& boolCellUnrefine,
	vector<int>& cLevel_recv, 
	vector<int>& cUnrefine_recv){

	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			if(mesh.cells[face.iL].level < mesh.cells[face.iR].level){
				boolCellUnrefine[face.iL] = false;
			}
			if(mesh.cells[face.iL].level > mesh.cells[face.iR].level){
				boolCellUnrefine[face.iR] = false;
			}
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			if(mesh.cells[face.iL].level < cLevel_recv[proc_num]){
				boolCellUnrefine[face.iL] = false;
			}
			++proc_num;
		}
	}
}



void MASCH_Poly_AMR_Builder::extractGroupCellListsCanUnrefine(
	MASCH_Mesh& mesh, 
	vector<vector<int>>& groupCellListsCanUnrefine
){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	vector<vector<int>> groupCellListsLevelZero;
	vector<vector<int>> groupCellLevelListsLevelZero;
	for(int i=0; i<mesh.cells.size(); ++i){
		
		groupCellListsLevelZero.push_back(vector<int>(1,i));
		groupCellLevelListsLevelZero.push_back(vector<int>(1,mesh.cells[i].level));
		
		if(mesh.cells[i].level<=0) continue;
		
		int groupOrg = mesh.cells[i].group;
		
		++i;
		while(groupOrg==mesh.cells[i].group && i<mesh.cells.size()){
			groupCellListsLevelZero.back().push_back(i);
			groupCellLevelListsLevelZero.back().push_back(mesh.cells[i].level);
			++i;
		}
		--i;
	}
	
	
	int testNum = 0;
	
	
	
	vector<int> test;
	for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		if(groupCellLevelListsLevelZero[i].size()==1) continue;
		
		bool boolAllLevelOne = true;
		for(auto& j : groupCellLevelListsLevelZero[i]){
			if(j != 1) {
				boolAllLevelOne = false;
				break;
			}
		}
		
		
		
		if(boolAllLevelOne==true){
			test.push_back(i);
			groupCellListsCanUnrefine.push_back(groupCellListsLevelZero[i]);
		}
		else{
			vector<vector<int>> tmpGroupLevel;
			vector<vector<int>> tmpGroupNumber;
			int tmpNum2 = 0;
			
			sortCellCanUnrefine(
				groupCellLevelListsLevelZero[i],
				0,
				tmpNum2,
				tmpGroupLevel,
				tmpGroupNumber
				);
		
			int startCellNum = groupCellListsLevelZero[i][0];
			for(auto& j : tmpGroupNumber){
				for(auto& k : j){
					k += startCellNum;
				}
			}
			
			for(auto& j : tmpGroupNumber){
				test.push_back(i);
				groupCellListsCanUnrefine.push_back(j);
			}
		}
		
	}
	
	
	
	
	for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){

		auto& group = groupCellListsCanUnrefine[i];
		
		auto& target_Cell = mesh.cells[group[0]];
		int cellCenterPoint = -1;
		for(auto& targetPoint : target_Cell.ipoints){
			bool boolCellCenterPoint = true;
			for(int j=1; j<group.size(); ++j){
				auto& cell = mesh.cells[group[j]];
				if(
				std::find(cell.ipoints.begin(),cell.ipoints.end(),targetPoint)
				== cell.ipoints.end()
				){
					boolCellCenterPoint = false;
					break;
				}
			}
			if(boolCellCenterPoint==true){
				cellCenterPoint = targetPoint;
				break;
			}
		}
		
		
		if(cellCenterPoint==-1){
			
			
			
			
			for(int j=0; j<group.size(); ++j){
				auto& cell = mesh.cells[group[j]];
				if(rank==0) cout << cell.level << endl;
				for(auto& k : cell.ipoints){
					if(rank==0) cout << rank << " " << group[j] << " " << k << endl;
				}
				
			}
			cout << endl;
			
			if(rank==0) cout << "cellLevel" << endl;
			
			
			testNum = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << groupCellListsLevelZero[test[i]][testNum] << " " << groupCellLevelListsLevelZero[test[i]][testNum] << " ";
				++testNum;
			}
			
			cout << endl;
			
			testNum = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << j << " ";
				++testNum;
			}
			cout << endl;
			
			testNum = 0;
			for(auto& j : groupCellListsLevelZero[0]){
				if(rank==0) cout << j << " ";
				++testNum;
			}	
			cout << endl;
			
			testNum = 0;
			int saveFirst = groupCellListsLevelZero[test[i]][0];
			int saveLast = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << mesh.cells[j].group << " ";
				saveLast = j;
				++testNum;
			}
			
			
			cout << endl;
			if(rank==0) cout << mesh.cells[saveFirst-1].group << endl;
			if(rank==0) cout << mesh.cells[saveLast+1].group << endl;
			if(rank==0) cout << mesh.cells.size() << endl;
			
			for(int j=0; j<mesh.cells.size(); ++j){
				if(rank==0) cout << mesh.cells[j].level << " ";
			}
			
			if(rank==0) MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	}
	
}





void MASCH_Poly_AMR_Builder::extractGroupUnrefineCells(
	MASCH_Mesh& mesh, 
	vector<bool>& boolCellUnrefine,
	vector<vector<int>>& groupCellListsCanUnrefine,
	vector<groupCells_Unrefine>& groupCellsUnrefine
	){

	for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		bool boolAllUnrefine = true;
		for(auto& j : groupCellListsCanUnrefine[i]){
			if(boolCellUnrefine[j] == false) {
				boolAllUnrefine = false;
				break;
			}
		}
		
		if(boolAllUnrefine){
			groupCellsUnrefine.push_back(groupCells_Unrefine());
			for(auto& j : groupCellListsCanUnrefine[i]){
				groupCellsUnrefine.back().ichild.push_back(j);
			}
		}
	}
}


