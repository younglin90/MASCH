
#include "./mesh.h"


void MASCH_Mesh::getFaceOrders(
int maxLevel, vector<MASCH_Point>& points, vector<MASCH_Face>& faces, 
vector<int>& procFaces, vector<int>& reorderFaceIds) {
	
	int faceSize = procFaces.size();

	
	reorderFaceIds.clear();

	vector<vector<int>> masterFaces(maxLevel+1);
	for (int level = 1; level <= maxLevel; ++level) {
		for (int i = 0; i < faceSize; ++i) {
			bool boolStartFace = false;
			for (auto& ipoint : faces[procFaces[i]].ipoints) {
				if (points[ipoint].level <= level-1) {
					boolStartFace = true;
					break;
				}
			}
			if (boolStartFace == true) masterFaces[level].push_back( procFaces[i] );
		}
	}

	{
		reverse(masterFaces[1].begin() + 1, masterFaces[1].end());
		for (auto& item : masterFaces[1]) {
			reorderFaceIds.push_back(item);
		}
		reverse(masterFaces[1].begin() + 1, masterFaces[1].end());

		for (int level = 1; level < maxLevel; ++level) {
			int nextLv = level + 1;
			for (int i = 0; i < masterFaces[level].size() - 1; ++i) {
				int str_id = masterFaces[level][i];
				int end_id = masterFaces[level][i + 1];
				int str = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), str_id) - masterFaces[nextLv].begin();
				int end = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), end_id) - masterFaces[nextLv].begin();
				vector<int> tmp_vec;
				for (int j = str; j < end; ++j) {
					tmp_vec.push_back(masterFaces[nextLv][j]);
				}
				if (tmp_vec.size() > 0) {
					reverse(tmp_vec.begin() + 1, tmp_vec.end());
					int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						tmp_vec.begin() + 1, tmp_vec.end());
				}
			}
			{
				int str_id = masterFaces[level].back();
				int str = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), str_id) - masterFaces[nextLv].begin();
				vector<int> tmp_vec;
				for (int j = str; j < masterFaces[nextLv].size(); ++j) {
					tmp_vec.push_back(masterFaces[nextLv][j]);
				}
				if (tmp_vec.size() > 0) {
					reverse(tmp_vec.begin() + 1, tmp_vec.end());
					int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						tmp_vec.begin() + 1, tmp_vec.end());
				}
			}
		}
	}
	
	// int faceSize = procFaces.size();

	
	// reorderFaceIds.clear();

	// vector<vector<int>> masterFaces(maxLevel+1);
	// for (int level = 1; level <= maxLevel; ++level) {
		// for (int i = 0; i < faceSize; ++i) {
			// bool boolStartFace = false;
			// for (auto& ipoint : faces[procFaces[i]].ipoints) {
				// if (points[ipoint].level <= level-1) {
					// boolStartFace = true;
					// break;
				// }
			// }
			// if (boolStartFace == true) masterFaces[level].push_back( procFaces[i] );
		// }
	// }

	// {
		// // if(masterFaces.size()<=1) cout << "AAAAAAAAA" << endl;
		// // if(masterFaces[1].size()<1) cout << "AAAAAAAAA" << endl;
		// // reverse(masterFaces[1].begin() + 1, masterFaces[1].end());
		// for (auto& item : masterFaces[1]) {
			// reorderFaceIds.push_back(item);
		// }
		// // reverse(masterFaces[1].begin() + 1, masterFaces[1].end());

		// for (int level = 1; level < maxLevel; ++level) {
			// int nextLv = level + 1;
			// // if(masterFaces.size()<=level) cout << "AAAAAAA" << endl;
			// auto& masterFaces_lv = masterFaces[level];
			// auto& masterFaces_nextLv = masterFaces[nextLv];
			// for (int i = 0; i < masterFaces_lv.size() - 1; ++i) {
				// int str_id = masterFaces_lv[i];
				// int end_id = masterFaces_lv[i + 1];
				// auto str_ptr = find(masterFaces_nextLv.begin(), masterFaces_nextLv.end(), str_id);
				// auto end_ptr = find(masterFaces_nextLv.begin(), masterFaces_nextLv.end(), end_id);
				// int str = 0;
				// // if(str_ptr==masterFaces_lv.end()) {
					// // str = masterFaces_lv.size()-1;
				// // }
				// // else{
					// // str = str_ptr - masterFaces_nextLv.begin();
					// str = masterFaces_nextLv.begin() - masterFaces_nextLv.begin();
				// // }
				// int end = 0;
				// // if(end_ptr==masterFaces_nextLv.end()) {
					// // end = masterFaces_nextLv.size()-1;
				// // }
				// // else{
					// // end = end_ptr - masterFaces_nextLv.begin();
					// end = masterFaces_nextLv.begin() - masterFaces_nextLv.begin();
				// // }
				// vector<int> tmp_vec2;
				// // cout << str << " " << end << endl;
				// for (int jj = str; jj < end; ++jj) {
					// // tmp_vec2.push_back(masterFaces_nextLv.at(j));
				// }
				// // // if (tmp_vec.size() > 0) {
					// // // // reverse(tmp_vec.begin() + 1, tmp_vec.end());
					// // // int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					// // // // reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						// // // // tmp_vec.begin() + 1, tmp_vec.end());
				// // // }
			// }
			// // {
				// // int str_id = masterFaces[level].back();
				// // auto str_ptr = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), str_id);
				// // int str = 0;
				// // if(str_ptr==masterFaces[nextLv].end()) {
					// // str = masterFaces[nextLv].size()-1;
				// // }
				// // else{
					// // str = str_ptr - masterFaces[nextLv].begin();
				// // }
				// // vector<int> tmp_vec;
				// // for (int j = str; j < masterFaces[nextLv].size(); ++j) {
					// // tmp_vec.push_back(masterFaces[nextLv][j]);
				// // }
				// // // if (tmp_vec.size() > 0) {
					// // // // reverse(tmp_vec.begin() + 1, tmp_vec.end());
					// // // int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					// // // // reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						// // // // tmp_vec.begin() + 1, tmp_vec.end());
				// // // }
			// // }
		// }
	// }

};


// void getFaceOrders(
// int maxLevel, vector<vector<int>>& facePoints, 
// vector<int>& pointLevels, vector<int>& reorderFaceIds) {
	// int faceSize = facePoints.size();

	// reorderFaceIds.clear();
	// //reorderFaceIds.resize(faceSize,-1);

	// vector<vector<int>> masterFaces(maxLevel+1);
	// for (int level = 1; level <= maxLevel; ++level) {
		// for (int i = 0; i < faceSize; ++i) {
			// bool boolStartFace = false;
			// for (auto& ipoint : facePoints[i]) {
				// if (pointLevels[ipoint] <= level-1) {
					// boolStartFace = true;
					// break;
				// }
			// }
			// if (boolStartFace == true) masterFaces[level].push_back(i);
		// }
	// }

	// {
		// reverse(masterFaces[1].begin() + 1, masterFaces[1].end());
		// for (auto& item : masterFaces[1]) {
			// reorderFaceIds.push_back(item);
		// }
		// reverse(masterFaces[1].begin() + 1, masterFaces[1].end());

		// for (int level = 1; level < maxLevel; ++level) {
			// int nextLv = level + 1;
			// for (int i = 0; i < masterFaces[level].size() - 1; ++i) {
				// int str_id = masterFaces[level][i];
				// int end_id = masterFaces[level][i + 1];
				// int str = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), str_id) - masterFaces[nextLv].begin();
				// int end = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), end_id) - masterFaces[nextLv].begin();
				// vector<int> tmp_vec;
				// for (int j = str; j < end; ++j) {
					// tmp_vec.push_back(masterFaces[nextLv][j]);
				// }
				// if (tmp_vec.size() > 0) {
					// reverse(tmp_vec.begin() + 1, tmp_vec.end());
					// int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					// reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						// tmp_vec.begin() + 1, tmp_vec.end());
				// }
			// }
			// {
				// int str_id = masterFaces[level].back();
				// int str = find(masterFaces[nextLv].begin(), masterFaces[nextLv].end(), str_id) - masterFaces[nextLv].begin();
				// vector<int> tmp_vec;
				// for (int j = str; j < masterFaces[nextLv].size(); ++j) {
					// tmp_vec.push_back(masterFaces[nextLv][j]);
				// }
				// if (tmp_vec.size() > 0) {
					// reverse(tmp_vec.begin() + 1, tmp_vec.end());
					// int insert_str = find(reorderFaceIds.begin(), reorderFaceIds.end(), tmp_vec[0]) - reorderFaceIds.begin();
					// reorderFaceIds.insert(reorderFaceIds.begin() + insert_str + 1,
						// tmp_vec.begin() + 1, tmp_vec.end());
				// }
			// }
		// }
	// }

	// //for (auto& item : reorderFaceIds) {
	// //		cout << item << endl;
	// //}


// };


void MASCH_Mesh::repartitioning(
vector<int>& idBlockCell, 
int maxLevel, 
vector<int>& to_new_cell_id,
vector<int>& parcel_ip
){
	
	
	int rank = (MPI::COMM_WORLD.Get_rank()); 
	int size = (MPI::COMM_WORLD.Get_size()); 
	
	auto& mesh = (*this);
	
		// if(rank==0){
			
			// int leng = 10;
			// int leng_glo;
			// MPI_Scatter(&leng, 1, MPI_INT, &leng_glo, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// }
		// else{
			
			// int leng = 10;
			// int leng_glo;
			// MPI_Scatter(NULL, 1, MPI_INT, &leng_glo, 1, MPI_INT, 0, MPI_COMM_WORLD);
			// cout << leng_glo << endl;
		// }
	// to_new_cell_id.resize(mesh.cells.size());
	
	
    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute Load Balancing ... ";
	}
	
	MASCH_MPI mpi;
	
	vector<MASCH_Mesh> meshNew;
	MASCH_Mesh meshComb;
	
	
	vector<int> recv_idBlockCell;
	vector<int> recv_rank;
	{
		vector<int> send_idBlockCell;
		vector<int> send_rank;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType()==MASCH_Face_Types::PROCESSOR){
				send_idBlockCell.push_back(idBlockCell[face.iL]);
				send_rank.push_back(rank);
			}
		}
		recv_idBlockCell.resize(send_idBlockCell.size());
		MPI_Alltoallv( send_idBlockCell.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_idBlockCell.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
		recv_rank.resize(send_rank.size());
		MPI_Alltoallv( send_rank.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_rank.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	
	
	
	int PR2IN = -1;
	int BC2BC = -2;
	int IN2PR = -4;
	int PR2PR = -3;
	vector<int> face_state;
	
	// ************************
	vector<vector<int>> recv_localCell_group;
	vector<vector<int>> recv_localCell_level;
	vector<int> nCells_local(size+1,0);
	// ************************
	//*************************
	vector<bool> procFace_boolReorder;
	//*************************
	
	//=======================================	
	// parcel 위치 정보 넘기기
	vector<vector<int>> recv_localParcel_size;
	vector<vector<int>> recv_localParcel_id;
	vector<vector<int>> recv_localParcel_icell;
	//=======================================	
		
	// 포인트
	{
		
		// 넘길 포인트 process 정립
		vector<vector<int>> send_localPoint_proc(mesh.points.size());
		for(int i=0, ip=0; i<mesh.cells.size(); ++i){
			int proc = idBlockCell[i];
			for(auto& ipoint : mesh.cells[i].ipoints){
				auto& point = send_localPoint_proc[ipoint];
				if(find(point.begin(),point.end(),proc)==point.end()){
					point.push_back(proc);
				}
			}
		}
		
		
		// 옮겨질 포인트 정보들 정리
		vector<vector<double>> send_localPoint_xyz(size);
		vector<vector<int>> send_localPoint_level(size);
		vector<vector<pair<int,int>>> send_localPoint_proc_id(mesh.points.size());
		vector<int> send_localPoint_n(size,0);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& proc : send_localPoint_proc[i]){
				send_localPoint_proc_id[i].push_back(make_pair(proc,send_localPoint_n[proc]++));
				send_localPoint_xyz[proc].push_back(point.x);
				send_localPoint_xyz[proc].push_back(point.y);
				send_localPoint_xyz[proc].push_back(point.z);
				send_localPoint_level[proc].push_back(point.level);
			}
		}
		
		// 포인트 x,y,z 옮기기
		vector<vector<double>> recv_localPoint_xyz;
		mpi.Alltoallv(send_localPoint_xyz, recv_localPoint_xyz);
		vector<vector<int>> recv_localPoint_level;
		mpi.Alltoallv(send_localPoint_level, recv_localPoint_level);
		
		// 옮겨진 포인트 x,y,z 넣기
		vector<int> str_points_glo(size+1,0);
		for(int ip=0, iter_glob=0; ip<size; ++ip){
			int tmp_size = recv_localPoint_xyz[ip].size()/3;
			str_points_glo[ip+1] = tmp_size;
			for(int i=0, iter=0; i<tmp_size; ++i){
				meshComb.addPoint();
				meshComb.points.back().x = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().y = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().z = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().level = recv_localPoint_level[ip][i];
			}
		}
		// 각 proc에 대한 포인트 시작지점
		for(int ip=0; ip<size; ++ip){
			str_points_glo[ip+1] = str_points_glo[ip+1] + str_points_glo[ip];
		}
		
		
		
		//========================================================
		// connPoints 대한 정보들 
		vector<vector<int>> send_localPoint_connId(size);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [proc, id] : point.connPoints){
				send_localPoint_connId[proc].push_back(id);
				send_localPoint_connId[proc].push_back(i);
				send_localPoint_connId[proc].push_back(send_localPoint_proc_id[i].size());
				for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
					send_localPoint_connId[proc].push_back(tmp_proc);
					send_localPoint_connId[proc].push_back(tmp_id_local);
				}
			}
		}
		
		// connPoints 대한 정보들 넘기기 
		vector<vector<int>> recv_localPoint_connId;
		mpi.Alltoallv(send_localPoint_connId, recv_localPoint_connId);
		
		// 받은 connPoints 대한 정보들에 대해서, 같은 proc으로 이동하는 포인트들 정보 입력 
		vector<vector<int>> send_localPoint_toProc_toId(size);
		// vector<vector<int>> send_globalPoint_toProc_toId(size);
		for(int ip=0; ip<size; ++ip){
			for(int i=0; i<recv_localPoint_connId[ip].size(); ++i){
				int ipoint = recv_localPoint_connId[ip][i];
				int right_i = recv_localPoint_connId[ip][++i];
				int tmp_size = recv_localPoint_connId[ip][++i];
				for(int j=0; j<tmp_size; ++j){
					int tmp_proc = recv_localPoint_connId[ip][++i];
					int tmp_id_local = recv_localPoint_connId[ip][++i];
					for(auto& [proc, id] : send_localPoint_proc_id[ipoint]){
						if(proc==tmp_proc){
							// if(n_tmp_nnnn2[ipoint] == 1){
								// if(rank==0) cout << ipoint << " " << proc << " " << id << " " << ip << " " << tmp_id_local << endl;
								send_localPoint_toProc_toId[proc].push_back(id); // proc 블록 에서 자신의 local id  
								send_localPoint_toProc_toId[proc].push_back(ip); // proc 블록 에서 자신과 중복한 포인트의 local proc
								send_localPoint_toProc_toId[proc].push_back(tmp_id_local); // proc 블록 에서 자신과 중복한 포인트의 local id
							// }
						}
					}
				}
			}
		}
		
		// 같은 proc으로 이동하는 포인트들 정보 넘기기
		vector<vector<int>> recv_localPoint_toProc_toId;
		mpi.Alltoallv(send_localPoint_toProc_toId, recv_localPoint_toProc_toId);
		
		
		// 삭제되는 중복 포인트들 찾기
		vector<bool> deletePoints(str_points_glo[size],false);
		vector<vector<int>> reorder_my_id_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightProcNo_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightId_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightId_global(str_points_glo[size]);
		{
			for(int ip=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=0, iter=0; i<recv_localPoint_toProc_toId[ip].size()/3; ++i){
					int my_id_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int my_id_global = str_points_glo[ip] + my_id_local;
					int rightProc_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_global = str_points_glo[rightProc_local] + rightId_local;
					
					reorder_my_id_local[my_id_global].push_back(my_id_local);
					reorder_rightProcNo_local[my_id_global].push_back(rightProc_local);
					reorder_rightId_local[my_id_global].push_back(rightId_local);
					reorder_rightId_global[my_id_global].push_back(rightId_global);
					
					// 중복 및 삭제
					if(rightId_global<my_id_global){
						deletePoints.at(my_id_global) = true;
					}
				}
			}
		}
		
		
		// 로컬 포인트 id -> 글로벌 포인트 id
		vector<vector<int>> points_id_local2global(size);
		{
			for(int ip=0, iter=0, before_nd=0, tmp_glob=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=str; i<end; ++i){
					points_id_local2global[ip].push_back(tmp_glob);
					if(deletePoints[i]==false) {
						++tmp_glob;
					}
					else{
						++before_nd;
					}
				}
			}
			for(int ip=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=str, iter=0; i<end; ++i){
					int iter = 0;
					int my_id_global = i;
					int min_rightId_global = str_points_glo[size];
					
					
					for(auto& my_id_local : reorder_my_id_local[i]){
						
						int rightProc_local = reorder_rightProcNo_local[i][iter];
						int rightId_local = reorder_rightId_local[i][iter];
						int rightId_global = reorder_rightId_global[i][iter];
						
						int copy_tmp = points_id_local2global[ip][my_id_local];
						if(rightId_global>i){
							points_id_local2global[rightProc_local][rightId_local] = copy_tmp;
						}
						
						++iter;
					}
				}
			}
		}
		
		
		
		//========================================================
		
		
		// connPoints 대한 정보들 
		vector<vector<int>> send_conn_proc2(size);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [proc, id] : point.connPoints){
				send_conn_proc2[proc].push_back(id);
				send_conn_proc2[proc].push_back(send_localPoint_proc_id[i].size());
				for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
					send_conn_proc2[proc].push_back(tmp_proc);
					send_conn_proc2[proc].push_back(rank);
					send_conn_proc2[proc].push_back(tmp_id_local);
				}
			}
		}
		
		// connPoints 대한 정보들 넘기기
		vector<vector<int>> recv_conn_proc2;
		mpi.Alltoallv(send_conn_proc2, recv_conn_proc2);
		
		
		// 각 포인트에 대한, 넘겨진 포인트의 옆 proc_glo, proc_loc, id_loc 저장
		vector<vector<int>> connPoints_all_send_procNo_glo(mesh.points.size());
		vector<vector<int>> connPoints_all_send_procNo_loc(mesh.points.size());
		vector<vector<int>> connPoints_all_send_id_loc(mesh.points.size());
		for(int ip=0, id=0; ip<size; ++ip){
			int iter=0;
			while(iter<recv_conn_proc2[ip].size()){
				int id_glo = recv_conn_proc2[ip][iter++];
				int tmp_size = recv_conn_proc2[ip][iter++];
				for(int i=0; i<tmp_size; ++i){
					int rightProcNo_glo = recv_conn_proc2[ip][iter++];
					int rightProcNo_loc = recv_conn_proc2[ip][iter++];
					int rightId_loc = recv_conn_proc2[ip][iter++];
					
					connPoints_all_send_procNo_glo[id_glo].push_back(rightProcNo_glo);
					connPoints_all_send_procNo_loc[id_glo].push_back(rightProcNo_loc);
					connPoints_all_send_id_loc[id_glo].push_back(rightId_loc);
					
				}
			}
		}
		
		// 각 포인트에 대한, 자기 자신 포인트의 옆 proc_glo, proc_loc, id_loc 저장
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
				connPoints_all_send_procNo_glo[i].push_back(tmp_proc);
				connPoints_all_send_procNo_loc[i].push_back(rank);
				connPoints_all_send_id_loc[i].push_back(tmp_id_local);
			}
			
		}
		
		
		
		
		
		vector<vector<int>> connPoints_all_send_id_glo(mesh.points.size());
		{
			vector<vector<int>> send_debug_procNo_loc(size);
			vector<vector<int>> send_debug_id_loc(size);
			vector<vector<double>> send_debug_8x(size);
			vector<vector<double>> send_debug_8y(size);
			vector<vector<double>> send_debug_8z(size);
			for(int i=0, ip=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				int iter=0;
				for(auto& proc : connPoints_all_send_procNo_glo[i]){
					int procNo_loc = connPoints_all_send_procNo_loc[i][iter];
					int id_loc = connPoints_all_send_id_loc[i][iter];
					
					send_debug_procNo_loc[proc].push_back(procNo_loc);
					send_debug_id_loc[proc].push_back(id_loc);
					
					send_debug_8x[proc].push_back(mesh.points[i].x);
					send_debug_8y[proc].push_back(mesh.points[i].y);
					send_debug_8z[proc].push_back(mesh.points[i].z);
					
					++iter;
				}
				
			}
			vector<vector<int>> recv_debug_procNo_loc;
			mpi.Alltoallv(send_debug_procNo_loc, recv_debug_procNo_loc);
			vector<vector<int>> recv_debug_id_loc;
			mpi.Alltoallv(send_debug_id_loc, recv_debug_id_loc);
			// vector<vector<int>> recv_debug_id_glo;
			// mpi.Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			vector<vector<double>> recv_debug_8x;
			mpi.Alltoallv(send_debug_8x, recv_debug_8x);
			vector<vector<double>> recv_debug_8y;
			mpi.Alltoallv(send_debug_8y, recv_debug_8y);
			vector<vector<double>> recv_debug_8z;
			mpi.Alltoallv(send_debug_8z, recv_debug_8z);
			
				
			vector<vector<int>> send_debug_id_glo(size);
			for(int ip=0; ip<size; ++ip){
				for(int i=0; i<recv_debug_procNo_loc[ip].size(); ++i){
					int procNo_loc = recv_debug_procNo_loc[ip][i];
					int id_loc = recv_debug_id_loc[ip][i];
					
					
					int id_glo = points_id_local2global[procNo_loc][id_loc];
					send_debug_id_glo[ip].push_back(id_glo);
					
					auto& point = meshComb.points[str_points_glo[procNo_loc]+id_loc];
					
				}
			}
			vector<vector<int>> recv_debug_id_glo;
			mpi.Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			
			vector<int> iter_tmp_n(size,0);
			for(int i=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				int iter=0;
				for(auto& proc : connPoints_all_send_procNo_glo[i]){
					int procNo_loc = connPoints_all_send_procNo_loc[i][iter];
					int id_loc = connPoints_all_send_id_loc[i][iter];
					
					int id_glo = recv_debug_id_glo[proc].at(iter_tmp_n[proc]);
					connPoints_all_send_id_glo[i].push_back(id_glo);
					
					++iter_tmp_n[proc];
					++iter;
				}
				
			}
				
		}
		
	
			
		// 포인트 삭제
		{
			int numN = 0;
			meshComb.points.erase( std::remove_if( meshComb.points.begin(), meshComb.points.end(), 
				[&deletePoints, &numN](MASCH_Point const& v) { 
				return deletePoints[numN++]; 
				}), meshComb.points.end());
		}
			
			
		{
			vector<vector<int>> send_debug_tmp_size(size);
			vector<vector<int>> send_debug_procNo_glo(size);
			vector<vector<int>> send_debug_procNo_loc(size);
			vector<vector<int>> send_debug_id_loc(size);
			vector<vector<int>> send_debug_id_glo(size);
			for(int i=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				for(auto& proc : send_localPoint_proc[i]){
					int tmp_size = connPoints_all_send_procNo_glo[i].size();
					send_debug_tmp_size[proc].push_back(tmp_size);
					for(int j=0; j<tmp_size; ++j){
						send_debug_procNo_glo[proc].push_back(connPoints_all_send_procNo_glo[i].at(j));
						send_debug_procNo_loc[proc].push_back(connPoints_all_send_procNo_loc[i].at(j));
						send_debug_id_loc[proc].push_back(connPoints_all_send_id_loc[i].at(j));
						send_debug_id_glo[proc].push_back(connPoints_all_send_id_glo[i].at(j));
						
					}
				}
			}
			vector<vector<int>> recv_debug_tmp_size;
			mpi.Alltoallv(send_debug_tmp_size, recv_debug_tmp_size);
			vector<vector<int>> recv_debug_procNo_glo;
			mpi.Alltoallv(send_debug_procNo_glo, recv_debug_procNo_glo);
			vector<vector<int>> recv_debug_procNo_loc;
			mpi.Alltoallv(send_debug_procNo_loc, recv_debug_procNo_loc);
			vector<vector<int>> recv_debug_id_loc;
			mpi.Alltoallv(send_debug_id_loc, recv_debug_id_loc);
			vector<vector<int>> recv_debug_id_glo;
			mpi.Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			
			
			for(int ip=0; ip<size; ++ip){
				int tmp_size = recv_debug_tmp_size[ip].size();
				for(int i=0, iter=0, iter2=0; i<tmp_size; ++i){
					
					int tmp2_size = recv_debug_tmp_size[ip][i];
					int llll=0;
					int iter_glob_tmp = -1;
					for(int j=0; j<tmp2_size; ++j){
						int procNo_glo = recv_debug_procNo_glo[ip][iter];
						int procNo_loc = recv_debug_procNo_loc[ip][iter];
						int id_loc = recv_debug_id_loc[ip][iter];
						int id_glo = recv_debug_id_glo[ip][iter];
						if(rank==procNo_glo){
							iter_glob_tmp = id_glo;
							// if(rank==0) cout << procNo_glo << " " << procNo_loc << " " << id_loc << " " << id_glo << " " << iter_glob << endl;
							// ++llll;
						}
						++iter;
					}
					for(int j=0; j<tmp2_size; ++j){
						int procNo_glo = recv_debug_procNo_glo[ip][iter2];
						int procNo_loc = recv_debug_procNo_loc[ip][iter2];
						int id_loc = recv_debug_id_loc[ip][iter2];
						int id_glo = recv_debug_id_glo[ip][iter2];
						if(rank!=procNo_glo){
							auto& connPoints = meshComb.points.at(iter_glob_tmp).connPoints;
							bool thereProc = false;
							for(auto& [proc, id] : connPoints){
								if(proc==procNo_glo) thereProc = true;
							}
							if(thereProc==false) 
								connPoints.push_back(make_pair(procNo_glo,id_glo));
						}
						++iter2;
					}
					
					
					
				}
				
			}
		}
		
		//======================================
		
		// ********************************************
		vector<vector<int>> send_localCell_group(size);
		vector<vector<int>> send_localCell_level(size);
		// ********************************************

		
		vector<pair<int,int>> send_localCell_proc_id;
		vector<vector<int>> send_localCell_id(size);
		vector<int> send_localCell_n(size,0);
		for(int i=0, ip=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			int proc = idBlockCell[i];
			int tmp_nCell = send_localCell_n.at(proc)++;
			send_localCell_proc_id.push_back(make_pair(proc,tmp_nCell));
			send_localCell_id[proc].push_back(tmp_nCell);
			send_localCell_group[proc].push_back(cell.group);
			send_localCell_level[proc].push_back(cell.level);
		}
		
		// //=======================================
		// // 디버그
		// {
			// int max_group_id = -1;
			// int min_group_id = 1e8;
			// for(int i=0, ip=0; i<mesh.cells.size(); ++i){
				// auto& cell = mesh.cells[i];
				// max_group_id = max(cell.group,max_group_id);
				// min_group_id = min(cell.group,min_group_id);
			// }
			
			// vector<vector<int>> debug_localCell_group(max_group_id-min_group_id+1);
			// for(int i=0, ip=0; i<mesh.cells.size(); ++i){
				// auto& cell = mesh.cells[i];
				// int proc = idBlockCell[i];
				// debug_localCell_group.at(cell.group-min_group_id).push_back(proc);
			// }
			
			// for(auto& item : debug_localCell_group){
				// int proc0 = item[0];
				// for(auto& item2 : item){
					// if(proc0 != item2){
						// cout << "proc0!=item2 " << proc0 << " " << item2 << endl;
					// }
				// }
				// // cout << endl;
			// }
			
		// }
		// //=======================================
		
		
		vector<vector<int>> recv_localCell_id;
		mpi.Alltoallv(send_localCell_id, recv_localCell_id);
		
		
		mpi.Alltoallv(send_localCell_group, recv_localCell_group);
		mpi.Alltoallv(send_localCell_level, recv_localCell_level);
		
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = recv_localCell_id[ip].size();
		}
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = nCells_local[ip+1] + nCells_local[ip];
		}
		
		
		
		
		
		vector<int> boundary_type(mesh.faces.size(),0);
		int nbc=0;
		for(int i=0; i<mesh.boundaries.size(); ++i){
			auto& boundary = mesh.boundaries[i];
			if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
				++nbc;
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int j=str; j<end; ++j){
					boundary_type[j] = i;
				}
			}
		}
		
		
		
		// processor 페이스 고유 번호 부여
		vector<int> proc_numbering;
		{
			int maxProcFaceNum = 0;
			for(int ip=rank; ip<size; ++ip){
				maxProcFaceNum += mesh.countsSendProcFaces[ip];
			}
			// if(rank==0) cout << maxProcFaceNum << endl;
			vector<int> tmp_maxProcFaceNum(size);
			MPI_Allgather(&maxProcFaceNum, 1, MPI_INT, tmp_maxProcFaceNum.data(), 1, MPI_INT, MPI_COMM_WORLD);
			
			vector<int> str_proc_numbering(size+1,0);
			for(int ip=0; ip<size; ++ip){
				str_proc_numbering[ip+1] = str_proc_numbering[ip] + tmp_maxProcFaceNum[ip];
			}
			
			
			// vector<int> send_proc_numbering;
			// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				// if(face.getType()==MASCH_Face_Types::PROCESSOR){
					// int tmp_num = str_proc_numbering[rank] + ip;
					// send_proc_numbering.push_back(tmp_num);
					// ++ip;
				// }
			// }
			
			
			
			
			vector<int> send_proc_numbering;
			for(int ip=0, numbering=0; ip<rank; ++ip){
				int str=mesh.displsSendProcFaces[ip];
				int end=str+mesh.countsSendProcFaces[ip];
				for(int i=str; i<end; ++i){
					send_proc_numbering.push_back(-1);
				}
			}
			for(int ip=rank, numbering=0; ip<size; ++ip){
				int str=mesh.displsSendProcFaces[ip];
				int end=str+mesh.countsSendProcFaces[ip];
				for(int i=str; i<end; ++i){
					send_proc_numbering.push_back(str_proc_numbering[rank] + i - mesh.displsSendProcFaces[rank]);
				}
			}
			vector<int> recv_proc_numbering;
			recv_proc_numbering.resize(send_proc_numbering.size());
			MPI_Alltoallv( send_proc_numbering.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
						   recv_proc_numbering.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
						   MPI_COMM_WORLD);
				
			
			for(int ip=0, numbering=0; ip<rank; ++ip){
				int str=mesh.displsRecvProcFaces[ip];
				int end=str+mesh.countsRecvProcFaces[ip];
				for(int i=str; i<end; ++i){
					proc_numbering.push_back(recv_proc_numbering[i]);
				}
			}
			for(int ip=rank, numbering=0; ip<size; ++ip){
				int str=mesh.displsRecvProcFaces[ip];
				int end=str+mesh.countsRecvProcFaces[ip];
				for(int i=str; i<end; ++i){
					proc_numbering.push_back(send_proc_numbering[i]);
				}
			}
			 


		
	// MPI_Barrier(MPI_COMM_WORLD);
			// for(int ip=0, numbering=0; ip<size; ++ip){
				// int str=mesh.displsProcFaces[ip];
				// int end=str+mesh.countsProcFaces[ip];
				// for(int i=str; i<end; ++i){
					// // if(rank==1) cout << rank << " " << proc_numbering[i] << endl;
					// cout << rank << " " << proc_numbering[i] << endl;
				// }
			// }
					
			// // mesh.displsProcFaces[rank]
		}
	
		
	
	
		
		
		
		
	
		
		vector<vector<int>> send_localFace_IN2IN_iL(size);
		vector<vector<int>> send_localFace_IN2IN_iR(size);
		vector<vector<int>> send_localFace_IN2IN_ipoints(size);
		
		vector<vector<int>> send_localFace_IN2PR_iL(size);
		vector<vector<int>> send_globalFace_IN2PR_toProc(size);
		vector<vector<int>> send_localFace_IN2PR_ipoints(size);
		vector<vector<int>> send_localFace_IN2PR_BoolReorder(size);
		
		vector<vector<int>> send_localFace_BC2BC_iL(size);
		vector<vector<int>> send_localFace_BC2BC_BCType(size);
		vector<vector<int>> send_localFace_BC2BC_ipoints(size);
		
		vector<vector<int>> send_localFace_PR2IN_iL(size);
		vector<vector<int>> send_localFace_PR2IN_toProc(size);
		vector<vector<int>> send_localFace_PR2IN_ipoints(size);
		vector<vector<int>> send_localFace_PR2IN_BoolReorder(size);
		
		vector<vector<int>> send_localFace_PR2PR_procFaceNum(size);
		vector<vector<int>> send_localFace_PR2PR_iL(size);
		vector<vector<int>> send_globalFace_PR2PR_toProc(size);
		vector<vector<int>> send_localFace_PR2PR_ipoints(size);
		vector<vector<int>> send_localFace_PR2PR_BoolReorder(size);
		
		
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL;
			int iR = face.iR;
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				auto& [ipR, iR_local] = send_localCell_proc_id[iR];
				if(ipL==ipR){
				// IN2IN
					send_localFace_IN2IN_iL[ipL].push_back(iL_local);
					send_localFace_IN2IN_iR[ipL].push_back(iR_local);
					send_localFace_IN2IN_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l = 0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) ipoint_local = ipoint;
							if(proc==ipL) ++l;
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						send_localFace_IN2IN_ipoints[ipL].push_back(ipoint_local);
					}
					
				}
				else{
				// IN2PR
					send_localFace_IN2PR_iL[ipL].push_back(iL_local);
					send_globalFace_IN2PR_toProc[ipL].push_back(ipR);
					send_localFace_IN2PR_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 53" << endl;
						send_localFace_IN2PR_ipoints[ipL].push_back(ipoint_local);
					}
					
					send_localFace_IN2PR_iL[ipR].push_back(iR_local);
					send_globalFace_IN2PR_toProc[ipR].push_back(ipL);
					vector<int> tmp_ipoints = face.ipoints;
					std::reverse(tmp_ipoints.begin()+1, tmp_ipoints.end());
					send_localFace_IN2PR_ipoints[ipR].push_back(face.ipoints.size());
					for(auto& item : tmp_ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipR) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 54" << endl;
						send_localFace_IN2PR_ipoints[ipR].push_back(ipoint_local);
					}
					
					
					
					//*************************
					if(ipL<ipR){
						// ipL, ipR 모두 페이스 리오더
						send_localFace_IN2PR_BoolReorder[ipL].push_back(1);
						send_localFace_IN2PR_BoolReorder[ipR].push_back(1);
					}
					else{
						send_localFace_IN2PR_BoolReorder[ipL].push_back(0);
						send_localFace_IN2PR_BoolReorder[ipR].push_back(0);
					}
					//*************************
					
					
					
				}
			}
			else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				// BC2BC
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				
				send_localFace_BC2BC_iL[ipL].push_back(iL_local);
				send_localFace_BC2BC_BCType[ipL].push_back(boundary_type[i]);
				send_localFace_BC2BC_ipoints[ipL].push_back(face.ipoints.size());
				for(auto& item : face.ipoints) {
					int ipoint_local = -1;
					for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
						if(proc==ipL) ipoint_local = ipoint;
					}
					if(ipoint_local==-1) cout << "ERROR" << endl;
					send_localFace_BC2BC_ipoints[ipL].push_back(ipoint_local);
				}
				
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				int ipR = recv_idBlockCell[ip];
				int rightProcNo = recv_rank[ip];
				
				if(ipL==ipR){
				// PR2IN
					send_localFace_PR2IN_iL[ipL].push_back(iL_local);
					send_localFace_PR2IN_toProc[ipL].push_back(rightProcNo);
					send_localFace_PR2IN_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) ipoint_local = ipoint;
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						send_localFace_PR2IN_ipoints[ipL].push_back(ipoint_local);
					}
					
					
					//************************
					if( rank<rightProcNo ){
						// 페이스 리오더
						send_localFace_PR2IN_BoolReorder[ipL].push_back(1);
					}
					else{
						send_localFace_PR2IN_BoolReorder[ipL].push_back(0);
					}
					//************************
					
				}
				else{
				// PR2PR
					send_localFace_PR2PR_procFaceNum[ipL].push_back(proc_numbering[ip]);
					
					send_localFace_PR2PR_iL[ipL].push_back(iL_local);
					send_globalFace_PR2PR_toProc[ipL].push_back(ipR);
					send_localFace_PR2PR_ipoints[ipL].push_back(face.ipoints.size());
					
					double avgx = 0.0;
					double avgy = 0.0;
					double avgz = 0.0;
					int tmp_size = face.ipoints.size();
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 55" << endl;
						// cout << l << endl;
						send_localFace_PR2PR_ipoints[ipL].push_back(ipoint_local);
						
						avgx += mesh.points[item].x/(double)tmp_size;
						avgy += mesh.points[item].y/(double)tmp_size;
						avgz += mesh.points[item].z/(double)tmp_size;
					}
					
					
					
					
					//************************
					if(
					(rank<rightProcNo && ipL<ipR) ||
					(rank>rightProcNo && ipL>ipR) ){
						send_localFace_PR2PR_BoolReorder[ipL].push_back(0);
					}
					else{
						// ipL, ipR 모두 페이스 리오더
						send_localFace_PR2PR_BoolReorder[ipL].push_back(1);
					}
					//************************
					
				}
				++ip;
			}
		}
		
		
		
		// //=========================================
		
		vector<vector<int>> recv_localFace_IN2IN_iL;
		vector<vector<int>> recv_localFace_IN2IN_iR;
		vector<vector<int>> recv_localFace_IN2IN_ipoints;
		vector<vector<int>> recv_localFace_IN2PR_iL;
		vector<vector<int>> recv_globalFace_IN2PR_toProc;
		vector<vector<int>> recv_localFace_IN2PR_ipoints;
		vector<vector<int>> recv_localFace_BC2BC_iL;
		vector<vector<int>> recv_localFace_BC2BC_BCType;
		vector<vector<int>> recv_localFace_BC2BC_ipoints;
		vector<vector<int>> recv_localFace_PR2IN_iL;
		vector<vector<int>> recv_localFace_PR2IN_toProc;
		vector<vector<int>> recv_localFace_PR2IN_ipoints;
		vector<vector<int>> recv_localFace_PR2PR_procFaceNum;
		vector<vector<int>> recv_localFace_PR2PR_iL;
		vector<vector<int>> recv_globalFace_PR2PR_toProc;
		vector<vector<int>> recv_localFace_PR2PR_ipoints;
		vector<vector<int>> recv_localFace_IN2PR_BoolReorder;
		vector<vector<int>> recv_localFace_PR2PR_BoolReorder;
		vector<vector<int>> recv_localFace_PR2IN_BoolReorder;
		
		
		mpi.Alltoallv(send_localFace_IN2IN_iL, recv_localFace_IN2IN_iL);
		mpi.Alltoallv(send_localFace_IN2IN_iR, recv_localFace_IN2IN_iR);
		mpi.Alltoallv(send_localFace_IN2IN_ipoints, recv_localFace_IN2IN_ipoints);
		
		mpi.Alltoallv(send_localFace_IN2PR_iL, recv_localFace_IN2PR_iL);
		mpi.Alltoallv(send_globalFace_IN2PR_toProc, recv_globalFace_IN2PR_toProc);
		mpi.Alltoallv(send_localFace_IN2PR_ipoints, recv_localFace_IN2PR_ipoints);
		mpi.Alltoallv(send_localFace_IN2PR_BoolReorder, recv_localFace_IN2PR_BoolReorder);
		
		mpi.Alltoallv(send_localFace_BC2BC_iL, recv_localFace_BC2BC_iL);
		mpi.Alltoallv(send_localFace_BC2BC_BCType, recv_localFace_BC2BC_BCType);
		mpi.Alltoallv(send_localFace_BC2BC_ipoints, recv_localFace_BC2BC_ipoints);
		
		mpi.Alltoallv(send_localFace_PR2IN_iL, recv_localFace_PR2IN_iL);
		mpi.Alltoallv(send_localFace_PR2IN_toProc, recv_localFace_PR2IN_toProc);
		mpi.Alltoallv(send_localFace_PR2IN_ipoints, recv_localFace_PR2IN_ipoints);
		mpi.Alltoallv(send_localFace_PR2IN_BoolReorder, recv_localFace_PR2IN_BoolReorder);
		
		mpi.Alltoallv(send_localFace_PR2PR_procFaceNum, recv_localFace_PR2PR_procFaceNum);
		mpi.Alltoallv(send_localFace_PR2PR_iL, recv_localFace_PR2PR_iL);
		mpi.Alltoallv(send_globalFace_PR2PR_toProc, recv_globalFace_PR2PR_toProc);
		mpi.Alltoallv(send_localFace_PR2PR_ipoints, recv_localFace_PR2PR_ipoints);
		mpi.Alltoallv(send_localFace_PR2PR_BoolReorder, recv_localFace_PR2PR_BoolReorder);
		
		
		
		
		
		
		// IN2IN
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_localFace_IN2IN_iL[ip].size(); ++i){
				
				face_state.push_back(meshComb.faces.size());
				
				meshComb.addFace();
				meshComb.faces.back().setType(MASCH_Face_Types::INTERNAL);
				{
					int id_local = recv_localFace_IN2IN_iL[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iL = id_global;
				}
				{
					int id_local = recv_localFace_IN2IN_iR[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iR = id_global;
				}
				int tmp_size = recv_localFace_IN2IN_ipoints[ip].at(iter++);
				for(int j=0; j<tmp_size; ++j){
					int ipoint_local = recv_localFace_IN2IN_ipoints[ip].at(iter++);
					int ipoint_global = points_id_local2global[ip].at(ipoint_local);
					meshComb.faces.back().ipoints.push_back(ipoint_global);
				}
			}
		}
		
		
		// PR2IN
		{
			int PR2IN_str = meshComb.faces.size();
			vector<int> tmp_nFaces_local(size,0);
			vector<int> tmp_strFaces_local(size+1,0);
			for(int ip=0; ip<size; ++ip){
				tmp_strFaces_local[ip+1] = tmp_strFaces_local[ip];
				for(int i=0, iter=0; i<recv_localFace_PR2IN_iL[ip].size(); ++i){
					int rightProc_local = recv_localFace_PR2IN_toProc[ip].at(i);
					
					
					bool boolReorder = (recv_localFace_PR2IN_BoolReorder[ip].at(i)==1);
					
					
					
					// 중복
					if(rightProc_local<ip){
						int tmp_id_global = PR2IN_str + 
											tmp_strFaces_local[rightProc_local] +
											tmp_nFaces_local.at(rightProc_local);
						++tmp_nFaces_local.at(rightProc_local);
						// ++tmp_nFaces_local.at(ip);
						auto& rightFace = meshComb.faces.at(tmp_id_global);
						int id_local = recv_localFace_PR2IN_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						
						if(boolReorder==true){
							rightFace.iR = id_global;
						}
						else{
							rightFace.iL = id_global;
						}
						
						int tmp_size = recv_localFace_PR2IN_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j) iter++;
					}
					else{
				
				face_state.push_back(PR2IN);
						
						++tmp_strFaces_local[ip+1];
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::INTERNAL);
						int id_local = recv_localFace_PR2IN_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						
						if(boolReorder==true){
							meshComb.faces.back().iR = id_global;
						}
						else{
							meshComb.faces.back().iL = id_global;
						}
						
						int tmp_size = recv_localFace_PR2IN_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j){
							int ipoint = recv_localFace_PR2IN_ipoints[ip].at(iter++);
							int tmp_id_global = points_id_local2global[ip].at(ipoint);
							meshComb.faces.back().ipoints.push_back(tmp_id_global);
						}
						
						if(boolReorder==true){
							auto& ipoint_ptr = meshComb.faces.back().ipoints;
							reverse(ipoint_ptr.begin()+1,ipoint_ptr.end());
						}
						
						
						
					}
				}
			}
		}
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
	
	
		
		// BC2BC
		if(size>1){
			int nbc_glo;
			MPI_Allreduce(&nbc, &nbc_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			nbc = nbc_glo;
		}
		
		
		vector<int> nFaces_boundary(nbc,0);
		
		vector<vector<pair<int,int>>> reorder_procFace_BC2BC_proc_id(nbc);
		vector<vector<int>> reorder_procFace_BC2BC_strId(nbc);
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_localFace_BC2BC_BCType[ip].size(); ++i){
				int ibc = recv_localFace_BC2BC_BCType[ip][i];
				// if(ibc>=nbc) cout << ibc << " " << nbc << endl;
				reorder_procFace_BC2BC_proc_id.at(ibc).push_back(make_pair(ip,i));
				reorder_procFace_BC2BC_strId[ibc].push_back(iter);
				int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(iter++);
				for(int j=0; j<tmp_size; ++j) iter++;
			}
		}
		
	
			
		vector<int> iter_BC2BC(size,0);
		for(int ibc=0; ibc<nbc; ++ibc){
			int str_size = meshComb.faces.size();
			int iter=0;
			for(auto& [ip, i] : reorder_procFace_BC2BC_proc_id[ibc]){
				
				
				face_state.push_back(BC2BC);
				
				meshComb.addFace();
				meshComb.faces.back().setType(MASCH_Face_Types::BOUNDARY);
				{
					int id_local = recv_localFace_BC2BC_iL[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iL = id_global;
				}
				int strId = reorder_procFace_BC2BC_strId[ibc][iter];
				int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(strId++);
				for(int j=0; j<tmp_size; ++j){
					int ipoint = recv_localFace_BC2BC_ipoints[ip].at(strId++);
					int id_global = points_id_local2global[ip].at(ipoint);
					meshComb.faces.back().ipoints.push_back(id_global);
				}
						
				// int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
				// for(int j=0; j<tmp_size; ++j){
					// int ipoint = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
					// int id_global = points_id_local2global[ip].at(ipoint);
					// meshComb.faces.back().ipoints.push_back(id_global);
				// }
				++iter;
			}
			nFaces_boundary[ibc] = meshComb.faces.size() - str_size;
		}
		
		
		
		
		// //************************************
		// // 프로세서 페이스 그룹핑 하기
		// vector<vector<int>> groupProcFaces;
		// {
			
			
			
			
			
			
			
		// }
		// //************************************
		
		
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
		
		// 프로세서
		vector<int> nFaces_processor(size,0);
		{
			// IN2PR
			vector<vector<pair<int,int>>> reorder_procFace_IN2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_IN2PR_strId(size);
			// for(int ip=size-1; ip>=0; --ip){
			for(int ip=0; ip<size; ++ip){
				for(int i=0, iter=0; i<recv_localFace_IN2PR_iL[ip].size(); ++i){
					int rightProc_global = recv_globalFace_IN2PR_toProc[ip][i];
					
					int boolReorder = recv_localFace_IN2PR_BoolReorder[ip][i];
					
					// if(rightProc_global<rank)
					{
						reorder_procFace_IN2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
						reorder_procFace_IN2PR_strId[rightProc_global].push_back(iter);
						int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j) iter++;
					}
				}
			}
			
			// PR2PR
			vector<vector<int>> reorder_procFace_PR2PR_procFaceNum(size);
			vector<vector<pair<int,int>>> reorder_procFace_PR2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_PR2PR_strId(size);
			vector<vector<vector<int>>> reorder_procFace_PR2PR_ipoints(size);
			for(int ip=0; ip<size; ++ip){
				for(int i=0, iter=0, iter2=0; i<recv_localFace_PR2PR_iL[ip].size(); ++i){
					// if(rank==0) cout << recv_localFace_PR2PR_iL[ip][i] << " " << i << endl;
					int procFaceNum = recv_localFace_PR2PR_procFaceNum[ip].at(i);
					
					int rightProc_global = recv_globalFace_PR2PR_toProc[ip].at(i);
					
					int boolReorder = recv_localFace_PR2PR_BoolReorder[ip][i];
					
					// if(rightProc_global!=rank)
					{
						reorder_procFace_PR2PR_procFaceNum[rightProc_global].push_back(procFaceNum);
						reorder_procFace_PR2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
						reorder_procFace_PR2PR_strId[rightProc_global].push_back(iter);
						int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(iter++);
						double avgx = 0.0;
						double avgy = 0.0;
						double avgz = 0.0;
						vector<int> tmp;
						for(int j=0; j<tmp_size; ++j) {
							int id_local = recv_localFace_PR2PR_ipoints[ip].at(iter);
							int id_global = points_id_local2global[ip].at(id_local);
							tmp.push_back(id_global);
							avgx += meshComb.points[id_global].x/(double)tmp_size;
							avgy += meshComb.points[id_global].y/(double)tmp_size;
							avgz += meshComb.points[id_global].z/(double)tmp_size;
							iter++;
						}
						reorder_procFace_PR2PR_ipoints[rightProc_global].push_back(tmp);
						
						// if(rank==0 && rightProc_global==1){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
						// if(rank==1 && rightProc_global==0){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
					}
				}
			}
			
	
			
			// face 번호로 리오더
			for(int ip=0; ip<size; ++ip){
				int tmp_size = reorder_procFace_PR2PR_procFaceNum[ip].size();
				vector<int> v_procFaceNum = reorder_procFace_PR2PR_procFaceNum[ip];
				sort(v_procFaceNum.begin(), v_procFaceNum.end());
				// reverse(v_procFaceNum.begin(), v_procFaceNum.end());
				
				vector<pair<int,int>> tmp_v1(tmp_size);
				vector<int> tmp_v2(tmp_size);
				vector<vector<int>> tmp_v3(tmp_size);
				for(int i=0; i<tmp_size; ++i){
					int procFaceNum = reorder_procFace_PR2PR_procFaceNum[ip].at(i);
					// if(rank==1) cout << rank << " " << ip << " " << v_procFaceNum[i] << endl;
					// if(rank==1) cout << ip << " " << v_procFaceNum[i] << " " << procFaceNum << endl;
					auto it = find (v_procFaceNum.begin(), v_procFaceNum.end(), procFaceNum);
					
					int order = (it - v_procFaceNum.begin());
					// if(rank==1) cout << ip << " " << order << endl;
					
					tmp_v1[order] = (reorder_procFace_PR2PR_proc_id[ip][i]);
					tmp_v2[order] = (reorder_procFace_PR2PR_strId[ip][i]);
					tmp_v3[order] = (reorder_procFace_PR2PR_ipoints[ip][i]);
					
				}
				
				reorder_procFace_PR2PR_procFaceNum[ip].clear();
				reorder_procFace_PR2PR_proc_id[ip].clear();
				reorder_procFace_PR2PR_strId[ip].clear();
				reorder_procFace_PR2PR_ipoints[ip].clear();
				for(int i=0; i<tmp_size; ++i){
					reorder_procFace_PR2PR_proc_id[ip].push_back(tmp_v1[i]);
					reorder_procFace_PR2PR_strId[ip].push_back(tmp_v2[i]);
					reorder_procFace_PR2PR_procFaceNum[ip].push_back(v_procFaceNum[i]);
					reorder_procFace_PR2PR_ipoints[ip].push_back(tmp_v3[i]);
				}
				
				
			}
			
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
			
			
			
			
			
		
			vector<int> iter_IN2PR(size,0);
			vector<int> iter_PR2PR(size,0);
			for(int proc=0; proc<size; ++proc){
				int str_size = meshComb.faces.size();
				// IN2PR
				{
					int iter=0;
					// cout << reorder_procFace_IN2PR_proc_id[proc].size() << 
					for(auto& [ip, i] : reorder_procFace_IN2PR_proc_id[proc]){
				
				face_state.push_back(IN2PR);
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						{
							int id_local = recv_localFace_IN2PR_iL[ip].at(i);
							int id_global = nCells_local[ip] + id_local;
							meshComb.faces.back().iL = id_global;
						}
						
						
						
						//*************************
						int boolReorder = recv_localFace_IN2PR_BoolReorder[ip][i];
						procFace_boolReorder.push_back( (boolReorder==1) );
						//*************************
						
						
						
						// cout << recv_localFace_IN2PR_ipoints[ip].size() << " " << iter_IN2PR[ip] << endl;
						// int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(iter_IN2PR[ip]++);
						int strId = reorder_procFace_IN2PR_strId[proc].at(iter);
						int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(strId++);
						for(int j=0; j<tmp_size; ++j){
							int ipoint = recv_localFace_IN2PR_ipoints[ip].at(strId++);
							int id_global = points_id_local2global[ip].at(ipoint);
							meshComb.faces.back().ipoints.push_back(id_global);
						}
						++iter;
					}
				}
				// PR2PR
				{
					
					int iter=0;
					int tmp2_size = reorder_procFace_PR2PR_proc_id[proc].size();
					for(int ii=0; ii<tmp2_size; ++ii){
						
				face_state.push_back(PR2PR);
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						
						auto& [ip, i] = reorder_procFace_PR2PR_proc_id[proc][ii];
						int id_local = recv_localFace_PR2PR_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						meshComb.faces.back().iL = id_global;
						
						
						
						//*************************
						int boolReorder = recv_localFace_PR2PR_BoolReorder[ip][i];
						procFace_boolReorder.push_back( (boolReorder==1) );
						//*************************
						
						
						
						int strId = reorder_procFace_PR2PR_strId[proc][ii];
						int procFaceNum = reorder_procFace_PR2PR_procFaceNum[proc][ii];
						vector<int> ipoints = reorder_procFace_PR2PR_ipoints[proc][ii];
				
						double avgx = 0.0;
						double avgy = 0.0;
						double avgz = 0.0;
						// int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						int tmp_size = ipoints.size();
						for(auto& ipoint : ipoints){
							// int id_loc = recv_localFace_PR2PR_ipoints[ip].at(strId++);
							// int id_glo = points_id_local2global[ip].at(id_loc);
							meshComb.faces.back().ipoints.push_back(ipoint);
							avgx += meshComb.points[ipoint].x/(double)tmp_size;
							avgy += meshComb.points[ipoint].y/(double)tmp_size;
							avgz += meshComb.points[ipoint].z/(double)tmp_size;
						}
						
						// if(rank==0 && proc==1){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
						// if(rank==1 && proc==0){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
					}
					
					
					// int iter=0;
					// for(auto& [ip, i] : reorder_procFace_PR2PR_proc_id[proc]){
				
				// face_state.push_back(PR2PR);
				
						// meshComb.addFace();
						// meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						// {
							// int id_local = recv_localFace_PR2PR_iL[ip].at(i);
							// int id_global = nCells_local[ip] + id_local;
							// meshComb.faces.back().iL = id_global;
						// }
						// int procFaceNum = reorder_procFace_PR2PR_procFaceNum[proc].at(iter);
						// double avgx = 0.0;
						// double avgy = 0.0;
						// double avgz = 0.0;
						// int strId = reorder_procFace_PR2PR_strId[proc].at(iter);
						// int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						// for(int j=0; j<tmp_size; ++j){
							// int ipoint = recv_localFace_PR2PR_ipoints[ip].at(strId++);
							// int id_global = points_id_local2global[ip].at(ipoint);
							// meshComb.faces.back().ipoints.push_back(id_global);
							// avgx += meshComb.points[id_global].x/(double)tmp_size;
							// avgy += meshComb.points[id_global].y/(double)tmp_size;
							// avgz += meshComb.points[id_global].z/(double)tmp_size;
						// }
						// ++iter;
						
						
						
						// // if(rank==0 && proc==1){
							// // cout << rank << " " << procFaceNum << endl;
							// // cout << avgx << " " << avgy << " " << avgz << endl;
							
						// // }
						// // if(rank==1 && proc==0){
							// // cout << rank << " " << procFaceNum << endl;
							// // cout << avgx << " " << avgy << " " << avgz << endl;
							
						// // }
					// }
				}
				nFaces_processor[proc] = meshComb.faces.size() - str_size;
			}
			
		}
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		// 바운더리 네임 뿌려주기
		vector<string> bc_names;
		if(rank == 0)
		{
			vector<string> send_bc_names;
			for(int ibc=0; ibc<nbc; ++ibc){
				auto& boundary = mesh.boundaries[ibc];
				string tmp_name = boundary.name;
				int leng = tmp_name.length();
				MPI_Bcast(&leng, 1, MPI_INT, 0, MPI_COMM_WORLD);
				char *buf = new char[leng];
				strcpy(buf,tmp_name.c_str());
				MPI_Bcast(buf, leng, MPI_CHAR, 0, MPI_COMM_WORLD);
				// MPI_Scatter(tmp_name.c_str(), leng_glo, MPI_CHAR, buf, leng_glo, MPI_CHAR, 0, MPI_COMM_WORLD);
				// string bla1(buf, leng_glo);
				// delete [] buf;
				bc_names.push_back(tmp_name);
			}
		}
		else
		{
			for(int ibc=0; ibc<nbc; ++ibc){
				int leng;
				MPI_Bcast(&leng, 1, MPI_INT, 0, MPI_COMM_WORLD);
				char *buf = new char[leng];
				MPI_Bcast(buf, leng, MPI_CHAR, 0, MPI_COMM_WORLD);
				string bla1(buf, leng);
				delete [] buf;
				bc_names.push_back(bla1);
			}
		}
		
		
		
		// 바운더리
		if(mesh.boundaries.size()<nbc){
			// for(auto& item : bc_names){
				// cout << item << endl;
			// }
			mesh.boundaries.resize(nbc);
		}
		meshComb.boundaries.clear();
		for(int ibc=0; ibc<nbc; ++ibc){
			// if(ibc>=mesh.boundaries.size()) break;
			auto& boundary = mesh.boundaries[ibc];
			meshComb.addBoundary();
			
			string bcName = bc_names[ibc];
			// string bcName = boundary.name;
			
			bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			meshComb.boundaries.back().name = bcName;
			meshComb.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			meshComb.boundaries.back().nFaces = nFaces_boundary[ibc];
			
			// if(rank==60) cout << "GG " << nFaces_boundary[ibc] << endl;
		}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		for(int ip=0; ip<size; ++ip){
			if(nFaces_processor[ip]>0){
				string bcnames = "procBoundary" + to_string(rank) + "to" + to_string(ip);
				
				meshComb.addBoundary();
				meshComb.boundaries.back().name = bcnames;
				meshComb.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
				meshComb.boundaries.back().nFaces = nFaces_processor[ip];
				meshComb.boundaries.back().myProcNo = rank;
				meshComb.boundaries.back().rightProcNo = ip;
			}
		}
		
// // int nSizeOrg, int nSizeTar, 
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		// if(rank==60) cout << meshComb.boundaries.size() << endl;
	
		// if(meshComb.boundaries.size()>0){
			int maxBCNum = meshComb.boundaries.size()-1;
			meshComb.boundaries[maxBCNum].startFace = meshComb.faces.size()-meshComb.boundaries[maxBCNum].nFaces;
			for(int i=maxBCNum-1; i>=0; --i){
				meshComb.boundaries[i].startFace = meshComb.boundaries[i+1].startFace-meshComb.boundaries[i].nFaces;
			}
		// }
		
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// variables 넘기기 위한 재료
		{
			vector<vector<int>> send_nCells_local(size);
			for(int ip=0; ip<size; ++ip) send_nCells_local[ip].push_back(nCells_local[ip]);
			vector<vector<int>> recv_nCells_local(size);
			mpi.Alltoallv(send_nCells_local, recv_nCells_local);
			int iter = 0;
			for(auto& [proc, id] : send_localCell_proc_id){
				to_new_cell_id.at(iter++) = recv_nCells_local[proc][0] + id;
				// to_new_cell_id.at(iter++);
			}
		}
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
		
		//=======================================	
		// parcel 위치 정보 넘기기
		// vector<vector<int>> recv_localParcel_size;
		// vector<vector<int>> recv_localParcel_id;
		// vector<vector<int>> recv_localParcel_icell;
		{
			// vector<vector<int>> send_localParcel_size(size);
			vector<vector<int>> send_localParcel_id(size);
			vector<vector<int>> send_localParcel_icell(size);
			// int iter = 0;
			// idBlockCell
			for(auto& parcel : mesh.parcels){
				auto& [proc, id] = send_localCell_proc_id[parcel.icell];
				send_localParcel_id[proc].push_back(parcel.id);
				send_localParcel_icell[proc].push_back(id);
			}
			// for(auto& [proc, id] : send_localCell_proc_id){
				// auto& cell = mesh.cells[iter++];
				// send_localParcel_size[proc].push_back(cell.iparcels.size());
				// for(auto& iparcel : cell.iparcels){
					// auto& parcel = mesh.parcels[iparcel];
					// send_localParcel_id[proc].push_back(parcel.id);
					// send_localParcel_icell[proc].push_back(id);
				// }
			// }
			// mpi.Alltoallv(send_localParcel_size, recv_localParcel_size);
			mpi.Alltoallv(send_localParcel_id, recv_localParcel_id);
			mpi.Alltoallv(send_localParcel_icell, recv_localParcel_icell);
			
			
			// for(int ip=0, tmp_id=0; ip<size; ++ip){
				// for(auto& id : recv_localParcel_icell[ip]){
					// int id_glo = nCells_local[ip] + id;
					// mesh.cells[id_glo].
				// }
				// int str = nCells_local[ip];
				// int end = nCells_local[ip+1];
				// for(int i=str; i<end; ++i){
					// auto& cell = mesh.cells[i];
					// // cell.group = recv_localCell_group[ip][i-str];
					// cell.level = recv_localCell_level[ip][i-str];
				// }
			// }
		
		}
		//=======================================	
		
		
		
		
		
		
		
		
		
		
		//=======================================	
		
		
		
		
	
		
	}
	
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
	
	// 원래 메쉬에 넣기
	{
		mesh.cells.clear();
		mesh.cells.shrink_to_fit();
		
		mesh.faces.clear();
		mesh.faces.shrink_to_fit();
		
		mesh.points.clear();
		mesh.points.shrink_to_fit();
		
		mesh.boundaries.clear();
		mesh.boundaries.shrink_to_fit();

		mesh.cells.reserve(meshComb.cells.size());
		mesh.faces.reserve(meshComb.faces.size());
		mesh.points.reserve(meshComb.points.size());
		
		meshComb.cells.clear();
		meshComb.cells.shrink_to_fit();
		
		for(auto& point : meshComb.points){
			mesh.addPoint();
			mesh.points.back().x = point.x;
			mesh.points.back().y = point.y;
			mesh.points.back().z = point.z;
			mesh.points.back().level = point.level;
			for(auto& [proc, id] : point.connPoints){
				mesh.points.back().connPoints.push_back(make_pair(proc, id));
			}
		}
		meshComb.points.clear();
		meshComb.points.shrink_to_fit();
		
		for(auto& face : meshComb.faces){
			mesh.addFace();
			mesh.faces.back().iL = face.iL;
			mesh.faces.back().iR = face.iR;
			mesh.faces.back().setType(face.getType());
			for(auto& ipoint : face.ipoints){
				mesh.faces.back().ipoints.push_back(ipoint);
			}
		}
		meshComb.faces.clear();
		meshComb.faces.shrink_to_fit();
		
		for(auto& boundary : meshComb.boundaries){
			mesh.addBoundary();
			mesh.boundaries.back().name = boundary.name;
			mesh.boundaries.back().startFace = boundary.startFace;
			mesh.boundaries.back().nFaces = boundary.nFaces;
			mesh.boundaries.back().rightProcNo = boundary.rightProcNo;
			mesh.boundaries.back().myProcNo = boundary.myProcNo;
			mesh.boundaries.back().setType(boundary.getType());
		}
		meshComb.boundaries.clear();
		meshComb.boundaries.shrink_to_fit();
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.buildCells();















	
		// {
			// double resi_normals = 0.7;
			
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
			// // vector<vector<vector<double>>> reorder_procFace_x(size);
			// // vector<vector<vector<double>>> reorder_procFace_y(size);
			// // vector<vector<vector<double>>> reorder_procFace_z(size);
			// for(int ip=0; ip<size; ++ip){
				// int size_procFace = tmp_procFace_group[ip].size();
					
				// vector<vector<int>> tmp_face_group;
				// // vector<vector<double>> tmp_face_group_x;
				// // vector<vector<double>> tmp_face_group_y;
				// // vector<vector<double>> tmp_face_group_z;
					
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
							// // vector<double> tmp_x;
							// // vector<double> tmp_y;
							// // vector<double> tmp_z;
							// tmp_face.push_back(tmp_id[j]);
							// // tmp_x.push_back(normal0x);
							// // tmp_y.push_back(normal0y);
							// // tmp_z.push_back(normal0z);
							
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
								// if(resi>resi_normals){
									// tmp_face.push_back(tmp_id[j]);
									// // tmp_x.push_back(normalx);
									// // tmp_y.push_back(normaly);
									// // tmp_z.push_back(normalz);
								// }
								// else{
									// break;
								// }
								
							// }
							// --j;
							
							// tmp_face_group.push_back(tmp_face);
							// // tmp_face_group_x.push_back(tmp_x);
							// // tmp_face_group_y.push_back(tmp_y);
							// // tmp_face_group_z.push_back(tmp_z);
						// }
					// }
				// }
					
				// reorder_procFace_id[ip] = (tmp_face_group);
				// // reorder_procFace_x[ip] = (tmp_face_group_x);
				// // reorder_procFace_y[ip] = (tmp_face_group_y);
				// // reorder_procFace_z[ip] = (tmp_face_group_z);
			// }
			
		// // MPI_Barrier(MPI_COMM_WORLD);
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
			
			// int str_proFace_id = -1;
			// for(auto& boundary : mesh.boundaries){
				// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
				// if(str_proFace_id==-1 && boundary.nFaces>0){
					// str_proFace_id = boundary.startFace;
					// break;
				// }
			// }
			
		// // MPI_Barrier(MPI_COMM_WORLD);
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			
			// for(int ip=0; ip<size; ++ip){
				// int size1 = reorder_procFace_id[ip].size();
				// for(int i=0; i<size1; ++i){
					// auto& procFace = reorder_procFace_id[ip][i];
					// int size2 = procFace.size();
					// if(size2<=1) continue;
					
					// int id0 = procFace[0];
					// if(procFace_boolReorder[id0-str_proFace_id]==false) continue;
					
					// vector<int> reorderFaceIds;
					// mesh.getFaceOrders(maxLevel, mesh.points, mesh.faces, procFace, reorderFaceIds);
					
					// vector<int> tmp_iL;
					// vector<vector<int>> tmp_ipoints;
					// // for(int j=size2-1; j>=0; --j){
					// for(int j=0; j<size2; ++j){
						// int id = reorderFaceIds[j];
						// tmp_iL.push_back(mesh.faces[id].iL);
						// tmp_ipoints.push_back(vector<int>());
						// for(auto& ipoint : mesh.faces[id].ipoints){
							// tmp_ipoints.back().push_back(ipoint);
						// }
					// }
						
					// for(int j=0; j<size2; ++j){
						// int id = procFace[j];
						// mesh.faces[id].iL = tmp_iL[j];
						// mesh.faces[id].ipoints.clear();
						// for(auto& ipoint : tmp_ipoints[j]){
							// mesh.faces[id].ipoints.push_back(ipoint);
						// }
					// }
				// }
				
			// }
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		
		
		
		
		
		
		
		
		
		
		
	// **************************
	// 셀 값
	{
		for(int ip=0, tmp_id=0; ip<size; ++ip){
			int str = nCells_local[ip];
			int end = nCells_local[ip+1];
			for(int i=str; i<end; ++i){
				auto& cell = mesh.cells[i];
				// cell.group = recv_localCell_group[ip][i-str];
				cell.level = recv_localCell_level[ip][i-str];
			}
		}
	}
	
	
	
	
	
	//=======================================	
	// parcel 위치 정보 넘기기
	// 파슬값 전달
	{
		// vector<int> parcel_ip(mesh.parcels.size());
		parcel_ip.resize(mesh.parcels.size());
		for(int i=0; i<mesh.parcels.size(); ++i){
			parcel_ip.at(i) = idBlockCell.at(mesh.parcels[i].icell);
		}
		
		mesh.parcels.clear();
		for(int ip=0; ip<size; ++ip){
			int iter=0;
			for(auto& id : recv_localParcel_icell[ip]){
				int id_glo = nCells_local[ip] + id;
				int parcel_id = recv_localParcel_id[ip][iter];
				mesh.addParcel(id_glo, parcel_id);
				++iter;
			}
		}
		
		// int org_nParcels = var.parcels.size();
		// int nParcels = mesh.parcels.size();
		// int nParcelsMax = max(nParcels,org_nParcels);
		// var.parcels.resize(nParcelsMax);
		// for(auto& item : var.parcels) item.resize(controls.nParcelVar);
		
		// int varSize = controls.nParcelVar;
		// for(int iprim=0; iprim<varSize; ++iprim){
			// vector<vector<double>> send_value(size);
			// for(int i=0; i<org_nParcels; ++i){
				// double org_value = var.parcels.at(i).at(iprim);
				// send_value.at(parcel_ip.at(i)).push_back(org_value);
			// }
			// vector<vector<double>> recv_value(size);
			// mpi.Alltoallv(send_value, recv_value);
			
			// for(int ip=0, iter=0; ip<size; ++ip){
				// for(auto& item : recv_value[ip]){
					// var.parcels.at(iter).at(iprim) = item;
					// ++iter;
				// }
			// }
		// }
		// var.parcels.resize(nParcels);
	}
	
	//=======================================
	
	
	
	
	// 그룹 번호 재정립 (그룹은 뭉쳐있음)
	int total_nGroup = 0;
	vector<int> str_icells(size+1,0);
	{
		for(int ip=0, tmp_id=0; ip<size; ++ip){
			int str = nCells_local[ip];
			int end = nCells_local[ip+1];
			for(int i=str; i<end; ++i){
				auto& cell0 = mesh.cells[i];
				int group0 = recv_localCell_group[ip][i-str];
				cell0.group = total_nGroup;
				while(1){
					++i;
					if(i==end) break;
					auto& cell = mesh.cells[i];
					int group = recv_localCell_group[ip][i-str];
					if(group0!=group) break;
					cell.group = total_nGroup;
				}
				--i;
				++total_nGroup;
			}
		}
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

		
		
		// int total_nGroup = 0;
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell0 = mesh.cells[i];
			// int group0 = cell0.group;
			// cell0.group = total_nGroup;
			// while(1){
				// auto& cell = mesh.cells[++i];
				// int group = cell.group;
				// if(group0!=group) break;
				// cell.group = total_nGroup;
			// }
			// --i;
			// ++total_nGroup;
		// }
		
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
			// for(int i=0, SIZE0=groupCellsProcFaces.size(); i<SIZE0; ++i){
				// int org_Rip0=-1;
				// if(groupCellsRightProcNo[i].size()!=0) 
					// org_Rip0 = groupCellsRightProcNo[i][0];
				// for(int j=0, SIZE1=groupCellsProcFaces[i].size(); j<SIZE1; ++j){
					// int org_id = groupCellsProcFaces[i][j];
					// int org_Rip = groupCellsRightProcNo[i][j];
					// if(org_Rip0==org_Rip) cout << i << " " << org_Rip << " " << org_id << endl;
				// }
			// }
		// }
		// //=======================================

	
		
		
		
		
	}
	// **************************
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
	mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.setFaceLevels();
	
	// mesh.setCellStencils();
	// mesh.setNumberOfFaces();
	// mesh.cellsGlobal();
	
	// mesh.debug_group_procFaces(0.99);

		
	//************************************
	// 프로세서 페이스 그룹핑 하기
	{
		// vector<vector<int>> groupCellsProcFaces(total_nGroup);
		// int cell_size = mesh.cells.size();
		// for(int i=0; i<cell_size; ++i){
			// auto& cell0 = mesh.cells[i];
			// int group0 = cell0.group;
		
			// for(auto& iface : cell0.ifaces){
				// auto& face = mesh.faces[iface];
				// if(face.getType()!=MASCH_Face_Types::PROCESSOR) continue;
				// int iGroup = group0-str_icells[rank];
				// // if(find(groupCellsProcFaces[iGroup].begin(),groupCellsProcFaces[iGroup].end(),
				// // iface)!=groupCellsProcFaces[iGroup].end()){
					// // cout << "ERROR123" << endl;
				// // }
				// groupCellsProcFaces[iGroup].push_back(iface);
			// }
		
			// while(1){
				// ++i;
				// if(i==cell_size) break;
				// auto& cell = mesh.cells[i];
				// int group = cell.group;
				// if(group0!=group) break;
				
				// for(auto& iface : cell.ifaces){
					// auto& face = mesh.faces[iface];
					// if(face.getType()!=MASCH_Face_Types::PROCESSOR) continue;
					// int iGroup = group-str_icells[rank];
					// // if(find(groupCellsProcFaces[iGroup].begin(),groupCellsProcFaces[iGroup].end(),
					// // iface)!=groupCellsProcFaces[iGroup].end()){
						// // cout << "ERROR123" << endl;
					// // }
					// groupCellsProcFaces[iGroup].push_back(iface);
				// }
			// }
			// --i;
		// }
		
		// vector<vector<int>> groupCellsProcFaces(total_nGroup);
		// vector<vector<int>> groupCellsRightProcNo(total_nGroup);
		// {
			// for(auto& boundary : mesh.boundaries){
				// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				// int rightProcNo = boundary.rightProcNo;
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					// int iL = face.iL;
					// int igroup = mesh.cells[iL].group-str_icells[rank];
					// groupCellsProcFaces[igroup].push_back(i);
					// groupCellsRightProcNo[igroup].push_back(rightProcNo);
					// // if(rank==0) cout << rightProcNo << " " << igroup << endl;
				// }
			// }
		// }
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		
		
		
		
		double resi_normals = 0.9;
		
		vector<vector<int>> tmp_procFace_id(size);
		vector<vector<int>> tmp_procFace_group(size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			int rightProcNo = boundary.rightProcNo;
			for(int i=str, iter=0; i<end; ++i){
				auto& face = mesh.faces[i];
				int iL = face.iL;
				// int igroup = mesh.cells[iL].group-min_group_id;
				int igroup = mesh.cells[iL].group;
				tmp_procFace_id.at(rightProcNo).push_back(i);
				tmp_procFace_group.at(rightProcNo).push_back(igroup);
			}
		}
	
		MASCH_Math math;
	
		
		vector<vector<vector<int>>> reorder_procFace_id(size);
		// vector<vector<vector<double>>> reorder_procFace_x(size);
		// vector<vector<vector<double>>> reorder_procFace_y(size);
		// vector<vector<vector<double>>> reorder_procFace_z(size);
		for(int ip=0; ip<size; ++ip){
			int size_procFace = tmp_procFace_group[ip].size();
				
			vector<vector<int>> tmp_face_group;
			for(int i=0; i<size_procFace; ++i){
				int id0 = tmp_procFace_id[ip].at(i);
				int igroup0 = tmp_procFace_group[ip].at(i);
				vector<int> tmp_id;
				tmp_id.push_back(id0);
				while(1){
					++i;
					if(i==size_procFace) break;
					int id = tmp_procFace_id[ip].at(i);
					int igroup = tmp_procFace_group[ip].at(i);
					if(igroup0!=igroup) break;
					tmp_id.push_back(id);
				}
				--i;
				{
					vector<vector<double>> procFace_unitNormals;
					for(auto& iface : tmp_id){
						auto& face = mesh.faces[iface];
						vector<double> Vx, Vy, Vz;
						for(auto& ipoint : face.ipoints){
							Vx.push_back(mesh.points.at(ipoint).x);
							Vy.push_back(mesh.points.at(ipoint).y);
							Vz.push_back(mesh.points.at(ipoint).z);
						}
						
						double VSn=0.0;
						vector<double> cellCentroid;
						vector<double> tmp_unitNormals(3,0.0);
						double area;
						math.calcUnitNormals_Area3dPolygon(
							face.ipoints.size(), Vx,Vy,Vz,
							tmp_unitNormals, area,
							face.x, face.y, face.z,
							VSn, cellCentroid);
							
						procFace_unitNormals.push_back(tmp_unitNormals);
					}
					
					int tmp_faceSize = procFace_unitNormals.size();
					for(int j=0; j<tmp_faceSize; ++j){
						double normal0x = procFace_unitNormals.at(j).at(0);
						double normal0y = procFace_unitNormals.at(j).at(1);
						double normal0z = procFace_unitNormals.at(j).at(2);
						
						vector<int> tmp_face;
						tmp_face.push_back(tmp_id.at(j));
						
						while(1){
							++j;
							if(j==tmp_faceSize) break;
							double normalx = procFace_unitNormals.at(j)[0];
							double normaly = procFace_unitNormals.at(j)[1];
							double normalz = procFace_unitNormals.at(j)[2];
							
							double resi = 0.0;
							resi += normalx*normal0x;
							resi += normaly*normal0y;
							resi += normalz*normal0z;
							if(resi>resi_normals){
								tmp_face.push_back(tmp_id.at(j));
								// tmp_x.push_back(normalx);
								// tmp_y.push_back(normaly);
								// tmp_z.push_back(normalz);
							}
							else{
								break;
							}
							
						}
						--j;
						
						tmp_face_group.push_back(tmp_face);
						// tmp_face_group_x.push_back(tmp_x);
						// tmp_face_group_y.push_back(tmp_y);
						// tmp_face_group_z.push_back(tmp_z);
					}
				}
			}
				
			reorder_procFace_id[ip] = (tmp_face_group);
			// reorder_procFace_x[ip] = (tmp_face_group_x);
			// reorder_procFace_y[ip] = (tmp_face_group_y);
			// reorder_procFace_z[ip] = (tmp_face_group_z);
		}
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		
		int str_proFace_id = -1;
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			if(str_proFace_id==-1 && boundary.nFaces>0){
				str_proFace_id = boundary.startFace;
				break;
			}
		}
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		// if(rank==5) cout << "AAAAAAAA " << rank << endl;
		
		for(int ip=0; ip<size; ++ip){
			int size1 = reorder_procFace_id[ip].size();
			for(int i=0; i<size1; ++i){
				auto& procFace = reorder_procFace_id[ip].at(i);
				int size2 = procFace.size();
				if(size2<=1) continue;
				
				int id0 = procFace.at(0);
				if(procFace_boolReorder.at(id0-str_proFace_id)==false) continue;
				
				vector<int> reorderFaceIds;
				// if(rank==5) cout << procFace.size() << endl;
				mesh.getFaceOrders(maxLevel, mesh.points, mesh.faces, procFace, reorderFaceIds);
				// if(rank==5) cout << procFace.size() << endl;
				
				vector<int> tmp_iL;
				vector<vector<int>> tmp_ipoints;
				// for(int j=size2-1; j>=0; --j){
				for(int j=0; j<size2; ++j){
					int id = reorderFaceIds.at(j);
					tmp_iL.push_back(mesh.faces.at(id).iL);
					tmp_ipoints.push_back(vector<int>());
					for(auto& ipoint : mesh.faces[id].ipoints){
						tmp_ipoints.back().push_back(ipoint);
					}
				}
					
				for(int j=0; j<size2; ++j){
					int id = procFace[j];
					mesh.faces[id].iL = tmp_iL[j];
					mesh.faces[id].ipoints.clear();
					for(auto& ipoint : tmp_ipoints[j]){
						mesh.faces[id].ipoints.push_back(ipoint);
					}
				}
			}
			
		}
		
		// if(rank==50 || rank==5) cout << "BBBBB " << rank << endl;
		
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.setFaceLevels();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	mesh.cellsGlobal();
	
	//************************************

		
	
	// proc 페이스 매칭 디버깅
	vector<vector<int>> send_procNFaces(size,vector<int>(1,0));
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
		int rightProcNo = boundary.rightProcNo;
		send_procNFaces[rightProcNo][0] = boundary.nFaces;
	}
	
	vector<vector<int>> recv_procNFaces;
	mpi.Alltoallv(send_procNFaces, recv_procNFaces);
	
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
		int rightProcNo = boundary.rightProcNo;
		if(recv_procNFaces[rightProcNo][0] != boundary.nFaces){
			cout << "#WARNING : not match proc faces" << endl;
		}
	}
	
	// // mesh.informations();
	
	// ===========================================
	// connPoints 디버깅
	mesh.debug_connPoints(1.e-5);
	// proc face points 디버깅
	mesh.debug_procFacePoints(1.e-5);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// ===========================================
	
	
	// // 디버깅
	// {
		// vector<int> send_idBlockCell;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// send_idBlockCell.push_back(idBlockCell[face.iL]);
			// }
		// }
		// recv_idBlockCell.resize(send_idBlockCell.size());
		// MPI_Alltoallv( send_idBlockCell.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   // recv_idBlockCell.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   // MPI_COMM_WORLD);
	// }
	
	
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
}

