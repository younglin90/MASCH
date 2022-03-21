
#include "./solvers.h"


void MASCH_Solver::calcCurvature(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
vector<string>& inp_cell){
	
	auto& solver = (*this);
	
	// cout << inp_cell.size() << endl;
	
	if(inp_cell.size()==0) return;
	
	
	// controls.log.push("test1");
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int inp_size = inp_cell.size();
	int proc_size = var.procRightCells.size();
	
	int id_area = controls.getId_faceVar("area");

	int id_inp[inp_size], id_curv_oup[inp_size], id_levelSet_inp[inp_size], 
	id_gradx_inp[inp_size], id_grady_inp[inp_size], 
	id_gradz_inp[inp_size], id_bc_face[inp_size],
	unitNormalId_nx[inp_size], unitNormalId_ny[inp_size], unitNormalId_nz[inp_size],
	gradx_unitNormalId_x[inp_size], grady_unitNormalId_y[inp_size], gradz_unitNormalId_z[inp_size];
	vector<string> unitNormalString, levelSetString;
	{
		int tmp_iter = 0;
		for(auto& item : inp_cell){
			id_inp[tmp_iter] = (controls.getId_cellVar("volume-fraction-"+item));
			id_curv_oup[tmp_iter] = (controls.getId_cellVar("curvature-"+item));
			id_levelSet_inp[tmp_iter] = (controls.getId_cellVar("level-set-"+item));
			levelSetString.push_back("level-set-"+item);
			
			// cout << levelSetString.back() << endl;
			
			id_gradx_inp[tmp_iter] = (controls.getId_cellVar("x-gradient level-set-"+item));
			id_grady_inp[tmp_iter] = (controls.getId_cellVar("y-gradient level-set-"+item));
			id_gradz_inp[tmp_iter] = (controls.getId_cellVar("z-gradient level-set-"+item));
			
			unitNormalId_nx[tmp_iter] = (controls.getId_cellVar("x-unit-normal-"+item));
			unitNormalId_ny[tmp_iter] = (controls.getId_cellVar("y-unit-normal-"+item));
			unitNormalId_nz[tmp_iter] = (controls.getId_cellVar("z-unit-normal-"+item));
			
			unitNormalString.push_back("x-unit-normal-"+item);
			unitNormalString.push_back("y-unit-normal-"+item);
			unitNormalString.push_back("z-unit-normal-"+item);
			
			gradx_unitNormalId_x[tmp_iter] = (controls.getId_cellVar("x-gradient x-unit-normal-"+item));
			// grady_unitNormalId_x[tmp_iter] = (controls.getId_cellVar("y-gradient x-unit-normal-"+item));
			// gradz_unitNormalId_x[tmp_iter] = (controls.getId_cellVar("z-gradient x-unit-normal-"+item));
			
			// gradx_unitNormalId_y[tmp_iter] = (controls.getId_cellVar("x-gradient y-unit-normal-"+item));
			grady_unitNormalId_y[tmp_iter] = (controls.getId_cellVar("y-gradient y-unit-normal-"+item));
			// gradz_unitNormalId_y[tmp_iter] = (controls.getId_cellVar("z-gradient y-unit-normal-"+item));
			
			// gradx_unitNormalId_z[tmp_iter] = (controls.getId_cellVar("x-gradient z-unit-normal-"+item));
			// grady_unitNormalId_z[tmp_iter] = (controls.getId_cellVar("y-gradient z-unit-normal-"+item));
			gradz_unitNormalId_z[tmp_iter] = (controls.getId_cellVar("z-gradient z-unit-normal-"+item));
			
			
		}
	}
	
	
	// 체적분율 to 레벨셋
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		
		for(int j=0; j<inp_size; ++j){
			int id_ls = id_levelSet_inp[j];
			int id_phi = id_inp[j];
			cellVar_i[id_ls] = cellVar_i[id_phi];
		}
	}

	
	// smoothing Ai
	int iterVFSmoothingMax = 6;
	for(int iter=0; iter<iterVFSmoothingMax; ++iter){
		
		for(int j=0; j<inp_size; ++j){
			int id_ls = id_levelSet_inp[j];

			vector<double> AiUp(mesh.cells.size(),0.0);
			vector<double> AiDown(mesh.cells.size(),0.0);
			auto AiUp_ptr = AiUp.data();
			auto AiDown_ptr = AiDown.data();
			
			
			for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
				auto& face = mesh.faces[i];
				auto faceVar_i = faceVar[i].data();
				
				double wCL = 0.5; double wCR = 1.0-wCL;
				
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					auto cellVar_iL = cellVar[face.iL].data();
					auto cellVar_iR = cellVar[face.iR].data();
					double area = faceVar_i[id_area];
				
					double AiF = wCL*cellVar_iL[id_ls] + wCR*cellVar_iR[id_ls];
						
					AiUp_ptr[face.iL] += AiF*area;
					AiUp_ptr[face.iR] += AiF*area;
					
					AiDown_ptr[face.iL] += area;
					AiDown_ptr[face.iR] += area;
				}
			}

			// boundary
			for(auto& boundary : mesh.boundaries){
				if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
					int str = boundary.startFace;
					int end = str + boundary.nFaces;
					for(int i=str; i<end; ++i){
						auto& face = mesh.faces[i];
						auto faceVar_i = faceVar[i].data();
						auto cellVar_iL = cellVar[face.iL].data();
						double AiF = cellVar_iL[id_ls];
						double area = faceVar_i[id_area];
						
						AiUp_ptr[face.iL] += AiF*area;
						
						AiDown_ptr[face.iL] += area;
					}
				}
			}
			
			if(size>1){
				vector<double> send_value;
				send_value.reserve(proc_size);
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto cellVar_i = cellVar[faces[i].iL].data();
							send_value.push_back(cellVar_i[id_ls]);
						}
					}
				}
				vector<double> recv_value(proc_size);
				MPI_Alltoallv( send_value.data(), mesh.countsProcFaces.data(), 
								mesh.displsProcFaces.data(), MPI_DOUBLE, 
								recv_value.data(), mesh.countsProcFaces.data(), 
								mesh.displsProcFaces.data(), MPI_DOUBLE, 
							   MPI_COMM_WORLD);
				auto recv_value_ptr = recv_value.data();
				int ip=0;
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto& face = mesh.faces[i];
							auto faceVar_i = faceVar[i].data();
							// auto cellVar_iL = cellVar[face.iL].data();
							double area = faceVar_i[id_area];
							
							AiUp_ptr[face.iL] += recv_value_ptr[ip]*area;
							
							AiDown_ptr[face.iL] += area;
							
							++ip;
						}
					}
				}
			}
		
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				auto cellVar_i = cellVar[i].data();
				cellVar_i[id_ls] = AiUp_ptr[i]/AiDown_ptr[i];
			}
		}
		
	}
	
	
	
	
	
	// 레벨셋 구배 계산
	solver.calcGradient.leastSquare_zeroGradient(
		mesh, controls, var, levelSetString);
	
	
	
	// 수직벡터 계산
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		
		for(int j=0; j<inp_size; ++j){
			int id_x = unitNormalId_nx[j];
			int id_y = unitNormalId_ny[j];
			int id_z = unitNormalId_nz[j];
			int id_gradx = id_gradx_inp[j];
			int id_grady = id_grady_inp[j];
			int id_gradz = id_gradz_inp[j];
			
			double magGrad = 0.0;
			magGrad += cellVar_i[id_gradx]*cellVar_i[id_gradx];
			magGrad += cellVar_i[id_grady]*cellVar_i[id_grady];
			magGrad += cellVar_i[id_gradz]*cellVar_i[id_gradz];
			magGrad = sqrt(magGrad);
			
			cellVar_i[id_x] = 0.0;
			cellVar_i[id_y] = 0.0;
			cellVar_i[id_z] = 0.0;
				
			// if( magGrad > std::numeric_limits<double>::min() ){
			if( magGrad > 1.e-15 ){
				cellVar_i[id_x] = cellVar_i[id_gradx]/magGrad;
				cellVar_i[id_y] = cellVar_i[id_grady]/magGrad;
				cellVar_i[id_z] = cellVar_i[id_gradz]/magGrad;
			}
		}
		
	}
	
	// if(size>1){
		// auto id_oup_ptr = id_oup.data();
		
		// vector<double> send_value;
		// send_value.reserve(proc_size*3*inp_size);
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				// for(int i=str; i<end; ++i){
					// auto cellVar_i = cellVar[faces[i].iL].data();
					// auto faceVar_i = faceVar[i].data();
					// for(int j=0; j<inp_size; ++j){
						// int id_x = unitNormalId_x[j];
						// int id_y = unitNormalId_y[j];
						// int id_z = unitNormalId_z[j];
						// send_value.push_back(cellVar_i[id_x]);
						// send_value.push_back(cellVar_i[id_y]);
						// send_value.push_back(cellVar_i[id_z]);
					// }
				// }
			// }
		// }
		
		// vector<int> tmp_counts(size,0);
		// vector<int> tmp_displs(size+1,0);
		// for(int ip=0; ip<size; ++ip){
			// tmp_counts[ip] = mesh.countsProcFaces[ip]*3*inp_size;
		// }
		// for(int ip=0; ip<size; ++ip){
			// tmp_displs[ip+1] = tmp_displs[ip] + tmp_counts[ip];
		// }
		
		// vector<double> recv_value(proc_size*3*inp_size);
		// MPI_Alltoallv( send_value.data(), tmp_counts.data(), 
						// tmp_displs.data(), MPI_DOUBLE, 
						// recv_value.data(), tmp_counts.data(), 
						// tmp_displs.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// auto recv_value_ptr = recv_value.data();
		// int iter=0;
		// for(auto& cells : var.procRightCells){
			// auto cellVar_i = cells.data();
			// for(int j=0; j<inp_size; ++j){
				// int id_x = unitNormalId_x[j];
				// int id_y = unitNormalId_y[j];
				// int id_z = unitNormalId_z[j];
				// cellVar_i[id_x] = recv_value_ptr[iter++];
				// cellVar_i[id_y] = recv_value_ptr[iter++];
				// cellVar_i[id_z] = recv_value_ptr[iter++];
			// }
		// }
	// }
	
	
	
	
	// 수직벡터 구배 계산
	solver.calcGradient.leastSquare_zeroGradient(
		mesh, controls, var, unitNormalString);
	
	
	// 곡률 구하기
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0; j<inp_size; ++j){
			int id = id_curv_oup[j];
			int id_gradx = gradx_unitNormalId_x[j];
			int id_grady = grady_unitNormalId_y[j];
			int id_gradz = gradz_unitNormalId_z[j];
			cellVar_i[id] = (cellVar_i[id_gradx]+cellVar_i[id_grady]+cellVar_i[id_gradz]);
		}
	}
	
	
	// 곡률 전달
	if(size>1){
		// auto id_oup_ptr = id_oup.data();
		
		vector<double> send_value;
		send_value.reserve(proc_size*inp_size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					auto cellVar_i = cellVar[faces[i].iL].data();
					auto faceVar_i = faceVar[i].data();
					for(int j=0; j<inp_size; ++j){
						int id = id_curv_oup[j];
						send_value.push_back(cellVar_i[id]);
					}
				}
			}
		}
		
		vector<int> tmp_counts(size,0);
		vector<int> tmp_displs(size+1,0);
		for(int ip=0; ip<size; ++ip){
			tmp_counts[ip] = mesh.countsProcFaces[ip]*inp_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_displs[ip+1] = tmp_displs[ip] + tmp_counts[ip];
		}
		
		vector<double> recv_value(proc_size*inp_size);
		MPI_Alltoallv( send_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
						recv_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		auto recv_value_ptr = recv_value.data();
		int iter=0;
		for(auto& cells : var.procRightCells){
			auto cellVar_i = cells.data();
			for(int j=0; j<inp_size; ++j){
				int id = id_curv_oup[j];
				cellVar_i[id] = recv_value_ptr[iter++];
			}
		}
	}
	
	
	
	
	
}