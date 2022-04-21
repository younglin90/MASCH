
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
	
	int id_vol = controls.getId_cellVar("volume");
	int id_area = controls.getId_faceVar("area");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_dLR = controls.getId_faceVar("distance of between left and right cell");
	int id_maxA = controls.getId_cellVar("maximum area");
	int id_minA = controls.getId_cellVar("minimum area");

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
	
	
	
	
	vector<bool> interfaceRegion(mesh.cells.size(),false);
    // bool* interfaceRegion_ptr = interfaceRegion.data();
	{
		for(int j=0; j<inp_size; ++j){
			int ip = 0;
			int id_phi = id_inp[j];
			for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
				int iL = faces[i].iL;
				int iR = faces[i].iR;
				auto& face = faces[i];
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					double alphaL = cellVar[iL][id_phi];
					double alphaR = cellVar[iR][id_phi];
					if(
					(alphaL>1.e-8 && alphaL<1.0-1.e-8) ||
					(alphaR>1.e-8 && alphaR<1.0-1.e-8) ||
					(alphaL>=1.0-1.e-8 && alphaR<=1.e-8) ||
					(alphaL<=1.e-8 && alphaR>=1.0-1.e-8)
					){
						interfaceRegion[iL] = true;
						interfaceRegion[iR] = true;
					}
				}
				else if(face.getType() == MASCH_Face_Types::BOUNDARY){
					double alphaL = cellVar[iL][id_phi];
					if(
					(alphaL>1.e-8 && alphaL<1.0-1.e-8)
					){
						interfaceRegion[iL] = true;
					}
				}
				else if(face.getType() == MASCH_Face_Types::PROCESSOR){
					auto cellVar_i = var.procRightCells[ip].data();
					double alphaL = cellVar_i[id_phi];
					if(
					(alphaL>1.e-8 && alphaL<1.0-1.e-8)
					){
						interfaceRegion[iL] = true;
					}
					++ip;
				}
			}
		}
	}
	
	double hmin = 1.e8;
	vector<double> hmin_loc(mesh.cells.size());
    auto hmin_loc_ptr = hmin_loc.data();
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		if(interfaceRegion[i]==false) continue;
		auto cellVar_i = cellVar[i].data();
		hmin_loc_ptr[i] = pow(cellVar_i[id_vol],0.3);
		hmin = min(hmin, hmin_loc_ptr[i]);
	}
	double hmin_glob;
	MPI_Allreduce(&hmin, &hmin_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	hmin = hmin_glob;
	
	// cout << hmin << endl;
	
	
	// 체적분율 to 레벨셋
	// vector<vector<bool>> interfacePresentCell(inp_size,vector<bool>(mesh.cells.size(),false));
	double Ccap = 0.99;
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		
		for(int j=0; j<inp_size; ++j){
			int id_ls = id_levelSet_inp[j];
			int id_phi = id_inp[j];
			double phi = cellVar_i[id_phi];
			
			// hmin[i] = pow(cellVar_i[id_vol],0.3);
				// cellVar_i[id_ls] = phi;
			// if(interfaceRegion[i]==true){
				// cellVar_i[id_ls] = hmin * phi;
				// cellVar_i[id_ls] = phi;
				cellVar_i[id_ls] = 1.0/(1.0-Ccap)*
					(min(max(phi,0.5*Ccap),1.0-0.5*Ccap)-0.5*Ccap);
				// cellVar_i[id_ls] = phi;
			// }
			// else{
				// cellVar_i[id_ls] = 0.0;
			// }
				
			// if(phi>1.e-8 && phi<1.0-1.e-8){
				// interfacePresentCell[j][i] = true;
				// // cellVar_i[id_ls] = phi;
			// }
		}
	}
	
	
	// smoothing 레벨셋
	int iterVFSmoothingMax = 8;
	for(int iter=0; iter<iterVFSmoothingMax; ++iter){
		for(int j=0; j<inp_size; ++j){
			int id_ls = id_levelSet_inp[j];
			vector<double> flux_diff(mesh.cells.size(),0.0);
            auto flux_diff_ptr = flux_diff.data();
			for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
				auto& face = mesh.faces[i];
				auto faceVar_i = faceVar[i].data();
				double area = faceVar_i[id_area];
				double dLR = faceVar_i[id_dLR];
				
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					auto cellVar_iL = cellVar[face.iL].data();
					auto cellVar_iR = cellVar[face.iR].data();
				
					double alphaL = cellVar_iL[id_ls];
					double alphaR = cellVar_iR[id_ls];
					flux_diff_ptr[face.iL] += (alphaR-alphaL)/dLR*area;
					flux_diff_ptr[face.iR] -= (alphaR-alphaL)/dLR*area;
				}
				else if(face.getType() == MASCH_Face_Types::BOUNDARY){
					auto cellVar_iL = cellVar[face.iL].data();
					double alphaL = cellVar_iL[id_ls];
					// flux_diff_ptr[face.iL] += (0.0-alphaL)/dLR*area;
				
				}
			}

			if(size>1){
				vector<double> send_value3;
				send_value3.reserve(proc_size);
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto cellVar_i = cellVar[faces[i].iL].data();
							send_value3.push_back(cellVar_i[id_ls]);
						}
					}
				}
				vector<double> recv_value3(proc_size);
				MPI_Alltoallv( send_value3.data(), mesh.countsSendProcFaces.data(), 
								mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								recv_value3.data(), mesh.countsRecvProcFaces.data(), 
								mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   MPI_COMM_WORLD);
				auto recv_value3_ptr = recv_value3.data();
				int ip=0;
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto& face = mesh.faces[i];
							auto faceVar_i = faceVar[i].data();
							auto cellVar_iL = cellVar[face.iL].data();
							double area = faceVar_i[id_area];
							double dLR = faceVar_i[id_dLR];
							double alphaL = cellVar_iL[id_ls];
							double alphaR = recv_value3_ptr[ip];
							flux_diff_ptr[face.iL] += (alphaR-alphaL)/dLR*area;
							
							++ip;
						}
					}
				}
			}
		
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				auto cellVar_i = cellVar[i].data();
				double vol = cellVar_i[id_vol];
				double maxA = cellVar_i[id_maxA];
				double minA = cellVar_i[id_minA];
				
				// cellVar_i[id_ls] += 0.0001*hmin*flux_diff_ptr[i];
				// cellVar_i[id_ls] += 1.e-9*flux_diff_ptr[i]/vol;
				// if(interfaceRegion[i]==false){
					// cellVar_i[id_ls] += hmin*flux_diff_ptr[i]/vol;
				// }
				double Dcoeff = hmin;
				double CFL = 0.3;
				double dt = minA/Dcoeff/6.0;
					cellVar_i[id_ls] += CFL * dt * (Dcoeff*flux_diff_ptr[i]/vol);
					
				if(abs(cellVar_i[id_ls])>=1.e3){
					cout << "#WARNING : level-set smoothing values >= 1.e3 in curvature" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				// cellVar_i[id_ls] = min(1.e5,max(-1.e5,cellVar_i[id_ls]));
			}
		}
	
	
	
	}

	
	// // smoothing Ai
	// int iterVFSmoothingMax = 1;
	// for(int iter=0; iter<iterVFSmoothingMax; ++iter){
		
		// for(int j=0; j<inp_size; ++j){
			// int id_ls = id_levelSet_inp[j];
			// int id_phi = id_inp[j];

			// vector<double> AiUp(mesh.cells.size(),0.0);
			// vector<double> AiDown(mesh.cells.size(),0.0);
			// auto AiUp_ptr = AiUp.data();
			// auto AiDown_ptr = AiDown.data();
			
			
			// for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
				// auto& face = mesh.faces[i];
				// auto faceVar_i = faceVar[i].data();
				
				// double wCL = 0.5; double wCR = 1.0-wCL;
				// double area = faceVar_i[id_area];
				
				// if(face.getType() == MASCH_Face_Types::INTERNAL){
					// auto cellVar_iL = cellVar[face.iL].data();
					// auto cellVar_iR = cellVar[face.iR].data();
				
					// double AiF = wCL*cellVar_iL[id_ls] + wCR*cellVar_iR[id_ls];
						
					// AiUp_ptr[face.iL] += AiF*area;
					// AiUp_ptr[face.iR] += AiF*area;
					
					// AiDown_ptr[face.iL] += area;
					// AiDown_ptr[face.iR] += area;
				// }
				// else if(face.getType() == MASCH_Face_Types::BOUNDARY){
					// auto cellVar_iL = cellVar[face.iL].data();
				
					// double AiF = cellVar_iL[id_ls];
						
					// AiUp_ptr[face.iL] += AiF*area;
					
					// AiDown_ptr[face.iL] += area;
				// }
			// }
			
			// if(size>1){
				// vector<double> send_value;
				// send_value.reserve(proc_size);
				// for(auto& boundary : mesh.boundaries){
					// if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						// int str = boundary.startFace;
						// int end = str + boundary.nFaces;
						// for(int i=str; i<end; ++i){
							// auto cellVar_i = cellVar[faces[i].iL].data();
							// send_value.push_back(cellVar_i[id_ls]);
						// }
					// }
				// }
				// vector<double> recv_value(proc_size);
				// MPI_Alltoallv( send_value.data(), mesh.countsProcFaces.data(), 
								// mesh.displsProcFaces.data(), MPI_DOUBLE, 
								// recv_value.data(), mesh.countsProcFaces.data(), 
								// mesh.displsProcFaces.data(), MPI_DOUBLE, 
							   // MPI_COMM_WORLD);
				// auto recv_value_ptr = recv_value.data();
				// int ip=0;
				// for(auto& boundary : mesh.boundaries){
					// if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						// int str = boundary.startFace;
						// int end = str + boundary.nFaces;
						// for(int i=str; i<end; ++i){
							// auto& face = mesh.faces[i];
							// auto faceVar_i = faceVar[i].data();
							// // auto cellVar_iL = cellVar[face.iL].data();
							// double area = faceVar_i[id_area];
							
							// AiUp_ptr[face.iL] += recv_value_ptr[ip]*area;
							
							// AiDown_ptr[face.iL] += area;
							
							// ++ip;
						// }
					// }
				// }
			// }
		
			// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				// auto cellVar_i = cellVar[i].data();
				// double phi = cellVar_i[id_phi];
				// double ls_phi = AiUp_ptr[i]/AiDown_ptr[i];
				// if(ls_phi>1.e-8 && ls_phi<1.0-1.e-8){
					// interfacePresentCell[j][i] = true;
					// cellVar_i[id_ls] = ls_phi;
				// }
			// }
		// }
		
	// }
	
	
	// // 레벨셋 구배 계산
	// solver.calcGradient.leastSquare_zeroGradient(
		// mesh, controls, var, levelSetString);
	
	
	// // Hamilton-Jacobi 방정식 풀기 , 레벨셋 계산
	// int iterHJ_Max = 6;
	// double CFL_LS = 0.04;
	// for(int j=0; j<inp_size; ++j){
		// int id_ls = id_levelSet_inp[j];
		// int id_phi = id_inp[j];
		// int id_gradx = id_gradx_inp[j];
		// int id_grady = id_grady_inp[j];
		// int id_gradz = id_gradz_inp[j];
		// vector<double> velocity0(mesh.cells.size(),0.0);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto cellVar_i = cellVar[i].data();
			// double phi = cellVar_i[id_phi];
			// double vol = cellVar_i[id_vol];
			// double tmp_value = cellVar_i[id_ls] - 0.5;
			// if(interfacePresentCell[j][i]){
				// double gradx_value = cellVar_i[id_gradx];
				// double grady_value = cellVar_i[id_grady];
				// double gradz_value = cellVar_i[id_gradz];
				// double magGrad2 = (
					// gradx_value*gradx_value+
					// grady_value*grady_value+
					// gradz_value*gradz_value);
				// if(magGrad2>1.e-2) cellVar_i[id_ls] = tmp_value / magGrad2;
			// }
			// velocity0[i] = (tmp_value > 0.0 ? 1.0 : -1.0) / sqrt(1.0 + pow(vol,0.666));
		// }
		// // 로컬 시간스템 구하기
		// int id_maxA = controls.getId_cellVar("maximum area");
		// vector<double> tau_ls_local(mesh.cells.size(),1.e8);
		// double tau_ls = 1.e15;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto cellVar_i = cellVar[i].data();
			// auto& cell = mesh.cells[i];
			// double maxA = cellVar_i[id_maxA];
			// double vol = cellVar_i[id_vol];
			// // tau_ls_local[i] = min(tau_ls_local[i],min(vol / maxA / 1.0, pow(cell.volume,0.3)));
			// tau_ls_local[i] = min(tau_ls_local[i],pow(vol,0.3));
			// tau_ls = min(tau_ls,tau_ls_local[i]);
		// }
		// double tau_ls_glob;
		// MPI_Allreduce(&tau_ls, &tau_ls_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		// tau_ls = tau_ls_glob;
		// for(int iter=0; iter<iterHJ_Max; ++iter){
			// // 다음 값 계산
			// double residual_ls = 0.0;
			// for(int i=0; i<mesh.cells.size(); ++i){
				// auto cellVar_i = cellVar[i].data();

				// if(!interfacePresentCell[j][i])
				// {
					// double gradx_value = cellVar_i[id_gradx];
					// double grady_value = cellVar_i[id_grady];
					// double gradz_value = cellVar_i[id_gradz];
					// double magGrad = sqrt(
						// gradx_value*gradx_value+
						// grady_value*grady_value+
						// gradz_value*gradz_value);
						
					// double resi_tmp = CFL_LS * tau_ls_local[i]*velocity0[i]*(1.0-magGrad);
					
					// // double interThick = abs(tanh(weightIterThick*smoothAi[i]));
					// // double resi_tmp = interThick * CFL_LS * tau_ls_local[i]*velocity0[i]*(1.0-magGrad);
					
					// cellVar_i[id_ls] += resi_tmp;
					// residual_ls += resi_tmp*resi_tmp;
				// }
			// }
			// // 레벨셋 구배 계산
			// solver.calcGradient.leastSquare_zeroGradient(
				// mesh, controls, var, levelSetString);
		// }
	// }
		
	
	// 레벨셋 구배 계산
	solver.calcGradient.leastSquare_zeroGradient(
		mesh, controls, var, levelSetString);
	
	
	
	// 수직벡터 계산
	vector<vector<double>> magGradAlpha(inp_size,vector<double>(mesh.cells.size()));
    auto magGradAlpha_ptr = magGradAlpha.data();
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		
		for(int j=0; j<inp_size; ++j){
            
            auto magGradAlpha_ptr_j = magGradAlpha_ptr[j].data();
            
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
			
			magGradAlpha_ptr_j[i] = magGrad;
		
			//magGrad += 1.e-8*pow(cellVar_i[id_vol],0.3333);
			cellVar_i[id_x] = 0.0;
			cellVar_i[id_y] = 0.0;
			cellVar_i[id_z] = 0.0;
			if(magGrad!=0.0){
				cellVar_i[id_x] = cellVar_i[id_gradx]/magGrad;
				cellVar_i[id_y] = cellVar_i[id_grady]/magGrad;
				cellVar_i[id_z] = cellVar_i[id_gradz]/magGrad;
			}
				
			// if( magGrad > std::numeric_limits<double>::min() ){
			// if( magGrad > 1.e-200 ){
			// }
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
	
	
	
	
	// // 수직벡터 구배 계산
	// solver.calcGradient.leastSquare_zeroGradient(
		// mesh, controls, var, unitNormalString);
		
	
	
	// 곡률구하기
	//  the quotient rule (Saha et al. and Martinez et al.)
	{
		for(int j=0; j<inp_size; ++j){
			int id = id_curv_oup[j];
			int id_ls = id_levelSet_inp[j];
			int id_x = unitNormalId_nx[j];
			int id_y = unitNormalId_ny[j];
			int id_z = unitNormalId_nz[j];
			// int id_gradx_nx = gradx_unitNormalId_x[j];
			// int id_grady_ny = grady_unitNormalId_y[j];
			// int id_gradz_nz = gradz_unitNormalId_z[j];
			int id_gradx = id_gradx_inp[j];
			int id_grady = id_grady_inp[j];
			int id_gradz = id_gradz_inp[j];

			vector<double> flux(mesh.cells.size(),0.0);
            auto flux_ptr = flux.data();
			// vector<double> flux_diff(mesh.cells.size(),0.0);
			// vector<double> gradx_magGradAlp(mesh.cells.size(),0.0);
			// vector<double> grady_magGradAlp(mesh.cells.size(),0.0);
			// vector<double> gradz_magGradAlp(mesh.cells.size(),0.0);
			for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
				auto& face = mesh.faces[i];
				auto faceVar_i = faceVar[i].data();
				
				double wCL = 0.5; double wCR = 1.0-wCL;
				double area = faceVar_i[id_area];
				double dLR = faceVar_i[id_dLR];
				double nvec[3];
				nvec[0] = faceVar_i[id_nx];
				nvec[1] = faceVar_i[id_ny];
				nvec[2] = faceVar_i[id_nz];
				
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					auto cellVar_iL = cellVar[face.iL].data();
					auto cellVar_iR = cellVar[face.iR].data();
				
					double n_n = (wCL*cellVar_iL[id_x] + wCR*cellVar_iR[id_x])*nvec[0];
					n_n += (wCL*cellVar_iL[id_y] + wCR*cellVar_iR[id_y])*nvec[1];
					n_n += (wCL*cellVar_iL[id_z] + wCR*cellVar_iR[id_z])*nvec[2];
						
					flux_ptr[face.iL] += n_n*area;
					flux_ptr[face.iR] -= n_n*area;
					
					// double magGradAlphaF = 0.5*(magGradAlpha[j][face.iL]+magGradAlpha[j][face.iR]);
					// gradx_magGradAlp[face.iL] += magGradAlphaF*nvec[0]*area;
					// grady_magGradAlp[face.iL] += magGradAlphaF*nvec[1]*area;
					// gradz_magGradAlp[face.iL] += magGradAlphaF*nvec[2]*area;
					
					// gradx_magGradAlp[face.iR] -= magGradAlphaF*nvec[0]*area;
					// grady_magGradAlp[face.iR] -= magGradAlphaF*nvec[1]*area;
					// gradz_magGradAlp[face.iR] -= magGradAlphaF*nvec[2]*area;
					
					// double alphaL = cellVar_iL[id_ls];
					// double alphaR = cellVar_iR[id_ls];
					// flux_diff[face.iL] += (alphaR-alphaL)/dLR*area;
					// flux_diff[face.iR] -= (alphaR-alphaL)/dLR*area;
				}
				else if(face.getType() == MASCH_Face_Types::INTERNAL){
					auto cellVar_iL = cellVar[face.iL].data();
				
					double n_n = cellVar_iL[id_x]*nvec[0];
					n_n += cellVar_iL[id_y]*nvec[1];
					n_n += cellVar_iL[id_z]*nvec[2];
                    
					flux_ptr[face.iL] += n_n*area;
					
					// double magGradAlphaF = magGradAlpha[j][face.iL];
					// gradx_magGradAlp[face.iL] += magGradAlphaF*nvec[0]*area;
					// grady_magGradAlp[face.iL] += magGradAlphaF*nvec[1]*area;
					// gradz_magGradAlp[face.iL] += magGradAlphaF*nvec[2]*area;
						
					// double alphaL = cellVar_iL[id_ls];
					// flux_diff[face.iL] += (alphaR-alphaL)/dLR*area;
				}
			}

			if(size>1){
				vector<double> send_value0,send_value1,send_value2,send_value3,send_value4;
				send_value0.reserve(proc_size);
				send_value1.reserve(proc_size);
				send_value2.reserve(proc_size);
				// send_value3.reserve(proc_size);
				// send_value4.reserve(proc_size);
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto cellVar_i = cellVar[faces[i].iL].data();
							send_value0.push_back(cellVar_i[id_x]);
							send_value1.push_back(cellVar_i[id_y]);
							send_value2.push_back(cellVar_i[id_z]);
							// send_value3.push_back(cellVar_i[id_ls]);
							// send_value4.push_back(magGradAlpha[j][faces[i].iL]);
						}
					}
				}
				vector<double> recv_value0(proc_size);
				vector<double> recv_value1(proc_size);
				vector<double> recv_value2(proc_size);
				// vector<double> recv_value3(proc_size);
				// vector<double> recv_value4(proc_size);
				MPI_Alltoallv( send_value0.data(), mesh.countsSendProcFaces.data(), 
								mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								recv_value0.data(), mesh.countsRecvProcFaces.data(), 
								mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   MPI_COMM_WORLD);
				MPI_Alltoallv( send_value1.data(), mesh.countsSendProcFaces.data(), 
								mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								recv_value1.data(), mesh.countsRecvProcFaces.data(), 
								mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   MPI_COMM_WORLD);
				MPI_Alltoallv( send_value2.data(), mesh.countsSendProcFaces.data(), 
								mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								recv_value2.data(), mesh.countsRecvProcFaces.data(), 
								mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   MPI_COMM_WORLD);
				// MPI_Alltoallv( send_value3.data(), mesh.countsSendProcFaces.data(), 
								// mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								// recv_value3.data(), mesh.countsRecvProcFaces.data(), 
								// mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   // MPI_COMM_WORLD);
				// MPI_Alltoallv( send_value4.data(), mesh.countsSendProcFaces.data(), 
								// mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								// recv_value4.data(), mesh.countsRecvProcFaces.data(), 
								// mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
							   // MPI_COMM_WORLD);
				auto recv_value0_ptr = recv_value0.data();
				auto recv_value1_ptr = recv_value1.data();
				auto recv_value2_ptr = recv_value2.data();
				// auto recv_value3_ptr = recv_value3.data();
				// auto recv_value4_ptr = recv_value4.data();
				int ip=0;
				for(auto& boundary : mesh.boundaries){
					if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						for(int i=str; i<end; ++i){
							auto& face = mesh.faces[i];
							auto faceVar_i = faceVar[i].data();
							auto cellVar_iL = cellVar[face.iL].data();
							double area = faceVar_i[id_area];
							double dLR = faceVar_i[id_dLR];
							double nvec[3];
							nvec[0] = faceVar_i[id_nx];
							nvec[1] = faceVar_i[id_ny];
							nvec[2] = faceVar_i[id_nz];
							
							double wCL = 0.5; double wCR = 1.0-wCL;
						
							double n_n = (wCL*cellVar_iL[id_x] + wCR*recv_value0_ptr[ip])*nvec[0];
							n_n += (wCL*cellVar_iL[id_y] + wCR*recv_value1_ptr[ip])*nvec[1];
							n_n += (wCL*cellVar_iL[id_z] + wCR*recv_value2_ptr[ip])*nvec[2];
							
							flux_ptr[face.iL] += n_n*area;
									
							// double magGradAlphaF = 0.5*(magGradAlpha[j][face.iL]+recv_value4_ptr[ip]);
							// gradx_magGradAlp[face.iL] += magGradAlphaF*nvec[0]*area;
							// grady_magGradAlp[face.iL] += magGradAlphaF*nvec[1]*area;
							// gradz_magGradAlp[face.iL] += magGradAlphaF*nvec[2]*area;
							
							// double alphaL = cellVar_iL[id_ls];
							// double alphaR = recv_value3_ptr[ip];
							// flux_diff[face.iL] += (alphaR-alphaL)/dLR*area;
							
							++ip;
						}
					}
				}
			}
		
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				auto cellVar_i = cellVar[i].data();
				double vol = cellVar_i[id_vol];
				// // cellVar_i[id] = -flux[i]/vol;
				// // double magGrad = 0.0;
				// // magGrad += cellVar_i[id_gradx]*cellVar_i[id_gradx];
				// // magGrad += cellVar_i[id_grady]*cellVar_i[id_grady];
				// // magGrad += cellVar_i[id_gradz]*cellVar_i[id_gradz];
				// // magGrad = sqrt(magGrad);
				
				// flux_diff[i] /= vol;
				// gradx_magGradAlp[i] /= vol;
				// grady_magGradAlp[i] /= vol;
				// gradz_magGradAlp[i] /= vol;
				
				// double first_term = 0.0;
				// first_term += cellVar_i[id_x]*gradx_magGradAlp[i];
				// first_term += cellVar_i[id_y]*grady_magGradAlp[i];
				// first_term += cellVar_i[id_z]*gradz_magGradAlp[i];
				
				// double tmp_kappa = (first_term - flux_diff[i]);
				
				// double magGrad = magGradAlpha[j][i];
				
				// cellVar_i[id] = 0.0;
				// if(magGrad!=0){
					// cellVar_i[id] = tmp_kappa/magGrad;
				// }
				
				
				cellVar_i[id] = -flux_ptr[i]/vol;
			}
		}
	}		
	
		
	
	
	// // 곡률 구하기
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto cellVar_i = cellVar[i].data();
		// for(int j=0; j<inp_size; ++j){
			// int id = id_curv_oup[j];
			// int id_gradx = gradx_unitNormalId_x[j];
			// int id_grady = grady_unitNormalId_y[j];
			// int id_gradz = gradz_unitNormalId_z[j];
			// cellVar_i[id] = -(cellVar_i[id_gradx]+cellVar_i[id_grady]+cellVar_i[id_gradz]);
		// }
	// }
	
	
	
	
	//================================================
	// smoothing process of Curvature (The Sharp Surface Force (SSF) Model)
	{
		int iterFirstSmoothingMax = 2;
		// first smoothing
		for(int i_id=0; i_id<inp_size; ++i_id){
			int id_phi = id_inp[i_id];
			// int id_ls = id_levelSet_inp[i_id];
			// int id_phi = id_levelSet_inp[i_id];
			int id_kappa = id_curv_oup[i_id];
			int id_x = unitNormalId_nx[i_id];
			int id_y = unitNormalId_ny[i_id];
			int id_z = unitNormalId_nz[i_id];
			vector<double> kappa(mesh.cells.size(),0.0);
            auto kappa_ptr = kappa.data();
			for(int i=0; i<mesh.cells.size(); ++i){
				auto cellVar_i = cellVar[i].data();
				kappa_ptr[i] = cellVar_i[id_kappa];
			}
			
			
			for(int iter=0; iter<iterFirstSmoothingMax; ++iter)
			{
				vector<double> smoothUp(mesh.cells.size(),0.0);
				vector<double> smoothDown(mesh.cells.size(),0.0);
                auto smoothUp_ptr = smoothUp.data();
                auto smoothDown_ptr = smoothDown.data();
				
				for(int i=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					int iR = face.iR;
					
					if(face.getType() == MASCH_Face_Types::INTERNAL){
						auto cellVar_L = cellVar[iL].data();
						auto cellVar_R = cellVar[iR].data();
						
						double alpha_org_L = cellVar_L[id_phi];
						double alpha_org_R = cellVar_R[id_phi];
						double alpha_org_F = 0.5*(alpha_org_L+alpha_org_R);
						double weight = sqrt(alpha_org_F*(1.0-alpha_org_F)+1.e-3);
						double kappa_F = 0.5*(kappa_ptr[iL]+kappa_ptr[iR]);
						
						smoothUp_ptr[iL] += weight*kappa_F;
						smoothDown_ptr[iL] += weight;
						
						smoothUp_ptr[iR] += weight*kappa_F;
						smoothDown_ptr[iR] += weight;
					}
					else if(face.getType() == MASCH_Face_Types::BOUNDARY){
						auto cellVar_L = cellVar[iL].data();
						
						double alpha_org_L = cellVar_L[id_phi];
						double alpha_org_F = alpha_org_L;
						double weight = sqrt(alpha_org_F*(1.0-alpha_org_F)+1.e-6);
						
						smoothUp_ptr[iL] += weight*kappa_ptr[iL];
						smoothDown_ptr[iL] += weight;
					}
					
				}
				
				if(size>1){
					// processor faces
					vector<double> sendValues, sendValues2;
                    sendValues.reserve(proc_size);
                    sendValues2.reserve(proc_size);
					for(int i=0; i<mesh.faces.size(); ++i){
						auto& face = mesh.faces[i];
						int iL = face.iL;
						
						if(face.getType() == MASCH_Face_Types::PROCESSOR){
							auto cellVar_L = cellVar[iL].data();
							sendValues.push_back(kappa_ptr[iL]);
							sendValues2.push_back(cellVar_L[id_phi]);
						}
					}
					vector<double> recvValues(sendValues.size()), smoothAi_recv(sendValues2.size());
					MPI_Alltoallv( sendValues.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								   recvValues.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
								   MPI_COMM_WORLD);
					MPI_Alltoallv( sendValues2.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_DOUBLE, 
								   smoothAi_recv.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_DOUBLE, 
								   MPI_COMM_WORLD);
					for(int i=0, ip=0; i<mesh.faces.size(); ++i){
						auto& face = mesh.faces[i];
						int iL = face.iL;
					
						if(face.getType() == MASCH_Face_Types::PROCESSOR){ 
							auto cellVar_L = cellVar[iL].data();
							
							double alpha_org_L = cellVar_L[id_phi];
							double alpha_org_R = smoothAi_recv[ip];
							double alpha_org_F = 0.5*(alpha_org_L+alpha_org_R);
							double weight = sqrt(alpha_org_F*(1.0-alpha_org_F)+1.e-3);
							double kappa_F = 0.5*(kappa_ptr[iL]+recvValues[ip]);
							
							smoothUp_ptr[iL] += weight*kappa_F;
							smoothDown_ptr[iL] += weight;
							
							++ip;
						}
					}
					
				}
				for(int i=0; i<mesh.cells.size(); ++i){
					auto cellVar_i = cellVar[i].data();
					double alpha_i = cellVar_i[id_phi];
					double kappa_i = smoothUp_ptr[i]/smoothDown_ptr[i];
					double weight = 2.0*sqrt(alpha_i*(1.0-alpha_i));
					kappa_ptr[i] = weight*kappa_ptr[i] + (1.0-weight)*kappa_i;
				}
			}
			
			for(int i=0; i<mesh.cells.size(); ++i){
				auto cellVar_i = cellVar[i].data();
				cellVar_i[id_kappa] = kappa_ptr[i];
			}
		}
	}
	//================================================

	
	
	// //================================================
	// // smoothing process of Curvature (denner)
	// {
		// int iterFirstSmoothingMax = 1;
		// int iterSecondSmoothingMax = 1;
		// double tanhWeight = 2.0;
		// // first smoothing
		// for(int i_id=0; i_id<inp_size; ++i_id){
			// int id_kappa = id_curv_oup[i_id];
			// int id_ls = id_levelSet_inp[i_id];
			// int id_x = unitNormalId_nx[i_id];
			// int id_y = unitNormalId_ny[i_id];
			// int id_z = unitNormalId_nz[i_id];
			// vector<double> kappa(mesh.cells.size(),0.0);
			// vector<double> smoothAi(mesh.cells.size(),0.0);
			// for(int i=0; i<mesh.cells.size(); ++i){
				// auto cellVar_i = cellVar[i].data();
				// kappa[i] = cellVar_i[id_kappa];
				// smoothAi[i] = cellVar_i[id_ls];
			// }
			
			
			// for(int iter=0; iter<iterFirstSmoothingMax; ++iter)
			// {
				// vector<double> smoothUp(mesh.cells.size(),0.0);
				// vector<double> smoothDown(mesh.cells.size(),0.0);
				// for(int i=0; i<mesh.cells.size(); ++i){
					// // double weight = pow(1.0-2.0*abs(0.5-orgAi[i]), 8.0);
					// double weight = 1.0-abs(tanh(tanhWeight*smoothAi[i]));
					// smoothUp[i] = kappa[i]*weight;
					// smoothDown[i] = weight;
				// }
				
				// for(int i=0; i<mesh.faces.size(); ++i){
					// auto& face = mesh.faces[i];
					// int iL = face.iL;
					// int iR = face.iR;
					
					// // double wCL = face.wC; double wCR = 1.0 - wCL;
					
					// if(face.getType() == MASCH_Face_Types::INTERNAL){
						
						// // double weightL = pow(1.0-2.0*abs(0.5-orgAi[iL]), 8.0);
						// // double weightR = pow(1.0-2.0*abs(0.5-orgAi[iR]), 8.0);
						// double weightL = 1.0-abs(tanh(tanhWeight*smoothAi[iL]));
						// double weightR = 1.0-abs(tanh(tanhWeight*smoothAi[iR]));
					
						// smoothUp[iL] += weightR*kappa[iR];
						// smoothDown[iL] += weightR;
						
						// smoothUp[iR] += weightL*kappa[iL];
						// smoothDown[iR] += weightL;
						
						// // double avgVar = wCL*smoothAi[iL] + wCR*smoothAi[iR];
						// // double weight = pow(1.0-2.0*avgVar, 8.0) * face.area;
						// // double avgKappa = wCL*kappa[iL] + wCR*kappa[iR];
					
						// // smoothUp[iL] += weight*avgKappa;
						// // smoothDown[iL] += weight;
						
						// // smoothUp[iR] += weight*avgKappa;
						// // smoothDown[iR] += weight;
					
					// }
					
				// }
				
				// if(size>1){
					// // processor faces
					// vector<double> sendValues, sendValues2;
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
						
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){
							// sendValues.push_back(kappa[iL]);
							// sendValues2.push_back(smoothAi[iL]);
							// // orgAi_send.push_back(orgAi[iL]);
						// }
					// }
					// vector<double> recvValues(sendValues.size()), smoothAi_recv(sendValues2.size());
					// MPI_Alltoallv( sendValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // recvValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // MPI_COMM_WORLD);
					// MPI_Alltoallv( sendValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // smoothAi_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // MPI_COMM_WORLD);
					// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
					
						// // double wCL = face.wC; double wCR = 1.0 - wCL;
						
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){ 
							// // double weightR = pow(1.0-2.0*abs(0.5-orgAi_recv[ip]), 8.0);
							// double weightR = 1.0-abs(tanh(tanhWeight*smoothAi_recv[ip]));
							
							// smoothUp[iL] += weightR*recvValues[ip];
							// smoothDown[iL] += weightR;
							
							// // double avgVar = wCL*smoothAi[iL] + wCR*smoothAi_recv[ip];
							// // double weight = pow(1.0-2.0*avgVar, 8.0) * face.area;
							// // double avgKappa = wCL*kappa[iL] + wCR*recvValues[ip];
							
							// // smoothUp[iL] += weight*avgKappa;
							// // smoothDown[iL] += weight;
							
							// ++ip;
						// }
					// }
					
				// }
			
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
						// kappa[i] = smoothUp[i]/smoothDown[i];
					// }
					// // else{
						// // kappa[i] = 0.0;
					// // }
				// }
				
			// }
			
			
			
			
			// // second smoothing
			// for(int iter=0; iter<iterSecondSmoothingMax; ++iter)
			// {
				// vector<double> smoothUp(mesh.cells.size(),0.0);
				// vector<double> smoothDown(mesh.cells.size(),0.0);
				// for(int i=0; i<mesh.cells.size(); ++i){
					// // double weight = pow(1.0-2.0*abs(0.5-orgAi[i]), 8.0);
					// double weight = 1.0-abs(tanh(tanhWeight*smoothAi[i]));
					// smoothUp[i] = kappa[i]*weight;
					// smoothDown[i] = weight;
				// }
				
				// for(int i=0; i<mesh.faces.size(); ++i){
					// auto& face = mesh.faces[i];
					// int iL = face.iL;
					// int iR = face.iR;
					// // double wCL = face.wC; double wCR = 1.0 - wCL;
					
					// if(face.getType() == MASCH_Face_Types::INTERNAL){
						// auto cellVarL = cellVar[iL].data();
						// auto cellVarR = cellVar[iR].data();
						
						// // double weightL = pow(1.0-2.0*abs(0.5-orgAi[iL]), 8.0);
						// // double weightR = pow(1.0-2.0*abs(0.5-orgAi[iR]), 8.0);
						// double weightL = 1.0-abs(tanh(tanhWeight*smoothAi[iL]));
						// double weightR = 1.0-abs(tanh(tanhWeight*smoothAi[iR]));
						
						// double weightL_m = 0.0;
						// weightL_m += cellVarR[id_x]*faceVar[i][id_nx];
						// weightL_m += cellVarR[id_y]*faceVar[i][id_ny];
						// weightL_m += cellVarR[id_z]*faceVar[i][id_nz];
						// weightL_m = abs(weightL_m);
						// weightL_m = pow(weightL_m, 8.0);
						// double weightR_m = 0.0;
						// weightR_m += cellVarL[id_x]*faceVar[i][id_nx];
						// weightR_m += cellVarL[id_y]*faceVar[i][id_ny];
						// weightR_m += cellVarL[id_z]*faceVar[i][id_nz];
						// weightR_m = abs(weightR_m);
						// weightR_m = pow(weightR_m, 8.0);
					
						// smoothUp[iL] += weightR*weightR_m*kappa[iR];
						// smoothUp[iR] += weightL*weightL_m*kappa[iL];
						
						// smoothDown[iL] += weightR*weightR_m;
						// smoothDown[iR] += weightL*weightL_m;
					
					// }
					
				// }
				
				// if(size>1){
					// // processor faces
					// vector<double> sendValues, sendValues2;
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
						
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){
							// sendValues.push_back(kappa[iL]);
							// sendValues2.push_back(smoothAi[iL]);
							// // orgAi_send.push_back(orgAi[iL]);
						// }
					// }
					// vector<double> recvValues(sendValues.size()), smoothAi_recv(sendValues2.size());
					// MPI_Alltoallv( sendValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // recvValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // MPI_COMM_WORLD);
					// MPI_Alltoallv( sendValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // smoothAi_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
								   // MPI_COMM_WORLD);
					// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// int iL = face.iL;
					
						// // double wCL = face.wC; double wCR = 1.0 - wCL;
						
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){ 
							// auto cellVarL = cellVar[iL].data();
							// // double weightR = pow(1.0-2.0*abs(0.5-orgAi_recv[ip]), 8.0);
							// double weightR = 1.0-abs(tanh(tanhWeight*smoothAi_recv[ip]));
							
							// double weightR_m = 0.0;
							// // weightR_m += surfNormalVec_recv[0][ip]*face.unitNomalsPN[0];
							// // weightR_m += surfNormalVec_recv[1][ip]*face.unitNomalsPN[1];
							// // weightR_m += surfNormalVec_recv[2][ip]*face.unitNomalsPN[2];
							// weightR_m += cellVarL[id_x]*faceVar[i][id_nx];
							// weightR_m += cellVarL[id_y]*faceVar[i][id_ny];
							// weightR_m += cellVarL[id_z]*faceVar[i][id_nz];
							// weightR_m = abs(weightR_m);
							// weightR_m = pow(weightR_m, 8.0);
							
							// smoothUp[iL] += weightR*weightR_m*recvValues[ip];
							// smoothDown[iL] += weightR*weightR_m;
							
							// ++ip;
						// }
					// }
					
				// }
			
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
						// kappa[i] = smoothUp[i]/smoothDown[i];
					// }
				// }
				
			// }
			// for(int i=0; i<mesh.cells.size(); ++i){
				// auto cellVar_i = cellVar[i].data();
				// cellVar_i[id_kappa] = kappa[i];
			// }
			
			
		// }
	// }
	// //================================================
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
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
			tmp_counts[ip] = mesh.countsSendProcFaces[ip]*inp_size;
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