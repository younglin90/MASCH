
#include "./mesh.h"
#include "./polyAMR.h"
// #include "geometric.h" 
#include "./mpi.h"
#include "./solvers.h" 
#include "./save.h"  

// void gradientTerms_AMR(MASCH_Mesh& mesh, MASCH_Control& controls, 
	// MASCH_Solver& solver, MASCH_Variables& var);

void MASCH_Poly_AMR_Builder::calcIndicators(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var,
int maxBuffer,
int maxLevel,
int maxCells,
double minVolume_AMR,
vector<vector<double>>& indicatorCriterion,
vector<vector<int>>& indicatorAMR_id, 
vector<bool>& boolCellRefine, 
vector<bool>& boolCellUnrefine,
vector<bool>& boolCellPreserved
){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int id_vol = controls.getId_cellVar("volume");
	
	// 인디케이터 값 초기화
	int nIndi = indicatorAMR_id.size();
	vector<vector<double>> indicatorScalValues(nIndi);
	vector<vector<double>> indicatorGradValues(nIndi);
	for(int j=0, iter=0; j<nIndi; ++j){
		indicatorScalValues[j].clear();
		indicatorScalValues[j].resize(mesh.cells.size(),0.0);
		indicatorGradValues[j].clear();
		indicatorGradValues[j].resize(mesh.cells.size(),0.0);
	}
	for(int i=0; i<mesh.cells.size(); ++i){
		for(int j=0, iter=0; j<nIndi; ++j){
			int id_scal = indicatorAMR_id[j][0];
			int id_gradx = indicatorAMR_id[j][1];
			int id_grady = indicatorAMR_id[j][2];
			int id_gradz = indicatorAMR_id[j][3];
			indicatorScalValues[j][i] = var.cells[i][id_scal];
			double mag_grad = sqrt(
				var.cells[i][id_gradx]*var.cells[i][id_gradx]+
				var.cells[i][id_grady]*var.cells[i][id_grady]+
				var.cells[i][id_gradz]*var.cells[i][id_gradz]);
			indicatorGradValues[j][i] = mag_grad;
		}
	}
	
	
	

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		for(int j=0; j<nIndi; ++j){
			double indicatorRefine_scal = indicatorCriterion[j][0];
			double indicatorRefine_grad = indicatorCriterion[j][1];
			double indiScalVal = indicatorScalValues[j][i];
			double indiGradVal = indicatorGradValues[j][i];
			if( 
			indiScalVal > indicatorRefine_scal &&
			indiGradVal > indicatorRefine_grad
			){
				boolCellPreserved[i] = true;
				boolCellRefine[i] = true;
				boolCellUnrefine[i] = false;
			}
		}
		// if(mesh.cells[i].volume < minVolume_AMR) boolCellRefine[i] = false;
		if(var.cells[i][id_vol] < minVolume_AMR) boolCellRefine[i] = false;
		if(cell.level >= maxLevel) boolCellRefine[i] = false;
		if(cell.level < 0) boolCellRefine[i] = false;
		// if(boolCellPreserved[i] == true) boolCellRefine[i] = false;
		
		
		
		
		// if(cell.level < 0) {
			
			// cout << "GGGGGGGGGGG" << endl;
		// }
		
	} 
	
	
	
	{
		// amr 정보 mpi 교환
		vector<int> cLevel_recv;
		vector<int> cRefine_recv;
		this->mpiLevelRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	
		// // for(int iLevel=maxLevel; iLevel>=0; --iLevel)
		// {
			// vector<bool> tmp_boolCellRefine(mesh.cells.size());
			// for(int i=0; i<mesh.cells.size(); ++i){
				// tmp_boolCellRefine[i] = boolCellRefine[i];
			// }
			// int iLevel=maxLevel;
			// for(int iBuffer=0; iBuffer<maxBuffer-1; ++iBuffer){
				// vector<int> cPreserved_recv;
				// vector<int> send_value2;
				// if(size>1){
					// vector<int> send_value;
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						// if(face.getType() == MASCH_Face_Types::PROCESSOR){
							// send_value.push_back(boolCellPreserved[face.iL]);
							// send_value2.push_back(boolCellRefine[face.iL]);
						// }
					// }
					// cPreserved_recv.resize(send_value.size());
					// MPI_Alltoallv( send_value.data(), mesh.countsProcFaces.data(), 
									// mesh.displsProcFaces.data(), MPI_INT, 
									// cPreserved_recv.data(), mesh.countsProcFaces.data(), 
									// mesh.displsProcFaces.data(), MPI_INT, 
								   // MPI_COMM_WORLD);
					// MPI_Alltoallv( send_value2.data(), mesh.countsProcFaces.data(), 
									// mesh.displsProcFaces.data(), MPI_INT, 
									// cRefine_recv.data(), mesh.countsProcFaces.data(), 
									// mesh.displsProcFaces.data(), MPI_INT, 
								   // MPI_COMM_WORLD);
					
				// }
				
				
				// // for(int i=0; i<mesh.nInternalFaces; ++i){
					// // auto& face = mesh.faces[i];
					// // int iL = face.iL;
					// // int iR = face.iR;
					// // auto& cellL = mesh.cells[iL];
					// // auto& cellR = mesh.cells[iR];
					// // if(boolCellPreserved[iR]==true){
						// // tmp_boolCellRefine[iL]=true;
					// // }
					// // if(boolCellPreserved[iL]==true){
						// // tmp_boolCellRefine[iR]=true;
					// // }
				// // }
				// // int ip=0;
				// // for(auto& boundary : mesh.boundaries){
					// // int str = boundary.startFace;
					// // int end = str + boundary.nFaces;
					// // if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
						// // for(int i=str; i<end; ++i){
							// // auto& face = mesh.faces[i];
							// // int iL = face.iL;
							// // auto& cellL = mesh.cells[iL];
							// // if(cPreserved_recv[ip]==1){
								// // tmp_boolCellRefine[iL]=true;
							// // }
							
							// // ++ip;
						// // }
					// // }
				// // }	
				
				// // for(int i=0; i<mesh.cells.size(); ++i){
					// // boolCellRefine[i]=tmp_boolCellRefine[i];
					// // if(boolCellRefine[i]==true){
						// // boolCellPreserved[i] = true;
						// // boolCellUnrefine[i] = false;
					// // }
				// // }
			// // }
		// }
	
	
		
		{
			vector<bool> tmp_boolCellPreserved(mesh.cells.size(),false);
			vector<bool> tmp_boolCellRefine(mesh.cells.size());
			vector<bool> boolCellRefine0(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i){
				tmp_boolCellPreserved[i] = boolCellPreserved[i];
				tmp_boolCellRefine[i] = boolCellRefine[i];
				boolCellRefine0[i] = boolCellRefine[i];
			}
			for(int iLevel=maxLevel; iLevel>=0; --iLevel)
			{
				int nBuffer = maxBuffer*(1+maxLevel-iLevel);
				for(int iBuffer=0; iBuffer<nBuffer; ++iBuffer){
					vector<int> cPreserved_recv;
					vector<int> send_value2;
					if(size>1){
						vector<int> send_value;
						for(int i=0; i<mesh.faces.size(); ++i){
							auto& face = mesh.faces[i];
							if(face.getType() == MASCH_Face_Types::PROCESSOR){
								send_value.push_back(boolCellPreserved[face.iL]);
								send_value2.push_back(boolCellRefine[face.iL]);
							}
						}
						cPreserved_recv.resize(send_value.size());
						MPI_Alltoallv( send_value.data(), mesh.countsProcFaces.data(), 
										mesh.displsProcFaces.data(), MPI_INT, 
										cPreserved_recv.data(), mesh.countsProcFaces.data(), 
										mesh.displsProcFaces.data(), MPI_INT, 
									   MPI_COMM_WORLD);
						MPI_Alltoallv( send_value2.data(), mesh.countsProcFaces.data(), 
										mesh.displsProcFaces.data(), MPI_INT, 
										cRefine_recv.data(), mesh.countsProcFaces.data(), 
										mesh.displsProcFaces.data(), MPI_INT, 
									   MPI_COMM_WORLD);
						
					}
					for(int i=0; i<mesh.nInternalFaces; ++i){
						auto& face = mesh.faces[i];
						int iL = face.iL;
						int iR = face.iR;
						auto& cellL = mesh.cells[iL];
						auto& cellR = mesh.cells[iR];
						if(
						boolCellPreserved[iR] == true &&
						boolCellPreserved[iL] == false
						){
							if(
							cellL.level == iLevel
							){
								tmp_boolCellPreserved[iL] = true;
								tmp_boolCellRefine[iL] = false;
								boolCellUnrefine[iL] = false;
							}
							else if(
							cellL.level == iLevel-1
							){
								tmp_boolCellPreserved[iL] = true;
								tmp_boolCellRefine[iL] = true;
								boolCellUnrefine[iL] = false;
							}
							else if(
							cellL.level == iLevel+1
							){
								tmp_boolCellPreserved[iL] = true;
								tmp_boolCellRefine[iL] = false;
								boolCellUnrefine[iL] = true;
							}
						}
						// ======================
						if(
						boolCellPreserved[iL] == true &&
						boolCellPreserved[iR] == false
						){
							if(
							cellR.level == iLevel
							){
								tmp_boolCellPreserved[iR] = true;
								tmp_boolCellRefine[iR] = false;
								boolCellUnrefine[iR] = false;
							}
							else if(
							cellR.level == iLevel-1
							){
								tmp_boolCellPreserved[iR] = true;
								tmp_boolCellRefine[iR] = true;
								boolCellUnrefine[iR] = false;
							}
							else if(
							cellR.level == iLevel+1
							){
								tmp_boolCellPreserved[iR] = true;
								tmp_boolCellRefine[iR] = false;
								boolCellUnrefine[iR] = true;
							}
						}
						
					}
					int ip=0;
					for(auto& boundary : mesh.boundaries){
						int str = boundary.startFace;
						int end = str + boundary.nFaces;
						if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
							for(int i=str; i<end; ++i){
								auto& face = mesh.faces[i];
								int iL = face.iL;
								auto& cellL = mesh.cells[iL];
								if(
								cPreserved_recv[ip] == true &&
								boolCellPreserved[iL] == false
								){
									if(
									cellL.level == iLevel
									){
										tmp_boolCellPreserved[iL] = true;
										tmp_boolCellRefine[iL] = false;
										boolCellUnrefine[iL] = false;
									}
									else if(
									cellL.level == iLevel-1
									){
										tmp_boolCellPreserved[iL] = true;
										tmp_boolCellRefine[iL] = true;
										boolCellUnrefine[iL] = false;
									}
									else if(
									cellL.level == iLevel+1
									){
										tmp_boolCellPreserved[iL] = true;
										tmp_boolCellRefine[iL] = false;
										boolCellUnrefine[iL] = true;
									}
								}
								++ip;
							}
						}
					}	
					
					for(int i=0; i<mesh.cells.size(); ++i){
						boolCellPreserved[i]=tmp_boolCellPreserved[i];
						boolCellRefine[i]=tmp_boolCellRefine[i];
					}
				}
			}
			// {
				// int iLevel = maxLevel-1;
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if(mesh.cells[i].level==iLevel) boolCellRefine[i] = boolCellRefine0[i];
				// }
			// }
			// {
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if(boolCellRefine[i]==true) boolCellPreserved[i] = true;
					// // if(boolCellRefine[i]==false) boolCellPreserved[i] = false;
				// }
			// }
		}
		
		// 대각선 방향에서, 레벨 차이 2 이상이면 리파인
		{
			// processor faces
			vector<int> recv_value;
			if(size>1){

				vector<int> send_value;
				send_value.reserve(mesh.send_StencilCellsId.size());
				for(auto& icell : mesh.send_StencilCellsId){
					send_value.push_back(mesh.cells[icell].level);
				}
				recv_value.resize(mesh.recv_displsStencilCells[size]);
				MPI_Alltoallv( send_value.data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_INT, 
							   recv_value.data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_INT, 
							   MPI_COMM_WORLD);
				
			}
			vector<int> tmp_maxLevel(mesh.cells.size());
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				auto& cell = mesh.cells[i];
				int maxInd = -100;
				for(auto& icell : cell.iStencils){
					maxInd = max(mesh.cells[icell].level,maxInd);
				}
				for(auto& icell : cell.recv_iStencils){
					maxInd = max(recv_value[icell],maxInd);
				}
				tmp_maxLevel[i] = maxInd;
			}
			// int inp_size = indicatorValues.size();
			for(int i=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				int my_level = cell.level;
				if(my_level==maxLevel) continue;
				if(my_level<0) continue;
				// if(boolCellPreserved[i] == true) continue;
				if(my_level<tmp_maxLevel[i]-2){
					boolCellRefine[i] = true;
					boolCellPreserved[i]=true;
					boolCellUnrefine[i]=false;
				}
			}
		}
			
	}
	
	
	
	
}


void MASCH_Poly_AMR_Builder::polyAMR_inline(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	MASCH_Solver& solver,
	MASCH_Variables& var,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Load load;
	
	// cout << iter+1 % 10 << endl;
	// for(int ii=0; ii<39; ++ii)
	if( (iter+1) %
		stoi(controls.dynamicMeshMap["AMR.interval"]) == 0)
	{
		int maxLevel = stoi(controls.dynamicMeshMap["AMR.maxLevel"]);
		int maxBuffer = stoi(controls.dynamicMeshMap["AMR.maxBufferLayer"]);
		double minVolume = stod(controls.dynamicMeshMap["AMR.minVolume"]);
		int maxCells = stoi(controls.dynamicMeshMap["AMR.maxCells"]);
		
		vector<string> sScalIndi = load.extractVector(controls.dynamicMeshMap["AMR.indicatorCellNames"]);
		vector<string> sGradIndi = load.extractVector(controls.dynamicMeshMap["AMR.indicatorGradientNames"]);
		vector<string> sScalIndiValues = load.extractVector(controls.dynamicMeshMap["AMR.indicatorCellValues"]);
		vector<string> sGradIndiValues = load.extractVector(controls.dynamicMeshMap["AMR.indicatorGradientValues"]);
		
		int nCriterion = sScalIndi.size();
		vector<vector<double>> indicatorCriterion(nCriterion);
		vector<vector<int>> indicatorAMR_id(nCriterion);
		vector<vector<double>> indicatorValues(nCriterion);
		for(int i=0; i<nCriterion; ++i){
			indicatorAMR_id[i].push_back(controls.getId_cellVar(sScalIndi[i]));
			indicatorCriterion[i].push_back(stod(sScalIndiValues[i]));
			indicatorAMR_id[i].push_back(controls.getId_cellVar("x-gradient "+sGradIndi[i]));
			indicatorAMR_id[i].push_back(controls.getId_cellVar("y-gradient "+sGradIndi[i]));
			indicatorAMR_id[i].push_back(controls.getId_cellVar("z-gradient "+sGradIndi[i]));
			indicatorCriterion[i].push_back(stod(sGradIndiValues[i]));
		}
		
		
		vector<bool> boolCellRefine(mesh.cells.size(),false);
		vector<bool> boolCellUnrefine(mesh.cells.size(),true);
		vector<bool> boolCellPreserved(mesh.cells.size(),false);
		calcIndicators(mesh, controls, var, 
			maxBuffer, maxLevel, maxCells, minVolume, 
			indicatorCriterion, indicatorAMR_id, 
			boolCellRefine, boolCellUnrefine, boolCellPreserved);
		
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		// vector<string> interpolRefine_s = load.extractVector(controls.dynamicMeshMap["AMR.interpolationNames"]);
		// vector<vector<int>> interpolRefine_id(interpolRefine_s.size());
		// for(int i=0; i<interpolRefine_s.size(); ++i){
			// interpolRefine_id[i].push_back(controls.getId_cellVar(interpolRefine_s[i]));
			// interpolRefine_id[i].push_back(controls.getId_cellVar("x-gradient "+interpolRefine_s[i]));
			// interpolRefine_id[i].push_back(controls.getId_cellVar("y-gradient "+interpolRefine_s[i]));
			// interpolRefine_id[i].push_back(controls.getId_cellVar("z-gradient "+interpolRefine_s[i]));
		// }
		
		// 리파인
		{
			// 그레디언트 계산
			// solver.gradientTerms_inline(mesh, controls, var);
			// gradientTerms_AMR(mesh, controls, solver, var);
			// calc_indicatorValues(mesh, var, maxBuffer, maxLevel,
				// indicatorCriterion, indicatorAMR_id, indicatorValues);
			
			// 원래 셀의 x,y,z 저장
			vector<vector<double>> org_xyz(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i) {
				auto& cell = mesh.cells[i];
				org_xyz[i].push_back(cell.x); org_xyz[i].push_back(cell.y); org_xyz[i].push_back(cell.z);
			}
			// if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
			
			vector<vector<int>> child_new_cell_id_of_org;
			polyRefine(mesh, controls, 
				maxLevel, maxCells, minVolume, 
				indicatorCriterion, indicatorValues, 
				child_new_cell_id_of_org, 
				boolCellPreserved, boolCellRefine, boolCellUnrefine,
				0);
			
			controls.setGeometricOnlyCell_xyz(mesh);
			
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START11" << endl;
			controls.resetVariableArray(mesh, var, org_xyz, child_new_cell_id_of_org, "refine");
			// controls.resetVariableArray(mesh, var, org_xyz, child_new_cell_id_of_org, interpolRefine_id, "refine");
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START12" << endl;
			// controls.setGeometric(mesh, var);
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START13" << endl;
			// var.setSparCSR(mesh, controls);
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14" << endl;
			// solver.calcGradient.init(mesh, controls, var);
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START15" << endl;
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		// 언리파인
		{
			// 그레디언트 계산
			// solver.gradientTerms_inline(mesh, controls, var);
			// gradientTerms_AMR(mesh, controls, solver, var);
			// calc_indicatorValues(mesh, var, maxBuffer, maxLevel,
				// indicatorCriterion, indicatorAMR_id, indicatorValues);
			
			vector<vector<double>> dummy;
			// if(rank==0) cout << "| exe. Poly AMR Unrefinement" << endl;
			vector<vector<int>> child_org_cell_id_of_new;
			polyUnrefine(mesh, controls, 
				maxLevel,
				indicatorCriterion, indicatorValues, 
				child_org_cell_id_of_new, 
				boolCellPreserved, boolCellRefine, boolCellUnrefine,
				0);
			
			controls.resetVariableArray(mesh, var, dummy, child_org_cell_id_of_new, "unrefine");
			// controls.setGeometric(mesh, var);
			// var.setSparCSR(mesh, controls);
			// solver.calcGradient.init(mesh, controls, var);
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		
		// 리파티셔닝
		// if( (iter+1) %
			// stoi(controls.dynamicMeshMap["AMR.intervalRepart"]) == 0)
		{
			// vector<vector<double>> dummy;
			// if(rank==0) cout << "| exe. Dynamic Load Balancing" << endl;
			vector<int> cell_ip(mesh.cells.size(),rank);
			mesh.repartParMETIS(size, cell_ip, mesh);
			vector<int> to_new_cell_id;
			// if(rank==0) cout << "AA" << endl;
			mesh.repartitioning(cell_ip, maxLevel, to_new_cell_id);
			// if(rank==0) cout << "BB" << endl;
			
			controls.setGeometricOnlyCell_xyz(mesh);
			
			
			controls.resetVariableArray(mesh, var, cell_ip, to_new_cell_id, "repart");
			// controls.setGeometric(mesh, var);
			// var.setSparCSR(mesh, controls);
			// solver.calcGradient.init(mesh, controls, var);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		
		controls.setGeometric(mesh, var);
		var.setSparCSR(mesh, controls);
		solver.calcGradient.init(mesh, controls, var);
		
		// solver.updateProcRightCellValues_All(mesh, controls, var);
		// solver.updateCellAddiValues_All(mesh, controls, var);
		// solver.updateBoundaryFacePrimValues_All(mesh, controls, var);
		// solver.gradientTerms_All(mesh, controls, var);
		// solver.curvatureTerms_All(mesh, controls, var);
		// solver.updateProcRightCellValues_All(mesh, controls, var);
		
		solver.updateProcRightCellPrimValues(mesh, controls, var, 0);
		solver.updateCellAddiValues(mesh, controls, var, 0);
		solver.updateProcRightCellAddiValues(mesh, controls, var, 0);
		solver.updateBoundaryFacePrimValues(mesh, controls, var, 0);
		solver.gradientTerms(mesh, controls, var, 0);
		solver.updateProcRightCellGradValues(mesh, controls, var, 0);
		
	}
	
	
	
}






// void gradientTerms_AMR(MASCH_Mesh& mesh, MASCH_Control& controls, 
	// MASCH_Solver& solver, MASCH_Variables& var){
	
	
	
	

	// {
		
		// int nSp = controls.faceVar["left mass fraction"].sub_name.size();
		
		// int id_dt = controls.getId_fieldVar("time-step");
		// int id_vol = controls.getId_cellVar("volume");
		// int id_p = controls.getId_cellVar("pressure");
		// int id_u = controls.getId_cellVar("x-velocity");
		// int id_v = controls.getId_cellVar("y-velocity");
		// int id_w = controls.getId_cellVar("z-velocity");
		// int id_T = controls.getId_cellVar("temperature");
		// int id_c = controls.getId_cellVar("speed of sound");
		// int id_rho = controls.getId_cellVar("density");
		// int id_Ht = controls.getId_cellVar("total enthalpy");
		// int id_Y[2];
		// for(int i=0; i<nSp; ++i){
			// string tmp_name = controls.cellVar["mass fraction"].sub_name[i];
			// id_Y[i] = controls.getId_cellVar(tmp_name);
		// }
		// int id_drhodp = controls.getId_cellVar("density diff with pressure");
		// int id_drhodT = controls.getId_cellVar("density diff with temperature");
		// int id_dHtdp = controls.getId_cellVar("total enthalpy diff with pressure");
		// int id_dHtdT = controls.getId_cellVar("total enthalpy diff with temperature");
		// int id_drhodY[nSp];
		// int id_dHtdY[nSp];
		// for(int i=0; i<nSp; ++i){
			// string tmp_name1 = controls.cellVar["density diff with mass fraction"].sub_name[i];
			// string tmp_name2 = controls.cellVar["total enthalpy diff with mass fraction"].sub_name[i];
			// id_drhodY[i] = controls.cellVar[tmp_name1].id;
			// id_dHtdY[i] = controls.cellVar[tmp_name2].id;
		// }
		
		// solver.updateProcRightCellPrimValues(mesh, controls, var, 0);
		// // solver.updateProcRightCellAddiValues(mesh, controls, var);
		
		// auto cellVar = var.cells.data();
		// auto faceVar = var.faces.data();
		// auto procRightCellVar = var.procRightCells.data();
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto cellVar_i = cellVar[i].data();
			// double pF = cellVar_i[id_p];
			// double TF = cellVar_i[id_T];
			// double uF = cellVar_i[id_u];
			// double vF = cellVar_i[id_v];
			// double wF = cellVar_i[id_w];
			// double YF[2]; 
			// YF[0] = cellVar_i[id_Y[0]];
			// YF[1] = 1.0-YF[0];

			// double rhoi[2], ci[2], Hti[2];
			// double drhodpi[2], drhodTi[2], dHtdpi[2], dHtdTi[2];
			// {
				// solver.eosNASG(
				// 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
				// pF,uF,vF,wF,TF,rhoi[0],ci[0],Hti[0], 
				// drhodpi[0], drhodTi[0], dHtdpi[0], dHtdTi[0]);
			// }
			// {
				// solver.eosIdeal(
				// 717.0,1.4,
				// pF,uF,vF,wF,TF,rhoi[1],ci[1],Hti[1], 
				// drhodpi[1], drhodTi[1], dHtdpi[1], dHtdTi[1]);
			// }
			// double rhoF = 0.0;
			// for(int ns=0; ns<2; ++ns){
				// rhoF += YF[ns]/rhoi[ns];
			// }
			// rhoF = 1.0/rhoF;
			
			// double HtF=0.0;
			// double drhodpF=0.0;
			// double drhodTF=0.0;
			// double dHtdpF=0.0;
			// double dHtdTF=0.0;
			// double drhodYF[nSp], dHtdYF[nSp];
			// for(int ns=0; ns<2; ++ns){
				// HtF += YF[ns]*Hti[ns];
				// drhodpF += rhoF*rhoF*(YF[ns]/rhoi[ns]/rhoi[ns]*drhodpi[ns]);
				// drhodTF += rhoF*rhoF*(YF[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
				// dHtdpF += YF[ns]*dHtdpi[ns];
				// dHtdTF += YF[ns]*dHtdTi[ns];
				// drhodYF[ns] = ( -rhoF*rhoF*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]) );
				// dHtdYF[ns] = ( Hti[ns]-Hti[nSp-1] );
			// }
			
			// double cF = drhodpF + 1.0/rhoF*drhodTF/dHtdTF*(1.0-rhoF*dHtdpF);
			// cF = sqrt( 1.0 / cF );


			// cellVar_i[id_rho] = rhoF;
			// cellVar_i[id_c] = cF;
			// cellVar_i[id_Ht] = HtF;
			
			// cellVar_i[id_drhodp] = drhodpF;
			// cellVar_i[id_drhodT] = drhodTF;
			// cellVar_i[id_dHtdp] = dHtdpF;
			// cellVar_i[id_dHtdT] = dHtdTF;
			// for(int ns=0; ns<nSp; ++ns){
				// cellVar_i[id_drhodY[ns]] = drhodYF[ns];
				// cellVar_i[id_dHtdY[ns]] = dHtdYF[ns];
			// }
		// }
		// for(auto& cellVar_i : var.procRightCells){
			// double pF = cellVar_i[id_p];
			// double TF = cellVar_i[id_T];
			// double uF = cellVar_i[id_u];
			// double vF = cellVar_i[id_v];
			// double wF = cellVar_i[id_w];
			// double YF[2]; 
			// YF[0] = cellVar_i[id_Y[0]];
			// YF[1] = 1.0-YF[0];

			// double rhoi[2], ci[2], Hti[2];
			// double drhodpi[2], drhodTi[2], dHtdpi[2], dHtdTi[2];
			// {
				// solver.eosNASG(
				// 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
				// pF,uF,vF,wF,TF,rhoi[0],ci[0],Hti[0], 
				// drhodpi[0], drhodTi[0], dHtdpi[0], dHtdTi[0]);
			// }
			// {
				// solver.eosIdeal(
				// 717.0,1.4,
				// pF,uF,vF,wF,TF,rhoi[1],ci[1],Hti[1], 
				// drhodpi[1], drhodTi[1], dHtdpi[1], dHtdTi[1]);
			// }
			// double rhoF = 0.0;
			// for(int ns=0; ns<2; ++ns){
				// rhoF += YF[ns]/rhoi[ns];
			// }
			// rhoF = 1.0/rhoF;
			
			// double HtF=0.0;
			// double drhodpF=0.0;
			// double drhodTF=0.0;
			// double dHtdpF=0.0;
			// double dHtdTF=0.0;
			// double drhodYF[nSp], dHtdYF[nSp];
			// for(int ns=0; ns<2; ++ns){
				// HtF += YF[ns]*Hti[ns];
				// drhodpF += rhoF*rhoF*(YF[ns]/rhoi[ns]/rhoi[ns]*drhodpi[ns]);
				// drhodTF += rhoF*rhoF*(YF[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
				// dHtdpF += YF[ns]*dHtdpi[ns];
				// dHtdTF += YF[ns]*dHtdTi[ns];
				// drhodYF[ns] = ( -rhoF*rhoF*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]) );
				// dHtdYF[ns] = ( Hti[ns]-Hti[nSp-1] );
			// }
			
			// double cF = drhodpF + 1.0/rhoF*drhodTF/dHtdTF*(1.0-rhoF*dHtdpF);
			// cF = sqrt( 1.0 / cF );


			// cellVar_i[id_rho] = rhoF;
			// cellVar_i[id_c] = cF;
			// cellVar_i[id_Ht] = HtF;
			
			// cellVar_i[id_drhodp] = drhodpF;
			// cellVar_i[id_drhodT] = drhodTF;
			// cellVar_i[id_dHtdp] = dHtdpF;
			// cellVar_i[id_dHtdT] = dHtdTF;
			// for(int ns=0; ns<nSp; ++ns){
				// cellVar_i[id_drhodY[ns]] = drhodYF[ns];
				// cellVar_i[id_dHtdY[ns]] = dHtdYF[ns];
			// }
		// }
		
	// }
	
	
	
	
	// {
		
		// int nSp = controls.faceVar["left mass fraction"].sub_name.size();
		
		// int id_p = controls.getId_cellVar("pressure");
		// int id_rho = controls.getId_cellVar("density");
		// int id_pL = controls.getId_faceVar("left pressure");
		// int id_rhoL = controls.getId_faceVar("left density");
		// int id_pR = controls.getId_faceVar("right pressure");
		// int id_rhoR = controls.getId_faceVar("right density");
		
		// int id_area = controls.getId_faceVar("area");
		// int id_nx = controls.getId_faceVar("x unit normal");
		// int id_ny = controls.getId_faceVar("y unit normal");
		// int id_nz = controls.getId_faceVar("z unit normal");
		
		// // int nEq = controls.nEq;
		
		// auto cellVar = var.cells.data();
		// auto faceVar = var.faces.data();
		// auto procRightCellVar = var.procRightCells.data();
		
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				// int iL = face.iL;
				// auto faceVar_i = faceVar[i].data();
				// auto cellVar_iL = cellVar[iL].data();
				
				// double pF = cellVar_iL[id_p];
				// double rhoF = cellVar_iL[id_rho];
				
				// faceVar_i[id_pL] = pF;
				// faceVar_i[id_rhoL] = rhoF;
				
				// faceVar_i[id_pR]=faceVar_i[id_pL];
				// faceVar_i[id_rhoR] = faceVar_i[id_rhoL];
			// }
		// }
	// }
			
	// vector<string> tmp_names;
	// tmp_names.push_back("pressure");
	// tmp_names.push_back("density");
	
	// solver.calcGradient.leastSquare(mesh, controls, var, tmp_names, tmp_names);
	
// }
