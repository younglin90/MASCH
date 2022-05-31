
#include "./solvers.h"
#include <amgcl/profiler.hpp>

void MASCH_Solver::fvm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& solver = (*this);
	
	bool debug_bool = false;
	bool debug_AMGCL_bool = false;
	
    amgcl::profiler<> prof("fvm solver");
	
	
	
	// // 타임스텝 구하기
	// if(debug_bool) controls.log.push("calcTempSteps");
	// if(debug_AMGCL_bool) prof.tic("calcTempSteps");
	// solver.calcTempSteps(mesh, controls, var, iSegEq);
	// if(debug_bool) controls.log.pop();
	// if(debug_AMGCL_bool) prof.toc("calcTempSteps");
	
	// 고차 reconstruction
	if(debug_bool) controls.log.push("highOrderTerms");
	if(debug_AMGCL_bool) prof.tic("highOrderTerms");
	solver.highOrderTerms(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("highOrderTerms");
	
	// 페이스 루프 텀, 컨벡티브 텀 + 디퓨젼 텀
	if(debug_bool) controls.log.push("convective & diffusion Terms");
	if(debug_AMGCL_bool) prof.tic("convective & diffusion Terms");
	solver.faceLoopTerms(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("convective & diffusion Terms");
	
	// if(iSegEq==1){
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	// 셀 루프 텀, temporal 텀 + 소스텀
	if(debug_bool) controls.log.push("temporal & source Terms");
	if(debug_AMGCL_bool) prof.tic("temporal & source Terms");
	solver.cellLoopTerms(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("temporal & source Terms");
	
	// 선형 솔버
	if(debug_bool) controls.log.push("linearSystem");
	if(debug_AMGCL_bool) prof.tic("linearSystem");
	solver.linearSystem(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("linearSystem");
	
	// 셀 원시변수 업데이트
	if(debug_bool) controls.log.push("updateCellPrimValues");
	if(debug_AMGCL_bool) prof.tic("updateCellPrimValues");
	solver.updateCellPrimValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateCellPrimValues");
	
	// 선형 시스템 0.0으로 초기화
	if(debug_bool) controls.log.push("clearLinearSystems");
	if(debug_AMGCL_bool) prof.tic("clearLinearSystems");
	var.clearLinearSystems(iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("clearLinearSystems");
		
		
	
	
	// proc right cell 로 셀의 원시변수 넘기기
	if(debug_bool) controls.log.push("updateProcRightCellPrimValues");
	if(debug_AMGCL_bool) prof.tic("updateProcRightCellPrimValues");
	solver.updateProcRightCellPrimValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateProcRightCellPrimValues");
	
	// 원시변수 제외한 나머지 셀값 업데이트
	if(debug_bool) controls.log.push("updateCellAddiValues");
	if(debug_AMGCL_bool) prof.tic("updateCellAddiValues");
	solver.updateCellAddiValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateCellAddiValues");
	
	// 셀 원시변수 제외한 나머지 proc 셀값 업데이트
	if(debug_bool) controls.log.push("updateProcRightCellAddiValues");
	if(debug_AMGCL_bool) prof.tic("updateProcRightCellAddiValues");
	solver.updateProcRightCellAddiValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateProcRightCellAddiValues");
	
	// B.C. 원시변수 업데이트
	if(debug_bool) controls.log.push("updateBoundaryFacePrimValues");
	if(debug_AMGCL_bool) prof.tic("updateBoundaryFacePrimValues");
	solver.updateBoundaryFacePrimValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateBoundaryFacePrimValues");
	
	// 셀 그레디언트
	if(debug_bool) controls.log.push("gradientTerms");
	if(debug_AMGCL_bool) prof.tic("gradientTerms");
	solver.gradientTerms(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("gradientTerms");
	
	// 셀 그레디언트 proc right cell 로 넘기기
	if(debug_bool) controls.log.push("updateProcRightCellGradValues");
	if(debug_AMGCL_bool) prof.tic("updateProcRightCellGradValues");
	solver.updateProcRightCellGradValues(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("updateProcRightCellGradValues");
	
	// 셀 곡률
	if(debug_bool) controls.log.push("curvatureTerms");
	if(debug_AMGCL_bool) prof.tic("curvatureTerms");
	solver.curvatureTerms(mesh, controls, var, iSegEq);
	if(debug_bool) controls.log.pop();
	if(debug_AMGCL_bool) prof.toc("curvatureTerms");
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	if (rank == 0 && debug_AMGCL_bool) cout << prof << endl;
	
	
	
	
	
	
	// // 시간 업데이트
	// var.fields[controls.fieldVar["time"].id] += 
	// var.fields[controls.fieldVar["time-step"].id];
	
	// // controls.log.push("linear solver ZERO");
	// fill(var.Avalues.begin(),var.Avalues.end(),0.0);
	// fill(var.Bvalues.begin(),var.Bvalues.end(),0.0);
	// fill(var.Xvalues.begin(),var.Xvalues.end(),0.0);
	// // controls.log.pop();
	
	// // 셀 루프 텀, temporal 텀 + 소스텀
	// // controls.log.push("temporal & source Terms");
	// solver.cellLoopTerms(mesh, controls, var);
	// // controls.log.pop();
	// // 페이스 루프 텀, 컨벡티브 텀 + 디퓨젼 텀
	// // controls.log.push("convective & diffusion Terms");
	// solver.faceLoopTerms(mesh, controls, var);
	// // controls.log.pop();
	
	// // 선형 솔버
	// // controls.log.push("linearSystem");
	// solver.linearSystem(mesh, controls, var);
	// // controls.log.pop();
	
	// // 셀 원시변수 업데이트
	// // controls.log.push("updateCellPrimValues");
	// solver.updateCellPrimValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // 선형 솔버 값들 0.0 으로 만들어주기
	// // controls.log.push("clearLinearSystems");
	// var.clearLinearSystems();
	// // controls.log.pop();
	
	// // 타임스텝 구하기
	// // controls.log.push("calcTempSteps");
	// solver.calcTempSteps(mesh, controls, var);
	// // controls.log.pop();
	
	// // 셀 그레디언트
	// // controls.log.push("gradientTerms");
	// solver.gradientTerms(mesh, controls, var);
	// // controls.log.pop();
	
	// // 고차 reconstruction
	// // controls.log.push("highOrderTerms");
	// solver.highOrderTerms(mesh, controls, var);
	// // controls.log.pop();
	
	// // 원시변수 제외한 나머지 셀값 업데이트
	// // controls.log.push("updateCellAddiValues");
	// solver.updateCellAddiValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // proc 셀의 원시변수 mpi 넘기기
	// // controls.log.push("updateProcRightCellPrimValues");
	// solver.updateProcRightCellPrimValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // 셀 원시변수 제외한 나머지 proc 셀값 업데이트
	// // controls.log.push("updateProcRightCellAddiValues");
	// solver.updateProcRightCellAddiValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // B.C. 원시변수 업데이트
	// // controls.log.push("updateBoundaryFacePrimValues");
	// solver.updateBoundaryFacePrimValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // B.C. 원시변수 제외한 나머지 값 업데이트
	// // controls.log.push("updateBoundaryFaceAddiValues");
	// solver.updateBoundaryFaceAddiValues(mesh, controls, var);
	// // controls.log.pop();
	
	// // old 값 업데이트
	// // controls.log.push("updateOldValues");
	// solver.updateOldValues(mesh, controls, var);
	// // controls.log.pop();
	
	
}


void MASCH_Solver::gradientTerms(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	auto& solver = (*this);
	
	if(solver.gradLSIds_cell_name.size()==0) return;
	
	// solver.calcGradient.leastSquare(mesh, controls, var, 
		// solver.gradLSIds_cell_name[iSegEq], solver.gradLSIds_bcFace_name[iSegEq]);
	// solver.calcGradient.leastSquare(mesh, controls, var, 
		// solver.gradLSIds_cell_name[iSegEq], solver.gradLSIds_bcFace_name[iSegEq], 
		// solver.minmaxInp_cell_name[iSegEq], solver.maxOut_cell_name[iSegEq], solver.minOut_cell_name[iSegEq]);
	solver.calcGradient.leastSquare(mesh, controls, var, 
		solver.gradLSIds_cell_name[iSegEq], solver.gradLSIds_bcFace_name[iSegEq], 
		solver.minmaxInp_cell_name[iSegEq], solver.maxOut_cell_name[iSegEq], solver.minOut_cell_name[iSegEq], 
		solver.minmaxInp_point_name[iSegEq], solver.maxOut_point_name[iSegEq], solver.minOut_point_name[iSegEq]);
	
}

void MASCH_Solver::curvatureTerms(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	auto& solver = (*this);
	
	if(solver.curvatureIds_cell_name.size()==0) return;
	
	solver.calcCurvature(mesh, controls, var, 
		solver.curvatureIds_cell_name[iSegEq]);
	
}

void MASCH_Solver::highOrderTerms(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	auto& solver = (*this);
	
	

	// 하이오더 리컨스트럭션 + 추가적 변수들
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto pointVar = var.points.data();
	auto procRightCellVar = var.procRightCells.data();
	auto calcAddSolPtr = solver.calcFaceAddiVal.data();
	int calcAddSolSize = solver.calcFaceAddiVal.size();
	auto sol = solver.calcHO_FaceVal[iSegEq];
	auto solAddi = solver.calcFaceAddiVal[iSegEq];
    
    int point_max_min_size = solver.minmaxInp_point_name[iSegEq].size();
	vector<int> id_oup_point_max, id_oup_point_min;
	for(auto& item : solver.minmaxInp_point_name[iSegEq]){
		id_oup_point_max.push_back(controls.getId_pointVar("maximum "+item));
		id_oup_point_min.push_back(controls.getId_pointVar("minimum "+item));
	}
    
	for(int i=0, ip=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType()==MASCH_Face_Types::BOUNDARY) continue;
		
		int iL = face.iL; int iR = face.iR;
		auto faceVar_i = faceVar[i].data();
		auto cellVar_iL = cellVar[iL].data();
        
        int point_size = face.ipoints.size();
        vector<vector<double>> point_xyz(point_size,vector<double>(3,0.0));
        vector<vector<double>> point_max(point_size,vector<double>(point_max_min_size,0.0));
        vector<vector<double>> point_min(point_size,vector<double>(point_max_min_size,0.0));
        int tmp_iter = 0;
        for(auto& ipoint : face.ipoints){
            point_xyz[tmp_iter][0] = mesh.points[ipoint].x;
            point_xyz[tmp_iter][1] = mesh.points[ipoint].y;
            point_xyz[tmp_iter][2] = mesh.points[ipoint].z;
            
            auto pointVar_i = pointVar[ipoint].data();
            
            for(int j=0; j<point_max_min_size; ++j){
                point_max[tmp_iter][j] = pointVar_i[id_oup_point_max[j]];
                point_min[tmp_iter][j] = pointVar_i[id_oup_point_min[j]];
            }
            ++tmp_iter;
        }
        
		// 하이 오더 리컨스트럭션
		if(face.getType()==MASCH_Face_Types::INTERNAL){
			// for(auto& sol : solver.calcHO_FaceVal){
				sol(var.fields.data(),cellVar_iL,cellVar[iR].data(),faceVar_i,
                    point_xyz, point_max, point_min);
			// }
		}
		else if(face.getType()==MASCH_Face_Types::PROCESSOR){
			// for(auto& sol : solver.calcHO_FaceVal){
				sol(var.fields.data(),cellVar_iL,procRightCellVar[ip++].data(),faceVar_i,
                    point_xyz, point_max, point_min);
			// }
		}
		
		// 추가적인 변수들
		// for(auto& sol : solver.calcFaceAddiVal){
			solAddi(faceVar[i].data());
		// }
		// for(int iSol=0; iSol<calcAddSolSize; ++iSol){
			// calcAddSolPtr[iSol](faceVar[i].data());
		// }
		
		
	}
	
	// auto& solver = (*this);
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto procRightCellsVar = var.procRightCells.data();
	
	// int funct_size = calcHO_FaceVal.size();
	// auto funct_ptr = calcHO_FaceVal.data();
	
	// int functAddi_size = calcFaceAddiVal.size();
	// auto functAddi_ptr = calcFaceAddiVal.data();
	
	// // for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		// // int iL = faces[i].iL;
		// // int iR = faces[i].iR;
		// // for(int iFunct=0; iFunct<funct_size; ++iFunct){
			// // funct_ptr[iFunct](cellVar[iL].data(), cellVar[iR].data(), faceVar[i].data());
		// // }
		// // for(int j=0; j<functAddi_size; ++j){
			// // functAddi_ptr[j](faceVar[i].data());
		// // }
	// // }
	
	// // // MPI_Barrier(MPI_COMM_WORLD);
	// // // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// // for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		// // auto& boundary = mesh.boundaries[i];
		// // int str = boundary.startFace;
		// // int end = str + boundary.nFaces;
		// // if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			// // for(int i=str; i<end; ++i, ++ip){
				// // auto& face = faces[i];
				// // int iL = face.iL;
				// // for(int iFunct=0; iFunct<funct_size; ++iFunct){
					// // funct_ptr[iFunct](cellVar[iL].data(), 
							// // procRightCellsVar[ip].data(), faceVar[i].data());
				// // }
				// // for(int j=0; j<functAddi_size; ++j){
					// // functAddi_ptr[j](faceVar[i].data());
				// // }
			// // }
		// // }
	// // }
	
}


void MASCH_Solver::cellLoopTerms(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	auto& solver = (*this);
	

	// 시간텀
	int nEq = controls.nEq[iSegEq];
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellVar = var.procRightCells.data();
	auto sol = solver.calcTemporal[iSegEq];
	for(int i=0; i<mesh.cells.size(); ++i){
		auto cellVar_i = cellVar[i].data();
		
		for(int iEq=0; iEq<nEq; ++iEq){
			for(int jEq=0; jEq<nEq; ++jEq){
				fluxA[iEq*nEq+jEq] = 0.0;
			}
			fluxB[iEq] = 0.0;
		}
		
		// for(auto item : solver.calcTemporal){
			sol(cellVar_i, fieldVar, fluxA, fluxB);
		// }
		for(int iEq=0; iEq<nEq; ++iEq){
			for(int jEq=0; jEq<nEq; ++jEq){
				var.accumSparD( iSegEq, i, iEq, jEq, fluxA[iEq*nEq+jEq] );
			}
			var.accumB( iSegEq, i, iEq, fluxB[iEq] );
		}
	}


	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto calcTemporal = solver.calcTemporal.data();
	// auto calcSource = solver.calcSource.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto fieldVar = var.fields.data();
	
	// int nEq = controls.nEq;
	// double fluxA[nEq*nEq];
	// double fluxB[nEq];
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto& cell = cells[i];
		
		// // 시간텀
		// for(auto item : solver.calcTemporal){
			// item(cellVar[i].data(), fieldVar, fluxA, fluxB);
		// }
		// for(int jEq=0; jEq<nEq; ++jEq){
			// for(int iEq=0; iEq<nEq; ++iEq){
				// // var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
				// var.Avalues[nEq*nEq*i+nEq*jEq+iEq] += fluxA[jEq*nEq+iEq];
			// }
			// // var.accumB( i, jEq, fluxB[jEq] );
			// var.Bvalues[nEq*i+jEq] += fluxB[jEq];
		// }
		
		// // // 소스텀
		// // for(auto item : solver.calcSource){
			// // item(cellVar[i].data(), fluxA, fluxB);
		// // }
		// // for(int jEq=0; jEq<nEq; ++jEq){
			// // for(int iEq=0; iEq<nEq; ++iEq){
				// // // var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
				// // var.Avalues[nEq*nEq*i+nEq*jEq+iEq] += fluxA[jEq*nEq+iEq];
			// // }
			// // // var.accumB( i, jEq, fluxB[jEq] );
			// // var.Bvalues[nEq*i+jEq] += fluxB[jEq];
		// // }
	// }
}




void MASCH_Solver::faceLoopTerms(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	auto& solver = (*this);

	// 컨벡티브 + 디퓨젼 풀기
	int nEq = controls.nEq[iSegEq];
	// double fluxAL_conv[nEq*nEq], fluxAL_lapl[nEq*nEq], fluxAL_nLapl[nEq*nEq];
	// double fluxAR_conv[nEq*nEq], fluxAR_lapl[nEq*nEq], fluxAR_nLapl[nEq*nEq];
	// double fluxB_conv[nEq], fluxB_lapl[nEq], fluxB_nLapl[nEq];
	double fluxA_LL[nEq*nEq];
	double fluxA_RR[nEq*nEq];
	double fluxA_LR[nEq*nEq];
	double fluxA_RL[nEq*nEq];
	double fluxB_LL[nEq];
	double fluxB_RR[nEq];
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellVar = var.procRightCells.data();
	auto sol_conv = solver.calcConvFlux[iSegEq];
	// auto sol_lapl = solver.calcLaplFlux[iSegEq];
	// auto sol_nLapl = solver.calcNLaplFlux[iSegEq];
	auto solImplicit = checkImplicit[iSegEq];
	
	// if(iSegEq==1){
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	// if(iSegEq==0) solImplicit = false;
	// if(iSegEq==1) solImplicit = true;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		auto& face = mesh.faces[i];
		int iL = face.iL;
		int iR = face.iR;
	
		for(int iEq=0; iEq<nEq; ++iEq){
			for(int jEq=0; jEq<nEq; ++jEq){
				fluxA_LL[iEq*nEq+jEq] = 0.0;
				fluxA_RR[iEq*nEq+jEq] = 0.0;
				fluxA_LR[iEq*nEq+jEq] = 0.0;
				fluxA_RL[iEq*nEq+jEq] = 0.0;
			}
			fluxB_LL[iEq] = 0.0;
			fluxB_RR[iEq] = 0.0;
		}
		
		// 컨벡티브 텀
		sol_conv(fieldVar, cellVar[iL].data(), cellVar[iR].data(), 
			faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB_LL, fluxB_RR);
		// // 디퓨젼 텀 (라플라스)
		// sol_lapl(cellVar[iL].data(), cellVar[iR].data(), 
			// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
		// // 디퓨젼 텀 (논-라플라스)
		// sol_nLapl(cellVar[iL].data(), cellVar[iR].data(), 
			// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
		
		// A sparse matrix 에 넣기
		if(solImplicit==true){
			for(int iEq=0; iEq<nEq; ++iEq){
				for(int jEq=0; jEq<nEq; ++jEq){
					var.accumSparD( iSegEq, iL, iEq, jEq, fluxA_LL[iEq*nEq+jEq] );
					var.accumSparD( iSegEq, iR, iEq, jEq, fluxA_RR[iEq*nEq+jEq] );
					var.accumSparLR( iSegEq, i, iL, iEq, jEq, fluxA_LR[iEq*nEq+jEq] );
					var.accumSparRL( iSegEq, i, iR, iEq, jEq, fluxA_RL[iEq*nEq+jEq] );
				}
			}
		}
		
		for(int iEq=0; iEq<nEq; ++iEq){
			var.accumB( iSegEq, iL, iEq, fluxB_LL[iEq] );
			var.accumB( iSegEq, iR, iEq, fluxB_RR[iEq] );
		}
	}
	for(int ibc=0, iter=0, ip=0, SIZE=mesh.boundaries.size(); ibc<SIZE; ++ibc){
		auto& boundary = mesh.boundaries[ibc];
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		// double* cellVar_iR;
		
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			auto& sol_conv_bc_prim = solver.calcConvFlux_BC[iSegEq][iter];
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				int iL = face.iL;
				
				for(int iEq=0; iEq<nEq; ++iEq){
					for(int jEq=0; jEq<nEq; ++jEq){
						fluxA_LL[iEq*nEq+jEq] = 0.0;
					}
					fluxB_LL[iEq] = 0.0;
				}
				
				// 컨벡티브 텀
				for(auto& sol_conv_bc : sol_conv_bc_prim){
					sol_conv_bc(fieldVar, cellVar[iL].data(), 
						faceVar[i].data(), fluxA_LL, fluxB_LL);
				}
				// sol_conv(cellVar[iL].data(), nullptr, 
					// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
				// // 디퓨젼 텀 (라플라스)
				// sol_lapl(cellVar[iL].data(), nullptr, 
					// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
				// // 디퓨젼 텀 (논-라플라스)
				// sol_nLapl(cellVar[iL].data(), nullptr, 
					// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
				
				// A sparse matrix 에 넣기
				if(solImplicit==true){
					for(int iEq=0; iEq<nEq; ++iEq){
						for(int jEq=0; jEq<nEq; ++jEq){
							var.accumSparD( iSegEq, iL, iEq, jEq, fluxA_LL[iEq*nEq+jEq] );
						}
					}
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iSegEq, iL, iEq, fluxB_LL[iEq] );
				}
			}
			++iter;
		}
		else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				int iL = face.iL;
				
				for(int iEq=0; iEq<nEq; ++iEq){
					for(int jEq=0; jEq<nEq; ++jEq){
						fluxA_LL[iEq*nEq+jEq] = 0.0;
						fluxA_RR[iEq*nEq+jEq] = 0.0;
						fluxA_LR[iEq*nEq+jEq] = 0.0;
						fluxA_RL[iEq*nEq+jEq] = 0.0;
					}
					fluxB_LL[iEq] = 0.0;
				}
				
				// 컨벡티브 텀
				sol_conv(fieldVar, cellVar[iL].data(), procRightCellVar[ip].data(), 
					faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB_LL, fluxB_RR);
				// // 디퓨젼 텀 (라플라스)
				// sol_lapl(cellVar[iL].data(), procRightCellVar[ip].data(), 
					// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
				// // 디퓨젼 텀 (논-라플라스)
				// sol_nLapl(cellVar[iL].data(), procRightCellVar[ip].data(), 
					// faceVar[i].data(), fluxA_LL, fluxA_RR, fluxA_LR, fluxA_RL, fluxB);
				
				// A sparse matrix 에 넣기
				if(solImplicit==true){
					for(int iEq=0; iEq<nEq; ++iEq){
						for(int jEq=0; jEq<nEq; ++jEq){
							var.accumSparD( iSegEq, iL, iEq, jEq, fluxA_LL[iEq*nEq+jEq] );
							var.accumSparLR( iSegEq, i, iL, iEq, jEq, fluxA_LR[iEq*nEq+jEq] );
						}
					}
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iSegEq, iL, iEq, fluxB_LL[iEq] );
				}
				
				++ip;
			}
		}
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// auto& solver = (*this);
	
	// int cellSize = mesh.cells.size();
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto procRightCellsVar = var.procRightCells.data();
	// int nEq = controls.nEq;
	
	// auto convFunc = solver.calcConvFlux.data();
	// auto calcLaplFlux = solver.calcLaplFlux.data();
	// auto calcNLaplFlux = solver.calcNLaplFlux.data();
	
	// double fluxA[nEq*nEq];
	// double fluxB[nEq];
	// double fluxLaplA[nEq*nEq];
	// double fluxLaplB[nEq];
	// double fluxNLaplA[nEq*nEq];
	// double fluxNLaplB[nEq];
	
	
	// for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		// auto& face = faces[i];
		// int iL = face.iL;
		// int iR = face.iR;
		
		// // 컨벡티브 텀
		// convFunc[0](cellVar[iL].data(), cellVar[iR].data(), 
					// faceVar[i].data(), fluxA, fluxB);
		
		// if(checkImplicitConvFlux[0]==true){
			// // A sparse matrix 에 넣기
		// }
		
		// for(int iEq=0; iEq<nEq; ++iEq){
			// // var.accumB( iL, iEq, +fluxB[iEq] );
			// // var.accumB( iR, iEq, -fluxB[iEq] );
			// var.Bvalues[nEq*iL+iEq] += (+fluxB[iEq]);
			// var.Bvalues[nEq*iR+iEq] += (-fluxB[iEq]);
		// }
		
		// // // 디퓨젼 텀
		// // calcLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					// // faceVar[i].data(), fluxLaplA, fluxLaplB);
		// // calcNLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					// // faceVar[i].data(), fluxNLaplA, fluxNLaplB);
		
		// // if(checkImplicitConvFlux[0]==true){
			// // // A sparse matrix 에 넣기
		// // }
		
		// // for(int iEq=0; iEq<nEq; ++iEq){
			// // // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
			// // // var.accumB( iR, iEq, -fluxLaplB[iEq]-fluxNLaplB[iEq] );
			// // var.Bvalues[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
			// // var.Bvalues[nEq*iR+iEq] += (-fluxLaplB[iEq]-fluxNLaplB[iEq]);
		// // }
		
		
	// }
	
	
	// for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		// auto& boundary = mesh.boundaries[i];
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// int iL = face.iL;
				
				// // 컨벡티브 텀
				// convFunc[0](cellVar[iL].data(), nullptr, 
							// faceVar[i].data(), fluxA, fluxB);
							
				// if(checkImplicitConvFlux[0]==true){
					// // A sparse matrix 에 넣기
				// }
				
				// double tmp_value=0.0;
				// for(int iEq=0; iEq<nEq; ++iEq){
					// // var.accumB( iL, iEq, +fluxB[iEq] );
					// var.Bvalues[nEq*iL+iEq] += (+fluxB[iEq]);
					// tmp_value+=( +fluxB[iEq]);
				// }
				
			// // double test1 = tmp_value;
			// // if(isnan(test1) || test1<-1.e12 || test1>1.e12){
				// // cout << test1 << endl;
				// // cout << fluxB[0] << endl;
				// // cout << fluxB[1] << endl;
				// // cout << fluxB[2] << endl;
				// // cout << fluxB[3] << endl;
				// // cout << fluxB[4] << endl;
			// // }
				// // // 디퓨젼 텀
				// // calcLaplFlux[0](cellVar[iL].data(), nullptr, 
							// // faceVar[i].data(), fluxLaplA, fluxLaplB);
				// // calcNLaplFlux[0](cellVar[iL].data(), nullptr, 
							// // faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				// // if(checkImplicitConvFlux[0]==true){
					// // // A sparse matrix 에 넣기
				// // }
				
				// // for(int iEq=0; iEq<nEq; ++iEq){
					// // // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
					// // var.Bvalues[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
				// // }
				
				
			// }
		// }
		// else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// int iL = face.iL;
				// int iR = i-str;
				
				// // 컨벡티브 텀
				// convFunc[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							// faceVar[i].data(), fluxA, fluxB);
							
				// if(checkImplicitConvFlux[0]==true){
					// // A sparse matrix 에 넣기
				// }
				// for(int iEq=0; iEq<nEq; ++iEq){
					// // var.accumB( iL, iEq, +fluxB[iEq] );
					// var.Bvalues[nEq*iL+iEq] += ( +fluxB[iEq]);
				// }
				
				// // // 디퓨젼 텀
				// // calcLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							// // faceVar[i].data(), fluxLaplA, fluxLaplB);
				// // calcNLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							// // faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				// // if(checkImplicitConvFlux[0]==true){
					// // // A sparse matrix 에 넣기
				// // }
				
				// // for(int iEq=0; iEq<nEq; ++iEq){
					// // // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
					// // var.Bvalues[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
				// // }
				
				// ++ip;
			// }
		// }
	// }
	
	
}



void MASCH_Solver::linearSystem(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	// controls.log.push("solve ddd0");
	auto& solver = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();

	int nEq = controls.nEq[iSegEq];
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellVar = var.procRightCells.data();
	auto solImplicit = checkImplicit[iSegEq];
	
	if(solImplicit==true){
		
		// AMGCL
		solveAMGCL(
			var.i_str_CSR[iSegEq], var.j_displ_CSR[iSegEq], var.Avalues[iSegEq], 
			var.Bvalues[iSegEq], var.Xvalues[iSegEq]);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	}
	else{
		
		// cout << "AAAAAAAAAAAA" << endl;
		for(int i=0; i<mesh.cells.size(); ++i){
			auto cellVar_i = cellVar[i].data();
			
			vector<vector<double>> matA(nEq,vector<double>(nEq,0.0));
			
			for(int iEq=0; iEq<nEq; ++iEq){
				for(int jEq=0; jEq<nEq; ++jEq){
					matA[iEq][jEq] = var.getSparD( iSegEq, i, iEq, jEq );
				}
			}
			
			math.GaussSeidelSOR(matA);
			
			for(int iEq=0; iEq<nEq; ++iEq){
				double tmp_val = 0.0;
				for(int jEq=0; jEq<nEq; ++jEq){
					// tmp_val += matA[iEq][jEq]*var.Bvalues[nEq*i+jEq];
					// tmp_val += matA[iEq][jEq]*var.Bvalues[mesh.cells.size()*jEq+i];
					tmp_val += matA[iEq][jEq]*var.getB( iSegEq, i, jEq );
				}
				// var.Xvalues[iSegEq][nEq*i+iEq] = tmp_val;
				var.setX(iSegEq, i, iEq, tmp_val);
			}
		}
	}
	
}

void MASCH_Solver::updateCellPrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	auto& solver = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();

	int nEq = controls.nEq[iSegEq];
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellVar = var.procRightCells.data();
	int id_resi = controls.fieldVar["residual"].id;
	int id_vol = controls.getId_cellVar("volume");
	auto sol = solver.calcUpdatePrim[iSegEq];
	
	double Xvalues[nEq];
	
	double tmp_volume = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		double volume = var.cells[i][id_vol];
		auto cellVar_i = cellVar[i].data();
	
		for(int iEq=0; iEq<nEq; ++iEq){
			Xvalues[iEq] = var.getX(iSegEq, i, iEq);
		}
		
		// for(auto& sol : solver.calcUpdatePrim){
			// sol(cellVar_i,&(var.Xvalues[iSegEq][nEq*i]));
			sol(fieldVar, cellVar_i, Xvalues);
		// }
		
		for(int iEq=0; iEq<nEq; ++iEq){
			double tmp_resi = Xvalues[iEq]*volume;
			fieldVar[id_resi] += tmp_resi*tmp_resi;
			tmp_volume += volume;
		}
	}
	if(size>1){
		vector<double> tmp_send(2,0.0);
		tmp_send[0] = (var.fields[id_resi]);
		tmp_send[1] = (tmp_volume);
		vector<double> tmp_recv(2,0.0);
		MPI_Allreduce(tmp_send.data(), tmp_recv.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		var.fields[id_resi] = tmp_recv[0];
		tmp_volume = tmp_recv[1];
	}
	var.fields[id_resi] = sqrt(var.fields[id_resi]/tmp_volume);
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// int cellSize = mesh.cells.size();
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto fieldsVar = var.fields.data();
	// int nEq = controls.nEq;
	// auto primVarIds = controls.primVarIds.data();
	// // function<void*(int iEq, double& value)> limitPrim;
	// // auto limitPrim = &MASCH_Control::limitPrim;
	// // auto limitPrim = &controls.limitPrim;
	// // void(*limitPrim)(int, double&) = &controls.limitPrim;
	
	// int id_resi = controls.fieldVar["residual"].id;
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto& cell = cells[i];
		// auto cellVar_i = cellVar[i].data();
		// for(int iEq=0; iEq<nEq; ++iEq){
			// // cellVar_i[iEq] += var.getX(i, iEq);
			// int id = primVarIds[iEq];
			// cellVar_i[id] += var.Xvalues[nEq*i+iEq];
			
			// fieldsVar[id_resi] += var.Xvalues[nEq*i+iEq]*var.Xvalues[nEq*i+iEq];
			
			// // if(iEq==3 && abs(cellVar_i[id])>1.e-16) {
				// // cout << cellVar_i[id] << " " << var.Xvalues[nEq*i+iEq] << endl;
			// // }		
			
			// // controls.limitPrim(iEq, cellVar_i[iEq]);
			// // (*&limitPrim)(iEq, cellVar_i[iEq]);
		// }
	// }
	
	// if(size>1){
		// double tmp_fieldVar = fieldsVar[id_resi];
		// double tmp_fieldVar_glo;
		// MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// fieldsVar[id_resi] = tmp_fieldVar_glo;
	// }
	
	
}


void MASCH_Solver::updateCellAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	auto& solver = (*this);
	
	auto sol = solver.calcCellAddiVal[iSegEq];
	
	// cell 추가적 변수
	auto cellVar = var.cells.data();
	for(int i=0; i<mesh.cells.size(); ++i){
		auto cellVar_i = cellVar[i].data();
		sol(cellVar_i);
	}
	
		// // cout << var.faces[0].size() << endl;
	// auto cells = mesh.cells.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto p_calcCellAddiVal = calcCellAddiVal.data();
	
	// int size_c = calcCellAddiVal.size();
	// int size_f = calcFaceAddiVal.size();
	
		// // int id_test = controls.cellVar["speed of sound"].id;
		
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto cellVar_i = cellVar[i].data();
		// for(int j=0; j<size_c; ++j){
			// p_calcCellAddiVal[j](cellVar_i);
		// }
		// // cout << cellVar_i[id_test] << endl;
	// }
	
	
}


void MASCH_Solver::updateProcRightCellPrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	if(size>1){
		int proc_size = var.procRightCells.size();
		int prim_size = controls.sendProcValueNames[iSegEq].size();
		vector<int> ids(prim_size);
		for(int i=0; i<prim_size; ++i){
			string name = controls.sendProcValueNames[iSegEq][i];
			ids[i] = controls.getId_cellVar(name);
		}
		
		vector<double> send_value;
		send_value.reserve(proc_size*prim_size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					auto cellVar_i = cellVar[faces[i].iL].data();
					auto faceVar_i = faceVar[i].data();
					for(int j=0; j<prim_size; ++j){
						send_value.push_back(cellVar_i[ids[j]]);
					}
				}
			}
		}
		
		vector<int> tmp_sendCounts(size,0);
		vector<int> tmp_recvCounts(size,0);
		vector<int> tmp_sendDispls(size+1,0);
		vector<int> tmp_recvDispls(size+1,0);
		for(int ip=0; ip<size; ++ip){
			tmp_sendCounts[ip] = mesh.countsSendProcFaces[ip]*prim_size;
			tmp_recvCounts[ip] = mesh.countsRecvProcFaces[ip]*prim_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_sendDispls[ip+1] = tmp_sendDispls[ip] + tmp_sendCounts[ip];
			tmp_recvDispls[ip+1] = tmp_recvDispls[ip] + tmp_recvCounts[ip];
		}
		
		vector<double> recv_value(proc_size*prim_size);
		MPI_Alltoallv( send_value.data(), tmp_sendCounts.data(), 
						tmp_sendDispls.data(), MPI_DOUBLE, 
						recv_value.data(), tmp_recvCounts.data(), 
						tmp_recvDispls.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		auto recv_value_ptr = recv_value.data();
		int iter=0;
		for(auto& cells : var.procRightCells){
			auto cellVar_i = cells.data();
			for(int j=0; j<prim_size; ++j){
				cellVar_i[ids[j]] = recv_value_ptr[iter++];
			}
		}
	}
	
	
	
}


void MASCH_Solver::updateProcRightCellGradValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& inp_s = (*this).gradLSIds_cell_name[iSegEq];
	int inp_size = inp_s.size();
	
	if(inp_size==0) return;
	
	vector<int> id_oup;
	{

		vector<string> vec_oup0, vec_oup1, vec_oup2, vec_bc_face;
		for(auto& item : inp_s){
			string tmp_oup0 = "x-gradient "; tmp_oup0+=item;
			string tmp_oup1 = "y-gradient "; tmp_oup1+=item;
			string tmp_oup2 = "z-gradient "; tmp_oup2+=item;
			// string Bvaluesc_face = "left "; Bvaluesc_face+=item;
			
			vec_oup0.push_back(tmp_oup0);
			vec_oup1.push_back(tmp_oup1);
			vec_oup2.push_back(tmp_oup2);
			// vec_bc_face.push_back(Bvaluesc_face);
		}
		
		
		vector<int> id_inp, id_oup0, id_oup1, id_oup2, id_bc_face;
		
		int iter_tmp = 0;
		for(auto& item : inp_s){
			int inp = controls.getId_cellVar(item);
			int oup0 = controls.getId_cellVar(vec_oup0[iter_tmp]);
			int oup1 = controls.getId_cellVar(vec_oup1[iter_tmp]);
			int oup2 = controls.getId_cellVar(vec_oup2[iter_tmp]);
			// int bc_face = controls.getId_faceVar(vec_bc_face[iter_tmp]);
			
			// id_inp.push_back(inp);
			id_oup.push_back(oup0);
			id_oup.push_back(oup1);
			id_oup.push_back(oup2);
			// id_bc_face.push_back(bc_face);
			
			++iter_tmp;
		}

		
		
	}
	
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	if(size>1){
		auto id_oup_ptr = id_oup.data();
		int proc_size = var.procRightCells.size();
		
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
						send_value.push_back(cellVar_i[id_oup_ptr[j]]);
					}
				}
			}
		}
		
		vector<int> tmp_sendCounts(size,0);
		vector<int> tmp_recvCounts(size,0);
		vector<int> tmp_sendDispls(size+1,0);
		vector<int> tmp_recvDispls(size+1,0);
		for(int ip=0; ip<size; ++ip){
			tmp_sendCounts[ip] = mesh.countsSendProcFaces[ip]*inp_size;
			tmp_recvCounts[ip] = mesh.countsRecvProcFaces[ip]*inp_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_sendDispls[ip+1] = tmp_sendDispls[ip] + tmp_sendCounts[ip];
			tmp_recvDispls[ip+1] = tmp_recvDispls[ip] + tmp_recvCounts[ip];
		}
		
		vector<double> recv_value(proc_size*inp_size);
		MPI_Alltoallv( send_value.data(), tmp_sendCounts.data(), 
						tmp_sendDispls.data(), MPI_DOUBLE, 
						recv_value.data(), tmp_recvCounts.data(), 
						tmp_recvDispls.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		auto recv_value_ptr = recv_value.data();
		int iter=0;
		for(auto& cells : var.procRightCells){
			auto cellVar_i = cells.data();
			for(int j=0; j<inp_size; ++j){
				cellVar_i[id_oup_ptr[j]] = recv_value_ptr[iter++];
			}
		}
	}
	
	
}


void MASCH_Solver::updateProcRightCellAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	auto& solver = (*this);
	
	auto sol = solver.calcCellAddiVal[iSegEq];
	// proc right cell 추가적 변수
	for(auto& cellVar_i : var.procRightCells){
		// for(auto& sol : solver.calcCellAddiVal){
			sol(cellVar_i.data());
		// }
	}
		
	// auto cells = mesh.cells.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto p_calcCellAddiVal = calcCellAddiVal.data();
	
	// int size_c = calcCellAddiVal.size();
	// int size_f = calcFaceAddiVal.size();
	
	// for(auto& cells : var.procRightCells){
		// auto cellVar_i = cells.data();
		// for(int j=0; j<size_c; ++j){
			// p_calcCellAddiVal[j](cellVar_i);
		// }
	// }
	
	
}

void MASCH_Solver::updateBoundaryFacePrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	auto& solver = (*this);
	// 바운더리 페이스 값들 정하기
	int id_t = controls.getId_fieldVar("time");
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellVar = var.procRightCells.data();
    
    
    
    // // update time varying boundary values
    // solver.updateTimeVaryingMappedFixedValue(mesh, controls, var);
    
    
    
	
	auto sol = solver.calcFaceAddiVal[iSegEq];
	
	for(int ibc=0, iter=0, SIZE=mesh.boundaries.size(); ibc<SIZE; ++ibc){
		auto& boundary = mesh.boundaries[ibc];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		auto calcSolPtr = solver.calcBoundFacePrimVal[iter].data();
		int calcSolSize = solver.calcBoundFacePrimVal[iter].size();
		// auto calcAddSolPtr = solver.calcFaceAddiVal.data();
		// int calcAddSolSize = solver.calcFaceAddiVal.size();
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL;
			// for(auto& sol : solver.calcBoundFacePrimVal[iter]){
			for(int iSol=0; iSol<calcSolSize; ++iSol){
				calcSolPtr[iSol](fieldVar[id_t],face.x,face.y,face.z,
					cellVar[iL].data(),faceVar[i].data());
			}
			// for(auto& sol : solver.calcFaceAddiVal){
			// for(int iSol=0; iSol<calcAddSolSize; ++iSol){
				sol(faceVar[i].data());
			// }
		}
		++iter;
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		// // cout << var.faces[0].size() << endl;
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	
	// int primSize = controls.primVarNames.size();
	// auto calcBoundFacePrimVal_ptr = calcBoundFacePrimVal.data();

	// // int iter=0;
	// // for(auto& boundary : mesh.boundaries){
		// // if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			// // for(int i=str; i<end; ++i){
				// // auto cellVar_i = cellVar[faces[i].iL].data();
				// // auto faceVar_i = faceVar[i].data();
				// // for(int j=0; j<primSize; ++j){
					// // calcBoundFacePrimVal_ptr[iter*primSize + j](cellVar_i, faceVar_i);
					
				// // // string type = boundary.types[j];
				// // // int id_inp = controls.primVarIds[j];
				// // // string name = controls.primVarNames[j];
				// // // string left_name = "left ";
				// // // string right_name = "right ";
				// // // left_name += name;
				// // // right_name += name;
				// // // int id_L_out = controls.faceVar[left_name].id;
				// // // int id_R_out = controls.faceVar[right_name].id;
				
				// // // if(boundary.name=="supinlet"){
				// // // cout << boundary.name << " " << name << " " << faceVar_i[id_L_out] << " " << faceVar_i[id_R_out] << " " << endl;
				// // // }
				// // }
				
				
			// // }
			// // ++iter;
		// // }
	// // }
}


void MASCH_Solver::updateBoundaryFaceAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
		// cout << var.faces[0].size() << endl;
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int size_c = calcCellAddiVal.size();
	int size_f = calcFaceAddiVal.size();
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				auto faceVar_i = faceVar[i].data();
				for(int j=0; j<size_f; ++j){
					calcFaceAddiVal[j](faceVar_i);
				}
			}
		}
	}
}

void MASCH_Solver::initOldValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	
	int currSize = solver.saveCurrValueNames.size();
	int old1Size = solver.saveOld1ValueNames.size();
	int old2Size = solver.saveOld2ValueNames.size();
	vector<int> currIds;
	vector<int> old1Ids;
	vector<int> old2Ids;
	for(int i=0; i<currSize; ++i) currIds.push_back(controls.getId_cellVar(solver.saveCurrValueNames[i]));
	for(int i=0; i<old1Size; ++i) old1Ids.push_back(controls.getId_cellVar(solver.saveOld1ValueNames[i]));
	for(int i=0; i<old2Size; ++i) old2Ids.push_back(controls.getId_cellVar(solver.saveOld2ValueNames[i]));
	
	// auto p_saveCurrIds = saveCurrIds.data();
	// auto p_saveOld1Ids = saveOld1Ids.data();
	// auto p_saveOld2Ids = saveOld2Ids.data();
	
	// auto size_saveCurrIds = saveCurrIds.size();
	// auto size_saveOld1Ids = saveOld1Ids.size();
	// auto size_saveOld2Ids = saveOld2Ids.size();
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0; j<old1Size; ++j){
			cellVar_i[old1Ids[j]] = cellVar_i[currIds[j]];
		}
		for(int j=0; j<old2Size; ++j){
			cellVar_i[old2Ids[j]] = cellVar_i[currIds[j]];
		}
	}
	
}



void MASCH_Solver::updateOldValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	
	int currSize = solver.saveCurrValueNames.size();
	int old1Size = solver.saveOld1ValueNames.size();
	int old2Size = solver.saveOld2ValueNames.size();
	vector<int> currIds;
	vector<int> old1Ids;
	vector<int> old2Ids;
	for(int i=0; i<currSize; ++i) currIds.push_back(controls.getId_cellVar(solver.saveCurrValueNames[i]));
	for(int i=0; i<old1Size; ++i) old1Ids.push_back(controls.getId_cellVar(solver.saveOld1ValueNames[i]));
	for(int i=0; i<old2Size; ++i) old2Ids.push_back(controls.getId_cellVar(solver.saveOld2ValueNames[i]));
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0; j<old2Size; ++j){
			cellVar_i[old2Ids[j]] = cellVar_i[old1Ids[j]];
		}
		for(int j=0; j<old1Size; ++j){
			cellVar_i[old1Ids[j]] = cellVar_i[currIds[j]];
		}
	}
	
}




void MASCH_Solver::calcTempSteps(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& solver = (*this);
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto procRightCellsVar = var.procRightCells.data();
	
	auto sol = solver.calcTempStepCell[iSegEq];
	
	
	int id_dt = controls.getId_fieldVar("time-step");
	fieldVar[id_dt] = 1.e12;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto cellVar_i = cellVar[i].data();
		// for(auto& sol : solver.calcTempStepCell){
			sol(cellVar_i, fieldVar);
		// }
		

		// for(auto& lim_phi_id : controls.limiterNamesForUnst){
			// cellVar_i[lim_phi_id] = 1.0;
		// }
		
		
	}
	if(size>1){
		double tmp_fieldVar = var.fields[id_dt];
		double tmp_fieldVar_glo;
		MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		var.fields[id_dt] = tmp_fieldVar_glo;
	}
	
	
	// for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		// int iL = faces[i].iL;
		// int iR = faces[i].iR;
		// for(int iFunct=0; iFunct<tempFunctFace_size; ++iFunct){
			// tempFunctFace_ptr[iFunct](cellVar[iL].data(), cellVar[iR].data(), 
				// faceVar[i].data(), fieldVar);
		// // cout << fieldVar[id_dt] << endl;
		// }
	// }
	// for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		// auto& boundary = mesh.boundaries[i];
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// int iL = face.iL;
				
				// for(int iFunct=0; iFunct<tempFunctFace_size; ++iFunct){
					// tempFunctFace_ptr[iFunct](
							// cellVar[iL].data(), nullptr, faceVar[i].data(), fieldVar);
				// }
			// }
		// }
		// else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// int iL = face.iL;
				// int iR = i-str;
				// for(int iFunct=0; iFunct<tempFunctFace_size; ++iFunct){
					// tempFunctFace_ptr[iFunct](
						// cellVar[iL].data(), nullptr, faceVar[i].data(), fieldVar);
				// }
				
				// ++ip;
			// }
		// }
	// }
	
	// if(size>1){
		// double tmp_fieldVar = fieldVar[id_dt];
		// double tmp_fieldVar_glo;
		// MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		// fieldVar[id_dt] = tmp_fieldVar_glo;
	// }
	
	
}