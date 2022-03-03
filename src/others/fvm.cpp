
#include "./solvers.h"

void MASCH_Solver::fvm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	
	// 시간 업데이트
	var.fields[controls.fieldVar["time"].id] += 
	var.fields[controls.fieldVar["time-step"].id];
	
	// controls.log.push("linear solver ZERO");
	fill(var.tmp_sparA.begin(),var.tmp_sparA.end(),0.0);
	fill(var.tmp_B.begin(),var.tmp_B.end(),0.0);
	fill(var.tmp_X.begin(),var.tmp_X.end(),0.0);
	// controls.log.pop();
	
	// 셀 루프 텀, temporal 텀 + 소스텀
	// controls.log.push("temporal & source Terms");
	solver.cellLoopTerms(mesh, controls, var);
	// controls.log.pop();
	// 페이스 루프 텀, 컨벡티브 텀 + 디퓨젼 텀
	// controls.log.push("convective & diffusion Terms");
	solver.faceLoopTerms(mesh, controls, var);
	// controls.log.pop();
	
	// 선형 솔버
	// controls.log.push("linearSystem");
	solver.linearSystem(mesh, controls, var);
	// controls.log.pop();
	
	// 셀 원시변수 업데이트
	// controls.log.push("updateCellPrimValues");
	solver.updateCellPrimValues(mesh, controls, var);
	// controls.log.pop();
	
	// 선형 솔버 값들 0.0 으로 만들어주기
	// controls.log.push("clearLinearSystems");
	var.clearLinearSystems();
	// controls.log.pop();
	
	// 타임스텝 구하기
	// controls.log.push("calcTempSteps");
	solver.calcTempSteps(mesh, controls, var);
	// controls.log.pop();
	
	// 셀 그레디언트
	// controls.log.push("gradientTerms");
	solver.gradientTerms(mesh, controls, var);
	// controls.log.pop();
	
	// 고차 reconstruction
	// controls.log.push("highOrderTerms");
	solver.highOrderTerms(mesh, controls, var);
	// controls.log.pop();
	
	// 원시변수 제외한 나머지 셀값 업데이트
	// controls.log.push("updateCellAddiValues");
	solver.updateCellAddiValues(mesh, controls, var);
	// controls.log.pop();
	
	// proc 셀의 원시변수 mpi 넘기기
	// controls.log.push("updateProcRightCellPrimValues");
	solver.updateProcRightCellPrimValues(mesh, controls, var);
	// controls.log.pop();
	
	// 셀 원시변수 제외한 나머지 proc 셀값 업데이트
	// controls.log.push("updateProcRightCellAddiValues");
	solver.updateProcRightCellAddiValues(mesh, controls, var);
	// controls.log.pop();
	
	// B.C. 원시변수 업데이트
	// controls.log.push("updateBoundaryFacePrimValues");
	solver.updateBoundaryFacePrimValues(mesh, controls, var);
	// controls.log.pop();
	
	// B.C. 원시변수 제외한 나머지 값 업데이트
	// controls.log.push("updateBoundaryFaceAddiValues");
	solver.updateBoundaryFaceAddiValues(mesh, controls, var);
	// controls.log.pop();
	
	// old 값 업데이트
	// controls.log.push("updateOldValues");
	solver.updateOldValues(mesh, controls, var);
	// controls.log.pop();
	
	
}


void MASCH_Solver::gradientTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	solver.calcGradient.leastSquare(mesh, controls, var, solver.gradLSIds_name);
	
}



void MASCH_Solver::highOrderTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto procRightCellsVar = var.procRightCells.data();
	
	int funct_size = calcHO_FaceVal.size();
	auto funct_ptr = calcHO_FaceVal.data();
	
	int functAddi_size = calcFaceAddiVal.size();
	auto functAddi_ptr = calcFaceAddiVal.data();
	
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		int iL = faces[i].iL;
		int iR = faces[i].iR;
		for(int iFunct=0; iFunct<funct_size; ++iFunct){
			funct_ptr[iFunct](cellVar[iL].data(), cellVar[iR].data(), faceVar[i].data());
		}
		for(int j=0; j<functAddi_size; ++j){
			functAddi_ptr[j](faceVar[i].data());
		}
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		auto& boundary = mesh.boundaries[i];
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=str; i<end; ++i, ++ip){
				auto& face = faces[i];
				int iL = face.iL;
				for(int iFunct=0; iFunct<funct_size; ++iFunct){
					funct_ptr[iFunct](cellVar[iL].data(), 
							procRightCellsVar[ip].data(), faceVar[i].data());
				}
				for(int j=0; j<functAddi_size; ++j){
					functAddi_ptr[j](faceVar[i].data());
				}
			}
		}
	}
	
}


void MASCH_Solver::cellLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto calcTemporal = solver.calcTemporal.data();
	auto calcSource = solver.calcSource.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	
	int nEq = controls.nEq;
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		
		// 시간텀
		for(auto item : solver.calcTemporal){
			item(cellVar[i].data(), fieldVar, fluxA, fluxB);
		}
		for(int jEq=0; jEq<nEq; ++jEq){
			for(int iEq=0; iEq<nEq; ++iEq){
				// var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
				var.tmp_sparA[nEq*nEq*i+nEq*jEq+iEq] += fluxA[jEq*nEq+iEq];
			}
			// var.accumB( i, jEq, fluxB[jEq] );
			var.tmp_B[nEq*i+jEq] += fluxB[jEq];
		}
		
		// // 소스텀
		// for(auto item : solver.calcSource){
			// item(cellVar[i].data(), fluxA, fluxB);
		// }
		// for(int jEq=0; jEq<nEq; ++jEq){
			// for(int iEq=0; iEq<nEq; ++iEq){
				// // var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
				// var.tmp_sparA[nEq*nEq*i+nEq*jEq+iEq] += fluxA[jEq*nEq+iEq];
			// }
			// // var.accumB( i, jEq, fluxB[jEq] );
			// var.tmp_B[nEq*i+jEq] += fluxB[jEq];
		// }
	}
}




void MASCH_Solver::faceLoopTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	
	auto& solver = (*this);
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto procRightCellsVar = var.procRightCells.data();
	int nEq = controls.nEq;
	
	auto convFunc = solver.calcConvFlux.data();
	auto calcLaplFlux = solver.calcLaplFlux.data();
	auto calcNLaplFlux = solver.calcNLaplFlux.data();
	
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	double fluxLaplA[nEq*nEq];
	double fluxLaplB[nEq];
	double fluxNLaplA[nEq*nEq];
	double fluxNLaplB[nEq];
	
	
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		auto& face = faces[i];
		int iL = face.iL;
		int iR = face.iR;
		
		// 컨벡티브 텀
		convFunc[0](cellVar[iL].data(), cellVar[iR].data(), 
					faceVar[i].data(), fluxA, fluxB);
		
		if(checkImplicitConvFlux[0]==true){
			// A sparse matrix 에 넣기
		}
		
		for(int iEq=0; iEq<nEq; ++iEq){
			// var.accumB( iL, iEq, +fluxB[iEq] );
			// var.accumB( iR, iEq, -fluxB[iEq] );
			var.tmp_B[nEq*iL+iEq] += (+fluxB[iEq]);
			var.tmp_B[nEq*iR+iEq] += (-fluxB[iEq]);
		}
		
		// // 디퓨젼 텀
		// calcLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					// faceVar[i].data(), fluxLaplA, fluxLaplB);
		// calcNLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					// faceVar[i].data(), fluxNLaplA, fluxNLaplB);
		
		// if(checkImplicitConvFlux[0]==true){
			// // A sparse matrix 에 넣기
		// }
		
		// for(int iEq=0; iEq<nEq; ++iEq){
			// // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
			// // var.accumB( iR, iEq, -fluxLaplB[iEq]-fluxNLaplB[iEq] );
			// var.tmp_B[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
			// var.tmp_B[nEq*iR+iEq] += (-fluxLaplB[iEq]-fluxNLaplB[iEq]);
		// }
		
		
	}
	
	
	for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		auto& boundary = mesh.boundaries[i];
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				
				// 컨벡티브 텀
				convFunc[0](cellVar[iL].data(), nullptr, 
							faceVar[i].data(), fluxA, fluxB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				
				double tmp_value=0.0;
				for(int iEq=0; iEq<nEq; ++iEq){
					// var.accumB( iL, iEq, +fluxB[iEq] );
					var.tmp_B[nEq*iL+iEq] += (+fluxB[iEq]);
					tmp_value+=( +fluxB[iEq]);
				}
				
			// double test1 = tmp_value;
			// if(isnan(test1) || test1<-1.e12 || test1>1.e12){
				// cout << test1 << endl;
				// cout << fluxB[0] << endl;
				// cout << fluxB[1] << endl;
				// cout << fluxB[2] << endl;
				// cout << fluxB[3] << endl;
				// cout << fluxB[4] << endl;
			// }
				// // 디퓨젼 텀
				// calcLaplFlux[0](cellVar[iL].data(), nullptr, 
							// faceVar[i].data(), fluxLaplA, fluxLaplB);
				// calcNLaplFlux[0](cellVar[iL].data(), nullptr, 
							// faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				// if(checkImplicitConvFlux[0]==true){
					// // A sparse matrix 에 넣기
				// }
				
				// for(int iEq=0; iEq<nEq; ++iEq){
					// // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
					// var.tmp_B[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
				// }
				
				
			}
		}
		else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				int iR = i-str;
				
				// 컨벡티브 텀
				convFunc[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							faceVar[i].data(), fluxA, fluxB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				for(int iEq=0; iEq<nEq; ++iEq){
					// var.accumB( iL, iEq, +fluxB[iEq] );
					var.tmp_B[nEq*iL+iEq] += ( +fluxB[iEq]);
				}
				
				// // 디퓨젼 텀
				// calcLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							// faceVar[i].data(), fluxLaplA, fluxLaplB);
				// calcNLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							// faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				// if(checkImplicitConvFlux[0]==true){
					// // A sparse matrix 에 넣기
				// }
				
				// for(int iEq=0; iEq<nEq; ++iEq){
					// // var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
					// var.tmp_B[nEq*iL+iEq] += (+fluxLaplB[iEq]+fluxNLaplB[iEq]);
				// }
				
				++ip;
			}
		}
	}
	
	
}



void MASCH_Solver::linearSystem(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	// controls.log.push("solve ddd0");
	auto& solver = (*this);
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int nEq = controls.nEq;
	
	vector<vector<double>> vec2Amat(nEq,vector<double>(nEq,0.0));
	auto vec2Amat_ptr = vec2Amat.data();
	
	vector<double> vecAmat(nEq*nEq,0.0);
	auto Amat = vecAmat.data();
	vector<double> vectestAmat(nEq*nEq,0.0);
	auto testAmat = vectestAmat.data();
		int id_w = controls.cellVar["z-velocity"].id;
	
	// controls.log.pop();
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
	// controls.log.push("solve ddd1");
		auto& cell = cells[i];
		
		vector<double> A(nEq*nEq), X(nEq), B(nEq);
		
		for(int jEq=0; jEq<nEq; ++jEq){
			for(int iEq=0; iEq<nEq; ++iEq){
				A[nEq*jEq+iEq] = var.tmp_sparA[nEq*nEq*i+nEq*jEq+iEq];
				// vec2Amat[jEq][iEq] = var.tmp_sparA[nEq*nEq*i+nEq*jEq+iEq];
			}
			// B[jEq] = var.tmp_B[nEq*i+jEq];
		}
		
		math.GaussSeidel(A.data(), nEq);
		// math.GaussSeidelSOR(vec2Amat);
		// math.GaussSeidelSOR(vec2Amat);
		// math.GaussSeidel(A.data(), B.data(), X.data(), nEq);
		
		for(int jEq=0; jEq<nEq; ++jEq){
			double tmp_val = 0.0;
			for(int iEq=0; iEq<nEq; ++iEq){
				tmp_val += A[nEq*jEq+iEq]*var.tmp_B[nEq*i+iEq];
				// tmp_val += vec2Amat[jEq][iEq]*var.tmp_B[nEq*i+iEq];
			}
			var.tmp_X[nEq*i+jEq] = tmp_val;
		}
		
	}
}

void MASCH_Solver::updateCellPrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldsVar = var.fields.data();
	int nEq = controls.nEq;
	auto primVarIds = controls.primVarIds.data();
	// function<void*(int iEq, double& value)> limitPrim;
	// auto limitPrim = &MASCH_Control::limitPrim;
	// auto limitPrim = &controls.limitPrim;
	// void(*limitPrim)(int, double&) = &controls.limitPrim;
	
	int id_resi = controls.fieldVar["residual"].id;
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();
		for(int iEq=0; iEq<nEq; ++iEq){
			// cellVar_i[iEq] += var.getX(i, iEq);
			int id = primVarIds[iEq];
			cellVar_i[id] += var.tmp_X[nEq*i+iEq];
			
			fieldsVar[id_resi] += var.tmp_X[nEq*i+iEq]*var.tmp_X[nEq*i+iEq];
			
			// if(iEq==3 && abs(cellVar_i[id])>1.e-16) {
				// cout << cellVar_i[id] << " " << var.tmp_X[nEq*i+iEq] << endl;
			// }		
			
			// controls.limitPrim(iEq, cellVar_i[iEq]);
			// (*&limitPrim)(iEq, cellVar_i[iEq]);
		}
	}
	
	if(size>1){
		double tmp_fieldVar = fieldsVar[id_resi];
		double tmp_fieldVar_glo;
		MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		fieldsVar[id_resi] = tmp_fieldVar_glo;
	}
	
	
}


void MASCH_Solver::updateCellAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
		// cout << var.faces[0].size() << endl;
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto p_calcCellAddiVal = calcCellAddiVal.data();
	
	int size_c = calcCellAddiVal.size();
	int size_f = calcFaceAddiVal.size();
	
		// int id_test = controls.cellVar["speed of sound"].id;
		
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0; j<size_c; ++j){
			p_calcCellAddiVal[j](cellVar_i);
		}
		// cout << cellVar_i[id_test] << endl;
	}
	
	
}


void MASCH_Solver::updateProcRightCellPrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	if(size>1){
		int proc_size = var.procRightCells.size();
		int prim_size = controls.primVarNames.size();
		vector<unsigned short> ids = controls.primVarIds;
		
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
		
		vector<int> tmp_counts(size,0);
		vector<int> tmp_displs(size+1,0);
		for(int ip=0; ip<size; ++ip){
			tmp_counts[ip] = mesh.countsProcFaces[ip]*prim_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_displs[ip+1] = tmp_displs[ip] + tmp_counts[ip];
		}
		
		vector<double> recv_value(proc_size*prim_size);
		MPI_Alltoallv( send_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
						recv_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
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
void MASCH_Solver::updateProcRightCellAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto p_calcCellAddiVal = calcCellAddiVal.data();
	
	int size_c = calcCellAddiVal.size();
	int size_f = calcFaceAddiVal.size();
	
	for(auto& cells : var.procRightCells){
		auto cellVar_i = cells.data();
		for(int j=0; j<size_c; ++j){
			p_calcCellAddiVal[j](cellVar_i);
		}
	}
	
	
}

void MASCH_Solver::updateBoundaryFacePrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
		// cout << var.faces[0].size() << endl;
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int primSize = controls.primVarNames.size();
	auto calcBoundFacePrimVal_ptr = calcBoundFacePrimVal.data();

	int iter=0;
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				auto cellVar_i = cellVar[faces[i].iL].data();
				auto faceVar_i = faceVar[i].data();
				for(int j=0; j<primSize; ++j){
					calcBoundFacePrimVal_ptr[iter*primSize + j](cellVar_i, faceVar_i);
					
				// string type = boundary.types[j];
				// int id_inp = controls.primVarIds[j];
				// string name = controls.primVarNames[j];
				// string left_name = "left ";
				// string right_name = "right ";
				// left_name += name;
				// right_name += name;
				// int id_L_out = controls.faceVar[left_name].id;
				// int id_R_out = controls.faceVar[right_name].id;
				
				// if(boundary.name=="supinlet"){
				// cout << boundary.name << " " << name << " " << faceVar_i[id_L_out] << " " << faceVar_i[id_R_out] << " " << endl;
				// }
				}
				
				
			}
			++iter;
		}
	}
}


void MASCH_Solver::updateBoundaryFaceAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
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

void MASCH_Solver::initOldValues(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	
	auto p_saveCurrIds = saveCurrIds.data();
	auto p_saveOld1Ids = saveOld1Ids.data();
	auto p_saveOld2Ids = saveOld2Ids.data();
	
	// cout << cellVar[0].size() << endl;
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0, jSIZE=saveCurrIds.size(); j<jSIZE; ++j){
			int curr_id = p_saveCurrIds[j];
			int old1_id = p_saveOld1Ids[j];
			int old2_id = p_saveOld2Ids[j];
	// cout << curr_id << " " << old1_id << " " << old2_id << endl;
			cellVar_i[old1_id] = cellVar_i[curr_id];
			cellVar_i[old2_id] = cellVar_i[curr_id];
		}
	}
	
	
	
}