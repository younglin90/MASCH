
#include "./solvers.h"

void MASCH_Solver::fvm(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	
	controls.log.push("gradientTerms");
	solver.gradientTerms(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("highOrderTerms");
	solver.highOrderTerms(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("temporalTerms");
	solver.temporalTerms(mesh, controls, var);
	controls.log.pop();
	controls.log.push("convectiveTerms");
	solver.convectiveTerms(mesh, controls, var);
	controls.log.pop();
	controls.log.push("diffusionTerms");
	solver.diffusionTerms(mesh, controls, var);
	controls.log.pop();
	controls.log.push("sourceTerms");
	solver.sourceTerms(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("linearSystem");
	solver.linearSystem(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("updateCellPrimValues");
	solver.updateCellPrimValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateCellAddiValues");
	solver.updateCellAddiValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateProcRightCellPrimValues");
	solver.updateProcRightCellPrimValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateProcRightCellAddiValues");
	solver.updateProcRightCellAddiValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateBoundaryFacePrimValues");
	solver.updateBoundaryFacePrimValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateBoundaryFaceAddiValues");
	solver.updateBoundaryFaceAddiValues(mesh, controls, var);
	controls.log.pop();
	
	controls.log.push("inner updateOldValues");
	solver.updateOldValues(mesh, controls, var);
	controls.log.pop();
	
	
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
			tmp_counts[ip] = mesh.send_countsStencilCells[ip]*prim_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_displs[ip+1] = tmp_displs[ip] + tmp_counts[ip];
		}
		
		vector<double> recv_value(proc_size*prim_size);
		MPI_Alltoallv( send_value.data(), tmp_counts.data(), tmp_displs.data(), MPI_DOUBLE, recv_value.data(), tmp_counts.data(), tmp_displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		int iter=0;
		for(auto& cells : var.procRightCells){
			auto cellVar_i = cells.data();
			for(int j=0; j<prim_size; ++j){
				cellVar_i[ids[j]] = recv_value[iter++];
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


void MASCH_Solver::updateCellAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
		// cout << var.faces[0].size() << endl;
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto p_calcCellAddiVal = calcCellAddiVal.data();
	
	int size_c = calcCellAddiVal.size();
	int size_f = calcFaceAddiVal.size();
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0; j<size_c; ++j){
			p_calcCellAddiVal[j](cellVar_i);
		}
	}
	
	
}

void MASCH_Solver::updateFaceAddiValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
		// cout << var.faces[0].size() << endl;
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int size_c = calcCellAddiVal.size();
	int size_f = calcFaceAddiVal.size();
	
	
	for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		auto faceVar_i = faceVar[i].data();
		for(int j=0; j<size_f; ++j){
			calcFaceAddiVal[j](faceVar_i);
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
	
	int size_f = calcBoundFacePrimVal.size();

	for(auto& boundary : mesh.boundaries){
		if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				auto cellVar_i = cellVar[faces[i].iL].data();
				auto faceVar_i = faceVar[i].data();
				for(int j=0; j<size_f; ++j){
					calcBoundFacePrimVal[j](cellVar_i, faceVar_i);
				}
			}
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
	
	for(auto& item : calcHO_FaceVal){
		for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
			int iL = faces[i].iL;
			int iR = faces[i].iR;
			item(cellVar[iL].data(), cellVar[iR].data(), faceVar[i].data());
		}
		
		for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
			auto& boundary = mesh.boundaries[i];
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
				for(int i=str; i<end; ++i, ++ip){
					auto& face = faces[i];
					int iL = face.iL;
					item(cellVar[iL].data(), 
					procRightCellsVar[ip].data(), faceVar[i].data());
				}
			}
		}
	
	}
	
	
	// int cellSize = mesh.cells.size();
	// auto cells = mesh.cells.data();
	// int faceSize = mesh.faces.size();
	// auto faces = mesh.faces.data();
	// int reconSize = solver.reconstruction.size();
	// auto reconstruction = solver.reconstruction.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	// auto procRightCellsVar = var.procRightCells.data();
	
	// int nEq = controls.nEq;
	
	// for(int i=0, SIZE=mesh.endInternalFace; i<SIZE; ++i){
		// int iL = faces[i].iL;
		// int iR = faces[i].iR;
		// for(int iRec=0; iRec<reconSize; ++iRec){
			// reconstruction[iRec](cellVar[iL].data(), cellVar[iR].data(), faceVar[i].data());
		// }
	// }
	
	// for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		// auto& boundary = mesh.boundaries[i];
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		// if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			// for(int i=str; i<end; ++i, ++ip){
				// auto& face = faces[i];
				// int iL = face.iL;
				// for(int iRec=0; iRec<reconSize; ++iRec){
					// reconstruction[iRec](cellVar[iL].data(), procRightCellsVar[ip].data(), faceVar[i].data());
				// }
			// }
		// }
	// }
	
	
}



void MASCH_Solver::temporalTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto tempFunc = solver.calcTemporal.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int nEq = controls.nEq;
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		
		for(auto item : solver.calcTemporal){
			item(cellVar[i].data(), fluxA, fluxB);
		}
		
		for(int jEq=0; jEq<nEq; ++jEq){
			for(int iEq=0; iEq<nEq; ++iEq){
				var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
			}
			var.accumB( i, jEq, fluxB[jEq] );
		}
	}
	
}

void MASCH_Solver::convectiveTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto convFunc = solver.calcConvFlux.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto procRightCellsVar = var.procRightCells.data();
	
	int nEq = controls.nEq;
	
	// vector<double> vecflux(nEq,0.0);
	// auto flux = vecflux.data();
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	
	
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		auto& face = faces[i];
		int iL = face.iL;
		int iR = face.iR;
		
		convFunc[0](cellVar[iL].data(), cellVar[iR].data(), 
					faceVar[i].data(), fluxA, fluxB);
		
		if(checkImplicitConvFlux[0]==true){
			// A sparse matrix 에 넣기
		}
		
		for(int iEq=0; iEq<nEq; ++iEq){
			var.accumB( iL, iEq, +fluxB[iEq] );
			var.accumB( iR, iEq, -fluxB[iEq] );
		}
		
	}
	
	
	for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		auto& boundary = mesh.boundaries[i];
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				
				convFunc[0](cellVar[iL].data(), nullptr, 
							faceVar[i].data(), fluxA, fluxB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iL, iEq, +fluxB[iEq] );
				}
			}
		}
		else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				int iR = i-str;
				
				convFunc[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							faceVar[i].data(), fluxA, fluxB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iL, iEq, +fluxB[iEq] );
				}
				++ip;
			}
		}
	}
	
}




void MASCH_Solver::diffusionTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	
	auto& solver = (*this);
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto calcLaplFlux = solver.calcLaplFlux.data();
	auto calcNLaplFlux = solver.calcNLaplFlux.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto procRightCellsVar = var.procRightCells.data();
	int nEq = controls.nEq;
	
	// vector<double> vecflux(nEq,0.0);
	// auto flux = vecflux.data();
	double fluxLaplA[nEq*nEq];
	double fluxLaplB[nEq];
	double fluxNLaplA[nEq*nEq];
	double fluxNLaplB[nEq];
	
	
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		auto& face = faces[i];
		int iL = face.iL;
		int iR = face.iR;
		
		calcLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					faceVar[i].data(), fluxLaplA, fluxLaplB);
		calcNLaplFlux[0](cellVar[iL].data(), cellVar[iR].data(), 
					faceVar[i].data(), fluxNLaplA, fluxNLaplB);
		
		if(checkImplicitConvFlux[0]==true){
			// A sparse matrix 에 넣기
		}
		
		for(int iEq=0; iEq<nEq; ++iEq){
			var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
			var.accumB( iR, iEq, -fluxLaplB[iEq]-fluxNLaplB[iEq] );
		}
		
	}
	
	for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
		auto& boundary = mesh.boundaries[i];
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				
				calcLaplFlux[0](cellVar[iL].data(), nullptr, 
							faceVar[i].data(), fluxLaplA, fluxLaplB);
				calcNLaplFlux[0](cellVar[iL].data(), nullptr, 
							faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
				}
			}
		}
		else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				int iL = face.iL;
				int iR = i-str;
				
				calcLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							faceVar[i].data(), fluxLaplA, fluxLaplB);
				calcNLaplFlux[0](cellVar[iL].data(), procRightCellsVar[ip].data(), 
							faceVar[i].data(), fluxNLaplA, fluxNLaplB);
							
				if(checkImplicitConvFlux[0]==true){
					// A sparse matrix 에 넣기
				}
				
				for(int iEq=0; iEq<nEq; ++iEq){
					var.accumB( iL, iEq, +fluxLaplB[iEq]+fluxNLaplB[iEq] );
				}
				++ip;
			}
		}
	}
	
	
}



void MASCH_Solver::sourceTerms(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto calcTemporal = solver.calcTemporal.data();
	auto calcSource = solver.calcSource.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	int nEq = controls.nEq;
	double fluxA[nEq*nEq];
	double fluxB[nEq];
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		
		for(auto item : solver.calcSource){
			item(cellVar[i].data(), fluxA, fluxB);
		}
		
		for(int jEq=0; jEq<nEq; ++jEq){
			for(int iEq=0; iEq<nEq; ++iEq){
				var.accumSparD( i, iEq, jEq, fluxA[jEq*nEq+iEq] );
			}
			var.accumB( i, jEq, fluxB[jEq] );
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
	
	vector<double> vecAmat(nEq*nEq,0.0);
	auto Amat = vecAmat.data();
	
	// controls.log.pop();
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
	// controls.log.push("solve ddd1");
		auto& cell = cells[i];
		
		for(int iEq=0; iEq<nEq; ++iEq){
			for(int jEq=0; jEq<nEq; ++jEq){
				Amat[iEq*nEq+jEq] = var.getSparD(i, iEq, jEq);
			}
		}
	// controls.log.pop();
	// controls.log.push("solve ddd2");
		math.GaussSeidelSOR(Amat, nEq);
	// controls.log.pop();
	// controls.log.push("solve ddd3");
		for(int iEq=0; iEq<nEq; ++iEq){
			double tmp_value = 0.0;
			for(int jEq=0; jEq<nEq; ++jEq){
				tmp_value += Amat[iEq*nEq+jEq]*var.getB(i, jEq);
			}
			var.setX(i, iEq, tmp_value);
		}
	// controls.log.pop();
	}
}

void MASCH_Solver::updateCellPrimValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	int cellSize = mesh.cells.size();
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	int nEq = controls.nEq;
	// function<void*(int iEq, double& value)> limitPrim;
	// auto limitPrim = &MASCH_Control::limitPrim;
	// auto limitPrim = &controls.limitPrim;
	// void(*limitPrim)(int, double&) = &controls.limitPrim;
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();
		for(int iEq=0; iEq<nEq; ++iEq){
			cellVar_i[iEq] += var.getX(i, iEq);
			
			controls.limitPrim(iEq, cellVar_i[iEq]);
			// (*&limitPrim)(iEq, cellVar_i[iEq]);
		}
	}
}