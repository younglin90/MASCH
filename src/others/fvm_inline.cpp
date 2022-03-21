
#include "./solvers.h"

void MASCH_Solver::fvm_inline(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	double CFL = controls.maxCFL;
	double coeffCellCFL = 1.0;
	{
		// proc right cell 로 셀의 원시변수 넘기기
		controls.log.push("calc2");
		solver.updateProcRightCellPrimValues(mesh, controls, var);
		// solver.updateProcRightCellAddiValues(mesh, controls, var);
		controls.log.pop();
		
		
		// cell 추가적 변수
		auto cellVar = var.cells.data();
		controls.log.push("calc3");
		for(int i=0; i<mesh.cells.size(); ++i){
			auto cellVar_i = cellVar[i].data();
			for(auto& sol : solver.calcCellAddiVal){
				sol(cellVar_i);
			}
		}
		controls.log.pop();
		
		// proc right cell 추가적 변수
		controls.log.push("calc4");
		for(auto& cellVar_i : var.procRightCells){
			for(auto& sol : solver.calcCellAddiVal){
				sol(cellVar_i.data());
			}
		}
		controls.log.pop();
		
	}
	{
		controls.log.push("calc5");
		
		
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_c = controls.getId_cellVar("speed-of-sound");
		
		auto fieldVar = var.fields.data();
		var.fields[id_dt] = 1.e12;
		// for(int i=0, ip=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
			// auto& face = mesh.faces[i];
			// auto cellVar_iL = cellVar[face.iL].data();
			// if(face.getType()==MASCH_Face_Types::INTERNAL){
				// auto cellVar_iR = cellVar[face.iR].data();
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,cellVar_iR,faceVar[i].data(),fieldVar);
				// }
			// }
			// else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// auto cellVar_iR = var.procRightCells[ip].data();
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,cellVar_iR,faceVar[i].data(),fieldVar);
				// }
				// ++ip;
			// }
			// else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,nullptr,faceVar[i].data(),fieldVar);
				// }
			// }
		// }
		
		
		
		for(int i=0; i<mesh.cells.size(); ++i){
			double mag_vel = sqrt(
				var.cells[i][id_u]*var.cells[i][id_u]+
				var.cells[i][id_v]*var.cells[i][id_v]+
				var.cells[i][id_w]*var.cells[i][id_w]);
			double c = var.cells[i][id_c];
			double dt_tmp = pow(var.cells[i][id_vol],0.3)/(mag_vel+c);
			var.fields[id_dt] = min(var.fields[id_dt],dt_tmp);
		}
		controls.log.pop();
		controls.log.push("calc6");
		if(size>1){
			double tmp_fieldVar = var.fields[id_dt];
			double tmp_fieldVar_glo;
			MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			var.fields[id_dt] = tmp_fieldVar_glo;
		}
		// MPI_Barrier(MPI_COMM_WORLD);
		controls.log.pop();
		var.fields[id_dt] *= CFL;
		
		
		
		// int id_CFL = controls.getId_cellVar("Courant-Friedrichs-Lewy number");
		// double tmp_dt = coeffCellCFL*var.fields[id_dt];
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double mag_vel = sqrt(
				// var.cells[i][id_u]*var.cells[i][id_u]+
				// var.cells[i][id_v]*var.cells[i][id_v]+
				// var.cells[i][id_w]*var.cells[i][id_w]);
			// double c = var.cells[i][id_c];
			// double dt_cell = pow(var.cells[i][id_vol],0.3)/(mag_vel+c);
			// var.cells[i][id_CFL] = tmp_dt/dt_cell;
		// }
		
	}
	{
		// controls.log.pop();
		controls.log.push("calc8");
		
		// 바운더리 페이스 값들 정하기
		int id_t = controls.getId_fieldVar("time");
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		for(int ibc=0, iter=0, SIZE=mesh.boundaries.size(); ibc<SIZE; ++ibc){
			auto& boundary = mesh.boundaries[ibc];
			if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
			auto calcSolPtr = solver.calcBoundFacePrimVal[iter].data();
			int calcSolSize = solver.calcBoundFacePrimVal[iter].size();
			auto calcAddSolPtr = solver.calcFaceAddiVal.data();
			int calcAddSolSize = solver.calcFaceAddiVal.size();
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
				for(int iSol=0; iSol<calcAddSolSize; ++iSol){
					calcAddSolPtr[iSol](faceVar[i].data());
				}
			}
			++iter;
		}
		controls.log.pop();
	}
		
	{
		controls.log.push("calc9.2");
		// 그레디언트 구하기
		solver.gradientTerms(mesh, controls, var);
		// 그레디언트 proc risght cell로 넘기기
		solver.updateProcRightCellGradValues(mesh, controls, var);
		controls.log.pop();
	}
	{
		controls.log.push("calc9.3");
		// 하이오더 리컨스트럭션 + 추가적 변수들
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		auto calcAddSolPtr = solver.calcFaceAddiVal.data();
		int calcAddSolSize = solver.calcFaceAddiVal.size();
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType()==MASCH_Face_Types::BOUNDARY) continue;
			
			int iL = face.iL; int iR = face.iR;
			auto faceVar_i = faceVar[i].data();
			auto cellVar_iL = cellVar[iL].data();
			
			// 하이 오더 리컨스트럭션
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				for(auto& sol : solver.calcHO_FaceVal){
					sol(var.fields.data(),cellVar_iL,cellVar[iR].data(),faceVar_i);
				}
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				for(auto& sol : solver.calcHO_FaceVal){
					sol(var.fields.data(),cellVar_iL,procRightCellVar[ip++].data(),faceVar_i);
				}
			}
			
			// 추가적인 변수들
			for(int iSol=0; iSol<calcAddSolSize; ++iSol){
				calcAddSolPtr[iSol](faceVar[i].data());
			}
			
			
		}
		controls.log.pop();
	}
	
	{
		controls.log.push("calc10");
		// 리니어 시스템 초기화
		var.clearLinearSystems();
		controls.log.pop();
		
	}
	{
		controls.log.push("calc11");
		// 컨벡티브 + 디퓨젼 풀기
		int nEq = controls.nEq;
		double fluxA[nEq*nEq];
		double fluxB[nEq];
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL;
			int iR = face.iR;
			
			// 컨벡티브 텀
			for(auto& sol : solver.calcConvFlux){
				sol(cellVar[iL].data(), cellVar[iR].data(), 
					faceVar[i].data(), fluxA, fluxB);
			}
			
			if(checkImplicitConvFlux[0]==true){
				// A sparse matrix 에 넣기
			}
			
			for(int iEq=0; iEq<nEq; ++iEq){
				var.accumB( iL, iEq, +fluxB[iEq] );
				var.accumB( iR, iEq, -fluxB[iEq] );
				// var.Bvalues[nEq*iL+iEq] += (+fluxB[iEq]);
				// var.Bvalues[nEq*iR+iEq] += (-fluxB[iEq]);
			}
		}
		for(int i=0, ip=0, SIZE=mesh.boundaries.size(); i<SIZE; ++i){
			auto& boundary = mesh.boundaries[i];
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			// double* cellVar_iR;
			if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					
					// 컨벡티브 텀
					for(auto& sol : solver.calcConvFlux){
						sol(cellVar[iL].data(), nullptr, 
							faceVar[i].data(), fluxA, fluxB);
					}
								
					if(checkImplicitConvFlux[0]==true){
						// A sparse matrix 에 넣기
					}
					
					for(int iEq=0; iEq<nEq; ++iEq){
						var.accumB( iL, iEq, +fluxB[iEq] );
						// var.Bvalues[nEq*iL+iEq] += (+fluxB[iEq]);
					}
				}
			}
			else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					
					// 컨벡티브 텀
					for(auto& sol : solver.calcConvFlux){
						sol(cellVar[iL].data(), procRightCellVar[ip++].data(), 
							faceVar[i].data(), fluxA, fluxB);
					}
								
					if(checkImplicitConvFlux[0]==true){
						// A sparse matrix 에 넣기
					}
					
					for(int iEq=0; iEq<nEq; ++iEq){
						var.accumB( iL, iEq, +fluxB[iEq] );
						// var.Bvalues[nEq*iL+iEq] += (+fluxB[iEq]);
					}
				}
			}
		}
		controls.log.pop();
	}
	
	
	{
		controls.log.push("calc12");
		// 시간텀
		int nEq = controls.nEq;
		double fluxA[nEq*nEq];
		double fluxB[nEq];
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		for(int i=0; i<mesh.cells.size(); ++i){
			auto cellVar_i = cellVar[i].data();
			for(auto item : solver.calcTemporal){
				item(cellVar_i, fieldVar, fluxA, fluxB);
			}
			for(int iEq=0; iEq<nEq; ++iEq){
				for(int jEq=0; jEq<nEq; ++jEq){
					var.accumSparD( i, iEq, jEq, fluxA[iEq*nEq+jEq] );
				}
				var.accumB( i, iEq, fluxB[iEq] );
			}
		}
		controls.log.pop();
	}
		
		
	{
		controls.log.push("calc14");
		int nEq = controls.nEq;
		int nSp = controls.nSp;
		double fluxA[nEq*nEq];
		double fluxB[nEq];
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		int id_resi = controls.fieldVar["residual"].id;
		int id_vol = controls.getId_cellVar("volume");
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		int id_Y[2];
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
			id_Y[i] = controls.getId_cellVar("mass-fraction-"+tmp_name);
		}
		
		fieldVar[id_resi] = 0.0;
		double tmp_volume = 0.0;
		for(int i=0; i<mesh.cells.size(); ++i){
			double volume = var.cells[i][id_vol];
			auto cellVar_i = cellVar[i].data();
			
			vector<vector<double>> matA(nEq,vector<double>(nEq,0.0));
			
			for(int iEq=0; iEq<nEq; ++iEq){
				for(int jEq=0; jEq<nEq; ++jEq){
					matA[iEq][jEq] = var.getSparD(i, iEq, jEq);
				}
			}
			
			math.GaussSeidelSOR(matA);
			
			for(int iEq=0; iEq<nEq; ++iEq){
				double tmp_val = 0.0;
				for(int jEq=0; jEq<nEq; ++jEq){
					// tmp_val += matA[iEq][jEq]*var.Bvalues[nEq*i+jEq];
					// tmp_val += matA[iEq][jEq]*var.Bvalues[mesh.cells.size()*jEq+i];
					tmp_val += matA[iEq][jEq]*var.getB(i, jEq);
				}
				var.Xvalues[nEq*i+iEq] = tmp_val;
			}
			cellVar_i[id_p] += var.Xvalues[nEq*i+0];
			cellVar_i[id_u] += var.Xvalues[nEq*i+1];
			cellVar_i[id_v] += var.Xvalues[nEq*i+2];
			cellVar_i[id_w] += var.Xvalues[nEq*i+3];
			cellVar_i[id_T] += var.Xvalues[nEq*i+4];
			cellVar_i[id_Y[0]] += var.Xvalues[nEq*i+5];
			
			// if(cellVar_i[id_p]<1000.0) cellVar_i[id_p]=1000.0;
			// if(cellVar_i[id_T]<100.0) cellVar_i[id_T]=100.0;
			// if(cellVar_i[id_T]>800.0) cellVar_i[id_T]=800.0;
			
			cellVar_i[id_Y[0]] = max(0.0,min(1.0,cellVar_i[id_Y[0]]));
			cellVar_i[id_Y[1]] = 1.0-cellVar_i[id_Y[0]];
			
			for(int jEq=0; jEq<nEq; ++jEq){
				double tmp_resi = var.Xvalues[nEq*i+jEq]*volume;
				fieldVar[id_resi] += tmp_resi*tmp_resi;
				tmp_volume += volume;
			}
		}
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
		
		var.fields[controls.fieldVar["time"].id] +=
		var.fields[controls.fieldVar["time-step"].id];
		
		controls.log.pop();
	}
	
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
}