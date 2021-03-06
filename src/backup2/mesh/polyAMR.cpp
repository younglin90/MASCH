#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>

#include "mesh.h"
#include "polyAMR.h"
#include "geometric.h" 
#include "../mpi/mpi.h"
#include "../solvers/solvers.h"  

void SEMO_Poly_AMR_Builder::polyAMR(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_Mesh_Geometric geometric;
	SEMO_Utility_Math math;
	SEMO_Solvers_Builder solvers;
	

	SEMO_MPI_Builder mpi;
	if(size>1){
		mpi.setCellDatasToFaceRight(mesh, 
					controls.VF[0], controls.fVF[0],
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	}
	// cout << "AAAAAAAA" << endl;
	vector<vector<double>> gradVF(mesh.cells.size(),vector<double>(3,0.0));
	// vector<double> dummy;
	// math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
		// controls.VF[0], controls.fVF[0], dummy, gradVF);
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	for(int i=0; i<mesh.cells.size(); ++i){
		mesh.cells[i].var[controls.indicatorAMR[0]] = 
			sqrt(gradVF[i][0]*gradVF[i][0]+
				 gradVF[i][1]*gradVF[i][1]+
				 gradVF[i][2]*gradVF[i][2]);
	}
	
	
	
	// ======================
	// add Buffer layer
	for(int iter=0; iter<controls.bufferLayer; ++iter){

		vector<double> newIndicatorAMR(mesh.cells.size(),0.0);
		vector<double> recvValues;
		if(size>1){
			// processor faces
			vector<double> sendValues;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					sendValues.push_back(mesh.cells[face.owner].var[controls.indicatorAMR[0]]);
				}
			}
			mpi.setProcsFaceDatas(
						sendValues, recvValues,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
		}
		
		int num_proc = 0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double maxInd = 0.0;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				maxInd = max(mesh.cells[face.owner].var[controls.indicatorAMR[0]],
							 mesh.cells[face.neighbour].var[controls.indicatorAMR[0]]);
				newIndicatorAMR[face.owner] = max(newIndicatorAMR[face.owner],maxInd);
				newIndicatorAMR[face.neighbour] = max(newIndicatorAMR[face.neighbour],maxInd);
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				maxInd = max(mesh.cells[face.owner].var[controls.indicatorAMR[0]],
							 recvValues[num_proc]);
				newIndicatorAMR[face.owner] = max(newIndicatorAMR[face.owner],maxInd);
				++num_proc;
			}
			else{
				maxInd = mesh.cells[face.owner].var[controls.indicatorAMR[0]];
				newIndicatorAMR[face.owner] = max(newIndicatorAMR[face.owner],maxInd);
			}
		}
		
		for(int i=0; i<mesh.cells.size(); ++i){
			mesh.cells[i].var[controls.indicatorAMR[0]] = newIndicatorAMR[i];
		}
	}
	// ======================
	
	
	
	// // ======================
	// double coeff_indi_refine = 5.0; // 2D
	// // double coeff_indi_refine = 100.0; // 3D
	// double coeff_indi_const = 1.0;
	// // ======================
	
	
	// for(auto& cell : mesh.cells){
		// if(cell.var[controls.indicatorAMR[0]] > coeff_indi_refine){
			// cell.var[controls.indicatorAMR[0]] = coeff_indi_const;  //coeff_mul_indi * controls.indicatorRefine;
		// }
		// else{
			// cell.var[controls.indicatorAMR[0]] = 0.0;
		// } 
	// }
	

	// vector<double> smoothAi;
	// vector<double> AiUp;
	// vector<double> AiDown;
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// smoothAi.push_back(mesh.cells[i].var[controls.indicatorAMR[0]]);
	// }
	
	// //================================================
	// // smoothing Ai
	// for(int iter=0; iter<5; ++iter){
		
		// AiUp.clear();
		// AiDown.clear();
		// for(int i=0; i<mesh.cells.size(); ++i){
			// AiUp.push_back(0.0);
			// AiDown.push_back(0.0);
		// }
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC;
			// // double wCL = 0.5;
			// double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
				// double AiF = 
					// wCL*smoothAi[face.owner] +
					// wCR*smoothAi[face.neighbour];
				
				// AiUp[face.owner] += AiF*face.area;
				// AiUp[face.neighbour] += AiF*face.area;
				
				// AiDown[face.owner] += face.area;
				// AiDown[face.neighbour] += face.area;
			
			// }
			
		// }

		// // boundary
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double AiF = smoothAi[face.owner];
					
					// AiUp[face.owner] += AiF*face.area;
				
					// AiDown[face.owner] += face.area;
					
				// }
			// }
		// }
		
		// if(size>1){
			// // processor faces
			// vector<double> sendValues;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// sendValues.push_back(smoothAi[face.owner]);
				// }
			// }
			// vector<double> recvValues;
			// mpi.setProcsFaceDatas(
						// sendValues, recvValues,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// int num=0;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC;
				// // double wCL = 0.5;
				// double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double AiF = wCL*smoothAi[face.owner]+wCR*recvValues[num];
					
					// AiUp[face.owner] += AiF*face.area;
				
					// AiDown[face.owner] += face.area;
					
					// ++num;
				// }
			// }
			
		// }
	
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// smoothAi[i] = AiUp[i]/AiDown[i];
			// if(mesh.cells[i].var[controls.indicatorAMR[0]] > coeff_indi_refine){
				// smoothAi[i] = coeff_indi_const; //coeff_mul_indi * controls.indicatorRefine;
			// }
			
		// }
		
	// }
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// mesh.cells[i].var[controls.indicatorAMR[0]] = smoothAi[i];
	// }
	
	
	
	
	// //=========================================================
	// for(auto& cell : mesh.cells){
		// if(cell.var[controls.indicatorAMR[0]] > controls.indicatorRefine){
			// cell.var[controls.indicatorAMR[0]] = 1.0;
		// }
		// else{
			// cell.var[controls.indicatorAMR[0]] = 0.0;
		// } 
	// }
	
	// for(int ipp=0; ipp<5; ++ipp){
		// vector<double> var_recv;
		// if(size>1){
			// vector<double> var_send;
			
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// var_send.push_back(mesh.cells[face.owner].var[controls.indicatorAMR[0]]);
				// }
			// }
			
			// var_recv.resize(var_send.size(),0);

			// MPI_Alltoallv( var_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
						   // var_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
						   // MPI_COMM_WORLD);
						   
		// }
		
		// vector<double> tmp_grad(mesh.cells.size(),0.0);
		// int proc_num=0;
		// for(auto& face : mesh.faces){
			// vector<double> distanceCells;
			// distanceCells.push_back(face.distCells[0]);
			// distanceCells.push_back(face.distCells[1]);
			// distanceCells.push_back(face.distCells[2]);
			// double dPN = face.magPN;
							  
			// double dtstepK = 0.6*0.5*dPN*dPN;//0.1*0.5*dPN*dPN*dPN*dPN;
			
			// if(face.getType()==SEMO_Types::INTERNAL_FACE){
				// double diff_flux = dtstepK*
					// (
					// mesh.cells[face.neighbour].var[controls.indicatorAMR[0]]-
					// mesh.cells[face.owner].var[controls.indicatorAMR[0]]
					// )
					// /dPN*face.area;
				// tmp_grad[face.owner] += diff_flux/mesh.cells[face.owner].volume;
				// tmp_grad[face.neighbour] -= diff_flux/mesh.cells[face.neighbour].volume;
			// }
			// else if(face.getType()==SEMO_Types::PROCESSOR_FACE){
				// double diff_flux = dtstepK*
					// (
					// var_recv[proc_num]-
					// mesh.cells[face.owner].var[controls.indicatorAMR[0]]
					// )
					// /dPN*face.area;
				// tmp_grad[face.owner] += diff_flux/mesh.cells[face.owner].volume;
				// ++proc_num;
			// }
		// }
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			// if(cell.var[controls.indicatorAMR[0]]<0.99999999999){
				// cell.var[controls.indicatorAMR[0]] += tmp_grad[i];
				// // if(cell.var[controls.indicatorAMR[0]] > 0.7){
					// // cell.var[controls.indicatorAMR[0]] = 1.0;
				// // }
			// }
			
		// }
	// }
	// //=========================================================
	
	// vector<double> maxGrad(controls.indicatorAMR.size(),-1.e10);
	// vector<double> minGrad(controls.indicatorAMR.size(),1.e10);
	// for(auto& cell : mesh.cells){
		// double tmp_var = cell.var[controls.indicatorAMR[0]];
		// maxGrad[0] = max(maxGrad[0],tmp_var);
		// minGrad[0] = min(minGrad[0],tmp_var);
	// }
	// for(auto& cell : mesh.cells){
		// cell.var[controls.indicatorAMR[0]] /= maxGrad[0];
			
		// // if(abs(maxGrad[0]-minGrad[0])<1.e-200){
			// // mesh.cells[i].var[controls.indicatorAMR[0]] = 0.0;
		// // }
		// // else{
			// // mesh.cells[i].var[controls.indicatorAMR[0]] =
				// // (tmp_var-minGrad[0])/(maxGrad[0]-minGrad[0]);
		// // }
	// }


	// for(int i=0; i<5; ++i){
		
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// mesh.cells[i].var[controls.indicatorAMR[0]] = 100000.0;
	// }
	
	// double strSolverClock = clock();
		
		
	if( (controls.iterReal+1) % controls.intervalRefine == 0){
		if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
		polyRefine(mesh, controls, 0);
	}
	
	// geometric.init(mesh);
	
	if( (controls.iterReal+1) % controls.intervalUnrefine == 0){
		if(rank==0) cout << "| exe. Poly AMR Un-refinement" << endl;
		polyUnrefine(mesh, controls, 0); 
	}
	
	
		// double dClock = clock() - strSolverClock;
		// dClock /= CLOCKS_PER_SEC;
		// cout << dClock << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
	// ???????????? ??? ??????
	mesh.cellsProcVar.resize(controls.nTotalCellVar,vector<double>());
	mesh.cellsProcGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	mesh.cellsGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	
	
	// geometric.init(mesh);
	
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "1111111" << endl;
	
	
	
	
	
	
	// mesh.informations();
	
	
	
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "22222" << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// SEMO_Mesh_Save save;
	// save.vtu("./save/1/", mesh, controls, species);
	// // save.vtu("./save/2/", mesh, controls, species);
	
	

	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	geometric.init(mesh);
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "33333" << endl;
	
	math.initLeastSquare(mesh); 
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << "44444" << endl;
	
	
	
	// solvers.calcIncomCellEOSVF(mesh, controls, species);
	// solvers.calcCellTransport(mesh, controls, species);
	
	

	// SEMO_Mesh_Save save;
	// // string tmpFile = "./Uf" + to_string(iter);
	// string tmpFile = "./";
	// save.vtu(tmpFile, rank, mesh);
	
	
	
	
}




void SEMO_Poly_AMR_Builder::mpiLevelRefine(
	SEMO_Mesh_Builder& mesh, 
	vector<bool>& boolCellRefine,
	vector<int>& cLevel_recv, 
	vector<int>& cRefine_recv){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size>1){
		vector<int> cLevel_send;
		vector<int> cRefine_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cLevel_send.push_back(mesh.cells[face.owner].level);
				
				if(boolCellRefine[face.owner]){
					cRefine_send.push_back(1);
				}
				else{
					cRefine_send.push_back(0);
				}
			}
		}
		
		cLevel_recv.resize(cLevel_send.size(),0);
		cRefine_recv.resize(cRefine_send.size(),0);

		MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
					   
		MPI_Alltoallv( cRefine_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cRefine_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
					   
	}
	



}




void SEMO_Poly_AMR_Builder::mpiRefines(
	SEMO_Mesh_Builder& mesh, 
	vector<bool>& boolCellRefine,
	vector<int>& cRefine_recv){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size>1){
		vector<int> cRefine_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				if(boolCellRefine[face.owner]){
					cRefine_send.push_back(1);
				}
				else{
					cRefine_send.push_back(0);
				}
			}
		}
		
		cRefine_recv.resize(cRefine_send.size(),0);
					   
		MPI_Alltoallv( cRefine_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cRefine_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
					   
	}
	



}

void SEMO_Poly_AMR_Builder::mpiLevels(
	SEMO_Mesh_Builder& mesh, 
	vector<int>& cLevel_recv){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size>1){
		vector<int> cLevel_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cLevel_send.push_back(mesh.cells[face.owner].level);
			}
		}
		
		// cout << "1!!!!!!" << endl;
		
		cLevel_recv.clear();
		cLevel_recv.resize(cLevel_send.size(),0);

		// cout << "2!!!!!!    : " << cLevel_send.size() << endl;
		MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
				
		// cout << "3!!!!!!" << endl;	   
	}
	
}




void SEMO_Poly_AMR_Builder::restrictCellRefine(
	SEMO_Mesh_Builder& mesh, 
	vector<bool>& boolCellRefine,
	vector<int>& cLevel_recv, 
	vector<int>& cRefine_recv){

	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			if(mesh.cells[face.owner].level > mesh.cells[face.neighbour].level){
				boolCellRefine[face.owner] = false;
			}
			if(mesh.cells[face.owner].level < mesh.cells[face.neighbour].level){
				boolCellRefine[face.neighbour] = false;
			}
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			if(mesh.cells[face.owner].level > cLevel_recv[proc_num]){
				boolCellRefine[face.owner] = false;
			}
			++proc_num;
		}
	}
	// cLevel_recv.clear();
	// cRefine_recv.clear();
}



void SEMO_Poly_AMR_Builder::createEdges(
	SEMO_Mesh_Builder& mesh, 
	vector<int>& edgesPoint0,
	vector<int>& edgesPoint1, 
	vector<vector<int>>& facesEdges,
	vector<vector<int>>& edgesFaces,
	vector<int>& edgeLevel){
		
		
	facesEdges.resize(mesh.faces.size(),vector<int>(0,0));
	vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		int pointSize = face.points.size();
		for(int j=0; j<pointSize; ++j){
			int ipoint0 = face.points[j];
			int ipoint1 = ( j+1 == pointSize ? face.points[0] : face.points[j+1] );
			vector<int> matchFaces;
			for(auto& k : pointsFaces[ipoint0]){
				for(auto& l : pointsFaces[ipoint1]){
					if( k == l ) {
						matchFaces.push_back(l);
					}
				}
			}
			
			if(matchFaces.size()==0){
				edgesPoint0.push_back(ipoint0);
				edgesPoint1.push_back(ipoint1);
				
				facesEdges[i].push_back(edgesPoint0.size()-1);
				
			}
			else{
				int iFace = matchFaces[0];
				int iEdgeSave = -1;
				for(auto& iEdge : facesEdges[iFace]){
					if(
					(edgesPoint0[iEdge]==ipoint0 && edgesPoint1[iEdge]==ipoint1) ||
					(edgesPoint1[iEdge]==ipoint0 && edgesPoint0[iEdge]==ipoint1) 
					){
						iEdgeSave = iEdge;
						break;
					}
				}
				
				facesEdges[i].push_back(iEdgeSave);
			}
			pointsFaces[ipoint0].push_back(i);
			
		}
	}
	pointsFaces.clear();
	
	edgesFaces.resize(edgesPoint0.size(),vector<int>(0,0));
	for(int i=0; i<mesh.faces.size(); ++i){
		for(auto& j : facesEdges[i]){
			edgesFaces[j].push_back(i);
		}
	}

	edgeLevel.resize(edgesPoint0.size(),0);
	for(int i=0; i<edgesPoint0.size(); ++i){
		int point0 = edgesPoint0[i];
		int point1 = edgesPoint1[i];
		
		edgeLevel[i] = max(
			mesh.points[point0].level, 
			mesh.points[point1].level);
	}
	
	
	
}



void SEMO_Poly_AMR_Builder::searchOriginalPoints(
	SEMO_Mesh_Builder& mesh, 
	vector<int>& points,
	int targetLevel, 
	vector<int>& originPoints){
		
	originPoints.clear();
	for(auto& i : points){
		if(mesh.points[i].level <= targetLevel){
			originPoints.push_back(i);
		}
	}
	
}