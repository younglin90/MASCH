
#include "./mesh.h"
#include "./polyAMR.h"
// #include "geometric.h" 
#include "./mpi.h"
#include "./solvers.h" 
#include "./save.h"  

void MASCH_Poly_AMR_Builder::polyAMR(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	MASCH_Variables& var,
	int iter){

	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// // SEMO_Mesh_Geometric geometric;
	// // SEMO_Utility_Math math;
	// // SEMO_Solvers_Builder solvers;
	
	// vector<vector<double>> indicatorValues;
	// // for(int i=0; i<5; ++i)
	// // {
		// // if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
		
		// // polyRefine(mesh, controls, indicatorValues, 0);
		
		// // controls.setVariableArray(mesh, var);
		// // controls.setGeometric(mesh, var);
		
		// // polyUnrefine(mesh, controls, indicatorValues, 0);
		
		// // controls.setVariableArray(mesh, var);
		// // controls.setGeometric(mesh, var);
		
		// // // // random
		// // // std::random_device rd;
		// // // std::default_random_engine eng(rd());
		// // // std::uniform_int_distribution<int> distr(0, size-1);
		// // // vector<int> cell_ip(mesh.cells.size());
		// // // for(int i=0; i<mesh.cells.size(); ++i){
			// // // cell_ip[i] = distr(eng);
		// // // }
		
		// // // mesh.informations();
		
		// // vector<int> cell_ip(mesh.cells.size(),0);
		// // mesh.repartParMETIS(size, cell_ip, mesh);
		// // mesh.repartitioning(cell_ip);
		
		// // controls.setVariableArray(mesh, var);
		// // controls.setGeometric(mesh, var);
	// // }
	
	// for(int ii=0; ii<39; ++ii)
	// {
		
		// if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
		// polyRefine(mesh, controls, indicatorValues, 0);
		// controls.setVariableArray(mesh, var);
		// controls.setGeometric(mesh, var);
		// mesh.debug_procFace_unitNomals(0.8);
		
		// int iterReal = ii+1;
		// int calcInter = 8;
		// if(iterReal % calcInter == 0){
			// if(rank==0) cout << "| exe. Dynamic Load Balancing" << endl;
			// vector<int> cell_ip(mesh.cells.size(),rank);
			// mesh.repartParMETIS(size, cell_ip, mesh);
			// mesh.repartitioning(cell_ip);
			// controls.setVariableArray(mesh, var);
			// controls.setGeometric(mesh, var);
			// mesh.debug_procFace_unitNomals(0.8);
		// }
		
		// string tttt = "calc";
		// tttt += to_string(ii);
		// controls.log.push(tttt);
		
		// vector<int> aaaa;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// aaaa.push_back(i);
			// // for(int j=0; j<mesh.cells.size(); ++j){
			// // }
		// }
		// // sort(aaaa.begin(),aaaa.end());
		// controls.log.pop();
		// controls.log.show();
		
		
		
		// if(rank==0) cout << "| exe. Poly AMR Unrefinement" << endl;
		// polyUnrefine(mesh, controls, indicatorValues, 0);
		// controls.setVariableArray(mesh, var);
		// controls.setGeometric(mesh, var);
		// mesh.debug_procFace_unitNomals(0.8);
		
		// // if(rank==0) cout << "| exe. Poly AMR Unrefinement" << endl;
		// // polyUnrefine(mesh, controls, indicatorValues, 0);
		// // controls.setVariableArray(mesh, var);
		// // controls.setGeometric(mesh, var);
		// // mesh.debug_procFace_unitNomals(0.8);
		
		// {
			// int ncells = mesh.cells.size();
			// vector<int> recv_ncells(size);
			// MPI_Allgather(&ncells, 1, MPI_INT, recv_ncells.data(), 1, MPI_INT, MPI_COMM_WORLD);
			
			// if(rank==0) {
				// for(auto& item : recv_ncells){
					// cout << " | " << item;
				// }
				// cout << endl;
			// }
		// }
		
		
		// if(iterReal % calcInter == 0){
			// string savestring = "./save/";
			// savestring += tttt;
			// MASCH_Mesh_Save save;
			// save.fvmFiles(savestring, rank, mesh, controls, var);
		// }
		
		
	// }
	
	
	// // polyRefine(mesh, controls, indicatorValues, 0);
		
	// // for(int i=0; i<2; ++i)
	// // {
		// // polyRefine(mesh, controls, indicatorValues, 0);
		
		
		// // // random
		// // std::random_device rd;
		// // std::default_random_engine eng(rd());
		// // std::uniform_int_distribution<int> distr(0, size-1);
		// // vector<int> idBlockCell(mesh.cells.size());
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // idBlockCell[i] = distr(eng);
		// // }
		// // mesh.repartitioning(idBlockCell);
	// // }
	

	// // SEMO_MPI_Builder mpi;
	// // if(size>1){
		// // mpi.setCellDatasToFaceRight(mesh, 
					// // controls.VF[0], controls.fVF[0],
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
	// // }
	// // vector<vector<double>> gradVF(mesh.cells.size(),vector<double>(3,0.0));
	// // math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // mesh.cells[i].var[controls.indicatorAMR[0]] = 
			// // sqrt(gradVF[i][0]*gradVF[i][0]+
				 // // gradVF[i][1]*gradVF[i][1]+
				 // // gradVF[i][2]*gradVF[i][2]);
	// // }
	
	
	// // // ======================
	// // // add Buffer layer
	// // for(int iter=0; iter<controls.bufferLayer; ++iter){

		// // vector<double> newIndicatorAMR(mesh.cells.size(),0.0);
		// // vector<double> recvValues;
		// // if(size>1){
			// // // processor faces
			// // vector<double> sendValues;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == MASCH_Face_Types::PROCESSOR){
					// // sendValues.push_back(mesh.cells[face.iL].var[controls.indicatorAMR[0]]);
				// // }
			// // }
			// // mpi.setProcsFaceDatas(
						// // sendValues, recvValues,
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // }
		
		// // int num_proc = 0;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // double maxInd = 0.0;
			// // if(face.getType() == MASCH_Face_Types::INTERNAL){
				// // maxInd = max(mesh.cells[face.iL].var[controls.indicatorAMR[0]],
							 // // mesh.cells[face.iR].var[controls.indicatorAMR[0]]);
				// // newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
				// // newIndicatorAMR[face.iR] = max(newIndicatorAMR[face.iR],maxInd);
			// // }
			// // else if(face.getType() == MASCH_Face_Types::PROCESSOR){
				// // maxInd = max(mesh.cells[face.iL].var[controls.indicatorAMR[0]],
							 // // recvValues[num_proc]);
				// // newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
				// // ++num_proc;
			// // }
			// // else{
				// // maxInd = mesh.cells[face.iL].var[controls.indicatorAMR[0]];
				// // newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
			// // }
		// // }
		
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // mesh.cells[i].var[controls.indicatorAMR[0]] = newIndicatorAMR[i];
		// // }
	// // }
	// // // ======================
	
	
	// // if( (controls.iterReal+1) % controls.intervalRefine == 0){
		// // if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
		// // polyRefine(mesh, controls, 0);
	// // }
	
	// // if( (controls.iterReal+1) % controls.intervalUnrefine == 0){
		// // if(rank==0) cout << "| exe. Poly AMR Un-refinement" << endl;
		// // polyUnrefine(mesh, controls, 0); 
	// // }
	
	// // // 추가적인 셀 값들
	// // mesh.cellsProcVar.resize(controls.nTotalCellVar,vector<double>());
	// // mesh.cellsProcGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	// // mesh.cellsGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	
	// // geometric.init(mesh);
	// // math.initLeastSquare(mesh); 
	
}




void MASCH_Poly_AMR_Builder::mpiLevelRefine(
	MASCH_Mesh& mesh, 
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
			
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				cLevel_send.push_back(mesh.cells[face.iL].level);
				
				if(boolCellRefine[face.iL]){
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




void MASCH_Poly_AMR_Builder::mpiRefines(
	MASCH_Mesh& mesh, 
	vector<bool>& boolCellRefine,
	vector<int>& cRefine_recv){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size>1){
		vector<int> cRefine_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				
				if(boolCellRefine[face.iL]){
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

void MASCH_Poly_AMR_Builder::mpiLevels(
	MASCH_Mesh& mesh, 
	vector<int>& cLevel_recv){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size>1){
		vector<int> cLevel_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				cLevel_send.push_back(mesh.cells[face.iL].level);
			}
		}
		
		cLevel_recv.clear();
		cLevel_recv.resize(cLevel_send.size(),0);

		MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
				   
	}
	
}




void MASCH_Poly_AMR_Builder::restrictCellRefine(
	MASCH_Mesh& mesh, 
	vector<bool>& boolCellRefine,
	vector<int>& cLevel_recv, 
	vector<int>& cRefine_recv){

	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			if(mesh.cells[face.iL].level > mesh.cells[face.iR].level){
				boolCellRefine[face.iL] = false;
			}
			if(mesh.cells[face.iL].level < mesh.cells[face.iR].level){
				boolCellRefine[face.iR] = false;
			}
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			if(mesh.cells[face.iL].level > cLevel_recv[proc_num]){
				boolCellRefine[face.iL] = false;
			}
			++proc_num;
		}
	}
}



void MASCH_Poly_AMR_Builder::createEdges(
	MASCH_Mesh& mesh, 
	vector<int>& edgesPoint0,
	vector<int>& edgesPoint1, 
	vector<vector<int>>& facesEdges,
	vector<vector<int>>& edgesFaces,
	vector<int>& edgeLevel){
		
		
	// facesEdges.resize(mesh.faces.size(),vector<int>(0,0));
	// vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// int pointSize = face.ipoints.size();
		// for(int j=0; j<pointSize; ++j){
			// int ipoint0 = face.ipoints[j];
			// int ipoint1 = ( j+1 == pointSize ? face.ipoints[0] : face.ipoints[j+1] );
			// vector<int> matchFaces;
			// for(auto& k : pointsFaces[ipoint0]){
				// for(auto& l : pointsFaces[ipoint1]){
					// if( k == l ) {
						// matchFaces.push_back(l);
					// }
				// }
			// }
			
			// if(matchFaces.size()==0){
				// edgesPoint0.push_back(ipoint0);
				// edgesPoint1.push_back(ipoint1);
				
				// facesEdges[i].push_back(edgesPoint0.size()-1);
				
			// }
			// else{
				// int iFace = matchFaces[0];
				// int iEdgeSave = -1;
				// for(auto& iEdge : facesEdges[iFace]){
					// if(
					// (edgesPoint0[iEdge]==ipoint0 && edgesPoint1[iEdge]==ipoint1) ||
					// (edgesPoint1[iEdge]==ipoint0 && edgesPoint0[iEdge]==ipoint1) 
					// ){
						// iEdgeSave = iEdge;
						// break;
					// }
				// }
				
				// facesEdges[i].push_back(iEdgeSave);
			// }
			// pointsFaces[ipoint0].push_back(i);
			
		// }
	// }
	// pointsFaces.clear();
	
	
	
	
	
	
	

	
	vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	// vector<vector<int>> facesEdges(mesh.faces.size());
	facesEdges.clear();
	facesEdges.resize(mesh.faces.size());
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		int pointSize = face.ipoints.size();
		for(int j=0; j<pointSize; ++j){
			int ipoint = face.ipoints[j];
			if(find(
			pointsFaces[ipoint].begin(),
			pointsFaces[ipoint].end(),
			i) == pointsFaces[ipoint].end()){
				pointsFaces[ipoint].push_back(i);
			}
			facesEdges[i].push_back(-100);
		}
	}
	edgesPoint0.clear();
	edgesPoint1.clear();
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		int pointSize = face.ipoints.size();
		for(int j=0; j<pointSize; ++j){
			int ipoint0 = face.ipoints[j];
			int ipoint1 = ( j+1 == pointSize ? face.ipoints[0] : face.ipoints[j+1] );
			
			if(facesEdges[i][j]!=-100) continue;
			
			int new_edge_id = edgesPoint0.size();
			facesEdges[i][j] = new_edge_id;
			
			for(auto& iface : pointsFaces[ipoint0]){
				auto& other_face = mesh.faces[iface];
				int other_pointSize = other_face.ipoints.size();
				for(int k=0; k<other_pointSize; ++k){
					int other_ipoint0 = other_face.ipoints[k];
					int other_ipoint1 = ( k+1 == other_pointSize ? other_face.ipoints[0] : other_face.ipoints[k+1] );
					if(
					(other_ipoint0==ipoint0 && other_ipoint1==ipoint1) ||
					(other_ipoint0==ipoint1 && other_ipoint1==ipoint0)){
						facesEdges[iface][k] = new_edge_id;
					}
				}
			}
			for(auto& iface : pointsFaces[ipoint1]){
				auto& other_face = mesh.faces[iface];
				int other_pointSize = other_face.ipoints.size();
				for(int k=0; k<other_pointSize; ++k){
					int other_ipoint0 = other_face.ipoints[k];
					int other_ipoint1 = ( k+1 == other_pointSize ? other_face.ipoints[0] : other_face.ipoints[k+1] );
					if(
					(other_ipoint0==ipoint0 && other_ipoint1==ipoint1) ||
					(other_ipoint0==ipoint1 && other_ipoint1==ipoint0)){
						facesEdges[iface][k] = new_edge_id;
					}
				}
			}
			
			edgesPoint0.push_back(ipoint0);
			edgesPoint1.push_back(ipoint1);
			
			
		}
	}
	// cout << edgesPoint0.size() << endl;
	
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



void MASCH_Poly_AMR_Builder::searchOriginalPoints(
	MASCH_Mesh& mesh, 
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