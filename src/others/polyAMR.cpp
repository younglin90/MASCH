
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
int maxRefineCellPerBlockAMR,
double minVolume_AMR,
vector<vector<double>>& indicatorCriterion,
vector<vector<int>>& indicatorAMR_id, 
vector<bool>& boolCellRefine, 
vector<bool>& boolCellUnrefine,
vector<bool>& boolCellPreserved
){
    
	MASCH_Load load;
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int id_vol = controls.getId_cellVar("volume");
	
	int maxUnrefineCellPerBlockAMR = stoi(controls.dynamicMeshMap["AMR.maxUnrefineCellPerBlockAMR"]);
	if(controls.dynamicMeshMap.find("AMR.maxUnrefineCellPerBlockAMR")==controls.dynamicMeshMap.end()){
		maxUnrefineCellPerBlockAMR = 1000000;
		if(rank==0) cout << "#WARNING : no seraching, AMR.maxUnrefineCellPerBlockAMR" << endl;
	}
    
    
    string sRegionType;
    vector<double> boxMin(3,-1.e12), boxMax(3,1.e12);
    if(controls.dynamicMeshMap.find("AMR.region.type")==controls.dynamicMeshMap.end()){
        sRegionType = "box";
        if(rank==0) cout << "#WARNING : no seraching, AMR.region.type" << endl;
    }
    else{
        if(controls.dynamicMeshMap["AMR.region.type"]=="box"){
            vector<string> sboxMin = load.extractVector(controls.dynamicMeshMap["AMR.region.min"]);
            vector<string> sboxMax = load.extractVector(controls.dynamicMeshMap["AMR.region.max"]);
            boxMin.clear(); boxMax.clear();
            for(auto& item : sboxMin) boxMin.push_back(stod(item));
            for(auto& item : sboxMax) boxMax.push_back(stod(item));
        }
        else{
            if(rank==0) cout << "#WARNING : no seraching, AMR.region.type" << endl;
        }
    }
	
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
	
	
	// // int nRefineCellPerBlockAMR = 0;
    // std::random_device rd;
    // std::default_random_engine eng(rd());
    // std::uniform_real_distribution<double> distr(0.0, 1.0);

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
				
				// cout << "AAAA" << endl;
				
				// if(j==0) ++nRefineCellPerBlockAMR;
				// if(nRefineCellPerBlockAMR > maxRefineCellPerBlockAMR) {
					// boolCellRefine[i] = false;
				// }
			}
		}
		// if(mesh.cells[i].volume < minVolume_AMR) boolCellRefine[i] = false;
		if(var.cells[i][id_vol] < minVolume_AMR) boolCellRefine[i] = false;
		if(cell.level >= maxLevel) boolCellRefine[i] = false;
		if(cell.level < 0) boolCellRefine[i] = false;
        
        if(cell.x<boxMin[0] || cell.x>boxMax[0]) boolCellRefine[i] = false;
        if(cell.y<boxMin[1] || cell.y>boxMax[1]) boolCellRefine[i] = false;
        if(cell.z<boxMin[2] || cell.z>boxMax[2]) boolCellRefine[i] = false;
        
        
        // boolCellPreserved[i] = false;
        // boolCellRefine[i] = false;
        // if(distr(eng) > 0.9){
            // boolCellRefine[i] = true;
            // boolCellUnrefine[i] = false;
            // boolCellPreserved[i] = true;
        // }
        
        // boolCellUnrefine[i] = false;
        // if(distr(eng) > 0.4){
            // boolCellRefine[i] = false;
            // boolCellUnrefine[i] = true;
            // boolCellPreserved[i] = true;
        // }
		
		// if(nRefineCellPerBlockAMR > maxRefineCellPerBlockAMR) {
			// boolCellRefine[i] = false;
		// }
		
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
						MPI_Alltoallv( send_value.data(), mesh.countsSendProcFaces.data(), 
										mesh.displsSendProcFaces.data(), MPI_INT, 
										cPreserved_recv.data(), mesh.countsRecvProcFaces.data(), 
										mesh.displsRecvProcFaces.data(), MPI_INT, 
									   MPI_COMM_WORLD);
						MPI_Alltoallv( send_value2.data(), mesh.countsSendProcFaces.data(), 
										mesh.displsSendProcFaces.data(), MPI_INT, 
										cRefine_recv.data(), mesh.countsRecvProcFaces.data(), 
										mesh.displsRecvProcFaces.data(), MPI_INT, 
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
				if(
				my_level<tmp_maxLevel[i]-2 
				){
					boolCellRefine[i] = true;
					boolCellPreserved[i]=true;
					boolCellUnrefine[i]=false;
				
					// ++nRefineCellPerBlockAMR;
				}
				// if(nRefineCellPerBlockAMR > maxRefineCellPerBlockAMR) {
					// boolCellRefine[i] = false;
				// }
			}
		}
			
	}
	
	

	int nRefineCellPerBlockAMR = 0;
	int nUnrefineCellPerBlockAMR = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
        auto& cell = mesh.cells[i];
		if(var.cells[i][id_vol] < minVolume_AMR) boolCellRefine[i] = false;
		// if(var.cells[i][id_vol] < minVolume_AMR) boolCellUnrefine[i] = false;
		
		if(boolCellUnrefine[i]==true) ++nUnrefineCellPerBlockAMR;
		if(nUnrefineCellPerBlockAMR > maxUnrefineCellPerBlockAMR) boolCellUnrefine[i] = false;
		if(boolCellRefine[i]==true) ++nRefineCellPerBlockAMR;
		if(nRefineCellPerBlockAMR > maxRefineCellPerBlockAMR) boolCellRefine[i] = false;
		
        
        if(cell.x<boxMin[0] || cell.x>boxMax[0]) boolCellRefine[i] = false;
        if(cell.y<boxMin[1] || cell.y>boxMax[1]) boolCellRefine[i] = false;
        if(cell.z<boxMin[2] || cell.z>boxMax[2]) boolCellRefine[i] = false;
        
	} 
	
	
	
}


void MASCH_Poly_AMR_Builder::polyAMR(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	MASCH_Solver& solver,
	MASCH_Variables& var,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Load load;
	
	if( (iter+1) %
		stoi(controls.dynamicMeshMap["AMR.interval"]) == 0)
	{
		int maxLevel = stoi(controls.dynamicMeshMap["AMR.maxLevel"]);
		if(controls.dynamicMeshMap.find("AMR.maxLevel")==controls.dynamicMeshMap.end()){
			maxLevel = 0;
			if(rank==0) cout << "#WARNING : no seraching, AMR.maxLevel" << endl;
		}
		int maxBuffer = stoi(controls.dynamicMeshMap["AMR.maxBufferLayer"]);
		if(controls.dynamicMeshMap.find("AMR.maxBufferLayer")==controls.dynamicMeshMap.end()){
			maxBuffer = 0;
			if(rank==0) cout << "#WARNING : no seraching, AMR.maxBufferLayer" << endl;
		}
		double minVolume = stod(controls.dynamicMeshMap["AMR.minVolume"]);
		if(controls.dynamicMeshMap.find("AMR.minVolume")==controls.dynamicMeshMap.end()){
			minVolume = 0.0;
			if(rank==0) cout << "#WARNING : no seraching, AMR.minVolume" << endl;
		}
		int maxCells = stoi(controls.dynamicMeshMap["AMR.maxCells"]);
		if(controls.dynamicMeshMap.find("AMR.maxCells")==controls.dynamicMeshMap.end()){
			maxCells = 1000000;
			if(rank==0) cout << "#WARNING : no seraching, AMR.maxCells" << endl;
		}
		int maxRefineCellPerBlockAMR = stoi(controls.dynamicMeshMap["AMR.maxRefineCellPerBlockAMR"]);
		if(controls.dynamicMeshMap.find("AMR.maxRefineCellPerBlockAMR")==controls.dynamicMeshMap.end()){
			maxRefineCellPerBlockAMR = 1000000;
			if(rank==0) cout << "#WARNING : no seraching, AMR.maxRefineCellPerBlockAMR" << endl;
		}
		
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
			maxBuffer, maxLevel, maxCells, maxRefineCellPerBlockAMR, minVolume, 
			indicatorCriterion, indicatorAMR_id, 
			boolCellRefine, boolCellUnrefine, boolCellPreserved);
		
		// 리파인
		{
			// 원래 셀의 x,y,z 저장
			vector<vector<double>> org_xyz(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i) {
				auto& cell = mesh.cells[i];
				org_xyz[i].push_back(cell.x); org_xyz[i].push_back(cell.y); org_xyz[i].push_back(cell.z);
			}
			// cout << maxCells << endl;
			vector<vector<int>> child_new_cell_id_of_org;
			polyRefine(mesh, controls, 
				maxLevel, maxCells, minVolume, 
				indicatorCriterion, indicatorValues, 
				child_new_cell_id_of_org, 
				boolCellPreserved, boolCellRefine, boolCellUnrefine,
				0);
			
			controls.setGeometricOnlyCell_xyz(mesh);
			
			controls.resetVariableArray(mesh, var, org_xyz, child_new_cell_id_of_org, "refine");
			// controls.resetParcelPosition(mesh, var, child_new_cell_id_of_org, "refine");
		}
		
		// 언리파인
		{
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
			// controls.resetParcelPosition(mesh, var, child_org_cell_id_of_new, "unrefine");
		}
		
		
	}
	
	
	// 리파티셔닝
	if( (iter+1) %
		stoi(controls.dynamicMeshMap["AMR.intervalRepart"]) == 0)
	{
		int maxLevel = stoi(controls.dynamicMeshMap["AMR.maxLevel"]);
		// vector<vector<double>> dummy;
		// if(rank==0) cout << "| exe. Dynamic Load Balancing" << endl;
		vector<int> cell_ip(mesh.cells.size(),rank);
		mesh.repartParMETIS(size, cell_ip, mesh);
		vector<int> to_new_cell_id(mesh.cells.size());
		// if(rank==0) cout << "AA" << endl;
		vector<int> parcel_ip;
		mesh.repartitioning(cell_ip, maxLevel, to_new_cell_id, parcel_ip);
		// if(rank==0) cout << "BB" << endl;
		
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start gemetric only cell xyz" << endl;
		controls.setGeometricOnlyCell_xyz(mesh);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| end gemetric only cell xyz" << endl;
		
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset variables" << endl;
		controls.resetVariableArray(mesh, var, cell_ip, to_new_cell_id, parcel_ip, "repart");
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| end reset variables" << endl;
	}
	
	
	
	// 값 다시 구하기
	if(
	( (iter+1) %
		stoi(controls.dynamicMeshMap["AMR.interval"]) == 0) ||
	( (iter+1) %
		stoi(controls.dynamicMeshMap["AMR.intervalRepart"]) == 0)
	){
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values1" << endl;
		controls.setGeometric(mesh, var);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values2" << endl;
		var.setSparCSR(mesh, controls);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values3" << endl;
		solver.calcGradient.init(mesh, controls, var);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values4" << endl;
		
		solver.updateProcRightCellPrimValues(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values5" << endl;
		solver.updateCellAddiValues(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values6" << endl;
		solver.updateProcRightCellAddiValues(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values7" << endl;
        
        // ======================================
        solver.calcBoundFacePrimVal.clear();
        
        // controls.timeVaryingMappedFixedValueNCycle.clear();
        controls.timeVaryingMappedFixedValueTimeCycle.clear();
        // controls.timeVaryingMappedFixedValueTime1.clear();
        // controls.timeVaryingMappedFixedValueTime2.clear();
        controls.timeVaryingMappedFixedValueTimeOrder1.clear();
        // controls.timeVaryingMappedFixedValueTimeOrder2.clear();
        controls.timeVaryingMappedFixedValueTime.clear();
        // controls.timeVaryingMappedFixedValueValue1.clear();
        // controls.timeVaryingMappedFixedValueValue2.clear();
        controls.timeVaryingMappedFixedValueFileName.clear();
        controls.timeVaryingMappedFixedValueValueIter.clear();
        
        solver.setBoundaryFunctions(mesh, controls, var);
        // ======================================
        
		solver.updateBoundaryFacePrimValues(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values8" << endl;
		solver.gradientTerms(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values9" << endl;
		solver.updateProcRightCellGradValues(mesh, controls, var, 0);
        // MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0) cout << "| start reset values10" << endl;
        
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// MPI::Finalize();
		// MPI::Init(); 
	}
	
	
	
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

		MPI_Alltoallv( cLevel_send.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
					   
		MPI_Alltoallv( cRefine_send.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   cRefine_recv.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
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
					   
		MPI_Alltoallv( cRefine_send.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   cRefine_recv.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
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

		MPI_Alltoallv( cLevel_send.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
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