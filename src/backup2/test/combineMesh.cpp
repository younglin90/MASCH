#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string>
#include <map>
// #include <cstring>
#include <random>

using namespace std;

#include "parmetis.h" 
#include "scotch.h" 
// #include "metis.h" 
// #include "scotchf.h" 

#include "../mesh/mesh.h" 
#include "../mesh/geometric.h" 

void parMETIS_Graph_Partition(int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh);
// void parMETIS_Mesh_Partition(int nBlocks, vector<int>& idBlockCell);
void partitionFromSerial(int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh, vector<MASCH_Mesh>& newMesh);
void loadOnlyMeshVtu(string folder, MASCH_Mesh &mesh);
void combineMesh(vector<int>& idBlockCell, MASCH_Mesh &mesh);
void MASCH_Gatherv(vector<int>& my_value, vector<int>& value, vector<int>& displs);
void MASCH_Gatherv(vector<double>& my_value, vector<double>& value, vector<int>& displs);

int iteeeer = 0;

int main(int argc, char* argv[]) {
	

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	// map<string,string> mapArgv;
	
	// for(int i=1; i<argc; i+=2){
		// string first = argv[i];
		// string second;
		// if(i+1==argc){
			// second = "nan";
		// }
		// else{
			// second = argv[i+1];
		// }
		// mapArgv.insert(make_pair(first, second));
	// }
	
	// if( 
	// mapArgv.find("-help") != mapArgv.end() ||
	// mapArgv.find("-h") != mapArgv.end()
	// ){
		// if(rank==0){
			// cout << endl;
			// cout << "┌─────── Partitioning helper ─────────────────────────────── " << endl;
			// cout << "| -n \"int num.\"   : # of partition" << endl;
			// cout << "| -s \"real num.\"  : scale of mesh" << endl;
			// cout << "└───────────────────────────────────────────────────────────────── " << endl;
			// cout << endl;
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI::Finalize();
		// return EXIT_SUCCESS;
	// }
	
	// int nBlocks = stoi(mapArgv["-n"]);
	// double scaleMesh = stod(mapArgv["-s"]);
	
	
	MASCH_Mesh mesh;
	
	
	// MASCH_Mesh_Load load;
	loadOnlyMeshVtu("./grid/0/", mesh);
	
	
	vector<int> idBlockCell(mesh.cells.size(),0);
	
    std::random_device rd;
    std::default_random_engine eng(rd());
    // std::uniform_real_distribution<double> distr(0.0, 1.0);
    std::uniform_int_distribution<int> distr_int(0, 3);
	
	// idBlockCell[0] = 1;

		
		// for(auto& item : idBlockCell){
			// item = rank;
		// }
		// // // idBlockCell[1] = 1;
		// // if(rank==0){
			// // for(auto& item : idBlockCell){
				// // item = 0;
			// // }
		// // }
		// // // idBlockCell[2] = 3;
		
		// if(rank==0) idBlockCell[0] = 1;
		// // if(rank==0) idBlockCell[1] = 1;
		// if(rank==1) idBlockCell[0] = 0;
		// // if(rank==1) idBlockCell[1] = 0;
		// combineMesh(idBlockCell, mesh);
		
		// idBlockCell.resize(mesh.cells.size());
		// for(auto& item : idBlockCell){
			// item = rank;
		// }
		// if(rank==0) idBlockCell[0] = 1;
		// // if(rank==0) idBlockCell[1] = 1;
		// if(rank==1) idBlockCell[0] = 0;
		// // if(rank==1) idBlockCell[1] = 0;
		// combineMesh(idBlockCell, mesh);
	
	for(int iter=0; iter<200; ++iter){
	
		idBlockCell.resize(mesh.cells.size());
		for(auto& item : idBlockCell){
			// item = rank;
			item = distr_int(eng);
		}
		
		combineMesh(idBlockCell, mesh);
	
	}
	
	
		// for(auto& item : idBlockCell){
			// item = distr_int(eng);
		// }
		
		// combineMesh(idBlockCell, mesh);
		
		// ++iteeeer;
		
		// idBlockCell.resize(mesh.cells.size());
		// for(auto& item : idBlockCell){
			// item = distr_int(eng);
		// }
		
		// combineMesh(idBlockCell, mesh);

	// vector<int> idBlockCell(static_cast<int>(mesh.cells.size()),0);

	// parMETIS_Graph_Partition(nBlocks, idBlockCell);
	// // // parMETIS_Mesh_Partition(nBlocks, idBlockCell);

	// vector<MASCH_Mesh> newMesh(nBlocks,MASCH_Mesh());
	
	// partitionFromSerial(nBlocks, idBlockCell, newMesh);
	
	// // 스케일
	// for(auto& item : newMesh) {
		// for(auto& point : item.points){
			// point.x *= scaleMesh;
			// point.y *= scaleMesh;
			// point.z *= scaleMesh;
		// }
	// }

	cout.precision(20);
	// for(int ip=0; ip<nBlocks; ++ip)
	{

		// SEMO_Utility_Math math;
		// SEMO_Mesh_Geometric geometric;
		// geometric.init(newMesh[ip]);
		
		MASCH_Mesh_Save save;
		save.vtu("./grid/test/", rank, mesh);
	}


	MPI::Finalize();
	return EXIT_SUCCESS;
}



// void libSCOTCH_PartGraph (
	// vector<int>& xadj,
	// vector<int>& adjncy,
	// int* nparts,
	// vector<int>& part) {
	// /* Scotch graph object to interface with libScotch */
	// SCOTCH_Graph grafdat;
	// SCOTCH_Strat stradat;
	// int vertnbr, edgenbr;

	// SCOTCH_graphInit (&grafdat);

	// vertnbr = xadj.size()-1;
	// edgenbr = xadj[vertnbr];
	
	
	// SCOTCH_graphBuild (&grafdat,
					   // 0,
					   // vertnbr, // 6, 
					   // xadj.data(), // tmpxadj, 
					   // NULL, 
					   // NULL, 
					   // NULL,
					   // edgenbr, // 14, 
					   // adjncy.data(), // tmpadjncy, 
					   // NULL);

	// SCOTCH_graphCheck(&grafdat);
	// SCOTCH_stratInit(&stradat);
	// SCOTCH_graphPart(&grafdat, *nparts, &stradat, part.data());
	// SCOTCH_stratExit(&stradat);
	// SCOTCH_graphExit(&grafdat);

// }




void parMETIS_Graph_Partition(int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh){

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute ParMETIS ... ";
	}
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	options[0] = 0;
	options[1] = 0;
	options[2] = 0;
	
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.02);
	
	real_t ubvec = 1.02;
	
	
	vector<int> temp_vtxdist(nBlocks);
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	// vector<int> vtxdist(nBlocks+1);
	vector<int> vtxdist(nBlocks+1);
	vtxdist[0] = 0;
	for(int i=1; i<nBlocks+1; ++i){
		vtxdist[i] = vtxdist[i-1] + static_cast<int>(temp_vtxdist[i-1]);
	}
	
	
	vector<int> xadj(ncells+1,0);
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		for(auto& face : cell.faces()){
			if(
			(*face).getType() == MASCH_Face_Types::INTERNAL ||
			(*face).getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		}
		++numt;
		xadj[numt] = xadj[numt-1] + static_cast<int>(numt2);
	}
	
	
	vector<int> adjncy(xadj[ncells],0);
	
	vector<int> num_ncell(ncells,0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL) {
			
			int tnum = xadj[face.iL] + num_ncell[face.iL];
			adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iR);
			++num_ncell[face.iL];
			
			tnum = xadj[face.iR] + num_ncell[face.iR];
			adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iL);
			++num_ncell[face.iR];
			
		}
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// int tnum = xadj[face.owner] + num_ncell[face.owner];
			// int rank_neighbour = iRank_Recv[proc_num];
			// // cout << proc_num << " " << iRank_Recv.size() << endl;
			// adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// ++num_ncell[face.owner];
			
			// ++proc_num;
		// }
	}
	

	ParMETIS_V3_PartKway(
		vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, &wgtflag, &numflag,
		&ncon, &nBlocks, tpwgts, &ubvec,
		options, &objval, idBlockCell.data(), &comm);
		
		
	// idx_t vsize = NULL;
	// real_t itr = 1000.0;
	// ParMETIS_V3_AdaptiveRepart(
		// vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, &wgtflag, &numflag,
		// &ncon, &nBlocks, tpwgts, &ubvec, &itr, 
		// options, &objval, idBlockCell.data(), &comm);

		
	// libSCOTCH_PartGraph(xadj, adjncy, &nBlocks, idBlockCell);


	// for(int i=0; i<ncells; ++i){
		// idBlockCell[i] = 1;
	// }
	// for(int i=0; i<ncells/2; ++i){
		// idBlockCell[i] = 0;
	// }

	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}





// void parMETIS_Mesh_Partition(int nBlocks, vector<int>& idBlockCell){

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();

    // if(rank == 0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute ParMETIS ... ";
	// }
		
	// int ncells = mesh.cells.size();
	// int npoints = mesh.points.size();
	// int nfaces = mesh.faces.size();
	// int ncon=1;
	// int ncommon=3;
	// int objval;
	
	// int options[METIS_NOPTIONS];
	// METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	

	// MPI_Comm comm;
	// MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	
	// int wgtflag=0;
	// int numflag=0;
	// real_t tpwgts[nBlocks*ncon];
	// for(int i=0;i<nBlocks*ncon;++i) 
		// tpwgts[i]=1.0/nBlocks;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.05);
	
	
	
	// int temp_vtxdist[nBlocks];
    // MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	// int vtxdist[nBlocks+1];
	// vtxdist[0] = 0;
	// for(int i=1; i<nBlocks+1; ++i){
		// vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	// }
	
	
	
	
	// vector<int> xadj(ncells+1,0);
	// xadj[0] = 0;
	// for(int i=1; i<ncells+1; ++i){
		// xadj[i] = xadj[i-1] + mesh.cells[i-1].points.size();
	// }
	
	
	
	// // vector<int> num_ncell(ncells,0);
	// // vector<int> adjncy(xadj[ncells],0);
	// vector<int> adjncy;
	// int proc_num = 0;
	// for(int i=0; i<ncells; ++i){
		// for(auto& point : mesh.cells[i].points){
			// adjncy.push_back(point);
		// }
	// }
		
	// // for(auto& face : mesh.points){
		
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			// // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // adjncy[tnum] = vtxdist[rank] + face.neighbour;
			// // ++num_ncell[face.owner];
			
			// // // cout << face.neighbour << endl;
			
			// // tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			// // adjncy[tnum] = vtxdist[rank] + face.owner;
			// // ++num_ncell[face.neighbour];
			
		// // }
		// // // else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// // // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // // int rank_neighbour = iRank_Recv[proc_num];
			// // // // cout << proc_num << " " << iRank_Recv.size() << endl;
			// // // adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// // // ++num_ncell[face.owner];
			
			// // // ++proc_num;
		// // // }
	// // }
	

	// int ncommonnodes = 3;

	// // ParMETIS_V3_PartMeshKway(
		// // vtxdist, xadj.data(), adjncy.data(), NULL, &wgtflag, &numflag,
		// // &ncon, &ncommonnodes, &nBlocks, tpwgts, ubvec,
		// // options, &objval, idBlockCell.data(), &comm);

	// // objval = adjncy.size()+1;

	// vector<int> npart(npoints,0);

	// METIS_PartMeshDual(
		// &ncells, &npoints, xadj.data(), adjncy.data(), NULL, NULL,
		// &ncommonnodes, &nBlocks, NULL, NULL,
		// &objval, idBlockCell.data(), npart.data());

	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	
	
// }







void partitionFromSerial(int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh, vector<MASCH_Mesh>& newMesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int ncells = static_cast<int>(mesh.cells.size());
	int nfaces = static_cast<int>(mesh.faces.size());
	int npoints = static_cast<int>(mesh.points.size());
	int nboundary = static_cast<int>(mesh.boundaries.size());
	
	vector<vector<int>> idBlockPoint(npoints,vector<int>(0,0)); // point block id (copies)
	vector<int> nCellsLocal(nBlocks,0); // local total cells
	vector<int> idCellLocal(ncells,0); // local cell id
	// std::fill_n(nCellsLocal, nBlocks, 0);
	for(int i=0; i<ncells; ++i) {
		for(auto& j : mesh.cells[i].ipoints){
			if(find(
			idBlockPoint[j].begin(),idBlockPoint[j].end(),idBlockCell[i])==
			idBlockPoint[j].end()){
				idBlockPoint[j].push_back( idBlockCell[i] );
			}
		}
		idCellLocal[i] = nCellsLocal[ idBlockCell[i] ]++;
	}
	
	int nTotalLocalPointSize = 0; // total all # of local point
	for(auto& item : idBlockPoint) 
		nTotalLocalPointSize += item.size();

	// point distribution (CSR format)
	//
	// / gP : global points / lP : local points / BL : block id /
	//
	//     gP[0]     gP[1]      gP[2]      gP[3]   gP[4] ...
	// - - - - - - - - - - - - - - - - - - - - - - - - - ....
	// |  lP[0,5]  |lP[6,9]| lP[10,15] |  ....   |     |
	// |  BL[0,5]  |BL[6,9]| BL[10,15] |  ....   |     |
	// |           |       |           |         |     |
	// strPoints[0]|  strPoints[2]  strPoints[3] |  strPoints[5]  ...
	//         strPoints[1]                  strPoints[4]
	//
	vector<int> nPointsLocal(nBlocks,0); // local total points
	vector<int> idPointLocal(nTotalLocalPointSize,0); // local point id
	vector<int> idBlockPointLocal(nTotalLocalPointSize,0); // local point block id
	vector<int> strPoints(npoints+1,0); // start of each global point
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(auto& item : idBlockPoint[i]){
			int idBlock = item;
			idPointLocal[nIndex] = nPointsLocal[ idBlock ]++;
			idBlockPointLocal[nIndex] = idBlock;
			++nIndex;
		}
	}
	strPoints[mesh.points.size()] = nIndex;
	
	nPointsLocal.clear();
	nPointsLocal.resize(nBlocks,0);
	for(int i=0; i<nTotalLocalPointSize; ++i)
		++nPointsLocal[ idBlockPointLocal[i] ];


	vector<int> nDisplPoint(nBlocks,0);
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			int idLocal = idPointLocal[j];
			
			newMesh[idBlock].addPoint();
			newMesh[idBlock].points.back().x = mesh.points[i].x;
			newMesh[idBlock].points.back().y = mesh.points[i].y;
			newMesh[idBlock].points.back().z = mesh.points[i].z;
			
			// MPI 컨넥트 포인트들 저장 (procNo, local point id)
			vector<pair<int,int>> tmp_idBlock_idLocal;
			for(int k=strPoints[i]; k<strPoints[i+1]; ++k){
				if(j==k) continue;
				tmp_idBlock_idLocal.push_back(
				make_pair(idBlockPointLocal[k],idPointLocal[k]));
			}
			newMesh[idBlock].points.back().connPoints =
			tmp_idBlock_idLocal;
			
			
			
			++nDisplPoint[idBlock];
		}
		
		
	}
	
	
	// for(auto& item : newMesh[2].points){
		// if(item.connPoints.size()>0){
		// cout << "========" << endl;
		// }
		// for(auto& i : item.connPoints){
			// cout << i.first << " " << i.second << endl;
		// }
	// }
	
	// face setting
	vector<int> nFacesLocal(nBlocks,0); // total local faces
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.iL];
		
		++nFacesLocal[idBlockOwner];
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			int idBlockNeighbour = idBlockCell[face.iR];
			if(idBlockOwner != idBlockNeighbour){
				face.setType(MASCH_Face_Types::PROCESSOR); // set processor face
				++nFacesLocal[idBlockNeighbour];
			}
		}
	}
	
	
	vector<int> idFaceLocal(nBlocks,0);  // temporary local face id
	
	// internal faces
	vector<vector<int>> idFacePoint_Send(nBlocks,vector<int>(0,0));
	
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			int idBlock = idBlockCell[face.iL];
			++idFaceLocal[idBlock];
			
			vector<int> idFacePoint; // face's points id
			for(auto& i : face.ipoints){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( static_cast<int>(idPointLocal[j]) );
						break;
					}
				}
			}
			// save face's points size , local points id
			newMesh[idBlock].addFace();
			for(int i=0; i<idFacePoint.size(); ++i){
				newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
			}
		}
	}
	
	
	// 페이스 타입 BC
	vector<int> faceTypeBC(nfaces,-1);
	for(int ibcs=0; ibcs<mesh.boundaries.size(); ++ibcs){
		int str = mesh.boundaries[ibcs].startFace;
		int end = str + mesh.boundaries[ibcs].nFaces;
		for(int i=str; i<end; ++i){
			faceTypeBC[i] = ibcs;
		}
	}
	
	
	// boundary faces
	// local boundary face size
	vector<int> nFacesBoundaryLocal(nBlocks,0); 
	// local each boundary face size
	vector<vector<int>> nFacesEachBoundaryLocal(nBlocks,vector<int>(nboundary,0)); 
	// local each boundary face start
	vector<vector<int>> nStartFaceEachBoundaryLocal(nBlocks,vector<int>(nboundary,0)); 
	vector<vector<bool>> boolFaceEachBoundaryLocal(nBlocks,vector<bool>(nboundary,false)); 
	for(int id=0; id<mesh.faces.size(); ++id){
		auto& face = mesh.faces[id];
		
		if(face.getType() == MASCH_Face_Types::BOUNDARY){
			int idBlock = idBlockCell[face.iL];
			++nFacesBoundaryLocal[ idBlock ];
			
			++nFacesEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ];
			if(boolFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] == false){
				nStartFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] = idFaceLocal[ idBlock ];
				boolFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] = true;
			}
			++idFaceLocal[ idBlock ];

			// wirte bc
			vector<int> idFacePoint;
			for(auto& i : face.ipoints){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			// save face's points size , local points id
			// idFacePoint_Send[idBlock].push_back(face.points.size());
			newMesh[idBlock].addFace();
			for(int i=0; i<idFacePoint.size(); ++i){
				newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
				// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
		
	}
	
	
	// processor faces
	// local processor face size
	vector<int> nFacesProcessorLocal(nBlocks,0); 
	 // sending size from each processor to each processor 
	vector<vector<int>> sendCountsProcs(nBlocks,vector<int>(nBlocks,0));
	vector<int> idFacesProcessor; // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::PROCESSOR){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.iL];
			int idBlockNeighbour = idBlockCell[face.iR];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			++sendCountsProcs[idBlockOwner][idBlockNeighbour];
			++sendCountsProcs[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			int k = idFacesProcessor[jp];
			int m = idBlockCell[ mesh.faces[k].iL ];
			int n = idBlockCell[ mesh.faces[k].iR ];
			if(n==ip) {
				int idBlock = m;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].ipoints){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				// idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				newMesh[idBlock].addFace();
				for(int i=0; i<idFacePoint.size(); ++i){
					// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
					
				}
				
			}
			else if(m==ip) {
				int idBlock = n;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].ipoints){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				// idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				newMesh[idBlock].addFace();
				std::reverse(idFacePoint.begin()+1,idFacePoint.end());
				for(int i=0; i<idFacePoint.size(); ++i){
					// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
				}
					
			}
		}
	}
	
	
	
	
	// write owner of internal face
	vector<vector<int>> idOwnerLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			int j = face.iL;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::BOUNDARY){
			int j = face.iL;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(int i=0; i<nBlocks; ++i){
		for(int j=0; j<idFacesProcessor.size(); ++j){
			int k = idFacesProcessor[j];
			int m = idBlockCell[ mesh.faces[k].iL ];
			int n = idBlockCell[ mesh.faces[k].iR ];
			
			if(n==i) {
				idOwnerLocal[m].push_back( idCellLocal[ mesh.faces[k].iL ] );
			}
			else if(m==i) {
				idOwnerLocal[n].push_back( idCellLocal[ mesh.faces[k].iR ] );
			}
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		int i=0;
		for(auto& iL : idOwnerLocal[ip]){
			newMesh[ip].faces[i].iL = iL;
			++i;
		}
	}
	
	// write neighbour of internal face
	vector<vector<int>> idNeighbourLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(
		face.getType() == MASCH_Face_Types::INTERNAL
		){
			int j = face.iR;
			int i = idBlockCell[j];
			idNeighbourLocal[i].push_back(idCellLocal[j]);
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		int i=0;
		for(auto& iR : idNeighbourLocal[ip]){
			newMesh[ip].faces[i].iR = iR;
			// newMesh[ip].faces[i].thereR = true;
			newMesh[ip].faces[i].setType(MASCH_Face_Types::INTERNAL);
			++i;
		}
	}
	

	// write of boundary faces
	// vector<int> nFaceBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
	// vector<int> nStartBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
	// int temp_bound=0;
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int ibcs=0; ibcs<mesh.boundaries.size(); ++ibcs){
			newMesh[ip].addBoundary();
			
			string bcnames = mesh.boundaries[ibcs].name;
			
			bcnames.erase(
			std::find_if(bcnames.rbegin(), 
			bcnames.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), 
			bcnames.end());
			
			bcnames.erase(
			bcnames.begin(), 
			std::find_if(bcnames.begin(), bcnames.end(), 
			std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			newMesh[ip].boundaries.back().name = bcnames;
			
			if( nFacesEachBoundaryLocal[ip][ibcs] > 0) {
				
				newMesh[ip].boundaries.back().nFaces = nFacesEachBoundaryLocal[ip][ibcs];
				newMesh[ip].boundaries.back().startFace = nStartFaceEachBoundaryLocal[ip][ibcs];
				// newMesh[ip].boundaries.back().thereR = false;
				newMesh[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			}
			else{
				newMesh[ip].boundaries.back().nFaces = 0;
				newMesh[ip].boundaries.back().startFace = 0;
				// newMesh[ip].boundaries.back().thereR = false;
				newMesh[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			}
		} 
	}
	
	
	
	
	// write of processor faces
	for(int ip=0; ip<nBlocks; ++ip){
		int n = nFacesLocal[ip];
		for(int jp=0; jp<nBlocks; ++jp){
			n -= sendCountsProcs[ip][jp];
		}
		
		for(int jp=0; jp<nBlocks; ++jp){
			if( sendCountsProcs[ip][jp] > 0) {
				newMesh[ip].addBoundary();
				
				string bcnames = "procBoundary" + to_string(ip) + "to" + to_string(jp);
				
				newMesh[ip].boundaries.back().name = bcnames;
				
				newMesh[ip].boundaries.back().nFaces = sendCountsProcs[ip][jp];
				newMesh[ip].boundaries.back().startFace = n;
				newMesh[ip].boundaries.back().myProcNo = ip;
				newMesh[ip].boundaries.back().rightProcNo = jp;
				// newMesh[ip].boundaries.back().thereR = true;
				newMesh[ip].boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
				
				n += sendCountsProcs[ip][jp];
			}
		}
	}
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		int maxBCNum = newMesh[ip].boundaries.size()-1;
		newMesh[ip].boundaries[maxBCNum].startFace = newMesh[ip].faces.size()-newMesh[ip].boundaries[maxBCNum].nFaces;
		for(int i=maxBCNum-1; i>=0; --i){
			newMesh[ip].boundaries[i].startFace = newMesh[ip].boundaries[i+1].startFace-newMesh[ip].boundaries[i].nFaces;
		}
	}
	
	// //==========================================
	
	vector<int> nPointsInf;
	vector<int> nFacesInf;
	vector<int> nCellsInf;
	vector<int> nFacesIntInf;
	vector<int> nFacesBCInf;
	vector<int> nFacesProcInf;
	vector<int> nTriangleInf;
	vector<int> nQuadrangleInf;
	vector<int> nPolygonInf;
	vector<int> nTetrahedronInf;
	vector<int> nHexahedronInf;
	vector<int> nPrismInf;
	vector<int> nPyramidInf;
	vector<int> nPolyhedronInf;
	
	int nbcs=0;
	
	for(int ip=0; ip<nBlocks; ++ip){

		newMesh[ip].check();
		newMesh[ip].setFaceTypes();
		newMesh[ip].buildCells();
		newMesh[ip].connectFacetoPointsCells();
		newMesh[ip].connectCelltoFaces();
		newMesh[ip].connectCelltoPoints();
		newMesh[ip].setCountsProcFaces();
		newMesh[ip].setDisplsProcFaces();
		
		nPointsInf.push_back(newMesh[ip].points.size());
		nFacesInf.push_back(newMesh[ip].faces.size());
		nCellsInf.push_back(newMesh[ip].cells.size());
	
		nbcs=0;
		int nFacesBC = 0;
		int nFacesProc = 0;
		for(auto& boundary : newMesh[ip].boundaries){
			
			if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				nFacesBC += boundary.nFaces;
				++nbcs;
			}
			else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
				nFacesProc += boundary.nFaces;
			}
			else{
				cout << "#ERROR : 868" << endl;
			}
		}
		
		
		// faces type
		int nFacesInt = 0;
		int nTriangle = 0;
		int nQuadrangle = 0;
		int nPolygon = 0;
		for(auto& face : newMesh[ip].faces){
			
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				++nFacesInt;
			}
			
			if(face.ipoints.size() == 3){
				++nTriangle;
			}
			else if(face.ipoints.size() == 4){
				++nQuadrangle;
			}
			else{
				++nPolygon;
			}
		}
		
		
		
		// cells type
		int nTetrahedron = 0;
		int nHexahedron = 0;
		int nPrism = 0;
		int nPyramid = 0;
		int nPolyhedron = 0;
		for(auto& cell : newMesh[ip].cells){
			if( cell.ipoints.size() == 4 && cell.ifaces.size() == 4 ){
				++nTetrahedron;
			}
			else if( cell.ipoints.size() == 5 && cell.ifaces.size() == 5 ){
				++nPyramid;
			}
			else if( cell.ipoints.size() == 6 && cell.ifaces.size() == 5 ){
				++nPrism;
			}
			else if( cell.ipoints.size() == 8 && cell.ifaces.size() == 6 ){
				++nHexahedron;
			}
			else {
				++nPolyhedron;
			}
		}
		
		nFacesIntInf.push_back(nFacesInt);
		nFacesBCInf.push_back(nFacesBC);
		nFacesProcInf.push_back(nFacesProc);
		nTriangleInf.push_back(nTriangle);
		nQuadrangleInf.push_back(nQuadrangle);
		nPolygonInf.push_back(nPolygon);
		nTetrahedronInf.push_back(nTetrahedron);
		nHexahedronInf.push_back(nHexahedron);
		nPrismInf.push_back(nPrism);
		nPyramidInf.push_back(nPyramid);
		nPolyhedronInf.push_back(nPolyhedron);
	}
	

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| partitioned MPI size : " << nBlocks << endl;
		cout << "| points size : ";
		for(auto& i : nPointsInf) cout << i << " | ";
		cout << endl;
		// gatherValue = mesh.faces.size();
		cout << "| faces size : ";
		for(auto& i : nFacesInf) cout << i << " | ";
		cout << endl;
		// gatherValue = mesh.cells.size();
		cout << "| cells size : ";
		for(auto& i : nCellsInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nFacesInt;
		cout << "| internal faces size : ";
		for(auto& i : nFacesIntInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nFacesBC;
		cout << "| boundary faces size : ";
		for(auto& i : nFacesBCInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nFacesProc;
		cout << "| processor faces size : ";
		for(auto& i : nFacesProcInf) cout << i << " | ";
		cout << endl;
		cout << "| boundary types : " << nbcs << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nTriangle;
		cout << "| Triangle faces : ";
		for(auto& i : nTriangleInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nQuadrangle;
		cout << "| Quadrangle faces : ";
		for(auto& i : nQuadrangleInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPolygon;
		cout << "| Polygon faces : ";
		for(auto& i : nPolygonInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nTetrahedron;
		cout << "| Tetrahedron cells : ";
		for(auto& i : nTetrahedronInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nHexahedron;
		cout << "| Hexahedron cells : ";
		for(auto& i : nHexahedronInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPrism;
		cout << "| Prism cells : ";
		for(auto& i : nPrismInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPyramid;
		cout << "| Pyramid cells : ";
		for(auto& i : nPyramidInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPolyhedron;
		cout << "| Polyhedron cells : ";
		for(auto& i : nPolyhedronInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
    }
	
	
}










void loadOnlyMeshVtu(string folder, MASCH_Mesh &mesh){
	
	

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load vtu data files ... ";
	}
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	// points
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	string nextToken;
	bool boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				boolCompress = true;
			}
			break;
		}
	}
	inputFile.clear();
	
	
	
	
	// vector<string> volFracName;
	// volFracName.push_back(controls.name[controls.VF[0]]);
	
	// vector<string> massFracName;
	// massFracName.push_back(controls.name[controls.MF[0]]);
	
	SEMO_Mesh_Load load;
	
	vector<double> NodeCoordinates;
	load.loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	// cout << NodeCoordinates.size() << endl;
	
	vector<int> connectivity;
	load.loadDatasAtVTU(inputFile, "connectivity", connectivity);
	
	vector<int> offsets;
	load.loadDatasAtVTU(inputFile, "offsets", offsets);
	
	vector<int> faces;
	load.loadDatasAtVTU(inputFile, "faces", faces);
	
	vector<int> faceoffsets;
	load.loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
	
	vector<int> owner;
	load.loadDatasAtVTU(inputFile, "owner", owner);
	
	vector<int> neighbour;
	load.loadDatasAtVTU(inputFile, "neighbour", neighbour);
	
	vector<string> bcName;
	load.loadDatasAtVTU(inputFile, "bcName", bcName);
	
	vector<int> bcStartFace;
	load.loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
	
	vector<int> bcNFaces;
	load.loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
	
	vector<int> bcNeighbProcNo;
	load.loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);
	
	vector<int> connPoints;
	load.loadDatasAtVTU(inputFile, "connPoints", connPoints);


	inputFile.close();
	
	
	
	
	

	
	for(int i=0; i<NodeCoordinates.size()/3; ++i){
		mesh.addPoint();
		mesh.points.back().x = NodeCoordinates[i*3+0];
		mesh.points.back().y = NodeCoordinates[i*3+1];
		mesh.points.back().z = NodeCoordinates[i*3+2];
	}
	NodeCoordinates.clear();
	
	
	int n=0;
	int ncells=0;
	mesh.faces.clear();
	for(auto& i : owner){
		mesh.addFace();
		mesh.faces.back().iL = i;
		ncells = max(ncells, mesh.faces.back().iL);
	}
	owner.clear();
	
	// cout << ncells << endl;
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();
	}
	
	
	
	n=0;
	for(auto& i : neighbour){
		mesh.faces[n].iR = i;
		mesh.faces[n].setType(MASCH_Face_Types::INTERNAL);
		++n;
	}
	neighbour.clear();
	
	
	int m=0;
	n=0;
	for(auto& i : offsets){
		for(int j=n; j<i; ++j){
			int point = connectivity[j];
			// cout << m << endl;
			mesh.cells[m].ipoints.push_back( static_cast<int>(point) );
			// if(rank==1) cout << point << endl;
		}
		
		// if(rank==0 && mesh.cells[m].ipoints.size()!=8) cout << mesh.cells[m].ipoints.size() << endl;
		
		n=i;
		++m;
	}
	
	
	
	n=0;
	int nFacesInt=0;
	for(auto& face : mesh.faces){
		// if(face.iR != -1){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			mesh.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
			++nFacesInt;
		}
		else{
			mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
		}
		++n;
	}
	
	
	m=0;
	n=0;
	for(auto& i : faceoffsets){
		// if(faces[n]>5) cout << faces[n] << endl;
		int N=0;
		int face_size = faces[m+N];
		for(int j=0; j<face_size; ++j){
			int face = mesh.cells[n].ifaces[j];
			++N;
			int point_size = faces[m+N];
			for(int k=0; k<point_size; ++k){
				++N;
				int point = faces[m+N];
				if(mesh.faces[ face ].ipoints.size() == point_size) continue;
				mesh.faces[ face ].ipoints.push_back( static_cast<int>(point) );
				// if(rank==1) cout << point << endl;
			}
		}
		m=i;
		++n;
	}
	faces.clear();
	faceoffsets.clear();
	
	
	
	n=0;
	for(auto& startFace : bcStartFace){
		
		mesh.addBoundary();
		
		string tmp_bcName = bcName[n];
		tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
		tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
		tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
		tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
		mesh.boundaries.back().name = tmp_bcName;
		mesh.boundaries.back().startFace = bcStartFace[n];
		mesh.boundaries.back().nFaces = bcNFaces[n];
		if(bcNeighbProcNo[n]<0){
			mesh.boundaries.back().rightProcNo = 0;
			mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
		}
		else{
			mesh.boundaries.back().rightProcNo = bcNeighbProcNo[n];
			mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
		}
		mesh.boundaries.back().myProcNo = rank;
		// if(bcNeighbProcNo[n] < 0){
			// mesh.boundaries.back().myProcNo = -1;
		// }
		
		++n;
		
	}
	
	int maxBCNum = mesh.boundaries.size()-1;
	mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
	for(int i=maxBCNum-1; i>=0; --i){
		mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
	}
	
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
	for(int i=0; i<connPoints.size(); ++i){
		// cout << mesh.points.size() << " " << connPoints[i] << " " << connPoints[i+1] << " " << connPoints[i+2] << " " << endl;
		mesh.points[connPoints[i]].connPoints.push_back(
		make_pair(connPoints[i+1],connPoints[i+2])
		);
		++i; ++i;
	}
	
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	for(auto& i : mesh.cells){
		
		if(rank==0 && i.ipoints.size()!=8) cout << i.ipoints.size() << endl;
		if(rank==0 && i.ifaces.size()!=6) cout << i.ifaces.size() << endl;
		
	}
	
	mesh.check();
	mesh.setFaceTypes();
	// mesh.buildCells();
	// mesh.connectFacetoPointsCells();
	// mesh.connectCelltoFaces();
	// mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	
	
	mesh.informations();

		
	
	
	
}










void MASCH_MPI_Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs);
void MASCH_MPI_Alltoallv(vector<vector<int>>& inp_send_value, vector<int>& recv_value, vector<int>& displs);
void MASCH_MPI_Alltoallv(vector<vector<double>>& inp_send_value, vector<vector<double>>& recv_value);
void MASCH_MPI_Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value);



class MASCH_Load_Balancing_Inform{
public:
	
	using vec1 = vector<int>;
	using vec2 = vector<vector<int>>;
	using vec3 = vector<vector<vector<int>>>;
	
	vec2 new_ifacesIN;
	vec2 new_faceIN_iL, new_faceIN_iR;
	vec3 new_faceBC_iL;
	vec3 new_facePR_iL;
	vec2 new_faceIN_ipoints;
	vec3 new_faceBC_ipoints;
	vec3 new_facePR_ipoints;
	vec2 new_ipoints;
	vector<vector<double>> new_point_x;
	vector<vector<double>> new_point_y;
	vector<vector<double>> new_point_z;
	
	vec2 new_icells;
	
	vector<string> new_boundaries_name;
	vec2 new_boundaries_nFaces;
	
};



void combineMesh(vector<int>& idBlockCell, MASCH_Mesh &mesh){
	
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	
	
	// {
		// vector<double> send_x0,send_y0,send_z0;
		// vector<double> send_x1,send_y1,send_z1;
		// vector<double> send_x2,send_y2,send_z2;
		// vector<double> send_x3,send_y3,send_z3;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// send_x0.push_back(mesh.points[face.ipoints[0]].x);
				// send_y0.push_back(mesh.points[face.ipoints[0]].y);
				// send_z0.push_back(mesh.points[face.ipoints[0]].z);
				// send_x1.push_back(mesh.points[face.ipoints[1]].x);
				// send_y1.push_back(mesh.points[face.ipoints[1]].y);
				// send_z1.push_back(mesh.points[face.ipoints[1]].z);
				// send_x2.push_back(mesh.points[face.ipoints[2]].x);
				// send_y2.push_back(mesh.points[face.ipoints[2]].y);
				// send_z2.push_back(mesh.points[face.ipoints[2]].z);
				// send_x3.push_back(mesh.points[face.ipoints[3]].x);
				// send_y3.push_back(mesh.points[face.ipoints[3]].y);
				// send_z3.push_back(mesh.points[face.ipoints[3]].z);
			// }
		// }
		// vector<double> recv_x0,recv_y0,recv_z0;
		// vector<double> recv_x1,recv_y1,recv_z1;
		// vector<double> recv_x2,recv_y2,recv_z2;
		// vector<double> recv_x3,recv_y3,recv_z3;
		// recv_x0.resize(send_x0.size());
		// recv_y0.resize(send_y0.size());
		// recv_z0.resize(send_z0.size());
		// recv_x1.resize(send_x1.size());
		// recv_y1.resize(send_y1.size());
		// recv_z1.resize(send_z1.size());
		// recv_x2.resize(send_x2.size());
		// recv_y2.resize(send_y2.size());
		// recv_z2.resize(send_z2.size());
		// recv_x3.resize(send_x3.size());
		// recv_y3.resize(send_y3.size());
		// recv_z3.resize(send_z3.size());
		// MPI_Alltoallv( send_x0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_x0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_y0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_y0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_z0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_z0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
					   
		// MPI_Alltoallv( send_x1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_x1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_y1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_y1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_z1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_z1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
					   
		// MPI_Alltoallv( send_x2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_x2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_y2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_y2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_z2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_z2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
					   
		// MPI_Alltoallv( send_x3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_x3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_y3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_y3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_z3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_z3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
					   

		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// if(rank==0){
				// cout << send_x0[ip] << " " << send_y0[ip] << " " << send_z0[ip] << endl;
				// cout << recv_x0[ip] << " " << recv_y0[ip] << " " << recv_z0[ip] << endl;
				// }
				
				// ++ip;
			// }
		// }
		
		// cout << endl;
		// cout << endl;
		// cout << endl;
		
	// }
	
	
	
	
	
	
	vector<MASCH_Mesh> meshNew;
	MASCH_Mesh meshComb;
	
	
	
	vector<int> recv_idBlockCell;
	vector<int> recv_rank;
	{
		vector<int> send_idBlockCell;
		vector<int> send_rank;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType()==MASCH_Face_Types::PROCESSOR){
				send_idBlockCell.push_back(idBlockCell[face.iL]);
				send_rank.push_back(rank);
			}
		}
		recv_idBlockCell.resize(send_idBlockCell.size());
		MPI_Alltoallv( send_idBlockCell.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   recv_idBlockCell.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
		recv_rank.resize(send_rank.size());
		MPI_Alltoallv( send_rank.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   recv_rank.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	
	
	
	int PR2IN = -1;
	int BC2BC = -2;
	int IN2PR = -4;
	int PR2PR = -3;
	vector<int> face_state;
	
	
	
	// 포인트
	{
		
		// 넘길 포인트 process 정립
		vector<vector<int>> send_localPoint_proc(mesh.points.size());
		for(int i=0, ip=0; i<mesh.cells.size(); ++i){
			int proc = idBlockCell[i];
			for(auto& ipoint : mesh.cells[i].ipoints){
				auto& point = send_localPoint_proc[ipoint];
				if(find(point.begin(),point.end(),proc)==point.end()){
					point.push_back(proc);
				}
			}
		}
		
		
		// 옮겨질 포인트 정보들 정리
		vector<vector<double>> send_localPoint_xyz(size);
		vector<vector<pair<int,int>>> send_localPoint_proc_id(mesh.points.size());
		vector<int> send_localPoint_n(size,0);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& proc : send_localPoint_proc[i]){
				send_localPoint_proc_id[i].push_back(make_pair(proc,send_localPoint_n[proc]++));
				send_localPoint_xyz[proc].push_back(point.x);
				send_localPoint_xyz[proc].push_back(point.y);
				send_localPoint_xyz[proc].push_back(point.z);
			}
		}
		
		// 포인트 x,y,z 옮기기
		vector<vector<double>> recv_localPoint_xyz;
		MASCH_MPI_Alltoallv(send_localPoint_xyz, recv_localPoint_xyz);
		
		// 옮겨진 포인트 x,y,z 넣기
		vector<int> str_points_glo(size+1,0);
		for(int ip=0, iter_glob=0; ip<size; ++ip){
			int tmp_size = recv_localPoint_xyz[ip].size()/3;
			str_points_glo[ip+1] = tmp_size;
			for(int i=0, iter=0, iter2=0; i<tmp_size; ++i){
				meshComb.addPoint();
				meshComb.points.back().x = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().y = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().z = recv_localPoint_xyz[ip][iter++];
				
			}
		}
		// 각 proc에 대한 포인트 시작지점
		for(int ip=0; ip<size; ++ip){
			str_points_glo[ip+1] = str_points_glo[ip+1] + str_points_glo[ip];
		}
		
		
		
		//========================================================
		// connPoints 대한 정보들 
		vector<vector<int>> send_localPoint_connId(size);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [proc, id] : point.connPoints){
				send_localPoint_connId[proc].push_back(id);
				send_localPoint_connId[proc].push_back(i);
				send_localPoint_connId[proc].push_back(send_localPoint_proc_id[i].size());
				for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
					send_localPoint_connId[proc].push_back(tmp_proc);
					send_localPoint_connId[proc].push_back(tmp_id_local);
				}
			}
		}
		
		// connPoints 대한 정보들 넘기기 
		vector<vector<int>> recv_localPoint_connId;
		MASCH_MPI_Alltoallv(send_localPoint_connId, recv_localPoint_connId);
		
		// 받은 connPoints 대한 정보들에 대해서, 같은 proc으로 이동하는 포인트들 정보 입력 
		vector<vector<int>> send_localPoint_toProc_toId(size);
		// vector<vector<int>> send_globalPoint_toProc_toId(size);
		for(int ip=0; ip<size; ++ip){
			for(int i=0; i<recv_localPoint_connId[ip].size(); ++i){
				int ipoint = recv_localPoint_connId[ip][i];
				int right_i = recv_localPoint_connId[ip][++i];
				int tmp_size = recv_localPoint_connId[ip][++i];
				for(int j=0; j<tmp_size; ++j){
					int tmp_proc = recv_localPoint_connId[ip][++i];
					int tmp_id_local = recv_localPoint_connId[ip][++i];
					for(auto& [proc, id] : send_localPoint_proc_id[ipoint]){
						if(proc==tmp_proc){
							// if(n_tmp_nnnn2[ipoint] == 1){
								// if(rank==0) cout << ipoint << " " << proc << " " << id << " " << ip << " " << tmp_id_local << endl;
								send_localPoint_toProc_toId[proc].push_back(id); // proc 블록 에서 자신의 local id  
								send_localPoint_toProc_toId[proc].push_back(ip); // proc 블록 에서 자신과 중복한 포인트의 local proc
								send_localPoint_toProc_toId[proc].push_back(tmp_id_local); // proc 블록 에서 자신과 중복한 포인트의 local id
							// }
						}
					}
				}
			}
		}
		
		// 같은 proc으로 이동하는 포인트들 정보 넘기기
		vector<vector<int>> recv_localPoint_toProc_toId;
		MASCH_MPI_Alltoallv(send_localPoint_toProc_toId, recv_localPoint_toProc_toId);
		
		
		// 삭제되는 중복 포인트들 찾기
		vector<bool> deletePoints(str_points_glo[size],false);
		vector<vector<int>> reorder_my_id_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightProcNo_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightId_local(str_points_glo[size]);
		vector<vector<int>> reorder_rightId_global(str_points_glo[size]);
		{
			for(int ip=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=0, iter=0; i<recv_localPoint_toProc_toId[ip].size()/3; ++i){
					int my_id_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int my_id_global = str_points_glo[ip] + my_id_local;
					int rightProc_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_global = str_points_glo[rightProc_local] + rightId_local;
					
					reorder_my_id_local[my_id_global].push_back(my_id_local);
					reorder_rightProcNo_local[my_id_global].push_back(rightProc_local);
					reorder_rightId_local[my_id_global].push_back(rightId_local);
					reorder_rightId_global[my_id_global].push_back(rightId_global);
					
					// 중복 및 삭제
					if(rightId_global<my_id_global){
						deletePoints.at(my_id_global) = true;
					}
				}
			}
		}
		
		
		// 로컬 포인트 id -> 글로벌 포인트 id
		vector<vector<int>> points_id_local2global(size);
		{
			for(int ip=0, iter=0, before_nd=0, tmp_glob=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=str; i<end; ++i){
					points_id_local2global[ip].push_back(tmp_glob);
					if(deletePoints[i]==false) {
						++tmp_glob;
					}
					else{
						++before_nd;
					}
				}
			}
			for(int ip=0; ip<size; ++ip){
				int str = str_points_glo[ip];
				int end = str_points_glo[ip+1];
				for(int i=str, iter=0; i<end; ++i){
					int iter = 0;
					int my_id_global = i;
					int min_rightId_global = str_points_glo[size];
					
					// if(rank==0){
						// cout << i << endl;
						// cout << 
						// meshComb.points[i].x << " " <<
						// meshComb.points[i].y << " " <<
						// meshComb.points[i].z << " " <<
						// endl;
					// }
					
					for(auto& my_id_local : reorder_my_id_local[i]){
						
						int rightProc_local = reorder_rightProcNo_local[i][iter];
						int rightId_local = reorder_rightId_local[i][iter];
						int rightId_global = reorder_rightId_global[i][iter];
						
						int copy_tmp = points_id_local2global[ip][my_id_local];
						// if(rank==0) cout << i << " " << rightId_global << " " << copy_tmp << endl;
						// if(rank==0) cout << "Delete = " << deletePoints[i] << endl;
						if(rightId_global>i){
							// if(rank==0){
								// cout << 
								// meshComb.points[rightId_global].x << " " <<
								// meshComb.points[rightId_global].y << " " <<
								// meshComb.points[rightId_global].z << " " <<
								// endl;
							// }
							// if(rank==0) cout << i << " " << rightId_global << " " << copy_tmp << endl;
							points_id_local2global[rightProc_local][rightId_local] = copy_tmp;
						}
						
						++iter;
					}
				}
			}
		}
		
		
		// // 디버그
		// {
		
			// for(int ip=0; ip<size; ++ip){
				// int str = str_points_glo[ip];
				// int end = str_points_glo[ip+1];
				// for(int i=str, iter=0; i<end; ++i){
					// int id_glo = points_id_local2global[ip][i-str];
					// if(rank==0){
						// cout << 
						// i << " " <<
						// id_glo << " " <<
						// meshComb.points[i].x << " " <<
						// meshComb.points[i].y << " " <<
						// meshComb.points[i].z << " " <<
						// endl;
					// }
					
				// }
			// }
		
		
		
		// }
		
		//========================================================
		
		
		// vector<vector<int>> send_localPoint_connId(size);
		// for(int i=0, ip=0; i<mesh.points.size(); ++i){
			// auto& point = mesh.points[i];
			// for(auto& [proc, id] : point.connPoints){
				// send_localPoint_connId[proc].push_back(id);
				// send_localPoint_connId[proc].push_back(i);
				// send_localPoint_connId[proc].push_back(send_localPoint_proc_id[i].size());
				// for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
					// send_localPoint_connId[proc].push_back(tmp_proc);
					// send_localPoint_connId[proc].push_back(tmp_id_local);
				// }
			// }
		// }
		
		// connPoints 대한 정보들 
		vector<vector<int>> send_conn_proc2(size);
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [proc, id] : point.connPoints){
				send_conn_proc2[proc].push_back(id);
				send_conn_proc2[proc].push_back(send_localPoint_proc_id[i].size());
				for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
					send_conn_proc2[proc].push_back(tmp_proc);
					send_conn_proc2[proc].push_back(rank);
					send_conn_proc2[proc].push_back(tmp_id_local);
				}
			}
		}
		
		// connPoints 대한 정보들 넘기기
		vector<vector<int>> recv_conn_proc2;
		MASCH_MPI_Alltoallv(send_conn_proc2, recv_conn_proc2);
		
		
		// 각 포인트에 대한, 넘겨진 포인트의 옆 proc_glo, proc_loc, id_loc 저장
		vector<vector<int>> connPoints_all_send_procNo_glo(mesh.points.size());
		vector<vector<int>> connPoints_all_send_procNo_loc(mesh.points.size());
		vector<vector<int>> connPoints_all_send_id_loc(mesh.points.size());
		for(int ip=0, id=0; ip<size; ++ip){
			int iter=0;
			while(iter<recv_conn_proc2[ip].size()){
				int id_glo = recv_conn_proc2[ip][iter++];
				int tmp_size = recv_conn_proc2[ip][iter++];
				for(int i=0; i<tmp_size; ++i){
					int rightProcNo_glo = recv_conn_proc2[ip][iter++];
					int rightProcNo_loc = recv_conn_proc2[ip][iter++];
					int rightId_loc = recv_conn_proc2[ip][iter++];
					
					connPoints_all_send_procNo_glo[id_glo].push_back(rightProcNo_glo);
					connPoints_all_send_procNo_loc[id_glo].push_back(rightProcNo_loc);
					connPoints_all_send_id_loc[id_glo].push_back(rightId_loc);
					
				}
			}
		}
		
		// 각 포인트에 대한, 자기 자신 포인트의 옆 proc_glo, proc_loc, id_loc 저장
		for(int i=0, ip=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [tmp_proc, tmp_id_local] : send_localPoint_proc_id[i]){
				connPoints_all_send_procNo_glo[i].push_back(tmp_proc);
				connPoints_all_send_procNo_loc[i].push_back(rank);
				connPoints_all_send_id_loc[i].push_back(tmp_id_local);
			}
			
		}
		
		// for(int i=0, ip=0; i<mesh.points.size(); ++i){
			// if(rank==0){
				// int iter=0;
				// for(auto& item : connPoints_all_send_procNo_glo[i]){
					// cout << 
					// connPoints_all_send_procNo_glo[i][iter] << " " <<
					// connPoints_all_send_procNo_loc[i][iter] << " " <<
					// connPoints_all_send_id_loc[i][iter] << " " << 
					// endl;
					// ++iter;
				// }
			// }
		// }
		
		// cout << endl;
		// cout << endl;
		// cout << endl;
		// cout << endl;
		
		// for(int i=0, ip=0; i<mesh.points.size(); ++i){
			// if(rank==1){
				// int iter=0;
				// for(auto& item : connPoints_all_send_procNo_glo[i]){
					// cout << 
					// connPoints_all_send_procNo_glo[i][iter] << " " <<
					// connPoints_all_send_procNo_loc[i][iter] << " " <<
					// connPoints_all_send_id_loc[i][iter] << " " << 
					// endl;
					// ++iter;
				// }
			// }
		// }
		
		
		
		
		
		
		vector<vector<int>> connPoints_all_send_id_glo(mesh.points.size());
		{
			vector<vector<int>> send_debug_procNo_loc(size);
			vector<vector<int>> send_debug_id_loc(size);
			vector<vector<double>> send_debug_8x(size);
			vector<vector<double>> send_debug_8y(size);
			vector<vector<double>> send_debug_8z(size);
			for(int i=0, ip=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				int iter=0;
				for(auto& proc : connPoints_all_send_procNo_glo[i]){
					int procNo_loc = connPoints_all_send_procNo_loc[i][iter];
					int id_loc = connPoints_all_send_id_loc[i][iter];
					
					send_debug_procNo_loc[proc].push_back(procNo_loc);
					send_debug_id_loc[proc].push_back(id_loc);
					
					send_debug_8x[proc].push_back(mesh.points[i].x);
					send_debug_8y[proc].push_back(mesh.points[i].y);
					send_debug_8z[proc].push_back(mesh.points[i].z);
					
					// if(rank==0 && proc==1){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
					
					++iter;
				}
				
			}
			vector<vector<int>> recv_debug_procNo_loc;
			MASCH_MPI_Alltoallv(send_debug_procNo_loc, recv_debug_procNo_loc);
			vector<vector<int>> recv_debug_id_loc;
			MASCH_MPI_Alltoallv(send_debug_id_loc, recv_debug_id_loc);
			// vector<vector<int>> recv_debug_id_glo;
			// MASCH_MPI_Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			vector<vector<double>> recv_debug_8x;
			MASCH_MPI_Alltoallv(send_debug_8x, recv_debug_8x);
			vector<vector<double>> recv_debug_8y;
			MASCH_MPI_Alltoallv(send_debug_8y, recv_debug_8y);
			vector<vector<double>> recv_debug_8z;
			MASCH_MPI_Alltoallv(send_debug_8z, recv_debug_8z);
			
			// cout << endl;
			// cout << endl;
			// cout << endl;
	// MPI_Barrier(MPI_COMM_WORLD);
				
			vector<vector<int>> send_debug_id_glo(size);
			for(int ip=0; ip<size; ++ip){
				for(int i=0; i<recv_debug_procNo_loc[ip].size(); ++i){
					int procNo_loc = recv_debug_procNo_loc[ip][i];
					int id_loc = recv_debug_id_loc[ip][i];
					
					if(rank==1 && ip==0){
						// if(rank==1){
						// cout <<
						// recv_debug_8x[ip][i] << " " <<
						// recv_debug_8y[ip][i] << " " <<
						// recv_debug_8z[ip][i] << " " <<
						// endl;
						
						// cout <<
						// recv_localPoint_xyz[procNo_loc][3*id_loc] << " " <<
						// recv_localPoint_xyz[procNo_loc][3*id_loc+1] << " " <<
						// recv_localPoint_xyz[procNo_loc][3*id_loc+2] << " " <<
						// endl;
					}
					
					int id_glo = points_id_local2global[procNo_loc][id_loc];
					send_debug_id_glo[ip].push_back(id_glo);
					
					auto& point = meshComb.points[str_points_glo[procNo_loc]+id_loc];
					// if(rank==1 && ip==0){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
					
				}
			}
			vector<vector<int>> recv_debug_id_glo;
			MASCH_MPI_Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			// cout << endl;
			// cout << endl;
			// cout << endl;
	// MPI_Barrier(MPI_COMM_WORLD);
			
			vector<int> iter_tmp_n(size,0);
			for(int i=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				int iter=0;
				for(auto& proc : connPoints_all_send_procNo_glo[i]){
					int procNo_loc = connPoints_all_send_procNo_loc[i][iter];
					int id_loc = connPoints_all_send_id_loc[i][iter];
					
					int id_glo = recv_debug_id_glo[proc].at(iter_tmp_n[proc]);
					connPoints_all_send_id_glo[i].push_back(id_glo);
					
					// if(rank==0 && proc==1){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
					
					++iter_tmp_n[proc];
					++iter;
				}
				
			}
			// for(int i=0, ip=0; i<mesh.points.size(); ++i){
				// auto& point = mesh.points[i];
				// int iter=0;
				// for(auto& proc : connPoints_all_send_procNo_glo[i]){
					// int procNo_loc = connPoints_all_send_procNo_loc[i][iter];
					// int id_loc = connPoints_all_send_id_loc[i][iter];
					
					// // send_debug_procNo_loc[proc].push_back(procNo_loc);
					// // send_debug_id_loc[proc].push_back(id_loc);
					
					// // send_debug_8x[proc].push_back(mesh.points[i].x);
					// // send_debug_8y[proc].push_back(mesh.points[i].y);
					// // send_debug_8z[proc].push_back(mesh.points[i].z);
					
					// if(rank==0 && proc==1){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
					
					// ++iter;
				// }
				
			// }
				
		}
		
	
		
		
		// // 디버그
		// {
			// vector<vector<int>> send_debug_tmp_size(size);
			// vector<vector<int>> send_debug_procNo_glo(size);
			// vector<vector<int>> send_debug_procNo_loc(size);
			// vector<vector<int>> send_debug_id_loc(size);
			// vector<vector<int>> send_debug_id_glo(size);
			// for(int i=0, ip=0; i<mesh.points.size(); ++i){
				// auto& point = mesh.points[i];
				// for(auto& proc : send_localPoint_proc[i]){
					// int tmp_size = connPoints_all_send_procNo_glo[i].size();
					// send_debug_tmp_size[proc].push_back(tmp_size);
					// for(int j=0; j<tmp_size; ++j){
						// send_debug_procNo_glo[proc].push_back(connPoints_all_send_procNo_glo[i].at(j));
						// send_debug_procNo_loc[proc].push_back(connPoints_all_send_procNo_loc[i].at(j));
						// send_debug_id_loc[proc].push_back(connPoints_all_send_id_loc[i].at(j));
						// send_debug_id_glo[proc].push_back(connPoints_all_send_id_glo[i].at(j));
						
						// // if(rank==0 && proc==1 && connPoints_all_send_procNo_glo[i].at(j) == 1){
							// // cout <<
							// // point.x << " " <<
							// // point.y << " " <<
							// // point.z << " " <<
							// // endl;
						// // }
					// }
				// }
			// }
			// vector<vector<int>> recv_debug_tmp_size;
			// MASCH_MPI_Alltoallv(send_debug_tmp_size, recv_debug_tmp_size);
			// vector<vector<int>> recv_debug_procNo_glo;
			// MASCH_MPI_Alltoallv(send_debug_procNo_glo, recv_debug_procNo_glo);
			// vector<vector<int>> recv_debug_procNo_loc;
			// MASCH_MPI_Alltoallv(send_debug_procNo_loc, recv_debug_procNo_loc);
			// vector<vector<int>> recv_debug_id_loc;
			// MASCH_MPI_Alltoallv(send_debug_id_loc, recv_debug_id_loc);
			// vector<vector<int>> recv_debug_id_glo;
			// MASCH_MPI_Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			// cout << endl;
			// cout << endl;
			// cout << endl;
			
			
			// for(int ip=0, iter_glob=0; ip<size; ++ip){
				// int tmp_size = recv_debug_tmp_size[ip].size();
				// for(int i=0, iter=0, iter2=0; i<tmp_size; ++i){
					// int tmp2_size = recv_debug_tmp_size[ip][i];
					// for(int j=0; j<tmp2_size; ++j){
						// int procNo_glo = recv_debug_procNo_glo[ip][iter];
						// int procNo_loc = recv_debug_procNo_loc[ip][iter];
						// int id_loc = recv_debug_id_loc[ip][iter];
						// int id_glo = recv_debug_id_glo[ip][iter];
						
						// // if(rank==1 && ip==0 && procNo_glo==1){
							// // cout << "AA = " <<
							// // meshComb.points[id_glo].x << " " <<
							// // meshComb.points[id_glo].y << " " <<
							// // meshComb.points[id_glo].z << " " <<
							// // endl;
							// // // cout << "BB = " <<
							// // // recv_localPoint_xyz[procNo_loc][3*id_loc] << " " <<
							// // // recv_localPoint_xyz[procNo_loc][3*id_loc+1] << " " <<
							// // // recv_localPoint_xyz[procNo_loc][3*id_loc+2] << " " <<
							// // // endl;
						// // }
						
						// if(rank!=procNo_glo){
							// auto& connPoints = meshComb.points.at(iter_glob).connPoints;
							// bool thereProc = false;
							// for(auto& [proc, id] : connPoints){
								// if(proc==procNo_glo) thereProc = true;
							// }
							// if(thereProc==false) 
								// connPoints.push_back(make_pair(procNo_glo,id_glo));
						// }
						
						
						
						// ++iter;
					// }
					
					
					// ++iter_glob;
					
				// }
			// }
			
			
			
			
		// }
		
		
			
			// // 포인트 삭제
			// {
				// int numN = 0;
				// meshComb.points.erase( std::remove_if( meshComb.points.begin(), meshComb.points.end(), 
					// [&deletePoints, &numN](MASCH_Point const& v) { 
					// return deletePoints[numN++]; 
					// }), meshComb.points.end());
			// }
			
		
		
		// // 디버그
		// {
			// vector<vector<int>> send_debug_tmp_size(size);
			// for(int i=0, ip=0; i<mesh.points.size(); ++i){
				// auto& point = mesh.points[i];
				// int iter=0;
				// for(auto& proc : connPoints_all_send_procNo_glo[i]){
					// // if(rank==0) cout << "proc1 = " << proc << endl;
					// if(rank==0 && proc==1){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
					// send_debug_tmp_size[proc].push_back(connPoints_all_send_id_glo[i][iter]);
					
					// ++iter;
				// }
					
				// // for(auto& proc : send_localPoint_proc[i]){
					// // if(rank==0) cout << "proc2 = " << proc << endl;
				// // }
			// }
			
			// vector<vector<int>> recv_debug_tmp_size;
			// MASCH_MPI_Alltoallv(send_debug_tmp_size, recv_debug_tmp_size);
			
			// cout << endl;
			// cout << endl;
			// cout << endl;
			
			
			// for(int ip=0, iter_glob=0; ip<size; ++ip){
				// int tmp_size = recv_debug_tmp_size[ip].size();
				// for(int i=0, iter=0, iter2=0; i<tmp_size; ++i){
					// int id_glo = recv_debug_tmp_size[ip][i];
					// auto& point = meshComb.points[id_glo];
					// if(rank==1 && ip==0){
						// cout << 
						// point.x << " " <<
						// point.y << " " <<
						// point.z << " " <<
						// endl;
					// }
				// }
			// }
		// }
		
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		
			
		// 포인트 삭제
		{
			int numN = 0;
			meshComb.points.erase( std::remove_if( meshComb.points.begin(), meshComb.points.end(), 
				[&deletePoints, &numN](MASCH_Point const& v) { 
				return deletePoints[numN++]; 
				}), meshComb.points.end());
		}
			
		

		
		{
			vector<vector<int>> send_debug_tmp_size(size);
			vector<vector<int>> send_debug_procNo_glo(size);
			vector<vector<int>> send_debug_procNo_loc(size);
			vector<vector<int>> send_debug_id_loc(size);
			vector<vector<int>> send_debug_id_glo(size);
			for(int i=0; i<mesh.points.size(); ++i){
				auto& point = mesh.points[i];
				for(auto& proc : send_localPoint_proc[i]){
					int tmp_size = connPoints_all_send_procNo_glo[i].size();
					send_debug_tmp_size[proc].push_back(tmp_size);
					for(int j=0; j<tmp_size; ++j){
						send_debug_procNo_glo[proc].push_back(connPoints_all_send_procNo_glo[i].at(j));
						send_debug_procNo_loc[proc].push_back(connPoints_all_send_procNo_loc[i].at(j));
						send_debug_id_loc[proc].push_back(connPoints_all_send_id_loc[i].at(j));
						send_debug_id_glo[proc].push_back(connPoints_all_send_id_glo[i].at(j));
						
					}
				}
			}
			vector<vector<int>> recv_debug_tmp_size;
			MASCH_MPI_Alltoallv(send_debug_tmp_size, recv_debug_tmp_size);
			vector<vector<int>> recv_debug_procNo_glo;
			MASCH_MPI_Alltoallv(send_debug_procNo_glo, recv_debug_procNo_glo);
			vector<vector<int>> recv_debug_procNo_loc;
			MASCH_MPI_Alltoallv(send_debug_procNo_loc, recv_debug_procNo_loc);
			vector<vector<int>> recv_debug_id_loc;
			MASCH_MPI_Alltoallv(send_debug_id_loc, recv_debug_id_loc);
			vector<vector<int>> recv_debug_id_glo;
			MASCH_MPI_Alltoallv(send_debug_id_glo, recv_debug_id_glo);
			
			
			
			for(int ip=0; ip<size; ++ip){
				int tmp_size = recv_debug_tmp_size[ip].size();
				for(int i=0, iter=0, iter2=0; i<tmp_size; ++i){
					
					int tmp2_size = recv_debug_tmp_size[ip][i];
					int llll=0;
					int iter_glob_tmp = -1;
					for(int j=0; j<tmp2_size; ++j){
						int procNo_glo = recv_debug_procNo_glo[ip][iter];
						int procNo_loc = recv_debug_procNo_loc[ip][iter];
						int id_loc = recv_debug_id_loc[ip][iter];
						int id_glo = recv_debug_id_glo[ip][iter];
						if(rank==procNo_glo){
							iter_glob_tmp = id_glo;
							// if(rank==0) cout << procNo_glo << " " << procNo_loc << " " << id_loc << " " << id_glo << " " << iter_glob << endl;
							// ++llll;
						}
						++iter;
					}
					for(int j=0; j<tmp2_size; ++j){
						int procNo_glo = recv_debug_procNo_glo[ip][iter2];
						int procNo_loc = recv_debug_procNo_loc[ip][iter2];
						int id_loc = recv_debug_id_loc[ip][iter2];
						int id_glo = recv_debug_id_glo[ip][iter2];
						if(rank!=procNo_glo){
							auto& connPoints = meshComb.points.at(iter_glob_tmp).connPoints;
							bool thereProc = false;
							for(auto& [proc, id] : connPoints){
								if(proc==procNo_glo) thereProc = true;
							}
							if(thereProc==false) 
								connPoints.push_back(make_pair(procNo_glo,id_glo));
						}
						++iter2;
					}
					
					
					
				}
				
			}
			
			
			
			
			
				
			
			
			// 디;버그
			{
		
				MPI_Barrier(MPI_COMM_WORLD);
				vector<vector<int>> send_test1(size);
				vector<vector<double>> send_test2(size);
				int iter=0;
				for(auto& point : meshComb.points){
					auto& connPoints = point.connPoints;
					for(auto& [proc, id] : connPoints){
						// if(rank==0){
							// cout << 
							// iter << " " <<
							// proc << " " <<
							// id << " " <<
							// endl;
						// }
						// if(rank==0 && proc==1){
							// cout << 
							// point.x << " " <<
							// point.y << " " <<
							// point.z << " " <<
							// endl;
						// }
						send_test1[proc].push_back(id);
						send_test2[proc].push_back(point.x);
						send_test2[proc].push_back(point.y);
						send_test2[proc].push_back(point.z);
					}
					++iter;
				}
				
				vector<vector<int>> recv_test1;
				MASCH_MPI_Alltoallv(send_test1, recv_test1);
				vector<vector<double>> recv_test2;
				MASCH_MPI_Alltoallv(send_test2, recv_test2);
				
				// cout << endl;
				// cout << endl;
				// cout << endl;
				
				for(int ip=0; ip<size; ++ip){
					int iter2 = 0;
					for(auto& ipoint : recv_test1[ip]){
						auto& point = meshComb.points[ipoint];
						
						// if(rank==1 && ip==0){
							// cout << 
							// point.x << " " <<
							// point.y << " " <<
							// point.z << " " <<
							// endl;
						// }
						
						// if(rank==0) cout << ipoint << endl;
						double resi = 0.0;
						resi += abs(point.x-recv_test2[ip][iter2++]);
						resi += abs(point.y-recv_test2[ip][iter2++]);
						resi += abs(point.z-recv_test2[ip][iter2++]);
						if(rank==0) {
							if(resi>1.e-16){
								cout << "NONON " << resi << endl;
							}
						}
					}
				}
				
				
			}
			
		
			
			
			// {
				// // conn 융합
				// for(int ip=0; ip<size; ++ip){
					// int str = str_points_glo[ip];
					// int end = str_points_glo[ip+1];
					// for(int i=str; i<end; ++i){
						// int id_loc = i-str;
						// auto& overlab_connPoints = meshComb.points.at(i).connPoints;
						
						// if(deletePoints[i]==true){
							// int original_id_glo = points_id_local2global[ip][id_loc];
							// auto& original_connPoints = meshComb.points.at(original_id_glo).connPoints;
							
							// for(auto& [over_proc, over_id] : overlab_connPoints){
								// bool thereProc = false;
								// for(auto& [proc, id] : original_connPoints){
									// if(proc==over_proc) thereProc = true;
								// }
								// if(thereProc==false) {
									// original_connPoints.push_back(make_pair(over_proc,over_id));
								// }
							// }
						// }
					// }
				// }
						
			// }
		
		
			
			// // 포인트 삭제
			// {
				// int numN = 0;
				// meshComb.points.erase( std::remove_if( meshComb.points.begin(), meshComb.points.end(), 
					// [&deletePoints, &numN](MASCH_Point const& v) { 
					// return deletePoints[numN++]; 
					// }), meshComb.points.end());
			// }
			
				
			
			
			
		}

		
			
		
		// {
			// int iter = 0;
			// for(auto& point : meshComb.points){
				
				// auto& connPoints = point.connPoints;
				
				// for(auto& [proc, id] : connPoints){
					
					// // if(rank==1) cout << iter << " " << proc << " " << id << endl;
				// }
				// ++iter;
			// }
		// }
	
	
		
		// {
			
			// vector<vector<double>> send_test1(size);
			// int iter = 0;
			// for(auto& point : meshComb.points){
				
				// auto& connPoints = point.connPoints;
				// if(connPoints.size()==0){
					// if(rank==0){
						// send_test1[1].push_back(point.x);
						// send_test1[1].push_back(point.y);
						// send_test1[1].push_back(point.z);
					// }
					// else{
						// send_test1[0].push_back(point.x);
						// send_test1[0].push_back(point.y);
						// send_test1[0].push_back(point.z);
					// }
				// }
			// }
			// vector<vector<double>> recv_test1;
			// MASCH_MPI_Alltoallv(send_test1, recv_test1);
			
			// for(int ip=0, iter2=0; ip<size; ++ip){
				// int tmp_size = recv_test1[ip].size();
				// for(int i=0, iter3=0; i<tmp_size/3; ++i){
					// double x_tmp = recv_test1[ip][iter3++];
					// double y_tmp = recv_test1[ip][iter3++];
					// double z_tmp = recv_test1[ip][iter3++];
					
					// for(auto& point : meshComb.points){
						// double resi = 0.0;
						// resi += abs(point.x-x_tmp);
						// resi += abs(point.y-y_tmp);
						// resi += abs(point.z-z_tmp);
						
						// if(resi<1.e-8){
							// cout << resi << endl;
						// }
					// }
					
					
				// }
				
			// }
			
		// }
			
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// // 디;버그
		// {
			// vector<vector<int>> send_test1(size);
			// vector<vector<double>> send_test2(size);
			// int iter=0;
			// for(auto& point : meshComb.points){
				// auto& connPoints = point.connPoints;
				
				// // if(rank==0) cout << iter << endl;
				// // if(rank==0) cout << point.x << " " << point.y << " " << point.z << endl;
				// for(auto& [proc, id] : connPoints){
					// send_test1[proc].push_back(id);
					// send_test2[proc].push_back(point.x);
					// send_test2[proc].push_back(point.y);
					// send_test2[proc].push_back(point.z);
				// }
				// ++iter;
			// }
			
			// vector<vector<int>> recv_test1;
			// MASCH_MPI_Alltoallv(send_test1, recv_test1);
			// vector<vector<double>> recv_test2;
			// MASCH_MPI_Alltoallv(send_test2, recv_test2);
			
			// for(int ip=0, iter2=0; ip<size; ++ip){
				// for(auto& ipoint : recv_test1[ip]){
					// auto& point = meshComb.points[ipoint];
					// // if(rank==0) cout << ipoint << endl;
					// double resi = 0.0;
					// resi += abs(point.x-recv_test2[ip][iter2++]);
					// resi += abs(point.y-recv_test2[ip][iter2++]);
					// resi += abs(point.z-recv_test2[ip][iter2++]);
					// if(rank==0) {
						// if(resi>1.e-8){
							// cout << "NONON " << resi << endl;
						// }
					// }
				// }
			// }
			
			
		// }
		
		
		
		//======================================
		

		
		vector<pair<int,int>> send_localCell_proc_id;
		vector<vector<int>> send_localCell_id(size);
		vector<int> send_localCell_n(size,0);
		for(int i=0, ip=0; i<mesh.cells.size(); ++i){
			int proc = idBlockCell[i];
			int tmp_nCell = send_localCell_n.at(proc)++;
			send_localCell_proc_id.push_back(make_pair(proc,tmp_nCell));
			send_localCell_id[proc].push_back(tmp_nCell);
		}
		
		
		vector<vector<int>> recv_localCell_id;
		MASCH_MPI_Alltoallv(send_localCell_id, recv_localCell_id);
		vector<int> nCells_local(size+1,0);
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = recv_localCell_id[ip].size();
		}
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = nCells_local[ip+1] + nCells_local[ip];
		}
		
		
		
		
		
		vector<int> boundary_type(mesh.faces.size(),0);
		int nbc=0;
		for(int i=0; i<mesh.boundaries.size(); ++i){
			auto& boundary = mesh.boundaries[i];
			if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
				++nbc;
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int j=str; j<end; ++j){
					boundary_type[j] = i;
				}
			}
		}
		
		
		
		// processor 페이스 고유 번호 부여
		vector<int> proc_numbering;
		{
			int maxProcFaceNum = 0;
			for(int ip=rank; ip<size; ++ip){
				maxProcFaceNum += mesh.countsProcFaces[ip];
			}
			// if(rank==0) cout << maxProcFaceNum << endl;
			vector<int> tmp_maxProcFaceNum(size);
			MPI_Allgather(&maxProcFaceNum, 1, MPI_INT, tmp_maxProcFaceNum.data(), 1, MPI_INT, MPI_COMM_WORLD);
			
			vector<int> str_proc_numbering(size+1,0);
			for(int ip=0; ip<size; ++ip){
				str_proc_numbering[ip+1] = str_proc_numbering[ip] + tmp_maxProcFaceNum[ip];
			}
			
			
			// vector<int> send_proc_numbering;
			// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				// if(face.getType()==MASCH_Face_Types::PROCESSOR){
					// int tmp_num = str_proc_numbering[rank] + ip;
					// send_proc_numbering.push_back(tmp_num);
					// ++ip;
				// }
			// }
			
			
			
			
			vector<int> send_proc_numbering;
			for(int ip=0, numbering=0; ip<rank; ++ip){
				int str=mesh.displsProcFaces[ip];
				int end=str+mesh.countsProcFaces[ip];
				for(int i=str; i<end; ++i){
					send_proc_numbering.push_back(-1);
				}
			}
			for(int ip=rank, numbering=0; ip<size; ++ip){
				int str=mesh.displsProcFaces[ip];
				int end=str+mesh.countsProcFaces[ip];
				for(int i=str; i<end; ++i){
					send_proc_numbering.push_back(str_proc_numbering[rank] + i - mesh.displsProcFaces[rank]);
				}
			}
			vector<int> recv_proc_numbering;
			recv_proc_numbering.resize(send_proc_numbering.size());
			MPI_Alltoallv( send_proc_numbering.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
						   recv_proc_numbering.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
						   MPI_COMM_WORLD);
				
			
			for(int ip=0, numbering=0; ip<rank; ++ip){
				int str=mesh.displsProcFaces[ip];
				int end=str+mesh.countsProcFaces[ip];
				for(int i=str; i<end; ++i){
					proc_numbering.push_back(recv_proc_numbering[i]);
				}
			}
			for(int ip=rank, numbering=0; ip<size; ++ip){
				int str=mesh.displsProcFaces[ip];
				int end=str+mesh.countsProcFaces[ip];
				for(int i=str; i<end; ++i){
					proc_numbering.push_back(send_proc_numbering[i]);
				}
			}
			 


		
	// MPI_Barrier(MPI_COMM_WORLD);
			// for(int ip=0, numbering=0; ip<size; ++ip){
				// int str=mesh.displsProcFaces[ip];
				// int end=str+mesh.countsProcFaces[ip];
				// for(int i=str; i<end; ++i){
					// // if(rank==1) cout << rank << " " << proc_numbering[i] << endl;
					// cout << rank << " " << proc_numbering[i] << endl;
				// }
			// }
					
			// // mesh.displsProcFaces[rank]
		}
	
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
		
		
		
		
	
		
		vector<vector<int>> send_localFace_IN2IN_iL(size);
		vector<vector<int>> send_localFace_IN2IN_iR(size);
		vector<vector<int>> send_localFace_IN2IN_ipoints(size);
		
		vector<vector<int>> send_localFace_IN2PR_iL(size);
		vector<vector<int>> send_globalFace_IN2PR_toProc(size);
		vector<vector<int>> send_localFace_IN2PR_ipoints(size);
		
		vector<vector<int>> send_localFace_BC2BC_iL(size);
		vector<vector<int>> send_localFace_BC2BC_BCType(size);
		vector<vector<int>> send_localFace_BC2BC_ipoints(size);
		
		vector<vector<int>> send_localFace_PR2IN_iL(size);
		vector<vector<int>> send_localFace_PR2IN_toProc(size);
		vector<vector<int>> send_localFace_PR2IN_ipoints(size);
		
		vector<vector<int>> send_localFace_PR2PR_procFaceNum(size);
		vector<vector<int>> send_localFace_PR2PR_iL(size);
		vector<vector<int>> send_globalFace_PR2PR_toProc(size);
		vector<vector<int>> send_localFace_PR2PR_ipoints(size);
		
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL;
			int iR = face.iR;
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				auto& [ipR, iR_local] = send_localCell_proc_id[iR];
				if(ipL==ipR){
					
					send_localFace_IN2IN_iL[ipL].push_back(iL_local);
					send_localFace_IN2IN_iR[ipL].push_back(iR_local);
					send_localFace_IN2IN_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l = 0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) ipoint_local = ipoint;
							if(proc==ipL) ++l;
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						//if(l>1) cout << "ERROR2" << endl;
						send_localFace_IN2IN_ipoints[ipL].push_back(ipoint_local);
					}
					
				}
				else{
					// cout << "AAAAAAAAAA" << endl;
					// cout << ipL << " " << ipR << endl;
					
					send_localFace_IN2PR_iL[ipL].push_back(iL_local);
					send_globalFace_IN2PR_toProc[ipL].push_back(ipR);
					send_localFace_IN2PR_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 53" << endl;
						send_localFace_IN2PR_ipoints[ipL].push_back(ipoint_local);
					}
					
					send_localFace_IN2PR_iL[ipR].push_back(iR_local);
					send_globalFace_IN2PR_toProc[ipR].push_back(ipL);
					vector<int> tmp_ipoints = face.ipoints;
					std::reverse(tmp_ipoints.begin()+1, tmp_ipoints.end());
					send_localFace_IN2PR_ipoints[ipR].push_back(face.ipoints.size());
					for(auto& item : tmp_ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipR) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 54" << endl;
						send_localFace_IN2PR_ipoints[ipR].push_back(ipoint_local);
					}
				}
			}
			else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				
				send_localFace_BC2BC_iL[ipL].push_back(iL_local);
				send_localFace_BC2BC_BCType[ipL].push_back(boundary_type[i]);
				send_localFace_BC2BC_ipoints[ipL].push_back(face.ipoints.size());
				for(auto& item : face.ipoints) {
					int ipoint_local = -1;
					for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
						if(proc==ipL) ipoint_local = ipoint;
					}
					if(ipoint_local==-1) cout << "ERROR" << endl;
					send_localFace_BC2BC_ipoints[ipL].push_back(ipoint_local);
				}
				
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				auto& [ipL, iL_local] = send_localCell_proc_id[iL];
				int ipR = recv_idBlockCell[ip];
				int rightProcNo = recv_rank[ip];
				
				// // if(rank==0) cout << ipR << " " << rightProcNo << endl;
				// if(rank==0){
					// if(ipL==ipR){
						// // for(auto& iface : face.ifaces)
						// {
							// // auto& face = mesh.faces[iface];
							// // cout << face_state[iface] << endl;
							// cout << "(" << mesh.points[face.ipoints[0]].x << " " <<
							// mesh.points[face.ipoints[0]].y << " " <<
							// mesh.points[face.ipoints[0]].z << ") ";
							// cout << "(" << mesh.points[face.ipoints[1]].x << " " <<
							// mesh.points[face.ipoints[1]].y << " " <<
							// mesh.points[face.ipoints[1]].z << ") ";
							// cout << "(" << mesh.points[face.ipoints[2]].x << " " <<
							// mesh.points[face.ipoints[2]].y << " " <<
							// mesh.points[face.ipoints[2]].z << ") ";
							// cout << "(" << mesh.points[face.ipoints[3]].x << " " <<
							// mesh.points[face.ipoints[3]].y << " " <<
							// mesh.points[face.ipoints[3]].z << ") " << endl;
						// }
						
					// }
				// }
				
				if(ipL==ipR){
					
					send_localFace_PR2IN_iL[ipL].push_back(iL_local);
					send_localFace_PR2IN_toProc[ipL].push_back(rightProcNo);
					send_localFace_PR2IN_ipoints[ipL].push_back(face.ipoints.size());
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) ipoint_local = ipoint;
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						send_localFace_PR2IN_ipoints[ipL].push_back(ipoint_local);
					}
				}
				else{
					
					send_localFace_PR2PR_procFaceNum[ipL].push_back(proc_numbering[ip]);
					
					send_localFace_PR2PR_iL[ipL].push_back(iL_local);
					send_globalFace_PR2PR_toProc[ipL].push_back(ipR);
					send_localFace_PR2PR_ipoints[ipL].push_back(face.ipoints.size());						
					double avgx = 0.0;
					double avgy = 0.0;
					double avgz = 0.0;
					int tmp_size = face.ipoints.size();
					for(auto& item : face.ipoints) {
						int ipoint_local = -1;
						int l=0;
						for(auto& [proc, ipoint] : send_localPoint_proc_id[item]){
							if(proc==ipL) {
								ipoint_local = ipoint;
								++l;
							}
						}
						if(ipoint_local==-1) cout << "ERROR" << endl;
						if(l>1) cout << "ERROR 55" << endl;
						// cout << l << endl;
						send_localFace_PR2PR_ipoints[ipL].push_back(ipoint_local);
						
						avgx += mesh.points[item].x/(double)tmp_size;
						avgy += mesh.points[item].y/(double)tmp_size;
						avgz += mesh.points[item].z/(double)tmp_size;
					}
					
					// if(ipL==0 && ipR==1){
						// cout << rank << " " << proc_numbering[ip] << endl;
						// cout << avgx << " " << avgy << " " << avgz << endl;
						
					// }
					// if(ipL==1 && ipR==0){
						// cout << rank << " " << proc_numbering[ip] << endl;
						// cout << avgx << " " << avgy << " " << avgz << endl;
						
					// }
				}
				++ip;
			}
		}
		
		
		
		vector<vector<int>> recv_localFace_IN2IN_iL;
		vector<vector<int>> recv_localFace_IN2IN_iR;
		vector<vector<int>> recv_localFace_IN2IN_ipoints;
		vector<vector<int>> recv_localFace_IN2PR_iL;
		vector<vector<int>> recv_globalFace_IN2PR_toProc;
		vector<vector<int>> recv_localFace_IN2PR_ipoints;
		vector<vector<int>> recv_localFace_BC2BC_iL;
		vector<vector<int>> recv_localFace_BC2BC_BCType;
		vector<vector<int>> recv_localFace_BC2BC_ipoints;
		vector<vector<int>> recv_localFace_PR2IN_iL;
		vector<vector<int>> recv_localFace_PR2IN_toProc;
		vector<vector<int>> recv_localFace_PR2IN_ipoints;
		vector<vector<int>> recv_localFace_PR2PR_procFaceNum;
		vector<vector<int>> recv_localFace_PR2PR_iL;
		vector<vector<int>> recv_globalFace_PR2PR_toProc;
		vector<vector<int>> recv_localFace_PR2PR_ipoints;
		
		
		MASCH_MPI_Alltoallv(send_localFace_IN2IN_iL, recv_localFace_IN2IN_iL);
		MASCH_MPI_Alltoallv(send_localFace_IN2IN_iR, recv_localFace_IN2IN_iR);
		MASCH_MPI_Alltoallv(send_localFace_IN2IN_ipoints, recv_localFace_IN2IN_ipoints);
		
		MASCH_MPI_Alltoallv(send_localFace_IN2PR_iL, recv_localFace_IN2PR_iL);
		MASCH_MPI_Alltoallv(send_globalFace_IN2PR_toProc, recv_globalFace_IN2PR_toProc);
		MASCH_MPI_Alltoallv(send_localFace_IN2PR_ipoints, recv_localFace_IN2PR_ipoints);
		
		MASCH_MPI_Alltoallv(send_localFace_BC2BC_iL, recv_localFace_BC2BC_iL);
		MASCH_MPI_Alltoallv(send_localFace_BC2BC_BCType, recv_localFace_BC2BC_BCType);
		MASCH_MPI_Alltoallv(send_localFace_BC2BC_ipoints, recv_localFace_BC2BC_ipoints);
		
		MASCH_MPI_Alltoallv(send_localFace_PR2IN_iL, recv_localFace_PR2IN_iL);
		MASCH_MPI_Alltoallv(send_localFace_PR2IN_toProc, recv_localFace_PR2IN_toProc);
		MASCH_MPI_Alltoallv(send_localFace_PR2IN_ipoints, recv_localFace_PR2IN_ipoints);
		
		MASCH_MPI_Alltoallv(send_localFace_PR2PR_procFaceNum, recv_localFace_PR2PR_procFaceNum);
		MASCH_MPI_Alltoallv(send_localFace_PR2PR_iL, recv_localFace_PR2PR_iL);
		MASCH_MPI_Alltoallv(send_globalFace_PR2PR_toProc, recv_globalFace_PR2PR_toProc);
		MASCH_MPI_Alltoallv(send_localFace_PR2PR_ipoints, recv_localFace_PR2PR_ipoints);
		
		
		
		
		
		
		// IN2IN
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_localFace_IN2IN_iL[ip].size(); ++i){
				
				face_state.push_back(meshComb.faces.size());
				
				meshComb.addFace();
				meshComb.faces.back().setType(MASCH_Face_Types::INTERNAL);
				{
					int id_local = recv_localFace_IN2IN_iL[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iL = id_global;
				}
				{
					int id_local = recv_localFace_IN2IN_iR[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iR = id_global;
				}
				int tmp_size = recv_localFace_IN2IN_ipoints[ip].at(iter++);
				for(int j=0; j<tmp_size; ++j){
					int ipoint_local = recv_localFace_IN2IN_ipoints[ip].at(iter++);
					int ipoint_global = points_id_local2global[ip].at(ipoint_local);
					meshComb.faces.back().ipoints.push_back(ipoint_global);
				}
			}
		}
		
		
		// PR2IN
		{
			int PR2IN_str = meshComb.faces.size();
			vector<int> tmp_nFaces_local(size,0);
			vector<int> tmp_strFaces_local(size+1,0);
			for(int ip=0; ip<size; ++ip){
				tmp_strFaces_local[ip+1] = tmp_strFaces_local[ip];
				for(int i=0, iter=0; i<recv_localFace_PR2IN_iL[ip].size(); ++i){
					int rightProc_local = recv_localFace_PR2IN_toProc[ip].at(i);
					// 중복
					if(rightProc_local<ip){
						int tmp_id_global = PR2IN_str + 
											tmp_strFaces_local[rightProc_local] +
											tmp_nFaces_local.at(rightProc_local);
						++tmp_nFaces_local.at(rightProc_local);
						// ++tmp_nFaces_local.at(ip);
						auto& rightFace = meshComb.faces.at(tmp_id_global);
						int id_local = recv_localFace_PR2IN_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						// if(rightFace.iR!=-1) cout << rightFace.iR << " " << "ERROR" << endl;
						rightFace.iR = id_global;
						
						int tmp_size = recv_localFace_PR2IN_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j) iter++;
					}
					else{
				
				face_state.push_back(PR2IN);
						
						++tmp_strFaces_local[ip+1];
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::INTERNAL);
						int id_local = recv_localFace_PR2IN_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						meshComb.faces.back().iL = id_global;
						int tmp_size = recv_localFace_PR2IN_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j){
							int ipoint = recv_localFace_PR2IN_ipoints[ip].at(iter++);
							int id_global = points_id_local2global[ip].at(ipoint);
							meshComb.faces.back().ipoints.push_back(id_global);
						}
					}
				}
			}
		}
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
	
	
		
		// BC2BC
		vector<int> nFaces_boundary(nbc,0);
		
		vector<vector<pair<int,int>>> reorder_procFace_BC2BC_proc_id(nbc);
		vector<vector<int>> reorder_procFace_BC2BC_strId(nbc);
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_localFace_BC2BC_BCType[ip].size(); ++i){
				int ibc = recv_localFace_BC2BC_BCType[ip][i];
				reorder_procFace_BC2BC_proc_id[ibc].push_back(make_pair(ip,i));
				reorder_procFace_BC2BC_strId[ibc].push_back(iter);
				int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(iter++);
				for(int j=0; j<tmp_size; ++j) iter++;
			}
		}
			
		vector<int> iter_BC2BC(size,0);
		for(int ibc=0; ibc<nbc; ++ibc){
			int str_size = meshComb.faces.size();
			int iter=0;
			for(auto& [ip, i] : reorder_procFace_BC2BC_proc_id[ibc]){
				
				face_state.push_back(BC2BC);
				
				meshComb.addFace();
				meshComb.faces.back().setType(MASCH_Face_Types::BOUNDARY);
				{
					int id_local = recv_localFace_BC2BC_iL[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iL = id_global;
				}
				int strId = reorder_procFace_BC2BC_strId[ibc][iter];
				int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(strId++);
				for(int j=0; j<tmp_size; ++j){
					int ipoint = recv_localFace_BC2BC_ipoints[ip].at(strId++);
					int id_global = points_id_local2global[ip].at(ipoint);
					meshComb.faces.back().ipoints.push_back(id_global);
				}
						
				// int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
				// for(int j=0; j<tmp_size; ++j){
					// int ipoint = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
					// int id_global = points_id_local2global[ip].at(ipoint);
					// meshComb.faces.back().ipoints.push_back(id_global);
				// }
				++iter;
			}
			nFaces_boundary[ibc] = meshComb.faces.size() - str_size;
		}
	
		
		// 프로세서
		vector<int> nFaces_processor(size,0);
		{
			// IN2PR
			vector<vector<pair<int,int>>> reorder_procFace_IN2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_IN2PR_strId(size);
			// for(int ip=size-1; ip>=0; --ip){
			for(int ip=0; ip<size; ++ip){
				for(int i=0, iter=0; i<recv_localFace_IN2PR_iL[ip].size(); ++i){
					int rightProc_global = recv_globalFace_IN2PR_toProc[ip][i];
					// if(rightProc_global<rank)
					{
						reorder_procFace_IN2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
						reorder_procFace_IN2PR_strId[rightProc_global].push_back(iter);
						int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(iter++);
						for(int j=0; j<tmp_size; ++j) iter++;
					}
				}
			}
			
			// PR2PR
			vector<vector<int>> reorder_procFace_PR2PR_procFaceNum(size);
			vector<vector<pair<int,int>>> reorder_procFace_PR2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_PR2PR_strId(size);
			vector<vector<vector<int>>> reorder_procFace_PR2PR_ipoints(size);
			for(int ip=0; ip<size; ++ip){
				for(int i=0, iter=0, iter2=0; i<recv_localFace_PR2PR_iL[ip].size(); ++i){
					// if(rank==0) cout << recv_localFace_PR2PR_iL[ip][i] << " " << i << endl;
					int procFaceNum = recv_localFace_PR2PR_procFaceNum[ip].at(i);
					
					int rightProc_global = recv_globalFace_PR2PR_toProc[ip].at(i);
					// if(rightProc_global!=rank)
					{
						reorder_procFace_PR2PR_procFaceNum[rightProc_global].push_back(procFaceNum);
						reorder_procFace_PR2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
						reorder_procFace_PR2PR_strId[rightProc_global].push_back(iter);
						int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(iter++);
						double avgx = 0.0;
						double avgy = 0.0;
						double avgz = 0.0;
						vector<int> tmp;
						for(int j=0; j<tmp_size; ++j) {
							int id_local = recv_localFace_PR2PR_ipoints[ip].at(iter);
							int id_global = points_id_local2global[ip].at(id_local);
							tmp.push_back(id_global);
							avgx += meshComb.points[id_global].x/(double)tmp_size;
							avgy += meshComb.points[id_global].y/(double)tmp_size;
							avgz += meshComb.points[id_global].z/(double)tmp_size;
							iter++;
						}
						reorder_procFace_PR2PR_ipoints[rightProc_global].push_back(tmp);
						
						// if(rank==0 && rightProc_global==1){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
						// if(rank==1 && rightProc_global==0){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
					}
				}
			}
			
	
			
			// face 번호로 리오더
			for(int ip=0; ip<size; ++ip){
				int tmp_size = reorder_procFace_PR2PR_procFaceNum[ip].size();
				vector<int> v_procFaceNum = reorder_procFace_PR2PR_procFaceNum[ip];
				sort(v_procFaceNum.begin(), v_procFaceNum.end());
				// reverse(v_procFaceNum.begin(), v_procFaceNum.end());
				
				vector<pair<int,int>> tmp_v1(tmp_size);
				vector<int> tmp_v2(tmp_size);
				vector<vector<int>> tmp_v3(tmp_size);
				for(int i=0; i<tmp_size; ++i){
					int procFaceNum = reorder_procFace_PR2PR_procFaceNum[ip].at(i);
					// if(rank==1) cout << rank << " " << ip << " " << v_procFaceNum[i] << endl;
					// if(rank==1) cout << ip << " " << v_procFaceNum[i] << " " << procFaceNum << endl;
					auto it = find (v_procFaceNum.begin(), v_procFaceNum.end(), procFaceNum);
					
					int order = (it - v_procFaceNum.begin());
					// if(rank==1) cout << ip << " " << order << endl;
					
					tmp_v1[order] = (reorder_procFace_PR2PR_proc_id[ip][i]);
					tmp_v2[order] = (reorder_procFace_PR2PR_strId[ip][i]);
					tmp_v3[order] = (reorder_procFace_PR2PR_ipoints[ip][i]);
					
				}
				
				reorder_procFace_PR2PR_procFaceNum[ip].clear();
				reorder_procFace_PR2PR_proc_id[ip].clear();
				reorder_procFace_PR2PR_strId[ip].clear();
				reorder_procFace_PR2PR_ipoints[ip].clear();
				for(int i=0; i<tmp_size; ++i){
					reorder_procFace_PR2PR_proc_id[ip].push_back(tmp_v1[i]);
					reorder_procFace_PR2PR_strId[ip].push_back(tmp_v2[i]);
					reorder_procFace_PR2PR_procFaceNum[ip].push_back(v_procFaceNum[i]);
					reorder_procFace_PR2PR_ipoints[ip].push_back(tmp_v3[i]);
				}
				
				
			}
			
			
			// // 디버그
			// for(int proc=0; proc<size; ++proc){
				// int iter=0;
				// int tmp2_size = reorder_procFace_PR2PR_proc_id[proc].size();
				// for(int ii=0; ii<tmp2_size; ++ii){
				// // for(auto& [ip, i] : reorder_procFace_PR2PR_proc_id[proc]){
					// auto& [ip, i] = reorder_procFace_PR2PR_proc_id[proc][ii];
					// int strId = reorder_procFace_PR2PR_strId[proc][ii];
					// int procFaceNum = reorder_procFace_PR2PR_procFaceNum[proc][ii];
					// vector<int> ipoints = reorder_procFace_PR2PR_ipoints[proc][ii];
			
					// double avgx = 0.0;
					// double avgy = 0.0;
					// double avgz = 0.0;
					// // int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
					// int tmp_size = ipoints.size();
					// for(auto& ipoint : ipoints){
						// // int id_loc = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						// // int id_glo = points_id_local2global[ip].at(id_loc);
						// avgx += meshComb.points[ipoint].x/(double)tmp_size;
						// avgy += meshComb.points[ipoint].y/(double)tmp_size;
						// avgz += meshComb.points[ipoint].z/(double)tmp_size;
					// }
					
					// if(rank==0 && proc==1){
						// cout << rank << " " << procFaceNum << endl;
						// cout << avgx << " " << avgy << " " << avgz << endl;
						
					// }
					// if(rank==1 && proc==0){
						// cout << rank << " " << procFaceNum << endl;
						// cout << avgx << " " << avgy << " " << avgz << endl;
						
					// }
				// }
			// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
			
			
		
			vector<int> iter_IN2PR(size,0);
			vector<int> iter_PR2PR(size,0);
			for(int proc=0; proc<size; ++proc){
				int str_size = meshComb.faces.size();
				// IN2PR
				{
					int iter=0;
					// cout << reorder_procFace_IN2PR_proc_id[proc].size() << 
					for(auto& [ip, i] : reorder_procFace_IN2PR_proc_id[proc]){
				
				face_state.push_back(IN2PR);
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						{
							int id_local = recv_localFace_IN2PR_iL[ip].at(i);
							int id_global = nCells_local[ip] + id_local;
							meshComb.faces.back().iL = id_global;
						}
						// cout << recv_localFace_IN2PR_ipoints[ip].size() << " " << iter_IN2PR[ip] << endl;
						// int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(iter_IN2PR[ip]++);
						int strId = reorder_procFace_IN2PR_strId[proc].at(iter);
						int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(strId++);
						for(int j=0; j<tmp_size; ++j){
							int ipoint = recv_localFace_IN2PR_ipoints[ip].at(strId++);
							int id_global = points_id_local2global[ip].at(ipoint);
							meshComb.faces.back().ipoints.push_back(id_global);
						}
						++iter;
					}
				}
				// PR2PR
				{
					
					int iter=0;
					int tmp2_size = reorder_procFace_PR2PR_proc_id[proc].size();
					for(int ii=0; ii<tmp2_size; ++ii){
						
				face_state.push_back(PR2PR);
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						
						auto& [ip, i] = reorder_procFace_PR2PR_proc_id[proc][ii];
						int id_local = recv_localFace_PR2PR_iL[ip].at(i);
						int id_global = nCells_local[ip] + id_local;
						meshComb.faces.back().iL = id_global;
						
						int strId = reorder_procFace_PR2PR_strId[proc][ii];
						int procFaceNum = reorder_procFace_PR2PR_procFaceNum[proc][ii];
						vector<int> ipoints = reorder_procFace_PR2PR_ipoints[proc][ii];
				
						double avgx = 0.0;
						double avgy = 0.0;
						double avgz = 0.0;
						// int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						int tmp_size = ipoints.size();
						for(auto& ipoint : ipoints){
							// int id_loc = recv_localFace_PR2PR_ipoints[ip].at(strId++);
							// int id_glo = points_id_local2global[ip].at(id_loc);
							meshComb.faces.back().ipoints.push_back(ipoint);
							avgx += meshComb.points[ipoint].x/(double)tmp_size;
							avgy += meshComb.points[ipoint].y/(double)tmp_size;
							avgz += meshComb.points[ipoint].z/(double)tmp_size;
						}
						
						// if(rank==0 && proc==1){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
						// if(rank==1 && proc==0){
							// cout << rank << " " << procFaceNum << endl;
							// cout << avgx << " " << avgy << " " << avgz << endl;
							
						// }
					}
					
					
					// int iter=0;
					// for(auto& [ip, i] : reorder_procFace_PR2PR_proc_id[proc]){
				
				// face_state.push_back(PR2PR);
				
						// meshComb.addFace();
						// meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						// {
							// int id_local = recv_localFace_PR2PR_iL[ip].at(i);
							// int id_global = nCells_local[ip] + id_local;
							// meshComb.faces.back().iL = id_global;
						// }
						// int procFaceNum = reorder_procFace_PR2PR_procFaceNum[proc].at(iter);
						// double avgx = 0.0;
						// double avgy = 0.0;
						// double avgz = 0.0;
						// int strId = reorder_procFace_PR2PR_strId[proc].at(iter);
						// int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						// for(int j=0; j<tmp_size; ++j){
							// int ipoint = recv_localFace_PR2PR_ipoints[ip].at(strId++);
							// int id_global = points_id_local2global[ip].at(ipoint);
							// meshComb.faces.back().ipoints.push_back(id_global);
							// avgx += meshComb.points[id_global].x/(double)tmp_size;
							// avgy += meshComb.points[id_global].y/(double)tmp_size;
							// avgz += meshComb.points[id_global].z/(double)tmp_size;
						// }
						// ++iter;
						
						
						
						// // if(rank==0 && proc==1){
							// // cout << rank << " " << procFaceNum << endl;
							// // cout << avgx << " " << avgy << " " << avgz << endl;
							
						// // }
						// // if(rank==1 && proc==0){
							// // cout << rank << " " << procFaceNum << endl;
							// // cout << avgx << " " << avgy << " " << avgz << endl;
							
						// // }
					// }
				}
				nFaces_processor[proc] = meshComb.faces.size() - str_size;
			}
			
		}
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// 바운더리
		meshComb.boundaries.clear();
		for(int ibc=0; ibc<nbc; ++ibc){
			auto& boundary = mesh.boundaries[ibc];
			meshComb.addBoundary();
			meshComb.boundaries.back().name = boundary.name;
			meshComb.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			meshComb.boundaries.back().nFaces = nFaces_boundary[ibc];
		}
		
		for(int ip=0; ip<size; ++ip){
			if(nFaces_processor[ip]>0){
				string bcnames = "procBoundary" + to_string(rank) + "to" + to_string(ip);
				
				meshComb.addBoundary();
				meshComb.boundaries.back().name = bcnames;
				meshComb.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
				meshComb.boundaries.back().nFaces = nFaces_processor[ip];
				meshComb.boundaries.back().myProcNo = rank;
				meshComb.boundaries.back().rightProcNo = ip;
			}
		}
		
	
		int maxBCNum = meshComb.boundaries.size()-1;
		meshComb.boundaries[maxBCNum].startFace = meshComb.faces.size()-meshComb.boundaries[maxBCNum].nFaces;
		for(int i=maxBCNum-1; i>=0; --i){
			meshComb.boundaries[i].startFace = meshComb.boundaries[i+1].startFace-meshComb.boundaries[i].nFaces;
		}
		 
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		// // connPoints (접촉하는 포인트 옆 포인트가 뭔지 저장)
		// vector<vector<int>> send_rightProcNo_rightId_loc(size);
		// for(int ip=0; ip<size; ++ip){
			// // cout << recv_globalPoint_toProc_toId[ip].size()/4 << endl;
			// for(int i=0, iter=0; i<recv_globalPoint_toProc_toId[ip].size()/4; ++i){
				// int my_id_loc = recv_globalPoint_toProc_toId[ip][iter++]; // 자신의 local id  
				// int rightProcNo_glo = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 global proc 넘버
				// int rightProcNo_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록 에서 자신과 중복한 포인트의 local proc
				// int rightId_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 자신과 중복한 포인트의 local id 넘버
				
				// // cout << rank << " " << rightProcNo_glo << endl;
				// // cout << rightProcNo_glo << endl;
				// send_rightProcNo_rightId_loc[rightProcNo_glo].push_back(rightProcNo_loc);
				// send_rightProcNo_rightId_loc[rightProcNo_glo].push_back(rightId_loc);
				
				// // int new_my_id_glo = points_id_local2global[ip].at(my_id_loc); // 새로운 포인트 id
			// }
				
		// }
		// vector<vector<int>> recv_rightProcNo_rightId_loc;
		// MASCH_MPI_Alltoallv(send_rightProcNo_rightId_loc, recv_rightProcNo_rightId_loc);
		// vector<vector<int>> send_new_rightId_loc(size);
		// for(int ip=0; ip<size; ++ip){
			// for(int i=0, iter=0; i<recv_rightProcNo_rightId_loc[ip].size()/2; ++i){
				// int rightProcNo_loc = recv_rightProcNo_rightId_loc[ip][iter++];
				// int rightId_loc = recv_rightProcNo_rightId_loc[ip][iter++];
				
				// int new_my_id_glo = points_id_local2global[rightProcNo_loc].at(rightId_loc); // 새로운 포인트 id
				
				// send_new_rightId_loc[ip].push_back(new_my_id_glo);
			// }
		// }
		// vector<vector<int>> recv_new_rightId_loc;
		// MASCH_MPI_Alltoallv(send_new_rightId_loc, recv_new_rightId_loc);
		// // if(rank==0) cout << recv_new_rightId_loc[1].size() << endl;
		// vector<int> ttt_size(size,0);
		// for(int ip=0; ip<size; ++ip){
			// for(int i=0, iter=0, iter2=0; i<recv_globalPoint_toProc_toId[ip].size()/4; ++i){
				// int my_id_loc = recv_globalPoint_toProc_toId[ip][iter++]; // 자신의 local id  
				// int rightProcNo_glo = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 global proc 넘버
				// int rightProcNo_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록 에서 자신과 중복한 포인트의 local proc
				// int rightId_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 자신과 중복한 포인트의 local id 넘버
				
				// // cout << rank << " " << rightProcNo_loc << endl;
				// int new_rightId_glo = recv_new_rightId_loc[rightProcNo_glo].at(ttt_size[rightProcNo_glo]++);
				
				// int new_my_id_glo = points_id_local2global[ip].at(my_id_loc); // 새로운 포인트 id
				
				// auto& connPoints = meshComb.points[new_my_id_glo].connPoints;
				// vector<int> tmp_first;
				// for(auto& [proc, id] : connPoints){
					// tmp_first.push_back(proc);
				// }
				
				// if(find(
				// tmp_first.begin(),tmp_first.end(),rightProcNo_glo)==
				// tmp_first.end()){
					// connPoints.push_back(make_pair(rightProcNo_glo,new_rightId_glo));
				// }
			// }
		// }
		
		
		
		
		
		
	}
	
	
	
	
	
	
	// 원래 메쉬에 넣기
	{
		mesh.cells.clear();
		mesh.cells.shrink_to_fit();
		
		mesh.faces.clear();
		mesh.faces.shrink_to_fit();
		
		mesh.points.clear();
		mesh.points.shrink_to_fit();
		
		mesh.boundaries.clear();
		mesh.boundaries.shrink_to_fit();

		mesh.cells.reserve(meshComb.cells.size());
		mesh.faces.reserve(meshComb.faces.size());
		mesh.points.reserve(meshComb.points.size());
		
		for(auto& point : meshComb.points){
			mesh.addPoint();
			mesh.points.back().x = point.x;
			mesh.points.back().y = point.y;
			mesh.points.back().z = point.z;
			for(auto& [proc, id] : point.connPoints){
				mesh.points.back().connPoints.push_back(make_pair(proc, id));
			}
		}
		
		for(auto& face : meshComb.faces){
			mesh.addFace();
			mesh.faces.back().iL = face.iL;
			mesh.faces.back().iR = face.iR;
			mesh.faces.back().setType(face.getType());
			for(auto& ipoint : face.ipoints){
				mesh.faces.back().ipoints.push_back(ipoint);
			}
		}
		
		for(auto& boundary : meshComb.boundaries){
			mesh.addBoundary();
			mesh.boundaries.back().name = boundary.name;
			mesh.boundaries.back().startFace = boundary.startFace;
			mesh.boundaries.back().nFaces = boundary.nFaces;
			mesh.boundaries.back().rightProcNo = boundary.rightProcNo;
			mesh.boundaries.back().myProcNo = boundary.myProcNo;
			mesh.boundaries.back().setType(boundary.getType());
		}
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.buildCells();
	mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// int tmptmtp=0;
	// for(auto& face : mesh.faces){
		// double tmp0 = pow(mesh.points[face.ipoints[0]].x-mesh.points[face.ipoints[2]].x,2.0);
		// tmp0 += pow(mesh.points[face.ipoints[0]].y-mesh.points[face.ipoints[2]].y,2.0);
		// tmp0 += pow(mesh.points[face.ipoints[0]].z-mesh.points[face.ipoints[2]].z,2.0);
		
		// double tmp1 = pow(mesh.points[face.ipoints[1]].x-mesh.points[face.ipoints[3]].x,2.0);
		// tmp1 += pow(mesh.points[face.ipoints[1]].y-mesh.points[face.ipoints[3]].y,2.0);
		// tmp1 += pow(mesh.points[face.ipoints[1]].z-mesh.points[face.ipoints[3]].z,2.0);
		
		// if(abs(tmp1-tmp0)>1.e-6){
			// cout << face_state[tmptmtp] << endl;
			// cout << "(" << mesh.points[face.ipoints[0]].x << " " <<
			// mesh.points[face.ipoints[0]].y << " " <<
			// mesh.points[face.ipoints[0]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[1]].x << " " <<
			// mesh.points[face.ipoints[1]].y << " " <<
			// mesh.points[face.ipoints[1]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[2]].x << " " <<
			// mesh.points[face.ipoints[2]].y << " " <<
			// mesh.points[face.ipoints[2]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[3]].x << " " <<
			// mesh.points[face.ipoints[3]].y << " " <<
			// mesh.points[face.ipoints[3]].z << ") " << endl;
		// }
		// if(face.ipoints.size()!=4){
			// cout << "EEEE" << endl;
		// }
		
		// ++tmptmtp;
	// }
	
	// for(auto& cell : mesh.cells){
		// vector<int> tmp_ipoints;
	
		// for(auto& iface : cell.ifaces){
			// auto& face = mesh.faces[iface];
			// for(auto& ipoint : face.ipoints){
				// if(find(
				// tmp_ipoints.begin(),tmp_ipoints.end(),ipoint)==
				// tmp_ipoints.end()){
					// tmp_ipoints.push_back( ipoint );
				// }
			// }
		// }
		// if(tmp_ipoints.size()!=8){
			// cout << "NON " << cell.ifaces.size() << endl;
			// cout << "NON " << tmp_ipoints.size() << endl;
			// for(auto& iface : cell.ifaces){
				// auto& face = mesh.faces[iface];
				// cout << face_state[iface] << endl;
				// cout << "(" << mesh.points[face.ipoints[0]].x << " " <<
				// mesh.points[face.ipoints[0]].y << " " <<
				// mesh.points[face.ipoints[0]].z << ") ";
				// cout << "(" << mesh.points[face.ipoints[1]].x << " " <<
				// mesh.points[face.ipoints[1]].y << " " <<
				// mesh.points[face.ipoints[1]].z << ") ";
				// cout << "(" << mesh.points[face.ipoints[2]].x << " " <<
				// mesh.points[face.ipoints[2]].y << " " <<
				// mesh.points[face.ipoints[2]].z << ") ";
				// cout << "(" << mesh.points[face.ipoints[3]].x << " " <<
				// mesh.points[face.ipoints[3]].y << " " <<
				// mesh.points[face.ipoints[3]].z << ") " << endl;
			// }
			
			
			// // for(auto& tmp : tmp_ipoints){
				// // cout << tmp << endl;
			// // }
		// }
		// else{
			// // cout << "NON " << cell.ifaces.size() << endl;
			// // for(auto& iface : cell.ifaces){
				// // cout << face_state[iface] << endl;
			// // }
		// }
		
	// }
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// {
		// MPI_Barrier(MPI_COMM_WORLD);
		// vector<double> send_x0,send_y0,send_z0;
		// vector<double> send_x1,send_y1,send_z1;
		// vector<double> send_x2,send_y2,send_z2;
		// vector<double> send_x3,send_y3,send_z3;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// double avgx = 0.0;
				// double avgy = 0.0;
				// double avgz = 0.0;
				// for(auto& ipoint : face.ipoints){
					// avgx += mesh.points[ipoint].x;
					// avgy += mesh.points[ipoint].y;
					// avgz += mesh.points[ipoint].z;
				// }
				// avgx /= (double)face.ipoints.size();
				// avgy /= (double)face.ipoints.size();
				// avgz /= (double)face.ipoints.size();
				// send_x0.push_back(avgx);
				// send_y0.push_back(avgy);
				// send_z0.push_back(avgz);
				// // send_x0.push_back(mesh.points[face.ipoints[0]].x);
				// // send_y0.push_back(mesh.points[face.ipoints[0]].y);
				// // send_z0.push_back(mesh.points[face.ipoints[0]].z);
				// // send_x1.push_back(mesh.points[face.ipoints[1]].x);
				// // send_y1.push_back(mesh.points[face.ipoints[1]].y);
				// // send_z1.push_back(mesh.points[face.ipoints[1]].z);
				// // send_x2.push_back(mesh.points[face.ipoints[2]].x);
				// // send_y2.push_back(mesh.points[face.ipoints[2]].y);
				// // send_z2.push_back(mesh.points[face.ipoints[2]].z);
				// // send_x3.push_back(mesh.points[face.ipoints[3]].x);
				// // send_y3.push_back(mesh.points[face.ipoints[3]].y);
				// // send_z3.push_back(mesh.points[face.ipoints[3]].z);
			// }
		// }
		// vector<double> recv_x0,recv_y0,recv_z0;
		// // vector<double> recv_x1,recv_y1,recv_z1;
		// // vector<double> recv_x2,recv_y2,recv_z2;
		// // vector<double> recv_x3,recv_y3,recv_z3;
		// recv_x0.resize(send_x0.size());
		// recv_y0.resize(send_y0.size());
		// recv_z0.resize(send_z0.size());
		// // recv_x1.resize(send_x1.size());
		// // recv_y1.resize(send_y1.size());
		// // recv_z1.resize(send_z1.size());
		// // recv_x2.resize(send_x2.size());
		// // recv_y2.resize(send_y2.size());
		// // recv_z2.resize(send_z2.size());
		// // recv_x3.resize(send_x3.size());
		// // recv_y3.resize(send_y3.size());
		// // recv_z3.resize(send_z3.size());
		// MPI_Alltoallv( send_x0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_x0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_y0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_y0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		// MPI_Alltoallv( send_z0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // recv_z0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
					   
		// // MPI_Alltoallv( send_x1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_x1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_y1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_y1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_z1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_z1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
					   
		// // MPI_Alltoallv( send_x2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_x2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_y2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_y2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_z2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_z2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
					   
		// // MPI_Alltoallv( send_x3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_x3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_y3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_y3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
		// // MPI_Alltoallv( send_z3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // recv_z3.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   // // MPI_COMM_WORLD);
					   

		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// if(rank==1){
					// cout << face_state[i] << endl;
					
					// double avgx = 0.0;
					// double avgy = 0.0;
					// double avgz = 0.0;
					// for(auto& ipoint : face.ipoints){
						// avgx += mesh.points[ipoint].x;
						// avgy += mesh.points[ipoint].y;
						// avgz += mesh.points[ipoint].z;
					// }
					// avgx /= (double)face.ipoints.size();
					// avgy /= (double)face.ipoints.size();
					// avgz /= (double)face.ipoints.size();
					// cout << recv_x0[ip] << " " << recv_y0[ip] << " " << recv_z0[ip] << endl;
					// cout << avgx << " " << avgy << " " << avgz << endl;
					
					
					
					// // cout << send_x0[ip] << " " << send_y0[ip] << " " << send_z0[ip] << endl;
					// // cout << recv_x0[ip] << " " << recv_y0[ip] << " " << recv_z0[ip] << endl;
					// // cout << send_x1[ip] << " " << send_y1[ip] << " " << send_z1[ip] << endl;
					// // cout << recv_x1[ip] << " " << recv_y1[ip] << " " << recv_z1[ip] << endl;
					// // cout << send_x2[ip] << " " << send_y2[ip] << " " << send_z2[ip] << endl;
					// // cout << recv_x2[ip] << " " << recv_y2[ip] << " " << recv_z2[ip] << endl;
					// // cout << send_x3[ip] << " " << send_y3[ip] << " " << send_z3[ip] << endl;
					// // cout << recv_x3[ip] << " " << recv_y3[ip] << " " << recv_z3[ip] << endl;
				// }
				
				// ++ip;
			// }
		// }
		
	// }
	
	
	
	
	
	
	// for(auto& cell : meshComb.cells){
		// if(cell.ifaces.size()!=6){
			// cout << "NO " << cell.ifaces.size() << endl;
			// for(auto& iface : cell.ifaces){
				// cout << face_state[iface] << endl;
			// }
		// }
		
	// }

	
	// int tmptmtp=0;
	// for(auto& face : mesh.faces){
		// double tmp0 = pow(mesh.points[face.ipoints[0]].x-mesh.points[face.ipoints[2]].x,2.0);
		// tmp0 += pow(mesh.points[face.ipoints[0]].y-mesh.points[face.ipoints[2]].y,2.0);
		// tmp0 += pow(mesh.points[face.ipoints[0]].z-mesh.points[face.ipoints[2]].z,2.0);
		
		// double tmp1 = pow(mesh.points[face.ipoints[1]].x-mesh.points[face.ipoints[3]].x,2.0);
		// tmp1 += pow(mesh.points[face.ipoints[1]].y-mesh.points[face.ipoints[3]].y,2.0);
		// tmp1 += pow(mesh.points[face.ipoints[1]].z-mesh.points[face.ipoints[3]].z,2.0);
		
		// if(abs(tmp1-tmp0)>0.1){
			// // cout << face_state[tmptmtp] << endl;
			// cout << "(" << mesh.points[face.ipoints[0]].x << " " <<
			// mesh.points[face.ipoints[0]].y << " " <<
			// mesh.points[face.ipoints[0]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[1]].x << " " <<
			// mesh.points[face.ipoints[1]].y << " " <<
			// mesh.points[face.ipoints[1]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[2]].x << " " <<
			// mesh.points[face.ipoints[2]].y << " " <<
			// mesh.points[face.ipoints[2]].z << ") ";
			// cout << "(" << mesh.points[face.ipoints[3]].x << " " <<
			// mesh.points[face.ipoints[3]].y << " " <<
			// mesh.points[face.ipoints[3]].z << ") " << endl;
		// }
		// ++tmptmtp;
	// }
	
	
	// for(auto& point : mesh.points){
		// int iii=0;
		// for(auto& point2 : mesh.points){
			// double resi = 0.0;
			// resi += abs(point.x-point2.x);
			// resi += abs(point.y-point2.y);
			// resi += abs(point.z-point2.z);
			// if(resi<1.e-8) ++iii;
		// }
		
		// if(iii!=1) cout << "ERRORRRRR" << endl;
		
	// }
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	// // 디버그 processor face
	// {
		// vector<double> send_value;
		// vector<double> recv_value;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// send_value.push_back(mesh.points[face.ipoints[0]].x);
			// }
		// }
	// }
	
	
	
	
	
	
	
	
	
	mesh.informations();
	
	// // if(rank==0)
	// {
		
		// cout.precision(20);
		// // for(int ip=0; ip<size; ++ip)
		// {

			// // SEMO_Utility_Math math;
			// // SEMO_Mesh_Geometric geometric;
			// // geometric.init(newMesh[ip]);
			
			// MASCH_Mesh_Save save;
			// save.vtu("./grid/test/", rank, meshComb);
		// }

	
	// }
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	
}




void MASCH_MPI_Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	recv_value.clear();
	displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<double> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	recv_value.resize(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, 
				   recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	displs.resize(size+1,0);
	displs = recv_displs;
	


}



void MASCH_MPI_Alltoallv(vector<vector<int>>& inp_send_value, vector<int>& recv_value, vector<int>& displs){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<int> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	recv_value.resize(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_INT, 
				   recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
	displs = recv_displs;
	


}




void MASCH_MPI_Alltoallv(vector<vector<double>>& inp_send_value, vector<vector<double>>& recv_value){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<double> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	vector<double> tmp_recv_value(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, 
				   tmp_recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	recv_value.clear();
	recv_value.resize(size);
	for(int ip=0; ip<size; ++ip){
		int str = recv_displs[ip];
		int end = recv_displs[ip+1];
		for(int i=str; i<end; ++i){
			recv_value[ip].push_back(tmp_recv_value[i]);
		}
	}


}



void MASCH_MPI_Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<int> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	vector<int> tmp_recv_value(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_INT, 
				   tmp_recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
	recv_value.clear();
	recv_value.resize(size);
	for(int ip=0; ip<size; ++ip){
		int str = recv_displs[ip];
		int end = recv_displs[ip+1];
		for(int i=str; i<end; ++i){
			recv_value[ip].push_back(tmp_recv_value[i]);
		}
	}


}




void MASCH_Gatherv(vector<int>& my_value, vector<int>& value, vector<int>& displs){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	int my_value_size = my_value.size();
	if(rank==0){
		vector<int> counts(size,0);
        MPI_Gather(&my_value_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		displs.clear(); displs.resize(size+1,0);
		for(int ip=0; ip<size; ++ip) displs[ip+1] = displs[ip]+counts[ip];
            
        value.resize(displs[size]);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_INT, value.data(), counts.data(), displs.data(), 
		MPI_INT, 0, MPI_COMM_WORLD);
	}
	else{
        MPI_Gather(&my_value_size, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_INT, NULL, NULL, NULL, 
		MPI_INT, 0, MPI_COMM_WORLD);
	}
	
}



void MASCH_Gatherv(vector<double>& my_value, vector<double>& value, vector<int>& displs){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	int my_value_size = my_value.size();
	if(rank==0){
		vector<int> counts(size,0);
        MPI_Gather(&my_value_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		displs.clear(); displs.resize(size+1,0);
		for(int ip=0; ip<size; ++ip) displs[ip+1] = displs[ip]+counts[ip];
            
        value.resize(displs[size]);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_DOUBLE, value.data(), counts.data(), displs.data(), 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else{
        MPI_Gather(&my_value_size, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_DOUBLE, NULL, NULL, NULL, 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	
}