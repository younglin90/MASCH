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
	
	
	
	combineMesh(idBlockCell, mesh);
	
	 

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

	// cout.precision(20);
	// for(int ip=0; ip<nBlocks; ++ip){

		// // SEMO_Utility_Math math;
		// // SEMO_Mesh_Geometric geometric;
		// // geometric.init(newMesh[ip]);
		
		// MASCH_Mesh_Save save;
		// save.vtu("./grid/0/", ip, newMesh[ip]);
	// }

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




void combineMesh(vector<int>& idBlockCell, MASCH_Mesh &mesh){
	
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	vector<MASCH_Mesh> meshNew;
	MASCH_Mesh meshComb;
	
	
	
	
	// 격자들 각 processor 로 옮기기
	// 포인트 xyz
	vector<double> recv_point_xyz;
	vector<int> recv_point_conn;
	vector<int> displs_point_xyz, displs_point_conn;
	vector<vector<pair<int,int>>> localPointsId(mesh.points.size());
	{
		// 옮길 포인트 정하기
		vector<vector<int>> moveProcNoPoints(mesh.points.size(),vector<int>());
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			int moveProcNo = idBlockCell[i];
			if(rank!=moveProcNo){
				for(auto& ipoint : cell.ipoints){
					if(find(
					moveProcNoPoints[ipoint].begin(),moveProcNoPoints[ipoint].end(),moveProcNo)==
					moveProcNoPoints[ipoint].end()){
						moveProcNoPoints[ipoint].push_back( moveProcNo );
					}
					
				}
			}
		}
		// 포인트 xyz, conn포인트 정보
		vector<vector<double>> point_xyz(size,vector<double>());
		vector<vector<int>> point_local_id_conn(size,vector<int>());
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& iproc : moveProcNoPoints[i]){
				if(rank==iproc) continue;
				
				point_xyz[iproc].push_back(point.x);
				point_xyz[iproc].push_back(point.y);
				point_xyz[iproc].push_back(point.z);
				
				// point_conn[iproc].push_back(point.connPoints.size());
				for(auto& [rightProcNo, rightId] : point.connPoints){
					// point_conn[iproc].push_back(rightProcNo);
					// point_conn[iproc].push_back(rightId);
					
					point_local_id_conn[rightProcNo].push_back(rightId);
				}
				
			}
		}
		// point_local_id_conn 넘기기
		vector<int> recv_point_local_id_conn, displs_point_local_id_conn;
		MASCH_MPI_Alltoallv(point_local_id_conn, recv_point_local_id_conn, displs_point_local_id_conn);
		
		// 로컬 포인트 번호 생성
		vector<int> proc_npoints(size,0);
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& iproc : moveProcNoPoints[i]){
				localPointsId[i].push_back(make_pair(iproc,proc_npoints[iproc]++));
				
			}
		}
		
		// conn id 로컬 변수로 저장
		vector<vector<int>> point_local_id(size,vector<int>());
		{
			for(int ip=0; ip<size; ++ip){
				int str = displs_point_local_id_conn[ip];
				int end = displs_point_local_id_conn[ip+1];
				for(int i=str; i<end; ++i){
					int ipoint = recv_point_local_id_conn[i];
					int local_ipoint = -1;
					for(auto& [rightProcNo, localId] : localPointsId[ipoint]){
						if(rightProcNo==ip){
							local_ipoint = localId;
							break;
						}
					}
					point_local_id[ip].push_back(local_ipoint);
				}
			}
		}
		// conn id 넘기기
		vector<int> recv_point_local_id, displs_point_local_id;
		MASCH_MPI_Alltoallv(point_local_id, recv_point_local_id, displs_point_local_id);
		
		// conn포인트 정보 로컬 변수로 저장
		vector<vector<int>> point_conn(size,vector<int>());
		vector<int> tmp_nnn(size,0);
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& iproc : moveProcNoPoints[i]){
				if(rank==iproc) continue;
				
				point_conn[iproc].push_back(point.connPoints.size());
				for(auto& [rightProcNo, rightId] : point.connPoints){
					point_conn[iproc].push_back(rightProcNo);
					// point_conn[iproc].push_back(rightId);
					int str = displs_point_local_id[rightProcNo];
					int tmp_id = recv_point_local_id[str+tmp_nnn[rightProcNo]];
					++tmp_nnn[rightProcNo];
					point_conn[iproc].push_back(tmp_id);
				}
				
			}
		}
		
		
		
		// 넘기기
		MASCH_MPI_Alltoallv(point_xyz, recv_point_xyz, displs_point_xyz);
		MASCH_MPI_Alltoallv(point_conn, recv_point_conn, displs_point_conn);
		
	}
	
	
	// 페이스의 포인트들
	vector<int> recv_face_ipoints;
	vector<int> displs_face_ipoints;
	{
		// 옮길 페이스 정하기
		vector<vector<int>> moveProcNoFaces(mesh.faces.size(),vector<int>());
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			int moveProcNo = idBlockCell[i];
			if(rank!=moveProcNo){
				for(auto& ifaces : cell.ifaces){
					if(find(
					moveProcNoFaces[ifaces].begin(),moveProcNoFaces[ifaces].end(),moveProcNo)==
					moveProcNoFaces[ifaces].end()){
						moveProcNoFaces[ifaces].push_back( moveProcNo );
					}
					
				}
			}
		}
		
		// 페이스의 포인트들 정보
		vector<vector<int>> face_ipoints(size,vector<int>());
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			for(auto& iproc : moveProcNoFaces[i]){
				if(rank==iproc) continue;
				face_ipoints[iproc].push_back(face.ipoints.size());
				for(auto& ipoint : face.ipoints){
					int local_ipoint = -1;
					for(auto& [rightProcNo, localId] : localPointsId[ipoint]){
						if(rightProcNo==iproc){
							local_ipoint = localId;
							break;
						}
					}
					face_ipoints[iproc].push_back(local_ipoint);
				}
			}
		}
		
		
		// 넘기기
		MASCH_MPI_Alltoallv(face_ipoints, recv_face_ipoints, displs_face_ipoints);
		
		
	}
	
	
	// 페이스의 L, R
	vector<int> recv_face_INTERNAL_iL, recv_face_INTERNAL_iR;
	vector<int> displs_face_INTERNAL_iL, displs_face_INTERNAL_iR;
	vector<int> recv_face_PROCESSOR_iL, recv_face_PROCESSOR_before_iproc;
	vector<int> recv_face_PROCESSOR_after_iproc, recv_face_BOUNDARY_iL;
	vector<int> recv_face_BOUNDARY_typeid, displs_face_PROCESSOR_before_iproc;
	vector<int> displs_face_PROCESSOR_after_iproc, displs_face_BOUNDARY_iL;
	vector<int> displs_face_BOUNDARY_typeid, displs_face_PROCESSOR_iL;
	{
		// 옮길 페이스 정하기
		vector<vector<int>> moveProcNoFaces(mesh.faces.size(),vector<int>());
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			int moveProcNo = idBlockCell[i];
			if(rank!=moveProcNo){
				for(auto& ifaces : cell.ifaces){
					if(find(
					moveProcNoFaces[ifaces].begin(),moveProcNoFaces[ifaces].end(),moveProcNo)==
					moveProcNoFaces[ifaces].end()){
						moveProcNoFaces[ifaces].push_back( moveProcNo );
					}
					
				}
			}
		}
		// processor 옆에 셀 프로세서 (현재, 미래)
		vector<int> recv_idBlockCell;
		{
			vector<int> send_idBlockCell;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				if(face.getType()==MASCH_Face_Types::PROCESSOR){
					send_idBlockCell.push_back(idBlockCell[face.iL]);
				}
			}
			recv_idBlockCell.resize(send_idBlockCell.size());
			MPI_Alltoallv( send_idBlockCell.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
						   recv_idBlockCell.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
						   MPI_COMM_WORLD);
		}
		vector<int> faceL_After_IProc(mesh.faces.size(),-1);
		vector<int> faceR_After_IProc(mesh.faces.size(),-1);
		int iproc_num = 0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			faceL_After_IProc[i] = idBlockCell[face.iL];
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				faceR_After_IProc[i] = idBlockCell[face.iR];
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				faceR_After_IProc[i] = recv_idBlockCell[iproc_num++];
			}
		}
		vector<int> faceR_Before_IProc(mesh.faces.size(),-1);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					faceR_Before_IProc[i] = boundary.rightProcNo;
				}
			}
		}
		// 바운더리 정보
		vector<int> faceR_Before_typeid(mesh.faces.size(),-1);
		int id =0;
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					faceR_Before_typeid[i] = id;
				}
			}
			++id;
		}
		// 로컬 셀 번호 생성
		vector<vector<pair<int,int>>> localCellsId(mesh.cells.size());
		vector<int> proc_ncells(size,0);
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			int iproc = idBlockCell[i];
			if(rank!=iproc){
				localCellsId[i].push_back(make_pair(iproc,proc_ncells[iproc]++));
			}
		}
		// 내부 페이스 부터,...
		vector<vector<int>> face_INTERNAL_iL(size,vector<int>());
		vector<vector<int>> face_INTERNAL_iR(size,vector<int>());
		vector<vector<int>> face_PROCESSOR_iL(size,vector<int>());
		vector<vector<int>> face_PROCESSOR_before_iproc(size,vector<int>());
		vector<vector<int>> face_PROCESSOR_after_iproc(size,vector<int>());
		vector<vector<int>> face_BOUNDARY_iL(size,vector<int>());
		vector<vector<int>> face_BOUNDARY_typeid(size,vector<int>());
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				if(idBlockCell[face.iL] == idBlockCell[face.iR]){
					for(auto& iproc : moveProcNoFaces[i]){
						int local_ipoint = -1;
						for(auto& [rightProcNo, localId] : localCellsId[face.iL]){
							if(rightProcNo==iproc){
								local_ipoint = localId;
								break;
							}
						}
						face_INTERNAL_iL[iproc].push_back(local_ipoint);
						
						local_ipoint = -1;
						for(auto& [rightProcNo, localId] : localCellsId[face.iR]){
							if(rightProcNo==iproc){
								local_ipoint = localId;
								break;
							}
						}
						face_INTERNAL_iR[iproc].push_back(local_ipoint);
					}
				}
				else{
					for(auto& iproc : moveProcNoFaces[i]){
						int local_ipoint = -1;
						for(auto& [rightProcNo, localId] : localCellsId[face.iL]){
							if(rightProcNo==iproc){
								local_ipoint = localId;
								break;
							}
						}
						face_PROCESSOR_iL[iproc].push_back(local_ipoint);
						face_PROCESSOR_before_iproc[iproc].push_back(rank);
						if(idBlockCell[face.iL] == iproc){
							face_PROCESSOR_after_iproc[iproc].push_back(faceR_After_IProc[i]);
						}
						else{
							face_PROCESSOR_after_iproc[iproc].push_back(faceL_After_IProc[i]);
						}
					}
				}
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				for(auto& iproc : moveProcNoFaces[i]){
					int local_ipoint = -1;
					for(auto& [rightProcNo, localId] : localCellsId[face.iL]){
						if(rightProcNo==iproc){
							local_ipoint = localId;
							break;
						}
					}
					face_PROCESSOR_iL[iproc].push_back(local_ipoint);
					face_PROCESSOR_before_iproc[iproc].push_back(faceR_Before_IProc[i]);
					face_PROCESSOR_after_iproc[iproc].push_back(faceR_After_IProc[i]);
				}
			}
			else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				for(auto& iproc : moveProcNoFaces[i]){
					int local_ipoint = -1;
					for(auto& [rightProcNo, localId] : localCellsId[face.iL]){
						if(rightProcNo==iproc){
							local_ipoint = localId;
							break;
						}
					}
					face_BOUNDARY_iL[iproc].push_back(local_ipoint);
					face_BOUNDARY_typeid[iproc].push_back(faceR_Before_typeid[i]);
				}
			}
		}
		
		
		// 넘기기
		MASCH_MPI_Alltoallv(face_INTERNAL_iL, recv_face_INTERNAL_iL, displs_face_INTERNAL_iL);
		MASCH_MPI_Alltoallv(face_INTERNAL_iR, recv_face_INTERNAL_iR, displs_face_INTERNAL_iR);
		MASCH_MPI_Alltoallv(face_PROCESSOR_iL, recv_face_PROCESSOR_iL, displs_face_PROCESSOR_iL);
		MASCH_MPI_Alltoallv(face_PROCESSOR_before_iproc, recv_face_PROCESSOR_before_iproc, displs_face_PROCESSOR_before_iproc);
		MASCH_MPI_Alltoallv(face_PROCESSOR_after_iproc, recv_face_PROCESSOR_after_iproc, displs_face_PROCESSOR_after_iproc);
		MASCH_MPI_Alltoallv(face_BOUNDARY_iL, recv_face_BOUNDARY_iL, displs_face_BOUNDARY_iL);
		MASCH_MPI_Alltoallv(face_BOUNDARY_typeid, recv_face_BOUNDARY_typeid, displs_face_BOUNDARY_typeid);
		
	}
	
	
	
	// 옮겨진 격자 생성
	meshNew.resize(size);
	meshNew[rank] = mesh;
	
	{
		for(int ip=0; ip<size; ++ip){
			if(ip==rank) continue;
			
			// 포인트 제작
			{
				int str = displs_point_xyz[ip];
				int end = displs_point_xyz[ip+1];
				for(int i=str; i<end; ++i){
					meshNew[ip].addPoint();
					meshNew[ip].points.back().x = recv_point_xyz[i++];
					meshNew[ip].points.back().y = recv_point_xyz[i++];
					meshNew[ip].points.back().z = recv_point_xyz[i];
					
				}
				
				str = displs_point_conn[ip];
				end = displs_point_conn[ip+1];
				int iter = 0;
				for(int i=str; i<end; ++i){
					int tmp_size = recv_point_conn[i];
					for(int j=0; j<tmp_size; ++j){
						meshNew[ip].points[iter].connPoints.push_back(
						make_pair(recv_point_conn[i+1],recv_point_conn[i+2]));
						++i;
						++i;
					}
					++iter;
				}
			}
			
			// 페이스 제작
			int startBCFace = 0;
			{
				int n=0;
				meshNew[ip].faces.clear();
				
				int str = displs_face_INTERNAL_iL[ip];
				int end = displs_face_INTERNAL_iL[ip+1];
				startBCFace = end-str;
				for(int i=str; i<end; ++i){
					meshNew[ip].addFace();
					meshNew[ip].faces.back().iL = recv_face_INTERNAL_iL[i];
					meshNew[ip].faces.back().iR = recv_face_INTERNAL_iR[i];
					meshNew[ip].faces.back().setType(MASCH_Face_Types::INTERNAL);
				}
				str = displs_face_BOUNDARY_iL[ip];
				end = displs_face_BOUNDARY_iL[ip+1];
				for(int i=str; i<end; ++i){
					meshNew[ip].addFace();
					meshNew[ip].faces.back().iL = recv_face_BOUNDARY_iL[i];
					meshNew[ip].faces.back().setType(MASCH_Face_Types::BOUNDARY);
				}
				str = displs_face_PROCESSOR_iL[ip];
				end = displs_face_PROCESSOR_iL[ip+1];
				for(int i=str; i<end; ++i){
					meshNew[ip].addFace();
					meshNew[ip].faces.back().iL = recv_face_PROCESSOR_iL[i];
					meshNew[ip].faces.back().setType(MASCH_Face_Types::PROCESSOR);
				}
			}
			
			// 셀 제작
			{
				int ncells=-1;
				for(auto& face : meshNew[ip].faces){
					ncells = max(ncells,face.iL);
				}
				for(int i=0; i<ncells+1; ++i){
					meshNew[ip].addCell();
				}
			}
			
			// 페이스 포인트 컨넥트
			{
				int str = displs_face_ipoints[ip];
				int end = displs_face_ipoints[ip+1];
				int iter=0;
				for(int i=str; i<end; ++i){
					int tmp_size = recv_face_ipoints[i];
					for(int j=0; j<tmp_size; ++j){
						meshNew[ip].faces[iter].ipoints.push_back(recv_face_ipoints[++i]);
					}
					++iter;
				}
				
				
			}
			
			// 바운더리
			{
				int nnnn = startBCFace;
				
				int str = displs_face_BOUNDARY_typeid[ip];
				int end = displs_face_BOUNDARY_typeid[ip+1];
				vector<int> boundary_nfaces(mesh.boundaries.size(),0);
				for(int i=str; i<end; ++i){
					int type = recv_face_BOUNDARY_typeid[i];
					++boundary_nfaces[type];
				}
				
				for(int ibcs=0; ibcs<mesh.boundaries.size(); ++ibcs){
					if(mesh.boundaries[ibcs].getType() != MASCH_Face_Types::BOUNDARY) continue;
					meshNew[ip].addBoundary();
					meshNew[ip].boundaries.back().name = mesh.boundaries[ibcs].name;
					
					if( boundary_nfaces[ibcs] > 0) {
						meshNew[ip].boundaries.back().nFaces = boundary_nfaces[ibcs];
						meshNew[ip].boundaries.back().startFace = nnnn;
						meshNew[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
						
						nnnn += boundary_nfaces[ibcs];
					}
					else{
						meshNew[ip].boundaries.back().nFaces = 0;
						meshNew[ip].boundaries.back().startFace = 0;
						meshNew[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
					}
				} 
				
				
				
				str = displs_face_PROCESSOR_after_iproc[ip];
				end = displs_face_PROCESSOR_after_iproc[ip+1];
				vector<int> proc_nfaces(size,0);
				for(int i=str; i<end; ++i){
					int rightProcNo = recv_face_PROCESSOR_after_iproc[i];
					++proc_nfaces[rightProcNo];
				}
				
				for(int j=0; j<size; ++j){
					meshNew[ip].addBoundary();
					string bcnames = "procBoundary" + to_string(rank) + "to" + to_string(j);
					
					meshNew[ip].boundaries.back().name = bcnames;
					meshNew[ip].boundaries.back().nFaces = proc_nfaces[j];
					meshNew[ip].boundaries.back().startFace = nnnn;
					meshNew[ip].boundaries.back().myProcNo = rank;
					meshNew[ip].boundaries.back().rightProcNo = j;
					meshNew[ip].boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
					
					nnnn += proc_nfaces[j];
					
				}
			}	
			
			
			{
				int maxBCNum = meshNew[ip].boundaries.size()-1;
				// if(rank==0) cout << meshNew[ip].boundaries[maxBCNum].startFace << endl;
				meshNew[ip].boundaries[maxBCNum].startFace = meshNew[ip].faces.size()-meshNew[ip].boundaries[maxBCNum].nFaces;
				for(int i=maxBCNum-1; i>=0; --i){
					meshNew[ip].boundaries[i].startFace = meshNew[ip].boundaries[i+1].startFace-meshNew[ip].boundaries[i].nFaces;
				}
				// if(rank==0) cout << meshNew[ip].boundaries[maxBCNum].startFace << endl;
			}
			
			meshNew[ip].check();
			meshNew[ip].setFaceTypes();
			// newMesh[ip].buildCells();
			meshNew[ip].connectFacetoPointsCells();
			meshNew[ip].connectCelltoFaces();
			meshNew[ip].connectCelltoPoints();
			meshNew[ip].setCountsProcFaces();
			meshNew[ip].setDisplsProcFaces();
		
		// cout << rank << endl;
			meshNew[ip].informations();
			
	
		}
		
		
		
		
	}
	
	
	
	
	
	
	// 현재 격자에서 떨어져나간 격자 제거
	{
		
		
		
		
		
	}
	
	
	
	
	
	
	
	
	
	// 옮겨진 격자들 합치기
	{
		// 삭제될 포인트들, 새로운 로컬 아이디
		vector<vector<bool>> delete_points(size,vector<bool>());
		vector<vector<int>> new_points_local_id(size,vector<int>());
		for(int ip=0, iter=0; ip<size; ++ip){
			for(int i=0; i<meshNew[ip].points.size(); ++i){
				delete_points[ip].push_back(true);
				if(delete_points[ip].back()==true){
					new_points_local_id[ip].push_back();
				}
			}
		}
		
		// 새로운 포인트들 아이디
		vector<int> new_points_id;
		vector<int> new_points_id_displs(size+1,0);
		for(int ip=0, iter=0; ip<size; ++ip){
			for(int i=0, iter2=0; i<meshNew[ip].points.size(); ++i){
				if(delete_points[ip][i]==false){
					new_points_id.push_back(iter++);
				}
				else{
					int new_local_id = new_points_local_id[ip][iter2++];
					new_points_id.push_back(new_local_id);
				}
			}
			new_points_id_displs[ip+1] = static_cast<int>(new_points_id.size());
		}
		
		
		
		
		
		// 중복 포인트 확인
		vector<vector<int>> overlapPoints(size,vector<int>());
		for(int ip=0; ip<size; ++ip){
			if(ip==rank) continue;
			for(auto& point : meshNew[ip].points){
				for(auto& [rightProcNo, rightId] : point.connPoints){
					if(rightProcNo<ip){
						overlapPoints[rightProcNo].push_back(rightId);
					}
				}
			}
		}
		
		
		for(int ip=0; ip<size; ++ip){
			if(ip==rank) continue;
			for(auto& point : meshNew[ip].points){
				for(auto& ipoint : overlapPoints[ip]){
					meshNew[ip].points[ipoint].
				}
			}
		}
		
		
		// 중복 포인트 제거
		// for(int ip=0; ip<size; ++ip){
			// if(ip==rank) continue;
			
			
			
			
			// meshNew[ip]
		// }
		
		
		// 포인트 합치기
		// for(int ip=0; ip<size; ++ip){
			// if(ip==rank) continue;
			
			
			
			
			// meshNew[ip]
		// }
		
		
		
		// 페이스 합치기
		// 페이스 포인트 합치기
		
		
		
		
		
		// 셀 합치기
		
		
		
	}
	
	
	
	
	
	
	
	
	
	if(rank==0){
		
		cout.precision(20);
		for(int ip=0; ip<size; ++ip){

			// SEMO_Utility_Math math;
			// SEMO_Mesh_Geometric geometric;
			// geometric.init(newMesh[ip]);
			
			MASCH_Mesh_Save save;
			save.vtu("./grid/test/", ip, meshNew[ip]);
		}

	
	}
	
	
	
	// // cout << mesh.startCellGlobal << endl;
	// // cout << mesh.ncellsTotal << endl;
	// // cout << mesh.startProcCellGlobal.size() << endl;
	
	// // int my_nCells = mesh.cells.size();
	// // int total_nCells;
    // // MPI_Reduce(&my_nCells, &total_nCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	// // int my_nFaces = mesh.faces.size();
	// // int total_nFaces;
    // // MPI_Reduce(&my_nFaces, &total_nFaces, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	// // int my_nPoints = mesh.points.size();
	// // int total_nPoints;
    // // MPI_Reduce(&my_nPoints, &total_nPoints, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	
	
	// // vector<int> my_iLs;
	// // for(auto& face : mesh.faces) {
		// // my_iLs.push_back(face.iL);
	// // }
	
	
	
	// // vector<int> my_iRs;
	// // for(auto& face : mesh.faces) my_iRs.push_back(face.iR);
	
	
	// // vector<int> my_face_ipoints;
	// // for(auto& face : mesh.faces) {
		// // for(auto& point : face.ipoints){
			// // my_face_ipoints.push_back(point);
		// // }
	// // }
	
	
	
	// // vector<double> my_point_x;
	// // for(auto& point : mesh.points) my_point_x.push_back(point.x);
	// // vector<double> my_point_y;
	// // for(auto& point : mesh.points) my_point_y.push_back(point.y);
	// // vector<double> my_point_z;
	// // for(auto& point : mesh.points) my_point_z.push_back(point.z);
		
		
		
	// vector<int> my_nCells;
	// my_nCells.push_back(static_cast<int>(mesh.cells.size()));
	// vector<int> nCells;
	// vector<int> nCells_disls;
	// MASCH_Gatherv(my_nCells, nCells, nCells_disls);
	
	// vector<int> my_nFaces;
	// my_nFaces.push_back(static_cast<int>(mesh.faces.size()));
	// vector<int> nFaces;
	// vector<int> nFaces_disls;
	// MASCH_Gatherv(my_nFaces, nFaces, nFaces_disls);
	
	// vector<int> my_nPoints;
	// my_nPoints.push_back(static_cast<int>(mesh.points.size()));
	// vector<int> nPoints;
	// vector<int> nPoints_disls;
	// MASCH_Gatherv(my_nPoints, nPoints, nPoints_disls);
	
	
	
	// vector<double> my_NodeCoordinates;
	// for(auto& point : mesh.points){
		// my_NodeCoordinates.push_back(point.x);
		// my_NodeCoordinates.push_back(point.y);
		// my_NodeCoordinates.push_back(point.z);
	// }
	// vector<double> NodeCoordinates;
	// vector<int> NodeCoordinates_disls;
	// MASCH_Gatherv(my_NodeCoordinates, NodeCoordinates, NodeCoordinates_disls);
	
		

	// vector<int> my_iLs_INTERNAL;
	// for(auto& face : mesh.faces) {
		// if(face.getType()==MASCH_Face_Types::INTERNAL){
			// my_iLs_INTERNAL.push_back(face.iL);
		// }
	// }
	// vector<int> iLs_INTERNAL;
	// vector<int> iLs_INTERNAL_disls;
	// MASCH_Gatherv(my_iLs_INTERNAL, iLs_INTERNAL, iLs_INTERNAL_disls);
	
	
	// vector<int> my_iLs_BOUNDARY;
	// for(auto& face : mesh.faces) {
		// if(face.getType()==MASCH_Face_Types::BOUNDARY){
			// my_iLs_BOUNDARY.push_back(face.iL);
		// }
	// }
	// vector<int> iLs_BOUNDARY;
	// vector<int> iLs_disls_BOUNDARY;
	// MASCH_Gatherv(my_iLs_BOUNDARY, iLs_BOUNDARY, iLs_disls_BOUNDARY);
	
	
	// vector<int> my_iLs_PROCESSOR;
	// for(auto& face : mesh.faces) {
		// if(face.getType()==MASCH_Face_Types::PROCESSOR){
			// my_iLs_PROCESSOR.push_back(face.iL);
		// }
	// }
	// vector<int> iLs_PROCESSOR;
	// vector<int> iLs_disls_PROCESSOR;
	// MASCH_Gatherv(my_iLs_PROCESSOR, iLs_PROCESSOR, iLs_disls_PROCESSOR);
	
	
	
	
	

	// vector<int> my_iRs_INTERNAL;
	// for(auto& face : mesh.faces) {
		// if(face.getType()==MASCH_Face_Types::INTERNAL){
			// my_iRs_INTERNAL.push_back(face.iR);
		// }
	// }
	// vector<int> iRs_INTERNAL;
	// vector<int> iRs_INTERNAL_disls;
	// MASCH_Gatherv(my_iRs_INTERNAL, iRs_INTERNAL, iRs_INTERNAL_disls);
	
	
	
	
	
	// vector<int> my_face_ipoints_INTERNAL;
	// vector<int> my_face_ipoints_BOUNDARY;
	// vector<int> my_face_ipoints_PROCESSOR;
	// for(auto& face : mesh.faces) {
		// if(face.getType()==MASCH_Face_Types::INTERNAL){
			// my_face_ipoints_INTERNAL.push_back(face.ipoints.size());
			// for(auto& point : face.ipoints){
				// my_face_ipoints_INTERNAL.push_back(point);
			// }
		// }
		// else if(face.getType()==MASCH_Face_Types::BOUNDARY){
			// my_face_ipoints_BOUNDARY.push_back(face.ipoints.size());
			// for(auto& point : face.ipoints){
				// my_face_ipoints_BOUNDARY.push_back(point);
			// }
		// }
		// else if(face.getType()==MASCH_Face_Types::PROCESSOR){
			// my_face_ipoints_PROCESSOR.push_back(face.ipoints.size());
			// for(auto& point : face.ipoints){
				// my_face_ipoints_PROCESSOR.push_back(point);
			// }
		// }
	// }
	// vector<int> face_ipoints_INTERNAL;
	// vector<int> face_ipoints_INTERNAL_disls;
	// MASCH_Gatherv(my_face_ipoints_INTERNAL, face_ipoints_INTERNAL, face_ipoints_INTERNAL_disls);
	// vector<int> face_ipoints_BOUNDARY;
	// vector<int> face_ipoints_BOUNDARY_disls;
	// MASCH_Gatherv(my_face_ipoints_BOUNDARY, face_ipoints_BOUNDARY, face_ipoints_BOUNDARY_disls);
	// vector<int> face_ipoints_PROCESSOR;
	// vector<int> face_ipoints_PROCESSOR_disls;
	// MASCH_Gatherv(my_face_ipoints_PROCESSOR, face_ipoints_PROCESSOR, face_ipoints_PROCESSOR_disls);
	
	
	
	
	// int iter = 0;
	// vector<int> my_point_connPoints;
	// for(auto& point : mesh.points) {
		// if(point.connPoints.size()>0){
			// my_point_connPoints.push_back(iter);
			// my_point_connPoints.push_back(point.connPoints.size());
			// for(auto& [neighb, id] : point.connPoints){
				// my_point_connPoints.push_back(neighb);
				// my_point_connPoints.push_back(id);
			// }
		// }
		// ++iter;
	// }
	// vector<int> point_connPoints;
	// vector<int> point_connPoints_disls;
	// MASCH_Gatherv(my_point_connPoints, point_connPoints, point_connPoints_disls);
	
	
	
	
	
	
	// if(rank==0){
		
		
		// vector<bool> boolOverlapPoints(point_connPoints_disls[size],false);
		// vector<int> new_ipoints(point_connPoints_disls[size],0);
		// iter = 0;
		// for(int ip=0; ip<size; ++ip){
			// int str = point_connPoints_disls[ip];
			// int end = point_connPoints_disls[ip+1];
			// for(int i=str; i<end; ++i){
				// int localId = point_connPoints[i++];
				// int connSize = point_connPoints[i++];
				// int rightProcNo = point_connPoints[i++];
				// int rightLocalId = point_connPoints[i];
				
				// if(rightProcNo<ip){
					// boolOverlapPoints[i] = true;
					// new_ipoints[i] = localId;
				// }
				// else{
					// new_ipoints[i] = iter++;
				// }
			// }
		// }
		
		
		
		
		
		
		// // for(int ip=0; ip<nPoints_disls.size()-1; ++ip){
			// // int str = nPoints_disls[ip];
			// // int end = nPoints_disls[ip+1];
			// // for(int i=str; i<end; ++i){
				// // nPoints[i]
			// // }
		// // }


		
		// for(int ip=0; ip<size; ++ip){
			// int str = NodeCoordinates_disls[ip];
			// int end = NodeCoordinates_disls[ip+1];
			// for(int i=str; i<end; ++i){
				// meshComb.addPoint();
				// meshComb.points.back().x = NodeCoordinates[i++];
				// meshComb.points.back().y = NodeCoordinates[i++];
				// meshComb.points.back().z = NodeCoordinates[i];
			// }
		// }
		
		

		
		// // for(int i=0; i<NodeCoordinates.size()/3; ++i){
			// // meshComb.addPoint();
			// // meshComb.points.back().x = NodeCoordinates[i*3+0];
			// // meshComb.points.back().y = NodeCoordinates[i*3+1];
			// // meshComb.points.back().z = NodeCoordinates[i*3+2];
		// // }
		// // NodeCoordinates.clear();
		
		
		// // int n=0;
		// // int ncells=0;
		// // meshComb.faces.clear();
		// // for(auto& i : owner){
			// // meshComb.addFace();
			// // meshComb.faces.back().iL = i;
			// // ncells = max(ncells, meshComb.faces.back().iL);
		// // }
		// // owner.clear();
		
		// // // cout << ncells << endl;
		
		// // meshComb.cells.clear();
		// // for(int i=0; i<ncells+1; ++i){
			// // meshComb.addCell();
		// // }
		
		
		
		// // n=0;
		// // for(auto& i : neighbour){
			// // meshComb.faces[n].iR = i;
			// // meshComb.faces[n].setType(MASCH_Face_Types::INTERNAL);
			// // ++n;
		// // }
		// // neighbour.clear();
		
		
		// // int m=0;
		// // n=0;
		// // for(auto& i : offsets){
			// // for(int j=n; j<i; ++j){
				// // int point = connectivity[j];
				// // // cout << m << endl;
				// // meshComb.cells[m].ipoints.push_back( static_cast<int>(point) );
				// // // if(rank==1) cout << point << endl;
			// // }
			
			// // // if(rank==0 && mesh.cells[m].ipoints.size()!=8) cout << mesh.cells[m].ipoints.size() << endl;
			
			// // n=i;
			// // ++m;
		// // }
		
		
		
		// // n=0;
		// // int nFacesInt=0;
		// // for(auto& face : meshComb.faces){
			// // // if(face.iR != -1){
			// // if(face.getType() == MASCH_Face_Types::INTERNAL){
				// // meshComb.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
				// // meshComb.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
				// // ++nFacesInt;
			// // }
			// // else{
				// // meshComb.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			// // }
			// // ++n;
		// // }
		
		
		// // m=0;
		// // n=0;
		// // for(auto& i : faceoffsets){
			// // // if(faces[n]>5) cout << faces[n] << endl;
			// // int N=0;
			// // int face_size = faces[m+N];
			// // for(int j=0; j<face_size; ++j){
				// // int face = meshComb.cells[n].ifaces[j];
				// // ++N;
				// // int point_size = faces[m+N];
				// // for(int k=0; k<point_size; ++k){
					// // ++N;
					// // int point = faces[m+N];
					// // if(meshComb.faces[ face ].ipoints.size() == point_size) continue;
					// // meshComb.faces[ face ].ipoints.push_back( static_cast<int>(point) );
					// // // if(rank==1) cout << point << endl;
				// // }
			// // }
			// // m=i;
			// // ++n;
		// // }
		// // faces.clear();
		// // faceoffsets.clear();
		
		
		
		// // n=0;
		// // for(auto& startFace : bcStartFace){
			
			// // meshComb.addBoundary();
			
			// // string tmp_bcName = bcName[n];
			// // tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
			// // tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
			// // tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
			// // tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
				
			// // meshComb.boundaries.back().name = tmp_bcName;
			// // meshComb.boundaries.back().startFace = bcStartFace[n];
			// // meshComb.boundaries.back().nFaces = bcNFaces[n];
			// // if(bcNeighbProcNo[n]<0){
				// // meshComb.boundaries.back().rightProcNo = 0;
				// // meshComb.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			// // }
			// // else{
				// // meshComb.boundaries.back().rightProcNo = bcNeighbProcNo[n];
				// // meshComb.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
			// // }
			// // meshComb.boundaries.back().myProcNo = rank;
			// // // if(bcNeighbProcNo[n] < 0){
				// // // mesh.boundaries.back().myProcNo = -1;
			// // // }
			
			// // ++n;
			
		// // }
		
		// // int maxBCNum = meshComb.boundaries.size()-1;
		// // meshComb.boundaries[maxBCNum].startFace = meshComb.faces.size()-meshComb.boundaries[maxBCNum].nFaces;
		// // for(int i=maxBCNum-1; i>=0; --i){
			// // meshComb.boundaries[i].startFace = meshComb.boundaries[i+1].startFace-meshComb.boundaries[i].nFaces;
		// // }
		
		// // bcName.clear();
		// // bcStartFace.clear();
		// // bcNFaces.clear();
		// // bcNeighbProcNo.clear();
		
			
		// // for(int i=0; i<connPoints.size(); ++i){
			// // // cout << mesh.points.size() << " " << connPoints[i] << " " << connPoints[i+1] << " " << connPoints[i+2] << " " << endl;
			// // meshComb.points[connPoints[i]].connPoints.push_back(
			// // make_pair(connPoints[i+1],connPoints[i+2])
			// // );
			// // ++i; ++i;
		// // }
		
		
		
		// // if(rank==0){
			// // cout << "-> completed" << endl;
			// // cout << "└────────────────────────────────────────────────────" << endl;
		// // }
		
		// // // for(auto& i : mesh.cells){
			
			// // // if(rank==0 && i.ipoints.size()!=8) cout << i.ipoints.size() << endl;
			// // // if(rank==0 && i.ifaces.size()!=6) cout << i.ifaces.size() << endl;
			
		// // // }
		
		// // meshComb.check();
		// // meshComb.setFaceTypes();
		// // // mesh.buildCells();
		// // // mesh.connectFacetoPointsCells();
		// // // mesh.connectCelltoFaces();
		// // // mesh.connectCelltoPoints();
		// // meshComb.setCountsProcFaces();
		// // meshComb.setDisplsProcFaces();
		// // meshComb.cellsGlobal();
		
		
		// // meshComb.informations();
		
					
	// }
	
		
	// // if(rank==0){
		// // // cout << total_nCells << endl;
		// // // cout << total_nFaces << endl;
		// // // cout << total_nPoints << endl;
		
		// // vector<int> cell_counts(size,0);
        // // MPI_Gather(&my_nCells, 1, MPI_INT, cell_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// // vector<int> cell_displs(size+1,0);
		// // for(int ip=0; ip<size; ++ip) cell_displs[ip+1] = cell_displs[ip]+cell_counts[ip];
		
		// // vector<int> face_counts(size,0);
        // // MPI_Gather(&my_nFaces, 1, MPI_INT, face_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// // vector<int> face_displs(size+1,0);
		// // for(int ip=0; ip<size; ++ip) face_displs[ip+1] = face_displs[ip]+face_counts[ip];
		
		// // vector<int> point_counts(size,0);
        // // MPI_Gather(&my_nPoints, 1, MPI_INT, point_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// // vector<int> point_displs(size+1,0);
		// // for(int ip=0; ip<size; ++ip) point_displs[ip+1] = point_displs[ip]+point_counts[ip];
		
		
		
		
		
		
        // // vector<int> iLs(total_nFaces,0);
		// // MPI_Gatherv(my_iLs.data(), my_nFaces, MPI_INT, iLs.data(), face_counts.data(), face_displs.data(), 
		// // MPI_INT, 0, MPI_COMM_WORLD);
            
        // // vector<int> iRs(total_nFaces,0);
		// // MPI_Gatherv(my_iRs.data(), my_nFaces, MPI_INT, iRs.data(), face_counts.data(), face_displs.data(), 
		// // MPI_INT, 0, MPI_COMM_WORLD);
            
		
		
		
		
		// // vector<int> face_ipoints_counts(size,0);
		// // int size_my_face_ipoints = my_face_ipoints.size();
        // // MPI_Gather(&size_my_face_ipoints, 1, MPI_INT, face_ipoints_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		// // vector<int> face_ipoints_displs(size+1,0);
		// // for(int ip=0; ip<size; ++ip) face_ipoints_displs[ip+1] = face_ipoints_displs[ip]+face_ipoints_counts[ip];
            
        // // vector<int> face_ipoints(face_ipoints_displs[size]);
		// // MPI_Gatherv(my_face_ipoints.data(), size_my_face_ipoints, MPI_INT, face_ipoints.data(), face_ipoints_counts.data(), face_ipoints_displs.data(), 
		// // MPI_INT, 0, MPI_COMM_WORLD);
            
            
			
			
			
			
        // // vector<double> point_x(total_nPoints);
		// // MPI_Gatherv(my_point_x.data(), my_nPoints, MPI_DOUBLE, point_x.data(), point_counts.data(), point_displs.data(), 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // // vector<double> point_y(total_nPoints);
		// // MPI_Gatherv(my_point_y.data(), my_nPoints, MPI_DOUBLE, point_y.data(), point_counts.data(), point_displs.data(), 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // // vector<double> point_z(total_nPoints);
		// // MPI_Gatherv(my_point_z.data(), my_nPoints, MPI_DOUBLE, point_z.data(), point_counts.data(), point_displs.data(), 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			
			
			
			
        // // // MPI_Gatherv(&gatherValue, 1, MPI_INT, Ls.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		
        // // // vector<int> buffer(mesh.ncellsTotal,0);
		
		
		
		
	// // }
	// // else{
		
        // // MPI_Gather(&my_nCells, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        // // MPI_Gather(&my_nFaces, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        // // MPI_Gather(&my_nPoints, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		// // MPI_Gatherv(my_iLs.data(), my_nFaces, MPI_INT, NULL, NULL, NULL, 
		// // MPI_INT, 0, MPI_COMM_WORLD);
		
		// // MPI_Gatherv(my_iRs.data(), my_nFaces, MPI_INT, NULL, NULL, NULL, 
		// // MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
		
		
		// // int size_my_face_ipoints = my_face_ipoints.size();
        // // MPI_Gather(&size_my_face_ipoints, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		// // MPI_Gatherv(my_face_ipoints.data(), size_my_face_ipoints, MPI_INT, NULL, NULL, NULL, 
		// // MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
		
		
		// // MPI_Gatherv(my_point_x.data(), my_nPoints, MPI_DOUBLE, NULL, NULL, NULL, 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// // MPI_Gatherv(my_point_y.data(), my_nPoints, MPI_DOUBLE, NULL, NULL, NULL, 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// // MPI_Gatherv(my_point_z.data(), my_nPoints, MPI_DOUBLE, NULL, NULL, NULL, 
		// // MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		// // // int gatherValue = mesh.points.size();
        // // // MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
		
		
		
		
		
		
		
	// // }
	
	
	
	
	
	
	
	
	
	
	
	
}




void MASCH_MPI_Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs){

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