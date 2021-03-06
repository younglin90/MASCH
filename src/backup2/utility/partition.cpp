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
#include "../mesh/partition/partition.h" 

// void parMETIS_Graph_Partition(int nSize, vector<int>& procNoCell);
// // void parMETIS_Mesh_Partition(int nSize, vector<int>& procNoCell);
// void partitionFromSerial(int nSize, vector<int>& procNoCell, vector<MASCH_Mesh>& newMesh);


int main(int argc, char* argv[]) {
	

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	map<string,string> mapArgv;
	
	for(int i=1; i<argc; i+=2){
		string first = argv[i];
		string second;
		if(i+1==argc){
			second = "nan";
		}
		else{
			second = argv[i+1];
		}
		mapArgv.insert(make_pair(first, second));
	}
	
	if( 
	mapArgv.find("-help") != mapArgv.end() ||
	mapArgv.find("-h") != mapArgv.end()
	){
		if(rank==0){
			cout << endl;
			cout << "┌─────── Partitioning helper ─────────────────────────────── " << endl;
			cout << "| -n \"int num.\"   : # of partition" << endl;
			cout << "| -s \"real num.\"  : scale of mesh" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	
	int nSize = stoi(mapArgv["-n"]);
	double scaleMesh = stod(mapArgv["-s"]);
	
	
	MASCH_Mesh mesh;
	
	// if(rank==0)
	{
		// cout << "AAAAAAAA" << endl;
		MASCH_Mesh_Load load;
		load.OpenFoam("./grid/", mesh);
	}
	// else{
		
	// }
	

	// cout << mesh.cells.size() << endl;
	vector<int> procNoCell(static_cast<int>(mesh.cells.size()),0);

	// {
		// int nIter = mesh.cells.size()/nSize;
		// int iter=0, inpRank=0;
		// for(auto& item : procNoCell){
			// item = inpRank;
			// ++iter;
			// if(iter % nIter == 0) ++inpRank;
		// }
	// }

	MASCH_Mesh_Partition partition;
	// partition.combine(procNoCell, mesh);
	
	partition.parMETIS_Graph_Partition(nSize, procNoCell, mesh);

	
	vector<MASCH_Mesh> newMesh(nSize,MASCH_Mesh());
	partition.partitionFromSerial(nSize, procNoCell, mesh, newMesh);
	
	// partition.parMETIS_Graph_Partition(nSize, procNoCell, mesh);
	// // parMETIS_Mesh_Partition(nSize, procNoCell);
	
	// partition.combineMesh(procNoCell, mesh)

	// // vector<MASCH_Mesh> newMesh(nSize,MASCH_Mesh());
	// // partitionFromSerial(nSize, procNoCell, newMesh);
	
	// 스케일
	for(auto& item : newMesh) {
		for(auto& point : item.points){
			point.x *= scaleMesh;
			point.y *= scaleMesh;
			point.z *= scaleMesh;
		}
	}

	{
		cout.precision(20);
		for(int ip=0; ip<nSize; ++ip){

			// SEMO_Utility_Math math;
			// SEMO_Mesh_Geometric geometric;
			// geometric.init(newMesh[ip]);
			
			MASCH_Mesh_Save save;
			save.vtu("./grid/0/", ip, newMesh[ip]);
		}
	}
	// if(rank==0)
		
	
	// {
		// cout.precision(20);
		// // for(int ip=0; ip<size; ++ip)
		// {
			// MASCH_Mesh_Save save;
			// save.vtu("./grid/0/", rank, meshComb);
		// }
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




// void parMETIS_Graph_Partition(int nSize, vector<int>& procNoCell){

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
	// // options[METIS_OPTION_DBGLVL]=1;
	// // options[0] = 0;
	// options[0] = 0;
	// options[1] = 0;
	// options[2] = 0;
	
	// MPI_Comm comm;
	// MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	// int wgtflag=0;
	// int numflag=0;
	// real_t tpwgts[nSize*ncon];
	// for(int i=0;i<nSize*ncon;++i) 
		// tpwgts[i]=1.0/nSize;

	// // real_t ubvec[ncon];
	// // std::fill_n(ubvec, ncon, 1.02);
	
	// real_t ubvec = 1.02;
	
	
	// vector<int> temp_vtxdist(nSize);
    // MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	// // vector<int> vtxdist(nSize+1);
	// vector<int> vtxdist(nSize+1);
	// vtxdist[0] = 0;
	// for(int i=1; i<nSize+1; ++i){
		// vtxdist[i] = vtxdist[i-1] + static_cast<int>(temp_vtxdist[i-1]);
	// }
	
	
	// vector<int> xadj(ncells+1,0);
	// xadj[0] = 0;
	// int numt = 0;
	// for(auto& cell : mesh.cells){
		// int numt2 = 0;
		// for(auto& face : cell.faces()){
			// if(
			// (*face).getType() == MASCH_Face_Types::INTERNAL ||
			// (*face).getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		// }
		// ++numt;
		// xadj[numt] = xadj[numt-1] + static_cast<int>(numt2);
	// }
	
	
	// vector<int> adjncy(xadj[ncells],0);
	
	// vector<int> num_ncell(ncells,0);
	// int proc_num = 0;
	// for(auto& face : mesh.faces){
		// if(face.getType() == MASCH_Face_Types::INTERNAL) {
			
			// int tnum = xadj[face.iL] + num_ncell[face.iL];
			// adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iR);
			// ++num_ncell[face.iL];
			
			// tnum = xadj[face.iR] + num_ncell[face.iR];
			// adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iL);
			// ++num_ncell[face.iR];
			
		// }
		// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // int rank_neighbour = iRank_Recv[proc_num];
			// // // cout << proc_num << " " << iRank_Recv.size() << endl;
			// // adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// // ++num_ncell[face.owner];
			
			// // ++proc_num;
		// // }
	// }
	

	// ParMETIS_V3_PartKway(
		// vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, &wgtflag, &numflag,
		// &ncon, &nSize, tpwgts, &ubvec,
		// options, &objval, procNoCell.data(), &comm);
		
		
	// // idx_t vsize = NULL;
	// // real_t itr = 1000.0;
	// // ParMETIS_V3_AdaptiveRepart(
		// // vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, &wgtflag, &numflag,
		// // &ncon, &nSize, tpwgts, &ubvec, &itr, 
		// // options, &objval, procNoCell.data(), &comm);

		
	// // libSCOTCH_PartGraph(xadj, adjncy, &nSize, procNoCell);


	// // for(int i=0; i<ncells; ++i){
		// // procNoCell[i] = 1;
	// // }
	// // for(int i=0; i<ncells/2; ++i){
		// // procNoCell[i] = 0;
	// // }

	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	
	
// }





// // void parMETIS_Mesh_Partition(int nSize, vector<int>& procNoCell){

    // // int rank = MPI::COMM_WORLD.Get_rank(); 
    // // int size = MPI::COMM_WORLD.Get_size();

    // // if(rank == 0){
		// // cout << "┌────────────────────────────────────────────────────" << endl;
		// // cout << "| execute ParMETIS ... ";
	// // }
		
	// // int ncells = mesh.cells.size();
	// // int npoints = mesh.points.size();
	// // int nfaces = mesh.faces.size();
	// // int ncon=1;
	// // int ncommon=3;
	// // int objval;
	
	// // int options[METIS_NOPTIONS];
	// // METIS_SetDefaultOptions(options);
	// // options[METIS_OPTION_DBGLVL]=1;
	// // options[0] = 0;
	

	// // MPI_Comm comm;
	// // MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	
	// // int wgtflag=0;
	// // int numflag=0;
	// // real_t tpwgts[nSize*ncon];
	// // for(int i=0;i<nSize*ncon;++i) 
		// // tpwgts[i]=1.0/nSize;

	// // real_t ubvec[ncon];
	// // std::fill_n(ubvec, ncon, 1.05);
	
	
	
	// // int temp_vtxdist[nSize];
    // // MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	// // int vtxdist[nSize+1];
	// // vtxdist[0] = 0;
	// // for(int i=1; i<nSize+1; ++i){
		// // vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	// // }
	
	
	
	
	// // vector<int> xadj(ncells+1,0);
	// // xadj[0] = 0;
	// // for(int i=1; i<ncells+1; ++i){
		// // xadj[i] = xadj[i-1] + mesh.cells[i-1].points.size();
	// // }
	
	
	
	// // // vector<int> num_ncell(ncells,0);
	// // // vector<int> adjncy(xadj[ncells],0);
	// // vector<int> adjncy;
	// // int proc_num = 0;
	// // for(int i=0; i<ncells; ++i){
		// // for(auto& point : mesh.cells[i].points){
			// // adjncy.push_back(point);
		// // }
	// // }
		
	// // // for(auto& face : mesh.points){
		
		// // // if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			// // // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // // adjncy[tnum] = vtxdist[rank] + face.neighbour;
			// // // ++num_ncell[face.owner];
			
			// // // // cout << face.neighbour << endl;
			
			// // // tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			// // // adjncy[tnum] = vtxdist[rank] + face.owner;
			// // // ++num_ncell[face.neighbour];
			
		// // // }
		// // // // else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// // // // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // // // int rank_neighbour = iRank_Recv[proc_num];
			// // // // // cout << proc_num << " " << iRank_Recv.size() << endl;
			// // // // adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// // // // ++num_ncell[face.owner];
			
			// // // // ++proc_num;
		// // // // }
	// // // }
	

	// // int ncommonnodes = 3;

	// // // ParMETIS_V3_PartMeshKway(
		// // // vtxdist, xadj.data(), adjncy.data(), NULL, &wgtflag, &numflag,
		// // // &ncon, &ncommonnodes, &nSize, tpwgts, ubvec,
		// // // options, &objval, procNoCell.data(), &comm);

	// // // objval = adjncy.size()+1;

	// // vector<int> npart(npoints,0);

	// // METIS_PartMeshDual(
		// // &ncells, &npoints, xadj.data(), adjncy.data(), NULL, NULL,
		// // &ncommonnodes, &nSize, NULL, NULL,
		// // &objval, procNoCell.data(), npart.data());

	
	// // if(rank==0){
		// // cout << "-> completed" << endl;
		// // cout << "└────────────────────────────────────────────────────" << endl;
	// // }
	
	
	
	
	
// // }







// void partitionFromSerial(int nSize, vector<int>& procNoCell, vector<MASCH_Mesh>& newMesh){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// int ncells = static_cast<int>(mesh.cells.size());
	// int nfaces = static_cast<int>(mesh.faces.size());
	// int npoints = static_cast<int>(mesh.points.size());
	// int nboundary = static_cast<int>(mesh.boundaries.size());
	
	// vector<vector<int>> idBlockPoint(npoints,vector<int>(0,0)); // point block id (copies)
	// vector<int> nCellsLocal(nSize,0); // local total cells
	// vector<int> idCellLocal(ncells,0); // local cell id
	// // std::fill_n(nCellsLocal, nSize, 0);
	// for(int i=0; i<ncells; ++i) {
		// for(auto& j : mesh.cells[i].ipoints){
			// if(find(
			// idBlockPoint[j].begin(),idBlockPoint[j].end(),procNoCell[i])==
			// idBlockPoint[j].end()){
				// idBlockPoint[j].push_back( procNoCell[i] );
			// }
		// }
		// idCellLocal[i] = nCellsLocal[ procNoCell[i] ]++;
	// }
	
	// int nTotalLocalPointSize = 0; // total all # of local point
	// for(auto& item : idBlockPoint) 
		// nTotalLocalPointSize += item.size();

	// // point distribution (CSR format)
	// //
	// // / gP : global points / lP : local points / BL : block id /
	// //
	// //     gP[0]     gP[1]      gP[2]      gP[3]   gP[4] ...
	// // - - - - - - - - - - - - - - - - - - - - - - - - - ....
	// // |  lP[0,5]  |lP[6,9]| lP[10,15] |  ....   |     |
	// // |  BL[0,5]  |BL[6,9]| BL[10,15] |  ....   |     |
	// // |           |       |           |         |     |
	// // strPoints[0]|  strPoints[2]  strPoints[3] |  strPoints[5]  ...
	// //         strPoints[1]                  strPoints[4]
	// //
	// vector<int> nPointsLocal(nSize,0); // local total points
	// vector<int> idPointLocal(nTotalLocalPointSize,0); // local point id
	// vector<int> idBlockPointLocal(nTotalLocalPointSize,0); // local point block id
	// vector<int> strPoints(npoints+1,0); // start of each global point
	// int nIndex = 0;
	// for(int i=0; i<mesh.points.size(); ++i){
		// strPoints[i] = nIndex;
		// for(auto& item : idBlockPoint[i]){
			// int idBlock = item;
			// idPointLocal[nIndex] = nPointsLocal[ idBlock ]++;
			// idBlockPointLocal[nIndex] = idBlock;
			// ++nIndex;
		// }
	// }
	// strPoints[mesh.points.size()] = nIndex;
	
	// nPointsLocal.clear();
	// nPointsLocal.resize(nSize,0);
	// for(int i=0; i<nTotalLocalPointSize; ++i)
		// ++nPointsLocal[ idBlockPointLocal[i] ];


	// vector<int> nDisplPoint(nSize,0);
	// for(int i=0; i<mesh.points.size(); ++i){
		// for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			// int idBlock = idBlockPointLocal[j];
			// int idLocal = idPointLocal[j];
			
			// newMesh[idBlock].addPoint();
			// newMesh[idBlock].points.back().x = mesh.points[i].x;
			// newMesh[idBlock].points.back().y = mesh.points[i].y;
			// newMesh[idBlock].points.back().z = mesh.points[i].z;
			
			// // MPI 컨넥트 포인트들 저장 (procNo, local point id)
			// vector<pair<int,int>> tmp_idBlock_idLocal;
			// for(int k=strPoints[i]; k<strPoints[i+1]; ++k){
				// if(j==k) continue;
				// tmp_idBlock_idLocal.push_back(
				// make_pair(idBlockPointLocal[k],idPointLocal[k]));
			// }
			// newMesh[idBlock].points.back().connPoints =
			// tmp_idBlock_idLocal;
			
			
			
			// ++nDisplPoint[idBlock];
		// }
		
		
	// }
	
	
	// // for(auto& item : newMesh[2].points){
		// // if(item.connPoints.size()>0){
		// // cout << "========" << endl;
		// // }
		// // for(auto& i : item.connPoints){
			// // cout << i.first << " " << i.second << endl;
		// // }
	// // }
	
	// // face setting
	// vector<int> nFacesLocal(nSize,0); // total local faces
	// for(auto& face : mesh.faces){
		// int idBlockOwner = procNoCell[face.iL];
		
		// ++nFacesLocal[idBlockOwner];
		
		// if(face.getType() == MASCH_Face_Types::INTERNAL){
			// int idBlockNeighbour = procNoCell[face.iR];
			// if(idBlockOwner != idBlockNeighbour){
				// face.setType(MASCH_Face_Types::PROCESSOR); // set processor face
				// ++nFacesLocal[idBlockNeighbour];
			// }
		// }
	// }
	
	
	// vector<int> idFaceLocal(nSize,0);  // temporary local face id
	
	// // internal faces
	// vector<vector<int>> idFacePoint_Send(nSize,vector<int>(0,0));
	
	// for(auto& face : mesh.faces){
		// if(face.getType() == MASCH_Face_Types::INTERNAL){
			// int idBlock = procNoCell[face.iL];
			// ++idFaceLocal[idBlock];
			
			// vector<int> idFacePoint; // face's points id
			// for(auto& i : face.ipoints){
				// for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					// if(idBlockPointLocal[j] == idBlock){
						// idFacePoint.push_back( static_cast<int>(idPointLocal[j]) );
						// break;
					// }
				// }
			// }
			// // save face's points size , local points id
			// newMesh[idBlock].addFace();
			// for(int i=0; i<idFacePoint.size(); ++i){
				// newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
			// }
		// }
	// }
	
	
	// // 페이스 타입 BC
	// vector<int> faceTypeBC(nfaces,-1);
	// for(int ibcs=0; ibcs<mesh.boundaries.size(); ++ibcs){
		// int str = mesh.boundaries[ibcs].startFace;
		// int end = str + mesh.boundaries[ibcs].nFaces;
		// for(int i=str; i<end; ++i){
			// faceTypeBC[i] = ibcs;
		// }
	// }
	
	
	// // boundary faces
	// // local boundary face size
	// vector<int> nFacesBoundaryLocal(nSize,0); 
	// // local each boundary face size
	// vector<vector<int>> nFacesEachBoundaryLocal(nSize,vector<int>(nboundary,0)); 
	// // local each boundary face start
	// vector<vector<int>> nStartFaceEachBoundaryLocal(nSize,vector<int>(nboundary,0)); 
	// vector<vector<bool>> boolFaceEachBoundaryLocal(nSize,vector<bool>(nboundary,false)); 
	// for(int id=0; id<mesh.faces.size(); ++id){
		// auto& face = mesh.faces[id];
		
		// if(face.getType() == MASCH_Face_Types::BOUNDARY){
			// int idBlock = procNoCell[face.iL];
			// ++nFacesBoundaryLocal[ idBlock ];
			
			// ++nFacesEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ];
			// if(boolFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] == false){
				// nStartFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] = idFaceLocal[ idBlock ];
				// boolFaceEachBoundaryLocal[ idBlock ][ faceTypeBC[id] ] = true;
			// }
			// ++idFaceLocal[ idBlock ];

			// // wirte bc
			// vector<int> idFacePoint;
			// for(auto& i : face.ipoints){
				// for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					// if(idBlockPointLocal[j] == idBlock){
						// idFacePoint.push_back( idPointLocal[j] );
						// break;
					// }
				// }
			// }
			// // save face's points size , local points id
			// // idFacePoint_Send[idBlock].push_back(face.points.size());
			// newMesh[idBlock].addFace();
			// for(int i=0; i<idFacePoint.size(); ++i){
				// newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
				// // idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			// }
			
		// }
		
	// }
	
	
	// // processor faces
	// // local processor face size
	// vector<int> nFacesProcessorLocal(nSize,0); 
	 // // sending size from each processor to each processor 
	// vector<vector<int>> sendCountsProcs(nSize,vector<int>(nSize,0));
	// vector<int> idFacesProcessor; // local processor face id
	// int temp_num_proc_face = 0;
	// for(auto& face : mesh.faces){
		// if(face.getType() == MASCH_Face_Types::PROCESSOR){
			// idFacesProcessor.push_back(temp_num_proc_face);
			
			// int idBlockOwner = procNoCell[face.iL];
			// int idBlockNeighbour = procNoCell[face.iR];
			
			// ++nFacesProcessorLocal[idBlockOwner];
			// ++nFacesProcessorLocal[idBlockNeighbour];
			// ++sendCountsProcs[idBlockOwner][idBlockNeighbour];
			// ++sendCountsProcs[idBlockNeighbour][idBlockOwner];
			
		// }
		// ++temp_num_proc_face;
	// }
	
	
	// for(int ip=0; ip<nSize; ++ip){
		// for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			// int k = idFacesProcessor[jp];
			// int m = procNoCell[ mesh.faces[k].iL ];
			// int n = procNoCell[ mesh.faces[k].iR ];
			// if(n==ip) {
				// int idBlock = m;
				// vector<int> idFacePoint;
				// for(auto& i : mesh.faces[k].ipoints){
					// for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						// if(idBlockPointLocal[j] == idBlock){
							// idFacePoint.push_back( idPointLocal[j] );
							// break;
						// }
					// }
				// }
			
				// // idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				// newMesh[idBlock].addFace();
				// for(int i=0; i<idFacePoint.size(); ++i){
					// // idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					// newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
					
				// }
				
			// }
			// else if(m==ip) {
				// int idBlock = n;
				// vector<int> idFacePoint;
				// for(auto& i : mesh.faces[k].ipoints){
					// for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						// if(idBlockPointLocal[j] == idBlock){
							// idFacePoint.push_back( idPointLocal[j] );
							// break;
						// }
					// }
				// }
			
				// // idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				// newMesh[idBlock].addFace();
				// std::reverse(idFacePoint.begin()+1,idFacePoint.end());
				// for(int i=0; i<idFacePoint.size(); ++i){
					// // idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					// newMesh[idBlock].faces.back().ipoints.push_back( static_cast<int>(idFacePoint[i]) );
				// }
					
			// }
		// }
	// }
	
	
	
	
	// // write owner of internal face
	// vector<vector<int>> idOwnerLocal(nSize,vector<int>(0,0));
	// for(auto& face : mesh.faces){
		// if(face.getType() == MASCH_Face_Types::INTERNAL){
			// int j = face.iL;
			// int i = procNoCell[j];
			// idOwnerLocal[i].push_back(idCellLocal[j]);
		// }
	// }
	// for(auto& face : mesh.faces){
		// if(face.getType() == MASCH_Face_Types::BOUNDARY){
			// int j = face.iL;
			// int i = procNoCell[j];
			// idOwnerLocal[i].push_back(idCellLocal[j]);
		// }
	// }
	// for(int i=0; i<nSize; ++i){
		// for(int j=0; j<idFacesProcessor.size(); ++j){
			// int k = idFacesProcessor[j];
			// int m = procNoCell[ mesh.faces[k].iL ];
			// int n = procNoCell[ mesh.faces[k].iR ];
			
			// if(n==i) {
				// idOwnerLocal[m].push_back( idCellLocal[ mesh.faces[k].iL ] );
			// }
			// else if(m==i) {
				// idOwnerLocal[n].push_back( idCellLocal[ mesh.faces[k].iR ] );
			// }
		// }
	// }
	
	// for(int ip=0; ip<nSize; ++ip){
		// int i=0;
		// for(auto& iL : idOwnerLocal[ip]){
			// newMesh[ip].faces[i].iL = iL;
			// ++i;
		// }
	// }
	
	// // write neighbour of internal face
	// vector<vector<int>> idNeighbourLocal(nSize,vector<int>(0,0));
	// for(auto& face : mesh.faces){
		// if(
		// face.getType() == MASCH_Face_Types::INTERNAL
		// ){
			// int j = face.iR;
			// int i = procNoCell[j];
			// idNeighbourLocal[i].push_back(idCellLocal[j]);
		// }
	// }
	
	// for(int ip=0; ip<nSize; ++ip){
		// int i=0;
		// for(auto& iR : idNeighbourLocal[ip]){
			// newMesh[ip].faces[i].iR = iR;
			// // newMesh[ip].faces[i].thereR = true;
			// newMesh[ip].faces[i].setType(MASCH_Face_Types::INTERNAL);
			// ++i;
		// }
	// }
	

	// // write of boundary faces
	// // vector<int> nFaceBoundaryFaceLocal_Send(nSize*mesh.boundary.size(),0);
	// // vector<int> nStartBoundaryFaceLocal_Send(nSize*mesh.boundary.size(),0);
	// // int temp_bound=0;
	
	
	// for(int ip=0; ip<nSize; ++ip){
		// for(int ibcs=0; ibcs<mesh.boundaries.size(); ++ibcs){
			// newMesh[ip].addBoundary();
			
			// string bcnames = mesh.boundaries[ibcs].name;
			
			// bcnames.erase(
			// std::find_if(bcnames.rbegin(), 
			// bcnames.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), 
			// bcnames.end());
			
			// bcnames.erase(
			// bcnames.begin(), 
			// std::find_if(bcnames.begin(), bcnames.end(), 
			// std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			// newMesh[ip].boundaries.back().name = bcnames;
			
			// if( nFacesEachBoundaryLocal[ip][ibcs] > 0) {
				
				// newMesh[ip].boundaries.back().nFaces = nFacesEachBoundaryLocal[ip][ibcs];
				// newMesh[ip].boundaries.back().startFace = nStartFaceEachBoundaryLocal[ip][ibcs];
				// // newMesh[ip].boundaries.back().thereR = false;
				// newMesh[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			// }
			// else{
				// newMesh[ip].boundaries.back().nFaces = 0;
				// newMesh[ip].boundaries.back().startFace = 0;
				// // newMesh[ip].boundaries.back().thereR = false;
				// newMesh[ip].boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
			// }
		// } 
	// }
	
	
	
	
	// // write of processor faces
	// for(int ip=0; ip<nSize; ++ip){
		// int n = nFacesLocal[ip];
		// for(int jp=0; jp<nSize; ++jp){
			// n -= sendCountsProcs[ip][jp];
		// }
		
		// for(int jp=0; jp<nSize; ++jp){
			// if( sendCountsProcs[ip][jp] > 0) {
				// newMesh[ip].addBoundary();
				
				// string bcnames = "procBoundary" + to_string(ip) + "to" + to_string(jp);
				
				// newMesh[ip].boundaries.back().name = bcnames;
				
				// newMesh[ip].boundaries.back().nFaces = sendCountsProcs[ip][jp];
				// newMesh[ip].boundaries.back().startFace = n;
				// newMesh[ip].boundaries.back().myProcNo = ip;
				// newMesh[ip].boundaries.back().rightProcNo = jp;
				// // newMesh[ip].boundaries.back().thereR = true;
				// newMesh[ip].boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
				
				// n += sendCountsProcs[ip][jp];
			// }
		// }
	// }
	
	
	// for(int ip=0; ip<nSize; ++ip){
		// int maxBCNum = newMesh[ip].boundaries.size()-1;
		// newMesh[ip].boundaries[maxBCNum].startFace = newMesh[ip].faces.size()-newMesh[ip].boundaries[maxBCNum].nFaces;
		// for(int i=maxBCNum-1; i>=0; --i){
			// newMesh[ip].boundaries[i].startFace = newMesh[ip].boundaries[i+1].startFace-newMesh[ip].boundaries[i].nFaces;
		// }
	// }
	
	// // //==========================================
	
	// vector<int> nPointsInf;
	// vector<int> nFacesInf;
	// vector<int> nCellsInf;
	// vector<int> nFacesIntInf;
	// vector<int> nFacesBCInf;
	// vector<int> nFacesProcInf;
	// vector<int> nTriangleInf;
	// vector<int> nQuadrangleInf;
	// vector<int> nPolygonInf;
	// vector<int> nTetrahedronInf;
	// vector<int> nHexahedronInf;
	// vector<int> nPrismInf;
	// vector<int> nPyramidInf;
	// vector<int> nPolyhedronInf;
	
	// int nbcs=0;
	
	// for(int ip=0; ip<nSize; ++ip){

		// newMesh[ip].check();
		// newMesh[ip].setFaceTypes();
		// newMesh[ip].buildCells();
		// newMesh[ip].connectFacetoPointsCells();
		// newMesh[ip].connectCelltoFaces();
		// newMesh[ip].connectCelltoPoints();
		// newMesh[ip].setCountsProcFaces();
		// newMesh[ip].setDisplsProcFaces();
		
		// nPointsInf.push_back(newMesh[ip].points.size());
		// nFacesInf.push_back(newMesh[ip].faces.size());
		// nCellsInf.push_back(newMesh[ip].cells.size());
	
		// nbcs=0;
		// int nFacesBC = 0;
		// int nFacesProc = 0;
		// for(auto& boundary : newMesh[ip].boundaries){
			
			// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// nFacesBC += boundary.nFaces;
				// ++nbcs;
			// }
			// else if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
				// nFacesProc += boundary.nFaces;
			// }
			// else{
				// cout << "#ERROR : 868" << endl;
			// }
		// }
		
		
		// // faces type
		// int nFacesInt = 0;
		// int nTriangle = 0;
		// int nQuadrangle = 0;
		// int nPolygon = 0;
		// for(auto& face : newMesh[ip].faces){
			
			// if(face.getType() == MASCH_Face_Types::INTERNAL){
				// ++nFacesInt;
			// }
			
			// if(face.ipoints.size() == 3){
				// ++nTriangle;
			// }
			// else if(face.ipoints.size() == 4){
				// ++nQuadrangle;
			// }
			// else{
				// ++nPolygon;
			// }
		// }
		
		
		
		// // cells type
		// int nTetrahedron = 0;
		// int nHexahedron = 0;
		// int nPrism = 0;
		// int nPyramid = 0;
		// int nPolyhedron = 0;
		// for(auto& cell : newMesh[ip].cells){
			// if( cell.ipoints.size() == 4 && cell.ifaces.size() == 4 ){
				// ++nTetrahedron;
			// }
			// else if( cell.ipoints.size() == 5 && cell.ifaces.size() == 5 ){
				// ++nPyramid;
			// }
			// else if( cell.ipoints.size() == 6 && cell.ifaces.size() == 5 ){
				// ++nPrism;
			// }
			// else if( cell.ipoints.size() == 8 && cell.ifaces.size() == 6 ){
				// ++nHexahedron;
			// }
			// else {
				// ++nPolyhedron;
			// }
		// }
		
		// nFacesIntInf.push_back(nFacesInt);
		// nFacesBCInf.push_back(nFacesBC);
		// nFacesProcInf.push_back(nFacesProc);
		// nTriangleInf.push_back(nTriangle);
		// nQuadrangleInf.push_back(nQuadrangle);
		// nPolygonInf.push_back(nPolygon);
		// nTetrahedronInf.push_back(nTetrahedron);
		// nHexahedronInf.push_back(nHexahedron);
		// nPrismInf.push_back(nPrism);
		// nPyramidInf.push_back(nPyramid);
		// nPolyhedronInf.push_back(nPolyhedron);
	// }
	

    // if(rank == 0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| partitioned MPI size : " << nSize << endl;
		// cout << "| points size : ";
		// for(auto& i : nPointsInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = mesh.faces.size();
		// cout << "| faces size : ";
		// for(auto& i : nFacesInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = mesh.cells.size();
		// cout << "| cells size : ";
		// for(auto& i : nCellsInf) cout << i << " | ";
		// cout << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// // gatherValue = nFacesInt;
		// cout << "| internal faces size : ";
		// for(auto& i : nFacesIntInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nFacesBC;
		// cout << "| boundary faces size : ";
		// for(auto& i : nFacesBCInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nFacesProc;
		// cout << "| processor faces size : ";
		// for(auto& i : nFacesProcInf) cout << i << " | ";
		// cout << endl;
		// cout << "| boundary types : " << nbcs << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// // gatherValue = nTriangle;
		// cout << "| Triangle faces : ";
		// for(auto& i : nTriangleInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nQuadrangle;
		// cout << "| Quadrangle faces : ";
		// for(auto& i : nQuadrangleInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nPolygon;
		// cout << "| Polygon faces : ";
		// for(auto& i : nPolygonInf) cout << i << " | ";
		// cout << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// // gatherValue = nTetrahedron;
		// cout << "| Tetrahedron cells : ";
		// for(auto& i : nTetrahedronInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nHexahedron;
		// cout << "| Hexahedron cells : ";
		// for(auto& i : nHexahedronInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nPrism;
		// cout << "| Prism cells : ";
		// for(auto& i : nPrismInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nPyramid;
		// cout << "| Pyramid cells : ";
		// for(auto& i : nPyramidInf) cout << i << " | ";
		// cout << endl;
		// // gatherValue = nPolyhedron;
		// cout << "| Polyhedron cells : ";
		// for(auto& i : nPolyhedronInf) cout << i << " | ";
		// cout << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
    // }
	
	
// }








