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

// #include "../mesh.h" 
// #include "../geometric.h" 
#include "./partition.h" 
#include "../../test1/mpi.h"
// #include "../../mpi/mpi.h" 

// void parMETIS_Mesh_Partition(int nBlocks, vector<int>& idBlockCell);
// void loadOnlyMeshVtu(string folder, MASCH_Mesh &mesh);
// void MASCH_Gatherv(vector<int>& my_value, vector<int>& value, vector<int>& displs);
// void MASCH_Gatherv(vector<double>& my_value, vector<double>& value, vector<int>& displs);



void MASCH_Mesh_Partition::partitionFromSerial(
int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh, vector<MASCH_Mesh>& newMesh){
	
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










// void MASCH_Mesh_Partition::loadOnlyMeshVtu(
// string folder, MASCH_Mesh &mesh){
	
	

	// int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	// int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load vtu data files ... ";
	// }
		
	// string saveFolderName = folder;
	// string saveFileName = "plot";
	// string saveRankName = to_string(rank);
	
	// ifstream inputFile;
	// string openFileName;
	
	// // points
	// openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	// inputFile.open(openFileName);
	// if(inputFile.fail()){
		// cerr << "Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	
	// string nextToken;
	// bool boolCompress = false;
	// while(getline(inputFile, nextToken)){
		// if( nextToken.find("VTKFile") != string::npos ){
			// if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				// boolCompress = true;
			// }
			// break;
		// }
	// }
	// inputFile.clear();
	
	
	
	
	// // vector<string> volFracName;
	// // volFracName.push_back(controls.name[controls.VF[0]]);
	
	// // vector<string> massFracName;
	// // massFracName.push_back(controls.name[controls.MF[0]]);
	
	// SEMO_Mesh_Load load;
	
	// vector<double> NodeCoordinates;
	// load.loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	// // cout << NodeCoordinates.size() << endl;
	
	// vector<int> connectivity;
	// load.loadDatasAtVTU(inputFile, "connectivity", connectivity);
	
	// vector<int> offsets;
	// load.loadDatasAtVTU(inputFile, "offsets", offsets);
	
	// vector<int> faces;
	// load.loadDatasAtVTU(inputFile, "faces", faces);
	
	// vector<int> faceoffsets;
	// load.loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
	
	// vector<int> owner;
	// load.loadDatasAtVTU(inputFile, "owner", owner);
	
	// vector<int> neighbour;
	// load.loadDatasAtVTU(inputFile, "neighbour", neighbour);
	
	// vector<string> bcName;
	// load.loadDatasAtVTU(inputFile, "bcName", bcName);
	
	// vector<int> bcStartFace;
	// load.loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
	
	// vector<int> bcNFaces;
	// load.loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
	
	// vector<int> bcNeighbProcNo;
	// load.loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);
	
	// vector<int> connPoints;
	// load.loadDatasAtVTU(inputFile, "connPoints", connPoints);


	// inputFile.close();
	
	
	
	
	

	
	// for(int i=0; i<NodeCoordinates.size()/3; ++i){
		// mesh.addPoint();
		// mesh.points.back().x = NodeCoordinates[i*3+0];
		// mesh.points.back().y = NodeCoordinates[i*3+1];
		// mesh.points.back().z = NodeCoordinates[i*3+2];
	// }
	// NodeCoordinates.clear();
	
	
	// int n=0;
	// int ncells=0;
	// mesh.faces.clear();
	// for(auto& i : owner){
		// mesh.addFace();
		// mesh.faces.back().iL = i;
		// ncells = max(ncells, mesh.faces.back().iL);
	// }
	// owner.clear();
	
	// // cout << ncells << endl;
	
	// mesh.cells.clear();
	// for(int i=0; i<ncells+1; ++i){
		// mesh.addCell();
	// }
	
	
	
	// n=0;
	// for(auto& i : neighbour){
		// mesh.faces[n].iR = i;
		// mesh.faces[n].setType(MASCH_Face_Types::INTERNAL);
		// ++n;
	// }
	// neighbour.clear();
	
	
	// int m=0;
	// n=0;
	// for(auto& i : offsets){
		// for(int j=n; j<i; ++j){
			// int point = connectivity[j];
			// // cout << m << endl;
			// mesh.cells[m].ipoints.push_back( static_cast<int>(point) );
			// // if(rank==1) cout << point << endl;
		// }
		
		// // if(rank==0 && mesh.cells[m].ipoints.size()!=8) cout << mesh.cells[m].ipoints.size() << endl;
		
		// n=i;
		// ++m;
	// }
	
	
	
	// n=0;
	// int nFacesInt=0;
	// for(auto& face : mesh.faces){
		// // if(face.iR != -1){
		// if(face.getType() == MASCH_Face_Types::INTERNAL){
			// mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			// mesh.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
			// ++nFacesInt;
		// }
		// else{
			// mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
		// }
		// ++n;
	// }
	
	
	// m=0;
	// n=0;
	// for(auto& i : faceoffsets){
		// // if(faces[n]>5) cout << faces[n] << endl;
		// int N=0;
		// int face_size = faces[m+N];
		// for(int j=0; j<face_size; ++j){
			// int face = mesh.cells[n].ifaces[j];
			// ++N;
			// int point_size = faces[m+N];
			// for(int k=0; k<point_size; ++k){
				// ++N;
				// int point = faces[m+N];
				// if(mesh.faces[ face ].ipoints.size() == point_size) continue;
				// mesh.faces[ face ].ipoints.push_back( static_cast<int>(point) );
				// // if(rank==1) cout << point << endl;
			// }
		// }
		// m=i;
		// ++n;
	// }
	// faces.clear();
	// faceoffsets.clear();
	
	
	
	// n=0;
	// for(auto& startFace : bcStartFace){
		
		// mesh.addBoundary();
		
		// string tmp_bcName = bcName[n];
		// tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
		// tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
		// tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
		// tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
		// mesh.boundaries.back().name = tmp_bcName;
		// mesh.boundaries.back().startFace = bcStartFace[n];
		// mesh.boundaries.back().nFaces = bcNFaces[n];
		// if(bcNeighbProcNo[n]<0){
			// mesh.boundaries.back().rightProcNo = 0;
			// mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
		// }
		// else{
			// mesh.boundaries.back().rightProcNo = bcNeighbProcNo[n];
			// mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
		// }
		// mesh.boundaries.back().myProcNo = rank;
		// // if(bcNeighbProcNo[n] < 0){
			// // mesh.boundaries.back().myProcNo = -1;
		// // }
		
		// ++n;
		
	// }
	
	// int maxBCNum = mesh.boundaries.size()-1;
	// mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
	// for(int i=maxBCNum-1; i>=0; --i){
		// mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
	// }
	
	// bcName.clear();
	// bcStartFace.clear();
	// bcNFaces.clear();
	// bcNeighbProcNo.clear();
	
		
	// for(int i=0; i<connPoints.size(); ++i){
		// // cout << mesh.points.size() << " " << connPoints[i] << " " << connPoints[i+1] << " " << connPoints[i+2] << " " << endl;
		// mesh.points[connPoints[i]].connPoints.push_back(
		// make_pair(connPoints[i+1],connPoints[i+2])
		// );
		// ++i; ++i;
	// }
	
	
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	// for(auto& i : mesh.cells){
		
		// if(rank==0 && i.ipoints.size()!=8) cout << i.ipoints.size() << endl;
		// if(rank==0 && i.ifaces.size()!=6) cout << i.ifaces.size() << endl;
		
	// }
	
	// mesh.check();
	// mesh.setFaceTypes();
	// // mesh.buildCells();
	// // mesh.connectFacetoPointsCells();
	// // mesh.connectCelltoFaces();
	// // mesh.connectCelltoPoints();
	// mesh.setCountsProcFaces();
	// mesh.setDisplsProcFaces();
	// mesh.cellsGlobal();
	
	
	// mesh.informations();

		
	
	
	
// }










// void mpi.Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs);
// void mpi.Alltoallv(vector<vector<int>>& inp_send_value, vector<int>& recv_value, vector<int>& displs);
// void mpi.Alltoallv(vector<vector<double>>& inp_send_value, vector<vector<double>>& recv_value);
// void mpi.Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value);




void MASCH_Mesh_Partition::combine(vector<int>& idBlockCell, MASCH_Mesh &mesh){
	
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	MASCH_MPI mpi;
	
	
	// vector<MASCH_Mesh> meshNew;
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
		
		// send_idBlockCell.size();
		
		// if(send_idBlockCell.size()==0){
			// send_idBlockCell.resize(1);
			// send_rank.resize(1);
		// }
		
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
	
	
	
	{
		
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
		
		
		vector<pair<int,int>> send_localCell_proc_id;
		vector<vector<int>> send_localCell_id(size);
		vector<int> send_localCell_n(size,0);
		for(int i=0, ip=0; i<mesh.cells.size(); ++i){
			int proc = idBlockCell[i];
			int tmp_nCell = send_localCell_n[proc]++;
			send_localCell_proc_id.push_back(make_pair(proc,tmp_nCell));
			send_localCell_id[proc].push_back(tmp_nCell);
		}
		vector<vector<int>> recv_localCell_id;
		mpi.Alltoallv(send_localCell_id, recv_localCell_id);
		vector<int> nCells_local(size+1,0);
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = recv_localCell_id[ip].size();
		}
		for(int ip=0; ip<size; ++ip){
			nCells_local[ip+1] = nCells_local[ip+1] + nCells_local[ip];
		}
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
		
		
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
		vector<vector<int>> recv_localPoint_connId;
		mpi.Alltoallv(send_localPoint_connId, recv_localPoint_connId);
		vector<vector<int>> send_localPoint_toProc_toId(size);
		vector<vector<int>> send_globalPoint_toProc_toId(size);
		for(int ip=0; ip<size; ++ip){
			for(int i=0; i<recv_localPoint_connId[ip].size(); ++i){
				int ipoint = recv_localPoint_connId[ip][i];
				int right_i = recv_localPoint_connId[ip][++i];
				int tmp_size = recv_localPoint_connId[ip][++i];
				auto& point = send_localPoint_proc[ipoint];
				for(int j=0; j<tmp_size; ++j){
					int tmp_proc = recv_localPoint_connId[ip][++i];
					int tmp_id_local = recv_localPoint_connId[ip][++i];
					for(auto& [proc, id] : send_localPoint_proc_id[ipoint]){
						if(proc==tmp_proc){
							send_localPoint_toProc_toId[proc].push_back(id); // proc 블록 에서 자신의 local id  
							send_localPoint_toProc_toId[proc].push_back(ip); // proc 블록 에서 자신과 중복한 포인트의 local proc
							send_localPoint_toProc_toId[proc].push_back(tmp_id_local); // proc 블록 에서 자신과 중복한 포인트의 local id
						}
						// connPoints 에 저장할 것
						else{
							send_globalPoint_toProc_toId[proc].push_back(id); // proc 블록 에서 자신의 local id  
							send_globalPoint_toProc_toId[proc].push_back(tmp_proc); // right 블록의 global proc 넘버
							send_globalPoint_toProc_toId[proc].push_back(ip); // right 블록 에서 자신과 중복한 포인트의 local proc
							send_globalPoint_toProc_toId[proc].push_back(tmp_id_local); // right 블록의 자신과 중복한 포인트의 local id 넘버
						}
					}
				}
			}
		}
		vector<vector<int>> recv_localPoint_toProc_toId;
		mpi.Alltoallv(send_localPoint_toProc_toId, recv_localPoint_toProc_toId);
		vector<vector<int>> recv_globalPoint_toProc_toId;
		mpi.Alltoallv(send_globalPoint_toProc_toId, recv_globalPoint_toProc_toId);
		
		vector<vector<double>> recv_localPoint_xyz;
		mpi.Alltoallv(send_localPoint_xyz, recv_localPoint_xyz);
		
		
		vector<int> nPoints_local(size+1,0);
		for(int ip=0, iter_glob=0; ip<size; ++ip){
			int tmp_size = recv_localPoint_xyz[ip].size()/3;
			nPoints_local[ip+1] = tmp_size;
			for(int i=0, iter=0; i<tmp_size; ++i){
				meshComb.addPoint();
				meshComb.points.back().x = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().y = recv_localPoint_xyz[ip][iter++];
				meshComb.points.back().z = recv_localPoint_xyz[ip][iter++];
				
				// points_id_local2global[ip].push_back(iter_glob++);
			}
			// int str = iter_glob-tmp_size;
			// for(int i=0, iter=0; i<recv_globalPoint_toProc_toId[ip].size()/4; ++i){
				// int my_id_local = recv_globalPoint_toProc_toId[ip][iter++];
				// int rightProc_global = recv_globalPoint_toProc_toId[ip][iter++];
				// int rightProc_local = recv_globalPoint_toProc_toId[ip][iter++];
				// int rightId_local = recv_globalPoint_toProc_toId[ip][iter++];
				// int id_global = str + my_id_local;
				// // cout << rightProc_global <<endl;
				// // meshComb.points[id_global].connPoints.push_back(make_pair(rightProc_global,));
			// }
		}
		for(int ip=0; ip<size; ++ip){
			nPoints_local[ip+1] = nPoints_local[ip+1] + nPoints_local[ip];
		}
		vector<bool> deletePoints(nPoints_local[size],false);
		vector<vector<int>> vec_rightProcNo_local(size);
		vector<vector<int>> vec_rightId_local(size);
		{
			for(int ip=0; ip<size; ++ip){
				int str = nPoints_local[ip];
				int end = nPoints_local[ip+1];
				for(int i=str; i<end; ++i){
					vec_rightProcNo_local[ip].push_back(-1);
					vec_rightId_local[ip].push_back(-1);
				}
				for(int i=0, iter=0; i<recv_localPoint_toProc_toId[ip].size()/3; ++i){
					int my_id_local = recv_localPoint_toProc_toId[ip][iter++];
					int my_id_global = nPoints_local[ip] + my_id_local;
					int rightProc_local = recv_localPoint_toProc_toId[ip][iter++];
					int rightId_local = recv_localPoint_toProc_toId[ip][iter++];
					int rightId_global = nPoints_local[rightProc_local] + rightId_local;
					
					vec_rightProcNo_local[ip][my_id_local] = rightProc_local;
					vec_rightId_local[ip][my_id_local] = rightId_local;
					// 중복 및 삭제
					if(rightProc_local<ip){
						deletePoints[my_id_global] = true;
					}
				}
			}
		}
		vector<int> tmp_point_global_id(size+1,0);
		vector<vector<int>> points_id_local2global(size);
		vector<vector<int>> delete_total_front_points(size);
		for(int ip=0, iter=0, before_nd=0, tmp_glob=0; ip<size; ++ip){
			int str = nPoints_local[ip];
			int end = nPoints_local[ip+1];
			for(int i=str; i<end; ++i){
				points_id_local2global[ip].push_back(tmp_glob);
				delete_total_front_points[ip].push_back(before_nd);
				if(deletePoints[i]==false) {
					++tmp_glob;
				}
				else{
					++before_nd;
				}
			}
			tmp_point_global_id[ip+1] = tmp_glob;
		}
		
		{
			for(int ip=0; ip<size; ++ip){
				int str = nPoints_local[ip];
				int end = nPoints_local[ip+1];
				// for(int i=str, iter=0; i<end; ++i){
				for(int i=0, iter=0; i<recv_localPoint_toProc_toId[ip].size()/3; ++i){
					int my_id_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int my_id_global = nPoints_local[ip] + my_id_local;
					int rightProc_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_local = recv_localPoint_toProc_toId[ip].at(iter++);
					int rightId_global = nPoints_local[rightProc_local] + rightId_local;
					
					// int new_id_global = tmp_point_global_id[rightProc_local] + rightId_local - 
											// delete_total_front_points[rightProc_local][rightId_local];
					// 중복 및 삭제
					if(rightProc_local<ip){
						bool continueCalc = true;
						int new_id_global = -1;
						while(continueCalc){
							new_id_global = points_id_local2global[rightProc_local][rightId_local];
							my_id_global = nPoints_local[rightProc_local] + rightId_local;
							int tmp_rightProc_local = vec_rightProcNo_local[rightProc_local][rightId_local];
							int tmp_rightId_local = vec_rightId_local[rightProc_local][rightId_local];
							rightProc_local = tmp_rightProc_local;
							rightId_local = tmp_rightId_local;
							if(deletePoints[my_id_global]==false) continueCalc=false;
						}
											
						// if(deletePoints[my_id_global]!=true) cout << "ERROR33" << endl;
						// if(deletePoints[rightId_global]==true) cout << "ERROR34" << endl;
						// points_id_local2global[ip][rightId_local] = rightId_global;
						points_id_local2global[ip][my_id_local] = new_id_global;
						// deletePoints[my_id_global] = true;
					}
				}
			}			
		}
		for(int ip=0; ip<size; ++ip){
			int str = nPoints_local[ip];
			int end = nPoints_local[ip+1];
			// for(int i=str, iter=0; i<end; ++i){
			for(int i=0, iter=0; i<recv_localPoint_toProc_toId[ip].size()/3; ++i){
				int my_id_local = recv_localPoint_toProc_toId[ip].at(iter++);
				int my_id_global = nPoints_local[ip] + my_id_local;
				int rightProc_local = recv_localPoint_toProc_toId[ip].at(iter++);
				int rightId_local = recv_localPoint_toProc_toId[ip].at(iter++);
				int rightId_global = nPoints_local[rightProc_local] + rightId_local;
				
				// int new_id_global = tmp_point_global_id[rightProc_local] + rightId_local - 
										// delete_total_front_points[rightProc_local][rightId_local];
				// 중복 및 삭제
				if(rightProc_local<ip){
					int new_id_global = points_id_local2global[rightProc_local][rightId_local];
										
					if(deletePoints[my_id_global]!=true) cout << "ERROR33" << endl;
					if(deletePoints[rightId_global]==true) cout << "ERROR34" << endl;
					// points_id_local2global[ip][rightId_local] = rightId_global;
					points_id_local2global[ip][my_id_local] = new_id_global;
					// deletePoints[my_id_global] = true;
				}
			}
		}
		
		// for(int ip=0, iter=0, nDelete=0; ip<size; ++ip){
			// for(int i=0; i<points_id_local2global[ip].size(); ++i){
				// if(deletePoints[iter++]==true) {
					// points_id_local2global[ip][i] -= (nDelete++);
				// }
			// }
		// }
		
		// cout << meshComb.points.size() << endl;
		
		// 포인트 삭제
		{
			int numN = 0;
			meshComb.points.erase( std::remove_if( meshComb.points.begin(), meshComb.points.end(), 
				[&deletePoints, &numN](MASCH_Point const& v) { 
				return deletePoints[numN++]; 
				}), meshComb.points.end());
		}
		
		
		//======================================
		
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
					
					send_localFace_PR2PR_iL[ipL].push_back(iL_local);
					send_globalFace_PR2PR_toProc[ipL].push_back(ipR);
					send_localFace_PR2PR_ipoints[ipL].push_back(face.ipoints.size());
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
					}
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
		vector<vector<int>> recv_localFace_PR2PR_iL;
		vector<vector<int>> recv_globalFace_PR2PR_toProc;
		vector<vector<int>> recv_localFace_PR2PR_ipoints;
		
		
		mpi.Alltoallv(send_localFace_IN2IN_iL, recv_localFace_IN2IN_iL);
		mpi.Alltoallv(send_localFace_IN2IN_iR, recv_localFace_IN2IN_iR);
		mpi.Alltoallv(send_localFace_IN2IN_ipoints, recv_localFace_IN2IN_ipoints);
		
		mpi.Alltoallv(send_localFace_IN2PR_iL, recv_localFace_IN2PR_iL);
		mpi.Alltoallv(send_globalFace_IN2PR_toProc, recv_globalFace_IN2PR_toProc);
		mpi.Alltoallv(send_localFace_IN2PR_ipoints, recv_localFace_IN2PR_ipoints);
		
		mpi.Alltoallv(send_localFace_BC2BC_iL, recv_localFace_BC2BC_iL);
		mpi.Alltoallv(send_localFace_BC2BC_BCType, recv_localFace_BC2BC_BCType);
		mpi.Alltoallv(send_localFace_BC2BC_ipoints, recv_localFace_BC2BC_ipoints);
		
		mpi.Alltoallv(send_localFace_PR2IN_iL, recv_localFace_PR2IN_iL);
		mpi.Alltoallv(send_localFace_PR2IN_toProc, recv_localFace_PR2IN_toProc);
		mpi.Alltoallv(send_localFace_PR2IN_ipoints, recv_localFace_PR2IN_ipoints);
		
		mpi.Alltoallv(send_localFace_PR2PR_iL, recv_localFace_PR2PR_iL);
		mpi.Alltoallv(send_globalFace_PR2PR_toProc, recv_globalFace_PR2PR_toProc);
		mpi.Alltoallv(send_localFace_PR2PR_ipoints, recv_localFace_PR2PR_ipoints);
		
		
		
		
		
		
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
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_localFace_BC2BC_BCType[ip].size(); ++i){
				int ibc = recv_localFace_BC2BC_BCType[ip][i];
				reorder_procFace_BC2BC_proc_id[ibc].push_back(make_pair(ip,i));
			}
		}
			
		vector<int> iter_BC2BC(size,0);
		for(int ibc=0; ibc<nbc; ++ibc){
			int str_size = meshComb.faces.size();
			for(auto& [ip, i] : reorder_procFace_BC2BC_proc_id[ibc]){
				
				face_state.push_back(BC2BC);
				
				meshComb.addFace();
				meshComb.faces.back().setType(MASCH_Face_Types::BOUNDARY);
				{
					int id_local = recv_localFace_BC2BC_iL[ip].at(i);
					int id_global = nCells_local[ip] + id_local;
					meshComb.faces.back().iL = id_global;
				}
				int tmp_size = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
				for(int j=0; j<tmp_size; ++j){
					int ipoint = recv_localFace_BC2BC_ipoints[ip].at(iter_BC2BC[ip]++);
					int id_global = points_id_local2global[ip].at(ipoint);
					meshComb.faces.back().ipoints.push_back(id_global);
				}
			}
			nFaces_boundary[ibc] = meshComb.faces.size() - str_size;
		}
	
		
		// 프로세서
		vector<int> nFaces_processor(size,0);
		{
			// IN2PR
			vector<vector<pair<int,int>>> reorder_procFace_IN2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_IN2PR_strId(size);
			for(int ip=0; ip<size; ++ip){
				// cout << recv_localFace_IN2PR_iL[ip].size() << " " <<  recv_globalFace_IN2PR_toProc[ip].size() << endl;
				for(int i=0, iter=0; i<recv_localFace_IN2PR_iL[ip].size(); ++i){
					int rightProc_global = recv_globalFace_IN2PR_toProc[ip][i];
					reorder_procFace_IN2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
					reorder_procFace_IN2PR_strId[rightProc_global].push_back(iter);
					int tmp_size = recv_localFace_IN2PR_ipoints[ip].at(iter++);
					for(int j=0; j<tmp_size; ++j) iter++;
				}
			}
			
			// PR2PR
			vector<vector<pair<int,int>>> reorder_procFace_PR2PR_proc_id(size);
			vector<vector<int>> reorder_procFace_PR2PR_strId(size);
			for(int ip=0; ip<size; ++ip){
				for(int i=0, iter=0; i<recv_localFace_PR2PR_iL[ip].size(); ++i){
					int rightProc_global = recv_globalFace_PR2PR_toProc[ip][i];
					reorder_procFace_PR2PR_proc_id[rightProc_global].push_back(make_pair(ip,i));
					reorder_procFace_PR2PR_strId[rightProc_global].push_back(iter);
					int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(iter++);
					for(int j=0; j<tmp_size; ++j) iter++;
				}
			}
			
		
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
					for(auto& [ip, i] : reorder_procFace_PR2PR_proc_id[proc]){
				
				face_state.push_back(PR2PR);
				
						meshComb.addFace();
						meshComb.faces.back().setType(MASCH_Face_Types::PROCESSOR);
						{
							int id_local = recv_localFace_PR2PR_iL[ip].at(i);
							int id_global = nCells_local[ip] + id_local;
							meshComb.faces.back().iL = id_global;
						}
						int strId = reorder_procFace_PR2PR_strId[proc].at(iter);
						int tmp_size = recv_localFace_PR2PR_ipoints[ip].at(strId++);
						for(int j=0; j<tmp_size; ++j){
							int ipoint = recv_localFace_PR2PR_ipoints[ip].at(strId++);
							int id_global = points_id_local2global[ip].at(ipoint);
							meshComb.faces.back().ipoints.push_back(id_global);
						}
						++iter;
					}
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
			}
		}
		
	
		int maxBCNum = meshComb.boundaries.size()-1;
		meshComb.boundaries[maxBCNum].startFace = meshComb.faces.size()-meshComb.boundaries[maxBCNum].nFaces;
		for(int i=maxBCNum-1; i>=0; --i){
			meshComb.boundaries[i].startFace = meshComb.boundaries[i+1].startFace-meshComb.boundaries[i].nFaces;
		}
		 
		
		
		// connPoints (접촉하는 포인트 옆 포인트가 뭔지 저장)
		vector<vector<int>> send_rightProcNo_rightId_loc(size);
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_globalPoint_toProc_toId[ip].size()/4; ++i){
				int my_id_loc = recv_globalPoint_toProc_toId[ip][iter++]; // 자신의 local id  
				int rightProcNo_glo = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 global proc 넘버
				int rightProcNo_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록 에서 자신과 중복한 포인트의 local proc
				int rightId_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 자신과 중복한 포인트의 local id 넘버
				
				// cout << rightProcNo_glo << endl;
				send_rightProcNo_rightId_loc[rightProcNo_glo].push_back(rightProcNo_loc);
				send_rightProcNo_rightId_loc[rightProcNo_glo].push_back(rightId_loc);
				
				// int new_my_id_glo = points_id_local2global[ip].at(my_id_loc); // 새로운 포인트 id
			}
				
		}
		vector<vector<int>> recv_rightProcNo_rightId_loc;
		mpi.Alltoallv(send_rightProcNo_rightId_loc, recv_rightProcNo_rightId_loc);
		vector<vector<int>> send_new_rightId_loc(size);
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0; i<recv_rightProcNo_rightId_loc[ip].size()/2; ++i){
				int rightProcNo_loc = recv_rightProcNo_rightId_loc[ip][iter++];
				int rightId_loc = recv_rightProcNo_rightId_loc[ip][iter++];
				
				int new_my_id_glo = points_id_local2global[rightProcNo_loc].at(rightId_loc); // 새로운 포인트 id
				
				send_new_rightId_loc[ip].push_back(new_my_id_glo);
			}
		}
		vector<vector<int>> recv_new_rightId_loc;
		mpi.Alltoallv(send_new_rightId_loc, recv_new_rightId_loc);
		for(int ip=0; ip<size; ++ip){
			for(int i=0, iter=0, iter2=0; i<recv_globalPoint_toProc_toId[ip].size()/4; ++i){
				int my_id_loc = recv_globalPoint_toProc_toId[ip][iter++]; // 자신의 local id  
				int rightProcNo_glo = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 global proc 넘버
				int rightProcNo_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록 에서 자신과 중복한 포인트의 local proc
				int rightId_loc = recv_globalPoint_toProc_toId[ip][iter++]; // right 블록의 자신과 중복한 포인트의 local id 넘버
				
				int new_rightId_glo = recv_globalPoint_toProc_toId[ip].at(iter2++);
				
				int new_my_id_glo = points_id_local2global[ip].at(my_id_loc); // 새로운 포인트 id
				
				auto& connPoints = meshComb.points[new_my_id_glo].connPoints;
				vector<int> tmp_first;
				for(auto& [proc, id] : connPoints){
					tmp_first.push_back(proc);
				}
				
				if(find(
				tmp_first.begin(),tmp_first.end(),rightProcNo_glo)==
				tmp_first.end()){
					connPoints.push_back(make_pair(rightProcNo_glo,new_rightId_glo));
				}
			}
		}
		
		
		
		
		
		
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	meshComb.check();
	meshComb.setFaceTypes();
	meshComb.buildCells();
	meshComb.connectFacetoPointsCells();
	meshComb.connectCelltoFaces();
	meshComb.connectCelltoPoints();
	meshComb.setCountsProcFaces();
	meshComb.setDisplsProcFaces();
	
	
	// int tmptmtp=0;
	// for(auto& face : meshComb.faces){
		// double tmp0 = pow(meshComb.points[face.ipoints[0]].x-meshComb.points[face.ipoints[2]].x,2.0);
		// tmp0 += pow(meshComb.points[face.ipoints[0]].y-meshComb.points[face.ipoints[2]].y,2.0);
		// tmp0 += pow(meshComb.points[face.ipoints[0]].z-meshComb.points[face.ipoints[2]].z,2.0);
		
		// double tmp1 = pow(meshComb.points[face.ipoints[1]].x-meshComb.points[face.ipoints[3]].x,2.0);
		// tmp1 += pow(meshComb.points[face.ipoints[1]].y-meshComb.points[face.ipoints[3]].y,2.0);
		// tmp1 += pow(meshComb.points[face.ipoints[1]].z-meshComb.points[face.ipoints[3]].z,2.0);
		
		// if(abs(tmp1-tmp0)>1.e-6){
			// cout << face_state[tmptmtp] << endl;
			// cout << "(" << meshComb.points[face.ipoints[0]].x << " " <<
			// meshComb.points[face.ipoints[0]].y << " " <<
			// meshComb.points[face.ipoints[0]].z << ") ";
			// cout << "(" << meshComb.points[face.ipoints[1]].x << " " <<
			// meshComb.points[face.ipoints[1]].y << " " <<
			// meshComb.points[face.ipoints[1]].z << ") ";
			// cout << "(" << meshComb.points[face.ipoints[2]].x << " " <<
			// meshComb.points[face.ipoints[2]].y << " " <<
			// meshComb.points[face.ipoints[2]].z << ") ";
			// cout << "(" << meshComb.points[face.ipoints[3]].x << " " <<
			// meshComb.points[face.ipoints[3]].y << " " <<
			// meshComb.points[face.ipoints[3]].z << ") " << endl;
		// }
		// if(face.ipoints.size()!=4){
			// cout << "EEEE" << endl;
		// }
		
		// ++tmptmtp;
	// }
	
	// for(auto& cell : meshComb.cells){
		// vector<int> tmp_ipoints;
	
		// for(auto& iface : cell.ifaces){
			// auto& face = meshComb.faces[iface];
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
			// for(auto& iface : cell.ifaces){
				// auto& face = meshComb.faces[iface];
				// cout << face_state[iface] << endl;
				// cout << "(" << meshComb.points[face.ipoints[0]].x << " " <<
				// meshComb.points[face.ipoints[0]].y << " " <<
				// meshComb.points[face.ipoints[0]].z << ") ";
				// cout << "(" << meshComb.points[face.ipoints[1]].x << " " <<
				// meshComb.points[face.ipoints[1]].y << " " <<
				// meshComb.points[face.ipoints[1]].z << ") ";
				// cout << "(" << meshComb.points[face.ipoints[2]].x << " " <<
				// meshComb.points[face.ipoints[2]].y << " " <<
				// meshComb.points[face.ipoints[2]].z << ") ";
				// cout << "(" << meshComb.points[face.ipoints[3]].x << " " <<
				// meshComb.points[face.ipoints[3]].y << " " <<
				// meshComb.points[face.ipoints[3]].z << ") " << endl;
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
	
	
	
	meshComb.informations();
	
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
	
	
	
	
	
	
	
	
}



