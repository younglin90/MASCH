#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>

#include "parmetis.h" 
#include "scotch.h" 

#include "../../mesh/mesh.h" 
#include "../../load/load.h" 
#include "../../mesh/geometric.h" 
#include "../../mesh/polyAMR.h"

#include "../../controls/controls.h" 
#include "../../solvers/solvers.h" 


int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	vector<SEMO_Species> species;
	
	SEMO_Controls_Builder controls;
		
	controls.readSpecies(species);
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	SEMO_Solvers_Builder solvers;
	
	SEMO_Mesh_Builder mesh;
	
	SEMO_MPI_Builder mpi;
	
	
	bool boolLoad = true;
	bool boolPartitioning = false;
	bool boolAMR = true;
	bool boolGeometric = true;
	
	if(boolLoad){
		
		SEMO_Mesh_Load load;
		double starttime = stod(controls.startFrom);
		string foldername;
		std::ostringstream streamObj;
		streamObj << starttime;
		foldername = "./save/" + streamObj.str() + "/";
		
		if(starttime == 0.0){
			foldername = "./save/0/";
		}
		
		load.vtu(foldername, mesh, controls, species);
		
		// solvers.calcCellEOSVF(mesh, controls, species);
		
		solvers.calcCellEOSMF(mesh, controls, species);
		
		for(auto& cell : mesh.cells){
			cell.var[controls.oldP] = cell.var[controls.P];
			cell.var[controls.oldU] = cell.var[controls.U];
			cell.var[controls.oldV] = cell.var[controls.V];
			cell.var[controls.oldW] = cell.var[controls.W];
			cell.var[controls.oldT] = cell.var[controls.T];
			cell.var[controls.oldRho] = cell.var[controls.Rho];
			cell.var[controls.oldHt] = cell.var[controls.Ht];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.oldMF[i]] = cell.var[controls.MF[i]];
			}
		}
			
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				face.var[controls.Un] = 0.5 * mesh.cells[face.owner].var[controls.U] * face.unitNormals[0];
				face.var[controls.Un] += 0.5 * mesh.cells[face.owner].var[controls.V] * face.unitNormals[1];
				face.var[controls.Un] += 0.5 * mesh.cells[face.owner].var[controls.W] * face.unitNormals[2];
				face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.U] * face.unitNormals[0];
				face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.V] * face.unitNormals[1];
				face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.W] * face.unitNormals[2];
				
				face.var[controls.oldUn] = face.var[controls.Un];
			}
		}
	
	}
	

	
	
	// geometric
	if(boolGeometric){
		SEMO_Utility_Math math;
		
		SEMO_Mesh_Geometric geometric;
		
		geometric.init(mesh);
		
		math.initLeastSquare(mesh);
	}

	

	bool initFlow = false;
	
	// flow initialization
	if(initFlow){
		
		solvers.setInitValues(mesh, controls);
		
		mesh.saveFile("vtu", "./save/0/", controls);
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	

	bool calcFlow = true;
	
	// flow calculation
	if(calcFlow){
		
		solvers.calcCellTransport(mesh, controls, species);
		
		
		SEMO_Mesh_Save save;
		
		while(
		controls.iterReal<controls.iterRealMax ||
		controls.time<stod(controls.stopAt) ){
			
	
			if(rank==0) {
				cout << "| real-time step = " << controls.iterReal 
				<< " | time = " << controls.time;
			}
		
			solvers.compressibleCoupled(mesh, controls, species);

			//==============================
			// AMR
			if(
			// controls.iterReal != 0 && 
			( (controls.iterReal+1) % controls.intervalRefine == 0 ||
			  (controls.iterReal+1) % controls.intervalUnrefine == 0)
			){ 
				SEMO_Poly_AMR_Builder AMR;
				AMR.polyAMR(mesh, controls, species, 0);
				mesh.cellsGlobal();
				solvers.calcCellEOSMF(mesh, controls, species);
				solvers.calcCellTransport(mesh, controls, species);
			} 
			//==============================
				
			controls.time += controls.timeStep;
			
			++controls.iterReal;
			
			if(controls.saveControl == "timeStep"){
				if(controls.iterReal % (int)controls.saveInterval == 0){
					string foldername;
					std::ostringstream streamObj;
					streamObj << controls.time;
					foldername = "./save/" + streamObj.str() + "/";
					
					solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					save.vtu(foldername, mesh, controls, species);
					
					// save.particles(foldername, mesh, controls, species);
					
				}
			}
			else if(controls.saveControl == "runTime"){
				int jung = controls.time / controls.saveInterval;
				double namuji = controls.time - (double)jung * controls.saveInterval;
				if(
				namuji < controls.timeStep &&
				namuji >= 0.0
				){
					string foldername;
					std::ostringstream streamObj;
					streamObj << controls.time;
					foldername = "./save/" + streamObj.str() + "/";
					
					solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					save.vtu(foldername, mesh, controls, species);
					
					// save.particles(foldername, mesh, controls, species);
					
				}
			}
		}
	}
	
	
	if(rank==0) cout << "| End Program" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	return EXIT_SUCCESS;
}