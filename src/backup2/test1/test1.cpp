#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <iomanip>
#include "parmetis.h" 
#include "scotch.h" 

#include "./mesh.h"  
#include "./mpi.h"
#include "./load.h" 
#include "./variables.h"
#include "./solvers.h"
#include "./controls.h"
#include "./save.h"



int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// MASCH_Log log;
	MASCH_Mesh_Save save;
	
	MASCH_Control control;
	control.setVariablesBasic();
	
	MASCH_Load load;
	// control = load.settingFiles("./setting/");
	load.settingFiles("./setting/", control);
	
	MASCH_Mesh mesh;
	// load.meshFiles(control.getFolderName(), control, mesh);
	
	// 초기조건 대입하고 저장 or 로딩하고 계쏙 계산
	load.meshFiles("./grid/0/", control, mesh);
	
	control.setBoundaryConditions(mesh);
	
	MASCH_Variables var;
	control.setMinMaxPrim();
	control.setVariableArray(mesh, var);
	
	var.setSparCSR(mesh, control);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	// var.setInitialization(mesh, control);
	control.setGeometric(mesh, var);
	// // load.fvmFiles(foldername, mesh, controls);
	// // load.dpmFiles(foldername, mesh, controls);
	
	MASCH_Solver solver;
	solver.setBoundaryFunctions(mesh, control, var);
	solver.setFunctions(mesh, control);
	solver.calcGradient.init(mesh, control, var);
	
	
	control.log.push("cell values");
	{
		int id = control.cellVar["pressure"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 101325.0;
		}
	}
	{
		int id = control.cellVar["x-velocity"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 0.0;
		}
	}
	{
		int id = control.cellVar["y-velocity"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 0.0;
		}
	}
	{
		int id = control.cellVar["z-velocity"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 0.0;
		}
	}
	{
		int id = control.cellVar["temperature"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 300.0;
		}
	}
	{
		int id = control.cellVar["water"].id;
		auto cells = mesh.cells.data();
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			cellVar_i[id] = 0.0;
		}
	}
	control.log.pop();
	
	control.log.push("updateCellAddiValues");
	solver.updateCellAddiValues(mesh, control, var);
	control.log.pop();
	control.log.push("updateProcRightCellPrimValues");
	solver.updateProcRightCellPrimValues(mesh, control, var);
	control.log.pop();
	control.log.push("updateProcRightCellAddiValues");
	solver.updateProcRightCellAddiValues(mesh, control, var);
	control.log.pop();
	control.log.push("updateBoundaryFacePrimValues");
	solver.updateBoundaryFacePrimValues(mesh, control, var);
	control.log.pop();
	control.log.push("updateBoundaryFaceAddiValues");
	solver.updateBoundaryFaceAddiValues(mesh, control, var);
	control.log.pop();
	control.log.push("initOldValues");
	solver.initOldValues(mesh, control, var);
	control.log.pop();
	
	
	
	// while( control.checkContinueRealTimeStep(var) ){
	// while( var.fields[control.fieldVar["time"].id]<control.endTime )
	{
		
		solver.fvm(mesh, control, var);
		// solver.dpm(mesh, control, var.cells, var.particles);
		
		// var.updateRealTime(control);
		
		var.fields[control.fieldVar["time"].id] += 
		var.fields[control.fieldVar["time-step"].id];
		
		// mesh.amr();
	
		// if( control.checkSaveFiles() ){
			// save.fvmFiles("./save/1/", rank, mesh, control, var);
			// // save.dpmFiles();
			// // save.udfFiles();
		// }
	}
	
	control.log.show();
	
	save.fvmFiles("./save/1/", rank, mesh, control, var);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	
	
	MPI::Finalize();
	return EXIT_SUCCESS;
}








void MASCH_Control::setMinMaxPrim(){
	
	limitMaxPrim.resize(primVarNames.size(),1.0);
	limitMinPrim.resize(primVarNames.size(),0.0);
	
	limitMaxPrim[0] = 1.e12;
	limitMinPrim[0] = 10.0;
	
	limitMaxPrim[1] = 1.e12;
	limitMinPrim[1] = -1.e12;
	
	limitMaxPrim[2] = 1.e12;
	limitMinPrim[2] = -1.e12;
	
	limitMaxPrim[3] = 1.e12;
	limitMinPrim[3] = -1.e12;
	
	limitMaxPrim[4] = 5000.0;
	limitMinPrim[4] = 10.0;
	
}




