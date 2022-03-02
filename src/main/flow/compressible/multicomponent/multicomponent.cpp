#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <iomanip>
#include "parmetis.h" 
#include "scotch.h" 

#include "../../../../others/mesh.h"  
#include "../../../../others/mpi.h"
#include "../../../../others/load.h" 
#include "../../../../others/variables.h"
#include "../../../../others/solvers.h"
#include "../../../../others/controls.h"
#include "../../../../others/save.h"
#include "../../../../others/polyAMR.h"

void print_help();

int main(int argc, char* argv[]) {

	// MPI 초기화
	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// 옵션 읽어들이기
	map<string,string> mapArgv;
	for(int i=1; i<argc; i+=2){
		string first = argv[i]; string second;
		if(i+1==argc){ second = "nan"; }
		else{ second = argv[i+1]; }
		mapArgv.insert(make_pair(first, second));
	}
	
	if(mapArgv.find("-help") != mapArgv.end() ||
	   mapArgv.find("-h") != mapArgv.end() ){
		if(rank==0) print_help();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	
	// 필요한 객체들 생성
	MASCH_Mesh_Save save;
	MASCH_Control control;
	MASCH_Load load;
	MASCH_Mesh mesh;
	MASCH_Variables var;
	MASCH_Solver solver;
	
	// 기본 variables 셋팅
	control.setVariablesBasic();
	// 셋팅 파일 로드
	load.settingFiles("./setting/", control);
	// 메쉬 파일 로드
	load.meshFiles("./grid/0/", control, mesh);
	// B.C. 셋팅
	control.setBoundaryConditions(mesh);
	// primitive 값 limit 셋팅
	control.setMinMaxPrim();
	// variable들 어레이 생성
	control.setVariableArray(mesh, var);
	
	// 초기조건 대입하고 저장 후 종료
	if(mapArgv.find("-init") != mapArgv.end() ){
		save.fvmFiles("./save/0/", rank, mesh, control, var);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	// sparse matrix의 CSR 포맷 셋팅
	var.setSparCSR(mesh, control);
	// 메쉬 지오메트릭 셋팅
	control.setGeometric(mesh, var);
	// B.C. 펑션 셋팅
	solver.setBoundaryFunctions(mesh, control, var);
	// 솔버 펑션 셋팅
	solver.setFunctions(mesh, control);
	// 그레디언트 계산시 필요한 값 셋팅
	solver.calcGradient.init(mesh, control, var);
	
	MASCH_Poly_AMR_Builder amr;
	amr.polyAMR(mesh, control, var, 0);
	save.fvmFiles("./save/1/", rank, mesh, control, var);
	
	// primitive 값 로드
	// save.fvmFiles(control.getFolderName(), rank, mesh, control, var);
	// // load.dpmFiles(foldername, mesh, controls);
	
	
	
	
	
	MPI::Finalize();
	return EXIT_SUCCESS;
	
	
	// control.log.push("cell values");
	// {
		// int id = control.cellVar["pressure"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 101325.0;
		// }
	// }
	// {
		// int id = control.cellVar["x-velocity"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 0.0;
		// }
	// }
	// {
		// int id = control.cellVar["y-velocity"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 0.0;
		// }
	// }
	// {
		// int id = control.cellVar["z-velocity"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 0.0;
		// }
	// }
	// {
		// int id = control.cellVar["temperature"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 300.0;
		// }
	// }
	// {
		// int id = control.cellVar["water"].id;
		// auto cells = mesh.cells.data();
		// auto cellVar = var.cells.data();
		// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// auto& cell = cells[i];
			// auto cellVar_i = cellVar[i].data();
			// cellVar_i[id] = 0.0;
		// }
	// }
	// control.log.pop();
	
	// control.log.push("updateCellAddiValues");
	// solver.updateCellAddiValues(mesh, control, var);
	// control.log.pop();
	// control.log.push("updateProcRightCellPrimValues");
	// solver.updateProcRightCellPrimValues(mesh, control, var);
	// control.log.pop();
	// control.log.push("updateProcRightCellAddiValues");
	// solver.updateProcRightCellAddiValues(mesh, control, var);
	// control.log.pop();
	// control.log.push("updateBoundaryFacePrimValues");
	// solver.updateBoundaryFacePrimValues(mesh, control, var);
	// control.log.pop();
	// control.log.push("updateBoundaryFaceAddiValues");
	// solver.updateBoundaryFaceAddiValues(mesh, control, var);
	// control.log.pop();
	// control.log.push("initOldValues");
	// solver.initOldValues(mesh, control, var);
	// control.log.pop();
	
	
	// while( control.checkContinueRealTimeStep(var) ){
	// while( var.fields[control.fieldVar["time"].id]<control.endTime )
	{
		
		// solver.fvm(mesh, control, var);
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





void print_help(){

	cout << endl;
	cout << "┌─────── Comp. Multicomponent helper ───────────────────────────── " << endl;
	cout << "| -init : initialization" << endl;
	cout << "└───────────────────────────────────────────────────────────────── " << endl;
	cout << endl;
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




