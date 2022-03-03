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
	MASCH_Control controls;
	MASCH_Load load;
	MASCH_Mesh mesh;
	MASCH_Solver solver;
	MASCH_Variables var;
	MASCH_Poly_AMR_Builder amr;
	
	// 기본 variables 셋팅
	controls.setVariablesBasic();
	// 셋팅 파일 로드
	load.settingFiles("./setting/", controls);
	
	// 초기조건 대입하고 저장 후 종료
	if(mapArgv.find("-init") != mapArgv.end()){
		controls.saveAfterInitial(mesh);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	
	// 메쉬 파일 로드
	load.meshFiles(controls.getLoadFolderName(), controls, mesh);
	// B.C. 셋팅
	controls.setBoundaryConditions(mesh);
	// primitive 값 limit 셋팅
	controls.setMinMaxPrim();
	// variable들 어레이 생성
	controls.setVariableArray(mesh, var);
	// cout << mesh.cells.size() << endl;
	
	// sparse matrix의 CSR 포맷 셋팅
	var.setSparCSR(mesh, controls);
	// 메쉬 지오메트릭 셋팅
	controls.setGeometric(mesh, var);
	// B.C. 펑션 셋팅
	solver.setBoundaryFunctions(mesh, controls, var);
	// 솔버 펑션 셋팅
	solver.setFunctions(mesh, controls);
	// 그레디언트 계산시 필요한 값 셋팅
	solver.calcGradient.init(mesh, controls, var);
	
	// primitive 값 로드
	var.fields[controls.fieldVar["time"].id] = 0.0;
	load.fvmFiles(controls.getLoadFolderName(), rank, mesh, controls, var);
	// DPM 파일 로드
	// load.dpmFiles(foldername, mesh, controls);
	
	// 선형 솔버 값들 0.0 으로 만들어주기
	var.clearLinearSystems();
	// 원시변수 제외한 나머지 셀값 업데이트
	solver.updateCellAddiValues(mesh, controls, var);
	// 올드 값 초기화
	solver.initOldValues(mesh, controls, var);
	// proc 셀의 원시변수 mpi 넘기기
	solver.updateProcRightCellPrimValues(mesh, controls, var);
	// 셀 원시변수 제외한 나머지 proc 셀값 업데이트
	solver.updateProcRightCellAddiValues(mesh, controls, var);
	// 셀 그레디언트
	solver.gradientTerms(mesh, controls, var);
	// 고차 reconstruction
	solver.highOrderTerms(mesh, controls, var);
	// B.C. 원시변수 업데이트
	solver.updateBoundaryFacePrimValues(mesh, controls, var);
	// B.C. 원시변수 제외한 나머지 값 업데이트
	solver.updateBoundaryFaceAddiValues(mesh, controls, var);
	// 타임스텝 구하기
	solver.calcTempSteps(mesh, controls, var);
	
	
			// cout <<  
			// var.fields[controls.fieldVar["time-step"].id] << endl;
	int iter=0;
	while( var.fields[controls.fieldVar["time"].id] < 4.0 )
	{
		// var.fields[controls.fieldVar["time-step"].id] *= 0.001;
		
		amr.polyAMR_inline(mesh, controls, solver, var, iter);
		
		solver.fvm_inline(mesh, controls, var);
		// solver.fvm(mesh, controls, var);
		// solver.dpm(mesh, control, var.cells, var.particles);
		
		// var.updateRealTime(control);
		
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		if( (iter+1) % 10 == 0 ){
			if(rank==0) cout << iter+1 << " " << var.fields[controls.fieldVar["time"].id] << 
			" " << var.fields[controls.fieldVar["time-step"].id] << 
			" " << var.fields[controls.fieldVar["residual"].id] << 
			endl;
			controls.log.show();
		}
		
		if(isnan(var.fields[controls.fieldVar["residual"].id]) ||
		var.fields[controls.fieldVar["residual"].id] < -1.e12 ||
		var.fields[controls.fieldVar["residual"].id] > 1.e12) break;
	
		// if( controls.checkSaveFiles() )
		if( (iter+1) % 100 == 0 )
		{
			string foldername;
			std::ostringstream streamObj;
			streamObj << var.fields[controls.fieldVar["time"].id];
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			save.fvmFiles(foldername, rank, mesh, controls, var);
			// save.dpmFiles();
			// save.udfFiles();
		}
		++iter;
	// controls.log.show();
	}
	// save.fvmFiles("./save/test/", rank, mesh, controls, var);
	
	
	
	
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





		// cout << var.fields[controls.fieldVar["time-step"].id] << endl;
		
		// {
					
			// int id_dt = controls.fieldVar["time-step"].id;
			// int id_vol = controls.cellVar["volume"].id;
			// int id_x = controls.cellVar["x-velocity"].id;
			// int id_y = controls.cellVar["y-velocity"].id;
			// int id_z = controls.cellVar["z-velocity"].id;
			// int id_c = controls.cellVar["speed of sound"].id;
			// var.fields[id_dt] = 1.e12;
			
			// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				// double inp_val = pow(var.cells[i][id_vol],0.3)/(
				// sqrt(
				// var.cells[i][id_x]*var.cells[i][id_x]+
				// var.cells[i][id_y]*var.cells[i][id_y]+
				// var.cells[i][id_z]*var.cells[i][id_z])+
				// var.cells[i][id_c]);
				// var.fields[id_dt] = min(var.fields[id_dt],inp_val*0.0001);
			// }
			
			// if(size>1){
				// double tmp_fieldVar = var.fields[id_dt];
				// double tmp_fieldVar_glo;
				// MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
				// var.fields[id_dt] = tmp_fieldVar_glo;
			// }
			
		// }
		