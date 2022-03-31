#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <iomanip>
#include "parmetis.h" 
#include "scotch.h" 

#include "../../../others/mesh.h"  
#include "../../../others/mpi.h"
#include "../../../others/load.h" 
#include "../../../others/variables.h"
#include "../../../others/solvers.h"
#include "../../../others/controls.h"
#include "../../../others/save.h"
#include "../../../others/polyAMR.h"

void print_help();

int main(int argc, char* argv[]) {

	// MPI 초기화
	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// 필요한 객체들 생성
	MASCH_Mesh_Save save;
	MASCH_Control controls;
	MASCH_Load load;
	MASCH_Mesh mesh;
	MASCH_Solver solver;
	MASCH_Variables var;
	MASCH_Poly_AMR_Builder amr;
	
	// 옵션 읽어들이기
	// map<string,string> mapArgv;
	for(int i=1; i<argc; i+=2){
		string first = argv[i]; string second;
		if(i+1==argc){ second = "nan"; }
		else{ second = argv[i+1]; }
		controls.mapArgv.insert(make_pair(first, second));
	}
	
	if(controls.mapArgv.find("-help") != controls.mapArgv.end() ||
	   controls.mapArgv.find("-h") != controls.mapArgv.end() ){
		if(rank==0) print_help();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	
	
	// 기본 variables 셋팅
	controls.setVariablesBasic();
	// 셋팅 파일 로드
	load.settingFiles("./setting/", controls);
	
	// 초기조건 대입하고 저장 후 종료
	if(controls.mapArgv.find("-init") != controls.mapArgv.end()){
		load.meshFiles("./grid/0/", controls, mesh);
		controls.setGeometricOnlyCell_xyz(mesh);
		controls.saveAfterInitial(mesh);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	
	// 메쉬 파일 로드
	load.meshFiles(controls.getLoadFolderName(), controls, mesh);
	
	// variable들 어레이 생성
	controls.setVariableArray(mesh, var);
	
	// 메쉬 지오메트릭 셋팅
	controls.setGeometric(mesh, var);
	// B.C. 펑션 셋팅
	solver.setBoundaryFunctions(mesh, controls, var);
	
	// 솔버 펑션 셋팅
	solver.setFunctions(mesh, controls);
	// 그레디언트 계산시 필요한 값 셋팅
	solver.calcGradient.init(mesh, controls, var);
	
	// primitive 값 로드
	var.fields[controls.fieldVar["time"].id] = stod(controls.startFrom);
	
	load.fvmFiles(controls.getLoadFolderName(), rank, mesh, controls, var);
	
	// sparse matrix의 CSR 포맷 셋팅
	var.setSparCSR(mesh, controls);
	
	// DPM 파일 로드
	// load.dpmFiles(foldername, mesh, controls);
	
	// 초기 셋팅
	solver.updateProcRightCellValues_All(mesh, controls, var);
	solver.updateCellAddiValues_All(mesh, controls, var);
	solver.updateBoundaryFacePrimValues_All(mesh, controls, var);
	solver.gradientTerms_All(mesh, controls, var);
	solver.curvatureTerms_All(mesh, controls, var);
	solver.updateProcRightCellValues_All(mesh, controls, var);
	solver.initOldValues(mesh, controls, var);
	
	
	var.fields[controls.fieldVar["parcel-injection-accum-time"].id] = 0.0;
	controls.iterReal=0;
	var.fields[controls.fieldVar["time-step"].id] = 0.1;
	// while( controls.iterReal < 1000 )
	{
		for(int ii=0; ii<10; ++ii){
			solver.dpm(mesh, controls, var);
			
			// 시간 업데이트
			var.fields[controls.fieldVar["time"].id] +=
			var.fields[controls.fieldVar["time-step"].id];
			
			// if(controls.iterReal % 100 == 0)
			{
				double time = var.fields[controls.fieldVar["time"].id];
				save.fvmFiles("./save/"+to_string(time)+"/", rank, mesh, controls, var);
				save.parcels("./save/"+to_string(time)+"/", rank, mesh, controls, var);
				// save.fvmFiles_boundary("./save/nan/", rank, mesh, controls, var);
			}
		}
		
		++controls.iterReal;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI::Finalize();
	return EXIT_SUCCESS;
	
	
	// 메인 솔버
	controls.iterReal=0; 
	bool bool_resi_isnan = false;
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["coupled.nCorrectors"]);
	int maxIterOuter_fvm = stoi(s_iter[0]);
	// int maxIterOuter_dpm = stoi(s_iter[0]);
	while( 
	var.fields[controls.getId_fieldVar("time")] < controls.stopAt &&
	!bool_resi_isnan
	){
		// amr.polyAMR(mesh, controls, solver, var, controls.iterReal);
		
		// old 값 업데이트
		solver.updateOldValues(mesh, controls, var);
		
		
		controls.iterPseudo=0;
		for(int iterOuter=0; iterOuter<maxIterOuter_fvm; ++iterOuter){
			// 레지듀얼 초기화
			var.fields[controls.getId_fieldVar("residual")] = 0.0;
			for(int iSegEq=0, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
				
				// solver.fvm(mesh, controls, var, iSegEq);
			}
			
			controls.show_residual(var);
			if(controls.check_isnan(var.fields[controls.getId_fieldVar("residual")])) 
				bool_resi_isnan=true;
			
			++controls.iterPseudo;
		}
		
		// solver.fvm_to_dpm(mesh, controls, var);
		
		// controls.iterPseudo=0;
		// for(int iterOuter=0; iterOuter<maxIterOuter_dpm; ++iterOuter){
			
			// solver.dpm(mesh, controls, var);
			
			// controls.show_dpm_state(var);
			
			// ++controls.iterPseudo;
		// }
		
		
		// 시간 업데이트
		var.fields[controls.fieldVar["time"].id] +=
		var.fields[controls.fieldVar["time-step"].id];
		
		controls.save_fvmFiles(mesh, var);
		
		++controls.iterReal;
		
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	{
		save.fvmFiles("./save/nan/", rank, mesh, controls, var);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	// controls.show_residual(var);
	
	
	MPI::Finalize();
	return EXIT_SUCCESS;
	
}





void print_help(){

	cout << endl;
	cout << "┌─────── Comp. Potential helper ───────────────────────────── " << endl;
	cout << "| -n : # of iteration" << endl;
	cout << "| -r : relaxation factor" << endl;
	cout << "└───────────────────────────────────────────────────────────────── " << endl;
	cout << endl;
}





// void MASCH_Control::setMinMaxPrim(){
	
	// limitMaxPrim.resize(primVarNames.size(),1.0);
	// limitMinPrim.resize(primVarNames.size(),0.0);
	
	// limitMaxPrim[0] = 1.e12;
	// limitMinPrim[0] = 10.0;
	
	// limitMaxPrim[1] = 1.e12;
	// limitMinPrim[1] = -1.e12;
	
	// limitMaxPrim[2] = 1.e12;
	// limitMinPrim[2] = -1.e12;
	
	// limitMaxPrim[3] = 1.e12;
	// limitMinPrim[3] = -1.e12;
	
	// limitMaxPrim[4] = 5000.0;
	// limitMinPrim[4] = 10.0;
	
// }





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
		