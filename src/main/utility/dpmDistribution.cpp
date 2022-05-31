#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <iomanip>
using namespace std;
#include "parmetis.h" 
#include "scotch.h" 

#include "../../others/mesh.h"  
#include "../../others/mpi.h"
#include "../../others/load.h" 
#include "../../others/variables.h"
#include "../../others/solvers.h"
#include "../../others/controls.h"
#include "../../others/save.h"
#include "../../others/polyAMR.h"

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
	// MASCH_Poly_AMR_Builder amr;
	
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
	
	// // 초기조건 대입하고 저장 후 종료
	// if(controls.mapArgv.find("-init") != controls.mapArgv.end()){
		// load.meshFiles("./grid/0/", controls, mesh);
		// controls.setGeometricOnlyCell_xyz(mesh);
		// controls.saveAfterInitial(mesh);
		// MPI::Finalize();
		// return EXIT_SUCCESS;
	// }
	
	// // 메쉬 파일 로드
	// load.meshFiles(controls.getLoadFolderName(), controls, mesh);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// dpm 갯수 & 필수 정보 로드
	load.dpmSizeFiles(controls.getLoadFolderName(), controls, mesh);
    
    
	// parcels
	controls.setVarible({"parcel"},"position","xyz","","","vector3",
						{"x-position","y-position","z-position"},{"u","v","w"},
						{"","",""});
	controls.setVarible({"parcel"},"diameter","","","","scalar");
	controls.setVarible({"parcel"},"density","","","","scalar");
	controls.setVarible({"parcel"},"velocity","U","","","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{"","",""});
	controls.setVarible({"parcel"},"temperature","","","","scalar");
	controls.setVarible({"parcel"},"number-of-parcel","","","","scalar");
	controls.setVarible({"field"},"parcel-injection-accum-time","","","","scalar");
	controls.setVarible({"field"},"time-step-parcels","","","","scalar");
    
    var.parcels.resize(mesh.parcels.size());
    for(auto& parcel : var.parcels){
        parcel.resize(controls.parcelVar.size());
    }
	
	load.dpmFiles(controls.getLoadFolderName(), rank, mesh, controls, var);
    
    int id_dia = controls.getId_parcelVar("diameter");
    // int dist_N = 10;
    // double dist_min = 1.e-5;
    // double dist_max = 2.e-4;
    int dist_N = stoi(controls.mapArgv["-n"]);
    double dist_min = stod(controls.mapArgv["-min"]);
    double dist_max = stod(controls.mapArgv["-max"]);
    vector<int> dpm_distribution(dist_N,0);
    vector<double> dpm_distribution_inter(dist_N+1);
    dpm_distribution_inter[0] = dist_min;
    dpm_distribution_inter[dist_N] = dist_max;
    double inter = (dist_max-dist_min)/(double)dist_N;
    for(int i=1; i<dist_N; ++i){
        dpm_distribution_inter[i] = dpm_distribution_inter[i-1] + inter;
    }
    for(auto& parcel : var.parcels){
        double diameter = parcel[id_dia];
        for(int i=0; i<dist_N; ++i){
            if(dpm_distribution_inter[i] < diameter && diameter < dpm_distribution_inter[i+1]){
                ++dpm_distribution[i];
                break;
            }
        }
    }
    
    
	MPI_Barrier(MPI_COMM_WORLD);
    vector<int> dpm_distributionTot(dist_N);
    MPI_Allreduce(dpm_distribution.data(), dpm_distributionTot.data(),dist_N,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
    if(rank==0){
        for(int i=0; i<dist_N; ++i){
            double mean_dist = 0.5*(dpm_distribution_inter[i]+dpm_distribution_inter[i+1]);
            cout << mean_dist << ", " << dpm_distributionTot[i] << endl;
        }
    }
    
    
    // cout << mesh.parcels.size() << endl;
	
	// // variable들 어레이 생성
	// controls.setVariableArray(mesh, var);
	
	// // 메쉬 지오메트릭 셋팅
	// controls.setGeometric(mesh, var);
	// // B.C. 펑션 셋팅
	// solver.setBoundaryFunctions(mesh, controls, var);
	
	// // 솔버 펑션 셋팅
	// solver.setFunctions(mesh, controls);
	// // // 그레디언트 계산시 필요한 값 셋팅
	// // solver.calcGradient.init(mesh, controls, var);
	
	// // // primitive 값 로드
	// // var.fields[controls.fieldVar["time"].id] = stod(controls.startFrom);
	
	// // load.fvmFiles(controls.getLoadFolderName(), rank, mesh, controls, var);
	// load.dpmFiles(controls.getLoadFolderName(), rank, mesh, controls, var);
    
    
	
	// // sparse matrix의 CSR 포맷 셋팅
	// var.setSparCSR(mesh, controls);
	
	// // DPM 파일 로드
	// // load.dpmFiles(foldername, mesh, controls);
	
	// // 초기 셋팅
	// solver.updateProcRightCellValues_All(mesh, controls, var);
	// solver.updateCellAddiValues_All(mesh, controls, var);
	// solver.updateBoundaryFacePrimValues_All(mesh, controls, var);
	// solver.gradientTerms_All(mesh, controls, var);
	// solver.curvatureTerms_All(mesh, controls, var); 
	// solver.updateProcRightCellValues_All(mesh, controls, var);
	// solver.initOldValues(mesh, controls, var);
	
	
	// // AMR 후 초기조건 대입하고 저장 후 종료
	// if(controls.mapArgv.find("-initAMR") != controls.mapArgv.end()){
		// controls.iterReal = stoi(controls.dynamicMeshMap["AMR.interval"])-1;
		// int maxIterAMR = stoi(controls.mapArgv["-initAMR"]);
		// for(int i=0; i<maxIterAMR; ++i){
			// amr.polyAMR(mesh, controls, solver, var, controls.iterReal);
			// controls.saveAfterInitialAMR(mesh, var);
		// }
		// save.fvmFiles("./save/0_AMR/", rank, mesh, controls, var);
		// save.fvmFiles_boundary("./save/0_AMR/", rank, mesh, controls, var);
		// MPI::Finalize();
		// return EXIT_SUCCESS;
	// }
	
	
	
	// // 메인 솔버
	// controls.iterReal=0; 
	// bool bool_resi_isnan = false;
	// vector<string> s_iter = load.extractVector(controls.fvSolutionMap["coupled.nCorrectors"]);
	// int maxIterOuter_fvm = stoi(s_iter[0]);
	// // int maxIterOuter_dpm = stoi(s_iter[0]);
	// while( 
	// var.fields[controls.getId_fieldVar("time")] < controls.stopAt &&
	// !bool_resi_isnan
	// ){
		
		// {
			// amr.polyAMR(mesh, controls, solver, var, controls.iterReal);
			// // save.fvmFiles("./save/afterAMR/", rank, mesh, controls, var);
		// }
		
		// // old 값 업데이트
		// solver.updateOldValues(mesh, controls, var);
		
		// // time-step 계산
		// solver.calcTempSteps(mesh, controls, var, 0); //solver.calcTempSteps(mesh, controls, var, iSegEq);
        
        // // cout.precision(15);
        // // cout << var.fields[controls.fieldVar["time-step"].id] << endl;
		
		// controls.iterPseudo=0;
		// for(int iterOuter=0; iterOuter<maxIterOuter_fvm; ++iterOuter){
			// // 레지듀얼 초기화
			// var.fields[controls.getId_fieldVar("residual")] = 0.0;
            // // var.fields[controls.getId_fieldVar("residual-volume")] = 0.0;
			// for(int iSegEq=0, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
				
				// solver.fvm(mesh, controls, var, iSegEq);
			// }
			
			// controls.show_residual(var);
			// if(controls.check_isnan(var.fields[controls.getId_fieldVar("residual")])) 
				// bool_resi_isnan=true;
			
			// ++controls.iterPseudo;
		// }
			
		// // if(controls.nameParcels.size()!=0)
		// {
			// solver.dpm(mesh, controls, var);
			// controls.show_dpm_information();
		// }
        
        // // mean 값 계산
        // solver.calcMeanValues(mesh, controls, var);
        
		
		// // 시간 업데이트
		// var.fields[controls.fieldVar["time"].id] +=
		// var.fields[controls.fieldVar["time-step"].id];
		
		// controls.save_fvmFiles(mesh, var);
		// controls.save_dpmFiles(mesh, var);
		// controls.save_pvdFile(mesh, var);
		
		// ++controls.iterReal;
		
	// }
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // {
		// // save.fvmFiles("./save/nan/", rank, mesh, controls, var);
	// // }
	// MPI_Barrier(MPI_COMM_WORLD); 
    // if(rank==0) cout << endl << "| Program END |" << endl << endl;
	
	// // controls.show_residual(var);
	
	
	MPI::Finalize();
	return EXIT_SUCCESS;
	
}


void MASCH_Control::setVariablesUDF(vector<string>& a){};
// void MASCH_Solver::setFunctions(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setBoundaryFunctions(MASCH_Mesh& a, MASCH_Control& b, MASCH_Variables& c){};
void MASCH_Solver::setOldVFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setSegEqUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setTimeStepFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setGradFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setCurvatureFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setTermsCellLoopFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setTermsFaceLoopFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setUpdatePrimFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setDPMFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b, MASCH_Variables& c){};
void MASCH_Solver::setMinMaxCellValuesFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
void MASCH_Solver::setMeanCellValuesFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};

void print_help(){

	cout << endl;
	cout << "┌─────── Comp. Coupled solve helper ───────────────────────────── " << endl;
	cout << "| -n : number of distribution" << endl;
	cout << "| -min : min value" << endl;
	cout << "| -max : max value" << endl;
	cout << "└───────────────────────────────────────────────────────────────── " << endl;
	cout << endl;
}


