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
    
    
    
    
    
    
    
    
    
	
	// 체적분율 압축 방정식 솔버, OpenFoam compression term 참고
	{
		
		
		int gradIterMax_LG = 1;
		int gradIterMax_GG = 1;
		

		// 로컬 시간스템 구하기
        int id_vol = controls.getId_cellVar("volume");
        int id_area = controls.getId_faceVar("area");
		double tau_ls = 1.e15;
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			double maxA=0.0;
			double minA=1000000.0;
			for(auto& j : cell.faces){
				maxA = max(maxA, var.faces[j][id_area]);
				minA = min(minA, var.faces[j][id_area]);
			}
			double minX = cell.volume / maxA;
			tau_ls = min(tau_ls,min(pow(var.cells[i][id_vol],0.3), minX));
		}
		double tau_ls_glob;
		MPI_Allreduce(&tau_ls, &tau_ls_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		tau_ls = tau_ls_glob;
		
	
		
		
		for(int outerIter=0; outerIter<iterMax; ++outerIter){
			
			vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
			for(int iter=0; iter<gradIterMax_LG; ++iter){
				vector<double> dummyVec;
				math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
					controls.VF[0], controls.fVF[0], dummyVec, gradAi);
			}
			// for(int iter=0; iter<gradIterMax_GG; ++iter){
				// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
			// }
			
			
			vector<double> UN(mesh.cells.size());
			vector<double> VN(mesh.cells.size());
			vector<double> WN(mesh.cells.size());
			vector<double> alphaN(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				
				vector<double> sufNorVec(3,0.0);
				double magGrad = sqrt(
					pow(gradAi[i][0],2.0)+
					pow(gradAi[i][1],2.0)+
					pow(gradAi[i][2],2.0));
				if(magGrad != 0.0){
					for(int ii=0; ii<3; ++ii){
						sufNorVec[ii] = gradAi[i][ii]/magGrad;
					}
				}
				
				// double magVel = sqrt(
					// pow(cell.var[controls.U],2.0)+
					// pow(cell.var[controls.V],2.0)+
					// pow(cell.var[controls.W],2.0));
				double magVel = 1.0;
				UN[i] = magVel*sufNorVec[0];
				VN[i] = magVel*sufNorVec[1];
				WN[i] = magVel*sufNorVec[2];
				
				double alpha = cell.var[controls.VF[0]];
				alphaN[i] = alpha*(1.0-alpha);
				
			}
			// processor faces
			vector<double> UN_recv, VN_recv, WN_recv, alphaN_recv;
			if(size>1){
				vector<double> UN_send, VN_send, WN_send, alphaN_send;
				for(int i=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						UN_send.push_back(UN[face.owner]);
						VN_send.push_back(VN[face.owner]);
						WN_send.push_back(WN[face.owner]);
						alphaN_send.push_back(alphaN[face.owner]);
					}
				}
				mpi.setProcsFaceDatas(
							UN_send, UN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							VN_send, VN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							WN_send, WN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							alphaN_send, alphaN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
			}
			
			
			
			vector<double> resi(mesh.cells.size(),0.0);
			for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				vector<double> nvec(3,0.0);
				nvec[0] = face.unitNormals[0];
				nvec[1] = face.unitNormals[1];
				nvec[2] = face.unitNormals[2];
			
				double wCL = face.wC;
				double wCR = 1.0-wCL;
				
				double UnF = wCL*(UN[face.owner]*nvec[0]+VN[face.owner]*nvec[1]+WN[face.owner]*nvec[2]);
				double alphaL = alphaN[face.owner];
				double alphaR = 0.0;
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					UnF += wCR*(UN[face.neighbour]*nvec[0]+VN[face.neighbour]*nvec[1]+WN[face.neighbour]*nvec[2]);
					alphaR = alphaN[face.neighbour];
				}
				else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					UnF += wCR*(UN_recv[ip]*nvec[0]+VN_recv[ip]*nvec[1]+WN_recv[ip]*nvec[2]);
					alphaR = alphaN_recv[ip];
				}
				else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
					UnF = 0.0;
				}
				
				double weightL = (UnF > 0.0) ? 1.0 : 0.0;
				double weightR = 1.0 - weightL;
			
				double alphaF = weightL*alphaL + weightR*alphaR;
				
				if(alphaL*(1.0-alphaL) < 1.e-15) alphaF = 0.0;
				if(alphaR*(1.0-alphaR) < 1.e-15) alphaF = 0.0;
				
				double flux = alphaF*UnF*face.area;
				
				resi[face.owner] -= flux;
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					resi[face.neighbour] += flux;
				}
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
				
			}
			
			double residual =0.0;
			for(int i=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				
				double updateVar = CFL * tau_ls * resi[i] / cell.volume;
				
				cell.var[controls.VF[0]] += updateVar;
				
				residual += (updateVar*updateVar);
				
				cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
				
			}
			residual = sqrt(residual);
			double residualReduced;
			MPI_Allreduce(&residual, &residualReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			controls.residual = residualReduced;
			
			if(rank==0) cout << 
			" | iteration = " << outerIter << 
			" | residual = " << residualReduced << endl;
		
		}
		
		
	}
		
		
	// 저장
	{
		string foldername;
		std::ostringstream streamObj;
		streamObj << controls.time;
		foldername = "./save/" + streamObj.str() + "_CompressVF/";
		
		// solvers.setCompValuesLeftRightFace(mesh, controls, species);
		
		save.vtu(foldername, mesh, controls, species);
		
		save.cellDataOverTime(foldername, controls);
		
		// save.particles(foldername, mesh, controls, species);
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
void MASCH_Solver::setDPMFunctionsUDF(MASCH_Mesh& a, MASCH_Control& b){};
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


