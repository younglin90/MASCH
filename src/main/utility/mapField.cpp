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
	
	
	for(int i=0; i<1020; ++i)
	{
		MASCH_Mesh tmp_mesh;
		MASCH_Variables tmp_var;
		
		load.meshFiles(controls.getLoadFolderName(), controls, tmp_mesh);
		controls.setBoundaryConditions(tmp_mesh);
		controls.setMinMaxPrim();
		controls.setVariableArray(tmp_mesh, tmp_var);
		tmp_var.setSparCSR(tmp_mesh, controls);
		controls.setGeometric(tmp_mesh, tmp_var);
		solver.setBoundaryFunctions(tmp_mesh, controls, tmp_var);
		solver.setFunctions(tmp_mesh, controls);
		solver.calcGradient.init(tmp_mesh, controls, tmp_var);
		
		load.fvmFiles(controls.getLoadFolderName(), rank, tmp_mesh, controls, tmp_var);
		
		
		
	}

	save.fvmFiles("./save/mapField/", rank, mesh, controls, var);
	
	MPI::Finalize();
	return EXIT_SUCCESS;
	
}



void mapField_vtuPrimitive(
	string folderName, int rank,
	MASCH_Control& controls, 
	MASCH_Mesh& mesh,
	MASCH_Variables& var){
		
		
	string saveFolderName = folderName;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	string nextToken;
	// boolBinary = false;
	boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				boolCompress = true;
			}
			break;
		}
	}
	
	vector<string> tmp_primVarNames;
	tmp_primVarNames.push_back("pressure");
	tmp_primVarNames.push_back("temperature");
	tmp_primVarNames.push_back("mass-fraction-water");
	tmp_primVarNames.push_back("mass-fraction-water");
	
	
	for(auto& name : controls.primVarNames){
		vector<double> tmp_cellVars;
		loadDatasAtVTU(inputFile, name, tmp_cellVars);
		int id = controls.cellVar[name].id;
		int iter = 0;
		for(auto& cell : var.cells){
			cell[id] = tmp_cellVars[iter];
			++iter;
		}
	}
	
	
	inputFile.close();
	
	
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


