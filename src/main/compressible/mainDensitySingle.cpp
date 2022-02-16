#include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <sstream>
// #include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
// #include <cstddef>
// #include <locale>
// #include <numeric>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <stack>
#include <cstring>

#include "parmetis.h" 
#include "scotch.h" 

#include "../../mesh/mesh.h"  
#include "../../mesh/geometric.h" 
// #include "../../mesh/polyAMR.h"
#include "../../load/load.h" 
#include "../../controls/controls.h" 
#include "../../mpi/mpi.h"
// #include "../../solvers/solvers.h"  

// wchar_t wStr[] = L"€áa¢cée£";
// int iStr = sizeof(wStr) / sizeof(wStr[0]);        // length of the string
// wchar_t *pStr = 0;

// using namespace chrono; // std 내에 chrono 가 존재  

#define MASCH_FFL __FILE__,__FUNCTION__,__LINE__

// class MASCH_MPI {
// private:
// public:
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	// MASCH_MPI(){}
	// // ~MASCH_MPI(){ MPI::Finalize(); }
	
	// vector<string> gatherv(vector<string>& inp){
		// if(size>1){
			// vector<int> counts(size,0);
			// vector<int> disp(size+1,0);
			
			// std::vector<char> cstrings;
			// for(std::string s: inp)
			// {
				// for(int i = 0; i < strlen(s.c_str()); ++i)
				// {
					// cstrings.push_back(s.c_str()[i]);
				// }
			// }
			
			// int my_count = cstrings.size();
			// MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
			// for(int i=1; i<size+1; ++i) disp[i] = disp[i-1] + counts[i-1];
			
			// std::vector<char> cstrings_recv(disp[size]);
			// if(rank==0){
				// MPI_Gatherv(
				// cstrings.data(), cstrings.size(), MPI_CHAR, 
				// cstrings_recv.data(), counts.data(), disp.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
				
			// }
			// else{
				// MPI_Gatherv(
				// cstrings.data(), cstrings.size(), MPI_CHAR, 
				// NULL, NULL, NULL, MPI_CHAR, 0, MPI_COMM_WORLD);
			// }
			
			// vector<string> results;
			// if(rank==0){
				// for(int ip=0; ip<size; ++ip){
					// std::stringstream ss;
					// for(int i=disp[ip]; i<disp[ip+1]; ++i){
						// ss << cstrings_recv[i];
					// }
					// results.push_back(ss.str());
				// }
			// }
			// return results;
		// }
	// }

// };


class MASCH_FileSystem {
private:

public:

	// 파일들 pretty하게 표시
	string showPrettySize(int fileSize){
		string logging = "";
		if(fileSize/1024 < 1){
			logging = to_string(fileSize) + " Byte";
		}
		else if(fileSize/1024 > 1 && fileSize/1024 <= 1024){
			logging = to_string(fileSize/1024) + " KB";
		}
		else if(fileSize/1024/1024 > 1 && fileSize/1024/1024 <= 1024){
			logging = to_string(fileSize/1024/1024) + " MB";
		}
		else{
			logging = to_string(fileSize/1024/1024/1024) + " GB";
		}
		return logging;
	}

	// 폴더 안에 있는 파일들 크기
	int calcFileSize(string folder){
		int fileSize = 0;
		for (const filesystem::directory_entry& entry :
			filesystem::recursive_directory_iterator(filesystem::current_path() / folder)) {
			// std::cout << entry.path() << std::endl;
			if(filesystem::is_regular_file(entry.path())){
				fileSize += std::filesystem::file_size(entry.path());
			}
		}
		return fileSize;
	}

};


class MASCH_TimeSystem {
private:
	using time = chrono::system_clock::time_point;
	stack<time> start;
	stack<string> name;
	stack<int> level;
	int levelNow = 0;
	vector<string> logCalcTime;
	MASCH_MPI mpi;
public:
	MASCH_TimeSystem& push(string inp_name){
		name.push(inp_name);
		start.push(chrono::system_clock::now());
		level.push(levelNow++);
		return (*this);
	}
	MASCH_TimeSystem& pop(){
		chrono::microseconds calcTime = 
		chrono::duration_cast<chrono::microseconds>(
		chrono::system_clock::now() - start.top()); 
		start.pop();
		std::stringstream ss;
		ss << setw(level.top()*2) << setfill('_') << 
		" " << name.top() << " : " << 
		calcTime.count() << " us";;
		logCalcTime.push_back(ss.str());
		name.pop();
		level.pop();
		return (*this);
	}
	void show(){
		if(mpi.rank==0){
			for(auto& item : logCalcTime){
				cout << "| " << item << endl;
			}
		}
		logCalcTime.clear();
	}

	chrono::microseconds showCalcTime(){
		chrono::microseconds calcTime = 
		chrono::duration_cast<chrono::microseconds>(
		chrono::system_clock::now() - start.top()); 
		name.pop();
		start.pop();
		return calcTime;
	}
	
	string now(){
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
		ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
		return ss.str();
	}

};



class MASCH_Error {
private:
	MASCH_MPI mpi;
public:
	void stop(string inp,string inp0, string inp1, long inp2){
		if(mpi.rank==0){
			std::cout << "| #Error | " << 
			inp << ", file(" << inp0 << "), func(" << inp1 << "), line(" << inp2 << ")" << endl;;
		}
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		MPI::Finalize();
		exit(0);
	}

};


class MASCH_Warning {
private:
	vector<string> warn_;
	MASCH_MPI mpi;
public:
	void push(string inp,string inp0, string inp1, long inp2){
		std::stringstream ss;
		ss << inp << ", file(" << inp0 << "), func(" << inp1 << "), line(" << inp2 << ")";
		warn_.push_back(ss.str());
	}
	void pop(){
		warn_.pop_back();
	}
	void show(){
		if(mpi.size>1){
			vector<string> warn_glo = mpi.gatherv(warn_);
			if(mpi.rank==0) {
				std::cout << "| #Warning |" << endl;
				int ip=0;
				for(auto& item : warn_glo){
					if(!item.empty())
						std::cout << "| proc(" << ip++ << ") : " << item << endl;
				}
			}
		}
		else{
			std::cout << "| #Warning : " << std::endl;
			for(auto& item : warn_){
				std::cout << item << std::endl;
			}
		}
		warn_.clear();
	}

};

// ==============================
// 로그 클래스
class MASCH_Log {
private:
	string language = "eng";
	string state = "";
	string fileName = "";
	string functionName = "";
	string error_print = "";
	string log_print = "";
	long lineName;
	int level = 0;
public:

	MASCH_TimeSystem calcTime;
	MASCH_FileSystem file;
	MASCH_Error error;
	MASCH_Warning warning;

	MASCH_Log& setLevel(int inp_level) {
		level = inp_level;
		return (*this);
	}

};



int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Log log;
	
	MASCH_Control controls;
	
	
	MASCH_Load load;


	


	// 컨트롤 파일 읽어서 map 에 넣기
	load.extractFile("./setting/controlDict", controls.controlDictMap);
	
	// 기본적인 컨트롤 목록
	controls.language = controls.controlDictMap["language"];
	controls.startFrom = (controls.controlDictMap["startFrom"]);
	controls.stopAt = stod(controls.controlDictMap["stopAt"]);
	controls.timeStep = stod(controls.controlDictMap["timeStep"]);
	controls.adjustTimeStep = controls.controlDictMap["adjustTimeStep"];
	controls.maxCFL = stod(controls.controlDictMap["maxCFL"]);
	controls.maxTimeStep = stod(controls.controlDictMap["maxTimeStep"]);
	controls.againWhenDiverge = (controls.controlDictMap["againWhenDiverge"]);
	controls.multiCFL = stod(controls.controlDictMap["multiCFL"]);
	controls.minTimeStep = stod(controls.controlDictMap["minTimeStep"]);
	
	// 세이브 컨트롤
	controls.saveControl = (controls.controlDictMap["saveControl"]);
	controls.saveInTimeStep = 1e8;
	controls.saveInRunTime = 1.e8;
	if(controls.saveControl=="timeStep"){
		controls.saveInTimeStep = stoi(controls.controlDictMap["saveInterval"]);
	}
	else if(controls.saveControl=="runTime"){
		controls.saveInRunTime = stod(controls.controlDictMap["saveInterval"]);
	}
	else{
		log.warning.push("no defined saveControl",MASCH_FFL);
	}
	
	// 세이브 포맷
	controls.saveFormat = (controls.controlDictMap["saveFormat"]);
	controls.saveCompression = stoi(controls.controlDictMap["saveCompression"]);
	controls.writePrecision = stoi(controls.controlDictMap["writePrecision"]);
	
	// 세이브 시킬 데이터 목록
	controls.saveMeshData = load.extractVector(controls.controlDictMap["saveMeshData"]);
	controls.saveGradientData = load.extractVector(controls.controlDictMap["saveGradientData"]);
	controls.saveThermodynamicData = load.extractVector(controls.controlDictMap["saveThermodynamicData"]);
	controls.saveBodyForceData = load.extractVector(controls.controlDictMap["saveBodyForceData"]);
	

	// 화학종 파일 읽어서 map 에 넣기
	load.extractFile("./setting/species", controls.speciesMap);
	vector<string> species = load.extractVector(controls.speciesMap["name"]);
	// 화학종 갯수
	controls.nSp = species.size();
	
	
	// 열역학 파라미터 파일 읽어서 map 에 넣기
	load.extractFile("./setting/physics/thermophysicalProperties", controls.thermophysicalProperties);
	
	
	// 난류 파라미터 파일 읽어서 map 에 넣기
	load.extractFile("./setting/physics/turbulenceProperties", controls.turbulenceProperties);
	controls.turbType = controls.turbulenceProperties["type"];
	controls.LESModel = controls.turbulenceProperties["LES.LESModel"];
	controls.RANSModel = controls.turbulenceProperties["RANS.RANSModel"];
	// if(controls.turbType=="RANS"){
		// if(controls.RANSModel=="kEpsilon"){
			// load.extractFile("./setting/boundary/k", controls.boundary_k);
			// load.extractFile("./setting/boundary/epsilon", controls.boundary_epsilon);
		// }
		// else if(controls.RANSModel=="kOmega"){
			// load.extractFile("./setting/boundary/k", controls.boundary_k);
			// load.extractFile("./setting/boundary/omega", controls.boundary_omega);
		// }
	// }
	
	
	// 다이나믹메쉬(AMR...) 파일 읽어서 map 에 넣기
	load.extractFile("./setting/dynamicMesh", controls.dynamicMeshMap);
	
	// 시간에 따른 셀값 추출 파일 읽어서 map 에 넣기
	load.extractFile("./setting/extractDatasOverTime", controls.extractDatasOverTimeMap);
	
	// 수치 스킴 파일 읽어서 map 에 넣기
	load.extractFile("./setting/fvScheme", controls.fvSchemeMap);
	
	// 솔루션 컨트롤 파일 읽어서 map 에 넣기
	load.extractFile("./setting/fvSolution", controls.fvSolutionMap);
	
	
	// primitive 변수 셋팅
	controls.primitive.push_back(make_pair("pressure",0));
	controls.primitive.push_back(make_pair("x_velocity",1));
	controls.primitive.push_back(make_pair("y_velocity",2));
	controls.primitive.push_back(make_pair("z_velocity",3));
	controls.primitive.push_back(make_pair("temperature",4));
	for(int i=0; i<controls.nSp-1; ++i){
		controls.primitive.push_back(make_pair(species[i],5+i));
	}
	
	// primitive 약어 셋팅
	controls.primitive_abb.push_back("p");
	controls.primitive_abb.push_back("U");
	controls.primitive_abb.push_back("T");
	controls.primitive_abb.push_back("Y");
	// for(int i=0; i<controls.nSp-1; ++i){
		// string speabb = "Y";
		// speabb += to_string(i);
		// controls.primitive_abb.push_back(speabb);
	// }
	
	// primitive 역할 셋팅
	controls.primitive_role.push_back("scalar");
	controls.primitive_role.push_back("transport");
	controls.primitive_role.push_back("scalar");
	controls.primitive_role.push_back("vector");
	
	
	
	
	// 바운더리 컨디션 파일 읽어서 map 에 넣기
	// for(auto& [name, value] : controls.primitive){
	for(auto& name : controls.primitive_abb){
		string boundaryName = "./setting/boundary/";
		boundaryName += name;
		map<string, string> tmp_map;
		load.extractFile(boundaryName, tmp_map);
		controls.boundaryMap.push_back(tmp_map);
	}
	
	
	
	
	  
	string foldername;
	double starttime = stod(controls.startFrom);
	std::ostringstream streamObj;
	streamObj << starttime;
	foldername = "./save/" + streamObj.str() + "/";
	
	
	MASCH_Mesh mesh;
	MASCH_Mesh_Load mesh_load;
	mesh_load.vtu(foldername, mesh, controls);
	
	

	// vector<SEMO_Species> species;
	
	// SEMO_Controls_Builder controls;
		
	// controls.readSpecies(species);
	// controls.readConfigures();
	// // controls.setValues(species);
	
	// // // SEMO_Solvers_Builder solvers;
	
	// // SEMO_Mesh_Builder mesh;
	
	// // SEMO_MPI_Builder mpi;
	
	
	
	
	
	
	// bool boolLoad = true;
	// bool boolPartitioning = false;
	// bool boolAMR = true;
	// bool boolGeometric = true;
	
	// if(boolLoad){
		
		// SEMO_Mesh_Load load;
		// double starttime = stod(controls.startFrom);
		// string foldername;
		// std::ostringstream streamObj;
		// streamObj << starttime;
		// foldername = "./save/" + streamObj.str() + "/";
		
		// if(starttime == 0.0){
			// foldername = "./save/0/";
		// }
		
		// load.vtu(foldername, mesh, controls, species);
		
		
		// // //solvers.calcCellEOSVF(mesh, controls, species);
		// // solvers.calcCellEOSMF(mesh, controls, species);
		
		// // solvers.calcCellTransport(mesh, controls, species);
		
		// // for(auto& cell : mesh.cells){
			// // vector<double> Q(controls.nEq,0.0);
			// // Q[0] = cell.var[controls.Rho];
			// // Q[1] = cell.var[controls.Rho] * cell.var[controls.U];
			// // Q[2] = cell.var[controls.Rho] * cell.var[controls.V];
			// // Q[3] = cell.var[controls.Rho] * cell.var[controls.W];
			// // Q[4] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			// // for(int ns=0; ns<controls.nSp-1; ++ns){
				// // Q[5+ns] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
				
				// // // cout << cell.var[controls.MF[ns]] << endl;
			// // }
			// // for(int i=0; i<controls.nEq; ++i){
				// // cell.var[controls.Qm[i]] = Q[i];
				// // cell.var[controls.Qn[i]] = Q[i];
			// // }
		// // }
			
		// // for(auto& face : mesh.faces){
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// // face.var[controls.Un] = 0.5 * mesh.cells[face.owner].var[controls.U] * face.unitNormals[0];
				// // face.var[controls.Un] += 0.5 * mesh.cells[face.owner].var[controls.V] * face.unitNormals[1];
				// // face.var[controls.Un] += 0.5 * mesh.cells[face.owner].var[controls.W] * face.unitNormals[2];
				// // face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.U] * face.unitNormals[0];
				// // face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.V] * face.unitNormals[1];
				// // face.var[controls.Un] += 0.5 * mesh.cells[face.neighbour].var[controls.W] * face.unitNormals[2];
				
				// // face.var[controls.oldUn] = face.var[controls.Un];
			// // }
		// // }
	
	// }
	// else {
		
		// // load serial OpenFoam files
		// mesh.loadFile("OpenFoam", "./grid/");
	// }
	
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
	
	// // // partitioning
	// // if(boolPartitioning){
		
		// // mesh.distributeOneToAll("EVENLY");
		
	// // }
	
	
	
	// // if(boolAMR){
		
		// // mesh.hexaOctAMR();
		
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // }
	

	
	// // geometric
	// if(boolGeometric){
		
		// SEMO_Utility_Math math;
		
		// SEMO_Mesh_Geometric geometric;
		
		// geometric.init(mesh);
		
		// math.initLeastSquare(mesh);
	
	// }
	
	

	// bool initFlow = false;
	
	// // flow initialization
	// if(initFlow){
		
		// solvers.setInitValues(mesh, controls);
		
		// mesh.saveFile("vtu", "./save/0/", controls);
		
		// MPI_Barrier(MPI_COMM_WORLD);
	// }
	
	
	


	// bool calcFlow = true;
	
	// // flow calculation
	// if(calcFlow){
		
		// solvers.calcCellTransport(mesh, controls, species);
		
		
		// SEMO_Mesh_Save save;
		
		
		// while(
		// controls.iterReal<controls.iterRealMax ||
		// controls.time<stod(controls.stopAt) ){
			
	
			// if(rank==0) {
				// cout << "| real-time step = " << controls.iterReal 
				// << " | time = " << controls.time;
			// }
		
			// solvers.compressibleDensityBasedSingleTime(mesh, controls, species);
			
			// //==============================
			// // AMR
			// if(
			// // controls.iterReal != 0 && 
			// ( (controls.iterReal+1) % controls.intervalRefine == 0 ||
			  // (controls.iterReal+1) % controls.intervalUnrefine == 0)
			// ){ 
				// SEMO_Poly_AMR_Builder AMR;
				// AMR.polyAMR(mesh, controls, species, 0);
				// mesh.cellsGlobal();
				// solvers.calcCellEOSMF(mesh, controls, species);
				// solvers.calcCellTransport(mesh, controls, species);
			// } 
			// //==============================
				
			// controls.time += controls.timeStep;
			
			// ++controls.iterReal;
			
			// if(controls.saveControl == "timeStep"){
				// if(controls.iterReal % (int)controls.saveInterval == 0){
					// string foldername;
					// std::ostringstream streamObj;
					// streamObj << controls.time;
					// foldername = "./save/" + streamObj.str() + "/";
					
					// solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					// save.vtu(foldername, mesh, controls, species);
					
					// // save.particles(foldername, mesh, controls, species);
					
				// }
			// }
			// else if(controls.saveControl == "runTime"){
				// int jung = controls.time / controls.saveInterval;
				// double namuji = controls.time - (double)jung * controls.saveInterval;
				// if(
				// namuji < controls.timeStep &&
				// namuji >= 0.0
				// ){
					// string foldername;
					// std::ostringstream streamObj;
					// streamObj << controls.time;
					// foldername = "./save/" + streamObj.str() + "/";
					
					// solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					// save.vtu(foldername, mesh, controls, species);
					
					// // save.particles(foldername, mesh, controls, species);
					
				// }
			// }
		// }
		
	// }
	
	
	// if(rank==0) cout << "| End Program" << endl;
	// MPI_Barrier(MPI_COMM_WORLD);
	// //MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	MPI::Finalize();
	return EXIT_SUCCESS;
}














