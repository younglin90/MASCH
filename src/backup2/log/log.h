#pragma once
#include <vector>
#include <chrono>
#include <filesystem>
#include <stack>
#include <cstring>
using namespace std;

// #include "../mpi/mpi.h" 
#include "../test1/mpi.h"
// #include "../mesh/mesh.h"  
// #include "../controls/controls.h" 


class MASCH_MPI;

class MASCH_FileSystem {
private:

public:

	// 파일들 pretty하게 표시
	string showPrettySize(int fileSize);
	// 폴더 안에 있는 파일들 크기
	int calcFileSize(string folder);

};



class MASCH_TimeSystem {
private:
	using time = chrono::system_clock::time_point;
	stack<time> start;
	stack<string> name;
	stack<int> level;
	int levelNow = 0;
	vector<string> logCalcTime;
	// MASCH_MPI mpi;
public:
	MASCH_TimeSystem& push(string inp_name);
	MASCH_TimeSystem& pop();
	void show();

	chrono::microseconds showCalcTime();
	
	string now();

};




class MASCH_Error {
private:
	// MASCH_MPI mpi;
public:
	void stop(string inp,string inp0, string inp1, long inp2);

};


class MASCH_Warning {
private:
	vector<string> warn_;
	// MASCH_MPI mpi;
public:
	void push(string inp,string inp0, string inp1, long inp2);
	void pop();
	void show();

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

	MASCH_Log& setLevel(int inp_level);

};

