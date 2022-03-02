#pragma once
#include <iostream>
#include <filesystem>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <mpi.h>
#include <zlib.h>
// #include "../species/species.h"
// #include "../test1/mesh.h" 
#include "controls.h" 

using namespace std;

class MASCH_Control;

// class MASCH_Mesh;
// class MASCH_Control;

class MASCH_Mesh_Save {
private:

public:
    MASCH_Mesh_Save() {
    }
    ~MASCH_Mesh_Save() {
    }
	
	
	int Base64encode_len(int len);
	int Base64encode(char *encoded, const char *string, int len);
	
	void vtu(MASCH_Mesh &in);
	void vtu(string folder, int rank, MASCH_Mesh &in);
	// void gnuplot(int iter, double dtime, vector<double>& norm);
	
	template<typename T>
	void writeDatasAtVTU(
		MASCH_Control &controls, ofstream& out, vector<T>& vecInp);

	template<typename T>
	void writeAscii(ofstream& out, vector<T>& vecInp);
		
	template<typename T>
	void writeBinary(ofstream& out, vector<T>& vecInp);
	
	template<typename T>
	void writeCompress(ofstream& out, vector<T>& vecInp, int compressSize);
	
	
	void vtu(
		string folder, 
		MASCH_Mesh &in, 
		MASCH_Control &controls);
		
	void fvmFiles(string folder, int myRank, MASCH_Mesh& in, 
		MASCH_Control& controls, MASCH_Variables& var);
		
	// void particles(
		// string folder, 
		// MASCH_Mesh &in, 
		// MASCH_Control &controls);
		
		
	
	// // 포인트 추출
	// void cellDataToMemoryOverTime(
		// SEMO_Mesh_Builder &mesh, 
		// SEMO_Controls_Builder &controls);
	// void cellDataOverTime(
		// string folder, 
		// SEMO_Controls_Builder &controls);

	
	
};


