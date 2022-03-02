#pragma once
#include <iostream>

#include <vector>
#include "../species/species.h"
#include "../mesh/mesh.h" 
#include "../controls/controls.h" 

class SEMO_Mesh_Builder;
class SEMO_Controls_Builder;

// template<typename T>
class SEMO_Mesh_Save {
private:

public:
    SEMO_Mesh_Save() {
    }
    ~SEMO_Mesh_Save() {
    }
	
	void mkdirs(char *dir_path);
	int Base64encode_len(int len);
	int Base64encode(char *encoded, const char *string, int len);
	
	void vtu(SEMO_Mesh_Builder &in);
		
	void vtu(string folder, int rank, SEMO_Mesh_Builder &in);
	// void vtuZlib(SEMO_Mesh_Builder &in, SEMO_Controls_Builder &controls);
	void gnuplot(int iter, double dtime, vector<double>& norm);
	
	template<typename T>
	void writeDatasAtVTU(
		SEMO_Controls_Builder &controls, ofstream& out, vector<T>& vecInp);

	template<typename T>
	void writeAscii(ofstream& out, vector<T>& vecInp);
		
	template<typename T>
	void writeBinary(ofstream& out, vector<T>& vecInp);
	
	template<typename T>
	void writeCompress(ofstream& out, vector<T>& vecInp, int compressSize);
	
	
	void vtu(
		string folder, 
		SEMO_Mesh_Builder &in, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species);
		
	void particles(
		string folder, 
		SEMO_Mesh_Builder &in, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species);
		
		
	
	// 포인트 추출
	void cellDataToMemoryOverTime(
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls);
	void cellDataOverTime(
		string folder, 
		SEMO_Controls_Builder &controls);

};




#include <filesystem>

class MASCH_Mesh;
class MASCH_Control;

class MASCH_Mesh_Save {
private:

public:
    MASCH_Mesh_Save() {
    }
    ~MASCH_Mesh_Save() {
    }
	
	// void mkdirs(char *dir_path);
	// int Base64encode_len(int len);
	// int Base64encode(char *encoded, const char *string, int len);
	
	// void vtu(MASCH_Mesh &in);
		
	void vtu(string folder, int rank, MASCH_Mesh &in);
	// // void vtuZlib(SEMO_Mesh_Builder &in, SEMO_Controls_Builder &controls);
	// void gnuplot(int iter, double dtime, vector<double>& norm);
	
	// template<typename T>
	// void writeDatasAtVTU(
		// MASCH_Control &controls, ofstream& out, vector<T>& vecInp);

	// template<typename T>
	// void writeAscii(ofstream& out, vector<T>& vecInp);
		
	// template<typename T>
	// void writeBinary(ofstream& out, vector<T>& vecInp);
	
	// template<typename T>
	// void writeCompress(ofstream& out, vector<T>& vecInp, int compressSize);
	
	
	// void vtu(
		// string folder, 
		// MASCH_Mesh &in, 
		// SEMO_Controls_Builder &controls,
		// vector<SEMO_Species>& species);
		
	// void particles(
		// string folder, 
		// MASCH_Mesh &in, 
		// SEMO_Controls_Builder &controls,
		// vector<SEMO_Species>& species);
		
		
	
	// // 포인트 추출
	// void cellDataToMemoryOverTime(
		// MASCH_Mesh &mesh, 
		// SEMO_Controls_Builder &controls);
	// void cellDataOverTime(
		// string folder, 
		// SEMO_Controls_Builder &controls);

};


