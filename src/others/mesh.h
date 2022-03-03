#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
#include <cmath>
#include <sstream>
#include <mpi.h>
#include "parmetis.h" 
#include "scotch.h" 
// #include <dlfnc.h>
using namespace std;

// #include "./load.h"
// #include "./save.h"
#include "./mpi.h"
#include "./math.h"
// #include "./math.h"

enum class MASCH_Face_Types{
	INTERNAL,
	BOUNDARY,
	PROCESSOR,
	TO_BE_INTERNAL,
	TO_BE_PROCESSOR,
	TO_BE_DELETE,
	NOT_AMR,
	NONE
};


class MASCH_Point{
public:
	MASCH_Point(){
		level=0;
	}
	int level;
	double x, y, z;
	
	vector<double> var;
	
	vector<int> istencil;
	
	// MPI 연결 포인트
	vector<pair<int,int>> connPoints;
};


class MASCH_Cell;

class MASCH_Face{
private:
	MASCH_Face_Types type_;
	MASCH_Cell* left_;
	MASCH_Cell* right_;
	vector<MASCH_Point*> points_;
	
public:
	MASCH_Face(){
		unitNormals.resize(3,0.0);
		level=0;
		iR=0;
		// thereR = false;
	}
	
	void setType(MASCH_Face_Types in){
		type_ = in;
	}
	MASCH_Face_Types& getType(){
		return type_;
	}

	int level;
	int group;
	vector<double> areaNormals;
	vector<double> unitNormals;
	double area;
	double wL;
	double distLR;
	double alpha;
	vector<double> vecLF;
	vector<double> vecRF;
	vector<double> vecLR;
	vector<double> unitNomalsLR;
	vector<double> vecSkewness;
	vector<double> vecLdL;
	vector<double> vecRdR;
	
	
	void setL(MASCH_Cell* inp){
		left_ = inp;
	}
	void setR(MASCH_Cell* inp){
		right_ = inp;
	}
	MASCH_Cell& L(){
		return *left_;
	}
	MASCH_Cell& R(){
		return *right_;
	}
	
	vector<MASCH_Point*>& setPoints(){
		return points_;
	}
	MASCH_Point& points(int id){
		return *points_[id];
	}
	
	int iL, iR;
	// bool thereR;
	vector<int> ipoints;
	double x, y, z;
	
	vector<double> values;
	
	vector<int> istencils;
		
		
};

class MASCH_Cell{
private:
	vector<MASCH_Face*> faces_;
	vector<MASCH_Point*> points_;

public:
	MASCH_Cell(){
		level=0;
	}
	vector<MASCH_Face*>& setFaces(){
		return faces_;
	}
	vector<MASCH_Face*>& faces(){
		return faces_;
	}
	MASCH_Face& faces(int id){
		return *faces_[id];
	}
	vector<MASCH_Point*>& setPoints(){
		return points_;
	}
	MASCH_Point& points(int id){
		return *points_[id];
	}
	
	int level;
	int group;
	vector<int> ipoints;
	vector<int> ifaces;
	
	double x, y, z;
	double volume;
	
	vector<int> iStencils;
	vector<int> recv_iStencils;
	vector<double> vars;
	
	int RCM;
	int invRCM;
		
};

class MASCH_Boundary{
private:
	MASCH_Face_Types type_;
	
public:
	MASCH_Boundary(){
		rightProcNo = -1;
		// thereR = false;
	}
	
	void setType(MASCH_Face_Types in){
		type_ = in;
	}
	MASCH_Face_Types& getType(){
		return type_;
	}
	
	typedef void (*setFunc_t)(double, double, double, double, double&);
	vector<setFunc_t> setFuncVars;
	
	// bool thereR;
	string name;
	vector<string> types;
	vector<double> values;
	int nFaces;
	int startFace;
	int myProcNo;
	int rightProcNo = -1;
};




class MASCH_Mesh{
private:
public:
	MASCH_Mesh() {}
	MASCH_Mesh& addPoint(){
		points.push_back(MASCH_Point());
		return *this;
	}
	MASCH_Mesh& addFace(){
		faces.push_back(MASCH_Face());
		return *this;
	}
	MASCH_Mesh& addCell(){
		cells.push_back(MASCH_Cell());
		return *this;
	}
	MASCH_Mesh& addBoundary(){
		boundaries.push_back(MASCH_Boundary());
		return *this;
	}
	
	MASCH_Point& point(int inp){
		return points[inp];
	}
	MASCH_Face& face(int inp){
		return faces[inp];
	}
	MASCH_Cell& cell(int inp){
		return cells[inp];
	}
	MASCH_Boundary& boundary(int inp){
		return boundaries[inp];
	}
	
	
	void check();
	void buildCells();
	void setFaceTypes();
	void setCountsProcFaces();
	void setDisplsProcFaces();
	void cellsGlobal();
	void informations();
	void connectCelltoFaces();
	void connectCelltoPoints();
	void connectFacetoPointsCells();
	void setCellStencils();
	void setFaceLevels();
	void setNumberOfFaces();
	void set(
		vector<double>& NodeCoordinates, vector<int>& connectivity, vector<int>& offsets, 
		vector<int>& faces, vector<int>& faceoffsets,
		vector<int>& owner, vector<int>& neighbour, vector<string>& bcName, vector<int>& bcStartFace, 
		vector<int>& bcNFaces, vector<int>& bcNeighbProcNo, vector<int>& connPoints,
		vector<int>& pointLevels, vector<int>& cellLevels, vector<int>& cellGroups);
	
	
	// mesh datas
	vector<MASCH_Point> points;
	vector<MASCH_Face> faces;
	vector<MASCH_Cell> cells;
	vector<MASCH_Boundary> boundaries;
	
	// 페이스 갯수 관련
	int nInternalFaces;
	int nBoundaryFaces;
	int nProcessorFaces;
	
	// 프로세스 관련
	vector<int> countsProcFaces;
	vector<int> displsProcFaces;
	
	vector<int> send_StencilCellsId;
	vector<int> send_countsStencilCells;
	vector<int> send_displsStencilCells;
	vector<int> recv_countsStencilCells;
	vector<int> recv_displsStencilCells;
	vector<double> recv_x_StencilCells;
	vector<double> recv_y_StencilCells;
	vector<double> recv_z_StencilCells;
	
	// 전체 도메인 관련
	int startCellGlobal;
	int ncellsTotal;
	vector<int> startProcCellGlobal;
	
	// 디버깅
	void debug_connPoints(double resi);
	void debug_procFacePoints(double resi);
	void debug_group_procFaces(double inp_resi);
	void debug_procFace_unitNomals(double inp_resi);
	
	// 로드밸런싱 리파티셔닝
	void repartitioning(
		vector<int>& idBlockCell, int maxLevel, vector<int>& connCelll);
	void repartParMETIS(
		int nSize, vector<int>& cell_ip, MASCH_Mesh &mesh);
	
	// 메쉬 형상 관련
	int nBoundTypes = 0;
	int nBoundFaces = 0;
	int nProcFaces = 0;
	int nInterFaces = 0;
	int nTriangle = 0;
	int nQuadrangle = 0;
	int nPolygon = 0;
	int nTetrahedron = 0;
	int nHexahedron = 0;
	int nPrism = 0;
	int nPyramid = 0;
	int nPolyhedron = 0;
		
};

