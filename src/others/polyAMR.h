#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>
#include <set>
#include <tuple>

#include "./controls.h" 

// #include "build.h"
class MASCH_Mesh_Builder;

class faces_refind{

private:

public:
	int id;
	vector<int> ipoints;
	int iL;
	int iR;
	bool booliL = false;
	bool booliR = false;
};

class faces_proc_refind{

private:

public:
	vector<int> ipoints;
	vector<int> connPoints;
};



class groupMesh_Refine{

private:

public:
	int parent;
	vector<faces_refind> faces;
	vector<int> vertexPoints;
	vector<int> vertexCenterPoints;
	int centerPoint;
	int type;
	vector<vector<int>> subOutEdgePoints;
	vector<vector<int>> subIntEdgePoints;
};


class faces_Unrefind{

private:

public:
	int id;
	bool boolCombine;
	vector<int> ichild;
	vector<int> ipoints;
	int iL;
	int iR;
};



class groupCells_Unrefine{

private:

public:
	vector<int> ichild;
	// vector<int> outChildFaces;
	vector<int> intChildFaces;
	vector<int> groupOutFaces_id;
	vector<int> ifaces;
	vector<int> ipoints;
	vector<int> iL;
	int iR;
	int cellCenterPoint;
	bool boolCombine;
	vector<int> faceCenterPoints;
};












class MASCH_Poly_AMR_Builder {
private:

public:
    MASCH_Poly_AMR_Builder() {
    }
    ~MASCH_Poly_AMR_Builder() {
    }
	
	void polyAMR(
		MASCH_Mesh& mesh, 
		MASCH_Control& controls,
		MASCH_Solver& solver,
		MASCH_Variables& var,
		int iter);
	
	// void polyAMR_inline(
		// MASCH_Mesh& mesh, 
		// MASCH_Control& controls,
		// MASCH_Solver& solver,
		// MASCH_Variables& var,
		// int iter);
		

	void calcIndicators(
		MASCH_Mesh& mesh, 
		MASCH_Control& controls,
		MASCH_Variables& var,
		int maxBuffer,
		int maxLevel,
		int maxCells,
		int maxRefineCellPerBlockAMR,
		double minVolume_AMR,
		vector<vector<double>>& indicatorCriterion,
		vector<vector<int>>& indicatorAMR_id, 
		vector<bool>& boolCellRefine, 
		vector<bool>& boolCellUnrefine,
		vector<bool>& boolCellPreserved
		);
		
		
	// 리파인 관련
	void polyRefine(
		MASCH_Mesh& mesh, 
		MASCH_Control& controls,
		int maxLevel_AMR, int maxCells_AMR, double minVolume_AMR, 
		vector<vector<double>> indicatorCriterion,
		vector<vector<double>>& indicatorValues,
		vector<vector<int>>& child_new_cell_id_of_org,
		vector<bool>& boolCellPreserved,
		vector<bool>& boolCellRefine,
		vector<bool>& boolCellUnrefine,
		int iter);
		
	void mpiLevelRefine(
		MASCH_Mesh& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cLevel_recv, 
		vector<int>& cRefine_recv);
		
	void mpiRefines(
		MASCH_Mesh& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cRefine_recv);
		
	void mpiLevels(
		MASCH_Mesh& mesh, 
		vector<int>& cLevel_recv);
		
	void restrictCellRefine(
		MASCH_Mesh& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cLevel_recv, 
		vector<int>& cRefine_recv);
		
	void createEdges(
		MASCH_Mesh& mesh, 
		vector<int>& edgesPoint0,
		vector<int>& edgesPoint1, 
		vector<vector<int>>& facesEdges,
		vector<vector<int>>& edgesFaces,
		vector<int>& edgeLevel);
		
	void searchOriginalPoints(
		MASCH_Mesh& mesh, 
		vector<int>& points,
		int targetLevel, 
		vector<int>& originPoints);
		

	void addCenterPoint(
		MASCH_Mesh& mesh, 
		vector<int> vertex, 
		int level
		){
		mesh.addPoint();
		mesh.points.back().x = 0.0;
		mesh.points.back().y = 0.0;
		mesh.points.back().z = 0.0;
		double n = (double)vertex.size();
		for(auto& i : vertex){
			mesh.points.back().x += mesh.points[i].x/n;
			mesh.points.back().y += mesh.points[i].y/n;
			mesh.points.back().z += mesh.points[i].z/n;
		}
		mesh.points.back().level = level + 1;
	};
		
		

	void extractVertexPoints(
		int faceLevel,
		vector<int>& facePoints,
		vector<int>& facePointLevels,
		vector<int>& vertexPoints
	);


	void extractVertexCenterPoints(
		int faceLevel,
		vector<int>& facePoints,
		vector<int>& facePointLevels,
		vector<int>& vertexPoints,
		vector<int>& vertexCenterPoints,
		vector<bool>& addVertexCeterPoints
	);


	void extractSubEdgePoints(
		vector<int>& facePoints,
		vector<int>& vertexPoints,
		vector<int>& vertexCenterPoints,
		vector<vector<int>>& subEdgePoints
	) ;

	void addVertexCenterPoint(
		MASCH_Mesh& mesh, 
		int faceLevel,
		vector<int>& vertexPoints,
		vector<int>& vertexCenterPoints
	) ;
	void addFaceCenterPoint(
		MASCH_Mesh& mesh, 
		int faceLevel,
		vector<int>& vertexPoints,
		int& centerPoint
	);
	void addEdgeFacesPointsToVertexCenterPoint(
		MASCH_Mesh& mesh,
		int i,
		vector<vector<int>>& facesEdges,
		vector<vector<int>>& edgesFaces,
		vector<bool>& boolEdgeRefine,
		vector<bool>& addVertexCeterPoints,
		vector<int>& canRefineEdgeOrders,
		vector<int>& vertexPoints,
		vector<int>& vertexCenterPoints
	);
	void addCellCenterPoint(
		MASCH_Mesh& mesh, 
		int cellLevel,
		vector<int>& vertexPoints,
		int& centerPoint
	);
	void addSubOuterFacesPoints(
		vector<int>& vertexPoints,
		int& centerPoint,
		vector<vector<int>>& subOutEdgePoints,
		vector<faces_refind>& faces
	);
	void extractSubInternalEdgesPoints(
		vector<int>& vertexCenterPoints,
		int& centerPoint,
		vector<vector<int>>& subIntEdgePoints
	) ;
	void extractCellVertexOrder(
		vector<int>& cellVertexPoints,
		map<int,int>& cellVertexOrder
	) ;
	void extractInternalFacesVertexs(
		MASCH_Cell& cell,
		map<int,int>& cellInternalFaces,
		vector<int>& groupChildFaces_id,
		vector<groupMesh_Refine>& groupChildFaces,
		groupMesh_Refine& groupChildFace,
		int& cellCenterPoint,
		vector<vector<int>>& intFacesVertexs
	) ;
	void addInternalFacesiLiR(
		MASCH_Mesh& mesh, 
		int& cellLevel,
		int& totalCellNum,
		map<int,int>& cellVertexOrder,
		vector<vector<int>>& intFacesVertexs,
		vector<faces_refind>& faces
	);
	void addOuterFacesiLiR(
		MASCH_Mesh& mesh, 
		MASCH_Cell& cell,
		int& i,
		int& totalCellNum,
		map<int,int>& cellVertexOrder,
		vector<int>& groupChildFaces_id,
		vector<groupMesh_Refine>& groupChildFaces
	);
	void reorderOuterFacesiLiR(
		MASCH_Mesh& mesh, 
		MASCH_Cell& cell,
		int& i,
		int& totalCellNum,
		vector<int>& groupChildFaces_id,
		vector<bool>& boolInputFacesiL,
		vector<bool>& boolInputFacesiR, 
		vector<groupMesh_Refine>& groupChildFaces
	);
	
		
	// 언리파인 관련
	void polyUnrefine(
		MASCH_Mesh& mesh, 
		MASCH_Control& controls,
		int maxLevel_AMR, 
		vector<vector<double>> indicatorCriterion,
		vector<vector<double>>& indicatorValues,
		vector<vector<int>>& child_org_cell_id_of_new,
		vector<bool>& boolCellPreserved,
		vector<bool>& boolCellRefine,
		vector<bool>& boolCellUnrefine,
		int iter);
		
	void sortCellCanUnrefine(
		vector<int>& vLevel, 
		int saveLevel, 
		int& num, 
		vector<vector<int>>& vGroupLevel, 
		vector<vector<int>>& vGroupNumber);

	void restrictCellUnrefine(
		MASCH_Mesh& mesh, 
		vector<bool>& boolCellUnrefine,
		vector<int>& cLevel_recv, 
		vector<int>& cUnrefine_recv);


	void extractGroupCellListsCanUnrefine(
		MASCH_Mesh& mesh, 
		vector<vector<int>>& groupCellListsCanUnrefine);


	void extractGroupUnrefineCells(
		MASCH_Mesh& mesh, 
		vector<bool>& boolCellUnrefine,
		vector<vector<int>>& groupCellListsCanUnrefine,
		vector<groupCells_Unrefine>& groupCellsUnrefine);
		
	
		
};




class MASCH_Poly_Edge_Refine {
private:

public:
	int level;
	int point0;
	int point1;
	int centerPoint;
	vector<int> faces;
	
	
};




class MASCH_Poly_Mesh_Refine {
private:

public:
    MASCH_Poly_Mesh_Refine() {
    }
    ~MASCH_Poly_Mesh_Refine() {
    }
	
	MASCH_Poly_Mesh_Refine& addPoint(){
		MASCH_Point e;
		points.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addEdge(){
		MASCH_Poly_Edge_Refine e;
		edges.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addOutFace(int n){
		vector<MASCH_Face> e;
		for(int i=0; i<n; ++n){
			MASCH_Face a;
			e.push_back(a);
		}
		outFaces.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addIntFace(int n){
		vector<MASCH_Face> e;
		for(int i=0; i<n; ++n){
			MASCH_Face a;
			e.push_back(a);
		}
		intFaces.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addBdrFace(int n){
		vector<MASCH_Face> e;
		for(int i=0; i<n; ++n){
			MASCH_Face a;
			e.push_back(a);
		}
		bdrFaces.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addPrcFace(int n){
		vector<MASCH_Face> e;
		for(int i=0; i<n; ++n){
			MASCH_Face a;
			e.push_back(a);
		}
		prcFaces.push_back(e);
		return *this;
	}
	MASCH_Poly_Mesh_Refine& addCell(){
		MASCH_Cell e;
		cells.push_back(e);
		return *this;
	}
	
	void addCenterPoint(
		MASCH_Mesh& mesh, 
		vector<int> vertex, 
		int level
		){
		this->addPoint();
		this->points.back().x = 0.0;
		this->points.back().y = 0.0;
		this->points.back().z = 0.0;
		double n = (double)vertex.size();
		for(auto& i : vertex){
			this->points.back().x += mesh.points[i].x/n;
			this->points.back().y += mesh.points[i].y/n;
			this->points.back().z += mesh.points[i].z/n;
		}
		this->points.back().level = level + 1;
	}
	
	
	vector<MASCH_Point> points;
	vector<MASCH_Poly_Edge_Refine> edges;
	vector<MASCH_Cell> cells;
	
	vector<int> outFacesParentFace;
	vector<int> intFacesParentCell;
	vector<int> bdrFacesParentFace;
	vector<int> prcFacesParentFace;
	vector<vector<MASCH_Face>> outFaces;
	vector<vector<MASCH_Face>> intFaces;
	vector<vector<MASCH_Face>> bdrFaces;
	vector<vector<MASCH_Face>> prcFaces;
	
	
	void separateEdgesPoints(
		MASCH_Mesh& mesh, 
		int faceLevel,
		vector<int>& points, 
		vector<int>& edges,
		vector<bool>& boolEdgeRefine,
		vector<int>& edgesCenter,
		vector<int>& edgesLevel,
		vector<int>& faceVertex,
		vector<int>& edgeCenterPoints,
		vector<vector<int>>& subFaceEdgePoints,
		int& newPointNumber);
	
	
};



