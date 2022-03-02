#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
#include <functional>
using namespace std;

#include "mesh.h"
#include "load.h"
#include "variables.h"
#include "log.h"
// #include "controls.h"
// #include "../species/species.h"


class MASCH_Variables;
// class MASCH_Control;

class MASCH_Control_Variable {
public:
	vector<string> name;
	int size;
	int varSize;
};


class MASCH_Control_Variable_Set {
public:
	MASCH_Control_Variable_Set(){};
	MASCH_Control_Variable_Set(
		int inp_id, string inp_name, string inp_abb, string inp_unit, string inp_role,
		string inp_shape, vector<string> inp_sub_name, vector<string> inp_sub_abb, vector<string> inp_sub_role){
		
		id = inp_id;
		name = inp_name;
		abb = inp_abb;
		unit = inp_unit;
		role = inp_role;
		shape = inp_shape;
		sub_name = inp_sub_name;
		sub_abb = inp_sub_abb;
		sub_role = inp_sub_role;
	}
	// string name;
	
	int id=-1;
	string name;
	string abb;
	string unit;
	string role;
	string shape;
	vector<string> sub_name;
	vector<string> sub_abb;
	vector<string> sub_role;
	
};



class MASCH_Control {
private:
public:
	MASCH_TimeSystem log;

	int nFieldVar=0;
	int nCellVar=0;
	int nFaceVar=0;
	int nPointVar=0;
	
	// map<string, string> boundaryMap;
	map<string, string> physicsMap;
	map<string, string> controlDictMap;
	map<string, string> dynamicMeshMap;
	map<string, string> extractDatasOverTimeMap;
	map<string, string> fvSchemeMap;
	map<string, string> fvSolutionMap;
	map<string, string> speciesMap;
	map<string, string> thermophysicalProperties;
	map<string, string> turbulenceProperties;
	vector<map<string, string>> boundaryMap;
	vector<map<string, string>> initialMap;
	
	string language, adjustTimeStep, againWhenDiverge, saveControl;
	string startFrom,saveFormat, turbType, LESModel, RANSModel;
	double stopAt, timeStep, maxCFL, maxTimeStep, multiCFL, minTimeStep;
	int saveInTimeStep, saveCompression, writePrecision;
	int nSp;
	double saveInRunTime;
	vector<string> saveMeshData, saveGradientData, saveThermodynamicData, saveBodyForceData;
	
	
	map<string, int> varDict;
	vector<string> species;
	// vector<pair<string, int>> variables;
	// vector<pair<string, int>> primitive;
	vector<string> primitive_abb;
	vector<string> primitive_role;
	
	
	map<string, MASCH_Control_Variable_Set> fieldVar;
	map<string, MASCH_Control_Variable_Set> cellVar;
	map<string, MASCH_Control_Variable_Set> faceVar;
	map<string, MASCH_Control_Variable_Set> pointVar;
	
	MASCH_Control_Variable cell;
	MASCH_Control_Variable face;
	MASCH_Control_Variable point;
	MASCH_Control_Variable boundary;
	MASCH_Control_Variable parcel;
	MASCH_Control_Variable field;
	
	int nEq, nLinSys, nPrim, nSupPrim;
	vector<unsigned short> primVarIds;
	vector<string> primVarNames;
	vector<string> supPrimVarAbbs;
	vector<string> supPrimVarNames;
	vector<unsigned short> supPrimVarSizes;
	map<string,int> primitiveMap;
	
	// vector<function<double(double* inp)>> boundaryFunct;
	
	string getFolderName(){
		string foldername;
		double starttime = stod((*this).startFrom);
		std::ostringstream streamObj;
		streamObj << starttime;
		foldername = "./save/" + streamObj.str() + "/";
		if(starttime == 0.0) foldername = "./save/0/";
		
		return foldername;
	}
	
	
	// void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape);
	// void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape, vector<string> sub_name, vector<string> sub_abb);
	void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape, vector<string> sub_name={""}, vector<string> sub_abb={""}, vector<string> sub_role={""});
	
	void setVariablesBasic();
	void setVariablesUDF(vector<string>& species);
	
	
	
	double time=0.0;
	double endTime=1.e8;
	double saveInterval;
	int iterReal;
	bool checkContinueRealTimeStep(MASCH_Variables& var);
	bool updateRealTime();
	
	void updateRealTime(MASCH_Variables& var);
	
	void setPrimitiveValues();

	void setBoundaryConditions(MASCH_Mesh& mesh);
	// void setBoundaryFunctions(MASCH_Mesh& mesh, MASCH_Variables& var);
	void setVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var);
	void setGeometric(MASCH_Mesh& mesh, MASCH_Variables& var);
	
	vector<double> limitMaxPrim, limitMinPrim;
	void setMinMaxPrim();
	void limitPrim(int iEq, double& value);
	
	// AMR 관련
	int maxLevelRefine = 1;
	double minVolumeRefine = 1.e-16;
	int maxCellsRefine = 10000000;
	vector<double> indicatorCriterion;
};


	