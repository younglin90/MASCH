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
#include "save.h"
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
	int nParcelVar=0;
	
	// 실행옵션
	map<string,string> mapArgv;
	
	// 원시변수
	vector<string> supPrimNames;
	vector<string> primAllScalarNames;
	vector<int> primAllScalarIds;
	vector<string> primVector3Names;
	vector<string> primVectorNames;
	vector<string> primScalarNames;
	
	// map<string, string> boundaryMap;
	map<string, string> physicsMap;
	map<string, string> controlDictMap;
	map<string, string> controlParcelsMap;
	map<string, string> dynamicMeshMap;
	map<string, string> extractDatasOverTimeMap;
	map<string, string> fvSchemeMap;
	map<string, string> fvSolutionMap;
	map<string, string> speciesMap;
	map<string, string> bodyforceMap;
	map<string, string> thermophysicalProperties;
	map<string, string> turbulenceProperties;
	map<string,map<string, string>> boundaryMap;
	map<string,map<string, string>> initialMap;
	
	string language, adjustTimeStep, againWhenDiverge, saveControl;
	string startFrom,saveFormat, turbType, LESModel, RANSModel;
	double stopAt, timeStep, maxCFL, maxTimeStep, multiCFL, minTimeStep;
	int saveInTimeStep, saveCompression, writePrecision;
	int nSp;
	double saveInRunTime;
	// vector<string> saveMeshData, saveGradientData, saveThermodynamicData, saveBodyForceData;
	vector<string> saveCellValues, saveGradientValues;
	
	
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
	map<string, MASCH_Control_Variable_Set> parcelVar;
	
	MASCH_Control_Variable cell;
	MASCH_Control_Variable face;
	MASCH_Control_Variable point;
	MASCH_Control_Variable boundary;
	MASCH_Control_Variable parcel;
	MASCH_Control_Variable field;
	
	vector<int> nEq;
	vector<vector<string>> sendProcValueNames;
	int nLinSys, nPrim, nSupPrim;
	vector<unsigned short> primVarIds;
	vector<string> primVarNames;
	vector<string> supPrimVarAbbs;
	vector<string> supPrimVarNames;
	vector<unsigned short> supPrimVarSizes;
	map<string,int> primitiveMap;
	
	// vector<function<double(double* inp)>> boundaryFunct;
	
	string getLoadFolderName(){
		string foldername;
		double starttime = stod((*this).startFrom);
		std::ostringstream streamObj;
		streamObj << starttime;
		foldername = "./save/" + streamObj.str() + "/";
		if(starttime == 0.0) foldername = "./save/0/";
		
		return foldername;
	}
	string getSaveFolderName(double inp_time){
		string foldername;
		std::ostringstream streamObj;
		streamObj << inp_time;
		streamObj.precision(12);
		foldername = "./save/" + streamObj.str() + "/";
		
		return foldername;
	}
	int getId_fieldVar(string inp_name){
		if(fieldVar.find(inp_name) != fieldVar.end()){
			return fieldVar[inp_name].id;
		}
		else{
			cout << "#WARNING, fieldVar not there, " << inp_name << endl;
		}
	}
	int getId_cellVar(string inp_name){
		if(cellVar.find(inp_name) != cellVar.end()){
			return cellVar[inp_name].id;
		}
		else{
			cout << "#WARNING, cellVar not there, " << inp_name << endl;
		}
	}
	int getId_faceVar(string inp_name){
		if(faceVar.find(inp_name) != faceVar.end()){
			return faceVar[inp_name].id;
		}
		else{
			cout << "#WARNING, faceVar not there, " << inp_name << endl;
		}
	}
	int getId_pointVar(string inp_name){
		if(pointVar.find(inp_name) != pointVar.end()){
			return pointVar[inp_name].id;
		}
		else{
			cout << "#WARNING, pointVar not there, " << inp_name << endl;
		}
	}
	int getId_parcelVar(string inp_name){
		if(parcelVar.find(inp_name) != parcelVar.end()){
			return parcelVar[inp_name].id;
		}
		else{
			cout << "#WARNING, parcelVar not there, " << inp_name << endl;
		}
	}
	
	
	// void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape);
	// void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape, vector<string> sub_name, vector<string> sub_abb);
	void setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape, vector<string> sub_name={""}, vector<string> sub_abb={""}, vector<string> sub_role={""});
	
	void setVariablesBasic();
	void setVariablesUDF(vector<string>& species);
	
	
	
	double time=0.0;
	double endTime=1.e8;
	double saveInterval;
	int iterReal, iterPseudo;
	bool checkContinueRealTimeStep(MASCH_Variables& var);
	bool updateRealTime();
	
	void updateRealTime(MASCH_Variables& var);
	
	void setPrimitiveValues();

	void setBoundaryConditions(MASCH_Mesh& mesh);
	// void setBoundaryFunctions(MASCH_Mesh& mesh, MASCH_Variables& var);
	void setVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var);
	void resetVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var,
		vector<vector<double>>& org_xyz, vector<vector<int>>& cellConn,
		string inp_option);
	void resetVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var,
		vector<int>& cell_ip, vector<int>& cellConn,
		vector<int>& parcel_ip,
		string inp_option);
	
	// 지오메트릭 관련
	void setGeometric(MASCH_Mesh& mesh, MASCH_Variables& var);
	void setGeometricOnlyCell_xyz(MASCH_Mesh& mesh);
	
	
	vector<double> limitMaxPrim, limitMinPrim;
	// void setMinMaxPrim();
	void limitPrim(int iEq, double& value);
	
	// AMR 관련
	int maxLevelRefine = 1;
	double minVolumeRefine = 1.e-16;
	int maxCellsRefine = 10000000;
	vector<double> indicatorCriterion;
	
	// species 관련
	vector<string> spName;
	
	// 세이브 관련
	void save_fvmFiles(MASCH_Mesh& mesh, MASCH_Variables& var);
	void save_dpmFiles(MASCH_Mesh& mesh, MASCH_Variables& var);
	void save_pvdFile(MASCH_Mesh& mesh, MASCH_Variables& var);
	
	// 체크 관련
	bool check_isnan(double value);
	
	// 초기화
	void saveAfterInitial(MASCH_Mesh& mesh);
	void saveAfterInitialAMR(MASCH_Mesh& mesh, MASCH_Variables& var);
	
	// 로그 관련
	void show_residual(MASCH_Variables& var);
	
	// parcel 관련
	int nIterDPM;
	vector<string> nameParcels;
	vector<int> calcDPM_iSeg;
	vector<pair<int,int>> idSetLagrangianEulerian;
	int nChangeParcelsE2L, nToProcsRishtParcels, nInsideParcels;
	int nReflectParcels, nEscapeParcels, nDeleteParcels;
	void show_dpm_information();
	
	// 리미터 관련
	vector<int> limiterNamesForUnst;
    
    
    // patch 바운더리 관련
    vector<double> timeVaryingMappedFixedValueNCycle;
    vector<double> timeVaryingMappedFixedValueTimeCycle;
    vector<double> timeVaryingMappedFixedValueTime1;
    vector<double> timeVaryingMappedFixedValueTime2;
    vector<int> timeVaryingMappedFixedValueTimeOrder;
    vector<vector<double>> timeVaryingMappedFixedValueTime;
    vector<vector<double>> timeVaryingMappedFixedValueValue1;
    vector<vector<double>> timeVaryingMappedFixedValueValue2;
    vector<int> timeVaryingMappedFixedValueValueIter;
    vector<vector<string>> timeVaryingMappedFixedValueFileName;
    
     
    // 평균값 관련
    vector<string> saveMeanCellValues;
    vector<string> saveSMDValues;
    vector<int> meanSurfInp_cell_id;
    vector<int> meanSurfOut_cell_id;
    vector<int> meanVolInp_cell_id;
    vector<int> meanVolOut_cell_id;
    vector<int> meanInp_cell_id;
    vector<int> meanInp_cell_totalTime_id;
    vector<int> meanOut_cell_id;
    vector<int> meanInp_parcel_id;
    vector<int> meanInp_parcel_totalTime_id;
    vector<int> meanOut_parcel_id;
    
    
};


	