
#include "./load.h" 

void MASCH_Load::settingFiles(string folderName, MASCH_Control& controls){

	auto& load = (*this);
	
	// 컨트롤 파일 읽어서 map 에 넣기
	string controlDictFile = folderName + "/controlDict";
	load.extractFile(controlDictFile, controls.controlDictMap);
	
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
		// controls.log.warning.push("no defined saveControl",MASCH_FFL);
	}
	
	// 세이브 포맷
	controls.saveFormat = (controls.controlDictMap["saveFormat"]);
	controls.saveCompression = stoi(controls.controlDictMap["saveCompression"]);
	controls.writePrecision = stoi(controls.controlDictMap["writePrecision"]);
	
	// 세이브 시킬 데이터 목록
	// controls.saveMeshData = load.extractVector(controls.controlDictMap["saveMeshData"]);
	// controls.saveGradientData = load.extractVector(controls.controlDictMap["saveGradientData"]);
	// controls.saveThermodynamicData = load.extractVector(controls.controlDictMap["saveThermodynamicData"]);
	// controls.saveBodyForceData = load.extractVector(controls.controlDictMap["saveBodyForceData"]);
	// controls.saveFaceValues = load.extractVector(controls.controlDictMap["saveFaceValues"]);
	controls.saveCellValues = load.extractVector(controls.controlDictMap["saveCellValues"]);
	controls.saveGradientValues = load.extractVector(controls.controlDictMap["saveGradientValues"]);
    
    
	controls.saveSMDValues = load.extractVector(controls.controlDictMap["saveSMDValues"]);
	controls.saveMeanCellValues = load.extractVector(controls.controlDictMap["saveMeanCellValues"]);

	// 화학종 파일 읽어서 map 에 넣기
	string speciesFile = folderName + "/species";
	load.extractFile(speciesFile, controls.speciesMap);
	vector<string> species = load.extractVector(controls.speciesMap["name"]);
	// 화학종 갯수
	controls.nSp = species.size();
	controls.spName = species;
	
	
	// 바디포스
	string bodyforceFile = folderName + "/physics/bodyforce";
	load.extractFile(bodyforceFile, controls.bodyforceMap);
	// 열역학 파라미터 파일 읽어서 map 에 넣기
	string thermophysicalPropertiesFile = folderName + "/physics/thermophysicalProperties";
	load.extractFile(thermophysicalPropertiesFile, controls.thermophysicalProperties);
	
	
	// 난류 파라미터 파일 읽어서 map 에 넣기
	string turbulencePropertiesFile = folderName + "/physics/turbulenceProperties";
	load.extractFile(turbulencePropertiesFile, controls.turbulenceProperties);
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
	string dynamicMeshfile = folderName + "/dynamicMesh";
	load.extractFile(dynamicMeshfile, controls.dynamicMeshMap);
	
	// 시간에 따른 셀값 추출 파일 읽어서 map 에 넣기
	string extractDatasOverTimefile = folderName + "/extractDatasOverTime";
	load.extractFile(extractDatasOverTimefile, controls.extractDatasOverTimeMap);
	
	// 수치 스킴 파일 읽어서 map 에 넣기
	string fvSchemefile = folderName + "/fvScheme";
	load.extractFile(fvSchemefile, controls.fvSchemeMap);
	
	// 솔루션 컨트롤 파일 읽어서 map 에 넣기
	string fvSolutionFile = folderName + "/fvSolution";
	load.extractFile(fvSolutionFile, controls.fvSolutionMap);
	
	// 컨트롤 파일 읽어서 map 에 넣기
	string controlParcelsFile = folderName + "/controlParcels";
	load.extractFile(controlParcelsFile, controls.controlParcelsMap);
	controls.nameParcels = load.extractVector(controls.controlParcelsMap["name"]);
	
	
	//=====================================
	// using var = MASCH_Control_Variable_Set;
	
	controls.setVariablesUDF(species);
	
	controls.setPrimitiveValues();
	
	
	// 바운더리 컨디션 파일 읽어서 map 에 넣기
	for(auto& inp : controls.supPrimNames){
		string boundaryName = folderName + "/boundary/";
		boundaryName += controls.cellVar[inp].abb;
		map<string, string> tmp_map;
		load.extractFile(boundaryName, tmp_map);
		controls.boundaryMap.insert(make_pair(inp,tmp_map));
	}
	
	
	// 초기화 파일 
	for(auto& inp : controls.supPrimNames){
		string initialName = folderName + "/initial/";
		initialName += controls.cellVar[inp].abb;
		map<string, string> tmp_map;
		load.extractFile(initialName, tmp_map);
		controls.initialMap.insert(make_pair(inp,tmp_map));
		// controls.initialMap.push_back(tmp_map);
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}









