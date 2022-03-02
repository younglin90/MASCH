

using namespace std;

#include "./load.h" 
#include <dlfcn.h>


/* aaaack but it's fast and const should make it shared text page. */
static const unsigned char pr2six[256] =
{
    /* ASCII table */
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 62, 64, 64, 64, 63,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 64, 64, 64, 64, 64, 64,
    64,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 64, 64, 64, 64, 64,
    64, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
};


//앞에 있는 개행 문자 제거 
string &MASCH_Mesh_Load::ltrim(std::string &s) { 
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
	return s; 
}

//뒤에 있는 개행 문자 제거 
string &MASCH_Mesh_Load::rtrim(std::string &s) { 
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end()); 
	return s; 
}

//양쪽 끝의 개행 문자 제거 
string &MASCH_Mesh_Load::trim(std::string &s) { 
	return ltrim(rtrim(s)); 
}
	





void MASCH_Load::settingFiles(string folderName, MASCH_Control& controls){

	auto& load = (*this);
	
	// MASCH_Log log;
	// // MASCH_Control controls;
	
	
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
	controls.saveMeshData = load.extractVector(controls.controlDictMap["saveMeshData"]);
	controls.saveGradientData = load.extractVector(controls.controlDictMap["saveGradientData"]);
	controls.saveThermodynamicData = load.extractVector(controls.controlDictMap["saveThermodynamicData"]);
	controls.saveBodyForceData = load.extractVector(controls.controlDictMap["saveBodyForceData"]);

	// 화학종 파일 읽어서 map 에 넣기
	string speciesFile = folderName + "/species";
	load.extractFile(speciesFile, controls.speciesMap);
	vector<string> species = load.extractVector(controls.speciesMap["name"]);
	// 화학종 갯수
	controls.nSp = species.size();
	
	
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
	
	
	//=====================================
	// using var = MASCH_Control_Variable_Set;
	
	controls.setVariablesUDF(species);
	
	controls.setPrimitiveValues();
	
	
	// 바운더리 컨디션 파일 읽어서 map 에 넣기
	for(auto& inp_abb : controls.supPrimVarAbbs){
		string abb = inp_abb;
		string boundaryName = folderName + "/boundary/";
		boundaryName += abb;
		map<string, string> tmp_map;
		load.extractFile(boundaryName, tmp_map);
		controls.boundaryMap.push_back(tmp_map);
	}
	
	
	// 초기화 파일 
	for(auto& inp_abb : controls.supPrimVarAbbs){
		string abb = inp_abb;
		string initialName = folderName + "/initial/";
		initialName += abb;
		map<string, string> tmp_map;
		load.extractFile(initialName, tmp_map);
		controls.initialMap.push_back(tmp_map);
	}
	
	
}












template void MASCH_Mesh_Load::loadDatasAtVTU<int>(
	ifstream& inputFile, string dataName, vector<int>& outData);
template void MASCH_Mesh_Load::loadDatasAtVTU<double>(
	ifstream& inputFile, string dataName, vector<double>& outData);
template void MASCH_Mesh_Load::loadDatasAtVTU<string>(
	ifstream& inputFile, string dataName, vector<string>& outData);
template<typename T>
void MASCH_Mesh_Load::loadDatasAtVTU(
	ifstream& inputFile, string dataName, vector<T>& outData){
	
	string nextToken;
	string combDataName = "\"" + dataName + "\"";
	trim(combDataName);
	
	inputFile.clear();
	inputFile.seekg( 0, std::ios_base::beg);
	outData.clear();
	
	bool startValue = false;
	bool boolBinary = false;
	while(getline(inputFile, nextToken)){
		
		string asignToken;

		if(startValue){
			if(nextToken.find("</DataArray>") != string::npos){
				startValue=false;
				boolBinary=false;
				break;
			}
			else{
				// istringstream iss(nextToken);
				// T tempint;
				// while(iss >> tempint){
					// outData.push_back(tempint);
				// }
				stringstream sstream(nextToken);
				string word;
				
				
				char del = ' ';
				if(boolBinary==false){
					while (getline(sstream, word, del)){
						istringstream iss(word);
						T tempint;
						iss >> tempint;
						outData.push_back(tempint);
					}
				}
				else{
					if(boolCompress==false){
						while (getline(sstream, word, del)){
							readBinary(word, outData);
						}
					}
					else{
						while (getline(sstream, word, del)){
							readCompress(word, outData);
						}
					}
				}
				
			}
		}
		else{
			if( nextToken.find(combDataName) != string::npos ){
				startValue=true;
				if( nextToken.find("format=\"binary\"") != string::npos ){
					boolBinary=true;
				}
			}
		}
		
	}
	
}

void MASCH_Mesh_Load::vtu(
	string folder, 
	MASCH_Control &controls, 
	MASCH_Mesh &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load vtu data files ... ";
	}
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	// points
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
	
	vector<int> pointLevels;
	loadDatasAtVTU(inputFile, "pointLevels", pointLevels);
	
	vector<int> cellLevels;
	loadDatasAtVTU(inputFile, "cellLevels", cellLevels);
	
	vector<int> cellGroups;
	loadDatasAtVTU(inputFile, "cellGroups", cellGroups);
	
	vector<double> NodeCoordinates;
	loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	
	vector<int> connectivity;
	loadDatasAtVTU(inputFile, "connectivity", connectivity);
	
	vector<int> offsets;
	loadDatasAtVTU(inputFile, "offsets", offsets);
	
	vector<int> faces;
	loadDatasAtVTU(inputFile, "faces", faces);
	
	vector<int> faceoffsets;
	loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
	
	vector<int> owner;
	loadDatasAtVTU(inputFile, "owner", owner);
	
	vector<int> neighbour;
	loadDatasAtVTU(inputFile, "neighbour", neighbour);
	
	vector<string> bcName;
	loadDatasAtVTU(inputFile, "bcName", bcName);
	
	vector<int> bcStartFace;
	loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
	
	vector<int> bcNFaces;
	loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
	
	vector<int> bcNeighbProcNo;
	loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);

	vector<int> connPoints;
	loadDatasAtVTU(inputFile, "connPoints", connPoints);

	inputFile.close();
	
	

	
	
	for(int i=0; i<NodeCoordinates.size()/3; ++i){
		mesh.addPoint();
		mesh.points.back().x = NodeCoordinates[i*3+0];
		mesh.points.back().y = NodeCoordinates[i*3+1];
		mesh.points.back().z = NodeCoordinates[i*3+2];
	}
	NodeCoordinates.clear();
	
	
	int n=0;
	int ncells=0;
	mesh.faces.clear();
	for(auto& i : owner){
		mesh.addFace();
		mesh.faces.back().iL = i;
		ncells = max(ncells, mesh.faces.back().iL);
	}
	owner.clear();
	
	// cout << ncells << endl;
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();
	}
	
	
	
	n=0;
	for(auto& i : neighbour){
		mesh.faces[n].iR = i;
		mesh.faces[n].setType(MASCH_Face_Types::INTERNAL);
		++n;
	}
	neighbour.clear();
	
	
	int m=0;
	n=0;
	for(auto& i : offsets){
		for(int j=n; j<i; ++j){
			int point = connectivity[j];
			mesh.cells[m].ipoints.push_back( static_cast<int>(point) );
		}
		n=i;
		++m;
	}
	
	
	
	n=0;
	int nFacesInt=0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
			mesh.cells[ face.iR ].ifaces.push_back( static_cast<int>(n) );
			++nFacesInt;
		}
		else{
			mesh.cells[ face.iL ].ifaces.push_back( static_cast<int>(n) );
		}
		++n;
	}
	
	
	m=0;
	n=0;
	for(auto& i : faceoffsets){
		int N=0;
		int face_size = faces[m+N];
		for(int j=0; j<face_size; ++j){
			int face = mesh.cells[n].ifaces[j];
			++N;
			int point_size = faces[m+N];
			for(int k=0; k<point_size; ++k){
				++N;
				int point = faces[m+N];
				if(mesh.faces[ face ].ipoints.size() == point_size) continue;
				mesh.faces[ face ].ipoints.push_back( static_cast<int>(point) );
			}
		}
		m=i;
		++n;
	}
	faces.clear();
	faceoffsets.clear();
	
	
	
	n=0;
	for(auto& startFace : bcStartFace){
		
		mesh.addBoundary();
		
		string tmp_bcName = bcName[n];
		tmp_bcName.erase(std::find_if(tmp_bcName.rbegin(), 
		tmp_bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), tmp_bcName.end());
		tmp_bcName.erase(tmp_bcName.begin(), std::find_if(
		tmp_bcName.begin(), tmp_bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
		mesh.boundaries.back().name = tmp_bcName;
		mesh.boundaries.back().startFace = bcStartFace[n];
		mesh.boundaries.back().nFaces = bcNFaces[n];
		if(bcNeighbProcNo[n]<0){
			mesh.boundaries.back().rightProcNo = 0;
			mesh.boundaries.back().setType(MASCH_Face_Types::BOUNDARY);
		}
		else{
			mesh.boundaries.back().rightProcNo = bcNeighbProcNo[n];
			mesh.boundaries.back().setType(MASCH_Face_Types::PROCESSOR);
		}
		mesh.boundaries.back().myProcNo = rank;
		
		++n;
		
	}
	
	int maxBCNum = mesh.boundaries.size()-1;
	mesh.boundaries[maxBCNum].startFace = mesh.faces.size()-mesh.boundaries[maxBCNum].nFaces;
	for(int i=maxBCNum-1; i>=0; --i){
		mesh.boundaries[i].startFace = mesh.boundaries[i+1].startFace-mesh.boundaries[i].nFaces;
	}
	
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
	for(int i=0; i<connPoints.size(); ++i){
		mesh.points[connPoints[i]].connPoints.push_back(
		make_pair(connPoints[i+1],connPoints[i+2])
		);
		++i; ++i;
	}
	
	
	
	// mesh.set(NodeCoordinates, connectivity, offsets, faces, faceoffsets,
		// owner, neighbour, bcName, bcStartFace, bcNFaces, bcNeighbProcNo, connPoints,
		// pointLevels, cellLevels, cellGroups);
	
		
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	mesh.check();
	mesh.setFaceTypes();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	mesh.setFaceLevels();
	mesh.setCellStencils();
	mesh.setNumberOfFaces();
	mesh.informations();

		
}



















void MASCH_Load::meshFiles(string folderName, MASCH_Control& controls, MASCH_Mesh& mesh){
	
	MASCH_Mesh_Load mesh_load;
	
	
	
	mesh_load.vtu(folderName, controls, mesh);
	
	
	
	
}




















int MASCH_Mesh_Load::Base64decode_len(const char *bufcoded)
{
    int nbytesdecoded;
    register const unsigned char *bufin;
    register int nprbytes;

    bufin = (const unsigned char *) bufcoded;
    while (pr2six[*(bufin++)] <= 63);

    nprbytes = (bufin - (const unsigned char *) bufcoded) - 1;
    nbytesdecoded = ((nprbytes + 3) / 4) * 3;

    return nbytesdecoded + 1;
}

int MASCH_Mesh_Load::Base64decode(char *bufplain, const char *bufcoded)
{
    int nbytesdecoded;
    register const unsigned char *bufin;
    register unsigned char *bufout;
    register int nprbytes;

    bufin = (const unsigned char *) bufcoded;
    while (pr2six[*(bufin++)] <= 63);
    nprbytes = (bufin - (const unsigned char *) bufcoded) - 1;
    nbytesdecoded = ((nprbytes + 3) / 4) * 3;

    bufout = (unsigned char *) bufplain;
    bufin = (const unsigned char *) bufcoded;

    while (nprbytes > 4) {
    *(bufout++) =
        (unsigned char) (pr2six[*bufin] << 2 | pr2six[bufin[1]] >> 4);
    *(bufout++) =
        (unsigned char) (pr2six[bufin[1]] << 4 | pr2six[bufin[2]] >> 2);
    *(bufout++) =
        (unsigned char) (pr2six[bufin[2]] << 6 | pr2six[bufin[3]]);
    bufin += 4;
    nprbytes -= 4;
    }

    /* Note: (nprbytes == 1) would be an error, so just ingore that case */
    if (nprbytes > 1) {
    *(bufout++) =
        (unsigned char) (pr2six[*bufin] << 2 | pr2six[bufin[1]] >> 4);
    }
    if (nprbytes > 2) {
    *(bufout++) =
        (unsigned char) (pr2six[bufin[1]] << 4 | pr2six[bufin[2]] >> 2);
    }
    if (nprbytes > 3) {
    *(bufout++) =
        (unsigned char) (pr2six[bufin[2]] << 6 | pr2six[bufin[3]]);
    }

    *(bufout++) = '\0';
    nbytesdecoded -= (4 - nprbytes) & 3;
    return nbytesdecoded;
}















template void MASCH_Mesh_Load::readBinary<int>(
	string& word, vector<int>& outData);
template void MASCH_Mesh_Load::readBinary<double>(
	string& word, vector<double>& outData);
template void MASCH_Mesh_Load::readBinary<string>(
	string& word, vector<string>& outData);
template<typename T>
void MASCH_Mesh_Load::readBinary(
	string& word, vector<T>& outData){
		
	outData.clear();
		
	int datasize = 8;
	int dataByteSize = sizeof(T);
		
	const char *data_in = word.c_str();
		
	int decoded_data_length = Base64decode_len(data_in);
	char* data_out = (char*)malloc(decoded_data_length);
	
	Base64decode(data_out, data_in);
	// printf("The string\n[%s]\ndecodes from base64 as:\n[%s]\n", data_in, data_out);
	
		// cout << sizeof(data_in) << endl;
		// cout << decoded_data_length << " " << sizeof(data_out) << endl;
	
	char buffer1[datasize];
	
	int pointer_end = 0;
	pointer_end += datasize;
	
	std::copy(data_out, data_out + pointer_end, buffer1);
	long long total_byte;
	std::memcpy( &total_byte, buffer1, datasize );
	int data_length = total_byte/8;
	
	// cout << data_length << endl;
	
	for(int i=0; i<data_length; ++i){
		int pointer_str = pointer_end;
		pointer_end += dataByteSize;
		
		char buffer2[dataByteSize];
		std::copy(data_out + pointer_str, data_out + pointer_end, buffer2);
		T data;
		std::memcpy( &data, buffer2, dataByteSize );
		
		// cout << data << " ";
		
		outData.push_back(data);
	}
			
			
	free(data_out);
	
	
}



















template void MASCH_Mesh_Load::readCompress<int>(
	string& word, vector<int>& outData);
template void MASCH_Mesh_Load::readCompress<double>(
	string& word, vector<double>& outData);
template void MASCH_Mesh_Load::readCompress<string>(
	string& word, vector<string>& outData);
template<typename T>
void MASCH_Mesh_Load::readCompress(
	string& word, vector<T>& outData){
		
	// 데이터 클리어
	outData.clear();
	
	// string 나누기, 44 는 header 사이즈
	string header_string = word.substr(0, 44);
	string value_string = word.substr(44, word.size()-44);
	
	// 인풋 데이터 캐릭터형으로 복사
	const char *header_in = header_string.c_str();
	const char *value_in = value_string.c_str();
	
	// 헤더 base64 디코딩 후 저장
	int header_byte_size = 8;
	int header_encode_length = header_string.size();
	char header_in_char[header_encode_length];
	std::copy(header_in, header_in + header_encode_length, header_in_char);
	
	int header_decoded_length = Base64decode_len(header_in_char);
	char* header_decoded = (char*)malloc(header_decoded_length);
	Base64decode(header_decoded, header_in_char);
	
	// 헤더 부분 저장
	int pointer_end = 0;
	vector<long long> header;
	for(int i=0; i<4; ++i){
		int pointer_str = pointer_end;
		pointer_end += header_byte_size;
		
		char header_buffer[header_byte_size];
		std::copy(header_decoded + pointer_str, header_decoded + pointer_end, header_buffer);
		long long header_tmp;
		std::memcpy( &header_tmp, header_buffer, header_byte_size );
		header.push_back(header_tmp);
	}
	free(header_decoded);
	
	// 데이터 부분만 base64 디코딩 후 저장
	int value_encode_length = value_string.size();
	char* value_in_char = (char*)malloc(value_encode_length);
	std::copy(value_in, value_in + value_encode_length, value_in_char);
	int value_decoded_length = Base64decode_len(value_in_char);
	char* value_decoded = (char*)malloc(value_decoded_length);
	Base64decode(value_decoded, value_in_char);
	free(value_in_char);

	// zlib 로 압축해제
	Bytef* value_decoded_byte = (Bytef*)malloc(value_decoded_length);
	uLong value_sizee = value_decoded_length;
	char* value_buffer = (char*)malloc(value_decoded_length);
	std::copy(value_decoded, value_decoded + value_decoded_length, value_buffer);
	std::memcpy( value_decoded_byte, value_buffer, value_decoded_length );
	uLong inflate_size = header[1];
	Bytef *inflate_data = (Bytef*)malloc(inflate_size);
	
	int error = uncompress2(inflate_data, &inflate_size, value_decoded_byte, &value_sizee);
	if(error != Z_OK){
		cout << "zlib error !!!!!!!!!!!!" << endl;
		cout << error << endl;
	}
	free(value_decoded);
	free(value_buffer);
	
	// 데이터 부분 저장
	int size_value_byte = sizeof(T);
	pointer_end = 0;
	for(int i=0; i<inflate_size/size_value_byte; ++i){
		int pointer_str = pointer_end;
		pointer_end += size_value_byte;
		
		char header_buffer[size_value_byte];
		std::copy(inflate_data + pointer_str, inflate_data + pointer_end, header_buffer);
		T header_tmp;
		std::memcpy( &header_tmp, header_buffer, size_value_byte );
		outData.push_back(header_tmp);
	}
	free(inflate_data);
	
	
	
}

















void MASCH_Mesh_Load::OpenFoam(string folder, MASCH_Mesh &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute load OpenFoam files ... ";
	
	// if(rank == 0)
	{
		
		string gridFolderName = folder;
		string pointsName = "points";
		string facesName = "faces";
		string ownerName = "owner";
		string neighbourName = "neighbour";
		string boundaryName = "boundary";
		
		ifstream inputFile;
		string openFileName;
		
		// points 읽기
		openFileName = gridFolderName + "/" + pointsName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		string nextToken;
		bool startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					vector<double> xyz(3,0.0);
					nextToken.erase(nextToken.find("("),1); 
					nextToken.erase(nextToken.find(")"),1); 
					stringstream sstream(nextToken);
					string word;
					char del = ' ';
					int num=0;
					while (getline(sstream, word, del)){
						xyz[num] = stold(word);
						++num;
					}
					mesh.addPoint();
					mesh.points.back().x = xyz[0];
					mesh.points.back().y = xyz[1];
					mesh.points.back().z = xyz[2];
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		

		// faces 읽기
		openFileName = gridFolderName + "/" + facesName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		startInput=false;
		bool continueInput=false;
		string saveToken;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" && !continueInput ){
					break;
				}
				else{
					if(nextToken.size()==1) continue;
					
					if(nextToken.find(")") == string::npos){
						saveToken.append(" ");
						rtrim(nextToken);
						saveToken.append(nextToken);
						continueInput = true;
						continue;
					}
					saveToken.append(" ");
					rtrim(nextToken);
					saveToken.append(nextToken);
					
					saveToken.replace(saveToken.find("("),1," ");
					saveToken.replace(saveToken.find(")"),1," ");
					istringstream iss(saveToken);
					int tempint;
					iss >> tempint;
					
					mesh.addFace();
					
					while(iss >> tempint){
						mesh.faces.back().ipoints.push_back(
						static_cast<int>(tempint));
					}
					saveToken.clear();
					continueInput = false;
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		// owner
		openFileName = gridFolderName + "/" + ownerName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		int temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					while(iss >> tempint){
						mesh.faces[temp_num].iL = static_cast<int>(tempint);
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		
		// neighbour
		openFileName = gridFolderName + "/" + neighbourName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					while(iss >> tempint){
						if(tempint<0) break;
						mesh.faces[temp_num].iR = static_cast<int>(tempint);
						// mesh.faces[temp_num].thereR = true;
						mesh.faces[temp_num].setType(MASCH_Face_Types::INTERNAL);
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		inputFile.close();
		
		
		
		// boundary
		openFileName = gridFolderName + "/" + boundaryName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		vector<string> boundary_name;
		vector<char> boundary_type;
		vector<int> boundary_nFaces;
		vector<int> boundary_startFace;
		vector<int> boundary_myProcNo;
		vector<int> boundary_neighbProcNo;
		
		string backToken;
		startInput=false;
		vector<string> setToken;
		while(getline(inputFile, nextToken)){
			string asignToken;
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				setToken.push_back(nextToken.c_str());
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}
			backToken = nextToken;
		}
		
		string names;
		vector<string> setToken2;
		startInput=false;
		for (auto item : setToken) {
			string asignToken;
			if(startInput){
				if( item.find("}") != string::npos ){
					trim(names);
					
					boundary_name.push_back(names);
					int l=0;
					for (auto item2 : setToken2) {
						if( item2.find("nFaces") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_nFaces.push_back(temptempint);
						}
						if( item2.find("startFace") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_startFace.push_back(temptempint);
						}
						if( item2.find("myProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_myProcNo.push_back(temptempint);
							++l;
						}
						if( item2.find("neighbProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_neighbProcNo.push_back(temptempint);
							++l;
						}
					}
					if(l==0) {
						boundary_myProcNo.push_back(-1);
						boundary_neighbProcNo.push_back(-1);
					}
					
					startInput=false;
					setToken2.clear();
				}
				setToken2.push_back(item.c_str());
			}
			else{ 
				if( item.find("{") != string::npos ){
					names.clear();
					names = backToken;
					startInput=true;
				}
			}
			backToken = item;
		}
		
		
		inputFile.close();
		mesh.boundaries.clear();
		int nbcs=boundary_name.size();
		for (int i=0; i<boundary_name.size(); ++i) {
			mesh.addBoundary();
		}
		for (int i=0; i<mesh.boundaries.size(); ++i) {
			mesh.boundaries[i].name = trim(boundary_name[i]);
			mesh.boundaries[i].nFaces = static_cast<int>(boundary_nFaces[i]);
			mesh.boundaries[i].startFace = static_cast<int>(boundary_startFace[i]);
			mesh.boundaries[i].myProcNo = static_cast<int>(boundary_myProcNo[i]);
			mesh.boundaries[i].rightProcNo = static_cast<int>(boundary_neighbProcNo[i]);
			// mesh.boundaries[i].thereR = false;
			mesh.boundaries[i].setType(MASCH_Face_Types::BOUNDARY);
		}
	}
	
	if(rank!=0){
		mesh.points.clear();
		mesh.faces.clear();
		mesh.cells.clear();
		for (int i=0; i<mesh.boundaries.size(); ++i) {
			mesh.boundaries[i].nFaces = 0;
			mesh.boundaries[i].startFace = 0;
		}
	}
	
	
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	

	mesh.check();
	mesh.setFaceTypes();
	mesh.buildCells();
	mesh.connectFacetoPointsCells();
	mesh.connectCelltoFaces();
	mesh.connectCelltoPoints();
	mesh.setCountsProcFaces();
	mesh.setDisplsProcFaces();
	mesh.cellsGlobal();
	mesh.informations();
	
	
}

