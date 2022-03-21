
#include "./save.h"

void MASCH_Mesh_Save::mkdirs(char *dir_path){
	char buff[1024];
	char *p_dir = buff;
	
	strcpy(buff, dir_path);
	buff[1024-1] = '\0';
	
	while(*p_dir){
		if('/'==*p_dir){
			*p_dir = '\0';
			mkdir(buff,0777);
			*p_dir = '/';
		}
		p_dir ++;
	}
	
}


template void MASCH_Mesh_Save::writeAscii<int>(ofstream& out, vector<int>& vecInp);
template void MASCH_Mesh_Save::writeAscii<double>(ofstream& out, vector<double>& vecInp);
template<typename T>
void MASCH_Mesh_Save::writeAscii(ofstream& out, vector<T>& vecInp)
{
	for(auto& inp : vecInp){
		out << inp << " ";
	}
	out << endl;
}


/* aaaack but it's fast and const should make it shared text page. */
static const char basis_64[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


int MASCH_Mesh_Save::Base64encode_len(int len)
{
    return ((len + 2) / 3 * 4) + 1;
}

int MASCH_Mesh_Save::Base64encode(char *encoded, const char *string, int len)
{
    int i;
    char *p;

    p = encoded;
    for (i = 0; i < len - 2; i += 3) {
    *p++ = basis_64[(string[i] >> 2) & 0x3F];
    *p++ = basis_64[((string[i] & 0x3) << 4) |
                    ((int) (string[i + 1] & 0xF0) >> 4)];
    *p++ = basis_64[((string[i + 1] & 0xF) << 2) |
                    ((int) (string[i + 2] & 0xC0) >> 6)];
    *p++ = basis_64[string[i + 2] & 0x3F];
    }
    if (i < len) {
    *p++ = basis_64[(string[i] >> 2) & 0x3F];
    if (i == (len - 1)) {
        *p++ = basis_64[((string[i] & 0x3) << 4)];
        *p++ = '=';
    }
    else {
        *p++ = basis_64[((string[i] & 0x3) << 4) |
                        ((int) (string[i + 1] & 0xF0) >> 4)];
        *p++ = basis_64[((string[i + 1] & 0xF) << 2)];
    }
    *p++ = '=';
    }

    *p++ = '\0';
    return p - encoded;
}
	
	
	
template void MASCH_Mesh_Save::writeBinary<int>(ofstream& out, vector<int>& vecInp);
template void MASCH_Mesh_Save::writeBinary<double>(ofstream& out, vector<double>& vecInp);
template<typename T>
void MASCH_Mesh_Save::writeBinary(ofstream& out, vector<T>& vecInp)
{
	int datasize = 8;
	int dataByteSize = sizeof(T);
	
	// cout << dataByteSize << endl;
	
	int data_length = datasize + dataByteSize * vecInp.size();
	int encoded_data_length = Base64encode_len(data_length);
	char* base64_string = (char*)malloc(encoded_data_length);
	char* data = (char*)malloc(data_length);
	{
		long long byte_size = datasize * vecInp.size();
		char* data1 = (char*)&byte_size;
		std::copy(data1, data1 + datasize, data);
	}
	int numm = 0;
	for(auto& inp : vecInp) {
		char* data1 = (char*)&inp;
		std::copy(data1, data1 + dataByteSize, data + datasize + numm*dataByteSize);
		++numm;
	}
	Base64encode(base64_string, data, data_length);
	out << base64_string << endl;
	free(data);
	free(base64_string);
	
}



template void MASCH_Mesh_Save::writeDatasAtVTU<int>(
	MASCH_Control &controls, ofstream& out, vector<int>& vecInp);
template void MASCH_Mesh_Save::writeDatasAtVTU<double>(
	MASCH_Control &controls, ofstream& out, vector<double>& vecInp);
template<typename T>
void MASCH_Mesh_Save::writeDatasAtVTU(
	MASCH_Control &controls, ofstream& out, vector<T>& vecInp){
	
	if(controls.saveFormat == "ascii"){
		writeAscii(out, vecInp);
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			writeBinary(out, vecInp);
		}
		else{
			writeCompress(out, vecInp, controls.saveCompression);
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
}



template void MASCH_Mesh_Save::writeCompress<int>(ofstream& out, vector<int>& vecInp, int compressSize);
template void MASCH_Mesh_Save::writeCompress<double>(ofstream& out, vector<double>& vecInp, int compressSize);
template<typename T>
void MASCH_Mesh_Save::writeCompress(ofstream& out, vector<T>& vecInp, int compressSize)
{
	int value_data_size = vecInp.size();
	
	int headerByteSize = 8;
	int dataByteSize = sizeof(T);
	
	
	// cout << dataByteSize << endl;
	
	int header_size = headerByteSize*4;
	int data_size = dataByteSize*value_data_size;
	long long header[4] = {0};
	header[0] = 1;
	header[1] = data_size;
	header[2] = data_size;
	header[3] = 2*data_size;
	Bytef *deflate_data = (Bytef*)malloc(header[3]);
	uLong deflate_size = header[3];
	Bytef* raw_data = (Bytef*)malloc(data_size);
	uLong raw_size = data_size;
	
	int numm = 0;
	for(auto& inp : vecInp){
		char* header_char0 = (char*)&inp;
		std::copy(header_char0, header_char0 + dataByteSize, raw_data + numm*dataByteSize);
		++numm;
	}
	if(compress2(deflate_data, &deflate_size, raw_data, raw_size, compressSize) != Z_OK){
		cout << "zlib error !!!!!!!!!!!!" << endl;
	}
	// cout << endl;
	// cout << deflate_size <<  endl;
	header[3] = deflate_size;
	int header_encoded_length = Base64encode_len(header_size);
	char* base64_string_header = (char*)malloc(header_encoded_length);
	char* header_data = (char*)malloc(header_size);
	
	char* header_char1 = (char*)&header[0];
	char* header_char2 = (char*)&header[1];
	char* header_char3 = (char*)&header[2];
	char* header_char4 = (char*)&header[3];
	std::copy(header_char1, header_char1 + headerByteSize, header_data+headerByteSize*0);
	std::copy(header_char2, header_char2 + headerByteSize, header_data+headerByteSize*1);
	std::copy(header_char3, header_char3 + headerByteSize, header_data+headerByteSize*2);
	std::copy(header_char4, header_char4 + headerByteSize, header_data+headerByteSize*3);
	Base64encode(base64_string_header, header_data, header_size);
	
	int data_encoded_length = Base64encode_len(deflate_size);
	char* base64_string_data = (char*)malloc(data_encoded_length);
	Base64encode(base64_string_data, (char*)deflate_data, deflate_size);
	
	// cout << base64_string_header << " " << header_encoded_length << endl;
	out << base64_string_header << base64_string_data << endl; 
	
	free(base64_string_header);
	free(base64_string_data);
	free(deflate_data);
	free(header_data);
	free(raw_data); 
}





void MASCH_Mesh_Save::vtu(
	string folder, 
	MASCH_Mesh &mesh, 
	MASCH_Control &controls){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_Utility_Math math;
	
	// int compressSize = controls.saveCompression;
	// // int compressSize = 6;
	
	// char folder_name[1000];
	// strcpy(folder_name, folder.c_str());
	// mkdirs(folder_name);

	// ofstream outputFile;
	// string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute save file (" << folder_name << "plot...) ... ";
	// }
	
	// outputFile.open(filenamePlot);
	
	// outputFile.precision( controls.writePrecision );
	
	// if(outputFile.fail()){
		// cerr << "Unable to write file for writing." << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	// // string out_line;
	// outputFile << "<?xml version=\"1.0\"?>" << endl;

	// string saveFormat = controls.saveFormat;
	// if(controls.saveFormat == "ascii"){
		// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	// }
	// else if(controls.saveFormat == "binary"){
		// if(controls.saveCompression==0){
			// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		// }
		// else{
			// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		// }
	// }
	// else{
		// cout << endl;
		// cout << endl;
		// cout << "| warning : not defined saveFormat at controlDic file" << endl;
		// cout << endl;
		// cout << endl;
	// }
	
	// outputFile << "  <UnstructuredGrid>" << endl;
	

	// // Field data
	// outputFile << "    <FieldData>" << endl;
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// values.push_back(controls.time);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// outputFile << "    </FieldData>" << endl;
	
	// outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// // Points data
	// outputFile << "    <PointData>" << endl;
	// {
		// outputFile << "     <DataArray type=\"Int32\" Name=\"pointLevels\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& point : mesh.points) values.push_back(point.level);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// outputFile << "    </PointData>" << endl;
	
	
	
	// // Cells data
	// outputFile << "    <CellData>" << endl;
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.P]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells){
			// values.push_back(cell.var[controls.U]); 
			// values.push_back(cell.var[controls.V]); 
			// values.push_back(cell.var[controls.W]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.T]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << controls.name[controls.VF[0]] << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.VF[0]]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << controls.name[controls.MF[0]] << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.MF[0]]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // cell levels
	// {
		// outputFile << "     <DataArray type=\"Int32\" Name=\"cellLevels\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.level);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // cell groups
	// {
		// outputFile << "     <DataArray type=\"Int32\" Name=\"cellGroups\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.group);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // 추가적인 데이터 저장
	// // mesh data
	// if(controls.saveMeshData["nonOrthogonality"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "nonOrthogonality" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> nonOrthogonality;
		// mesh.calcNonOrthogonality(nonOrthogonality);
		// writeDatasAtVTU(controls, outputFile, nonOrthogonality);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// if(controls.saveMeshData["uniformity"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "uniformity" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> uniformity;
		// mesh.calcUniformity(uniformity);
		// writeDatasAtVTU(controls, outputFile, uniformity);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveMeshData["skewness"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "skewness" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> skewness;
		// mesh.calcSkewness(skewness);
		// writeDatasAtVTU(controls, outputFile, skewness);
		// outputFile << "     </DataArray>" << endl;
	// }
	// // cout << "B1" << endl;
	// // gradient
	// if(controls.saveGradientData["pressure"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-pressure" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.P, controls.fP, gradient);
		// math.calcGaussGreen(mesh, controls.P, controls.fP, gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveGradientData["x-velocity"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-x-velocity" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.U, controls.fU, gradient);
		// math.calcGaussGreen(mesh, controls.U, controls.fU, gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveGradientData["y-velocity"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-y-velocity" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.V, controls.fV, gradient);
		// math.calcGaussGreen(mesh, controls.V, controls.fV, gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveGradientData["z-velocity"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-z-velocity" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.W, controls.fW, gradient);
		// math.calcGaussGreen(mesh, controls.W, controls.fW, gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// // // if(controls.saveGradientData["velocityMagnitude"]){
		// // // outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-pressure" << "\" format=\"" << saveFormat << "\">" << endl;
		// // // vector<double> gradient(mesh.cells.size(),0.0);
		// // // math.calcGaussGreen(mesh, controls.P, controls.fP, gradient);
		// // // math.calcGaussGreen(mesh, controls.P, controls.fP, gradient);
		// // // writeDatasAtVTU(controls, outputFile, gradient);
		// // // outputFile << "     </DataArray>" << endl;
	// // // }
	// if(controls.saveGradientData["temperature"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-temperature" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.T, controls.fT, gradient);
		// math.calcGaussGreen(mesh, controls.T, controls.fT, gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveGradientData[controls.name[controls.VF[0]]]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-" << controls.name[controls.VF[0]] << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradient);
		// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradient);
	// // cout << "B2" << endl;
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveGradientData[controls.name[controls.MF[0]]]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gradient-" << controls.name[controls.MF[0]] << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<vector<double>> gradient(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.MF[0], controls.fMF[0], gradient);
		// math.calcGaussGreen(mesh, controls.MF[0], controls.fMF[0], gradient);
		// vector<double> values;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// values.push_back(gradient[i][0]); 
			// values.push_back(gradient[i][1]); 
			// values.push_back(gradient[i][2]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// // thermodynamic
	// if(controls.saveThermodynamicData["density"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "density" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.Rho]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveThermodynamicData["speedOfSound"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "speedOfSound" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.C]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveThermodynamicData["enthalpy"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "enthalpy" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(
			// cell.var[controls.Ht]-0.5*(
			// cell.var[controls.U]*cell.var[controls.U]+
			// cell.var[controls.V]*cell.var[controls.V]+
			// cell.var[controls.W]*cell.var[controls.W]));
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveThermodynamicData["totalEnthalpy"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "totalEnthalpy" << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.Ht]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// // body-force
	// if(controls.saveBodyForceData["gravity"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "gravity" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells){
			// values.push_back(cell.var[controls.sourceGravity[0]]); 
			// values.push_back(cell.var[controls.sourceGravity[1]]);  
			// values.push_back(cell.var[controls.sourceGravity[2]]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// if(controls.saveBodyForceData["surfaceTension"]){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << "surface-tension" << "\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells){
			// values.push_back(cell.var[controls.sourceSurfaceTension[0]]); 
			// values.push_back(cell.var[controls.sourceSurfaceTension[1]]);  
			// values.push_back(cell.var[controls.sourceSurfaceTension[2]]); 
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // UDV
	// for(int i=0; i<controls.UDV.size(); ++i){
		// outputFile << "     <DataArray type=\"Float64\" Name=\"" << controls.name[controls.UDV[i]] << "\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& cell : mesh.cells) values.push_back(cell.var[controls.UDV[i]]);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	
	// // // level-set
	// // {
		// // outputFile << "     <DataArray type=\"Float64\" Name=\"" << "level-set" << "\" format=\"" << saveFormat << "\">" << endl;
		// // vector<double> values;
		// // for(auto& cell : mesh.cells) values.push_back(cell.var[controls.LS]);
		// // writeDatasAtVTU(controls, outputFile, values);
		// // outputFile << "     </DataArray>" << endl;
	// // }
	
	
	
	
	
	// outputFile << "    </CellData>" << endl;
	
	
	// // Points
	// outputFile << "    <Points>" << endl;
	// {
		// outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		// vector<double> values;
		// for(auto& point : mesh.points) {
			// values.push_back(point.x);
			// values.push_back(point.y);
			// values.push_back(point.z);
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	// outputFile << "   </Points>" << endl;
	
	// // cells
	// outputFile << "   <Cells>" << endl; 
	// // connectivity (cell's points)
	// {
		// outputFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& cell : mesh.cells){
			// for(auto i : cell.points){
				// values.push_back(i);
			// }
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // offsets (cell's points offset)
	// {
		// outputFile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// int cellFaceOffset = 0;
		// for(auto& cell : mesh.cells){
			// cellFaceOffset += cell.points.size();
			// values.push_back(cellFaceOffset);
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // types (cell's type, 42 = polyhedron)
	// {
		// outputFile << "    <DataArray type=\"Int32\" Name=\"types\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values(mesh.cells.size(),42);
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // faces (cell's faces number, each face's point number, cell's faces's points)
	// {
		// outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faces\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& cell : mesh.cells){
			// values.push_back(cell.faces.size());
			// for(auto& i : cell.faces){
				// values.push_back(mesh.faces[i].points.size());
				// for(auto& j : mesh.faces[i].points){
					// values.push_back(j);
				// }
			// }
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// // faceoffsets (cell's face offset)
	// {
		// outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faceoffsets\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// int cellFacePointOffset = 0;
		// for(auto& cell : mesh.cells){
			// int numbering = 1 + cell.faces.size();
			// for(auto& i : cell.faces){
				// numbering += mesh.faces[i].points.size();
			// }
			// cellFacePointOffset += numbering;
			// values.push_back(cellFacePointOffset);
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << "     </DataArray>" << endl;
	// }
	
	// outputFile << "   </Cells>" << endl;
	
	// outputFile << "  </Piece>" << endl;
	// outputFile << " </UnstructuredGrid>" << endl;
	

	// // additional informations
	// {
		// outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& face : mesh.faces){
			// values.push_back(face.owner);
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << " </DataArray>" << endl;
	// }
	// {
		// outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		// for(auto& face : mesh.faces){
			// values.push_back(face.neighbour);
		// }
		// writeDatasAtVTU(controls, outputFile, values);
		// outputFile << " </DataArray>" << endl;
	// }
	
	

	// // boundary informations
	// {
		// outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// // cout << boundary.name << endl;
			// // trim;
			// string bcName = boundary.name;
			
			// bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			// bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// // outputFile << boundary.name << " ";
			// outputFile << bcName << " ";
		// }
		// outputFile << endl;
		// outputFile << " </DataArray>" << endl;
		
		// outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// outputFile << boundary.startFace << " ";
		// }
		// outputFile << endl;
		// outputFile << " </DataArray>" << endl;
		
		// outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// outputFile << boundary.nFaces << " ";
		// }
		// outputFile << endl;
		// outputFile << " </DataArray>" << endl;
		
		// outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// outputFile << boundary.neighbProcNo << " ";
		// }
		// outputFile << endl;
		// outputFile << " </DataArray>" << endl;
		
	// }
	
	// outputFile << "</VTKFile>" << endl;
	
	// outputFile.close();
	
	
	
	// // ==========================================
	// // pvtu file
	// if(rank==0){
		// string filenamePvtu = "./save/plot.";
		// string stime = folder;
		// stime.erase(stime.find("./save/"),7);
		// stime.erase(stime.find("/"),1);
		// filenamePvtu += stime;
		// filenamePvtu += ".pvtu";
		
		// outputFile.open(filenamePvtu);
		// if(outputFile.fail()){
			// cerr << "Unable to write file for writing." << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		// // string out_line;
		// outputFile << "<?xml version=\"1.0\"?>" << endl;
		// outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		// outputFile << "  <PUnstructuredGrid>" << endl;

		// outputFile << "   <PFieldData>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"TimeValue\"/>" << endl;
		// outputFile << "   </PFieldData>" << endl;
		
		// outputFile << "   <PPoints>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		// outputFile << "   </PPoints>" << endl;
		// for(int ip=0; ip<size; ++ip){
			// string filenamevtus = "./" + stime;
			// filenamevtus += "/plot.";
			// filenamevtus += to_string(ip);
			// filenamevtus += ".vtu";
			// outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		// }
		// outputFile << "   <PPointData>" << endl;

		// outputFile << "    <PDataArray type=\"Int32\" Name=\"pointLevels\"/>" << endl;
		
		// outputFile << "   </PPointData>" << endl;
		// outputFile << "   <PCellData>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		// // outputFile << "    <PDataArray type=\"Float64\" Name=\"density\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"" << controls.name[controls.VF[0]] << "\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"" << controls.name[controls.MF[0]] << "\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Int32\" Name=\"cellLevels\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Int32\" Name=\"cellGroups\"/>" << endl;
		// // 추가적인 데이터 저장
		// // mesh data
		// if(controls.saveMeshData["nonOrthogonality"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"nonOrthogonality\"/>" << endl;
		// }
		// if(controls.saveMeshData["uniformity"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"uniformity\"/>" << endl;
		// }
		// if(controls.saveMeshData["skewness"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"skewness\"/>" << endl;
		// }
		// // gradient
		// if(controls.saveGradientData["pressure"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-pressure\"/ NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveGradientData["x-velocity"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-x-velocity\"/ NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveGradientData["y-velocity"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-y-velocity\"/ NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveGradientData["z-velocity"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-z-velocity\"/ NumberOfComponents=\"3\"/>" << endl;
		// }
		// // if(controls.saveGradientData["velocityMagnitude"]){
			// // outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient\"/>" << endl;
		// // }
		// if(controls.saveGradientData["temperature"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-temperature\"/ NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveGradientData[controls.name[controls.VF[0]]]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-" << controls.name[controls.VF[0]] << "\" NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveGradientData[controls.name[controls.MF[0]]]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gradient-" << controls.name[controls.MF[0]] << "\" NumberOfComponents=\"3\"/>" << endl;
		// }
		// // thermodynamic
		// if(controls.saveThermodynamicData["density"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"density\"/>" << endl;
		// }
		// if(controls.saveThermodynamicData["speedOfSound"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"speedOfSound\"/>" << endl;
		// }
		// if(controls.saveThermodynamicData["enthalpy"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"enthalpy\"/>" << endl;
		// }
		// if(controls.saveThermodynamicData["totalEnthalpy"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"totalEnthalpy\"/>" << endl;
		// }
		// // body-force
		// if(controls.saveBodyForceData["gravity"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"gravity\" NumberOfComponents=\"3\"/>" << endl;
		// }
		// if(controls.saveBodyForceData["surfaceTension"]){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"surface-tension\" NumberOfComponents=\"3\"/>" << endl;
		// }
		// // UDV
		// for(int i=0; i<controls.UDV.size(); ++i){
			// outputFile << "    <PDataArray type=\"Float64\" Name=\"" << controls.name[controls.UDV[i]] << "\"/>" << endl;
		// }
		// // // level-set
		// // {
			// // outputFile << "    <PDataArray type=\"Float64\" Name=\"" << "level-set" << "\"/>" << endl;
		// // }
		
		
		// outputFile << "   </PCellData>" << endl;
		// outputFile << "  </PUnstructuredGrid>" << endl;
		// outputFile << "</VTKFile>" << endl;
		
		
		// outputFile.close();
		
	// }
	
	
	
	
	
	// // ==========================================
	// // PVD file
	// if(rank==0){
		// string filenamePvtu = "./save/plot.pvd";
		// string stime = folder;
		// stime.erase(stime.find("./save/"),7);
		// stime.erase(stime.find("/"),1);
		
		// ifstream inputFile;
		// inputFile.open(filenamePvtu);
		// if(inputFile && stod(stime)-controls.timeStep != 0.0){

			// string nextToken;
			// int lineStart = 0;
			// while(getline(inputFile, nextToken)){
				// if( nextToken.find("</Collection>") != string::npos ){
					// break;
				// }
				// ++lineStart;
			// }
			// inputFile.close();
			
			// outputFile.open(filenamePvtu, ios::in);
			// outputFile.seekp(-26, ios::end);
			
			// // outputFile << saveLines;
			// outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			// outputFile << "  </Collection>" << endl;
			// outputFile << "</VTKFile>";
			// outputFile.close();
			
		// }
		// else{
			// inputFile.close();
			
			// outputFile.open(filenamePvtu);
			
			// // string out_line;
			// outputFile << "<?xml version=\"1.0\"?>" << endl;
			// outputFile << " <VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
			// outputFile << "  <Collection>" << endl;
			// outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			// outputFile << "  </Collection>" << endl;
			// outputFile << "</VTKFile>";
			
			
			// outputFile.close();
			
		// }
		
		
	// }
	
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
}







void MASCH_Mesh_Save::vtu(string folder, int myRank, MASCH_Mesh &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	// 폴더 만들기
    // auto ret = filesystem::create_directories(folder);
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(myRank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << filenamePlot << ") ... ";
	}
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	outputFile << "    </PointData>" << endl;
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;  
	outputFile.precision(20);  
	// for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.ipoints){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.ipoints.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	{
		// faces (cell's faces number, each face's point number, cell's faces's points)
		outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
		
		for(auto& cell : mesh.cells){
			outputFile << cell.ifaces.size() << endl;
			for(auto& i : cell.ifaces){
				outputFile << mesh.faces[i].ipoints.size() << " ";
				int tmp_space = 0;
				for(auto& j : mesh.faces[i].ipoints){
					if(tmp_space!=mesh.faces[i].ipoints.size()-1){
						outputFile << j << " ";
					}
					else{
						outputFile << j;
					}
					++tmp_space;
				}
				outputFile << endl;
			}
		}
		// for(auto& cell : mesh.cells){
			// outputFile << cell.ifaces.size();
			// for(auto& i : cell.ifaces){
				// outputFile << " " << mesh.faces[i].ipoints.size();
				// for(auto& j : mesh.faces[i].ipoints){
					// outputFile << " " << j;
				// }
				// // outputFile << endl;
			// }
			// outputFile << " ";
		// }
		// outputFile << endl;
		
		outputFile << "    </DataArray>" << endl;
	
		// outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faces\" format=\"" << "ascii" << "\">" << endl;
		// vector<int> values;
		// for(auto& cell : mesh.cells){
			// values.push_back(cell.ifaces.size());
			// for(auto& i : cell.ifaces){
				// values.push_back(mesh.faces[i].ipoints.size());
				// for(auto& j : mesh.faces[i].ipoints){
					// values.push_back(j);
				// }
			// }
		// }
		// for(auto& item : values){
			// outputFile << item << endl;
		// }
		// outputFile << "     </DataArray>" << endl;
	
	}
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.ifaces.size();
		for(auto& i : cell.ifaces){
			numbering += mesh.faces[i].ipoints.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	
	

	// additional informations
	{
		string saveFormat = "ascii";
		outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		for(auto& face : mesh.faces){
			// values.push_back(face.owner);
			outputFile << face.iL << " ";
		}
		// writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	{
		string saveFormat = "ascii";
		outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		for(auto& face : mesh.faces){
			// values.push_back(face.neighbour);
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				outputFile << face.iR << " ";
			}
		}
		// writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	
	

	// boundary informations
	{
		outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// // cout << boundary.name << endl;
			// // trim;
			// string bcName = boundary.name;
			
			// bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			// bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// // outputFile << boundary.name << " ";
			// outputFile << bcName << " ";
		// }
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.name << endl;
			outputFile << boundary.name << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.startFace << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.nFaces << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.rightProcNo << endl;
			outputFile << boundary.rightProcNo << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"connPoints\" format=\"" << "ascii" << "\">" << endl;
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [item, item2] : point.connPoints){
				outputFile << i << " " << item << " " << item2 << " ";
			}
			// if(!point.connPoints.empty()) outputFile << endl;
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		
	}
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}










void MASCH_Mesh_Save::fvmFiles(
string folder, int rank, MASCH_Mesh& mesh, 
MASCH_Control& controls, MASCH_Variables& var){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	
	// 폴더 만들기
    // auto ret = filesystem::create_directories(folder);
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	if(rank==0){
		string printPlot = folder + "plot.proc.vtu";
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << printPlot << ") ... ";
	}
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	string saveFormat = controls.saveFormat;
	
	
	
	
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;

	if(controls.saveFormat == "ascii"){
		outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		}
		else{
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
	
	outputFile << "  <UnstructuredGrid>" << endl;
	





	// Field data
	outputFile << "    <FieldData>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.push_back(var.fields[controls.fieldVar["time"].id]);
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"pointLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.points.size());
		for(auto& point : mesh.points) values.push_back(point.level);
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </PointData>" << endl;
	
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.cells.size());
		for(auto& cell : mesh.cells) values.push_back(cell.level);
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellGroups\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		values.reserve(mesh.cells.size());
		for(auto& cell : mesh.cells) values.push_back(cell.group);
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	
	
	
	
	
	
	
	// 원시변수들 네이밍
	vector<string> scal_cell_save_name;
	for(auto& item : controls.primScalarNames) scal_cell_save_name.push_back(item);
	vector<string> vec_cell_save_sup_name;
	vector<vector<string>> vec_cell_save_name;
	for(auto& item : controls.primVector3Names) {
		vec_cell_save_sup_name.push_back(item);
		vec_cell_save_name.push_back(vector<string>());
		vector<string> sub_names = controls.cellVar[item].sub_name;
		for(int i=0; i<sub_names.size(); ++i){
			vec_cell_save_name.back().push_back(sub_names[i]);
		}
	}
	for(auto& item : controls.primVectorNames) {
		// vec_cell_save_name.push_back(vector<string>());
		vector<string> sub_roles = controls.cellVar[item].sub_role;
		vector<string> sub_names = controls.cellVar[item].sub_name;
		for(int i=0; i<sub_names.size(); ++i){
			if(sub_roles[i]!="primitive") continue;
			string name = item;
			name += "-"+sub_names[i];
			scal_cell_save_name.push_back(name);
		}
	}
	// 추가적 데이터 네이밍
	// for(auto& item : controls.saveFaceValues) scal_cell_save_name.push_back(item);
	for(auto& item : controls.saveCellValues) scal_cell_save_name.push_back(item);
	for(auto& item : controls.saveGradientValues) {
		string xGrad = "x-gradient "; string yGrad = "y-gradient "; string zGrad = "z-gradient ";
		xGrad += item; yGrad += item; zGrad += item;
		vec_cell_save_sup_name.push_back("gradient "+item);
		vec_cell_save_name.push_back(vector<string>());
		vec_cell_save_name.back().push_back(xGrad);
		vec_cell_save_name.back().push_back(yGrad);
		vec_cell_save_name.back().push_back(zGrad);
	}
	
	// 스칼라 형식 데이터 저장
	for(auto& name : scal_cell_save_name)
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
		name << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.reserve(mesh.cells.size());
		int id_phi = controls.getId_cellVar(name);
		auto cellVar = var.cells.data();
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i) {
			auto cellVar_i = cellVar[i];
			values.push_back(cellVar_i[id_phi]);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	// 벡터형식 데이터 저장
	{
		int tmp_iter = 0;
		for(auto& sup_name : vec_cell_save_sup_name)
		{
			outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
			sup_name << "\" NumberOfComponents=\"" << vec_cell_save_name[tmp_iter].size() << 
			"\" format=\"" << saveFormat << "\">" << endl;
			vector<double> values;
			values.reserve(3*mesh.cells.size());
			vector<int> id_phi;
			for(auto& item : vec_cell_save_name[tmp_iter]){
				id_phi.push_back(controls.getId_cellVar(item));
			}
			auto cellVar = var.cells.data();
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i) {
				auto cellVar_i = cellVar[i];
				for(auto& item : id_phi){
					values.push_back(cellVar_i[item]);
				}
			}
			writeDatasAtVTU(controls, outputFile, values);
			outputFile << "     </DataArray>" << endl;
			
			++tmp_iter;
		}
	}
	
	outputFile << "    </CellData>" << endl;
	
	
	
	
	
	
	
	
	
	
	
	
	// Points
	outputFile << "    <Points>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& point : mesh.points) {
			values.push_back(point.x);
			values.push_back(point.y);
			values.push_back(point.z);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "   </Points>" << endl;
	
	
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			for(auto i : cell.ipoints){
				values.push_back(i);
			}
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// offsets (cell's points offset)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFaceOffset = 0;
		for(auto& cell : mesh.cells){
			cellFaceOffset += cell.ipoints.size();
			values.push_back(cellFaceOffset);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// types (cell's type, 42 = polyhedron)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"types\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values(mesh.cells.size(),42);
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faces\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			values.push_back(cell.ifaces.size());
			for(auto& i : cell.ifaces){
				values.push_back(mesh.faces[i].ipoints.size());
				for(auto& j : mesh.faces[i].ipoints){
					values.push_back(j);
				}
			}
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// faceoffsets (cell's face offset)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faceoffsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFacePointOffset = 0;
		for(auto& cell : mesh.cells){
			int numbering = 1 + cell.ifaces.size();
			for(auto& i : cell.ifaces){
				numbering += mesh.faces[i].ipoints.size();
			}
			cellFacePointOffset += numbering;
			values.push_back(cellFacePointOffset);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	// additional informations
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			values.push_back(face.iL);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << " </DataArray>" << endl;
	}
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				values.push_back(face.iR);
			}
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	
	
	
	// boundary informations
	{
		outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.name << endl;
			// trim;
			string bcName = boundary.name;
			
			bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// outputFile << boundary.name << " ";
			outputFile << bcName << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.startFace << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.nFaces << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.rightProcNo << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"connPoints\" format=\"" << "ascii" << "\">" << endl;
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [item, item2] : point.connPoints){
				outputFile << i << " " << item << " " << item2 << " ";
			}
			// if(!point.connPoints.empty()) outputFile << endl;
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		
	}
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	
	
	// ==========================================
	// pvtu file
	if(rank==0){
		string filenamePvtu = "./save/plot.";
		string stime = folder;
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		filenamePvtu += stime;
		filenamePvtu += ".pvtu";
		
		outputFile.open(filenamePvtu);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// string out_line;
		outputFile << "<?xml version=\"1.0\"?>" << endl;
		outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		outputFile << "  <PUnstructuredGrid>" << endl;

		outputFile << "   <PFieldData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"TimeValue\"/>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0, SIZE=MPI::COMM_WORLD.Get_size(); ip<SIZE; ++ip){
			string filenamevtus = "./" + stime;
			filenamevtus += "/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;

		outputFile << "    <PDataArray type=\"Int32\" Name=\"pointLevels\"/>" << endl;
		
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		
		for(auto& name : scal_cell_save_name){
			outputFile << "    <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << endl;
		}
		{
			int tmp_iter=0;
			for(auto& sup_name : vec_cell_save_sup_name){
				outputFile << "    <PDataArray type=\"Float64\" Name=\"" <<
				sup_name << "\" NumberOfComponents=\"" << 3 << "\"/>" << endl;
				// outputFile << "    <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << endl;
				++tmp_iter;
			}
		}
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	
	
	
	// ==========================================
	// PVD file
	if(rank==0){
		string filenamePvtu = "./save/plot.pvd";
		string stime = folder;
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		
		ifstream inputFile;
		inputFile.open(filenamePvtu);
		// if(inputFile && stod(stime)-controls.timeStep != 0.0){
		if(inputFile){

			string nextToken;
			int lineStart = 0;
			while(getline(inputFile, nextToken)){
				if( nextToken.find("</Collection>") != string::npos ){
					break;
				}
				++lineStart;
			}
			inputFile.close();
			
			outputFile.open(filenamePvtu, ios::in);
			outputFile.seekp(-26, ios::end);
			
			// outputFile << saveLines;
			outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			outputFile << "  </Collection>" << endl;
			outputFile << "</VTKFile>";
			outputFile.close();
			
		}
		else{
			inputFile.close();
			
			outputFile.open(filenamePvtu);
			
			// string out_line;
			outputFile << "<?xml version=\"1.0\"?>" << endl;
			outputFile << " <VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
			outputFile << "  <Collection>" << endl;
			outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			outputFile << "  </Collection>" << endl;
			outputFile << "</VTKFile>";
			
			
			outputFile.close();
			
		}
		
		
	}
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
	
	
	
}









void MASCH_Mesh_Save::fvmFiles_boundary(
string folder, int rank, MASCH_Mesh& mesh, 
MASCH_Control& controls, MASCH_Variables& var){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	
	// 폴더 만들기
    // auto ret = filesystem::create_directories(folder);
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save vtp boundary face files ... ";
	}
	
	for(auto& boundary : mesh.boundaries){
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			ofstream outputFile;
			string filenamePlot = folder + "plot_bc_" + boundary.name + "." + to_string(rank) + ".vtp";
			
			
			outputFile.open(filenamePlot);
			if(outputFile.fail()){
				cerr << "Unable to write file for writing." << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			controls.saveFormat = "ascii";
			string saveFormat = controls.saveFormat;
			
			if(controls.saveFormat == "ascii"){
				outputFile << " <VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
			}
			else if(controls.saveFormat == "binary"){
				if(controls.saveCompression==0){
					outputFile << " <VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
				}
				else{
					outputFile << " <VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
				}
			}
			else{
				cout << endl;
				cout << endl;
				cout << "| warning : not defined saveFormat at controlDic file" << endl;
				cout << endl;
				cout << endl;
			}
			
			outputFile << "  <PolyData>" << endl;
			

			int new_nPoints = 0;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			vector<int> point_id(mesh.points.size(),-10);
			vector<double> point_x, point_y, point_z;
			vector<vector<int>> face_new_ipoints;
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				face_new_ipoints.push_back(vector<int>());
				for(auto& ipoint : face.ipoints){
					if(point_id[ipoint]==-10) {
						point_id[ipoint] = new_nPoints++;
						point_x.push_back(mesh.points[ipoint].x);
						point_y.push_back(mesh.points[ipoint].y);
						point_z.push_back(mesh.points[ipoint].z);
						face_new_ipoints.back().push_back(point_id[ipoint]);
					}
					else{
						face_new_ipoints.back().push_back(point_id[ipoint]);
					}
				}
			}


			outputFile << "   <Piece NumberOfPoints=\"" << point_x.size() << "\" NumberOfPolys=\"" << boundary.nFaces << "\">" << endl;
			
			
			// Points
			outputFile << "    <Points>" << endl;
			{
				outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
				vector<double> values;
				for(int i=0; i<point_x.size(); ++i) {
					values.push_back(point_x[i]);
					values.push_back(point_y[i]);
					values.push_back(point_z[i]);
				}
				writeDatasAtVTU(controls, outputFile, values);
				outputFile << "     </DataArray>" << endl;
			}
			outputFile << "   </Points>" << endl;
			
			
			
			// Polys
			outputFile << "   <Polys>" << endl; 
			// connectivity (cell's points)
			{
				outputFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << saveFormat << "\">" << endl;
				vector<int> values;
				for(auto& faces : face_new_ipoints){
					for(auto ipoint : faces){
						values.push_back(ipoint);
					}
				}
				writeDatasAtVTU(controls, outputFile, values);
				outputFile << "     </DataArray>" << endl;
			}
			
			// offsets (cell's points offset)
			{
				outputFile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"" << saveFormat << "\">" << endl;
				vector<int> values;
				int cellFaceOffset = 0;
				for(auto& faces : face_new_ipoints){
					cellFaceOffset += faces.size();
					values.push_back(cellFaceOffset);
				}
				writeDatasAtVTU(controls, outputFile, values);
				outputFile << "     </DataArray>" << endl;
			}
			
			
			outputFile << "   </Polys>" << endl;
			
			
			outputFile << "  </Piece>" << endl;
			outputFile << " </PolyData>" << endl;
			
			
			outputFile << "</VTKFile>" << endl;
			
			outputFile.close();
			
			
			
			// ==========================================
			// pvtu file
			if(rank==0){
				string filenamePvtu = "./save/plot_bc_";
				filenamePvtu += boundary.name;
				filenamePvtu += ".";
				string stime = folder;
				stime.erase(stime.find("./save/"),7);
				stime.erase(stime.find("/"),1);
				filenamePvtu += stime;
				filenamePvtu += ".pvtp";
				
				outputFile.open(filenamePvtu);
				if(outputFile.fail()){
					cerr << "Unable to write file for writing." << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				// string out_line;
				outputFile << "<?xml version=\"1.0\"?>" << endl;
				outputFile << " <VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
				outputFile << "  <PPolyData>" << endl;

				outputFile << "   <PPoints>" << endl;
				outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
				outputFile << "   </PPoints>" << endl;
				for(int ip=0, SIZE=MPI::COMM_WORLD.Get_size(); ip<SIZE; ++ip){
					string filenamevtus = "./" + stime;
					filenamevtus += "/plot_bc_";
					filenamevtus += boundary.name;
					filenamevtus += ".";
					filenamevtus += to_string(ip);
					filenamevtus += ".vtp";
					outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
				}
				outputFile << "  </PPolyData>" << endl;
				outputFile << "</VTKFile>" << endl;
				
				
				outputFile.close();
				
			}
			
			
		}
	}
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
	
}

















void MASCH_Mesh_Save::parcels(
string folder, int rank, MASCH_Mesh& mesh, 
MASCH_Control& controls, MASCH_Variables& var){
	
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save parcel files ... ";
	}

	ofstream outputFile;
	string filenamePlot = folder + "parcels." + to_string(rank) + ".vtu";
	
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	string saveFormat = controls.saveFormat;
	
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	if(controls.saveFormat == "ascii"){
		outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		}
		else{
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
	
	outputFile << "  <UnstructuredGrid>" << endl;

	outputFile << "    <FieldData>" << endl;
	outputFile << "    </FieldData>" << endl;

	outputFile << "   <Piece NumberOfPoints=\"" << mesh.parcels.size() << "\" NumberOfCells=\"" << 0 << "\">" << endl;
	
	outputFile << "    <PointData>" << endl;
	
	
	
	vector<string> scal_save_name;
	scal_save_name.push_back("diameter");
	scal_save_name.push_back("temperature");
	vector<string> vec_save_sup_name;
	vec_save_sup_name.push_back("velocity");
	
	// 스칼라 형식 데이터 저장
	for(auto& name : scal_save_name)
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
		name << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.reserve(mesh.parcels.size());
		int id_phi = controls.getId_parcelVar(name);
		auto parcelVar = var.parcels.data();
		for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
			auto parcelVar_i = parcelVar[i].data();
			values.push_back(parcelVar_i[id_phi]);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	// 벡터형식 데이터 저장
	{
		int tmp_iter = 0;
		for(auto& sup_name : vec_save_sup_name)
		{
			outputFile << "     <DataArray type=\"Float64\" Name=\"" <<
			sup_name << "\" NumberOfComponents=\"" << 3 << 
			"\" format=\"" << saveFormat << "\">" << endl;
			vector<double> values;
			values.reserve(3*mesh.parcels.size());
			vector<int> id_phi;
			id_phi.push_back(controls.getId_parcelVar("x-"+sup_name));
			id_phi.push_back(controls.getId_parcelVar("y-"+sup_name));
			id_phi.push_back(controls.getId_parcelVar("z-"+sup_name));
			auto parcelVar = var.parcels.data();
			for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
				auto parcelVar_i = parcelVar[i].data();
				for(auto& item : id_phi){
					values.push_back(parcelVar_i[item]);
				}
			}
			writeDatasAtVTU(controls, outputFile, values);
			outputFile << "     </DataArray>" << endl;
			
			++tmp_iter;
		}
	}
	
	outputFile << "    </PointData>" << endl;

	
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	
	// Points
	outputFile << "    <Points>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		auto parcelVar = var.parcels.data();
		int id_x = controls.getId_parcelVar("x-location");
		int id_y = controls.getId_parcelVar("y-location");
		int id_z = controls.getId_parcelVar("z-location");
		for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i) {
			auto parcelVar_i = parcelVar[i].data();
			values.push_back(parcelVar_i[id_x]);
			values.push_back(parcelVar_i[id_y]);
			values.push_back(parcelVar_i[id_z]);
		}
		writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "   </Points>" << endl;
	
	
	outputFile << "   <Cells>" << endl; 
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	outputFile << "</VTKFile>" << endl;
	outputFile.close();
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
}
