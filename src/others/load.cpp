
#include "./load.h" 

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
		
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
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
		if(rank==0) cout << "zlib error !!!!!!!!!!!!" << endl;
		if(rank==0) cout << error << endl;
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

