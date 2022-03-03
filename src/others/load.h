#pragma once
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <mpi.h>
#include <cstring>
#include <zlib.h>
#include <dlfcn.h>

using namespace std;

#include "./mesh.h"  
#include "./controls.h"  
#include "./solvers.h" 
#include "./log.h"  
#include "./variables.h"  

#define MASCH_FFL __FILE__,__FUNCTION__,__LINE__

class MASCH_Control;
class MASCH_Mesh;


class MASCH_Mesh_Load {
private:

public:
	template<typename T>
	void loadDatasAtVTU(
		ifstream& inputFile, string dataName, vector<T>& outData);
		
	template<typename T>
	void readBinary(
		string& word, vector<T>& outData);
		
	template<typename T>
	void readCompress(
		string& word, vector<T>& outData);
		
	int Base64decode_len(const char *bufcoded);
	int Base64decode(char *bufplain, const char *bufcoded);
		
	string &ltrim(std::string &s);
	string &rtrim(std::string &s);
	string &trim(std::string &s);
	
	void vtu(string folder, MASCH_Control &controls, MASCH_Mesh &mesh);
	void OpenFoam(string folder, MASCH_Mesh &in);
	void vtuPrimitive(
		string folderName, int rank,
		MASCH_Control& controls, 
		MASCH_Mesh& mesh,
		MASCH_Variables& var);
	
	bool boolCompress;
};


class MASCH_Load {
public:
	string fileName;
	ifstream inputFile;

	void cuttingCommets(string& inp) {
		inp = erase(inp, "//");
	};
	string erase(string& inp, string er) {
		if (inp.find(er.c_str()) != std::string::npos) {
			inp.erase(inp.find(er.c_str()));
		}
		trim(inp);
		return inp;
	};

	vector<string> strtok(string inp) {
		istringstream ss(inp);
		string word;
		vector<string> out;
		while (getline(ss, word, ' ')) {
			out.push_back(word);
		}
		return out;
	};

	int extract(vector<string>& saveToken, map<string,string>& output, string front="", int str = 0) {

		bool boolSub = false;
		string saveSubName;
		for (int i = str; i < saveToken.size(); ++i) {
			if (saveToken[i].find('}') != std::string::npos) {
				return i;
			}

			if (saveToken[i].find('{') != std::string::npos) {
				string front2 = front + (erase(saveToken[i], "{") + ".");
				i = extract(saveToken, output, front2, i+1);
			}
			else{
				if (saveToken[i].find('=') != std::string::npos) {
					vector<string> tmp; 
					erase(saveToken[i], ";");
					tmp = strtok(saveToken[i]);
					string tmp_value = "";
					for (int j = 2; j < tmp.size(); ++j) {
						tmp_value += (tmp[j] + " ");
					}
					output.insert(make_pair(front+tmp[0], trim(tmp_value)));
					//cout << front << tmp[0] << " : " << tmp[2] << endl;
				}
			}
		}
		return saveToken.size();

	};


	vector<string> extractVector(string inp){
		istringstream ss(inp);
		string word;
		vector<string> out;
		while (getline(ss, word, ' ')) {
			out.push_back(word);
		}

		vector<string> out2;
		vector<int> level;
		recursiveExtractVector(out, out2, level);
		
		int maxLevel = *max_element(level.begin(), level.end());
		// if(maxLevel==1){
			
		// }
		// if(maxLevel==2){
			
		// }
		
		
		return out2;

	}


	int recursiveExtractVector(vector<string>& inp, vector<string>& out, vector<int>& level, int lv=0, int pos=0) {
		vector<string> sub;
		for (int i = pos; i < inp.size(); ++i) {

			if (inp[i].find("(") != string::npos && inp[i].rfind(")") != string::npos) {
				inp[i].insert(inp[i].find("("), " ");
				inp[i].erase(inp[i].find("("), 1);
				inp[i].insert(inp[i].rfind(")"), " ");
				inp[i].erase(inp[i].rfind(")"), 1);
				trim(inp[i]);
				if (inp[i].empty()) continue;
				out.push_back(inp[i]);
				inp[i].erase(inp[i].find(inp[i]), 1);
				trim(inp[i]);
				level.push_back(lv);
				
			}
			else if (inp[i].find("(") != string::npos) {
				inp[i].insert(inp[i].find("("), " ");
				inp[i].erase(inp[i].find("("), 1);
				trim(inp[i]);
				istringstream ss(inp[i]);
				string word;
				vector<string> subOut;
				while (getline(ss, word, ' ')) {
					subOut.push_back(trim(word));
				}
				if (subOut.size() == 1) {
					out.push_back(subOut[0]);
					level.push_back(lv+1);
					inp[i].erase(inp[i].find(subOut[0]), 1);
					trim(inp[i]);
				}
				else {
					out.push_back(subOut[0]);
					level.push_back(lv);
					inp[i].erase(inp[i].find(subOut[0]), 1);
					trim(inp[i]);
					out.push_back(subOut[1]);
					level.push_back(lv+1);
					inp[i].erase(inp[i].find(subOut[1]), 1);
					trim(inp[i]);
				}

				i = recursiveExtractVector(inp, out, level, lv+1, i+1);
			}
			else if (inp[i].find(")") != string::npos) {
				inp[i].insert(inp[i].find(")"), " ");
				inp[i].erase(inp[i].find(")"), 1);
				trim(inp[i]);
				if (inp[i].empty()) return i;
				istringstream ss(inp[i]);
				string word;
				vector<string> subOut;
				while (getline(ss, word, ' ')) {
					subOut.push_back(trim(word));
				}
				out.push_back(subOut[0]);
				level.push_back(lv);
				inp[i].erase(inp[i].find(subOut[0]), 1);
				trim(inp[i]);
				return i;
			}
			else {
				if (inp[i].empty()) continue;
				out.push_back(inp[i]);
				level.push_back(lv);
				inp[i].erase(inp[i].find(inp[i]), 1);
				trim(inp[i]);

			}
		}
		return 0;

	};

	
	
	
	
	
	
	


	void extractFile(string filename, map<string, string>& outMaps) {
		inputFile.open(filename);
		if (inputFile.fail()) {
			cerr << "Unable to open file for reading : " << filename << endl;
		}
		string nextToken;
		vector<string> saveToken;
		while (getline(inputFile, nextToken)) {
			cuttingCommets(nextToken);
			saveToken.push_back(nextToken);
		}
		extract(saveToken, outMaps);

		inputFile.close();

	}


	// trim from left 
	std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
	{
		s.erase(0, s.find_first_not_of(t));
		return s;
	}
	// trim from right 
	std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
	{
		s.erase(s.find_last_not_of(t) + 1);
		return s;
	}
	// trim from left & right 
	std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
	{
		return ltrim(rtrim(s, t), t);
	}

	
	void settingFiles(string folderName, MASCH_Control& controls);
	void meshFiles(string folderName, MASCH_Control& controls, MASCH_Mesh& mesh);
	
	// fvm file 로드
	void fvmFiles(
		string folderName, 
		int rank,
		MASCH_Mesh& mesh,
		MASCH_Control& controls, 
		MASCH_Variables& var);
};


