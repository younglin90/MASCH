#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

#include "../mesh/mesh.h"
#include "../load/load.h"
#include "../species/species.h"

class SEMO_Controls_Builder{
	public:
		void readSpecies(vector<SEMO_Species>& species);
		void readConfigures();
		void setValues(vector<SEMO_Species>& species);
		
		map<string,int> dict;
		
		string application;
		string startFrom;
		string stopAt;
		double saveInterval;
		string saveControl;
		string saveFormat;
		string timeFormat;
		string adjustTimeStep;
		int saveCompression;
		double maxCo;
		double maxVfCo;
		double maxTimeStep;
		double orgTimeStep;
		int writePrecision;
		
		// user save datas
		map<string,bool> saveMeshData;
		map<string,bool> saveGradientData;
		map<string,bool> saveThermodynamicData;
		map<string,bool> saveBodyForceData;
		
		// extract 관련
		double residual;
		vector<string> extractFieldDatas;
		vector<string> extractAverageDatas;
		vector<string> extractNames;
		vector<string> extractCellValueTargets;
		vector<vector<double>> extractCenterPoints;
		vector<double> extractRadii;
		vector<vector<double>> extractDatas;
		
		string fluxScheme;
		string gradScheme;
		string divScheme;
	
	
		
		vector<double> gravityAcceleration;
		
		int nTotalCellVar;
		int nEq;
		int nSp;
		int P;
		int U;
		int V;
		int W;
		int T;
		vector<int> VF;
		vector<int> MF;
		int Rho;
		int C;
		int Ht;
		vector<int> indicatorAMR;
		
		// level-set
		int LS;

		// vector<int> Q;
		vector<int> Qn;
		vector<int> Qm;
		
		int oldP;
		int oldU;
		int oldV;
		int oldW;
		int oldT;
		vector<int> oldVF;
		vector<int> oldMF;
		int oldRho;
		int oldHt;
		
		int dRhoDP;
		int dHtDP;
		int dRhoDT;
		int dHtDT;
		vector<int> dRhoDVF;
		vector<int> dHtDVF;
		vector<int> dRhoDMF;
		vector<int> dHtDMF;
		
		int dtPseudo;
		int Ur;
		
		vector<int> UDV;
		
		int Un;
		int oldUn;
		// int wC;
		
		
		// int RCM;
		// int invRCM;
		
		double Uco;
		double Lch;
		
		// // bc
		// double subinletU;
		// double subinletV;
		// double subinletW;
		// vector<double> subinletVF;
		
		// double suboutletP;
		
		// double supinletP;
		// double supinletU;
		// double supinletV;
		// double supinletW;
		// vector<double> supinletVF;
		
		
		vector<string> name;
		
		// face LR var
		int nTotalFaceLRVar;
		int fP;
		int fU;
		int fV;
		int fW;
		int fT;
		vector<int> fVF;
		vector<int> fMF;
		int fRho;
		int fC;
		int fHt;
		int fmu;
		
		int fdRhoDP;
		int fdRhoDT;
		int fdHtDP;
		int fdHtDT;
		vector<int> fdRhoDMF;
		vector<int> fdHtDMF;
		
		vector<int> fMF_HO;
		int fRho_HO;
		int fHt_HO;
		
		int fVF_NVD;
	
		
		
		// face var
		int nTotalFaceVar;
		// vector<int> fDistCells;
		
		int dUdx;
		int dUdy;
		int dUdz;
		int dVdx;
		int dVdy;
		int dVdz;
		int dWdx;
		int dWdy;
		int dWdz;
		
		// 최대 최소 값
		int maximumP;
		int minimumP;
		int maximumU;
		int minimumU;
		int maximumV;
		int minimumV;
		int maximumW;
		int minimumW;
		int maximumT;
		int minimumT;
		vector<int> maximumVF;
		vector<int> minimumVF;
		vector<int> maximumMF;
		vector<int> minimumMF;
		
		// turbulence
		string turbType;
		string turbLESModel;
		string turbRANSModel;
		
		
		// solvers
		string solverP;
		double toleranceP;
		double relTolP;
		string preconditionerP;
		int maxIterP;
		
		string solverFinalP;
		double toleranceFinalP;
		double relTolFinalP;
		string preconditionerFinalP;
		int maxIterFinalP;
		
		string solverU;
		double toleranceU;
		double relTolU;
		string preconditionerU;
		int maxIterU;
		
		string solverFinalU;
		double toleranceFinalU;
		double relTolFinalU;
		string preconditionerFinalU;
		int maxIterFinalU;
		
		vector<string> solverVF;
		vector<double> toleranceVF;
		vector<double> relTolVF;
		vector<string> preconditionerVF;
		vector<int> maxIterVF;
		
		vector<string> solverFinalVF;
		vector<double> toleranceFinalVF;
		vector<double> relTolFinalVF;
		vector<string> preconditionerFinalVF;
		vector<int> maxIterFinalVF;
		
		
		// URF
		double momVelURF;
		string momVelAdjustRF;
		vector<int> momVelAdjustSteps;
		vector<double> momVelAdjustValues;
		
		double prePreURF;
		string prePreAdjustRF;
		vector<int> prePreAdjustSteps;
		vector<double> prePreAdjustValues;
		
		double preVelURF;
		string preVelAdjustRF;
		vector<int> preVelAdjustSteps;
		vector<double> preVelAdjustValues;
		
		double vofVofURF;
		string vofVofAdjustRF;
		vector<int> vofVofAdjustSteps;
		vector<double> vofVofAdjustValues;
		
		// URF
		double dualTimeURF_P;
		string dualTimeAdjustRF_P;
		vector<int> dualTimeAdjustSteps_P;
		vector<double> dualTimeAdjustValues_P;
		
		double dualTimeURF_U;
		string dualTimeAdjustRF_U;
		vector<int> dualTimeAdjustSteps_U;
		vector<double> dualTimeAdjustValues_U;
		
		double dualTimeURF_T;
		string dualTimeAdjustRF_T;
		vector<int> dualTimeAdjustSteps_T;
		vector<double> dualTimeAdjustValues_T;
		
		double dualTimeURF_MF;
		string dualTimeAdjustRF_MF;
		vector<int> dualTimeAdjustSteps_MF;
		vector<double> dualTimeAdjustValues_MF;
		
		// limiter
		double maxP;
		double minP;
		double maxU;
		double minU;
		double maxV;
		double minV;
		double maxW;
		double minW;
		double maxT;
		double minT;
		
		
		// scheme
		double time;
		double timeStep;
		double oldTimeStep;
		double old2TimeStep;
		double pseudoCo;
		double specifiedCFL;
		double allowableCFL;
		
		// iterator
		int iterRealMax;
		int iterReal;
		
		int iterPBs;
		int iterPBsMax;
		int iterMom;
		int iterMomMax;
		int iterPre;
		int iterPreMax;
		int iterVof;
		int iterVofMax;
	
		
		int iterPseudo;
		int iterPseudoMax;
		
		int iterTotal;
		double startClock;
		
		// transport & turbulence
		int mu;
		int muT;
		int muEffective;
		int k;
		int kEffective;
		int D;
		int DEffective;
		int cv;
		int cp;
		int kSGS;
		
		// source terms
		vector<int> sourceGravity;
		vector<int> sourceSurfaceTension;
		

		double PrT;
		double ScT;
		
		int kappa;
		
		// dynamic Mesh
		vector<double> indicatorCriterion;
		int intervalRefine;
		// double indicatorRefine;
		int maxLevelRefine;
		int maxCellsRefine;
		double minVolumeRefine;
		
		int bufferLayer;
		
		int intervalUnrefine;
		// double indicatorUnrefine;
		
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


};







class MASCH_Control {
public:
	// map<string, string> boundaryMap;
	map<string, string> initialMap;
	map<string, string> physicsMap;
	map<string, string> controlDictMap;
	map<string, string> dynamicMeshMap;
	map<string, string> extractDatasOverTimeMap;
	map<string, string> fvSchemeMap;
	map<string, string> fvSolutionMap;
	map<string, string> speciesMap;
	map<string, string> thermophysicalProperties;
	map<string, string> turbulenceProperties;
	
	// map<string, string> boundary_p;
	// map<string, string> boundary_u;
	// map<string, string> boundary_v;
	// map<string, string> boundary_w;
	// map<string, string> boundary_T;
	// map<string, string> boundary_k;
	// map<string, string> boundary_epsilon;
	// map<string, string> boundary_omega;
	// map<string, string> boundary_Y;
	vector<map<string, string>> boundaryMap;
	
	string language, adjustTimeStep, againWhenDiverge, saveControl;
	string startFrom,saveFormat, turbType, LESModel, RANSModel;
	double stopAt, timeStep, maxCFL, maxTimeStep, multiCFL, minTimeStep;
	int saveInTimeStep, saveCompression, writePrecision;
	int nSp;
	double saveInRunTime;
	vector<string> saveMeshData, saveGradientData, saveThermodynamicData, saveBodyForceData;
	
	
	map<string, int> varDict;
	vector<string> species;
	vector<pair<string, int>> variables;
	vector<pair<string, int>> primitive;
	vector<string> primitive_abb;
	vector<string> primitive_role;
};


	