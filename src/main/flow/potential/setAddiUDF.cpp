
#include "../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	using Bound_Face_Funct_type = function<int(double* faces)>;
	using Bound_Cell_Funct_type = function<int(double* cells)>;
	
	Bound_Cell_Funct_type setCellFunction;
	Bound_Face_Funct_type setFaceFunction;
	
	

	// species 관련 재료들
	int nSp = controls.nSp;
	auto& thermoMap = controls.thermophysicalProperties;
	// vector<int> id_spec_type;
	vector<vector<double>> spInf;
	vector<int> id_alpha;
	for(int i=0; i<nSp; ++i){
		spInf.push_back(vector<double>());
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.value";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		// id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+name));
	}
	int id_rho = controls.getId_cellVar("density");
	
	int id_rhoL = controls.getId_faceVar("left density");
	int id_rhoR = controls.getId_faceVar("right density");
	
	
	
	{
		calcCellAddiVal.push_back(
			// [&solver] (
			// double* cells) ->int {
				// return 0;
			// }
			[&solver,nSp,spInf,id_rho] (
			double* cells) ->int {
				// double tmp_rho = 0.0;
				// double tmp_sum = 0.0;
				// for(int i=0; i<nSp-1; ++i){
					// tmp_sum += cells[id_alpha[i]];
					// tmp_rho += cells[id_alpha[i]]*spInf[i][0];
				// }
				// tmp_rho += (1.0-tmp_sum)*spInf[nSp-1][0];
				// cells[id_rho] = tmp_rho;
				cells[id_rho] = spInf[0][0];
				return 0;
			}
		);
		calcFaceAddiVal.push_back( 
			[&solver,spInf,id_rhoL,id_rhoR] (
			double* faces) ->int {
				faces[id_rhoL] = spInf[0][0];
				faces[id_rhoR] = spInf[0][0];
				return 0;
			}
		);
	}
	
	{
		calcCellAddiVal.push_back(
			[&solver] (
			double* cells) ->int {
				return 0;
			}
		);
		calcFaceAddiVal.push_back( 
			[&solver] (
			double* faces) ->int {
				return 0;
			}
		);
	}
	
}

