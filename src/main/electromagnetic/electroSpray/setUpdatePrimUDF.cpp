
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["segregated.relaxationFactors"]);
	{
		int id_p = controls.getId_cellVar("pressure");
		vector<string> p_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.p"]);
		double p_min = stod(p_lim[0]); double p_max = stod(p_lim[1]);
		int id_u = controls.getId_cellVar("x-velocity");
		vector<string> u_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.u"]);
		double u_min = stod(u_lim[0]); double u_max = stod(u_lim[1]);
		int id_v = controls.getId_cellVar("y-velocity");
		vector<string> v_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.v"]);
		double v_min = stod(v_lim[0]); double v_max = stod(v_lim[1]);
		int id_w = controls.getId_cellVar("z-velocity");
		vector<string> w_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.w"]);
		double w_min = stod(w_lim[0]); double w_max = stod(w_lim[1]);
		vector<int> id_alpha;
		vector<string> alpha_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.alpha"]);
		double alpha_min = stod(alpha_lim[0]); double alpha_max = stod(alpha_lim[1]);
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_alpha.push_back(controls.getId_cellVar(tmp_name));
		}
		int id_dp = controls.getId_cellVar("delta-pressure");
        
        
		int id_rhoe = controls.getId_cellVar("charge-density");
		int id_xFe = controls.getId_cellVar("x-electric-force");
		int id_yFe = controls.getId_cellVar("y-electric-force");
		int id_zFe = controls.getId_cellVar("z-electric-force");
        
		int nSpm1 = controls.spName.size()-1;
		
		double relaxation_factor0 = stod(s_iter[0]);
		double relaxation_factor1 = stod(s_iter[1]);
		double relaxation_factor2 = stod(s_iter[2]);
		double relaxation_factor3 = stod(s_iter[3]);
		double relaxation_factor4 = stod(s_iter[4]);
		double relaxation_factor5 = stod(s_iter[5]);
		double relaxation_factor6 = stod(s_iter[6]);
		double relaxation_factor7 = stod(s_iter[7]);
		
		// 1번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			for(int iSp=0; iSp<nSp-1; ++iSp){
				cells[id_alpha[iSp]] += relaxation_factor0*Xvalues[0];
			}
		}); 
		
		// 2번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_rhoe] += relaxation_factor1*Xvalues[0];
		}); 
		
		// 3번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_phie] += relaxation_factor2*Xvalues[0];
		}); 
		
		// 4번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_E] += relaxation_factor3*Xvalues[0];
		}); 
		
		// 5번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_Fe] += relaxation_factor4*Xvalues[0];
		}); 
		
		// 6번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_u] += relaxation_factor5*Xvalues[0];
			cells[id_v] += relaxation_factor5*Xvalues[1];
			cells[id_w] += relaxation_factor5*Xvalues[2];
		}); 
		
		// 7번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_dp] = relaxation_factor6*Xvalues[0];
			cells[id_p] += cells[id_dp];
		}); 
		
		// 8번째
		calcUpdatePrim.push_back(
		[](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_u] += relaxation_factor7*Xvalues[0];
			cells[id_v] += relaxation_factor7*Xvalues[1];
			cells[id_w] += relaxation_factor7*Xvalues[2];
		}); 
		
		
	}
	
}
