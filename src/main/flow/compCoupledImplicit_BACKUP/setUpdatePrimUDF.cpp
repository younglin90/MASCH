
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["coupled.relaxationFactors"]);
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
		int id_T = controls.getId_cellVar("temperature");
		vector<string> T_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.T"]);
		double T_min = stod(T_lim[0]); double T_max = stod(T_lim[1]);
		vector<int> id_Y;
		vector<string> Y_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.Y"]);
		double Y_min = stod(Y_lim[0]); double Y_max = stod(Y_lim[1]);
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_Y.push_back(controls.getId_cellVar(tmp_name));
		}
		int nSpm1 = controls.spName.size()-1;
		
		double relaxation_factor = stod(s_iter[0]);
		
		calcUpdatePrim.push_back(
		[id_p,id_u,id_v,id_w,id_T,id_Y,relaxation_factor,nSpm1,
		u_min,u_max,v_min,v_max,w_min,w_max,T_min,T_max,p_min,p_max,Y_min,Y_max](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_p] += relaxation_factor*Xvalues[0];
			cells[id_u] += relaxation_factor*Xvalues[1];
			cells[id_v] += relaxation_factor*Xvalues[2];
			cells[id_w] += relaxation_factor*Xvalues[3];
			cells[id_T] += relaxation_factor*Xvalues[4];
			
			cells[id_p] = max(p_min,min(p_max,cells[id_p]));
			cells[id_u] = max(u_min,min(u_max,cells[id_u]));
			cells[id_v] = max(v_min,min(v_max,cells[id_v]));
			cells[id_w] = max(w_min,min(w_max,cells[id_w]));
			cells[id_T] = max(T_min,min(T_max,cells[id_T]));
			
			for(int isp=0; isp<nSpm1; ++isp){
				cells[id_Y[isp]] += relaxation_factor*Xvalues[5+isp];
				cells[id_Y[isp]] = max(Y_min,min(Y_max,cells[id_Y[isp]]));
			}
			
		}); 
		
	}
	
}
