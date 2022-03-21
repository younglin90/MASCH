
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["coupled.relaxationFactors"]);
	{
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		vector<int> id_Y;
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_Y.push_back(controls.getId_cellVar(tmp_name));
		}
		int nSpm1 = controls.spName.size()-1;
		
		double relaxation_factor = stod(s_iter[0]);
		
		calcUpdatePrim.push_back(
		[id_p,id_u,id_v,id_w,id_T,id_Y,relaxation_factor,nSpm1](
		double* fields, double* cells, double* Xvalues) ->int {
			cells[id_p] += relaxation_factor*Xvalues[0];
			cells[id_u] += relaxation_factor*Xvalues[1];
			cells[id_v] += relaxation_factor*Xvalues[2];
			cells[id_w] += relaxation_factor*Xvalues[3];
			cells[id_T] += relaxation_factor*Xvalues[4];
			
			cells[id_p] = max(100.0,cells[id_p]);
			cells[id_T] = max(50.0,cells[id_T]);
			
			for(int isp=0; isp<nSpm1; ++isp){
				cells[id_Y[isp]] += relaxation_factor*Xvalues[5+isp];
				cells[id_Y[isp]] = max(0.0,min(1.0,cells[id_Y[isp]]));
			}
			
		}); 
		
	}
	
}
