
#include "../../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		int nSp = controls.nSp;
		
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
			id_Y[i] = controls.getId_cellVar("mass-fraction-"+tmp_name);
		}
		
		calcUpdatePrim.push_back(
		[id_p,id_u,id_v,id_w,id_T,id_Y,
		nSp](
		double* cells, double* Xvalues) ->int {
			cells[id_p] += Xvalues[0];
			cells[id_u] += Xvalues[1];
			cells[id_v] += Xvalues[2];
			cells[id_w] += Xvalues[3];
			cells[id_T] += Xvalues[4];
			double tmp_value = 0.0;
			for(int i=0; i<nSp-1; ++i){
				cells[id_Y[i]] += Xvalues[5+i];
				cells[id_Y[i]] = max(0.0,min(1.0,cells[id_Y[i]]));
				tmp_value += cells[id_Y[i]];
			}
			cells[id_Y[nSp-1]] = 1.0-tmp_value;
			
		}); 
		
	}
	
	
}
