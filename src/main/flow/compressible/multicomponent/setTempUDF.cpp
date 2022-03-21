
#include "../../../../others/solvers.h"

// temporal
void MASCH_Solver::setTempFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	{
		using US = unsigned short;
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		
		// int nSp = controls.cellVar["mass-fraction"].sub_name.size();
		int nSp = controls.nSp;
		
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass-fraction"].sub_name[i];
			id_Y[i] = controls.getId_cellVar("mass-fraction-"+tmp_name);
		}
		int id_c = controls.getId_cellVar("speed-of-sound");
		int id_rho = controls.getId_cellVar("density");
		int id_Ht = controls.getId_cellVar("total-enthalpy");
		int id_drhodp = controls.getId_cellVar("density-diff-with-pressure");
		int id_drhodT = controls.getId_cellVar("density-diff-with-temperature");
		int id_dHtdp = controls.getId_cellVar("total-enthalpy-diff-with-pressure");
		int id_dHtdT = controls.getId_cellVar("total-enthalpy-diff-with-temperature");
		vector<int> id_drhodY(nSp);
		vector<int> id_dHtdY(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name1 = controls.cellVar["density-diff-with-mass-fraction"].sub_name[i];
			string tmp_name2 = controls.cellVar["total-enthalpy-diff-with-mass-fraction"].sub_name[i];
			id_drhodY[i] = controls.getId_cellVar("density-diff-with-mass-fraction-"+tmp_name1);
			id_dHtdY[i] = controls.getId_cellVar("total-enthalpy-diff-with-mass-fraction-"+tmp_name2);
		}
		
		
		// 메쉬관련
		int id_volume = controls.getId_cellVar("volume");
		
		// 필드값
		int id_dt = controls.getId_fieldVar("time-step");
		
		solver.calcTemporal.push_back(
		[id_volume,nSp,id_p,id_u,id_v,id_w,id_T,id_Y,
		id_rho,id_c,id_Ht,id_dt,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			int iter=0;
			double volume = cells[id_volume];
			double dt = fields[id_dt];
			
			fluxA[iter++] = cells[id_drhodp] * volume / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = cells[id_drhodT] * volume / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]] * volume / dt;
			}
			
			fluxA[iter++] = cells[id_drhodp]*cells[id_u] * volume / dt;
			fluxA[iter++] = cells[id_rho] * volume / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = cells[id_drhodT]*cells[id_u] * volume / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*cells[id_u] * volume / dt;
			}
			
			fluxA[iter++] = cells[id_drhodp]*cells[id_v] * volume / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = cells[id_rho] * volume / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = cells[id_drhodT]*cells[id_v] * volume / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*cells[id_v] * volume / dt;
			}
			
			fluxA[iter++] = cells[id_drhodp]*cells[id_w] * volume / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = cells[id_rho] * volume / dt;
			fluxA[iter++] = cells[id_drhodT]*cells[id_w] * volume / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*cells[id_w] * volume / dt;
			}
			
			fluxA[iter++] = (cells[id_drhodp]*cells[id_Ht]+
							cells[id_rho]*cells[id_dHtdp]-1.0) * volume / dt;
			fluxA[iter++] = cells[id_rho]*cells[id_u] * volume / dt;
			fluxA[iter++] = cells[id_rho]*cells[id_v] * volume / dt;
			fluxA[iter++] = cells[id_rho]*cells[id_w] * volume / dt;
			fluxA[iter++] = (cells[id_drhodT]*cells[id_Ht]+
							cells[id_rho]*cells[id_dHtdT]) * volume / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = (cells[id_drhodY[j]]*cells[id_Ht]+
								cells[id_rho]*cells[id_dHtdY[j]]) * volume / dt;
			}
			
			for(int i=0; i<nSp-1; ++i){
				int ii = id_drhodY[i];
				fluxA[iter++] = cells[id_drhodp]*cells[id_Y[i]] * volume / dt;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = cells[id_drhodT]*cells[id_Y[i]] * volume / dt;
				for(int j=0; j<nSp-1; ++j){
					fluxA[iter] = cells[id_drhodY[j]]*cells[id_Y[i]] * volume / dt;
					if(i==j) fluxA[iter] += cells[id_rho] * volume / dt;
					++iter;
				}
			}
			
			iter = 0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = 0.0;
			fluxB[iter++] = 0.0;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = 0.0;
			}
			
		}); 
	}
	
	
	
}