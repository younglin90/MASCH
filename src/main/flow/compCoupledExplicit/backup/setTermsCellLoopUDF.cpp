
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTermsCellLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		int nSp = controls.spName.size();
	
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		vector<int> id_Y, id_Y_old;
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_Y.push_back(controls.getId_cellVar(tmp_name));
			id_Y_old.push_back(controls.getId_cellVar("old "+tmp_name));
		}
		
		int id_p_old = controls.getId_cellVar("old pressure");
		int id_u_old = controls.getId_cellVar("old x-velocity");
		int id_v_old = controls.getId_cellVar("old y-velocity");
		int id_w_old = controls.getId_cellVar("old z-velocity");
		int id_rho_old = controls.getId_cellVar("old density");
		int id_Ht_old = controls.getId_cellVar("old total-enthalpy");
		
		int id_c = controls.getId_cellVar("speed-of-sound");
		int id_rho = controls.getId_cellVar("density");
		int id_Ht = controls.getId_cellVar("total-enthalpy");
		int id_drhodp = controls.getId_cellVar("partial-density-pressure");
		int id_drhodT = controls.getId_cellVar("partial-density-temperature");
		int id_dHtdp = controls.getId_cellVar("partial-total-enthalpy-pressure");
		int id_dHtdT = controls.getId_cellVar("partial-total-enthalpy-temperature");
		vector<int> id_drhodY, id_dHtdY;
		for(int i=0; i<controls.spName.size()-1; ++i){
			id_drhodY.push_back(controls.getId_cellVar("partial-density-mass-fraction-"+controls.spName[i]));
			id_dHtdY.push_back(controls.getId_cellVar("partial-total-enthalpy-mass-fraction-"+controls.spName[i]));
		}
		
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_u,id_v,id_w,
		id_u_old,id_v_old,id_w_old,nSp,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_Y,id_Y_old,id_rho,id_Ht,id_p,
		id_p_old,id_rho_old,id_Ht_old](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			int iter=0;
			double volume = cells[id_vol];
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
			
			
			return 0;
		}); 
	}
	
	
}