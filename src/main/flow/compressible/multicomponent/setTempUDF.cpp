
#include "../../../../others/solvers.h"

// temporal
void MASCH_Solver::setTempFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	{
		using US = unsigned short;
		int id_p = controls.cellVar["pressure"].id;
		int id_u = controls.cellVar["x-velocity"].id;
		int id_v = controls.cellVar["y-velocity"].id;
		int id_w = controls.cellVar["z-velocity"].id;
		int id_T = controls.cellVar["temperature"].id;
		
		int nSp = controls.cellVar["mass fraction"].sub_name.size();
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass fraction"].sub_name[i];
			id_Y[i] = controls.cellVar[tmp_name].id;
		}
		
		int id_rho = controls.cellVar["density"].id;
		int id_c = controls.cellVar["speed of sound"].id;
		int id_Ht = controls.cellVar["total enthalpy"].id;
		
		int id_drhodp = controls.cellVar["density diff with pressure"].id;
		int id_drhodT = controls.cellVar["density diff with temperature"].id;
		int id_dHtdp = controls.cellVar["total enthalpy diff with pressure"].id;
		int id_dHtdT = controls.cellVar["total enthalpy diff with temperature"].id;
		vector<int> id_drhodY(nSp);
		vector<int> id_dHtdY(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name1 = controls.cellVar["density diff with mass fraction"].sub_name[i];
			string tmp_name2 = controls.cellVar["total enthalpy diff with mass fraction"].sub_name[i];
			id_drhodY[i] = controls.cellVar[tmp_name1].id;
			id_dHtdY[i] = controls.cellVar[tmp_name2].id;
		}
		
		// 메쉬관련
		int id_volume = controls.cellVar["volume"].id;
		
		// 필드값
		int id_dt = controls.fieldVar["time-step"].id;
		
		solver.calcTemporal.push_back(
		[id_volume,nSp,id_p,id_u,id_v,id_w,id_T,id_Y,
		id_rho,id_c,id_Ht,id_dt,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			int iter=0;
			double volume = cells[id_volume];
			double dt = fields[id_dt];
			
			// cout << cells[id_volume] << " " <<
			// dt << " " <<
			// cells[id_drhodp] << " " <<
			// cells[id_drhodT] << " " <<
			// cells[id_u] << " " <<
			// cells[id_v] << " " <<
			// cells[id_w] << " " <<
			// cells[id_rho] << " " <<
			// cells[id_Ht] << " " <<
			// cells[id_dHtdp] << " " <<
			// cells[id_dHtdT] << " " <<
			// endl;
			
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