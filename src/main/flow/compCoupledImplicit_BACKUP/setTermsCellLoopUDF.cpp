
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTermsCellLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	{
		int nSp = controls.spName.size();
		
		vector<string> gVecS = load.extractVector(controls.bodyforceMap["g"]);
		double gVec_x = stod(gVecS[0]);
		double gVec_y = stod(gVecS[1]);
		double gVec_z = stod(gVecS[2]);
	
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
		id_p_old,id_rho_old,id_Ht_old,
		gVec_x,gVec_y,gVec_z](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double p = cells[id_p]; double p_old = cells[id_p_old]; 
			double u = cells[id_u]; double u_old = cells[id_u_old];
			double v = cells[id_v]; double v_old = cells[id_v_old]; 
			double w = cells[id_w]; double w_old = cells[id_w_old];
			double rho = cells[id_rho]; double rho_old = cells[id_rho_old];
			double Ht = cells[id_Ht]; double Ht_old = cells[id_Ht_old];
			double drhodp = cells[id_drhodp];
			double drhodT = cells[id_drhodT];
			double dHtdp = cells[id_dHtdp];
			double dHtdT = cells[id_dHtdT];
			double Y[nSp], Y_old[nSp];
			for(int i=0; i<nSp-1; ++i){
				Y[i] = cells[id_Y[i]];
				Y_old[i] = cells[id_Y_old[i]];
			}
			
			int iter=0;
			fluxA[iter++] = drhodp * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = drhodT * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]] * vol / dt;
			}
			
			fluxA[iter++] = drhodp*u * vol / dt;
			fluxA[iter++] = rho * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = drhodT*u * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*u * vol / dt;
			}
			
			fluxA[iter++] = drhodp*v * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = rho * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = drhodT*v * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*v * vol / dt;
			}
			
			fluxA[iter++] = drhodp*w * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = rho * vol / dt;
			fluxA[iter++] = drhodT*w * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = cells[id_drhodY[j]]*w * vol / dt;
			}
			
			fluxA[iter++] = (drhodp*Ht+rho*dHtdp-1.0) * vol / dt;
			fluxA[iter++] = rho*u * vol / dt;
			fluxA[iter++] = rho*v * vol / dt;
			fluxA[iter++] = rho*w * vol / dt;
			fluxA[iter++] = (drhodT*Ht+rho*dHtdT) * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = (cells[id_drhodY[j]]*Ht+
								rho*cells[id_dHtdY[j]]) * vol / dt;
			}
			
			for(int i=0; i<nSp-1; ++i){
				int ii = id_drhodY[i];
				fluxA[iter++] = drhodp*cells[id_Y[i]] * vol / dt;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = drhodT*cells[id_Y[i]] * vol / dt;
				for(int j=0; j<nSp-1; ++j){
					fluxA[iter] = cells[id_drhodY[j]]*cells[id_Y[i]] * vol / dt;
					if(i==j) fluxA[iter] += rho * vol / dt;
					++iter;
				}
			}
			
			iter = 0;
			fluxB[iter++] = -(rho-rho_old)*vol/dt;
			fluxB[iter++] = -(rho*u-rho_old*u_old)*vol/dt + gVec_x*rho*vol;
			fluxB[iter++] = -(rho*v-rho_old*v_old)*vol/dt + gVec_y*rho*vol;
			fluxB[iter++] = -(rho*w-rho_old*w_old)*vol/dt + gVec_z*rho*vol;
			fluxB[iter++] = -(rho*Ht-p-rho_old*Ht_old+p_old)*vol/dt;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = -(rho*Y[i]-rho_old*Y_old[i])*vol/dt;
			}
			
			
			
			return 0;
		}); 
	}
	
	
}