
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
		vector<int> id_Y, id_Y_old, id_Y_old2;
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_Y.push_back(controls.getId_cellVar(tmp_name));
			id_Y_old.push_back(controls.getId_cellVar("old "+tmp_name));
			id_Y_old2.push_back(controls.getId_cellVar("old2 "+tmp_name));
		}
		
		int id_p_old = controls.getId_cellVar("old pressure");
		int id_u_old = controls.getId_cellVar("old x-velocity");
		int id_v_old = controls.getId_cellVar("old y-velocity");
		int id_w_old = controls.getId_cellVar("old z-velocity");
		int id_rho_old = controls.getId_cellVar("old density");
		int id_Ht_old = controls.getId_cellVar("old total-enthalpy");
        
		int id_p_old2 = controls.getId_cellVar("old2 pressure");
		int id_u_old2 = controls.getId_cellVar("old2 x-velocity");
		int id_v_old2 = controls.getId_cellVar("old2 y-velocity");
		int id_w_old2 = controls.getId_cellVar("old2 z-velocity");
		int id_rho_old2 = controls.getId_cellVar("old2 density");
		int id_Ht_old2 = controls.getId_cellVar("old2 total-enthalpy");
		
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
        
        
        
        double C_dt_coeff1 = 1.0;
        double C_dt_coeff2 = 1.0;
        double C_dt_coeff3 = 0.0;
        if(controls.fvSchemeMap["temporalScheme"]=="2nd") C_dt_coeff1 = 1.5;
        if(controls.fvSchemeMap["temporalScheme"]=="2nd") C_dt_coeff2 = 2.0;
        if(controls.fvSchemeMap["temporalScheme"]=="2nd") C_dt_coeff3 = 0.5;
        
        
		
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_u,id_v,id_w,
		id_u_old,id_v_old,id_w_old,nSp,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_Y,id_Y_old,id_rho,id_Ht,id_p,
		id_p_old,id_rho_old,id_Ht_old,
        id_Y_old2,id_p_old2,id_u_old2,id_v_old2,id_w_old2,id_rho_old2,id_Ht_old2,
        C_dt_coeff1,C_dt_coeff2,C_dt_coeff3,
		gVec_x,gVec_y,gVec_z](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double p = cells[id_p]; double p_old = cells[id_p_old]; double p_old2 = cells[id_p_old2];
			double u = cells[id_u]; double u_old = cells[id_u_old]; double u_old2 = cells[id_u_old2];
			double v = cells[id_v]; double v_old = cells[id_v_old]; double v_old2 = cells[id_v_old2]; 
			double w = cells[id_w]; double w_old = cells[id_w_old]; double w_old2 = cells[id_w_old2];
			double rho = cells[id_rho]; double rho_old = cells[id_rho_old]; double rho_old2 = cells[id_rho_old2];
			double Ht = cells[id_Ht]; double Ht_old = cells[id_Ht_old]; double Ht_old2 = cells[id_Ht_old2];
			double drhodp = cells[id_drhodp];
			double drhodT = cells[id_drhodT];
			double dHtdp = cells[id_dHtdp];
			double dHtdT = cells[id_dHtdT];
			double Y[nSp], Y_old[nSp], Y_old2[nSp];
			for(int i=0; i<nSp-1; ++i){
				Y[i] = cells[id_Y[i]];
				Y_old[i] = cells[id_Y_old[i]];
				Y_old2[i] = cells[id_Y_old2[i]];
			}
            
            
            double C_dt = C_dt_coeff1;
            // double C_dt = 1.5;
            
			
			int iter=0;
			fluxA[iter++] = C_dt * drhodp * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = C_dt * drhodT * vol / dt;
			for(int j=0; j<nSp-1; ++j){
				fluxA[iter++] = C_dt * cells[id_drhodY[j]] * vol / dt;
			}
			
			fluxA[iter++] = C_dt * drhodp*u * vol / dt - drhodp*gVec_x*vol;
			fluxA[iter++] = C_dt * rho * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = C_dt * drhodT*u * vol / dt - drhodT*gVec_x*vol;
			for(int j=0; j<nSp-1; ++j){
				double drhodY = cells[id_drhodY[j]];
				fluxA[iter++] = C_dt * drhodY*u * vol / dt - drhodY*gVec_x*vol;
			}
			
			fluxA[iter++] = C_dt * drhodp*v * vol / dt - drhodp*gVec_y*vol;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = C_dt * rho * vol / dt;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = C_dt * drhodT*v * vol / dt - drhodT*gVec_y*vol;
			for(int j=0; j<nSp-1; ++j){
				double drhodY = cells[id_drhodY[j]];
				fluxA[iter++] = C_dt * drhodY*v * vol / dt - drhodY*gVec_y*vol;
			}
			
			fluxA[iter++] = C_dt * drhodp*w * vol / dt - drhodp*gVec_z*vol;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = 0.0;
			fluxA[iter++] = C_dt * rho * vol / dt;
			fluxA[iter++] = C_dt * drhodT*w * vol / dt - drhodT*gVec_z*vol;
			for(int j=0; j<nSp-1; ++j){
				double drhodY = cells[id_drhodY[j]];
				fluxA[iter++] = C_dt * drhodY*w * vol / dt - drhodY*gVec_z*vol;
			}
			
			fluxA[iter++] = C_dt * (drhodp*Ht+rho*dHtdp-1.0) * vol / dt - 
							drhodp*(gVec_x*u+gVec_y*v+gVec_z*w)*vol;
			fluxA[iter++] = C_dt * rho*u * vol / dt - rho*gVec_x*vol;
			fluxA[iter++] = C_dt * rho*v * vol / dt - rho*gVec_y*vol;
			fluxA[iter++] = C_dt * rho*w * vol / dt - rho*gVec_z*vol;
			fluxA[iter++] = C_dt * (drhodT*Ht+rho*dHtdT) * vol / dt - 
							drhodT*(gVec_x*u+gVec_y*v+gVec_z*w)*vol;
			for(int j=0; j<nSp-1; ++j){
				double drhodY = cells[id_drhodY[j]];
				fluxA[iter++] = C_dt * (drhodY*Ht+rho*cells[id_dHtdY[j]]) * vol / dt - 
								drhodY*(gVec_x*u+gVec_y*v+gVec_z*w)*vol;
			}
			
			for(int i=0; i<nSp-1; ++i){
				double Yi = cells[id_Y[i]];
				fluxA[iter++] = C_dt * drhodp*Yi * vol / dt;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = 0.0;
				fluxA[iter++] = C_dt * drhodT*Yi * vol / dt;
				for(int j=0; j<nSp-1; ++j){
					fluxA[iter] = C_dt * cells[id_drhodY[j]]*Yi * vol / dt;
					if(i==j) fluxA[iter] += C_dt * rho * vol / dt;
					++iter;
				}
			}
			
			iter = 0;
            
            
			// fluxB[iter++] = -(rho-rho_old)*vol/dt;
			// fluxB[iter++] = -(rho*u-rho_old*u_old)*vol/dt + gVec_x*rho*vol;
			// fluxB[iter++] = -(rho*v-rho_old*v_old)*vol/dt + gVec_y*rho*vol;
			// fluxB[iter++] = -(rho*w-rho_old*w_old)*vol/dt + gVec_z*rho*vol;
			// fluxB[iter++] = -(rho*Ht-p-rho_old*Ht_old+p_old)*vol/dt +
				// gVec_x*u*rho*vol + gVec_y*v*rho*vol + gVec_z*w*rho*vol;
			// for(int i=0; i<nSp-1; ++i){
				// fluxB[iter++] = -(rho*Y[i]-rho_old*Y_old[i])*vol/dt;
			// }
			
			fluxB[iter++] = -(C_dt_coeff1*rho-C_dt_coeff2*rho_old+C_dt_coeff3*rho_old2)*vol/dt;
			fluxB[iter++] = -(C_dt_coeff1*rho*u-C_dt_coeff2*rho_old*u_old+C_dt_coeff3*rho_old2*u_old2)*vol/dt + gVec_x*rho*vol;
			fluxB[iter++] = -(C_dt_coeff1*rho*v-C_dt_coeff2*rho_old*v_old+C_dt_coeff3*rho_old2*v_old2)*vol/dt + gVec_y*rho*vol;
			fluxB[iter++] = -(C_dt_coeff1*rho*w-C_dt_coeff2*rho_old*w_old+C_dt_coeff3*rho_old2*w_old2)*vol/dt + gVec_z*rho*vol;
			fluxB[iter++] = -(C_dt_coeff1*(rho*Ht-p)-C_dt_coeff2*(rho_old*Ht_old-p_old)+C_dt_coeff3*(rho_old2*Ht_old2-p_old2))*vol/dt +
				gVec_x*u*rho*vol + gVec_y*v*rho*vol + gVec_z*w*rho*vol;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = -(C_dt_coeff1*rho*Y[i]-C_dt_coeff2*rho_old*Y_old[i]+C_dt_coeff3*rho_old2*Y_old2[i])*vol/dt;
			}

			// cout << C_dt_coeff1 << " " << C_dt_coeff2 << " " << C_dt_coeff3 << endl;
			
			
			return 0;
		}); 
	}
	
	
}