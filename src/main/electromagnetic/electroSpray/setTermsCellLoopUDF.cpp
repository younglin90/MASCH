
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTermsCellLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	{
		int nSp = controls.spName.size();
		
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_u_old = controls.getId_cellVar("old x-velocity");
		int id_v_old = controls.getId_cellVar("old y-velocity");
		int id_w_old = controls.getId_cellVar("old z-velocity");
		vector<int> id_alpha, id_alpha_old;
		for(int i=0; i<controls.spName.size()-1; ++i){
			id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+controls.spName[i]));
			id_alpha_old.push_back(controls.getId_cellVar("old volume-fraction-"+controls.spName[i]));
		}
		int id_rhoe = controls.getId_cellVar("charge-density");
		int id_rhoe_old = controls.getId_cellVar("old charge-density");
		
		int id_xFe = controls.getId_cellVar("x-electric-force");
		int id_yFe = controls.getId_cellVar("y-electric-force");
		int id_zFe = controls.getId_cellVar("z-electric-force");
		
		int id_rho = controls.getId_cellVar("density");
		
		
		// 1번째
		solver.calcTemporal.push_back(
		[id_dt,id_vol,nSp,id_alpha,id_alpha_old](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double alpha[nSp], alpha_old[nSp];
			for(int i=0; i<nSp-1; ++i){
				alpha[i] = cells[id_alpha[i]];
				alpha_old[i] = cells[id_alpha_old[i]];
			}
			
			int iter=0;
			for(int i=0; i<nSp-1; ++i){
				for(int j=0; j<nSp-1; ++j){
					if(i==j) fluxA[iter] = vol / dt;
					++iter;
				}
			}
			
			iter = 0;
			for(int i=0; i<nSp-1; ++i){
				fluxB[iter++] = -(alpha[i]-alpha_old[i])*vol/dt;
			}
			
			return 0;
		}); 
		
		
		// 2번째
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_rhoe,id_rhoe_old](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double rhoe = cells[id_rhoe]; double rhoe_old = cells[id_rhoe_old];
			
			fluxA[0] = vol/dt;
			
			fluxB[0] = -(rhoe-rhoe_old)*vol/dt;
			
			return 0;
		}); 
		
		
		// 3번째
		solver.calcTemporal.push_back(
		[id_vol,id_rhoe](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double vol = cells[id_vol];
			double rhoe = cells[id_rhoe];
			
			fluxA[0] = 0.0;
			
			fluxB[0] = -rhoe*vol;
			
			return 0;
		}); 
		
		
		// 4번째
		solver.calcTemporal.push_back(
		[id_vol](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double vol = cells[id_vol];
			
			fluxA[0] = vol;
			
			fluxB[0] = 0.0;
			
			return 0;
		}); 
		
		
		// 5번째
		solver.calcTemporal.push_back(
		[id_vol](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double vol = cells[id_vol];
			
			fluxA[0] = vol;
			
			fluxB[0] = 0.0;
			
			return 0;
		}); 
		
		
		// 6번째
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_rho,id_u,id_u_old,id_v,id_v_old,id_w,id_w_old,
		id_xFe,id_yFe,id_zFe](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double rho = cells[id_rho];
			double u = cells[id_u]; double u_old = cells[id_u_old];
			double v = cells[id_v]; double v_old = cells[id_v_old];
			double w = cells[id_w]; double w_old = cells[id_w_old];
			
			double xFe = cells[id_xFe];
			double yFe = cells[id_yFe];
			double zFe = cells[id_zFe];
			
			int iter=0;
			
			int nEq = 3;
			
			iter = nEq*0+0;
			fluxA[iter] = rho*vol/dt;
			
			iter = nEq*1+1;
			fluxA[iter] = rho*vol/dt;
			
			iter = nEq*2+2;
			fluxA[iter] = rho*vol/dt;
			
			fluxB[0] = -rho*(u-u_old)*vol/dt + (xFe)*vol;
			fluxB[1] = -rho*(v-v_old)*vol/dt + (yFe)*vol;
			fluxB[2] = -rho*(w-w_old)*vol/dt + (zFe)*vol;
			
			return 0;
		}); 
		
		
		// 7번째
		solver.calcTemporal.push_back(
		[](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			fluxA[0] = 0.0;
			
			fluxB[0] = 0.0;
			
			return 0;
		}); 
		
		
		// 8번째
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_rho](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double rho = cells[id_rho];
			
			int iter=0;
			
			int nEq = 3;
			
			iter = nEq*0+0;
			fluxA[iter] = rho*vol/dt;
			
			iter = nEq*1+1;
			fluxA[iter] = rho*vol/dt;
			
			iter = nEq*2+2;
			fluxA[iter] = rho*vol/dt;
			
			fluxB[0] = 0.0;
			fluxB[1] = 0.0;
			fluxB[2] = 0.0;
			
			return 0;
		}); 
	}
	
	
}