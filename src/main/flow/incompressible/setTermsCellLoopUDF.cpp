
#include "../../../others/solvers.h"

// temporal
void MASCH_Solver::setTermsFaceLoopFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	{
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		
		int id_u_old = controls.getId_cellVar("old x-velocity");
		int id_v_old = controls.getId_cellVar("old y-velocity");
		int id_w_old = controls.getId_cellVar("old z-velocity");
		
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_u,id_v,id_w,
		id_u_old,id_v_old,id_w_old](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			double u = cells[id_u]; double u_old = cells[id_u_old];
			double v = cells[id_v]; double v_old = cells[id_v_old]; 
			double w = cells[id_w]; double w_old = cells[id_w_old];
			
			fluxA[0] = vol/dt; fluxA[1] = 0.0; fluxA[2] = 0.0;
			fluxA[3] = 0.0; fluxA[4] = vol/dt; fluxA[5] = 0.0;
			fluxA[6] = 0.0; fluxA[7] = 0.0; fluxA[8] = vol/dt;
			
			fluxB[0] = -(u-u_old)*vol/dt;
			fluxB[1] = -(v-v_old)*vol/dt;
			fluxB[2] = -(w-w_old)*vol/dt;
			return 0;
		}); 
	}
	
	{
		solver.calcTemporal.push_back(
		[](double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			fluxA[0] = 0.0;
			fluxB[0] = 0.0;
			return 0;
		}); 
	}
	
	{
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		
		solver.calcTemporal.push_back(
		[id_dt,id_vol,id_u,id_v,id_w](
		double* cells, double* fields, double* fluxA, double* fluxB) ->int {
			double dt = fields[id_dt];
			double vol = cells[id_vol];
			
			fluxA[0] = -vol/dt; fluxA[1] = 0.0; fluxA[2] = 0.0;
			fluxA[3] = 0.0; fluxA[4] = -vol/dt; fluxA[5] = 0.0;
			fluxA[6] = 0.0; fluxA[7] = 0.0; fluxA[8] = -vol/dt;
			
			// fluxA[0] = 1.0; fluxA[1] = 0.0; fluxA[2] = 0.0;
			// fluxA[3] = 0.0; fluxA[4] = 1.0; fluxA[5] = 0.0;
			// fluxA[6] = 0.0; fluxA[7] = 0.0; fluxA[8] = 1.0;
			
			fluxB[0] = 0.0;
			fluxB[1] = 0.0;
			fluxB[2] = 0.0;
			return 0;
		}); 
	}
	
}