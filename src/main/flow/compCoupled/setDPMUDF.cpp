
#include "../../../others/solvers.h"

void MASCH_Solver::setDPMFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Load load;
	
	int id_t = controls.getId_fieldVar("time");
	int id_dt = controls.getId_fieldVar("time-step");
	int id_dt_p = controls.getId_fieldVar("time-step-parcels");
	
	// double mdot_inj = stod(controls.controlParcelsMap["water.injection.mdot"]);
	// double d_inj = stod(controls.controlParcelsMap["water.injection.d"]);
	// double theta_inj = stod(controls.controlParcelsMap["water.injection.theta"]);
	// vector<string> xyz_inj = load.extractVector(controls.controlParcelsMap["water.injection.position"]);
	// double x_inj = stod(xyz_inj[0]);
	// double y_inj = stod(xyz_inj[1]);
	// double z_inj = stod(xyz_inj[2]);
	double rho_drop = stod(controls.controlParcelsMap["water.thermodynamics.rho"]);
	double T_drop = stod(controls.controlParcelsMap["water.thermodynamics.T"]);
	int id_rho_p = controls.getId_parcelVar("density");
	int id_u_p = controls.getId_parcelVar("x-velocity");
	int id_v_p = controls.getId_parcelVar("y-velocity");
	int id_w_p = controls.getId_parcelVar("z-velocity");
	int id_T_p = controls.getId_parcelVar("temperature");
	int id_nparcel_p = controls.getId_parcelVar("number-of-parcel");
	int id_d_p = controls.getId_parcelVar("diameter");
	int id_x_p = controls.getId_parcelVar("x-position");
	int id_y_p = controls.getId_parcelVar("y-position");
	int id_z_p = controls.getId_parcelVar("z-position");
	
	
	
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_T = controls.getId_cellVar("temperature");
	int id_rho = controls.getId_cellVar("density");
	
	
	
	controls.calcDPM_iSeg.push_back(0);
	
	
	controls.idSetLagrangianEulerian.push_back(make_pair(id_u_p, id_u));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_v_p, id_v));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_w_p, id_w));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_T_p, id_T));
	
	
	solver.calcDPM_parcelLoop.push_back(
		[id_u,id_v,id_w,id_u_p,id_v_p,id_w_p,id_rho_p,id_rho,id_d_p,id_dt,id_dt_p,
		id_x_p,id_y_p,id_z_p](
		double* cells, double* fields, double* parcels, double* fluxB) ->int {
			double dt_p = fields[id_dt_p];
			
			double u = cells[id_u];
			double v = cells[id_v];
			double w = cells[id_w];
			double u_p = parcels[id_u_p];
			double v_p = parcels[id_v_p];
			double w_p = parcels[id_w_p];
			double rho_p = parcels[id_rho_p];
			// relative velocity
			double urel = u - u_p;
			double vrel = v - v_p;
			double wrel = w - w_p;
			double urelAbs = sqrt(urel*urel+vrel*vrel+wrel*wrel);

			// Reynolds Number
			double mu = 0.001;
			double rho = cells[id_rho];
			double d_p = parcels[id_d_p];
			double Red = rho * d_p * urelAbs / mu;

			// Putnam's drag model
			double Cd = 0.424;
			if(Red<1000.0) {
				Cd = 24.0 / (Red+1.e-16) * (1.0 + pow(Red,0.666666666)/6.0);
			}
			if(Red<1.e-5) Cd = 0.0;

			// // vaporization effect, Han et al. eq.15
			// if(parcel(i)%BM < 0.78) then
				// Cd = Cd / (R1 + parcel(i)%BM)
			// else
				// Cd = Cd / (R1 + parcel(i)%BM)**(0.75)   ! Sazhin et al. Int.Jou.of Thermal Sciences 44,610,2005
			// end if

			// Drag is negative force
			double coeff = Cd * rho * 3.141592 * d_p*d_p * urelAbs / 8.0;
			double force_x = coeff * urel;
			double force_y = coeff * vrel;
			double force_z = coeff * wrel;

			parcels[id_x_p] += 0.5*u_p*dt_p;
			parcels[id_y_p] += 0.5*v_p*dt_p;
			parcels[id_z_p] += 0.5*w_p*dt_p;
			
			double mass = 4.0/3.0*3.141592*0.125*d_p*d_p*d_p * rho_p;
			parcels[id_u_p] += force_x / mass *dt_p;
			parcels[id_v_p] += force_y / mass *dt_p;
			parcels[id_w_p] += force_z / mass *dt_p;
			
			// double energy = h*3.141592*d^2*(Tg-Td);
			
			fluxB[1] -= force_x;
			fluxB[2] -= force_y;
			fluxB[3] -= force_z;
			// fluxB[4] -= energy;
			

			
			return 0;
			
		}
	);
	
	
	
	
}



// void MASCH_DPM::calcMechanicalForceModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }
	
// void MASCH_DPM::calcBreakupModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }
	
// void MASCH_DPM::calcHeatAndMassTransferModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }