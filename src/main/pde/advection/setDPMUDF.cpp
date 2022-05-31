
#include "../../../others/solvers.h"

void MASCH_Solver::setDPMFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
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
    
    
    
    
    // 수정해야됨 지금은 한 종의 parcel밖에 안됨
    vector<string> names = load.extractVector(controls.controlParcelsMap["name"]);
    
    
    
    
    
	double rho_drop = 0.0;
    if(names.size()>0) rho_drop = stod(controls.controlParcelsMap[names[0]+".thermodynamics.rho"]);
	double T_drop = 0.0;
    if(names.size()>0) T_drop = stod(controls.controlParcelsMap[names[0]+".thermodynamics.T"]);
    
    
    double minBreakupDiameter = stod(controls.controlParcelsMap["minBreakupDiameter"]);
    double minBreakupChildDiameter = stod(controls.controlParcelsMap["minBreakupChildDiameter"]);
    double minBreakupChildVolume = stod(controls.controlParcelsMap["minBreakupChildVolume"]);
    double maxBreakupNumberOfParcel = stod(controls.controlParcelsMap["maxBreakupNumberOfParcel"]);

    
    
    
    
    
	int id_rho_p = controls.getId_parcelVar("density");
	int id_u_p = controls.getId_parcelVar("x-velocity");
	int id_v_p = controls.getId_parcelVar("y-velocity");
	int id_w_p = controls.getId_parcelVar("z-velocity");
	int id_T_p = controls.getId_parcelVar("temperature");
	int id_N_p = controls.getId_parcelVar("number-of-parcel");
	int id_time_p = controls.getId_parcelVar("time");
	int id_d_p = controls.getId_parcelVar("diameter");
	int id_x_p = controls.getId_parcelVar("x-position");
	int id_y_p = controls.getId_parcelVar("y-position");
	int id_z_p = controls.getId_parcelVar("z-position");
	
	
	
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_T = controls.getId_cellVar("temperature");
	int id_rho = controls.getId_cellVar("density");
	int id_mu = controls.getId_cellVar("viscosity");
	
	
	
	controls.calcDPM_iSeg.push_back(0);
	
	
	controls.idSetLagrangianEulerian.push_back(make_pair(id_u_p, id_u));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_v_p, id_v));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_w_p, id_w));
	controls.idSetLagrangianEulerian.push_back(make_pair(id_T_p, id_T));
	
	
	solver.calcDPM_parcelLoop.push_back(
		[&controls, &mesh, &var, id_u,id_v,id_w,id_u_p,id_v_p,id_w_p,id_rho_p,id_rho,id_d_p,id_dt,id_dt_p,
		id_x_p,id_y_p,id_z_p,id_mu,id_N_p,id_time_p,id_T_p,
        minBreakupDiameter,minBreakupChildDiameter,
        minBreakupChildVolume,maxBreakupNumberOfParcel](
		int icell, double* cells, double* fields, double* parcels, double* fluxA, double* fluxB) ->int {
			double dt_p = fields[id_dt_p];
			
			double u = cells[id_u];
			double v = cells[id_v];
			double w = cells[id_w];
			double u_p = parcels[id_u_p];
			double v_p = parcels[id_v_p];
			double w_p = parcels[id_w_p];
			double rho_p = parcels[id_rho_p];
			double N_p = parcels[id_N_p];
			double time_p = parcels[id_time_p];
			// relative velocity
			double urel = u - u_p;
			double vrel = v - v_p;
			double wrel = w - w_p;
			double urelAbs = sqrt(urel*urel+vrel*vrel+wrel*wrel);

			// Reynolds Number
			double mu = cells[id_mu];//0.001;
			double rho = cells[id_rho];
			double d_p = parcels[id_d_p];
            double d_p0 = d_p;
			double Red = rho * d_p * urelAbs / mu;

			// // Putnam's drag model
			// double Cd = 0.424;
			// if(Red<1000.0) {
				// Cd = 24.0 / (Red+1.e-16) * (1.0 + pow(Red,0.666666666)/6.0);
			// }
			// if(Red<1.e-5) Cd = 0.0;
            
            // Drag coefficient, Morrison (2013), Reza Barati, 2014
            double Cd = 
                 24.0/Red + 2.6*(Red/5.0)/(1.0+pow(Red/5.0,1.52)) + 
                 0.411*pow(Red/263000.0,-7.94)/(1.0+pow(Red/263000.0,-8.0)) + 
                 pow(Red,0.8)/461000.0;
                 
            // double Cd = 
                 // 24.0/Red + 2.6*(Red/5.0)/(1.0+pow(Red/5.0,1.52)) + 
                 // 0.411*pow(Red/263000.0,-7.94)/(1.0+pow(Red/263000.0,-8.0)) + 
                 // 0.25*(Red/1.e6)/(1.0+(Red/1.e6));
                 



			// // vaporization effect, Han et al. eq.15
			// if(parcel(i)%BM < 0.78) then
				// Cd = Cd / (R1 + parcel(i)%BM)
			// else
				// Cd = Cd / (R1 + parcel(i)%BM)**(0.75)   ! Sazhin et al. Int.Jou.of Thermal Sciences 44,610,2005
			// end if
            
            double force_x = 0.0;
            double force_y = 0.0;
            double force_z = 0.0;
            
			if(Red>1.0) {

                // Drag is negative force
                double coeff = Cd * rho * 3.141592 * d_p*d_p * urelAbs / 8.0;
                force_x = coeff * urel;
                force_y = coeff * vrel;
                force_z = coeff * wrel;

                parcels[id_x_p] += 0.5*u_p*dt_p;
                parcels[id_y_p] += 0.5*v_p*dt_p;
                parcels[id_z_p] += 0.5*w_p*dt_p;
                
                double mass = 4.0/3.0*3.141592*0.125*d_p*d_p*d_p * rho_p;
                parcels[id_u_p] += force_x / mass *dt_p;
                parcels[id_v_p] += force_y / mass *dt_p;
                parcels[id_w_p] += force_z / mass *dt_p;
                
                // 2nd order location
                parcels[id_x_p] += 0.5*parcels[id_u_p]*dt_p;
                parcels[id_y_p] += 0.5*parcels[id_v_p]*dt_p;
                parcels[id_z_p] += 0.5*parcels[id_w_p]*dt_p;

                
                double energy = force_x*parcels[id_u_p]+force_y*parcels[id_v_p]+force_z*parcels[id_w_p];
                // double energy = h*3.141592*d^2*(Tg-Td);
                
                
                
                fluxB[1] -= force_x * N_p;
                fluxB[2] -= force_y * N_p;
                fluxB[3] -= force_z * N_p;
                fluxB[4] -= energy * N_p;
                
                
                // if(N_p>10000.0) { cout.precision(10); cout << N_p << " " << d_p << endl; }
                // if(N_p>10000.0)  mesh.parcels.back().setType(MASCH_Parcel_Types::INSIDE);
                
            }
            else{
                parcels[id_u_p] = u;
                parcels[id_v_p] = v;
                parcels[id_w_p] = w;
                
                parcels[id_x_p] += u*dt_p;
                parcels[id_y_p] += v*dt_p;
                parcels[id_z_p] += w*dt_p;
            }
            
            
            // time update
            parcels[id_time_p] += dt_p;
            
            
            
            
    // double minBreakupDiameter = stod(controls.controlParcelsMap["minBreakupDiameter"]);
    // double minBreakupChildDiameter = stod(controls.controlParcelsMap["minBreakupChildDiameter"]);
    // double minBreakupChildVolume = stod(controls.controlParcelsMap["minBreakupChildVolume"]);
    // double maxBreakupNumberOfParcel = stod(controls.controlParcelsMap["maxBreakupNumberOfParcel"]);
            
            // double workingMinBreakupDiam = 1.e-8;
            // double workingMinAfterBreakupDiam = 1.e-8;
            // double minChildDropVol = 1.e-18;
            // double maxBreakupNParcel = 1000.0;
            
            
            if(d_p>minBreakupDiameter && N_p < maxBreakupNumberOfParcel){
                //============================
                // breakup model (KH-RT)
                double C_tao = 1.0;
                double sigma = 0.0728;
                double mu_p = 0.001;
                
                double B0 = 0.61;
                double B1 = 10.0;
                double A1 = 0.188;
                double C3 = 0.2;
                
                double r_p = 0.5*d_p;
                double We = urelAbs*urelAbs * r_p / sigma;
                double We_p = rho_p * We;
                double We_c = rho * We;
                double Re_p = rho_p * urelAbs * r_p/ mu_p;
                
                double Zo = sqrt(We_p) / Re_p;
                double Ta = Zo * sqrt(We_c);
                double du_v[3];
                du_v[0] = -urel;
                du_v[1] = -vrel;
                du_v[2] = -wrel;
                double N_p0 = N_p;
                
                double double13 = 1.0/3.0;
                
                //######################## Kelvin-Helmholtz ###########################
                double lemtaKH = 9.02 * (1.0 + 0.45 * sqrt(Zo))*(1.0 + 0.4 * pow(Ta,0.7)) /  
                               pow(1.0 + 0.87 * pow(We_c,1.67), 0.6);
                double omegaKH = (0.34 + 0.38 * We_c * sqrt(We_c)) / 
                               (( 1.0 + Zo ) * (1.0 + 1.4 * pow(Ta,0.6)));
                
                lemtaKH = lemtaKH * r_p;
                omegaKH = omegaKH * sqrt(sigma /(rho_p * r_p*r_p*r_p));
                
                double r_child_KH = B0 * lemtaKH;

                if (r_child_KH > r_p) {
                   double r1 = 3.141592*urelAbs/(2.0*omegaKH);
                   double r2 = lemtaKH/4.0;
                   r_child_KH = pow(min(r1,r2) * 3.0 * r_p*r_p, double13);
                }
                //######################## Rayleigh-Taylor ##########################

                // call spray_coupling(mp_dot, f_p, q_dot,cdd,y_at)

                double m_p = 3.141592 / 6.0 * d_p*d_p*d_p * rho_p;
                double f_p_sc = sqrt(force_x*force_x+force_y*force_y+force_z*force_z);
                double g_c = f_p_sc / m_p; 

                double lemtaRT = 3.141592*sqrt( (3.0*sigma) / abs(g_c*(rho_p-rho)) );         
                double omegaRT = sqrt( ( pow(abs((g_c)*(rho_p-rho)), 1.5) ) / (rho_p + rho) );
                omegaRT = omegaRT*sqrt(2.0/(3.0*sqrt(3.0*sigma)));

                double tao_RT = C_tao / omegaRT;
                double r_child_RT = C3 * lemtaRT;     // in the Mixture formation in IC Engine

                // ======================
                // KH-RT breakup
                double somega = 0.0;
                double lemta = 0.0;
                double d_child = 0.0;
                double v0 = 0.0;
                bool boolBreakupWork = false;
                // RT model
                if (time_p > tao_RT) 
                {
                    if (r_child_RT <= r_p) {
                        v0  = N_p * d_p*d_p*d_p;  // mass conservation, 
                        // D_p is the original droplet diameter
                        r_p = r_child_RT;
                        d_child = 2.0 * r_child_RT;
                        d_p = d_child;           // D_p is changed to child diameter
                        // double v0_parent = d_p*d_p*d_p;        // v0_parent = v0_child
                        // double v0_child  = v0_parent;

                        somega = omegaRT;
                        lemta = lemtaRT;
                        boolBreakupWork = true;
                    }
                }
                // KH model
                else {         

                    if(r_child_KH <= r_p) {
                        v0  = N_p0 * d_p*d_p*d_p;  // mass conservation

                        double tao_KH  = 3.726 * B1 * r_p / (lemtaKH * omegaKH);

                        // explicit
                        double dadt =  -(r_p - r_child_KH) / tao_KH;
                        // D_p = D_p + 2 * dadt * dt_p
                        //    implicit for 
                        double dtodtao =  dt_p / tao_KH;
                        r_p = (r_p + r_child_KH * dtodtao)/(1.0 + dtodtao);
                        d_child = 2.0 * r_child_KH;
                        d_p = 2.0 * r_p;

                        // double v0_parent = d_p*d_p*d_p;
                        // double v0_child  = d_child*d_child*d_child;

                        somega = omegaKH;
                        lemta = lemtaKH;
                        boolBreakupWork = true;
                    }
                }
                
        
                if(boolBreakupWork==true){
                        // cout << "B1" << endl;
                    double dm = 0.0;
                    double n_child = 0.0;
                    if ( d_p <= minBreakupChildDiameter ) { 
                        d_p = d_p0;//pow(v0 / N_p, 0.333333333333);
                        N_p = N_p0;                      // n_p0 is original number of particle
                        n_child = 0.0;
                    }
                    else{
                        N_p = v0 / (d_p*d_p*d_p);           // original volume / child diameter          
                        dm  = (N_p - N_p0)/N_p * v0;     // mass will be shedded  
                        n_child = dm / (d_child*d_child*d_child);
                    }
                        // cout << "B2" << endl;
                    
                    bool boolChildDropAdd = false;
                    int nParcelVar = controls.nParcelVar;
                    vector<double> vars(nParcelVar,0.0);
                    if (0.5236*dm > minBreakupChildVolume && n_child > N_p0) 
                    { 
                        boolChildDropAdd = true;
                    
                        // cout << "B3" << endl;
                        // random
                        std::random_device rd;
                        std::default_random_engine eng(rd());
                        std::uniform_real_distribution<double> distr(0.0, 1.0);
                        double ran0 = distr(eng);
                        double ran1 = distr(eng);
                        
                        // cout << "B4" << endl;

                        // new child droplet/parcel
                        // id = new_drop()
                        // child parcel 생성
                        
                        // child droplet velocity in KH-RT model
                        double theta = atan(A1 * lemta * somega / urelAbs)*2.0;
                        double ran_theta = ran0 * theta;
                        double u_r = tan(ran_theta) * urelAbs;    // du is relative velocity magitude
                        
                        // calculate child drop velocity.  v1 is a unit vector orthogonal
                        //   to du_v,  v2 is a 2nd vector orthogonal to du_v and v1
                        //  cf) du_v is relative velocity, vector
                        ran_theta = ran1;
                        double phi = ran_theta * 2.0 * 3.141592; // double_pi = 2 * pi
                        double r1 = du_v[0]/(urelAbs*urelAbs); // du2 = du**2
                        double v1[3];
                        v1[0] = 1.0 - r1 * du_v[0];
                        v1[1] = - r1 * du_v[1];
                        v1[2] = - r1 * du_v[2];

                        // if v1 is parallel to du_v
                        double sum_v1 = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
                        if (sum_v1<=1.e-5) v1[1] = 1.0;

                        r1 = 1.0 / sqrt(sum_v1);
                        // normalized
                        v1[0] = v1[0] * r1; v1[1] = v1[1] * r1; v1[2] = v1[2] * r1;

                        double v2[3];
                        double r2 = 0.0;
                        if (urelAbs != 0.0){
                            // v2 = cshift(du_v,1)*cshift(v1,-1)-cshift(du_v,-1)*cshift(v1,1)
                            v2[0] = du_v[1]*v1[2]-du_v[2]*v1[1];
                            v2[1] = du_v[2]*v1[0]-du_v[0]*v1[2];
                            v2[2] = du_v[0]*v1[1]-du_v[1]*v1[0];
                            
                            double sum_v2 = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2];
                            r2 = 1.0 / sqrt(sum_v2);
                            // normalized
                            v2[0] = v2[0] * r2; v2[1] = v2[1] * r2; v2[2] = v2[2] * r2;
                        }
                        else{
                            v2[0] = 0.0; v2[1] = 0.0; v2[2] = 1.0;
                        }

                        r1 = u_r * sin(phi);
                        r2 = u_r * cos(phi);
                        double u_child[3];
                        // u_child[0] = u_p + v1[0]*r1 + v2[0]*r2;
                        // u_child[1] = v_p + v1[1]*r1 + v2[1]*r2;
                        // u_child[2] = w_p + v1[2]*r1 + v2[2]*r2; 
                        u_child[0] = u_p;
                        u_child[1] = v_p;
                        u_child[2] = w_p; 
                        
                        // call set_parcel(Par(id),node=nd,n=n_child,n0=n_child,i=ijk, & !IJK_RESET, &
                        // state=SPRAY_Normal,u=u_child,x=x_p,d=d_child,&
                        // t=T_p,yat=0.0_dk,vat=0.0_dk,tp_at=0.0_dk,ttp_at=total_timep,injstat=injstat_p,&
                        // DROPLET=2.d0)
                        // call plocation(Par(id:id))
                        
                        // cout << "A1" << endl;
                        
                        // int pres_icell = icell;
                        // int id = 1;
                        // int nParcelVar = controls.nParcelVar;
                        // vector<double> vars(nParcelVar,0.0);
                        vars[id_x_p] = parcels[id_x_p]; 
                        vars[id_y_p] = parcels[id_y_p]; 
                        vars[id_z_p] = parcels[id_z_p];
                        vars[id_u_p] = u_child[0]; 
                        vars[id_v_p] = u_child[1]; 
                        vars[id_w_p] = u_child[2];
                        vars[id_rho_p] = parcels[id_rho_p];
                        vars[id_d_p] = d_child;
                        vars[id_N_p] = n_child;
                        vars[id_time_p] = 0.0;
                        
                        // reset parent droplet
                        N_p = N_p0;
                        d_p = pow((v0 - d_child*d_child*d_child * n_child)/N_p, double13);
                        
                        // cout << "A3.1" << endl;

                    }
                    
                        // cout << "A4" << " " << N_p << " " << d_p << endl;
                    parcels[id_N_p] = N_p;
                    parcels[id_d_p] = d_p;
                    parcels[id_time_p] = 0.0;
                    
                    if(boolChildDropAdd==true){
                        int id = 1;
                        mesh.addParcel(icell, id);
                        var.addParcel(nParcelVar, vars);
                        
                    }
                    
                        // cout << "A5" << endl;
                }
            }

			
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