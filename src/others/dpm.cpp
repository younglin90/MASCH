
#include "./solvers.h"
#include <amgcl/profiler.hpp>
#include <random>


void MASCH_Solver::dpm(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& solver = (*this);
	
	bool debug_bool = false;
	bool debug_AMGCL_bool = false;
		
	amgcl::profiler<> prof("dpm solver");
	
	if(controls.nameParcels.size()==0) return;
	
	// parcel 추가 루틴
	// for(auto& name : controls.nameParcels)
	{
		// eulerian to lagrangian 
		if(debug_bool) controls.log.push("eulerianToLagrangian");
		if(debug_AMGCL_bool) prof.tic("eulerianToLagrangian");
		// if( (controls.iterReal+1) %
			// stoi(controls.controlParcelsMap[name+".eulerianToLagrangian.interval"]) == 0){
				
			solver.eulerianToLagrangian(mesh, controls, var);
			
		// }
		if(debug_bool) controls.log.pop();
		if(debug_AMGCL_bool) prof.toc("eulerianToLagrangian");
		
		// // parcel injection
		// if(debug_bool) controls.log.push("addParcelModel");
		// if(debug_AMGCL_bool) prof.tic("addParcelModel");
		// solver.addParcelModel(mesh, controls, var);
		// if(debug_bool) controls.log.pop();
		// if(debug_AMGCL_bool) prof.toc("addParcelModel");
	}
        
    // if(mesh.parcels.size()==0) return;
	
	// int nstep = solver.calcTimeStepParcel(mesh, controls, var);
	solver.calcTimeStepParcel(mesh, controls, var);
	
	double resi_dt = var.fields[controls.getId_fieldVar("time-step")];
	// double dt_p = var.fields[controls.getId_fieldVar("time-step-parcels")];
	controls.nIterDPM = 0;
	while(resi_dt>0.0)
    {
		
		resi_dt -= var.fields[controls.getId_fieldVar("time-step-parcels")];
		
		// parcel 루프 계산
		if(debug_bool) controls.log.push("parcelLoop");
		if(debug_AMGCL_bool) prof.tic("parcelLoop");
		solver.parcelLoop(mesh, controls, var);
		if(debug_bool) controls.log.pop();
		if(debug_AMGCL_bool) prof.toc("parcelLoop");
		
		// 이동한 parcel 찾기
		if(debug_bool) controls.log.push("searchLocationParcelsToOutside");
		if(debug_AMGCL_bool) prof.tic("searchLocationParcelsToOutside");
		solver.searchLocationParcelsToOutside(mesh, controls, var);
		if(debug_bool) controls.log.pop();
		if(debug_AMGCL_bool) prof.toc("searchLocationParcelsToOutside");
		
		// 다른 processor로 이동
		if(debug_bool) controls.log.push("updateProcRightParcels");
		if(debug_AMGCL_bool) prof.tic("updateProcRightParcels");
		solver.updateProcRightParcels(mesh, controls, var);
		if(debug_bool) controls.log.pop();
		if(debug_AMGCL_bool) prof.toc("updateProcRightParcels");
		
		// parcel 재정립
		if(debug_bool) controls.log.push("refreshParcels");
		if(debug_AMGCL_bool) prof.tic("refreshParcels");
		solver.refreshParcels(mesh, controls, var);
		if(debug_bool) controls.log.pop();
		if(debug_AMGCL_bool) prof.toc("refreshParcels");
		
		solver.calcTimeStepParcel(mesh, controls, var);
		
		if(resi_dt<var.fields[controls.getId_fieldVar("time-step-parcels")]){
			var.fields[controls.getId_fieldVar("time-step-parcels")] = resi_dt;
		}
		
		++controls.nIterDPM;
		
		// if(rank==0) cout << precision(9) << resi_dt << endl;
		
			// cout << scientific; cout.precision(2);
			// cout << fixed; cout.precision(0);
	}
		
		
	if (rank == 0 && debug_AMGCL_bool) cout << prof << endl;
	
}


int MASCH_Solver::calcTimeStepParcel(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	double cfl = stod(controls.controlParcelsMap["dpmCFL"]);
	double dt = var.fields[controls.getId_fieldVar("time-step")];
	
	int id_dt_p = controls.getId_fieldVar("time-step-parcels");
	
	int id_u_p = controls.getId_parcelVar("x-velocity");
	int id_v_p = controls.getId_parcelVar("y-velocity");
	int id_w_p = controls.getId_parcelVar("z-velocity");
	// int id_x_p = controls.getId_parcelVar("x-location");
	
	int id_vol = controls.getId_cellVar("volume");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	
	auto cellVar = var.cells.data();
	auto parcels_i = mesh.parcels.data();
	auto parcelVar = var.parcels.data();
	
	double dt_p = stod(controls.controlParcelsMap["dpmMaxTimeStep"]);
	for(int ipar=0, SIZE=mesh.parcels.size(); ipar<SIZE; ++ipar){
		auto& parcel = parcels_i[ipar];
		auto parcelVar_i = parcelVar[ipar].data();
		auto cellVar_i = cellVar[parcel.icell].data();
		double u = cellVar_i[id_u];
		double v = cellVar_i[id_v];
		double w = cellVar_i[id_w];
		double U = sqrt(u*u+v*v+w*w);
		double u_p = parcelVar_i[id_u_p];
		double v_p = parcelVar_i[id_v_p];
		double w_p = parcelVar_i[id_w_p];
		double U_p = sqrt(u_p*u_p+v_p*v_p+w_p*w_p);
		double vol = cellVar_i[id_vol];
		dt_p = min(cfl*pow(vol,0.3333)/(U_p+1.e-200), dt_p);
		// dt_p = min(cfl*pow(vol,0.3333)/(U+1.e-200), dt_p);
	}
	if(size>1){
		double tmp_fieldVar = dt_p;
		double tmp_fieldVar_glo;
		MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		var.fields[id_dt_p] = tmp_fieldVar_glo;
	}
	
	int nstep = (int)(dt/var.fields[id_dt_p]);
	
	if(nstep==0) nstep = 1;
	var.fields[id_dt_p] = dt/((double)nstep);
	
	return nstep;
	
	
	
}



void MASCH_Solver::addParcelModel(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var){
	
	
	if(controls.controlParcelsMap.find("water.injection.type")==
	controls.controlParcelsMap.end()) {
		cout << "#WARNING : not defined, water.injection.type" << endl;
		return;
	}
	if(controls.controlParcelsMap["water.injection.type"] == "none"){
		return;
	}
	
	
	
	
	
	double eps = -1.e-16;
	// double eps = -1.e-4;
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	int id_t = controls.getId_fieldVar("time");
	int id_dt = controls.getId_fieldVar("time-step");
	int id_dt_p = controls.getId_fieldVar("time-step-parcels");
	int id_time_accum = controls.getId_fieldVar("parcel-injection-accum-time");
	
	double time = var.fields[id_t];
	double dt = var.fields[id_dt];
	double dt_p = var.fields[id_dt_p];
	
	double mdot_inj = stod(controls.controlParcelsMap["water.injection.mdot"]);
	double d_inj = stod(controls.controlParcelsMap["water.injection.d"]);
	double theta_inj = stod(controls.controlParcelsMap["water.injection.theta"]);
	vector<string> xyz_inj = load.extractVector(controls.controlParcelsMap["water.injection.position"]);
	double x_inj = stod(xyz_inj[0]);
	double y_inj = stod(xyz_inj[1]);
	double z_inj = stod(xyz_inj[2]);
	double rho_drop = stod(controls.controlParcelsMap["water.thermodynamics.rho"]);
	double T_drop = stod(controls.controlParcelsMap["water.thermodynamics.T"]);
	int id_rho = controls.getId_parcelVar("density");
	int id_u = controls.getId_parcelVar("x-velocity");
	int id_v = controls.getId_parcelVar("y-velocity");
	int id_w = controls.getId_parcelVar("z-velocity");
	int id_T = controls.getId_parcelVar("temperature");
	int id_nparcel = controls.getId_parcelVar("number-of-parcel");
	int id_d = controls.getId_parcelVar("diameter");
	int id_x = controls.getId_parcelVar("x-position");
	int id_y = controls.getId_parcelVar("y-position");
	int id_z = controls.getId_parcelVar("z-position");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	
	auto cell_i = mesh.cells.data();
	auto face_i = mesh.faces.data();
	auto faceVar = var.faces.data();
	
	double pi = 3.141592;
	double r_drop = 0.5*d_inj;
	double mass_drop = 4.0/3.0*pi*r_drop*r_drop*r_drop*rho_drop;
	double dt_inj = mass_drop / mdot_inj;
	
	
	if(var.fields[id_time_accum]<dt_inj){
		var.fields[id_time_accum] += dt;
		return;
	}
	
	var.fields[id_time_accum] -= dt_inj;
	
	// cout << "INJECTION " << endl;
	
	double area_hole = pi*d_inj*d_inj/4.0;
	double U_inj = mdot_inj / (area_hole*rho_drop);
	double mass_inj = mdot_inj * dt_inj;//dt_inj;
	double n_drop_per_parcel = mass_inj / mass_drop;
	
	// int icell = solver.get_icell(mesh, controls, var, x_inj,y_inj,z_inj);
	
	
	
	
	
	// for(int i=0; i<100; ++i)
	{
		
		int icell = -1;
	
		// std::random_device rd;
		// std::default_random_engine eng(rd());
		// std::uniform_real_distribution<double> distr(0.0, 1.0);
		// x_inj = distr(eng);
		// y_inj = distr(eng);
		// z_inj = distr(eng);
		// // x_inj = 0.49;
		// // y_inj = 0.47;
		// // z_inj = 0.1;
		
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = mesh.cells[i];
			bool boolInside = true;
			for(auto& iface : cell.ifaces){
				auto& face = mesh.faces[iface];
				auto faceVar_i = faceVar[iface].data();
				double x_pF = x_inj - face.x;
				double y_pF = y_inj - face.y;
				double z_pF = z_inj - face.z;
				double nx = -faceVar_i[id_nx];
				double ny = -faceVar_i[id_ny];
				double nz = -faceVar_i[id_nz];
				if(face.iL!=i) {
					nx = -nx; ny = -ny; nz = -nz;
				}
				double Lcond = x_pF*nx + y_pF*ny + z_pF*nz;
				if(Lcond < eps){
					boolInside = false;
					break;
				}
			}
			if(boolInside==false) continue;
			
			icell = i;
			break;
			
		}
		
		
		if(icell!=-1){
			mesh.cells[icell].iparcels.push_back(mesh.parcels.size());
			mesh.addParcel();
			mesh.parcels.back().icell = icell;
			mesh.parcels.back().setType(MASCH_Parcel_Types::INSIDE);
			
			var.parcels.push_back(vector<double>());
			var.parcels.back().resize(controls.nParcelVar);
			var.parcels.back()[id_x] = x_inj;
			var.parcels.back()[id_y] = y_inj;
			var.parcels.back()[id_z] = z_inj;
			var.parcels.back()[id_rho] = rho_drop;
			var.parcels.back()[id_u] = 0.0;
			var.parcels.back()[id_v] = 0.0;
			var.parcels.back()[id_w] = 0.0;
			var.parcels.back()[id_T] = T_drop;
			var.parcels.back()[id_nparcel] = n_drop_per_parcel;
			var.parcels.back()[id_d] = d_inj;
			
			// cout << n_drop_per_parcel << endl;
		}
	}
			
	
}


void MASCH_Solver::parcelLoop(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
	
	int nSeg_DPM = controls.calcDPM_iSeg.size();
	
	
	for(int iSeg=0; iSeg<nSeg_DPM; ++iSeg){
		
		int iSegEq = controls.calcDPM_iSeg[iSeg];
		
		auto sol = solver.calcDPM_parcelLoop[iSegEq];
		
		auto parcelVar = var.parcels.data();
		auto parcels_i = mesh.parcels.data();
		auto cell_i = mesh.cells.data();
		auto fieldVar = var.fields.data();
		auto faceVar = var.faces.data();
		auto cellVar = var.cells.data();
		
		int nEq = controls.nEq[iSegEq];
        double fluxA[nEq*nEq];
		double fluxB[nEq];
		
		for(int ipar=0, SIZE=mesh.parcels.size(); ipar<SIZE; ++ipar){
			auto& parcel = mesh.parcels[ipar];
			int i = parcel.icell;
			// auto parcelVar_i = parcelVar[ipar].data();
			auto parcelVar_i = var.parcels[ipar].data();
			auto cellVar_i = cellVar[parcel.icell].data();
			
			for(int iEq=0; iEq<nEq; ++iEq){
                for(int jEq=0; jEq<nEq; ++jEq){
                    fluxA[iEq*nEq+jEq] = 0.0;
                }
				fluxB[iEq] = 0.0;
			}
			
			sol(i, cellVar_i, fieldVar, parcelVar_i, fluxA, fluxB);
			
			for(int iEq=0; iEq<nEq; ++iEq){
                for(int jEq=0; jEq<nEq; ++jEq){
                    var.accumSparD( iSegEq, i, iEq, jEq, fluxA[iEq*nEq+jEq] );
                }
				var.accumB( iSegEq, i, iEq, fluxB[iEq] );
			}
			
		}
	
	}

}




void MASCH_Solver::searchLocationParcelsToOutside(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var){
	
	double eps = -1.e-16;
	// double eps = -1.e-4;
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	
	int id_x = controls.getId_parcelVar("x-position");
	int id_y = controls.getId_parcelVar("y-position");
	int id_z = controls.getId_parcelVar("z-position");
	
	auto parcels_i = mesh.parcels.data();
	auto cell_i = mesh.cells.data();
	auto face_i = mesh.faces.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto parcelVar = var.parcels.data();
	
	for(int i=0, SIZE=mesh.nInternalFaces; i<SIZE; ++i){
		auto& face = face_i[i];
		auto faceVar_i = faceVar[i].data();
		int iL = face.iL;
		int iR = face.iR;
		int tmp_iL[2], tmp_iR[2];
		tmp_iL[0] = iL; tmp_iR[0] = iR;
		tmp_iL[1] = iR; tmp_iR[1] = iL;
		
		double nx[2], ny[2], nz[2];
		nx[0] = -faceVar_i[id_nx]; nx[1] = faceVar_i[id_nx];
		ny[0] = -faceVar_i[id_ny]; ny[1] = faceVar_i[id_ny];
		nz[0] = -faceVar_i[id_nz]; nz[1] = faceVar_i[id_nz];
		
		for(int ii=0; ii<2; ++ii){
			auto& cell = cell_i[tmp_iL[ii]];
			for(auto& iparcel : cell.iparcels){
				auto& parcel = parcels_i[iparcel];
				auto parcelVar_i = parcelVar[iparcel].data();
				if(parcel.getType() != MASCH_Parcel_Types::INSIDE) continue;
				double x_pF = parcelVar_i[id_x] - face.x;
				double y_pF = parcelVar_i[id_y] - face.y;
				double z_pF = parcelVar_i[id_z] - face.z;
				double Lcond = x_pF*nx[ii] + y_pF*ny[ii] + z_pF*nz[ii];
				if(Lcond <= eps){
					parcel.icell = tmp_iR[ii];
					parcel.setType(MASCH_Parcel_Types::INSIDE);
				}
			}
		}
	}
	
	for(auto& boundary : mesh.boundaries){
		
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
		
			string type = controls.controlParcelsMap["boundary."+boundary.name+".type"];
			MASCH_Parcel_Types type_p;
			if(type=="reflect"){
				type_p = MASCH_Parcel_Types::REFLECT;
			}
			else if(type=="escape"){
				type_p = MASCH_Parcel_Types::ESCAPE;
			}
			else{
				cout << "#WARNING, not defined DPM boundary type = " << type << endl;
			}
			
			for(int i=str; i<end; ++i){
				auto& face = face_i[i];
				auto faceVar_i = faceVar[i].data();
				int iL = face.iL;
				auto& cell = cell_i[iL];
				
				for(auto& iparcel : cell.iparcels){
					auto& parcel = parcels_i[iparcel];
					auto parcelVar_i = parcelVar[iparcel].data();
					if(parcel.getType() != MASCH_Parcel_Types::INSIDE) continue;
					double x_pF = parcelVar_i[id_x] - face.x;
					double y_pF = parcelVar_i[id_y] - face.y;
					double z_pF = parcelVar_i[id_z] - face.z;
					double nx = -faceVar_i[id_nx];
					double ny = -faceVar_i[id_ny];
					double nz = -faceVar_i[id_nz];
					double Lcond = x_pF*nx + y_pF*ny + z_pF*nz;
					if(Lcond <= eps){
						parcel.setType(type_p);
					}
				}
			}
		}
	}
	
	
}
	
void MASCH_Solver::updateProcRightParcels(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	if(size==1) return;
	
	double eps = -1.e-16;
	// double eps = -1.e-4;
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	
	int id_x = controls.getId_parcelVar("x-position");
	int id_y = controls.getId_parcelVar("y-position");
	int id_z = controls.getId_parcelVar("z-position");
	
	auto parcels_i = mesh.parcels.data();
	auto cell_i = mesh.cells.data();
	auto face_i = mesh.faces.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto parcelVar = var.parcels.data();
	int proc_size = var.procRightCells.size();
	int send_total_values = controls.nParcelVar;
	
	vector<int> send_each_count_parcels;
	send_each_count_parcels.reserve(proc_size);
	vector<double> send_parcel_values;
	vector<vector<double>> inp_send_value(size);
	int tmp_nToProcsRishtParcels = 0;
	for(auto& boundary : mesh.boundaries){
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			int rightProcNo = boundary.rightProcNo;
			for(int i=str; i<end; ++i){
				auto& face = face_i[i];
				auto faceVar_i = faceVar[i].data();
				int iL = face.iL;
				auto& cell = cell_i[iL];
				
				int tmp_num = 0;
				for(auto& iparcel : cell.iparcels){
					auto& parcel = parcels_i[iparcel];
					auto parcelVar_i = parcelVar[iparcel].data();
					if(parcel.getType() != MASCH_Parcel_Types::INSIDE) continue;
					double x_pF = parcelVar_i[id_x] - face.x;
					double y_pF = parcelVar_i[id_y] - face.y;
					double z_pF = parcelVar_i[id_z] - face.z;
					double nx = -faceVar_i[id_nx];
					double ny = -faceVar_i[id_ny];
					double nz = -faceVar_i[id_nz];
					double Lcond = x_pF*nx + y_pF*ny + z_pF*nz;
					if(Lcond < eps){
						parcel.setType(MASCH_Parcel_Types::TO_PROCS_RIGHT);
						++tmp_nToProcsRishtParcels;
						for(int j=0; j<send_total_values; ++j){
							send_parcel_values.push_back(parcelVar_i[j]);
							inp_send_value[rightProcNo].push_back(parcelVar_i[j]);
						}
						++tmp_num;
					}
				}
				send_each_count_parcels.push_back(tmp_num);
			}
		}
	}
	
	int nToProcsRishtParcels_glo = tmp_nToProcsRishtParcels;
	if(size>1){
		MPI_Allreduce(&tmp_nToProcsRishtParcels, &nToProcsRishtParcels_glo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	controls.nToProcsRishtParcels = nToProcsRishtParcels_glo;
	
	MASCH_MPI mpi;
	
	vector<int> displs;
	vector<double> recv_parcel_values;
	mpi.Alltoallv(inp_send_value, recv_parcel_values, displs);
	
	vector<int> recv_each_count_parcels;
	mpi.procFace_Alltoallv(
		send_each_count_parcels, recv_each_count_parcels,
		mesh.countsSendProcFaces, mesh.countsRecvProcFaces, 
		mesh.displsSendProcFaces, mesh.displsRecvProcFaces);

	int tmp_parcel_size = mesh.parcels.size();
	
	for(int i=0, iter=0, SIZE=recv_parcel_values.size()/send_total_values; i<SIZE; ++i){
		mesh.parcels.push_back(MASCH_Parcel());
		var.parcels.push_back(vector<double>());
		var.parcels.back().resize(send_total_values);
		
		for(int j=0; j<send_total_values; ++j){
			var.parcels.back()[j] = recv_parcel_values[iter++];
		}
	}
	
	int ip = 0;
	for(auto& boundary : mesh.boundaries){
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		if(boundary.getType() == MASCH_Face_Types::PROCESSOR){
			int rightProcNo = boundary.rightProcNo;
			for(int i=str; i<end; ++i){
				auto& face = face_i[i];
				auto faceVar_i = faceVar[i].data();
				int iL = face.iL;
				auto& cell = cell_i[iL];
				
				int count_parcel = recv_each_count_parcels[ip++];
				
				for(int j=0; j<count_parcel; ++j){
					mesh.parcels[tmp_parcel_size].setType(MASCH_Parcel_Types::INSIDE);
					mesh.parcels[tmp_parcel_size++].icell = iL;
				}
			}
		}
	}
	
}
	
	
void MASCH_Solver::refreshParcels(
MASCH_Mesh& mesh, 
MASCH_Control& controls,
MASCH_Variables& var){
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto parcelVar = var.parcels.data();
	auto cell_i = mesh.cells.data();
	auto parcel_i = mesh.parcels.data();
    
	int id_rho = controls.getId_parcelVar("density");
	int id_u = controls.getId_parcelVar("x-velocity");
	int id_v = controls.getId_parcelVar("y-velocity");
	int id_w = controls.getId_parcelVar("z-velocity");
	int id_T = controls.getId_parcelVar("temperature");
	int id_nparcel = controls.getId_parcelVar("number-of-parcel");
	int id_d = controls.getId_parcelVar("diameter");
	int id_x = controls.getId_parcelVar("x-position");
	int id_y = controls.getId_parcelVar("y-position");
	int id_z = controls.getId_parcelVar("z-position");
    
    
    double maxNParcel = stod(controls.controlParcelsMap["maxNumberOfParcel"]);
    double minDia = stod(controls.controlParcelsMap["minDiameter"]);
	
	int tmp_nInsideParcels=0;
	int tmp_nReflectParcels=0;
	int tmp_nEscapeParcels=0;
	int tmp_nDeleteParcels=0;
	for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i){
		auto& parcel = parcel_i[i];
        auto parcelVar_i = parcelVar[i];
        
        double N_p = parcelVar_i[id_nparcel];
        double d_p = parcelVar_i[id_d];
        if(N_p>maxNParcel)  parcel.setType(MASCH_Parcel_Types::TO_BE_DELETE);
        if(d_p<minDia)  parcel.setType(MASCH_Parcel_Types::TO_BE_DELETE);
        
		if(parcel.getType() == MASCH_Parcel_Types::INSIDE){
			++tmp_nInsideParcels;
		}
		else if(parcel.getType() == MASCH_Parcel_Types::REFLECT){
			++tmp_nReflectParcels;
		}
		else if(parcel.getType() == MASCH_Parcel_Types::ESCAPE){
			++tmp_nEscapeParcels;
		}
		else if(parcel.getType() == MASCH_Parcel_Types::TO_BE_DELETE){
			++tmp_nDeleteParcels;
		}
	}
	int nInsideParcels_glo = tmp_nInsideParcels;
	int nReflectParcels_glo = tmp_nReflectParcels;
	int nEscapeParcels_glo = tmp_nEscapeParcels;
	int nDeleteParcels_glo = tmp_nDeleteParcels;
	if(size>1){
		MPI_Allreduce(&tmp_nInsideParcels, &nInsideParcels_glo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tmp_nReflectParcels, &nReflectParcels_glo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tmp_nEscapeParcels, &nEscapeParcels_glo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&tmp_nDeleteParcels, &nDeleteParcels_glo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	controls.nInsideParcels = nInsideParcels_glo;
	controls.nReflectParcels = nReflectParcels_glo;
	controls.nEscapeParcels = nEscapeParcels_glo;
	controls.nDeleteParcels = nDeleteParcels_glo;
	
	// force eavporation
	{
		int num=0;
		auto& vars = var.parcels;
		var.parcels.erase( std::remove_if( var.parcels.begin(), var.parcels.end(), 
			[&mesh,&num](auto const& parcel) { 
			return (mesh.parcels[num++].getType() != MASCH_Parcel_Types::INSIDE);
			}), var.parcels.end());
	}
	{
		int num=0;
		auto& vars = mesh.parcels;
		mesh.parcels.erase( std::remove_if( mesh.parcels.begin(), mesh.parcels.end(), 
			[&mesh,&num](MASCH_Parcel const& parcel) { 
			return (mesh.parcels[num++].getType() != MASCH_Parcel_Types::INSIDE);
			}), mesh.parcels.end());
	}
	
	// cell to iparcel 컨넥션
	for(auto& cell : mesh.cells){
		cell.iparcels.clear();
	}
	for(int i=0, SIZE=mesh.parcels.size(); i<SIZE; ++i){
		auto& parcel = parcel_i[i];
		// auto parcelVar_i = parcelVar[i].data();
		
		auto& cell = cell_i[parcel.icell];
		if(find(cell.iparcels.begin(),cell.iparcels.end(),i)==cell.iparcels.end()){
			cell.iparcels.push_back(i);
		}
	}
	
}
	
	
// int MASCH_Solver::get_icell(
// MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var,
// double parcel_x, double parcel_y, double parcel_z){
	
	// // int id_nx = controls.getId_faceVar("x unit normal");
	// // int id_ny = controls.getId_faceVar("y unit normal");
	// // int id_nz = controls.getId_faceVar("z unit normal");
	
	// // auto cell_i = mesh.cells.data();
	// // auto face_i = mesh.faces.data();
	// // auto faceVar = var.faces.data();
	
	// // for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// // auto& cell = mesh.cells[i];
		
		// // // bool boolInside = checkInsideCell(cell, var.faces.data(), i, x, y, z);
		// // bool boolInside = true;
		// // for(auto& iface : cell.ifaces){
			// // auto& face = mesh.faces[iface];
			// // auto faceVar_i = faceVar[iface].data();
			// // double x_pF = parcel_x - face.x;
			// // double y_pF = parcel_y - face.y;
			// // double z_pF = parcel_z - face.z;
			// // double nx = -faceVar_i[id_nx];
			// // double ny = -faceVar_i[id_ny];
			// // double nz = -faceVar_i[id_nz];
			// // double Lcond = x_pF*nx + y_pF*ny + z_pF*nz;
			// // if(Lcond < 0.0){
				// // boolInside = false;
				// // break;
			// // }
		// // }
		// // if(boolInside==false) continue;
		
		// // return i;
		
	// // }
	
	// // return -1;
	
	
// }


