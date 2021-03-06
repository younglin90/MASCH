
#include "./solvers.h"

void MASCH_Solver::eulerianToLagrangian(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& solver = (*this);
		
	double eps = -1.e-4;
	// double eps = -1.e-16;
	
	
	int id_par_rho = controls.getId_parcelVar("density");
	int id_par_u = controls.getId_parcelVar("x-velocity");
	int id_par_v = controls.getId_parcelVar("y-velocity");
	int id_par_w = controls.getId_parcelVar("z-velocity");
	int id_par_T = controls.getId_parcelVar("temperature");
	int id_par_nparcel = controls.getId_parcelVar("number-of-parcel");
	int id_par_time = controls.getId_parcelVar("time");
	int id_par_d = controls.getId_parcelVar("diameter");
	int id_par_x = controls.getId_parcelVar("x-position");
	int id_par_y = controls.getId_parcelVar("y-position");
	int id_par_z = controls.getId_parcelVar("z-position");
	int id_vol = controls.getId_cellVar("volume");
	int id_rho = controls.getId_cellVar("density");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_p = controls.getId_cellVar("pressure");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
    
    int id_c = controls.getId_cellVar("speed-of-sound");
    int id_Ht = controls.getId_cellVar("total-enthalpy");
    int id_drhodp = controls.getId_cellVar("partial-density-pressure");
    int id_drhodT = controls.getId_cellVar("partial-density-temperature");
    int id_dHtdp = controls.getId_cellVar("partial-total-enthalpy-pressure");
    int id_dHtdT = controls.getId_cellVar("partial-total-enthalpy-temperature");
    
    int nSp = controls.spName.size();
    int id_Y[nSp], id_drhodY[nSp], id_dHtdY[nSp];
    for(int i=0; i<nSp-1; ++i){
        id_Y[i] = controls.getId_cellVar("mass-fraction-"+controls.spName[i]);
        id_drhodY[i] = controls.getId_cellVar("partial-density-mass-fraction-"+controls.spName[i]);
        id_dHtdY[i] = controls.getId_cellVar("partial-total-enthalpy-mass-fraction-"+controls.spName[i]);
    }
	
	double dt = var.fields[controls.getId_fieldVar("time-step")];
	
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	auto idSetLagrangianEulerian = controls.idSetLagrangianEulerian;


    bool boolDebug = false;
    
    
    


	for(int id=0; id<controls.nameParcels.size(); ++id){
        
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START1" << endl;
    
		string name = controls.nameParcels[id];
		
		int iSp = -1;
		for(int i=0; i<controls.spName.size(); ++i){
			if(name == controls.spName[i]){
				iSp = i;
				break;
			}
		}
		if(iSp==-1) cout << "#WARNING : can not fined spname, " << name << endl;
		
		if(controls.controlParcelsMap.find(name+".eulerianToLagrangian.interval")==
		controls.controlParcelsMap.end()) {
			cout << "#WARNING : not defined, " << name << ".eulerianToLagrangian.interval" << endl;
			continue;
		}
		if( (controls.iterReal+1) %
			stoi(controls.controlParcelsMap[name+".eulerianToLagrangian.interval"]) != 0){
			continue;	
		}
			
		if(controls.controlParcelsMap.find(name+".eulerianToLagrangian.type")==
		controls.controlParcelsMap.end()) {
			cout << "#WARNING : not defined, " << name << ".eulerianToLagrangian.type" << endl;
			continue;
		}
		if(controls.controlParcelsMap[name+".eulerianToLagrangian.type"] == "none"){
			continue;
		}
		
		double eps_max = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.criterionVolumeFraction"]);
		double R_cri = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.criterionRadius"]);
		double m_cri = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.criterionRatioSize"]);
		double alpha_cri = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.criterionRatioRadius"]);
		double dx_UDV = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.minRadius"]);
			
		// parcel ?????? ?????? (MPI ?????? ??????. ?????? ??????????????? ??????)
		// int id_Ht = controls.getId_cellVar("total-enthalpy");
		int id_Y_target = controls.getId_cellVar("mass-fraction-"+name);
		int id_alpha_target = controls.getId_cellVar("volume-fraction-"+name);
		
		double rho_drop = stod(controls.controlParcelsMap[name+".thermodynamics.rho"]);
        
        
        
        
        string sRegionType;
        vector<double> boxMin(3,-1.e12), boxMax(3,1.e12);
        if(controls.controlParcelsMap.find(name+".eulerianToLagrangian.region.type")==controls.controlParcelsMap.end()){
            sRegionType = "box";
            if(rank==0) cout << "#WARNING : no seraching, eulerianToLagrangian.region.type" << endl;
        }
        else{
            if(controls.controlParcelsMap[name+".eulerianToLagrangian.region.type"]=="box"){
                MASCH_Load load;
                vector<string> sboxMin = load.extractVector(controls.controlParcelsMap[name+".eulerianToLagrangian.region.min"]);
                vector<string> sboxMax = load.extractVector(controls.controlParcelsMap[name+".eulerianToLagrangian.region.max"]);
                boxMin.clear(); boxMax.clear();
                for(auto& item : sboxMin) boxMin.push_back(stod(item));
                for(auto& item : sboxMax) boxMax.push_back(stod(item));
            }
            else{
                if(rank==0) cout << "#WARNING : no seraching, eulerianToLagrangian.region.type" << endl;
            }
        }
        
        
        int max_n_E2L = 10000000;
        if(controls.controlParcelsMap.find(name+".eulerianToLagrangian.maxNumberPerIter")==controls.controlParcelsMap.end()){
            if(rank==0) cout << "#WARNING : no seraching, eulerianToLagrangian.maxNumberPerIter" << endl;
        }
        else{
            max_n_E2L = stoi(controls.controlParcelsMap[name+".eulerianToLagrangian.maxNumberPerIter"]);
        }
        
        
		
		vector<bool> choicedCellGroups(mesh.cells.size(),false);
		vector<int> id_parcel(mesh.cells.size(),-1);
		vector<vector<int>> parcelGroups_icells;
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START2" << endl;
    
		int n_parcel = 0;
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// if(choicedCellGroups[i]==true) continue;
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			
			if(cellVar_i[id_alpha_target]>=eps_max){
				
				vector<int> tmp_parcelGroup_id;
				for(auto& icell : cell.iStencils){
					if(choicedCellGroups.at(icell)==true){
						tmp_parcelGroup_id.push_back(id_parcel.at(icell));
						break;
					}
				}
				
				// ?????? id ??????
				if(tmp_parcelGroup_id.size()==0){
					parcelGroups_icells.push_back(vector<int>());
					parcelGroups_icells.back().push_back(i);
					id_parcel.at(i) = n_parcel++;
				}
				else{
					int tmp_id = tmp_parcelGroup_id.at(0);
					// cout << tmp_id << endl;
					parcelGroups_icells.at(tmp_id).push_back(i);
					id_parcel.at(i) = tmp_id;
				}
				choicedCellGroups[i]=true;
			
			}
		}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START3" << endl;
		
		// ??? ??????
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			int i_id = id_parcel[i];
			if(i_id==-1) continue;
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			
			for(auto& icell : cell.iStencils){
				auto& cellSten = cells[icell];
				auto cellStenVar_i = cellVar[icell].data();
				int iRight_id = id_parcel[icell];
				if(i_id!=iRight_id && iRight_id!=-1){
					for(auto& icell_tar : parcelGroups_icells[iRight_id]){
						id_parcel[icell_tar] = i_id;
					}
					parcelGroups_icells[i_id].insert(parcelGroups_icells[i_id].end(),
					parcelGroups_icells[iRight_id].begin(),parcelGroups_icells[iRight_id].end());
					parcelGroups_icells[iRight_id].clear();
				}
			}
		}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START4" << endl;
		
		// proc ????????? ???????????? ?????? ????????? ?????? ?????? (MPI ?????? ???????????? ?????????)
		for(auto& boundary : mesh.boundaries){
			// if(boundary.getType()!=MASCH_Face_Types::PROCESSOR) continue;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				int iL = mesh.faces[i].iL;
				int i_id = id_parcel[iL];
				if(i_id==-1) continue;
				parcelGroups_icells[i_id].clear();
			}
		}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START5" << endl;
		
		
		// ??? ?????? , ?????? id ??? ??????, ?????? ?????????
		vector<double> target_volume_inCell(mesh.cells.size());
		vector<double> target_mass_inCell(mesh.cells.size());
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			double tmp_rho = cellVar_i[id_rho];
			double tmp_Y = cellVar_i[id_Y_target];
			double tmp_alpha = cellVar_i[id_alpha_target];
			double cell_vol = cellVar_i[id_vol];
			
			target_volume_inCell[i] = cell_vol;//tmp_alpha*cell_vol;
			target_mass_inCell[i] = tmp_rho*tmp_Y*cell_vol;
		}
		
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START6" << endl;
		
		
		// parcel ??? ??????.. ?????? ????????? ??????
		int nParcelVar = controls.nParcelVar;
		vector<vector<double>> groupParcelVars(parcelGroups_icells.size(),vector<double>(nParcelVar));
		vector<double> groupParcelTotalVolume(parcelGroups_icells.size());
		{
			int iter = 0;
			for(auto& icells : parcelGroups_icells){
				if(icells.size()==0) {
					++iter;
					continue;
				}
				vector<double> avg_values(nParcelVar,0.0); 
				double total_mass = 0.0;
				double total_volume = 0.0;
				for(auto& icell : icells){
					auto& cellSten = cells[icell];
					auto cellStenVar_i = cellVar[icell].data();
					
					double target_volume = target_volume_inCell[icell];
					double target_mass = target_volume;//target_mass_inCell[icell];
					
					total_volume += target_volume;
					total_mass += target_mass;
					
					// ?????? ??????
					avg_values[id_par_x] += cellSten.x * target_mass;
					avg_values[id_par_y] += cellSten.y * target_mass;
					avg_values[id_par_z] += cellSten.z * target_mass;
					
					// ???????????? ??? ??????
					for(auto& [id_lag, id_eul] : idSetLagrangianEulerian){
						// int id_eu = id_matchEuler[j++];
						avg_values[id_lag] += cellStenVar_i[id_eul] * target_mass;
					}
					
				}
				int iter2=0;
				for(auto& avg_val : avg_values) {
					avg_val /= total_mass;
					groupParcelVars[iter][iter2++] = avg_val;
				}
				// cout << scientific; cout.precision(2);
				// cout << avg_values[id_par_x] << " " << avg_values[id_par_y] << " " << avg_values[id_par_z] << endl;
				// cout << fixed; cout.precision(0);
				groupParcelTotalVolume[iter] = total_volume;
				++iter;
			}
		}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START7" << endl;
		
		// ?????? ?????? ??????
		{
			vector<bool> boolCanChangeE2P(parcelGroups_icells.size(),false);
			for(int i=0, SIZE=parcelGroups_icells.size(); i<SIZE; ++i){
				auto& icells = parcelGroups_icells[i];
				if(icells.size()==0) continue;
				
				double total_volume = groupParcelTotalVolume[i];
				
				double mass_center_x = groupParcelVars[i][id_par_x];
				double mass_center_y = groupParcelVars[i][id_par_y];
				double mass_center_z = groupParcelVars[i][id_par_z];
				
				// ?????? ?????? ??????
				double R_max = -1.e8;
                double tmp_min_vol = 1.e12;
				for(auto& icell : icells){
					auto& cellSten = cells[icell];
					auto cellStenVar_i = cellVar[icell].data();
					
					double dx = cellSten.x-mass_center_x;
					double dy = cellSten.y-mass_center_y;
					double dz = cellSten.z-mass_center_z;
					double dist = sqrt(dx*dx+dy*dy+dz*dz);
					R_max = max(dist,R_max);
					tmp_min_vol = min(tmp_min_vol,cellStenVar_i[id_vol]);
				}
				
				double tmp_R = pow(3.0/4.0/3.141592*total_volume,1.0/3.0);
				groupParcelVars[i][id_par_d] = 2.0*tmp_R;
				
				groupParcelVars[i][id_par_nparcel] = 1.0;
				groupParcelVars[i][id_par_time] = 0.0;
				groupParcelVars[i][id_par_rho] = rho_drop;
				
					// cout << scientific; cout.precision(2);
					// cout << groupParcelVars[i][id_par_x] << " " << groupParcelVars[i][id_par_y] << " " << groupParcelVars[i][id_par_z] << endl;
					// cout << total_volume << " " << R_cri << " " << R_max << endl;
					// cout << fixed; cout.precision(0);
				// if(
				// total_volume <= 4.0/3.0*3.141592*R_cri*R_cri*R_cri &&
				// R_max/(max(dx_UDV,tmp_R) <= alpha_cri)
				// ){
				// cout << scientific; cout.precision(2);
				// cout << total_volume << " " << 4.0/3.0*3.141592*R_cri*R_cri*R_cri << endl;
				// cout << R_max << " " << tmp_R << " " << alpha_cri << endl;
				// cout << fixed; cout.precision(0);
                
                double R_cri_ratio = min( R_cri, m_cri*pow(tmp_min_vol,1.0/3.0) );
                double R_cri_ratio3 = R_cri_ratio*R_cri_ratio*R_cri_ratio;
                double vol_cri = 4.0/3.0*3.141592*R_cri_ratio3;
				if(
				(total_volume <= vol_cri) &&
				(R_max/(max(dx_UDV,tmp_R)) <= alpha_cri) 
				){
					boolCanChangeE2P[i] = true;
				}
				
				
			}
			int numN = 0;
            // if(groupParcelVars.size()>0){
                groupParcelVars.erase( std::remove_if( groupParcelVars.begin(), groupParcelVars.end(), 
                    [&boolCanChangeE2P, &numN](vector<double> const& v) { 
                    return !boolCanChangeE2P[numN++]; 
                    }), groupParcelVars.end());
            // }
			numN = 0;
            // if(parcelGroups_icells.size()>0){
                parcelGroups_icells.erase( std::remove_if( parcelGroups_icells.begin(), parcelGroups_icells.end(), 
                    [&boolCanChangeE2P, &numN](vector<int> const& v) { 
                    return !boolCanChangeE2P[numN++]; 
                    }), parcelGroups_icells.end());
            // }
		}
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START8" << endl;
    
		// int tmp_nChangeParcels = groupParcelVars.size();
		// int nChangeParcels_glb = tmp_nChangeParcels;
		// if(size>1){
			// MPI_Allreduce(&tmp_nChangeParcels, &nChangeParcels_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		// }
		// controls.nChangeParcelsE2L = nChangeParcels_glb;
		// // MPI_Barrier(MPI_COMM_WORLD);
		// // if(rank==0){
			// // // cout << "-> completed" << endl;
			// // cout << 
			// // "| parcel(EtoL) = +" << nChangeParcels_glb << endl;
			// // // " | face = +" << addedFaceSize_glb <<
			// // // " | point = +" << addedPointSize_glb << endl;
			// // // cout << "???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????" << endl;
		// // }
		
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START9" << endl;
	
    
		// int tmp_nChangeParcels = 0;
		
		// ??????
		{
			int iSegEq = 0;
			int nEq = controls.nEq[iSegEq];
			
            int numE2L = 0;
            
			for(int iPar=0, SIZE=groupParcelVars.size(); iPar<SIZE; ++iPar){
				auto& vars = groupParcelVars[iPar];
				auto& icells = parcelGroups_icells[iPar];
                
                
                if(numE2L > max_n_E2L) break;
                
                
                bool boolRegion = true;
				
				// parcel ??? ?????? ??? ????????? ??????
				double pres_x = vars[id_par_x];
				double pres_y = vars[id_par_y];
				double pres_z = vars[id_par_z];
				int pres_icell = -1;
				for(auto& i : icells){
					auto& cell = mesh.cells[i];
                    
                    

                    if(
                    cell.x<boxMin[0] || cell.x>boxMax[0] ||
                    cell.y<boxMin[1] || cell.y>boxMax[1] ||
                    cell.z<boxMin[2] || cell.z>boxMax[2]
                    ) {
                        boolRegion = false;
                        break;
                    }
                    
                    
					bool boolInside = true;
					for(auto& iface : cell.ifaces){
						auto& face = mesh.faces[iface];
						auto faceVar_i = faceVar[iface].data();
						double x_pF = pres_x - face.x;
						double y_pF = pres_y - face.y;
						double z_pF = pres_z - face.z;
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
					
					pres_icell = i;
					break;
				}
                
                
                if(boolRegion==false) continue;
                
				
				if(pres_icell==-1) {
					cout << scientific; cout.precision(8);
					cout << "#WARNING : not find parcel location -> deleted" << endl;
					continue;
					cout << pres_x << " " << pres_y << " " << pres_z << endl;
					cout << icells.size() << endl;
					int iter3 = 0;
					for(auto& i : icells){
						auto& cell = mesh.cells[i];
						cout << iter3 << " | " << cell.x << " " << cell.y << " " << cell.z << endl;
						++iter3;
					}
					cout << fixed; cout.precision(0);
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				// parcel ??????
                // vars[id_N_p] = 1.0;
                // vars[id_time_p] = 0.0;
                
                // cout.precision(6);
                // cout << vars[id_par_d] << " " << vars[id_par_nparcel] << " " <<  pres_icell << endl;
                
                // ++tmp_nChangeParcels;
                
				mesh.addParcel(pres_icell, id);
				var.addParcel(controls.nParcelVar, vars);
				
				// euler ??? ?????????
				for(auto& icell : icells){
					auto& cellSten = cells[icell];
					auto cellStenVar_i = cellVar[icell].data();
					
					double vol = cellStenVar_i[id_vol];
					double Y_org = cellStenVar_i[id_Y_target];
                    
                    double p = cellStenVar_i[id_p];
                    double u = cellStenVar_i[id_u];
                    double v = cellStenVar_i[id_v]; 
                    double w = cellStenVar_i[id_w];
                    double rho = cellStenVar_i[id_rho];
                    double Ht = cellStenVar_i[id_Ht];
                    double drhodp = cellStenVar_i[id_drhodp];
                    double drhodT = cellStenVar_i[id_drhodT];
                    double dHtdp = cellStenVar_i[id_dHtdp];
                    double dHtdT = cellStenVar_i[id_dHtdT];
                    double Y[nSp], drhodY[nSp], dHtdY[nSp];
                    for(int i=0; i<nSp-1; ++i){
                        Y[i] = cellStenVar_i[id_Y[i]];
                        drhodY[i] = cellStenVar_i[id_drhodY[i]];
                        dHtdY[i] = cellStenVar_i[id_dHtdY[i]];
                    }
					
					double fluxB[nEq];
					double fluxA[nEq*nEq];
					for(int iEq=0; iEq<nEq; ++iEq){
                        for(int jEq=0; jEq<nEq; ++jEq){
                            fluxA[iEq*nEq+jEq] = 0.0;
                        }
						fluxB[iEq] = 0.0;
					}
                    
                    
                    double Y_org_diss = (Y_org-1.e-8);
                    // double Y_org_diss = Y_org;
                    
                    
					fluxB[0] -= Y_org_diss*rho/dt*vol;
					fluxB[1] -= Y_org_diss*rho*u/dt*vol;
					fluxB[2] -= Y_org_diss*rho*v/dt*vol;
					fluxB[3] -= Y_org_diss*rho*w/dt*vol;
					fluxB[4] -= Y_org_diss*(rho*Ht-p)/dt*vol;
					fluxB[5+iSp] -= Y_org_diss*rho/dt*vol;
                    
                    
                    
                    
                    // ????????? ????????????
                    int iter = 0;
                    fluxA[iter++] += Y_org_diss * drhodp * vol / dt;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += Y_org_diss * drhodT * vol / dt;
                    for(int j=0; j<nSp-1; ++j){
                        // if(iSp == j) fluxA[iter] += rho * vol / dt;
                        fluxA[iter++] += Y_org_diss * drhodY[j] * vol / dt;
                    }
                    
                    fluxA[iter++] += Y_org_diss * drhodp*u * vol / dt;
                    fluxA[iter++] += Y_org_diss * rho * vol / dt;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += Y_org_diss * drhodT*u * vol / dt;
                    for(int j=0; j<nSp-1; ++j){
                        // if(iSp == j) fluxA[iter] += rho*u * vol / dt;
                        fluxA[iter++] += Y_org_diss * drhodY[j]*u * vol / dt;
                    }
                    
                    fluxA[iter++] += Y_org_diss * drhodp*v * vol / dt;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += Y_org_diss * rho * vol / dt;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += Y_org_diss * drhodT*v * vol / dt;
                    for(int j=0; j<nSp-1; ++j){
                        // if(iSp == j) fluxA[iter] += rho*v * vol / dt;
                        fluxA[iter++] += Y_org_diss * drhodY[j]*v * vol / dt;
                    }
                    
                    fluxA[iter++] += Y_org_diss * drhodp*w * vol / dt;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += 0.0;
                    fluxA[iter++] += Y_org_diss * rho * vol / dt;
                    fluxA[iter++] += Y_org_diss * drhodT*w * vol / dt;
                    for(int j=0; j<nSp-1; ++j){
                        // if(iSp == j) fluxA[iter] += rho*w * vol / dt;
                        fluxA[iter++] += Y_org_diss * drhodY[j]*w * vol / dt;
                    }
                    
                    fluxA[iter++] += Y_org_diss * (drhodp*Ht+rho*dHtdp-1.0) * vol / dt;
                    fluxA[iter++] += Y_org_diss * rho*u * vol / dt;
                    fluxA[iter++] += Y_org_diss * rho*v * vol / dt;
                    fluxA[iter++] += Y_org_diss * rho*w * vol / dt;
                    fluxA[iter++] += Y_org_diss * (drhodT*Ht+rho*dHtdT) * vol / dt;
                    for(int j=0; j<nSp-1; ++j){
                        // if(iSp == j) fluxA[iter] += (rho*Ht-p) * vol / dt;
                        fluxA[iter++] += Y_org_diss * (drhodY[j]*Ht+rho*dHtdY[j]) * vol / dt;
                    }
                    
                    for(int i=0; i<nSp-1; ++i){
                        fluxA[iter++] += Y_org_diss * drhodp * vol / dt;
                        fluxA[iter++] += 0.0;
                        fluxA[iter++] += 0.0;
                        fluxA[iter++] += 0.0;
                        fluxA[iter++] += Y_org_diss * drhodT * vol / dt;
                        for(int j=0; j<nSp-1; ++j){
                            // if(iSp == j) fluxA[iter] += rho * vol / dt;
                            fluxA[iter++] += Y_org_diss * drhodY[j] * vol / dt;
                        }
                    }
                    
                    
					for(int iEq=0; iEq<nEq; ++iEq){
                        for(int jEq=0; jEq<nEq; ++jEq){
                            var.accumSparD( iSegEq, icell, iEq, jEq, fluxA[iEq*nEq+jEq] );
                        }
						var.accumB( iSegEq, icell, iEq, fluxB[iEq] );
					}
                    
                    
                    
                    
					
					// cellStenVar_i[id_Y_target] = 1.e-16;
					
				}
                
                
                ++numE2L;
                
			}
            
            
            
            int nChangeParcels_glb = numE2L;
            if(size>1){
                MPI_Allreduce(&numE2L, &nChangeParcels_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            }
            controls.nChangeParcelsE2L = nChangeParcels_glb;
		}
        
    
		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && boolDebug) cout << "START10" << endl;
    
    
    
		// MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0){
			// // cout << "-> completed" << endl;
			// cout << 
			// "| parcel(EtoL) = +" << nChangeParcels_glb << endl;
			// // " | face = +" << addedFaceSize_glb <<
			// // " | point = +" << addedPointSize_glb << endl;
			// // cout << "???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????" << endl;
		// }
	
	// if(controls.nChangeParcelsE2L!=0){
		
		// cout << controls.nChangeParcelsE2L << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	// solver.searchParcelGroups(mesh, controls, var);
	
	// ????????? parcel ??? ?????????
	// solver.changeToParcels(mesh, var);
	
}



// void MASCH_Solver::searchParcelGroups(
// MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	
	
	
	
	
// }



// void MASCH_Solver::changeToParcels(
// MASCH_Mesh& mesh, MASCH_Variables& var){
	
	
// }