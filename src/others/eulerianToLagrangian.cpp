
#include "./solvers.h"

void MASCH_Solver::eulerianToLagrangian(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto& solver = (*this);
		
	double eps = -1.e-16;
	
	
	int id_par_rho = controls.getId_parcelVar("density");
	int id_par_u = controls.getId_parcelVar("x-velocity");
	int id_par_v = controls.getId_parcelVar("y-velocity");
	int id_par_w = controls.getId_parcelVar("z-velocity");
	int id_par_T = controls.getId_parcelVar("temperature");
	int id_par_nparcel = controls.getId_parcelVar("number-of-parcel");
	int id_par_d = controls.getId_parcelVar("diameter");
	int id_par_x = controls.getId_parcelVar("x-position");
	int id_par_y = controls.getId_parcelVar("y-position");
	int id_par_z = controls.getId_parcelVar("z-position");
	int id_vol = controls.getId_cellVar("volume");
	int id_rho = controls.getId_cellVar("density");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	auto idSetLagrangianEulerian = controls.idSetLagrangianEulerian;


	for(int id=0; id<controls.nameParcels.size(); ++id){
		string name = controls.nameParcels[id];
			
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
		double alpha_cri = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.criterionRatioRadius"]);
		double dx_UDV = stod(controls.controlParcelsMap[name+".eulerianToLagrangian.minRadius"]);
			
		// parcel 그룹 찾기 (MPI 사용 안함. 블럭 내부에서만 변경)
		int id_Y = controls.getId_cellVar("mass-fraction-"+name);
		int id_alpha = controls.getId_cellVar("volume-fraction-"+name);
		
		double rho_drop = stod(controls.controlParcelsMap[name+".thermodynamics.rho"]);
		
		vector<bool> choicedCellGroups(mesh.cells.size(),false);
		vector<int> id_parcel(mesh.cells.size(),-1);
		vector<vector<int>> parcelGroups_icells;
		
		int n_parcel = 0;
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			// if(choicedCellGroups[i]==true) continue;
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			
			if(cellVar_i[id_alpha]>=eps_max){
				
				vector<int> tmp_parcelGroup_id;
				for(auto& icell : cell.iStencils){
					if(choicedCellGroups.at(icell)==true){
						tmp_parcelGroup_id.push_back(id_parcel.at(icell));
						break;
					}
				}
				
				// 새로 id 추가
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
		
		
		// 재 합침
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
		
		
		// proc 면하고 접촉하고 있는 그룹은 모두 취소 (MPI 사용 안할거기 떄문에)
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
		
		
		
		// 각 셀의 , 타겟 id 의 부피, 질량 구하기
		vector<double> target_volume_inCell(mesh.cells.size());
		vector<double> target_mass_inCell(mesh.cells.size());
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = cells[i];
			auto cellVar_i = cellVar[i].data();
			double tmp_rho = cellVar_i[id_rho];
			double tmp_Y = cellVar_i[id_Y];
			double tmp_alpha = cellVar_i[id_alpha];
			double cell_vol = cellVar_i[id_vol];
			
			target_volume_inCell[i] = tmp_alpha*cell_vol;
			target_mass_inCell[i] = tmp_rho*tmp_Y*cell_vol;
		}
		
		
		
		
		// parcel 의 값들.. 질량 평균값 넣기
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
					double target_mass = target_mass_inCell[icell];
					total_volume += target_volume;
					total_mass += target_mass;
					
					// 위치 평균
					avg_values[id_par_x] += cellSten.x * target_mass;
					avg_values[id_par_y] += cellSten.y * target_mass;
					avg_values[id_par_z] += cellSten.z * target_mass;
					
					// 추가적인 값 평균
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
		
		
		// 변환 기준 계산
		{
			vector<bool> boolCanChangeE2P(parcelGroups_icells.size(),false);
			for(int i=0, SIZE=parcelGroups_icells.size(); i<SIZE; ++i){
				auto& icells = parcelGroups_icells[i];
				if(icells.size()==0) continue;
				
				double total_volume = groupParcelTotalVolume[i];
				
				double mass_center_x = groupParcelVars[i][id_par_x];
				double mass_center_y = groupParcelVars[i][id_par_y];
				double mass_center_z = groupParcelVars[i][id_par_z];
				
				// 최대 거리 계산
				double R_max = -1.e8;
				for(auto& icell : icells){
					auto& cellSten = cells[icell];
					auto cellStenVar_i = cellVar[icell].data();
					
					double dx = cellSten.x-mass_center_x;
					double dy = cellSten.y-mass_center_y;
					double dz = cellSten.z-mass_center_z;
					double dist = sqrt(dx*dx+dy*dy+dz*dz);
					R_max = max(dist,R_max);
				}
				
				double tmp_R = pow(3.0/4.0/3.141592*total_volume,0.3333);
				groupParcelVars[i][id_par_d] = 2.0*tmp_R;
				
				groupParcelVars[i][id_par_nparcel] = 1.0;
				groupParcelVars[i][id_par_rho] = rho_drop;
				
					// cout << scientific; cout.precision(2);
					// cout << groupParcelVars[i][id_par_x] << " " << groupParcelVars[i][id_par_y] << " " << groupParcelVars[i][id_par_z] << endl;
					// cout << total_volume << " " << R_cri << " " << R_max << endl;
					// cout << fixed; cout.precision(0);
				// if(
				// total_volume <= 4.0/3.0*3.141592*R_cri*R_cri*R_cri &&
				// R_max/(max(dx_UDV,tmp_R) <= alpha_cri)
				// ){
				if(
				(total_volume <= 4.0/3.0*3.141592*R_cri*R_cri*R_cri) &&
				(R_max/(max(dx_UDV,tmp_R)) <= alpha_cri)
				){
					boolCanChangeE2P[i] = true;
				}
				
				
			}
			int numN = 0;
			groupParcelVars.erase( std::remove_if( groupParcelVars.begin(), groupParcelVars.end(), 
				[&boolCanChangeE2P, &numN](vector<double> const& v) { 
				return !boolCanChangeE2P[numN++]; 
				}), groupParcelVars.end());
			numN = 0;
			parcelGroups_icells.erase( std::remove_if( parcelGroups_icells.begin(), parcelGroups_icells.end(), 
				[&boolCanChangeE2P, &numN](vector<int> const& v) { 
				return !boolCanChangeE2P[numN++]; 
				}), parcelGroups_icells.end());
		}
		
		int tmp_nChangeParcels = groupParcelVars.size();
		int nChangeParcels_glb = tmp_nChangeParcels;
		if(size>1){
			MPI_Allreduce(&tmp_nChangeParcels, &nChangeParcels_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
		controls.nChangeParcelsE2L = nChangeParcels_glb;
		// MPI_Barrier(MPI_COMM_WORLD);
		// if(rank==0){
			// // cout << "-> completed" << endl;
			// cout << 
			// "| parcel(EtoL) = +" << nChangeParcels_glb << endl;
			// // " | face = +" << addedFaceSize_glb <<
			// // " | point = +" << addedPointSize_glb << endl;
			// // cout << "└────────────────────────────────────────────────────" << endl;
		// }
		
		
	
		
		// 변환
		{
			for(int iPar=0, SIZE=groupParcelVars.size(); iPar<SIZE; ++iPar){
				auto& vars = groupParcelVars[iPar];
				auto& icells = parcelGroups_icells[iPar];
				
				// parcel 이 있는 셀 아이디 찾기
				double pres_x = vars[id_par_x];
				double pres_y = vars[id_par_y];
				double pres_z = vars[id_par_z];
				int pres_icell = -1;
				for(auto& i : icells){
					auto& cell = mesh.cells[i];
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
				
				if(pres_icell==-1) {
					cout << scientific; cout.precision(2);
					cout << "#WARNING : not find parcel location" << endl;
					cout << pres_x << " " << pres_y << " " << pres_z << endl;
					cout << icells.size() << endl;
					cout << fixed; cout.precision(0);
				}
				
				// parcel 생성
				mesh.addParcel(pres_icell, id);
				var.addParcel(controls.nParcelVar, vars);
				
				// euler 값 지우기
				for(auto& icell : icells){
					auto& cellSten = cells[icell];
					auto cellStenVar_i = cellVar[icell].data();
					
					cellStenVar_i[id_Y] = 0.0;
					
				}
			}
		}
		
	}
	
	
	// if(controls.nChangeParcelsE2L!=0){
		
		// cout << controls.nChangeParcelsE2L << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	// solver.searchParcelGroups(mesh, controls, var);
	
	// 그룹을 parcel 로 바꿔줌
	// solver.changeToParcels(mesh, var);
	
}



// void MASCH_Solver::searchParcelGroups(
// MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	
	
	
	
	
// }



// void MASCH_Solver::changeToParcels(
// MASCH_Mesh& mesh, MASCH_Variables& var){
	
	
// }