
#include "./variables.h"

void MASCH_Variables::clearLinearSystems(int iSegEq){
	std::fill((*this).Avalues[iSegEq].begin(), (*this).Avalues[iSegEq].end(), 0.0);
	std::fill((*this).Xvalues[iSegEq].begin(), (*this).Xvalues[iSegEq].end(), 0.0);
	std::fill((*this).Bvalues[iSegEq].begin(), (*this).Bvalues[iSegEq].end(), 0.0);
}



void MASCH_Variables::accumSparD(int iSegEq, int i, int iEq, int jEq, double inp){
	int str = i_str_CSR[iSegEq][nCells*iEq + i];
	int tmp_i = str + (cellnFaces[i]+1)*jEq + cellcell_displ[i];
	(*this).Avalues.data()[iSegEq].data()[tmp_i] += inp;
	// (*this).Avalues[iSegEq][tmp_i] += inp;
	// if(tmp_i>=Avalues.at(iSegEq).size()) {
		// cout << "BB1" << endl;
		// cout << nCells << " " << i << " " << cellnFaces[i] << " " << cellcell_displ[i] << " " << str << endl;
	// }
	// (*this).Avalues.at(iSegEq).at(tmp_i) += inp;
	
}

double MASCH_Variables::getSparD(int iSegEq, int i, int iEq, int jEq){
	int str = i_str_CSR[iSegEq][nCells*iEq + i];
	int tmp_i = str + (cellnFaces[i]+1)*jEq + cellcell_displ[i];
	return (*this).Avalues.data()[iSegEq].data()[tmp_i];
	// return (*this).Avalues[iSegEq][tmp_i];
	// if(tmp_i>=Avalues.at(iSegEq).size()) {
		// cout << "BB2" << endl;
		// cout << nCells << " " << i << " " << cellnFaces[i] << " " << cellcell_displ[i] << " " << str << endl;
	// }
	// return (*this).Avalues.at(iSegEq).at(tmp_i);
}


void MASCH_Variables::accumSparLR(int iSegEq, int i, int iL, int iEq, int jEq, double inp){
	int str = i_str_CSR[iSegEq][nCells*iEq + iL];
	int tmp_i = str + (cellnFaces[iL]+1)*jEq + face_LR_displ[i];
	(*this).Avalues.data()[iSegEq].data()[tmp_i] += inp;
	// (*this).Avalues[iSegEq][tmp_i] += inp;
	// if(tmp_i>=Avalues.at(iSegEq).size()) {
		// cout << "BB3" << endl;
		// cout << nCells << " " << i << " " << cellnFaces[i] << " " << face_LR_displ[i] << " " << str << endl;
	// }
	// (*this).Avalues.at(iSegEq).at(tmp_i) += inp;
	
}

void MASCH_Variables::accumSparRL(int iSegEq, int i, int iR, int iEq, int jEq, double inp){
	int str = i_str_CSR[iSegEq][nCells*iEq + iR];
	int tmp_i = str + (cellnFaces[iR]+1)*jEq + face_RL_displ[i];
	(*this).Avalues.data()[iSegEq].data()[tmp_i] += inp;
	// (*this).Avalues[iSegEq][tmp_i] += inp;
	// if(tmp_i>=Avalues.at(iSegEq).size()) {
		// cout << "BB4" << endl;
		// cout << nCells << " " << i << " " << cellnFaces[i] << " " << face_RL_displ[i] << " " << str << endl;
	// }
	// (*this).Avalues.at(iSegEq).at(tmp_i) += inp;
	
}

void MASCH_Variables::accumB(int iSegEq, int i, int iEq, double inp){
	(*this).Bvalues.data()[iSegEq].data()[nCells*iEq+i] += inp;
	// (*this).Bvalues[iSegEq][nCells*iEq+i] += inp;
	// (*this).Bvalues.at(iSegEq).at(nCells*iEq+i) += inp;
}
double MASCH_Variables::getB(int iSegEq, int i, int iEq){
	return (*this).Bvalues.data()[iSegEq].data()[nCells*iEq+i];
	// return (*this).Bvalues[iSegEq][nCells*iEq+i];
	// return (*this).Bvalues.at(iSegEq).at(nCells*iEq+i);
	
}
void MASCH_Variables::setX(int iSegEq, int i, int iEq, double inp){
	(*this).Xvalues.data()[iSegEq].data()[nCells*iEq+i] += inp;
	// (*this).Xvalues[iSegEq][nCells*iEq+i] += inp;
	// (*this).Xvalues.at(iSegEq).at(nCells*iEq+i) += inp;
	
}
double MASCH_Variables::getX(int iSegEq, int i, int iEq){
	return (*this).Xvalues.data()[iSegEq].data()[nCells*iEq+i];
	// return (*this).Xvalues[iSegEq][nCells*iEq+i];
	// return (*this).Xvalues.at(iSegEq).at(nCells*iEq+i);
	
}

void MASCH_Variables::setSparCSR(MASCH_Mesh& mesh, MASCH_Control& controls){
	// auto& var = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	int cellSize = mesh.cells.size();
	this->nCells = cellSize;
	
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14" << endl;

	vector<int> recv_rank;
	vector<int> recv_iL;
	if(size>1){
		int proc_size = (*this).procRightCells.size();
		vector<int> send_value0;
		vector<int> send_value1;
		// send_value0.reserve(proc_size);
		// send_value1.reserve(proc_size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					send_value0.push_back(rank);
					send_value1.push_back(mesh.faces[i].iL);
				}
			}
		}
		recv_rank.resize(proc_size);
		recv_iL.resize(proc_size);
		MPI_Alltoallv( send_value0.data(), mesh.countsProcFaces.data(), 
						mesh.displsProcFaces.data(), MPI_INT, 
						recv_rank.data(), mesh.countsProcFaces.data(), 
						mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( send_value1.data(), mesh.countsProcFaces.data(), 
						mesh.displsProcFaces.data(), MPI_INT, 
						recv_iL.data(), mesh.countsProcFaces.data(), 
						mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.5" << endl;
	
	int nSegEq = controls.nEq.size();
	
	this->cellcell_displ.clear(); this->cellcell_displ.resize(mesh.cells.size());
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.6" << endl;
	this->face_LR_displ.clear(); this->face_LR_displ.resize(mesh.faces.size());
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.7" << endl;
	this->face_RL_displ.clear(); this->face_RL_displ.resize(mesh.faces.size());
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.8" << endl;
	this->cellnFaces.clear(); this->cellnFaces.resize(mesh.cells.size());
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.9" << endl;
	
	this->i_str_CSR.clear(); this->i_str_CSR.resize(nSegEq); 
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.1" << endl;
	// this->Avalues[0].clear();
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.2" << endl;
	this->Avalues.clear();
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.dddd" << endl; 
	this->Avalues.resize(nSegEq); 
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.2" << endl;
	this->Xvalues.clear(); this->Xvalues.resize(nSegEq); 
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.3" << endl;
	this->Bvalues.clear(); this->Bvalues.resize(nSegEq);
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.4" << endl;
	this->j_displ_CSR.clear(); this->j_displ_CSR.resize(nSegEq);
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START14.5" << endl;
	
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START15" << endl;
	
	for(int iSegEq=0; iSegEq<nSegEq; ++iSegEq){
			
		int nEq = controls.nEq[iSegEq];
		
		i_str_CSR[iSegEq].resize(nEq*mesh.cells.size()+1,0);
		
		for(int i=0, str = 0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			auto& cell = mesh.cells[i];
			int tmp_size = 1;
			for(auto& iface : cell.ifaces){
				auto& face = mesh.faces[iface];
				if(
				face.getType() == MASCH_Face_Types::INTERNAL ||
				face.getType() == MASCH_Face_Types::PROCESSOR
				){
					++tmp_size;
				}
			}
			cellnFaces[i] = tmp_size-1;
			str += tmp_size*nEq;
			i_str_CSR[iSegEq][i+1] = str;
		}
		for(int ii=1; ii<nEq; ++ii){
			for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
				i_str_CSR[iSegEq].at(ii*cellSize+i+1) = i_str_CSR[iSegEq].at(ii*cellSize+i) + (cellnFaces[i]+1)*nEq;
			}
		}
		
		// 변수 사이즈 조정
		this->Avalues[iSegEq].resize(i_str_CSR[iSegEq][nEq*cellSize]);
		this->Xvalues[iSegEq].resize(nEq*cellSize);
		this->Bvalues[iSegEq].resize(nEq*cellSize);
		
		
		this->j_displ_CSR[iSegEq].resize(i_str_CSR[iSegEq][nEq*cellSize],-10);
		
		
		// vector<int> proc_face_iRs(mesh.faces.size(),-10);
		vector<int> str_glo_R_proc(mesh.faces.size(),-10);
		vector<int> step_loc_R_proc(mesh.faces.size(),-10);
		vector<int> i_loc_R_proc(mesh.faces.size(),-10);
		for(int i=0, ip=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
			if(mesh.faces[i].getType() == MASCH_Face_Types::PROCESSOR){
				str_glo_R_proc[i] = mesh.startProcCellGlobal[recv_rank[ip]];
				step_loc_R_proc[i] = mesh.startProcCellGlobal[recv_rank[ip]+1] -
									 mesh.startProcCellGlobal[recv_rank[ip]];
				i_loc_R_proc[i] = recv_iL[ip];
				
				// proc_face_iRs[i] = mesh.startProcCellGlobal[recv_rank[ip]]*nEq + recv_iL[ip];
				++ip;
			}
		}
		
		if(mesh.startProcCellGlobal[size] != mesh.ncellsTotal){
			cout << "#WARNING, " << "mesh.startProcCellGlobal[rank] != mesh.ncellsTotal" << endl;
		}


		// int str_glo_L = mesh.startCellGlobal*B_n;
		// int str_glo_R = mesh.startCellGlobal*B_n;
		// int step_loc_L = mesh.cells.size();
		// int step_loc_R = mesh.cells.size();
		// int i_loc_L = face.owner;
		// int i_loc_R = face.neighbour;
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// str_glo_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]]*B_n;
			// step_loc_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]+1] - 
			             // mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
			// i_loc_R = mesh.procNeighbCellNo[proc_num];
		// }
		
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START16" << endl;
		
		int max_test = -10;
		for(int i=0, ip=0; i<cellSize; ++i){
			auto& cell = mesh.cells[i];
			vector<int> face_iRs;
			vector<bool> faceL;
			vector<string> face_stype;
			vector<int> proc_str, proc_iR;
			for(auto& iface : cell.ifaces){
				auto& face = mesh.faces[iface];
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					if(face.iL==i){
						faceL.push_back(true);
						face_stype.push_back("internal face iR");
						proc_str.push_back(-10); proc_iR.push_back(-10);
						face_iRs.push_back(mesh.startProcCellGlobal[rank]+face.iR);
					}
					else{
						faceL.push_back(false);
						face_stype.push_back("internal face iL");
						proc_str.push_back(-10); proc_iR.push_back(-10);
						face_iRs.push_back(mesh.startProcCellGlobal[rank]+face.iL);
					}
				}
				else if(face.getType() == MASCH_Face_Types::PROCESSOR){
					faceL.push_back(true);
					face_stype.push_back("procs face");
					proc_str.push_back(str_glo_R_proc[iface]); proc_iR.push_back(i_loc_R_proc[iface]);
					face_iRs.push_back(str_glo_R_proc[iface] + i_loc_R_proc[iface]);
					if(str_glo_R_proc[iface]<0){
						cout << "#WARNING : str_glo_R_proc[iface]<0" << endl;
					}
					if(i_loc_R_proc[iface]<0){
						cout << "#WARNING : i_loc_R_proc[iface]<0" << endl;
					}
				}
			}
			
			
			
			face_stype.push_back("cell");
			proc_str.push_back(-10); proc_iR.push_back(-10);
			face_iRs.push_back(mesh.startProcCellGlobal[rank]+i);
			int step = face_iRs.size();
			
			vector<int> sorted_face_iRs = face_iRs;
			sort(sorted_face_iRs.begin(), sorted_face_iRs.end());
			vector<int> tmp_order(step,-1);
			for(int j=0; j<step; ++j){
				int num = sorted_face_iRs.at(j);
				int order = find(face_iRs.begin(), face_iRs.end(), num) - face_iRs.begin();
				tmp_order.at(order) = j;
			}
			for(int j=0; j<step; ++j){
				if(tmp_order[j]==-1){
					cout << "#WARNING !! : tmp_order[j]==-1" << endl;
					int tmpiter2=0;
					for(auto& item : face_iRs){
						cout << item << endl;
						cout << face_stype[tmpiter2] << endl;
						cout << proc_str[tmpiter2] << endl;
						cout << proc_iR[tmpiter2] << endl;
						++tmpiter2;
					}
					break;
				}
			}
			
			int iter = 0;
			for(auto& iface : cell.ifaces){
				auto& face = mesh.faces[iface];
				if(face.getType() == MASCH_Face_Types::INTERNAL){
					if(face.iL==i){
						// face_LR_displ[iface] = i_str_CSR[i]+tmp_order[iter];
						this->face_LR_displ[iface] = tmp_order[iter];
						for(int ii=0; ii<nEq; ++ii){
							int str = i_str_CSR[iSegEq][ii*cellSize+i];
							for(int jj=0; jj<nEq; ++jj){
								
								
								this->j_displ_CSR[iSegEq].at(str+step*jj+tmp_order[iter]) = 
									mesh.startProcCellGlobal[rank]*nEq + 
									jj*cellSize + 
									face.iR;
									
								// max_test = max(max_test,j_displ_CSR[iSegEq].at(str+step*jj+tmp_order[iter]));
							}
						}
					}
					else{
						// face_RL_displ[iface] = i_str_CSR[i]+tmp_order[iter];
						this->face_RL_displ[iface] = tmp_order[iter];
						for(int ii=0; ii<nEq; ++ii){
							int str = i_str_CSR[iSegEq][ii*cellSize+i];
							for(int jj=0; jj<nEq; ++jj){
								
								
								this->j_displ_CSR[iSegEq].at(str+step*jj+tmp_order[iter]) = 
									mesh.startProcCellGlobal[rank]*nEq + 
									jj*cellSize + 
									face.iL;
									
									
								// max_test = max(max_test,j_displ_CSR[iSegEq].at(str+step*jj+tmp_order[iter]));
							}
						}
					}
					++iter;
				}
				else if(face.getType() == MASCH_Face_Types::PROCESSOR){
					// face_LR_displ[iface] = i_str_CSR[i]+tmp_order[iter];
					this->face_LR_displ.at(iface) = tmp_order.at(iter);
					for(int ii=0; ii<nEq; ++ii){
						int str = i_str_CSR[iSegEq].at(ii*cellSize+i);
						for(int jj=0; jj<nEq; ++jj){
							
							
							this->j_displ_CSR[iSegEq].at(str+step*jj+tmp_order.at(iter)) = 
								str_glo_R_proc[iface]*nEq +
								jj*step_loc_R_proc[iface] + 
								i_loc_R_proc[iface];
								
								
						}
					}
					++iter;
				}
			}
			
			
			this->cellcell_displ[i] = tmp_order.back();
			for(int ii=0; ii<nEq; ++ii){
				int str = i_str_CSR[iSegEq][ii*cellSize+i];
				for(int jj=0; jj<nEq; ++jj){
					
					
					this->j_displ_CSR[iSegEq][str+step*jj+tmp_order.back()] = 
						mesh.startProcCellGlobal[rank]*nEq + 
						jj*cellSize + 
						i;
						
					// max_test = max(max_test,j_displ_CSR[iSegEq].at(str+step*jj+tmp_order.back()));
				}
			}
			
		}
		
			// MPI_Barrier(MPI_COMM_WORLD);
			// if(rank==0) cout << "START17" << endl;

		// cout << max_test << " " << cellSize << endl;
		
	}
        // id = step_loc_L*(B_n*0+1) + i_loc_L; 
		// i_glo = str_glo_L + step_loc_L*0 + i_loc_L; 
		// j_glo = str_glo_R + step_loc_R*1 + i_loc_R;
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// 임시 변수
	// tmp_sparA.resize(nEq*nEq*cellSize);
	// tmp_B.resize(nEq*cellSize);
	// tmp_X.resize(nEq*cellSize);
	
	
}

