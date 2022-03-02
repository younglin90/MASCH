
#include "./variables.h"

void MASCH_Variables::setSparCSR(MASCH_Mesh& mesh, MASCH_Control& controls){
	// auto& var = (*this);
	
	int nEq = controls.nEq;
	int cellSize = mesh.cells.size();
	
	// var.cellcell_displ = new int[mesh.cells.size()];
	// var.face_LR_displ = new int[mesh.faces.size()];
	// var.face_RL_displ = new int[mesh.faces.size()];
	// cellnFaces = new int[mesh.cells.size()];
	cellcell_displ.resize(mesh.cells.size());
	face_LR_displ.resize(mesh.faces.size());
	face_RL_displ.resize(mesh.faces.size());
	cellnFaces.resize(mesh.cells.size());
	
	i_str_CSR.clear();
	i_str_CSR.resize(nEq*mesh.cells.size()+1,0);
	
	for(int i=0, str = 0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		int tmp_size = 1;
		for(auto& iface : cell.ifaces){
			auto& face = mesh.faces[iface];
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				++tmp_size;
			}
		}
		cellnFaces[i] = tmp_size-1;
		str += tmp_size*nEq;
		i_str_CSR[i+1] = str;
	}
	for(int ii=1; ii<nEq; ++ii){
		for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
			i_str_CSR[ii*cellSize+i+1] = i_str_CSR[ii*cellSize+i] + 
			i_str_CSR[(ii-1)*cellSize+i+1] - i_str_CSR[(ii-1)*cellSize+i];
		}
	}
	
	j_displ_CSR.clear();
	j_displ_CSR.resize(i_str_CSR[nEq*cellSize]);
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		vector<int> face_iRs;
		vector<bool> faceL;
		for(auto& iface : cell.ifaces){
			auto& face = mesh.faces[iface];
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				if(face.iL==i){
					faceL.push_back(true);
					face_iRs.push_back(face.iR);
				}
				else{
					faceL.push_back(false);
					face_iRs.push_back(face.iL);
				}
			}
		}
		face_iRs.push_back(i);
		
		vector<int> sorted_face_iRs = face_iRs;
		sort(sorted_face_iRs.begin(), sorted_face_iRs.end());
		vector<int> tmp_order(face_iRs.size(),-1);
		for(int j=0; j<sorted_face_iRs.size(); ++j){
			int num = sorted_face_iRs[j];
			auto it = find(face_iRs.begin(), face_iRs.end(), num);
			
			int order = (it - face_iRs.begin());
			tmp_order[order] = j;
		}
		
		int step = face_iRs.size();
		
		int iter = 0;
		for(auto& iface : cell.ifaces){
			auto& face = mesh.faces[iface];
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				if(face.iL==i){
					// face_LR_displ[iface] = i_str_CSR[i]+tmp_order[iter];
					face_LR_displ[iface] = tmp_order[iter];
					for(int ii=0; ii<nEq; ++ii){
						int str = i_str_CSR[ii*cellSize+i];
						for(int jj=0; jj<nEq; ++jj){
							j_displ_CSR[str+step*jj+tmp_order[iter]] = jj*cellSize + face.iR;
						}
					}
				}
				else{
					// face_RL_displ[iface] = i_str_CSR[i]+tmp_order[iter];
					face_RL_displ[iface] = tmp_order[iter];
					for(int ii=0; ii<nEq; ++ii){
						int str = i_str_CSR[ii*cellSize+i];
						for(int jj=0; jj<nEq; ++jj){
							j_displ_CSR[str+step*jj+tmp_order[iter]] = jj*cellSize + face.iL;
						}
					}
				}
				++iter;
			}
		}
		
		cellcell_displ[i] = i_str_CSR[i]+tmp_order.back();
		for(int ii=0; ii<nEq; ++ii){
			int str = i_str_CSR[ii*cellSize+i];
			for(int jj=0; jj<nEq; ++jj){
				j_displ_CSR[str+step*jj+tmp_order.back()] = jj*cellSize + i;
			}
		}
		
	}
	
	// 변수 사이즈 조정
	Avar.resize(i_str_CSR[nEq*cellSize]);
	Xvar.resize(nEq*mesh.cells.size());
	Bvar.resize(nEq*mesh.cells.size());
	
}





void MASCH_Variables::accumSparD(int i, int iEq, int jEq, double inp){
	auto& var = (*this);
	int str = i_str_CSR[nCells*iEq + i];
	int tmp_i = str + (cellnFaces[i]+1)*jEq + cellcell_displ[i];
	var.Avar[tmp_i] += inp;
	
}

double MASCH_Variables::getSparD(int i, int iEq, int jEq){
	auto& var = (*this);
	int str = i_str_CSR[nCells*iEq + i];
	int tmp_i = str + (cellnFaces[i]+1)*jEq + cellcell_displ[i];
	return var.Avar[tmp_i];
}


void MASCH_Variables::accumSparLR(int i, int iL, int iEq, int jEq, double inp){
	auto& var = (*this);
	int str = i_str_CSR[nCells*iEq + iL];
	int tmp_i = str + (cellnFaces[iL]+1)*jEq + face_LR_displ[i];
	var.Avar[tmp_i] += inp;
	
}

void MASCH_Variables::accumSparRL(int i, int iR, int iEq, int jEq, double inp){
	auto& var = (*this);
	int str = i_str_CSR[nCells*iEq + iR];
	int tmp_i = str + (cellnFaces[iR]+1)*jEq + face_RL_displ[i];
	var.Avar[tmp_i] += inp;
	
}


void MASCH_Variables::accumB(int i, int iEq, double inp){
	auto& var = (*this);
	int str = nCells*iEq;
	int tmp_i = str + i;
	var.Bvar[tmp_i] += inp;
	
}
double MASCH_Variables::getB(int i, int iEq){
	auto& var = (*this);
	int str = nCells*iEq;
	int tmp_i = str + i;
	return var.Bvar[tmp_i];
	
}

void MASCH_Variables::setX(int i, int iEq, double inp){
	auto& var = (*this);
	int str = nCells*iEq;
	int tmp_i = str + i;
	var.Xvar[tmp_i] += inp;
	
}
double MASCH_Variables::getX(int i, int iEq){
	auto& var = (*this);
	int str = nCells*iEq;
	int tmp_i = str + i;
	return var.Xvar[tmp_i];
	
}