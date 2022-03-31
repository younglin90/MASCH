
#include "./mesh.h"

void MASCH_Mesh::repartParMETIS(
int nSize, vector<int>& cell_ip, MASCH_Mesh &mesh){

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute ParMETIS ... ";
	}
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	options[0] = 0;
	options[1] = 0;
	options[2] = 0;
	
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nSize*ncon];
	for(int i=0;i<nSize*ncon;++i) 
		tpwgts[i]=1.0/nSize;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.02);
	
	real_t ubvec = 1.02;
	
	
	vector<int> temp_vtxdist(nSize);
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	// vector<int> vtxdist(nSize+1);
	vector<int> vtxdist(nSize+1);
	vtxdist[0] = 0;
	for(int i=1; i<nSize+1; ++i){
		vtxdist.at(i) = vtxdist.at(i-1) + (temp_vtxdist.at(i-1));
	}
	
	
	vector<int> xadj(ncells+1,0);
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		// for(auto& face : cell.faces()){
			// if(
			// (*face).getType() == MASCH_Face_Types::INTERNAL ||
			// (*face).getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		// }
		for(auto& iface : cell.ifaces){
			auto& face = mesh.faces[iface];
			if(
			face.getType() == MASCH_Face_Types::INTERNAL ||
			face.getType() == MASCH_Face_Types::PROCESSOR) ++numt2;
		}
		++numt;
		xadj.at(numt) = xadj.at(numt-1) + (numt2);
	}
	
	
	vector<int> adjncy(xadj.at(ncells),0);
	
	vector<int> recv_right_rank;
	vector<int> recv_right_icell;
	if(size>1){
		vector<int> send_right_rank;
		vector<int> send_right_icell;
		for(auto& face : mesh.faces){
			if(face.getType() == MASCH_Face_Types::PROCESSOR) {
				send_right_rank.push_back(rank);
				send_right_icell.push_back(face.iL);
			}
		}
		recv_right_rank.resize(send_right_rank.size());
		recv_right_icell.resize(send_right_icell.size());
		MPI_Alltoallv( send_right_rank.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_right_rank.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( send_right_icell.data(), mesh.countsSendProcFaces.data(), mesh.displsSendProcFaces.data(), MPI_INT, 
					   recv_right_icell.data(), mesh.countsRecvProcFaces.data(), mesh.displsRecvProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	vector<int> num_ncell(ncells,0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL) {
			
			int tnum = xadj.at(face.iL) + num_ncell.at(face.iL);
			adjncy.at(tnum) = vtxdist.at(rank) + face.iR;
			++num_ncell.at(face.iL);
			
			tnum = xadj.at(face.iR) + num_ncell.at(face.iR);
			adjncy.at(tnum) = vtxdist.at(rank) + face.iL;
			++num_ncell.at(face.iR);
			
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR) {
			
			int tnum = xadj.at(face.iL) + num_ncell.at(face.iL);
			int right_rank = recv_right_rank.at(proc_num);
			adjncy.at(tnum) = vtxdist.at(right_rank) + recv_right_icell.at(proc_num);
			++num_ncell.at(face.iL);
			
			++proc_num;
		}
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// for(auto& item : adjncy){
		// cout << item << endl;
	// }

	// ParMETIS_V3_PartKway(
		// vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, &wgtflag, &numflag,
		// &ncon, &nSize, tpwgts, &ubvec,
		// options, &objval, cell_ip.data(), &comm);
		
	idx_t vsize = NULL;
	real_t itr = 1000.0;
	ParMETIS_V3_AdaptiveRepart(
		vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, &wgtflag, &numflag,
		&ncon, &nSize, tpwgts, &ubvec, &itr, 
		options, &objval, cell_ip.data(), &comm);
		
	
	// 그룹 재정립
	vector<int> nIps(nSize,0);
	for(int i=0; i<mesh.cells.size(); ++i){
		int group0 = mesh.cells[i].group;
		int cell_ip0 = cell_ip[i];
		int group = group0;
		vector<int> tmp_icell;
		tmp_icell.push_back(i);
		vector<int> tmp_ip;
		tmp_ip.push_back(cell_ip[i]); ++nIps[cell_ip[i]];
		int ip_max = cell_ip[i];
		while(1){
			++i;
			if(i==mesh.cells.size()) break;
			group = mesh.cells[i].group;
			if(group0!=group) break;
			tmp_icell.push_back(i);
			tmp_ip.push_back(cell_ip[i]); ++nIps[cell_ip[i]];
			if(nIps[cell_ip[i-1]]<nIps[cell_ip[i]]) ip_max = cell_ip[i];
		}
		--i;
		
		for(auto& icell : tmp_icell){
			cell_ip[icell] = ip_max;
		}
	}
	
	// // 디버깅
	// {
		// for(int i=0; i<mesh.cells.size(); ++i){
			// int group0 = mesh.cells[i].group;
			// int cell_ip0 = cell_ip[i];
			// int group = group0;
			// while(1){
				// group = mesh.cells[++i].group;
				// if(group0!=group) break;
				// if(cell_ip0!=cell_ip[i]){
					// cout << "not equal group's cell_ip, group0 ip = " << cell_ip0 <<
					// " another group ip = " << cell_ip[i] << endl;
				// }
			// }
			// --i;
		// }
	// }
	
	
	// for(auto& item : cell_ip){
		// item = rank;
	// }
	
		
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}


