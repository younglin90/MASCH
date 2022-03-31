
#include "./partition.h" 

void MASCH_Mesh_Partition::parMETIS_Graph_Partition(
int nSize, vector<int>& procNoCell, MASCH_Mesh &mesh){

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
		vtxdist[i] = vtxdist[i-1] + static_cast<int>(temp_vtxdist[i-1]);
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
		xadj[numt] = xadj[numt-1] + static_cast<int>(numt2);
	}
	
	
	vector<int> adjncy(xadj[ncells],0);
	
	vector<int> num_ncell(ncells,0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == MASCH_Face_Types::INTERNAL) {
			
			int tnum = xadj[face.iL] + num_ncell[face.iL];
			adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iR);
			++num_ncell[face.iL];
			
			tnum = xadj[face.iR] + num_ncell[face.iR];
			adjncy[tnum] = vtxdist[rank] + static_cast<int>(face.iL);
			++num_ncell[face.iR];
			
		}
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// int tnum = xadj[face.owner] + num_ncell[face.owner];
			// int rank_neighbour = iRank_Recv[proc_num];
			// // cout << proc_num << " " << iRank_Recv.size() << endl;
			// adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// ++num_ncell[face.owner];
			
			// ++proc_num;
		// }
	}
	

	ParMETIS_V3_PartKway(
		vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, &wgtflag, &numflag,
		&ncon, &nSize, tpwgts, &ubvec,
		options, &objval, procNoCell.data(), &comm);
		
		
	// idx_t vsize = NULL;
	// real_t itr = 1000.0;
	// ParMETIS_V3_AdaptiveRepart(
		// vtxdist.data(), xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, &wgtflag, &numflag,
		// &ncon, &nSize, tpwgts, &ubvec, &itr, 
		// options, &objval, procNoCell.data(), &comm);

		
	// libSCOTCH_PartGraph(xadj, adjncy, &nSize, procNoCell);


	// for(int i=0; i<ncells; ++i){
		// procNoCell[i] = 1;
	// }
	// for(int i=0; i<ncells/2; ++i){
		// procNoCell[i] = 0;
	// }

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}





// void parMETIS_Mesh_Partition(int nSize, vector<int>& procNoCell){

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();

    // if(rank == 0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute ParMETIS ... ";
	// }
		
	// int ncells = mesh.cells.size();
	// int npoints = mesh.points.size();
	// int nfaces = mesh.faces.size();
	// int ncon=1;
	// int ncommon=3;
	// int objval;
	
	// int options[METIS_NOPTIONS];
	// METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	

	// MPI_Comm comm;
	// MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	
	// int wgtflag=0;
	// int numflag=0;
	// real_t tpwgts[nSize*ncon];
	// for(int i=0;i<nSize*ncon;++i) 
		// tpwgts[i]=1.0/nSize;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.05);
	
	
	
	// int temp_vtxdist[nSize];
    // MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	// int vtxdist[nSize+1];
	// vtxdist[0] = 0;
	// for(int i=1; i<nSize+1; ++i){
		// vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	// }
	
	
	
	
	// vector<int> xadj(ncells+1,0);
	// xadj[0] = 0;
	// for(int i=1; i<ncells+1; ++i){
		// xadj[i] = xadj[i-1] + mesh.cells[i-1].points.size();
	// }
	
	
	
	// // vector<int> num_ncell(ncells,0);
	// // vector<int> adjncy(xadj[ncells],0);
	// vector<int> adjncy;
	// int proc_num = 0;
	// for(int i=0; i<ncells; ++i){
		// for(auto& point : mesh.cells[i].points){
			// adjncy.push_back(point);
		// }
	// }
		
	// // for(auto& face : mesh.points){
		
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			// // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // adjncy[tnum] = vtxdist[rank] + face.neighbour;
			// // ++num_ncell[face.owner];
			
			// // // cout << face.neighbour << endl;
			
			// // tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			// // adjncy[tnum] = vtxdist[rank] + face.owner;
			// // ++num_ncell[face.neighbour];
			
		// // }
		// // // else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// // // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // // int rank_neighbour = iRank_Recv[proc_num];
			// // // // cout << proc_num << " " << iRank_Recv.size() << endl;
			// // // adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// // // ++num_ncell[face.owner];
			
			// // // ++proc_num;
		// // // }
	// // }
	

	// int ncommonnodes = 3;

	// // ParMETIS_V3_PartMeshKway(
		// // vtxdist, xadj.data(), adjncy.data(), NULL, &wgtflag, &numflag,
		// // &ncon, &ncommonnodes, &nSize, tpwgts, ubvec,
		// // options, &objval, procNoCell.data(), &comm);

	// // objval = adjncy.size()+1;

	// vector<int> npart(npoints,0);

	// METIS_PartMeshDual(
		// &ncells, &npoints, xadj.data(), adjncy.data(), NULL, NULL,
		// &ncommonnodes, &nSize, NULL, NULL,
		// &objval, procNoCell.data(), npart.data());

	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	
	
// }





