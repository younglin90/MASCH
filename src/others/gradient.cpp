
#include "./solvers.h"
void MASCH_Gradient::gaussGreen(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int inp, string inp_bc, int oup0, int oup1, int oup2){
	
	
	
	
	
	
}


void MASCH_Gradient::init(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	int id_coeff1 = controls.getId_cellVar("coeff of least square-(1,1)");
	int id_coeff2 = controls.getId_cellVar("coeff of least square-(1,2)");
	int id_coeff3 = controls.getId_cellVar("coeff of least square-(1,3)");
	int id_coeff4 = controls.getId_cellVar("coeff of least square-(2,2)");
	int id_coeff5 = controls.getId_cellVar("coeff of least square-(2,3)");
	int id_coeff6 = controls.getId_cellVar("coeff of least square-(3,3)");
	
	// cout << id_coeff1 << endl;
	
	// if(rank==0){
		// cout << id_coeff1 << " " << id_coeff6 << endl;
	// }

	// SEMO_MPI_Builder mpi;
	
	// processor faces
	if(size>1){

		for(int i=0, proc=0; i<mesh.faces.size(); ++i){
			if(mesh.faces[i].getType() == MASCH_Face_Types::PROCESSOR) ++proc;
		}
		
		vector<double> send_x, send_y, send_z;
		for(auto& icell : mesh.send_StencilCellsId){
			send_x.push_back(mesh.cells.at(icell).x);
			send_y.push_back(mesh.cells.at(icell).y);
			send_z.push_back(mesh.cells.at(icell).z);
		}
		mesh.recv_x_StencilCells.resize(mesh.recv_displsStencilCells[size]);
		mesh.recv_y_StencilCells.resize(mesh.recv_displsStencilCells[size]);
		mesh.recv_z_StencilCells.resize(mesh.recv_displsStencilCells[size]);
		
		// if(rank==1){
			// cout << mesh.send_countsStencilCells[0] << endl;
		// }
		
		MPI_Alltoallv( send_x.data(), mesh.send_countsStencilCells.data(), 
					   mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
					   mesh.recv_x_StencilCells.data(), mesh.recv_countsStencilCells.data(), 
					   mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( send_y.data(), mesh.send_countsStencilCells.data(), 
					   mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
					   mesh.recv_y_StencilCells.data(), mesh.recv_countsStencilCells.data(), 
					   mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( send_z.data(), mesh.send_countsStencilCells.data(), 
					   mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
					   mesh.recv_z_StencilCells.data(), mesh.recv_countsStencilCells.data(), 
					   mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		
		
		
	}

	vector<vector<double>> vsum(mesh.cells.size(),vector<double>(6,0.0));

	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		
		for(auto& icell : cell.iStencils){
			auto& cellSten = mesh.cells[icell];
			
			double distX = cellSten.x - cell.x;
			double distY = cellSten.y - cell.y;
			double distZ = cellSten.z - cell.z;
			
			
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// wk = pow(wk,n_weight);
			
			// if(isnan(wk) || wk<-10000000.0 || wk>10000000.0){
				// cout << i << " " << icell << endl;
			// }
			// vector<double> vari(9,0.0);;
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			// 대칭행렬
			for(int ii=0, num=0; ii<3; ++ii){
				for(int jj=0; jj<3; ++jj){
					if(jj>=ii){
						vsum[i][num++] += wk * vari[ii]*vari[jj];
					}
				}
			}
			
		}
		
		for(auto& icell : cell.recv_iStencils){
			// auto& cellSten = mesh.cells[icell];
			// double distX = recv_x[icell] - cell.x;
			// double distY = recv_y[icell] - cell.y;
			// double distZ = recv_z[icell] - cell.z;
			double distX = mesh.recv_x_StencilCells.at(icell) - cell.x;
			double distY = mesh.recv_y_StencilCells.at(icell) - cell.y;
			double distZ = mesh.recv_z_StencilCells.at(icell) - cell.z;
		
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// // wk = pow(wk,n_weight);
			// if(rank==1){
			// if(isnan(wk) || wk<-10000000.0 || wk>10000000.0){
				// cout << i << " " << icell << endl;
				// cout << mesh.recv_x_StencilCells[icell] << " " << 
				// cell.x << endl;
				// cout << mesh.recv_y_StencilCells[icell] << " " << 
				// cell.y << endl;
				// cout << mesh.recv_z_StencilCells[icell] << " " << 
				// cell.z << endl;
			// }
			// }
			// // vector<double> vari(9,0.0);;
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			// 대칭행렬
			for(int ii=0, num=0; ii<3; ++ii){
				for(int jj=0; jj<3; ++jj){
					if(jj>=ii){
						vsum[i][num++] += wk * vari[ii]*vari[jj];
					}
				}
			}
		}
	}
		
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundaries){
		
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.iL];
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
				// wk = pow(wk,n_weight);
			
				// vector<double> vari(9,0.0);;
				double vari[3];
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
			
				// 대칭행렬
				for(int ii=0, num=0; ii<3; ++ii){
					for(int jj=0; jj<3; ++jj){
						if(jj>=ii){
							vsum[face.iL][num++] += wk * vari[ii]*vari[jj];
						}
					}
				}

			}
		}
	}
	
	// 역행렬 계산 뒤 저장
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		
		// // 1차 least-square 저장
		double detA = 
			vsum[i][0]*vsum[i][3]*vsum[i][5] + vsum[i][1]*vsum[i][4]*vsum[i][2]
		  + vsum[i][2]*vsum[i][1]*vsum[i][4] - vsum[i][0]*vsum[i][4]*vsum[i][4]
		  - vsum[i][2]*vsum[i][3]*vsum[i][2] - vsum[i][1]*vsum[i][1]*vsum[i][5];
		// cout << vsum[i][0] << " " << vsum[i][5] << endl;
		if(detA==0.0 || std::isnan(detA)){
			cerr << "| #Error, detA=0.0 at leat-sqare" << endl;
			// cout << vsum[i][0] << endl;
			// cout << vsum[i][1] << endl;
			// cout << vsum[i][2] << endl;
			// cout << vsum[i][3] << endl;
			// cout << vsum[i][4] << endl;
			// cout << vsum[i][5] << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		var.cells[i][id_coeff1] = 
			(vsum[i][3] * vsum[i][5] - vsum[i][4] * vsum[i][4]) / detA;    // inv_A(1,1)
		var.cells[i][id_coeff2] = 
			(vsum[i][2] * vsum[i][4] - vsum[i][1] * vsum[i][5]) / detA;    // inv_A(1,2) = (2,1)
		var.cells[i][id_coeff3] = 
			(vsum[i][1] * vsum[i][4] - vsum[i][2] * vsum[i][3]) / detA;    // inv_A(1,3) = (3,1)
		var.cells[i][id_coeff4] = 
			(vsum[i][0] * vsum[i][5] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
		var.cells[i][id_coeff5] = 
			(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][4]) / detA;    // inv_A(2,3) = (3,2)
		var.cells[i][id_coeff6] = 
			(vsum[i][0] * vsum[i][3] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
		
	}

	
}

// void MASCH_Gradient::leastSquare(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, string inp_s){
	
	// // controls.log.push("test1");
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size();
	
	// double tmp_n_weight = n_weight;
	
	// int id_coeff1 = controls.getId_cellVar("coeff of least square-(1,1)");
	// int id_coeff2 = controls.getId_cellVar("coeff of least square-(1,2)");
	// int id_coeff3 = controls.getId_cellVar("coeff of least square-(1,3)");
	// int id_coeff4 = controls.getId_cellVar("coeff of least square-(2,2)");
	// int id_coeff5 = controls.getId_cellVar("coeff of least square-(2,3)");
	// int id_coeff6 = controls.getId_cellVar("coeff of least square-(3,3)");
	
	// int id_primOrder = controls.primitiveMap[inp_s];
	// bool primId_there = true;
    // if (controls.primitiveMap.find(inp_s) == controls.primitiveMap.end()) {
        // primId_there = false;
    // }
	// string type_bc;
// // primId_there=false;
	// if(primId_there==false) type_bc = "zeroGradient";

	// string tmp_oup0 = "x-gradient "; tmp_oup0+=inp_s;
	// string tmp_oup1 = "y-gradient "; tmp_oup1+=inp_s;
	// string tmp_oup2 = "z-gradient "; tmp_oup2+=inp_s;
	// string tmp_bc_face = "left "; tmp_bc_face+=inp_s;
	
	// int inp = controls.cellVar[inp_s].id;
	// int oup0 = controls.cellVar[tmp_oup0].id;
	// int oup1 = controls.cellVar[tmp_oup1].id;
	// int oup2 = controls.cellVar[tmp_oup2].id;
	// int bc_face = controls.faceVar[tmp_bc_face].id;

	// // cout << inp << endl;
	// // cout << oup0 << endl;
	// // cout << oup1 << endl;
	// // cout << oup2 << endl;
	// // cout << bc_face << endl;
	// // cout << var.cells[0].size() << endl;
	// // cout << var.faces[0].size() << endl;

	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	
	// // processor faces
	// vector<double> recv_value;
	// if(size>1){

		// vector<double> send_value;
		// send_value.reserve(mesh.send_StencilCellsId.size());
		// for(auto& icell : mesh.send_StencilCellsId){
			// auto cellVar_i = cellVar[icell].data();
			// send_value.push_back(cellVar_i[inp]);
		// }
		// recv_value.resize(mesh.recv_displsStencilCells[size]);
		// MPI_Alltoallv( send_value.data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
					   // recv_value.data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
					   // MPI_COMM_WORLD);
		
	// }

	
	
	// // controls.log.pop();
	// // controls.log.push("test2");
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto cellVar_i = cellVar[i].data();
		// cellVar_i[oup0] = 0.0;
		// cellVar_i[oup1] = 0.0;
		// cellVar_i[oup2] = 0.0;
	// }
	
	// // controls.log.pop();
	// // controls.log.push("test3");
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto& cell = cells[i];
		// auto cellVar_i = cellVar[i].data();
		// double cell_var = cellVar_i[inp];
		// double cell_x = cell.x;
		// double cell_y = cell.y;
		// double cell_z = cell.z;
		// auto& p_oup0 = cellVar_i[oup0];
		// auto& p_oup1 = cellVar_i[oup1];
		// auto& p_oup2 = cellVar_i[oup2];
		// // cout << "asdf " << cell.iStencils.size() << endl;
		// for(auto& icell : cell.iStencils){
		// // cout << icell << endl;
			// auto& cellSten = cells[icell];
			// auto cellStenVar_i = cellVar[icell].data();
			
			// double distX = cellSten.x - cell_x;
			// double distY = cellSten.y - cell_y;
			// double distZ = cellSten.z - cell_z;
			
			// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// // int tmptmp = tmp_n_weight;
			// // wk = pow(wk,tmptmp);
			
			// // vector<double> vari(9,0.0);
			// double vari[3];
			// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			// double DVar = cellStenVar_i[inp] - cell_var;
			// p_oup0 += wk * vari[0] * DVar;
			// p_oup1 += wk * vari[1] * DVar;
			// p_oup2 += wk * vari[2] * DVar;
			
		// }
		
		// for(auto& icell : cell.recv_iStencils){
			// auto cellSten_x = mesh.recv_x_StencilCells.data();
			// auto cellSten_y = mesh.recv_y_StencilCells.data();
			// auto cellSten_z = mesh.recv_z_StencilCells.data();
			// // double distX = recv_x[icell] - cell.x;
			// // double distY = recv_y[icell] - cell.y;
			// // double distZ = recv_z[icell] - cell.z;
			// double distX = cellSten_x[icell] - cell.x;
			// double distY = cellSten_y[icell] - cell.y;
			// double distZ = cellSten_z[icell] - cell.z;
		
			// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// // wk = pow(wk,n_weight);
			
			// // vector<double> vari(9,0.0);;
			// double vari[3];
			// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			// double DVar = recv_value[icell] - cellVar_i[inp];
			// cellVar_i[oup0] += wk * vari[0] * DVar;
			// cellVar_i[oup1] += wk * vari[1] * DVar;
			// cellVar_i[oup2] += wk * vari[2] * DVar;
		// }
	// }
		
	// // controls.log.pop();
	// // controls.log.push("test4");
	
	// // boundary face's nodes
	// for(auto& boundary : mesh.boundaries){
		
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// auto& cell = cells[face.iL];
				// auto cellVar_i = cellVar[face.iL].data();
				// auto faceVar_i = faceVar[i].data();
					
				// double distX = face.x - cell.x;
				// double distY = face.y - cell.y;
				// double distZ = face.z - cell.z;
				
				// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
				// // wk = pow(wk,n_weight);
			
				// // vector<double> vari(9,0.0);;
				// double vari[3];
				// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				
				
				// double DVar = faceVar_i[bc_face] - cellVar_i[inp];
				// cellVar_i[oup0] += wk * vari[0] * DVar;
				// cellVar_i[oup1] += wk * vari[1] * DVar;
				// cellVar_i[oup2] += wk * vari[2] * DVar;

			// }
		// }
	// }
	
	// // controls.log.pop();
	// // controls.log.push("test5");
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto& cell = cells[i];
		// auto cellVar_i = cellVar[i].data();

		// double tmp0, tmp1, tmp2;
		// tmp0 = 0.0;
		// tmp0 += cellVar_i[id_coeff1] * cellVar_i[oup0];
		// tmp0 += cellVar_i[id_coeff2] * cellVar_i[oup1];
		// tmp0 += cellVar_i[id_coeff3] * cellVar_i[oup2];
		// tmp1 = 0.0;
		// tmp1 += cellVar_i[id_coeff2] * cellVar_i[oup0];
		// tmp1 += cellVar_i[id_coeff4] * cellVar_i[oup1];
		// tmp1 += cellVar_i[id_coeff5] * cellVar_i[oup2];
		// tmp2 = 0.0;
		// tmp2 += cellVar_i[id_coeff3] * cellVar_i[oup0];
		// tmp2 += cellVar_i[id_coeff5] * cellVar_i[oup1];
		// tmp2 += cellVar_i[id_coeff6] * cellVar_i[oup2];
		
		// cellVar_i[oup0] = tmp0;
		// cellVar_i[oup1] = tmp1;
		// cellVar_i[oup2] = tmp2;
		
	// }
	
	// // controls.log.pop();
	
	// // // 바운더리 커렉션
	// // int maxIbc = 2;
	
	// // int nBoundary = 0;
	// // for(auto& boundary : mesh.boundaries){
		// // if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// // nBoundary += boundary.nFaces;
		// // }
	// // }
	
	// // vector<vector<double>> old_cellVar_bc(nBoundary,vector<double>(3,0.0));
	// // auto p_old_cellVar_bc = old_cellVar_bc.data();
	// // int i_bc = 0;
	// // for(auto& boundary : mesh.boundaries){
		// // if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // for(int i=str; i<end; ++i){
				// // auto& face = faces[i];
				// // auto& cell = cells[face.iL];
				// // auto cellVar_i = cellVar[face.iL].data();
				// // auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
				// // p_old_cellVar_bc_i[0] = cellVar_i[oup0];
				// // p_old_cellVar_bc_i[1] = cellVar_i[oup1];
				// // p_old_cellVar_bc_i[2] = cellVar_i[oup2];
				// // ++i_bc;
			// // }
		// // }
	// // }
	
	// // for(int ibc=0; ibc<maxIbc; ++ibc)
	// // {
		// // vector<vector<double>> cellVar_bc(nBoundary,vector<double>(3,0.0));
		// // auto p_cellVar_bc = cellVar_bc.data();
		// // i_bc = 0;
		// // for(auto& boundary : mesh.boundaries){
			// // if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// // if(primId_there==true){
					// // type_bc.clear();
					// // type_bc = boundary.types[id_primOrder];
				// // }
				
				// // if(type_bc!="zeroGradient") continue;
		
				// // int str = boundary.startFace;
				// // int end = str + boundary.nFaces;
				
				// // for(int i=str; i<end; ++i){
					// // auto& face = faces[i];
					// // auto& cell = cells[face.iL];
					// // auto cellVar_i = cellVar[face.iL].data();
						
					// // double distX = face.x - cell.x;
					// // double distY = face.y - cell.y;
					// // double distZ = face.z - cell.z;
							
					// // double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
					// // // wk = pow(wk,n_weight);
				
					// // // vector<double> vari(9,0.0);;
					// // double vari[3];
					// // vari[0] = distX; vari[1] = distY; vari[2] = distZ;
					
					// // // for(int iEq=0; iEq<nEq; ++iEq){
					// // // }
					
					// // double DVar = 0.0;
					// // DVar += cellVar_i[oup0]*distX;
					// // DVar += cellVar_i[oup1]*distY;
					// // DVar += cellVar_i[oup2]*distZ;
					
					// // // double resi[3];
					// // // resi[0] = wk * vari[0] * DVar;
					// // // resi[1] = wk * vari[1] * DVar;
					// // // resi[2] = wk * vari[2] * DVar;
							
					// // auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// // p_cellVar_bc_i[0] += wk * vari[0] * DVar;
					// // p_cellVar_bc_i[1] += wk * vari[1] * DVar;
					// // p_cellVar_bc_i[2] += wk * vari[2] * DVar;
					// // ++i_bc;
				// // }
			// // }
		// // }
		
		
		// // i_bc = 0;
		// // for(auto& boundary : mesh.boundaries){
			// // if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// // if(primId_there==true){
					// // type_bc.clear();
					// // type_bc = boundary.types[id_primOrder];
				// // }
				
				// // if(type_bc!="zeroGradient") continue;
		
				// // int str = boundary.startFace;
				// // int end = str + boundary.nFaces;
				
				// // for(int i=str; i<end; ++i){
					// // auto& face = faces[i];
					// // auto& cell = cells[face.iL];
					// // auto cellVar_i = cellVar[face.iL].data();
					
					// // auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// // double tmp0, tmp1, tmp2;
					// // tmp0 = 0.0;
					// // tmp0 += cellVar_i[id_coeff1] * p_cellVar_bc_i[0];
					// // tmp0 += cellVar_i[id_coeff2] * p_cellVar_bc_i[1];
					// // tmp0 += cellVar_i[id_coeff3] * p_cellVar_bc_i[2];
					// // tmp1 = 0.0;
					// // tmp1 += cellVar_i[id_coeff2] * p_cellVar_bc_i[0];
					// // tmp1 += cellVar_i[id_coeff4] * p_cellVar_bc_i[1];
					// // tmp1 += cellVar_i[id_coeff5] * p_cellVar_bc_i[2];
					// // tmp2 = 0.0;
					// // tmp2 += cellVar_i[id_coeff3] * p_cellVar_bc_i[0];
					// // tmp2 += cellVar_i[id_coeff5] * p_cellVar_bc_i[1];
					// // tmp2 += cellVar_i[id_coeff6] * p_cellVar_bc_i[2];
				
					// // auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
					// // cellVar_i[oup0] = p_old_cellVar_bc_i[0] + tmp0;
					// // cellVar_i[oup1] = p_old_cellVar_bc_i[1] + tmp1;
					// // cellVar_i[oup2] = p_old_cellVar_bc_i[2] + tmp2;
					
					// // // cout << tmp0 << " " << tmp1 << " " << tmp2 << endl;

					// // ++i_bc;
				// // }
			// // }
		// // }
	// // }
	
	
	
	
// }









void MASCH_Gradient::leastSquare(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
vector<string>& inp_cell, vector<string>& inp_bcFace){
	
	// controls.log.push("test1");
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	double tmp_n_weight = n_weight;
	
	int id_coeff1 = controls.getId_cellVar("coeff of least square-(1,1)");
	int id_coeff2 = controls.getId_cellVar("coeff of least square-(1,2)");
	int id_coeff3 = controls.getId_cellVar("coeff of least square-(1,3)");
	int id_coeff4 = controls.getId_cellVar("coeff of least square-(2,2)");
	int id_coeff5 = controls.getId_cellVar("coeff of least square-(2,3)");
	int id_coeff6 = controls.getId_cellVar("coeff of least square-(3,3)");

	int inp_size = inp_cell.size();

	vector<string> vec_oup0, vec_oup1, vec_oup2, vec_bc_face;
	{
		int tmp_iter = 0;
		for(auto& item : inp_cell){
			string tmp_oup0 = "x-gradient "; tmp_oup0+=item;
			string tmp_oup1 = "y-gradient "; tmp_oup1+=item;
			string tmp_oup2 = "z-gradient "; tmp_oup2+=item;
			string tmp_bc_face = inp_bcFace[tmp_iter];
			
			vec_oup0.push_back(tmp_oup0);
			vec_oup1.push_back(tmp_oup1);
			vec_oup2.push_back(tmp_oup2);
			vec_bc_face.push_back(tmp_bc_face);
			++tmp_iter;
		}
	}
	
	
	vector<int> id_inp, id_oup0, id_oup1, id_oup2, id_bc_face;
	
	int iter_tmp = 0;
	for(auto& item : inp_cell){
		int inp = controls.getId_cellVar(item);
		int oup0 = controls.getId_cellVar(vec_oup0[iter_tmp]);
		int oup1 = controls.getId_cellVar(vec_oup1[iter_tmp]);
		int oup2 = controls.getId_cellVar(vec_oup2[iter_tmp]);
		int bc_face = controls.getId_faceVar(vec_bc_face[iter_tmp]);
		
		id_inp.push_back(inp);
		id_oup0.push_back(oup0);
		id_oup1.push_back(oup1);
		id_oup2.push_back(oup2);
		id_bc_face.push_back(bc_face);
		
		++iter_tmp;
	}

	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// cout << inp << endl;
	// cout << oup0 << endl;
	// cout << oup1 << endl;
	// cout << oup2 << endl;
	// cout << bc_face << endl;
	// cout << var.cells[0].size() << endl;
	// cout << var.faces[0].size() << endl;

	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	// processor faces
	vector<vector<double>> recv_value(inp_size);
	if(size>1){

		vector<vector<double>> send_value(inp_size);
		for(int j=0; j<inp_size; ++j){
			send_value[j].reserve(mesh.send_StencilCellsId.size());
		}
		for(auto& icell : mesh.send_StencilCellsId){
			auto cellVar_i = cellVar[icell].data();
			for(int j=0; j<inp_size; ++j){
				send_value[j].push_back(cellVar_i[id_inp[j]]);
			}
		}
		for(int j=0; j<inp_size; ++j){
			recv_value[j].resize(mesh.recv_displsStencilCells[size]);
			MPI_Alltoallv( send_value[j].data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
						   recv_value[j].data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		
	}

	
	
	// controls.log.pop();
	// controls.log.push("test2");
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(auto& item : id_oup0){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup1){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup2){
			cellVar_i[item] = 0.0;
		}
	}
	
	// controls.log.pop();
	// controls.log.push("test3");
	
	// vector<double> cell_var(id_oup0.size(),0.0);
	// vector<double> cell_var(id_oup0.size(),0.0);
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();
		double cell_x = cell.x;
		double cell_y = cell.y;
		double cell_z = cell.z;
		double cell_var[inp_size];
		double* p_oup0[inp_size];
		double* p_oup1[inp_size];
		double* p_oup2[inp_size];
	// controls.log.push("test2");
		for(int j=0; j<inp_size; ++j){
			cell_var[j] = cellVar_i[id_inp[j]];
			p_oup0[j] = &cellVar_i[id_oup0[j]];
			p_oup1[j] = &cellVar_i[id_oup1[j]];
			p_oup2[j] = &cellVar_i[id_oup2[j]];
		}
		for(auto& icell : cell.iStencils){
		// cout << icell << endl;
			auto& cellSten = cells[icell];
			auto cellStenVar_i = cellVar[icell].data();
			
			double distX = cellSten.x - cell_x;
			double distY = cellSten.y - cell_y;
			double distZ = cellSten.z - cell_z;
			
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// int tmptmp = tmp_n_weight;
			// wk = pow(wk,tmptmp);
			
			// vector<double> vari(9,0.0);
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			double DVar[inp_size];
			for(int j=0; j<inp_size; ++j){
				double DVar = cellStenVar_i[id_inp[j]] - cell_var[j];
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
			
		}
		
	// controls.log.pop();
		for(auto& icell : cell.recv_iStencils){
			auto cellSten_x = mesh.recv_x_StencilCells.data();
			auto cellSten_y = mesh.recv_y_StencilCells.data();
			auto cellSten_z = mesh.recv_z_StencilCells.data();
			// double distX = recv_x[icell] - cell.x;
			// double distY = recv_y[icell] - cell.y;
			// double distZ = recv_z[icell] - cell.z;
			double distX = cellSten_x[icell] - cell_x;
			double distY = cellSten_y[icell] - cell_y;
			double distZ = cellSten_z[icell] - cell_z;
		
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// wk = pow(wk,n_weight);
			
			// vector<double> vari(9,0.0);;
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			for(int j=0; j<inp_size; ++j){
				double DVar = recv_value[j][icell] - cell_var[j];
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
		}
	}
		
	// controls.log.pop();
	// controls.log.push("test4");
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundaries){
		
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				auto& cell = cells[face.iL];
				auto cellVar_i = cellVar[face.iL].data();
				auto faceVar_i = faceVar[i].data();
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
				// wk = pow(wk,n_weight);
			
				// vector<double> vari(9,0.0);;
				double vari[3];
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				
				for(int j=0; j<inp_size; ++j){
					double DVar = faceVar_i[id_bc_face[j]] - cellVar_i[id_inp[j]];
					cellVar_i[id_oup0[j]] += wk * vari[0] * DVar;
					cellVar_i[id_oup1[j]] += wk * vari[1] * DVar;
					cellVar_i[id_oup2[j]] += wk * vari[2] * DVar;
				}
			}
		}
	}
	
	// controls.log.pop();
	// controls.log.push("test5");
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();

		double tmp0, tmp1, tmp2;
		for(int j=0; j<inp_size; ++j){
			tmp0 = 0.0;
			tmp0 += cellVar_i[id_coeff1] * cellVar_i[id_oup0[j]];
			tmp0 += cellVar_i[id_coeff2] * cellVar_i[id_oup1[j]];
			tmp0 += cellVar_i[id_coeff3] * cellVar_i[id_oup2[j]];
			tmp1 = 0.0;
			tmp1 += cellVar_i[id_coeff2] * cellVar_i[id_oup0[j]];
			tmp1 += cellVar_i[id_coeff4] * cellVar_i[id_oup1[j]];
			tmp1 += cellVar_i[id_coeff5] * cellVar_i[id_oup2[j]];
			tmp2 = 0.0;
			tmp2 += cellVar_i[id_coeff3] * cellVar_i[id_oup0[j]];
			tmp2 += cellVar_i[id_coeff5] * cellVar_i[id_oup1[j]];
			tmp2 += cellVar_i[id_coeff6] * cellVar_i[id_oup2[j]];
			
			cellVar_i[id_oup0[j]] = tmp0;
			cellVar_i[id_oup1[j]] = tmp1;
			cellVar_i[id_oup2[j]] = tmp2;
		}
		
	}
	
	// controls.log.pop();
	
	// // 바운더리 커렉션
	// int maxIbc = 2;
	
	// int nBoundary = 0;
	// for(auto& boundary : mesh.boundaries){
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// nBoundary += boundary.nFaces;
		// }
	// }
	
	// vector<vector<double>> old_cellVar_bc(nBoundary,vector<double>(3,0.0));
	// auto p_old_cellVar_bc = old_cellVar_bc.data();
	// int i_bc = 0;
	// for(auto& boundary : mesh.boundaries){
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// auto& cell = cells[face.iL];
				// auto cellVar_i = cellVar[face.iL].data();
				// auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
				// p_old_cellVar_bc_i[0] = cellVar_i[oup0];
				// p_old_cellVar_bc_i[1] = cellVar_i[oup1];
				// p_old_cellVar_bc_i[2] = cellVar_i[oup2];
				// ++i_bc;
			// }
		// }
	// }
	
	// for(int ibc=0; ibc<maxIbc; ++ibc)
	// {
		// vector<vector<double>> cellVar_bc(nBoundary,vector<double>(3,0.0));
		// auto p_cellVar_bc = cellVar_bc.data();
		// i_bc = 0;
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// if(primId_there==true){
					// type_bc.clear();
					// type_bc = boundary.types[id_primOrder];
				// }
				
				// if(type_bc!="zeroGradient") continue;
		
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = faces[i];
					// auto& cell = cells[face.iL];
					// auto cellVar_i = cellVar[face.iL].data();
						
					// double distX = face.x - cell.x;
					// double distY = face.y - cell.y;
					// double distZ = face.z - cell.z;
							
					// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
					// // wk = pow(wk,n_weight);
				
					// // vector<double> vari(9,0.0);;
					// double vari[3];
					// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
					
					// // for(int iEq=0; iEq<nEq; ++iEq){
					// // }
					
					// double DVar = 0.0;
					// DVar += cellVar_i[oup0]*distX;
					// DVar += cellVar_i[oup1]*distY;
					// DVar += cellVar_i[oup2]*distZ;
					
					// // double resi[3];
					// // resi[0] = wk * vari[0] * DVar;
					// // resi[1] = wk * vari[1] * DVar;
					// // resi[2] = wk * vari[2] * DVar;
							
					// auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// p_cellVar_bc_i[0] += wk * vari[0] * DVar;
					// p_cellVar_bc_i[1] += wk * vari[1] * DVar;
					// p_cellVar_bc_i[2] += wk * vari[2] * DVar;
					// ++i_bc;
				// }
			// }
		// }
		
		
		// i_bc = 0;
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// if(primId_there==true){
					// type_bc.clear();
					// type_bc = boundary.types[id_primOrder];
				// }
				
				// if(type_bc!="zeroGradient") continue;
		
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = faces[i];
					// auto& cell = cells[face.iL];
					// auto cellVar_i = cellVar[face.iL].data();
					
					// auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// double tmp0, tmp1, tmp2;
					// tmp0 = 0.0;
					// tmp0 += cellVar_i[id_coeff1] * p_cellVar_bc_i[0];
					// tmp0 += cellVar_i[id_coeff2] * p_cellVar_bc_i[1];
					// tmp0 += cellVar_i[id_coeff3] * p_cellVar_bc_i[2];
					// tmp1 = 0.0;
					// tmp1 += cellVar_i[id_coeff2] * p_cellVar_bc_i[0];
					// tmp1 += cellVar_i[id_coeff4] * p_cellVar_bc_i[1];
					// tmp1 += cellVar_i[id_coeff5] * p_cellVar_bc_i[2];
					// tmp2 = 0.0;
					// tmp2 += cellVar_i[id_coeff3] * p_cellVar_bc_i[0];
					// tmp2 += cellVar_i[id_coeff5] * p_cellVar_bc_i[1];
					// tmp2 += cellVar_i[id_coeff6] * p_cellVar_bc_i[2];
				
					// auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
					// cellVar_i[oup0] = p_old_cellVar_bc_i[0] + tmp0;
					// cellVar_i[oup1] = p_old_cellVar_bc_i[1] + tmp1;
					// cellVar_i[oup2] = p_old_cellVar_bc_i[2] + tmp2;
					
					// // cout << tmp0 << " " << tmp1 << " " << tmp2 << endl;

					// ++i_bc;
				// }
			// }
		// }
	// }
	
	
	
	
}






void MASCH_Gradient::leastSquare(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
vector<string>& inp_cell, vector<string>& inp_bcFace, 
vector<string>& minmaxInp_cell_name, vector<string>& maxOut_cell_name, vector<string>& minOut_cell_name,
vector<string>& minmaxInp_point_name, vector<string>& maxOut_point_name, vector<string>& minOut_point_name){
	
	// controls.log.push("test1");
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	double tmp_n_weight = n_weight;
	
	int id_coeff1 = controls.getId_cellVar("coeff of least square-(1,1)");
	int id_coeff2 = controls.getId_cellVar("coeff of least square-(1,2)");
	int id_coeff3 = controls.getId_cellVar("coeff of least square-(1,3)");
	int id_coeff4 = controls.getId_cellVar("coeff of least square-(2,2)");
	int id_coeff5 = controls.getId_cellVar("coeff of least square-(2,3)");
	int id_coeff6 = controls.getId_cellVar("coeff of least square-(3,3)");

	int inp_size = inp_cell.size();
	int inp_minmax_size = minmaxInp_cell_name.size();

	vector<string> vec_oup0, vec_oup1, vec_oup2, vec_bc_face;
	{
		int tmp_iter = 0;
		for(auto& item : inp_cell){
			string tmp_oup0 = "x-gradient "; tmp_oup0+=item;
			string tmp_oup1 = "y-gradient "; tmp_oup1+=item;
			string tmp_oup2 = "z-gradient "; tmp_oup2+=item;
			string tmp_bc_face = inp_bcFace[tmp_iter];
			
			vec_oup0.push_back(tmp_oup0);
			vec_oup1.push_back(tmp_oup1);
			vec_oup2.push_back(tmp_oup2);
			vec_bc_face.push_back(tmp_bc_face);
			++tmp_iter;
		}
	}
	
	
	vector<int> id_inp, id_oup0, id_oup1, id_oup2, id_bc_face;
	
	int iter_tmp = 0;
	for(auto& item : inp_cell){
		int inp = controls.getId_cellVar(item);
		int oup0 = controls.getId_cellVar(vec_oup0[iter_tmp]);
		int oup1 = controls.getId_cellVar(vec_oup1[iter_tmp]);
		int oup2 = controls.getId_cellVar(vec_oup2[iter_tmp]);
		int bc_face = controls.getId_faceVar(vec_bc_face[iter_tmp]);
		
		id_inp.push_back(inp);
		id_oup0.push_back(oup0);
		id_oup1.push_back(oup1);
		id_oup2.push_back(oup2);
		id_bc_face.push_back(bc_face);
		
		++iter_tmp;
	}
	
	int iter_tmp2 = 0;
	vector<int> id_inp_minmax, id_oup_max, id_oup_min;
	for(auto& item : minmaxInp_cell_name){
		id_inp_minmax.push_back(controls.getId_cellVar(item));
		id_oup_max.push_back(controls.getId_cellVar(maxOut_cell_name[iter_tmp2]));
		id_oup_min.push_back(controls.getId_cellVar(minOut_cell_name[iter_tmp2]));
		
		++iter_tmp2;
	}
	



	int inp_minmax_point_size = minmaxInp_point_name.size();
	vector<int> id_inp_point_minmax, id_oup_point_max, id_oup_point_min;
    iter_tmp2 = 0;
	for(auto& item : minmaxInp_point_name){
		id_inp_point_minmax.push_back(controls.getId_cellVar(item));
		id_oup_point_max.push_back(controls.getId_pointVar("maximum "+item));
		id_oup_point_min.push_back(controls.getId_pointVar("minimum "+item));
		
		++iter_tmp2;
	}



	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// cout << inp << endl;
	// cout << oup0 << endl;
	// cout << oup1 << endl;
	// cout << oup2 << endl;
	// cout << bc_face << endl;
	// cout << var.cells[0].size() << endl;
	// cout << var.faces[0].size() << endl;

	auto cells = mesh.cells.data();
	auto points = mesh.points.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto pointVar = var.points.data();
	auto faceVar = var.faces.data();
	
	// processor faces
	vector<vector<double>> recv_value(inp_size);
	vector<vector<double>> recv_value2(inp_minmax_size);
	vector<vector<double>> recv_value3(inp_minmax_point_size);
	if(size>1){

		vector<vector<double>> send_value(inp_size);
		vector<vector<double>> send_value2(inp_minmax_size);
		vector<vector<double>> send_value3(inp_minmax_point_size);
		for(int j=0; j<inp_size; ++j){
			send_value[j].reserve(mesh.send_StencilCellsId.size());
		}
		for(int j=0; j<inp_minmax_size; ++j){
			send_value2[j].reserve(mesh.send_StencilCellsId.size());
		}
		for(int j=0; j<inp_minmax_point_size; ++j){
			send_value3[j].reserve(mesh.send_StencilCellsId.size());
		}
		for(auto& icell : mesh.send_StencilCellsId){
			auto cellVar_i = cellVar[icell].data();
			for(int j=0; j<inp_size; ++j){
				send_value[j].push_back(cellVar_i[id_inp[j]]);
			}
			for(int j=0; j<inp_minmax_size; ++j){
				send_value2[j].push_back(cellVar_i[id_inp_minmax[j]]);
			}
			for(int j=0; j<inp_minmax_point_size; ++j){
				send_value3[j].push_back(cellVar_i[id_inp_point_minmax[j]]);
			}
		}
		for(int j=0; j<inp_size; ++j){
			recv_value[j].resize(mesh.recv_displsStencilCells[size]);
			MPI_Alltoallv( send_value[j].data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
						   recv_value[j].data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		for(int j=0; j<inp_minmax_size; ++j){
			recv_value2[j].resize(mesh.recv_displsStencilCells[size]);
			MPI_Alltoallv( send_value2[j].data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
						   recv_value2[j].data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		for(int j=0; j<inp_minmax_point_size; ++j){
			recv_value3[j].resize(mesh.recv_displsStencilCells[size]);
			MPI_Alltoallv( send_value3[j].data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
						   recv_value3[j].data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		
	}

	
	
    //=================================
    // 포인트 min max
	for(int i=0, SIZE=mesh.points.size(); i<SIZE; ++i){
		auto& point = points[i];
		auto pointVar_i = pointVar[i].data();
		
		double maxVal[inp_minmax_point_size];
		double minVal[inp_minmax_point_size];
		for(int j=0; j<inp_minmax_point_size; ++j){
			maxVal[j] = -1.e200;
			minVal[j] = 1.e200;
		}
		
		for(auto& icell : point.iStencils){
			auto& cellSten = cells[icell];
			auto cellStenVar_i = cellVar[icell].data();
			
			// min & max
			for(int j=0; j<inp_minmax_point_size; ++j){
				double tmp_val = cellStenVar_i[id_inp_point_minmax[j]];
				maxVal[j] = max(tmp_val,maxVal[j]);
				minVal[j] = min(tmp_val,minVal[j]);
			}
		}
		
		for(auto& icell : point.recv_iStencils){
			// min & max
			for(int j=0; j<inp_minmax_point_size; ++j){
				double tmp_val = recv_value3[j][icell];
				maxVal[j] = max(tmp_val,maxVal[j]);
				minVal[j] = min(tmp_val,minVal[j]);
			}
		}

		// min & max
		for(int j=0; j<inp_minmax_point_size; ++j){
			pointVar_i[id_oup_point_max[j]] = maxVal[j];
			pointVar_i[id_oup_point_min[j]] = minVal[j];
		}
			
		
	}
    //=================================
		
    
    
    
    
	// controls.log.pop();
	// controls.log.push("test2");
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(auto& item : id_oup0){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup1){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup2){
			cellVar_i[item] = 0.0;
		}
		
		// 리미터 초기화
		for(auto& lim_phi_id : controls.limiterNamesForUnst){
			cellVar_i[lim_phi_id] = 1.0;
		}
	}
    
    
    
    
    
	
	// controls.log.pop();
	// controls.log.push("test3");
	
	// vector<double> cell_var(id_oup0.size(),0.0);
	// vector<double> cell_var(id_oup0.size(),0.0);
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();
		
		double maxVal[inp_minmax_size];
		double minVal[inp_minmax_size];
		for(int j=0; j<inp_minmax_size; ++j){
			maxVal[j] = cellVar_i[id_inp_minmax[j]];
			minVal[j] = cellVar_i[id_inp_minmax[j]];
		}
		
		double cell_x = cell.x;
		double cell_y = cell.y;
		double cell_z = cell.z;
		double cell_var[inp_size];
		double* p_oup0[inp_size];
		double* p_oup1[inp_size];
		double* p_oup2[inp_size];
	// controls.log.push("test2");
		for(int j=0; j<inp_size; ++j){
			cell_var[j] = cellVar_i[id_inp[j]];
			p_oup0[j] = &cellVar_i[id_oup0[j]];
			p_oup1[j] = &cellVar_i[id_oup1[j]];
			p_oup2[j] = &cellVar_i[id_oup2[j]];
		}
		for(auto& icell : cell.iStencils){
		// cout << icell << endl;
			auto& cellSten = cells[icell];
			auto cellStenVar_i = cellVar[icell].data();
			
			// min & max
			for(int j=0; j<inp_minmax_size; ++j){
				double tmp_val = cellStenVar_i[id_inp_minmax[j]];
				maxVal[j] = max(tmp_val,maxVal[j]);
				minVal[j] = min(tmp_val,minVal[j]);
			}
			
			double distX = cellSten.x - cell_x;
			double distY = cellSten.y - cell_y;
			double distZ = cellSten.z - cell_z;
			
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// int tmptmp = tmp_n_weight;
			// wk = pow(wk,tmptmp);
			
			// vector<double> vari(9,0.0);
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			double DVar[inp_size];
			for(int j=0; j<inp_size; ++j){
				double DVar = cellStenVar_i[id_inp[j]] - cell_var[j];
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
			
		}
		
	// controls.log.pop();
		for(auto& icell : cell.recv_iStencils){
			auto cellSten_x = mesh.recv_x_StencilCells.data();
			auto cellSten_y = mesh.recv_y_StencilCells.data();
			auto cellSten_z = mesh.recv_z_StencilCells.data();
			
			// min & max
			for(int j=0; j<inp_minmax_size; ++j){
				double tmp_val = recv_value2[j][icell];
				maxVal[j] = max(tmp_val,maxVal[j]);
				minVal[j] = min(tmp_val,minVal[j]);
			}
			
			// double distX = recv_x[icell] - cell.x;
			// double distY = recv_y[icell] - cell.y;
			// double distZ = recv_z[icell] - cell.z;
			double distX = cellSten_x[icell] - cell_x;
			double distY = cellSten_y[icell] - cell_y;
			double distZ = cellSten_z[icell] - cell_z;
		
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// wk = pow(wk,n_weight);
			
			// vector<double> vari(9,0.0);;
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			for(int j=0; j<inp_size; ++j){
				double DVar = recv_value[j][icell] - cell_var[j];
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
		}

		// min & max
		for(int j=0; j<inp_minmax_size; ++j){
			cellVar_i[id_oup_max[j]] = maxVal[j];
			cellVar_i[id_oup_min[j]] = minVal[j];
		}
			
		
	}
		
	// controls.log.pop();
	// controls.log.push("test4");
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundaries){
		
		if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = faces[i];
				auto& cell = cells[face.iL];
				auto cellVar_i = cellVar[face.iL].data();
				auto faceVar_i = faceVar[i].data();
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
				// wk = pow(wk,n_weight);
			
				// vector<double> vari(9,0.0);;
				double vari[3];
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				
				for(int j=0; j<inp_size; ++j){
					double DVar = faceVar_i[id_bc_face[j]] - cellVar_i[id_inp[j]];
					cellVar_i[id_oup0[j]] += wk * vari[0] * DVar;
					cellVar_i[id_oup1[j]] += wk * vari[1] * DVar;
					cellVar_i[id_oup2[j]] += wk * vari[2] * DVar;
				}
			}
		}
	}
	
	// controls.log.pop();
	// controls.log.push("test5");
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();

		double tmp0, tmp1, tmp2;
		for(int j=0; j<inp_size; ++j){
			tmp0 = 0.0;
			tmp0 += cellVar_i[id_coeff1] * cellVar_i[id_oup0[j]];
			tmp0 += cellVar_i[id_coeff2] * cellVar_i[id_oup1[j]];
			tmp0 += cellVar_i[id_coeff3] * cellVar_i[id_oup2[j]];
			tmp1 = 0.0;
			tmp1 += cellVar_i[id_coeff2] * cellVar_i[id_oup0[j]];
			tmp1 += cellVar_i[id_coeff4] * cellVar_i[id_oup1[j]];
			tmp1 += cellVar_i[id_coeff5] * cellVar_i[id_oup2[j]];
			tmp2 = 0.0;
			tmp2 += cellVar_i[id_coeff3] * cellVar_i[id_oup0[j]];
			tmp2 += cellVar_i[id_coeff5] * cellVar_i[id_oup1[j]];
			tmp2 += cellVar_i[id_coeff6] * cellVar_i[id_oup2[j]];
			
			cellVar_i[id_oup0[j]] = tmp0;
			cellVar_i[id_oup1[j]] = tmp1;
			cellVar_i[id_oup2[j]] = tmp2;
		}
		
	}
	
	// controls.log.pop();
	
	// // 바운더리 커렉션
	// int maxIbc = 2;
	
	// int nBoundary = 0;
	// for(auto& boundary : mesh.boundaries){
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// nBoundary += boundary.nFaces;
		// }
	// }
	
	// vector<vector<double>> old_cellVar_bc(nBoundary,vector<double>(3,0.0));
	// auto p_old_cellVar_bc = old_cellVar_bc.data();
	// int i_bc = 0;
	// for(auto& boundary : mesh.boundaries){
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// auto& cell = cells[face.iL];
				// auto cellVar_i = cellVar[face.iL].data();
				// auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
				// p_old_cellVar_bc_i[0] = cellVar_i[oup0];
				// p_old_cellVar_bc_i[1] = cellVar_i[oup1];
				// p_old_cellVar_bc_i[2] = cellVar_i[oup2];
				// ++i_bc;
			// }
		// }
	// }
	
	// for(int ibc=0; ibc<maxIbc; ++ibc)
	// {
		// vector<vector<double>> cellVar_bc(nBoundary,vector<double>(3,0.0));
		// auto p_cellVar_bc = cellVar_bc.data();
		// i_bc = 0;
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// if(primId_there==true){
					// type_bc.clear();
					// type_bc = boundary.types[id_primOrder];
				// }
				
				// if(type_bc!="zeroGradient") continue;
		
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = faces[i];
					// auto& cell = cells[face.iL];
					// auto cellVar_i = cellVar[face.iL].data();
						
					// double distX = face.x - cell.x;
					// double distY = face.y - cell.y;
					// double distZ = face.z - cell.z;
							
					// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
					// // wk = pow(wk,n_weight);
				
					// // vector<double> vari(9,0.0);;
					// double vari[3];
					// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
					
					// // for(int iEq=0; iEq<nEq; ++iEq){
					// // }
					
					// double DVar = 0.0;
					// DVar += cellVar_i[oup0]*distX;
					// DVar += cellVar_i[oup1]*distY;
					// DVar += cellVar_i[oup2]*distZ;
					
					// // double resi[3];
					// // resi[0] = wk * vari[0] * DVar;
					// // resi[1] = wk * vari[1] * DVar;
					// // resi[2] = wk * vari[2] * DVar;
							
					// auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// p_cellVar_bc_i[0] += wk * vari[0] * DVar;
					// p_cellVar_bc_i[1] += wk * vari[1] * DVar;
					// p_cellVar_bc_i[2] += wk * vari[2] * DVar;
					// ++i_bc;
				// }
			// }
		// }
		
		
		// i_bc = 0;
		// for(auto& boundary : mesh.boundaries){
			// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
				// if(primId_there==true){
					// type_bc.clear();
					// type_bc = boundary.types[id_primOrder];
				// }
				
				// if(type_bc!="zeroGradient") continue;
		
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = faces[i];
					// auto& cell = cells[face.iL];
					// auto cellVar_i = cellVar[face.iL].data();
					
					// auto p_cellVar_bc_i = p_cellVar_bc[i_bc].data();
					// double tmp0, tmp1, tmp2;
					// tmp0 = 0.0;
					// tmp0 += cellVar_i[id_coeff1] * p_cellVar_bc_i[0];
					// tmp0 += cellVar_i[id_coeff2] * p_cellVar_bc_i[1];
					// tmp0 += cellVar_i[id_coeff3] * p_cellVar_bc_i[2];
					// tmp1 = 0.0;
					// tmp1 += cellVar_i[id_coeff2] * p_cellVar_bc_i[0];
					// tmp1 += cellVar_i[id_coeff4] * p_cellVar_bc_i[1];
					// tmp1 += cellVar_i[id_coeff5] * p_cellVar_bc_i[2];
					// tmp2 = 0.0;
					// tmp2 += cellVar_i[id_coeff3] * p_cellVar_bc_i[0];
					// tmp2 += cellVar_i[id_coeff5] * p_cellVar_bc_i[1];
					// tmp2 += cellVar_i[id_coeff6] * p_cellVar_bc_i[2];
				
					// auto p_old_cellVar_bc_i = p_old_cellVar_bc[i_bc].data();
					// cellVar_i[oup0] = p_old_cellVar_bc_i[0] + tmp0;
					// cellVar_i[oup1] = p_old_cellVar_bc_i[1] + tmp1;
					// cellVar_i[oup2] = p_old_cellVar_bc_i[2] + tmp2;
					
					// // cout << tmp0 << " " << tmp1 << " " << tmp2 << endl;

					// ++i_bc;
				// }
			// }
		// }
	// }
	
	
	
	
}
















void MASCH_Gradient::leastSquare_zeroGradient(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, 
vector<string>& inp_cell){
	
	// controls.log.push("test1");
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	double tmp_n_weight = n_weight;
	
	int id_coeff1 = controls.getId_cellVar("coeff of least square-(1,1)");
	int id_coeff2 = controls.getId_cellVar("coeff of least square-(1,2)");
	int id_coeff3 = controls.getId_cellVar("coeff of least square-(1,3)");
	int id_coeff4 = controls.getId_cellVar("coeff of least square-(2,2)");
	int id_coeff5 = controls.getId_cellVar("coeff of least square-(2,3)");
	int id_coeff6 = controls.getId_cellVar("coeff of least square-(3,3)");
	
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// if(var.cells[i][id_coeff1]!=0.0) cout << var.cells[i][id_coeff1] << endl;
		// if(var.cells[i][id_coeff2]!=0.0) cout << var.cells[i][id_coeff2] << endl;
		// if(var.cells[i][id_coeff3]!=0.0) cout << var.cells[i][id_coeff3] << endl;
		// if(var.cells[i][id_coeff4]!=0.0) cout << var.cells[i][id_coeff4] << endl;
		// if(var.cells[i][id_coeff5]!=0.0) cout << var.cells[i][id_coeff5] << endl;
		// if(var.cells[i][id_coeff6]!=0.0) cout << var.cells[i][id_coeff6] << endl;
		
	// }
	

	int inp_size = inp_cell.size();

	vector<string> vec_oup0, vec_oup1, vec_oup2, vec_bc_face;
	{
		int tmp_iter = 0;
		for(auto& item : inp_cell){
			string tmp_oup0 = "x-gradient "; tmp_oup0+=item;
			string tmp_oup1 = "y-gradient "; tmp_oup1+=item;
			string tmp_oup2 = "z-gradient "; tmp_oup2+=item;
			// string tmp_bc_face = inp_bcFace[tmp_iter];
			
			vec_oup0.push_back(tmp_oup0);
			vec_oup1.push_back(tmp_oup1);
			vec_oup2.push_back(tmp_oup2);
			// vec_bc_face.push_back(tmp_bc_face);
			++tmp_iter;
		}
	}
	
	
	vector<int> id_inp, id_oup0, id_oup1, id_oup2, id_bc_face;
	
	int iter_tmp = 0;
	for(auto& item : inp_cell){
		int inp = controls.getId_cellVar(item);
		int oup0 = controls.getId_cellVar(vec_oup0[iter_tmp]);
		int oup1 = controls.getId_cellVar(vec_oup1[iter_tmp]);
		int oup2 = controls.getId_cellVar(vec_oup2[iter_tmp]);
		// int bc_face = controls.getId_faceVar(vec_bc_face[iter_tmp]);
		
		// cout << item << endl;
		id_inp.push_back(inp);
		id_oup0.push_back(oup0);
		id_oup1.push_back(oup1);
		id_oup2.push_back(oup2);
		// id_bc_face.push_back(bc_face);
		
		++iter_tmp;
	}

	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// cout << inp << endl;
	// cout << oup0 << endl;
	// cout << oup1 << endl;
	// cout << oup2 << endl;
	// cout << bc_face << endl;
	// cout << var.cells[0].size() << endl;
	// cout << var.faces[0].size() << endl;

	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	// processor faces
	vector<vector<double>> recv_value(inp_size);
	if(size>1){

		vector<vector<double>> send_value(inp_size);
		for(int j=0; j<inp_size; ++j){
			send_value[j].reserve(mesh.send_StencilCellsId.size());
		}
		for(auto& icell : mesh.send_StencilCellsId){
			auto cellVar_i = cellVar[icell].data();
			for(int j=0; j<inp_size; ++j){
				send_value[j].push_back(cellVar_i[id_inp[j]]);
			}
		}
		for(int j=0; j<inp_size; ++j){
			recv_value[j].resize(mesh.recv_displsStencilCells[size]);
			MPI_Alltoallv( send_value[j].data(), mesh.send_countsStencilCells.data(), mesh.send_displsStencilCells.data(), MPI_DOUBLE, 
						   recv_value[j].data(), mesh.recv_countsStencilCells.data(), mesh.recv_displsStencilCells.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		
	}

	
	
	// controls.log.pop();
	// controls.log.push("test2");
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(auto& item : id_oup0){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup1){
			cellVar_i[item] = 0.0;
		}
		for(auto& item : id_oup2){
			cellVar_i[item] = 0.0;
		}
	}
	
	// controls.log.pop();
	// controls.log.push("test3");
	
	// vector<double> cell_var(id_oup0.size(),0.0);
	// vector<double> cell_var(id_oup0.size(),0.0);
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();
		double cell_x = cell.x;
		double cell_y = cell.y;
		double cell_z = cell.z;
		double cell_var[inp_size];
		double* p_oup0[inp_size];
		double* p_oup1[inp_size];
		double* p_oup2[inp_size];
	// controls.log.push("test2");
		for(int j=0; j<inp_size; ++j){
			cell_var[j] = cellVar_i[id_inp[j]];
			p_oup0[j] = &cellVar_i[id_oup0[j]];
			p_oup1[j] = &cellVar_i[id_oup1[j]];
			p_oup2[j] = &cellVar_i[id_oup2[j]];
			// cout << cellVar_i[id_inp[j]] << endl;
			// if(cell_var[j]!=0.0) cout << cell_var[j] << endl;
		}
		for(auto& icell : cell.iStencils){
		// cout << icell << endl;
			auto& cellSten = cells[icell];
			auto cellStenVar_i = cellVar[icell].data();
			
			double distX = cellSten.x - cell_x;
			double distY = cellSten.y - cell_y;
			double distZ = cellSten.z - cell_z;
			
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// int tmptmp = tmp_n_weight;
			// wk = pow(wk,tmptmp);
			
			// vector<double> vari(9,0.0);
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			double DVar[inp_size];
			for(int j=0; j<inp_size; ++j){
				double DVar = cellStenVar_i[id_inp[j]] - cell_var[j];
				// if(DVar!=0.0) cout << DVar << endl;
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
			
		}
		
	// controls.log.pop();
		for(auto& icell : cell.recv_iStencils){
			auto cellSten_x = mesh.recv_x_StencilCells.data();
			auto cellSten_y = mesh.recv_y_StencilCells.data();
			auto cellSten_z = mesh.recv_z_StencilCells.data();
			// double distX = recv_x[icell] - cell.x;
			// double distY = recv_y[icell] - cell.y;
			// double distZ = recv_z[icell] - cell.z;
			double distX = cellSten_x[icell] - cell_x;
			double distY = cellSten_y[icell] - cell_y;
			double distZ = cellSten_z[icell] - cell_z;
		
			double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
			// wk = pow(wk,n_weight);
			
			// vector<double> vari(9,0.0);;
			double vari[3];
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
		
			for(int j=0; j<inp_size; ++j){
				double DVar = recv_value[j][icell] - cell_var[j];
				*p_oup0[j] += wk * vari[0] * DVar;
				*p_oup1[j] += wk * vari[1] * DVar;
				*p_oup2[j] += wk * vari[2] * DVar;
			}
		}
	}
		
	// controls.log.pop();
	// controls.log.push("test4");
	
	// // boundary face's nodes
	// for(auto& boundary : mesh.boundaries){
		
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// auto& cell = cells[face.iL];
				// auto cellVar_i = cellVar[face.iL].data();
				// auto faceVar_i = faceVar[i].data();
					
				// double distX = face.x - cell.x;
				// double distY = face.y - cell.y;
				// double distZ = face.z - cell.z;
				
				// double wk = 1.0 / (distX*distX+distY*distY+distZ*distZ);
				// // wk = pow(wk,n_weight);
			
				// // vector<double> vari(9,0.0);;
				// double vari[3];
				// vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				
				// for(int j=0; j<inp_size; ++j){
					// double DVar = faceVar_i[id_bc_face[j]] - cellVar_i[id_inp[j]];
					// cellVar_i[id_oup0[j]] += wk * vari[0] * DVar;
					// cellVar_i[id_oup1[j]] += wk * vari[1] * DVar;
					// cellVar_i[id_oup2[j]] += wk * vari[2] * DVar;
				// }
			// }
		// }
	// }
	
	// controls.log.pop();
	// controls.log.push("test5");
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = cells[i];
		auto cellVar_i = cellVar[i].data();

		double tmp0, tmp1, tmp2;
		for(int j=0; j<inp_size; ++j){
			tmp0 = 0.0;
			tmp0 += cellVar_i[id_coeff1] * cellVar_i[id_oup0[j]];
			tmp0 += cellVar_i[id_coeff2] * cellVar_i[id_oup1[j]];
			tmp0 += cellVar_i[id_coeff3] * cellVar_i[id_oup2[j]];
			tmp1 = 0.0;
			tmp1 += cellVar_i[id_coeff2] * cellVar_i[id_oup0[j]];
			tmp1 += cellVar_i[id_coeff4] * cellVar_i[id_oup1[j]];
			tmp1 += cellVar_i[id_coeff5] * cellVar_i[id_oup2[j]];
			tmp2 = 0.0;
			tmp2 += cellVar_i[id_coeff3] * cellVar_i[id_oup0[j]];
			tmp2 += cellVar_i[id_coeff5] * cellVar_i[id_oup1[j]];
			tmp2 += cellVar_i[id_coeff6] * cellVar_i[id_oup2[j]];
			
			// if(cellVar_i[id_coeff3]!=0.0) cout << cellVar_i[id_coeff1] << endl;
			
			cellVar_i[id_oup0[j]] = tmp0;
			cellVar_i[id_oup1[j]] = tmp1;
			cellVar_i[id_oup2[j]] = tmp2;
			
		// if(var.cells[i][id_coeff1]!=0.0) cout << var.cells[i][id_coeff1] << endl;
		// if(var.cells[i][id_coeff2]!=0.0) cout << var.cells[i][id_coeff2] << endl;
		// if(var.cells[i][id_coeff3]!=0.0) cout << var.cells[i][id_coeff3] << endl;
		// if(var.cells[i][id_coeff4]!=0.0) cout << var.cells[i][id_coeff4] << endl;
		// if(var.cells[i][id_coeff5]!=0.0) cout << var.cells[i][id_coeff5] << endl;
		// if(var.cells[i][id_coeff6]!=0.0) cout << var.cells[i][id_coeff6] << endl;
		}
		
	}
	
	
	
	
}