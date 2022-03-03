
#include "./controls.h"

void MASCH_Control::setGeometric(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	}
	
	int id_volume = controls.cellVar["volume"].id;
	int id_area = controls.faceVar["area"].id;
	int id_nx = controls.faceVar["x unit normal"].id;
	int id_ny = controls.faceVar["y unit normal"].id;
	int id_nz = controls.faceVar["z unit normal"].id;
	int id_xLR = controls.faceVar["x distance of between left and right cell"].id;
	int id_yLR = controls.faceVar["y distance of between left and right cell"].id;
	int id_zLR = controls.faceVar["z distance of between left and right cell"].id;
	int id_xSkew = controls.faceVar["x skewness"].id;
	int id_ySkew = controls.faceVar["y skewness"].id;
	int id_zSkew = controls.faceVar["z skewness"].id;
	int id_xCN = controls.faceVar["x distance project to face normal"].id;
	int id_yCN = controls.faceVar["y distance project to face normal"].id;
	int id_zCN = controls.faceVar["z distance project to face normal"].id;
	int id_Wc = controls.faceVar["distance weight"].id;
	int id_dLR = controls.faceVar["distance of between left and right cell"].id;
	int id_alpha = controls.faceVar["cosine angle of between face normal and cells"].id;
	
	MASCH_Math math;
	
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		
		var.cells[i][id_volume] = 0.0;
		cell.x = 0.0;
		cell.y = 0.0;
		cell.z = 0.0;
	}
	
	// polygon face normal vectors & polygon face area
	// polygon face center x,y,z
	// 3D Polygon Areas
	for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		auto& face = mesh.faces[i];
		vector<double> Vx, Vy, Vz;
		for(auto& ipoint : face.ipoints){
			Vx.push_back(mesh.points[ipoint].x);
			Vy.push_back(mesh.points[ipoint].y);
			Vz.push_back(mesh.points[ipoint].z);
		}
		
		double VSn=0.0;
		vector<double> cellCentroid;
		vector<double> unitNormals(3,0.0);
		double area;
		math.calcUnitNormals_Area3dPolygon(
			face.ipoints.size(), Vx,Vy,Vz,
			unitNormals, area,
			face.x, face.y, face.z,
			VSn, cellCentroid);
			
		var.faces[i][id_nx] = unitNormals[0];
		var.faces[i][id_ny] = unitNormals[1];
		var.faces[i][id_nz] = unitNormals[2];
		var.faces[i][id_area] = area;
			
		int iL = face.iL;
		int iR = face.iR;
			
		// mesh.cells[iL].volume += VSn / 3.0;
		var.cells[iL][id_volume] += VSn / 3.0;
		mesh.cells[iL].x += cellCentroid[0];
		mesh.cells[iL].y += cellCentroid[1];
		mesh.cells[iL].z += cellCentroid[2];
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			// mesh.cells[iR].volume -= VSn / 3.0;
			var.cells[iR][id_volume] -= VSn / 3.0;
			mesh.cells[iR].x -= cellCentroid[0];
			mesh.cells[iR].y -= cellCentroid[1];
			mesh.cells[iR].z -= cellCentroid[2];
		}
			
		
	}
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		cell.x *= 1.0/(24.0*2.0*var.cells[i][id_volume]);
		cell.y *= 1.0/(24.0*2.0*var.cells[i][id_volume]);
		cell.z *= 1.0/(24.0*2.0*var.cells[i][id_volume]);
		
		
		// if(var.cells[i][id_volume]<1.e-30) cout << var.cells[i][id_volume] << endl;
	}
	
	
	// for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		// auto& cell = mesh.cells[i];
		// double avgx = 0.0;
		// double avgy = 0.0;
		// double avgz = 0.0;
		// for(auto& ipoint : cell.ipoints){
			// avgx += mesh.points[ipoint].x;
			// avgy += mesh.points[ipoint].y;
			// avgz += mesh.points[ipoint].z;
		// }
		// avgx /= (double)cell.ipoints.size();
		// avgy /= (double)cell.ipoints.size();
		// avgz /= (double)cell.ipoints.size();
		
		// double resi = 0.0;
		// resi += abs(avgx-cell.x);
		// resi += abs(avgy-cell.y);
		// resi += abs(avgz-cell.z);
		
		// // if(rank==1){
			// if(resi>1.e-12){
				// cout << resi << endl;
			// }
		// // }
	// }
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		auto& face = mesh.faces[i];
		int iL = face.iL;
		int iR = face.iR;
		double nx = var.faces[i][id_nx];
		double ny = var.faces[i][id_ny];
		double nz = var.faces[i][id_nz];
		
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			var.faces[i][id_xLR] = mesh.cells[iR].x - mesh.cells[iL].x;
			var.faces[i][id_yLR] = mesh.cells[iR].y - mesh.cells[iL].y;
			var.faces[i][id_zLR] = mesh.cells[iR].z - mesh.cells[iL].z;
			
			double dxFP = face.x-mesh.cells[iL].x;
			double dyFP = face.y-mesh.cells[iL].y;
			double dzFP = face.z-mesh.cells[iL].z;
			double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
			
			double dxFN = face.x-mesh.cells[iR].x;
			double dyFN = face.y-mesh.cells[iR].y;
			double dzFN = face.z-mesh.cells[iR].z;
			double dFN = sqrt(pow(dxFN,2.0)+pow(dyFN,2.0)+pow(dzFN,2.0));
			
			double dxPN = var.faces[i][id_xLR];
			double dyPN = var.faces[i][id_yLR];
			double dzPN = var.faces[i][id_zLR];
			double dPN = sqrt(pow(dxPN,2.0)+pow(dyPN,2.0)+pow(dzPN,2.0));
			
			// at openfoam 
			double dFP_of = abs(dxFP*nx+dyFP*ny+dzFP*nz);
			double dFN_of = abs(dxFN*nx+dyFN*ny+dzFN*nz);
				
			// at openfoam
			var.faces[i][id_Wc] = dFN_of/(dFP_of+dFN_of);
			// 절단오차 걸러내기 위한 과정
			var.faces[i][id_Wc] = (float)(round(var.faces[i][id_Wc] * 10000000) / 10000000);
			
			var.faces[i][id_dLR] = dPN;
			
			// original
			vector<double> unitNomalsPN(3,0.0);
			unitNomalsPN[0] = dxPN/dPN;
			unitNomalsPN[1] = dyPN/dPN;
			unitNomalsPN[2] = dzPN/dPN;
			double alphaF = 0.0;
			alphaF += nx*unitNomalsPN[0];
			alphaF += ny*unitNomalsPN[1];
			alphaF += nz*unitNomalsPN[2];
			var.faces[i][id_alpha] = 1.0/abs(alphaF);
			
			// skewness
			double D_plane = -(nx*face.x+ny*face.y+nz*face.z);
			double u_line = 
				(nx*mesh.cells[iL].x+
				 ny*mesh.cells[iL].y+
				 nz*mesh.cells[iL].z+
				 D_plane) /
				(nx*(mesh.cells[iL].x-mesh.cells[iR].x)+
				 ny*(mesh.cells[iL].y-mesh.cells[iR].y)+
				 nz*(mesh.cells[iL].z-mesh.cells[iR].z));
			vector<double> vecSkewness(3,0.0);
			vecSkewness[0] = mesh.cells[iL].x-u_line*(mesh.cells[iL].x-mesh.cells[iR].x);
			vecSkewness[1] = mesh.cells[iL].y-u_line*(mesh.cells[iL].y-mesh.cells[iR].y);
			vecSkewness[2] = mesh.cells[iL].z-u_line*(mesh.cells[iL].z-mesh.cells[iR].z);
			vecSkewness[0] = face.x-vecSkewness[0];
			vecSkewness[1] = face.y-vecSkewness[1];
			vecSkewness[2] = face.z-vecSkewness[2];
			
			var.faces[i][id_xSkew] = vecSkewness[0];
			var.faces[i][id_ySkew] = vecSkewness[1];
			var.faces[i][id_zSkew] = vecSkewness[2];
			
			
			
		}
		else if(face.getType() == MASCH_Face_Types::BOUNDARY){
			var.faces[i][id_xLR] = face.x - mesh.cells[iL].x;
			var.faces[i][id_yLR] = face.y - mesh.cells[iL].y;
			var.faces[i][id_zLR] = face.z - mesh.cells[iL].z;
			
			double dxFP = face.x-mesh.cells[iL].x;
			double dyFP = face.y-mesh.cells[iL].y;
			double dzFP = face.z-mesh.cells[iL].z;
			double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
			
			// at openfoam 
			double dFP_of = abs(dxFP*nx+dyFP*ny+dzFP*nz);
				
			// at openfoam
			var.faces[i][id_Wc] = 0.5;
			
			var.faces[i][id_dLR] = dFP;
			
			// original
			vector<double> unitNomalsPN(3,0.0);
			unitNomalsPN[0] = dxFP/dFP;
			unitNomalsPN[1] = dyFP/dFP;
			unitNomalsPN[2] = dzFP/dFP;
			double alphaF = 0.0;
			alphaF += nx*unitNomalsPN[0];
			alphaF += ny*unitNomalsPN[1];
			alphaF += nz*unitNomalsPN[2];
			var.faces[i][id_alpha] = 1.0/abs(alphaF);
			
			// skewness
			var.faces[i][id_xSkew] = 0.0;
			var.faces[i][id_ySkew] = 0.0;
			var.faces[i][id_zSkew] = 0.0;
			
		}
		else if(face.getType() == MASCH_Face_Types::PROCESSOR){
			var.faces[i][id_xLR] = face.x - mesh.cells[iL].x;
			var.faces[i][id_yLR] = face.y - mesh.cells[iL].y;
			var.faces[i][id_zLR] = face.z - mesh.cells[iL].z;
		}
	}
	
	
	
	if(size>1){
		MASCH_MPI mpi;
		
		vector<vector<double>> sendValues(3);
		vector<vector<double>> recvValues(3);
		for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
			auto& face = mesh.faces[i];
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				sendValues[0].push_back(var.faces[i][id_xLR]);
				sendValues[1].push_back(var.faces[i][id_yLR]);
				sendValues[2].push_back(var.faces[i][id_zLR]);
			}
		}
		for(int i=0; i<3; ++i){
			mpi.procFace_Alltoallv(sendValues[i], recvValues[i],
				mesh.countsProcFaces, mesh.countsProcFaces,
				mesh.displsProcFaces, mesh.displsProcFaces);
		}	
		for(int i=0, ip=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL;
			int iR = face.iR;
			double nx = var.faces[i][id_nx];
			double ny = var.faces[i][id_ny];
			double nz = var.faces[i][id_nz];
		
			if(face.getType() == MASCH_Face_Types::PROCESSOR){
				var.faces[i][id_xLR] -= recvValues[0][ip];
				var.faces[i][id_yLR] -= recvValues[1][ip];
				var.faces[i][id_zLR] -= recvValues[2][ip];
				
				double dxFP = face.x-mesh.cells[iL].x;
				double dyFP = face.y-mesh.cells[iL].y;
				double dzFP = face.z-mesh.cells[iL].z;
				double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
				
				double dxFN = recvValues[0][ip];
				double dyFN = recvValues[1][ip];
				double dzFN = recvValues[2][ip];
				double dFN = sqrt(pow(dxFN,2.0)+pow(dyFN,2.0)+pow(dzFN,2.0));
				
				double dxPN = var.faces[i][id_xLR];
				double dyPN = var.faces[i][id_yLR];
				double dzPN = var.faces[i][id_zLR];
				double dPN = sqrt(pow(dxPN,2.0)+pow(dyPN,2.0)+pow(dzPN,2.0));
				
				// at openfoam 
				double dFP_of = abs(dxFP*nx+dyFP*ny+dzFP*nz);
				double dFN_of = abs(dxFN*nx+dyFN*ny+dzFN*nz);
					
				// at openfoam
				var.faces[i][id_Wc] = dFN_of/(dFP_of+dFN_of);
				// 절단오차 걸러내기 위한 과정
				var.faces[i][id_Wc] = (float)(round(var.faces[i][id_Wc] * 10000000) / 10000000);
				
				var.faces[i][id_dLR] = dPN;
				
				// original
				vector<double> unitNomalsPN(3,0.0);
				unitNomalsPN[0] = dxPN/dPN;
				unitNomalsPN[1] = dyPN/dPN;
				unitNomalsPN[2] = dzPN/dPN;
				double alphaF = 0.0;
				alphaF += nx*unitNomalsPN[0];
				alphaF += ny*unitNomalsPN[1];
				alphaF += nz*unitNomalsPN[2];
				var.faces[i][id_alpha] = 1.0/abs(alphaF);
				
				// skewness
				double D_plane = -(nx*face.x+ny*face.y+nz*face.z);
				double u_line = 
					(nx*mesh.cells[iL].x+ny*mesh.cells[iL].y+nz*mesh.cells[iL].z+
					 D_plane) / (nx*(-dxPN)+ny*(-dyPN)+nz*(-dzPN));
				vector<double> vecSkewness(3,0.0);
				vecSkewness[0] = mesh.cells[iL].x-u_line*(-dxPN);
				vecSkewness[1] = mesh.cells[iL].y-u_line*(-dyPN);
				vecSkewness[2] = mesh.cells[iL].z-u_line*(-dzPN);
				vecSkewness[0] = face.x-vecSkewness[0];
				vecSkewness[1] = face.y-vecSkewness[1];
				vecSkewness[2] = face.z-vecSkewness[2];
				
				var.faces[i][id_xSkew] = vecSkewness[0];
				var.faces[i][id_ySkew] = vecSkewness[1];
				var.faces[i][id_zSkew] = vecSkewness[2];
				
				++ip;
			}
		}
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	}
		
}



void MASCH_Control::setGeometricOnlyCell_xyz(MASCH_Mesh& mesh){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	}
	
	MASCH_Math math;
	
	vector<double> tmp_volume(mesh.cells.size(),0.0);
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		cell.x = 0.0;
		cell.y = 0.0;
		cell.z = 0.0;
	}
	
	// polygon face normal vectors & polygon face area
	// polygon face center x,y,z
	// 3D Polygon Areas
	for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
		auto& face = mesh.faces[i];
		vector<double> Vx, Vy, Vz;
		for(auto& ipoint : face.ipoints){
			Vx.push_back(mesh.points[ipoint].x);
			Vy.push_back(mesh.points[ipoint].y);
			Vz.push_back(mesh.points[ipoint].z);
		}
		
		double VSn=0.0;
		vector<double> cellCentroid;
		vector<double> unitNormals(3,0.0);
		double area;
		math.calcUnitNormals_Area3dPolygon(
			face.ipoints.size(), Vx,Vy,Vz,
			unitNormals, area,
			face.x, face.y, face.z,
			VSn, cellCentroid);
			
		int iL = face.iL;
		int iR = face.iR;
			
		// mesh.cells[iL].volume += VSn / 3.0;
		tmp_volume[iL] += VSn / 3.0;
		mesh.cells[iL].x += cellCentroid[0];
		mesh.cells[iL].y += cellCentroid[1];
		mesh.cells[iL].z += cellCentroid[2];
		if(face.getType() == MASCH_Face_Types::INTERNAL){
			// mesh.cells[iR].volume -= VSn / 3.0;
			tmp_volume[iR] -= VSn / 3.0;
			mesh.cells[iR].x -= cellCentroid[0];
			mesh.cells[iR].y -= cellCentroid[1];
			mesh.cells[iR].z -= cellCentroid[2];
		}
			
		
	}
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto& cell = mesh.cells[i];
		cell.x *= 1.0/(24.0*2.0*tmp_volume[i]);
		cell.y *= 1.0/(24.0*2.0*tmp_volume[i]);
		cell.z *= 1.0/(24.0*2.0*tmp_volume[i]);
		
		
		// if(var.cells[i][id_volume]<1.e-30) cout << var.cells[i][id_volume] << endl;
	}
	
}


