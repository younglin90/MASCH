
#include "./controls.h"

void MASCH_Control::setVariablesBasic(){
	
	string scal = "scalar";
	string vec = "vector";
	string cell = "cell";
	string face = "face";
	string point = "point";
	string field = "field";
	
	(*this).setVarible({field},"time","t","sec.","physics",scal);
	(*this).setVarible({field},"time-step","dt","sec.","physics",scal);
	
	(*this).setVarible({face},"area","dA","m^2","mesh",scal);
	(*this).setVarible({face},"unit normal vector","nvec","","mesh",vec,
						{"x unit normal","y unit normal","z unit normal"},{"nx","ny","nz"},
							{"mesh","mesh","mesh"});
	(*this).setVarible({face},"distance vector of between left and right cell",
						"vLR","m","mesh",vec,{
						"x distance of between left and right cell",
						"y distance of between left and right cell",
						"z distance of between left and right cell"},
						{"xLR","yLR","zLR"},
							{"mesh","mesh","mesh"});
	(*this).setVarible({face},"distance weight","Wc","","mesh",scal);
	(*this).setVarible({face},"distance of between left and right cell","dLR","","mesh",scal);
	(*this).setVarible({face},"cosine angle of between face normal and cells","alpha","","mesh",scal);
	(*this).setVarible({face},"skewness vector","vSkew",
						"m","mesh",vec,{
							"x skewness",
							"y skewness",
							"z skewness"},{"xSkew","ySkew","zSkew"},
							{"mesh","mesh","mesh"});
	
	(*this).setVarible({cell},"Courant-Friedrichs-Lewy number","CFL","","numeric",scal);
	(*this).setVarible({cell},"pseudo time-step","dtau","sec.","numeric",scal);
	(*this).setVarible({cell},"volume","V","m^3","mesh",scal);
	(*this).setVarible({cell},"distance vector project to face normal","vCN",
						"m","mesh",vec,{
							"x distance project to face normal",
							"y distance project to face normal",
							"z distance project to face normal"},{"xCN","yCN","zCN"},
							{"mesh","mesh","mesh"});
	
	(*this).setVarible({cell},"coeff of least square","coeffLS",
						"","mesh",vec,{
							"coeff(1,1) of least square",
							"coeff(1,2) of least square",
							"coeff(1,3) of least square",
							"coeff(2,2) of least square",
							"coeff(2,3) of least square",
							"coeff(3,3) of least square"},
							{"coeffLS1","coeffLS2","coeffLS3","coeffLS4","coeffLS5","coeffLS6"},
							{"mesh","mesh","mesh","mesh","mesh","mesh"});
	
}


void MASCH_Control::limitPrim(int iEq, double& value){
	
	double maxPrim = limitMaxPrim[iEq];
	double minPrim = limitMinPrim[iEq];
	if(value>maxPrim) value = maxPrim;
	if(value<minPrim) value = minPrim;
	
}

void MASCH_Control::setPrimitiveValues(){
	auto& controls = (*this);
	
	// primitive 에서 중복되지 않게 벡터와 스칼라만 빼오기
	vector<string> vectorPrim;
	vector<string> tmp_primitive;
	for(auto& [name, tmp_var] : controls.cellVar){
		if(tmp_var.role=="primitive" &&
		tmp_var.shape=="vector"
		){
			vectorPrim.push_back(name);
		}
		if(tmp_var.role=="primitive"){
			tmp_primitive.push_back(name);
		}
	}
	for(auto& prim : vectorPrim){
		for(auto& sub_name : controls.cellVar[prim].sub_name){
			if(find(tmp_primitive.begin(),tmp_primitive.end(),sub_name)!=tmp_primitive.end()){
				tmp_primitive.erase(remove(tmp_primitive.begin(), tmp_primitive.end(), sub_name), tmp_primitive.end());
			}
		}
	}
	
	int iter=0;
	int iter_nPrim=0;
	for(auto& item : tmp_primitive){
		controls.primitiveMap.insert(make_pair(item,iter));
		controls.supPrimVarNames.push_back(item);
		controls.supPrimVarAbbs.push_back(controls.cellVar[item].abb);
		vector<string> tmp_string;
		vector<unsigned short> tmp_unshort;
		if(controls.cellVar[item].shape != "vector") {
			tmp_string.push_back(item);
			tmp_unshort.push_back(controls.cellVar[item].id);
			++iter_nPrim;
		}
		else{
			int iter2 = 0;
			for(auto& item2 : controls.cellVar[item].sub_name){
				if(controls.cellVar[item].sub_role[iter2++] != "primitive")
					continue;
			// cout << item2 << endl;
				tmp_string.push_back(item2);
				tmp_unshort.push_back(controls.cellVar[item2].id);
				++iter_nPrim;
			}
		}
		
		for(auto& item : tmp_string){
			controls.primVarNames.push_back(item);
		}
		for(auto& item : tmp_unshort){
			controls.primVarIds.push_back(item);
		}
		controls.supPrimVarSizes.push_back(tmp_string.size());
		
		++iter;
	}
	controls.nPrim = iter_nPrim;
	controls.nSupPrim = controls.supPrimVarNames.size();
	
	controls.nEq = iter_nPrim;
	
}




void MASCH_Control::setVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int nCells = mesh.cells.size();
	int nFaces = mesh.faces.size();
	int nProcFaces = mesh.nProcessorFaces;
	int nPoints = mesh.points.size();
	int nBoundaries = mesh.boundaries.size();
	// int nParticles = controls.parcel.size();
	int nParticles = 0;
	
	// var.cells = new double*[nCells];
	// var.faces = new double*[nFaces];
	// var.points = new double*[nPoints];
	// // var.edges = new double*[nPoints];
	// var.boundaries = new double*[nBoundaries];
	
	var.cells.resize(nCells);
	var.faces.resize(nFaces);
	var.procRightCells.resize(nProcFaces);
	var.points.resize(nPoints);
	var.boundaries.resize(nBoundaries);
	
	var.parcel.resize(nParticles);
	
	// int nVarCell = controls.cell.varSize;
	// int nVarFace = controls.face.varSize;
	// int nVarPoint = controls.point.varSize;
	// int nVarBoundary = controls.boundary.varSize;
	// int nVarField = controls.field.varSize;
	// int nVarParticle = controls.parcel.varSize;
	int nVarBoundary = 0;
	for(auto& [name, tmp_var] : controls.cellVar){
		if(tmp_var.id>=0 && tmp_var.role=="primitive"){
			// cout << name << endl;
			++nVarBoundary;
		}
	}
	int nVarParticle = 100;
	
	// cout << nVarBoundary << endl;
	
	// var.fields = new double[controls.nFieldVar];
	var.fields.resize(controls.nFieldVar);
	if(controls.nCellVar>0){
		// for(int i=0; i<nCells; ++i) var.cells[i] = new double[controls.nCellVar];
		for(int i=0; i<nCells; ++i) var.cells[i].resize(controls.nCellVar);
		for(int i=0; i<nProcFaces; ++i) var.procRightCells[i].resize(controls.nCellVar);
	}
	if(controls.nFaceVar>0){
		// for(int i=0; i<nFaces; ++i) var.faces[i] = new double[controls.nFaceVar];
		for(int i=0; i<nFaces; ++i) var.faces[i].resize(controls.nFaceVar);
	}
	if(controls.nPointVar>0){
		// for(int i=0; i<nPoints; ++i) var.points[i] = new double[controls.nPointVar];
		for(int i=0; i<nPoints; ++i) var.points[i].resize(controls.nPointVar);
	}
	// for(int i=0; i<nBoundaries; ++i) var.boundaries[i] = new double[nVarBoundary];
	for(int i=0; i<nBoundaries; ++i) var.boundaries[i].resize(nVarBoundary);
	for(int i=0; i<nParticles; ++i) var.parcel[i].resize(nVarParticle);
	
	
	
}


void MASCH_Control::setBoundaryConditions(MASCH_Mesh& mesh){
	auto& controls = (*this);

	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load boundary property files ... ";
	// }
	
	MASCH_Load load;
	
	
	int nPrimitive = controls.nPrim;
	
	// cout << nPrimitive << endl;
	// cout << controls.nSupPrim << endl;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		
		if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			boundary.types.resize(nPrimitive,"");
			boundary.values.resize(nPrimitive,0.0);
			string name = boundary.name;
			{
				int iter=0;
				int iter2=0;
				for(auto& item : controls.boundaryMap){
					string tmp_nametype = name;
					tmp_nametype += ".type";
					string tmp_namevalue = name;
					tmp_namevalue += ".value";
					string type_name = item[tmp_nametype];
					
					// int tmp_size = controls.primVarNames[iter].size();
					int tmp_size = controls.supPrimVarSizes[iter];
					
					if(tmp_size==1){
						boundary.types[iter2] = type_name;
						if(type_name=="fixedValue"){
							boundary.values[iter2] = stod(item[tmp_namevalue]);
						}
						++iter2;
					}
					else{
						vector<string> output = 
						load.extractVector(item[tmp_namevalue]);
						for(int ii=0; ii<tmp_size; ++ii){
							// cout << controls.primVarNames[iter].size() << endl;
							string sub_primName = controls.primVarNames[iter2];
							boundary.types[iter2] = type_name;
							if(type_name=="fixedValue"){
								boundary.values[iter2] = stod(output[ii]);
								// cout << item[tmp_namevalue] << endl;
								// boundary.values[iter2] = stod(item[tmp_namevalue])
							}
							else if(type_name=="slip"){
								
							}
							else if(type_name=="noSlip"){
								
							}
							++iter2;
						}
					}
					++iter;
				}
			}
		}
	}
	
}





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
		for(auto ipoint : face.ipoints){
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
	}
	
	
	
	
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







void MASCH_Control::updateRealTime(MASCH_Variables& var){
	var.fields[fieldVar["time"].id] +=
	var.fields[fieldVar["time-step"].id];
}
bool MASCH_Control::checkContinueRealTimeStep(MASCH_Variables& var){
	time = var.fields[(*this).fieldVar["time"].id];
	if(time>endTime){
		return false;
	}
	return true;
}
bool MASCH_Control::updateRealTime(){
	if(saveControl == "timeStep"){
		if(iterReal % (int)saveInterval == 0){
			return true;
		}
		else{
			return false;
		}
	}
	else if(saveControl == "runTime"){
		int jung = time / (*this).saveInterval;
		double namuji = time - (double)jung * saveInterval;
		if(namuji < timeStep && namuji >= 0.0){
			return true;
		}
		else{
			return false;
		}
	}
	
}

void MASCH_Control::setVarible(vector<string> save_where, string name, string abb, string unit, string role, string shape, vector<string> sub_name, vector<string> sub_abb, vector<string> sub_role){
	
	for(auto& item : save_where){
		if(item=="field"){
			if(shape=="scalar"){
				fieldVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				nFieldVar,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				++nFieldVar;
			}
			else if(shape=="vector"){
				fieldVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					fieldVar.insert(
					make_pair(sub_name[j], MASCH_Control_Variable_Set(
					nFieldVar,sub_name[j],sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nFieldVar;
				}
			}
		}
		if(item=="cell"){
			if(shape=="scalar"){
				cellVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				nCellVar,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				++nCellVar;
			}
			else if(shape=="vector"){
				cellVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					cellVar.insert(
					make_pair(sub_name[j], MASCH_Control_Variable_Set(
					nCellVar,sub_name[j],sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nCellVar;
				}
			}
			
		}
		if(item=="face"){
			if(shape=="scalar"){
				faceVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				nFaceVar,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				++nFaceVar;
			}
			else if(shape=="vector"){
				faceVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					faceVar.insert(
					make_pair(sub_name[j], MASCH_Control_Variable_Set(
					nFaceVar,sub_name[j],sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nFaceVar;
				}
			}
		}
		if(item=="point"){
			if(shape=="scalar"){
				pointVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				nPointVar,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				++nPointVar;
			}
			else if(shape=="vector"){
				pointVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					pointVar.insert(
					make_pair(sub_name[j], MASCH_Control_Variable_Set(
					nPointVar,sub_name[j],sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nPointVar;
				}
			}
		}
	}
	
	
}



