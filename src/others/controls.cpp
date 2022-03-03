
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
	(*this).setVarible({field},"residual","","","physics",scal);
	
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
	(*this).setVarible({face},"Courant-Friedrichs-Lewy number","CFL","","numeric",scal);
	
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
	
	// var.nCells = nCells;
	
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

void MASCH_Control::resetVariableArray(
MASCH_Mesh& mesh, MASCH_Variables& var,
vector<vector<double>>& org_xyz, vector<vector<int>>& cellConn,
string inp_option){
	
	auto& controls = (*this);
	
	int tmp_nPrimitive = 0;
	vector<string> tmp_name_prim;
	vector<string> tmp_name_x_grad_prim;
	vector<string> tmp_name_y_grad_prim;
	vector<string> tmp_name_z_grad_prim;
	for(auto& [name, tmp_var] : controls.cellVar){
		if(tmp_var.id>=0 && tmp_var.role=="primitive"){
			tmp_name_prim.push_back(name);
			string x_grad_name = "x-gradient ";
			string y_grad_name = "y-gradient ";
			string z_grad_name = "z-gradient ";
			x_grad_name += name;
			y_grad_name += name;
			z_grad_name += name;
			tmp_name_x_grad_prim.push_back(x_grad_name);
			tmp_name_y_grad_prim.push_back(y_grad_name);
			tmp_name_z_grad_prim.push_back(z_grad_name);
			++tmp_nPrimitive;
		}
	}
	
	if(inp_option=="refine"){
		
		
		int org_nCells = cellConn.size();
		
		vector<vector<double>> send_val(tmp_nPrimitive);
		vector<vector<double>> send_x_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_y_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_z_grad_val(tmp_nPrimitive);
		for(int i=0; i<org_nCells; ++i){
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				int id_x_grad_prim = controls.getId_cellVar(tmp_name_x_grad_prim[iprim]);
				int id_y_grad_prim = controls.getId_cellVar(tmp_name_y_grad_prim[iprim]);
				int id_z_grad_prim = controls.getId_cellVar(tmp_name_z_grad_prim[iprim]);
				double value = var.cells[i][id_prim];
				send_val[iprim].push_back(value);
				send_x_grad_val[iprim].push_back(var.cells[i][id_x_grad_prim]);
				send_y_grad_val[iprim].push_back(var.cells[i][id_y_grad_prim]);
				send_z_grad_val[iprim].push_back(var.cells[i][id_z_grad_prim]);
			}
		}
		
			
		int nCells = mesh.cells.size();
		int nFaces = mesh.faces.size();
		int nProcFaces = mesh.nProcessorFaces;
		int nPoints = mesh.points.size();
		int nBoundaries = mesh.boundaries.size();
		// int nParticles = controls.parcel.size();
		int nParticles = 0;
		
		var.cells.resize(nCells);
		var.faces.resize(nFaces);
		var.procRightCells.resize(nProcFaces);
		var.points.resize(nPoints);
		var.boundaries.resize(nBoundaries);
		
		var.parcel.resize(nParticles);
		
		int nVarBoundary = 0;
		for(auto& [name, tmp_var] : controls.cellVar){
			if(tmp_var.id>=0 && tmp_var.role=="primitive"){
				// cout << name << endl;
				++nVarBoundary;
			}
		}
		int nVarParticle = 100;
		
		var.fields.resize(controls.nFieldVar);
		if(controls.nCellVar>0){
			for(int i=0; i<nCells; ++i) var.cells[i].resize(controls.nCellVar);
			for(int i=0; i<nProcFaces; ++i) var.procRightCells[i].resize(controls.nCellVar);
		}
		if(controls.nFaceVar>0){
			for(int i=0; i<nFaces; ++i) var.faces[i].resize(controls.nFaceVar);
		}
		if(controls.nPointVar>0){
			for(int i=0; i<nPoints; ++i) var.points[i].resize(controls.nPointVar);
		}
		for(int i=0; i<nBoundaries; ++i) var.boundaries[i].resize(nVarBoundary);
		for(int i=0; i<nParticles; ++i) var.parcel[i].resize(nVarParticle);
	
		
		
		for(int i=0; i<org_nCells; ++i){
			int sub_size = cellConn[i].size();
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				double org_value = send_val[iprim][i];
				double cell_x = org_xyz[i][0];
				double cell_y = org_xyz[i][1];
				double cell_z = org_xyz[i][2];
				double x_grad = send_x_grad_val.at(iprim).at(i);
				double y_grad = send_y_grad_val.at(iprim).at(i);
				double z_grad = send_z_grad_val.at(iprim).at(i);
				double avg_values = 0.0;
				for(int j=0; j<sub_size; ++j){
					int id_cell = cellConn.at(i).at(j);
					double value_interpol = org_value;
					// value_interpol += x_grad*(mesh.cells.at(id_cell).x - cell_x);
					// value_interpol += y_grad*(mesh.cells.at(id_cell).y - cell_y);
					// value_interpol += z_grad*(mesh.cells.at(id_cell).z - cell_z);
					
					var.cells.at(id_cell).at(id_prim) = value_interpol;
					
					avg_values += value_interpol/(double)sub_size;
				}
				// if(avg_values!=org_value){
					// double coeff = sub_size*org_value/(avg_values+1.e-12);
					// for(int j=0; j<sub_size; ++j){
						// int id_cell = cellConn[i][j];
						// var.cells.at(id_cell).at(id_prim) *= coeff;
					// }
				// }
			}
		}
		

			
		
		
	}
	else if(inp_option=="repart"){
		
	}
	else if(inp_option=="unrefine"){
		
		// cout << mesh.cells.size() << " " <<
		// cellConn.size() << " " <<
		// cellConn.back().back() << " " << 
		// var.cells.size() <<
		// endl;
		
		vector<vector<double>> send_val(tmp_nPrimitive);
		vector<vector<double>> send_x_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_y_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_z_grad_val(tmp_nPrimitive);
		for(auto& item : cellConn){
			for(auto& i : item){
				for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
					int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
					int id_x_grad_prim = controls.getId_cellVar(tmp_name_x_grad_prim[iprim]);
					int id_y_grad_prim = controls.getId_cellVar(tmp_name_y_grad_prim[iprim]);
					int id_z_grad_prim = controls.getId_cellVar(tmp_name_z_grad_prim[iprim]);
					double value = var.cells[i][id_prim];
					send_val[iprim].push_back(value);
				}
			}
		}
		
			
		int nCells = mesh.cells.size();
		int nFaces = mesh.faces.size();
		int nProcFaces = mesh.nProcessorFaces;
		int nPoints = mesh.points.size();
		int nBoundaries = mesh.boundaries.size();
		// int nParticles = controls.parcel.size();
		int nParticles = 0;
		
		var.cells.resize(nCells);
		var.faces.resize(nFaces);
		var.procRightCells.resize(nProcFaces);
		var.points.resize(nPoints);
		var.boundaries.resize(nBoundaries);
		
		var.parcel.resize(nParticles);
		
		int nVarBoundary = 0;
		for(auto& [name, tmp_var] : controls.cellVar){
			if(tmp_var.id>=0 && tmp_var.role=="primitive"){
				// cout << name << endl;
				++nVarBoundary;
			}
		}
		int nVarParticle = 100;
		
		var.fields.resize(controls.nFieldVar);
		if(controls.nCellVar>0){
			for(int i=0; i<nCells; ++i) var.cells[i].resize(controls.nCellVar);
			for(int i=0; i<nProcFaces; ++i) var.procRightCells[i].resize(controls.nCellVar);
		}
		if(controls.nFaceVar>0){
			for(int i=0; i<nFaces; ++i) var.faces[i].resize(controls.nFaceVar);
		}
		if(controls.nPointVar>0){
			for(int i=0; i<nPoints; ++i) var.points[i].resize(controls.nPointVar);
		}
		for(int i=0; i<nBoundaries; ++i) var.boundaries[i].resize(nVarBoundary);
		for(int i=0; i<nParticles; ++i) var.parcel[i].resize(nVarParticle);
	
		// MPI_Barrier(MPI_COMM_WORLD);
		// cout << "44" << endl;
		
		
		for(int i=0; i<cellConn.size(); ++i){
			int sub_size = cellConn[i].size();
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				
				double new_value = 0.0;
				for(int j=0; j<sub_size; ++j){
					int sub_id = cellConn[i][j];
					double org_value = send_val.at(iprim).at(sub_id);
					new_value += org_value/(double)sub_size;
				}
				
				var.cells.at(i).at(id_prim) = new_value;
			}
		}
		
		
		
	}
	else{
		cout << "#WARNING" << endl;
	}
			
	
	
}


void MASCH_Control::resetVariableArray(
MASCH_Mesh& mesh, MASCH_Variables& var,
vector<int>& cell_ip, vector<int>& cellConn,
string inp_option){
	
	auto& controls = (*this);
	
	
	if(inp_option=="repart"){
		
		int tmp_nPrimitive = 0;
		vector<string> tmp_name_prim;
		vector<string> tmp_name_x_grad_prim;
		vector<string> tmp_name_y_grad_prim;
		vector<string> tmp_name_z_grad_prim;
		for(auto& [name, tmp_var] : controls.cellVar){
			if(tmp_var.id>=0 && tmp_var.role=="primitive"){
				tmp_name_prim.push_back(name);
				string x_grad_name = "x-gradient ";
				string y_grad_name = "y-gradient ";
				string z_grad_name = "z-gradient ";
				x_grad_name += name;
				y_grad_name += name;
				z_grad_name += name;
				tmp_name_x_grad_prim.push_back(x_grad_name);
				tmp_name_y_grad_prim.push_back(y_grad_name);
				tmp_name_z_grad_prim.push_back(z_grad_name);
				++tmp_nPrimitive;
			}
		}
		
		int org_nCells = cellConn.size();
		
		vector<vector<double>> send_val(tmp_nPrimitive);
		vector<vector<double>> send_x_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_y_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_z_grad_val(tmp_nPrimitive);
		for(int i=0; i<org_nCells; ++i){
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				int id_x_grad_prim = controls.getId_cellVar(tmp_name_x_grad_prim[iprim]);
				int id_y_grad_prim = controls.getId_cellVar(tmp_name_y_grad_prim[iprim]);
				int id_z_grad_prim = controls.getId_cellVar(tmp_name_z_grad_prim[iprim]);
				double value = var.cells[i][id_prim];
				send_val[iprim].push_back(value);
				send_x_grad_val[iprim].push_back(var.cells[i][id_x_grad_prim]);
				send_y_grad_val[iprim].push_back(var.cells[i][id_y_grad_prim]);
				send_z_grad_val[iprim].push_back(var.cells[i][id_z_grad_prim]);
			}
		}
		
			
		int nCells = mesh.cells.size();
		int nFaces = mesh.faces.size();
		int nProcFaces = mesh.nProcessorFaces;
		int nPoints = mesh.points.size();
		int nBoundaries = mesh.boundaries.size();
		// int nParticles = controls.parcel.size();
		int nParticles = 0;
		
		var.cells.resize(nCells);
		var.faces.resize(nFaces);
		var.procRightCells.resize(nProcFaces);
		var.points.resize(nPoints);
		var.boundaries.resize(nBoundaries);
		
		var.parcel.resize(nParticles);
		
		int nVarBoundary = 0;
		for(auto& [name, tmp_var] : controls.cellVar){
			if(tmp_var.id>=0 && tmp_var.role=="primitive"){
				// cout << name << endl;
				++nVarBoundary;
			}
		}
		int nVarParticle = 100;
		
		var.fields.resize(controls.nFieldVar);
		if(controls.nCellVar>0){
			for(int i=0; i<nCells; ++i) var.cells[i].resize(controls.nCellVar);
			for(int i=0; i<nProcFaces; ++i) var.procRightCells[i].resize(controls.nCellVar);
		}
		if(controls.nFaceVar>0){
			for(int i=0; i<nFaces; ++i) var.faces[i].resize(controls.nFaceVar);
		}
		if(controls.nPointVar>0){
			for(int i=0; i<nPoints; ++i) var.points[i].resize(controls.nPointVar);
		}
		for(int i=0; i<nBoundaries; ++i) var.boundaries[i].resize(nVarBoundary);
		for(int i=0; i<nParticles; ++i) var.parcel[i].resize(nVarParticle);
	
		
		int rank = MPI::COMM_WORLD.Get_rank();
		int size = MPI::COMM_WORLD.Get_size();
	
		vector<vector<double>> send_value(size);
		for(int i=0; i<org_nCells; ++i){
			int send_id = cellConn[i];
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				double org_value = send_val[iprim][i];
				double x_grad = send_x_grad_val.at(iprim).at(i);
				double y_grad = send_y_grad_val.at(iprim).at(i);
				double z_grad = send_z_grad_val.at(iprim).at(i);
				double avg_values = 0.0;
				
				double value_interpol = org_value;
				
				send_value[cell_ip[i]].push_back(value_interpol);
					
			}
		}
		MASCH_MPI mpi;
		vector<vector<double>> recv_value(size);
		mpi.Alltoallv(send_value, recv_value);
		for(int ip=0, id=0; ip<size; ++ip){
			int tmp_size = recv_value[ip].size()/tmp_nPrimitive;
			for(int i=0, iter=0; i<tmp_size; ++i){
				for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
					int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
					double tmp_value = recv_value[ip][iter++];
					var.cells[id][id_prim] = tmp_value;
				}
				++id;
			}
		}
		
		
	}
			
	
	
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


void MASCH_Control::saveAfterInitial(MASCH_Mesh& mesh){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// 메쉬 파일 로드
	MASCH_Load load;
	load.meshFiles("./grid/0/", controls, mesh);
	// variable들 어레이 생성
	MASCH_Variables var;
	controls.setVariableArray(mesh, var);
	
	// 값 초기화
	int ii=0;
	for(auto& inp_abb : controls.supPrimVarAbbs){
		if(inp_abb=="p") {
			int id = controls.cellVar["pressure"].id;
			double value = stod(controls.initialMap[ii]["value"]);
			int iter = 0;
			for(auto& cell : var.cells){
				cell[id] = value;
				++iter;
			}
		}
		if(inp_abb=="T") {
			int id = controls.cellVar["temperature"].id;
			double value = stod(controls.initialMap[ii]["value"]);
			int iter = 0;
			for(auto& cell : var.cells){
				cell[id] = value;
				++iter;
			}
		}
		if(inp_abb=="U") {
			int id0 = controls.cellVar["x-velocity"].id;
			int id1 = controls.cellVar["y-velocity"].id;
			int id2 = controls.cellVar["z-velocity"].id;
			vector<string> output = load.extractVector(controls.initialMap[ii]["value"]);
			double value0 = stod(output[0]);
			double value1 = stod(output[1]);
			double value2 = stod(output[2]);
			int iter = 0;
			for(auto& cell : var.cells){
				cell[id0] = value0;
				cell[id1] = value1;
				cell[id2] = value2;
				++iter;
			}
		}
		
		++ii;
	}
	// if(controls.initialMap[0]["type"] == "fixedValue"){
		// int id = controls.cellVar["pressure"].id;
		// double value = stod(controls.initialMap[0]["value"]);
		// int iter = 0;
		// for(auto& cell : var.cells){
			// cell[id] = value;
			// ++iter;
		// }
	// }
	
	// if(controls.initialMap[2]["type"] == "fixedValue"){
		// int id0 = controls.cellVar["x-velocity"].id;
		// int id1 = controls.cellVar["y-velocity"].id;
		// int id2 = controls.cellVar["z-velocity"].id;
		// vector<string> output = load.extractVector(controls.initialMap[2]["value"]);
		// double value0 = stod(output[0]);
		// double value1 = stod(output[1]);
		// double value2 = stod(output[2]);
		// int iter = 0;
		// for(auto& cell : var.cells){
			// cell[id0] = value0;
			// cell[id1] = value1;
			// cell[id2] = value2;
			// ++iter;
		// }
	// }
	
	// if(controls.initialMap[1]["type"] == "fixedValue"){
		// int id = controls.cellVar["temperature"].id;
		// double value = stod(controls.initialMap[1]["value"]);
		// int iter = 0;
		// for(auto& cell : var.cells){
			// cell[id] = value;
			// ++iter;
		// }
	// }
	
	// // if(controls.initialMap[0]["water.type"] == "fixedValue"){
		// // int id = controls.cellVar["water"].id;
		// // double value = stod(controls.initialMap[0]["water.value"]);
		// // int iter = 0;
		// // for(auto& cell : var.cells){
			// // cell[id] = value;
			// // ++iter;
		// // }
	// // }
	
	
	
	var.fields[controls.fieldVar["time"].id] = 0.0;
	
	
	// save 하기
	MASCH_Mesh_Save save;
	save.fvmFiles("./save/0/", rank, mesh, controls, var);
	// save.vtu("./save/0/", rank, mesh);

}