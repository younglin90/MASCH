
#include "./controls.h"

void MASCH_Control::setVariablesBasic(){
	
	string scal = "scalar";
	string vec3 = "vector3";
	string vec = "vector";
	string cell = "cell";
	string face = "face";
	string point = "point";
	string field = "field";
	
	// 필드 기본값
	(*this).setVarible({field},"time","t","sec.","",scal);
	(*this).setVarible({field},"time-step","dt","sec.","",scal);
	(*this).setVarible({field},"residual","","","",scal);
	
	// 페이스 기본값
	(*this).setVarible({face},"area","dA","m^2","",scal);
	(*this).setVarible({face},"unit normal vector","nvec","","",vec3,
						{"x unit normal","y unit normal","z unit normal"},{"nx","ny","nz"},{"","",""});
	(*this).setVarible({face},"distance vector of between left and right cell",
						"vLR","m","",vec3,{
						"x distance of between left and right cell",
						"y distance of between left and right cell",
						"z distance of between left and right cell"},
						{"xLR","yLR","zLR"},{"","",""});
	(*this).setVarible({face},"distance vector of between cell and face",
						"vLF","m","",vec3,{
						"x distance of between left cell and face",
						"y distance of between left cell and face",
						"z distance of between left cell and face"},
						{"xLF","yLF","zLF"},{"","",""});
	(*this).setVarible({face},"distance vector of between cell and face",
						"vRF","m","",vec3,{
						"x distance of between right cell and face",
						"y distance of between right cell and face",
						"z distance of between right cell and face"},
						{"xRF","yRF","zRF"},{"","",""});
	(*this).setVarible({face},"distance weight","Wc","","",scal);
	(*this).setVarible({face},"distance of between left and right cell","dLR","","",scal);
	(*this).setVarible({face},"cosine angle of between face normal and cells","alpha","","",scal);
	(*this).setVarible({face},"skewness vector","vSkew","m","",vec3,{"x skewness","y skewness","z skewness"},
						{"xSkew","ySkew","zSkew"},{"","",""});
	(*this).setVarible({face},"Courant-Friedrichs-Lewy number","CFL","","",scal);
	
	// 셀 기본값
	(*this).setVarible({cell},"Courant-Friedrichs-Lewy number","CFL","","",scal);
	(*this).setVarible({cell},"pseudo time-step","dtau","sec.","",scal);
	(*this).setVarible({cell},"volume","V","m^3","mesh",scal);
	(*this).setVarible({cell},"distance vector project to face normal","vCN",
						"m","mesh",vec3,{
						"x distance project to face normal",
						"y distance project to face normal",
						"z distance project to face normal"},{"xCN","yCN","zCN"},{"","",""});
	
	(*this).setVarible({cell},"coeff of least square","coeffLS",
						"","mesh",vec,{
						"(1,1)","(1,2)","(1,3)","(2,2)","(2,3)","(3,3)"},
						{"coeffLS1","coeffLS2","coeffLS3","coeffLS4","coeffLS5","coeffLS6"},
						{"","","","","",""});
	(*this).setVarible({cell},"maximum area","maxA","m^2","mesh",scal);
	(*this).setVarible({cell},"minimum area","minA","m^2","mesh",scal);
	
	(*this).setVarible({face},"left cell to face normal vector shortest distant","","m","",vec3,
		{"left cell to face normal vector shortest x distance",
		 "left cell to face normal vector shortest y distance",
		 "left cell to face normal vector shortest z distance"},{"","",""},{"","",""});
	(*this).setVarible({face},"right cell to face normal vector shortest distant","","m","",vec3,
		{"right cell to face normal vector shortest x distance",
		 "right cell to face normal vector shortest y distance",
		 "right cell to face normal vector shortest z distance"},{"","",""},{"","",""});
	
	
	(*this).setVarible({face},"unit normal of between left and right cell",
						"vLR","m","",vec3,{
						"x unit normal of between left and right cell",
						"y unit normal of between left and right cell",
						"z unit normal of between left and right cell"},{"","",""},{"","",""});
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
	vector<string> vector3Prim;
	vector<string> tmp_primitive;
	for(auto& [name, tmp_var] : controls.cellVar){
		if(tmp_var.role=="primitive" &&
		tmp_var.shape=="vector"
		){
			vectorPrim.push_back(name);
			// cout << name << endl;
		}
		else if(tmp_var.role=="primitive" &&
		tmp_var.shape=="vector3"
		){
			vector3Prim.push_back(name);
		}
		else if(tmp_var.role=="primitive"){
			tmp_primitive.push_back(name);
			// cout << name << endl;
		}
	}
	
	controls.primAllScalarNames = tmp_primitive;
	for(auto& item : tmp_primitive){
		controls.primAllScalarIds.push_back(
		controls.getId_cellVar(item)
		);
	}
	controls.primVector3Names = vector3Prim;
	controls.primVectorNames = vectorPrim;
	
	// controls.nEq = tmp_primitive.size();
	// cout << controls.nEq << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	for(auto& prim : vectorPrim){
		for(auto& sub_name : controls.cellVar[prim].sub_name){
			string tmp_name = prim;
			tmp_name += ("-" + sub_name);
			if(find(tmp_primitive.begin(),tmp_primitive.end(),tmp_name)!=tmp_primitive.end()){
				tmp_primitive.erase(remove(tmp_primitive.begin(), tmp_primitive.end(), tmp_name), tmp_primitive.end());
			}
		}
	}
	for(auto& prim : vector3Prim){
		for(auto& sub_name : controls.cellVar[prim].sub_name){
			if(find(tmp_primitive.begin(),tmp_primitive.end(),sub_name)!=tmp_primitive.end()){
				tmp_primitive.erase(remove(tmp_primitive.begin(), tmp_primitive.end(), sub_name), tmp_primitive.end());
			}
		}
	}
	
	controls.primScalarNames = tmp_primitive;
	
	for(auto& item : controls.primScalarNames){
		controls.supPrimNames.push_back(item);
	}
	for(auto& item : controls.primVector3Names){
		controls.supPrimNames.push_back(item);
	}
	for(auto& item : controls.primVectorNames){
		controls.supPrimNames.push_back(item);
	}
	
	
	// int iter=0;
	// int iter_nPrim=0;
	// for(auto& item : tmp_primitive){
		// controls.primitiveMap.insert(make_pair(item,iter));
		// controls.supPrimVarNames.push_back(item);
		// controls.supPrimVarAbbs.push_back(controls.cellVar[item].abb);
		// vector<string> tmp_string;
		// vector<unsigned short> tmp_unshort;
		// if(controls.cellVar[item].shape != "vector") {
			// tmp_string.push_back(item);
			// tmp_unshort.push_back(controls.cellVar[item].id);
			// ++iter_nPrim;
		// }
		// else{
			// int iter2 = 0;
			// for(auto& item2 : controls.cellVar[item].sub_name){
				// if(controls.cellVar[item].sub_role[iter2++] != "primitive")
					// continue;
			// // cout << item2 << endl;
				// tmp_string.push_back(item2);
				// tmp_unshort.push_back(controls.cellVar[item2].id);
				// ++iter_nPrim;
			// }
		// }
		
		// for(auto& item : tmp_string){
			// controls.primVarNames.push_back(item);
		// }
		// for(auto& item : tmp_unshort){
			// controls.primVarIds.push_back(item);
		// }
		// controls.supPrimVarSizes.push_back(tmp_string.size());
		
		// ++iter;
	// }
	// controls.nPrim = iter_nPrim;
	// controls.nSupPrim = controls.supPrimVarNames.size();
	
	// controls.nEq = iter_nPrim;
	
}




void MASCH_Control::setVariableArray(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int nCells = mesh.cells.size();
	int nFaces = mesh.faces.size();
	int nProcFaces = mesh.nProcessorFaces;
	int nPoints = mesh.points.size();
	int nBoundaries = mesh.boundaries.size();
	int nParcels = mesh.parcels.size();
	
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
	var.parcels.resize(nParcels);
	
	
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
	// int nVarParticle = 100;
	
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
	
	for(int i=0; i<nParcels; ++i) var.parcels[i].resize(controls.nParcelVar);
	
	
}

void MASCH_Control::resetVariableArray(
MASCH_Mesh& mesh, MASCH_Variables& var,
vector<vector<double>>& org_xyz, vector<vector<int>>& cellConn,
// vector<int>& interpolRefine_id,
string inp_option){
	
	auto& controls = (*this);
	
	// MASCH_Solver solver;
	
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
	
	bool boolInterpolRefineValues = false;
	if(controls.dynamicMeshMap["AMR.interpolationRefineValues"]=="yes") 
		boolInterpolRefineValues = true;
	
	
	if(inp_option=="refine"){
		
		
		int org_nCells = cellConn.size();
		
		vector<vector<double>> send_val(tmp_nPrimitive);
		vector<vector<double>> send_max_val(tmp_nPrimitive);
		vector<vector<double>> send_min_val(tmp_nPrimitive);
		vector<vector<double>> send_x_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_y_grad_val(tmp_nPrimitive);
		vector<vector<double>> send_z_grad_val(tmp_nPrimitive);
		for(int i=0; i<org_nCells; ++i){
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				int id_max_prim = controls.getId_cellVar("maximum "+tmp_name_prim[iprim]);
				int id_min_prim = controls.getId_cellVar("minimum "+tmp_name_prim[iprim]);
				int id_x_grad_prim = controls.getId_cellVar(tmp_name_x_grad_prim[iprim]);
				int id_y_grad_prim = controls.getId_cellVar(tmp_name_y_grad_prim[iprim]);
				int id_z_grad_prim = controls.getId_cellVar(tmp_name_z_grad_prim[iprim]);
				double value = var.cells[i][id_prim];
				send_val[iprim].push_back(value);
				send_max_val[iprim].push_back(var.cells[i][id_max_prim]);
				send_min_val[iprim].push_back(var.cells[i][id_min_prim]);
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
		int nParcels = mesh.parcels.size();
		
		var.cells.resize(nCells);
		var.faces.resize(nFaces);
		var.procRightCells.resize(nProcFaces);
		var.points.resize(nPoints);
		var.boundaries.resize(nBoundaries);
		var.parcels.resize(nParcels);
		
		// var.parcels.resize(nParticles);
		
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
		for(int i=0; i<nParcels; ++i) var.parcels[i].resize(controls.nParcelVar);
	
		
		
		for(int i=0; i<org_nCells; ++i){
			int sub_size = cellConn[i].size();
			for(int iprim=0; iprim<tmp_nPrimitive; ++iprim){
				int id_prim = controls.getId_cellVar(tmp_name_prim[iprim]);
				double org_value = send_val[iprim][i];
				double org_max_value = send_max_val[iprim][i];
				double org_min_value = send_min_val[iprim][i];
				double cell_x = org_xyz[i][0];
				double cell_y = org_xyz[i][1];
				double cell_z = org_xyz[i][2];
				double x_grad = send_x_grad_val.at(iprim).at(i);
				double y_grad = send_y_grad_val.at(iprim).at(i);
				double z_grad = send_z_grad_val.at(iprim).at(i);
				double avg_values = 0.0;
				for(int j=0; j<sub_size; ++j){
					int id_cell = cellConn.at(i).at(j);
					double Delta_minus = 0.0;
					Delta_minus += x_grad*(mesh.cells.at(id_cell).x - cell_x);
					Delta_minus += y_grad*(mesh.cells.at(id_cell).y - cell_y);
					Delta_minus += z_grad*(mesh.cells.at(id_cell).z - cell_z);
					double eta = 1.e-16;
					// double limGrad = solver.limiter_MLP(org_value,org_max_value,org_min_value,gradInter_val, eta);
					double phi = org_value;
					double phi_max = org_max_value;
					double phi_min = org_min_value;
					double limGrad = 1.0;
					// if(Delta_minus>0.0){
						// if(phi+Delta_minus>=phi_max) limGrad = (phi_max-phi)/Delta_minus;
					// }
					// else if(Delta_minus<0.0){
						// if(phi+Delta_minus<=phi_min) limGrad = (phi_min-phi)/Delta_minus;
					// }
					// limGrad = max(0.0,min(1.0,limGrad));
					if(Delta_minus>0.0){
						double Delta_plus = phi_max - phi;
						limGrad= 1.0/Delta_minus*
							((Delta_plus*Delta_plus+eta*eta)*
							Delta_minus+2.0*Delta_minus*Delta_minus*Delta_plus)/
							(Delta_plus*Delta_plus+2.0*Delta_minus*Delta_minus+
							Delta_minus*Delta_plus+eta*eta);
					}
					else if(Delta_minus<0.0){
						double Delta_plus = phi_min - phi;
						limGrad= 1.0/Delta_minus*
							((Delta_plus*Delta_plus+eta*eta)*
							Delta_minus+2.0*Delta_minus*Delta_minus*Delta_plus)/
							(Delta_plus*Delta_plus+2.0*Delta_minus*Delta_minus+
							Delta_minus*Delta_plus+eta*eta);
					}
					
					double value_interpol = org_value;
					if(boolInterpolRefineValues==true) value_interpol += limGrad*Delta_minus;
					var.cells.at(id_cell).at(id_prim) = value_interpol;
					
					avg_values += value_interpol/(double)sub_size;
				}
				// if(avg_values!=org_value){
					// double coeff = 1.0;
					// if(abs(avg_values) > 1.e-16) coeff = org_value/avg_values;
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
		int nParcels = mesh.parcels.size();
		
		var.cells.resize(nCells);
		var.faces.resize(nFaces);
		var.procRightCells.resize(nProcFaces);
		var.points.resize(nPoints);
		var.boundaries.resize(nBoundaries);
		var.parcels.resize(nParcels);
		
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
		for(int i=0; i<nParcels; ++i) var.parcels[i].resize(controls.nParcelVar);
	
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
		
		// var.parcel.resize(nParticles);
		
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
		// for(int i=0; i<nParticles; ++i) var.parcel[i].resize(nVarParticle);
	
		
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
	
	
	// int nPrimitive = controls.nPrim;
	
	

	
	
	// for(int i=0; i<mesh.boundaries.size(); ++i){
		// auto& boundary = mesh.boundaries[i];
		
		// if(boundary.getType()==MASCH_Face_Types::BOUNDARY){
			// boundary.types.resize(nPrimitive,"");
			// boundary.values.resize(nPrimitive,0.0);
			// string name = boundary.name;
			// {
				// for(auto& prim_name : controls.primScalarNames){
					// string type = controls.boundaryMap[prim_name]["type"];
					// int id = controls.getId_cellVar(prim_name);
					
					// vec_inpId.push_back(vector<int>());
					// vec_inpId.back().push_back(id);
					
					// if(type=="fixedValue"){
						// double value = stod(controls.boundaryMap[prim_name]["value"]);
						// calcInitial.push_back(
						// [value, id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
							// cells[id] = value;
							// return 0;
						// });
						// // calcBoundFacePrimVal
					// }
					// else if(type=="function"){
						// string inp_file = controls.boundaryMap[prim_name]["file"];
						// string inp_funct_name = controls.boundaryMap[prim_name]["name"];
						// void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
						// char *error = nullptr;
						// if (handle) {
							// using setFunc_t = int(*)(double, double, double, double, int*, double*);
							// setFunc_t setFunction;
							// *(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
							// calcInitial.push_back(*setFunction);
						// }
						// else{
							// cout << "#WARNING : file not there, " << inp_file << endl;
						// }
					// }
				// }
				
				
				// // int iter=0;
				// // int iter2=0;
				// // for(auto& item : controls.boundaryMap){
					// // string tmp_nametype = name;
					// // tmp_nametype += ".type";
					// // string tmp_namevalue = name;
					// // tmp_namevalue += ".value";
					// // string type_name = item[tmp_nametype];
					
					// // // int tmp_size = controls.primVarNames[iter].size();
					// // int tmp_size = controls.supPrimVarSizes[iter];
					
					// // if(tmp_size==1){
						// // boundary.types[iter2] = type_name;
						// // if(type_name=="fixedValue"){
							// // boundary.values[iter2] = stod(item[tmp_namevalue]);
						// // }
						// // ++iter2;
					// // }
					// // else{
						// // vector<string> output = 
						// // load.extractVector(item[tmp_namevalue]);
						// // for(int ii=0; ii<tmp_size; ++ii){
							// // // cout << controls.primVarNames[iter].size() << endl;
							// // string sub_primName = controls.primVarNames[iter2];
							// // boundary.types[iter2] = type_name;
							// // if(type_name=="fixedValue"){
								// // boundary.values[iter2] = stod(output[ii]);
								// // // cout << item[tmp_namevalue] << endl;
								// // // boundary.values[iter2] = stod(item[tmp_namevalue])
							// // }
							// // else if(type_name=="slip"){
								
							// // }
							// // else if(type_name=="noSlip"){
								
							// // }
							// // ++iter2;
						// // }
					// // }
					// // ++iter;
				// // }
			// }
		// }
	// }
	
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
					string tmp_name = name;
					tmp_name += ("-" + sub_name[j]);
					fieldVar.insert(
					make_pair(tmp_name, MASCH_Control_Variable_Set(
					nFieldVar,tmp_name,sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nFieldVar;
				}
			}
			else{
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
					string tmp_name = name;
					tmp_name += ("-" + sub_name[j]);
					cellVar.insert(
					make_pair(tmp_name, MASCH_Control_Variable_Set(
					nCellVar,tmp_name,sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nCellVar;
				}
			}
			else{
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
					string tmp_name = name;
					tmp_name += ("-" + sub_name[j]);
					faceVar.insert(
					make_pair(tmp_name, MASCH_Control_Variable_Set(
					nFaceVar,tmp_name,sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nFaceVar;
				}
			}
			else{
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
					string tmp_name = name;
					tmp_name += ("-" + sub_name[j]);
					pointVar.insert(
					make_pair(tmp_name, MASCH_Control_Variable_Set(
					nPointVar,tmp_name,sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nPointVar;
				}
			}
			else{
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
		if(item=="parcel"){
			if(shape=="scalar"){
				parcelVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				nParcelVar,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				++nParcelVar;
			}
			else if(shape=="vector"){
				parcelVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					string tmp_name = name;
					tmp_name += ("-" + sub_name[j]);
					parcelVar.insert(
					make_pair(tmp_name, MASCH_Control_Variable_Set(
					nParcelVar,tmp_name,sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nParcelVar;
				}
			}
			else{
				parcelVar.insert(
				make_pair(name, MASCH_Control_Variable_Set(
				-5,name,abb,unit,role,shape,sub_name,sub_abb,sub_role
				)));
				for(int j=0, SIZE=sub_name.size(); j<SIZE; ++j){
					parcelVar.insert(
					make_pair(sub_name[j], MASCH_Control_Variable_Set(
					nParcelVar,sub_name[j],sub_abb[j],unit,sub_role[j],"scalar",{""},{""},{""}
					)));
					++nParcelVar;
				}
			}
		}
	}
	
	
}


void MASCH_Control::saveAfterInitial(MASCH_Mesh& mesh){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	{
		for(auto& point : mesh.points){
			point.level = 0;
		}
		// 중복 포인트, 면 > 3 찾기
		int noAMR_Cells = 0;
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			
			for(auto& p0 : cell.ipoints){
				int p0Num = 0;
				for(auto& iface : cell.ifaces){
					auto& face = mesh.faces[iface];
					auto& facePoints = face.ipoints;
					if( std::find(facePoints.begin(),facePoints.end(),p0)!=facePoints.end() ){
						++p0Num;
					}
				}
				if(p0Num>3){
					cell.level = -1;
				}
			}
			if(cell.level==-1) ++noAMR_Cells;
			
		}
		int noAMR_CellsTot;
		MPI_Allreduce(&noAMR_Cells, &noAMR_CellsTot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		if(rank==0) cout << "| #WARNING : " << noAMR_CellsTot << " cells can NOT AMR" << endl;
	}
	
	
	// 메쉬 파일 로드
	MASCH_Load load;
	// variable들 어레이 생성
	MASCH_Variables var;
	controls.setVariableArray(mesh, var);
	
	// controls.setGeometric(mesh, var);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	vector<vector<int>> vec_inpId;
	
	// 초기화 함수
	using funct_type = 
	function<int(double time, double x, double y, double z, int* inp_id, double* cells)>;
	vector<funct_type> calcInitial;
	
	for(auto& name : controls.primScalarNames){
		string type = controls.initialMap[name]["type"];
		int id = controls.getId_cellVar(name);
		
		vec_inpId.push_back(vector<int>());
		vec_inpId.back().push_back(id);
		
		if(type=="fixedValue"){
			double value = stod(controls.initialMap[name]["value"]);
			calcInitial.push_back(
			[value, id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
				cells[id] = value;
				return 0;
			});
		}
		else if(type=="function"){
			string inp_file = controls.initialMap[name]["file"];
			string inp_funct_name = controls.initialMap[name]["name"];
			void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
			char *error = nullptr;
			if (handle) {
				using setFunc_t = int(*)(double, double, double, double, int*, double*);
				setFunc_t setFunction;
				*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
				calcInitial.push_back(*setFunction);
			}
			else{
				cout << "#WARNING : file not there, " << inp_file << endl;
			}
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(auto& name : controls.primVector3Names){
		string type = controls.initialMap[name]["type"];
		vector<string> sub_names = controls.cellVar[name].sub_name;
		vector<int> sub_id;
		for(auto& item : sub_names){
			sub_id.push_back(controls.getId_cellVar(item));
		}
		
		vec_inpId.push_back(sub_id);
		
		if(type=="fixedValue"){
			vector<string> s_value = load.extractVector(controls.initialMap[name]["value"]);
			vector<double> value;
			for(auto& item : s_value){
				value.push_back(stod(item));
			}
			calcInitial.push_back(
			[value, sub_id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
				cells[sub_id[0]] = value[0];
				cells[sub_id[1]] = value[1];
				cells[sub_id[2]] = value[2];
				return 0;
			});
		}
		else if(type=="function"){
			string inp_file = controls.initialMap[name]["file"];
			string inp_funct_name = controls.initialMap[name]["name"];
			void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
			char *error = nullptr;
			if (handle) {
				using setFunc_t = int(*)(double, double, double, double, int*, double*);
				setFunc_t setFunction;
				*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
				calcInitial.push_back(*setFunction);
			}
			else{
				cout << "#WARNING : file not there, " << inp_file << endl;
			}
		}
	}
	for(auto& sup_name : controls.primVectorNames){
		vector<string> sub_names = controls.cellVar[sup_name].sub_name;
		vector<string> sub_roles = controls.cellVar[sup_name].sub_role;
		int iter=0;
		for(auto& name : sub_names){
			if(sub_roles[iter]!="primitive") continue;
			string type = controls.initialMap[sup_name][name+".type"];
			int id = controls.getId_cellVar(sup_name+"-"+name);
			
			vec_inpId.push_back(vector<int>());
			vec_inpId.back().push_back(id);
		
			if(type=="fixedValue"){
				double value = stod(controls.initialMap[sup_name][name+".value"]);
				calcInitial.push_back(
				[value, id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
					cells[id] = value;
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.initialMap[sup_name][name+".file"];
				string inp_funct_name = controls.initialMap[sup_name][name+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, int*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcInitial.push_back(*setFunction);
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			++iter;
		}
	}
	
	
	var.fields[controls.fieldVar["time"].id] = 0.0;
	{
		int iter=0;
		for(auto& cellVar : var.cells){
			auto& cell = mesh.cells[iter];
			auto cellVar_ptr = cellVar.data();
			int iter2 = 0;
			for(auto& funct : calcInitial){
				funct(0.0,cell.x,cell.y,cell.z,vec_inpId[iter2].data(),cellVar_ptr);
				++iter2;
			}
			++iter;
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// save 하기
	MASCH_Mesh_Save save;
	save.fvmFiles("./save/0/", rank, mesh, controls, var);
	save.fvmFiles_boundary("./save/0/", rank, mesh, controls, var);

}




void MASCH_Control::saveAfterInitialAMR(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
	// 메쉬 파일 로드
	MASCH_Load load;
	vector<vector<int>> vec_inpId;
	
	// 초기화 함수
	using funct_type = 
	function<int(double time, double x, double y, double z, int* inp_id, double* cells)>;
	vector<funct_type> calcInitial;
	
	for(auto& name : controls.primScalarNames){
		string type = controls.initialMap[name]["type"];
		int id = controls.getId_cellVar(name);
		
		vec_inpId.push_back(vector<int>());
		vec_inpId.back().push_back(id);
		
		if(type=="fixedValue"){
			double value = stod(controls.initialMap[name]["value"]);
			calcInitial.push_back(
			[value, id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
				cells[id] = value;
				return 0;
			});
		}
		else if(type=="function"){
			string inp_file = controls.initialMap[name]["file"];
			string inp_funct_name = controls.initialMap[name]["name"];
			void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
			char *error = nullptr;
			if (handle) {
				using setFunc_t = int(*)(double, double, double, double, int*, double*);
				setFunc_t setFunction;
				*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
				calcInitial.push_back(*setFunction);
			}
			else{
				cout << "#WARNING : file not there, " << inp_file << endl;
			}
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(auto& name : controls.primVector3Names){
		string type = controls.initialMap[name]["type"];
		vector<string> sub_names = controls.cellVar[name].sub_name;
		vector<int> sub_id;
		for(auto& item : sub_names){
			sub_id.push_back(controls.getId_cellVar(item));
		}
		
		vec_inpId.push_back(sub_id);
		
		if(type=="fixedValue"){
			vector<string> s_value = load.extractVector(controls.initialMap[name]["value"]);
			vector<double> value;
			for(auto& item : s_value){
				value.push_back(stod(item));
			}
			calcInitial.push_back(
			[value, sub_id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
				cells[sub_id[0]] = value[0];
				cells[sub_id[1]] = value[1];
				cells[sub_id[2]] = value[2];
				return 0;
			});
		}
		else if(type=="function"){
			string inp_file = controls.initialMap[name]["file"];
			string inp_funct_name = controls.initialMap[name]["name"];
			void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
			char *error = nullptr;
			if (handle) {
				using setFunc_t = int(*)(double, double, double, double, int*, double*);
				setFunc_t setFunction;
				*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
				calcInitial.push_back(*setFunction);
			}
			else{
				cout << "#WARNING : file not there, " << inp_file << endl;
			}
		}
	}
	for(auto& sup_name : controls.primVectorNames){
		vector<string> sub_names = controls.cellVar[sup_name].sub_name;
		vector<string> sub_roles = controls.cellVar[sup_name].sub_role;
		int iter=0;
		for(auto& name : sub_names){
			if(sub_roles[iter]!="primitive") continue;
			string type = controls.initialMap[sup_name][name+".type"];
			int id = controls.getId_cellVar(sup_name+"-"+name);
			
			vec_inpId.push_back(vector<int>());
			vec_inpId.back().push_back(id);
		
			if(type=="fixedValue"){
				double value = stod(controls.initialMap[sup_name][name+".value"]);
				calcInitial.push_back(
				[value, id](double time, double x, double y, double z, int* inp_id, double* cells) ->int {
					cells[id] = value;
					return 0;
				});
			}
			else if(type=="function"){
				string inp_file = controls.initialMap[sup_name][name+".file"];
				string inp_funct_name = controls.initialMap[sup_name][name+".name"];
				void *handle = dlopen(inp_file.c_str(), RTLD_NOW);
				char *error = nullptr;
				if (handle) {
					using setFunc_t = int(*)(double, double, double, double, int*, double*);
					setFunc_t setFunction;
					*(void **) (&setFunction) = dlsym(handle, inp_funct_name.c_str());
					calcInitial.push_back(*setFunction);
				}
				else{
					cout << "#WARNING : file not there, " << inp_file << endl;
				}
			}
			++iter;
		}
	}
	
	
	var.fields[controls.fieldVar["time"].id] = 0.0;
	{
		int iter=0;
		for(auto& cellVar : var.cells){
			auto& cell = mesh.cells[iter];
			auto cellVar_ptr = cellVar.data();
			int iter2 = 0;
			for(auto& funct : calcInitial){
				funct(0.0,cell.x,cell.y,cell.z,vec_inpId[iter2].data(),cellVar_ptr);
				++iter2;
			}
			++iter;
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// // save 하기
	// MASCH_Mesh_Save save;
	// save.fvmFiles("./save/0_AMR/", rank, mesh, controls, var);
	// save.fvmFiles_boundary("./save/0_AMR/", rank, mesh, controls, var);

}




void MASCH_Control::save_fvmFiles(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Mesh_Save save;
		
	double tmp_time = var.fields[controls.fieldVar["time"].id];
	double tmp_timestep = var.fields[controls.fieldVar["time-step"].id];
	if(controls.saveControl == "timeStep"){
		if((controls.iterReal+1) % controls.saveInTimeStep == 0){
			string foldername;
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			save.fvmFiles(foldername, rank, mesh, controls, var);
			// save.fvmFiles_boundary(foldername, rank, mesh, controls, var);
		}
	}
	else if(controls.saveControl == "runTime"){
		int jung = tmp_time / controls.saveInRunTime;
		double namuji = tmp_time - (double)jung * controls.saveInRunTime;
		if(namuji < tmp_timestep && namuji >= 0.0){
			string foldername;
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			save.fvmFiles(foldername, rank, mesh, controls, var);
			// save.fvmFiles_boundary(foldername, rank, mesh, controls, var);
		}
	}
	
}




void MASCH_Control::save_dpmFiles(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	MASCH_Mesh_Save save;
		
	double tmp_time = var.fields[controls.fieldVar["time"].id];
	double tmp_timestep = var.fields[controls.fieldVar["time-step"].id];
	if(controls.saveControl == "timeStep"){
		if((controls.iterReal+1) % controls.saveInTimeStep == 0){
			string foldername;
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			// save.fvmFiles(foldername, rank, mesh, controls, var);
			save.parcels(foldername, rank, mesh, controls, var);
			// save.fvmFiles_boundary(foldername, rank, mesh, controls, var);
		}
	}
	else if(controls.saveControl == "runTime"){
		int jung = tmp_time / controls.saveInRunTime;
		double namuji = tmp_time - (double)jung * controls.saveInRunTime;
		if(namuji < tmp_timestep && namuji >= 0.0){
			string foldername;
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			// save.fvmFiles(foldername, rank, mesh, controls, var);
			save.parcels(foldername, rank, mesh, controls, var);
			// save.fvmFiles_boundary(foldername, rank, mesh, controls, var);
		}
	}
	
}


void MASCH_Control::save_pvdFile(MASCH_Mesh& mesh, MASCH_Variables& var){
	auto& controls = (*this);
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// MASCH_Mesh_Save save;
		
	ofstream outputFile;
	
	string foldername;
	bool startSave = false;
	double tmp_time = var.fields[controls.fieldVar["time"].id];
	double tmp_timestep = var.fields[controls.fieldVar["time-step"].id];
	if(controls.saveControl == "timeStep"){
		if((controls.iterReal+1) % controls.saveInTimeStep == 0){
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			startSave = true;
		}
	}
	else if(controls.saveControl == "runTime"){
		int jung = tmp_time / controls.saveInRunTime;
		double namuji = tmp_time - (double)jung * controls.saveInRunTime;
		if(namuji < tmp_timestep && namuji >= 0.0){
			// string foldername;
			std::ostringstream streamObj;
			streamObj << tmp_time;
			streamObj.precision(12);
			foldername = "./save/" + streamObj.str() + "/";
			startSave = true;
		}
	}
	// ==========================================
	// PVD file
	if(rank==0 && startSave==true){
// cout << foldername << endl;
		string filenamePvtu = "./save/plot.pvd";
		string stime = foldername;
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		
		ifstream inputFile;
		inputFile.open(filenamePvtu);
		// if(inputFile && stod(stime)-controls.timeStep != 0.0){
		if(inputFile){

			string nextToken;
			int lineStart = 0;
			while(getline(inputFile, nextToken)){
				if( nextToken.find("</Collection>") != string::npos ){
					break;
				}
				++lineStart;
			}
			inputFile.close();
			
			outputFile.open(filenamePvtu, ios::in);
			outputFile.seekp(-26, ios::end);
			
			// outputFile << saveLines;
			outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			if(controls.nameParcels.size()!=0){
				outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"parcels." << stime << ".pvtu\"/>" << endl;
			}
			outputFile << "  </Collection>" << endl;
			outputFile << "</VTKFile>";
			outputFile.close();
			
		}
		else{
			inputFile.close();
			
			outputFile.open(filenamePvtu);
			
			// string out_line;
			outputFile << "<?xml version=\"1.0\"?>" << endl;
			outputFile << " <VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
			outputFile << "  <Collection>" << endl;
			outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"plot." << stime << ".pvtu\"/>" << endl;
			if(controls.nameParcels.size()!=0){
				outputFile << "    <DataSet timestep=\"" << stime << "\" group=\"\" part=\"0\" file=\"parcels." << stime << ".pvtu\"/>" << endl;
			}
			outputFile << "  </Collection>" << endl;
			outputFile << "</VTKFile>";
			
			
			outputFile.close();
			
		}
		
		
	}
}


bool MASCH_Control::check_isnan(double value){
	if(std::isnan(value) || value < -1.e12 || value > 1.e12){
		return true;
	}
	else{
		return false;
	}
}



void MASCH_Control::show_residual(MASCH_Variables& var){
	auto& controls = (*this);
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	if( (controls.iterReal+1) % stoi(controls.controlDictMap["printLog"]) == 0 ){
		cout << scientific; cout.precision(2);
		if(rank==0) cout << 
		"| iReal = " << controls.iterReal+1 << 
		", " << "iPseudo = " << controls.iterPseudo+1 << 
		", " << "t = " << var.fields[controls.getId_fieldVar("time")] << 
		", " << "dt = " << var.fields[controls.getId_fieldVar("time-step")] << 
		", " << "resi = " << var.fields[controls.getId_fieldVar("residual")] << 
		" |" << endl;
		cout << fixed; cout.precision(0);
		controls.log.show();
	}
}



void MASCH_Control::show_dpm_information(){
	auto& controls = (*this);
	
	if(controls.nameParcels.size()==0) return;
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	int nChangeParcelsE2L, nToProcsRishtParcels, nInsideParcels;
	int nReflectParcels, nEscapeParcels, nDeleteParcels;
	if(rank==0){
		cout << 
		"| nIterDPM = " << controls.nIterDPM <<
		" | EtoL = +" << controls.nChangeParcelsE2L <<
		" | inside = " << controls.nInsideParcels <<
		" | toProc = " << controls.nToProcsRishtParcels <<
		" | reflect = " << controls.nReflectParcels <<
		" | escape = -" << controls.nEscapeParcels <<
		" | delete = -" << controls.nDeleteParcels << 
		" |" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	}
	// if( (controls.iterReal+1) % stoi(controls.controlDictMap["printLog"]) == 0 ){
		// cout << scientific; cout.precision(2);
		// if(rank==0) cout << 
		// "| iReal = " << controls.iterReal+1 << 
		// ", " << "iPseudo = " << controls.iterPseudo+1 << 
		// ", " << "t = " << var.fields[controls.getId_fieldVar("time")] << 
		// ", " << "dt = " << var.fields[controls.getId_fieldVar("time-step")] << 
		// ", " << "resi = " << var.fields[controls.getId_fieldVar("residual")] << 
		// " |" << endl;
		// cout << fixed; cout.precision(0);
		// controls.log.show();
	// }
}
