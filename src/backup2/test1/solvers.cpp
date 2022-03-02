
#include "./solvers.h"



void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	using Bound_Funct_type = function<int(double* cells, double* faces)>;
	
	vector<string> tmp_type;
	vector<int> tmp_cell_id;
	vector<int> tmp_L_id;
	vector<int> tmp_R_id;
	vector<int> tmp_vel_cell_id;
	vector<int> tmp_vel_L_id;
	vector<int> tmp_vel_R_id;
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		int tmp_size = boundary.types.size();
		int iter=0;
		for(int i=0; i<boundary.types.size(); ++i){
			string type = boundary.types[i];
			int id_inp = controls.primVarIds[i];
			string name = controls.primVarNames[i];
			string left_name = "left ";
			string right_name = "right ";
			left_name += name;
			right_name += name;
			int id_L_out = controls.faceVar[left_name].id;
			int id_R_out = controls.faceVar[right_name].id;
			
			tmp_type.push_back(type);
			tmp_cell_id.push_back(id_inp);
			tmp_L_id.push_back(id_L_out);
			tmp_R_id.push_back(id_R_out);
			
			if(type=="slip" || type=="noSlip"){
				tmp_vel_cell_id.push_back(id_inp);
				tmp_vel_L_id.push_back(id_L_out);
				tmp_vel_R_id.push_back(id_R_out);
			}
				
		}
	}
	
	
	for(int i=0; i<tmp_type.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		string type = tmp_type[i];
		int id_inp = tmp_cell_id[i];
		int id_L_out = tmp_L_id[i];
		int id_R_out = tmp_R_id[i];
		
		Bound_Funct_type setFunction;
		// cout << type << endl;
		if(type=="zeroGradient"){
			setFunction = [id_inp, id_L_out, id_R_out](
			double* cells, double* faces) ->int {
				faces[id_L_out] = cells[id_inp];
				faces[id_R_out] = faces[id_L_out];
				// return inp[0];
				return 0;
			};
		}
		else if(type=="fixedValue"){
			double value = boundary.values[i];
			setFunction = [id_inp, id_L_out, id_R_out, value](
			double* cells, double* faces) ->int {
				faces[id_L_out] = value;
				faces[id_R_out] = value;
				return 0;
			};
		}
		else if(type=="slip"){
			int id_u = tmp_vel_cell_id[0]; 
			int id_v = tmp_vel_cell_id[1]; 
			int id_w = tmp_vel_cell_id[2]; 
			int id_u_L = tmp_vel_L_id[0]; 
			int id_v_L = tmp_vel_L_id[1]; 
			int id_w_L = tmp_vel_L_id[2];
			int id_u_R = tmp_vel_R_id[0]; 
			int id_v_R = tmp_vel_R_id[1]; 
			int id_w_R = tmp_vel_R_id[2];
			int id_nx = controls.faceVar["x unit normal"].id;
			int id_ny = controls.faceVar["y unit normal"].id;
			int id_nz = controls.faceVar["z unit normal"].id;
			setFunction = [id_u, id_v, id_w, 
			id_u_L, id_v_L, id_w_L, id_u_R, id_v_R, id_w_R,
			id_nx,id_ny,id_nz](
			double* cells, double* faces) ->int {
				double u_vel = cells[id_u];
				double v_vel = cells[id_v];
				double w_vel = cells[id_w];
				double norVel = u_vel * faces[id_nx] + 
								v_vel * faces[id_ny] + 
								w_vel * faces[id_nz];
				double invU = u_vel - norVel * faces[id_nx];
				double invV = v_vel - norVel * faces[id_ny];
				double invW = w_vel - norVel * faces[id_nz];
				
				faces[id_u_L] = invU; faces[id_u_R] = invU;
				faces[id_v_L] = invU; faces[id_v_R] = invV;
				faces[id_w_L] = invU; faces[id_w_R] = invW;
				return 0;
			};
			// ++i;
			// ++i;
		}
		else if(type=="noSlip"){
			setFunction = [id_inp, id_L_out, id_R_out](
			double* cells, double* faces) ->int {
				faces[id_L_out] = 0.0;
				faces[id_R_out] = 0.0;
				return 0;
			};
		}
		else if(type=="function"){
			string libraryName = "./constant/boundaryFunctions.so";
			void *handle = dlopen(libraryName.c_str(), RTLD_NOW);
			// typedef void (*setFunc_t)(double, double, double, double, double&);
			*(void **) (&setFunction) = dlsym(handle, "setFunction");
		}
		solver.calcBoundFacePrimVal.push_back(setFunction);
			
	}
	
	
	// cout << mesh.boundaries.size() << endl;
	// cout << controls.boundaryFunct.size() << endl;
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
}




void MASCH_Solver::setFunctions(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	/*
	cell : p u v w T Y1 Y2 rho c Ht drhodp drhodT drhodY1 drhodY2 dHtdp dHtdT dHtdY1 dH2dY2
	face : p u v w T Y1 Y2 rho c Ht
	field : dt
	
	*/
	
	solver.setOldVFunctionsUDF(mesh, controls);
	
	solver.setGradFunctionsUDF(mesh, controls);
	
	solver.setRecoFunctionsUDF(mesh, controls);
	solver.setFValFunctionsUDF(mesh, controls);
	
	solver.setTempFunctionsUDF(mesh, controls);
	solver.setConvFunctionsUDF(mesh, controls);
	solver.setDiffFunctionsUDF(mesh, controls);
	solver.setSourFunctionsUDF(mesh, controls);
	
	solver.setAddiFunctionsUDF(mesh, controls);
	
	
}







double MASCH_NVD::Minmod(double phiUU, double phiU, double phiD) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5){
		double tildeCf = 1.5*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<1.0){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	return ( gamma_f * phiD + (1.0-gamma_f) * phiU );

}






void MASCH_Solver::updateOldValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto cells = mesh.cells.data();
	auto cellVar = var.cells.data();
	
	auto p_saveCurrIds = saveCurrIds.data();
	auto p_saveOld1Ids = saveOld1Ids.data();
	auto p_saveOld2Ids = saveOld2Ids.data();
	
	// cout << cellVar[0].size() << endl;
	
	for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
		auto cellVar_i = cellVar[i].data();
		for(int j=0, jSIZE=saveCurrIds.size(); j<jSIZE; ++j){
			int curr_id = p_saveCurrIds[j];
			int old1_id = p_saveOld1Ids[j];
			int old2_id = p_saveOld2Ids[j];
			cellVar_i[old2_id] = cellVar_i[old1_id];
			cellVar_i[old1_id] = cellVar_i[curr_id];
		}
	}
	
	
}
void MASCH_Solver::updateBoundaryFaceValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	// auto& solver = (*this);
	
	// vector<int> prim_local_id;
	// vector<int> prim_face_left_id;
	// vector<int> prim_face_right_id;
	// for(auto [key, value] : controls.primitiveMap){
		// prim_local_id.push_back(value);
		// string left_s, right_s;
		// left_s = "left " + key;
		// right_s = "right " + key;
		// prim_face_left_id.push_back(controls.faceVar[left_s].id);
		// prim_face_right_id.push_back(controls.faceVar[right_s].id);
	// }
	// int tmp_prim_size = prim_local_id.size();
	
	// auto cells = mesh.cells.data();
	// auto faces = mesh.faces.data();
	// auto cellVar = var.cells.data();
	// auto faceVar = var.faces.data();
	
	// for(auto& boundary : mesh.boundaries){
		
		// if(boundary.getType() == MASCH_Face_Types::BOUNDARY){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = faces[i];
				// auto& cell = cells[face.iL];
				// auto faceVar_i = faceVar[i].data();
				// auto cellVar_i = cellVar[face.iL].data();
					
				// // primitive values
				// for(int j=0; j<tmp_prim_size; ++j){
					// int tmp_prim_local_id = prim_local_id[j];
					// int tmp_prim_face_left_id = prim_face_left_id[j];
					// int tmp_prim_face_right_id = prim_face_right_id[j];
					// faceVar_i[tmp_prim_face_left_id] = 
					// controls.boundaryFunct[tmp_prim_local_id](cellVar_i);
					// faceVar_i[tmp_prim_face_right_id] = 
					// faceVar_i[tmp_prim_face_left_id];
				// }
				// // additional values
				// solver.calcAddiFaceValues(faceVar_i);
			// }
		// }
	// }
	
	
	
}










void MASCH_Solver::eosIdeal(
	double cv, double gam,
	double& P, double& U, double& V, double& W, double& T,
	double& rho, double& C, double& Ht,
	double& drhodP, double& drhodT, double& dhdP, double& dhdT){
		
		
	double usqrt = U*U+V*V+W*W;
	
		
	// double cv = species.cv;
    // double gam = species.gamma;
	
	double cp = gam*cv;
		
	// density of each phase
	rho = 1.0/( (gam-1.0)*cv*T/(P) );
	C = sqrt( gam/(rho*rho)*(P)/(1.0/rho) );

	// d(rho)/d(p)
	drhodP = rho*rho*(1.0/rho)/(P); 
	
	// d(rho)/d(T)
	drhodT = -rho*rho*(1.0/rho)/T;

	// d(h)/d(p)
	dhdP = 0.0;
	// d(h)/d(T)
	dhdT = gam*cv;

	// internal energy of each phase
	// internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	// enthalpy of each phase
	double enthalpy = gam*cv*T;

		   
	// eti = internal_energy + 0.5*usqrt
	Ht = enthalpy + 0.5*usqrt;

	// cvi = cv
	// cpi = cp	


		
		
}


void MASCH_Solver::eosNASG(
	double pinf, double cv, double gam, double bNASG, double q,
	double& P, double& U, double& V, double& W, double& T,
	double& rho, double& C, double& Ht,
	double& drhodP, double& drhodT, double& dhdP, double& dhdT){
		
		
	double usqrt = U*U+V*V+W*W;
	
		
    // double pinf = species.Pinf; 
	// double cv = species.cv;
    // double gam = species.gamma;
    // double bNASG = species.b; 
	// double q = species.q;
	
	double cp = gam*cv;
		
	// density of each phase
	rho = 1.0/( (gam-1.0)*cv*T/(P+pinf)+bNASG );
	C = sqrt( gam/(rho*rho)*(P+pinf)/(1.0/rho-bNASG) );

	// d(rho)/d(p)
	drhodP = rho*rho*(1.0/rho-bNASG)/(P+pinf); 
	
	// d(rho)/d(T)
	drhodT = -rho*rho*(1.0/rho-bNASG)/T;

	// d(h)/d(p)
	dhdP = bNASG;
	// d(h)/d(T)
	dhdT = gam*cv;

	// internal energy of each phase
	// internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	// enthalpy of each phase
	double enthalpy = gam*cv*T + bNASG*P + q;

		   
	// eti = internal_energy + 0.5*usqrt
	Ht = enthalpy + 0.5*usqrt;

	// cvi = cv
	// cpi = cp	


		
		
}