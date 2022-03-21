
#include "../../../../others/solvers.h"


void MASCH_Solver::setTempStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	using tmp_dtstep_Cell_Funct_type = 
	function<int(double* faces, double* fields)>;
	using tmp_dtstep_Face_Funct_type = 
	function<int(double* cellsL, double* cellsR, double* faces, double* fields)>;
	auto& tempFunctCell = solver.calcTempStepCell;
	auto& tempFunctFace = solver.calcTempStepFace;
	// functTypeCell setFunctionCell;
	// functTypeFace setFunctionFace;
	
	{
		// using US = unsigned short;
		
		int id_dt = controls.getId_fieldVar("time-step");
		
		int id_vol = controls.getId_cellVar("volume");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_c = controls.getId_cellVar("speed-of-sound");
		
		double CFL = controls.maxCFL;
		
		solver.calcTempStepCell.push_back(
		[id_dt,id_vol,
		id_u,id_v,id_w,id_c,
		CFL](
		double* cells, double* fields) ->int {
			double mag_vel = sqrt(
				cells[id_u]*cells[id_u]+
				cells[id_v]*cells[id_v]+
				cells[id_w]*cells[id_w]);
			double c = cells[id_c];
			double dt_tmp = CFL * pow(cells[id_vol],0.3)/(mag_vel+c);
			fields[id_dt] = min(fields[id_dt],dt_tmp);
			
		});
		
		
		// for(int i=0, ip=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
			// auto& face = mesh.faces[i];
			// auto cellVar_iL = cellVar[face.iL].data();
			// if(face.getType()==MASCH_Face_Types::INTERNAL){
				// auto cellVar_iR = cellVar[face.iR].data();
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,cellVar_iR,faceVar[i].data(),fieldVar);
				// }
			// }
			// else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				// auto cellVar_iR = var.procRightCells[ip].data();
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,cellVar_iR,faceVar[i].data(),fieldVar);
				// }
				// ++ip;
			// }
			// else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				// for(auto& sol : solver.calcTempStepFace){
					// sol(cellVar_iL,nullptr,faceVar[i].data(),fieldVar);
				// }
			// }
		// }
		
		
		
		
		// int id_CFL = controls.getId_cellVar("Courant-Friedrichs-Lewy number");
		// double tmp_dt = coeffCellCFL*var.fields[id_dt];
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double mag_vel = sqrt(
				// var.cells[i][id_u]*var.cells[i][id_u]+
				// var.cells[i][id_v]*var.cells[i][id_v]+
				// var.cells[i][id_w]*var.cells[i][id_w]);
			// double c = var.cells[i][id_c];
			// double dt_cell = pow(var.cells[i][id_vol],0.3)/(mag_vel+c);
			// var.cells[i][id_CFL] = tmp_dt/dt_cell;
		// }
		
		
		// tempFunctCell.push_back(
		// [id_dt,id_area,id_volume,
		// id_u,id_v,id_w,id_c,
		// id_u_L,id_v_L,id_w_L,id_c_L,
		// id_u_R,id_v_R,id_w_R,id_c_R](
		// double* cellsL, double* cellsR, 
		// double* faces, double* fields) ->int {
			// double mag_vel = sqrt(
				// var.cells[i][id_u]*var.cells[i][id_u]+
				// var.cells[i][id_v]*var.cells[i][id_v]+
				// var.cells[i][id_w]*var.cells[i][id_w]);
			// double c = var.cells[i][id_c];
			// double dt_tmp = pow(var.cells[i][id_vol],0.3)/(mag_vel+c);
			// var.fields[id_dt] = min(var.fields[id_dt],dt_tmp);
			
			// double uF = 0.0;//0.5*(cellsL[id_u]+cellsR[id_u]);
			// double vF = 0.0;//0.5*(cellsL[id_v]+cellsR[id_v]);
			// double wF = 0.0;//0.5*(cellsL[id_w]+cellsR[id_w]);
			// double cF = 0.0;//0.5*(cellsL[id_c]+cellsR[id_c]);
			// double volume = 0.0;
			// double magVel = 0.0;
			// if(cellsR!=nullptr){
				// // uF = 0.5*(cellsL[id_u]+cellsR[id_u]);
				// // vF = 0.5*(cellsL[id_v]+cellsR[id_v]);
				// // wF = 0.5*(cellsL[id_w]+cellsR[id_w]);
				// // cF = 0.5*(cellsL[id_c]+cellsR[id_c]);
				// // volume = 0.5*(cellsL[id_volume]+cellsR[id_volume]);
				// double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					// cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				// double magVel_R = sqrt(cellsR[id_u]*cellsR[id_u]+
					// cellsR[id_v]*cellsR[id_v]+cellsR[id_w]*cellsR[id_w]);
				// magVel = max(magVel_L,magVel_R);
				// cF = max(cellsL[id_c],cellsR[id_c]);
				// // cF = cellsL[id_c];
				// volume = min(cellsL[id_volume],cellsR[id_volume]);
				// // volume = cellsL[id_volume];
				// // if(volume<1.e-30) cout << volume << endl;
			// }
			// else{
				// double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					// cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				// magVel=magVel_L;
				// cF = cellsL[id_c];
				// volume = cellsL[id_volume];
			// }
			
			// fields[id_dt] = min(fields[id_dt], volume/( (magVel+cF)*faces[id_area] ));
			
		// }); 
		
		
		
		// tempFunctFace.push_back(
		// [id_dt,id_area,id_volume,
		// id_u,id_v,id_w,id_c,
		// id_u_L,id_v_L,id_w_L,id_c_L,
		// id_u_R,id_v_R,id_w_R,id_c_R](
		// double* cellsL, double* cellsR, 
		// double* faces, double* fields) ->int {
			// double uF = 0.0;//0.5*(cellsL[id_u]+cellsR[id_u]);
			// double vF = 0.0;//0.5*(cellsL[id_v]+cellsR[id_v]);
			// double wF = 0.0;//0.5*(cellsL[id_w]+cellsR[id_w]);
			// double cF = 0.0;//0.5*(cellsL[id_c]+cellsR[id_c]);
			// double volume = 0.0;
			// double magVel = 0.0;
			// if(cellsR!=nullptr){
				// // uF = 0.5*(cellsL[id_u]+cellsR[id_u]);
				// // vF = 0.5*(cellsL[id_v]+cellsR[id_v]);
				// // wF = 0.5*(cellsL[id_w]+cellsR[id_w]);
				// // cF = 0.5*(cellsL[id_c]+cellsR[id_c]);
				// // volume = 0.5*(cellsL[id_volume]+cellsR[id_volume]);
				// double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					// cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				// double magVel_R = sqrt(cellsR[id_u]*cellsR[id_u]+
					// cellsR[id_v]*cellsR[id_v]+cellsR[id_w]*cellsR[id_w]);
				// magVel = max(magVel_L,magVel_R);
				// cF = max(cellsL[id_c],cellsR[id_c]);
				// // cF = cellsL[id_c];
				// volume = min(cellsL[id_volume],cellsR[id_volume]);
				// // volume = cellsL[id_volume];
				// // if(volume<1.e-30) cout << volume << endl;
			// }
			// else{
				// double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					// cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				// magVel=magVel_L;
				// cF = cellsL[id_c];
				// volume = cellsL[id_volume];
			// }
			
			// fields[id_dt] = min(fields[id_dt], volume/( (magVel+cF)*faces[id_area] ));
			
		// }); 
		
	}
	
	
}
