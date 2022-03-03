
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
		using US = unsigned short;
		
		// 필드 값
		int id_dt = controls.fieldVar["time-step"].id;
		
		// 페이스 값
		// int id_F_CFL = controls.faceVar["Courant-Friedrichs-Lewy number"].id;
		
		int id_u_L = controls.faceVar["left x-velocity"].id;
		int id_v_L = controls.faceVar["left y-velocity"].id;
		int id_w_L = controls.faceVar["left z-velocity"].id;
		int id_c_L = controls.faceVar["left speed of sound"].id;
		
		int id_u_R = controls.faceVar["right x-velocity"].id;
		int id_v_R = controls.faceVar["right y-velocity"].id;
		int id_w_R = controls.faceVar["right z-velocity"].id;
		int id_c_R = controls.faceVar["right speed of sound"].id;
		
		// 셀 값
		// int id_C_CFL = controls.cellVar["Courant-Friedrichs-Lewy number"].id;
		int id_u = controls.cellVar["x-velocity"].id;
		int id_v = controls.cellVar["y-velocity"].id;
		int id_w = controls.cellVar["z-velocity"].id;
		int id_c = controls.cellVar["speed of sound"].id;
		
		// cout << id_u << " " << id_v << " " << id_w << " " << id_c << endl;
		
		// 메쉬관련
		int id_area = controls.faceVar["area"].id;
		int id_volume = controls.cellVar["volume"].id;
		
		tempFunctFace.push_back(
		[id_dt,id_area,id_volume,
		id_u,id_v,id_w,id_c,
		id_u_L,id_v_L,id_w_L,id_c_L,
		id_u_R,id_v_R,id_w_R,id_c_R](
		double* cellsL, double* cellsR, 
		double* faces, double* fields) ->int {
			double uF = 0.0;//0.5*(cellsL[id_u]+cellsR[id_u]);
			double vF = 0.0;//0.5*(cellsL[id_v]+cellsR[id_v]);
			double wF = 0.0;//0.5*(cellsL[id_w]+cellsR[id_w]);
			double cF = 0.0;//0.5*(cellsL[id_c]+cellsR[id_c]);
			double volume = 0.0;
			double magVel = 0.0;
			if(cellsR!=nullptr){
				// uF = 0.5*(cellsL[id_u]+cellsR[id_u]);
				// vF = 0.5*(cellsL[id_v]+cellsR[id_v]);
				// wF = 0.5*(cellsL[id_w]+cellsR[id_w]);
				// cF = 0.5*(cellsL[id_c]+cellsR[id_c]);
				// volume = 0.5*(cellsL[id_volume]+cellsR[id_volume]);
				double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				double magVel_R = sqrt(cellsR[id_u]*cellsR[id_u]+
					cellsR[id_v]*cellsR[id_v]+cellsR[id_w]*cellsR[id_w]);
				magVel = max(magVel_L,magVel_R);
				cF = max(cellsL[id_c],cellsR[id_c]);
				// cF = cellsL[id_c];
				volume = min(cellsL[id_volume],cellsR[id_volume]);
				// volume = cellsL[id_volume];
				// if(volume<1.e-30) cout << volume << endl;
			}
			else{
				double magVel_L = sqrt(cellsL[id_u]*cellsL[id_u]+
					cellsL[id_v]*cellsL[id_v]+cellsL[id_w]*cellsL[id_w]);
				magVel=magVel_L;
				cF = cellsL[id_c];
				volume = cellsL[id_volume];
			}
			
			fields[id_dt] = min(fields[id_dt], volume/( (magVel+cF)*faces[id_area] ));
			
		}); 
		
	}
	
	
}
