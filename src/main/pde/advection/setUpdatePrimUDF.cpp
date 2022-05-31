
#include "../../../others/solvers.h"


void MASCH_Solver::setUpdatePrimFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	vector<string> s_iter = load.extractVector(controls.fvSolutionMap["coupled.relaxationFactors"]);
	{
		int id_p = controls.getId_cellVar("pressure");
		vector<string> p_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.p"]);
		double p_min = stod(p_lim[0]); double p_max = stod(p_lim[1]);
		int id_u = controls.getId_cellVar("x-velocity");
		vector<string> u_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.u"]);
		double u_min = stod(u_lim[0]); double u_max = stod(u_lim[1]);
		int id_v = controls.getId_cellVar("y-velocity");
		vector<string> v_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.v"]);
		double v_min = stod(v_lim[0]); double v_max = stod(v_lim[1]);
		int id_w = controls.getId_cellVar("z-velocity");
		vector<string> w_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.w"]);
		double w_min = stod(w_lim[0]); double w_max = stod(w_lim[1]);
		int id_T = controls.getId_cellVar("temperature");
		vector<string> T_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.T"]);
		double T_min = stod(T_lim[0]); double T_max = stod(T_lim[1]);
		vector<int> id_Y;
		vector<string> Y_lim = load.extractVector(controls.fvSolutionMap["limiterCellValues.Y"]);
		double Y_min = stod(Y_lim[0]); double Y_max = stod(Y_lim[1]);
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("mass-fraction-"+controls.spName[i]);
			id_Y.push_back(controls.getId_cellVar(tmp_name));
		}
		int nSpm1 = controls.spName.size()-1;
		
		vector<double> relaxation_factor;
        for(auto& item : s_iter){
            relaxation_factor.push_back(stod(item));
        }
        
        
        
        
		int id_p_old = controls.getId_cellVar("old pressure");
		int id_u_old = controls.getId_cellVar("old x-velocity");
		int id_v_old = controls.getId_cellVar("old y-velocity");
		int id_w_old = controls.getId_cellVar("old z-velocity");
		int id_T_old = controls.getId_cellVar("old temperature");
		vector<int> id_Y_old;
		for(int i=0; i<controls.spName.size()-1; ++i){
			string tmp_name = ("old mass-fraction-"+controls.spName[i]);
			id_Y_old.push_back(controls.getId_cellVar(tmp_name));
		}
        
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		
		calcUpdatePrim.push_back(
		[&controls,id_dt,id_vol,id_p,id_u,id_v,id_w,id_T,id_Y,relaxation_factor,nSpm1,
		u_min,u_max,v_min,v_max,w_min,w_max,T_min,T_max,p_min,p_max,Y_min,Y_max,
        id_p_old,id_u_old,id_v_old,id_w_old,id_T_old,id_Y_old](
		double* fields, double* cells, double* Xvalues) ->int {
            int ii = controls.iterPseudo;
            double relax = relaxation_factor[ii];
			// cells[id_p] += Xvalues[0];
			// cells[id_u] += Xvalues[1];
			// cells[id_v] += Xvalues[2];
			// cells[id_w] += Xvalues[3];
			// cells[id_T] += Xvalues[4];
            
            // cells[id_p] = (1.0-relax)*cells[id_p_old] + relax*cells[id_p];
            // cells[id_u] = (1.0-relax)*cells[id_u_old] + relax*cells[id_u];
            // cells[id_v] = (1.0-relax)*cells[id_v_old] + relax*cells[id_v];
            // cells[id_w] = (1.0-relax)*cells[id_w_old] + relax*cells[id_w];
            // cells[id_T] = (1.0-relax)*cells[id_T_old] + relax*cells[id_T];
			
			// cells[id_p] = max(p_min,min(p_max,cells[id_p]));
			// cells[id_u] = max(u_min,min(u_max,cells[id_u]));
			// cells[id_v] = max(v_min,min(v_max,cells[id_v]));
			// cells[id_w] = max(w_min,min(w_max,cells[id_w]));
			// cells[id_T] = max(T_min,min(T_max,cells[id_T]));
            
			double vol = cells[id_vol];
			double dt = fields[id_dt];
			
			for(int isp=0; isp<nSpm1; ++isp){
				cells[id_Y[isp]] += Xvalues[5+isp];
                
                
                // cells[id_Y[isp]] = (1.0-relax)*cells[id_Y_old[isp]] + relax*cells[id_Y[isp]];
                
                
				cells[id_Y[isp]] = max(Y_min,min(Y_max,cells[id_Y[isp]]));
			}
            
            // if((cells[id_p])>1.e12 || (cells[id_p])<-1.e12) cout << "p " << cells[id_p] << Xvalues[0] << endl;
            // if((cells[id_u])) cout << "u " << cells[id_u] << endl;
            // if((cells[id_v])) cout << "v " << cells[id_v] << endl;
            // if((cells[id_w])) cout << "w " << cells[id_w] << endl;
            // if((cells[id_T])) cout << "T " << cells[id_T] << endl;
			
		}); 
		
	}
	
}
