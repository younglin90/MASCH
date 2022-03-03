
#include "./solvers.h"

void MASCH_Solver::fvm_inline(MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	
	auto& solver = (*this);
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	double CFL = 0.1;
	{
		controls.log.push("calc1");
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		int id_c = controls.getId_cellVar("speed of sound");
		int id_rho = controls.getId_cellVar("density");
		int id_Ht = controls.getId_cellVar("total enthalpy");
		int id_drhodp = controls.getId_cellVar("density diff with pressure");
		int id_drhodT = controls.getId_cellVar("density diff with temperature");
		int id_dHtdp = controls.getId_cellVar("total enthalpy diff with pressure");
		int id_dHtdT = controls.getId_cellVar("total enthalpy diff with temperature");
		controls.log.pop();
		controls.log.push("calc2");
		solver.updateProcRightCellPrimValues(mesh, controls, var);
		// solver.updateProcRightCellAddiValues(mesh, controls, var);
		
		controls.log.pop();
		controls.log.push("calc3");
		for(int i=0; i<mesh.cells.size(); ++i){
			var.cells[i][id_rho] = var.cells[i][id_p]/287.0/var.cells[i][id_T];
			var.cells[i][id_c] = sqrt(1.4*287.0*var.cells[i][id_T]);
			var.cells[i][id_Ht] = 1.4*717.0*var.cells[i][id_T] + 0.5*(
				var.cells[i][id_u]*var.cells[i][id_u]+
				var.cells[i][id_v]*var.cells[i][id_v]+
				var.cells[i][id_w]*var.cells[i][id_w]);
			
			var.cells[i][id_drhodp] = var.cells[i][id_rho]/var.cells[i][id_p];
			var.cells[i][id_drhodT] = -var.cells[i][id_rho]/var.cells[i][id_T];
			var.cells[i][id_dHtdp] = 0.0;
			var.cells[i][id_dHtdT] = 1.4*717.0;
		}
		controls.log.pop();
		controls.log.push("calc4");
		for(auto& cells_i : var.procRightCells){
			cells_i[id_rho] = cells_i[id_p]/287.0/cells_i[id_T];
			cells_i[id_c] = sqrt(1.4*287.0*cells_i[id_T]);
			cells_i[id_Ht] = 1.4*717.0*cells_i[id_T] + 0.5*(
				cells_i[id_u]*cells_i[id_u]+
				cells_i[id_v]*cells_i[id_v]+
				cells_i[id_w]*cells_i[id_w]);
			
			cells_i[id_drhodp] = cells_i[id_rho]/cells_i[id_p];
			cells_i[id_drhodT] = -cells_i[id_rho]/cells_i[id_T];
			cells_i[id_dHtdp] = 0.0;
			cells_i[id_dHtdT] = 1.4*717.0;
		}
		controls.log.pop();
		controls.log.push("calc5");
		var.fields[id_dt] = 1.e12;
		for(int i=0; i<mesh.cells.size(); ++i){
			double mag_vel = sqrt(
				var.cells[i][id_u]*var.cells[i][id_u]+
				var.cells[i][id_v]*var.cells[i][id_v]+
				var.cells[i][id_w]*var.cells[i][id_w]);
			double c = var.cells[i][id_c];
			double dt_tmp = pow(var.cells[i][id_vol],0.3)/(mag_vel+c);
			var.fields[id_dt] = min(var.fields[id_dt],dt_tmp);
		}
		controls.log.pop();
		controls.log.push("calc6");
		if(size>1){
			double tmp_fieldVar = var.fields[id_dt];
			double tmp_fieldVar_glo;
			MPI_Allreduce(&tmp_fieldVar, &tmp_fieldVar_glo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			var.fields[id_dt] = tmp_fieldVar_glo;
		}
		var.fields[id_dt] *= CFL;
		
	}
	{
		controls.log.pop();
		controls.log.push("calc8");
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		int id_c = controls.getId_cellVar("speed of sound");
		int id_rho = controls.getId_cellVar("density");
		int id_Ht = controls.getId_cellVar("total enthalpy");
		
		int id_pL = controls.getId_faceVar("left pressure");
		int id_uL = controls.getId_faceVar("left x-velocity");
		int id_vL = controls.getId_faceVar("left y-velocity");
		int id_wL = controls.getId_faceVar("left z-velocity");
		int id_TL = controls.getId_faceVar("left temperature");
		int id_cL = controls.getId_faceVar("left speed of sound");
		int id_rhoL = controls.getId_faceVar("left density");
		int id_HtL = controls.getId_faceVar("left total enthalpy");
		
		int id_pR = controls.getId_faceVar("right pressure");
		int id_uR = controls.getId_faceVar("right x-velocity");
		int id_vR = controls.getId_faceVar("right y-velocity");
		int id_wR = controls.getId_faceVar("right z-velocity");
		int id_TR = controls.getId_faceVar("right temperature");
		int id_cR = controls.getId_faceVar("right speed of sound");
		int id_rhoR = controls.getId_faceVar("right density");
		int id_HtR = controls.getId_faceVar("right total enthalpy");
		
		int id_area = controls.getId_faceVar("area");
		int id_nx = controls.getId_faceVar("x unit normal");
		int id_ny = controls.getId_faceVar("y unit normal");
		int id_nz = controls.getId_faceVar("z unit normal");
		
		int nEq = controls.nEq;
		
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto procRightCellVar = var.procRightCells.data();
		
		controls.log.pop();
		controls.log.push("calc9");
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			if(boundary.name=="invwall"){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					auto faceVar_i = faceVar[i].data();
					auto cellVar_iL = cellVar[iL].data();
					
					faceVar_i[id_pL] = cellVar_iL[id_p];
					faceVar_i[id_TL] = cellVar_iL[id_T];
					
					double u_vel = cellVar_iL[id_u];
					double v_vel = cellVar_iL[id_v];
					double w_vel = cellVar_iL[id_w];
					double norVel = u_vel * faceVar_i[id_nx] + 
									v_vel * faceVar_i[id_ny] + 
									w_vel * faceVar_i[id_nz];
					double invU = u_vel - norVel * faceVar_i[id_nx];
					double invV = v_vel - norVel * faceVar_i[id_ny];
					double invW = w_vel - norVel * faceVar_i[id_nz];
					
					faceVar_i[id_uL] = invU;
					faceVar_i[id_vL] = invV;
					faceVar_i[id_wL] = invW;
					
					faceVar_i[id_rhoL] = faceVar_i[id_pL]/287.0/faceVar_i[id_TL];
					faceVar_i[id_cL] = sqrt(1.4*287.0*faceVar_i[id_TL]);
					faceVar_i[id_HtL] = 1.4*717.0*faceVar_i[id_TL]+0.5*(
						faceVar_i[id_uL]*faceVar_i[id_uL]+
						faceVar_i[id_vL]*faceVar_i[id_vL]+
						faceVar_i[id_wL]*faceVar_i[id_wL]);
						
					faceVar_i[id_pR]=faceVar_i[id_pL];
					faceVar_i[id_uR]=faceVar_i[id_uL];
					faceVar_i[id_vR]=faceVar_i[id_vL];
					faceVar_i[id_wR]=faceVar_i[id_wL];
					faceVar_i[id_TR]=faceVar_i[id_TL];
					faceVar_i[id_rhoR]=faceVar_i[id_rhoL];
					faceVar_i[id_cR]=faceVar_i[id_cL];
					faceVar_i[id_HtR]=faceVar_i[id_HtL];
				}
			}
			else if(boundary.name=="viswall"){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					auto faceVar_i = faceVar[i].data();
					auto cellVar_iL = cellVar[iL].data();
					faceVar_i[id_pL] = cellVar_iL[id_p];
					faceVar_i[id_TL] = cellVar_iL[id_T];
					faceVar_i[id_uL] = 0.0;
					faceVar_i[id_vL] = 0.0;
					faceVar_i[id_wL] = 0.0;
					
					faceVar_i[id_rhoL] = faceVar_i[id_pL]/287.0/faceVar_i[id_TL];
					faceVar_i[id_cL] = sqrt(1.4*287.0*faceVar_i[id_TL]);
					faceVar_i[id_HtL] = 1.4*717.0*faceVar_i[id_TL]+0.5*(
						faceVar_i[id_uL]*faceVar_i[id_uL]+
						faceVar_i[id_vL]*faceVar_i[id_vL]+
						faceVar_i[id_wL]*faceVar_i[id_wL]);
						
					faceVar_i[id_pR]=faceVar_i[id_pL];
					faceVar_i[id_uR]=faceVar_i[id_uL];
					faceVar_i[id_vR]=faceVar_i[id_vL];
					faceVar_i[id_wR]=faceVar_i[id_wL];
					faceVar_i[id_TR]=faceVar_i[id_TL];
					faceVar_i[id_rhoR]=faceVar_i[id_rhoL];
					faceVar_i[id_cR]=faceVar_i[id_cL];
					faceVar_i[id_HtR]=faceVar_i[id_HtL];
				}
			}
			else if(boundary.name=="supinlet"){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					auto faceVar_i = faceVar[i].data();
					auto cellVar_iL = cellVar[iL].data();
					faceVar_i[id_pL] = 1.0;
					faceVar_i[id_TL] = 0.0035;
					faceVar_i[id_uL] = 3.6;
					faceVar_i[id_vL] = 0.0;
					faceVar_i[id_wL] = 0.0;
					
					faceVar_i[id_rhoL] = faceVar_i[id_pL]/287.0/faceVar_i[id_TL];
					faceVar_i[id_cL] = sqrt(1.4*287.0*faceVar_i[id_TL]);
					faceVar_i[id_HtL] = 1.4*717.0*faceVar_i[id_TL]+0.5*(
						faceVar_i[id_uL]*faceVar_i[id_uL]+
						faceVar_i[id_vL]*faceVar_i[id_vL]+
						faceVar_i[id_wL]*faceVar_i[id_wL]);
						
					faceVar_i[id_pR]=faceVar_i[id_pL];
					faceVar_i[id_uR]=faceVar_i[id_uL];
					faceVar_i[id_vR]=faceVar_i[id_vL];
					faceVar_i[id_wR]=faceVar_i[id_wL];
					faceVar_i[id_TR]=faceVar_i[id_TL];
					faceVar_i[id_rhoR]=faceVar_i[id_rhoL];
					faceVar_i[id_cR]=faceVar_i[id_cL];
					faceVar_i[id_HtR]=faceVar_i[id_HtL];
				}
			}
			else if(boundary.name=="outlet"){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					int iL = face.iL;
					auto faceVar_i = faceVar[i].data();
					auto cellVar_iL = cellVar[iL].data();
					faceVar_i[id_pL] = cellVar_iL[id_p];
					faceVar_i[id_TL] = cellVar_iL[id_T];
					faceVar_i[id_uL] = cellVar_iL[id_u];
					faceVar_i[id_vL] = cellVar_iL[id_v];
					faceVar_i[id_wL] = cellVar_iL[id_w];
					
					faceVar_i[id_rhoL] = faceVar_i[id_pL]/287.0/faceVar_i[id_TL];
					faceVar_i[id_cL] = sqrt(1.4*287.0*faceVar_i[id_TL]);
					faceVar_i[id_HtL] = 1.4*717.0*faceVar_i[id_TL]+0.5*(
						faceVar_i[id_uL]*faceVar_i[id_uL]+
						faceVar_i[id_vL]*faceVar_i[id_vL]+
						faceVar_i[id_wL]*faceVar_i[id_wL]);
						
					faceVar_i[id_pR]=faceVar_i[id_pL];
					faceVar_i[id_uR]=faceVar_i[id_uL];
					faceVar_i[id_vR]=faceVar_i[id_vL];
					faceVar_i[id_wR]=faceVar_i[id_wL];
					faceVar_i[id_TR]=faceVar_i[id_TL];
					faceVar_i[id_rhoR]=faceVar_i[id_rhoL];
					faceVar_i[id_cR]=faceVar_i[id_cL];
					faceVar_i[id_HtR]=faceVar_i[id_HtL];
				}
			}
			else{
				cout << "#WARNING, NO boundary name" << endl;
			}
		}
		
		controls.log.pop();
		controls.log.push("calc10");
		
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<nEq; ++j){
				var.tmp_B[i*nEq+j]=0.0;
			}
		}
		controls.log.pop();
		controls.log.push("calc11");
		
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			int iL = face.iL; int iR = face.iR;
			auto faceVar_i = faceVar[i].data();
			auto cellVar_iL = cellVar[iL].data();
		
			double pL = cellVar_iL[id_p]; double pR = 0.0;
			double uL = cellVar_iL[id_u]; double uR = 0.0;
			double vL = cellVar_iL[id_v]; double vR = 0.0;
			double wL = cellVar_iL[id_w]; double wR = 0.0;
			double TL = cellVar_iL[id_T]; double TR = 0.0;
			double rhoL = cellVar_iL[id_rho]; double rhoR = 0.0;
			double cL = cellVar_iL[id_c]; double cR = 0.0;
			double HtL = cellVar_iL[id_Ht]; double HtR = 0.0;
			
			if(face.getType()==MASCH_Face_Types::INTERNAL){
				pR = cellVar[iR][id_p];
				uR = cellVar[iR][id_u];
				vR = cellVar[iR][id_v];
				wR = cellVar[iR][id_w];
				TR = cellVar[iR][id_T];
				rhoR = cellVar[iR][id_rho];
				cR = cellVar[iR][id_c];
				HtR = cellVar[iR][id_Ht];
			}
			else if(face.getType()==MASCH_Face_Types::BOUNDARY){
				pL = faceVar_i[id_pL];
				uL = faceVar_i[id_uL];
				vL = faceVar_i[id_vL];
				wL = faceVar_i[id_wL];
				TL = faceVar_i[id_TL];
				rhoL = faceVar_i[id_rhoL];
				cL = faceVar_i[id_cL];
				HtL = faceVar_i[id_HtL];
				
				pR = faceVar_i[id_pR];
				uR = faceVar_i[id_uR];
				vR = faceVar_i[id_vR];
				wR = faceVar_i[id_wR];
				TR = faceVar_i[id_TR];
				rhoR = faceVar_i[id_rhoR];
				cR = faceVar_i[id_cR];
				HtR = faceVar_i[id_HtR];
			}
			else if(face.getType()==MASCH_Face_Types::PROCESSOR){
				pR = procRightCellVar[ip][id_p];
				uR = procRightCellVar[ip][id_u];
				vR = procRightCellVar[ip][id_v];
				wR = procRightCellVar[ip][id_w];
				TR = procRightCellVar[ip][id_T];
				rhoR = procRightCellVar[ip][id_rho];
				cR = procRightCellVar[ip][id_c];
				HtR = procRightCellVar[ip][id_Ht];
				
				++ip;
			}
			
			{
				double nvec[3];
				nvec[0] = faceVar_i[id_nx];
				nvec[1] = faceVar_i[id_ny];
				nvec[2] = faceVar_i[id_nz];
				double UnL = uL*nvec[0] + vL*nvec[1] + wL*nvec[2];
				double UnR = uR*nvec[0] + vR*nvec[1] + wR*nvec[2];
				
				double Unhat = 0.5*(UnL+UnR);
				double rhohat = 0.5*(rhoL+rhoR);
				double chat= 0.5*(cL+cR);
				
				double ML = UnL/chat; 
				double MR = UnR/chat;
				double KLR = sqrt(0.5*(uL*uL+vL*vL+wL*wL+uR*uR+vR*vR+wR*wR));
				double MLP = 0.5*(ML+abs(ML));
				if( abs(ML) < 1.0 ) {
					MLP = 0.25*(ML + 1.0)*(ML + 1.0);
				}
				double MRM = 0.5*(MR-abs(MR));
				if( abs(MR) < 1.0 ) {
					MRM = -0.25*(MR - 1.0)*(MR - 1.0);
				}
				double preP = 0.5*(1.0 + ( ML>0.0 ? 1.0 : -1.0 ) );
				if( abs(ML) < 1.0 ) {
					preP = 0.25*(ML+1.0)*(ML+1.0)*(2.0-ML);
				} 
				double preM = 0.5*(1.0 - ( MR>0.0 ? 1.0 : -1.0 ) );
				if( abs(MR) < 1.0 ) {
					preM = 0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);
				} 
				
				double Mbar = ( rhoL*abs(ML)+rhoR*abs(MR) ) / ( rhoL + rhoR );

				// SLAU
				double Mcy = min(1.0,KLR/chat);
				double phi_c = (1.0-Mcy)*(1.0-Mcy);
				double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );
				double D_L = ML+(1.0-g_c)*abs(ML);
				double D_R = MR-(1.0-g_c)*abs(MR);
				double D_rho = Mbar*g_c;
				double MLPL = 0.5*(D_L+D_rho);
				double MRMR = 0.5*(D_R-D_rho);

				double mdot = rhoL*chat*MLPL + rhoR*chat*MRMR - 0.5*phi_c/chat*(pR-pL);

				double f1L = mdot;
				double f1R = 0.0;
				if( mdot < 0.0 ) {
					f1L = 0.0; f1R = mdot;
				}

				double PLR = 0.5*(pL+pR) - 
							 0.5*Mcy*preP*preM*0.5*(pL+pR)/chat*(UnR-UnL) + 
							 Mcy*0.5*(pL+pR)*(preP+preM-1.0) - 
							 0.5*(preP-preM)*(pR-pL);
				double area = var.faces[i][id_area];
				
				double flux[nEq];
				flux[0] = (f1L+f1R)*area;
				flux[1] = (f1L*uL + f1R*uR + PLR*nvec[0])*area;
				flux[2] = (f1L*vL + f1R*vR + PLR*nvec[1])*area;
				flux[3] = (f1L*wL + f1R*wR + PLR*nvec[2])*area;
				flux[4] = (f1L*HtL+ f1R*HtR)*area;
			
				for(int j=0; j<nEq; ++j){
					var.tmp_B[iL*nEq+j] -= (flux[j]);
				}
				if(face.getType()==MASCH_Face_Types::INTERNAL){
					for(int j=0; j<nEq; ++j){
						var.tmp_B[iR*nEq+j] += (flux[j]);
					}
				}
			}
		}
		
	}
	
	
	{
		controls.log.pop();
		controls.log.push("calc12");
		
		int id_dt = controls.getId_fieldVar("time-step");
		int id_vol = controls.getId_cellVar("volume");
		int id_p = controls.getId_cellVar("pressure");
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_T = controls.getId_cellVar("temperature");
		int id_c = controls.getId_cellVar("speed of sound");
		int id_rho = controls.getId_cellVar("density");
		int id_Ht = controls.getId_cellVar("total enthalpy");
		int id_drhodp = controls.getId_cellVar("density diff with pressure");
		int id_drhodT = controls.getId_cellVar("density diff with temperature");
		int id_dHtdp = controls.getId_cellVar("total enthalpy diff with pressure");
		int id_dHtdT = controls.getId_cellVar("total enthalpy diff with temperature");
		
		int nEq = controls.nEq;
		
		double dt = var.fields[id_dt];
		int id_resi = controls.fieldVar["residual"].id;
		
		double tmp_volume = 0.0;
		
		auto cellVar = var.cells.data();
		auto faceVar = var.faces.data();
		auto fieldVar = var.fields.data();
		auto procRightCellVar = var.procRightCells.data();
		
		controls.log.pop();
		controls.log.push("calc13");
		for(int i=0; i<mesh.cells.size(); ++i){
			
			double volume = var.cells[i][id_vol];
			auto cellVar_i = cellVar[i].data();
			
			vector<vector<double>> matA(nEq,vector<double>(nEq,0.0));
			
			matA[0][0] = cellVar_i[id_drhodp] * volume / dt;
			matA[0][4] = cellVar_i[id_drhodT] * volume / dt;
			
			matA[1][0] = cellVar_i[id_drhodp]*cellVar_i[id_u] * volume / dt;
			matA[1][1] = cellVar_i[id_rho] * volume / dt;
			matA[1][4] = cellVar_i[id_drhodT]*cellVar_i[id_u] * volume / dt;
			
			matA[2][0] = cellVar_i[id_drhodp]*cellVar_i[id_v] * volume / dt;
			matA[2][2] = cellVar_i[id_rho] * volume / dt;
			matA[2][4] = cellVar_i[id_drhodT]*cellVar_i[id_v] * volume / dt;
			
			matA[3][0] = cellVar_i[id_drhodp]*cellVar_i[id_w] * volume / dt;
			matA[3][3] = cellVar_i[id_rho] * volume / dt;
			matA[3][4] = cellVar_i[id_drhodT]*cellVar_i[id_w] * volume / dt;
			
			matA[4][0] = (cellVar_i[id_drhodp]*cellVar_i[id_Ht]+
							cellVar_i[id_rho]*cellVar_i[id_dHtdp]-1.0) * volume / dt;
			matA[4][1] = cellVar_i[id_rho]*cellVar_i[id_u] * volume / dt;
			matA[4][2] = cellVar_i[id_rho]*cellVar_i[id_v] * volume / dt;
			matA[4][3] = cellVar_i[id_rho]*cellVar_i[id_w] * volume / dt;
			matA[4][4] = (cellVar_i[id_drhodT]*cellVar_i[id_Ht]+
							cellVar_i[id_rho]*cellVar_i[id_dHtdT]) * volume / dt;
				
			// math.GaussSeidel(matA.data(), nEq);
			math.GaussSeidelSOR(matA);
			// math.GaussSeidel(A.data(), B.data(), X.data(), nEq);
			
			
			for(int jEq=0; jEq<nEq; ++jEq){
				double tmp_val = 0.0;
				for(int iEq=0; iEq<nEq; ++iEq){
					// tmp_val += A[nEq*jEq+iEq]*var.tmp_B[nEq*i+iEq];
					tmp_val += matA[jEq][iEq]*var.tmp_B[nEq*i+iEq];
				}
				var.tmp_X[nEq*i+jEq] = tmp_val;
				
			}
			
			// vector<double> matA(nEq*nEq);
			
			// matA[0*nEq+0] = cellVar_i[id_drhodp] * volume / dt;
			// matA[0*nEq+4] = cellVar_i[id_drhodT] * volume / dt;
			
			// matA[1*nEq+0] = cellVar_i[id_drhodp]*cellVar_i[id_u] * volume / dt;
			// matA[1*nEq+1] = cellVar_i[id_rho] * volume / dt;
			// matA[1*nEq+4] = cellVar_i[id_drhodT]*cellVar_i[id_u] * volume / dt;
			
			// matA[2*nEq+0] = cellVar_i[id_drhodp]*cellVar_i[id_v] * volume / dt;
			// matA[2*nEq+2] = cellVar_i[id_rho] * volume / dt;
			// matA[2*nEq+4] = cellVar_i[id_drhodT]*cellVar_i[id_v] * volume / dt;
			
			// matA[3*nEq+0] = cellVar_i[id_drhodp]*cellVar_i[id_w] * volume / dt;
			// matA[3*nEq+3] = cellVar_i[id_rho] * volume / dt;
			// matA[3*nEq+4] = cellVar_i[id_drhodT]*cellVar_i[id_w] * volume / dt;
			
			// matA[4*nEq+0] = (cellVar_i[id_drhodp]*cellVar_i[id_Ht]+
							// cellVar_i[id_rho]*cellVar_i[id_dHtdp]-1.0) * volume / dt;
			// matA[4*nEq+1] = cellVar_i[id_rho]*cellVar_i[id_u] * volume / dt;
			// matA[4*nEq+2] = cellVar_i[id_rho]*cellVar_i[id_v] * volume / dt;
			// matA[4*nEq+3] = cellVar_i[id_rho]*cellVar_i[id_w] * volume / dt;
			// matA[4*nEq+4] = (cellVar_i[id_drhodT]*cellVar_i[id_Ht]+
							// cellVar_i[id_rho]*cellVar_i[id_dHtdT]) * volume / dt;
				
			// math.GaussSeidel(matA.data(), nEq);
			
			// for(int jEq=0; jEq<nEq; ++jEq){
				// double tmp_val = 0.0;
				// for(int iEq=0; iEq<nEq; ++iEq){
					// // tmp_val += A[nEq*jEq+iEq]*var.tmp_B[nEq*i+iEq];
					// tmp_val += matA[jEq*nEq+iEq]*var.tmp_B[nEq*i+iEq];
				// }
				// var.tmp_X[nEq*i+jEq] = tmp_val;
				
			// }
			
			cellVar_i[id_p] += var.tmp_X[nEq*i+0];
			cellVar_i[id_u] += var.tmp_X[nEq*i+1];
			cellVar_i[id_v] += var.tmp_X[nEq*i+2];
			cellVar_i[id_w] += var.tmp_X[nEq*i+3];
			cellVar_i[id_T] += var.tmp_X[nEq*i+4];
			
			for(int jEq=0; jEq<nEq; ++jEq){
				double tmp_resi = var.tmp_X[nEq*i+jEq]*volume;
				fieldVar[id_resi] += tmp_resi*tmp_resi;
				tmp_volume += volume;
			}
		}
		controls.log.pop();
		controls.log.push("calc14");
		if(size>1){
			vector<double> tmp_send(2,0.0);
			tmp_send[0] = (var.fields[id_resi]);
			tmp_send[1] = (tmp_volume);
			vector<double> tmp_recv(2,0.0);
			MPI_Allreduce(tmp_send.data(), tmp_recv.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			var.fields[id_resi] = tmp_recv[0];
			tmp_volume = tmp_recv[1];
		}
		var.fields[id_resi] = sqrt(var.fields[id_resi]/tmp_volume);
		
		var.fields[controls.fieldVar["time"].id] +=
		var.fields[controls.fieldVar["time-step"].id];
		
		controls.log.pop();
	}
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}