
#include "../../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	using Bound_Face_Funct_type = function<int(double* faces)>;
	using Bound_Cell_Funct_type = function<int(double* cells)>;
	
	Bound_Cell_Funct_type setCellFunction;
	Bound_Face_Funct_type setFaceFunction;
	
	{
		int id_p = controls.cellVar["pressure"].id;
		int id_u = controls.cellVar["x-velocity"].id;
		int id_v = controls.cellVar["y-velocity"].id;
		int id_w = controls.cellVar["z-velocity"].id;
		int id_T = controls.cellVar["temperature"].id;
		
		int nSp = controls.cellVar["mass fraction"].sub_name.size();
		vector<int> id_Y(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name = controls.cellVar["mass fraction"].sub_name[i];
			id_Y[i] = controls.cellVar[tmp_name].id;
		}
		
		int id_rho = controls.cellVar["density"].id;
		int id_c = controls.cellVar["speed of sound"].id;
		int id_Ht = controls.cellVar["total enthalpy"].id;
		
		int id_drhodp = controls.cellVar["density diff with pressure"].id;
		int id_drhodT = controls.cellVar["density diff with temperature"].id;
		int id_dHtdp = controls.cellVar["total enthalpy diff with pressure"].id;
		int id_dHtdT = controls.cellVar["total enthalpy diff with temperature"].id;
		vector<int> id_drhodY(nSp);
		vector<int> id_dHtdY(nSp);
		for(int i=0; i<nSp; ++i){
			string tmp_name1 = controls.cellVar["density diff with mass fraction"].sub_name[i];
			string tmp_name2 = controls.cellVar["total enthalpy diff with mass fraction"].sub_name[i];
			id_drhodY[i] = controls.cellVar[tmp_name1].id;
			id_dHtdY[i] = controls.cellVar[tmp_name2].id;
		}
		
		int id_mu = controls.cellVar["viscosity"].id;
		
		setCellFunction = [&solver,id_p,id_u,id_v,id_w,id_T,id_Y,nSp,
		id_rho,id_c,id_Ht,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_mu] (
		double* cells) ->int {
			double p = cells[id_p];
			double u = cells[id_u];
			double v = cells[id_v];
			double w = cells[id_w];
			double T = cells[id_T];
			double Y[nSp], alpha[nSp];
			for(int i=0; i<nSp; ++i){
				Y[i] = cells[id_Y[i]];
			}
			double Y_sum = 0.0;
			for(int i=0; i<nSp-1; ++i){
				Y_sum += Y[i];
			}
			Y[nSp-1] = 1.0-Y_sum;
			cells[id_Y[nSp-1]] = Y[nSp-1];
			
			double rhoi[nSp], ci[nSp], Hti[nSp], 
			drhodpi[nSp], drhodTi[nSp], dHtdpi[nSp], 
			dHtdTi[nSp];
			solver.eosNASG(
			621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
			p,u,v,w,T,rhoi[0],ci[0],Hti[0], 
			drhodpi[0], drhodTi[0], dHtdpi[0], dHtdTi[0]);
			solver.eosIdeal(
			717.0,1.4,
			p,u,v,w,T,rhoi[1],ci[1],Hti[1], 
			drhodpi[1], drhodTi[1], dHtdpi[1], dHtdTi[1]);
			
		
			double rho = 0.0;
			for(int ns=0; ns<nSp; ++ns){
				rho += Y[ns]/rhoi[ns];
			}
			rho = 1.0/rho;
			
			for(int ns=0; ns<nSp; ++ns){
				alpha[ns] = rho*Y[ns]/rhoi[ns];
			}
			
			double Ht=0.0;
			double drhodp=0.0;
			double drhodT=0.0;
			double dHtdp=0.0;
			double dHtdT=0.0;
			double drhodY[nSp];
			double dHtdY[nSp];
			for(int ns=0; ns<nSp; ++ns){
				Ht += Y[ns]*Hti[ns];
				drhodp += rho*rho*(Y[ns]/rhoi[ns]/rhoi[ns]*drhodpi[ns]);
				drhodT += rho*rho*(Y[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
				dHtdp += Y[ns]*dHtdpi[ns];
				dHtdT += Y[ns]*dHtdTi[ns];
				drhodY[ns] = ( -rho*rho*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]) );
				dHtdY[ns] = ( Hti[ns]-Hti[nSp-1] );
			}
			
			double c = drhodp + 1.0/rho*drhodT/dHtdT*(1.0-rho*dHtdp);
			c = sqrt( 1.0 / c );
			
			cells[id_rho] = rho;
			cells[id_c] = c;
			cells[id_Ht] = Ht;
			cells[id_drhodp] = drhodp;
			cells[id_drhodT] = drhodT;
			cells[id_dHtdp] = dHtdp;
			cells[id_dHtdT] = dHtdT;
			for(int ns=0; ns<nSp; ++ns){
				cells[id_drhodY[ns]] = drhodY[ns];
				cells[id_dHtdY[ns]] = dHtdY[ns];
			}
			
			
			
			// 점성계수 관련
			cells[id_mu] = 0.0;
			cells[id_mu] += alpha[0]*1.e-3;
			cells[id_mu] += alpha[1]*1.e-5;
			// for(int ns=0; ns<nSp; ++ns){
				// cells[id_mu] += alpha[ns]*
			// }
			
			
			return 0;
		};
		calcCellAddiVal.push_back(setCellFunction);
	}
	
	// setFaceFunction = [id_inp, id_L_out, id_R_out, value](
	// double* faces) ->int {
		// faces[id_L_out] = value;
		// faces[id_R_out] = value;
		// return 0;
	// };
	for(int ii=0; ii<2; ++ii)
	{
		int nSp;
		
		int id_p;
		int id_u;
		int id_v;
		int id_w;
		int id_T;
		
		vector<int> id_Y(nSp);
		
		int id_rho;
		int id_c;
		int id_Ht;
		
		if(ii==0){
				
			nSp = controls.faceVar["left mass fraction"].sub_name.size();
			
			id_p = controls.faceVar["left pressure"].id;
			id_u = controls.faceVar["left x-velocity"].id;
			id_v = controls.faceVar["left y-velocity"].id;
			id_w = controls.faceVar["left z-velocity"].id;
			id_T = controls.faceVar["left temperature"].id;
			
			for(int i=0; i<nSp; ++i){
				string tmp_name = controls.faceVar["left mass fraction"].sub_name[i];
				id_Y[i] = controls.faceVar[tmp_name].id;
			}
			
			id_rho = controls.faceVar["left density"].id;
			id_c = controls.faceVar["left speed of sound"].id;
			id_Ht = controls.faceVar["left total enthalpy"].id;
			
		}
		else{
			nSp = controls.faceVar["right mass fraction"].sub_name.size();
			
			id_p = controls.faceVar["right pressure"].id;
			id_u = controls.faceVar["right x-velocity"].id;
			id_v = controls.faceVar["right y-velocity"].id;
			id_w = controls.faceVar["right z-velocity"].id;
			id_T = controls.faceVar["right temperature"].id;
			
			for(int i=0; i<nSp; ++i){
				string tmp_name = controls.faceVar["right mass fraction"].sub_name[i];
				id_Y[i] = controls.faceVar[tmp_name].id;
			}
			
			id_rho = controls.faceVar["right density"].id;
			id_c = controls.faceVar["right speed of sound"].id;
			id_Ht = controls.faceVar["right total enthalpy"].id;
			
		}
		
		setFaceFunction = [&solver,nSp,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht] (
		double* cells) ->int {
			double p = cells[id_p];
			double u = cells[id_u];
			double v = cells[id_v];
			double w = cells[id_w];
			double T = cells[id_T];
			double Y[nSp];
			for(int i=0; i<nSp; ++i){
				Y[i] = cells[id_Y[i]];
			}
			double Y_sum = 0.0;
			for(int i=0; i<nSp-1; ++i){
				Y_sum += Y[i];
			}
			Y[nSp-1] = 1.0-Y_sum;
			cells[id_Y[nSp-1]] = Y[nSp-1];
			
			double rhoi[nSp], ci[nSp], Hti[nSp], 
			drhodpi[nSp], drhodTi[nSp], dHtdpi[nSp], 
			dHtdTi[nSp];
			solver.eosNASG(
			621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
			p,u,v,w,T,rhoi[0],ci[0],Hti[0], 
			drhodpi[0], drhodTi[0], dHtdpi[0], dHtdTi[0]);
			solver.eosIdeal(
			717.0,1.4,
			p,u,v,w,T,rhoi[1],ci[1],Hti[1], 
			drhodpi[1], drhodTi[1], dHtdpi[1], dHtdTi[1]);
			
		
			double rho = 0.0;
			for(int ns=0; ns<nSp; ++ns){
				rho += Y[ns]/rhoi[ns];
			}
			rho = 1.0/rho;
			
			double Ht=0.0;
			double drhodp=0.0;
			double drhodT=0.0;
			double dHtdp=0.0;
			double dHtdT=0.0;
			for(int ns=0; ns<nSp; ++ns){
				Ht += Y[ns]*Hti[ns];
				drhodp += rho*rho*(Y[ns]/rhoi[ns]/rhoi[ns]*drhodpi[ns]);
				drhodT += rho*rho*(Y[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
				dHtdp += Y[ns]*dHtdpi[ns];
				dHtdT += Y[ns]*dHtdTi[ns];
			}
			
			double c = drhodp + 1.0/rho*drhodT/dHtdT*(1.0-rho*dHtdp);
			c = sqrt( 1.0 / c );
			
			cells[id_rho] = rho;
			cells[id_c] = c;
			cells[id_Ht] = Ht;
			
			return 0;
		};
		calcFaceAddiVal.push_back(setFaceFunction);
	}
	
}

