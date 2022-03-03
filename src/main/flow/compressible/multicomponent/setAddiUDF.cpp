
#include "../../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	using Bound_Face_Funct_type = function<int(double* faces)>;
	using Bound_Cell_Funct_type = function<int(double* cells)>;
	
	Bound_Cell_Funct_type setCellFunction;
	Bound_Face_Funct_type setFaceFunction;
	
	

	// species 관련 재료들
	auto& thermoMap = controls.thermophysicalProperties;
	vector<int> id_spec_type;
	vector<vector<double>> spInf;
	for(int i=0; i<controls.nSp; ++i){
		spInf.push_back(vector<double>());
		// cout << controls.spName[i] << endl;
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.pinf";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.cv";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.gamma";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.b";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.q";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		{
			string name = controls.spName[i];
			name += ".transport.mu.value";
			if (thermoMap.find(name) != thermoMap.end()) {
				spInf.back().push_back(stod(thermoMap[name]));
			}
			else{
				spInf.back().push_back(0.0);
			}
		}
		
		{
			string name = controls.spName[i];
			name += ".thermodynamics.rho.type";
			if(
			controls.thermophysicalProperties[name] ==
			"ideal"){
				id_spec_type.push_back(0);
			}
			else if(
			controls.thermophysicalProperties[name] ==
			"NASG"){
				id_spec_type.push_back(1);
			}
		}
	}
	
	
	
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
		id_rho,id_c,id_Ht,spInf,id_spec_type,
		id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
		id_mu] (
		double* cells) ->int {
			auto spInf_ptr = spInf.data();
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
			for(int i=0; i<nSp; ++i){
				auto spInf_ptr_i = spInf_ptr[i].data();
				if(id_spec_type[i]==0){
					solver.eosIdeal(
					// 717.0,1.4,
					spInf_ptr_i[1],spInf_ptr_i[2],
					p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
					drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
				}
				else if(id_spec_type[i]==1){
					solver.eosNASG(
					// 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
					spInf_ptr_i[0],spInf_ptr_i[1],spInf_ptr_i[2],spInf_ptr_i[3],spInf_ptr_i[4],
					p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
					drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
				}
			}
			
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
			for(int i=0; i<nSp; ++i){
				cells[id_mu] += alpha[0]*spInf[i][5];
				// cout << spInf[i][5] << endl;
			}
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
		
		setFaceFunction = 
		[&solver,nSp,spInf,id_spec_type,
		id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht] (
		double* faces) ->int {
			auto spInf_ptr = spInf.data();
			double p = faces[id_p];
			double u = faces[id_u];
			double v = faces[id_v];
			double w = faces[id_w];
			double T = faces[id_T];
			double Y[nSp];
			for(int i=0; i<nSp; ++i){
				Y[i] = faces[id_Y[i]];
			}
			double Y_sum = 0.0;
			for(int i=0; i<nSp-1; ++i){
				Y_sum += Y[i];
			}
			Y[nSp-1] = 1.0-Y_sum;
			faces[id_Y[nSp-1]] = Y[nSp-1];
			
			double rhoi[nSp], ci[nSp], Hti[nSp], 
			drhodpi[nSp], drhodTi[nSp], dHtdpi[nSp], 
			dHtdTi[nSp];

			for(int i=0; i<nSp; ++i){
				auto spInf_ptr_i = spInf_ptr[i].data();
				if(id_spec_type[i]==0){
					solver.eosIdeal(
					// 717.0,1.4,
					spInf_ptr_i[1],spInf_ptr_i[2],
					p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
					drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
				}
				else if(id_spec_type[i]==1){
					solver.eosNASG(
					// 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
					spInf_ptr_i[0],spInf_ptr_i[1],spInf_ptr_i[2],spInf_ptr_i[3],spInf_ptr_i[4],
					p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
					drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
				}
			}
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
			
			faces[id_rho] = rho;
			faces[id_c] = c;
			faces[id_Ht] = Ht;
			
			// if(c<=1.e-8){
			// cout << id_rho << " " << id_c << " " << id_Ht << endl;
			// }
			
			return 0;
		};
		calcFaceAddiVal.push_back(setFaceFunction);
	}
	
}

