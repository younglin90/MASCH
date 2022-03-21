
#include "../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int nSp = controls.spName.size();
	
	int id_p = controls.getId_cellVar("pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	int id_T = controls.getId_cellVar("temperature");
	vector<int> id_Y;
	for(int i=0; i<controls.spName.size(); ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		id_Y.push_back(controls.getId_cellVar(tmp_name));
	}
	int id_drhodp = controls.getId_cellVar("partial-density-pressure");
	int id_drhodT = controls.getId_cellVar("partial-density-temperature");
	int id_dHtdp = controls.getId_cellVar("partial-total-enthalpy-pressure");
	int id_dHtdT = controls.getId_cellVar("partial-total-enthalpy-temperature");
	vector<int> id_drhodY, id_dHtdY;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_drhodY.push_back(controls.getId_cellVar("partial-density-mass-fraction-"+controls.spName[i]));
		id_dHtdY.push_back(controls.getId_cellVar("partial-total-enthalpy-mass-fraction-"+controls.spName[i]));
	}
	int id_rho = controls.getId_cellVar("density");
	int id_c = controls.getId_cellVar("speed-of-sound");
	int id_Ht = controls.getId_cellVar("total-enthalpy");
	// int id_mu = controls.getId_cellVar("viscosity");

	int id_pL = controls.getId_faceVar("left pressure"); int id_pR = controls.getId_faceVar("right pressure");
	int id_uL = controls.getId_faceVar("left x-velocity"); int id_uR = controls.getId_faceVar("right x-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity"); int id_vR = controls.getId_faceVar("right y-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity"); int id_wR = controls.getId_faceVar("right z-velocity");
	int id_TL = controls.getId_faceVar("left temperature"); int id_TR = controls.getId_faceVar("right temperature");
	vector<int> id_YL, id_YR;
	for(int i=0; i<controls.spName.size()-1; ++i){
		id_YL.push_back(controls.getId_faceVar("left mass-fraction-"+controls.spName[i]));
		id_YR.push_back(controls.getId_faceVar("right mass-fraction-"+controls.spName[i]));
	}
	int id_rhoL = controls.getId_faceVar("left density"); int id_rhoR = controls.getId_faceVar("right density");
	int id_cL = controls.getId_faceVar("left speed-of-sound"); int id_cR = controls.getId_faceVar("right speed-of-sound");
	int id_HtL = controls.getId_faceVar("left total-enthalpy"); int id_HtR = controls.getId_faceVar("right total-enthalpy");
	// int id_muF = controls.getId_faceVar("viscosity");
	
	vector<int> eos_type(nSp);
	vector<vector<double>> spInform(nSp,vector<double>(6,0.0));
	for(int i=0; i<controls.spName.size(); ++i){
		string tmp_name = controls.spName[i];
		if(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.type"] == "ideal"){
			eos_type[i] = 0;
			spInform[i][0] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.cv"]);
			spInform[i][1] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.gamma"]);
		}
		else if(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.type"] == "NASG"){
			eos_type[i] = 1;
			spInform[i][0] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.cv"]);
			spInform[i][1] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.gamma"]);
			spInform[i][2] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.pinf"]);
			spInform[i][3] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.b"]);
			spInform[i][4] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.q"]);
		}
		else{
			cout << "#WARNING : not defined EOS type" << endl;
		}
		
		
		if(controls.thermophysicalProperties[tmp_name+".transport.mu.type"] == "constant"){
			spInform[i][5] = stod(controls.thermophysicalProperties[tmp_name+".transport.mu.value"]);
		}
		else{
			cout << "#WARNING : not defined transport type" << endl;
		}
		
	}
	
		
	
	{
		calcCellAddiVal.push_back(
			[&solver,spInform,nSp,
			id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
			id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
			id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR,
			id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
			eos_type] (
			double* cells) ->int {
					
				auto spInf_ptr = spInform.data();
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
					if(eos_type[i]==0){
						solver.eosIdeal(
						// 717.0,1.4,
						spInf_ptr_i[0],spInf_ptr_i[1],
						p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
						drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
					}
					else if(eos_type[i]==1){
						// solver.eosNASG(
						// // 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
						// spInf_ptr_i[0],spInf_ptr_i[1],spInf_ptr_i[2],spInf_ptr_i[3],spInf_ptr_i[4],
						// p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
						// drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
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
				for(int ns=0; ns<nSp-1; ++ns){
					cells[id_drhodY[ns]] = drhodY[ns];
					cells[id_dHtdY[ns]] = dHtdY[ns];
				}
				
				
				return 0;
			}
		);
		
		
		
		
		
		
		
		
		
		calcFaceAddiVal.push_back( 
			[&solver,spInform,nSp,
			id_p,id_u,id_v,id_w,id_T,id_Y,id_rho,id_c,id_Ht,
			id_pL,id_uL,id_vL,id_wL,id_TL,id_YL,id_rhoL,id_cL,id_HtL,
			id_pR,id_uR,id_vR,id_wR,id_TR,id_YR,id_rhoR,id_cR,id_HtR,
			id_drhodp,id_drhodT,id_dHtdp,id_dHtdT,id_drhodY,id_dHtdY,
			eos_type] (
			double* faces) ->int {
				
				for(int jj=0; jj<2; ++jj)
				{
					auto spInf_ptr = spInform.data();
					double p = faces[id_pL];
					double u = faces[id_uL];
					double v = faces[id_vL];
					double w = faces[id_wL];
					double T = faces[id_TL];
					double Y[nSp], alpha[nSp];
					double Y_sum = 0.0;
					for(int i=0; i<nSp-1; ++i){
						Y[i] = faces[id_YL[i]];
						Y_sum += Y[i];
					}
					Y[nSp-1] = 1.0-Y_sum;
					
					if(jj==1){
						p = faces[id_pR];
						u = faces[id_uR];
						v = faces[id_vR];
						w = faces[id_wR];
						T = faces[id_TR];
						Y_sum = 0.0;
						for(int i=0; i<nSp-1; ++i){
							Y[i] = faces[id_YR[i]];
							Y_sum += Y[i];
						}
						Y[nSp-1] = 1.0-Y_sum;
					}
					
					double rhoi[nSp], ci[nSp], Hti[nSp], 
					drhodpi[nSp], drhodTi[nSp], dHtdpi[nSp], 
					dHtdTi[nSp];
					for(int i=0; i<nSp; ++i){
						auto spInf_ptr_i = spInf_ptr[i].data();
						if(eos_type[i]==0){
							solver.eosIdeal(
							// 717.0,1.4,
							spInf_ptr_i[0],spInf_ptr_i[1],
							p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
							drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
						}
						else if(eos_type[i]==1){
							// solver.eosNASG(
							// // 621780000.0,3610.0,1.19,6.7212e-4,-1177788.0,
							// spInf_ptr_i[2],spInf_ptr_i[0],spInf_ptr_i[1],spInf_ptr_i[3],spInf_ptr_i[4],
							// p,u,v,w,T,rhoi[i],ci[i],Hti[i], 
							// drhodpi[i], drhodTi[i], dHtdpi[i], dHtdTi[i]);
						}
					}
					
					double rho = 0.0;
					for(int i=0; i<nSp; ++i){
						rho += Y[i]/rhoi[i];
					}
					rho = 1.0/rho;
					
					for(int i=0; i<nSp; ++i){
						alpha[i] = rho*Y[i]/rhoi[i];
					}
					
					double Ht=0.0;
					double drhodp=0.0;
					double drhodT=0.0;
					double dHtdp=0.0;
					double dHtdT=0.0;
					double drhodY[nSp];
					double dHtdY[nSp];
					for(int i=0; i<nSp; ++i){
						Ht += Y[i]*Hti[i];
						drhodp += rho*rho*(Y[i]/rhoi[i]/rhoi[i]*drhodpi[i]);
						drhodT += rho*rho*(Y[i]/rhoi[i]/rhoi[i]*drhodTi[i]);
						dHtdp += Y[i]*dHtdpi[i];
						dHtdT += Y[i]*dHtdTi[i];
						drhodY[i] = ( -rho*rho*(1.0/rhoi[i]-1.0/rhoi[nSp-1]) );
						dHtdY[i] = ( Hti[i]-Hti[nSp-1] );
					}
					
					double c = drhodp + 1.0/rho*drhodT/dHtdT*(1.0-rho*dHtdp);
					c = sqrt( 1.0 / c );
					
					
					if(jj==0){
						faces[id_rhoL] = rho;
						faces[id_cL] = c;
						faces[id_HtL] = Ht;
					}
					else{
						faces[id_rhoR] = rho;
						faces[id_cR] = c;
						faces[id_HtR] = Ht;
					}
				}
				
				
				return 0;
			}
		);
	}
	
	
}

