

#include "../../../../others/solvers.h"

void MASCH_Solver::setHOReconFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	MASCH_Load load;
	
	int nSp = controls.spName.size();
	
	int id_c = controls.getId_cellVar("speed-of-sound");
	
	int id_p = controls.getId_cellVar("pressure");
	int id_pL = controls.getId_faceVar("left pressure");
	int id_pR = controls.getId_faceVar("right pressure");
	int id_dpdx = controls.getId_cellVar("x-gradient pressure");
	int id_dpdy = controls.getId_cellVar("y-gradient pressure");
	int id_dpdz = controls.getId_cellVar("z-gradient pressure");
	int id_pMax = controls.getId_cellVar("maximum pressure");
	int id_pMin = controls.getId_cellVar("minimum pressure");

	int id_u = controls.getId_cellVar("x-velocity");
	int id_uL = controls.getId_faceVar("left x-velocity");
	int id_uR = controls.getId_faceVar("right x-velocity");
	int id_dudx = controls.getId_cellVar("x-gradient x-velocity");
	int id_dudy = controls.getId_cellVar("y-gradient x-velocity");
	int id_dudz = controls.getId_cellVar("z-gradient x-velocity");
	int id_uMax = controls.getId_cellVar("maximum x-velocity");
	int id_uMin = controls.getId_cellVar("minimum x-velocity");
	
	int id_v = controls.getId_cellVar("y-velocity");
	int id_vL = controls.getId_faceVar("left y-velocity");
	int id_vR = controls.getId_faceVar("right y-velocity");
	int id_dvdx = controls.getId_cellVar("x-gradient y-velocity");
	int id_dvdy = controls.getId_cellVar("y-gradient y-velocity");
	int id_dvdz = controls.getId_cellVar("z-gradient y-velocity");
	int id_vMax = controls.getId_cellVar("maximum y-velocity");
	int id_vMin = controls.getId_cellVar("minimum y-velocity");
	
	int id_w = controls.getId_cellVar("z-velocity");
	int id_wL = controls.getId_faceVar("left z-velocity");
	int id_wR = controls.getId_faceVar("right z-velocity");
	int id_dwdx = controls.getId_cellVar("x-gradient z-velocity");
	int id_dwdy = controls.getId_cellVar("y-gradient z-velocity");
	int id_dwdz = controls.getId_cellVar("z-gradient z-velocity");
	int id_wMax = controls.getId_cellVar("maximum z-velocity");
	int id_wMin = controls.getId_cellVar("minimum z-velocity");
	
	int id_T = controls.getId_cellVar("temperature");
	int id_TL = controls.getId_faceVar("left temperature");
	int id_TR = controls.getId_faceVar("right temperature");
	int id_dTdx = controls.getId_cellVar("x-gradient temperature");
	int id_dTdy = controls.getId_cellVar("y-gradient temperature");
	int id_dTdz = controls.getId_cellVar("z-gradient temperature");
	int id_TMax = controls.getId_cellVar("maximum temperature");
	int id_TMin = controls.getId_cellVar("minimum temperature");
	
	vector<int> id_Y, id_YL, id_YR, id_dYdx, id_dYdy, id_dYdz, id_YMax, id_YMin;
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+controls.spName[i]);
		id_Y.push_back(controls.getId_cellVar(tmp_name));
		id_dYdx.push_back(controls.getId_cellVar("x-gradient "+tmp_name));
		id_dYdy.push_back(controls.getId_cellVar("y-gradient "+tmp_name));
		id_dYdz.push_back(controls.getId_cellVar("z-gradient "+tmp_name));
		id_YMax.push_back(controls.getId_cellVar("maximum "+tmp_name));
		id_YMin.push_back(controls.getId_cellVar("minimum "+tmp_name));
	}
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("left mass-fraction-"+controls.spName[i]);
		id_YL.push_back(controls.getId_faceVar(tmp_name));
	}
	for(int i=0; i<controls.spName.size()-1; ++i){
		string tmp_name = ("right mass-fraction-"+controls.spName[i]);
		id_YR.push_back(controls.getId_faceVar(tmp_name));
	}
	
	int id_rho = controls.getId_cellVar("density");
	int id_rhoL = controls.getId_faceVar("left density");
	int id_rhoR = controls.getId_faceVar("right density");
	
	int id_Ht = controls.getId_cellVar("total-enthalpy");
	int id_HtL = controls.getId_faceVar("left total-enthalpy");
	int id_HtR = controls.getId_faceVar("right total-enthalpy");
	
	
	// int id_dtrho = controls.getId_faceVar("time-step-density");
	
	// int id_UnF = controls.getId_faceVar("contravariant-velocity");
	
	int id_dt = controls.getId_fieldVar("time-step");
	
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
	
	int id_dLR = controls.getId_faceVar("distance of between left and right cell"); 
	int id_wd = controls.getId_faceVar("distance weight"); 
	int id_xLR = controls.getId_faceVar("x distance of between left and right cell");
	int id_yLR = controls.getId_faceVar("y distance of between left and right cell");
	int id_zLR = controls.getId_faceVar("z distance of between left and right cell");
	
	int id_xLF = controls.getId_faceVar("x distance of between left cell and face");
	int id_yLF = controls.getId_faceVar("y distance of between left cell and face");
	int id_zLF = controls.getId_faceVar("z distance of between left cell and face");
	int id_xRF = controls.getId_faceVar("x distance of between right cell and face");
	int id_yRF = controls.getId_faceVar("y distance of between right cell and face");
	int id_zRF = controls.getId_faceVar("z distance of between right cell and face");
	
	int id_xLNv = controls.getId_faceVar("left cell to face normal vector shortest x distance");
	int id_yLNv = controls.getId_faceVar("left cell to face normal vector shortest y distance");
	int id_zLNv = controls.getId_faceVar("left cell to face normal vector shortest z distance");
	int id_xRNv = controls.getId_faceVar("right cell to face normal vector shortest x distance");
	int id_yRNv = controls.getId_faceVar("right cell to face normal vector shortest y distance");
	int id_zRNv = controls.getId_faceVar("right cell to face normal vector shortest z distance");

	
	int id_alpha = controls.getId_faceVar("cosine angle of between face normal and cells");
	
	int id_nLRx = controls.getId_faceVar("x unit normal of between left and right cell");
	int id_nLRy = controls.getId_faceVar("y unit normal of between left and right cell");
	int id_nLRz = controls.getId_faceVar("z unit normal of between left and right cell");
	
	
	using funcType = double(MASCH_NVD::*)(double,double,double,double,double);
	vector<funcType> NVD_scheme;
	{
		string type = controls.fvSchemeMap["highOrderScheme.p"];
		{
			if(type=="minmod") { NVD_scheme.push_back(&MASCH_NVD::minmod); }
			else if(type=="vanLeer") { NVD_scheme.push_back(&MASCH_NVD::vanLeer); }
			else if(type=="QUICK") { NVD_scheme.push_back(&MASCH_NVD::QUICK); }
			else if(type=="boundedCD") { NVD_scheme.push_back(&MASCH_NVD::boundedCD); }
			else if(type=="OSHER") { NVD_scheme.push_back(&MASCH_NVD::OSHER); }
			else if(type=="SMART") { NVD_scheme.push_back(&MASCH_NVD::SMART); }
			else if(type=="modifiedSMART") { NVD_scheme.push_back(&MASCH_NVD::modifiedSMART); }
			else if(type=="STOIC") { NVD_scheme.push_back(&MASCH_NVD::STOIC); }
			else if(type=="modifiedSTOIC") { NVD_scheme.push_back(&MASCH_NVD::modifiedSTOIC); }
			else if(type=="MUSCL") { NVD_scheme.push_back(&MASCH_NVD::MUSCL); }
			else if(type=="SUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::SUPERBEE); }
			else if(type=="modifiedSUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::modifiedSUPERBEE); }
			
			// else if(type=="HRIC") { NVD_scheme.push_back(&MASCH_NVD::HRIC); }
			// else if(type=="CICSAM") { NVD_scheme.push_back(&MASCH_NVD::CICSAM); }
			// else if(type=="STACS") { NVD_scheme.push_back(&MASCH_NVD::STACS); }
			// else if(type=="FBICS") { NVD_scheme.push_back(&MASCH_NVD::FBICS); }
			// else if(type=="SAISH") { NVD_scheme.push_back(&MASCH_NVD::SAISH); }
			// else if(type=="MSTACS") { NVD_scheme.push_back(&MASCH_NVD::MSTACS); }
			else{ NVD_scheme.push_back(&MASCH_NVD::none); }
		}
	}
	{
		vector<string> types = load.extractVector(controls.fvSchemeMap["highOrderScheme.U"]);
		for(int i=0; i<types.size(); ++i){
			string type = types[i];
			if(type=="minmod") { NVD_scheme.push_back(&MASCH_NVD::minmod); }
			else if(type=="vanLeer") { NVD_scheme.push_back(&MASCH_NVD::vanLeer); }
			else if(type=="QUICK") { NVD_scheme.push_back(&MASCH_NVD::QUICK); }
			else if(type=="boundedCD") { NVD_scheme.push_back(&MASCH_NVD::boundedCD); }
			else if(type=="OSHER") { NVD_scheme.push_back(&MASCH_NVD::OSHER); }
			else if(type=="SMART") { NVD_scheme.push_back(&MASCH_NVD::SMART); }
			else if(type=="modifiedSMART") { NVD_scheme.push_back(&MASCH_NVD::modifiedSMART); }
			else if(type=="STOIC") { NVD_scheme.push_back(&MASCH_NVD::STOIC); }
			else if(type=="modifiedSTOIC") { NVD_scheme.push_back(&MASCH_NVD::modifiedSTOIC); }
			else if(type=="MUSCL") { NVD_scheme.push_back(&MASCH_NVD::MUSCL); }
			else if(type=="SUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::SUPERBEE); }
			else if(type=="modifiedSUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::modifiedSUPERBEE); }
			
			// else if(type=="HRIC") { NVD_scheme.push_back(&MASCH_NVD::HRIC); }
			// else if(type=="CICSAM") { NVD_scheme.push_back(&MASCH_NVD::CICSAM); }
			// else if(type=="STACS") { NVD_scheme.push_back(&MASCH_NVD::STACS); }
			// else if(type=="FBICS") { NVD_scheme.push_back(&MASCH_NVD::FBICS); }
			// else if(type=="SAISH") { NVD_scheme.push_back(&MASCH_NVD::SAISH); }
			// else if(type=="MSTACS") { NVD_scheme.push_back(&MASCH_NVD::MSTACS); }
			else{ NVD_scheme.push_back(&MASCH_NVD::none); }
		}
	}
	{
		string type = controls.fvSchemeMap["highOrderScheme.T"];
		{
			if(type=="minmod") { NVD_scheme.push_back(&MASCH_NVD::minmod); }
			else if(type=="vanLeer") { NVD_scheme.push_back(&MASCH_NVD::vanLeer); }
			else if(type=="QUICK") { NVD_scheme.push_back(&MASCH_NVD::QUICK); }
			else if(type=="boundedCD") { NVD_scheme.push_back(&MASCH_NVD::boundedCD); }
			else if(type=="OSHER") { NVD_scheme.push_back(&MASCH_NVD::OSHER); }
			else if(type=="SMART") { NVD_scheme.push_back(&MASCH_NVD::SMART); }
			else if(type=="modifiedSMART") { NVD_scheme.push_back(&MASCH_NVD::modifiedSMART); }
			else if(type=="STOIC") { NVD_scheme.push_back(&MASCH_NVD::STOIC); }
			else if(type=="modifiedSTOIC") { NVD_scheme.push_back(&MASCH_NVD::modifiedSTOIC); }
			else if(type=="MUSCL") { NVD_scheme.push_back(&MASCH_NVD::MUSCL); }
			else if(type=="SUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::SUPERBEE); }
			else if(type=="modifiedSUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::modifiedSUPERBEE); }
			
			// else if(type=="HRIC") { NVD_scheme.push_back(&MASCH_NVD::HRIC); }
			// else if(type=="CICSAM") { NVD_scheme.push_back(&MASCH_NVD::CICSAM); }
			// else if(type=="STACS") { NVD_scheme.push_back(&MASCH_NVD::STACS); }
			// else if(type=="FBICS") { NVD_scheme.push_back(&MASCH_NVD::FBICS); }
			// else if(type=="SAISH") { NVD_scheme.push_back(&MASCH_NVD::SAISH); }
			// else if(type=="MSTACS") { NVD_scheme.push_back(&MASCH_NVD::MSTACS); }
			else{ NVD_scheme.push_back(&MASCH_NVD::none); }
		}
	}
	{
		vector<string> types = load.extractVector(controls.fvSchemeMap["highOrderScheme.Y"]);
		for(int i=0; i<types.size(); ++i){
			string type = types[i];
			if(type=="minmod") { NVD_scheme.push_back(&MASCH_NVD::minmod); }
			else if(type=="vanLeer") { NVD_scheme.push_back(&MASCH_NVD::vanLeer); }
			else if(type=="QUICK") { NVD_scheme.push_back(&MASCH_NVD::QUICK); }
			else if(type=="boundedCD") { NVD_scheme.push_back(&MASCH_NVD::boundedCD); }
			else if(type=="OSHER") { NVD_scheme.push_back(&MASCH_NVD::OSHER); }
			else if(type=="SMART") { NVD_scheme.push_back(&MASCH_NVD::SMART); }
			else if(type=="modifiedSMART") { NVD_scheme.push_back(&MASCH_NVD::modifiedSMART); }
			else if(type=="STOIC") { NVD_scheme.push_back(&MASCH_NVD::STOIC); }
			else if(type=="modifiedSTOIC") { NVD_scheme.push_back(&MASCH_NVD::modifiedSTOIC); }
			else if(type=="MUSCL") { NVD_scheme.push_back(&MASCH_NVD::MUSCL); }
			else if(type=="SUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::SUPERBEE); }
			else if(type=="modifiedSUPERBEE") { NVD_scheme.push_back(&MASCH_NVD::modifiedSUPERBEE); }
			
			else if(type=="HRIC") { NVD_scheme.push_back(&MASCH_NVD::HRIC); }
			else if(type=="CICSAM") { NVD_scheme.push_back(&MASCH_NVD::CICSAM); }
			else if(type=="STACS") { NVD_scheme.push_back(&MASCH_NVD::STACS); }
			else if(type=="FBICS") { NVD_scheme.push_back(&MASCH_NVD::FBICS); }
			else if(type=="SAISH") { NVD_scheme.push_back(&MASCH_NVD::SAISH); }
			else if(type=="MSTACS") { NVD_scheme.push_back(&MASCH_NVD::MSTACS); }
			else{ NVD_scheme.push_back(&MASCH_NVD::none); }
		}
	}
	
	if(NVD_scheme.size() != 5+nSp-1) cout << "#WARNING : NVD_scheme not matching nEq" << endl;
	
	
	
	
	int id_lim_p = controls.getId_cellVar("limiter-unstructured pressure");
	int id_lim_u = controls.getId_cellVar("limiter-unstructured x-velocity");
	int id_lim_v = controls.getId_cellVar("limiter-unstructured y-velocity");
	int id_lim_w = controls.getId_cellVar("limiter-unstructured z-velocity");
	controls.limiterNamesForUnst.push_back(id_lim_p);
	controls.limiterNamesForUnst.push_back(id_lim_u);
	controls.limiterNamesForUnst.push_back(id_lim_v);
	controls.limiterNamesForUnst.push_back(id_lim_w);
	
	
	
	{
		calcHO_FaceVal.push_back(
			[&solver,NVD_scheme,
			id_pMax,id_pMin,id_uMax,id_uMin,id_vMax,id_vMin,id_wMax,id_wMin,
			id_TMax,id_TMin,id_YMax,id_YMin,
			id_dt,id_nx,id_ny,id_nz,id_area,id_dLR,id_xLR,id_yLR,id_zLR,
			id_p,id_pL,id_pR,id_dpdx,id_dpdy,id_dpdz,
			id_u,id_uL,id_uR,id_dudx,id_dudy,id_dudz,
			id_v,id_vL,id_vR,id_dvdx,id_dvdy,id_dvdz,
			id_w,id_wL,id_wR,id_dwdx,id_dwdy,id_dwdz,
			id_T,id_TL,id_TR,id_dTdx,id_dTdy,id_dTdz,
			id_Y,id_YL,id_YR,id_dYdx,id_dYdy,id_dYdz,
			id_wd,id_rho,id_rhoL,id_rhoR,id_Ht,id_HtL,id_HtR,
			id_c,nSp,id_alpha,id_nLRx,id_nLRy,id_nLRz,
			id_xLNv,id_yLNv,id_zLNv,id_xRNv,id_yRNv,id_zRNv,
			id_lim_p,id_lim_u,id_lim_v,id_lim_w,
			id_xLF,id_yLF,id_zLF,id_xRF,id_yRF,id_zRF](
			double* fields, double* cellsL, double* cellsR, double* faces)->int{
				auto NVD_scheme_i = NVD_scheme.data();
				double nvec[3];
				nvec[0] = faces[id_nx]; nvec[1] = faces[id_ny]; nvec[2] = faces[id_nz];
				double area = faces[id_area];
				double dLR = faces[id_dLR]; 
				double dAlpha = faces[id_alpha];
				double nLR[3];
				nLR[0] = faces[id_nLRx]; nLR[1] = faces[id_nLRy]; nLR[2] = faces[id_nLRz];
				double LNv[3];
				LNv[0] = faces[id_xLNv]; LNv[1] = faces[id_yLNv]; LNv[2] = faces[id_zLNv];
				double RNv[3];
				RNv[0] = faces[id_xRNv]; RNv[1] = faces[id_yRNv]; RNv[2] = faces[id_zRNv];
				double xyzLR[3];
				xyzLR[0] = faces[id_xLR]; xyzLR[1] = faces[id_yLR]; xyzLR[2] = faces[id_zLR];
				
				double xyzLF[3],xyzRF[3];
				xyzLF[0] = faces[id_xLF]; xyzLF[1] = faces[id_yLF]; xyzLF[2] = faces[id_zLF];
				xyzRF[0] = faces[id_xRF]; xyzRF[1] = faces[id_yRF]; xyzRF[2] = faces[id_zRF];
				
				
				double dt = fields[id_dt];
				// double rhoL = cellsL[id_rho]; double rhoR = cellsR[id_rho];
				
				double cL = cellsL[id_c]; double cR = cellsR[id_c];
				double pL = cellsL[id_p]; double pR = cellsR[id_p];
				double dpdxL = cellsL[id_dpdx]; double dpdxR = cellsR[id_dpdx];
				double dpdyL = cellsL[id_dpdy]; double dpdyR = cellsR[id_dpdy];
				double dpdzL = cellsL[id_dpdz]; double dpdzR = cellsR[id_dpdz];
				double pL_max = cellsL[id_pMax]; double pR_max = cellsR[id_pMax];
				double pL_min = cellsL[id_pMin]; double pR_min = cellsR[id_pMin];
				
				double uL = cellsL[id_u]; double uR = cellsR[id_u];
				double dudxL = cellsL[id_dudx]; double dudxR = cellsR[id_dudx];
				double dudyL = cellsL[id_dudy]; double dudyR = cellsR[id_dudy];
				double dudzL = cellsL[id_dudz]; double dudzR = cellsR[id_dudz];
				double uL_max = cellsL[id_uMax]; double uR_max = cellsR[id_uMax];
				double uL_min = cellsL[id_uMin]; double uR_min = cellsR[id_uMin];
				
				double vL = cellsL[id_v]; double vR = cellsR[id_v];
				double dvdxL = cellsL[id_dvdx]; double dvdxR = cellsR[id_dvdx];
				double dvdyL = cellsL[id_dvdy]; double dvdyR = cellsR[id_dvdy];
				double dvdzL = cellsL[id_dvdz]; double dvdzR = cellsR[id_dvdz];
				double vL_max = cellsL[id_vMax]; double vR_max = cellsR[id_vMax];
				double vL_min = cellsL[id_vMin]; double vR_min = cellsR[id_vMin];
				
				double wL = cellsL[id_w]; double wR = cellsR[id_w];
				double dwdxL = cellsL[id_dwdx]; double dwdxR = cellsR[id_dwdx];
				double dwdyL = cellsL[id_dwdy]; double dwdyR = cellsR[id_dwdy];
				double dwdzL = cellsL[id_dwdz]; double dwdzR = cellsR[id_dwdz];
				double wL_max = cellsL[id_wMax]; double wR_max = cellsR[id_wMax];
				double wL_min = cellsL[id_wMin]; double wR_min = cellsR[id_wMin];
				
				double TL = cellsL[id_T]; double TR = cellsR[id_T];
				double dTdxL = cellsL[id_dTdx]; double dTdxR = cellsR[id_dTdx];
				double dTdyL = cellsL[id_dTdy]; double dTdyR = cellsR[id_dTdy];
				double dTdzL = cellsL[id_dTdz]; double dTdzR = cellsR[id_dTdz];
				double TL_max = cellsL[id_TMax]; double TR_max = cellsR[id_TMax];
				double TL_min = cellsL[id_TMin]; double TR_min = cellsR[id_TMin];
				
				double YL[nSp]; double YR[nSp];
				double dYdxL[nSp]; double dYdxR[nSp];
				double dYdyL[nSp]; double dYdyR[nSp];
				double dYdzL[nSp]; double dYdzR[nSp];
				double YL_max[nSp]; double YR_max[nSp];
				double YL_min[nSp]; double YR_min[nSp];
				for(int i=0; i<nSp-1; ++i){
					YL[i] = cellsL[id_Y[i]]; YR[i] = cellsR[id_Y[i]];
					dYdxL[i] = cellsL[id_dYdx[i]]; dYdxR[i] = cellsR[id_dYdx[i]];
					dYdyL[i] = cellsL[id_dYdy[i]]; dYdyR[i] = cellsR[id_dYdy[i]];
					dYdzL[i] = cellsL[id_dYdz[i]]; dYdzR[i] = cellsR[id_dYdz[i]];
					YL_max[i] = cellsL[id_YMax[i]]; YR_max[i] = cellsR[id_YMax[i]];
					YL_min[i] = cellsL[id_YMin[i]]; YR_min[i] = cellsR[id_YMin[i]];
				}
				
				// double UnL = uL*nvec[0]+vL*nvec[1]+wL*nvec[2];
				// double UnR = uR*nvec[0]+vR*nvec[1]+wR*nvec[2];
				
				// double delphiL[5+nSp-1],delphiR[5+nSp-1];
				// double limGrad_phiL[5+nSp-1],limGrad_phiR[5+nSp-1];
				double dx = sqrt(faces[id_area]);
				double eta = 1000.0*dx*dx*dx;
				
				// delphiL[0] = dpdxL*LNv[0] + dpdyL*LNv[1] + dpdzL*LNv[2];
				// delphiR[0] = dpdxR*RNv[0] + dpdyR*RNv[1] + dpdzR*RNv[2];
				// limGrad_phiL[0] = solver.limiter_MLP(pL,pL_max,pL_min,delphiL[0], eta);
				// limGrad_phiR[0] = solver.limiter_MLP(pR,pR_max,pR_min,delphiR[0], eta);
				
				// delphiL[1] = dudxL*LNv[0] + dudyL*LNv[1] + dudzL*LNv[2];
				// delphiR[1] = dudxR*RNv[0] + dudyR*RNv[1] + dudzR*RNv[2];
				// limGrad_phiL[1] = solver.limiter_MLP(uL,uL_max,uL_min,delphiL[1], eta);
				// limGrad_phiR[1] = solver.limiter_MLP(uR,uR_max,uR_min,delphiR[1], eta);
				
				// delphiL[2] = dvdxL*LNv[0] + dvdyL*LNv[1] + dvdzL*LNv[2];
				// delphiR[2] = dvdxR*RNv[0] + dvdyR*RNv[1] + dvdzR*RNv[2];
				// limGrad_phiL[2] = solver.limiter_MLP(vL,vL_max,vL_min,delphiL[2], eta);
				// limGrad_phiR[2] = solver.limiter_MLP(vR,vR_max,vR_min,delphiR[2], eta);
				
				// delphiL[3] = dwdxL*LNv[0] + dwdyL*LNv[1] + dwdzL*LNv[2];
				// delphiR[3] = dwdxR*RNv[0] + dwdyR*RNv[1] + dwdzR*RNv[2];
				// limGrad_phiL[3] = solver.limiter_MLP(wL,wL_max,wL_min,delphiL[3], eta);
				// limGrad_phiR[3] = solver.limiter_MLP(wR,wR_max,wR_min,delphiR[3], eta);
				
				// delphiL[4] = dTdxL*LNv[0] + dTdyL*LNv[1] + dTdzL*LNv[2];
				// delphiR[4] = dTdxR*RNv[0] + dTdyR*RNv[1] + dTdzR*RNv[2];
				// limGrad_phiL[4] = solver.limiter_MLP(TL,TL_max,TL_min,delphiL[4], eta);
				// limGrad_phiR[4] = solver.limiter_MLP(TR,TR_max,TR_min,delphiR[4], eta);
				
				// for(int i=0; i<nSp-1; ++i){
					// delphiL[5+i] = dYdxL[i]*LNv[0] + dYdyL[i]*LNv[1] + dYdzL[i]*LNv[2];
					// delphiR[5+i] = dYdxR[i]*RNv[0] + dYdyR[i]*RNv[1] + dYdzR[i]*RNv[2];
					// limGrad_phiL[5+i] = solver.limiter_MLP(YL[i],YL_max[i],YL_min[i],delphiL[5+i], eta);
					// limGrad_phiR[5+i] = solver.limiter_MLP(YR[i],YR_max[i],YR_min[i],delphiR[5+i], eta);
				// }
				
				
				
				cellsL[id_lim_p] = min(cellsL[id_lim_p], 
					solver.limiter_MLP(pL,pL_max,pL_min, 
						dpdxL*xyzLF[0] + dpdyL*xyzLF[1] + dpdzL*xyzLF[2], eta));
				cellsR[id_lim_p] = min(cellsR[id_lim_p], 
					solver.limiter_MLP(pR,pR_max,pR_min, 
						dpdxR*xyzRF[0] + dpdyR*xyzRF[1] + dpdzR*xyzRF[2], eta));
						
				cellsL[id_lim_u] = min(cellsL[id_lim_u], 
					solver.limiter_MLP(uL,uL_max,uL_min, 
						dudxL*xyzLF[0] + dudyL*xyzLF[1] + dudzL*xyzLF[2], eta));
				cellsR[id_lim_u] = min(cellsR[id_lim_u], 
					solver.limiter_MLP(uR,uR_max,uR_min, 
						dudxR*xyzRF[0] + dudyR*xyzRF[1] + dudzR*xyzRF[2], eta));
						
				cellsL[id_lim_v] = min(cellsL[id_lim_v], 
					solver.limiter_MLP(vL,vL_max,vL_min, 
						dvdxL*xyzLF[0] + dvdyL*xyzLF[1] + dvdzL*xyzLF[2], eta));
				cellsR[id_lim_v] = min(cellsR[id_lim_v], 
					solver.limiter_MLP(vR,vR_max,vR_min, 
						dvdxR*xyzRF[0] + dvdyR*xyzRF[1] + dvdzR*xyzRF[2], eta));
						
				cellsL[id_lim_w] = min(cellsL[id_lim_w], 
					solver.limiter_MLP(wL,wL_max,wL_min, 
						dwdxL*xyzLF[0] + dwdyL*xyzLF[1] + dwdzL*xyzLF[2], eta));
				cellsR[id_lim_w] = min(cellsR[id_lim_w], 
					solver.limiter_MLP(wR,wR_max,wR_min, 
						dwdxR*xyzRF[0] + dwdyR*xyzRF[1] + dwdzR*xyzRF[2], eta));
						
				
		
				for(int ii=0; ii<2; ++ii){
				
					double phiL2[5+nSp], phiL1[5+nSp], phiR1[5+nSp];
					double tmp_limGrad = 1.0;
					if(ii==0){
						phiL1[0] = pL; phiR1[0] = pR;
						phiL2[0] = 0.0;
						phiL2[0] -= 2.0*dpdxL*xyzLR[0]; 
						phiL2[0] -= 2.0*dpdyL*xyzLR[1];
						phiL2[0] -= 2.0*dpdzL*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[0],pL_max,pL_min,phiL2[0], eta);
						phiL2[0] = phiR1[0] + tmp_limGrad*phiL2[0];
						
						phiL1[1] = uL; phiR1[1] = uR;
						phiL2[1] = 0.0;
						phiL2[1] -= 2.0*dudxL*xyzLR[0];
						phiL2[1] -= 2.0*dudyL*xyzLR[1];
						phiL2[1] -= 2.0*dudzL*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[1],uL_max,uL_min,phiL2[1], eta);
						phiL2[1] = phiR1[1] + tmp_limGrad*phiL2[1];
						
						phiL1[2] = vL; phiR1[2] = vR;
						phiL2[2] = 0.0;
						phiL2[2] -= 2.0*dvdxL*xyzLR[0];
						phiL2[2] -= 2.0*dvdyL*xyzLR[1];
						phiL2[2] -= 2.0*dvdzL*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[2],vL_max,vL_min,phiL2[2], eta);
						phiL2[2] = phiR1[2] + tmp_limGrad*phiL2[2];
						
						phiL1[3] = wL; phiR1[3] = wR;
						phiL2[3] = 0.0;
						phiL2[3] -= 2.0*dwdxL*xyzLR[0];
						phiL2[3] -= 2.0*dwdyL*xyzLR[1];
						phiL2[3] -= 2.0*dwdzL*xyzLR[2];
						tmp_limGrad = solver.limiter_MLP(phiL1[3],wL_max,wL_min,phiL2[3], eta);
						phiL2[3] = phiR1[3] + tmp_limGrad*phiL2[3];
						
						phiL1[4] = TL; phiR1[4] = TR;
						phiL2[4] = 0.0;
						phiL2[4] -= 2.0*dTdxL*xyzLR[0];
						phiL2[4] -= 2.0*dTdyL*xyzLR[1];
						phiL2[4] -= 2.0*dTdzL*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[4],TL_max,TL_min,phiL2[4], eta);
						phiL2[4] = phiR1[4] + tmp_limGrad*phiL2[4];
						
						for(int i=0; i<nSp-1; ++i){
							int tmp_i = 5+i;
							phiL1[tmp_i] = YL[i]; phiR1[tmp_i] = YR[i];
							phiL2[tmp_i] = 0.0;
							phiL2[tmp_i] -= 2.0*dYdxL[i]*xyzLR[0];
							phiL2[tmp_i] -= 2.0*dYdyL[i]*xyzLR[1];
							phiL2[tmp_i] -= 2.0*dYdzL[i]*xyzLR[2];
							// tmp_limGrad = solver.limiter_MLP(phiL1[tmp_i],YL_max[i],YL_min[i],phiL2[tmp_i], eta);
							phiL2[tmp_i] = phiR1[tmp_i] + tmp_limGrad*phiL2[tmp_i];
							
							phiL2[tmp_i] = max(0.0,min(1.0,phiL2[tmp_i]));
						}
						
						// // skewness
						// phiL1[0] += limGrad_phiL[0]*delphiL[0];
						// phiR1[0] += limGrad_phiR[0]*delphiR[0];
						// phiL2[0] += 2.0*limGrad_phiL[0]*delphiL[0];
						// phiL1[1] += limGrad_phiL[1]*delphiL[1];
						// phiR1[1] += limGrad_phiR[1]*delphiR[1];
						// phiL2[1] += 2.0*limGrad_phiL[1]*delphiL[1];
						// phiL1[2] += limGrad_phiL[2]*delphiL[2];
						// phiR1[2] += limGrad_phiR[2]*delphiR[2];
						// phiL2[2] += 2.0*limGrad_phiL[2]*delphiL[2];
						// phiL1[3] += limGrad_phiL[3]*delphiL[3];
						// phiR1[3] += limGrad_phiR[3]*delphiR[3];
						// phiL2[3] += 2.0*limGrad_phiL[3]*delphiL[3];
						// phiL1[4] += limGrad_phiL[4]*delphiL[4];
						// phiR1[4] += limGrad_phiR[4]*delphiR[4];
						// phiL2[4] += 2.0*limGrad_phiL[4]*delphiL[4];
						// for(int i=0; i<nSp-1; ++i){
							// int tmp_i = 5+i;
							// phiL1[tmp_i] += limGrad_phiL[tmp_i]*delphiL[tmp_i];
							// phiR1[tmp_i] += limGrad_phiR[tmp_i]*delphiR[tmp_i];
							// phiL2[tmp_i] += 2.0*limGrad_phiL[tmp_i]*delphiL[tmp_i];
							
							// phiL2[tmp_i] = max(0.0,min(1.0,phiL2[tmp_i]));
							// phiL1[tmp_i] = max(0.0,min(1.0,phiL1[tmp_i]));
							// phiR1[tmp_i] = max(0.0,min(1.0,phiR1[tmp_i]));
						// }
						
					}
					else{
						phiL1[0] = pR; phiR1[0] = pL;
						phiL2[0] = 0.0;
						phiL2[0] += 2.0*dpdxR*xyzLR[0];
						phiL2[0] += 2.0*dpdyR*xyzLR[1];
						phiL2[0] += 2.0*dpdzR*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[0],pR_max,pR_min,phiL2[0], eta);
						phiL2[0] = phiR1[0] + tmp_limGrad*phiL2[0];
						
						phiL1[1] = uR; phiR1[1] = uL;
						phiL2[1] = 0.0;
						phiL2[1] += 2.0*dudxR*xyzLR[0];
						phiL2[1] += 2.0*dudyR*xyzLR[1];
						phiL2[1] += 2.0*dudzR*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[1],uR_max,uR_min,phiL2[1], eta);
						phiL2[1] = phiR1[1] + tmp_limGrad*phiL2[1];
						
						phiL1[2] = vR; phiR1[2] = vL;
						phiL2[2] = 0.0;
						phiL2[2] += 2.0*dvdxR*xyzLR[0];
						phiL2[2] += 2.0*dvdyR*xyzLR[1];
						phiL2[2] += 2.0*dvdzR*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[2],vR_max,vR_min,phiL2[2], eta);
						phiL2[2] = phiR1[2] + tmp_limGrad*phiL2[2];
						
						phiL1[3] = wR; phiR1[3] = wL;
						phiL2[3] = 0.0;
						phiL2[3] += 2.0*dwdxR*xyzLR[0];
						phiL2[3] += 2.0*dwdyR*xyzLR[1];
						phiL2[3] += 2.0*dwdzR*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[3],wR_max,wR_min,phiL2[3], eta);
						phiL2[3] = phiR1[3] + tmp_limGrad*phiL2[3];
						
						phiL1[4] = TR; phiR1[4] = TL;
						phiL2[4] = 0.0;
						phiL2[4] += 2.0*dTdxR*xyzLR[0];
						phiL2[4] += 2.0*dTdyR*xyzLR[1];
						phiL2[4] += 2.0*dTdzR*xyzLR[2];
						// tmp_limGrad = solver.limiter_MLP(phiL1[4],TR_max,TR_min,phiL2[4], eta);
						phiL2[4] = phiR1[4] + tmp_limGrad*phiL2[4];
						
						for(int i=0; i<nSp-1; ++i){
							int tmp_i = 5+i;
							phiL1[tmp_i] = YR[i]; phiR1[tmp_i] = YL[i];
							phiL2[tmp_i] = 0.0;
							phiL2[tmp_i] += 2.0*dYdxR[i]*xyzLR[0];
							phiL2[tmp_i] += 2.0*dYdyR[i]*xyzLR[1];
							phiL2[tmp_i] += 2.0*dYdzR[i]*xyzLR[2];
							// tmp_limGrad = solver.limiter_MLP(phiL1[tmp_i],YR_max[i],YR_min[i],phiL2[tmp_i], eta);
							phiL2[tmp_i] = phiR1[tmp_i] + tmp_limGrad*phiL2[tmp_i];
							
							phiL2[tmp_i] = max(0.0,min(1.0,phiL2[tmp_i]));
						}
						
						// // skewness
						// phiL1[0] += limGrad_phiR[0]*delphiR[0];
						// phiR1[0] += limGrad_phiL[0]*delphiL[0];
						// phiL2[0] += 2.0*limGrad_phiR[0]*delphiR[0];
						// phiL1[2] += limGrad_phiR[2]*delphiR[2];
						// phiR1[2] += limGrad_phiL[2]*delphiL[2];
						// phiL2[2] += 2.0*limGrad_phiR[2]*delphiR[2];
						// phiL1[3] += limGrad_phiR[3]*delphiR[3];
						// phiR1[3] += limGrad_phiL[3]*delphiL[3];
						// phiL2[3] += 2.0*limGrad_phiR[3]*delphiR[3];
						// phiL1[4] += limGrad_phiR[4]*delphiR[4];
						// phiR1[4] += limGrad_phiL[4]*delphiL[4];
						// phiL2[4] += 2.0*limGrad_phiR[4]*delphiR[4];
						// for(int i=0; i<nSp-1; ++i){;
							// int tmp_i = 5+i;
							// phiL1[tmp_i] += limGrad_phiR[tmp_i]*delphiR[tmp_i];
							// phiR1[tmp_i] += limGrad_phiL[tmp_i]*delphiL[tmp_i];
							// phiL2[tmp_i] += 2.0*limGrad_phiR[tmp_i]*delphiR[tmp_i];
							
							// phiL2[tmp_i] = max(0.0,min(1.0,phiL2[tmp_i]));
							// phiL1[tmp_i] = max(0.0,min(1.0,phiL1[tmp_i]));
							// phiR1[tmp_i] = max(0.0,min(1.0,phiR1[tmp_i]));
						// }
						
						
					}
					
					// cout << NVD_scheme.size() << endl;
					for(int i=0; i<5; ++i) phiL1[i] = (solver.NVD.*NVD_scheme_i[i])(phiL2[i],phiL1[i],phiR1[i],0.0,0.0);
					
					{
						// double dt_tmp_L = dx/(sqrt(uL*uL+vL*vL+wL*wL)+cL);
						// double dt_tmp_R = dx/(sqrt(uR*uR+vR*vR+wR*wR)+cR);
						// double dt_tmp = dx/(sqrt(phiL1[1]*phiL1[1]+phiL1[2]*phiL1[2]+phiL1[3]*phiL1[3])+cL);
						// if(ii==1){
							// dt_tmp = dx/(sqrt(phiL1[1]*phiL1[1]+phiL1[2]*phiL1[2]+phiL1[3]*phiL1[3])+cR);
						// }		
						// double max_u = max(uL,uR);
						// double max_v = max(vL,vR);
						// double max_w = max(wL,wR);
						// double max_c = max(cL,cR);
						// double dt_tmp = dx/(sqrt(max_u*max_u+max_v*max_v+max_w*max_w)+max_c);
						double max_u = 0.5*(uL+uR);
						double max_v = 0.5*(vL+vR);
						double max_w = 0.5*(wL+wR);
						double max_c = 0.5*(cL+cR);
						double dt_tmp = dx/(sqrt(max_u*max_u+max_v*max_v+max_w*max_w)+1.e-200);
						double coDD = 1.0*dt/dt_tmp;
						
						for(int i=0; i<nSp-1; ++i){
							double mfLR[3];
							mfLR[0] = 0.5*dYdxL[i]+0.5*dYdxR[i];
							mfLR[1] = 0.5*dYdyL[i]+0.5*dYdyR[i];
							mfLR[2] = 0.5*dYdzL[i]+0.5*dYdzR[i];
							double magMfLR = 0.0;
							for(int j=0; j<3; ++j) magMfLR += mfLR[j]*mfLR[j];
							magMfLR = sqrt(magMfLR);
							for(int j=0; j<3; ++j) mfLR[j] = mfLR[j]/(magMfLR+1.e-200);
							double cosTheta = 0.0;
							cosTheta += mfLR[0]*faces[id_nx];
							cosTheta += mfLR[1]*faces[id_ny];
							cosTheta += mfLR[2]*faces[id_nz];
							cosTheta = abs(cosTheta);
							double gamF = min(cosTheta*cosTheta*cosTheta*cosTheta,1.0);
							
							phiL1[5+i] = (solver.NVD.*NVD_scheme_i[5+i])(phiL2[5+i],phiL1[5+i],phiR1[5+i],coDD,gamF);
						}
					}
					
					if(ii==0){
						faces[id_pL] = phiL1[0];
						faces[id_uL] = phiL1[1];
						faces[id_vL] = phiL1[2];
						faces[id_wL] = phiL1[3]; 
						faces[id_TL] = phiL1[4];
						for(int i=0; i<nSp-1; ++i){
							faces[id_YL[i]] = phiL1[5+i];
						}
					}
					else{
						faces[id_pR] = phiL1[0];
						faces[id_uR] = phiL1[1]; 
						faces[id_vR] = phiL1[2];
						faces[id_wR] = phiL1[3];
						faces[id_TR] = phiL1[4];
						for(int i=0; i<nSp-1; ++i){
							faces[id_YR[i]] = phiL1[5+i];
						}
					}
						
				
				}
				
				
				return 0;
			}
		);
	}
	// cout << "BBBBBBBB" << endl;
	
}


