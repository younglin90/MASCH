
#include "../../../../others/solvers.h"


void MASCH_Solver::setTimeStepFunctionsUDF(
MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	// for(int i=0; i<3; ++i)
	{
		int id_u = controls.getId_cellVar("x-velocity");
		int id_v = controls.getId_cellVar("y-velocity");
		int id_w = controls.getId_cellVar("z-velocity");
		int id_maxA = controls.getId_cellVar("maximum area");
		int id_minA = controls.getId_cellVar("minimum area");
		int id_vol = controls.getId_cellVar("volume");
		int id_mu = controls.getId_cellVar("viscosity");
		int id_rho = controls.getId_cellVar("density");
		int id_dt = controls.getId_fieldVar("time-step");
		double maxCFL = stod(controls.controlDictMap["maxCFL"]);
		double maxdt = stod(controls.controlDictMap["maxTimeStep"]);
		double dt = controls.timeStep;
		
		// 곡률관련
		vector<double> surf_sigma;
		for(int i=0; i<controls.spName.size(); ++i){
			string name = controls.spName[i];
			string type = controls.thermophysicalProperties[name+".transport.sigma.type"];
			if(type=="constant"){
				double value = stod(controls.thermophysicalProperties[name+".transport.sigma.value"]);
				if(value<1.e-200) continue;
				surf_sigma.push_back(value);
			}
		}
		int nCurv = surf_sigma.size();
	
	
		
		solver.calcTempStepCell.push_back(
		[id_dt,dt,id_vol,id_u,id_v,id_w,maxCFL,maxdt,id_c,id_mu,id_maxA,id_rho,
		id_minA,surf_sigma,nCurv](
		double* cells, double* fields) ->int {
			double maxA = cells[id_maxA];
			double minA = cells[id_minA];
			double vol = cells[id_vol];
			double u = cells[id_u]; double v = cells[id_v]; double w = cells[id_w];
			double char_speed = sqrt(u*u+v*v+w*w);
			double convTime = maxCFL*vol/maxA/(char_speed + 1.e-200);
			double mu = cells[id_mu];
			double rho = cells[id_rho];
			double viscTime = maxCFL* minA / ( mu + 1.e-200 ) / 6.0 * rho;
			double surfTime = 1.e8;
			for(int i=0; i<nCurv; ++i){
				double sigma = surf_sigma[i];
				surfTime = min(surfTime,maxCFL* sqrt(rho*vol/(2.0*3.141592*sigma)));
			}
			fields[id_dt] = min(fields[id_dt],min(convTime,min(surfTime,viscTime)));
			
			return 0;
		});
	}
	
}


