
#include "./solvers.h"


void MASCH_Solver::setFunctions(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	/*
	cell : p u v w T Y1 Y2 rho c Ht drhodp drhodT drhodY1 drhodY2 dHtdp dHtdT dHtdY1 dH2dY2
	face : p u v w T Y1 Y2 rho c Ht
	field : dt
	
	*/
	
	solver.setSegEqUDF(mesh, controls);
	
	solver.setTimeStepFunctionsUDF(mesh, controls);
	
	solver.setOldVFunctionsUDF(mesh, controls);
	
	solver.setGradFunctionsUDF(mesh, controls);
	solver.setCurvatureFunctionsUDF(mesh, controls);
	
	solver.setHOReconFunctionsUDF(mesh, controls);
	// solver.setFValFunctionsUDF(mesh, controls);
	
	solver.setTermsCellLoopFunctionsUDF(mesh, controls);
	solver.setTermsFaceLoopFunctionsUDF(mesh, controls);
	// solver.setDiffFunctionsUDF(mesh, controls);
	// solver.setSourFunctionsUDF(mesh, controls);
	
	solver.setAddiFunctionsUDF(mesh, controls);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	solver.setUpdatePrimFunctionsUDF(mesh, controls);
	
	solver.setDPMFunctionsUDF(mesh, controls);
	
	solver.setMinMaxCellValuesFunctionsUDF(mesh, controls);
	
	solver.setMeanCellValuesFunctionsUDF(mesh, controls);
	
}





double MASCH_NVD::none(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	// return phiU;
	return 1.0;
}


double MASCH_NVD::minmod(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
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
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::vanLeer(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		double tildeCf = tildeCd + (tildeCd*(1.0-tildeCd))
						/(2.0-4.0*tildeCd+4.0*tildeCd*tildeCd);
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::QUICK(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		double tildeCf = 0.375 + 0.75*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::boundedCD(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::OSHER(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.6666){
		double tildeCf = 1.5*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.6666 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::SMART(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.83333){
		double tildeCf = 0.75*tildeCd + 0.375;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.83333 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::modifiedSMART(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.166666){
		double tildeCf = 3.0*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.166666 && tildeCd<0.7){
		double tildeCf = 0.75*tildeCd + 0.375;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.7 && tildeCd<1.0){
		double tildeCf = 0.3333*tildeCd + 0.6666;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::STOIC(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.83333){
		double tildeCf = 0.75*tildeCd + 0.375;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.83333 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::modifiedSTOIC(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.2){
		double tildeCf = 3.0*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.2 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.7){
		double tildeCf = 0.75*tildeCd + 0.375;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.7 && tildeCd<1.0){
		double tildeCf = 0.3333*tildeCd + 0.6666;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::MUSCL(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.25){
		double tildeCf = 2.0*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.25 && tildeCd<0.75){
		double tildeCf = tildeCd + 0.25;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.75 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::SUPERBEE(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.6666){
		double tildeCf = 1.5*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.6666 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	
	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

double MASCH_NVD::modifiedSUPERBEE(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.3333){
		double tildeCf = 2.0*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.3333 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.6666){
		double tildeCf = 1.5*tildeCd;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.6666 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}

	// return ( gamma_f * phiD + (1.0-gamma_f) * phiU );
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}

// ---------- lowwers ----------------
//
// sharp-interface schemes
// 
// diifuse 잘되는 순서대로 배열함
// STACS>HRIC>MIG>EB>FBICS>SAISH>CICSAM=MSTACS

double MASCH_NVD::STACS(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	// compressive differencing scheme (CDS)
	// SUPERBEE
	double gamma_f0 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.3333){
		double tildeCf = 2.0*tildeCd;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.3333 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.6666){
		double tildeCf = 1.5*tildeCd;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.6666 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfCompressive = phiU + gamma_f0 * (phiD - phiU);
	
	
	// high-resolution
	// STOIC
	double gamma_f1 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.2){
		double tildeCf = 3.0*tildeCd;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.2 && tildeCd<0.5){
		double tildeCf = 0.5*tildeCd + 0.5;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.83333){
		double tildeCf = 0.75*tildeCd + 0.375;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.83333 && tildeCd<1.0){
		double tildeCf = 1.0;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfDiffusive = phiU + gamma_f1 * (phiD - phiU);
	
    // double phiF = wsd * phiU + (1.0-wsd) * ( gamF*vfCompressive + (1.0-gamF)*vfDiffusive );
	// return phiF;
    
    return wsd + (1.0-wsd) * (gamF*(1.0-gamma_f0) + (1.0-gamF)*(1.0-gamma_f1));


}

double MASCH_NVD::HRIC(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double tildeCf = tildeCd;
    double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5){
		tildeCf = 2.0*tildeCd;
        double tildeStarF = gamF*tildeCf + (1.0-gamF)*tildeCd;
        gamma_f = (tildeStarF - tildeCd) / (1.0 - tildeCd);
        
        if(coDD>=0.3 && coDD<0.7){
            double tildeCf = tildeCd+(tildeStarF-tildeCd)*(0.7-coDD)/(0.7-0.3);
            gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
        }
        else if(coDD>=0.7 && coDD<1.0){
            double tildeCf = tildeCd;
            gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
        }
	}
	else if(tildeCd>=0.5 && tildeCd<1.0){
		tildeCf = 1.0;
        double tildeStarF = gamF*tildeCf + (1.0-gamF)*tildeCd;
        gamma_f = (tildeStarF - tildeCd) / (1.0 - tildeCd);
        
        if(coDD>=0.3 && coDD<0.7){
            double tildeCf = tildeCd+(tildeStarF-tildeCd)*(0.7-coDD)/(0.7-0.3);
            gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
        }
        else if(coDD>=0.7 && coDD<1.0){
            double tildeCf = tildeCd;
            gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
        }
	}
    
    // double phiF = wsd * phiU + (1.0-wsd) * ( phiU + gamma_f * (phiD - phiU) );
	// return phiF;

    return wsd + (1.0-wsd) * (1.0-gamma_f);

}


double MASCH_NVD::MIG(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	double tildeCf = tildeCd;
    double gamma_f = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5){
		tildeCf = -2.0*tildeCd*tildeCd + 3.0*tildeCd;
        double tildeStarF = gamF*tildeCf + (1.0-gamF)*tildeCd;
        gamma_f = (tildeStarF - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<1.0){
		tildeCf = 1.0;
        double tildeStarF = gamF*tildeCf + (1.0-gamF)*tildeCd;
        gamma_f = (tildeStarF - tildeCd) / (1.0 - tildeCd);
	}
    
    // double phiF = wsd * phiU + (1.0-wsd) * ( phiU + gamma_f * (phiD - phiU) );
	// return phiF;
    
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}


// TVD-based scheme
double MASCH_NVD::EB(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
    
	double rF = (phiU-phiUU)/(abs(phiD-phiU)+1.e-200)*( phiD>phiU ? 1.0 : -1.0 );
    
    double gamma_f = 0.5 * max(0.0, min( min(2.0/(1.0-coDD), 2.0*rF/coDD), 2.0+1.5*(rF-1.0) ) );
	
    // double phiF = wsd * phiU + (1.0-wsd) * ( phiU + gamma_f * (phiD - phiU) );
	// return phiF;
    
    return wsd + (1.0-wsd) * (1.0-gamma_f);

}



double MASCH_NVD::FBICS(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	// compressive differencing scheme (CDS)
	// BD-FBICS
	double gamma_f0 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.3333333333333){
		double tildeCf = min(1.0,3.0*tildeCd);
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.3333333333333 && tildeCd<1.0) {
		double tildeCf = 1.0;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfCompressive = phiU + gamma_f0 * (phiD - phiU);
	
	
	// high-resolution
	// HR-FBICS
	double gamma_f1 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.125) {
		double tildeCf = min(1.0,3.0*tildeCd);
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.125 && tildeCd<0.75) {
		double tildeCf = tildeCd + 0.25;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.75 && tildeCd<1.0) {
		double tildeCf = 1.0;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfDiffusive = phiU + gamma_f1 * (phiD - phiU);
	
    // double phiF = wsd * phiU + (1.0-wsd) * ( gamF*vfCompressive + (1.0-gamF)*vfDiffusive );
	// return phiF;

    return wsd + (1.0-wsd) * (gamF*(1.0-gamma_f0) + (1.0-gamF)*(1.0-gamma_f1));

}





double MASCH_NVD::SAISH(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	// compressive differencing scheme (CDS)
	// BD-SAISH
	double gamma_f0 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.25){
		double tildeCf = min(1.0,4.0*tildeCd);
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.25 && tildeCd<1.0) {
		double tildeCf = 1.0;
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfCompressive = phiU + gamma_f0 * (phiD - phiU);
	
	
	// high-resolution
	// HR-SAISH
	double gamma_f1 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.5) {
		double tildeCf = min(1.0,tildeCd*(2.0-tildeCd));
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.75) {
		double tildeCf = (tildeCd + 0.25);
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.75 && tildeCd<1.0) {
		double tildeCf = 1.0;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfDiffusive = phiU + gamma_f1 * (phiD - phiU);
	
    // double phiF = wsd * phiU + (1.0-wsd) * ( gamF*vfCompressive + (1.0-gamF)*vfDiffusive );
	// return phiF;
    
    return wsd + (1.0-wsd) * (gamF*(1.0-gamma_f0) + (1.0-gamF)*(1.0-gamma_f1));

}


double MASCH_NVD::CICSAM(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	// compressive differencing scheme (CDS)
	// Hyper-C
	double gamma_f0 = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		double tildeCf = min(1.0,tildeCd/coDD);
		gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfCompressive = phiU + gamma_f0 * (phiD - phiU);
	
	// high-resolution
	// ULTIMATE QUICKEST (UQ)
	double gamma_f1 = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		double tildeCf = min((8.0*coDD*tildeCd+(1.0-coDD)*(6.0*tildeCd+3.0))/8.0,gamma_f0);
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	double vfDiffusive = phiU + gamma_f1 * (phiD - phiU);
    // double phiF = wsd * phiU + (1.0-wsd) * ( gamF*vfCompressive + (1.0-gamF)*vfDiffusive );
	// return phiF;
    
    return wsd + (1.0-wsd) * (gamF*(1.0-gamma_f0) + (1.0-gamF)*(1.0-gamma_f1));

}




double MASCH_NVD::MSTACS(double phiUU, double phiU, double phiD, double coDD, double gamF, double wsd) {
	
	double tildeCd = (phiU-phiUU)/(abs(phiD-phiUU)+1.e-200)*( phiD>phiUU ? 1.0 : -1.0 );
	
	// compressive differencing scheme (CDS)
	// CDS-MSTACS
	double gamma_f0 = 0.0;
	if(tildeCd>=0.0 && tildeCd<1.0){
		if(coDD>0.0 && coDD<=0.33){
			double tildeCf = min(1.0,tildeCd/coDD);
			gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
		}
		else if(coDD>0.33 && coDD<=1.0){
			double tildeCf = min(1.0,3.0*tildeCd);
			gamma_f0 = (tildeCf - tildeCd) / (1.0 - tildeCd);
		}
	}
	double vfCompressive = phiU + gamma_f0 * (phiD - phiU);
	
	
	// high-resolution
	// HR-STOIC
	double gamma_f1 = 0.0;
	if(tildeCd>=0.0 && tildeCd<0.2) {
		double tildeCf = min(1.0,3.0*tildeCd);
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.2 && tildeCd<0.5) {
		double tildeCf = 0.5 + 0.5*tildeCd;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.5 && tildeCd<0.8333) {
		double tildeCf = 0.375 + 0.75*tildeCd;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd);
	}
	else if(tildeCd>=0.8333 && tildeCd<1.0) {
		double tildeCf = 1.0;
		gamma_f1 = (tildeCf - tildeCd) / (1.0 - tildeCd); 
	}
	double vfDiffusive = phiU + gamma_f1 * (phiD - phiU);
	
    // double phiF = wsd * phiU + (1.0-wsd) * ( gamF*vfCompressive + (1.0-gamF)*vfDiffusive );
	// return phiF;
    
    return wsd + (1.0-wsd) * (gamF*(1.0-gamma_f0) + (1.0-gamF)*(1.0-gamma_f1));

}








// ===================================
// unstructured based limiter

double MASCH_Solver::limiter_MLP(double phi, double phi_max, double phi_min, double Delta_minus, double eta) {
	
	if(Delta_minus>0.0){
		double Delta_plus = phi_max - phi;
		return 1.0/Delta_minus*
			((Delta_plus*Delta_plus+eta*eta)*
			Delta_minus+2.0*Delta_minus*Delta_minus*Delta_plus)/
			(Delta_plus*Delta_plus+2.0*Delta_minus*Delta_minus+
			Delta_minus*Delta_plus+eta*eta);
	}
	else if(Delta_minus<0.0){
		double Delta_plus = phi_min - phi;
		return 1.0/Delta_minus*
			((Delta_plus*Delta_plus+eta*eta)*
			Delta_minus+2.0*Delta_minus*Delta_minus*Delta_plus)/
			(Delta_plus*Delta_plus+2.0*Delta_minus*Delta_minus+
			Delta_minus*Delta_plus+eta*eta);
	}
	return 1.0;

}





void MASCH_Solver::updateBoundaryFaceValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var, int iSegEq){
	
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
	
	// cout << rho << " " << gam << " " << cv << " " << T << " " << P << " " << pinf << " " << bNASG << endl;
	
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


void MASCH_Solver::updateCellAddiValues_All(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	for(int iSegEq=0, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
		updateCellAddiValues(mesh, controls, var, iSegEq);
	}
}

void MASCH_Solver::updateBoundaryFacePrimValues_All(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	for(int iSegEq=0, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
		updateBoundaryFacePrimValues(mesh, controls, var, iSegEq);
	}
}

void MASCH_Solver::gradientTerms_All(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	if(solver.gradLSIds_cell_name.size()==0) return;
		
	vector<string> name = solver.gradLSIds_cell_name[0];
	vector<string> bcname = solver.gradLSIds_bcFace_name[0];
		
	for(int iSegEq=1, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
		for(auto& item : solver.gradLSIds_cell_name[iSegEq]){
			if(find(name.begin(),name.end(),item)==name.end()){
				name.push_back(item);
			}
		}
		for(auto& item : solver.gradLSIds_bcFace_name[iSegEq]){
			if(find(bcname.begin(),bcname.end(),item)==bcname.end()){
				bcname.push_back(item);
			}
		}
	}
		
	solver.calcGradient.leastSquare(mesh, controls, var, name, bcname);
	
		
}

void MASCH_Solver::curvatureTerms_All(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
	
	// cout << curvatureIds_cell_name.size() << endl;
	
	if(solver.curvatureIds_cell_name.size()==0) return;
		
	vector<string> name = solver.curvatureIds_cell_name[0];
	// vector<string> bcname = solver.curvatureLSIds_bcFace_name[0];
		
	for(int iSegEq=1, nSegEq=controls.nEq.size(); iSegEq<nSegEq; ++iSegEq){
		for(auto& item : solver.curvatureIds_cell_name[iSegEq]){
			if(find(name.begin(),name.end(),item)==name.end()){
				name.push_back(item);
			}
		}
		// for(auto& item : solver.gradLSIds_bcFace_name[iSegEq]){
			// if(find(bcname.begin(),bcname.end(),item)==bcname.end()){
				// bcname.push_back(item);
			// }
		// }
	}
		
	solver.calcCurvature(mesh, controls, var, name);
	
		
}

void MASCH_Solver::updateProcRightCellValues_All(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
		
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	
	if(size>1){
		int proc_size = var.procRightCells.size();
		int prim_size = var.cells.back().size();
		// vector<int> ids = controls.primIdsSegEq[iSegEq];
		
		vector<double> send_value;
		// send_value.reserve(proc_size*prim_size);
		for(auto& boundary : mesh.boundaries){
			if(boundary.getType()==MASCH_Face_Types::PROCESSOR){
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				for(int i=str; i<end; ++i){
					auto cellVar_i = cellVar[faces[i].iL].data();
					auto faceVar_i = faceVar[i].data();
					for(int j=0; j<prim_size; ++j){
						send_value.push_back(cellVar_i[j]);
					}
				}
			}
		}
		
		vector<int> tmp_counts(size,0);
		vector<int> tmp_displs(size+1,0);
		for(int ip=0; ip<size; ++ip){
			tmp_counts[ip] = mesh.countsSendProcFaces[ip]*prim_size;
		}
		for(int ip=0; ip<size; ++ip){
			tmp_displs[ip+1] = tmp_displs[ip] + tmp_counts[ip];
		}
		
		vector<double> recv_value(proc_size*prim_size);
		MPI_Alltoallv( send_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
						recv_value.data(), tmp_counts.data(), 
						tmp_displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		auto recv_value_ptr = recv_value.data();
		int iter=0;
		for(auto& cells : var.procRightCells){
			auto cellVar_i = cells.data();
			for(int j=0; j<prim_size; ++j){
				// cout << recv_value_ptr[iter] << endl;
				cellVar_i[j] = recv_value_ptr[iter++];
			}
			
			// cout << cellVar_i[12] << " " << cellVar_i[16] << endl;
			// cout << prim_size << " " << cellVar_i[0] << endl;
		}
	}
	
	
}





void MASCH_Solver::calcMeanValues(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){
	auto& solver = (*this);
		
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	auto cells = mesh.cells.data();
	auto faces = mesh.faces.data();
	auto cellVar = var.cells.data();
	auto faceVar = var.faces.data();
	auto fieldVar = var.fields.data();
	auto parcelVar = var.parcels.data();
    
    if(controls.saveSMDValues.size()==0) return;
    
    string name = controls.saveSMDValues[0];
    
    
	int id_vol = controls.getId_cellVar("volume");
	int id_nx = controls.getId_faceVar("x unit normal");
	int id_ny = controls.getId_faceVar("y unit normal");
	int id_nz = controls.getId_faceVar("z unit normal");
	int id_area = controls.getId_faceVar("area");
    
    vector<vector<double>> gradx, grady, gradz; 
    gradx.resize(controls.meanSurfOut_cell_id.size(),vector<double>(mesh.cells.size(),0.0));
    grady.resize(controls.meanSurfOut_cell_id.size(),vector<double>(mesh.cells.size(),0.0));
    gradz.resize(controls.meanSurfOut_cell_id.size(),vector<double>(mesh.cells.size(),0.0));
    for(int i=0, SIZE=mesh.faces.size(); i<SIZE; ++i){
        auto& face = mesh.faces[i];
        int iL = face.iL;
        int iR = face.iR;
        double nvec[3];
        nvec[0] = faceVar[i][id_nx]; nvec[1] = faceVar[i][id_ny]; nvec[2] = faceVar[i][id_nz];
        double area = faceVar[i][id_area];
        if(face.getType()==MASCH_Face_Types::INTERNAL){
            int tmp_iter = 0;
            for(auto& item : controls.meanSurfInp_cell_id){
                double flux = 0.5*(cellVar[iL][item] + cellVar[iR][item]);
                gradx[tmp_iter][iL] += flux*nvec[0]*area;
                grady[tmp_iter][iL] += flux*nvec[1]*area;
                gradz[tmp_iter][iL] += flux*nvec[2]*area;
                
                gradx[tmp_iter][iR] -= flux*nvec[0]*area;
                grady[tmp_iter][iR] -= flux*nvec[1]*area;
                gradz[tmp_iter][iR] -= flux*nvec[2]*area;
                
                ++tmp_iter;
            }
        }
        else{
            int tmp_iter = 0;
            for(auto& item : controls.meanVolInp_cell_id){
                double flux = cellVar[iL][item];
                gradx[tmp_iter][iL] += flux*nvec[0]*area;
                grady[tmp_iter][iL] += flux*nvec[1]*area;
                gradz[tmp_iter][iL] += flux*nvec[2]*area;
                
                ++tmp_iter;
            }
        }
    }
    for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
        auto& cell = mesh.cells[i];
        int tmp_iter = 0;
        for(auto& item : controls.meanSurfOut_cell_id){
            double mag = 0.0;
            mag += gradx[tmp_iter][i]*gradx[tmp_iter][i];
            mag += grady[tmp_iter][i]*grady[tmp_iter][i];
            mag += gradz[tmp_iter][i]*gradz[tmp_iter][i];
            // cellVar[i][item] = sqrt(mag)*cellVar[i][id_vol];
            cellVar[i][item] = sqrt(mag);
            ++tmp_iter;
        }
        tmp_iter = 0;
        for(auto& item : controls.meanVolOut_cell_id){
            int tmp_j = controls.meanVolInp_cell_id[tmp_iter];
            cellVar[i][item] = cellVar[i][tmp_j]*cellVar[i][id_vol];
            ++tmp_iter;
        }
    }
    
    
    
    int id_values_size = controls.meanInp_cell_id.size();
    int id_parcel_values_size = controls.meanInp_parcel_id.size();
    vector<int> id_values = controls.meanOut_cell_id;
    vector<int> id_parcel_values = controls.meanOut_parcel_id;
    vector<int> id_targets = controls.meanInp_cell_id;
    vector<int> id_parcel_targets = controls.meanInp_parcel_id;
    vector<double> totalTimes;
    // vector<double> parcel_totalTimes;
    for(auto& item : controls.meanInp_cell_totalTime_id){
        totalTimes.push_back(fieldVar[item]);
    }
    // for(auto& item : controls.meanInp_parcel_totalTime_id){
        // parcel_totalTimes.push_back(fieldVar[item]);
    // }
    
    double dt = fieldVar[controls.getId_fieldVar("time-step")];
    int id_dia = controls.getId_parcelVar("diameter");
    double& totalTime_area = fieldVar[controls.getId_fieldVar("total-time-of-parcel-mean-surface-area-"+name)];
    double& totalTime_vol = fieldVar[controls.getId_fieldVar("total-time-of-parcel-mean-volume-"+name)];
    int id_d_area = controls.getId_cellVar("parcel-mean-surface-area-"+name);
    int id_d_vol = controls.getId_cellVar("parcel-mean-volume-"+name);
     
    
	auto parcels_i = mesh.parcels.data();
    for(int i=0, SIZE=mesh.cells.size(); i<SIZE; ++i){
        auto& cell = mesh.cells[i];
        auto cellVars_i = cellVar[i].data();
        // cell
        for(int j=0; j<id_values_size; ++j){
            int tmp_id1 = id_values[j];
            int tmp_id2 = id_targets[j];
            double tmp_tTime = totalTimes[j];
            
            cellVars_i[tmp_id1] = 
                (totalTimes[j]*cellVars_i[tmp_id1] + dt*cellVars_i[tmp_id2]) / (tmp_tTime+dt);
        }
        // parcel 
        {
            double mean_area = 0.0;
            double mean_vol = 0.0;
            double tmp_size = 0.0;
			for(auto& iparcel : cell.iparcels){
				auto& parcel = parcels_i[iparcel];
				auto parcelVar_i = parcelVar[iparcel].data();
				if(parcel.getType() != MASCH_Parcel_Types::INSIDE) continue;
                double d_p = parcelVar_i[id_dia];
                double d_vol = 4.0/3.0*3.141592*0.125*d_p*d_p*d_p;
                double d_area = 3.141592*d_p*d_p; //3.141592*0.25*d_p*d_p;
                mean_area += d_area;
                mean_vol += d_vol;
                tmp_size += 1.0;
			}
            if(cell.iparcels.size()>0){
                mean_area /= tmp_size;
                mean_vol /= tmp_size;
            }
            cellVars_i[id_d_area] = 
                (totalTime_area*cellVars_i[id_d_area] + dt*mean_area) / (totalTime_area+dt);
                // (totalTime_area*cellVars_i[id_d_area]) / (totalTime_area+dt + 1.e-200);
            cellVars_i[id_d_vol] = 
                (totalTime_vol*cellVars_i[id_d_vol] + dt*mean_vol) / (totalTime_vol+dt);
                // (dt*mean_vol) / (totalTime_vol+dt + 1.e-20);
                
        }
    }
    
    
    for(auto& item : totalTimes){
        item += dt;
    }
    totalTime_area += dt;
    totalTime_vol += dt;
    // for(auto& item : parcel_totalTimes){
        // item += dt;
    // }

    
    
    
    
    
    
    
}