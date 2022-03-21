
#include "../../../others/solvers.h"

void MASCH_Solver::setDPMFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	solver.calcDPM_parcelLoop.push_back(
		[](
		double* cells, double* fields, double* parcels, double* fluxB) ->int {
			return 0;
			
		}
	);
	
	
	
	
	
	
}



// void MASCH_DPM::calcMechanicalForceModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }
	
// void MASCH_DPM::calcBreakupModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }
	
// void MASCH_DPM::calcHeatAndMassTransferModel(
// MASCH_Mesh& mesh, 
// MASCH_Control& controls,
// MASCH_Variables& var){
	
	
	
// }