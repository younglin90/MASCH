
#include "../../../others/solvers.h"

void MASCH_Solver::setMinMaxCellValuesFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
    
    
    // 셀 min max
	minmaxInp_cell_name.resize(1);
	maxOut_cell_name.resize(1);
	minOut_cell_name.resize(1);
	
	minmaxInp_cell_name[0].push_back("pressure");
	maxOut_cell_name[0].push_back("maximum pressure");
	minOut_cell_name[0].push_back("minimum pressure");
	
	minmaxInp_cell_name[0].push_back("x-velocity");
	maxOut_cell_name[0].push_back("maximum x-velocity");
	minOut_cell_name[0].push_back("minimum x-velocity");
	
	minmaxInp_cell_name[0].push_back("y-velocity");
	maxOut_cell_name[0].push_back("maximum y-velocity");
	minOut_cell_name[0].push_back("minimum y-velocity");
	
	minmaxInp_cell_name[0].push_back("z-velocity");
	maxOut_cell_name[0].push_back("maximum z-velocity");
	minOut_cell_name[0].push_back("minimum z-velocity");
	
	minmaxInp_cell_name[0].push_back("temperature");
	maxOut_cell_name[0].push_back("maximum temperature");
	minOut_cell_name[0].push_back("minimum temperature");
	
	for(int i=0; i<controls.spName.size()-1; ++i){
		minmaxInp_cell_name[0].push_back("mass-fraction-"+controls.spName[i]);
		maxOut_cell_name[0].push_back("maximum mass-fraction-"+controls.spName[i]);
		minOut_cell_name[0].push_back("minimum mass-fraction-"+controls.spName[i]);
	}
    
    
    
    
    // 포인트 min max
	minmaxInp_point_name.resize(1);
	maxOut_point_name.resize(1);
	minOut_point_name.resize(1);
    
	minmaxInp_point_name[0].push_back("x-velocity");
	maxOut_point_name[0].push_back("maximum x-velocity");
	minOut_point_name[0].push_back("minimum x-velocity");
    
	minmaxInp_point_name[0].push_back("y-velocity");
	maxOut_point_name[0].push_back("maximum y-velocity");
	minOut_point_name[0].push_back("minimum y-velocity");
    
	minmaxInp_point_name[0].push_back("z-velocity");
	maxOut_point_name[0].push_back("maximum z-velocity");
	minOut_point_name[0].push_back("minimum z-velocity");
    
	for(int i=0; i<controls.spName.size()-1; ++i){
		minmaxInp_point_name[0].push_back("mass-fraction-"+controls.spName[i]);
		maxOut_point_name[0].push_back("maximum mass-fraction-"+controls.spName[i]);
		minOut_point_name[0].push_back("minimum mass-fraction-"+controls.spName[i]);
	}
    
    
	
}