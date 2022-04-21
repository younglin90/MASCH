
#include "../../../others/solvers.h"

void MASCH_Solver::setMeanCellValuesFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
    
	int nSp = controls.spName.size();
    
    // MASCH_Load load;
	// controls.saveSMDValues = load.extractVector(controls.controlDictMap["saveSMDValues"]);
	// controls.saveMeanCellValues = load.extractVector(controls.controlDictMap["saveMeanCellValues"]);
    
    for(auto& item : controls.saveSMDValues){
        controls.meanSurfInp_cell_id.push_back(controls.getId_cellVar("volume-fraction-"+item));
        controls.meanSurfOut_cell_id.push_back(controls.getId_cellVar("fvm-surface-area-"+item));
        
        controls.meanVolInp_cell_id.push_back(controls.getId_cellVar("volume-fraction-"+item));
        controls.meanVolOut_cell_id.push_back(controls.getId_cellVar("fvm-volume-"+item));
        
        controls.meanInp_cell_id.push_back(controls.getId_cellVar("fvm-surface-area-"+item));
        controls.meanInp_cell_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-fvm-mean-surface-area-"+item));
        controls.meanOut_cell_id.push_back(controls.getId_cellVar("fvm-mean-surface-area-"+item));
        
        controls.meanInp_cell_id.push_back(controls.getId_cellVar("fvm-volume-"+item));
        controls.meanInp_cell_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-fvm-mean-volume-"+item));
        controls.meanOut_cell_id.push_back(controls.getId_cellVar("fvm-mean-volume-"+item));
    }
     
    
    
    for(auto& item : controls.saveMeanCellValues){
        controls.meanInp_cell_id.push_back(controls.getId_cellVar(item));
        controls.meanInp_cell_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-mean-"+item));
        controls.meanOut_cell_id.push_back(controls.getId_cellVar("mean-"+item));
    }
    
    
     
	// controls.meanInp_parcel_id.push_back(controls.getId_cellVar("parcel-surface-area"));
	// controls.meanInp_parcel_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-parcel-mean-surface-area"));
	// controls.meanOut_parcel_id.push_back(controls.getId_cellVar("parcel-mean-surface-area"));
    
	// controls.meanInp_parcel_id.push_back(controls.getId_parcelVar("diameter"));
	// controls.meanInp_parcel_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-parcel-mean-volume"));
	// controls.meanOut_parcel_id.push_back(controls.getId_cellVar("parcel-mean-volume"));
    
	
	
	
}


