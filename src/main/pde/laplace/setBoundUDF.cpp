
#include "../../../others/solvers.h"

void MASCH_Solver::setBoundaryFunctions(
MASCH_Mesh& mesh, MASCH_Control& controls, MASCH_Variables& var){

	auto& solver = (*this);
	
	MASCH_Load load;
	
	using funct_type = 
	function<int(
	double time, double x, double y, double z, 
	double* cells, double* faces)>;
	
	for(int i=0; i<mesh.boundaries.size(); ++i){
		auto& boundary = mesh.boundaries[i];
		if(boundary.getType()!=MASCH_Face_Types::BOUNDARY) continue;
		
		string bcName = boundary.name;
		calcBoundFacePrimVal.push_back(vector<funct_type>());
		
		string name = "unknown";
		// string delta_name = "delta-unknown";
		{
			string type = controls.boundaryMap[name][bcName+".type"];
			
			
			int id_C = controls.getId_cellVar(name);
			int id_FL = controls.getId_faceVar("left "+name);
			int id_FR = controls.getId_faceVar("right "+name);
			
			// int id_dC = controls.getId_cellVar(delta_name);
			// int id_dFL = controls.getId_faceVar("left "+delta_name);
			// int id_dFR = controls.getId_faceVar("right "+delta_name);
		
			funct_type setFunct;
			
			if(type=="fixedValue"){
				double value = stod(controls.boundaryMap[name][bcName+".value"]);
				setFunct = (
				[value, id_C, id_FL, id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					faces[id_FL] = value;
					faces[id_FR] = value;
					return 0;
				});
			}
			else if(type=="zeroGradient"){
				setFunct = (
				[id_C, id_FL, id_FR](
				double time, double x, double y, double z, 
				double* cells, double* faces) ->int {
					double tmp_val = cells[id_C];
					faces[id_FL] = tmp_val;
					faces[id_FR] = tmp_val;
					return 0;
				});
			}
			else{
				cout << "#WARNING : not defiend " << name << " Boundary name = " << bcName << ", Boundary Type = " << type << endl;
			}
			
			calcBoundFacePrimVal.back().push_back(setFunct);
		}
		
		
		
	}
	
	// cout << calcBoundFacePrimVal.size() << endl;
	
}



