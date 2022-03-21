
#include "../../../others/solvers.h"

void MASCH_Control::setVariablesUDF(vector<string>& species){
	
	
	vector<string> vari;
	vector<string> abb;
	vector<string> role;
	vector<string> shape;
	vector<string> unit;
	vector<vector<string>> sub_vari;
	vector<vector<string>> sub_abb;
	
	string prim = "primitive";
	string thermo = "thermo property";
	string scal = "scalar";
	string vec = "vector";
	string vec3 = "vector3";
	string cell = "cell";
	string face = "face";
	string point = "point";
	string field = "field";
	
	// 셀 값 정의
	(*this).setVarible({cell},"pressure","p","Pa",prim,scal);
	(*this).setVarible({cell},"velocity","U","m/s",prim,vec3,
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{prim,prim,prim});
	;
	(*this).setVarible({cell},"delta-pressure","dp","Pa",prim,scal);
	(*this).setVarible({cell},"density","rho","kg/m^3",thermo,scal);
	
	(*this).setVarible({cell},"gradient delta-pressure","dpdX","","gradient",vec3,
			{"x-gradient delta-pressure","y-gradient delta-pressure","z-gradient delta-pressure"},
			{"dpdx","dpdy","dpdz"},{"","",""});
	
	// 페이스 값 정의
	(*this).setVarible({face},"left pressure","p","Pa","value",scal);
	(*this).setVarible({face},"left delta-pressure","p","Pa","value",scal);
	(*this).setVarible({face},"left velocity","U","m/s","value",vec3,
						{"left x-velocity","left y-velocity","left z-velocity"},{"u","v","w"},{"","",""});
	(*this).setVarible({face},"left density","rho","kg/m^3",thermo,scal);
	
	(*this).setVarible({face},"right pressure","p","Pa","value",scal);
	(*this).setVarible({face},"right delta-pressure","p","Pa","value",scal);
	(*this).setVarible({face},"right velocity","U","m/s","value",vec3,
						{"right x-velocity","right y-velocity","right z-velocity"},{"u","v","w"},{"","",""});
	(*this).setVarible({face},"right density","rho","kg/m^3",thermo,scal);
	
	
	
	
}

