
#include "../../../others/solvers.h"

void MASCH_Control::setVariablesUDF(vector<string>& species){
	
	// 셀 값 정의
	(*this).setVarible({"cell"},"pressure","p","","primitive","scalar");
	(*this).setVarible({"cell"},"velocity","U","","primitive","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{"primitive","primitive","primitive"});
	(*this).setVarible({"cell"},"temperature","T","","primitive","scalar");
	vector<string> sub_names, sub_rols, sub_abb;
	for(int i=0; i<species.size()-1; ++i){
		sub_names.push_back("mass-fraction-"+species[i]);
		sub_rols.push_back("primitive");
		sub_abb.push_back("");
	}
	// sub_names.push_back("mass-fraction-"+species.back());
	// sub_rols.push_back("");
	// sub_abb.push_back("");
	(*this).setVarible({"cell"},"mass-fraction","Y","","primitive","vector",
						sub_names,sub_abb,sub_rols);
	(*this).setVarible({"cell"},"mass-fraction-"+species.back(),"","","","scalar");

	// 추가적 셀값
	(*this).setVarible({"cell"},"partial-density-pressure","","","","scalar");
	(*this).setVarible({"cell"},"partial-density-temperature","","","","scalar");
	(*this).setVarible({"cell"},"partial-total-enthalpy-pressure","","","","scalar");
	(*this).setVarible({"cell"},"partial-total-enthalpy-temperature","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		(*this).setVarible({"cell"},"partial-density-mass-fraction-"+species[i],"","","","scalar");
		(*this).setVarible({"cell"},"partial-total-enthalpy-mass-fraction-"+species[i],"","","","scalar");
	}
	(*this).setVarible({"cell"},"density","","","","scalar");
	(*this).setVarible({"cell"},"speed-of-sound","","","","scalar");
	(*this).setVarible({"cell"},"total-enthalpy","","","","scalar");
	(*this).setVarible({"cell"},"viscosity","","","","scalar");
	(*this).setVarible({"cell"},"curvature","","","","scalar");
	

	// 그레디언트 값
	(*this).setVarible({"cell"},"gradient pressure","","","","vector3",
						{"x-gradient pressure","y-gradient pressure","z-gradient pressure"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"gradient x-velocity","","","","vector3",
						{"x-gradient x-velocity","y-gradient x-velocity","z-gradient x-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"gradient y-velocity","","","","vector3",
						{"x-gradient y-velocity","y-gradient y-velocity","z-gradient y-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"gradient z-velocity","","","","vector3",
						{"x-gradient z-velocity","y-gradient z-velocity","z-gradient z-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"gradient temperature","","","","vector3",
						{"x-gradient temperature","y-gradient temperature","z-gradient temperature"},{"","",""},{"","",""});
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+species[i]);
		(*this).setVarible({"cell"},"gradient "+tmp_name,"","","","vector3",
							{"x-gradient "+tmp_name,"y-gradient "+tmp_name,"z-gradient "+tmp_name},{"","",""},{"","",""});
	}
	
	
	// 올드 값
	(*this).setVarible({"cell"},"old pressure","","","","scalar");
	(*this).setVarible({"cell"},"old velocity","","","","vector3",
						{"old x-velocity","old y-velocity","old z-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"old density","","","","scalar");
	(*this).setVarible({"cell"},"old total-enthalpy","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("mass-fraction"+species[i]);
		(*this).setVarible({"cell"},"old "+tmp_name,"","","","scalar");
	}
	
	// 페이스 값 정의
	(*this).setVarible({"face"},"left pressure","","","","scalar");
	(*this).setVarible({"face"},"right pressure","","","","scalar");
	(*this).setVarible({"face"},"contravariant-velocity","","","","scalar");
	(*this).setVarible({"face"},"left x-velocity","","","","scalar");
	(*this).setVarible({"face"},"right x-velocity","","","","scalar");
	(*this).setVarible({"face"},"left y-velocity","","","","scalar");
	(*this).setVarible({"face"},"right y-velocity","","","","scalar");
	(*this).setVarible({"face"},"left z-velocity","","","","scalar");
	(*this).setVarible({"face"},"right z-velocity","","","","scalar");
	(*this).setVarible({"face"},"left temperature","","","","scalar");
	(*this).setVarible({"face"},"right temperature","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		(*this).setVarible({"face"},"left mass-fraction-"+species[i],"","","","scalar");
		(*this).setVarible({"face"},"right mass-fraction-"+species[i],"","","","scalar");
	}
	(*this).setVarible({"face"},"left density","","","","scalar");
	(*this).setVarible({"face"},"right density","","","","scalar");
	(*this).setVarible({"face"},"left speed-of-sound","","","","scalar");
	(*this).setVarible({"face"},"right speed-of-sound","","","","scalar");
	(*this).setVarible({"face"},"left total-enthalpy","","","","scalar");
	(*this).setVarible({"face"},"right total-enthalpy","","","","scalar");
	// (*this).setVarible({"face"},"viscosity","","","","scalar");
	
}

