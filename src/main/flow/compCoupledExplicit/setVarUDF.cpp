
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
		// sub_names.push_back("mass-fraction-"+species[i]);
		sub_names.push_back(species[i]);
		sub_rols.push_back("primitive");
		sub_abb.push_back("");
	}
	// sub_names.push_back("mass-fraction-"+species.back());
	// sub_rols.push_back("");
	// sub_abb.push_back("");
	(*this).setVarible({"cell"},"mass-fraction","Y","","primitive","vector",
						sub_names,sub_abb,sub_rols);
	(*this).setVarible({"cell"},"mass-fraction-"+species.back(),"","","","scalar");
	for(int i=0; i<species.size(); ++i){
		(*this).setVarible({"cell"},"volume-fraction-"+species[i],"","","","scalar");
	}
	

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
	
	// 최대 최소값
	(*this).setVarible({"cell"},"maximum pressure","","","","scalar");
	(*this).setVarible({"cell"},"minimum pressure","","","","scalar");
	
	(*this).setVarible({"cell"},"maximum x-velocity","","","","scalar");
	(*this).setVarible({"cell"},"minimum x-velocity","","","","scalar");
	
	(*this).setVarible({"cell"},"maximum y-velocity","","","","scalar");
	(*this).setVarible({"cell"},"minimum y-velocity","","","","scalar");
	
	(*this).setVarible({"cell"},"maximum z-velocity","","","","scalar");
	(*this).setVarible({"cell"},"minimum z-velocity","","","","scalar");
	
	(*this).setVarible({"cell"},"maximum temperature","","","","scalar");
	(*this).setVarible({"cell"},"minimum temperature","","","","scalar");
	
	for(int i=0; i<species.size()-1; ++i){
		(*this).setVarible({"cell"},"maximum mass-fraction-"+species[i],"","","","scalar");
		(*this).setVarible({"cell"},"minimum mass-fraction-"+species[i],"","","","scalar");
	}

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
	(*this).setVarible({"cell"},"gradient density","","","","vector3",
						{"x-gradient density","y-gradient density","z-gradient density"},{"","",""},{"","",""});
	
	
	// 올드 값
	(*this).setVarible({"cell"},"old pressure","","","","scalar");
	(*this).setVarible({"cell"},"old velocity","","","","vector3",
						{"old x-velocity","old y-velocity","old z-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"old density","","","","scalar");
	(*this).setVarible({"cell"},"old total-enthalpy","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+species[i]);
		(*this).setVarible({"cell"},"old "+tmp_name,"","","","scalar");
	}
	
	// 페이스 값 정의
	(*this).setVarible({"face"},"left pressure","","","","scalar");
	(*this).setVarible({"face"},"right pressure","","","","scalar");
	// (*this).setVarible({"face"},"pressure-org","","","","scalar");
	// (*this).setVarible({"face"},"contravariant-velocity","","","","scalar");
	(*this).setVarible({"face"},"left velocity","","","","vector3",
						{"left x-velocity","left y-velocity","left z-velocity"},{"u","v","w"},{"","",""});
	(*this).setVarible({"face"},"right velocity","","","","vector3",
						{"right x-velocity","right y-velocity","right z-velocity"},{"u","v","w"},{"","",""});
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
	(*this).setVarible({"face"},"left viscosity","","","","scalar");
	(*this).setVarible({"face"},"right viscosity","","","","scalar");
	
	
	// 곡률 관련
	for(int i=0; i<species.size(); ++i){
		string type = (*this).thermophysicalProperties[species[i]+".transport.sigma.type"];
		if(type=="constant"){
			double value = stod((*this).thermophysicalProperties[species[i]+".transport.sigma.value"]);
			if(value<1.e-200) continue;
			(*this).setVarible({"cell"},"curvature-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"level-set-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"x-gradient level-set-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"y-gradient level-set-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"z-gradient level-set-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"x-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"y-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"z-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"x-gradient x-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"x-gradient y-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"x-gradient z-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"y-gradient x-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"y-gradient y-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"y-gradient z-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"z-gradient x-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"z-gradient y-unit-normal-"+species[i],"","","","scalar");
			(*this).setVarible({"cell"},"z-gradient z-unit-normal-"+species[i],"","","","scalar");
		}
	}
	
	// 리미터 관련
	(*this).setVarible({"cell"},"limiter-unstructured pressure","","","","scalar");
	(*this).setVarible({"cell"},"limiter-unstructured x-velocity","","","","scalar");
	(*this).setVarible({"cell"},"limiter-unstructured y-velocity","","","","scalar");
	(*this).setVarible({"cell"},"limiter-unstructured z-velocity","","","","scalar");
	
	
	
	
	
	
	// parcels
	(*this).setVarible({"parcel"},"x-location","","","","scalar");
	(*this).setVarible({"parcel"},"y-location","","","","scalar");
	(*this).setVarible({"parcel"},"z-location","","","","scalar");
	(*this).setVarible({"parcel"},"diameter","","","","scalar");
	(*this).setVarible({"parcel"},"density","","","","scalar");
	(*this).setVarible({"parcel"},"velocity","U","","","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{"","",""});
	(*this).setVarible({"parcel"},"temperature","","","","scalar");
	(*this).setVarible({"parcel"},"number-of-parcel","","","","scalar");
	(*this).setVarible({"field"},"parcel-injection-accum-time","","","","scalar");
	(*this).setVarible({"field"},"time-step-parcels","","","","scalar");
	
	
	
	
}

