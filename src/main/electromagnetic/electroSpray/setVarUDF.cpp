
#include "../../../others/solvers.h"

void MASCH_Control::setVariablesUDF(vector<string>& species){
	
	// 셀 값 정의
	(*this).setVarible({"cell"},"pressure","p","","primitive","scalar");
	(*this).setVarible({"cell"},"velocity","U","","primitive","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{"primitive","primitive","primitive"});
	(*this).setVarible({"cell"},"charge-density","rhoe","","primitive","scalar");
	(*this).setVarible({"cell"},"electric-potential","phi","","primitive","scalar");
	vector<string> sub_names, sub_rols, sub_abb;
	for(int i=0; i<species.size()-1; ++i){
		sub_names.push_back(species[i]);
		sub_rols.push_back("primitive");
		sub_abb.push_back("");
	}
	(*this).setVarible({"cell"},"volume-fraction","alpha","","primitive","vector",
						sub_names,sub_abb,sub_rols);
	(*this).setVarible({"cell"},"volume-fraction-"+species.back(),"","","","scalar");
	for(int i=0; i<species.size(); ++i){
		(*this).setVarible({"cell"},"mass-fraction-"+species[i],"","","","scalar");
	}
	(*this).setVarible({"cell"},"electric-field","E","","","vector3",
						{"x-electric-field","y-electric-field","z-electric-field"},{"","",""},
						{"","",""});
	(*this).setVarible({"cell"},"electric-force","Fe","","","vector3",
						{"x-electric-force","y-electric-force","z-electric-force"},{"","",""},
						{"","",""});
	

	// 추가적 셀값
	(*this).setVarible({"cell"},"density","","","","scalar");
	(*this).setVarible({"cell"},"viscosity","","","","scalar");
	(*this).setVarible({"cell"},"conductivity","","","","scalar");
	(*this).setVarible({"cell"},"permittivity","","","","scalar");
	(*this).setVarible({"cell"},"curvature","","","","scalar");
	
	// 그레디언트 값
	(*this).setVarible({"cell"},"gradient pressure","","","","vector3",
						{"x-gradient pressure","y-gradient pressure","z-gradient pressure"},{"","",""},{"","",""});
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("volume-fraction-"+species[i]);
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
	(*this).setVarible({"cell"},"old charge-density","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("volume-fraction-"+species[i]);
		(*this).setVarible({"cell"},"old "+tmp_name,"","","","scalar");
	}
	
	// 페이스 값 정의
	(*this).setVarible({"face"},"pressure","","","","scalar");
	(*this).setVarible({"face"},"velocity","","","","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},{"","",""});
	for(int i=0; i<species.size()-1; ++i){
		(*this).setVarible({"face"},"volume-fraction-"+species[i],"","","","scalar");
		// (*this).setVarible({"face"},"left volume-fraction-"+species[i],"","","","scalar");
		// (*this).setVarible({"face"},"right volume-fraction-"+species[i],"","","","scalar");
	}
	(*this).setVarible({"face"},"density","","","","scalar");
	(*this).setVarible({"face"},"viscosity","","","","scalar");
	(*this).setVarible({"face"},"charge-density","","","","scalar");
	(*this).setVarible({"face"},"x-electric-field","","","","scalar");
	(*this).setVarible({"face"},"y-electric-field","","","","scalar");
	(*this).setVarible({"face"},"z-electric-field","","","","scalar");
	(*this).setVarible({"face"},"electric-potential","","","","scalar");
	(*this).setVarible({"face"},"conductivity","","","","scalar");
	(*this).setVarible({"face"},"permittivity","","","","scalar");
    
	(*this).setVarible({"face"},"contravariant-velocity","","","","scalar");
	
	
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
	
	
	
	// parcels
	(*this).setVarible({"parcel"},"position","xyz","","","vector3",
						{"x-position","y-position","z-position"},{"u","v","w"},
						{"","",""});
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

