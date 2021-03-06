
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
                            
		// string tmp_name2 = ("volume-fraction-"+species[i]);
		// (*this).setVarible({"cell"},"gradient "+tmp_name2,"","","","vector3",
							// {"x-gradient "+tmp_name2,"y-gradient "+tmp_name2,"z-gradient "+tmp_name2},{"","",""},{"","",""});
	}
	(*this).setVarible({"cell"},"gradient density","","","","vector3",
						{"x-gradient density","y-gradient density","z-gradient density"},{"","",""},{"","",""});
	
	
	// 올드 값
	(*this).setVarible({"cell"},"old pressure","","","","scalar");
	(*this).setVarible({"cell"},"old velocity","","","","vector3",
						{"old x-velocity","old y-velocity","old z-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"old temperature","","","","scalar");
	(*this).setVarible({"cell"},"old density","","","","scalar");
	(*this).setVarible({"cell"},"old total-enthalpy","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+species[i]);
		(*this).setVarible({"cell"},"old "+tmp_name,"","","","scalar");
	}
    
	(*this).setVarible({"cell"},"old2 pressure","","","","scalar");
	(*this).setVarible({"cell"},"old2 velocity","","","","vector3",
						{"old2 x-velocity","old2 y-velocity","old2 z-velocity"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"old2 temperature","","","","scalar");
	(*this).setVarible({"cell"},"old2 density","","","","scalar");
	(*this).setVarible({"cell"},"old2 total-enthalpy","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
		string tmp_name = ("mass-fraction-"+species[i]);
		(*this).setVarible({"cell"},"old2 "+tmp_name,"","","","scalar");
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
    
    // 페이스 리미터
    (*this).setVarible({"face"},"left limiter-pressure","","","","scalar");
    (*this).setVarible({"face"},"right limiter-pressure","","","","scalar");
    (*this).setVarible({"face"},"left limiter-x-velocity","","","","scalar");
    (*this).setVarible({"face"},"right limiter-x-velocity","","","","scalar");
    (*this).setVarible({"face"},"left limiter-y-velocity","","","","scalar");
    (*this).setVarible({"face"},"right limiter-y-velocity","","","","scalar");
    (*this).setVarible({"face"},"left limiter-z-velocity","","","","scalar");
    (*this).setVarible({"face"},"right limiter-z-velocity","","","","scalar");
    (*this).setVarible({"face"},"left limiter-temperature","","","","scalar");
    (*this).setVarible({"face"},"right limiter-temperature","","","","scalar");
	for(int i=0; i<species.size()-1; ++i){
        (*this).setVarible({"face"},"left limiter-mass-fraction-"+species[i],"","","","scalar");
        (*this).setVarible({"face"},"right limiter-mass-fraction-"+species[i],"","","","scalar");
    }
	
	
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
	
	
	
	// 셀 평균 플랏팅 관련
	(*this).setVarible({"cell"},"mean-pressure","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-pressure","","","","scalar");
	(*this).setVarible({"cell"},"mean-x-velocity","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-x-velocity","","","","scalar");
	(*this).setVarible({"cell"},"mean-y-velocity","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-y-velocity","","","","scalar");
	(*this).setVarible({"cell"},"mean-z-velocity","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-z-velocity","","","","scalar");
	(*this).setVarible({"cell"},"mean-temperature","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-temperature","","","","scalar");
	(*this).setVarible({"cell"},"mean-density","","","","scalar");
	(*this).setVarible({"field"},"total-time-of-mean-density","","","","scalar");
    
	// (*this).setVarible({"cell"},"mean-liquid-density","","","","scalar");
	// (*this).setVarible({"field"},"total-time-of-mean-liquid-density","","","","scalar");
	// (*this).setVarible({"cell"},"mean-liquid-velocity","","","","scalar");
	// (*this).setVarible({"field"},"total-time-of-mean-liquid-velocity","","","","scalar");
    
    
	for(int i=0; i<species.size()-1; ++i){
        
        (*this).setVarible({"cell"},"mean-mass-fraction-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-mean-mass-fraction-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"mean-volume-fraction-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-mean-volume-fraction-"+species[i],"","","","scalar");
        
        (*this).setVarible({"cell"},"fvm-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"fvm-volume-"+species[i],"","","","scalar");
        
        (*this).setVarible({"cell"},"fvm-mean-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-fvm-mean-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"parcel-mean-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-parcel-mean-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"mean-surface-area-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-mean-surface-area-"+species[i],"","","","scalar");
        
        (*this).setVarible({"cell"},"fvm-mean-volume-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-fvm-mean-volume-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"parcel-mean-volume-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-parcel-mean-volume-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"mean-volume-"+species[i],"","","","scalar");
        (*this).setVarible({"field"},"total-time-of-mean-volume-"+species[i],"","","","scalar");
        
        (*this).setVarible({"cell"},"fvm-sauter-mean-diameter-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"parcel-sauter-mean-diameter-"+species[i],"","","","scalar");
        (*this).setVarible({"cell"},"sauter-mean-diameter-"+species[i],"","","","scalar");
          
        // controls.meanInp_cell_id.push_back(controls.getId_cellVar("fvm-volume-"+item));
        // controls.meanInp_cell_totalTime_id.push_back(controls.getId_fieldVar("total-time-of-fvm-mean-volume-"+item));
        // controls.meanOut_cell_id.push_back(controls.getId_cellVar("fvm-mean-volume-"+item));
        // (*this).setVarible({"field"},"total-time-of-mean-pressure","","","","scalar");
        // (*this).setVarible({"field"},"total-time-of-mean-x-velocity","","","","scalar");
        // (*this).setVarible({"field"},"total-time-of-mean-y-velocity","","","","scalar");
        // (*this).setVarible({"field"},"total-time-of-mean-z-velocity","","","","scalar");
        // (*this).setVarible({"field"},"total-time-of-mean-temperature","","","","scalar");
        
        
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
	(*this).setVarible({"parcel"},"time","","","","scalar");
	(*this).setVarible({"field"},"parcel-injection-accum-time","","","","scalar");
	(*this).setVarible({"field"},"time-step-parcels","","","","scalar");
	
    
    
    
    // 포인트 min max
	(*this).setVarible({"point"},"maximum x-velocity","","","","scalar");
	(*this).setVarible({"point"},"minimum x-velocity","","","","scalar");
	(*this).setVarible({"point"},"maximum y-velocity","","","","scalar");
	(*this).setVarible({"point"},"minimum y-velocity","","","","scalar");
	(*this).setVarible({"point"},"maximum z-velocity","","","","scalar");
	(*this).setVarible({"point"},"minimum z-velocity","","","","scalar");
    
	(*this).setVarible({"point"},"maximum temperature","","","","scalar");
	(*this).setVarible({"point"},"minimum temperature","","","","scalar");
    
	for(int i=0; i<species.size()-1; ++i){
        (*this).setVarible({"point"},"maximum mass-fraction-"+species[i],"","","","scalar");
        (*this).setVarible({"point"},"minimum mass-fraction-"+species[i],"","","","scalar");
        
    }
	
	
	
}

