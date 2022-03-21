
#include "../../../others/solvers.h"

void MASCH_Control::setVariablesUDF(vector<string>& species){
	
	// 셀 값 정의
	(*this).setVarible({"cell"},"pressure","p","","primitive","scalar");
	(*this).setVarible({"cell"},"delta-pressure","p","","","scalar");
	(*this).setVarible({"cell"},"velocity","U","","primitive","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{"primitive","primitive","primitive"});
	(*this).setVarible({"cell"},"gradient pressure","","","","vector3",
						{"x-gradient pressure","y-gradient pressure","z-gradient pressure"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"gradient delta-pressure","","","","vector3",
						{"x-gradient delta-pressure","y-gradient delta-pressure","z-gradient delta-pressure"},{"","",""},{"","",""});
	(*this).setVarible({"cell"},"old velocity","","","","vector3",
						{"old x-velocity","old y-velocity","old z-velocity"},{"","",""},{"","",""});
	// 페이스 값 정의
	(*this).setVarible({"face"},"pressure","","","","scalar");
	(*this).setVarible({"face"},"delta-pressure","","","","scalar");
	(*this).setVarible({"face"},"contravariant-velocity","","","","scalar");
	// (*this).setVarible({"face"},"time-step-density","","","","scalar");
	(*this).setVarible({"face"},"velocity","","","","vector3",
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},{"","",""});
	
}

