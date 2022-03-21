
#include "../../../others/solvers.h"

void MASCH_Control::setVariablesUDF(vector<string>& species){
	
	// 셀 값 정의
	(*this).setVarible({"cell"},"unknown","phi","","primitive","scalar");
	// (*this).setVarible({"cell"},"delta-unknown","dphi","","","scalar");
	// 페이스 값 정의
	(*this).setVarible({"face"},"left unknown","","","","scalar");
	(*this).setVarible({"face"},"right unknown","","","","scalar");
	// (*this).setVarible({"face"},"left delta-unknown","","","","scalar");
	// (*this).setVarible({"face"},"right delta-unknown","","","","scalar");
	
}

