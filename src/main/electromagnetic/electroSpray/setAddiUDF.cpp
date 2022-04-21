
#include "../../../others/solvers.h"


void MASCH_Solver::setAddiFunctionsUDF(MASCH_Mesh& mesh, MASCH_Control& controls){
	
	auto& solver = (*this);
	
	int nSp = controls.spName.size();
	
	int id_p = controls.getId_cellVar("pressure");
	int id_u = controls.getId_cellVar("x-velocity");
	int id_v = controls.getId_cellVar("y-velocity");
	int id_w = controls.getId_cellVar("z-velocity");
	vector<int> id_alpha;
	for(int i=0; i<controls.spName.size(); ++i){
		id_alpha.push_back(controls.getId_cellVar("volume-fraction-"+controls.spName[i]));
	}
	int id_rho = controls.getId_cellVar("density");
	int id_rhoe = controls.getId_cellVar("charge-density");
	int id_k = controls.getId_cellVar("conductivity");
	int id_epsilon = controls.getId_cellVar("permittivity");
    
    
	int id_pF = controls.getId_faceVar("pressure");
	int id_uF = controls.getId_faceVar("x-velocity");
	int id_vF = controls.getId_faceVar("y-velocity");
	int id_wF = controls.getId_faceVar("z-velocity");
	int id_rhoeF = controls.getId_faceVar("charge-density");
	int id_phiF = controls.getId_faceVar("electric-potential");
    vector<int> id_alphaF, id_alphaL, id_alphaR;
	for(int i=0; i<controls.spName.size()-1; ++i){
        id_alphaF.push_back(controls.getId_faceVar("volume-fraction-"+controls.spName[i]));
        // id_alphaL.push_back(controls.getId_faceVar("left volume-fraction-"+controls.spName[i]));
        // id_alphaR.push_back(controls.getId_faceVar("right volume-fraction-"+controls.spName[i]));
	}
    
	int id_rhoF = controls.getId_faceVar("density");
	int id_muF = controls.getId_faceVar("viscosity");
	int id_xEF = controls.getId_faceVar("x-electric-field");
	int id_yEF = controls.getId_faceVar("y-electric-field");
	int id_zEF = controls.getId_faceVar("z-electric-field");
	int id_kF = controls.getId_faceVar("conductivity");
	int id_epsilonF = controls.getId_faceVar("permittivity");
    
    
    
	vector<double> rhoi(nSp);
	vector<double> mui(nSp);
	vector<double> ki(nSp);
	vector<double> epsiloni(nSp);
	for(int i=0; i<controls.spName.size(); ++i){
		string tmp_name = controls.spName[i];
		if(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.type"] == "constant"){
			rhoi[i] = stod(controls.thermophysicalProperties[tmp_name+".thermodynamics.rho.value"]);
		}
		else{cout << "#WARNING : not defined thermodynamics rho type" << endl;}
		
		if(controls.thermophysicalProperties[tmp_name+".transport.mu.type"] == "constant"){
			mui[i] = stod(controls.thermophysicalProperties[tmp_name+".transport.mu.value"]);
		}
		else{cout << "#WARNING : not defined transport mu type" << endl;}
		
		if(controls.thermophysicalProperties[tmp_name+".transport.k.type"] == "constant"){
			ki[i] = stod(controls.thermophysicalProperties[tmp_name+".transport.k.value"]);
		}
		else{cout << "#WARNING : not defined transport k type" << endl;}
        
		if(controls.thermophysicalProperties[tmp_name+".transport.epsilon.type"] == "constant"){
			epsiloni[i] = stod(controls.thermophysicalProperties[tmp_name+".transport.epsilon.value"]);
		}
		else{cout << "#WARNING : not defined transport epsilon type" << endl;}
	}
	
    
    
    
	
    // 셀 값
	{
		// 1번째
		calcCellAddiVal.push_back(
			[nSp,id_alpha,id_rho,rhoi,mui,ki,epsiloni,
            id_k,id_epsilon,id_mu] (
			double* cells) ->int {
				
				// volume-fraction
                double alpha[nSp];
				double alpha_sum = 0.0;
				for(int i=0; i<nSp-1; ++i){
					alpha[i] = cells[id_alpha[i]];
					alpha_sum += alpha[i];
				}
				alpha[nSp-1] = 1.0 - alpha_sum;
				cells[id_alpha[nSp-1]] = alpha[nSp-1];
				
				// mixture density
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]*rhoi[i];
					}
					cells[id_rho] = tmp_val;
				}
				
				
				// conductivity & permittivity
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]/ki[i];
					}
					cells[id_k] = 1.0/tmp_val;
				}
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]/epsiloni[i];
					}
					cells[id_epsilon] = 1.0/tmp_val;
				}
				
				
				// viscosity
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]*mui[i];
					}
					cells[id_mu] = tmp_val;
				}
				
				return 0;
			}
		);
		
		// 2번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 3번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 4번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 5번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 6번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 7번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
		// 8번째
		calcCellAddiVal.push_back(
			[] (
			double* cells) ->int {
				return 0;
			}
		);
		
    }
		
		
		
	


    // 페이스 값
    {
		// 1번째
		calcFaceAddiVal.push_back( 
			[nSp,id_alphaF,id_rhoF,rhoi,mui,ki,epsiloni,
            id_kF,id_epsilonF,id_muF] (
			double* faces) ->int {
					
				// volume-fraction
				double alpha_sum = 0.0;
                double alpha[nSp];
				for(int i=0; i<nSp-1; ++i){
					alpha[i] = faces[id_alphaF[i]];
					alpha_sum += alpha[i];
				}
				alpha[nSp-1] = 1.0 - alpha_sum;
				faces[id_alphaF[nSp-1]] = alpha[nSp-1];
				
				// mixture density
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]*rhoi[i];
					}
					faces[id_rhoF] = tmp_val;
				}
				
				
				// conductivity & permittivity
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]/ki[i];
					}
					faces[id_kF] = 1.0/tmp_val;
				}
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]/epsiloni[i];
					}
					faces[id_epsilonF] = 1.0/tmp_val;
				}
				
				
				// viscosity
				{
					double tmp_val = 0.0;
					for(int i=0; i<nSp; ++i){
						tmp_val += alpha[i]*mui[i][5];
					}
					faces[id_muF] = tmp_val;
				}
				
				return 0;
			}
		);
		
		// 2번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 3번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 4번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 5번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 6번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 7번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
		
		// 8번째
		calcFaceAddiVal.push_back(
			[] (
			double* faces) ->int {
				return 0;
			}
		);
	}
	
	
}

