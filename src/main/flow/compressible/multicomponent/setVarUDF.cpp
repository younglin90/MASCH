
#include "../../../../others/solvers.h"

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
	vector<string> mass_frac_species;
	vector<string> species_abb;
	vector<string> species_diff_rho;
	vector<string> species_diff_Ht;
	vector<string> species_diff_rho_abb;
	vector<string> species_diff_Ht_abb;
	vector<string> species_shape_prim;
	vector<string> species_shape;
	vector<string> left_species;
	vector<string> right_species;
	for(int i=0, SIZE=species.size(); i<SIZE; ++i){
		string Y_abb = "Y";
		Y_abb += to_string(i);
		species_abb.push_back(Y_abb);
		
		// mass_frac_species.push_back("mass-fraction-" + species[i]);
		mass_frac_species.push_back(species[i]);
		
		string diff_rho = "density-diff-with ";
		string diff_rho_abb = "drhodY";
		string diff_Ht = "total-enthalpy-diff-with ";
		string diff_Ht_abb = "dHtdY";
		diff_rho += ("mass-fraction-" + species[i]);
		diff_Ht += ("mass-fraction-" + species[i]);
		diff_rho_abb += to_string(i);
		diff_Ht_abb += to_string(i);
		species_diff_rho.push_back(diff_rho);
		species_diff_Ht.push_back(diff_Ht);
		species_diff_rho_abb.push_back(diff_rho_abb);
		species_diff_Ht_abb.push_back(diff_Ht_abb);
		
		string tmp_left_species = "left ";
		string tmp_right_species = "right ";
		tmp_left_species += ("mass-fraction-" + species[i]);
		tmp_right_species += ("mass-fraction-" + species[i]);
		left_species.push_back(tmp_left_species);
		right_species.push_back(tmp_right_species);
	}
	for(int i=0, SIZE=species.size()-1; i<SIZE; ++i){
		species_shape_prim.push_back(prim);
		species_shape.push_back("thermo");
	}
	species_shape_prim.push_back("thermo");
	species_shape.push_back("thermo");
	
	// 셀 값 정의
	(*this).setVarible({cell},"pressure","p","Pa",prim,scal);
	(*this).setVarible({cell},"velocity","U","m/s",prim,vec3,
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{prim,prim,prim});
	(*this).setVarible({cell},"temperature","T","K",prim,scal);
	(*this).setVarible({cell},"mass-fraction","Y","",prim,vec,
						mass_frac_species,species_abb,species_shape_prim);
						
	(*this).setVarible({cell},"pressure-old","p","Pa","old",scal);
	(*this).setVarible({cell},"velocity-old","U","m/s","old",vec3,
						{"x-velocity-old","y-velocity-old","z-velocity-old"},{"u","v","w"},
						{"old","old","old"});
	(*this).setVarible({cell},"density-old","rho","kg/m^3","old",scal);
	(*this).setVarible({cell},"total-enthalpy-old","Ht","J/kg","old",scal);
						
	(*this).setVarible({cell},"pressure-old2","p","Pa","old",scal);
	(*this).setVarible({cell},"velocity-old2","U","m/s","old",vec3,
						{"x-velocity-old2","y-velocity-old2","z-velocity-old2"},{"u","v","w"},
						{"old","old","old"});
	(*this).setVarible({cell},"density-old2","rho","kg/m^3","old",scal);
	(*this).setVarible({cell},"total-enthalpy-old2","Ht","J/kg","old",scal);
	
	(*this).setVarible({cell},"density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({cell},"speed-of-sound","c","m/s",thermo,scal);
	(*this).setVarible({cell},"total-enthalpy","Ht","J/kg",thermo,scal);
	
	(*this).setVarible({cell},"density-diff-with-pressure","drhodp","kg/m^3/Pa",thermo,scal);
	(*this).setVarible({cell},"density-diff-with-temperature","drhodT","kg/m^3/K",thermo,scal);
	(*this).setVarible({cell},"density-diff-with-mass-fraction","drhodY","kg/m^3",thermo,vec,
						mass_frac_species,species_diff_rho_abb,species_shape);
	(*this).setVarible({cell},"total-enthalpy-diff-with-pressure","dHtdp","J/kg/Pa",thermo,scal);
	(*this).setVarible({cell},"total-enthalpy-diff-with-temperature","dHtdT","J/kg/K",thermo,scal);
	(*this).setVarible({cell},"total-enthalpy-diff-with-mass-fraction","dHtdY","J/kg",thermo,vec,
						mass_frac_species,species_diff_Ht_abb,species_shape);
	
	(*this).setVarible({cell},"gradient pressure","dpdX","","gradient",vec3,
			{"x-gradient pressure","y-gradient pressure","z-gradient pressure"},
			{"dpdx","dpdy","dpdz"},{"","",""});
			
	(*this).setVarible({cell},"gradient x-velocity","dudX","","gradient",vec3,
			{"x-gradient x-velocity","y-gradient x-velocity","z-gradient x-velocity"},
			{"dudx","dudy","dudz"},{"","",""});
			
	(*this).setVarible({cell},"gradient y-velocity","dvdX","","gradient",vec3,
			{"x-gradient y-velocity","y-gradient y-velocity","z-gradient y-velocity"},
			{"dvdx","dvdy","dvdz"},{"","",""});
			
	(*this).setVarible({cell},"gradient z-velocity","dwdX","","gradient",vec3,
			{"x-gradient z-velocity","y-gradient z-velocity","z-gradient z-velocity"},
			{"dwdx","dwdy","dwdz"},{"","",""});
			
	(*this).setVarible({cell},"gradient temperature","dTdX","","gradient",vec3,
			{"x-gradient temperature","y-gradient temperature","z-gradient temperature"},
			{"dwdx","dwdy","dwdz"},{"","",""});
			
	for(int i=0, SIZE=species.size(); i<SIZE; ++i){
		string GradName = "gradient ";
		string xGradName = "x-gradient ";
		string yGradName = "y-gradient ";
		string zGradName = "z-gradient ";
		GradName += "mass-fraction-"+species[i];
		xGradName += "mass-fraction-"+species[i];
		yGradName += "mass-fraction-"+species[i];
		zGradName += "mass-fraction-"+species[i];
			
		(*this).setVarible({cell},GradName,"dYdX","","gradient",vec3,
				{xGradName,yGradName,zGradName},
				{"dYdx","dYdy","dYdz"},{"","",""});
	}
	
	(*this).setVarible({cell},"gradient density","drhodX","","gradient",vec3,
			{"x-gradient density","y-gradient density","z-gradient density"},
			{"drhodx","drhody","drhodz"},{"","",""});
			
	(*this).setVarible({cell},"viscosity","mu","Pa*s",thermo,scal);
	
	// 페이스 값 정의
	(*this).setVarible({face},"left pressure","p","Pa","value",scal);
	(*this).setVarible({face},"left velocity","U","m/s","value",vec3,
						{"left x-velocity","left y-velocity","left z-velocity"},{"u","v","w"},{"","",""});
	(*this).setVarible({face},"left temperature","T","K","value",scal);
	(*this).setVarible({face},"left mass-fraction","Y","","value",vec,
						mass_frac_species,species_abb,species_shape);
	(*this).setVarible({face},"left density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({face},"left speed-of-sound","c","m/s",thermo,scal);
	(*this).setVarible({face},"left total-enthalpy","Ht","J/kg",thermo,scal);
	
	(*this).setVarible({face},"right pressure","p","Pa","value",scal);
	(*this).setVarible({face},"right velocity","U","m/s","value",vec3,
						{"right x-velocity","right y-velocity","right z-velocity"},{"u","v","w"},{"","",""});
	(*this).setVarible({face},"right temperature","T","K","value",scal);
	(*this).setVarible({face},"right mass-fraction","Y","","value",vec,
						mass_frac_species,species_abb,species_shape);
	(*this).setVarible({face},"right density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({face},"right speed-of-sound","c","m/s",thermo,scal);
	(*this).setVarible({face},"right total-enthalpy","Ht","J/kg",thermo,scal);
	
	
	// 필드 값 정의
	{
		vector<string> species_material_name;
		vector<string> species_material_abb;
		vector<string> species_material_unit;
		for(int i=0, SIZE=species.size(); i<SIZE; ++i){
			{
				string tmp_name = species[i];
				tmp_name += " pinf";
				(*this).setVarible({field},"-step","","","thermo",scal);
			}
			{
				string tmp_name = species[i];
				tmp_name += " cv";
				(*this).setVarible({field},"-step","","","thermo",scal);
			}
			{
				string tmp_name = species[i];
				tmp_name += " gamma";
				(*this).setVarible({field},"-step","","","thermo",scal);
			}
			{
				string tmp_name = species[i];
				tmp_name += " bNASG";
				(*this).setVarible({field},"-step","","","thermo",scal);
			}
			{
				string tmp_name = species[i];
				tmp_name += " q";
				(*this).setVarible({field},"-step","","","thermo",scal);
			}
		}
	}
}

