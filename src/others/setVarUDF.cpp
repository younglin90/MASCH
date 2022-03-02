
#include "./solvers.h"

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
	string cell = "cell";
	string face = "face";
	string point = "point";
	string field = "field";
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
		
		string diff_rho = "density diff with ";
		string diff_rho_abb = "drhodY";
		string diff_Ht = "total enthalpy diff with ";
		string diff_Ht_abb = "dHtdY";
		diff_rho += species[i];
		diff_Ht += species[i];
		diff_rho_abb += to_string(i);
		diff_Ht_abb += to_string(i);
		species_diff_rho.push_back(diff_rho);
		species_diff_Ht.push_back(diff_Ht);
		species_diff_rho_abb.push_back(diff_rho_abb);
		species_diff_Ht_abb.push_back(diff_Ht_abb);
		
		string tmp_left_species = "left ";
		string tmp_right_species = "left ";
		tmp_left_species += species[i];
		tmp_right_species += species[i];
		left_species.push_back(tmp_left_species);
		right_species.push_back(tmp_right_species);
	}
	for(int i=0, SIZE=species.size()-1; i<SIZE; ++i){
		species_shape_prim.push_back(prim);
		species_shape.push_back("thermo");
	}
	species_shape_prim.push_back("thermo");
	species_shape.push_back("thermo");
	
	(*this).setVarible({cell},"pressure","p","Pa",prim,scal);
	(*this).setVarible({cell},"velocity","U","m/s",prim,vec,
						{"x-velocity","y-velocity","z-velocity"},{"u","v","w"},
						{prim,prim,prim});
	(*this).setVarible({cell},"temperature","T","K",prim,scal);
	(*this).setVarible({cell},"mass fraction","Y","",prim,vec,
						species,species_abb,species_shape_prim);
						
	(*this).setVarible({cell},"old pressure","p","Pa","old",scal);
	(*this).setVarible({cell},"old velocity","U","m/s","old",vec,
						{"old x-velocity","old y-velocity","old z-velocity"},{"u","v","w"},
						{"old","old","old"});
	(*this).setVarible({cell},"old density","rho","kg/m^3","old",scal);
	(*this).setVarible({cell},"old total enthalpy","Ht","J/kg","old",scal);
						
	(*this).setVarible({cell},"old2 pressure","p","Pa","old",scal);
	(*this).setVarible({cell},"old2 velocity","U","m/s","old",vec,
						{"old2 x-velocity","old2 y-velocity","old2 z-velocity"},{"u","v","w"},
						{"old","old","old"});
	(*this).setVarible({cell},"old2 density","rho","kg/m^3","old",scal);
	(*this).setVarible({cell},"old2 total enthalpy","Ht","J/kg","old",scal);
	
	(*this).setVarible({cell},"density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({cell},"speed of sound","c","m/s",thermo,scal);
	(*this).setVarible({cell},"total enthalpy","Ht","J/kg",thermo,scal);
	
	(*this).setVarible({cell},"density diff with pressure","drhodp","kg/m^3/Pa",thermo,scal);
	(*this).setVarible({cell},"density diff with temperature","drhodT","kg/m^3/K",thermo,scal);
	(*this).setVarible({cell},"density diff with mass fraction","drhodY","kg/m^3",thermo,vec,
						species_diff_rho,species_diff_rho_abb,species_shape);
	(*this).setVarible({cell},"total enthalpy diff with pressure","dHtdp","J/kg/Pa",thermo,scal);
	(*this).setVarible({cell},"total enthalpy diff with temperature","dHtdT","J/kg/K",thermo,scal);
	(*this).setVarible({cell},"total enthalpy diff with mass fraction","dHtdY","J/kg",thermo,vec,
						species_diff_Ht,species_diff_Ht_abb,species_shape);
	
	(*this).setVarible({cell},"gradient pressure","dpdX","","gradient",vec,
			{"x-gradient pressure","y-gradient pressure","z-gradient pressure"},
			{"dpdx","dpdy","dpdz"},
			{"gradient","gradient","gradient"});
			
	(*this).setVarible({cell},"gradient x-velocity","dudX","","gradient",vec,
			{"x-gradient x-velocity","y-gradient x-velocity","z-gradient x-velocity"},
			{"dudx","dudy","dudz"},
			{"gradient","gradient","gradient"});
			
	(*this).setVarible({cell},"gradient y-velocity","dvdX","","gradient",vec,
			{"x-gradient y-velocity","y-gradient y-velocity","z-gradient y-velocity"},
			{"dvdx","dvdy","dvdz"},
			{"gradient","gradient","gradient"});
			
	(*this).setVarible({cell},"gradient z-velocity","dwdX","","gradient",vec,
			{"x-gradient z-velocity","y-gradient z-velocity","z-gradient z-velocity"},
			{"dwdx","dwdy","dwdz"},
			{"gradient","gradient","gradient"});
			
	(*this).setVarible({cell},"viscosity","mu","Pa*s",thermo,scal);
	
	(*this).setVarible({face},"left pressure","p","Pa","value",scal);
	(*this).setVarible({face},"left velocity","U","m/s","value",vec,
						{"left x-velocity","left y-velocity","left z-velocity"},{"u","v","w"},
						{"value","value","value"});
	(*this).setVarible({face},"left temperature","T","K","value",scal);
	(*this).setVarible({face},"left mass fraction","Y","","value",vec,
						left_species,species_abb,species_shape);
	(*this).setVarible({face},"left density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({face},"left speed of sound","c","m/s",thermo,scal);
	(*this).setVarible({face},"left total enthalpy","Ht","J/kg",thermo,scal);
	
	(*this).setVarible({face},"right pressure","p","Pa","value",scal);
	(*this).setVarible({face},"right velocity","U","m/s","value",vec,
						{"right x-velocity","right y-velocity","right z-velocity"},{"u","v","w"},
						{"value","value","value"});
	(*this).setVarible({face},"right temperature","T","K","value",scal);
	(*this).setVarible({face},"right mass fraction","Y","","value",vec,
						right_species,species_abb,species_shape);
	(*this).setVarible({face},"right density","rho","kg/m^3",thermo,scal);
	(*this).setVarible({face},"right speed of sound","c","m/s",thermo,scal);
	(*this).setVarible({face},"right total enthalpy","Ht","J/kg",thermo,scal);
	
	
}

