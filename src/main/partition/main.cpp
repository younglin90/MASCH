
#include "./partition.h" 

void print_help();


int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	
	map<string,string> mapArgv;
	for(int i=1; i<argc; i+=2){
		string first = argv[i]; string second;
		if(i+1==argc){ second = "nan"; }
		else{ second = argv[i+1]; }
		mapArgv.insert(make_pair(first, second));
	}
	
	if(mapArgv.find("-help") != mapArgv.end() ||
	   mapArgv.find("-h") != mapArgv.end() ){
		if(rank==0) print_help();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI::Finalize();
		return EXIT_SUCCESS;
	}
	// set_mapArgv(argc, argv, mapArgv);
	int nSize = stoi(mapArgv["-n"]);
	double scaleMesh = stod(mapArgv["-s"]);
	
	
	// MASCH_Load load;
	MASCH_Mesh_Load load;
	MASCH_Mesh_Save save;
	MASCH_Control control;
	MASCH_Mesh mesh;
	
	load.OpenFoam("./grid/", mesh);
	
	vector<int> procNoCell(static_cast<int>(mesh.cells.size()),0);
	MASCH_Mesh_Partition partition;
	partition.parMETIS_Graph_Partition(nSize, procNoCell, mesh);
	vector<MASCH_Mesh> newMesh(nSize,MASCH_Mesh());
	partition.partitionFromSerial(nSize, procNoCell, mesh, newMesh);
	
	// 스케일
	for(auto& item : newMesh) {
		for(auto& point : item.points){
			point.x *= scaleMesh;
			point.y *= scaleMesh;
			point.z *= scaleMesh;
		}
	}

	{
		cout.precision(20);
		for(int ip=0; ip<nSize; ++ip){

			// SEMO_Utility_Math math;
			// SEMO_Mesh_Geometric geometric;
			// geometric.init(newMesh[ip]);
			// cout << "AAA" << endl;
			MASCH_Mesh_Save save;
			save.vtu("./grid/0/", ip, newMesh[ip]);
		}
	}
	
	MPI::Finalize();
	return EXIT_SUCCESS;
}




void print_help(){

	cout << endl;
	cout << "┌─────── Partitioning helper ─────────────────────────────── " << endl;
	cout << "| -n \"int num.\"   : # of partition" << endl;
	cout << "| -s \"real num.\"  : scale of mesh" << endl;
	cout << "└───────────────────────────────────────────────────────────────── " << endl;
	cout << endl;
}