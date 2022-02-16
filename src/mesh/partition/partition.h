#pragma once

#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
// #include <dlfnc.h>
using namespace std;

#include "../mesh.h" 
#include "../geometric.h" 
#include "../../load/load.h"
#include "../../save/save.h"
#include "../../mpi/mpi.h"
#include "../../math/math.h"
// #include "../../mpi/mpi.h" 


enum class MASCH_Partition_Types{
	IN2IN,
	PR2IN, 
	BC2BC, 
	IN2PR, 
	PR2PR
};


class MASCH_Mesh_Partition{
public:
	void parMETIS_Graph_Partition(int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh);
	
	void partitionFromSerial(
		int nBlocks, vector<int>& idBlockCell, MASCH_Mesh &mesh, vector<MASCH_Mesh>& newMesh);
		
	void combine(vector<int>& idBlockCell, MASCH_Mesh &mesh);
};
