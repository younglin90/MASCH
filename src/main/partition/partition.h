#pragma once 
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string>
#include <map>
#include <random>

#include "parmetis.h" 
#include "scotch.h" 

#include "../../others/mesh.h"  
#include "../../others/load.h" 
#include "../../others/mpi.h"
#include "../../others/controls.h"
#include "../../others/save.h"
#include "../../others/math.h"
using namespace std;

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
