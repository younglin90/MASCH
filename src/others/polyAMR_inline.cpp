
#include "./mesh.h"
#include "./polyAMR.h"
// #include "geometric.h" 
#include "./mpi.h"
#include "./solvers.h" 
#include "./save.h"  


void calc_indicatorValues(
MASCH_Mesh& mesh, 
MASCH_Variables& var,
vector<int>& indicatorAMR_id, vector<vector<double>>& indicatorValues
){
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// 인디케이터 값 초기화
	indicatorValues[0].clear();
	indicatorValues[0].resize(mesh.cells.size(),0.0);
	for(int i=0; i<mesh.cells.size(); ++i){
		int id_gradx = indicatorAMR_id[0];
		int id_grady = indicatorAMR_id[1];
		int id_gradz = indicatorAMR_id[2];
		double mag_grad = sqrt(
			var.cells[i][id_gradx]*var.cells[i][id_gradx]+
			var.cells[i][id_grady]*var.cells[i][id_grady]+
			var.cells[i][id_gradz]*var.cells[i][id_gradz]);
		indicatorValues[0][i] = mag_grad;
	}
	
	// 인디케이터 값 Buffer layer
	for(int iter=0; iter<3; ++iter){
		
		vector<double> newIndicatorAMR(mesh.cells.size(),0.0);
		
		vector<double> recvValues;
		if(size>1){
			// processor faces
			vector<double> sendValues;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == MASCH_Face_Types::PROCESSOR){
					sendValues.push_back(indicatorValues[0][face.iL]);
				}
			}
			recvValues.resize(sendValues.size());
			MPI_Alltoallv( sendValues.data(), mesh.countsProcFaces.data(), 
							mesh.displsProcFaces.data(), MPI_DOUBLE, 
							recvValues.data(), mesh.countsProcFaces.data(), 
							mesh.displsProcFaces.data(), MPI_DOUBLE, 
						   MPI_COMM_WORLD);
		}
		
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double maxInd = 0.0;
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				maxInd = max(indicatorValues[0][face.iL],
							 indicatorValues[0][face.iR]);
				newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
				newIndicatorAMR[face.iR] = max(newIndicatorAMR[face.iR],maxInd);
			}
			else if(face.getType() == MASCH_Face_Types::PROCESSOR){
				maxInd = max(indicatorValues[0][face.iL], recvValues[ip]);
				newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
				++ip;
			}
			else{
				maxInd = indicatorValues[0][face.iL];
				newIndicatorAMR[face.iL] = max(newIndicatorAMR[face.iL],maxInd);
			}
		}
		
		for(int i=0; i<mesh.cells.size(); ++i){
			indicatorValues[0][i] = newIndicatorAMR[i];
		}
	}
}


void MASCH_Poly_AMR_Builder::polyAMR_inline(
	MASCH_Mesh& mesh, 
	MASCH_Control& controls,
	MASCH_Solver& solver,
	MASCH_Variables& var,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// cout << iter+1 % 10 << endl;
	// for(int ii=0; ii<39; ++ii)
	if( (iter+1) % 10 == 0)
	{
		int maxLevel = 3;
		
		vector<vector<double>> indicatorCriterion(1);
		indicatorCriterion[0].push_back(10.0);
		indicatorCriterion[0].push_back(50.0);
		indicatorCriterion[0].push_back(100.0);
		
		// 그레디언트 mag 계산
		vector<vector<double>> indicatorValues(1);
		vector<int> indicatorAMR_id;
		indicatorAMR_id.push_back(controls.getId_cellVar("x-gradient pressure"));
		indicatorAMR_id.push_back(controls.getId_cellVar("y-gradient pressure"));
		indicatorAMR_id.push_back(controls.getId_cellVar("z-gradient pressure"));
		
	
		// 리파인
		{
			// 그레디언트 계산
			solver.gradientTerms(mesh, controls, var);
			calc_indicatorValues(mesh, var, indicatorAMR_id, indicatorValues);
			
			// 원래 셀의 x,y,z 저장
			vector<vector<double>> org_xyz(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i) {
				auto& cell = mesh.cells[i];
				org_xyz[i].push_back(cell.x); org_xyz[i].push_back(cell.y); org_xyz[i].push_back(cell.z);
			}
			// if(rank==0) cout << "| exe. Poly AMR Refinement" << endl;
			
			vector<vector<int>> child_new_cell_id_of_org;
			polyRefine(mesh, controls, 
				maxLevel, 50000, 1.e-16, 
				indicatorCriterion, indicatorValues, 
				child_new_cell_id_of_org, 0);
			
			controls.setGeometricOnlyCell_xyz(mesh);
			
			controls.resetVariableArray(mesh, var, org_xyz, child_new_cell_id_of_org, "refine");
			controls.setGeometric(mesh, var);
			var.setSparCSR(mesh, controls);
			solver.calcGradient.init(mesh, controls, var);
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		// 리파티셔닝
		{
			// vector<vector<double>> dummy;
			// if(rank==0) cout << "| exe. Dynamic Load Balancing" << endl;
			vector<int> cell_ip(mesh.cells.size(),rank);
			mesh.repartParMETIS(size, cell_ip, mesh);
			vector<int> to_new_cell_id;
			mesh.repartitioning(cell_ip, maxLevel, to_new_cell_id);
			
			controls.setGeometricOnlyCell_xyz(mesh);
			
			controls.resetVariableArray(mesh, var, cell_ip, to_new_cell_id, "repart");
			controls.setGeometric(mesh, var);
			var.setSparCSR(mesh, controls);
			solver.calcGradient.init(mesh, controls, var);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		// 언리파인
		{
			// 그레디언트 계산
			solver.gradientTerms(mesh, controls, var);
			calc_indicatorValues(mesh, var, indicatorAMR_id, indicatorValues);
			
			vector<vector<double>> dummy;
			// if(rank==0) cout << "| exe. Poly AMR Unrefinement" << endl;
			vector<vector<int>> child_org_cell_id_of_new;
			polyUnrefine(mesh, controls, 
				indicatorCriterion, indicatorValues, 
				child_org_cell_id_of_new, 0);
			
			controls.resetVariableArray(mesh, var, dummy, child_org_cell_id_of_new, "unrefine");
			controls.setGeometric(mesh, var);
			var.setSparCSR(mesh, controls);
			solver.calcGradient.init(mesh, controls, var);
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// mesh.debug_procFace_unitNomals(0.8);
		}
		
		
		
	}
	
	
	
}

