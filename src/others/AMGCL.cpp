#include "./solvers.h"
// #include "../mpi/mpi.h"

// #define AMGCL_NO_BOOST

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/schur_pressure_correction.hpp>
#include <amgcl/mpi/block_preconditioner.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/mpi/cpr.hpp>


#include <amgcl/mpi/coarsening/runtime.hpp>
#include <amgcl/mpi/relaxation/runtime.hpp>
#include <amgcl/mpi/partition/runtime.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>
#include <amgcl/mpi/solver/runtime.hpp>


#include <amgcl/mpi/relaxation/as_preconditioner.hpp>

#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/amg.hpp>
#include <amgcl/mpi/relaxation/spai0.hpp>
#include <amgcl/mpi/relaxation/spai1.hpp>
#include <amgcl/mpi/coarsening/aggregation.hpp>
#include <amgcl/mpi/coarsening/smoothed_aggregation.hpp>
#include <amgcl/mpi/solver/bicgstab.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>


#include <amgcl/mpi/relaxation/ilu0.hpp>
#include <amgcl/mpi/relaxation/gauss_seidel.hpp>
#include <amgcl/mpi/solver/preonly.hpp>
#include <amgcl/mpi/solver/idrs.hpp>
#include <amgcl/mpi/relaxation/iluk.hpp>
#include <amgcl/mpi/solver/fgmres.hpp>
#include <amgcl/mpi/solver/lgmres.hpp>
#include <amgcl/mpi/solver/gmres.hpp>
#include <amgcl/mpi/solver/richardson.hpp>

// #include <amgcl/adapter/reorder.hpp>
#include <amgcl/mpi/direct_solver/skyline_lu.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/scope_exit.hpp>




void MASCH_Solver::solveAMGCL(
	vector<int>& istr_CSR, vector<int>& j_CSR, vector<double>& Aval, 
	vector<double>& Bval, vector<double>& Xval){
		
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, var.i_str_CSR[iSegEq], var.j_displ_CSR[iSegEq], var.Avalues[iSegEq]));
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, istr_CSR, j_CSR, Aval));
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	
	// // COO to CSR
	// int strRow = mesh.startCellGlobal * B_n;
	
    // //compute number of non-zero entries per row of A 
	// int n_row = B_vals.size();
	// int nnz = A_rows.size();
	
	// vector<double> val(nnz,0.0);
	// vector<int> col(nnz,0);
	// vector<int> ptr(n_row+1,0);

    // for (int n = 0; n < nnz; n++){            
        // ptr[A_rows[n]-strRow]++;
    // }

	// // cout << "AAAAAAAAAAAAAAA" << endl;
    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = ptr[i];
        // ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // ptr[n_row] = nnz; 

	// // cout << "BBBBBBBB" << endl;
    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = A_rows[n] - strRow;
        // int dest = ptr[row];
	// // cout << row << " " << dest << endl;

        // col[dest] = A_cols[n];
        // val[dest] = A_vals[n];
		
		// if(A_cols[n]>= mesh.ncellsTotal*B_n){
			// cout << "ERROR " << A_cols[n] << " " << mesh.ncellsTotal*B_n <<  endl;
		// }

        // ptr[row]++;
    // }
	// // MPI_Barrier(MPI_COMM_WORLD);

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = ptr[i];
        // ptr[i]  = last;
        // last   = temp;
    // }

	// // cout << "1" << endl;
    // //now ptr,col,val form a CSR representation (with possible duplicates)
	
	amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // The profiler:
    // amgcl::profiler<> prof("Flows Solve");
	
	
	// ptrdiff_t chunk = B_vals.size();

    // Compose the solver type
    typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    // typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    // typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend

	// cout << ptr.size() << " " << col.size() << " " << val.size() << endl;

	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, ptr, col, val));

	ptrdiff_t chunk = Bval.size();
	// cout << chunk << endl;
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, var.i_str_CSR[iSegEq], var.j_displ_CSR[iSegEq], var.Avalues[iSegEq]));
	auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     world, std::tie(chunk, istr_CSR, j_CSR, Aval));

	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// if(equation == "coupled")
	{
			
		// iterative solver
		typedef amgcl::mpi::make_solver<
			amgcl::mpi::relaxation::as_preconditioner<
			// amgcl::mpi::relaxation::spai0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			>,
			// amgcl::mpi::amg<
			// SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// amgcl::mpi::solver::bicgstab<SBackend>
			amgcl::mpi::solver::fgmres<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// amgcl::mpi::solver::preonly<SBackend>
			// amgcl::mpi::solver::lgmres<SBackend>
		> Solver; 
		
		
		Solver::params prm;
		 
		// prm.precond.direct_coarse = true;
		// prm.precond.npre = 1;
		// prm.precond.npost = 5;
		// prm.precond.pre_cycles = 5;
		
		// prm.precond.coarsening.over_interp = 1.0;
		
		// prm.precond.k = 0;
		// prm.precond.damping = 3.0;
		// // // prm.precond.relax.k = 2;
		// // // prm.solver.K = 5;
		// // // prm.solver.M = 50;
		// // // prm.solver.s = 8;
		// // // prm.solver.omega = 0.7;
		prm.solver.maxiter = 50;//1000;
		prm.solver.tol = 1.e-8;
		// // // prm.solver.ns_search = true;
		// // // prm.solver.replacement = true;
		// // // prm.solver.smoothing = true;
		
		
		
		// typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::relaxation::as_preconditioner<
			// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// amgcl::mpi::amg<
			// SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
			// // amgcl::mpi::solver::lgmres<SBackend>
			// // amgcl::mpi::solver::fgmres<SBackend>
		// > Solver; 
		// Solver::params prm;

		// // prm.precond.relax.k = 2;
		// // prm.precond.npre = 1;
		// // prm.precond.npost = 5;
		// // prm.precond.coarsening.over_interp = 1.2;
		// // prm.precond.coarsening.relax = 0.7;

		// // prm.solver.ns_search = true;
		// // prm.solver.verbose = false;
		// // if(world.rank == 0) prm.solver.verbose = true;
		// // prm.solver.tol = 1.e-6;//1e-9;
		// // prm.solver.replacement = true;
		// prm.solver.maxiter = 50;//1000;
		// // prm.solver.smoothing = true;
		
		
		
		// prof.tic("setup");
		Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		int iters; double error;
		std::tie(iters, error) = solve(*A, Bval, Xval);
		// prof.toc("solve");
		
		if (world.rank == 0) {
			cout << scientific; cout.precision(2);
			// cout << solve << endl;
			// cout << prof << endl;
			cout << "| nLinSol = " << iters << ", resiLinSol = " << error << endl;
			cout << fixed; cout.precision(0);
		}
		

	
	
	}
	// else if(equation == "volume_fraction"){
		
		// typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::relaxation::as_preconditioner<
			// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// amgcl::mpi::amg<
			// SBackend,
			// // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// amgcl::mpi::coarsening::aggregation<SBackend>,
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// // amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
			// amgcl::mpi::solver::fgmres<SBackend>
			// // amgcl::mpi::solver::lgmres<SBackend>
		// > Solver; 
		// Solver::params prm;

		// // prm.solver.tol = 1e-10;
		// // prm.solver.tol = 1e-12;
		// prm.precond.npre = 1;
		// prm.precond.npost = 5;
		// // prm.precond.ncycle = 1;
		// // prm.precond.relax.k = 2;
		// // // prm.precond.max_levels = 40;
		// prm.precond.coarsening.over_interp = 1.8;
		// // prm.solver.K = 4;
		// // prm.solver.M = 35;
		// prm.solver.tol = 1e-9;
		// // prm.solver.relTol = 1.e-200;
		// prm.solver.maxiter = 100;//1000;
		
		// // prm.solver.ns_search = true;
		// // prm.solver.verbose = false;
		// // if(world.rank == 0) prm.solver.verbose = true;
		// // prm.solver.maxiter = 50;
		// // prm.solver.replacement = true;
		// // prm.solver.smoothing = true;
		
		
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		
		
		// // if (world.rank == 0) {
			// // cout << solve << endl;
			// // cout << prof << endl;
			// // cout << "|--- volume frac. : " << iters << " " << error << endl;
		// // }
		

	// }
	// else if(equation == "pressure"){
		
		// typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::relaxation::as_preconditioner<
			// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// amgcl::mpi::amg<
			// SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
			// // amgcl::mpi::solver::lgmres<SBackend>
			// // amgcl::mpi::solver::fgmres<SBackend>
		// > Solver; 
		// Solver::params prm;

		// // prm.precond.relax.k = 2;
		// prm.precond.npre = 1;
		// prm.precond.npost = 5;
		// // prm.precond.coarsening.over_interp = 1.2;
		// // prm.precond.coarsening.relax = 0.7;

		// // prm.solver.ns_search = true;
		// // prm.solver.verbose = false;
		// // if(world.rank == 0) prm.solver.verbose = true;
		// prm.solver.tol = 1e-9;
		// // prm.solver.replacement = true;
		// prm.solver.maxiter = 1000;
		// // prm.solver.smoothing = true;
		

		
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		

		// // if (world.rank == 0) {
			// // cout << solve << endl;
			// // cout << prof << endl;
			// // cout << "|--- pressure : " << iters << " " << error << endl;
		// // }

	// }
	

	
	
	
}




// //+++++++++++++++++++++++++++
// // Momentum
// //+++++++++++++++++++++++++++
// void SEMO_Solvers_Builder::solveAMGCL(
	// string equation,
	// SEMO_Mesh_Builder& mesh,
	// int B_n, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B0_vals, vector<double>& B1_vals, vector<double>& B2_vals,
	// vector<double>& resiVar0, vector<double>& resiVar1, vector<double>& resiVar2){
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 
	
	
	// // COO to CSR
	// int strRow = mesh.startCellGlobal * B_n;
	
    // //compute number of non-zero entries per row of A 
	// int n_row = B0_vals.size();
	// int nnz = A_rows.size();
	
	// vector<double> val(nnz,0.0);
	// vector<int> col(nnz,0);
	// vector<int> ptr(n_row+1,0);

    // for (int n = 0; n < nnz; n++){            
        // ptr[A_rows[n]-strRow]++;
    // }

	// // cout << "AAAAAAAAAAAAAAA" << endl;
    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = ptr[i];
        // ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // ptr[n_row] = nnz; 

	// // cout << "BBBBBBBB" << endl;
    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = A_rows[n] - strRow;
        // int dest = ptr[row];
	// // cout << row << " " << dest << endl;

        // col[dest] = A_cols[n];
        // val[dest] = A_vals[n];
		
		// if(A_cols[n]>= mesh.ncellsTotal*B_n){
			// cout << "ERROR " << A_cols[n] << " " << mesh.ncellsTotal*B_n <<  endl;
		// }

        // ptr[row]++;
    // }
	// // MPI_Barrier(MPI_COMM_WORLD);

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = ptr[i];
        // ptr[i]  = last;
        // last   = temp;
    // }

	// // cout << "1" << endl;
    // //now ptr,col,val form a CSR representation (with possible duplicates)
	
	// amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // // The profiler:
    // amgcl::profiler<> prof("Flows Solve");
	
	
	// ptrdiff_t chunk = B0_vals.size();

    // // Compose the solver type
    // typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    // typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    // //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    // typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend

	// // cout << ptr.size() << " " << col.size() << " " << val.size() << endl;

	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, ptr, col, val));

	// typedef amgcl::mpi::make_solver<
		// // amgcl::mpi::relaxation::as_preconditioner<
		// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
		// // >,
		// amgcl::mpi::amg<
		// SBackend,
		// // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
		// amgcl::mpi::coarsening::aggregation<SBackend>,
		// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
		// >,
		// // amgcl::mpi::solver::bicgstab<SBackend>
		// // amgcl::mpi::solver::idrs<SBackend>
		// // amgcl::mpi::solver::preonly<SBackend>
		// // amgcl::mpi::solver::lgmres<SBackend>
		// amgcl::mpi::solver::fgmres<SBackend>
	// > Solver; 
	// Solver::params prm;
	
	
	// prm.precond.npre = 1;
	// prm.precond.npost = 5;
	// // prm.precond.ncycle = 2;
	// // prm.precond.pre_cycles = 2;
	// prm.precond.coarsening.over_interp = 1.0;
	// // prm.precond.coarsening.relax = 0.7;

	// // prm.solver.K = 4;
	// // prm.solver.M = 35;
	// // prm.solver.ns_search = true;
	// prm.solver.tol = 1e-9;
	// prm.solver.maxiter = 100;//1000;
	// // prm.solver.replacement = true;
	// // prm.solver.smoothing = true;
	
	
	// // bool reorder = true;
	
    // // if (reorder) {
        // // amgcl::adapter::reorder<> perm(A);

        // // Solver solve(world, perm(A), prm);
		// // int iters; double error;
		// // std::tie(iters, error) = solve(*A, perm(B0_vals), resiVar0);
		// // std::tie(iters, error) = solve(*A, perm(B1_vals), resiVar1);
		// // std::tie(iters, error) = solve(*A, perm(B2_vals), resiVar2);

        // // perm.inverse(resiVar0, resiVar0);
        // // perm.inverse(resiVar1, resiVar1);
        // // perm.inverse(resiVar2, resiVar2);
	// // }
	// // else{
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters0; double error0;
		// int iters1; double error1;
		// int iters2; double error2;
		// std::tie(iters0, error0) = solve(*A, B0_vals, resiVar0);
		// std::tie(iters1, error1) = solve(*A, B1_vals, resiVar1);
		// std::tie(iters2, error2) = solve(*A, B2_vals, resiVar2);
		// prof.toc("solve");
	// // }
	
	

	// // if (world.rank == 0) {
		// // cout << solve << endl;
		// // cout << prof << endl;
		// // cout << "|--- momentum : " << iters0+iters1+iters2 << " " << error0+error1+error2 << endl;
	// // }
	
	
// }


