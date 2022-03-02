
#include "./mpi.h"


void MASCH_MPI::procFace_Alltoallv(
	vector<int>& sendValues, vector<int>& recvValues,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	int size = sendCounts.size(); 
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
}

void MASCH_MPI::procFace_Alltoallv(
	vector<double>& sendValues, vector<double>& recvValues,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	int size = sendCounts.size(); 
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
	
}



void MASCH_MPI::Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	recv_value.clear();
	displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<double> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	recv_value.resize(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, 
				   recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	displs.resize(size+1,0);
	displs = recv_displs;
	


}



void MASCH_MPI::Alltoallv(vector<vector<int>>& inp_send_value, vector<int>& recv_value, vector<int>& displs){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<int> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	recv_value.resize(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_INT, 
				   recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
	displs = recv_displs;
	


}




void MASCH_MPI::Alltoallv(vector<vector<double>>& inp_send_value, vector<vector<double>>& recv_value){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<double> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	vector<double> tmp_recv_value(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, 
				   tmp_recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	recv_value.clear();
	recv_value.resize(size);
	for(int ip=0; ip<size; ++ip){
		int str = recv_displs[ip];
		int end = recv_displs[ip+1];
		for(int i=str; i<end; ++i){
			recv_value[ip].push_back(tmp_recv_value[i]);
		}
	}


}



void MASCH_MPI::Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value){

	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	// displs.clear();
	
	vector<int> send_counts(size,0);
	for(int ip=0; ip<size; ++ip){
		send_counts[ip] = inp_send_value[ip].size();
	}
	vector<int> send_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) send_displs[ip+1] = send_displs[ip] + send_counts[ip];
	
	vector<int> recv_counts(size,0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	
	vector<int> recv_displs(size+1,0);
	for(int ip=0; ip<size; ++ip) recv_displs[ip+1] = recv_displs[ip] + recv_counts[ip];
	
	
	vector<int> send_value;
	for(int ip=0; ip<size; ++ip){
		for(auto& item : inp_send_value[ip]){
			send_value.push_back(item);
		}
	}
	
	vector<int> tmp_recv_value(recv_displs[size]);
	MPI_Alltoallv( send_value.data(), send_counts.data(), send_displs.data(), MPI_INT, 
				   tmp_recv_value.data(), recv_counts.data(), recv_displs.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
	recv_value.clear();
	recv_value.resize(size);
	for(int ip=0; ip<size; ++ip){
		int str = recv_displs[ip];
		int end = recv_displs[ip+1];
		for(int i=str; i<end; ++i){
			recv_value[ip].push_back(tmp_recv_value[i]);
		}
	}


}




void MASCH_MPI::Gatherv(vector<int>& my_value, vector<int>& value, vector<int>& displs){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	int my_value_size = my_value.size();
	if(rank==0){
		vector<int> counts(size,0);
        MPI_Gather(&my_value_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		displs.clear(); displs.resize(size+1,0);
		for(int ip=0; ip<size; ++ip) displs[ip+1] = displs[ip]+counts[ip];
            
        value.resize(displs[size]);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_INT, value.data(), counts.data(), displs.data(), 
		MPI_INT, 0, MPI_COMM_WORLD);
	}
	else{
        MPI_Gather(&my_value_size, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_INT, NULL, NULL, NULL, 
		MPI_INT, 0, MPI_COMM_WORLD);
	}
	
}



void MASCH_MPI::Gatherv(vector<double>& my_value, vector<double>& value, vector<int>& displs){
	
	int rank = static_cast<int>(MPI::COMM_WORLD.Get_rank()); 
	int size = static_cast<int>(MPI::COMM_WORLD.Get_size()); 
	
	
	int my_value_size = my_value.size();
	if(rank==0){
		vector<int> counts(size,0);
        MPI_Gather(&my_value_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		displs.clear(); displs.resize(size+1,0);
		for(int ip=0; ip<size; ++ip) displs[ip+1] = displs[ip]+counts[ip];
            
        value.resize(displs[size]);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_DOUBLE, value.data(), counts.data(), displs.data(), 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else{
        MPI_Gather(&my_value_size, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(my_value.data(), my_value_size, MPI_DOUBLE, NULL, NULL, NULL, 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	
}




// #include "../../solvers/solvers.h"  

// wchar_t wStr[] = L"€áa¢cée£";
// int iStr = sizeof(wStr) / sizeof(wStr[0]);        // length of the string
// wchar_t *pStr = 0;

// using namespace chrono; // std 내에 chrono 가 존재  

// #define MASCH_FFL __FILE__,__FUNCTION__,__LINE__

// class MASCH_MPI {
// private:
// public:
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	// MASCH_MPI(){}
	// // ~MASCH_MPI(){ MPI::Finalize(); }
	
	// vector<string> gatherv(vector<string>& inp){
		// if(size>1){
			// vector<int> counts(size,0);
			// vector<int> disp(size+1,0);
			
			// std::vector<char> cstrings;
			// for(std::string s: inp)
			// {
				// for(int i = 0; i < strlen(s.c_str()); ++i)
				// {
					// cstrings.push_back(s.c_str()[i]);
				// }
			// }
			
			// int my_count = cstrings.size();
			// MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
			// for(int i=1; i<size+1; ++i) disp[i] = disp[i-1] + counts[i-1];
			
			// std::vector<char> cstrings_recv(disp[size]);
			// if(rank==0){
				// MPI_Gatherv(
				// cstrings.data(), cstrings.size(), MPI_CHAR, 
				// cstrings_recv.data(), counts.data(), disp.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
				
			// }
			// else{
				// MPI_Gatherv(
				// cstrings.data(), cstrings.size(), MPI_CHAR, 
				// NULL, NULL, NULL, MPI_CHAR, 0, MPI_COMM_WORLD);
			// }
			
			// vector<string> results;
			// if(rank==0){
				// for(int ip=0; ip<size; ++ip){
					// std::stringstream ss;
					// for(int i=disp[ip]; i<disp[ip+1]; ++i){
						// ss << cstrings_recv[i];
					// }
					// results.push_back(ss.str());
				// }
			// }
			// return results;
		// }
	// }

// };