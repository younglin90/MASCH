#include "mpi.h"
#include <vector>
// #include "mpi.h"
using namespace std;

void SEMO_MPI_Builder::exchangeDatas(
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendValues, vector<int>& recvValues){

	int size = sendCounts.size(); 
	
// cout << size<<endl;
	vector<int> sendDisps(size,0);
	vector<int> recvDisps(size,0);
		
	for(int i=1; i<size; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<size; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
}


void SEMO_MPI_Builder::exchangeDatas(
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<double>& sendValues, vector<double>& recvValues){
	int size = sendCounts.size(); 
	
	vector<int> sendDisps(size,0);
	vector<int> recvDisps(size,0);
		
	for(int i=1; i<size; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<size; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0.0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
	
}


void SEMO_MPI_Builder::setCellDatasToFaceRight(
	SEMO_Mesh_Builder& mesh, 
	int cin, int fin,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps
	){
		
	vector<double> sendValues;
	vector<double> recvValues;
		
	int size = sendCounts.size(); 
	
	sendValues.clear();
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			sendValues.push_back(mesh.cells[face.owner].var[cin]);
		}
	}
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0.0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	vector<double>::iterator iter = recvValues.begin();
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.varR[fin] = *iter;
			++iter;
		}
	}
	
}



void SEMO_MPI_Builder::setProcsFaceDatas(
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


// void SEMO_MPI_Builder::setProcsFaceDatas(
	// vector<bool>& sendValues, vector<bool>& recvValues,
	// vector<int>& sendCounts, vector<int>& recvCounts, 
	// vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	// int size = sendCounts.size(); 
	
	// recvValues.clear();
	// recvValues.resize(recvDisps[size-1] + recvCounts[size-1],false);
	
	// // MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_C_BOOL, 
				   // // recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_C_BOOL, 
				   // // MPI_COMM_WORLD);
	
// }


void SEMO_MPI_Builder::setProcsFaceDatas(
	vector<double>& sendValues, vector<double>& recvValues,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	int size = sendCounts.size(); 
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0.0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
	
}




// void SEMO_MPI_Builder::setProcsFaceDatasDouble(
	// vector<double>& sendValues, vector<double>& recvValues,
	// vector<int>& sendCounts, vector<int>& recvCounts, 
	// vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	// int size = sendCounts.size(); 
	
	// recvValues.clear();
	// recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	// MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   // recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   // MPI_COMM_WORLD);
	
// }


void SEMO_MPI_Builder::sendRecvTemporaryData(
	SEMO_Mesh_Builder& mesh,
	vector<double>& sendValues, vector<double>& recvValues){
		
		
	recvValues.clear();
		
	vector<double> recv;
	vector<double> send;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			send.push_back(sendValues[face.owner]);
		}
	}
	setProcsFaceDatas(
				send, recv,
				mesh.countsProcFaces, mesh.countsProcFaces, 
				mesh.displsProcFaces, mesh.displsProcFaces);
	int num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			recvValues.push_back(recv[num]);
			
			++num;
		}
	}
	
}



void SEMO_MPI_Builder::sendRecvTemporaryVectorData(
	SEMO_Mesh_Builder& mesh,
	vector<vector<double>>& sendValues, vector<vector<double>>& recvValues){
		
	recvValues.clear();
		
	vector<double> x_recv, y_recv, z_recv;
	vector<double> x_send, y_send, z_send;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			x_send.push_back(sendValues[face.owner][0]);
			y_send.push_back(sendValues[face.owner][1]);
			z_send.push_back(sendValues[face.owner][2]);
		}
	}
	setProcsFaceDatas(
				x_send, x_recv,
				mesh.countsProcFaces, mesh.countsProcFaces, 
				mesh.displsProcFaces, mesh.displsProcFaces);
	setProcsFaceDatas(
				y_send, y_recv,
				mesh.countsProcFaces, mesh.countsProcFaces, 
				mesh.displsProcFaces, mesh.displsProcFaces);
	setProcsFaceDatas(
				z_send, z_recv,
				mesh.countsProcFaces, mesh.countsProcFaces, 
				mesh.displsProcFaces, mesh.displsProcFaces);
	int num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			vector<double> tmp(3,0.0);
			tmp[0] = x_recv[num];
			tmp[1] = y_recv[num];
			tmp[2] = z_recv[num];
			recvValues.push_back(tmp);
			
			++num;
		}
	}
		
	
}



void SEMO_MPI_Builder::sendRecvTemporaryCellData(
	SEMO_Mesh_Builder& mesh,
	int& cnum,
	vector<double>& recvValues){
		
	recvValues.clear();
	// recvValues.resize(cnum.size(),vector<double>());
		
	vector<double> recv;
	vector<double> send;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			send.push_back(mesh.cells[face.owner].var[cnum]);
		}
	}
	setProcsFaceDatas(
				send, recv,
				mesh.countsProcFaces, mesh.countsProcFaces, 
				mesh.displsProcFaces, mesh.displsProcFaces);
	int num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			recvValues.push_back(recv[num]);
			
			++num;
		}
	}
		
	
}



void SEMO_MPI_Builder::sendRecvTemporaryCellVectorData(
	SEMO_Mesh_Builder& mesh,
	vector<int>& cnum,
	vector<vector<double>>& recvValues){
		
	recvValues.clear();
	// recvValues.resize(cnum.size(),vector<double>());
		
	vector<vector<double>> recv(cnum.size(),vector<double>());
	vector<vector<double>> send(cnum.size(),vector<double>());
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			for(int ii=0; ii<cnum.size(); ++ii){
				send[ii].push_back(mesh.cells[face.owner].var[cnum[ii]]);
			}
		}
	}
	for(int ii=0; ii<cnum.size(); ++ii){
		setProcsFaceDatas(
					send[ii], recv[ii],
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	}
	int num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			vector<double> tmp;
			for(int ii=0; ii<cnum.size(); ++ii){
				tmp.push_back(recv[ii][num]);
			}
			recvValues.push_back(tmp);
			
			++num;
		}
	}
		
	
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