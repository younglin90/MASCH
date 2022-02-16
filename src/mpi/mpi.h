#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <list>
#include <mpi.h>
#include "mpi.h"
using namespace std;

#include "../mesh/mesh.h"


class SEMO_MPI_Builder{
	public:
		SEMO_MPI_Builder() {
		}
		void init(){
			MPI::Init(); 
		}
		void setRank(){
			rank = MPI::COMM_WORLD.Get_rank(); 
		}
		void setSize(){
			size = MPI::COMM_WORLD.Get_size();
		}
		int getRank() {
			return rank;
		}
		int getSize() {
			return size;
		}
		void setSendCounts(int in){
			sendCounts = new int[in];
		}
		void setRecvCounts(int in){
			recvCounts = new int[in];
		}
		void setSendDispls(int in){
			sendDispls = new int[in];
		}
		void setRecvDispls(int in){
			recvDispls = new int[in];
		}
		
		void exchangeDatas(
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendValues, vector<int>& recvValues);
		void exchangeDatas(
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<double>& sendValues, vector<double>& recvValues);
		void setCellDatasToFaceRight(
				SEMO_Mesh_Builder& mesh, 
				int cin, int fin,
				vector<int>& sendCounts, vector<int>& recvCounts, 
				vector<int>& sendDisps, vector<int>& recvDisps
				);
		
		
		void setProcsFaceDatas(
			vector<int>& sendValues, vector<int>& recvValues,
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendDisps, vector<int>& recvDisps);
		
		void setProcsFaceDatas(
			vector<double>& sendValues, vector<double>& recvValues,
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendDisps, vector<int>& recvDisps);
			

		void sendRecvTemporaryData(
			SEMO_Mesh_Builder& mesh,
			vector<double>& sendValues, vector<double>& recvValues);

		void sendRecvTemporaryVectorData(
			SEMO_Mesh_Builder& mesh,
			vector<vector<double>>& sendValues, vector<vector<double>>& recvValues);
			
		void sendRecvTemporaryCellData(
			SEMO_Mesh_Builder& mesh,
			int& cnum,
			vector<double>& recvValues);
		
		void sendRecvTemporaryCellVectorData(
			SEMO_Mesh_Builder& mesh,
			vector<int>& cnum,
			vector<vector<double>>& recvValues);
	
		// void setProcsFaceDatas(
			// vector<bool>& sendValues, vector<bool>& recvValues,
			// vector<int>& sendCounts, vector<int>& recvCounts, 
			// vector<int>& sendDisps, vector<int>& recvDisps);
		
		
		// void setProcsFaceDatasDouble(
			// vector<double>& sendValues, vector<double>& recvValues,
			// vector<int>& sendCounts, vector<int>& recvCounts, 
			// vector<int>& sendDisps, vector<int>& recvDisps);
		
		int* sendCounts;
		int* recvCounts;
		int* sendDispls;
		int* recvDispls;
		int nSend;
		int nRecv;
		
	private:
		int rank;
		int size;
		MPI::Status status;
		
};



#include <cstring>

// namespace MASCH_MPI{
class MASCH_MPI{
public:
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	void Alltoallv(vector<vector<double>>& inp_send_value, vector<double>& recv_value, vector<int>& displs);
	void Alltoallv(vector<vector<int>>& inp_send_value, vector<int>& recv_value, vector<int>& displs);
	void Alltoallv(vector<vector<double>>& inp_send_value, vector<vector<double>>& recv_value);
	void Alltoallv(vector<vector<int>>& inp_send_value, vector<vector<int>>& recv_value);
	void Gatherv(vector<int>& my_value, vector<int>& value, vector<int>& displs);
	void Gatherv(vector<double>& my_value, vector<double>& value, vector<int>& displs);
	MASCH_MPI(){}
	// ~MASCH_MPI(){ MPI::Finalize(); }
	
	vector<string> gatherv(vector<string>& inp){
		if(size>1){
			vector<int> counts(size,0);
			vector<int> disp(size+1,0);
			
			std::vector<char> cstrings;
			for(std::string s: inp)
			{
				for(int i = 0; i < strlen(s.c_str()); ++i)
				{
					cstrings.push_back(s.c_str()[i]);
				}
			}
			
			int my_count = cstrings.size();
			MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
			for(int i=1; i<size+1; ++i) disp[i] = disp[i-1] + counts[i-1];
			
			std::vector<char> cstrings_recv(disp[size]);
			if(rank==0){
				MPI_Gatherv(
				cstrings.data(), cstrings.size(), MPI_CHAR, 
				cstrings_recv.data(), counts.data(), disp.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
				
			}
			else{
				MPI_Gatherv(
				cstrings.data(), cstrings.size(), MPI_CHAR, 
				NULL, NULL, NULL, MPI_CHAR, 0, MPI_COMM_WORLD);
			}
			
			vector<string> results;
			if(rank==0){
				for(int ip=0; ip<size; ++ip){
					std::stringstream ss;
					for(int i=disp[ip]; i<disp[ip+1]; ++i){
						ss << cstrings_recv[i];
					}
					results.push_back(ss.str());
				}
			}
			return results;
		}
	}
};

