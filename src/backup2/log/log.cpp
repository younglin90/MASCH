
#include "./log.h"

// 파일들 pretty하게 표시
string MASCH_FileSystem::showPrettySize(int fileSize){
	string logging = "";
	if(fileSize/1024 < 1){
		logging = to_string(fileSize) + " Byte";
	}
	else if(fileSize/1024 > 1 && fileSize/1024 <= 1024){
		logging = to_string(fileSize/1024) + " KB";
	}
	else if(fileSize/1024/1024 > 1 && fileSize/1024/1024 <= 1024){
		logging = to_string(fileSize/1024/1024) + " MB";
	}
	else{
		logging = to_string(fileSize/1024/1024/1024) + " GB";
	}
	return logging;
}

// 폴더 안에 있는 파일들 크기
int MASCH_FileSystem::calcFileSize(string folder){
	int fileSize = 0;
	for (const filesystem::directory_entry& entry :
		filesystem::recursive_directory_iterator(filesystem::current_path() / folder)) {
		// std::cout << entry.path() << std::endl;
		if(filesystem::is_regular_file(entry.path())){
			fileSize += std::filesystem::file_size(entry.path());
		}
	}
	return fileSize;
}





MASCH_TimeSystem& MASCH_TimeSystem::push(string inp_name){
	name.push(inp_name);
	start.push(chrono::system_clock::now());
	level.push(levelNow++);
	return (*this);
}
MASCH_TimeSystem& MASCH_TimeSystem::pop(){
	chrono::microseconds calcTime = 
	chrono::duration_cast<chrono::microseconds>(
	chrono::system_clock::now() - start.top()); 
	start.pop();
	std::stringstream ss;
	ss << setw(level.top()*2) << setfill('_') << 
	" " << name.top() << " : " << 
	calcTime.count() << " us";;
	(*this).logCalcTime.push_back(ss.str());
	name.pop();
	level.pop();
	return (*this);
}
void MASCH_TimeSystem::show(){
	if(MPI::COMM_WORLD.Get_rank()==0){
		for(auto& item : (*this).logCalcTime){
			cout << "| " << item << endl;
		}
	}
	(*this).logCalcTime.clear();
}

chrono::microseconds MASCH_TimeSystem::showCalcTime(){
	chrono::microseconds calcTime = 
	chrono::duration_cast<chrono::microseconds>(
	chrono::system_clock::now() - start.top()); 
	name.pop();
	start.pop();
	return calcTime;
}

string MASCH_TimeSystem::now(){
	auto now = std::chrono::system_clock::now();
	auto in_time_t = std::chrono::system_clock::to_time_t(now);
	std::stringstream ss;
	ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
	return ss.str();
}





void MASCH_Error::stop(string inp,string inp0, string inp1, long inp2){
	if(MPI::COMM_WORLD.Get_rank()==0){
		std::cout << "| #Error | " << 
		inp << ", file(" << inp0 << "), func(" << inp1 << "), line(" << inp2 << ")" << endl;;
	}
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	MPI::Finalize();
	exit(0);
}





void MASCH_Warning::push(string inp,string inp0, string inp1, long inp2){
	std::stringstream ss;
	ss << inp << ", file(" << inp0 << "), func(" << inp1 << "), line(" << inp2 << ")";
	warn_.push_back(ss.str());
}
void MASCH_Warning::pop(){
	warn_.pop_back();
}
void MASCH_Warning::show(){
	MASCH_MPI mpi;
	if(MPI::COMM_WORLD.Get_size()>1){
		vector<string> warn_glo = mpi.gatherv(warn_);
		if(MPI::COMM_WORLD.Get_rank()==0) {
			std::cout << "| #Warning |" << endl;
			int ip=0;
			for(auto& item : warn_glo){
				if(!item.empty())
					std::cout << "| proc(" << ip++ << ") : " << item << endl;
			}
		}
	}
	else{
		std::cout << "| #Warning : " << std::endl;
		for(auto& item : warn_){
			std::cout << item << std::endl;
		}
	}
	warn_.clear();
}




MASCH_Log& MASCH_Log::setLevel(int inp_level) {
	level = inp_level;
	return (*this);
}