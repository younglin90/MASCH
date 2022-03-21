cd ./setting/
mpiicpc -fPIC -c lib_initial.cpp 
mpiicpc -shared -o lib_initial.so lib_initial.o
cd ..
