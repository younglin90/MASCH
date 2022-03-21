cd ./setting/
mpiicpc -fPIC -c lib_boundary.cpp 
mpiicpc -shared -o lib_boundary.so lib_boundary.o
cd ..
