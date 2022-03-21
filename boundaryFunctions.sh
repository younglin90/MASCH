cd ./setting/
mpiicpc -fPIC -c boundaryFunctions.cpp 
mpiicpc -shared -o boundaryFunctions.so boundaryFunctions.o
cd ..
