SimplePhysics: *.cpp
	g++ *.cpp -fopenmp -std=c++17 -O3 -o SimplePhysics $(shell pkg-config --libs raylib)
#	g++-13 *.cpp -fopenmp -std=c++17 -O3 -o SimplePhysics $(shell pkg-config --libs raylib)
# g++-13 *.cpp -std=c++17 -g -o SimplePhysics $(shell pkg-config --libs raylib)
# clang++ *.cpp -std=c++17 -g -o SimplePhysics $(shell pkg-config --libs raylib)
