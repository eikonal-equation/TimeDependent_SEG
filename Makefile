ICXX = icpc
CXXFLAGS = -std=c++14 -O3
BOOSTPATH = /usr/local/include/boost
EIGENPATH = /usr/local/include/eigen3/Eigen
INCLUDE = -I $(BOOSTPATH) -I $(EIGENPATH) -I ./include/ -I ./external/
LIB = -L $(BOOSTPATH) -L $(EIGENPATH)
SRC = ./src/
TEST =

CPPFILES= $(SRC)main.cpp \
							$(SRC)CTimeGrid.cpp \
							$(SRC)CFMM.cpp \
							$(SRC)CTimeDependentHjbSolver.cpp \
							$(SRC)CTimeDependentTracer.cpp \
							$(SRC)CAdversarialPlan.cpp\
							$(SRC)CTimeDependentPlanner.cpp\
							$(SRC)CMovingObserver.cpp \
							$(SRC)CStationaryTerrain.cpp \
							$(SRC)CStationaryObstacle.cpp \
							$(SRC)OptimizationAlgorithm.cpp \
							$(SRC)CTerrain.cpp \
							./external/QuadProg++.cc \

main: $(CPPFILES)
	 $(CXX)  $(CXXFLAGS) $(INCLUDE) $(CPPFILES)  $(LIB) -o $@

run: main
		./main $(TEST)

debug: $(CPPFILES)
	g++  -o1 -g -std=c++14 $(INCLUDE) $^  $(LIB)
	gdb -ex=r --args 1 $(TEST)

memleak: $(CPPFILES)
	g++ -g -o0 -std=c++14  $(INCLUDE) $^  $(LIB)

profile: $(CPPFILES)
	icpc $(CXXFLAGS) -ofast -profile-functions -profile-loops=all -profile-loops-report=2 $(INCLUDE) $^  $(LIB)

clean:
	rm *.o 0 main
	rm loop_prof*
	rm -rf 0.dSYM

realclean:
	rm output/ex*
	clean
