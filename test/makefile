# FIXME: do it right.
# TODO: implement CMake

INC_DIR="../src/"
CXXFLAGS=-std=c++11 -g -Wall -E -I $(INC_DIR) -v


all: units.cc
	g++ $(CXXFLAGS) -o ../bin/units units.cc fevor_crystal.hh fevor_crystal.cc vector_tensor_opperations.hh vector_tensor_opperations.cc 
