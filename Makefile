
CXX = g++
CXXFLAGS = -O3 -std=c++11
INCLUDES = -IC:/boost/boost_1_49_0

linterp_test:	linterp_test.cpp
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDES) $^

all:	linterp_test