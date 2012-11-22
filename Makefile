
PYTHON_DIR = c:/Python27
PYTHON_INC_DIRS = -I$(PYTHON_DIR)/include -I$(PYTHON_DIR)/lib/site-packages/numpy/core/include
PYTHON_LIB_DIRS = $(PYTHON_DIR)/libs 
BOOST_PYTHON_LIB = -LC:/boost/boost_1_49_0/lib -lboost_python-mgw46-mt-1_49

CXX = g++
CXXFLAGS = -O2 -std=c++11
INCLUDES = -IC:/boost/boost_1_49_0 -I. $(PYTHON_INC_DIRS)

linterp_test: linterp_test.o
	$(CXX) -o $@ $^
	
%.o:	%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $^

linterp_python.dll: linterp_python.o
	$(CXX) -shared -o $@ $(BOOST_PYTHON_LIB) $^
	
all:	linterp_test
#linterp_python.dll