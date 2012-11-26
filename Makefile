
DEBUG = 0
COMPILER = gcc
	
BOOST_DIR = $(HOME)/local/src/boost_1_51_0
BOOST_INC_DIR = $(BOOST_DIR)
PYTHON_DIR = $(HOME)/local
PYTHON_INC_DIRS = $(HOME)/local/include/python2.7 $(HOME)/local/lib/python2.7/site-packages/numpy/core/include
PYUBLAS_INC_DIR = .
INCLUDES = -I$(BOOST_INC_DIR) -I$(PYUBLAS_INC_DIR) \
	$(foreach dir, $(PYTHON_INC_DIRS), -I$(dir))
PYTHON_LIB = python27
BOOST_LIB_DIR = $(BOOST_DIR)/lib
ifeq ($(COMPILER), gcc)	
	BOOST_PYTHON_LIB = boost_python
	LIB_PATHS = $(HOME)/local/lib
	LIBS = -l$(BOOST_PYTHON_LIB) -l$(PYTHON_LIB)
	LIB_DIRS = $(foreach dir, $(LIB_PATHS), -L$(dir))
else	
	BOOST_PYTHON_LIB = boost_python-vc100-mt-1_49
	LIB_PATHS = $(BOOST_LIB_DIR) $(PYTHON_LIB_DIRS)
	LIB_NAMES = $(BOOST_PYTHON_LIB) $(PYTHON_LIB)
	LIB_DIRS = $(foreach dir, $(LIB_PATHS), /LIBPATH:$(dir))
	LIBS = tbb.lib $(BOOST_LIB_DIR)/$(BOOST_PYTHON_LIB).lib $(PYTHON_LIB).lib
endif

BUILD_DIR = build

#############################
# compile/link flags
#############################

ifeq ($(COMPILER), gcc)			
	CC = gcc
	CXX = g++ 
	LINK = g++	
	LINKFLAGS = -shared $(LINKFLAGS_EXE)
	ifeq ($(DEBUG), 0)
		CXXFLAGS = -O2 -std=c++11
	else
		CXXFLAGS = -g -std=c++11
	endif
else ifeq ($(COMPILER), msvc)
	CC = cl.exe
	CXX = cl.exe
#--compiler-options /MD
	LINK = link.exe
	ifeq ($(DEBUG), 0)
		CXXFLAGS = /nologo /O2 /MD /W3 /GS- /Zi /EHsc /Zm1000 /DNDEBUG /WL
		LINKFLAGS_EXE = /nologo /INCREMENTAL:NO /DEBUG
	else
		CXXFLAGS = /nologo /Od /MD /W3 /GS- /Zi /EHsc /Zm1000 /WL 
		LINKFLAGS_EXE = /nologo /INCREMENTAL:NO /DEBUG
	endif
	LINKFLAGS = /DLL $(LINKFLAGS_EXE)
	ifeq ($(FP_STRICT), 1)
		CXXFLAGS += /fp:strict /DFP_STRICT
	endif
endif #msvc

#########################
## generic rules
#########################

ifeq ($(COMPILER), gcc)

$(BUILD_DIR)/%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

_%.pyd:	$(BUILD_DIR)/%.o
	$(LINK) -o $@ $< $(LINKFLAGS) $(LIB_DIRS) $(LIBS)

%:	$(BUILD_DIR)/%.o
	$(LINK) $(LINKFLAGS_EXE) $(LIB_DIRS) $(LIBS) $< -o $@
	
else ifeq ($(COMPILER), msvc)

$(BUILD_DIR)/%.obj: %.cpp %.h
	$(CXX) /c $(CXXFLAGS) $(INCLUDES) /Tp$< /Fo$@
	
$(BUILD_DIR)/%.obj: %.cpp
	$(CXX) /c $(CXXFLAGS) $(INCLUDES) /Tp$< /Fo$@

_%.pyd:	$(BUILD_DIR)/%.obj _myfuncs.lib
	$(LINK) $(LINKFLAGS) $(LIB_DIRS) $(LIBS) _myfuncs.lib $< /OUT:$@ IMPLIB:$*.lib MANIFESTFILE:$*.pyd.manifest

%.exe: $(BUILD_DIR)/%.obj
	$(LINK) $(LINKFLAGS_EXE) $(LIB_DIRS) $(LIBS) $< /OUT:$@

endif

#########################
## targets
#########################
	
all:	linterp_test _linterp_python.pyd

