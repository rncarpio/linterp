//
// Copyright (c) 2012 Ronaldo Carpio
//                                     
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and   
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.          
// It is provided "as is" without express or implied warranty.
//                                                            

#include "linterp.h"
#include <Python.h>
#include <numpy/arrayobject.h>
#include <pyublas/numpy.hpp>

namespace bpl = boost::python;

#define bpl_assert(exp, s)	\
	if (!(exp)) {				\
      PyErr_SetString(PyExc_ValueError, (s));	\
      boost::python::throw_error_already_set();	\
	}	

typedef pyublas::numpy_vector<double> dPyArr;
typedef vector<dPyArr> dPyArrVector;

template <class Interpolator>
class NDInterpolator_wrap {
public:
  typedef Interpolator interp_type;  
  typedef typename interp_type::value_type T;
  static const int N = interp_type::m_N;
  
  std::unique_ptr<interp_type> m_p_interp_obj;
  
  NDInterpolator_wrap(vector<dPyArr> const &gridList, dPyArr f) {
    const double *grids_begin[N];
	int grids_len[N];
	for (int i=0; i<N; i++)		{ 
	  grids_begin[i] = &(gridList[i][0]); 
	  grids_len[i] = gridList[i].size();
	}
	m_p_interp_obj.reset(new interp_type(grids_begin, grids_len, &f[0], &f[0] + f.size()));
  }
  
  // each dimension is a separate array
  dPyArr interp_vec(vector<dVec> const &coord_arrays) const {
    bpl_assert(coord_arrays.size() == N, "wrong number of coords");
	int n = coord_arrays[0].size();
	for (int i=0; i<N; i++) {
	  bpl_assert(coord_arrays[i].size() == n, "coord arrays have different sizes");
    }	  
	dPyArr result(n);

    clock_t t1 = clock();	
    m_p_interp_obj->interp_vec(n, coord_arrays.begin(), coord_arrays.end(), result.begin());
    clock_t t2 = clock();
	//printf("interp: %d points, %d clocks, %f sec\n", n, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);	
	return result;
  }
  
  double interp(dPyArr x) const {
    bpl_assert(x.size() == N, "wrong number of coords");
	return m_p_interp_obj->interp(x.begin());
  }
};

typedef NDInterpolator_wrap<InterpSimplex<1,double,false,true,dPyArr,dPyArr> > NDInterpolator_1_S_wrap;
typedef NDInterpolator_wrap<InterpSimplex<2,double,false,true,dPyArr,dPyArr> > NDInterpolator_2_S_wrap;
typedef NDInterpolator_wrap<InterpSimplex<3,double,false,true,dPyArr,dPyArr> > NDInterpolator_3_S_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<1,double,false,true,dPyArr,dPyArr> > NDInterpolator_1_ML_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<2,double,false,true,dPyArr,dPyArr> > NDInterpolator_2_ML_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<3,double,false,true,dPyArr,dPyArr> > NDInterpolator_3_ML_wrap;

//////////////////////////////////////////////////////////////////////////////////////////////
// Python glue code

// vector of T -> Python list
template <class VecT>
struct vec_to_list {
  static PyObject* convert(VecT const &vec) {
    bpl::list result;
	typename VecT::const_iterator iter;
	for (iter=vec.begin(); iter != vec.end(); iter++) {
	  result.append(*iter);
	}
    return boost::python::incref(result.ptr());
  }
};

// Python list of numpy arrays -> vector<dVec>
struct dVecVec_from_python_list {
  static void *check(PyObject *pObj) {    
    if (!PyList_Check(pObj)) {
	  return NULL;
	} else {
	  bpl::handle<> h(bpl::borrowed(pObj));
	  bpl::list l(h);
	  for (int i=0; i<bpl::len(l); i++) {
	    bpl::object obj(l[i]);
	    if (!PyArray_Check(obj.ptr())) {
		  return NULL;
		}
	  }
	}
	return pObj;
  }
  static void construct(PyObject *pObj, boost::python::converter::rvalue_from_python_stage1_data* data) {
	bpl::handle<> h(bpl::borrowed(pObj));
	bpl::list l(h);
	void *storage = ((boost::python::converter::rvalue_from_python_storage<vector<dVec>>*)data)->storage.bytes;
    new (storage) vector<dVec>(bpl::len(l));
	vector<dVec> *pVec = static_cast<vector<dVec>*>(storage);
	data->convertible = storage;
	// copy elements
	for (int i=0; i<bpl::len(l); i++) {	  
	  dPyArr a = bpl::extract<dPyArr>(l[i]);
	  (*pVec)[i].resize(a.size());
	  std::copy(a.begin(), a.end(), (*pVec)[i].begin());	  
	}
  }    
};

// Python list of numpy arrays -> vector<dPyArr>
struct dPyArrVec_from_python_list {
  static void *check(PyObject *pObj) {    
    if (!PyList_Check(pObj)) {
	  return NULL;
	} else {
	  bpl::handle<> h(bpl::borrowed(pObj));
	  bpl::list l(h);
	  for (int i=0; i<bpl::len(l); i++) {
	    bpl::object obj(l[i]);
	    if (!PyArray_Check(obj.ptr())) {
		  return NULL;
		}
	  }
	}
	return pObj;
  }
  static void construct(PyObject *pObj, boost::python::converter::rvalue_from_python_stage1_data* data) {
	bpl::handle<> h(bpl::borrowed(pObj));
	bpl::list l(h);
	void *storage = ((boost::python::converter::rvalue_from_python_storage<vector<dVec>>*)data)->storage.bytes;
    new (storage) vector<dPyArr>(bpl::len(l));
	vector<dPyArr> *pVec = static_cast<vector<dPyArr>*>(storage);
	data->convertible = storage;
	// copy elements
	for (int i=0; i<bpl::len(l); i++) {
      (*pVec)[i] = bpl::extract<dPyArr>(l[i]);
	}
  }    
};

BOOST_PYTHON_MODULE(_linterp_python)
{   
  bpl::to_python_converter<dVec, vec_to_list<dVec>>();		// vec of doubles to list  
  bpl::to_python_converter<iVec, vec_to_list<iVec>>();		// vec of ints to list
  bpl::to_python_converter<dPyArrVector, vec_to_list<dPyArrVector>>();		// vec of pyublas arrays
  bpl::to_python_converter<array<double,2>, vec_to_list<array<double,2>>>();
  bpl::to_python_converter<array<double,3>, vec_to_list<array<double,3>>>();
  bpl::to_python_converter<array<double,4>, vec_to_list<array<double,4>>>();
  bpl::to_python_converter<array<double,5>, vec_to_list<array<double,5>>>();
  bpl::to_python_converter<array<double,6>, vec_to_list<array<double,6>>>();
  bpl::to_python_converter<array<int,2>, vec_to_list<array<int,2>>>();
  bpl::to_python_converter<array<int,3>, vec_to_list<array<int,3>>>();  
  bpl::to_python_converter<array<int,4>, vec_to_list<array<int,4>>>();
  bpl::to_python_converter<array<int,5>, vec_to_list<array<int,5>>>();
  bpl::to_python_converter<array<int,6>, vec_to_list<array<int,6>>>();
  bpl::to_python_converter<array<int,7>, vec_to_list<array<int,7>>>();
  
  bpl::converter::registry::push_back(&dVecVec_from_python_list::check, &dVecVec_from_python_list::construct, bpl::type_id<vector<dVec> >());
  bpl::converter::registry::push_back(&dPyArrVec_from_python_list::check, &dPyArrVec_from_python_list::construct, bpl::type_id<vector<dPyArr> >());

  
    bpl::class_<NDInterpolator_1_S_wrap, boost::noncopyable>("Interp_1_S", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_1_S_wrap::interp_vec)
		.def("interp", &NDInterpolator_1_S_wrap::interp_vec)
	;
    bpl::class_<NDInterpolator_2_S_wrap, boost::noncopyable>("Interp_2_S", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_2_S_wrap::interp_vec)
		.def("interp", &NDInterpolator_2_S_wrap::interp_vec)
	;
    bpl::class_<NDInterpolator_3_S_wrap, boost::noncopyable>("Interp_3_S", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_3_S_wrap::interp_vec)
		.def("interp", &NDInterpolator_3_S_wrap::interp_vec)
	;

    bpl::class_<NDInterpolator_1_ML_wrap, boost::noncopyable>("Interp_1_ML", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_1_ML_wrap::interp_vec)
		.def("interp", &NDInterpolator_1_ML_wrap::interp_vec)
	;
    bpl::class_<NDInterpolator_2_ML_wrap, boost::noncopyable>("Interp_2_ML", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_2_ML_wrap::interp_vec)
		.def("interp", &NDInterpolator_2_ML_wrap::interp_vec)
	;
    bpl::class_<NDInterpolator_3_ML_wrap, boost::noncopyable>("Interp_3_ML", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_3_ML_wrap::interp_vec)
		.def("interp", &NDInterpolator_3_ML_wrap::interp_vec)
	;
	
}


