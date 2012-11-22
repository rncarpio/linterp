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
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>

namespace bpl = boost::python;

#define bpl_assert(exp, s)	\
	if (!(exp)) {				\
      PyErr_SetString(PyExc_ValueError, (s));	\
      boost::python::throw_error_already_set();	\
	}	

typedef pyublas::numpy_vector<double> dPyArr;
typedef vector<dPyArr> dPyArrVector;

// Boost::Python wrapper class for NDInterpolator
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
};

typedef NDInterpolator_wrap<InterpSimplex<1,double,false,true,dPyArr,dPyArr> > NDInterpolator_1_S_wrap;
typedef NDInterpolator_wrap<InterpSimplex<2,double,false,true,dPyArr,dPyArr> > NDInterpolator_2_S_wrap;
typedef NDInterpolator_wrap<InterpSimplex<3,double,false,true,dPyArr,dPyArr> > NDInterpolator_3_S_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<1,double,false,true,dPyArr,dPyArr> > NDInterpolator_1_ML_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<2,double,false,true,dPyArr,dPyArr> > NDInterpolator_2_ML_wrap;
typedef NDInterpolator_wrap<InterpMultilinear<3,double,false,true,dPyArr,dPyArr> > NDInterpolator_3_ML_wrap;

BOOST_PYTHON_MODULE(_linterp_python)
{   
  bpl::class_<NDInterpolator_1_S_wrap, boost::noncopyable>("Interp_1_S", bpl::init<vector<dPyArr>, dPyArr>())
		.def("interp_vec", &NDInterpolator_1_ML_wrap::interp_vec) 
	;
}