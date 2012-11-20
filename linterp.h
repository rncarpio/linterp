
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

/*
This is a C++ header-only library for N-dimensional linear interpolation on a rectangular grid. Implements two methods:
* Multilinear: Interpolate using the N-dimensional hypercube containing the point. Interpolation step is O(2^N) 
* Simplicial: Interpolate using the N-dimensional simplex containing the point. Interpolation step is O(N log N), but less accurate.
Requires boost/multi_array library.

For a description of the algorithms, see:
* Weiser & Zarantonello (1988), "A Note on Piecewise Linear and Multilinear Table Interpolation in Many Dimensions", _Mathematics of Computation_ 50 (181), p. 189-196
* Davies (1996), "Multidimensional Triangulation and Interpolation for Reinforcement Learning", _Proceedings of Neural Information Processing Systems 1996_
*/

#ifndef _linterp_h
#define _linterp_h

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <cstdarg>
#include <string>
#include <vector>
#include <array>
#include <functional>

#include <boost/multi_array.hpp>

using std::vector;
using std::array;
typedef unsigned int uint;
typedef vector<int> iVec;
typedef vector<double> dVec;


// TODO:
//  - specify behavior past grid boundaries.
//    1) clamp
//    2) return a pre-determined value (e.g. NaN)

// compile-time params:
//   1) number of dimensions
//   2) scalar type T
//   3) copy data or not (default: false). The grids will always be copied
//   4) ref count class (default: none)
//   5) continuous or not. If not continuous, each dimension will have twice as many data points

// run-time constructor params:
//   1) f
//   2) grids
//   3) behavior outside grid: default=clamp
//   4) value to return outside grid, defaut=nan

struct EmptyClass {};

template <int N, class T, bool CopyData = true, class RefCountT = EmptyClass, bool Continuous = true>
class NDInterpolator {
public:
  typedef T value_type;
  typedef RefCountT array_ref_count_type;
  
  static const int m_N = N;
  static const bool m_bCopyData = CopyData;
  static const bool m_bContinuous = Continuous;
    
  typedef vector<T> grid_type;  
  typedef boost::const_multi_array_ref<T, N> array_type; 
  typedef std::unique_ptr<array_type> array_type_ptr;
  
  array_type_ptr m_pF;
  RefCountT m_refF;							// reference count for m_pF
  vector<T> m_FCopy;
  
  array<vector<T>, N> m_gridList;
  
  // non ref-counted constructor.
  template <class IterT1, class IterT2>  														
  NDInterpolator(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end) {	
    init(begins_begin, begins_end, ends_begin, f_begin, f_end);  
  }
  // ref-counted constructor
  template <class IterT1, class IterT2>
  NDInterpolator(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end, RefCountT &refF) {
    init_refcount(begins_begin, begins_end, ends_begin, f_begin, f_end, refF);
  }	
  
  template <class IterT1, class IterT2>  														
  void init(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end) {
    set_grids(begins_begin, begins_end, ends_begin);
	set_f_array(f_begin, f_end, m_bCopyData);
  }  
  template <class IterT1, class IterT2>
  void init_refcount(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end, RefCountT &refF) {	
    set_grids(begins_begin, begins_end, ends_begin);
	set_f_array_ref(f_begin, f_end, refF);    	
  }	
  
  template <class IterT1>  
  void set_grids(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin) {
    int nGrids = begins_end - begins_begin;
    assert(nGrids == N);
	for (int i=0; i<N; i++) {
	  int gridLength = ends_begin[i] - begins_begin[i];
	  m_gridList[i].assign(begins_begin[i], ends_begin[i]);
	}
  }

  // assumes that [f_begin, f_end) is a contiguous array in C-order  
  template <class IterT>  
  void set_f_array(IterT f_begin, IterT f_end, bool bCopy) {
    unsigned int nGridPoints = 1;
	array<int,N> sizes;
	for (unsigned int i=0; i<m_gridList.size(); i++) {
	  sizes[i] = m_gridList[i].size();
	  nGridPoints *= sizes[i];
	}

	int f_len = f_end - f_begin;
	if ( (m_bContinuous && f_len != nGridPoints) || (!m_bContinuous && f_len != 2 * nGridPoints) ) {
	  throw std::invalid_argument("f has wrong size");
	}
	for (unsigned int i=0; i<m_gridList.size(); i++) {
	  if (!m_bContinuous) { sizes[i] *= 2; }	  
	}

    if (bCopy == false) {
	  m_pF.reset(new array_type(f_begin, sizes));
	} else {
	  m_FCopy = vector<T>(f_begin, f_end);
	  m_pF.reset(new array_type(&m_FCopy[0], sizes));
	}
  }
  template <class IterT>  
  void set_f_array_ref(IterT f_begin, IterT f_end, RefCountT &refF, bool bCopy=CopyData) {
    set_f_array(f_begin, f_end, bCopy);
    m_refF = refF;
  }
  
  // -1 is before the first grid point
  // N-1 (where grid.size() == N) is after the last grid point
  int find_cell(int dim, T x) const {  
    grid_type const &grid(m_gridList[dim]);
    if (x < *(grid.begin())) return -1;
    else if (x >= *(grid.end()-1)) return grid.size()-1;
    else {
      auto i_upper = std::upper_bound(grid.begin(), grid.end(), x);
	  return i_upper - grid.begin() - 1;
    }	
  }
  
  // return the value of f at the given cell and vertex
  T get_f_val(array<int,N> const &cell_index, array<int,N> const &v_index) const {
    array<int,N> f_index;
	
	if (m_bContinuous) {	  
	  for (int i=0; i<N; i++) {
	    if (cell_index[i] < 0) {
		  f_index[i] = 0;		  
		} else if (cell_index[i] >= m_gridList[i].size()-1) {
		  f_index[i] = m_gridList[i].size()-1;		  
		} else {
		  f_index[i] = cell_index[i] + v_index[i];		  
		}
	  }
	} else {
	  for (int i=0; i<N; i++) {
	    if (cell_index[i] < 0) {
		  f_index[i] = 0;
		} else if (cell_index[i] >= m_gridList[i].size()-1) {
		  f_index[i] = (2*m_gridList[i].size())-1;
		} else {
		  f_index[i] = 1 + (2*cell_index[i]) + v_index[i];
		}
	  }
	}
    return (*m_pF)(f_index);
  }
  
  T get_f_val(array<int,N> const &cell_index, int v) const {
    array<int,N> v_index;
    for (int dim=0; dim<N; dim++) {
	  v_index[dim] = (v >> (N-dim-1)) & 1;						// test if the i-th bit is set
    }
    return get_f_val(cell_index, v_index);
  }	
};

template <int N, class T, bool CopyData = true, class RefCountT = EmptyClass, bool Continuous = true>
class InterpSimplex : public NDInterpolator<N,T,CopyData,RefCountT,Continuous> {
public:
  typedef NDInterpolator<N,T,CopyData,RefCountT,Continuous> super;
  
  template <class IterT1, class IterT2>  
  InterpSimplex(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end)
    : super(begins_begin, begins_end, ends_begin, f_begin, f_end)
  {}
  template <class IterT1, class IterT2>  
  InterpSimplex(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end, RefCountT &refF)
    : super(begins_begin, begins_end, ends_begin, f_begin, f_end, refF)
  {}
  
  template <class IterT1, class IterT2>
  void interp_vec(int n, IterT1 coord_iter_begin, IterT1 coord_iter_end, IterT2 i_result) const {
    assert(N == coord_iter_end - coord_iter_begin);
	
	array<int,N> cell_index, v_index;
	array<std::pair<T, int>,N> xipair;	
	int c;
	T y, v0, v1;
	for (int i=0; i<n; i++) {			// for each point
	  for (int dim=0; dim<N; dim++) {
	    typename super::grid_type const &grid(super::m_gridList[dim]);
		c = find_cell(dim, coord_iter_begin[dim][i]);
		if (c == -1) {					// before first grid point
		  y = 1.0;
		} else if (c == grid.size()-1) {	// after last grid point
		  y = 0.0;
		} else {
		  y = (coord_iter_begin[dim][i] - grid[c]) / (grid[c + 1] - grid[c]);
		  if (y < 0.0) y=0.0;
		  else if (y > 1.0) y=1.0;
		}
        xipair[dim].first = y;
	    xipair[dim].second = dim;		
		cell_index[dim] = c;
	  }  	  
      // sort xi's and get the permutation    
      std::sort(xipair.begin(), xipair.end(), [](std::pair<T, int> const &a, std::pair<T, int> const &b) {
        return (a.first < b.first);
      });
      // walk the vertices of the simplex determined by the permutation  
      for (int j=0; j<N; j++) {
        v_index[j] = 1;
      }		      
	  v0 = get_f_val(cell_index, v_index);
	  y = v0;
      for (int j=0; j<N; j++) {
	    v_index[xipair[j].second]--;		
		v1 = get_f_val(cell_index, v_index);
	    y += (1.0 - xipair[j].first) * (v1-v0);		// interpolate
		v0 = v1;
	  }
	  *i_result++ = y;
    }    
  }  
};

template <int N, class T, bool CopyData = true, class RefCountT = EmptyClass, bool Continuous = true>
class InterpMultilinear : public NDInterpolator<N,T,CopyData,RefCountT,Continuous> {
public:
  typedef NDInterpolator<N,T,CopyData,RefCountT,Continuous> super;
  
  template <class IterT1, class IterT2>  
  InterpMultilinear(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end)
    : super(begins_begin, begins_end, ends_begin, f_begin, f_end)
  {}
  template <class IterT1, class IterT2>  
  InterpMultilinear(IterT1 begins_begin, IterT1 begins_end, IterT1 ends_begin, IterT2 f_begin, IterT2 f_end, RefCountT &refF)
    : super(begins_begin, begins_end, ends_begin, f_begin, f_end, refF)
  {}

  template <class IterT1, class IterT2>
  static T linterp_nd_unitcube(IterT1 f_begin, IterT1 f_end, IterT2 xi_begin, IterT2 xi_end) {
    int n = xi_end - xi_begin;
    int f_len = f_end - f_begin;
    assert(1 << n == f_len);
    T sub_lower, sub_upper;
    if (n == 1) {
      sub_lower = f_begin[0];
  	  sub_upper = f_begin[1];
    } else {
      sub_lower = linterp_nd_unitcube(f_begin, f_begin + (f_len/2), xi_begin + 1, xi_end);
	  sub_upper = linterp_nd_unitcube(f_begin + (f_len/2), f_end, xi_begin + 1, xi_end);
    }  
    T result = sub_lower + (*xi_begin)*(sub_upper - sub_lower);
    return result;
  }

  template <class IterT1, class IterT2>
  void interp_vec(int n, IterT1 coord_iter_begin, IterT1 coord_iter_end, IterT2 i_result) const {
    assert(N == coord_iter_end - coord_iter_begin);
	array<int,N> index;
	int c;
	T y, xi;
    vector<T> f(1 << N);
	array<T,N> x;
	
	for (int i=0; i<n; i++) {								// loop over each point
	  for (int dim=0; dim<N; dim++) {						// loop over each dimension
	    auto const &grid(super::m_gridList[dim]);		
		xi = coord_iter_begin[dim][i];
		c = find_cell(dim, coord_iter_begin[dim][i]);
		if (c == -1) {					// before first grid point
		  y = 1.0;
		} else if (c == grid.size()-1) {	// after last grid point
		  y = 0.0;
		} else {
		  y = (coord_iter_begin[dim][i] - grid[c]) / (grid[c + 1] - grid[c]);
		  if (y < 0.0) y=0.0;
		  else if (y > 1.0) y=1.0;
		}
		index[dim] = c;
		x[dim] = y;
	  }
	  // copy f values at vertices
	  for (int v=0; v < (1 << N); v++) {					// loop over each vertex of hypercube
        f[v] = get_f_val(index, v);
	  }
	  *i_result++ = linterp_nd_unitcube(f.begin(), f.end(), x.begin(), x.end());
	}
  }
};	

typedef InterpSimplex<1,double> NDInterpolator_1_S;
typedef InterpSimplex<2,double> NDInterpolator_2_S;
typedef InterpSimplex<3,double> NDInterpolator_3_S;
typedef InterpSimplex<4,double> NDInterpolator_4_S;
typedef InterpSimplex<5,double> NDInterpolator_5_S;
typedef InterpMultilinear<1,double> NDInterpolator_1_ML;
typedef InterpMultilinear<2,double> NDInterpolator_2_ML;
typedef InterpMultilinear<3,double> NDInterpolator_3_ML;
typedef InterpMultilinear<4,double> NDInterpolator_4_ML;
typedef InterpMultilinear<5,double> NDInterpolator_5_ML;


#endif //_linterp_h