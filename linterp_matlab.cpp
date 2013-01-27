
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <cstdarg>
#include <string>
#include <vector>
#include <array>
#include <functional>
#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"
#include "linterp.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
				  
  int n_dims = (nrhs-1) / 2;
  assert(n_dims > 0);
  if (nrhs < 3 || (nrhs%2) == 0) {
    mexErrMsgTxt("incorrect number of args. linterp(grid1,grid2,grid3,...,V,Y1,Y2,Y3,...)");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  for (int i=0; i<nrhs; i++) {
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i])) {
	  mexErrMsgTxt("Input must be noncomplex double array");
	}
  }
  
  // grid vectors
  vector<const double*> grid_list;
  vector<size_t> grid_len_list;
  size_t total_size = 1;
  for (int i=0; i<n_dims; i++) {				// vectors will be flattened
    grid_list.push_back(mxGetPr(prhs[i]));
	size_t grid_len = mxGetNumberOfElements(prhs[i]);
	grid_len_list.push_back(grid_len);
	total_size *= grid_len;
  }
  
  // F array
  if (total_size != mxGetNumberOfElements(prhs[n_dims])) {
    char pcTemp[1024];
	sprintf(pcTemp, "array sizes do not match. total_size=%d, mxGetNumberOfElements=%d", total_size, mxGetNumberOfElements(prhs[n_dims]));
    mexErrMsgTxt(pcTemp);
  }
  const double *p_F = mxGetPr(prhs[n_dims]);
  
  // xi vectors
  vector<const double*> xi_list;    
  size_t min_len = mxGetNumberOfElements(prhs[n_dims+1]);
  for (int i=n_dims+1; i<nrhs; i++) {			// vectors will be flattened
    xi_list.push_back(mxGetPr(prhs[i]));
	min_len = (mxGetNumberOfElements(prhs[i]) < min_len) ? mxGetNumberOfElements(prhs[i]) : min_len;
  }

  // create result
  plhs[0] = mxCreateDoubleMatrix(1, min_len, mxREAL);
  double *p_result = mxGetPr(plhs[0]);
  
  // call interp
  if (n_dims == 1) {
    const int N = 1;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 2) {
    const int N = 2;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 3) {
    const int N = 3;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else {
    mexErrMsgTxt("dimension not implemented");
  }
}


