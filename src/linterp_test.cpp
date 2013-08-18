
#include <stdio.h>
#include <ctime>
#include "linterp.h"

const int N=4;

dVec linspace(double first, double last, int len) {
  dVec result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) {
    result[i] = first + i*step;
  }
  return result;
}

array<double,N> get_grid_point(vector<dVec> const &gridList, array<int,N> const &index) {
  array<double,N> result;
  for (int i=0; i<N; i++) {
    result[i] = gridList[i][index[i]];
  }
  return result;
}

template <class IterT>
unsigned int IndexListTo1D(iVec const &lenArray, IterT indexArray_begin, IterT indexArray_end) {
  unsigned int i;  
  assert(lenArray.size() == indexArray_end - indexArray_begin);
  if (lenArray.size() == 0) {
    return 0;
  }
  unsigned int result = indexArray_begin[0];
  for (i=1; i<lenArray.size(); i++) {
    result *= lenArray[i];           
    result += indexArray_begin[i];
  }                         
  return result;
}               

template <class IterT1, class IterT2>
void Index1DToArray(unsigned int index, IterT1 lenArray_begin, IterT1 lenArray_end, IterT2 indexArray_begin, IterT2 indexArray_end) {
  int i;
  int len = lenArray_end - lenArray_begin;
  for (i=len-1; i>=0; i--) {
    indexArray_begin[i] = index % lenArray_begin[i];      
        index /= lenArray_begin[i];     
  }                          
  assert(index == 0);
  return;            
}        

template <class IterT>
void ndgrid(IterT gridList_begin, IterT gridList_end, vector<dVec> &result) {
  int total=1;
  int n = gridList_end - gridList_begin;
  iVec lenArray(n);
  for (int i=0; i<n; i++) {
    total *= gridList_begin[i].size();
	lenArray[i] = gridList_begin[i].size();
  }
  for (int i=0; i<n; i++) {
    result[i].resize(total);
  }
  iVec index(n);
  
  for (int i=0; i<total; i++) {
    Index1DToArray(i, lenArray.begin(), lenArray.end(), index.begin(), index.end());
	for (int j=0; j<n; j++) {
	  result[j][i] = gridList_begin[j][index[j]];
	}
  }
}
  
template<class T>
double nd_circle_sine(T begin, T end) {
  double sum = 0.0;
  for (; begin != end; begin++) {
    sum += (*begin)*(*begin);
  }
  return sqrt(sum);
}

template <class IterT>
std::pair<vector<typename IterT::value_type::const_iterator>, vector<typename IterT::value_type::const_iterator> > get_begins_ends(IterT iters_begin, IterT iters_end) {
  typedef typename IterT::value_type T;
  typedef vector<typename T::const_iterator> VecT;
  int N = iters_end - iters_begin;
  std::pair<VecT, VecT> result;
  result.first.resize(N);
  result.second.resize(N);
  for (int i=0; i<N; i++) {
    result.first[i] = iters_begin[i].begin();
	result.second[i] = iters_begin[i].end();
  }
  return result;
}

int main(int argc, char **argv) {
  //printf("press return to continue\n");
  //char pcTemp[1024];
  //fgets(pcTemp, 1024, stdin);
  
  int f_len = 30;
  int p_len = 40;
  dVec f_grid1 = linspace(-5, 5, f_len);
  dVec p_grid1 = linspace(-8, 8, p_len);
  
  vector<dVec> f_gridList(N);
  for (int i=0; i<N; i++) {
    f_gridList[i] = f_grid1;
  }
  
  array<int,N> f_sizes, p_sizes;
  for (int i=0; i<N; i++) {
    f_sizes[i] = f_len;
	p_sizes[i] = p_len;
  }
  boost::multi_array<double,N> f(f_sizes);
  array<int,N> index;
  array<double,N> arg;
  for (int i=0; i<f_sizes[0]; i++) {
    for (int j=0; j<f_sizes[1]; j++) {
	  for (int k=0; k<f_sizes[2]; k++) {
	    for (int l=0; l<f_sizes[3]; l++) {
          index[0] = i;
		  index[1] = j;
		  index[2] = k;
		  index[3] = l;		  
		  arg = get_grid_point(f_gridList, index);  
		  f(index) = nd_circle_sine(arg.begin(), arg.end());
		}
      }
    }
  }	

  vector<dVec> p_gridList(N);
  for (int i=0; i<N; i++) {
    p_gridList[i] = p_grid1;
  }
  dVec result;
  int total;
  vector<dVec> p_ndgrid(N);
  ndgrid(p_gridList.begin(), p_gridList.end(), p_ndgrid);
  total = p_ndgrid[0].size();
  result.resize(total);
  printf ("interpolating %d points\n", total);
  auto begins_ends = get_begins_ends(f_gridList.begin(), f_gridList.end());
  
  InterpSimplex<4, double, false> s1(begins_ends.first.begin(), f_sizes, f.data(), f.data() + f.num_elements());
  InterpMultilinear<4, double, false> ml1(begins_ends.first.begin(), f_sizes, f.data(), f.data() + f.num_elements());

  clock_t t1, t2;
  t1 = clock();	
  ml1.interp_vec(total, p_ndgrid.begin(), p_ndgrid.end(), result.begin());
  t2 = clock();
  printf("multilinear interp: %d clocks, %f sec\n", (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  t1 = clock();	
  s1.interp_vec(total, p_ndgrid.begin(), p_ndgrid.end(), result.begin());
  t2 = clock();
  printf("simplex interp: %d clocks, %f sec\n", (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  return 0;
}

