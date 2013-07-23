#include <ctime>
#include "linterp.h"

// return an evenly spaced 1-d grid of doubles.
std::vector<double> linspace(double first, double last, int len) {
  std::vector<double> result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) { result[i] = first + i*step; }
  return result;
}

// the function to interpolate.
double fn (double x1, double x2) { return sin(x1 + x2); }

int main (int argc, char **argv) {
  const int length = 10;
    
  // construct the grid in each dimension. 
  // note that we will pass in a sequence of iterators pointing to the beginning of each grid
  std::vector<double> grid1 = linspace(0.0, 3.0, length);
  std::vector<double> grid2 = linspace(0.0, 3.0, length);
  std::vector< std::vector<double>::iterator > grid_iter_list;
  grid_iter_list.push_back(grid1.begin());
  grid_iter_list.push_back(grid2.begin());
  
  // the size of the grid in each dimension
  array<int,2> grid_sizes;
  grid_sizes[0] = length;
  grid_sizes[1] = length;
  
  // total number of elements
  int num_elements = grid_sizes[0] * grid_sizes[1];
  
  // fill in the values of f(x) at the gridpoints. 
  // we will pass in a contiguous sequence, values are assumed to be laid out C-style
  std::vector<double> f_values(num_elements);
  for (int i=0; i<grid_sizes[0]; i++) {
    for (int j=0; j<grid_sizes[1]; j++) {
	  f_values[i*grid_sizes[0] + j] = fn(grid1[i], grid2[j]);
	}
  }
  
  // construct the interpolator. the last two arguments are pointers to the underlying data
  InterpMultilinear<2, double> interp_ML(grid_iter_list.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);
  
  // interpolate one value
  array<double,2> args = {1.5, 1.5};
  printf("%f, %f -> %f\n", args[0], args[1], interp_ML.interp(args.begin()));

  // interpolate multiple values: create sequences for each coordinate
  std::vector<double> interp_grid = linspace(0.0, 3.0, length*10);
  int num_interp_elements = interp_grid.size() * interp_grid.size();
  std::vector<double> interp_x1(num_interp_elements);
  std::vector<double> interp_x2(num_interp_elements);
  for (int i=0; i<interp_grid.size(); i++) {
    for (int j=0; j<interp_grid.size(); j++) {
	  interp_x1[i*interp_grid.size() + j] = interp_grid[i];
	  interp_x2[i*interp_grid.size() + j] = interp_grid[j];
	}
  }
  std::vector<double> result(num_interp_elements);
  
  // pass in a sequence of iterators, one for each coordinate
  std::vector< std::vector<double>::iterator > interp_x_list;
  interp_x_list.push_back(interp_x1.begin());
  interp_x_list.push_back(interp_x2.begin());
  
  // interpolate a sequence of values
  clock_t t1, t2;
  t1 = clock();	
  interp_ML.interp_vec(num_interp_elements, interp_x_list.begin(), interp_x_list.end(), result.begin());
  t2 = clock();
  printf("multilinear: %d interpolations, %d clocks, %f sec\n", num_interp_elements, (t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

  // calculate the squared errors
  std::vector<double> true_f_vals(num_interp_elements);
  double SSE = 0.0;
  for (int i=0; i<num_interp_elements; i++) {
    true_f_vals[i] = fn(interp_x1[i], interp_x2[i]);
	double diff = true_f_vals[i] - result[i];
	SSE += diff*diff;
  }
  printf("sum of squared errors: %f\n", SSE);
  return 0;
}