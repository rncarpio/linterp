##
## Copyright (c) 2011 Ronaldo Carpio
##                                     
## Permission to use, copy, modify, distribute and sell this software
## and its documentation for any purpose is hereby granted without fee,
## provided that the above copyright notice appear in all copies and   
## that both that copyright notice and this permission notice appear
## in supporting documentation.  The authors make no representations
## about the suitability of this software for any purpose.          
## It is provided "as is" without express or implied warranty.
##                                                            
  
import scipy, scipy.interpolate, types
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyublas, _linterp_python as _linterp
import itertools, time, operator
		
# test functions for interpolation
def nd_sin_sum(x):
	s = x[0]
	for xi in x[1:]:
		s += xi
	return scipy.sin(s)

def nd_circle_sine(veclist):
	sum = scipy.zeros(len(veclist[0]))
	for vec in veclist:
		sum += vec*vec
	r = scipy.sqrt(sum)	
	return scipy.sin(r)
	
def ndgrid(gridList):
	coord_list = zip(*itertools.product(*gridList))
	return [scipy.array(l) for l in coord_list]

def exp_linspace(first, last, n):
	x1 = scipy.linspace(0, 1, n/2)
	x2 = scipy.exp(x1)
	x3 = x2 - x2[0] + (x2[1] - x2[0])/2
	x4 = x3 / x3[-1]
	x5 = scipy.array(list(reversed(-x4)) + list(x4))
	len = last-first
	mid = (first + last)/2
	x6 = x5 * (len/2)
	x7 = x6 + mid
	return x7
	
def test_nd(n=2, fn=nd_circle_sine, plot=True, timing=False, use_griddata=True, reps=100, f_len=20, p_len=30):
	# classes that implement interpolation
	class InterpGriddata:
		def __init__(self, grid_list, f_vals):
			(self.grid_list, self.f_vals) = (grid_list, f_vals)
			self.ndgrid = tuple(ndgrid(grid_list))
		def interp_vec(self, xi):			
			result = scipy.interpolate.griddata(self.ndgrid, self.f_vals, tuple(xi), method='linear', fill_value=0.0)
			return result

	f_grid_list = [exp_linspace(-5, 5, f_len)] * n;
	f_ndgrid = ndgrid(f_grid_list)
	f_vals = fn(f_ndgrid);												# evaluate fn at each grid point
	
	# interpolation points grid
	p_grid_list = [scipy.linspace(-8, 8, p_len)] * n;
	p_ndgrid = ndgrid(p_grid_list)		

	# methods to compare
	if (n==2): interpObjList = [(_linterp.Interp_2_ML, "ND multilinear recursive"), (_linterp.Interp_2_S, "ND simplex")]
	elif (n==3): interpObjList = [(_linterp.Interp_3_ML, "ND multilinear recursive"), (_linterp.Interp_3_S, "ND simplex")]
	elif (n==4): interpObjList = [(_linterp.Interp_4_ML, "ND multilinear recursive"), (_linterp.Interp_4_S, "ND simplex")]
	else: assert(false)
	if (use_griddata):
		# griddata can't handle large numbers of points
		interpObjList.append((InterpGriddata, "griddata"))
		
	if (plot):
		if (n==2):
			plot3d("actual fn", f_ndgrid[0], f_ndgrid[1], f_vals);		# original fn		
			for (i, (objClass, title)) in enumerate(interpObjList):
				interp_obj = objClass(f_grid_list, f_vals)				
				interp_vals = interp_obj.interp_vec(p_ndgrid)			
				plot3d(title, p_ndgrid[0], p_ndgrid[1], interp_vals);			

	if (timing):
		for (i, (objClass, title)) in enumerate(interpObjList):
			print("executing %s" % title)
			t1 = time.time()
			for j in range(reps):
				interp_obj = objClass(f_grid_list, f_vals)
			t2 = time.time()
			for j in range(reps):
				interp_vals = interp_obj.interp_vec(p_ndgrid)
			t3 = time.time()
			print("%s: %d iters, setup time %f (%f per iter), interp %f (%f per iter)" % (title, reps, (t2-t1), (t2-t1)/reps, (t3-t2), (t3-t2)/reps))
	
	# check accuracy	
	def in_grid(grid_list, x):
		in_grid_list = [(xi >= grid[0] and xi <= grid[-1]) for (grid, xi) in zip(grid_list, x)]
		return (sum(in_grid_list) == len(x))		
	# remove points outside the grid
	in_grid_list = []
	for (f_grid, p_grid) in zip(f_grid_list, p_grid_list):
		in_points = [x for x in p_grid if (x >= f_grid[0] and x <= f_grid[-1])]
		in_grid_list.append(scipy.array(in_points))
	in_ndgrid = ndgrid(in_grid_list)
	for (i, (objClass, title)) in enumerate(interpObjList):		
		interp_obj = objClass(f_grid_list, f_vals)
		interp_vals = interp_obj.interp_vec(in_ndgrid)
		true_vals = fn(in_ndgrid)
		err = interp_vals - true_vals
		ss = scipy.sum(err*err)
		err_max_i = scipy.argmax(err*err)
		print("%s: error sum of squares %f" % (title, ss))
		print("largest error: %f" % err[err_max_i])		
		if (plot and n==2): plot3d(title + " errors", in_ndgrid[0], in_ndgrid[1], err);			
		
def plot3d(title, xlist, ylist, zlist):
	fig = plt.figure()	
	ax = Axes3D(fig)
	ax.scatter(xlist, ylist, zlist, s=2)
	plt.title(title)
	
	
if __name__ == "__main__":
	test_nd()
	
	