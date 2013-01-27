
Copyright (c) 2012 Ronaldo Carpio
									 
Permission to use, copy, modify, distribute and sell this software
and its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice appear in all copies and   
that both that copyright notice and this permission notice appear
in supporting documentation.  The authors make no representations
about the suitability of this software for any purpose.          
It is provided "as is" without express or implied warranty.

Project page: http://rncarpio.github.com/linterp
                                                            
This is a C++ header-only library for N-dimensional linear interpolation on a rectangular grid. Implements two methods:
* Multilinear: Interpolate using the N-dimensional hypercube containing the point. Interpolation step is O(2^N) 
* Simplicial: Interpolate using the N-dimensional simplex containing the point. Interpolation step is O(N log N), but less accurate.
Requires boost/multi_array library.

For a description of the algorithms, see:
* Weiser & Zarantonello (1988), "A Note on Piecewise Linear and Multilinear Table Interpolation in Many Dimensions", _Mathematics of Computation_ 50 (181), p. 189-196
* Davies (1996), "Multidimensional Triangulation and Interpolation for Reinforcement Learning", _Proceedings of Neural Information Processing Systems 1996_