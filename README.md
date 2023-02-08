# geometrytp3
Elan Biswas and Ignacio Villanueva

## 3.1 Uniform map
We set up and solve the linear system presented in the problem. As the layout of the code is very similar for the uniform mapping as it is for the harmonic mapping, we created a helper function that takes as input a weight function and the mesh and outputs the parameter coordinates.

We encountered an issue in which some of the boundary vertices stuck out in jagged spikes beyond the fixed circle. This was due to the fact that we added edge weights to matrix elements for all vertices rather than only non-boundary vertices. This was resolved by adding a boundary check for each vertex in the edge.  

## 3.2 Harmonic map
The solution to this problem was the exact same as for the uniform map save for the cotangent weight computation. We pass this function into the higher-order(ish) function calc_parameterization to compute the mapping.

## 3.3 Texture mapping
We attempted to solve this problem by scaling the UV coordinates to the range specified by iRepeats, but this resulted in the example texture mapping only green onto the bunny. Due to time constraints we were unable to resolve this issue.

## 3.4 Parameterization distortion
Implementing the function was a matter of implementing the formulas given in the spec. The results show that both angle and area distortion are consistently greater in the uniform mapping than in the harmonic mapping. 