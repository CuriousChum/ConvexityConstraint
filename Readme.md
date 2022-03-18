# Optimization subject to convexity constraints
## A geometric approach, with C++ routines

This repository gathers some numerical codes devoted to optimization problems posed on the set of convex functions, or convex bodies. The convexity constraint is imposed using a barrier function, which is related with the entropy of the image measure defined by the convex function.

The same approach is illustrated and commented in detail in the following python Python notebooks:
- Convex bodies, and the problems of Minkowski and Meissner. [Notebook](https://nbviewer.org/github/Mirebeau/AdaptiveGridDiscretizations_showcase/blob/master/Notebooks_Algo/Meissner.ipynb)
- Convex functions, and the principal agent or *monopolist* problem. [Notebook](https://nbviewer.org/github/Mirebeau/AdaptiveGridDiscretizations/blob/master/Notebooks_Algo/Monopolist.ipynb)

Unfortunately, the performance and convenience of the Python implementation do not match the (earlier) c++ implementation in a number of regards, which is why I make the latter available in this code repo. 


## Dependencies

The following open source libraries are needed to compile the provided code. They can be installed using e.g. the provided anaconda environment:
 - The [CGAL](http://www.cgal.org/) computational geometry algorithms library.
 - The [nlopt](http://ab-initio.mit.edu/nlopt) non-linear optimization routines.
 - The [Eigen](http://eigen.tuxfamily.org/) linear solver.

Optionally, the Mathematica(R) software can be used to visualize the results, see the provided notebooks.

## Installation


<!---
Original Meissner code : 
/Users/mirebeau/Dropbox/Programmes/2015/9_Septembre/Minkowski 

Original Principal agent C++ code : 
/Users/mirebeau/Dropbox/Programmes/2015/1_Janvier/OptimizationTests/CGalTest

Original Principal agent visualization code : 
/Users/mirebeau/Dropbox/Shared/Quentin/SubgradientBarrier/Mathematica
--->