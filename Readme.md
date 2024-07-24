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

Additional dependency for visualization below.

Optionally, the Mathematica(R) software can be used to visualize the results, see the provided notebooks.

## Installation and use

To use this software, proceed according to the following steps:
- Download the current repository.
- Download the [HamiltonFastMarching repository](https://github.com/Mirebeau/HamiltonFastMarching).
- Install conda or (preferably) [miniconda](https://docs.conda.io/en/latest/miniconda.html).
- Open a terminal and cd in the directory containing this `Readme.md` file.
- Create a conda environnement `conda install -f ConvexityConstraint.yaml`
- Open the file `common.cmake` and update the `Conda_ConvexityConstraint_Env_Dir` and `HamiltonFastMarching_Source_Dir` paths.
- Install the [cmake](https://cmake.org) tool.
- Generate the project of your choice : `Meissner`, `Monopolist`, or both, using cmake. 
- Compile the project.
- Run the generated executable. Use the the argument `--help`, to learn about the proper command line arguments.
- Open the `[Meissner|Monopolist]/Visualization/CommonInitialization.nb` notebook using [Mathematica](https://www.wolfram.com/mathematica/).
- Setup the binary path in the `CommonInitialization.nb` notebook, and execute it. 
- Open the visualization notebook of your choice. 

## Basic Visualization of Gradient and Hessian of Monopolist

Dependencies: NumPy, Matplotlib, SciPy
All dependencies are added to ConvexityConstraint.yaml. Make sure the virtual environment is activated with `conda activate ConvexityConstraint`.

Modify `vis.py` to show different visualizations. I am currently creating a different project with similar functionalities, but for the time being, this should do.


<!---
Original Meissner code : 
/Users/mirebeau/Dropbox/Programmes/2015/9_Septembre/Minkowski 

Original Principal agent C++ code : 
/Users/mirebeau/Dropbox/Programmes/2015/1_Janvier/OptimizationTests/CGalTest

Original Principal agent visualization code : 
/Users/mirebeau/Dropbox/Shared/Quentin/SubgradientBarrier/Mathematica
--->
