//
//  Meissner_2_Test.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <fstream>
#include "Headers/Meissner_2.h"
#include "nlopt.hpp"

namespace Meissner_2_Test {
    using namespace Meissner_2;
    typedef N::VectorType NV;
    
	/**
	 Numerically solve the two dimensional Meissner problem, also known as Reulaux problem,
	 of the convex body with unit width and minimum volume.
	 
	 Input :
	  - n : the number of facets of the convex body. The facet normals are regularly spaced unit vectors.
	  - delta : a small parameter defining the initial condition, which is a perturbation of the sphere of radius 1/2.
	 (Setting delta = 0 fails, since the sphere of radius 1/2 is a stationnary point.)
	 */
    void OptimizeMeissner2(int n=20, ScalarType delta=0.05){
        Minkowski_2::Functional mink;
        G::VectorList & pts = mink.pts;
        
        pts.resize(Dimension,n);
        for(int i=0; i<n; ++i) pts.col(i) = G::Vector{cos(2*i*M_PI/n),sin(2*i*M_PI/n)};

		/*
		 Note :
		 This code was originally developed for the 2.4.2 version of nlopt.
		 Recent tests show that the performance with the 2.7.1 version is downgraded.
		 
		 2.7.2 -> minimal area = 0.735
		 2.4.2 -> minimal area = 0.705678
		 exact value for continuous problem -> minimal area = 0.704771
		 (The obtained shape is also less sharp with 2.7.2.)
		 
		 */
		
		nlopt::opt algo=nlopt::opt(nlopt::LD_SLSQP , n/2);
//		nlopt::opt algo=nlopt::opt(nlopt::LD_MMA , n/2);
        algo.set_min_objective(EvaluateVolume, static_cast<void*>(&mink));
        algo.add_inequality_constraint(EvaluateConstraint, static_cast<void*>(&mink));
        algo.set_lower_bounds(0);
        algo.set_upper_bounds(1);
        algo.set_maxeval(1000);

//        algo.set_stopval(0.71);

        ScalarType minimum;
        std::vector<ScalarType> x(n/2,0.5),init;
        for(int i=0; i<n/2; ++i) x[i]+=delta*cos(3*2*i*M_PI/n);
        init=x;
        
        std::cout << "Opt begin\n";
        nlopt::result result = algo.optimize(x,minimum);
		std::cout << result << std::endl;
        std::cout << "Opt end\n";
        
        std::ofstream os; os.open("Meissner_2.txt");
        os << "{"
        ExportVarArrow(minimum)
        ExportVarArrow(result)
        ExportArrayArrow(x)
        ExportArrayArrow(init)
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(0))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(1))))
        << "}";
    }
}


int main(int argc, const char * argv[]) {
    Meissner_2_Test::OptimizeMeissner2(200,0.05);
    return 0;
}
