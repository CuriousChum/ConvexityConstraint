//
//  Meissner_2_Test.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Meissner_2_Test_h
#define Minkowski_Meissner_2_Test_h

#include <fstream>
#include "Meissner_2.h"
#include "nlopt.hpp"
namespace Meissner_2_Test {
    using namespace Meissner_2;
    typedef N::VectorType NV;
    
    void Test1(int n=20, ScalarType delta=0.05){
        Minkowski_2::Functional mink;
        G::VectorList & pts = mink.pts;
        
        pts.resize(Dimension,n);
        for(int i=0; i<n; ++i) pts.col(i) = G::Vector{cos(2*i*M_PI/n),sin(2*i*M_PI/n)};

        nlopt::opt algo=nlopt::opt(nlopt::LD_SLSQP , n/2); //LD_MMA
        algo.set_min_objective(EvaluateVolume, static_cast<void*>(&mink));
        algo.add_inequality_constraint(EvaluateConstraint, static_cast<void*>(&mink));
        algo.set_lower_bounds(0);
        algo.set_upper_bounds(1);
        algo.set_maxeval(1000);
        
        algo.set_stopval(0.71);
        
        ScalarType minimum;
        std::vector<ScalarType> x(n/2,0.5),init;
        for(int i=0; i<n/2; ++i) x[i]+=delta*cos(3*2*i*M_PI/n);
        init=x;
        
        std::cout << "Opt begin\n";
        nlopt::result result = algo.optimize(x,minimum);
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

#endif
