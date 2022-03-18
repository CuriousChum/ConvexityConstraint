//
//  Minkowski_2_Test.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 20/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Minkowski_2_Test_h
#define CGalTest_Minkowski_2_Test_h

#include "MinkowskiConstraint_2.h"
#include "nlopt.hpp"
#include "MinkowskiEnergy_2.h"

namespace Minkowski_2_Test {
    using namespace Minkowski_2;
    
    void Test0(){
        
        const int n=6;
        std::vector<Vector> pts;
        for(int i=0; i<n; ++i) pts.push_back({cos(2*i*M_PI/n),sin(2*i*M_PI/n)});
        
        ConvexityConstraint cstr(pts);
        
        std::vector<ScalarType> input;
        input.resize(n,1);
        cstr.SetValues(input);
        cstr.Compute(31);
        
        std::cout << cstr << "\n";
        //Hexagon edge length : 2/sqrt(3) = 1.1547...
        //Hexagon area : 2 sqrt(3) = 3.46...
    }
    
    
    // Solve a 2D Minkowski problem
    
    struct {
        ScalarType operator () (const Vector & v) const {
            const ScalarType x=v[0], y=v[1];
            return sqrt(x*x+2*y*y-2*x*y);
//            return std::max(x+2*y,0.)+sqrt(x*x+y*y);
        }
    } norm;
    
    
    void Test1(){ // Issue : involves the stockage and inversion of a non-sparse hessian.
        // Construct points and target
        const int n=20;
        std::vector<Vector> pts;
        for(int i=0; i<n; ++i) pts.push_back({cos(2*i*M_PI/n),sin(2*i*M_PI/n)});

        ConvexityConstraint cstr(pts);
        std::vector<ScalarType> input;
        input.reserve(n);
        for(auto v : pts) input.push_back(norm(v));
        cstr.SetValues(input);
        cstr.Compute(ConstraintType::RValues);
        const std::vector<ScalarType> target = cstr.values;
        
        // Now, we want to solve for target.
        
        MinkowskiEnergy energy(pts,target);
        
        PositivityConstraint pos;
        std::vector<NS::Functionnal*> constraints = {&pos};
        NS::NewtonConstrained newton;
        newton.maxIter = 50;
        newton.multiplier = 1./pts.size();
        newton.multiplierBound *= newton.multiplier;
        
        NS::VectorType xx(pts.size());
        xx.fill(1);
        newton.Solve(energy, xx, constraints);
        
        std::vector<ScalarType> x(xx.size()); for(int i=0; i<xx.size(); ++i) x[i]=xx[i];
        cstr.SetValues(x);
        cstr.Compute(ConstraintType::RValues);
        std::vector<ScalarType> diff(target.size());
        for(int i=0; i<target.size(); ++i) diff[i] = cstr.values[i]/cstr.area - target[i];
        auto & os = std::cout
        ExportArrayArrow(cstr.values)
        ExportArrayArrow(target)
        ExportArrayArrow(diff)
        << "\n";
    }
    
    void Test2(){
        // Construct points and target
        const int n=200;
        std::vector<Vector> pts;
        for(int i=0; i<n; ++i) pts.push_back({cos(2*i*M_PI/n),sin(2*i*M_PI/n)});

        ConvexityConstraint cstr(pts);
        std::vector<ScalarType> input;
        input.reserve(n);
        for(auto v : pts) input.push_back(norm(v));
        cstr.SetValues(input);
        cstr.Compute(ConstraintType::RValues);
        const std::vector<ScalarType> target = cstr.values;
        
        // Now, we want to solve for target.
        
        MinkowskiEnergy energy(pts,target);
        
        
        nlopt::opt algo(nlopt::LD_LBFGS, n);
        algo.set_min_objective(MinkowskiEnergy::Evaluate, &energy);
        algo.set_lower_bounds(0);
        algo.set_maxeval(1000);
        std::vector<ScalarType> x; x.resize(n,1);
        ScalarType minimum;
        algo.optimize(x,minimum);
        
        
        cstr.SetValues(x);
        cstr.Compute(ConstraintType::RValues);
        std::vector<ScalarType> diff(target.size());
        for(int i=0; i<target.size(); ++i) diff[i] = cstr.values[i]/cstr.area - target[i];
        auto & os = std::cout
        ExportArrayArrow(cstr.values)
        ExportArrayArrow(target)
        ExportArrayArrow(diff)
        ExportArrayArrow(x)
        ExportVarArrow(minimum)
        << "\n";
    }
    
    void Test3(){
        
    }
}

#endif
