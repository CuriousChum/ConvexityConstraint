//
//  Minkowski_2_Test.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include "Headers/Minkowski_2.h"

namespace Minkowski_2_Test {
    using namespace Minkowski_2;
    std::ostream & os = std::cout;
    typedef N::VectorType NV;
    typedef G::Vector GV;
    
    struct {
        ScalarType operator () (const G::Vector & v) const {
            const ScalarType x=v[0], y=v[1];
            return sqrt(x*x+2*y*y-2*x*y);
            //            return std::max(x+2*y,0.)+sqrt(x*x+y*y);
        }
    } norm;
    
    void Test1(){
        const int n=20;
        Functional op;
        G::VectorList & pts = op.pts;
        pts.resize(2,n);
        for(int i=0; i<n; ++i) pts.col(i) = G::Vector{cos(2*i*M_PI/n),sin(2*i*M_PI/n)};
        
        N::VectorType sol(n), target(n), sol_copy;
//        input.setConstant(1);
        for(int i=0; i<n; ++i) sol[i]=norm(pts.col(i));
        sol_copy=sol;
        op.Compute(sol, target);
        
        G::Vector normAvg = pts * target;
        
        os
        ExportArrayArrow(N::StdFromEigen(op.Barycenter())) << "\n"
        ExportArrayArrow(N::StdFromEigen(NV(sol-sol_copy))) << "\n"
        ExportArrayArrow(N::StdFromEigen(target)) << "\n"
        ExportArrayArrow(N::StdFromEigen(GV(pts*target)))
        << "\n";

        N::VectorType x(n), val(n); x.setConstant(1);
        op.Compute(x,val);
        NV diff=val-target; 
        op.DescentDirection(diff);
        
        os
        ExportArrayArrow(N::StdFromEigen(x))
        ExportArrayArrow(N::StdFromEigen(val))
        ExportArrayArrow(N::StdFromEigen(diff))
        << "\n";
        
        NV x2 = x+diff;
        op.Compute(x2, val);
        
        os
        ExportArrayArrow(N::StdFromEigen(val))
        ExportArrayArrow(N::StdFromEigen(x2))
        ExportArrayArrow(N::StdFromEigen(NV(val-target)))
        ExportArrayArrow(N::StdFromEigen(NV(sol-x2)))
        ExportVarArrow((sol-x2).norm())
        << "\n";
        
        // Indeed, the system is solved in one step
        
    }
    
    void Test2(){
        // Construct points and target
        const int n=20;
        Functional op;
        G::VectorList & pts = op.pts;
        pts.resize(2,n);
        for(int i=0; i<n; ++i) pts.col(i) = G::Vector{cos(2*i*M_PI/n),sin(2*i*M_PI/n)};
        
        N::VectorType input(n), target(n);
        for(int i=0; i<n; ++i) input[i]=norm(pts.col(i));
        
        op.Compute(input, target);
        
        // Now, we want to solve for target.
        N::NewtonScheme solver;
        N::VectorType x = N::VectorType::Constant(n, 1);
        
        solver.Solve(op,target,x);
        
        
        std::ostream & os = std::cout
        ExportVarArrow((input-x).norm())
        ExportArrayArrow(N::StdFromEigen(input))
        ExportArrayArrow(N::StdFromEigen(x))
        ExportArrayArrow(N::StdFromEigen(target))
        << "\n";
    }
    
}

int main(int argc, const char * argv[]) {

    Minkowski_2_Test::Test2(); 

}