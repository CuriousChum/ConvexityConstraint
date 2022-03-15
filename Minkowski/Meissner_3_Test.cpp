//
//  Meissner_3_Test.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <fstream>
#include "Headers/Meissner_3.h"
#include "nlopt.hpp"

namespace Meissner_3_Test {
    using namespace Meissner_3;
    typedef N::VectorType NV;
    
    
    void Test1(int subdivision, ScalarType stopval = 0.4197){
        const ScalarType delta=0.02;
        std::ifstream is;
        is.open("HalfPts_"+std::to_string(subdivision)+".txt");
        
        int n,checkDim;
        is >> n;
        is >> checkDim;
        assert(checkDim==Dimension);
        
        
        Minkowski_3::Functional mink;
        G::VectorList & pts = mink.pts;
        pts.resize(Dimension,2*n);
        for(int i=0; i<n; ++i)
            for(int j=0; j<Dimension; ++j)
                is >> pts(j,i);
            
        
        for(int i=0; i<n; ++i)
            pts.col(n+i) = -pts.col(i);
        
        {
            std::ostream & os = std::cout
            ExportVarArrow(n)
            ExportArrayArrow(N::StdFromEigen(NV(pts.col(0))))
            << "\n";
        }
        
        Eigen::Matrix<ScalarType,4,3> tetra;
        tetra <<
        0.,0.,-1.,
        -0.942809,0.,0.333333,
        0.471405,-0.816497,0.333333,
        0.471405,0.816497,0.333333;
        
        
        std::vector<ScalarType> x(n,0.5),init;
        for(int i=0; i<n; ++i){
            for(int j=0; j<4; ++j){
                ScalarType s = pts.col(i).dot(tetra.row(j));
                x[i] += delta*(s*s*s);
            }
        }
        init=x;
        
        
        {
//            for(int i=0; i<2*n; ++i)
//                std::cout << pts.col(i) << "\n\n";
            
            NV input(2*n), output(2*n);
            for(int i=0; i<n; ++i){
                input[i]=x[i];
                input[n+i]=1-x[i];
                }
            const int result = mink.Compute(input,output);
            std::ostream & os = std::cout
            ExportVarArrow(result)
            ExportArrayArrow(N::StdFromEigen(input))
            ExportArrayArrow(N::StdFromEigen(output))
            << "\n";
        }
        
        
        nlopt::opt algo=nlopt::opt(nlopt::LD_SLSQP, n); //LD_MMA
        algo.set_min_objective(EvaluateVolume, static_cast<void*>(&mink));
        algo.add_inequality_constraint(EvaluateConstraint, static_cast<void*>(&mink));
        algo.set_lower_bounds(0);
        algo.set_upper_bounds(1);
        algo.set_maxeval(10000);
        
        algo.set_stopval(stopval);
        
        ScalarType minimum;
        std::cout << "Opt begin\n";
        nlopt::result result = algo.optimize(x,minimum);
        std::cout << "Opt end\n";
        
        std::ofstream os; os.open("Meissner_3_"+std::to_string(subdivision)+".txt");
        os << "{"
        ExportVarArrow(minimum)
        ExportVarArrow(result)
        ExportArrayArrow(x)
        ExportArrayArrow(init)
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(0))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(1))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(2))))
        << "}";

    }
}

int main(int argc, const char * argv[]) {
    /*
    if(argc>1){
        const int i= atoi(argv[1]);
        std::cout << i << std::endl;
        
        double stopVal = 0.4197;
        if(argc>2) stopVal = atof(argv[2]);
        
        if(2<=i && i<=10)
            Meissner_3_Test::Test1(i,stopVal);
    }
 */
    
    for(int i=5; i<10; ++i)
        Meissner_3_Test::Test1(i);

    return 0;
}
