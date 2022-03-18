//
//  ConvexityConstraint_Test.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_ConvexityConstraint_Test_hxx
#define CGalTest_ConvexityConstraint_Test_hxx

#include "ConvexityConstraint.h"
#include <fstream>


namespace ConvexityConstraint_Test {
    using namespace Geometry_2;
    
    void Pos0(){
        int n=10;
        int M=100;
        std::vector<ScalarType> val;
        val.reserve(n);
        
        for(int i=0; i<n; ++i)
            val.push_back( (1+(std::rand()%M))/ScalarType(M));
        
        PositivityConstraint c;
        c.SetValues(val);
        
        typedef PositivityConstraint::Request PC;
        c.Compute(PC(PC::RLogSum | PC::RLogGrad | PC::RLogHessian
                           | PC::RValues | PC::RJacobian));
        
        std::ofstream os;
        os.open("Pos0.txt");
        os << "{"
        ExportVarArrow(c)
        << "}";
    }
    
    void Cvx0(){
        std::vector<PointType> pts;
        std::vector<ScalarType> values;
        int n=10;
        for(int i=0; i<n; ++i){
            for(int j=0; j<n; ++j){
                pts.push_back(PointType{ScalarType(i),ScalarType(j)});
                values.push_back(i*i+j*j);
            }
        }
        
        ConvexityConstraint c(pts);
        c.SetValues(values);
        
        typedef ConstraintType::Request PC;
        c.Compute(PC(PC::RLogSum | PC::RLogGrad | PC::RLogHessian
                     | PC::RValues | PC::RJacobian));
        
        std::ofstream os;
        os.open("Cvx0.txt");
        
        os << "{"
        ExportVarArrow(c)
        << "}";

    }
    
}

#endif
