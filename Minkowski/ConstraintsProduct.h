//
//  ConstraintsProduct.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_ConstraintsProduct_h
#define Minkowski_ConstraintsProduct_h

#include "QuotientedNewton.h"
#include <functional>

namespace QuotientedNewton {
    
    ScalarType ConstraintsProduct(const VectorType & x, const std::vector<Triplet> & jac,
                                  VectorType & grad){
        const int n=(int)x.size();
        grad = VectorType::Constant(n,0);
        if(x.minCoeff()<=0) return 0;
        
        // Taylor development for the sum of logarithms
        for(const Triplet & t : jac)
            grad[t.col()] += t.value()/x[t.row()];
//        ScalarType logSum = x.unaryExpr(std::ptr_fun(log)).sum(); // Removed C++17
		ScalarType logSum = x.unaryExpr([](ScalarType s){return std::log(s);}).sum();

        
        // Normalization
        grad/=n;
        logSum/=n;
        
        // Exponentiation
        const ScalarType result = exp(logSum);
        grad*=result;
        
        
/*        std::ostream & os = std::cout
        ExportVarArrow(result)
        << "\n";*/
        
        return result;
    }
    
    
    
}

#endif
