//
//  PrincipalAgent_3.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_3_h
#define CGalTest_PrincipalAgent_3_h

#include "Geometry_3.h"

namespace Geometry_3 {
    
    struct PrincipalAgent : NS::Functionnal {
        
        PrincipalAgent(const CGT::RT &, ScalarType degen=1e-8);
        
        // Glue code : Barrier for the constraint
        virtual ScalarType Value(const NS::VectorType &) override;
        virtual const NS::VectorType & Gradient() override {return grad_;} // At latest position.
        virtual const NS::SparseMatrixType & Hessian() override {return hess_;}
        
        
        void PrintSelf(std::ostream &) const override;
        friend std::ostream & operator << (std::ostream & os, const PrincipalAgent & a){
            a.PrintSelf(os); return os;}
        
    protected:
        NS::VectorType grad0_, grad_; // gradient at origin, and at current position
        NS::SparseMatrixType hess_;
        
        
    };

    
#include "PrincipalAgent_3.hxx"
    
}


#endif
