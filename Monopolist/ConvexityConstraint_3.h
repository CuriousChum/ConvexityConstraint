//
//  ConvexityConstraint_3.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_ConvexityConstraint_3_h
#define CGalTest_ConvexityConstraint_3_h

#include "Geometry_3.h"

namespace Geometry_3 {
    
    struct ConvexityConstraint : ConstraintType {
        
        // Put weights (! not heights) at zero.
        ConvexityConstraint(const std::vector<CGT::Full_point> &, ScalarType = 1e-8);
        
        virtual void SetValues(const std::vector<ScalarType> &);
        // To do : update regular triangulation with new heights by edge flipping ?
        
        virtual void ComputeValJacHess(FlagType);
        
        CGT::RT rt;
        virtual void PrintSelf(std::ostream & os) const;
        virtual std::string Name() const {return "ConvexityConstraint_3";}
    protected:
        std::vector<CGT::Full_point> pts;
        const ScalarType degen;
        std::vector<TensorCoef> hessian;
    };
    
#include "ConvexityConstraint_3.hxx"
}

#endif