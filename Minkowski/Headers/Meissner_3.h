//
//  Meissner_3.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

/*
 Numerical solution of the three dimensional Meissner problem. See Meissner_2 for comments.
 */

#ifndef Minkowski_Meissner_3_h
#define Minkowski_Meissner_3_h

#include "Minkowski_3.h"
#include "ConstraintsProduct.h"

namespace Meissner_3 {
    namespace N = QuotientedNewton;
    namespace G = Minkowski_3::G;
    const int Dimension = 3;
    typedef N::ScalarType ScalarType;
    
    
    // NloptInterface. Pass Minkowski_3::Functional as data, with pts(n/2+i) = -pts(i)
    ScalarType EvaluateVolume(const std::vector<ScalarType> &, std::vector<ScalarType> &, void*);
    ScalarType EvaluateConstraint(const std::vector<ScalarType> &, std::vector<ScalarType> &, void*);
    
    #include "Meissner_3.hxx"
}



#endif
