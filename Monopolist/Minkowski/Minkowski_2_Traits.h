//
//  Minkowski_2_Traits.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 20/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Minkowski_2_Traits_h
#define CGalTest_Minkowski_2_Traits_h

#include <array>
#include "Constraint.h"

namespace Minkowski_2 {
    using namespace Constraint;
    const int Dimension=2;
    
    typedef std::array<ScalarType,2> Vector;
    
    
    Vector operator / (const Vector & v,ScalarType a){return Vector{v[0]/a,v[1]/a};}
    ScalarType ScalarProduct(const Vector & u, const Vector & v){return u[0]*v[0]+u[1]*v[1];}
    ScalarType Determinant(const Vector & u, const Vector & v){return u[0]*v[1]-u[1]*v[0];}
    ScalarType Orientation(const Vector & u, const Vector & v, const Vector & w){
        //Det (v-u,w-u)
        return (v[0]-u[0])*(w[1]-u[1])-(v[1]-u[1])*(w[0]-u[0]);
    }

    
    
}


#endif
