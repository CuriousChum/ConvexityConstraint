//
//  Meissner_2.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

/*
 Meissner's problem asks for the convex body whose of a given width,
 and whose volume is minimal. In this file, we numerically solve the two dimensional
 instantiation of this problem, whose solution is known as Reulaux's triangle.
*/


#ifndef Minkowski_Meissner_2_h
#define Minkowski_Meissner_2_h

#include "ConstraintsProduct.h"
#include "Minkowski_2.h"

namespace Meissner_2 {
    namespace N = QuotientedNewton;
    namespace G = Minkowski_2::G;
    const int Dimension = 2;
    typedef N::ScalarType ScalarType;
    
	// See implementation file for comments
    static ScalarType EvaluateVolume(const std::vector<ScalarType> &,
									 std::vector<ScalarType> &, void*);
    static ScalarType EvaluateConstraint(const std::vector<ScalarType> &,
										 std::vector<ScalarType> &, void*);
#include "Meissner_2.hxx"
}

#endif
