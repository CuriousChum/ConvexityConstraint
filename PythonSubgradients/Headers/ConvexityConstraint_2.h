#pragma once
//
//  ConvexityConstraint_2.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015. Some from Q. Merigot.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <vector>

#include "Geometry_2.h"
#include "ConvexityConstraint_1.h"

namespace Geometry_2 {

/**
 This functional evaluates the area of the subgradient cells associated with all the interior discretization points.
 Taking their logarithms, we obtain a barrier function for convexity in the domain interior.
 
 The boundary points are ignored however.
 */
struct ConvexityConstraint : ConstraintType {
	using Point = CGT::Point;
	ConvexityConstraint(const std::vector<Point> &,
						const std::vector<ScalarType> &,
						const std::vector<bool> &, FlagType);

    CGT::RT rt;
    virtual void PrintSelf(std::ostream & os) const override;
    virtual std::string Name() const override {return "ConvexityConstraint_2";}
protected:
	std::vector<ScalarType> heights;
	ScalarType Height(IndexType index) const {assert(0<=index && index<heights.size());
		return heights[index];}
	virtual void ComputeValJacHess(FlagType) override;
};
    
#include "ConvexityConstraint_2.hxx"

}
