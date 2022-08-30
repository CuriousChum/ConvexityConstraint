#pragma once
//
//  ConvexityConstraint_3.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include "Geometry_3.h"
namespace Geometry_3 {
    
struct ConvexityConstraint : ConstraintType {
	using Point = CGT::Point;
	ConvexityConstraint(const std::vector<Point> &,
						const std::vector<ScalarType> &,
						const std::vector<bool> &,
						FlagType, ScalarType degen = 1e-8);
	
	CGT::RT rt;
	virtual void PrintSelf(std::ostream & os) const override;
	virtual std::string Name() const override {return "ConvexityConstraint_3";}
protected:
	std::vector<ScalarType> heights;
	ScalarType Height(IndexType index) const {assert(0<=index && index<heights.size());
		return heights[index];}
	const ScalarType degen;
	virtual void ComputeValJacHess(FlagType) override;
};

#include "ConvexityConstraint_3.hxx"
}
