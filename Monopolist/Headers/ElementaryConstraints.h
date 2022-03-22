#pragma once
//
//  ElementaryConstraints.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//

namespace Constraint {

/**
 Standard log barrier constraint for positivity.
 */
struct PositivityConstraint : ConstraintType {
	virtual void SetValues(const std::vector<ScalarType> &);
	virtual void ComputeValJacHess(FlagType);
	virtual std::string Name() const {return "PositivityConstraint";}
protected:
	std::vector<ScalarType> val;
};

/**
 Defines a (log barrier for a) constraint applied to a subset of indices.
 */
struct ResampledConstraint : ConstraintType {
	std::unique_ptr<ConstraintType> constraint;
	std::vector<IndexType> indices;
	virtual void SetValues(const std::vector<ScalarType> &);
	virtual void ComputeValJacHess(FlagType);
	virtual std::string Name() const {return "Resampled_"+constraint->Name();}
protected:
	std::vector<ScalarType> x;
};

#include "ElementaryConstraints.hxx"
}
