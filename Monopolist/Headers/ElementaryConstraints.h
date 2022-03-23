#pragma once
//
//  ElementaryConstraints.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//

#include "Constraint.h"

namespace Constraint {

/**
 Standard log barrier constraint for positivity.
 */
struct PositivityConstraint : ConstraintType {
	virtual void SetValues(const std::vector<ScalarType> &) override;
	virtual void ComputeValJacHess(FlagType) override;
	virtual std::string Name() const override {return "PositivityConstraint";}
protected:
	std::vector<ScalarType> val;
};

/**
 Defines a (log barrier for a) constraint applied to a subset of indices.
 F( x_{i_k}; 0<=k<K ) > 0,
 where the constraint F and the indices (i_k) are given.
 
 TODO : It might be better to subclass the constraint rather than have it as a field.
 But for that, less functions should be overridable.
 */
struct ResampledConstraint : ConstraintType {
	std::unique_ptr<ConstraintType> constraint;
	std::vector<IndexType> indices;
	ResampledConstraint(std::unique_ptr<ConstraintType>, const std::vector<IndexType> &);
	virtual void SetValues(const std::vector<ScalarType> &) override;
	virtual void ComputeValJacHess(FlagType) override;
	virtual std::string Name() const override {return "Resampled_"+constraint->Name();}
	virtual void PrintSelf(std::ostream & os) const override;
protected:
	std::vector<ScalarType> x;
};

#include "ElementaryConstraints.hxx"
}
