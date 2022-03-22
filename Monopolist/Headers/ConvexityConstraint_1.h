#pragma once
//
//  ConvexityConstraint_1.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//

#include "Constraint.h"

namespace Geometry_1 {
	using namespace Constraint;
    
	/**
	 One dimensional convexity constraint
	 */
    struct ConvexityConstraint : ConstraintType {
        
		// Input the base points, ordered in an increasing manner
        ConvexityConstraint(const std::vector<ScalarType> &);
        
        virtual void SetValues(const std::vector<ScalarType> &) override;
        virtual void ComputeValJacHess(FlagType) override;
        
        virtual void PrintSelf(std::ostream & os) const override;
        virtual std::string Name() const override {return "ConvexityConstraint_1";}
    protected:
        std::vector<ScalarType> IDpts;//Inverse difference between successive point positions
		std::vector<ScalarType> x;
    };
    
//#include "ConvexityConstraint_1.hxx"


ConvexityConstraint::ConvexityConstraint(const std::vector<ScalarType> & pts){
	assert(!pts.empty());
	numberOfConstraints = (IndexType)pts.size()-2;
	numberOfUnknowns = (IndexType)pts.size();
	IDpts.resize(pts.size()-1);
	for(int i=0, error=0; i<IDpts.size(); ++i){
		const ScalarType Dpt = pts[i+1]-pts[i];
		if(Dpt<=0) throw NS::DomainError("ConvexityConstraint 1D : non-increasing points"); 
		error += (Dpt<=0);
		IDpts[i] = 1/Dpt;
	}
}

void ConvexityConstraint::SetValues(const std::vector<ScalarType> & x_){
	assert(x_.size() == numberOfConstraints+2);
	x = x_;
}

void ConvexityConstraint::ComputeValJacHess(FlagType r){
	if(r & RValues) values.resize(numberOfConstraints);
	if(r & RJacobian) assert(jacobian.empty()); 

	for(int i=1; i<x.size()-1; ++i){
		const ScalarType Il =IDpts[i-1], Ir = IDpts[i], // Inverse length, left and right
		d2x = (x[i+1]-x[i])*Ir + (x[i-1]-x[i])*Il;
		if(d2x<=0) throw NS::DomainError("ConvexityConstraint 1D : non-convex data");
		const int iConstraint = i-1;
		if(r & RValues){values[iConstraint] = d2x;}
		if(r & RJacobian){
			jacobian.push_back({iConstraint,i+1, Ir});
			jacobian.push_back({iConstraint,i,  -Ir-Il});
			jacobian.push_back({iConstraint,i-1, Il});
		}
		// The constraints are linear, hence the hessian remains empty
	}
}

void ConvexityConstraint::PrintSelf(std::ostream & os) const {
	os << "{"
	ExportArrayArrow(IDpts)
	ExportArrayArrow(x)
	<< "\"superclass\"->";
	ConstraintType::PrintSelf(os);
	os << "}";
}

}
