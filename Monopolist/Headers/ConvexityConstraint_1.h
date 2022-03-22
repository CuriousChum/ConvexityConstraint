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
        
        virtual void SetValues(const std::vector<ScalarType> &);        
        virtual void ComputeValJacHess(FlagType);
        
        virtual void PrintSelf(std::ostream & os) const;
        virtual std::string Name() const {return "ConvexityConstraint_1";}
    protected:
        std::vector<ScalarType> IDpts;//Inverse difference between successive point positions
		std::vector<ScalarType> x;
    };
    
//#include "ConvexityConstraint_1.hxx"


ConvexityConstraint::ConvexityConstraint(const std::vector<ScalarType> & pts){
	assert(!pts.empty());
	IDpts.resize(pts.size()-1);
	for(int i=0; i<IDpts.size(); ++i){
		const ScalarType Dpt = pts[i+1]-pts[i];
		if(Dpt<=0) throw "Convexity_constraint 1D error : points must be increasing";
		IDpts[i] = 1/Dpt;
	}
}

void ConvexityConstraint::SetValues(const std::vector<ScalarType> & x_){
	assert(x_.size() == IDpts.size()+1);
	x = x_;
	// TODO ? report error here, or at computation ?
}

void ConvexityConstraint::ComputeValJacHess(FlagType r){
	error=false;
	if(r & RValues){
		values.resize(x.size());
		values[0] = Infinity;
		values[x.size()-1]=Infinity;
	}
	for(int i=1; i<x.size()-1; ++i){
		const ScalarType Il =IDpts[i-1], Ir = IDpts[i], // Inverse length, left and right
		d2x = (x[i+1]-x[i])*Ir + (x[i-1]-x[i])*Il;
		error = error || (d2x <= 0); 
		if(r & RValues){values[i] = d2x;}
		if(r & RJacobian){
			jacobian.push_back({i,i+1, Ir});
			jacobian.push_back({i,i,  -Ir-Il});
			jacobian.push_back({i,i-1, Il});
		}
		// The hessian remains empty
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
