#pragma once
//
//  ElementaryConstraints.hxx
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//


// %%%%%%%%%%%%%%%%%%%%%%%%%%%% Positivity %%%%%%%%%%%%%%%%%%%%%%%%%

void PositivityConstraint::SetValues(const std::vector<ScalarType> & val_){
	numberOfUnknowns = (IndexType)val_.size();
    val = val_;
	for(ScalarType x : val) {
		if(x<=0) throw NS::DomainError("PositivityConstraint : non-positive data");}
};


void PositivityConstraint::ComputeValJacHess(FlagType r){
    numberOfConstraints=(IndexType)val.size();
    if(r & RValues) values.reserve(numberOfConstraints);
    if(r & RJacobian) jacobian.reserve(numberOfConstraints);
    for(int i=0 ; i<numberOfConstraints; ++i) {
        if(r & RValues) values.push_back(val[i]);
        if(r & RJacobian) jacobian.push_back({i,i,1});
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%% Bounds %%%%%%%%%%%%%%%%%%%%%%%%%

void BoundConstraint::Init(){
	if(lb.size() != ub.size())
		throw NS::DomainError("Bound constraint : different number of bounds");
	numberOfUnknowns = (IndexType)lb.size();
	numberOfConstraints=0;
	for(const ScalarType & x:lb) numberOfConstraints += (x > -Infinity);
	for(const ScalarType & x:ub) numberOfConstraints += (x <  Infinity);
}

void BoundConstraint::SetValues(const std::vector<ScalarType> & val_){
    val = val_;
	for(IndexType i=0; i<numberOfUnknowns; ++i){
		if(!(lb[i]<val[i] && val[i]<ub[i])){
			std::cout
			ExportVarArrow(i)
			ExportVarArrow(val[i])
			ExportVarArrow(lb[i])
			ExportVarArrow(ub[i]) << std::endl;
			throw NS::DomainError("BoundConstraint : data must be strictly within bounds");}
	}
};


void BoundConstraint::ComputeValJacHess(FlagType r){
    if(r & RValues) values.reserve(numberOfConstraints);
    if(r & RJacobian) jacobian.reserve(numberOfConstraints);
    for(int i=0 ; i<numberOfUnknowns; ++i) {
		if(lb[i] > -Infinity){
			if(r & RValues) values.push_back(val[i]-lb[i]);
			if(r & RJacobian) jacobian.push_back({i,i,1});
		}
		if(ub[i] < Infinity){
			if(r & RValues) values.push_back(ub[i]-val[i]);
			if(r & RJacobian) jacobian.push_back({i,i,-1});
		}
    }
}

// %%%%%%%%%%%%%%%%%%%%%% Resampled constraint %%%%%%%%%%%%%%%%%%%

ResampledConstraint::
ResampledConstraint(std::unique_ptr<ConstraintType> constraint_,
					const std::vector<IndexType> & indices_):
					constraint(std::move(constraint_)),indices(indices_){
						x.resize(indices.size());}


void ResampledConstraint::SetValues(const std::vector<ScalarType> & x_){
	numberOfUnknowns = (IndexType)x_.size();
	assert(x.size()==indices.size());
	for(int i=0; i<indices.size(); ++i) {
		assert(0<=indices[i] && indices[i]<x_.size());
		x[i] = x_[indices[i]];}
	constraint->SetValues(x);
}

void ResampledConstraint::ComputeValJacHess(FlagType r_){
	FlagType r = r_ & ~(RLogSum | RLogGrad | RLogHessian);
	constraint->Compute(r);
	numberOfConstraints = constraint->numberOfConstraints;
	if(r & RValues){values = constraint->values;}
	if(r & RJacobian){
		for(const MatCoef & c : constraint->jacobian){
			assert(0<=c.j && c.j<indices.size());
			jacobian.push_back({c.i,indices[c.j],c.v});}
	}
	if(r & RHessian){
		for(const TensorCoef & c : constraint->hessian){
			assert(0<=c.j && c.j<indices.size());
			assert(0<=c.k && c.k<indices.size());
			hessian.push_back({c.i,indices[c.j],indices[c.k],c.v});}
	}
}

void ResampledConstraint::PrintSelf(std::ostream & os) const {
	os << "{"
	ExportVarArrow(*constraint)
	ExportArrayArrow(indices)
	<< "}";
}
