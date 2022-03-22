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
    error = std::accumulate(val.begin(), val.end(), 0, [](int k, ScalarType x) {return k+int(x<=0);});
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
	error = constraint->error;
}

void ResampledConstraint::ComputeValJacHess(FlagType r_){
	FlagType r = r_ & ~(RLogSum | RLogGrad | RLogHessian);
	constraint->Compute(r);
	error = constraint->error;
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
