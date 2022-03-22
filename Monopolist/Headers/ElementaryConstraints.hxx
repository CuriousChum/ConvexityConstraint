#pragma once
//
//  ElementaryConstraints.hxx
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//


// %%%%%%%%%%%%%%%%%%%%%%%%%%%% Positivity %%%%%%%%%%%%%%%%%%%%%%%%%

void PositivityConstraint::SetValues(const std::vector<ScalarType> & val_){
    val = val_;
    error = std::accumulate(val.begin(), val.end(), 0, [](int k, ScalarType x) {return k+int(x<=0);});
};


void PositivityConstraint::ComputeValJacHess(FlagType r){
    numberOfConstraints=(int)val.size();
    if(r & RValues) values.reserve(numberOfConstraints);
    if(r & RJacobian) jacobian.reserve(numberOfConstraints);
    for(int i=0 ; i<numberOfConstraints; ++i) {
        if(r & RValues) values.push_back(val[i]);
        if(r & RJacobian) jacobian.push_back({i,i,1});
    }
}

// %%%%%%%%%%%%%%%%%%%%%% Resampled constraint %%%%%%%%%%%%%%%%%%%



void SetValues(const std::vector<ScalarType> & x_){
	if(x.empty()){
		IndexType imin,imax;
		
	}
}
