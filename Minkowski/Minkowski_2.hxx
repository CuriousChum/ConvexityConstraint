//
//  Minkowski_2.hxx
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 02/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Minkowski_2_hxx
#define Minkowski_Minkowski_2_hxx


int Functional::Compute(N::VectorType & input, N::VectorType & output){
    ErrCode error;
    error=CheckInput(input); if(error!=ErrCode::NoError) return (int)error;
    Compute(input);
    output = facetLengths;
    Normalize(input);
    return (int)error;
}

Functional::ErrCode Functional::CheckInput(const N::VectorType & input) const {
    if(input.size()!=Size()){assert(false); return ErrCode::Size;}
    // Check strict convexity
    const int n=Size();
    for(int si=0; si<n; ++si){
        const int ri=(si+n-1)%n, ti=(si+1)%n;
        const G::Vector r = pts.col(ri)/input[ri], s=pts.col(si)/input[si], t=pts.col(ti)/input[ti];
        if( G::Orientation(r,s,t)<=0 ) return ErrCode::Orientation;
        if(input[si]<=0) return ErrCode::Sign;
    }
    return ErrCode::NoError;
}


void Functional::Compute(const N::VectorType & input){
    const int n=Size();
    area=0;
    jacobian.clear();
    facetLengths.resize(n);
    facetMoments.resize(2,n);
    moment=G::Vector{0,0};
    
    for(int si=0; si<n; ++si){
        const int ri=(si+n-1)%n, ti=(si+1)%n;
        
        const ScalarType ra=input[ri], sa=input[si], ta=input[ti];
        const G::Vector r = pts.col(ri), s=pts.col(si), t=pts.col(ti);
        
        const ScalarType Srs = G::ScalarProduct(r,s),   Sst = G::ScalarProduct(s,t);
        const ScalarType Drs = G::Determinant(r,s),     Dst = G::Determinant(s,t);
        
        const ScalarType Lst = (ta-sa*Sst)/Dst, Lrs = (ra-sa*Srs)/Drs;
        const ScalarType facetLength = Lst+Lrs;
        const G::Vector facetMoment = facetLength*(sa*s + ((Lst-Lrs)/2)*G::Perp(s));
        
        jacobian.push_back({si,ti, 1/Dst});
        jacobian.push_back({si,si, -Sst/Dst -Srs/Drs});
        jacobian.push_back({si,ri, 1/Drs});
        
        facetLengths[si]=facetLength;
        facetMoments.col(si)=facetMoment;
    }
    area = facetLengths.dot(input)/2.;
    moment = (facetMoments*input)/3.;
}
void Functional::Normalize(N::VectorType & input) {
    input -= pts.transpose() * Barycenter();
}

int Functional::DescentDirection(N::VectorType & diff) {
    if(referenceIndices[0]==N::BadIndex) SetReferenceIndices();
    for(int i: referenceIndices)
        jacobian.push_back({i,i,1});
    N::SparseMatrixType m;
    m.resize(Size(), Size());
    m.setFromTriplets(jacobian.begin(), jacobian.end());
    
    Eigen::SimplicialLDLT<N::SparseMatrixType> solver(m);
    diff = solver.solve(diff);
    N::VectorType & dir = diff;
    dir=-dir;
    
    // Correction
    const G::Vector dirMoment = facetMoments*dir;
    const G::Vector correction = (facetMoments*pts.transpose()).inverse()*dirMoment;
    dir -= pts.transpose()*correction;
    return (int)ErrCode::NoError;
}

void Functional::SetReferenceIndices(){
    const int maxIt = 100;
    int i;
    for(i=0; i<maxIt; ++i){
        for(int j=0; j<Dimension; ++j)
            referenceIndices[j] = std::rand()%Size();
        if(ReferenceIndicesConditionNumber()>referenceIndicesConditionNumberThreshold)
            break;
    }
    if(i==maxIt){ assert(false);
        throw "Could not set reference indices, check point set for coplanarity.";
    }
}

ScalarType Functional::ReferenceIndicesConditionNumber(){
    Eigen::Matrix<ScalarType,Dimension,Dimension> m;
    for(int i=0; i<Dimension; ++i)
        m.col(i) = pts.col(referenceIndices[i]);
    return fabs(m.determinant())/ std::pow(m.norm(),Dimension);
}


#endif
