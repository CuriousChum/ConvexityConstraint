//
//  MinkowskiConstraint_2.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 20/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_MinkowskiConstraint_2_hxx
#define CGalTest_MinkowskiConstraint_2_hxx


ConvexityConstraint::ConvexityConstraint(const std::vector<Vector> & pts_):
pts(pts_){numberOfConstraints=(int)pts_.size();}

void ConvexityConstraint::SetValues(const std::vector<ScalarType> & input_){
    assert(input_.size()==pts.size());
    input=input_;
    const int n=numberOfConstraints;
    
    // Check strict convexity
    error=0;
    for(int si=0; si<n; ++si){
        const int ri=(si+n-1)%n, ti=(si+1)%n;
        const Vector r = pts[ri]/input[ri], s=pts[si]/input[si], t=pts[ti]/input[ti];
        if( Orientation(r,s,t)<=0 ) ++error;
        if(input[si]<=0) ++error;
    }
}

void ConvexityConstraint::ComputeValJacHess(FlagType r){
    const int n=numberOfConstraints;
    values.resize(n);
    area=0;
    
    for(int si=0; si<n; ++si){
        const int ri=(si+n-1)%n, ti=(si+1)%n;
        
        const ScalarType ra=input[ri], sa=input[si], ta=input[ti];
        const Vector r = pts[ri], s=pts[si], t=pts[ti];
        
        const ScalarType Srs = ScalarProduct(r,s), Sst = ScalarProduct(s,t);
        const ScalarType Drs = Determinant(r,s), Dst = Determinant(s,t);
        
        values[si] = (ta-sa*Sst)/Dst + (ra-sa*Srs)/Drs;
        
        /*
         std::cout
         ExportVarArrow(Dst)
         ExportVarArrow(Sst)
         ExportVarArrow(Drs)
         ExportVarArrow(values[si])
         ExportVarArrow(ra)
         ExportVarArrow(sa)
         ExportVarArrow(ta)
         << "\n";*/
        
        jacobian.push_back({si,ti, 1/Dst});
        jacobian.push_back({si,si, -Sst/Dst -Srs/Drs});
        jacobian.push_back({si,ri, 1/Drs});
        
        area+=values[si]*sa;
    }
    area/=2;
}

void ConvexityConstraint::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayArrow(input)
    ExportVarArrow(area)
    << '"' << "superclass" << '"' << "->";
    ConstraintType::PrintSelf(os);
    os << "}";
}


#endif
